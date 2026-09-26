#!/usr/bin/env python3
"""Generate a 2D chtMultiRegion case for a stagnant molten-salt fuel slot between graphite walls.

Geometry (x across the slot, z vertical, 1 cell in y, 'empty' front/back):

    x = 0          w              w+d            2w+d
    | graphite wall | fuel salt (d) | graphite wall |
  coolant BC                                    coolant BC

Both faces of the fuel slot are cooled through a 0.5 cm graphite wall.  This matches the
repo's plate stacking (F | wall | C | wall | F ...): every fuel slot has a coolant row on each side.
Top and bottom of the fuel slot are closed, adiabatic, no-slip walls (zero net flow).

Usage: make_case.py <outdir> [--d_f 0.005] [--q_avg 1e7] [--gravity 1] [--coolant linear|uniform]
                    [--h 0 (0 = fixed wall temperature)] [--nxf 40] [--nz 815] [--solver steady|transient]
"""
import argparse, math, os, shutil, stat, textwrap
import numpy as np

HDR = """FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    location    "{loc}";
    object      {obj};
}}
"""

# ---------------------------------------------------------------- properties (SI)
# Fuel salt 7LiF-BeF2-ZrF4-UF4 61.83-29.17-5.0-4.0 mol%
#  density: MSRE fuel correlation rho = 2.575 - 5.13e-4*T[C] g/cc (ORNL-4658 Table 8.2, Cantor method),
#  rescaled x1.15648 to 2.593 g/cc at 922 K with Cantor's additive molar volumes (ORNL-TM-4308) as in
#  ../graphite_slab_core.  -> rho(T[K]) = 3140.0 - 0.59328*T  kg/m3, beta = 2.29e-4 1/K at 922 K.
RHO_A, RHO_B = 3140.0009568, -0.5932765258
#  viscosity: mu = 0.116 exp(3755/T) cP (MSRE fuel, Cantor ORNL-TM-4308 / ORNL-4658); 5th-order polynomial
#  fit over 800-1250 K (max error 0.35 %).
MU_COEFFS = [1.21431605e+00, -5.16254204e-03, 8.96105594e-06, -7.87497715e-09, 3.48874372e-12, -6.21636273e-16]
K_FUEL = 1.05      # W/m/K
CP_FUEL = 1970.0   # J/kg/K (0.47 cal/g/K, MSRE fuel)
# Graphite
K_GR, RHO_GR, CP_GR = 30.0, 1870.0, 1700.0

H = 1.63            # m, fuel slot height
W_GR = 0.005        # m, graphite shared wall
PEAK_AVG = 1.4      # axial peak/average of the cosine power shape
T_COOL_BOT, T_COOL_TOP = 550.0 + 273.15, 650.0 + 273.15


def extrapolated_height(h=H, ratio=PEAK_AVG):
    """H_e such that cos(pi*(z-H/2)/H_e) over [0,H] has peak/average = ratio."""
    lo, hi = 1e-6, math.pi / 2          # x = pi*H/(2*H_e); average/peak = sin(x)/x
    for _ in range(200):
        x = 0.5 * (lo + hi)
        if math.sin(x) / x > 1 / ratio:
            lo = x
        else:
            hi = x
    return math.pi * h / (2 * x)


def band_sources(q_avg, nb):
    He = extrapolated_height()
    qpk = PEAK_AVG * q_avg
    zb = np.linspace(0, H, nb + 1)
    s = np.sin(math.pi * (zb - H / 2) / He)
    return zb, qpk * He / math.pi * np.diff(s) / np.diff(zb), He


def coolant_T(z, mode):
    if mode == "uniform":
        return np.full_like(z, 600.0 + 273.15)
    return T_COOL_BOT + (T_COOL_TOP - T_COOL_BOT) * z / H


def w(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(text)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--d_f", type=float, default=0.005)
    ap.add_argument("--q_avg", type=float, default=1e7)
    ap.add_argument("--gravity", type=int, default=1)
    ap.add_argument("--coolant", default="linear")
    ap.add_argument("--h", type=float, default=0.0)
    ap.add_argument("--nxf", type=int, default=40)
    ap.add_argument("--nxg", type=int, default=8)
    ap.add_argument("--nz", type=int, default=815)
    ap.add_argument("--nbands", type=int, default=163)
    ap.add_argument("--solver", default="transient")
    ap.add_argument("--endTime", type=float, default=6000)
    ap.add_argument("--iters", type=int, default=20000)
    a = ap.parse_args()

    c = a.outdir
    if os.path.exists(c):
        shutil.rmtree(c)
    os.makedirs(c)
    d = a.d_f
    xs = [0.0, W_GR, W_GR + d, 2 * W_GR + d]
    Y = 1e-3
    steady = a.solver == "steady"
    app = "chtMultiRegionSimpleFoam" if steady else "chtMultiRegionFoam"

    # ------------------------------------------------------------ case parameters (read by scripts)
    zb, qb, He = band_sources(a.q_avg, a.nbands)
    w(f"{c}/caseParams.json", __import__("json").dumps(dict(
        d_f=d, q_avg=a.q_avg, q_peak=PEAK_AVG * a.q_avg, H=H, H_extrap=He, w_graphite=W_GR, gravity=a.gravity,
        coolant=a.coolant, h=a.h, nxf=a.nxf, nxg=a.nxg, nz=a.nz, nbands=a.nbands, solver=app,
        k_fuel=K_FUEL, cp_fuel=CP_FUEL, rho_A=RHO_A, rho_B=RHO_B, mu_coeffs=MU_COEFFS, k_graphite=K_GR), indent=1))

    # ------------------------------------------------------------ blockMesh
    V = []
    for k, zz in enumerate([0.0, H]):
        for j, yy in enumerate([0.0, Y]):
            for i, xx in enumerate(xs):
                V.append((xx, yy, zz))
    vid = lambda i, j, k: k * 8 + j * 4 + i
    def hexv(i):
        return " ".join(str(vid(*t)) for t in [(i, 0, 0), (i + 1, 0, 0), (i + 1, 1, 0), (i, 1, 0),
                                                (i, 0, 1), (i + 1, 0, 1), (i + 1, 1, 1), (i, 1, 1)])
    def face(i0, pts):
        return "(" + " ".join(str(vid(*p)) for p in pts) + ")"
    nx = [a.nxg, a.nxf, a.nxg]
    zone = ["graphite", "fuel", "graphite"]
    blocks = "\n".join(f"    hex ({hexv(i)}) {zone[i]} ({nx[i]} 1 {a.nz}) simpleGrading (1 1 1)" for i in range(3))
    bottom = "\n".join("            " + face(0, [(i, 0, 0), (i, 1, 0), (i + 1, 1, 0), (i + 1, 0, 0)]) for i in range(3))
    top = "\n".join("            " + face(0, [(i, 0, 1), (i + 1, 0, 1), (i + 1, 1, 1), (i, 1, 1)]) for i in range(3))
    fb = "\n".join("            " + face(0, [(i, 0, 0), (i + 1, 0, 0), (i + 1, 0, 1), (i, 0, 1)]) + "\n" +
                   "            " + face(0, [(i, 1, 0), (i, 1, 1), (i + 1, 1, 1), (i + 1, 1, 0)]) for i in range(3))
    w(f"{c}/system/blockMeshDict", HDR.format(cls="dictionary", loc="system", obj="blockMeshDict") + f"""
scale 1;
vertices
(
{chr(10).join(f'    ({x:.8g} {y:.8g} {z:.8g})' for x, y, z in V)}
);
blocks
(
{blocks}
);
edges ();
boundary
(
    coolantL {{ type wall; faces ( {face(0, [(0,0,0),(0,0,1),(0,1,1),(0,1,0)])} ); }}
    coolantR {{ type wall; faces ( {face(0, [(3,0,0),(3,1,0),(3,1,1),(3,0,1)])} ); }}
    bottom
    {{
        type wall;
        faces
        (
{bottom}
        );
    }}
    top
    {{
        type wall;
        faces
        (
{top}
        );
    }}
    frontAndBack
    {{
        type empty;
        faces
        (
{fb}
        );
    }}
);
mergePatchPairs ();
""")

    # ------------------------------------------------------------ global system files
    g = "(0 0 -9.81)" if a.gravity else "(0 0 0)"
    w(f"{c}/constant/g", HDR.format(cls="uniformDimensionedVectorField", loc="constant", obj="g") +
      f"dimensions [0 1 -2 0 0 0 0];\nvalue {g};\n")
    w(f"{c}/constant/regionProperties", HDR.format(cls="dictionary", loc="constant", obj="regionProperties") +
      "regions\n(\n    fluid (fuel)\n    solid (graphite)\n);\n")
    fo = """
functions
{
    fuelMinMax
    {
        type            fieldMinMax;
        libs            ("libfieldFunctionObjects.so");
        region          fuel;
        fields          (T U);
        mode            magnitude;
        location        false;
        writeControl    timeStep;
        writeInterval   %d;
    }
    graphiteMinMax
    {
        type            fieldMinMax;
        libs            ("libfieldFunctionObjects.so");
        region          graphite;
        fields          (T);
        location        false;
        writeControl    timeStep;
        writeInterval   %d;
    }
}
""" % ((20, 20) if steady else (10, 10))
    if steady:
        ctrl = f"""application {app};
startFrom latestTime;
startTime 0;
stopAt endTime;
endTime {a.iters};
deltaT 1;
writeControl timeStep;
writeInterval 1000;
purgeWrite 2;
writeFormat ascii;
writePrecision 8;
writeCompression off;
timeFormat general;
timePrecision 8;
runTimeModifiable yes;
"""
    else:
        ctrl = f"""application {app};
startFrom latestTime;
startTime 0;
stopAt endTime;
endTime {a.endTime};
deltaT 0.01;
writeControl adjustableRunTime;
writeInterval 250;
purgeWrite 0;
writeFormat ascii;
writePrecision 8;
writeCompression off;
timeFormat general;
timePrecision 8;
runTimeModifiable yes;
adjustTimeStep yes;
maxCo 2;
maxDi 1000;
maxDeltaT 2;
"""
    # NOTE: no function objects - the Debian OpenFOAM v1912 build aborts in functionObjectList::read()
    # (SHA1 digest IOstream bug); monitoring is done from the solver log (Min/max T) and written fields.
    w(f"{c}/system/controlDict", HDR.format(cls="dictionary", loc="system", obj="controlDict") + ctrl)
    w(f"{c}/system/fvSchemes", HDR.format(cls="dictionary", loc="system", obj="fvSchemes") +
      "ddtSchemes{} gradSchemes{} divSchemes{} laplacianSchemes{} interpolationSchemes{} snGradSchemes{}\n")
    w(f"{c}/system/fvSolution", HDR.format(cls="dictionary", loc="system", obj="fvSolution") +
      ("" if steady else "PIMPLE\n{\n    nOuterCorrectors 1;\n}\n"))

    # ------------------------------------------------------------ region: fuel (fluid)
    ddt = "steadyState" if steady else "Euler"
    w(f"{c}/constant/fuel/thermophysicalProperties",
      HDR.format(cls="dictionary", loc="constant/fuel", obj="thermophysicalProperties") + f"""
// Fuel salt 7LiF-BeF2-ZrF4-UF4 61.83-29.17-5.0-4.0 mol%. rho(T) linear (thermal expansion drives buoyancy,
// full variable-density formulation, no Boussinesq approximation needed); mu(T) polynomial fit of MSRE correlation.
thermoType
{{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       polynomial;
    thermo          hPolynomial;
    equationOfState icoPolynomial;
    specie          specie;
    energy          sensibleEnthalpy;
}}
mixture
{{
    specie          {{ molWeight 40.0; }}
    equationOfState {{ rhoCoeffs<8> ({RHO_A:.10g} {RHO_B:.10g} 0 0 0 0 0 0); }}
    thermodynamics
    {{
        Hf              0;
        Sf              0;
        CpCoeffs<8>     ({CP_FUEL} 0 0 0 0 0 0 0);
    }}
    transport
    {{
        muCoeffs<8>     ({' '.join(f'{m:.10g}' for m in MU_COEFFS)} 0 0);
        kappaCoeffs<8>  ({K_FUEL} 0 0 0 0 0 0 0);
    }}
}}
""")
    w(f"{c}/constant/fuel/turbulenceProperties",
      HDR.format(cls="dictionary", loc="constant/fuel", obj="turbulenceProperties") + "simulationType laminar;\n")
    w(f"{c}/constant/fuel/radiationProperties",
      HDR.format(cls="dictionary", loc="constant/fuel", obj="radiationProperties") + "radiation off;\nradiationModel none;\n")
    # heat source: nbands axial bands, each uniform at the band-average of the cosine shape
    fv = [HDR.format(cls="dictionary", loc="constant/fuel", obj="fvOptions"),
          f"// Volumetric fission heating, axial cosine: q(z) = q_peak cos(pi (z - H/2)/H_e), H_e = {He:.4f} m,\n"
          f"// q_avg = {a.q_avg:.6g} W/m3, q_peak = {PEAK_AVG*a.q_avg:.6g} W/m3, {a.nbands} bands of {H/a.nbands*100:.2f} cm.\n"]
    for i, q in enumerate(qb):
        fv.append(f"""heat{i:03d}
{{
    type            scalarSemiImplicitSource;
    active          yes;
    scalarSemiImplicitSourceCoeffs
    {{
        selectionMode   cellSet;
        cellSet         band{i:03d};
        volumeMode      specific;
        injectionRateSuSp {{ h ({q:.8g} 0); }}
    }}
}}
""")
    w(f"{c}/constant/fuel/fvOptions", "".join(fv))
    ts = [HDR.format(cls="dictionary", loc="system/fuel", obj="topoSetDict"), "actions\n(\n"]
    for i in range(a.nbands):
        ts.append(f"    {{ name band{i:03d}; type cellSet; action new; source boxToCell; sourceInfo {{ box (-1 -1 {zb[i]:.9g}) (1 1 {zb[i+1]:.9g}); }} }}\n")
    ts.append(");\n")
    w(f"{c}/system/fuel/topoSetDict", "".join(ts))

    if steady:
        fschemes = """
ddtSchemes      { default steadyState; }
gradSchemes     { default Gauss linear; }
divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss linearUpwind grad(U);
    div(phi,K)      bounded Gauss linear;
    div(phi,h)      bounded Gauss linearUpwind grad(h);
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes   { default corrected; }
"""
        fsol = """
solvers
{
    rho { solver PCG; preconditioner DIC; tolerance 1e-8; relTol 0; }
    p_rgh { solver GAMG; smoother GaussSeidel; tolerance 1e-9; relTol 0.001; }
    "(U|h)" { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0.05; }
}
SIMPLE
{
    frozenFlow      %s;   // yes for the conduction-only comparison (U = 0, energy only)
    momentumPredictor yes;
    nNonOrthogonalCorrectors 0;
    pRefCell        0;
    pRefValue       1e5;
    rhoMin          2000;
    rhoMax          3000;
}
relaxationFactors
{
    // U 0.3 / h 0.5 / tight p_rgh tolerance are needed: looser settings (U 0.5, h 0.9, p_rgh relTol 0.01)
    // stall with large continuity errors that act as spurious mixing and under-predict the peak temperature.
    fields { rho 1; p_rgh 0.7; }
    equations { U 0.3; h %s; }
}
""" % (("no", "0.5") if a.gravity else ("yes", "1"))
    else:
        fschemes = """
ddtSchemes      { default Euler; }
gradSchemes     { default Gauss linear; }
divSchemes
{
    default         none;
    div(phi,U)      Gauss linearUpwind grad(U);
    div(phi,K)      Gauss linear;
    div(phi,h)      Gauss linearUpwind grad(h);
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes   { default corrected; }
"""
        fsol = """
solvers
{
    "rho.*" { solver PCG; preconditioner DIC; tolerance 1e-8; relTol 0; }
    p_rgh { solver GAMG; smoother GaussSeidel; tolerance 1e-8; relTol 0.01; }
    p_rghFinal { $p_rgh; relTol 0; }
    "(U|h)" { solver PBiCGStab; preconditioner DILU; tolerance 1e-9; relTol 0.05; }
    "(U|h)Final" { $U; relTol 0; }
}
PIMPLE
{
    momentumPredictor yes;
    nCorrectors     2;
    nNonOrthogonalCorrectors 0;
    pRefCell        0;
    pRefValue       1e5;
}
relaxationFactors { equations { ".*" 1; } }
"""
    w(f"{c}/system/fuel/fvSchemes", HDR.format(cls="dictionary", loc="system/fuel", obj="fvSchemes") + fschemes)
    w(f"{c}/system/fuel/fvSolution", HDR.format(cls="dictionary", loc="system/fuel", obj="fvSolution") + fsol)

    # ------------------------------------------------------------ region: graphite (solid)
    w(f"{c}/constant/graphite/thermophysicalProperties",
      HDR.format(cls="dictionary", loc="constant/graphite", obj="thermophysicalProperties") + f"""
thermoType
{{
    type            heSolidThermo;
    mixture         pureMixture;
    transport       constIso;
    thermo          hConst;
    equationOfState rhoConst;
    specie          specie;
    energy          sensibleEnthalpy;
}}
mixture
{{
    specie          {{ molWeight 12; }}
    transport       {{ kappa {K_GR}; }}
    thermodynamics  {{ Hf 0; Cp {CP_GR}; }}
    equationOfState {{ rho {RHO_GR}; }}
}}
""")
    w(f"{c}/constant/graphite/radiationProperties",
      HDR.format(cls="dictionary", loc="constant/graphite", obj="radiationProperties") + "radiation off;\nradiationModel none;\n")
    w(f"{c}/system/graphite/fvSchemes", HDR.format(cls="dictionary", loc="system/graphite", obj="fvSchemes") + f"""
ddtSchemes      {{ default {ddt}; }}
gradSchemes     {{ default Gauss linear; }}
divSchemes      {{ default none; }}
laplacianSchemes {{ default Gauss linear corrected; }}
interpolationSchemes {{ default linear; }}
snGradSchemes   {{ default corrected; }}
""")
    if steady:
        gs = """solvers { h { solver PCG; preconditioner DIC; tolerance 1e-10; relTol 0.01; } }
SIMPLE { nNonOrthogonalCorrectors 0; }
relaxationFactors { equations { h 1; } }
"""
    else:
        gs = """solvers { h { solver PCG; preconditioner DIC; tolerance 1e-9; relTol 0.01; } hFinal { $h; relTol 0; } }
PIMPLE { nNonOrthogonalCorrectors 0; }
"""
    w(f"{c}/system/graphite/fvSolution", HDR.format(cls="dictionary", loc="system/graphite", obj="fvSolution") + gs)
    for r in ("fuel", "graphite"):
        w(f"{c}/system/{r}/decomposeParDict", HDR.format(cls="dictionary", loc=f"system/{r}", obj="decomposeParDict") +
          "numberOfSubdomains 1;\nmethod simple;\ncoeffs { n (1 1 1); }\n")
    w(f"{c}/system/decomposeParDict", HDR.format(cls="dictionary", loc="system", obj="decomposeParDict") +
      "numberOfSubdomains 1;\nmethod simple;\ncoeffs { n (1 1 1); }\n")

    # ------------------------------------------------------------ Allrun / Allclean
    w(f"{c}/Allrun", textwrap.dedent(f"""\
        #!/bin/bash
        # 2D fuel-slot conjugate heat transfer case: mesh, split regions, set BCs/sources, run {app}.
        cd "${{0%/*}}" || exit 1
        [ -n "$WM_PROJECT_DIR" ] || . /usr/share/openfoam/etc/bashrc > /dev/null 2>&1
        SCRIPTS=../../scripts
        set -e
        blockMesh > log.blockMesh 2>&1
        splitMeshRegions -cellZonesOnly -overwrite > log.splitMeshRegions 2>&1
        topoSet -region fuel > log.topoSet.fuel 2>&1
        mkdir -p 0/fuel 0/graphite
        python3 $SCRIPTS/set_fields.py .          # writes 0/fuel/* and 0/graphite/* (coolant T(z) BC, init T)
        if grep -q "frozenFlow      no" system/fuel/fvSolution; then
            # stage 1: conduction-only initial field (frozen flow, U = 0), then stage 2: buoyant flow
            sed -i 's/frozenFlow      no/frozenFlow      yes/' system/fuel/fvSolution
            sed -i 's/^endTime .*/endTime 1000;/' system/controlDict
            {app} > log.{app}.stage1 2>&1
            sed -i 's/frozenFlow      yes/frozenFlow      no/' system/fuel/fvSolution
            sed -i 's/^endTime .*/endTime {a.iters};/' system/controlDict
        fi
        {app} > log.{app} 2>&1
        """))
    w(f"{c}/Allclean", "#!/bin/bash\ncd \"${0%/*}\" || exit 1\nrm -rf 0 [1-9]* constant/polyMesh constant/fuel/polyMesh constant/graphite/polyMesh postProcessing log.* processor*\n")
    for f in ("Allrun", "Allclean"):
        os.chmod(f"{c}/{f}", 0o755)


if __name__ == "__main__":
    main()
