#!/usr/bin/env python3
"""Write initial/boundary fields for a fuel-slot case after blockMesh + splitMeshRegions.

Coolant boundary (graphite outer faces coolantL/coolantR):
  h == 0 : fixedValue T = T_cool(z)                      (wetted-surface temperature prescribed)
  h  > 0 : mixed (Robin) q'' = h (T_wall - T_cool(z)), valueFraction = h/(h + k/dn)
T_cool(z) = 550 C at z=0 -> 650 C at z=H (linear), or 600 C uniform (coolant = "uniform").
Internal T initialised to T_cool(z) in both regions.
"""
import json, os, re, sys
import numpy as np

HDR = """FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    location    "0/{reg}";
    object      {obj};
}}
dimensions      {dim};
"""


def _foam_list_body(path):
    txt = open(path).read()
    txt = re.sub(r"FoamFile\s*\{.*?\}", "", txt, count=1, flags=re.S)
    txt = re.sub(r"//.*", "", txt)
    m = re.search(r"(\d+)\s*\(", txt)
    return int(m.group(1)), txt[m.end():]


def read_mesh_centres(case, reg):
    """Cell centres (mean of face centres; exact for these rectangular cells) and patch face centres."""
    pm = f"{case}/constant/{reg}/polyMesh"
    n, body = _foam_list_body(f"{pm}/points")
    pts = np.array([[float(v) for v in s.split()] for s in re.findall(r"\(([^()]*)\)", body)[:n]])
    n, body = _foam_list_body(f"{pm}/faces")
    faces = np.array([list(map(int, f.split())) for f in re.findall(r"\d+\(([^()]*)\)", body)[:n]])
    fc = pts[faces].mean(axis=1)
    n, body = _foam_list_body(f"{pm}/owner")
    own = np.array(body.split(")")[0].split(), dtype=int)[:n]
    n, body = _foam_list_body(f"{pm}/neighbour")
    nei = np.array(body.split(")")[0].split(), dtype=int)[:n]
    ncell = own.max() + 1
    acc = np.zeros((ncell, 3)); cnt = np.zeros(ncell)
    np.add.at(acc, own, fc[:len(own)]); np.add.at(cnt, own, 1)
    np.add.at(acc, nei, fc[:len(nei)]); np.add.at(cnt, nei, 1)
    cells = acc / cnt[:, None]
    btxt = open(f"{pm}/boundary").read()
    patches = {}
    for m in re.finditer(r"(\w+)\s*\{[^{}]*?nFaces\s+(\d+);\s*startFace\s+(\d+);", btxt, flags=re.S):
        nf, sf = int(m.group(2)), int(m.group(3))
        patches[m.group(1)] = fc[sf:sf + nf]
    return cells, patches


def lst(vals):
    return f"nonuniform List<scalar> {len(vals)}\n(\n" + "\n".join(f"{v:.6f}" for v in vals) + "\n)"


def main(case):
    p = json.load(open(f"{case}/caseParams.json"))
    H = p["H"]

    def Tc(z):
        z = np.asarray(z)
        if p["coolant"] == "uniform":
            return np.full(z.shape, 873.15)
        return 823.15 + 100.0 * z / H

    os.makedirs(f"{case}/0/fuel", exist_ok=True)
    os.makedirs(f"{case}/0/graphite", exist_ok=True)
    # ---------------- graphite
    cells, patches = read_mesh_centres(case, "graphite")
    bcs = []
    for name in ("coolantL", "coolantR"):
        z = patches[name][:, 2]
        if p["h"] > 0:
            dn = p["w_graphite"] / p["nxg"] / 2.0
            f = p["h"] / (p["h"] + p["k_graphite"] / dn)
            bcs.append(f"""    {name}
    {{
        type            mixed;
        refValue        {lst(Tc(z))};
        refGradient     uniform 0;
        valueFraction   uniform {f:.8g};
        value           {lst(Tc(z))};
    }}""")
        else:
            bcs.append(f"""    {name}
    {{
        type            fixedValue;
        value           {lst(Tc(z))};
    }}""")
    bcs.append("""    "graphite_to_.*"
    {
        type            compressible::turbulentTemperatureCoupledBaffleMixed;
        Tnbr            T;
        kappaMethod     solidThermo;
        value           uniform 873.15;
    }
    "(top|bottom)"  { type zeroGradient; }
    frontAndBack    { type empty; }""")
    open(f"{case}/0/graphite/T", "w").write(
        HDR.format(cls="volScalarField", reg="graphite", obj="T", dim="[0 0 0 1 0 0 0]") +
        f"internalField   {lst(Tc(cells[:, 2]))};\nboundaryField\n{{\n" + "\n".join(bcs) + "\n}\n")
    open(f"{case}/0/graphite/p", "w").write(
        HDR.format(cls="volScalarField", reg="graphite", obj="p", dim="[1 -1 -2 0 0 0 0]") +
        "internalField uniform 1e5;\nboundaryField\n{\n    \".*\" { type calculated; value uniform 1e5; }\n    frontAndBack { type empty; }\n}\n")

    # ---------------- fuel
    cells, patches = read_mesh_centres(case, "fuel")
    open(f"{case}/0/fuel/T", "w").write(
        HDR.format(cls="volScalarField", reg="fuel", obj="T", dim="[0 0 0 1 0 0 0]") +
        f"internalField   {lst(Tc(cells[:, 2]))};\n" + """boundaryField
{
    "fuel_to_.*"
    {
        type            compressible::turbulentTemperatureCoupledBaffleMixed;
        Tnbr            T;
        kappaMethod     fluidThermo;
        value           uniform 873.15;
    }
    "(top|bottom)"  { type zeroGradient; }
    frontAndBack    { type empty; }
}
""")
    open(f"{case}/0/fuel/U", "w").write(
        HDR.format(cls="volVectorField", reg="fuel", obj="U", dim="[0 1 -1 0 0 0 0]") +
        "internalField uniform (0 0 0);\nboundaryField\n{\n    \".*\" { type noSlip; }\n    frontAndBack { type empty; }\n}\n")
    open(f"{case}/0/fuel/p", "w").write(
        HDR.format(cls="volScalarField", reg="fuel", obj="p", dim="[1 -1 -2 0 0 0 0]") +
        "internalField uniform 1e5;\nboundaryField\n{\n    \".*\" { type calculated; value uniform 1e5; }\n    frontAndBack { type empty; }\n}\n")
    open(f"{case}/0/fuel/p_rgh", "w").write(
        HDR.format(cls="volScalarField", reg="fuel", obj="p_rgh", dim="[1 -1 -2 0 0 0 0]") +
        "internalField uniform 1e5;\nboundaryField\n{\n    \".*\" { type fixedFluxPressure; value uniform 1e5; }\n    frontAndBack { type empty; }\n}\n")
    for reg in ("fuel", "graphite", "."):
        if os.path.exists(f"{case}/0/{reg}/cellToRegion"):
            os.remove(f"{case}/0/{reg}/cellToRegion")


if __name__ == "__main__":
    main(sys.argv[1])
