#!/bin/bash
# Continue a finished case from its latest time with (possibly edited) settings:
#   scripts/continue_case.sh <case> <newEndTime> [U_relax h_relax p_rgh_relax]
case_dir=$1; end=$2
[ -n "$WM_PROJECT_DIR" ] || . /usr/share/openfoam/etc/bashrc > /dev/null 2>&1
cd "$case_dir" || exit 1
sed -i "s/^endTime .*/endTime $end;/" system/controlDict
if [ -n "$3" ]; then
    sed -i "s/equations { U [0-9.]*; h [0-9.]*; }/equations { U $3; h $4; }/; s/fields { rho 1; p_rgh [0-9.]*; }/fields { rho 1; p_rgh $5; }/" system/fuel/fvSolution
fi
[ -d constant/fuel/polyMesh/sets ] || topoSet -region fuel -time 0 > log.topoSet.fuel 2>&1
app=$(sed -n 's/^application \(.*\);/\1/p' system/controlDict)
[ -f log.$app ] && grep -q "^End" log.$app && mv log.$app log.$app.part$(ls log.$app.part* 2>/dev/null | wc -l)
$app > log.$app 2>&1
