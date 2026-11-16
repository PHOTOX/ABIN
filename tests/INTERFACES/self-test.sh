#!/bin/bash
# Prove that fixture mismatches fail instead of being accepted or re-blessed.
set -euo pipefail
export LC_ALL=C
SUITE=$(cd "$(dirname "$0")" && pwd)
work=$(mktemp -d /tmp/abin-fixture-selftest.XXXXXX)
trap 'rm -rf -- "$work"' EXIT
fail() {
    echo "ERROR: Interface fixture self-test: $*" >&2
    if [[ -f $work/log ]]; then cat "$work/log" >&2; fi
    exit 1
}
compare() { awk -v atol=1e-10 -v rtol=1e-9 -f "$SUITE/compare-numbers.awk" "$work/reference" "$work/actual"; }
printf '1 -2 3\n' > "$work/reference"
printf '1D0 -2d0 3.0000000001\n' > "$work/actual"
compare
for bad in '1 -2 -3' '1 -2' '1 -2 NaN' '1 -2 1e999' ''; do
    printf '%s\n' "$bad" > "$work/actual"
    if compare > "$work/log" 2>&1; then echo "Comparator accepted invalid data: $bad" >&2; exit 1; fi
done
printf '1 -2 3\n4 5 6\n' > "$work/actual"
if compare > "$work/log" 2>&1; then echo 'Comparator accepted an extra row' >&2; exit 1; fi
# Reproduce the MNDO2005 GETGEO failure independently of the input snapshots.
cat > "$work/four-keyword-lines" <<'INPUT'
OM2 jop=-2 igeom=1 iform=1 nsav15=3 kharge=0 imult=0 +
ipubo=1 ktrial=11 immdp=-1 kci=5 ici1=8 ici2=8 iroot=5 +
nciref=1 mciref=1 levexc=2 ioutci=1 kitscf=50000 iuvcd=2 ifast=2 +
icross=1 ncigrd=5 lroot=1
Timestep 0

C 1.153944 1 -0.82293659 1 -0.17499567 1

 1 2 3 4 5
INPUT
if awk -f "$SUITE/validate-mndo-input.awk" "$work/four-keyword-lines" > "$work/log" 2>&1; then
    fail 'MNDO layout check accepted four MOPAC-style keyword lines'
fi
grep -q 'at most three keyword lines' "$work/log" || fail 'Wrong MNDO layout diagnostic'
# A continuation marker beyond column 80 is not read by MNDO2005. This
# reported input truncates immdp to an unrecognized standalone I at column 80.
cat > "$work/long-keyword-lines" <<'INPUT'
OM2 jop=-2 igeom=1 iform=1 nsav15=3 kharge=0 imult=0 lroot=1 ipubo=1 ktrial=11 immdp=-1 kci=5 +
 ici1=8 ici2=8 iroot=5 icross=1 ncigrd=5 nciref=1 mciref=1 levexc=2 ioutci=1 kitscf=50000 iuvcd=2 ifast=2
Timestep 0

C 1.153944 1 -0.82293659 1 -0.17499567 1

2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0
1 2 3 4 5
INPUT
if awk -f "$SUITE/validate-mndo-input.awk" "$work/long-keyword-lines" > "$work/log" 2>&1; then
    fail 'MNDO layout check accepted keyword lines longer than 80 columns'
fi
grep -q '80-column limit' "$work/log" || fail 'Wrong MNDO keyword-length diagnostic'
awk 'NF!=16' "$SUITE/cases/mndo-5-state/steps/01-gradients/inputs/MNDO-EH/input003" > "$work/missing-reference"
if awk -f "$SUITE/validate-mndo-input.awk" "$work/missing-reference" > "$work/log" 2>&1; then
    fail 'MNDO layout check accepted mciref=1 without reference occupations'
fi
grep -q 'reference occupation count' "$work/log" || fail 'Wrong missing-reference diagnostic'
# Automatic reference modes must not consume a row of explicit occupations.
# Otherwise its leading 2s become five requests for the same electronic state.
sed 's/mciref=1/mciref=3/' "$SUITE/cases/mndo-5-state/steps/01-gradients/inputs/MNDO-EH/input003" > "$work/extra-reference"
if awk -f "$SUITE/validate-mndo-input.awk" "$work/extra-reference" > "$work/log" 2>&1; then
    fail 'MNDO layout check accepted explicit occupations in automatic reference mode'
fi
grep -q 'invalid CI state label' "$work/log" || fail 'Wrong extra-reference diagnostic'
# Mutations must start from a passing fixture. Otherwise an unrelated adapter
# change can mask the intended mismatch and make the self-check misleading.
if ! INTERFACE_FIXTURES_DIR="$SUITE/cases" KEEP_INTERFACE_TESTS=0 \
        bash "$SUITE/run.sh" mndo-3-state > "$work/log" 2>&1; then
    fail 'Unmodified MNDO baseline failed; inspect the adapter/fixture mismatch below.'
fi
mkdir "$work/cases"
for change in input arguments values truncated energy-output gradient-output repeated-state-output; do
    cp -R "$SUITE/cases/mndo-3-state" "$work/cases/$change"
    case_dir="$work/cases/$change"
    step="$case_dir/steps/01-gradients"
    case $change in
        input)
            printf '\nDeliberate input mismatch\n' >> "$step/inputs/MNDO-EH/input003"
            ;;
        arguments) printf 'unexpected argument\n' > "$step/argv" ;;
        values)
            awk 'NR==1 {$1=$1+1} {print}' "$step/expected/engrad.dat.003" > "$work/changed"
            mv "$work/changed" "$step/expected/engrad.dat.003"
            ;;
        truncated)
            sed '$d' "$step/expected/engrad.dat.003" > "$work/changed"
            mv "$work/changed" "$step/expected/engrad.dat.003"
            ;;
        energy-output)
            # The mock still reports SCF TOTAL ENERGY and COMPUTATION TIME.
            sed '/State  2 /d' "$step/stdout" > "$work/changed"
            mv "$work/changed" "$step/stdout"
            ;;
        gradient-output)
            sed '$d' "$step/exports/fort.15" > "$work/changed"
            mv "$work/changed" "$step/exports/fort.15"
            ;;
        repeated-state-output)
            sed 's/FOR STATE[[:space:]]*[123]/FOR STATE 2/' "$step/exports/fort.15" > "$work/changed"
            mv "$work/changed" "$step/exports/fort.15"
            ;;
    esac
    if INTERFACE_FIXTURES_DIR="$work/cases" KEEP_INTERFACE_TESTS=0 \
            bash "$SUITE/run.sh" "$change" > "$work/log" 2>&1; then
        fail "Fixture runner accepted altered $change reference"
    fi
    # Ensure the intended comparison failed, rather than a setup error.
    case $change in
        input) diagnostic='Generated input differs' ;;
        arguments) diagnostic='Executable arguments differ' ;;
        values) diagnostic='row 1, column 1' ;;
        truncated) diagnostic='extra row' ;;
        energy-output) diagnostic='missing energy for state 2 of 3' ;;
        gradient-output) diagnostic='incomplete gradient for state 2: expected 6 components, got 5' ;;
        repeated-state-output) diagnostic='duplicate gradient for state 2' ;;
    esac
    grep -q "$diagnostic" "$work/log" || fail "Unexpected failure while checking $change reference"
done
echo 'Fixture self-check passed: input/argument differences, numeric differences, invalid numbers, wrong shapes and incomplete converged exports are rejected.'
