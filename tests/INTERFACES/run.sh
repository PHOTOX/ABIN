#!/bin/bash
# Run the real Ehrenfest adapters against static input/output fixtures.
set -euo pipefail
export LC_ALL=C
SUITE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$SUITE/../.." && pwd)
FIXTURES=${INTERFACE_FIXTURES_DIR:-$SUITE/cases}
selection=${1:-all}
if [[ $selection == --list ]]; then
    for case_dir in "$FIXTURES"/*; do basename "$case_dir"; done
    exit 0
fi
run_case() (
    set -euo pipefail
    case_dir=$1
    source "$case_dir/parameters.sh"
    work=$(mktemp -d /tmp/abin-interface.XXXXXX)
    trap 'if [[ ${KEEP_INTERFACE_TESTS:-0} == 1 ]]; then echo "Artifacts: $work"; else rm -rf -- "$work"; fi' EXIT
    cp -R "$ROOT/interfaces/$adapter" "$work/$adapter"
    cp -R "$case_dir/setup/." "$work/"
    cp "$SUITE/mock-program.sh" "$work/mock-program"
    chmod +x "$work/mock-program"
    export FIXTURE_WORK=$work
    export FIXTURE_MNDO_VALIDATOR="$SUITE/validate-mndo-input.awk"
    # Do not inherit a user's real executable, site setup or shell startup file.
    cat > "$work/SetEnvironment.sh" <<'ENV'
MNDOEXE="$FIXTURE_WORK/mock-program"
MOLPROEXE="$FIXTURE_WORK/mock-program"
BAGELEXE="$FIXTURE_WORK/mock-program"
export MNDOEXE MOLPROEXE BAGELEXE
ENV
    cd "$work"
    for step in "$case_dir"/steps/*; do
        expected_launches=1
        source "$step/parameters.sh"
        export FIXTURE_STEP=$step
        : > launches
        : > fixture-errors.log
        rm -f mock.failed
        # Seed stale public output: it must not contaminate the new result.
        printf 'stale output\n' > engrad.dat.003
        printf 'stale output\n' > nacm.dat
        status=0
        env -u BASH_ENV -u ENV -u MNDOEXE -u MOLPROEXE \
            BAGELEXE="$work/mock-program" NCORES=2 \
            bash "$adapter/$script" "${arguments[@]}" < "$step/request" > interface.log 2>&1 || status=$?
        if ((status!=0)) || [[ -f mock.failed ]]; then
            cat interface.log fixture-errors.log >&2
            echo "FAILED: $(basename "$case_dir")/$(basename "$step") (exit $status)" >&2
            exit 1
        fi
        [[ $(wc -l < launches) -eq $expected_launches ]] || {
            echo "Incorrect executable launch count in $step" >&2; exit 1;
        }
        compared=0
        while IFS= read -r -d '' reference; do
            relative=${reference#"$step/expected/"}
            [[ -f $relative ]] || { echo "Missing extracted output: $relative" >&2; exit 1; }
            awk -v atol=1e-10 -v rtol=1e-9 -f "$SUITE/compare-numbers.awk" "$reference" "$relative"
            compared=$((compared+1))
        done < <(find "$step/expected" -type f -print0)
        ((compared>0)) || { echo "No numeric references in $step" >&2; exit 1; }
        echo "PASS $(basename "$case_dir")/$(basename "$step")"
    done
)
passed=0
failed=0
for case_dir in "$FIXTURES"/*; do
    name=$(basename "$case_dir")
    [[ $selection == all || $name == "$selection" || $name == "$selection"-* ]] || continue
    # Start a separate Bash so errexit remains active inside each case.
    if ( set -e; export -f run_case; export SUITE ROOT; env -u BASH_ENV -u ENV bash -e -u -o pipefail -c 'run_case "$1"' -- "$case_dir" ); then
        passed=$((passed+1))
    else
        failed=$((failed+1))
    fi
done
((passed+failed>0)) || { echo "No matching interface case: $selection" >&2; exit 2; }
printf '%s cases passed, %s failed\n' "$passed" "$failed"
((failed==0))
