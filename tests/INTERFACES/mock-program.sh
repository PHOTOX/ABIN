#!/bin/bash
# Replace only the external chemistry executable, leaving the adapter intact.
set -euo pipefail
error() {
    echo "$*" >> "$FIXTURE_WORK/fixture-errors.log"
    touch "$FIXTURE_WORK/mock.failed"
    exit 90
}
normalize() { sed "s|$FIXTURE_WORK|@WORK@|g" "$1"; }
printf 'launch\n' >> "$FIXTURE_WORK/launches"
: > "$FIXTURE_WORK/actual.argv"
for arg in "$@"; do printf '%s\n' "$arg" >> "$FIXTURE_WORK/actual.argv"; done
normalize "$FIXTURE_WORK/actual.argv" > "$FIXTURE_WORK/normalized.argv"
diff -u "$FIXTURE_STEP/argv" "$FIXTURE_WORK/normalized.argv" >> "$FIXTURE_WORK/fixture-errors.log" 2>&1 || error 'Executable arguments differ'
if [[ -f $FIXTURE_STEP/stdin-input ]]; then
    cat > "$FIXTURE_WORK/actual.stdin"
    input=$(cat "$FIXTURE_STEP/stdin-input")
    cmp -s "$FIXTURE_WORK/actual.stdin" "$FIXTURE_WORK/$input" || error 'Executable stdin differs from generated input'
fi
[[ -d $FIXTURE_STEP/inputs ]] || error 'No input references for external launch'
checked=0
while IFS= read -r -d '' reference; do
    relative=${reference#"$FIXTURE_STEP/inputs/"}
    actual="$FIXTURE_WORK/$relative"
    [[ -f $actual ]] || error "Missing generated input: $relative"
    normalize "$actual" > "$FIXTURE_WORK/normalized.input"
    diff -u "$reference" "$FIXTURE_WORK/normalized.input" >> "$FIXTURE_WORK/fixture-errors.log" 2>&1 || error "Generated input differs: $relative"
    checked=$((checked+1))
done < <(find "$FIXTURE_STEP/inputs" -type f -print0)
((checked>0)) || error 'Empty input reference directory'
if [[ ${input:-} == MNDO-EH/* ]]; then
    awk -f "$FIXTURE_MNDO_VALIDATOR" "$FIXTURE_WORK/$input" >> "$FIXTURE_WORK/fixture-errors.log" 2>&1 ||
        error 'Invalid MNDO keyword/title/geometry layout'
fi
# These files are prerecorded fixtures; nothing is computed from the input.
if [[ -d $FIXTURE_STEP/exports ]]; then cp -R "$FIXTURE_STEP/exports/." .; fi
if [[ -f $FIXTURE_STEP/stdout ]]; then cat "$FIXTURE_STEP/stdout"; fi
