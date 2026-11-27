#!/bin/bash
# Sourced by the BAGEL-EH entry points. Requires Bash and standard Unix tools.
bagel_eh_error() {
    echo "ERROR: BAGEL-EH: $*" >&2
    exit 2
}

# Prefer .out when BAGEL provides both supported export formats.
bagel_eh_export() {
    local suffix
    for suffix in .out .txt; do
        if [[ -f $1$suffix ]]; then
            printf '%s\n' "$1$suffix"
            return
        fi
    done
    bagel_eh_error "BAGEL did not export $1.out or $1.txt"
}

# Validate and convert an energy or vector export for ABIN.
bagel_eh_read_export() {
    local mode=$1 expected=$2 file=$3
    awk -v mode="$mode" -v expected="$expected" '
    # Validate BAGEL exports before ABIN can consume them.
    function fail(message) {
        print "ERROR: BAGEL-EH: " FILENAME ": " message > "/dev/stderr"
        failed = 1
        exit 2
    }
    function number(value, text) {
        gsub(/[dD]/, "e", value)
        if (value !~ /^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/)
            fail("Invalid number: " value)
        text = sprintf("%.17e", value + 0)
        if (text ~ /[iInN]/) fail("Non-finite number: " value)
        return text
    }
    mode == "energy" {
        for (i = 1; i <= NF; i++) {
            print number($i)
            rows++
        }
        next
    }
    NF == 0 || /^[[:space:]]*#/ { next }
    rows == 0 && NF == 3 && $1 == "0" && $2 == "1" && $3 == "2" { next }
    {
        if (NF != 4 || $1 !~ /^[0-9]+$/ || $1 + 0 != rows)
            fail("Invalid or unordered atom row: " $0)
        x = number($2); y = number($3); z = number($4)
        print x, y, z
        rows++
    }
    END {
        if (!failed && rows != expected) fail("Expected " expected " rows, got " rows + 0)
    }
    ' "$file"
}

run_bagel_eh() (
    set -eE -o pipefail
    export LC_ALL=C
    local mode=${1:-}
    [[ $# == 4 && ( $mode == gradients || $mode == nac ) ]] ||
        bagel_eh_error 'Expected mode timestep bead accuracy'
    local timestep=$2 bead=$3 accuracy=$4
    [[ $timestep =~ ^[+-]?[0-9]+$ && $bead =~ ^[0-9]+$ && $accuracy =~ ^[0-9]+$ ]] ||
        bagel_eh_error 'Invalid timestep, bead index or accuracy'
    accuracy=$((10#$accuracy))
    (( accuracy > 0 )) || bagel_eh_error 'Accuracy must be positive'
    cd "$(dirname "${BASH_SOURCE[0]}")"
    local output=../engrad.dat.$bead temporary= nac_temporary=
    local cache=.bagel-eh.$bead
    [[ $mode == gradients ]] || output=../nacm.dat
    trap 'echo "ERROR: BAGEL-EH: Interface failed at line $LINENO" >&2' ERR
    trap 'if [[ -n $temporary ]]; then rm -f -- "$temporary"; fi;
          if [[ -n $nac_temporary ]]; then rm -f -- "$nac_temporary"; fi' EXIT
    rm -f -- "$output"
    if [[ $mode == gradients ]]; then
        rm -f -- "$cache.accuracy" ../nacm.dat
    fi
    [[ -f bagel.inp ]] || bagel_eh_error 'bagel.inp is missing'
    source ./bagel.inp
    case $method in
        sa_casscf|ms_caspt2|xms_caspt2) ;;
        *) bagel_eh_error "Unsupported BAGEL-EH method: $method" ;;
    esac
    local key
    for key in charge nspin nact nclosed maxiter_CASSCF maxiter_CASPT2 maxziter; do
        [[ ${!key} =~ ^-?[0-9]+$ ]] || bagel_eh_error "Invalid integer setting: $key"
    done
    [[ -n ${BAGELEXE//[[:space:]]/} ]] || bagel_eh_error 'BAGELEXE is empty'
    local nstates flags expected flag
    IFS= read -r nstates || bagel_eh_error 'Missing state count'
    [[ $nstates =~ ^[0-9]+$ ]] || bagel_eh_error 'Invalid state count'
    nstates=$((10#$nstates))
    expected=$nstates
    [[ $mode == gradients ]] || expected=$((nstates * (nstates - 1) / 2))
    (( nstates > 0 && expected > 0 )) || bagel_eh_error 'Empty state or NAC request'
    read -r -a flags || bagel_eh_error 'Missing state flags'
    [[ ${#flags[@]} == "$expected" ]] || bagel_eh_error 'Ehrenfest requires all states or all upper-triangular NAC pairs'
    for flag in "${flags[@]}"; do
        [[ $flag == 1 ]] || bagel_eh_error 'Ehrenfest requires all states or all upper-triangular NAC pairs'
    done
    # Cache identity includes the effective settings, including environment overrides.
    local request cached_accuracy
    request=$(printf '%s\n' "$timestep" "$bead" "$nstates" "$basis" "$df_basis" \
        "$hf_variant" "$method" "$charge" "$nspin" "$nact" "$nclosed" \
        "$maxiter_CASSCF" "$maxiter_CASPT2" "$maxziter" "$thresh_CASPT2" "$shift" "$BAGELEXE")
    if [[ $mode == nac ]]; then
        [[ -f $cache.accuracy && -f $cache.request && -f $cache.geom && -f $cache.nac ]] ||
            bagel_eh_error 'NAC request requires a successful combined gradient calculation'
        [[ $(cat -- "$cache.request") == "$request" ]] &&
            cmp -s -- ../geom.dat."$bead" "$cache.geom" ||
            bagel_eh_error 'NAC request does not match the combined gradient calculation'
        IFS= read -r cached_accuracy < "$cache.accuracy"
        [[ $cached_accuracy =~ ^[0-9]+$ ]] && (( cached_accuracy >= accuracy )) ||
            bagel_eh_error 'NAC accuracy exceeds the combined gradient calculation; rerun gradients at the required accuracy'
        temporary=$(mktemp "$output.tmp.XXXXXX")
        cp -- "$cache.nac" "$temporary"
        mv -f -- "$temporary" "$output"
        temporary=
        return
    fi
    local prefix=input$bead
    local input=$prefix.json log=$prefix.out reference=hf
    if [[ -f orbitals.archive ]]; then
        reference=load
    fi
    # AWK writes JSON directly; numeric tokens are validated before emitting them.
    awk -v ns="$nstates" -v accuracy="$accuracy" \
        -v basis="$basis" -v df_basis="$df_basis" -v hf="$hf_variant" \
        -v method="$method" -v charge="$charge" -v nspin="$nspin" \
        -v nact="$nact" -v nclosed="$nclosed" -v maxcas="$maxiter_CASSCF" \
        -v maxpt="$maxiter_CASPT2" -v maxz="$maxziter" \
        -v threshpt="$thresh_CASPT2" -v shift="$shift" -v reference="$reference" \
        '
        # Generate the all-state BAGEL gradient/NAC input from an ABIN geometry.
        function fail(message) {
            print "ERROR: BAGEL-EH: " message > "/dev/stderr"
            exit 2
        }
        function number(value, text) {
            gsub(/[dD]/, "e", value)
            if (value !~ /^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/)
                fail("Invalid number: " value)
            text = sprintf("%.17e", value + 0)
            if (text ~ /[iInN]/) fail("Non-finite number: " value)
            return text
        }
        function quote(value, result, i, c) {
            result = "\""
            for (i = 1; i <= length(value); i++) {
                c = substr(value, i, 1)
                if (c == "\\" || c == "\"") result = result "\\"
                if (c == "\n") c = "\\n"
                if (c == "\r") c = "\\r"
                if (c == "\t") c = "\\t"
                result = result c
            }
            return result "\""
        }
        {
            if (NF != 4 || $1 !~ /^[A-Za-z][A-Za-z]?$/) fail("Invalid geometry row: " $0)
            geometry[++atoms] = "{\"atom\":" quote($1) ",\"xyz\":[" number($2) "," number($3) "," number($4) "]}"
        }
        END {
            if (!atoms) fail("Empty BAGEL geometry")
            print "{\"bagel\":["
            print "{\"title\":\"molecule\",\"basis\":" quote(basis) ",\"df_basis\":" quote(df_basis) ",\"angstrom\":true,\"geometry\":["
            for (a = 1; a <= atoms; a++) print geometry[a] (a < atoms ? "," : "")
            print "]},"
            if (reference == "load")
                print "{\"title\":\"load_ref\",\"file\":\"orbitals\",\"continue_geom\":false},"
            else
                print "{\"title\":" quote(hf) ",\"charge\":" number(charge) ",\"nspin\":" number(nspin) "},"
            print "{\"title\":\"forces\",\"export\":true,\"export_single\":false,\"grads\":["
            separator = ""
            for (i = 0; i < ns; i++) {
                printf "%s{\"title\":\"force\",\"target\":%d,\"maxziter\":%d}", separator, i, maxz
                separator = ",\n"
            }
            for (i = 0; i < ns; i++) {
                for (j = i + 1; j < ns; j++) {
                    printf "%s{\"title\":\"nacme\",\"target\":%d,\"target2\":%d,\"nacmtype\":\"full\",\"maxziter\":%d}", separator, i, j, maxz
                    separator = ",\n"
                }
            }
            print "],\"method\":[{"
            print "\"title\":" quote(method == "sa_casscf" ? "casscf" : "caspt2") ","
            printf "\"charge\":%d,\"nspin\":%d,\"nstate\":%d,\"nact\":%d,\"nclosed\":%d,\"maxiter\":%d,\"thresh\":1e-%d", charge, nspin, ns, nact, nclosed, maxcas, accuracy + 1
            if (method != "sa_casscf") {
                print ",\"smith\":{\"method\":\"caspt2\",\"ms\":true,\"xms\":" (method == "xms_caspt2" ? "true" : "false") ","
                print "\"maxiter\":" number(maxpt) ",\"thresh\":" number(threshpt) ",\"shift\":" number(shift) ",\"shift_imag\":true}"
            }
            print "}]}"
            print ",{\"title\":\"save_ref\",\"file\":\"orbitals\"}"
            print "]}"
        }
        ' ../geom.dat."$bead" > "$input"
    # Export files have no timestep stamp: never reuse a previous calculation.
    rm -f -- ENERGY.out ENERGY.txt FORCE_*.out FORCE_*.txt NACME_*.out NACME_*.txt
    cp -- "$input" "$log"
    local status=0
    # BAGELEXE is a trusted shell command from bagel.inp or the environment.
    # This supports MPI launchers and quoted executable paths containing spaces.
    bash -c "exec $BAGELEXE \"\$@\"" -- "$input" >> "$log" 2>&1 || status=$?
    printf '\nTIMESTEP = %s\n' "$timestep" >> "$log.all"
    cat -- "$log" >> "$log.all"
    if (( status != 0 )) || grep -qi 'EXCEPTION RAISED' "$log"; then
        bagel_eh_error "BAGEL failed; inspect $log (exit status $status)"
    fi
    local natoms state i j file
    natoms=$(awk 'END {print NR}' ../geom.dat."$bead")
    temporary=$(mktemp "$output.tmp.XXXXXX")
    file=$(bagel_eh_export ENERGY)
    bagel_eh_read_export energy "$nstates" "$file" > "$temporary"
    for ((state=0; state<nstates; state++)); do
        file=$(bagel_eh_export FORCE_$state)
        bagel_eh_read_export vector "$natoms" "$file" >> "$temporary"
    done
    nac_temporary=$(mktemp "$cache.nac.tmp.XXXXXX")
    for ((i=0; i<nstates-1; i++)); do
        for ((j=i+1; j<nstates; j++)); do
            file=$(bagel_eh_export NACME_${i}_$j)
            bagel_eh_read_export vector "$natoms" "$file" >> "$nac_temporary"
        done
    done
    # Publish only after every energy, gradient and NAC export has been validated.
    mv -f -- "$nac_temporary" "$cache.nac"
    nac_temporary=
    cp -- ../geom.dat."$bead" "$cache.geom"
    printf '%s\n' "$request" > "$cache.request"
    mv -f -- "$temporary" "$output"
    temporary=
    printf '%s\n' "$accuracy" > "$cache.accuracy"
)
