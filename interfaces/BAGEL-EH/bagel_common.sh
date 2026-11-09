#!/bin/bash
# Sourced by the BAGEL-EH entry points. Keep configuration compatible with bagel.inp.
run_bagel_eh() {
    local mode=$1
    shift
    cd "$(dirname "${BASH_SOURCE[0]}")" || return 2
    if [[ ! -f bagel.inp ]]; then
        echo "ERROR: BAGEL-EH/bagel.inp is missing" >&2
        return 2
    fi
    source ./bagel.inp || return 2
    export basis df_basis hf_variant method charge nspin nact nclosed
    export maxiter_CASSCF maxiter_CASPT2 maxziter thresh_CASPT2 shift BAGELEXE
    python3 ./bagel_eh.py "$mode" "$@"
}
