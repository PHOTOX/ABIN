#!/bin/bash
# Common functions for BAGEL BASH interface
# NOTE: This file is meant to be sourced, not executed!
 
function exec_bagel {
    local input_file=$1
    local output_file=$2
    local timestep=$3
    # BAGEL does not print its input so we will prepend it.
    cp "$input_file" "$output_file"
    $BAGELEXE "$input_file" >> "$output_file" 2>&1
    local exitcode=$?

    # Concatenate all output files
    all_output="$output_file".all
    sep="\n############################\n"
    header="${sep}TIMESTEP = $timestep\n$(date)${sep}"
    echo -e "${header}" >> "$all_output"
    cat "$output_file" >> "$all_output"
    return $exitcode
}

function bagel_error {
    local errmsg="$1"
    local bagel_output="$2"
    local copy="$3"
    if [[ -f "$bagel_output" && -n "$copy" ]]; then
        cp "$bagel_output" "$copy"
    fi
    >&2 echo -e "${errmsg}"
    >&2 echo "Inspect file $(dirname "$0")/$copy"
    exit 2
}

function file_exists {
    if [[ ! -f $1 ]];then
        >&2 echo "File \"$1\" does not exist!"
        exit 2
    fi
}

### FUNCTIONS FOR GENERATING BAGEL INPUT
 
# This function must be called first to start the bagel input
function specify_molecule {
    # XYZ coordinates provided by ABIN in geom.dat
    local xyz=$1
    local basis=$2
    local df_basis=$3
    local inp=$4
    local natom
    natom=$(wc -l < "$xyz")

    # NOTE: $basis and $df_basis MUST be defined at the top of the file!
    cat > "$inp" << EOF
{"bagel": [
  {
    "title": "molecule",
    "basis": "$basis",
    "df_basis": "$df_basis",
    "angstrom": "true",
    "geometry": [
EOF

    # Geometry specification in json format
    awk -v natom="$natom" '{printf "\t{\"atom\": \"%s\", \"xyz\": [ %.10f, %.10f, %.10f ]}", $1, $2, $3, $4} NR != natom {print ","}END{print ""}'  "$xyz" >> "$inp"
    echo -e "    ]\n  }," >> "$inp"
}

function generate_hf_orbitals {
    local hf_variant=$1
    local charge=$2
    local nspin=$3
    local inp=$4
    cat >> "$inp" << EOF
  {
    "title" : "$hf_variant",
    "charge" : $charge,
    "nspin" : $nspin
  },
EOF
}

function load_orbitals {
    local orbfile=$1
    local inp=$2
    cat >> "$inp" << EOF
  {
    "title": "load_ref",
    "file": "$orbfile",
    "continue_geom": false
  },
EOF
}

function save_orbitals {
    local orbfile=$1
    local inp=$2
    cat >> "$inp" << EOF
  {
    "title": "save_ref",
    "file": "$orbfile"
  },
EOF
}

function print_molden {
    local moldenfile=$1
    local inp=$2
    cat >> "$inp" << EOF
  {
    "title": "print",
    "file": "$moldenfile",
    "orbitals": true
  }
EOF
}

# Specify for which state we calculate the forces
# (for now we only support one target state)
# This needs to be called before the CAS section.
function print_forces {
    local target_state=$1
    local inp=$2
    cat >> "$inp" << EOF
  {
    "title" : "forces",
    "export" : true,
    "grads" : [
      {
        "title": "force",
        "target": $target_state
      }
    ],
EOF
}

function print_cas {
    local method=$1
    local nspin=$2
    local charge=$3
    local nstate=$4
    local nact=$5
    local nclosed=$6
    local thresh_CASSCF=$7
    local thresh_CASPT2=$8
    local maxiter_CASSCF=$9
    local maxiter_CASPT2=${10}
    local shift=${11}
    local inp=${12}

    if [[ $method == "xms_caspt2" ]]; then

        print_casscf "caspt2" "$nspin" "$charge" "$nstate" "$nact" "$nclosed" "$thresh_CASSCF" "$maxiter_CASSCF" "$inp"
        print_caspt2 "true" "$thresh_CASPT2" "$maxiter_CASPT2" "$shift" "$inp"

    elif [[ $method == "ms_caspt2" ]]; then

        print_casscf "caspt2" "$nspin" "$charge" "$nstate" "$nact" "$nclosed" "$thresh_CASSCF" "$maxiter_CASSCF" "$inp"
        print_caspt2 "false" "$thresh_CASPT2" "$maxiter_CASPT2" "$shift" "$inp"

    elif [[ $method == "sa_casscf" ]]; then

        print_casscf "casscf" "$nspin" "$charge" "$nstate" "$nact" "$nclosed" "$thresh_CASSCF" "$maxiter_CASSCF" "$inp"
        # Remove extra dangling comma
        sed -i '$ s/,$//' "$inp"

    else

        >&2 echo "ERROR: Unknown method ($method). Specify one of \"xms_caspt2\", \"ms_caspt2\" or \"sa_casscf\" in bagel.inp"
        exit 2

    fi
}

function print_casscf {
    cat >> "$9" << EOF
    "method": [{
      "title": "$1",
      "nspin": $2,
      "charge": $3,
      "nstate": $4,
      "nact": $5,
      "nclosed": $6,
      "thresh": $7,
      "maxiter": $8,
EOF
# NOTE: We have to leave this JSON section unclosed since we might need to append CASPT2 input
}

function print_caspt2 {
    # TODO: Make it possible to use a real shift instead of imaginary
    cat >> "$5" << EOF
      "smith": {
        "method": "caspt2",
        "xms": $1,
        "thresh": $2,
        "maxiter": $3,
        "shift": $4,
        "shift_imag": true
      }
EOF
}
