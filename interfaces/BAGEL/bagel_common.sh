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

    echo "TIMESTEP = $timestep" >> "$output_file".all
    date >> "$output_file".all
    echo "####################" >> "$output_file".all
    cat "$output_file" >> "$output_file".all
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
    # Name of the bagel input file we are creating
    local inp=$2
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
    local inp=$1
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
    local inp=$2
    # CAS section
    if [[ $method == "xms_caspt2" ]]; then
        print_cas "caspt2" "$thresh_CASSCF" "$maxiter_CASSCF" "$input"
        print_caspt2 "true" "$thresh_CASPT2" "$maxiter_CASPT2" "$input"
    elif [[ $method == "ms_caspt2" ]]; then
        print_cas "caspt2" "$thresh_CASSCF" "$maxiter_CASSCF" "$input"
        print_caspt2 "false" "$thresh_CASPT2" "$maxiter_CASPT2" "$input"
    elif [[ $method == "sa_casscf" ]]; then
        print_cas "casscf" "$thresh_CASSCF" "$maxiter_CASSCF" "$input"
        # Remove extra dangling comma
        sed -i '$ s/,$//' "$input"
    else
        >&2 echo "ERROR: Unknown method ($method). Specify one of \"xms_caspt2\", \"ms_caspt2\" or \"sa_casscf\" in bagel.inp"
        exit 2
    fi
}

function print_casscf {
    local method=$1
    local thresh=$2
    local maxiter=$3
    local inp=$4
    cat >> "$inp" << EOF
    "method": [{
      "title": "$method",
      "nspin": $nspin,
      "charge": $charge,
      "maxiter": $maxiter,
      "thresh": $thresh,
      "nact": $nact,
      "nclosed": $nclosed,
      "nstate": $nstates,
EOF
# NOTE: We have to leave this JSON section unclosed since we might need to append CASPT2 input
}

function print_caspt2 {
    # XMS-CASPT2 - "true" or "false"
    local xms=$1
    local thresh=$2
    local maxiter=$3
    local inp=$4
    # TODO: Make it possible to use a real shift instead of imaginary
    cat >> "$inp" << EOF
      "smith": {
        "method": "caspt2",
        "xms": "$xms",
        "shift": $shift,
        "shift_imag": true,
        "maxiter": $maxiter,
        "thresh": $thresh
      }
EOF
}

function converge_casscf {
    local input=$1
    local output=$2
    local casscf_error='EXCEPTION RAISED:  Max iteration reached during the second-order optimization'
    generate_input "$input" "$method"

    # Execute BAGEL
    exec_bagel "$input" "$output" "$timestep"
    local returncode=$?
    if grep -q "$casscf_error" "$output"; then
        error=true
        local savefile=$output.casscf_error.$timestep
        cp "$output" "$savefile"
        >&2 echo "ERROR: CASSCF did not converge, see file $savefile"
	# Return if we reached the maximum threshold value,
	# further increase would lead to inaccurate results.
	if awk "BEGIN{exit !($thresh_CASSCF >= $max_thresh_CASSCF)}" ;then
	    return 1
	fi
	# Retry again 5*thresh
	thresh_CASSCF=$(awk "BEGIN{print $thresh_CASSCF*5}")
        maxiter_CASSCF=600
        >&2 echo "Trying again with thresh=$thresh_CASSCF and maxiter=$maxiter_CASSCF"
	converge_casscf "$input" "$output"
        returncode=$?
    fi
    return $returncode
}
