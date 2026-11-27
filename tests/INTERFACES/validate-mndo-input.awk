# Check the MOPAC-style record layout used by the MNDO-EH adapters.
# Unlike a snapshot comparison, this rejects a fourth keyword line even if
# somebody accidentally updates the reference to match the broken generator.
function fail(message) {
    print "Invalid MNDO input at line " NR ": " message > "/dev/stderr"
    failed=1
    exit 1
}
BEGIN {section="keywords"}
section=="keywords" {
    if (++keyword_lines>3) fail("MOPAC-style input permits at most three keyword lines")
    if (length($0)>80) fail("keyword line exceeds MNDO2005's 80-column limit")
    for (i=1; i<=NF; i++) {
        split(tolower($i), option, "=")
        if (option[1]=="iroot") roots=option[2]+0
        if (option[1]=="ncigrd") gradients=option[2]+0
        if (option[1]=="icross") crossing=option[2]+0
        if (option[1]=="mciref") reference_mode=option[2]+0
        if (option[1]=="nciref") references=option[2]+0
        if (option[1]=="ici1" || option[1]=="ici2") active+=option[2]
    }
    if ($0 !~ /[[:space:]][+&][[:space:]]*$/) {section="titles"; titles=2}
    next
}
section=="titles" {
    if (--titles==0) section="geometry"
    next
}
section=="geometry" {
    if (NF==0) {
        if (!atoms) fail("geometry is empty or title records are misaligned")
        section=((reference_mode==1 || reference_mode==4) && references>0 ? "references" : "states")
        next
    }
    if (NF!=7 || $1 !~ /^([A-Za-z][A-Za-z]?|[0-9]+)$/)
        fail("expected atom and three Cartesian coordinate/flag pairs")
    atoms++
    next
}
section=="references" {
    if (NF!=active) fail("reference occupation count does not match active orbitals")
    for (i=1; i<=NF; i++) if ($i !~ /^[012]$/) fail("reference occupations must be 0, 1 or 2")
    if (++reference_rows==references) section="states"
    next
}
section=="states" {
    for (i=1; i<=NF; i++) {
        if ($i !~ /^[1-9][0-9]*$/ || $i>roots || seen[$i]++) fail("invalid CI state label")
        states++
    }
}
END {
    if (failed) exit 1
    if (section!="states" || !atoms || roots<1 || (crossing!=1 && crossing!=7))
        fail("incomplete keywords, titles or geometry")
    if (states!=gradients) fail("CI state list does not match ncigrd")
}
