# Compare numeric files by row and column, with absolute + relative tolerance.
# Usage: awk -v atol=1e-10 -v rtol=1e-9 -f compare-numbers.awk reference actual
function fail(message) {
    print FILENAME ": " message > "/dev/stderr"
    failed=1
    exit 1
}
function number(text) {
    gsub(/[dD]/, "e", text)
    if (text !~ /^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/)
        fail("invalid number: " text)
    if (sprintf("%.17e", text+0) ~ /[iInN]/) fail("non-finite number: " text)
    return text+0
}
function abs(x) {return x<0 ? -x : x}
NF==0 {next}
FILENAME==ARGV[1] {
    columns[++expected_rows]=NF
    for (i=1; i<=NF; i++) expected[expected_rows,i]=number($i)
    next
}
{
    row++
    if (row>expected_rows) fail("extra row " row)
    if (NF!=columns[row]) fail("wrong column count at row " row)
    for (i=1; i<=NF; i++) {
        actual=number($i); reference=expected[row,i]
        if (abs(actual-reference)>atol+rtol*abs(reference))
            fail("row " row ", column " i ": expected " reference ", got " actual)
    }
}
END {
    if (!failed && (row!=expected_rows || !expected_rows))
        fail("expected " expected_rows+0 " rows, got " row+0)
}
