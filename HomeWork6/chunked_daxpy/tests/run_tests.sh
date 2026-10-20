#!/usr/bin/env bash
# Tests for chunked_daxpy: run with "make test" from the chunked_daxpy folder.
#
# Each test writes a small config file, runs the program and checks
#   - the exit code (0 = all checks passed, 1 = bad input, 2 = a check failed)
#   - optionally, something in the output or in the HDF5 file.

cd "$(dirname "$0")/.." || exit 1
PROG=./chunked_daxpy
TMP=tests/tmp
rm -rf "$TMP"   # start from a clean folder: never test old output files
mkdir -p "$TMP"
pass=0
fail=0

# run_case NAME EXPECTED_EXIT CONFIG_TEXT [GREP_PATTERN_IN_OUTPUT]
run_case() {
    local name=$1 expected=$2 config=$3 pattern=$4
    printf '%s\n' "$config" > "$TMP/$name.cfg"
    out=$("$PROG" "$TMP/$name.cfg" 2>&1)
    code=$?
    if [ "$code" -ne "$expected" ]; then
        echo "FAIL  $name: exit code $code, expected $expected"
        echo "$out" | sed 's/^/      /'
        fail=$((fail + 1))
    elif [ -n "$pattern" ] && ! grep -q -- "$pattern" <<< "$out"; then
        echo "FAIL  $name: output does not contain '$pattern'"
        echo "$out" | sed 's/^/      /'
        fail=$((fail + 1))
    else
        echo "PASS  $name"
        pass=$((pass + 1))
    fi
}

# number of rows (= chunks) of /chunks/d in an HDF5 file
count_chunks() { h5ls "$1/chunks/d" | sed -E 's/.*\{([0-9]+),.*/\1/'; }

echo "== Correct input files (expected exit code 0)"
run_case divisible        0 $'n = 32\nchunk_size = 8\na = 2\nx = 1\ny = 1\noutput = '$TMP'/t1.h5' "number_of_chunks = 4"
run_case not_divisible    0 $'n = 25\nchunk_size = 8\na = 2\nx = 1\ny = 1\noutput = '$TMP'/t2.h5' "last chunk: 1 elements"
run_case chunk_bigger_n   0 $'n = 5\nchunk_size = 100\na = 2\nx = 1\ny = 1\noutput = '$TMP'/t3.h5' "number_of_chunks = 1"
run_case chunk_size_one   0 $'n = 7\nchunk_size = 1\na = 2\nx = 1\ny = 1\noutput = '$TMP'/t4.h5' "number_of_chunks = 7"
run_case number_of_chunks 0 $'n = 10\nnumber_of_chunks = 3\na = 2\nx = 1\ny = 1\noutput = '$TMP'/t5.h5' "chunk_size = 4"
run_case steps_and_1e6    0 $'n = 1e6\nchunk_size = 1000\na = 1\nx = 0.1\nx_step = 1e-6\ny = 0.2\noutput = '$TMP'/t6.h5' "equal within rounding"
run_case no_extension     0 $'n = 4\nchunk_size = 2\na = 1\nx = 1\ny = 1\noutput = '$TMP'/t7' "t7.h5"

echo "== HDF5 content"
for t in "t1 4" "t2 4" "t3 1" "t4 7" "t5 3"; do
    set -- $t
    got=$(count_chunks "$TMP/$1.h5")
    if [ "$got" = "$2" ]; then echo "PASS  $1.h5 has $2 chunks"; pass=$((pass + 1))
    else echo "FAIL  $1.h5 has '$got' chunks, expected $2"; fail=$((fail + 1)); fi
done
# n = 25, chunk_size = 8: chunk lengths 8 8 8 1, 4 partial sums, last row padded with NaN
if h5dump -d /chunks/length "$TMP/t2.h5" | grep -q '8, 8, 8, 1' &&
   h5ls "$TMP/t2.h5" | grep -q 'partial_chunk_sum *Dataset {4}' &&
   h5dump -d /chunks/d "$TMP/t2.h5" | grep -q 'nan'; then
    echo "PASS  t2.h5 chunk lengths 8 8 8 1, 4 partial sums, NaN padding"; pass=$((pass + 1))
else
    echo "FAIL  t2.h5 chunk lengths"; fail=$((fail + 1))
fi
# the partial sums and the total sums are stored in the file
if h5dump -d /partial_chunk_sum -a /total_sum/chunked "$TMP/t1.h5" > /dev/null; then
    echo "PASS  t1.h5 contains partial_chunk_sum and total_sum"; pass=$((pass + 1))
else
    echo "FAIL  t1.h5 partial_chunk_sum / total_sum"; fail=$((fail + 1))
fi

echo "== Wrong input files (expected exit code 1)"
run_case missing_n        1 $'chunk_size = 8\na = 2\nx = 1\ny = 1'              "missing required variable 'n'"
run_case missing_a        1 $'n = 10\nchunk_size = 8\nx = 1\ny = 1'             "missing required variable 'a'"
run_case both_chunk_keys  1 $'n = 10\nchunk_size = 2\nnumber_of_chunks = 5\na = 1\nx = 1\ny = 1' "exactly one"
run_case no_chunk_key     1 $'n = 10\na = 1\nx = 1\ny = 1'                       "exactly one"
run_case chunk_zero       1 $'n = 10\nchunk_size = 0\na = 1\nx = 1\ny = 1'       "positive integer"
run_case n_not_integer    1 $'n = 2.5\nchunk_size = 1\na = 1\nx = 1\ny = 1'      "positive integer"
run_case a_not_number     1 $'n = 10\nchunk_size = 2\na = abc\nx = 1\ny = 1'     "must be a number"
run_case too_many_chunks  1 $'n = 3\nnumber_of_chunks = 5\na = 1\nx = 1\ny = 1'  "cannot be larger"
if "$PROG" "$TMP/does_not_exist.cfg" > /dev/null 2>&1; then
    echo "FAIL  missing file accepted"; fail=$((fail + 1))
else
    echo "PASS  missing_file"; pass=$((pass + 1))
fi

echo "== Warnings"
run_case typo_warning     0 $'n = 4\nchunk_size = 2\na = 1\nx = 1\ny = 1\nx_stp = 3\noutput = '$TMP'/t8.h5' "unknown variable 'x_stp'"

echo ""
echo "$pass passed, $fail failed"
[ "$fail" -eq 0 ]
