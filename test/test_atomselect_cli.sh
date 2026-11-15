#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Nathanael Nardoto Batista
#
# test/test_atomselect_cli.sh
#
# Command-line tests for the atomselect executable (docs/design/contracts.md §6).
# Run from the repository root:   bash test/test_atomselect_cli.sh
# Exits non-zero if any check fails.

set -u

cd "$(dirname "$0")/.." || exit 1

fpm build >/dev/null 2>&1 || { echo "fpm build failed"; exit 1; }
BIN=$(fpm run atomselect --runner echo 2>/dev/null | tail -1)
if [ -z "$BIN" ] || [ ! -x "$BIN" ]; then
    BIN=$(find build -path '*app/atomselect' -type f | head -1)
fi
if [ -z "$BIN" ] || [ ! -x "$BIN" ]; then
    echo "atomselect binary not found"; exit 1
fi

STRUCT=test/data/slab8.nml
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

npass=0
nfail=0

pass() { npass=$((npass + 1)); echo "  ok   $1"; }
fail() { nfail=$((nfail + 1)); echo "  FAIL $1"; }

# check NAME CONDITION...   (CONDITION is evaluated as a command)
check() {
    local name=$1; shift
    if "$@"; then pass "$name"; else fail "$name"; fi
}

# expect_error NAME CMD...   command must exit with code 1 and print the error prefix
expect_error() {
    local name=$1; shift
    local out rc
    out=$("$@" 2>&1); rc=$?
    if [ "$rc" -eq 1 ] && grep -q '^\[atomselect error\]' <<<"$out"; then
        pass "$name"
    else
        fail "$name (rc=$rc)"; echo "$out" | sed 's/^/       | /'
    fi
}

# ----------------------------------------------------------------------
echo "1. inspect"
out=$("$BIN" inspect "$STRUCT" --xyz "$TMP/s.xyz" --table "$TMP/s.txt"); rc=$?
check "inspect exit 0" test "$rc" -eq 0
nrows=$(awk '$1 ~ /^[0-9]+$/ && NF == 5' <<<"$out" | wc -l | tr -d ' ')
check "inspect prints 8 rows" test "$nrows" -eq 8
check "xyz has 10 lines" test "$(wc -l < "$TMP/s.xyz" | tr -d ' ')" -eq 10
check "xyz first line is nat" test "$(head -1 "$TMP/s.xyz")" = "8"
tblrows=$(awk '$1 ~ /^[0-9]+$/ && NF == 5' "$TMP/s.txt" | wc -l | tr -d ' ')
check "table file has 8 rows" test "$tblrows" -eq 8
expect_error "inspect unknown option" "$BIN" inspect "$STRUCT" --bogus
expect_error "inspect missing structure" "$BIN" inspect "$TMP/nonexistent.nml"

# ----------------------------------------------------------------------
echo "2. fingerprint"
fp1=$("$BIN" fingerprint "$STRUCT")
fp2=$("$BIN" fingerprint "$STRUCT")
check "fingerprint is 16 hex chars" grep -Eq '^[0-9a-f]{16}$' <<<"$fp1"
check "fingerprint stable" test "$fp1" = "$fp2"

# ----------------------------------------------------------------------
echo "3. select top/bottom -> partition"
REG=$TMP/regions.nml
"$BIN" select "$STRUCT" --regions "$REG" --name top --atoms "1-2,3 ,4" >/dev/null
check "select top exit 0" test $? -eq 0
check "regions file created" test -f "$REG"
"$BIN" select "$STRUCT" --name bottom --atoms "5-8" --regions "$REG" >/dev/null
check "select bottom exit 0 (options in any order)" test $? -eq 0
out=$("$BIN" validate "$STRUCT" "$REG"); rc=$?
check "validate exit 0" test "$rc" -eq 0
check "partition yes" grep -q 'partition *: yes' <<<"$out"
check "fingerprint stored" grep -q "fingerprint *= *'$fp1'" "$REG"

# ----------------------------------------------------------------------
echo "4. select overlapping mid"
"$BIN" select "$STRUCT" --regions "$REG" --name mid --atoms "3-6" >/dev/null
check "select mid exit 0" test $? -eq 0
out=$("$BIN" validate "$STRUCT" "$REG"); rc=$?
check "validate exit 0" test "$rc" -eq 0
check "overlap top/mid reported" grep -q "overlap: 'top' and 'mid' share 2 atoms" <<<"$out"
check "overlap bottom/mid reported" grep -q "overlap: 'bottom' and 'mid' share 2 atoms" <<<"$out"
check "partition no" grep -q 'partition *: no' <<<"$out"

# ----------------------------------------------------------------------
echo "5. errors"
expect_error "duplicate name"  "$BIN" select "$STRUCT" --regions "$REG" --name top --atoms "1"
expect_error "atoms 0"         "$BIN" select "$STRUCT" --regions "$REG" --name x1 --atoms "0"
expect_error "atoms 9"         "$BIN" select "$STRUCT" --regions "$REG" --name x2 --atoms "9"
expect_error "atoms 1,1"       "$BIN" select "$STRUCT" --regions "$REG" --name x3 --atoms "1,1"
expect_error "atoms 3-1"       "$BIN" select "$STRUCT" --regions "$REG" --name x4 --atoms "3-1"
expect_error "name top-bad"    "$BIN" select "$STRUCT" --regions "$REG" --name top-bad --atoms "1"
expect_error "name total"      "$BIN" select "$STRUCT" --regions "$REG" --name total --atoms "1"
expect_error "missing --atoms" "$BIN" select "$STRUCT" --regions "$REG" --name x5
expect_error "unknown subcommand" "$BIN" frobnicate "$STRUCT"
expect_error "unknown option"  "$BIN" validate "$STRUCT" "$REG" --bogus

MOVED=$TMP/slab8_moved.nml
sed 's/^Ga   0.000000   0.000000   9.500000/Ga   0.100000   0.000000   9.500000/' "$STRUCT" > "$MOVED"
check "moved structure differs" test "$("$BIN" fingerprint "$MOVED")" != "$fp1"
expect_error "fingerprint mismatch" "$BIN" validate "$MOVED" "$REG"

NOFP=$TMP/nofp.nml
cat > "$NOFP" <<'EOF'
&REGION_FILE
    format_version = 1
    natoms         = 8
/
&REGION
    name  = 'all'
    atoms = 1, 2, 3, 4, 5, 6, 7, 8
/
EOF
expect_error "no fingerprint rejected" "$BIN" validate "$STRUCT" "$NOFP"
out=$("$BIN" validate "$STRUCT" "$NOFP" --allow-no-fingerprint); rc=$?
check "no fingerprint accepted with --allow-no-fingerprint" test "$rc" -eq 0
check "partition yes (single region)" grep -q 'partition *: yes' <<<"$out"

# ----------------------------------------------------------------------
echo "6. layers"
TILT=test/data/slab_tilted.nml
out=$("$BIN" layers "$STRUCT" --normal 0 0 1 --table "$TMP/layers.txt" --regions "$TMP/layers.nml"); rc=$?
check "layers exit 0" test "$rc" -eq 0
nrows=$(awk '$1 ~ /^[0-9]+$/ && NF == 7 && $7 ~ /^[0-9]+$/' <<<"$out" | wc -l | tr -d ' ')
check "layers table has 8 rows" test "$nrows" -eq 8
check "layers header" grep -q 'Index *Element *X *Y *Z *Projected *Layer' <<<"$out"
check "atom 1 in layer 4" grep -Eq '^ *1 +Ga +.* 4$' <<<"$out"
check "4 layers" grep -q 'nlayer *: 4' <<<"$out"
check "summary layer 1 = 7-8" grep -Eq '^# +1 +2\.000000 +2 +7-8$' <<<"$out"
tblrows=$(awk '$1 ~ /^[0-9]+$/ && NF == 7 && $7 ~ /^[0-9]+$/' "$TMP/layers.txt" | wc -l | tr -d ' ')
check "layers table file has 8 rows" test "$tblrows" -eq 8
check "layer_001 region" grep -q "name  = 'layer_001'" "$TMP/layers.nml"
check "layer_004 region" grep -q "name  = 'layer_004'" "$TMP/layers.nml"
out=$("$BIN" validate "$STRUCT" "$TMP/layers.nml"); rc=$?
check "layer regions validate" test "$rc" -eq 0
check "layer regions are a partition" grep -q 'partition *: yes' <<<"$out"
expect_error "layers regions duplicate names" "$BIN" layers "$STRUCT" --normal 0 0 1 --regions "$TMP/layers.nml"
out=$("$BIN" layers "$STRUCT" --normal 0 0 1 --tol 3); rc=$?
check "layers tol 3 merges to 1 layer" grep -q 'nlayer *: 1' <<<"$out"
out=$("$BIN" layers "$STRUCT" --normal 0 0 1 --refs "9.5, 2.0,4.5,7"); rc=$?
check "layers --refs exit 0" test "$rc" -eq 0
check "layers --refs 4 layers" grep -q 'nlayer *: 4' <<<"$out"
expect_error "layers --refs atom far from refs" "$BIN" layers "$STRUCT" --normal 0 0 1 --refs "2,4.5,7"
out=$("$BIN" layers "$TILT" --normal 1 1 0); rc=$?
check "tilted layers exit 0" test "$rc" -eq 0
check "tilted: 4 layers along (1,1,0)" grep -q 'nlayer *: 4' <<<"$out"
check "tilted: no split warning" test "$(grep -c 'atomselect warning' <<<"$out")" -eq 0
out=$("$BIN" layers "$TILT" --normal 0 0 1)
check "tilted: 2 layers along z" grep -q 'nlayer *: 2' <<<"$out"
WRAP=$TMP/wrapped.nml
awk '/^&ATOMS/{inatoms=1} inatoms && $1=="Ga" && $4=="2.000000"{$4="19.000000"} {print}' "$STRUCT" > "$WRAP"
sed -i.bak 's/^Ga   0.000000   0.000000   9.500000/Ga   0.000000   0.000000   0.500000/' "$WRAP"
out=$("$BIN" layers "$WRAP" --normal 0 0 1 --tol 1.5); rc=$?
check "split slab: exit 0" test "$rc" -eq 0
check "split slab: warning printed" grep -q '^\[atomselect warning\].*periodic' <<<"$out"
expect_error "layers zero normal" "$BIN" layers "$STRUCT" --normal 0 0 0
expect_error "layers missing --normal" "$BIN" layers "$STRUCT"
expect_error "layers --normal two values" "$BIN" layers "$STRUCT" --normal 0 1
expect_error "layers --normal non-number" "$BIN" layers "$STRUCT" --normal 0 x 1
expect_error "layers tol 0" "$BIN" layers "$STRUCT" --normal 0 0 1 --tol 0

# ----------------------------------------------------------------------
echo "7. geometric select"
GREG=$TMP/geom.nml
out=$("$BIN" select "$STRUCT" --regions "$GREG" --name upper --normal 0 0 1 --smin 6 --smax 10); rc=$?
check "select --normal exit 0" test "$rc" -eq 0
check "select --normal atoms 1-4" grep -q "upper': 4 atoms \[1-4\]" <<<"$out"
out=$("$BIN" select "$STRUCT" --regions "$GREG" --name low --normal 0 0 -1 --smin -3); rc=$?
check "select -z normal smin only" grep -q "low': 2 atoms \[7-8\]" <<<"$out"
out=$("$BIN" select "$STRUCT" --regions "$GREG" --name boxed --box "1 * * * 4 8"); rc=$?
check "select --box exit 0" test "$rc" -eq 0
check "select --box atoms 3,6" grep -q "boxed': 2 atoms \[3,6\]" <<<"$out"
out=$("$BIN" select "$TILT" --regions "$TMP/tilt_sel.nml" --name t --normal 1 1 0 --smin 6.5 --smax 10); rc=$?
check "tilted select atoms 1-4" grep -q "t': 4 atoms \[1-4\]" <<<"$out"
expect_error "--atoms and --normal exclusive" "$BIN" select "$STRUCT" --regions "$GREG" --name x6 --atoms "1" --normal 0 0 1 --smin 0
expect_error "--normal and --box exclusive" "$BIN" select "$STRUCT" --regions "$GREG" --name x7 --box "* * * * * *" --normal 0 0 1 --smin 0
expect_error "--normal without range" "$BIN" select "$STRUCT" --regions "$GREG" --name x8 --normal 0 0 1
expect_error "--smin without --normal" "$BIN" select "$STRUCT" --regions "$GREG" --name x9 --atoms "1" --smin 0
expect_error "empty geometric selection" "$BIN" select "$STRUCT" --regions "$GREG" --name x10 --normal 0 0 1 --smin 50
expect_error "--box wrong count" "$BIN" select "$STRUCT" --regions "$GREG" --name x11 --box "0 1 0 1"

# ----------------------------------------------------------------------
echo "8. surface"
SREG=$TMP/surf.nml
out=$("$BIN" surface "$STRUCT" --regions "$SREG" --normal 0 0 1 --top 1 --bottom 1); rc=$?
check "surface exit 0" test "$rc" -eq 0
check "surface top = 1-2" grep -q "region 'top': layers 4-4, 2 atoms \[1-2\]" <<<"$out"
check "surface bottom = 7-8" grep -q "region 'bottom': layers 1-1, 2 atoms \[7-8\]" <<<"$out"
"$BIN" validate "$STRUCT" "$SREG" >/dev/null
check "surface regions validate" test $? -eq 0
out=$("$BIN" surface "$STRUCT" --regions "$SREG" --normal 0 0 1 --top 2 --top-name top2); rc=$?
check "surface --top only" grep -q "region 'top2': layers 3-4, 4 atoms \[1-4\]" <<<"$out"
out=$("$BIN" surface "$STRUCT" --regions "$SREG" --normal 0 0 -1 --bottom 1 --bottom-name zbot); rc=$?
check "surface -z normal: bottom = highest z" grep -q "region 'zbot': .*\[1-2\]" <<<"$out"
expect_error "surface 2N > nlayer" "$BIN" surface "$STRUCT" --regions "$SREG" --normal 0 0 1 --top 3 --bottom 3 --top-name a --bottom-name b
expect_error "surface no top/bottom" "$BIN" surface "$STRUCT" --regions "$SREG" --normal 0 0 1
expect_error "surface missing --regions" "$BIN" surface "$STRUCT" --normal 0 0 1 --top 1
expect_error "surface duplicate name" "$BIN" surface "$STRUCT" --regions "$SREG" --normal 0 0 1 --top 1
expect_error "surface top 0" "$BIN" surface "$STRUCT" --regions "$SREG" --normal 0 0 1 --top 0 --top-name z0
TREG=$TMP/surf_tilt.nml
out=$("$BIN" surface "$TILT" --regions "$TREG" --normal 1 1 0 --tol 0.5 --top 1 --bottom 2); rc=$?
check "tilted surface exit 0" test "$rc" -eq 0
check "tilted surface top = 1-2" grep -q "region 'top': layers 4-4, 2 atoms \[1-2\]" <<<"$out"
check "tilted surface bottom = 5-8" grep -q "region 'bottom': layers 1-2, 4 atoms \[5-8\]" <<<"$out"
out=$("$BIN" surface "$TILT" --regions "$TMP/surf_tilt2.nml" --normal 1 1 0 --refs "2,4.5,7,9.5" --top 1); rc=$?
check "tilted surface with --refs" grep -q "region 'top': layers 4-4, 2 atoms \[1-2\]" <<<"$out"
expect_error "tilted surface 2N > nlayer" "$BIN" surface "$TILT" --regions "$TMP/surf_tilt3.nml" --normal 1 1 0 --top 3 --bottom 2

# ----------------------------------------------------------------------
echo "9. help"
"$BIN" help | grep -q 'atomselect select'
check "help lists subcommands" test $? -eq 0
"$BIN" >/dev/null
check "no args exits 0" test $? -eq 0

echo
if [ "$nfail" -ne 0 ]; then
    echo "atomselect CLI tests: $npass passed, $nfail FAILED"
    exit 1
fi
echo "atomselect CLI tests: $npass passed"
