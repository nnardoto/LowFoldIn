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
echo "6. help"
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
