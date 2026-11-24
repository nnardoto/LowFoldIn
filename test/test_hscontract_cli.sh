#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Nathanael Nardoto Batista
#
# Testes de linha de comando do hscontract:
#  - caso sintético (gerado pelo test_basis_contract com 'gen <dir>'):
#    contração, saídas, resumo, entradas intactas e re-execução identidade
#    sobre o resultado (valida o HSIO reduzido no fluxo existente);
#  - colisão de saída com entrada;
#  - códigos de saída em erro (argumentos, arquivo ausente, base ausente, l).
# Uso: bash test/test_hscontract_cli.sh   (a partir da raiz do repositório)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
PROFILE="${PROFILE:-release}"
export OPENBLAS_NUM_THREADS=1
fpm build --profile "$PROFILE" >/dev/null
binpath() { echo "$ROOT/$(fpm run --profile "$PROFILE" "$@" --runner echo 2>/dev/null | tail -n 1)"; }
HSC="$(binpath hscontract)"
fpm test test_basis_contract --profile "$PROFILE" >/dev/null 2>&1 || true
TESTBIN="$(find build -type f -name test_basis_contract -printf '%T@ %p\n' 2>/dev/null | sort -rn | head -n 1 | cut -d' ' -f2-)"
if [ -z "$TESTBIN" ]; then echo "  FAIL: test_basis_contract binary not found"; exit 1; fi

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
"$TESTBIN" gen "$TMP" >/dev/null
cd "$TMP"
pass=0; fail=0
ok()  { pass=$((pass + 1)); }
bad() { fail=$((fail + 1)); echo "  FAIL: $*"; }
check() { if eval "$1"; then ok; else bad "$2"; fi; }
# roda $HSC ARGS; captura saída em OUT e código em RC (não dispara set -e)
run() { OUT="$("$HSC" "$@" 2>&1)" && RC=0 || RC=$?; }

cat > contract.nml <<'EOF'
&HSCONTRACT
    spin_mode = 'none'
    h_file = 'case_H.hsio'
    s_file = 'case_S.hsio'
    structure_file = 'case.nml'
    output_prefix = 'red'
    basis = 'Aa s1p1', 'Bb s2'
    check_kgrid = 2, 2, 2
    tol_s_min = 1.0e-10
/
&RADIAL_RULE
    species = 'Aa'
    l = 1
    n_out = 0
/
&RADIAL_RULE
    species = 'Bb'
    l = 0
    n_out = 1
    coefficients = 0.6, 0.3
/
EOF

# 1. contração do caso sintético (6 -> 2 orbitais: Aa 4->1, Bb 2->1)
SUM_IN="$(md5sum case_H.hsio case_S.hsio)"
run contract.nml
check "[ $RC -eq 0 ]" "hscontract exit 0 (got $RC)"
check "[ -s red_H.hsio ]" "red_H.hsio written"
check "[ -s red_S.hsio ]" "red_S.hsio written"
check "echo \"\$OUT\" | grep -q 'orbitals       : 6 -> 2'" "summary reports 6 -> 2"
check "echo \"\$OUT\" | grep -q 'lambda_min(S_r)'" "summary reports lambda_min"
check "[ \"\$(md5sum case_H.hsio case_S.hsio)\" = \"\$SUM_IN\" ]" "inputs not modified"

# 2. o resultado reduzido é um HSIO válido: re-execução identidade (2 -> 2)
cat > contract2.nml <<'EOF'
&HSCONTRACT
    spin_mode = 'none'
    h_file = 'red_H.hsio'
    s_file = 'red_S.hsio'
    structure_file = 'case.nml'
    output_prefix = 'id'
    basis = 'Aa s1', 'Bb s1'
/
EOF
run contract2.nml
check "[ $RC -eq 0 ]" "identity re-run on reduced files exit 0 (got $RC)"
check "echo \"\$OUT\" | grep -q 'orbitals       : 2 -> 2'" "identity re-run keeps 2 orbitals"

# 3. colisão de saída com entrada (prefixo 'case' gera case_H.hsio = entrada)
cat > contract3.nml <<'EOF'
&HSCONTRACT
    spin_mode = 'none'
    h_file = 'case_H.hsio'
    s_file = 'case_S.hsio'
    structure_file = 'case.nml'
    output_prefix = 'case'
    basis = 'Aa s1', 'Bb s1'
/
EOF
run contract3.nml
check "[ $RC -ne 0 ]" "output/input collision is an error"
check "echo \"\$OUT\" | grep -qi 'overwrite'" "collision message names the problem"

# 4. caso collinear (H_up/H_dn/S compartilham a mesma T espacial)
cat > contract4.nml <<'NMLEOF'
&HSCONTRACT
    spin_mode = 'collinear'
    h_up_file = 'coll_H_up.hsio'
    h_dn_file = 'coll_H_dn.hsio'
    s_file = 'coll_S.hsio'
    structure_file = 'case.nml'
    output_prefix = 'redc'
    basis = 'Aa s1p1', 'Bb s2'
    check_kgrid = 2, 2, 2
    tol_s_min = 1.0e-10
/
&RADIAL_RULE
    species = 'Aa'
    l = 1
    n_out = 0
/
&RADIAL_RULE
    species = 'Bb'
    l = 0
    n_out = 1
    coefficients = 0.6, 0.3
/
NMLEOF
run contract4.nml
check "[ $RC -eq 0 ]" "collinear contract exit 0 (got $RC)"
check "[ -s redc_H_up.hsio ]" "redc_H_up.hsio written"
check "[ -s redc_H_dn.hsio ]" "redc_H_dn.hsio written"
check "[ -s redc_S.hsio ]" "redc_S.hsio written"
check "echo \"\$OUT\" | grep -q 'orbitals       : 6 -> 2'" "collinear summary reports 6 -> 2"
check "echo \"\$OUT\" | grep -q 'collinear'" "summary names collinear mode"

# 5. codigos de saida em erro
if "$HSC" >/dev/null 2>&1; then bad "no arguments should fail"; else ok; fi
if "$HSC" missing.nml >/dev/null 2>&1; then bad "missing nml should fail"; else ok; fi
cat > bad.nml <<'EOF'
&HSCONTRACT
    spin_mode = 'none'
    h_file = 'case_H.hsio'
    s_file = 'case_S.hsio'
    structure_file = 'case.nml'
    output_prefix = 'zz'
    basis = 'Aa s1'
/
EOF
if "$HSC" bad.nml >/dev/null 2>&1; then bad "species without basis entry should fail"; else ok; fi
cat > badk.nml <<'EOF'
&HSCONTRACT
    spin_mode = 'none'
    h_file = 'case_H.hsio'
    s_file = 'case_S.hsio'
    structure_file = 'case.nml'
    output_prefix = 'zz'
    basis = 'Aa s1', 'Bb s1'
    tol_s_min = 1.0e-10
/
&RADIAL_RULE
    species = 'Aa'
    l = 4
    n_out = 0
/
EOF
if "$HSC" badk.nml >/dev/null 2>&1; then bad "invalid l should fail"; else ok; fi

echo "[test_hscontract_cli] $pass passed, $fail failed"
[ $fail -eq 0 ]
