#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Nathanael Nardoto Batista
#
# Validacao de bandas full vs reduzida perto de E_F (fora do executavel
# hscontract, como definido no contrato):
#   1. modelos sinteticos 1D NAO ORTOGONAIS, multi-zeta (Aa s2p1d1), shell d
#      totalmente desacoplada das shells s/p com energias onsite altas;
#   2. hscontract com T = I (sem regras): bandas identicas em precisao
#      numerica sobre TODAS as bandas;
#   3. hscontract com poda d -> 0 (contrao nao trivial, exata para as
#      bandas s/p): bandas na janela E_F = 0 +- 3 eV identicas;
#   3b. CONTRACAO APROXIMADA s2 -> 1 com pesos [0.9, 0.1] nos dois zetas
#      (nao normalizados) + poda d exata, no modelo com o zeta 2 do s em
#      +5 eV (estado removido fora da janela; 4 estados s/p comparaveis):
#      max |dE|, RMS e extremos com limite deterministico do modelo;
#   4. as bandas sao calculadas pelo executavel LowFoldIn EXISTENTE (sem
#      modificacao) sobre os .hsio completos e reduzidos, no MESMO caminho
#      k e na mesma referencia de energia (E_F = 0, absoluta);
#   5. comparacao via test_basis_contract bandcmp (conjuntos ordenados por
#      k na janela; tolerancia numerica deterministica 1e-9 eV, nao um
#      limite fisico arbitrario).
# Uso: bash test/test_hscontract_bands.sh   (a partir da raiz do repositorio)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
PROFILE="${PROFILE:-release}"
export OPENBLAS_NUM_THREADS=1
fpm build --profile "$PROFILE" >/dev/null
binpath() { echo "$ROOT/$(fpm run --profile "$PROFILE" "$@" --runner echo 2>/dev/null | tail -n 1)"; }
HSC="$(binpath hscontract)"
LFI="$(binpath LowFoldIn)"
fpm test test_basis_contract --profile "$PROFILE" >/dev/null 2>&1 || true
TESTBIN="$(find build -type f -name test_basis_contract -printf '%T@ %p\n' 2>/dev/null | sort -rn | head -n 1 | cut -d' ' -f2-)"
TESTBIN="$ROOT/$TESTBIN"
if [ -z "$TESTBIN" ]; then echo "  FAIL: test_basis_contract binary not found"; exit 1; fi

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
"$TESTBIN" genbands "$TMP" >/dev/null
cd "$TMP"
pass=0; fail=0
ok()  { pass=$((pass + 1)); }
bad() { fail=$((fail + 1)); echo "  FAIL: $*"; }
check() { if eval "$1"; then ok; else bad "$2"; fi; }

# caminho k identico para os tres calculos: G -> X -> G (42 pontos)
lfi_input() {  # $1 = H file, $2 = S file, $3 = project
    cat <<EOF
&SYSTEM
    ProjectName = '$3'
    Calculation = 'bands'
    nonOrtho = .true.
    HDataFile = '$1'
    SDataFile = '$2'
    OutputDir = '.'
/
&BANDS
3
0.0 0.0 0.0   21
0.5 0.0 0.0   21
0.0 0.0 0.0   1
EOF
}

# 1. contraes
run() { OUT="$("$HSC" "$@" 2>&1)" && RC=0 || RC=$?; }
run contract_id.nml
check "[ $RC -eq 0 ]" "hscontract identity exit 0 (got $RC)"
check "[ -s rid_H.hsio ]" "rid_H.hsio written"
run contract_prune.nml
check "[ $RC -eq 0 ]" "hscontract prune exit 0 (got $RC)"
check "[ -s rpr_H.hsio ]" "rpr_H.hsio written"

# 2. bandas pelo LowFoldIn (fluxo normal) nos tres conjuntos de arquivos
lfi_input model_H.hsio model_S.hsio full  | "$LFI" >/dev/null 2>&1
lfi_input rid_H.hsio  rid_S.hsio  rid  | "$LFI" >/dev/null 2>&1
lfi_input rpr_H.hsio  rpr_S.hsio  rpr  | "$LFI" >/dev/null 2>&1
check "[ -s full.bands.bin ]" "full.bands.bin written by LowFoldIn"
check "[ -s rid.bands.bin ]" "rid.bands.bin written by LowFoldIn"
check "[ -s rpr.bands.bin ]" "rpr.bands.bin written by LowFoldIn"

# 3. T = I: todas as bandas (janela cobre tudo), precisao numerica
OUTID="$("$TESTBIN" bandcmp full.bands.bin rid.bands.bin -100 100 1.0e-9 2>&1)" && RC=0 || RC=$?
check "[ $RC -eq 0 ]" "identity bands equal to numerical precision (got $RC)"
echo "$OUTID"

# 4. poda d -> 0: janela E_F +- 3 eV (apenas bandas s/p), precisao numerica
OUTPR="$("$TESTBIN" bandcmp full.bands.bin rpr.bands.bin -3.0 3.0 1.0e-9 2>&1)" && RC=0 || RC=$?
check "[ $RC -eq 0 ]" "pruned bands equal in the E_F window (got $RC)"
echo "$OUTPR"
check "echo \"\$OUTPR\" | grep -q 'bands in window per k: 5$'" "window holds exactly the 5 s/p bands"
check "echo \"\$OUTPR\" | grep -q 'PASS'" "bandcmp reports PASS"

# 5. contracao radial APROXIMADA s2 -> 1 com pesos [0.9, 0.1] nos DOIS zetas
#    (nao normalizados), combinada com a poda d exata, no modelo com o
#    zeta 2 do s em +5.0 eV: o estado removido (banda ~5 eV) fica fora da
#    janela E_F +- 3 eV e os 4 estados s/p relevantes sao comparados.
#    Limite de regressao DETERMINISTICO deste modelo (0.5 eV, margem ~2x
#    sobre o erro medido max |dE| ~ 0.25 eV): valida ESTE modelo sintetico,
#    nao estabelece tolerancia universal para DFT.
run contract_appr.nml
check "[ $RC -eq 0 ]" "hscontract approximate exit 0 (got $RC)"
check "[ -s rap_H.hsio ]" "rap_H.hsio written"
lfi_input appmodel_H.hsio appmodel_S.hsio full2 | "$LFI" >/dev/null 2>&1
lfi_input rap_H.hsio  rap_S.hsio  rap  | "$LFI" >/dev/null 2>&1
check "[ -s full2.bands.bin ]" "full2.bands.bin written by LowFoldIn"
check "[ -s rap.bands.bin ]" "rap.bands.bin written by LowFoldIn"
OUTAP="$("$TESTBIN" bandcmp full2.bands.bin rap.bands.bin -3.0 3.0 0.5 2>&1)" && RC=0 || RC=$?
check "[ $RC -eq 0 ]" "approximate bands within the model limit (got $RC)"
echo "$OUTAP"
check "echo \"\$OUTAP\" | grep -q 'bands in window per k: 4$'" "window holds the 4 relevant s/p states"
check "echo \"\$OUTAP\" | grep -q 'PASS'" "approximate bandcmp reports PASS"

echo "[test_hscontract_bands] $pass passed, $fail failed"
[ $fail -eq 0 ]
