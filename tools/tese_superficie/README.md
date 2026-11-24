# Scripts da superfície Te(10-10) plana (tese)

Scripts Python (e um auxiliar Fortran) que geraram os termos de face e o pós-processamento de campo
da superfície semi-infinita de Te(10-10) calculada com o `surfgreen` (modo `extracted` + `&SCC` + `&SITEPOT`).
O código Fortran do LowFoldIn só lê os termos de face por arquivo (`&SITEPOT file = '...sitepot'`,
formato `LFSITEPOT 1`, blocos on-site 29 × 29 em eV na base s3p3d2f1); estes scripts os produzem.
Nenhum dado grande (cubos, `.hsio`) fica no repositório.

Caminhos: o padrão é o layout original da campanha; para outro lugar use as variáveis de ambiente
`LFI_TESE` (raiz `lfi_tese/`, com `vizinhos/`, `campo_semiinf/termo/`, `barreira/a20/`),
`LFI_PAO` (padrão `test/data/pao/Te7.0.pao` deste repositório), `LFI_TDEN` (`TeBulk.tden.cube` do OpenMX),
`LFI_PAOEVAL` (binário `paoeval_bin`) e, no `surf_field.py`, `A8_ROOT` (padrão: diretório corrente),
`A8_RHOMAP` e `A8_VOL`.

## Termo (A): VNA + xc dos vizinhos ausentes (A11 → A12)

| script | faz | entradas | saídas |
|---|---|---|---|
| `a11_pbe.py` | PBE sem spin em malha (V_xc por passo complexo); biblioteca | — | — |
| `a11_vizinhos.py` | deconvolui os potenciais locais do volume (OpenMX) em potenciais atômicos v_at(r) | cubos `TeBulk.{vna,vhart,v0,vxc0}.cube`, `a11_blocos.npz` | `a11_vat_TeBulk.npz`, `a11_saida_*.txt` |
| `a11_r2.py` | ajuste em espaço real f_VNA(r) (B-splines, R_c 6 bohr) e δV_xc PBE; confere com ΔH da DFT | `a11_vat_TeBulk.npz`, `a11_blocos.npz`, `tden` | `a11_r2_fvna.npz`, `a11_r2_saida.txt` |
| `a12_termo.py film` | ΔH(i) = −Σ_aus⟨μ|f_VNA|ν⟩ + ⟨μ|δV_xc|ν⟩ nos 12 átomos centrais do filme L8 (colunas 4–7) | `a11_r2_fvna.npz`, `TeBulk.tden.cube` (fo2), `hs/TeFilm_x1_L8.structure.nml`, `paoeval_bin` | `termo_aus.sitepot` (caso A), `zero.sitepot`, `a12_blocos.npz` |
| `paoeval_bin.f90` | avalia os 29 PAO pelo próprio `pao_eval` do LowFoldIn (`mod_pao_basis`) | `in.bin` (n, posições em bohr), `.pao` | `out.bin` (n × 29) |

`paoeval_bin.f90` não está no `fpm.toml`; compile-o junto dos dois módulos de que depende (na raiz do repositório):
`gfortran -O2 -fopenmp src/mod_types.f90 src/mod_pao_basis.f90 tools/tese_superficie/paoeval_bin.f90 -o paoeval_bin`.
`a12_termo.py check` confere os blocos contra o A11 (átomos 0, 1, 30, 32 do slab N11).

## Termo (B): (A) + barreira FD no vácuo + L = 2 excedente do volume (A18 → A20)

| script | faz | entradas | saídas |
|---|---|---|---|
| `a18_barreira.py V0 d w out [h]` | barreira V_b(y) = V0/(1 + e^{−(y−y0)/w}), y0 = y_topo + d, em blocos on-site (malha h = 0,15 bohr), somada a `termo_aus.sitepot` | `termo_aus.sitepot`, estrutura L8, `paoeval_bin` | `.sitepot` |
| `a20_cf_vol.py` | campo cristalino L = 2 do bloco p1 do volume decomposto em termos uniaxiais por camada de vizinhos; resto na face além do VNA+xc | `a11_blocos.npz`, `termo_aus.sitepot` | `barreira/a20/c_vol/R_cf_{2,3}sh.npy` |
| `a20_sitepot.py V0 d w c out [--tensor arq.npy]` | junta (A) + barreira (chama `a18_barreira.py`) + tensor L = 2 no bloco p1 | idem | `termo.sitepot` |

Caso (B) da tese (`a20_v2_d05_w035_cfvol`): V0 = 2 eV, d = 0,5 Å, w = 0,35 Å, tensor só-volume de 2 camadas
restrito aos átomos 22 e 24 da face (`R_cf_2sh_22_24.npy` = `R_cf_2sh.npy` com as chaves 21 e 23, índices 0-based):

    python3 a20_cf_vol.py
    python3 -c "import numpy as np; d=np.load('c_vol/R_cf_2sh.npy',allow_pickle=True).item(); np.save('c_vol/R_cf_2sh_22_24.npy',{k:d[k] for k in (21,23)},allow_pickle=True)"
    python3 a20_sitepot.py 2.0 0.5 0.35 0 termo.sitepot --tensor c_vol/R_cf_2sh_22_24.npy

## Cadeia do surfgreen (por caso)

1. `surfgreen f3_iii.nml`: SCC (monopolos + dipolos) do filme L8 sobre o volume semi-infinito, com `&SITEPOT`.
2. `surfgreen iii_k4.nml`: one-shot no perfil convergido, `rho_out` na janela [−45, 0] eV (64 pontos de contorno, k 4 × 4) →
   `iii_k4_scc.rho_pairs.dat`.
3. Volume: `surfgreen vol_k4_fo2.nml` (e `vol_k4.nml`) → `vol/vol_k4{,_fo2}.rho_pairs.dat`.

## Campo local nos núcleos (A8)

| script | faz | entradas | saídas |
|---|---|---|---|
| `a8lib.py` | densidade atômica do `.pao` (como o OpenMX), Poisson aberto por G∥, gradientes nos núcleos; biblioteca | — | — |
| `surf_field.py dir prefixo P M` | δρ = ρ − Σρ_at (via `rhomap`), pilha M períodos de volume + sistema, E_loc = +(1/e)∇δV_H por sítio e E₀(j) por camada | `dir/prefixo.rho_pairs.dat`, `vol/$A8_VOL.rho_pairs.dat`, `bin/rhomap` | `dir/E_prefixo_P{P}_M{M}.dat` (E por sítio, V/Å), `.npz`, cubos |
| `vol_check.py P` | volume: E_loc periódico × aberto (validação; ⟨|E_loc|⟩ = 1,358 V/Å) | `vol/volslab_h{P}.cube.*` (rhomap de `vol_k4`) | `vol/E_vol_h{P}.npy` |
| `mkrhomap.py` | escreve um `.nml` do `rhomap` na malha cartesiana usada | — | stdout |

Exemplo (a partir da pasta `campo_semiinf/`, com `bin/rhomap` deste build):

    A8_VOL=vol_k4_fo2 python3 surf_field.py t1vx iii_k4_scc 40 6
    python3 vol_check.py 40
