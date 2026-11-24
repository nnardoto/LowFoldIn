# Trincheira de Te(10-10) só a partir do volume (apêndice da tese)

Cadeia removida numa região central de 3 células, acoplada ao volume semi-infinito abaixo, com vácuo acima e periódica em x
(11 colunas, 52 Å, k_x = Γ). Usa o mesmo SCC com dipolos e os termos de face de ../README.md. A barreira de Fermi–Dirac
plana é generalizada por uma distância assinada suave, s = d_λ(presentes) − d_λ(ausentes), com
d_λ = −λ ln Σ_j exp(−|r − R_j|/λ) e λ = 2,5 Å. (s₀, w_s) = (−0,771 Å, 0,396 Å) são calibrados uma vez nos blocos da
barreira plana (V₀ = 2 eV, d = 0,5 Å, w = 0,35 Å; resíduo 6 %).

| arquivo | o que faz |
|---|---|
| a25_face.py | termos de face da trincheira (VNA+xc dos vizinhos ausentes + barreira generalizada) → .sitepot; calibração de (s₀, w_s) |
| gen.slurm | job do scicom que roda a25_face.py |
| entradas/tr*.nml | SCC e mapa espectral da trincheira (surfgreen, modo extracted) |
| entradas/occ.nml, emp.nml, map_*.nml | densidades nas janelas [E_F−0,5; E_F] e [E_F; E_F+0,5] eV (one-shot) e rhomap 'integrate' |
| a25_tr_post.py | espectros e carga induzida por cadeia |
| fig_trincheira_vol_bordas.py | figura DFT × LowFoldIn dos estados em espaço real |

Caminhos: variáveis de ambiente A25_ROOT, LFI_TESE e TESE_FIGSCRIPTS (padrão: os caminhos do autor). Os dados grandes
(.hsio, cubos, densidades) não estão no repositório.
