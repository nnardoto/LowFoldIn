# LowFoldIn

[![CI](https://github.com/nnardoto/LowFoldIn/actions/workflows/ci.yml/badge.svg)](https://github.com/nnardoto/LowFoldIn/actions/workflows/ci.yml)
[![Licença: GPL v3+](https://img.shields.io/badge/licen%C3%A7a-GPL--3.0--or--later-blue.svg)](LICENSE)
[![Wiki](https://img.shields.io/badge/docs-wiki%20em%20portugu%C3%AAs-2a62c9.svg)](https://github.com/nnardoto/LowFoldIn/wiki)

**Pós-processamento de estrutura eletrônica em bases localizadas (LCAO), em Fortran moderno.**

Você rodou um cálculo DFT numa supercélula com defeito, liga ou interface, ou num slab, e agora precisa enxergar a física por trás dele. O LowFoldIn pega o Hamiltoniano `H(R)` e a sobreposição `S(R)` do seu cálculo e oferece:

- 🔁 **Unfolding de bandas:** traz as bandas da supercélula de volta para a zona de Brillouin da célula primitiva, para comparar com o material puro ou com ARPES.
- 🧭 **Análise espectral por região:** mostra em que parte do material (superfície de cima, de baixo, interior, uma camada, uma sub-rede…) cada estado vive, em mapas `A(k,E)` e em mapas `kx × ky` de energia constante.
- 📈 **Estrutura de bandas convencional**, com o problema de autovalores generalizado `H(k)c = εS(k)c` para bases não ortogonais.

Os dados entram direto do **OpenMX** (`.scfout`), por um conversor nativo.

<p align="center">
  <img src="docs/images/unfold_conceito_2x2.png" width="820" alt="Supercélula 2x2 de grafeno: bandas dobradas e resultado do unfolding">
  <br><sub>Supercélula 2×2 de grafeno: as 8 bandas dobradas (esquerda) e o unfolding recuperando as bandas da célula primitiva (direita).</sub>
</p>

---

## ✨ Destaques

- **Unfolding com qualquer supercélula:** a matriz `M` pode ser diagonal (rotina fechada, custo O(m)) ou geral, girada ou cisalhada (projetor de grupo, custo O(m²)). Por exemplo, √3×√3 R30°.
- **Mapeamento automático de átomos** entre a célula primitiva e a supercélula, a partir das posições.
- **Bases não ortogonais** de ponta a ponta, com ortogonalização de Löwdin e pesos de Mulliken ou Löwdin.
- **Spin:** sem spin, colinear e não colinear (espinores).
- **Separação entre cálculo e visualização:** o `spectral` diagonaliza uma única vez; o `specmap` refaz alargamentos, janelas de energia e combinações de regiões em segundos.
- **Paralelismo OpenMP** sobre os pontos k, com resultados idênticos bit a bit para qualquer número de threads.
- **Verificações embutidas** em cada execução: regra de soma do unfolding e normalização dos pesos.
- **Testado:** 12 suítes de testes unitários, mais testes de linha de comando e dos exemplos, rodando no CI.

<table>
<tr>
<td width="50%"><img src="docs/images/unfold_impureza.png" alt="Função espectral desdobrada com impureza"></td>
<td width="50%"><img src="docs/images/grhbn_camadas.png" alt="Mapa grafeno menos hBN e gap em K"></td>
</tr>
<tr>
<td><sub>Grafeno 3×3 com uma impureza: função espectral desdobrada, com gaps de hibridização e bandas-fantasma.</sub></td>
<td><sub>Grafeno sobre hBN: mapa por camada (vermelho = grafeno, azul = hBN) e o gap de 68 meV aberto em K.</sub></td>
</tr>
</table>

<sub>Figuras geradas com modelos tight-binding sintéticos (não é DFT). As entradas para reproduzi-las estão na wiki, em [Modelos Sintéticos](https://github.com/nnardoto/LowFoldIn/wiki/Modelos-Sinteticos).</sub>

---

## 🚀 Começando

### Requisitos

- GNU Fortran **≥ 14** (o código usa `-std=f2023`)
- [fpm](https://fpm.fortran-lang.org/) **≥ 0.13**
- BLAS e LAPACK (OpenBLAS, Accelerate, MKL…)
- OpenMP (opcional, recomendado) e gnuplot (opcional, para os gráficos)

```bash
# Debian/Ubuntu
sudo apt install gfortran-14 liblapack-dev libblas-dev
# Arch Linux
sudo pacman -S gcc-fortran blas lapack
# macOS (Homebrew; LAPACK vem do Accelerate)
brew install gcc fpm gnuplot
```

> O `pip install fpm` instala a versão 0.12, que é antiga demais. Use o binário oficial disponível nas [releases do fpm](https://github.com/fortran-lang/fpm/releases).

### Compilar e testar

```bash
git clone https://github.com/nnardoto/LowFoldIn.git
cd LowFoldIn
fpm build --profile release     # -O3 + OpenMP
fpm test  --profile release     # suíte de testes
```

### Ver tudo funcionando em um minuto

```bash
bash examples/surface_spectral/run_example.sh
```

Esse script roda o fluxo completo de análise de superfícies num slab sintético e gera os mapas em PNG/PDF.

---

## 🧰 As ferramentas

| Programa | Para que serve |
|---|---|
| `openmx2hsio` | Converte o `.scfout` do OpenMX nas matrizes `H(R)`/`S(R)` (formato `.hsio`) e, com `--write-structure`, extrai a estrutura |
| `LowFoldIn` | Estrutura de bandas ao longo de um caminho k |
| `unfold` | Unfolding das bandas da supercélula para a célula primitiva |
| `atomselect` | Define regiões atômicas nomeadas: por índices, camadas, superfícies, caixas… |
| `spectral` | Energias e pesos de projeção por região, num caminho k ou numa grade 2D |
| `specmap` | Mapas espectrais alargados, combinações como `top-bottom` e mapas de energia constante |

### Fluxos de trabalho

```text
OpenMX .scfout ──► openmx2hsio ──► H.hsio, S.hsio (+ estrutura)
                                        │
             ┌──────────────────────────┼──────────────────────────────┐
             ▼                          ▼                              ▼
        LowFoldIn                    unfold                  atomselect → spectral → specmap
     (bandas comuns)       (bandas desdobradas, A(k,E))      (mapas por região, contornos de Fermi)
```

### Um gostinho: unfolding

```fortran
&UNFOLD
    ProjectName = "sc"
    HDataFile   = "sc_H.hsio"
    SDataFile   = "sc_S.hsio"
    OutputDir   = "resultados"
    M(1,:) = 2, 0, 0            ! linhas de M: A_i = Σ_j M_ij a_j
    M(2,:) = 0, 2, 0
    M(3,:) = 0, 0, 1
    PCStructureFile = "pc.nml"  ! mapeamento automático de átomos
    SCStructureFile = "sc.structure.nml"
    WriteSpectralFunction = .true.
    Broadening = 0.05
    Emin = -5.0, Emax = 5.0, nE = 1000
/
&PCPATH
4
0.000000  0.000000  0.0   100   ! Γ
0.500000  0.000000  0.0    58   ! M
0.666667  0.333333  0.0   115   ! K
0.000000  0.000000  0.0     1   ! Γ
/
```

```bash
mkdir -p resultados
OMP_NUM_THREADS=8 fpm run unfold --profile release -- unfold.nml
```

Ao final, o programa informa o desvio da **regra de soma** (Σ W = número de orbitais da célula primitiva). Esse valor deve ficar na precisão de máquina.

---

## 📚 Documentação

A documentação completa, **em português**, está na **[wiki do projeto](https://github.com/nnardoto/LowFoldIn/wiki)**:

| Comece por | Aprofunde em | Consulte |
|---|---|---|
| [Instalação](https://github.com/nnardoto/LowFoldIn/wiki/Instalacao-e-Compilacao) | [Fundamentos teóricos](https://github.com/nnardoto/LowFoldIn/wiki/Fundamentos-Teoricos) | [Formatos de arquivo](https://github.com/nnardoto/LowFoldIn/wiki/Formatos-de-Arquivo) |
| [Guia rápido](https://github.com/nnardoto/LowFoldIn/wiki/Guia-Rapido) | [Metodologia do unfolding](https://github.com/nnardoto/LowFoldIn/wiki/Metodologia-do-Unfolding) | [Solução de problemas](https://github.com/nnardoto/LowFoldIn/wiki/Solucao-de-Problemas) |
| [Modelos sintéticos (grafeno, hBN)](https://github.com/nnardoto/LowFoldIn/wiki/Modelos-Sinteticos) | [Metodologia da análise espectral](https://github.com/nnardoto/LowFoldIn/wiki/Metodologia-da-Analise-Espectral) | [Desempenho e paralelismo](https://github.com/nnardoto/LowFoldIn/wiki/Desempenho-e-Paralelismo) |

A referência detalhada de cada ferramenta de superfície, em inglês, fica em [`docs/`](docs): [`atomselect`](docs/atomselect.md), [`spectral`](docs/spectral.md), [`specmap`](docs/specmap.md) e [desempenho](docs/performance.md).

---

## 🔬 Base científica

O unfolding segue o formalismo para orbitais centrados em átomos de:

> J. Quan, N. Rybin, M. Scheffler e C. Carbogno,
> *Efficient band structure unfolding with atom-centered orbitals: General theory and application*,
> **Phys. Rev. B 113**, 085112 (2026). DOI: [10.1103/7xym-7388](https://doi.org/10.1103/7xym-7388)

Se você usar o unfolding, cite esse trabalho. Uma citação própria do LowFoldIn (com DOI) será adicionada quando houver uma release arquivada.

> Os mapas espectrais são **pesos de projeção alargados**, não simulações de ARPES: não incluem elementos de matriz de fotoemissão, estados finais nem efeitos de muitos corpos.

---

## 🚧 Estado do projeto

Software de pesquisa em desenvolvimento ativo. As ferramentas foram validadas em modelos sintéticos, em pequenos casos de teste e no exemplo DFT do artigo de referência. Antes de publicar, confira sempre a regra de soma e as normalizações. As limitações atuais (por exemplo: só OpenMX 3.9/4.0 e coordenadas cartesianas em Å) estão listadas na [wiki](https://github.com/nnardoto/LowFoldIn/wiki/Limitacoes).

## 🤝 Contribuindo

Issues e pull requests são muito bem-vindos, e uma boa porta de entrada é um conversor para outro código de base local (SIESTA, FHI-aims, CP2K…). Em mudanças numéricas, informe:

- compilador e versão, implementação de BLAS/LAPACK e número de threads;
- uma entrada mínima que reproduza o resultado;
- os desvios da regra de soma esperados e obtidos.

## 📄 Licença

Software livre sob a **GNU General Public License, versão 3 ou posterior** (`GPL-3.0-or-later`). Veja [`LICENSE`](LICENSE).

**Autor:** Nathanael Nardoto Batista
