# LowFoldIn

[![CI](https://github.com/nnardoto/LowFoldIn/actions/workflows/ci.yml/badge.svg)](https://github.com/nnardoto/LowFoldIn/actions/workflows/ci.yml)
[![Licença: GPL v3+](https://img.shields.io/badge/licen%C3%A7a-GPL--3.0--or--later-blue.svg)](LICENSE)
[![Wiki](https://img.shields.io/badge/docs-wiki%20em%20portugu%C3%AAs-2a62c9.svg)](https://github.com/nnardoto/LowFoldIn/wiki)
![Fortran 2023](https://img.shields.io/badge/Fortran-2023-734f96.svg)

**Pós-processamento de estrutura eletrônica em bases localizadas (LCAO), em Fortran moderno.**

Você rodou um cálculo DFT numa supercélula com defeito, liga ou interface, ou num slab, e ficou com dezenas de bandas emaranhadas que não lembram em nada o material de partida. O LowFoldIn pega o Hamiltoniano `H(R)` e a sobreposição `S(R)` do seu cálculo e devolve a física que estava escondida ali:

- 🔁 **Unfolding de bandas:** traz as bandas da supercélula de volta para a zona de Brillouin da célula primitiva, para comparar com o material puro ou com ARPES.
- 🧭 **Análise espectral por região:** mostra em que parte do material (superfície de cima, de baixo, interior, uma camada, uma sub-rede…) cada estado vive, em mapas `A(k,E)` e em mapas `kx × ky` de energia constante.
- 📈 **Estrutura de bandas convencional**, com o problema de autovalores generalizado `H(k)c = εS(k)c` para bases não ortogonais.

Os dados entram direto do **OpenMX** (`.scfout`), por um conversor nativo. Não é preciso refazer o cálculo autoconsistente: tudo parte das matrizes que você já tem.

<p align="center">
  <img src="docs/images/unfold_conceito_2x2.png" width="820" alt="Supercélula 2x2 de grafeno: bandas dobradas e resultado do unfolding">
  <br><sub>Supercélula 2×2 de grafeno: as 8 bandas dobradas (esquerda) e o unfolding recuperando as bandas da célula primitiva (direita).</sub>
</p>

<p align="center">
  <a href="#-começando">Começando</a> ·
  <a href="#-como-funciona">Como funciona</a> ·
  <a href="#-as-ferramentas">Ferramentas</a> ·
  <a href="#-exemplos-incluídos">Exemplos</a> ·
  <a href="#-documentação">Documentação</a> ·
  <a href="#-contribuindo">Contribuindo</a>
</p>

---

## ✨ Destaques

- **Unfolding com qualquer supercélula:** a matriz `M` pode ser diagonal (rotina fechada, custo O(m)) ou geral, girada ou cisalhada (projetor de grupo, custo O(m²)). Por exemplo, √3×√3 R30°.
- **Mapeamento automático de átomos** entre a célula primitiva e a supercélula, a partir das posições.
- **Bases não ortogonais** de ponta a ponta, com ortogonalização de Löwdin e pesos de Mulliken ou Löwdin.
- **Spin:** sem spin, colinear e não colinear (espinores).
- **Separação entre cálculo e visualização:** o `spectral` diagonaliza uma única vez; o `specmap` refaz alargamentos, janelas de energia e combinações de regiões em segundos, sem tocar no Hamiltoniano.
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

> [!NOTE]
> O `pip install fpm` instala a versão 0.12, que é antiga demais. Use o binário oficial disponível nas [releases do fpm](https://github.com/fortran-lang/fpm/releases).

### Compilar e testar

```bash
git clone https://github.com/nnardoto/LowFoldIn.git
cd LowFoldIn
fpm build --profile release     # -O3 + OpenMP
fpm test  --profile release     # suíte de testes
```

Além do `release`, o [`fpm.toml`](fpm.toml) traz outros perfis prontos:

| Perfil | Quando usar |
|---|---|
| `release` | Produção portátil (`-O3`, OpenMP) — o padrão recomendado |
| `native` | Produção otimizada para a CPU da máquina (`-march=native`) |
| `development` | Dia a dia de desenvolvimento (`-O1`, símbolos de depuração) |
| `debug` | Caça a bugs: `-fcheck=all`, armadilhas de ponto flutuante, inicialização com NaN |
| `debug-sanitize` | Como `debug`, mais AddressSanitizer e UBSan |

### Ver tudo funcionando em um minuto

```bash
bash examples/surface_spectral/run_example.sh
```

Esse script gera um slab sintético com estados de superfície, roda o fluxo completo de análise por região e produz os mapas em PNG/PDF. É o jeito mais rápido de ver as ferramentas conversando entre si.

---

## 🧠 Como funciona

**Unfolding.** Cada autoestado $|KN\rangle$ da supercélula é levado à base ortogonal de Löwdin, $C'_{KN} = S(K)^{1/2}\,C_{KN}$, e projetado nos autoestados $F_{k,n}$ do operador de translação da célula primitiva. O peso espectral no ponto $k$ da célula primitiva é

$$
W^{k}_{KN} \;=\; \sum_{n=1}^{p} \big|\, F_{k,n}^{\dagger}\, C'_{KN} \big|^{2},
$$

onde a soma percorre os $p$ orbitais distintos da célula primitiva. A função espectral desdobrada é $A(k,E) = \sum_N W^{k}_{KN}\,\delta(E-\varepsilon_{KN})$, com a delta alargada por uma lorentziana ou gaussiana. Como $\sum_N W^{k}_{KN}$ deve ser igual ao número de orbitais da célula primitiva, essa **regra de soma** serve de verificação automática em toda execução.

**Análise por região.** Para uma região $P$ (um conjunto de átomos), o peso $w_n^P(k)$ é a fração do estado $n$ que vive em $P$ (Mulliken ou Löwdin; soma 1 sobre uma partição). O mapa projetado é

$$
A_P(k,E) \;=\; \sum_n w_n^P(k)\; K\big(E-\varepsilon_n(k)\big),
$$

com $K$ um núcleo de alargamento escolhido por você. Combinações como `top-bottom` destacam em vermelho e azul onde cada estado de superfície mora.

Os detalhes, inclusive o mapeamento de pontos k para `M` não diagonal, estão em [Metodologia do unfolding](https://github.com/nnardoto/LowFoldIn/wiki/Metodologia-do-Unfolding) e [Metodologia da análise espectral](https://github.com/nnardoto/LowFoldIn/wiki/Metodologia-da-Analise-Espectral).

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

Todas leem a entrada em **namelists Fortran**: arquivos de texto curtos, comentáveis e fáceis de versionar junto com o seu projeto.

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

### Outro gostinho: quem mora na superfície?

Primeiro, nomeie as regiões a partir da geometria (os comandos abaixo supõem os binários no `PATH`; sem isso, use `fpm run atomselect --profile release -- <argumentos>`):

```bash
atomselect surface slab.nml --regions regions.nml --normal 0 0 1 --top 1 --bottom 1
atomselect select  slab.nml --regions regions.nml --name interior --atoms "2-7"
atomselect validate slab.nml regions.nml       # confere se as regiões formam uma partição
```

Depois, diagonalize uma vez e gere quantos mapas quiser:

```fortran
&SPECTRAL                                 ! spectral.nml
    h_file = 'slab_H.hsio',  s_file = 'slab_S.hsio'
    structure_file = 'slab.nml',  regions_file = 'regions.nml'
    method = 'mulliken',  output_prefix = 'slab'
    energy_reference = 'fermi'            ! openmx2hsio já desloca E_F para 0
/
&KPATH
    nvertex = 3
    vertex(:,1) = 0.5, 0.0, 0.0
    vertex(:,2) = 0.0, 0.0, 0.0
    vertex(:,3) = 0.5, 0.5, 0.0
    label   = 'X', 'G', 'M'
    npoints = 60, 85
/
```

```fortran
&SPECMAP                                  ! specmap.nml — rode de novo à vontade
    spectral_file = 'slab.spectral.bin'
    kernel = 'lorentzian',  width_ev = 0.05
    emin = -3.5, emax = 3.5, ne = 701
    maps = 'total', 'top', 'bottom', 'top-bottom'
    output_prefix = 'slab'
/
```

Trocar o alargamento, a janela de energia ou as combinações de regiões é só editar o `specmap.nml`: nada é rediagonalizado. Trocando `&KPATH` por `&KGRID`, o mesmo fluxo produz mapas de energia constante e contornos de Fermi.

---

## 🗂️ Exemplos incluídos

| Pasta | O que mostra |
|---|---|
| [`examples/surface_spectral`](examples/surface_spectral) | Fluxo completo `atomselect → spectral → specmap` num slab sintético com estados de superfície |
| [`examples/UnfoldChain1D`](examples/UnfoldChain1D) | O caso mais simples possível: unfolding de uma cadeia 1D |
| [`examples/UnfoldChain1D_Auto`](examples/UnfoldChain1D_Auto) | A mesma cadeia, com mapeamento automático de átomos |
| [`examples/UnfoldChain1DShearedM`](examples/UnfoldChain1DShearedM) | Matriz `M` não diagonal (supercélula cisalhada) |
| [`examples/UnfoldDiatomic`](examples/UnfoldDiatomic) | Base com mais de um átomo por célula primitiva |
| [`examples/UnfoldSpinorNonCollinear`](examples/UnfoldSpinorNonCollinear) | Spin não colinear (espinores) |
| [`examples/UnfoldGaNPaperH`](examples/UnfoldGaNPaperH) | Validação contra o exemplo DFT de GaN do artigo de referência |

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

> [!IMPORTANT]
> Os mapas espectrais são **pesos de projeção alargados**, não simulações de ARPES: não incluem elementos de matriz de fotoemissão, estados finais nem efeitos de muitos corpos.

---

## 🚧 Estado do projeto

Software de pesquisa em desenvolvimento ativo. As ferramentas foram validadas em modelos sintéticos, em pequenos casos de teste e no exemplo DFT do artigo de referência. Antes de publicar, confira sempre a regra de soma e as normalizações. As limitações atuais (por exemplo: só OpenMX 3.9/4.0 e coordenadas cartesianas em Å) estão listadas na [wiki](https://github.com/nnardoto/LowFoldIn/wiki/Limitacoes).

## 🤝 Contribuindo

Issues e pull requests são muito bem-vindos, inclusive dúvidas de uso e relatos do tipo "rodei no meu sistema e deu isto". Uma ótima porta de entrada é um conversor para outro código de base local (SIESTA, FHI-aims, CP2K…): basta produzir `H(R)` e `S(R)` no formato `.hsio`, e todo o resto funciona.

Em mudanças numéricas, informe:

- compilador e versão, implementação de BLAS/LAPACK e número de threads;
- uma entrada mínima que reproduza o resultado;
- os desvios da regra de soma esperados e obtidos.

Antes de abrir o PR, rode `fpm test --profile debug`: o perfil de depuração liga as checagens de limites e as armadilhas de ponto flutuante.

## 📄 Licença

Software livre sob a **GNU General Public License, versão 3 ou posterior** (`GPL-3.0-or-later`). Veja [`LICENSE`](LICENSE).

**Autor:** Nathanael Nardoto Batista
