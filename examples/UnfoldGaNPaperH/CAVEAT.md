# UnfoldGaNPaperH -- validacao com o exemplo do proprio artigo (Fig. 1)

PC cubica de H (a = 3.5 Angstrom, 1 atomo, base H5.0-s1) tiled numa SC de 8
atomos via M = [[2,2,0],[2,-2,0],[0,0,1]] (Eq. 9 do artigo, det = -8). Essa
M e' SIMETRICA (M = M^T), entao serve como primeiro teste de M nao-diagonal
com dado DFT real sem ambiguidade de convencao linha/coluna.

## Como rodar

1. `cd openmx && openmx pc_H.dat && openmx sc_H.dat`
2. Converter os `.scfout` com `openmx2hsio`.
3. `fpm run unfold -- examples/UnfoldGaNPaperH/unfold.nml` (ajuste os
   caminhos de HDataFile/SDataFile no unfold.nml conforme o prefixo que o
   seu openmx2hsio gerou).

## Resultado esperado vs. observado (IMPORTANTE)

Rodando o pipeline completo com um SCF de fato (nao so' a parte geometrica),
os PESOS saem exatamente 0/1 ao longo de todo o caminho Gamma-X (regra de
soma em ~1e-15) -- ou seja, o MOTOR de unfolding (a parte nao-diagonal
generalizada, Sec. II D do artigo) esta' correta.

PORÉM: comparando a energia do peso-1 contra a banda da PC calculada
diretamente (mesma .hsio da PC, via LowFoldIn no modo "bands"), os valores
NAO batem -- e a discrepancia nao e' um deslocamento rigido, ela encolhe de
~13 eV em Gamma para ~1 eV em X. Ao inspecionar os elementos de matriz H(R)
diretamente, o termo on-site H(0,0,0) difere por ~20x entre o calculo da PC
e o da SC (deveria ser IDENTICO -- e' uma quantidade puramente local que nao
depende de como voce descreve a periodicidade). O hopping para o vizinho
mais proximo difere por ~33x. Como as duas razoes sao diferentes, nao e' um
erro de conversao de unidade (isso daria a MESMA razao para os dois); e' o
padrao esperado de dois SCFs que convergiram para densidades eletronicas
diferentes.

Hipotese mais provavel: este e' um sistema fisicamente estranho (H simples
cubico bem compacto, 3.5 Angstrom, raio de corte do PAO ~2.65 Angstrom --
orbitais de vizinhos se sobrepondo bastante), possivelmente com multiplas
solucoes autoconsistentes possiveis, ou que exige configuracoes de SCF mais
cuidadosas (malha de k mais densa, scf.energycutoff maior, scf.criterion
mais apertado, ou uma verificacao manual de que a energia total realmente
parou de oscilar) do que as usadas aqui (que foram escolhidas so' pela
leveza computacional, nao por conta de fisica realista).

CONCLUSAO: os arquivos deste diretorio validam a MATEMATICA do motor de
unfolding (quantizacao 0/1 perfeita, degenerescencia correta na fronteira
da BZ). Eles NAO devem ser tratados, no estado atual, como uma validacao
fisica completa (banda reproduzida corretamente) ate' que o SCF do PC e/ou
da SC seja reconvergido com configuracoes mais cuidadosas e os dois
resultados batam. Se alguem retomar este exemplo, comece checando os logs
de SCF (energia total convergida? quantas iteracoes? algum aviso?) antes de
mexer no codigo.
