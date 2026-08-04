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

## Resultado (validado)

Regra de soma em ~1e-15, pesos exatamente 0/1 ao longo de todo o caminho
Gamma-X. Comparando a energia do peso-1 contra a banda da PC calculada
diretamente (mesma .hsio da PC, via LowFoldIn no modo "bands"), o maior
desvio entre os 41 pontos do caminho e' de ~6e-5 eV -- ruido numerico de
formatacao em arquivo texto, nao discrepancia fisica.

Isso valida, com dado DFT real e independente (nao so' um exemplo
sintetico), tanto a matematica do motor de unfolding nao-diagonal (Sec.
II D do artigo) quanto a fisica (banda recuperada por unfolding == banda
direta da PC).

## Nota historica

Uma primeira rodada deste exemplo (guardada so' na conversa, nao neste
repositorio) tinha um erro na configuracao/execucao do calculo da PC que
inflava artificialmente o overlap S entre vizinhos mais proximos (0.351 em
vez do 0.00616 correto, facilmente identificavel comparando S da PC com o
S entre os mesmos pares fisicos de atomos na SC -- overlap e' puramente
geometrico, entao os dois calculos TEM que concordar exatamente pra
qualquer distancia real igual). Fica registrado aqui como lembrete: ao
depurar uma discrepancia PC-vs-SC, comparar o overlap (que nao depende de
autoconsistencia nem nivel de Fermi) e' um jeito rapido de isolar se o
problema esta' na geometria/base ou em outro lugar (SCF, Fermi shift,
etc.) antes de suspeitar do codigo de unfolding em si.

