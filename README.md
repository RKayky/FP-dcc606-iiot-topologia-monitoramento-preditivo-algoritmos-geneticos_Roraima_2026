[README.md](https://github.com/user-attachments/files/29226652/README.md)
<div align="center">

<img src="https://img.shields.io/badge/UFRR-Ciência%20da%20Computação-065A82?style=for-the-badge&logo=university&logoColor=white"/>
<img src="https://img.shields.io/badge/DCC606-Análise%20de%20Algoritmos-02C39A?style=for-the-badge"/>
<img src="https://img.shields.io/badge/2026-Projeto%20Final-1C7293?style=for-the-badge"/>

# Otimização Topológica de Redes IIoT
### para Monitoramento Preditivo em Ambientes Industriais
#### Uma Abordagem Baseada em Algoritmos Genéticos

[![Linguagem](https://img.shields.io/badge/Linguagem-C99-A8B9CC?style=flat-square&logo=c&logoColor=white)](https://en.wikipedia.org/wiki/C99)
[![Visualização](https://img.shields.io/badge/Visualização-Python%203-3776AB?style=flat-square&logo=python&logoColor=white)](https://www.python.org/)
[![Licença](https://img.shields.io/badge/Licença-MIT-green?style=flat-square)](LICENSE)
[![Tema](https://img.shields.io/badge/Tema-5%20%2F%20DCC606-orange?style=flat-square)]()
[![Status](https://img.shields.io/badge/Status-Completo-brightgreen?style=flat-square)]()

</div>

---

## Índice

- [Visão Geral](#-visão-geral)
- [Problema Abordado](#-problema-abordado)
- [Formulação Matemática](#-formulação-matemática)
- [Algoritmo Implementado](#-algoritmo-implementado)
- [Estrutura do Repositório](#-estrutura-do-repositório)
- [Dependências e Compilação](#-dependências-e-compilação)
- [Como Executar](#-como-executar)
- [Visualização dos Resultados](#-visualização-dos-resultados)
- [Análise de Complexidade](#-análise-de-complexidade)
- [Referências Bibliográficas](#-referências-bibliográficas)
- [Autores](#-autores)

---

## Visão Geral

No contexto da **Indústria 4.0**, redes de sensores IIoT (*Industrial Internet of Things*) são infraestrutura crítica para monitoramento preditivo de máquinas, automação de processos e coleta de dados em tempo real. Centenas de nós sensores — tipicamente baseados em microcontroladores de baixo consumo como o **ESP32** — precisam transmitir dados continuamente a gateways coordenadores.

O **planejamento topológico** dessas redes é um problema de **otimização combinatória NP-difícil**: estruturas metálicas e paredes de concreto criam zonas de sombra que bloqueiam sinais, sensores operam com bateria (eficiência energética é crítica) e o número de posições candidatas a gateway gera um espaço de busca exponencial de $2^m$ subconjuntos.

Este projeto implementa um **motor evolutivo completo em C99** baseado em **Algoritmos Genéticos** para sintetizar automaticamente a topologia ótima de uma rede IIoT industrial, determinando **quais** gateways candidatos ativar e **onde** posicioná-los em coordenadas tridimensionais, maximizando simultaneamente vida útil, cobertura e redundância de caminhos.

---

## Problema Abordado

O problema é classificado como variante NP-difícil do **Facility Location Problem** combinado com o **Maximum Coverage Problem** [Li et al., 2009]. A determinação exata do *Minimum Connected Dominating Set* (MCDS) de gateways que garante cobertura total é NP-difícil, inviabilizando busca exaustiva para instâncias industriais.

| Parâmetro | Valor na Instância Padrão |
|---|---|
| Dimensões da planta | 100 × 100 × 10 m |
| Sensores fixos (*n*) | 20 (grade 5×4, *z* = 2,0 m) |
| Candidatos a gateway (*m*) | 15 posições |
| Orçamento máximo de gateways | 5 |
| Obstáculos industriais | 4 (2 estruturas metálicas + 2 paredes) |
| Potência TX (ESP32) | 20 dBm |
| Limiar mínimo de RSSI | −90 dBm |
| Expoente de atenuação (γ) | 2,7 (ambiente industrial interno) |
| Capacidade da bateria | 5.000 mAh |
| Vida útil máxima estimada | 8.760 h (1 ano) |

**Espaço de busca:** $2^{15} = 32.768$ subconjuntos binários de ativação × espaço contínuo $\mathbb{R}^3$ de coordenadas → inviável por força bruta.

---

## Formulação Matemática

### Cromossomo Misto

Cada indivíduo na população do AG codifica duas classes de variáveis:

```
x = ( a₁, a₂, …, a₁₅ ,  (x₁,y₁,z₁), (x₂,y₂,z₂), …, (x₁₅,y₁₅,z₁₅) )
      ←── parte binária ──→  ←──────────── parte real ───────────────→

aⱼ ∈ {0,1}        →  ativar gateway candidato j?
(xⱼ, yⱼ, zⱼ) ∈ ℝ³  →  coordenadas refinadas do gateway j
```

Tamanho por cromossomo: `15×4 + 15×3×8 + 8 = 428 bytes`  
Memória total (2 populações × 120 indivíduos): ≈ **102,7 KB**

### Modelo de Propagação — Equação Log-Distância

$$P_{\text{loss}}(d) = P_{\text{loss}}(d_0) + 10 \cdot \gamma \cdot \log_{10}\!\left(\frac{d}{d_0}\right) + X_\sigma$$

$$\text{RSSI}(s_i, g_j) = P_{TX} - P_{\text{loss}}(d) - A_{\text{obs}}(s_i, g_j)$$

A atenuação adicional $A_{\text{obs}}$ é calculada pelo **algoritmo de Liang-Barsky** (O(1) por obstáculo), que detecta interseções entre o segmento sensor→gateway e os retângulos dos obstáculos industriais.

### Função de Aptidão Multiobjetivo

$$\max\ F(\mathbf{x}) = \underbrace{0{,}50 \cdot \hat{L}}_{\text{vida útil}} - \underbrace{0{,}30 \cdot \hat{D}}_{\text{latência}} + \underbrace{0{,}20 \cdot \hat{R}}_{\text{redundância}} - \underbrace{2{,}0 \cdot (1-\rho_{\text{cob}})}_{\phi_{\text{cob}}} - \underbrace{0{,}5 \cdot \max(0,|G^*|-5)}_{\phi_{\text{bud}}}$$

### Restrições

| # | Restrição | Tipo | Tratamento |
|---|---|---|---|
| C1 | $\forall s_i \in S,\ \exists g_j \in G^*: \text{RSSI}(s_i,g_j) \geq -90\ \text{dBm}$ | Hard | Penalidade $\phi_{\text{cob}}$ |
| C2 | $\|G^*\| \leq M_{\max} = 5$ | Hard | `repair_budget()` determinístico |
| C3 | $(x_j, y_j) \in [0,100]^2$, $z_j \in [2,10]$ m | Domínio | `clamp_d()` determinístico |

---

## Algoritmo Implementado

### Fluxo Geracional

```
┌─────────────────────────────────────────────────────────────┐
│                    CICLO EVOLUTIVO                          │
│                                                             │
│  [Inicialização]  →  [Avaliação Eq.5]  →  [Ordenação]      │
│        ↓                                                    │
│  [Elitismo k=6]  →  [Seleção]  →  [Cruzamento]             │
│                          ↓              ↓                   │
│                    [Torneio O(1)]  [1pt / 2pt aleatório]    │
│                    [Roleta  O(N)]       ↓                   │
│                                   [Mutação mista]           │
│                                    bit-flip + N(0,σ)        │
│                                         ↓                   │
│                                   [Reparo budget]           │
│                                         ↓                   │
│                            [Nova Geração]                   │
│        ↑_______________________________________↓            │
│                                                             │
│  [g = MAX_GENS?] → SIM → [Exportar CSV + Mapa ASCII]        │
└─────────────────────────────────────────────────────────────┘
```

### Parâmetros do AG

| Parâmetro | Símbolo | Valor |
|---|---|---|
| Tamanho da população | $N$ | 120 |
| Número máximo de gerações | $G_{\max}$ | 600 |
| Elitismo | $k$ | 6 |
| Tamanho do torneio | $k_T$ | 5 |
| Probabilidade de cruzamento | $p_c$ | 0,85 |
| Taxa de mutação binária | $p_{mb}$ | 0,04 |
| Taxa de mutação real | $p_{mr}$ | 0,10 |
| Desvio padrão gaussiano XY | $\sigma_{XY}$ | 5,0 m |
| Desvio padrão gaussiano Z | $\sigma_Z$ | 1,5 m |

### Operadores de Seleção Comparados

| Operador | Complexidade | Pressão | Característica |
|---|---|---|---|
| **Torneio Estocástico** ($k_T=5$) | $O(1)$ por seleção | Moderada e controlável | Independente da escala de aptidão. Maior diversidade. |
| **Roleta Viciada por Aptidão** | $O(N)$ por seleção | Proporcional à aptidão | Convergência inicial rápida. Risco de *crowding*. |

Experimentos com **sementes independentes** ($\Delta_{\text{seed}} = 7919$, primo) garantem comparação estatisticamente válida.

---

##  Estrutura do Repositório

```
FP_DCC606_Tema_5_RR_2026/
│
├── README.md                          ← Este arquivo
│
├── src/
│   └── iiot_ga_1_.c                     ← Código-fonte principal (C99)
│
├── analysis/
│   └── plot_resultados.py               ← Gerador de 4 gráficos de análise
│
├── results/                          ← Saídas geradas pelo programa
│   ├── hist_torneio.csv                 ← Histórico por geração — Torneio
│   └── hist_roleta.csv                  ← Histórico por geração — Roleta
│
├── docs/
│   ├── relatorio_DCC606_Tema5_RR_2026.md   ← Artigo científico (formato IEEE)
│   └── DCC606_Tema5_Apresentacao.pptx      ← Slides do seminário
│
└── LICENSE
```

---

## Dependências e Compilação

### Requisitos

| Ferramenta | Versão mínima | Uso |
|---|---|---|
| `gcc` | 7.0+ | Compilar o código C99 |
| `make` | qualquer | Atalho de compilação (opcional) |
| `python3` | 3.8+ | Gerador de gráficos |
| `pandas` | 1.0+ | Leitura dos CSVs |
| `matplotlib` | 3.3+ | Plotagem dos gráficos |
| `numpy` | 1.18+ | Suavização de curvas |

### Instalar dependências Python

```bash
pip install pandas matplotlib numpy
```

### Compilar o programa C

```bash
# Compilação padrão com otimização
gcc -O2 src/iiot_ga_1_.c -lm -o iiot_ga

# Compilação com informações de debug
gcc -O0 -g src/iiot_ga_1_.c -lm -o iiot_ga_debug

# Verificação de warnings (recomendada antes da submissão)
gcc -O2 -Wall -Wextra src/iiot_ga_1_.c -lm -o iiot_ga
```

> **Nota:** o programa utiliza apenas as bibliotecas padrão C (`stdio.h`, `stdlib.h`, `math.h`, `string.h`, `time.h`). Nenhuma biblioteca externa de otimização é necessária.

---

## Como Executar

### 1. Executar o motor evolutivo

```bash
./iiot_ga
```

O programa executa **dois experimentos consecutivos** (Torneio e Roleta) e produz automaticamente:

```
╔══════════════════════════════════════════════════╗
║  Experimento 1: Seleção por Torneio Estocástico  ║
╚══════════════════════════════════════════════════╝
Geração   1 | Melhor: X.XXXXX | Média: X.XXXXX
Geração  50 | Melhor: X.XXXXX | Média: X.XXXXX
...
Geração 600 | Melhor: X.XXXXX | Média: X.XXXXX

[Mapa ASCII da planta — 54×27 caracteres]
[Tabela de RSSI por par sensor-gateway]
[Métricas consolidadas da melhor topologia]

[Arquivo salvo: hist_torneio.csv]

╔══════════════════════════════════════════════════╗
║  Experimento 2: Seleção por Roleta Viciada       ║
╚══════════════════════════════════════════════════╝
...
[Arquivo salvo: hist_roleta.csv]
```

#### Legenda do mapa ASCII

| Símbolo | Significado |
|---|---|
| `[G]` | Gateway ativo na melhor solução |
| `[s]` | Sensor fixo |
| `[#]` | Obstáculo industrial |
| `.` | Área livre da planta |

### 2. Gerar os gráficos de análise

```bash
# Certifique-se de que hist_torneio.csv e hist_roleta.csv estão no mesmo diretório
python3 analysis/plot_resultados.py
```

O script gera e salva automaticamente **4 gráficos**:

| Arquivo | Conteúdo |
|---|---|
| `g1_convergencia.png` | Curva de convergência — melhor aptidão por geração |
| `g2_aptidao_media.png` | Evolução da aptidão média da população |
| `g3_comparativo.png` | Comparativo sobrepostos (torneio vs roleta) |
| `g4_ganho_marginal.png` | Ganho marginal ΔF por geração (janela de 30 gerações) |

E imprime no terminal uma **tabela de resumo comparativo**:

```
════════════════════════════════════════════════════════
  RESUMO COMPARATIVO DOS EXPERIMENTOS
════════════════════════════════════════════════════════
  Métrica                        Torneio     Roleta
  ────────────────────────────────────────────────────
  Melhor aptidão final           X.XXXXX    X.XXXXX
  Aptidão média final (suav.)    X.XXXXX    X.XXXXX
  Geração de convergência (99%)      XXX        XXX
  Maior ganho marginal (ΔF)      X.XXXXX    X.XXXXX
  Gerações sem melhoria (20%)        XXX        XXX
  ────────────────────────────────────────────────────
  Vencedor: Torneio / Roleta
════════════════════════════════════════════════════════
```

---

## Análise de Complexidade

Análise formal fundamentada em **Cormen et al. [4]** (Introduction to Algorithms, 3ª ed., MIT Press, 2009).

| Função | Complexidade de Tempo | Complexidade de Espaço |
|---|---|---|
| `segment_hits_obstacle()` | $O(1)$ | $O(1)$ |
| `compute_rssi()` | $O(N_{\text{obs}}) = O(1)$ | $O(1)$ |
| `eval_fitness()` | $O(n \cdot m)$ | $O(1)$ |
| `sel_tournament()` | $O(k_T) = O(1)$ | $O(1)$ |
| `sel_roulette()` | $O(N)$ | $O(N)$ |
| `crossover_1pt / 2pt()` | $O(m)$ | $O(m)$ |
| `mutate()` | $O(m)$ | $O(1)$ |
| `repair_budget()` | $O(m)$ | $O(m)$ |
| `evolve_one_generation()` (torneio) | $O(N \cdot n \cdot m)$ | $O(N \cdot m)$ |
| `evolve_one_generation()` (roleta) | $O(N^2 + N \cdot n \cdot m)$ | $O(N \cdot m)$ |
| **`run_experiment()` — TOTAL** | $\mathbf{O(G \cdot N \cdot n \cdot m)}$ | $\mathbf{O(N \cdot m)}$ |

### Comparação com Força Bruta

| Abordagem | Complexidade | Instância Padrão |
|---|---|---|
| **AG (este trabalho)** | $O(G \cdot N \cdot n \cdot m)$ | $600 \times 120 \times 20 \times 15 \approx$ **21,6 M ops** |
| Força Bruta | $O(2^m \cdot n)$ + busca contínua ℝ³ | $2^{15} \times 20 \approx$ **655K** nós × custo altíssimo por nó |

- **Complexidade espacial total:** $O(N \cdot m) \approx$ **102,7 KB** — viável para embarcados industriais.
- **Ordenação por aptidão** (`qsort`): $O(N \log N) = O(120 \log 120) \approx 840$ comparações/geração [4].

---

## Referências Bibliográficas

```bibtex
@article{janikow1995,
  author  = {Janikow, Cezary Z. and Clair, David},
  title   = {Simulating Nature's Methods of Evolving the Best Design Solution},
  journal = {IEEE Potentials},
  volume  = {14},
  number  = {1},
  pages   = {13--17},
  year    = {1995}
}

@article{li2009,
  author  = {Li, Ji and Andrew, Lachlan L. H. and Foh, Chuan Heng
             and Zukerman, Moshe and Chen, Hsiao-Hwa},
  title   = {Connectivity, Coverage and Placement in Wireless Sensor Networks},
  journal = {Sensors},
  volume  = {9},
  number  = {10},
  pages   = {7664--7693},
  year    = {2009},
  doi     = {10.3390/s91007664}
}

@article{marks2010,
  author  = {Marks, Michał},
  title   = {A Survey of Multi-Objective Deployment in Wireless Sensor Networks},
  journal = {Journal of Telecommunications and Information Technology},
  number  = {3},
  pages   = {36--41},
  year    = {2010}
}

@book{cormen2009,
  author    = {Cormen, Thomas H. and Leiserson, Charles E.
               and Rivest, Ronald L. and Stein, Clifford},
  title     = {Introduction to Algorithms},
  edition   = {3},
  publisher = {MIT Press},
  year      = {2009},
  isbn      = {978-0-262-03384-8}
}
```

---

## Autores

<div align="center">

| Campo | Informação |
|---|---|
| **Instituição** | Universidade Federal de Roraima — UFRR |
| **Departamento** | Departamento de Ciência da Computação — DCC |
| **Disciplina** | DCC606 — Análise de Algoritmos |
| **Projeto** | Projeto Final 2026 — Tema 5 |
| **Repositório** | `FP_DCC606_Tema_5_RR_2026` |

</div>

---

## Licença

Este projeto é disponibilizado sob a licença **MIT** para fins acadêmicos.  
Consulte o arquivo [LICENSE](LICENSE) para detalhes.

---

<div align="center">

**DCC606 · Análise de Algoritmos · UFRR 2026**

*"A evolução não é apenas um processo biológico — é um paradigma computacional."*  
— Janikow & Clair [1]

</div>
