/*
 * =============================================================================
 *  MOTOR DE OTIMIZAÇÃO EVOLUTIVA PARA REDES IIoT
 *  Algoritmos Genéticos para Síntese Topológica em Ambientes Industriais
 *  DCC606 — Projeto Final 2026
 * =============================================================================
 *
 *  Referência da Função Objetivo (Eq. 5):
 *    max F = α·VidaÚtil(S) − β·LatênciaMédia(G*) + δ·Redundância(S)
 *
 *  Modelo de Propagação Log-Distância (Eq. 4):
 *    Ploss(d) = Ploss(d0) + 10·γ·log10(d/d0) + Xσ
 *
 *  Restrições:
 *    1. RSSI(si, gj) ≥ −90 dBm  para todo sensor si  (conectividade)
 *    2. |G*| ≤ Mmax               (orçamento de infraestrutura)
 *
 *  Operadores implementados:
 *    - Codificação mista: binária (ativação) + real (coordenadas 3D)
 *    - Seleção: Torneio Estocástico  e  Roleta Viciada por Aptidão
 *    - Cruzamento: 1 ponto e 2 pontos (alternados aleatoriamente)
 *    - Mutação: bit-flip (binária) + Gaussiana (coordenadas contínuas)
 *    - Elitismo: top-k sobrevivem intactos a cada geração
 * =============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>

/* Macro auxiliar para converter constante numérica em string literal */
#define STRINGIFY_INNER(x) #x
#define STRINGIFY(x)       STRINGIFY_INNER(x)

/* ─────────────────────────────────────────────────────────────────────────────
 *  CONSTANTES DE CONFIGURAÇÃO
 * ───────────────────────────────────────────────────────────────────────────*/

/* Planta industrial */
#define PLANT_W        100.0    /* largura total (m)          */
#define PLANT_H        100.0    /* comprimento total (m)      */
#define PLANT_Z_MAX    10.0     /* altura do teto (m)         */
#define PLANT_Z_MIN     2.0     /* altura mínima dos GWs (m)  */

/* Topologia da rede */
#define N_SENSORS      20       /* sensores fixos             */
#define N_CANDIDATES   15       /* posições candidatas de GW  */
#define MAX_GW          5       /* Mmax: orçamento de GWs     */

/* Parâmetros do AG */
#define POP_SIZE       120      /* tamanho da população        */
#define MAX_GENS       600      /* número máximo de gerações   */
#define ELITE_K          6      /* k melhores preservados      */
#define TOURNAMENT_K     5      /* tamanho do torneio          */
#define P_CROSS         0.85    /* probabilidade de cruzamento */
#define P_MUT_BIN       0.04    /* taxa de mutação binária     */
#define P_MUT_REAL      0.10    /* taxa de mutação real        */
#define SIGMA_MUT_XY    5.0     /* desvio padrão Gaussiano XY  */
#define SIGMA_MUT_Z     1.5     /* desvio padrão Gaussiano Z   */

/* Modelo de Propagação Log-Distância (Eq. 4) */
#define PL_D0          1.0      /* distância de referência (m) */
#define PL_REF        40.0      /* Ploss em d0 (dB)            */
#define PL_GAMMA       2.7      /* expoente de atenuação γ     */
#define TX_POWER      20.0      /* potência de TX do ESP32 (dBm)*/
#define RSSI_THRESH  -90.0      /* limiar mínimo RSSI (dBm)    */
#define OBS_ATTEN     20.0      /* atenuação por obstáculo (dB)*/

/* Modelo de Energia */
#define BATT_CAPACITY 5000.0    /* capacidade da bateria (mAh) */
#define TX_CURR_BASE   100.0    /* corrente TX base (mA)       */
#define MAX_LIFETIME  8760.0    /* vida máxima = 1 ano (horas) */

/* Pesos da Função Multiobjetivo (Eq. 5) */
#define W_ALPHA   0.50          /* peso: VidaÚtil              */
#define W_BETA    0.30          /* peso: LatênciaMédia         */
#define W_DELTA   0.20          /* peso: Redundância           */

/* Constante de Pi */
#define PI 3.14159265358979323846

/* ─────────────────────────────────────────────────────────────────────────────
 *  ESTRUTURAS DE DADOS
 * ───────────────────────────────────────────────────────────────────────────*/

/* Ponto 3D */
typedef struct {
    double x, y, z;
} Vec3;

/* Obstáculo industrial (projeção retangular XY) */
typedef struct {
    double xmin, xmax;
    double ymin, ymax;
    double attn_db;     /* atenuação extra no sinal (dB) */
} Obstacle;

/*
 * Cromossomo de Representação Mista
 * ─────────────────────────────────
 *  Genótipo binário : active[i] ∈ {0,1}  — ativa/desativa o gateway i
 *  Genótipo real    : pos[i][3]  ∈ R³    — coordenadas (x,y,z) refinadas
 *
 *  Mapeamento genotípico→fenotípico (Janikow & Clair, 1995):
 *  - "active" codifica a DECISÃO estrutural (qual infraestrutura implantar)
 *  - "pos"    codifica os PARÂMETROS contínuos do arranjo físico
 */
typedef struct {
    int    active[N_CANDIDATES];      /* parte binária           */
    double pos[N_CANDIDATES][3];      /* parte real: x, y, z     */
    double fitness;                   /* aptidão avaliada        */
} Chromo;

/* Resultado consolidado de um experimento */
typedef struct {
    double best_fit;
    int    conv_gen;
    int    gw_count;
    int    covered;
    Chromo best_chromo;
} ExpResult;

/* ─────────────────────────────────────────────────────────────────────────────
 *  VARIÁVEIS GLOBAIS DO AMBIENTE
 * ───────────────────────────────────────────────────────────────────────────*/

static Vec3     g_sensors[N_SENSORS];
static Vec3     g_candidates[N_CANDIDATES];
static Obstacle g_obs[8];
static int      g_n_obs = 0;

static Chromo   g_pop[POP_SIZE];
static Chromo   g_new_pop[POP_SIZE];
static Chromo   g_best;

/* ─────────────────────────────────────────────────────────────────────────────
 *  FUNÇÕES UTILITÁRIAS
 * ───────────────────────────────────────────────────────────────────────────*/

/* Número aleatório ∈ [0, 1) */
static double rnd(void) {
    return (double)rand() / ((double)RAND_MAX + 1.0);
}

/* Número aleatório ∈ [lo, hi) */

/* Clamp de double */
static double clamp_d(double v, double lo, double hi) {
    return (v < lo) ? lo : (v > hi) ? hi : v;
}

/*
 * Gerador de variável aleatória Normal N(mean, sigma)
 * Método Box-Muller (aprovado para uso em AG)
 */
static double gauss(double mean, double sigma) {
    double u1 = rnd() + 1e-12;   /* evita log(0) */
    double u2 = rnd() + 1e-12;
    double z  = sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
    return mean + sigma * z;
}

/* Distância Euclidiana 3D */
static double dist3(Vec3 a, Vec3 b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    double dz = a.z - b.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

/* Comparador para qsort — ordem decrescente por aptidão */
static int cmp_desc(const void *a, const void *b) {
    double fa = ((const Chromo *)a)->fitness;
    double fb = ((const Chromo *)b)->fitness;
    return (fa < fb) - (fa > fb);
}

/* ─────────────────────────────────────────────────────────────────────────────
 *  INICIALIZAÇÃO DO AMBIENTE INDUSTRIAL
 * ───────────────────────────────────────────────────────────────────────────*/

static void setup_environment(void) {
    /*
     * Sensores fixos — pontos obrigatórios de monitoramento
     * Dispostos em grade 5x4 cobrindo a planta de 100x100m
     */
    static const double sd[N_SENSORS][3] = {
        { 8,  8, 2}, {22,  8, 2}, {38,  8, 2}, {54,  8, 2}, {70,  8, 2},
        { 8, 25, 2}, {22, 25, 2}, {38, 25, 2}, {54, 25, 2}, {70, 25, 2},
        { 8, 45, 2}, {22, 45, 2}, {38, 45, 2}, {54, 45, 2}, {70, 45, 2},
        { 8, 65, 2}, {22, 65, 2}, {38, 65, 2}, {54, 65, 2}, {70, 65, 2}
    };
    for (int i = 0; i < N_SENSORS; i++) {
        g_sensors[i].x = sd[i][0];
        g_sensors[i].y = sd[i][1];
        g_sensors[i].z = sd[i][2];
    }

    /*
     * Posições candidatas para gateways (G = {g1,...,gm})
     * O AG seleciona G* ⊆ G com |G*| ≤ Mmax
     */
    static const double cd[N_CANDIDATES][3] = {
        {18, 18, 5}, {50, 18, 5}, {82, 18, 5},
        {18, 50, 5}, {50, 50, 5}, {82, 50, 5},
        {18, 82, 5}, {50, 82, 5}, {82, 82, 5},
        {34, 34, 5}, {66, 34, 5}, {34, 66, 5},
        {66, 66, 5}, {50,  8, 5}, {50, 92, 5}
    };
    for (int i = 0; i < N_CANDIDATES; i++) {
        g_candidates[i].x = cd[i][0];
        g_candidates[i].y = cd[i][1];
        g_candidates[i].z = cd[i][2];
    }

    /*
     * Obstáculos fixos — estruturas metálicas e paredes industriais
     * Cada obstáculo atenua o sinal em 20 dB ao bloquear a LOS
     */
    g_n_obs = 4;
    g_obs[0] = (Obstacle){28, 45, 18, 40, OBS_ATTEN};  /* Estrutura metálica A */
    g_obs[1] = (Obstacle){58, 76, 18, 40, OBS_ATTEN};  /* Estrutura metálica B */
    g_obs[2] = (Obstacle){12, 32, 52, 74, OBS_ATTEN};  /* Parede norte-sul     */
    g_obs[3] = (Obstacle){52, 76, 52, 74, OBS_ATTEN};  /* Parede leste         */
}

/* ─────────────────────────────────────────────────────────────────────────────
 *  MODELO DE PROPAGAÇÃO — LOG-DISTÂNCIA (Eq. 4)
 * ───────────────────────────────────────────────────────────────────────────*/

/*
 * Verifica se o segmento de reta [A,B] (projeção XY) intersecta o obstáculo.
 * Algoritmo de interseção paramétrica Liang-Barsky simplificada.
 */
static int segment_hits_obstacle(Vec3 a, Vec3 b, const Obstacle *o) {
    double dx = b.x - a.x, dy = b.y - a.y;
    double tmin = 0.0, tmax = 1.0;

    /* eixo X */
    if (fabs(dx) < 1e-9) {
        if (a.x < o->xmin || a.x > o->xmax) return 0;
    } else {
        double t1 = (o->xmin - a.x) / dx;
        double t2 = (o->xmax - a.x) / dx;
        if (t1 > t2) { double tmp = t1; t1 = t2; t2 = tmp; }
        if (t1 > tmin) tmin = t1;
        if (t2 < tmax) tmax = t2;
        if (tmin > tmax) return 0;
    }

    /* eixo Y */
    if (fabs(dy) < 1e-9) {
        if (a.y < o->ymin || a.y > o->ymax) return 0;
    } else {
        double t1 = (o->ymin - a.y) / dy;
        double t2 = (o->ymax - a.y) / dy;
        if (t1 > t2) { double tmp = t1; t1 = t2; t2 = tmp; }
        if (t1 > tmin) tmin = t1;
        if (t2 < tmax) tmax = t2;
        if (tmin > tmax) return 0;
    }

    return 1; /* intersecção detectada */
}

/* Soma de atenuação extra por todos os obstáculos no caminho sensor→gateway */
static double extra_attenuation(Vec3 sensor, Vec3 gw) {
    double total = 0.0;
    for (int i = 0; i < g_n_obs; i++) {
        if (segment_hits_obstacle(sensor, gw, &g_obs[i]))
            total += g_obs[i].attn_db;
    }
    return total;
}

/*
 * Equação 4 — Modelo de Atenuação Log-Distância:
 *   Ploss(d) = Ploss(d0) + 10·γ·log10(d/d0) + Xσ
 *
 * RSSI(dBm) = Ptx(dBm) − Ploss(d) − AtenuaçãoObstáculos
 */
static double compute_rssi(Vec3 sensor, Vec3 gw) {
    double d = dist3(sensor, gw);
    if (d < PL_D0) d = PL_D0;

    double path_loss = PL_REF + 10.0 * PL_GAMMA * log10(d / PL_D0);
    double obs_loss  = extra_attenuation(sensor, gw);

    return TX_POWER - path_loss - obs_loss;
}

/* ─────────────────────────────────────────────────────────────────────────────
 *  OPERAÇÕES SOBRE CROMOSSOMOS
 * ───────────────────────────────────────────────────────────────────────────*/

/* Conta gateways ativos no cromossomo */
static int n_active(const Chromo *c) {
    int k = 0;
    for (int i = 0; i < N_CANDIDATES; i++) k += c->active[i];
    return k;
}

/*
 * Reparação de restrições de orçamento após mutação/cruzamento.
 * Garante: 1 ≤ |G*| ≤ Mmax  (restrição de budget)
 */
static void repair_budget(Chromo *c) {
    /* Garante ao menos 1 gateway ativo */
    if (n_active(c) == 0)
        c->active[rand() % N_CANDIDATES] = 1;

    /* Reduz excesso removendo aleatoriamente */
    int active_idx[N_CANDIDATES], cnt = 0;
    for (int i = 0; i < N_CANDIDATES; i++)
        if (c->active[i]) active_idx[cnt++] = i;

    while (cnt > MAX_GW) {
        int pick = rand() % cnt;
        c->active[active_idx[pick]] = 0;
        active_idx[pick] = active_idx[--cnt];
    }
}

/* Inicializa um cromossomo aleatório válido */
static void init_chromo(Chromo *c) {
    for (int i = 0; i < N_CANDIDATES; i++) {
        /* ativação aleatória com ~25% de chance por candidato */
        c->active[i] = (rand() % 4 == 0) ? 1 : 0;

        /*
         * Coordenadas iniciais próximas ao candidato + perturbação Gaussiana.
         * Mapeamento fenotípico: o gene real "refina" a posição estrutural.
         */
        c->pos[i][0] = clamp_d(g_candidates[i].x + gauss(0.0, 4.0),
                                0.0, PLANT_W);
        c->pos[i][1] = clamp_d(g_candidates[i].y + gauss(0.0, 4.0),
                                0.0, PLANT_H);
        c->pos[i][2] = clamp_d(g_candidates[i].z + gauss(0.0, 1.0),
                                PLANT_Z_MIN, PLANT_Z_MAX);
    }
    repair_budget(c);
    c->fitness = 0.0;
}

/* ─────────────────────────────────────────────────────────────────────────────
 *  FUNÇÃO DE APTIDÃO MULTIOBJETIVO (Eq. 5)
 * ───────────────────────────────────────────────────────────────────────────*/

/*
 * Componente 1 — VidaÚtil(S) normalizada ∈ [0, 1]
 *
 * Modelo de energia: quanto maior o RSSI recebido pelo sensor,
 * menor a potência TX necessária → menor corrente → maior vida útil.
 * Baseado no perfil de potência adaptativa do ESP32.
 */
static double compute_lifetime(const Chromo *c) {
    double total_life = 0.0;

    for (int s = 0; s < N_SENSORS; s++) {
        double best_rssi = -999.0;

        for (int g = 0; g < N_CANDIDATES; g++) {
            if (!c->active[g]) continue;
            Vec3 gp = {c->pos[g][0], c->pos[g][1], c->pos[g][2]};
            double r = compute_rssi(g_sensors[s], gp);
            if (r > best_rssi) best_rssi = r;
        }

        if (best_rssi >= RSSI_THRESH) {
            /* Margem acima do limiar → fração de corrente TX economizada */
            double margin   = best_rssi - RSSI_THRESH;         /* dB  */
            double tx_frac  = 1.0 / (1.0 + margin / 25.0);    /* ∈ (0,1] */
            double current  = TX_CURR_BASE * tx_frac;          /* mA  */
            double life_h   = BATT_CAPACITY / current;         /* h   */
            total_life += (life_h > MAX_LIFETIME) ? MAX_LIFETIME : life_h;
        }
        /* sensor não coberto contribui 0 para o somatório */
    }

    return (total_life / N_SENSORS) / MAX_LIFETIME; /* normalizado */
}

/*
 * Componente 2 — LatênciaMédia(G*) normalizada ∈ [0, 1]
 *
 * Proxy de latência: distância ao gateway mais próximo com RSSI válido.
 * Distância máxima possível na planta: diagonal = √2 × 100 m.
 */
static double compute_latency(const Chromo *c) {
    double total_dist = 0.0;
    int    connected  = 0;

    for (int s = 0; s < N_SENSORS; s++) {
        double best_dist = DBL_MAX;

        for (int g = 0; g < N_CANDIDATES; g++) {
            if (!c->active[g]) continue;
            Vec3 gp = {c->pos[g][0], c->pos[g][1], c->pos[g][2]};
            double r = compute_rssi(g_sensors[s], gp);
            if (r >= RSSI_THRESH) {
                double d = dist3(g_sensors[s], gp);
                if (d < best_dist) best_dist = d;
            }
        }
        if (best_dist < DBL_MAX) {
            total_dist += best_dist;
            connected++;
        }
    }

    if (connected == 0) return 1.0; /* penalidade máxima */

    double avg = total_dist / connected;
    double max_dist = sqrt(2.0) * (PLANT_W > PLANT_H ? PLANT_W : PLANT_H);
    return avg / max_dist; /* normalizado */
}

/*
 * Componente 3 — Redundância(S) normalizada ∈ [0, 1]
 *
 * Número médio de gateways alcançáveis por sensor,
 * normalizado pelo máximo possível (Mmax).
 * Mitigação de Pontos Únicos de Falha.
 */
static double compute_redundancy(const Chromo *c) {
    double total_paths = 0.0;

    for (int s = 0; s < N_SENSORS; s++) {
        for (int g = 0; g < N_CANDIDATES; g++) {
            if (!c->active[g]) continue;
            Vec3 gp = {c->pos[g][0], c->pos[g][1], c->pos[g][2]};
            if (compute_rssi(g_sensors[s], gp) >= RSSI_THRESH)
                total_paths += 1.0;
        }
    }

    return (total_paths / N_SENSORS) / (double)MAX_GW; /* normalizado */
}

/* Fração de sensores cobertos ∈ [0, 1] */
static double coverage_ratio(const Chromo *c) {
    int covered = 0;
    for (int s = 0; s < N_SENSORS; s++) {
        for (int g = 0; g < N_CANDIDATES; g++) {
            if (!c->active[g]) continue;
            Vec3 gp = {c->pos[g][0], c->pos[g][1], c->pos[g][2]};
            if (compute_rssi(g_sensors[s], gp) >= RSSI_THRESH) {
                covered++;
                break; /* basta 1 gateway alcançável */
            }
        }
    }
    return (double)covered / N_SENSORS;
}

/*
 * Avaliação da Aptidão — Equação 5:
 *   F = α·VidaÚtil − β·Latência + δ·Redundância − penalidades
 *
 * Penalidade de cobertura: desincentiva soluções com sensores "cegos".
 * Penalidade de orçamento: salvaguarda caso repair_budget falhe.
 */
static void eval_fitness(Chromo *c) {
    double cov  = coverage_ratio(c);
    double life = compute_lifetime(c);
    double lat  = compute_latency(c);
    double red  = compute_redundancy(c);

    /* Penalidade proporcional à cobertura perdida */
    double pen_cov = (1.0 - cov) * 2.0;

    /* Penalidade por excesso de gateways (restrição de budget) */
    int    act     = n_active(c);
    double pen_bud = (act > MAX_GW) ? (double)(act - MAX_GW) * 0.5 : 0.0;

    c->fitness = W_ALPHA * life
               - W_BETA  * lat
               + W_DELTA * red
               - pen_cov
               - pen_bud;
}

/* ─────────────────────────────────────────────────────────────────────────────
 *  OPERADORES DE SELEÇÃO
 * ───────────────────────────────────────────────────────────────────────────*/

/*
 * Seleção por Torneio Estocástico
 * ─────────────────────────────────
 * Seleciona TOURNAMENT_K indivíduos aleatoriamente e retorna o melhor.
 * Vantagem: pressão de seleção controlável; sem necessidade de ordenação global.
 * Complexidade: O(k)
 */
static int sel_tournament(void) {
    int best = rand() % POP_SIZE;
    for (int i = 1; i < TOURNAMENT_K; i++) {
        int c = rand() % POP_SIZE;
        if (g_pop[c].fitness > g_pop[best].fitness)
            best = c;
    }
    return best;
}

/*
 * Seleção por Roleta Viciada por Aptidão (Fitness-Proportionate)
 * ─────────────────────────────────────────────────────────────────
 * P(i) = f_shifted(i) / Σ f_shifted(j)
 *
 * Aptidões são deslocadas para garantir positividade antes da normalização.
 * Vantagem: favorece fortemente os melhores indivíduos.
 * Risco: convergência prematura se um indivíduo dominar muito.
 */
static int sel_roulette(void) {
    double min_f = g_pop[0].fitness;
    for (int i = 1; i < POP_SIZE; i++)
        if (g_pop[i].fitness < min_f) min_f = g_pop[i].fitness;

    double total = 0.0;
    double shifted[POP_SIZE];
    for (int i = 0; i < POP_SIZE; i++) {
        shifted[i] = g_pop[i].fitness - min_f + 1e-9; /* ε evita 0 */
        total += shifted[i];
    }

    double r = rnd() * total, acc = 0.0;
    for (int i = 0; i < POP_SIZE; i++) {
        acc += shifted[i];
        if (r <= acc) return i;
    }
    return POP_SIZE - 1; /* fallback numérico */
}

/* ─────────────────────────────────────────────────────────────────────────────
 *  OPERADORES DE CRUZAMENTO
 * ───────────────────────────────────────────────────────────────────────────*/

/*
 * Cruzamento de 1 Ponto
 * ──────────────────────
 * Divide os cromossomos em [0, pt) e [pt, N_CANDIDATES).
 * Opera simultaneamente sobre os genes binários e reais.
 */
static void crossover_1pt(const Chromo *p1, const Chromo *p2,
                           Chromo *c1, Chromo *c2)
{
    int pt = 1 + rand() % (N_CANDIDATES - 1); /* ponto de corte */

    for (int i = 0; i < N_CANDIDATES; i++) {
        const Chromo *src1 = (i < pt) ? p1 : p2;
        const Chromo *src2 = (i < pt) ? p2 : p1;

        c1->active[i] = src1->active[i];
        c2->active[i] = src2->active[i];

        for (int d = 0; d < 3; d++) {
            c1->pos[i][d] = src1->pos[i][d];
            c2->pos[i][d] = src2->pos[i][d];
        }
    }
}

/*
 * Cruzamento de 2 Pontos
 * ───────────────────────
 * Segmento interno [pt1, pt2) é trocado entre os pais.
 * Preserva melhor o contexto das extremidades do cromossomo.
 */
static void crossover_2pt(const Chromo *p1, const Chromo *p2,
                           Chromo *c1, Chromo *c2)
{
    int pt1 = rand() % N_CANDIDATES;
    int pt2 = rand() % N_CANDIDATES;
    if (pt1 > pt2) { int t = pt1; pt1 = pt2; pt2 = t; }

    for (int i = 0; i < N_CANDIDATES; i++) {
        int inner = (i >= pt1 && i < pt2); /* dentro do segmento trocado? */
        const Chromo *src1 = inner ? p2 : p1;
        const Chromo *src2 = inner ? p1 : p2;

        c1->active[i] = src1->active[i];
        c2->active[i] = src2->active[i];

        for (int d = 0; d < 3; d++) {
            c1->pos[i][d] = src1->pos[i][d];
            c2->pos[i][d] = src2->pos[i][d];
        }
    }
}

/* ─────────────────────────────────────────────────────────────────────────────
 *  OPERADOR DE MUTAÇÃO MISTA
 * ───────────────────────────────────────────────────────────────────────────*/

/*
 * Mutação de Bit-Flip (genes binários)
 *   active[i] ← 1 - active[i]  com probabilidade P_MUT_BIN
 *
 * Mutação Gaussiana (genes reais de coordenadas)
 *   pos[i][d] ← pos[i][d] + N(0, σ)  com probabilidade P_MUT_REAL
 *
 * Após mutação, repair_budget corrige violações de orçamento.
 */
static void mutate(Chromo *c) {
    /* Mutação binária: inversão de bit */
    for (int i = 0; i < N_CANDIDATES; i++)
        if (rnd() < P_MUT_BIN)
            c->active[i] ^= 1;

    /* Mutação Gaussiana nas coordenadas contínuas */
    for (int i = 0; i < N_CANDIDATES; i++) {
        if (rnd() < P_MUT_REAL) {
            c->pos[i][0] = clamp_d(c->pos[i][0] + gauss(0.0, SIGMA_MUT_XY),
                                   0.0, PLANT_W);
            c->pos[i][1] = clamp_d(c->pos[i][1] + gauss(0.0, SIGMA_MUT_XY),
                                   0.0, PLANT_H);
            c->pos[i][2] = clamp_d(c->pos[i][2] + gauss(0.0, SIGMA_MUT_Z),
                                   PLANT_Z_MIN, PLANT_Z_MAX);
        }
    }

    repair_budget(c);
}

/* ─────────────────────────────────────────────────────────────────────────────
 *  GERENCIAMENTO DA POPULAÇÃO
 * ───────────────────────────────────────────────────────────────────────────*/

/* Inicializa população aleatória e avalia aptidões iniciais */
static void init_population(unsigned seed) {
    srand(seed);
    for (int i = 0; i < POP_SIZE; i++) {
        init_chromo(&g_pop[i]);
        eval_fitness(&g_pop[i]);
    }
    g_best = g_pop[0];
    for (int i = 1; i < POP_SIZE; i++)
        if (g_pop[i].fitness > g_best.fitness)
            g_best = g_pop[i];
}

/*
 * Evolução de Uma Geração
 * ────────────────────────
 * 1. Ordena população (qsort descendente)
 * 2. Aplica elitismo: copia top-k sem alteração
 * 3. Preenche o restante com seleção → cruzamento → mutação
 * 4. Atualiza g_best
 *
 * sel_method: 0 = torneio, 1 = roleta
 */
static void evolve_one_generation(int sel_method) {
    int i;

    /* Ordena por aptidão descendente */
    qsort(g_pop, POP_SIZE, sizeof(Chromo), cmp_desc);

    /* ── Elitismo: os k melhores sobrevivem intactos ── */
    for (i = 0; i < ELITE_K; i++)
        g_new_pop[i] = g_pop[i];

    /* ── Geração de filhos ── */
    int idx = ELITE_K;
    while (idx < POP_SIZE) {
        /* Seleção dos pais */
        int pi1 = (sel_method == 0) ? sel_tournament() : sel_roulette();
        int pi2 = (sel_method == 0) ? sel_tournament() : sel_roulette();

        Chromo c1, c2;

        /* Cruzamento (alternância 1 ponto / 2 pontos para diversidade) */
        if (rnd() < P_CROSS) {
            if (rand() % 2 == 0)
                crossover_1pt(&g_pop[pi1], &g_pop[pi2], &c1, &c2);
            else
                crossover_2pt(&g_pop[pi1], &g_pop[pi2], &c1, &c2);
        } else {
            c1 = g_pop[pi1];
            c2 = g_pop[pi2];
        }

        /* Mutação e avaliação */
        mutate(&c1); eval_fitness(&c1);
        mutate(&c2); eval_fitness(&c2);

        g_new_pop[idx++] = c1;
        if (idx < POP_SIZE)
            g_new_pop[idx++] = c2;
    }

    memcpy(g_pop, g_new_pop, sizeof(Chromo) * POP_SIZE);

    /* Atualiza melhor global */
    for (i = 0; i < POP_SIZE; i++)
        if (g_pop[i].fitness > g_best.fitness)
            g_best = g_pop[i];
}

/* ─────────────────────────────────────────────────────────────────────────────
 *  SAÍDA: MAPA ASCII DA PLANTA
 * ───────────────────────────────────────────────────────────────────────────*/

#define MAP_W 54
#define MAP_H 27

static void print_map(const Chromo *c) {
    char map[MAP_H][MAP_W + 1];
    int  x, y;

    for (y = 0; y < MAP_H; y++) {
        for (x = 0; x < MAP_W; x++) map[y][x] = ' ';
        map[y][MAP_W] = '\0';
    }

    /* Obstáculos */
    for (int o = 0; o < g_n_obs; o++) {
        int x1 = (int)(g_obs[o].xmin / PLANT_W * MAP_W);
        int x2 = (int)(g_obs[o].xmax / PLANT_W * MAP_W);
        int y1 = (int)(g_obs[o].ymin / PLANT_H * MAP_H);
        int y2 = (int)(g_obs[o].ymax / PLANT_H * MAP_H);
        for (y = y1; y < y2 && y < MAP_H; y++)
            for (x = x1; x < x2 && x < MAP_W; x++)
                map[y][x] = '#';
    }

    /* Sensores */
    for (int s = 0; s < N_SENSORS; s++) {
        int mx = (int)(g_sensors[s].x / PLANT_W * MAP_W);
        int my = (int)(g_sensors[s].y / PLANT_H * MAP_H);
        if (mx >= MAP_W) mx = MAP_W - 1;
        if (my >= MAP_H) my = MAP_H - 1;
        if (map[my][mx] == ' ') map[my][mx] = 's';
    }

    /* Gateways ativos */
    for (int g = 0; g < N_CANDIDATES; g++) {
        if (!c->active[g]) continue;
        int mx = (int)(c->pos[g][0] / PLANT_W * MAP_W);
        int my = (int)(c->pos[g][1] / PLANT_H * MAP_H);
        if (mx >= MAP_W) mx = MAP_W - 1;
        if (my >= MAP_H) my = MAP_H - 1;
        map[my][mx] = 'G';
    }

    printf("\n  Mapa Topológico — Planta 100 × 100 m\n");
    printf("  Legenda: [G]=Gateway  [s]=Sensor  [#]=Obstáculo\n\n");
    printf("  +");
    for (x = 0; x < MAP_W; x++) putchar('-');
    puts("+");
    for (y = MAP_H - 1; y >= 0; y--)
        printf("  |%s|\n", map[y]);
    printf("  +");
    for (x = 0; x < MAP_W; x++) putchar('-');
    puts("+\n");
}

/* ─────────────────────────────────────────────────────────────────────────────
 *  SAÍDA: RELATÓRIO DETALHADO DA SOLUÇÃO
 * ───────────────────────────────────────────────────────────────────────────*/

static void print_solution(const Chromo *c) {
    int s, g;

    puts("\n=================================================================");
    puts("  MELHOR SOLUÇÃO ENCONTRADA");
    puts("=================================================================");
    printf("  Aptidão (F):      %.6f\n", c->fitness);
    printf("  Gateways Ativos:  %d / %d (orçamento)\n", n_active(c), MAX_GW);
    printf("  Cobertura:        %.1f %%\n", coverage_ratio(c) * 100.0);

    puts("\n  ┌─────────────────── Gateways Ativos ────────────────────────┐");
    printf("  │ %-4s %-9s %-9s %-9s %-10s %-9s │\n",
           "GW#", "X (m)", "Y (m)", "Z (m)", "RSSI_min", "Cobertos");
    puts("  ├─────────────────────────────────────────────────────────────┤");

    for (g = 0; g < N_CANDIDATES; g++) {
        if (!c->active[g]) continue;

        Vec3   gp  = {c->pos[g][0], c->pos[g][1], c->pos[g][2]};
        double wr  = 0.0;        /* worst RSSI para este GW    */
        int    cov = 0;          /* sensores cobertos por este GW */
        int    first = 1;

        for (s = 0; s < N_SENSORS; s++) {
            double r = compute_rssi(g_sensors[s], gp);
            if (first || r < wr) { wr = r; first = 0; }
            if (r >= RSSI_THRESH) cov++;
        }

        printf("  │ %-4d %-9.2f %-9.2f %-9.2f %-10.2f %-9d │\n",
               g + 1, gp.x, gp.y, gp.z, wr, cov);
    }
    puts("  └─────────────────────────────────────────────────────────────┘");

    puts("\n  ┌────────────── Cobertura por Sensor ─────────────────────────┐");
    printf("  │ %-6s %-12s %-10s %-8s                   │\n",
           "Sensor", "RSSI (dBm)", "Caminhos", "Status");
    puts("  ├─────────────────────────────────────────────────────────────┤");

    int all_ok = 1;
    for (s = 0; s < N_SENSORS; s++) {
        double br   = -999.0;
        int    paths = 0;

        for (g = 0; g < N_CANDIDATES; g++) {
            if (!c->active[g]) continue;
            Vec3   gp = {c->pos[g][0], c->pos[g][1], c->pos[g][2]};
            double r  = compute_rssi(g_sensors[s], gp);
            if (r > br) br = r;
            if (r >= RSSI_THRESH) paths++;
        }

        const char *status = (br >= RSSI_THRESH) ? "OK" : "SEM SINAL";
        if (br < RSSI_THRESH) all_ok = 0;

        printf("  │ %-6d %-12.2f %-10d %-8s                   │\n",
               s + 1, br, paths, status);
    }
    puts("  └─────────────────────────────────────────────────────────────┘");

    double life_h = compute_lifetime(c) * MAX_LIFETIME;
    double lat_m  = compute_latency(c)
                  * sqrt(2.0) * (PLANT_W > PLANT_H ? PLANT_W : PLANT_H);
    double red    = compute_redundancy(c) * MAX_GW;

    puts("\n  ┌────────────── Métricas de Rede ─────────────────────────────┐");
    printf("  │  Vida Útil Estimada  : %7.1f h  (%5.1f dias)            │\n",
           life_h, life_h / 24.0);
    printf("  │  Distância Média GW  : %7.2f m                           │\n",
           lat_m);
    printf("  │  Redundância Média   : %7.2f caminhos/sensor             │\n",
           red);
    printf("  │  Conectividade Total : %s                         │\n",
           all_ok ? "SIM — RESTRIÇÃO SATISFEITA   " : "NÃO — RESTRIÇÃO VIOLADA  !!!!");
    puts("  └─────────────────────────────────────────────────────────────┘\n");
}

/* ─────────────────────────────────────────────────────────────────────────────
 *  EXECUÇÃO DE UM EXPERIMENTO COMPLETO
 * ───────────────────────────────────────────────────────────────────────────*/

static ExpResult run_experiment(const char *label, int sel_method, unsigned seed) {
    int g;

    printf("\n─────────────────────────────────────────────────────────────────\n");
    printf("  Experimento: %s\n", label);
    printf("─────────────────────────────────────────────────────────────────\n");
    printf("  %-9s  %-16s  %-16s\n", "Geração", "Melhor Aptidão", "Aptidão Média");
    puts("  ─────────────────────────────────────────────────────────");

    init_population(seed);

    double hist_best[MAX_GENS];
    double hist_avg[MAX_GENS];
    int    conv_gen = MAX_GENS - 1;

    for (g = 0; g < MAX_GENS; g++) {
        evolve_one_generation(sel_method);

        double avg = 0.0;
        int    i;
        for (i = 0; i < POP_SIZE; i++) avg += g_pop[i].fitness;
        avg /= POP_SIZE;

        hist_best[g] = g_best.fitness;
        hist_avg[g]  = avg;

        if (g % 100 == 0 || g == MAX_GENS - 1)
            printf("  %-9d  %-16.6f  %-16.6f\n", g, g_best.fitness, avg);
    }

    /* Geração de convergência (quando atingiu 99% do melhor final) */
    double threshold = g_best.fitness * 0.99;
    for (g = 0; g < MAX_GENS; g++) {
        if (hist_best[g] >= threshold) { conv_gen = g; break; }
    }

    /* Exporta histórico de aptidão para análise externa */
    char fname[80];
    snprintf(fname, sizeof(fname), "hist_%s.csv",
             sel_method == 0 ? "torneio" : "roleta");
    FILE *fp = fopen(fname, "w");
    if (fp) {
        fprintf(fp, "geracao,melhor_aptidao,aptidao_media\n");
        for (g = 0; g < MAX_GENS; g++)
            fprintf(fp, "%d,%.6f,%.6f\n", g, hist_best[g], hist_avg[g]);
        fclose(fp);
        printf("\n  [INFO] Histórico exportado: %s\n", fname);
    }

    /* Conta sensores cobertos na melhor solução */
    int covered = 0, s;
    for (s = 0; s < N_SENSORS; s++) {
        for (g = 0; g < N_CANDIDATES; g++) {
            if (!g_best.active[g]) continue;
            Vec3 gp = {g_best.pos[g][0], g_best.pos[g][1], g_best.pos[g][2]};
            if (compute_rssi(g_sensors[s], gp) >= RSSI_THRESH) {
                covered++; break;
            }
        }
    }

    print_map(&g_best);
    print_solution(&g_best);

    ExpResult res;
    res.best_fit     = g_best.fitness;
    res.conv_gen     = conv_gen;
    res.gw_count     = n_active(&g_best);
    res.covered      = covered;
    res.best_chromo  = g_best;
    return res;
}

/* ─────────────────────────────────────────────────────────────────────────────
 *  FUNÇÃO PRINCIPAL
 * ───────────────────────────────────────────────────────────────────────────*/

int main(void) {
    unsigned seed = (unsigned)time(NULL);

    /* Banner */
    puts("=================================================================");
    puts("  Motor de Otimização Evolutiva — Redes IIoT");
    puts("  Algoritmos Genéticos para Síntese Topológica");
    puts("  DCC606 — Projeto Final 2026");
    puts("=================================================================");
    printf("  Planta    : %.0f x %.0f m  |  Sensores: %d  |  Obstáculos: 4\n",
           PLANT_W, PLANT_H, N_SENSORS);
    printf("  Candidatos: %d GWs        |  Orçamento: Mmax = %d\n",
           N_CANDIDATES, MAX_GW);
    printf("  População : %d            |  Gerações:  %d\n",
           POP_SIZE, MAX_GENS);
    printf("  Elitismo  : top-%d        |  Semente:   %u\n", ELITE_K, seed);
    printf("  Crossover : %.0f%%        |  Mut. bin.: %.0f%%  Mut. real: %.0f%%\n",
           P_CROSS * 100, P_MUT_BIN * 100, P_MUT_REAL * 100);
    puts("=================================================================\n");

    setup_environment();
    printf("  [OK] Ambiente industrial inicializado.\n");
    printf("       %d sensores | %d candidatos de gateway | %d obstáculos\n\n",
           N_SENSORS, N_CANDIDATES, g_n_obs);

    /* ── Experimento 1: Seleção por Torneio ── */
    ExpResult r_torneio = run_experiment(
        "Torneio Estocástico (k=" STRINGIFY(TOURNAMENT_K) ")",
        0, seed);

    /* ── Experimento 2: Seleção por Roleta ── */
    ExpResult r_roleta = run_experiment(
        "Roleta Viciada por Aptidão",
        1, seed + 7919 /* primo → independência estatística */);

    /* ── Tabela Comparativa Final ── */
    puts("=================================================================");
    puts("  ANÁLISE COMPARATIVA DOS OPERADORES DE SELEÇÃO");
    puts("=================================================================");
    printf("  %-30s %-10s %-12s %-8s %-10s\n",
           "Método", "Aptidão", "Convergência", "GWs", "Cobertura");
    puts("  ───────────────────────────────────────────────────────────");
    printf("  %-30s %-10.5f %-12d %-8d %d/%d\n",
           "Torneio Estocástico",
           r_torneio.best_fit, r_torneio.conv_gen,
           r_torneio.gw_count, r_torneio.covered, N_SENSORS);
    printf("  %-30s %-10.5f %-12d %-8d %d/%d\n",
           "Roleta Viciada",
           r_roleta.best_fit,  r_roleta.conv_gen,
           r_roleta.gw_count,  r_roleta.covered,  N_SENSORS);
    puts("  ───────────────────────────────────────────────────────────");

    if (r_torneio.best_fit >= r_roleta.best_fit)
        puts("  VENCEDOR ► Seleção por Torneio Estocástico");
    else
        puts("  VENCEDOR ► Seleção por Roleta Viciada por Aptidão");

    puts("\n  Arquivos gerados:");
    puts("    hist_torneio.csv — histórico de aptidão (torneio)");
    puts("    hist_roleta.csv  — histórico de aptidão (roleta)");
    puts("    (use plot.gnu para visualizar a convergência com gnuplot)");
    puts("=================================================================\n");

    return 0;
}

