"""
=============================================================================
 Visualização de Resultados — Motor de Otimização Evolutiva IIoT
 DCC606 — Projeto Final 2026
=============================================================================
 Gera 4 gráficos a partir das médias de múltiplas amostras (ex: 13 rodadas)
 produzidas pelo programa em C, de forma a garantir validade estatística:
   1. Curvas de convergência — Melhor Aptidão por geração
   2. Evolução da Aptidão Média por geração
   3. Comparativo sobrepostos (melhor + média dos dois métodos)
   4. Análise de ganho marginal (variação da melhor aptidão por geração)

 Uso (com ambiente virtual ativado):
   python3 plot_resultados.py
=============================================================================
"""

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
#  CONFIGURAÇÃO GLOBAL DE ESTILO
# ─────────────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "figure.dpi":       150,
    "figure.facecolor": "#F8F9FA",
    "axes.facecolor":   "#FFFFFF",
    "axes.grid":        True,
    "axes.spines.top":  False,
    "axes.spines.right":False,
    "axes.labelsize":   12,
    "axes.titlesize":   13,
    "axes.titleweight": "bold",
    "grid.linestyle":   "--",
    "grid.color":       "#CCCCCC",
    "legend.framealpha":0.9
})

COR_TORNEIO = "#1f77b4"
COR_ROLETA = "#d62728"

# ─────────────────────────────────────────────────────────────────────────────
#  FUNÇÕES DE LEITURA E PROCESSAMENTO
# ─────────────────────────────────────────────────────────────────────────────
def carregar_lote_csv(prefixo: str, num_amostras: int) -> pd.DataFrame:
    """Lê múltiplas amostras e devolve a média agrupada por geração."""
    dfs = []
    base = os.path.dirname(os.path.abspath(__file__))
    
    for i in range(1, num_amostras + 1):
        caminho = os.path.join(base, f"{prefixo}_amostra_{i}.csv")
        if os.path.exists(caminho):
            dfs.append(pd.read_csv(caminho))
        else:
            print(f"[AVISO] Ficheiro não encontrado: {caminho}")
    
    if not dfs:
        print(f"[ERRO] Nenhuma amostra encontrada com o prefixo '{prefixo}'.")
        print("Certifique-se de que executou o laço no terminal (bash) primeiro.")
        sys.exit(1)
        
    # Concatena todos os DataFrames e calcula a média exata para cada geração
    df_concat = pd.concat(dfs)
    df_medio = df_concat.groupby("geracao").mean().reset_index()
    return df_medio

def suavizar(serie: pd.Series, janela: int = 30) -> pd.Series:
    """Suaviza a curva usando média móvel para visualização de tendências."""
    return serie.rolling(window=janela, min_periods=1).mean()

def imprimir_resumo(df_t: pd.DataFrame, df_r: pd.DataFrame) -> None:
    """Imprime um resumo estatístico das duas abordagens no terminal."""
    # Pressupõe que a última linha contém o valor final atingido
    best_t = df_t["melhor_aptidao"].iloc[-1] if "melhor_aptidao" in df_t else 0
    best_r = df_r["melhor_aptidao"].iloc[-1] if "melhor_aptidao" in df_r else 0
    
    print("\n" + "="*60)
    print(" RESUMO ESTATÍSTICO (MÉDIA DAS AMOSTRAS)")
    print("="*60)
    print(f" Torneio Estocástico : Melhor Aptidão Média Final = {best_t:.4f}")
    print(f" Roleta Viciada      : Melhor Aptidão Média Final = {best_r:.4f}")
    print("="*60 + "\n")

# ─────────────────────────────────────────────────────────────────────────────
#  FUNÇÕES DE PLOTAGEM
# ─────────────────────────────────────────────────────────────────────────────
def criar_janela(titulo: str):
    """Cria uma figura e um eixo com o título e rodapé padrão."""
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.canvas.manager.set_window_title(titulo)
    ax.set_title(titulo, pad=15)
    ax.set_xlabel("Geração")
    
    # Rodapé padrão do Projeto
    fig.text(
        0.5, 0.02,
        "DCC606 · Projeto Final 2026 — Otimização Topológica de Redes IIoT",
        ha="center", fontsize=8, color="gray"
    )
    return fig, ax

def grafico_convergencia(df_t: pd.DataFrame, df_r: pd.DataFrame, ax: plt.Axes) -> None:
    ax.plot(df_t["geracao"], df_t["melhor_aptidao"], label="Torneio (Melhor)", color=COR_TORNEIO, lw=2)
    ax.plot(df_r["geracao"], df_r["melhor_aptidao"], label="Roleta (Melhor)", color=COR_ROLETA, lw=2)
    ax.set_ylabel("Aptidão (Fitness)")
    ax.legend()

def grafico_aptidao_media(df_t: pd.DataFrame, df_r: pd.DataFrame, ax: plt.Axes) -> None:
    if "aptidao_media" in df_t.columns and "aptidao_media" in df_r.columns:
        ax.plot(df_t["geracao"], df_t["aptidao_media"], label="Torneio (Média Pop.)", color=COR_TORNEIO, linestyle="--", lw=2)
        ax.plot(df_r["geracao"], df_r["aptidao_media"], label="Roleta (Média Pop.)", color=COR_ROLETA, linestyle="--", lw=2)
        ax.set_ylabel("Aptidão Média")
        ax.legend()

def grafico_comparativo(df_t: pd.DataFrame, df_r: pd.DataFrame, ax: plt.Axes) -> None:
    ax.plot(df_t["geracao"], df_t["melhor_aptidao"], label="Torneio (Melhor)", color=COR_TORNEIO, lw=2)
    if "aptidao_media" in df_t.columns:
        ax.plot(df_t["geracao"], df_t["aptidao_media"], label="Torneio (Média Pop.)", color=COR_TORNEIO, lw=1.5, linestyle=":", alpha=0.7)
    
    ax.plot(df_r["geracao"], df_r["melhor_aptidao"], label="Roleta (Melhor)", color=COR_ROLETA, lw=2)
    if "aptidao_media" in df_r.columns:
        ax.plot(df_r["geracao"], df_r["aptidao_media"], label="Roleta (Média Pop.)", color=COR_ROLETA, lw=1.5, linestyle=":", alpha=0.7)
        
    ax.set_ylabel("Aptidão")
    ax.legend()

def grafico_ganho_marginal(df_t: pd.DataFrame, df_r: pd.DataFrame, ax: plt.Axes) -> None:
    for df, cor, nome in [(df_t, COR_TORNEIO, "Torneio"),
                          (df_r, COR_ROLETA,  "Roleta")]:
        if "melhor_aptidao" in df.columns:
            # clip(lower=0) ignora quedas de aptidão, considerando apenas ganhos positivos
            delta = df["melhor_aptidao"].diff().clip(lower=0)
            delta_suave = suavizar(delta, janela=30)
            ax.plot(df["geracao"], delta_suave, color=cor, lw=1.5, label=nome)
            ax.fill_between(df["geracao"], 0, delta_suave, color=cor, alpha=0.15)
            
    ax.set_ylabel("Ganho (Δ Aptidão)")
    ax.legend()

# ─────────────────────────────────────────────────────────────────────────────
#  MAIN
# ─────────────────────────────────────────────────────────────────────────────
def main() -> None:
    # ⚠️ Ajuste este número para corresponder à quantidade de rodadas feitas no Bash
    NUM_AMOSTRAS = 13
    
    print(f"A iniciar processamento de {NUM_AMOSTRAS} amostras de ficheiros CSV...")
    
    # Carrega e agrupa as médias matemáticas
    df_t = carregar_lote_csv("hist_torneio", NUM_AMOSTRAS)
    df_r = carregar_lote_csv("hist_roleta", NUM_AMOSTRAS)

    imprimir_resumo(df_t, df_r)

    # Definição do pipeline de gráficos: (função_geradora, Título da Janela, nome_do_ficheiro)
    graficos = [
        (grafico_convergencia,    "Convergência — Melhor Aptidão",         "g1_convergencia.png"),
        (grafico_aptidao_media,   "Aptidão Média da População",            "g2_aptidao_media.png"),
        (grafico_comparativo,     "Comparativo Geral — Torneio vs Roleta", "g3_comparativo.png"),
        (grafico_ganho_marginal,  "Ganho Marginal por Geração",            "g4_ganho_marginal.png"),
    ]

    for func, titulo, arquivo in graficos:
        fig, ax = criar_janela(titulo)
        func(df_t, df_r, ax)
        
        # Ajuste para evitar que as margens cortem as etiquetas e o rodapé
        fig.tight_layout()
        fig.subplots_adjust(bottom=0.15) 
        
        # Guarda a imagem localmente
        fig.savefig(arquivo)
        print(f"Gráfico renderizado e guardado: {arquivo}")

    # Exibe as 4 janelas no ecrã no final
    plt.show()

if __name__ == "__main__":
    main()