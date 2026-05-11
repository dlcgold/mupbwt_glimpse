import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def parse_bench(filepath):
    filename = os.path.basename(filepath)
    name = filename.replace("bench_", "").replace(".tsv", "")
    df = pd.read_csv(filepath, sep='\t')
    
    if name == "std":
        return {'Mode': 'std', 'MAF': 0.0, 'Label': 'Default', 'Time_s': df['s'].iloc[0], 'RAM_MB': df['max_rss'].iloc[0]}
    else:
        maf_val = float(name.replace("common_maf", ""))
        return {'Mode': 'common', 'MAF': maf_val, 'Label': str(maf_val), 'Time_s': df['s'].iloc[0], 'RAM_MB': df['max_rss'].iloc[0]}

def parse_conc(filepath):
    filename = os.path.basename(filepath)
    name = filename.replace("concordance_", "").split(".rsquare")[0]
    
    if name == "std":
        mode = "std"
        maf_val = 0.0
        label = "Default"
    else:
        mode = "common"
        maf_val = float(name.replace("common_maf", ""))
        label = str(maf_val)
        
    df = pd.read_csv(filepath, header=None, sep=r'\s+')
    
    # Rinominiamo le colonne importanti per chiarezza
    # Col 0: Bin ID, Col 1: Numero varianti, Col 2: MAF media del bin
    df.rename(columns={0: 'Bin_ID', 1: 'N_variants', 2: 'Bin_Avg_MAF'}, inplace=True)
    
    # L'R2 è l'ultima colonna
    df['R2'] = df.iloc[:, -1]
    
    # Calcoliamo la media per il cruscotto principale
    mean_r2 = df['R2'].mean()
    summary = {'Mode': mode, 'MAF': maf_val, 'Label': label, 'Mean_R2': mean_r2}
    
    # Aggiungiamo i metadati al dataframe dei bin per il secondo plot
    df['Mode'] = mode
    df['Threshold'] = label
    
    return summary, df[['Mode', 'Threshold', 'Bin_ID', 'N_variants', 'Bin_Avg_MAF', 'R2']]


# --- 1. CARICAMENTO DATI ---
bench_data = [parse_bench(f) for f in snakemake.input.benchmarks]
conc_results = [parse_conc(f) for f in snakemake.input.concordance]

# Separiamo il riassunto (per il Pareto) dai dati granulari (per i bin)
conc_summary = [res[0] for res in conc_results]
df_bins = pd.concat([res[1] for res in conc_results], ignore_index=True)


# --- 2. GRAFICO 1: DASHBOARD PARETO (Medie complessive) ---
df = pd.merge(pd.DataFrame(bench_data), pd.DataFrame(conc_summary), on=['Mode', 'MAF', 'Label'])
df_common = df[df['Mode'] == 'common'].sort_values('MAF')
df_std = df[df['Mode'] == 'std']

sns.set_theme(style="whitegrid")
fig, axs = plt.subplots(2, 2, figsize=(16, 12))

sns.scatterplot(data=df, x='Time_s', y='Mean_R2', hue='Label', style='Mode', s=200, palette='tab10', ax=axs[0, 0])
axs[0, 0].set_title("Pareto: Accuracy vs Phasing Time", fontsize=14, fontweight='bold')
axs[0, 0].set_xlabel("Phasing Time (Seconds) \u2192 Faster is better")
axs[0, 0].set_ylabel("Mean R\u00B2 \u2192 Higher is better")
for _, row in df.iterrows():
    axs[0, 0].text(row['Time_s'], row['Mean_R2'] + 0.0005, f" {row['Mode']}-{row['Label']}", fontsize=9)

sns.scatterplot(data=df, x='RAM_MB', y='Mean_R2', hue='Label', style='Mode', s=200, palette='tab10', ax=axs[0, 1])
axs[0, 1].set_title("Pareto: Accuracy vs Memory Usage", fontsize=14, fontweight='bold')
axs[0, 1].set_xlabel("Peak RAM (MB) \u2192 Lower is better")
axs[0, 1].set_ylabel("Mean R\u00B2 \u2192 Higher is better")
for _, row in df.iterrows():
    axs[0, 1].text(row['RAM_MB'], row['Mean_R2'] + 0.0005, f" {row['Mode']}-{row['Label']}", fontsize=9)

if not df_common.empty:
    sns.lineplot(data=df_common, x='MAF', y='Time_s', marker='o', markersize=8, color='blue', label='common mode', ax=axs[1, 0])
    axs[1, 0].set_xscale('log')
if not df_std.empty:
    axs[1, 0].axhline(df_std['Time_s'].values[0], color='red', linestyle='--', label='std baseline (Default)')
axs[1, 0].set_title("Performance Trend: Time vs Sparse MAF", fontsize=14, fontweight='bold')
axs[1, 0].set_xlabel("Sparse MAF Threshold (Log Scale)")
axs[1, 0].set_ylabel("Phasing Time (Seconds)")
axs[1, 0].legend()

if not df_common.empty:
    sns.lineplot(data=df_common, x='MAF', y='Mean_R2', marker='o', markersize=8, color='green', label='common mode', ax=axs[1, 1])
    axs[1, 1].set_xscale('log')
if not df_std.empty:
    axs[1, 1].axhline(df_std['Mean_R2'].values[0], color='red', linestyle='--', label='std baseline (Default)')
axs[1, 1].set_title("Accuracy Trend: Mean R\u00B2 vs Sparse MAF", fontsize=14, fontweight='bold')
axs[1, 1].set_xlabel("Sparse MAF Threshold (Log Scale)")
axs[1, 1].set_ylabel("Mean R\u00B2")
axs[1, 1].legend()

plt.tight_layout()
plt.savefig(snakemake.output.plot, dpi=300)
plt.close()


# --- 3. GRAFICO 2: DETTAGLIO PER BIN DI FREQUENZA ---
plt.figure(figsize=(12, 8))

# Disegniamo l'andamento R2 rispetto alla MAF media del bin
sns.lineplot(
    data=df_bins, 
    x='Bin_Avg_MAF', 
    y='R2', 
    hue='Threshold', 
    style='Mode', 
    palette='tab10',
    markers=True, 
    dashes=False, 
    linewidth=2, 
    markersize=8
)

plt.xscale('log') # La scala logaritmica sull'asse X per espandere la zona delle varianti rare
plt.title("Imputation Accuracy (R²) across Allele Frequency Spectrum", fontsize=16, fontweight='bold')
plt.xlabel("Bin Average MAF (Log Scale)", fontsize=12)
plt.ylabel("Concordance R²", fontsize=12)
plt.legend(title='Model - MAF Threshold', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.savefig(snakemake.output.plot_bins, dpi=300)
plt.close()
