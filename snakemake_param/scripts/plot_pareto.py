import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def parse_filename(filepath):
    filename = os.path.basename(filepath)
    name = filename.replace("bench_", "").replace("concordance_", "").replace(".tsv", "").split(".rsquare")[0]
    parts = name.split("_")
    mode = parts[0]
    min_thr = parts[1].replace("min", "")
    med_thr = parts[2].replace("med", "")
    max_val = parts[3].replace("max", "")
    label = f"{mode} | min:{min_thr} med:{med_thr} max:{max_val}"
    return mode, min_thr, med_thr, max_val, label

def parse_bench(filepath):
    mode, min_thr, med_thr, max_val, label = parse_filename(filepath)
    df = pd.read_csv(filepath, sep='\t')
    return {
        'Mode': mode, 'Min_Thr': min_thr, 'Med_Thr': med_thr, 'Max_Val': max_val, 
        'Label': label, 'Time_s': df['s'].iloc[0], 'RAM_MB': df['max_rss'].iloc[0]
    }

def parse_conc(filepath):
    mode, min_thr, med_thr, max_val, label = parse_filename(filepath)
    df = pd.read_csv(filepath, header=None, sep=r'\s+')
    df.rename(columns={0: 'Bin_ID', 1: 'N_variants', 2: 'Bin_Avg_MAF'}, inplace=True)
    df['R2'] = df.iloc[:, -1]
    
    mean_r2 = df['R2'].mean()
    rare_r2 = df.loc[df['Bin_ID'] == df['Bin_ID'].min(), 'R2'].values[0]
    
    summary = {
        'Mode': mode, 'Min_Thr': min_thr, 'Med_Thr': med_thr, 'Max_Val': max_val, 
        'Label': label, 'Mean_R2': mean_r2, 'Rare_R2': rare_r2
    }
    
    df['Mode'] = mode
    df['Configuration'] = label
    return summary, df[['Mode', 'Configuration', 'Bin_ID', 'N_variants', 'Bin_Avg_MAF', 'R2']]

def get_pareto_front(df, x_col, y_col):
    df_sorted = df.sort_values(by=x_col, ascending=True)
    front_x, front_y = [], []
    max_y = -float('inf')
    for _, row in df_sorted.iterrows():
        if row[y_col] >= max_y:
            front_x.append(row[x_col])
            front_y.append(row[y_col])
            max_y = row[y_col]
    return front_x, front_y

os.makedirs(os.path.dirname(snakemake.output.main), exist_ok=True)

bench_data = [parse_bench(f) for f in snakemake.input.benchmarks]
conc_results = [parse_conc(f) for f in snakemake.input.concordance]

conc_summary = [res[0] for res in conc_results]
df_bins = pd.concat([res[1] for res in conc_results], ignore_index=True)

df = pd.merge(pd.DataFrame(bench_data), pd.DataFrame(conc_summary), on=['Mode', 'Min_Thr', 'Med_Thr', 'Max_Val', 'Label'])
df.sort_values(by=['Max_Val', 'Min_Thr', 'Med_Thr'], inplace=True)

sns.set_theme(style="whitegrid")

fig, axs = plt.subplots(1, 2, figsize=(18, 8))
sns.scatterplot(data=df, x='Time_s', y='Mean_R2', hue='Max_Val', style='Mode', s=150, palette='Set1', ax=axs[0])
px, py = get_pareto_front(df, 'Time_s', 'Mean_R2')
axs[0].plot(px, py, color='black', linestyle='--', linewidth=2, label='Pareto Frontier')
axs[0].set_title("Main Pareto: Mean R\u00B2 vs Phasing Time", fontsize=14, fontweight='bold')
axs[0].set_xlabel("Phasing Time (Seconds) \u2192 Faster is better")
axs[0].set_ylabel("Mean R\u00B2 \u2192 Higher is better")
axs[0].legend()

sns.scatterplot(data=df, x='RAM_MB', y='Mean_R2', hue='Max_Val', style='Mode', s=150, palette='Set1', ax=axs[1])
px, py = get_pareto_front(df, 'RAM_MB', 'Mean_R2')
axs[1].plot(px, py, color='black', linestyle='--', linewidth=2, label='Pareto Frontier')
axs[1].set_title("Main Pareto: Mean R\u00B2 vs Memory Usage", fontsize=14, fontweight='bold')
axs[1].set_xlabel("Peak RAM (MB) \u2192 Lower is better")
axs[1].set_ylabel("Mean R\u00B2 \u2192 Higher is better")
axs[1].legend()

plt.tight_layout()
plt.savefig(snakemake.output.main, dpi=300)
plt.close()

fig, ax = plt.subplots(figsize=(10, 8))
sns.scatterplot(data=df, x='Time_s', y='Rare_R2', hue='Max_Val', style='Mode', s=200, palette='Set1', ax=ax)
px, py = get_pareto_front(df, 'Time_s', 'Rare_R2')
ax.plot(px, py, color='black', linestyle='--', linewidth=2, label='Pareto Frontier')
ax.set_title("Rare Variants Pareto: Accuracy (Lowest MAF Bin) vs Time", fontsize=16, fontweight='bold')
ax.set_xlabel("Phasing Time (Seconds) \u2192 Faster is better")
ax.set_ylabel("R\u00B2 (Rarest Variants) \u2192 Higher is better")
ax.legend()
plt.tight_layout()
plt.savefig(snakemake.output.rare, dpi=300)
plt.close()

fig, axs = plt.subplots(2, 3, figsize=(18, 10))
sns.boxplot(data=df, x='Min_Thr', y='Time_s', ax=axs[0, 0], palette='Blues')
axs[0, 0].set_title("Time by Min_Thr")
sns.boxplot(data=df, x='Med_Thr', y='Time_s', ax=axs[0, 1], palette='Greens')
axs[0, 1].set_title("Time by Med_Thr")
sns.boxplot(data=df, x='Max_Val', y='Time_s', ax=axs[0, 2], palette='Reds')
axs[0, 2].set_title("Time by Max_Val")

sns.boxplot(data=df, x='Min_Thr', y='Mean_R2', ax=axs[1, 0], palette='Blues')
axs[1, 0].set_title("Mean R\u00B2 by Min_Thr")
sns.boxplot(data=df, x='Med_Thr', y='Mean_R2', ax=axs[1, 1], palette='Greens')
axs[1, 1].set_title("Mean R\u00B2 by Med_Thr")
sns.boxplot(data=df, x='Max_Val', y='Mean_R2', ax=axs[1, 2], palette='Reds')
axs[1, 2].set_title("Mean R\u00B2 by Max_Val")
plt.tight_layout()
plt.savefig(snakemake.output.impact, dpi=300)
plt.close()

plt.figure(figsize=(14, 8))
sns.lineplot(data=df_bins, x='Bin_Avg_MAF', y='R2', hue='Configuration', linewidth=1.2, alpha=0.6)
plt.xscale('log')
plt.title("Imputation Accuracy (R\u00B2) across Allele Frequency Spectrum", fontsize=16, fontweight='bold')
plt.xlabel("Bin Average MAF (Log Scale)")
plt.ylabel("Concordance R\u00B2")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='xx-small')
plt.tight_layout()
plt.savefig(snakemake.output.bins, dpi=300)
plt.close()
