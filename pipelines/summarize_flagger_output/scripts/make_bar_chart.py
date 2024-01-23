import pandas as pd

all_df = pd.DataFrame()

import matplotlib.pyplot as plt

for sample, flagger_bed in snakemake.input.items():
    df = pd.read_csv(flagger_bed, sep='\t', usecols=[0,1,2,3], names=['#chrom', 'start', 'end', 'asm_type'])
    df['len'] = (df['end'] - df['start']) / 1e6
    out_df = df.groupby(['asm_type']).sum().reset_index()[['len', 'asm_type']].T
    out_df = out_df.rename(columns={col : f'{out_df[col].iloc[1]}' for col in out_df.columns}).iloc[0:1]
    out_df['sample'] = sample
    out_df = out_df.set_index('sample')
    all_df = pd.concat([all_df, out_df])

ASM_COLORS = {
    "Hap": (0, 0.541, 0),
    "Col": (0.667, 0, 1),
    "Dup": (0.980, 0.407, 0),
    "Unk": (0.502, 0.502, 0.502),
    "Err": (0.635, 0, 0.145),
}


all_df.plot(kind='bar', stacked=True, color=[ASM_COLORS[x] for x in all_df.columns])
plt.xlabel('Sample')
plt.xticks(rotation=45)
plt.ylabel('Flagger Bases (Mbp)')
plt.savefig(snakemake.output.tsv_all, bbox_inches='tight')


all_df[[x for x in all_df.columns if x != 'Hap']].plot(kind='bar', stacked=True, color=[ASM_COLORS[x] for x in all_df[[x for x in all_df.columns if x != 'Hap']].columns])
plt.xlabel('Sample')
plt.xticks(rotation=45)
plt.ylabel('Flagger Bases (Mbp)')
plt.savefig(snakemake.output.tsv_nohap, bbox_inches='tight')
