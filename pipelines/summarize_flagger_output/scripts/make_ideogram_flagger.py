import pandas as pd
from pybedtools import BedTool
import pandas as pd
import os
import numpy as np
import sys
import argparse
import pylab as pl
import matplotlib as mpl
from matplotlib import collections  as mc
sys.path.append('/net/eichler/vol26/7200/software/pipelines/svpop/svpop-3.4.0/dep')
sys.path.append('/net/eichler/vol26/7200/software/pipelines/svpop/svpop-3.4.0')
import svpoplib

OUT_PLOTS_DIR = sys.argv[1]
FLAGGER_RESULTS_FILES = sys.argv[2:]

#FLAGGER_RESULTS_DIR = '/net/eichler/vol26/projects/flagger_tmp/nobackups/saffire/CHM13_v2.0/visualization_pipeline/converted_paf' # Beds/ PAFs with lifted over coordinates

BIN_SIZE = np.int32(1e6)

#SPACER_PROP = 0.325  # Shift lower bars by this amount to make space for ideo

#LABEL_SPACE = 0.325  # Add this proportion of the y range to the upper limit to make space for the chromosome label

FAI_FILE_NAME = '/net/eichler/vol28/eee_shared/assemblies/CHM13/T2T/v2.0/T2T-CHM13v2.fasta.fai'

BAND_FILE_NAME = '/net/eichler/vol28/eee_shared/assemblies/CHM13/T2T/v2.0/anno/cyto.bed'
GAP_FILE_NAME = '/net/eichler/vol28/eee_shared/assemblies/CHM13/T2T/v2.0/anno/T2T-CHM13v2_gap.bed'
SD_FILE_NAME = '/net/eichler/vol28/eee_shared/assemblies/CHM13/T2T/v2.0/flagger/chm13v2.0.sd_frac_match.bed'
TR_FILE_NAME = '/net/eichler/vol28/eee_shared/assemblies/CHM13/T2T/v2.0/anno/chm13_t2t_trf.bed'

#BINS_BED = BedTool('/net/eichler/vol26/projects/flagger_tmp/nobackups/saffire/CHM13_v2.0/visualization_pipeline/windows_chr22.bed')
BINS_BED = BedTool('/net/eichler/vol26/projects/flagger_tmp/nobackups/saffire/CHM13_v2.0/visualization_pipeline/windows.bed')

ASM_COLORS = {
    'Hap': (0, 0.541, 0),
    'Col': (0.667, 0, 1),
    'Dup': (0.980, 0.407, 0),
    'Unk': (0.502, 0.502, 0.502),
    'Err': (0.635, 0, 0.145)
}

df_band = pd.read_csv(BAND_FILE_NAME, sep='\t')
df_gap = pd.read_csv(GAP_FILE_NAME, sep='\t')
df_sd = pd.read_csv(SD_FILE_NAME, sep='\t')
df_tr = pd.read_csv(
    TR_FILE_NAME, sep='\t', header=None, names=('#CHROM', 'POS', 'END')
)

def ideo_cb(df, chrom, ax, fig):
    # SVPop removes all columns except #CHROM, POS, END,
    # and label_col from the df and divides bins in half, use original df
    mpl.rcParams.update({'font.size':8})
    asms = df['ASM'].unique()
    max_bar_height = len(samples) * len(asms)
    df = asm_dict[asms[0]] if (len(asms) == 1) else asm_dict_all

    # Subset to chromosome
    df_chrom = df.loc[df['#CHROM'] == chrom].copy()
    for start_pos in df_chrom['POS'].unique():
        start_height = 0
        for j, asm in enumerate(asms):
            n_samples = len(
                df_chrom[(df_chrom['POS'] == start_pos) & (df_chrom['ASM'] == asm)]
            )
            sample_height = n_samples / max_bar_height
            b = ax.bar(
                start_pos, width=BIN_SIZE * .8, height=sample_height,
                bottom=start_height, align='edge', label=asm,
                color=ASM_COLORS[asm]
            )
            start_height += sample_height # Stack bars for each Flag

    # Need small y_max in order for ideogram to show, but set y labels as
    # actual number of samples
    ax.set_ylim(0, 1)
    ax.set_ylabel('# Samples')
    ax.set_yticks([0, 1])
    ax.set_yticklabels([0, max_bar_height])

    # Show legend
    if len(asms) > 1 and chrom == 'chrY':
        h, l = ax.get_legend_handles_labels()
        d = {l[i]:h[i] for i in range(len(h))}
        asms_for_legend = list(asms)[::-1]
        ax.legend(
            labels=asms_for_legend, handles=[d[x] for x in asms_for_legend], 
            loc='center left', bbox_to_anchor=(-1.5, 0.35)
        )

# asm_dict: All binned blocks (#CHROM, POS, END) for each Flag (Err, Dup, Hap, Col, Unk).
# Bedtools intersect each Flag of each sample w coord windows
# and concat to the overall df for the Flag
asm_dict = {}
for infile in FLAGGER_RESULTS_FILES:
    print(f'Reading in {infile}')
    assert(infile.lower().endswith('.paf'))
    #if infile not in ['GM18989.paf', 'GM19129.paf', 'GM19331.paf']:
    #    continue
    df = pd.read_csv(
        infile, sep='\t', header=None,
        usecols=[0, 2, 3, 12], names=['#CHROM', 'POS', 'END', 'ASM']
    )
    df['ASM'] = df['ASM'].apply(lambda x: x.split(':')[-1].split('_')[0])
    for asm in df['ASM'].unique():
        if not asm:
            continue
        #if asm != 'Dup':
        #    continue
        if asm not in asm_dict:
            asm_dict[asm] = pd.DataFrame()
        asm_bed = BedTool.from_dataframe(df.loc[df['ASM'] == asm])
        asm_df = BINS_BED.intersect(asm_bed, wa=True, u=True).to_dataframe()
        if asm_df.empty:
            continue
        asm_df.columns=['#CHROM', 'POS', 'END']
        asm_df['ASM'] = asm
        asm_df['SAMPLE'] = os.path.basename(infile).split('.')[0]
        asm_dict[asm] = pd.concat([asm_dict[asm], asm_df])

print('Making per-Flag plots')
samples = sorted(asm_dict['Hap']['SAMPLE'].unique())

asm_dict_all = pd.DataFrame()
for asm in asm_dict:
    asm_dict_all = pd.concat([asm_dict_all, asm_dict[asm]])
    ideo_hist = svpoplib.plot.ideo.ideo_hist(
        asm_dict[asm], FAI_FILE_NAME, df_band, df_gap, df_sd,
        df_tr, cb_func=ideo_cb, label_col='ASM'
    )
    ideo_hist.fig.savefig(f'{OUT_PLOTS_DIR}/{asm}.png', bbox_inches='tight')
    ideo_hist.fig.savefig(f'{OUT_PLOTS_DIR}/{asm}.svg', bbox_inches='tight')

print('Making overall plot')
ideo_hist = svpoplib.plot.ideo.ideo_hist(
    asm_dict_all, FAI_FILE_NAME, df_band, df_gap, df_sd, df_tr,
    cb_func=ideo_cb, label_col='ASM'
)
ideo_hist.fig.savefig(f'{OUT_PLOTS_DIR}/all.png', bbox_inches='tight')
ideo_hist.fig.savefig(f'{OUT_PLOTS_DIR}/all.svg', bbox_inches='tight')



