import pandas as pd
from pybedtools import BedTool
import pandas as pd
import os
import numpy as np
import sys
import pylab as pl
import matplotlib as mpl
from matplotlib import collections as mc

sys.path.append("/net/eichler/vol26/7200/software/pipelines/svpop/svpop-3.4.0/dep")
sys.path.append("/net/eichler/vol26/7200/software/pipelines/svpop/svpop-3.4.0")
import svpoplib

BIN_SIZE = np.int32(1e6)

ASM_COLORS = {
    "Hap": (0, 0.541, 0),
    "Col": (0.667, 0, 1),
    "Dup": (0.980, 0.407, 0),
    "Unk": (0.502, 0.502, 0.502),
    "Err": (0.635, 0, 0.145),
}

fai_filename = snakemake.input.fai
df_band = pd.read_csv(snakemake.input.chrom_bands, sep="\t")
df_gap = pd.read_csv(snakemake.input.gaps, sep="\t")
df_sd = pd.read_csv(snakemake.input.sd, sep="\t")
df_tr = pd.read_csv(
    snakemake.input.tandem_repeats,
    sep="\t",
    header=None,
    names=("#CHROM", "POS", "END"),
)

flagger_category = snakemake.wildcards.flag
infiles = snakemake.input.paf
bins_bed = BedTool(snakemake.input.bed)
outfiles = snakemake.output
plot_bed = BedTool(snakemake.input.plot_regions) if (
    snakemake.input.plot_regions
) else None

BAR_HEIGHT_NUM_BASES = plot_bed is not None

def ideo_cb(df, chrom, ax, fig):
    # SVPop removes all columns except #CHROM, POS, END,
    # and label_col from the df and divides bins in half, use original df

    mpl.rcParams.update({"font.size": 8})
    asms = sorted(df["ASM"].unique(), reverse=True)
    df = all_samples_df

    if not BAR_HEIGHT_NUM_BASES:
        max_bar_height = len(samples) * len(asms)

    # Subset to chromosome
    df_chrom = df.loc[df["#CHROM"] == chrom].copy()
    if BAR_HEIGHT_NUM_BASES:
        assert len(samples) == 1
        max_bar_height = df_chrom.groupby(['ASM', 'POS']).agg({"OVERLAP":"sum"}).reset_index()["OVERLAP"].max()
    for start_pos in df_chrom["POS"].unique():
        start_height = 0
        for j, asm in enumerate(asms):
            curr_rows = df_chrom[(df_chrom["POS"] == start_pos) & (df_chrom["ASM"] == asm)]
            if BAR_HEIGHT_NUM_BASES:
                n_bases_overlap = curr_rows["OVERLAP"].max()
                sample_height = n_bases_overlap / max_bar_height
            else:
                n_samples = len(curr_rows)
                sample_height = n_samples / max_bar_height
            b = ax.bar(
                start_pos,
                width=BIN_SIZE * 0.8,
                height=sample_height,
                bottom=start_height,
                align="edge",
                label=asm,
                color=ASM_COLORS[asm],
            )
            start_height += sample_height  # Stack bars for each Flag

    # Need small y_max in order for ideogram to show, but set y labels as
    # actual number of samples
    ax.set_ylim(0, 1)
    ax.set_ylabel("# Kilobases" if (BAR_HEIGHT_NUM_BASES) else '# Samples')
    ax.set_yticks([0, 1])
    y_max = 'N/A'
    if max_bar_height and not pd.isna(max_bar_height):
        y_max = f'{int(max_bar_height/1000.0):,}' if (BAR_HEIGHT_NUM_BASES) else max_bar_height
    ax.set_yticklabels([0, y_max])

    # Show legend
    if len(asms) > 1:
        h, l = ax.get_legend_handles_labels()
        legend_d.update({l[i]: h[i] for i in range(len(h))})

        if chrom == "chrY":
            asms_for_legend = asms[::-1]
            ax.legend(
                labels=asms_for_legend,
                handles=[legend_d[x] for x in asms_for_legend],
                loc="center left",
                bbox_to_anchor=(-1.5, 0.35),
            )


# asm_dict: All binned blocks (#CHROM, POS, END) for each Flag (Err, Dup, Hap, Col, Unk).
# Bedtools intersect each Flag of each sample w coord windows
# and concat to the overall df for the Flag
all_samples_df = pd.DataFrame()
legend_d = {}

for infile in infiles:
    # print(f'Reading in {infile}')
    sample = os.path.basename(infile).split(".")[0]  # Since PAFs are named {sample}.paf
    print(f"Sample {sample}")

    df = pd.read_csv(
        infile,
        sep="\t",
        header=None,
        usecols=[0, 2, 3, 12],
        names=["#CHROM", "POS", "END", "ASM"],
    )
    df["ASM"] = df["ASM"].apply(lambda x: x.split(":")[-1].split("_")[0])
    if flagger_category == "all":
        iter_range = df["ASM"].unique()
    else:
        iter_range = [flagger_category]
    for asm in iter_range:
        if not asm:
            continue
        asm_bed = BedTool.from_dataframe(df.loc[df["ASM"] == asm])
        if BAR_HEIGHT_NUM_BASES:
            asm_bed_filtered = asm_bed.intersect(plot_bed)
            asm_df = bins_bed.intersect(asm_bed_filtered, wo=True).to_dataframe(names=[
                '#CHROM', 'POS', 'END', '#CHROM_B', 'POS_B', 'END_B',
                'ASM', 'OVERLAP'
            ])
            asm_df = asm_df.groupby([
                '#CHROM', 'POS', 'END'
            ]).agg({'OVERLAP':'sum'}).reset_index()
        else:
            asm_df = bins_bed.intersect(asm_bed, wa=True, u=True).to_dataframe()
            asm_df.columns = ["#CHROM", "POS", "END"]
        asm_df["ASM"] = asm
        asm_df["SAMPLE"] = sample

        if asm_df.empty:
            continue

        all_samples_df = pd.concat([all_samples_df, asm_df])
if all_samples_df.empty:
    sys.exit()

print(f"Making plot for Flagger category {flagger_category}")
samples = sorted(all_samples_df["SAMPLE"].unique())
ideo_hist = svpoplib.plot.ideo.ideo_hist(
    all_samples_df,
    fai_filename,
    df_band,
    df_gap,
    df_sd,
    df_tr,
    cb_func=ideo_cb,
    label_col="ASM",
)
ideo_hist.fig.savefig(outfiles.png, bbox_inches="tight")
ideo_hist.fig.savefig(outfiles.svg, bbox_inches="tight")
