from pybedtools import BedTool
import pandas as pd
import os
import numpy as np
import sys
import argparse
from matplotlib import collections  as mc
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib as mpl


def get_df_fai(fai_file_name, usecols=('CHROM', 'LEN'), index_col='CHROM', squeeze=True, index_type=str):
    """
    Read an FAI File name. By default, return a Series of chromosome lengths keyed by the chromosome (or contig) name.

    :param fai_file_name: File to read.
    :param usecols: Use these columns. Must include index column if set.
    :param index_col: Index column name.
    :param squeeze: Squeeze to Series if only one column is selected (not counting the index column).
    :param index_type: Set index data type (dtype) to this type if not None.

    :return: A DataFrame or Series of the FAI file.
    """

    df = pd.read_csv(
        fai_file_name,
        sep='\t',
        names=['CHROM', 'LEN', 'POS', 'LINE_BP', 'LINE_BYTES'],
        usecols=usecols,
        dtype={'CHROM': str, 'LEN': int, 'POS': int, 'LINE_BP': int, 'LINE_BYTES': int}
    )

    if index_col is not None:
        df.set_index(index_col, inplace=True)

    if squeeze:
        df = df.squeeze(axis=1)

    if index_type is not None:
        df.index = df.index.astype(index_type)

    return df

"""
Generate ideogram plots.
"""

IDEO_HIST_PARAMS = {
    'x_adjust': 0.1,  # Adjust x in by this proportion to make space for y-axis labels
    'y_pad': 0.3,     # Add this proportion to the total height and distribute between plots as padding
    'band_color': {   # Color for bands keyed by "gieStain" is chromosome BED (from UCSC).
        'gneg': '0.875',
        'gpos100': '0.00',
        'gpos75': '0.16',
        'gpos50': '0.32',
        'gpos25': '0.48',
        'acen': 'darkred',
        'gvar': '0.00',
        'stalk': '0.64'
    },
    'chroms': [  # Chromosome order. This order minimizes dead space in multi-faceted plot with two columns (hg38)
        'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chrY',
        'chrX', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22'
    ],
    'anno_fig_height': 0.25,  # Proportion of original figure to expand below for annotation space
    'band_prop': 0.70,        # Proportion of annotation space for chromosome bands
    'gap_prop': 0.20,         # Proportion of annotation space for gap lines through chromosome bands
    'sd_color': 'goldenrod',       # SD overlay color
    'sd_alpha_range': (0.2, 0.8),  # The range of alpha values by SD identity are scaled within this range
    'tr_color': 'darkcyan',  # Tandem repeat overlay color
    'tr_alpha': 0.8,         # Tandem repeat overlay alpha
    'tr_min': 5000           # Tandem repeat minimum length
}

class IdeoHistogram:
    """
    Object returned from ideo_hist. Contains the figure and a dictionary of axes with tuple keys with the same
    indexes as a matrix of chromosome names (matrix_chr_name). This object makes it possible to write the figure as
    well as modify it (e.g. add additional annotations).
    """

    def __init__(self, fig, ax_dict, matrix_chr_name):
        self.fig = fig
        self.ax_dict = ax_dict
        self.matrix_chr_name = matrix_chr_name

    def close(self):
        """
        Close figure.
        """
        return plt.close(self.fig)

    def savefig(self, *args, **kwargs):
        return self.fig.savefig(*args, **kwargs)


def ideo_hist(
        df, fai_file_name,
        df_band, df_gap, df_sd=None, df_tr=None,
        width=10, height=12, dpi=300,
        label_col='SVTYPE',
        label_order=['INS', 'DEL', 'INV', 'SNV'],
        label_color={'INS': 'blue', 'DEL': 'red', 'INV': 'green', 'SNV': 'black'},
        plot_params=None,
        ideo_patch_z_start=100,
        ylim_dict=None,
        cb_func=None,
        verbose=False
    ):
    """
    Create an ideogram with a panel for each chromosome.

    :param df: DataFrame of values to plot. Must have "#CHROM", "POS", "END", and `label_col'. May be `None` if all
        plotting is handled by a callback function.
    :param fai_file_name: FAI file name.
    :param df_band: DataFrame of chromosome bands.
    :param df_gap: DataFrame of assembly gaps.
    :param df_sd: DataFrame of SDs with identity.
    :param df_tr: DataFrame of tandem repeat locations.
    :param width: Figure width in inches.
    :param height: Figure height in inches.
    :param dpi: Figure DPI.
    :param label_col: Labels to plot against. Stacked bar colors are separated on this column in `df`.
    :param label_order: Order of labels. Labels not in this list are ignored.
    :param label_color: Label color (color of stacked bars). May be a list of the same length as `label_order` or a
        dict of values with `label_order` items as keys. Must be valid matplotlib colors.
    :param plot_params: Plot parameters as a dictionary. `IDEO_HIST_PARAMS` is the default set of parameters.
    :param ideo_patch_z_start: Relative order of patches (chromosome ideo pieces). Set to 100 by default should make
        it appear behind other plot elements if they overlap.
    :param ylim_dict: Optional y-limits as a dict keyed on chromosomes.
    :param cb_func: If not `None`, call this function for each chromosome. Function signature is
        `cb_func(df, chrom, ax, fig)` where `df` is the DataFrame given to this function (not subset for the
        chromosome), `chrom` is the chromosome name, 'ax` is the matplotlib axes object for this chromosome, and `fig`
        is the matplotlib figure object.
    :param verbose: If `True`, print verbose information while plotting.

    :return: A `IdeoHistogram` object with the figure, a dict of axes, and a matrix of chromosome names.
    """

    # Init parameters
    if plot_params is None:
        plot_params = IDEO_HIST_PARAMS

    plot_params = plot_params.copy()

    chroms = plot_params['chroms']

    if ylim_dict is None:
        ylim_dict = dict()

    if issubclass(label_color.__class__, dict):
        label_color = [label_color[val] for val in label_order]

    ### Read ###

    # Get variant midpoints
    if df is not None:
        df = df[['#CHROM', 'POS', 'END', label_col]].copy()

        df['POS'] = (df['POS'] + df['END']) // 2

    # Read FAI
    df_fai = get_df_fai(fai_file_name)

    if df_sd is not None:
        df_sd = df_sd.copy()

        # Scale identity to span SD_ALPHA_RANGE
        sd_min = np.min(df_sd['MATCH'])

        df_sd['ALPHA'] = (df_sd['MATCH'] - sd_min) / (1 - sd_min) * (plot_params['sd_alpha_range'][1] - plot_params['sd_alpha_range'][0]) + plot_params['sd_alpha_range'][0]

    # Read TR
    if df_tr is not None:
        df_tr = df_tr.loc[(df_tr['END'] - df_tr['POS']) >= plot_params['tr_min']].copy()

    ### Assign chrom matrix ###
    # Note: assumes len(chroms) is an even number (update to handle more general cases)

    # Get names
    matrix_chr_name = np.array(
        [
            chroms[:len(chroms) // 2],
            chroms[-1:len(chroms) // 2 - 1:-1]
        ]
    ).T

    # Get chromosome lengths
    matrix_chr_len = np.array(
        list(map(lambda val: df_fai[val], matrix_chr_name))
    )

    # Find max horizontal pair
    max_pair_bp = np.max(
        np.apply_along_axis(np.sum, 1, matrix_chr_len)
    ) * (
        (1 + plot_params['x_adjust'] * 2)
    )

    # Get height per subplot
    subplot_height = 1 / matrix_chr_name.shape[0]  # Proportion for each row
    subplot_height *= (1 - plot_params['y_pad'])  # Shrink by total pad height

    pad_height = plot_params['y_pad'] / (matrix_chr_name.shape[0])  # Pad between each row

    ### Plot ###

    fig = plt.figure(figsize=(width, height), dpi=dpi)

    ax_dict = dict()

    for i in range(matrix_chr_name.shape[0]):
        for j in range(matrix_chr_name.shape[1]):

            chrom = matrix_chr_name[i, j]

            # Make axes
            ax_len = matrix_chr_len[i, j] / max_pair_bp

            ax = fig.add_axes((
                plot_params['x_adjust'] if j == 0 else 1 - ax_len,  # Left align if column 0, right align if column 1
                1 - (subplot_height * (i + 1) + pad_height * i),
                ax_len,
                subplot_height
            ))

            ax_dict[(i, j)] = ax

            # Histogram
            if df is not None:
                if verbose:
                    print('SVTYPE hist: {}'.format(chrom))

                ax.hist(
                    [
                        df.loc[
                            (df['#CHROM'] == chrom) & (df[label_col] == label),
                            'POS'
                        ] for label in label_order
                    ],
                    histtype='bar',
                    stacked=True,
                    bins=int(matrix_chr_len[i, j] // 1e6),
                    color=label_color,
                    label=label_order
                )

            # Callback plot function
            if cb_func is not None:
                if verbose:
                    print('Callback: {}'.format(chrom))

                cb_func(df, chrom, ax, fig)

            if chrom in ylim_dict:
                ax.set_ylim(ylim_dict[chrom])

            # Scale x to Mbp
            ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x, pos: str(int(x // 1e6))))

            # Adjust x-axis (0 to chromosome max)
            ax.set_xlim((-(df_fai[chrom] * 0.01), df_fai[chrom] * 1.01))

            # Remove spines
            for loc in ('top', 'right'):
                ax.spines[loc].set_visible(False)

            # Title
            ax.text(
                0.5, 0.95,
                matrix_chr_name[i, j],
                horizontalalignment='center',
                verticalalignment='top',
                transform=ax.transAxes,
                bbox={'boxstyle': 'round', 'facecolor': 'white', 'alpha': 0.5}
            )

            # Aestetics
            ax.tick_params(axis='both', which='major', labelsize=8)

    # Add bands, gaps, SDs, and TRs
    if verbose:
        print('Adding ideo bands')

    patch_z = ideo_patch_z_start

    for i in range(matrix_chr_name.shape[0]):
        for j in range(matrix_chr_name.shape[1]):

            ax = ax_dict[(i, j)]
            chrom = matrix_chr_name[i, j]

            y_min, y_max = ax.get_ylim()

            # Create annotation space below histogram
            anno_height = (y_max - y_min) * plot_params['anno_fig_height']

            ax.set_ylim(
                np.min([y_min, -anno_height]),
                np.max([y_max, anno_height])
            )

            # Set positions for band annotation
            band_y_pos = - anno_height * ((1 - plot_params['band_prop']) / 2)
            band_y_height = anno_height * plot_params['band_prop']

            # Set positions for gap annotations
            gap_y_pos = - anno_height * ((1 - plot_params['gap_prop']) / 2)
            gap_y_height = anno_height * plot_params['gap_prop']

            # Add Bands
            for index, row in df_band.loc[df_band['#chrom'] == chrom].iterrows():

                if row['gieStain'] == 'acen':
                    continue  # Skip centromeres (plot over all other annotations)

                ax.add_patch(
                    mpl.patches.Rectangle(
                        (row['start'], band_y_pos),
                        row['end'] - row['start'],
                        -band_y_height,
                        facecolor=plot_params['band_color'][row['gieStain']],
                        zorder=patch_z
                    )
                )

                patch_z += 1

            # Add Gaps
            for index, row in df_gap.loc[df_gap['#CHROM'] == chrom].iterrows():
                ax.add_patch(
                    mpl.patches.Rectangle(
                        (row['START'], gap_y_pos),
                        row['END'] - row['START'],
                        -gap_y_height,
                        facecolor='black',
                        zorder=patch_z
                    )
                )

                patch_z += 1

            # Add TR
            if df_tr is not None:
                for index, row in df_tr.loc[df_tr['#CHROM'] == chrom].iterrows():
                    ax.add_patch(
                        mpl.patches.Rectangle(
                            (row['POS'], band_y_pos),
                            row['END'] - row['POS'],
                            -band_y_height,
                            facecolor=plot_params['tr_color'],
                            alpha=plot_params['tr_alpha'],
                            zorder=patch_z
                        )
                    )

                    patch_z += 1

            # Add SD
            if df_sd is not None:
                for index, row in df_sd.loc[df_sd['#CHROM'] == chrom].iterrows():
                    ax.add_patch(
                        mpl.patches.Rectangle(
                            (row['POS'], band_y_pos),
                            row['END'] - row['POS'],
                            -band_y_height,
                            facecolor=plot_params['sd_color'],
                            alpha=row['ALPHA'],
                            zorder=patch_z
                        )
                    )

                    patch_z += 1

            # Add centromere bands
            for index, row in df_band.loc[(df_band['#chrom'] == chrom) & (df_band['gieStain'] == 'acen')].iterrows():
                ax.add_patch(
                    mpl.patches.Rectangle(
                        (row['start'], band_y_pos),
                        row['end'] - row['start'],
                        -band_y_height,
                        facecolor=plot_params['band_color'][row['gieStain']],
                        zorder=patch_z
                    )
                )

                patch_z += 1

    # Return
    return IdeoHistogram(fig, ax_dict, matrix_chr_name)





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


def ideo_cb(df, chrom, ax, fig):
    # SVPop removes all columns except #CHROM, POS, END,
    # and label_col from the df and divides bins in half, use original df

    mpl.rcParams.update({"font.size": 8})
    asms = sorted(df["ASM"].unique(), reverse=True)

    max_bar_height = len(samples) * len(asms)
    df = all_samples_df

    # Subset to chromosome
    df_chrom = df.loc[df["#CHROM"] == chrom].copy()
    for start_pos in df_chrom["POS"].unique():
        start_height = 0
        for j, asm in enumerate(asms):
            n_samples = len(
                df_chrom[(df_chrom["POS"] == start_pos) & (df_chrom["ASM"] == asm)]
            )
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
    ax.set_ylabel("# Samples")
    ax.set_yticks([0, 1])
    ax.set_yticklabels([0, max_bar_height])

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
    df['LEN'] = df["END"] - df['POS']
    df = df.loc[df['LEN'] >= int(snakemake.wildcards.filt)].copy()
    if flagger_category == "all":
        iter_range = df["ASM"].unique()
    elif flagger_category == "nohap":
        iter_range = [x for x in df["ASM"].unique() if x != "Hap"]
    else:
        iter_range = [flagger_category]
    for asm in iter_range:
        if not asm:
            continue
        asm_bed = BedTool.from_dataframe(df.loc[df["ASM"] == asm])
        asm_df = bins_bed.intersect(asm_bed, wa=True, u=True).to_dataframe()

        if asm_df.empty:
            continue
        asm_df.columns = ["#CHROM", "POS", "END"]
        asm_df["ASM"] = asm
        asm_df["SAMPLE"] = sample
        all_samples_df = pd.concat([all_samples_df, asm_df])

print(f"Making plot for Flagger category {flagger_category}")
# If SAMPLE not in all_samples_df, converted_paf may be empty
samples = sorted(all_samples_df["SAMPLE"].unique())
ideo_hist = ideo_hist(
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
