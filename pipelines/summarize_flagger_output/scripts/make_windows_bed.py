import sys
import glob

out_fname = "windows.bed"
BIN_SIZE = 1000000

chrom_lengths_chm13_t2t = {
    "chr1": 248387328,
    "chr2": 242696752,
    "chr3": 201105948,
    "chr4": 193574945,
    "chr5": 182045439,
    "chr6": 172126628,
    "chr7": 160567428,
    "chr8": 146259331,
    "chr9": 150617247,
    "chr10": 134758134,
    "chr11": 135127769,
    "chr12": 133324548,
    "chr13": 113566686,
    "chr14": 101161492,
    "chr15": 99753195,
    "chr16": 96330374,
    "chr17": 84276897,
    "chr18": 80542538,
    "chr19": 61707364,
    "chr20": 66210255,
    "chr21": 45090682,
    "chr22": 51324926,
    "chrX": 154259566,
    "chrY": 62460029,
}

with open(out_fname, "w") as out_f:
    for chrom in chrom_lengths_chm13_t2t:
        out_intervals = list(
            range(0, chrom_lengths_chm13_t2t[chrom] + BIN_SIZE, BIN_SIZE)
        )
        for i in range(len(out_intervals) - 1):
            out_f.write(
                "\t".join([chrom, str(out_intervals[i]), str(out_intervals[i + 1])])
                + "\n"
            )
