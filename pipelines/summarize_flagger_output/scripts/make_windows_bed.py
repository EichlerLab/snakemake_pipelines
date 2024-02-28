import sys
import glob
import csv

BIN_SIZE = 1000000

out_fname = snakemake.output.bed
fai_filename = snakemake.input.fai

with open(out_fname, "w") as out_f:
    with open(fai_filename) as f:
        r = csv.DictReader(f, delimiter="\t", fieldnames=["#CHROM", "LENGTH"])
        for row in r:
            out_intervals = list(range(0, int(row["LENGTH"]) + BIN_SIZE, BIN_SIZE))
            for i in range(len(out_intervals) - 1):
                out_f.write(
                    "\t".join(
                        [
                            row["#CHROM"],
                            str(out_intervals[i]),
                            str(out_intervals[i + 1]),
                        ]
                    )
                    + "\n"
                )
