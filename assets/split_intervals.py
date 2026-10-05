#!/usr/bin/env python3
"""
Reimplementation of the getzlab split_intervals.py CLI contract, inferred
from how MuTect1_Scatter_Gather.wdl's `split_intervals` task invokes it:

    split_intervals.py -bam T.bam -bai T.bam.bai -interval_type picard \
        -N 200 [-target_list targets.bed] -chrs chr1,chr2,...,chrX,chrY

Produces <interval_type>/<NNNN>-scattered.interval_list -- one Picard-style
interval_list file per shard, each with a valid @HD/@SQ header copied from
the BAM's own sequence dictionary, followed by tab-separated interval rows,
covering a disjoint, roughly equal-sized (by total bp) slice of the requested
chromosomes (optionally restricted to -target_list first).

This is NOT a port of the lab's actual script -- nobody outside the lab has
its source. It's a from-scratch reimplementation of the same job (partition a
genome/BED into N roughly-equal interval shards), written to satisfy the same
CLI contract and output shape. Shard *boundaries* will not match the
original tool's byte-for-byte, which is fine: nothing downstream (MuTect1,
the gather step) depends on getting the same boundaries, only on the shards
being a correct, non-overlapping, complete partition of the requested region.
"""
import argparse
import os
import sys

try:
    import pysam
except ImportError:
    print("ERROR: this script requires pysam (pip install pysam)", file=sys.stderr)
    raise


def parse_bed_or_interval_list(path):
    """Minimal BED / Picard interval_list reader -> list of (chrom, start, end), 1-based inclusive."""
    intervals = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("@") or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0]
            # Picard interval_list is 1-based inclusive and has >=5 columns
            # (chrom, start, end, strand, name); BED is 0-based half-open.
            if len(fields) >= 5:
                start, end = int(fields[1]), int(fields[2])
            else:
                start, end = int(fields[1]) + 1, int(fields[2])
            intervals.append((chrom, start, end))
    return intervals


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-bam", required=True)
    p.add_argument("-bai", required=True)
    p.add_argument("-interval_type", required=True, choices=["picard"])
    p.add_argument("-N", type=int, required=True)
    p.add_argument("-target_list", default=None)
    p.add_argument("-chrs", required=True, help="comma-separated chromosome names")
    p.add_argument("-padding", type=int, default=0, help="base pairs to pad each interval on both sides")
    args = p.parse_args()

    selected_chrs = args.chrs.split(",")

    bam = pysam.AlignmentFile(args.bam, "rb")
    contig_lengths = dict(zip(bam.references, bam.lengths))
    contig_order = {sn: i for i, sn in enumerate(bam.references)}

    missing = [c for c in selected_chrs if c not in contig_lengths]
    if missing:
        print(f"WARNING: chromosomes not found in BAM header, skipping: {missing}", file=sys.stderr)
    selected_chrs = [c for c in selected_chrs if c in contig_lengths]

    header_lines = ["@HD\tVN:1.6\tSO:coordinate"] + [
        f"@SQ\tSN:{sn}\tLN:{contig_lengths[sn]}" for sn in selected_chrs
    ]

    if args.target_list:
        intervals = [iv for iv in parse_bed_or_interval_list(args.target_list) if iv[0] in selected_chrs]
        if not intervals:
            print("ERROR: -target_list produced zero intervals after restricting to -chrs", file=sys.stderr)
            sys.exit(1)
    else:
        intervals = [(chrom, 1, contig_lengths[chrom]) for chrom in selected_chrs]

    intervals.sort(key=lambda iv: (contig_order[iv[0]], iv[1]))

    total_bp = sum(end - start + 1 for _, start, end in intervals)
    if total_bp == 0:
        print("ERROR: zero total bp across selected intervals", file=sys.stderr)
        sys.exit(1)
    target_bp_per_shard = max(1, total_bp // args.N)

    outdir = args.interval_type
    os.makedirs(outdir, exist_ok=True)

    state = {"shard_idx": 0, "current_bp": 0, "current_rows": []}

    def flush():
        if not state["current_rows"]:
            return
        idx = state["shard_idx"]
        out_path = os.path.join(outdir, f"{idx:04d}-scattered.interval_list")
        with open(out_path, "w") as out:
            out.write("\n".join(header_lines) + "\n")
            for chrom, start, end in state["current_rows"]:
                padded_start = max(1, start - args.padding)
                padded_end = min(contig_lengths[chrom], end + args.padding)
                out.write(f"{chrom}\t{padded_start}\t{padded_end}\t+\tshard_{idx:04d}\n")
        state["shard_idx"] += 1
        state["current_bp"] = 0
        state["current_rows"] = []

    for chrom, start, end in intervals:
        pos = start
        while pos <= end:
            # once we're on the last allowed shard, absorb everything
            # remaining instead of trying to keep chunking further
            if state["shard_idx"] == args.N - 1:
                state["current_rows"].append((chrom, pos, end))
                state["current_bp"] += end - pos + 1
                pos = end + 1
                continue

            space_left = target_bp_per_shard - state["current_bp"]
            if space_left <= 0:
                flush()
                space_left = target_bp_per_shard

            chunk_end = min(end, pos + space_left - 1)
            state["current_rows"].append((chrom, pos, chunk_end))
            state["current_bp"] += chunk_end - pos + 1
            pos = chunk_end + 1

    flush()

    if state["shard_idx"] == 0:
        print("ERROR: produced zero shards", file=sys.stderr)
        sys.exit(1)

    print(f"Wrote {state['shard_idx']} interval_list file(s) to {outdir}/", flush=True)


if __name__ == "__main__":
    main()
