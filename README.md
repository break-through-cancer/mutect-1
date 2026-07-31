# Combined MuTect1 pipeline (NeoDisc + getzlab), packaged for Cirro

One Nextflow DSL2 workflow merging NeoDisc's `Mutect_v1_calling.sh` (bash,
on-prem/HPC) and the getzlab `MuTect1_Scatter_Gather.wdl` (Cromwell,
cloud-native), rebuilt to match the conventions of an existing production
Cirro Mutect2 pipeline rather than a generic Nextflow layout.

## NeoDisc vs. getzlab: similarities and differences

**Similar:** both run the same core tool (GATK3 MuTect1, tumor+germline BAM
in, VCF + call_stats out); both use a scatter/gather architecture; both
apply MuTect1's own PASS/REJECT call as the baseline filter on the merged
output (NeoDisc: "drop REJECT", getzlab: "keep PASS" — the same filter);
both validate defensively before spending compute.

**Different:** NeoDisc hand-rolls concurrency with a bash polling loop;
getzlab uses Cromwell's native scatter with cloud retries. NeoDisc splits by
a fixed 24-way per-chromosome list from a pre-built interval directory;
getzlab dynamically splits into N shards (default 200) from the BAM header.
NeoDisc never measures contamination; getzlab runs ContEst and feeds the
result into every shard. getzlab's gather step also re-sorts, stamps sample
names into the VCF header, and produces a DeTiN-prep call_stats subset that
NeoDisc doesn't have. NeoDisc treats dbSNP/COSMIC as mandatory; getzlab
treats them as optional. **Not included:** NeoDisc's multi-caller
`MinOverlap` ensemble filter — that's a separate downstream stage needing
HaplotypeCaller/Mutect2/Varscan2 output too, not part of "the MuTect1
workflow" in either source pipeline.

## Why this version looks different from the earlier draft

An earlier draft of this pipeline used a Nextflow-side sample-sheet CSV
read via `splitCsv()`, with a `bin/prepare_samplesheet.py` validator. That
doesn't match how Cirro actually works, and doesn't match your own working
Mutect2 pipeline: Cirro discovers samples from the dataset's file manifest
*before* Nextflow ever runs, via a `cirro.helpers.preprocess_dataset.
PreprocessDataset` hook that injects a run list straight into `params`.
This version follows that pattern instead — see `preprocess.py`.

Two decisions carried over directly from your Mutect2 pipeline (confirmed,
not assumed):

- **Shared-normal model.** One normal BAM, shared across every tumor sample
  in the dataset — not one normal per tumor pair. This removes any need for
  tumor-only/paired branching: every run has a normal by construction, so
  `normal_bam`/`normal_bai` are broadcast values like `refFasta`, never
  something requiring a per-row join.
- **Sample names derived from BAM header SM tags** at runtime (the same
  `samtools view -H | awk` idiom your `mutect_wrapper` process uses), not
  threaded through as separate params.

Also carried over: `label 'process_low'/'process_medium'/'process_high'`
resource tiers (mapped in `nextflow.config`) instead of literal `cpus`/
`memory` per process, `container "${params.x ?: 'default:tag'}"` fallback
style, `test -s <file>` assertions after every tool invocation, and
`params.containsKey(...)` guards on every default so a top-level assignment
never clobbers a CLI/`-params-file` override — that last one was a real bug
in an earlier draft of this file, fixed by matching your `mutect_runs`
preprocessing script's own guard pattern once I re-read it closely.

## What's different from your Mutect2 pipeline

MuTect1 shards are given the full tumor/normal BAM restricted with `-L
<interval_shard>`, rather than physically subsetting each BAM per-shard via
`samtools view` the way your `subset_tumor_all_shards`/
`subset_shared_normal_all_shards` processes do for Mutect2. That subsetting
step is a reasonable I/O optimization if this pipeline turns out to need it
at scale, but wasn't present in either source pipeline (NeoDisc or getzlab)
and adds real complexity, so it's left out for now — flagging it explicitly
rather than silently deciding for you.

## Layout

```
main.nf          all processes + the workflow block
nextflow.config  container/resource-label defaults, local/awsbatch profiles
preprocess.py    Cirro dataset preprocessing hook -- discovers tumor BAMs +
                 the one shared normal BAM, injects params.mutect1_runs
assets/NO_FILE   empty sentinel for optional file inputs (self-provisioned
                 at pipeline start, same as your NO_PON_VCF/NO_ALLELES_VCF)
```

**Obsolete, left in place:** `bin/prepare_samplesheet.py` and
`tests/samplesheet.example.csv` are leftover from the earlier CSV-sample-sheet
draft superseded by `preprocess.py` above. They aren't referenced by
`main.nf` and can be ignored (or deleted manually) — they couldn't be
removed automatically from this folder.

## Running it

Inside Cirro: register `preprocess.py` as this pipeline's preprocessing
script the same way it's registered for the Mutect2 process, point
`--refFasta`/`--refFastaIdx`/`--refFastaDict` (and any of the optional
dbsnp/cosmic/ContEst params) at your reference bundle, and run.

Outside Cirro (e.g. local testing), `params.mutect1_runs` has no file-based
equivalent to point `--input` at — you'd set it directly in a
`-params-file` JSON:

```json
{
  "mutect1_runs": [
    {
      "output_prefix": "SAMPLE001",
      "tumor_reads": "/data/SAMPLE001.tumor.bam",
      "tumor_reads_index": "/data/SAMPLE001.tumor.bam.bai",
      "normal_reads": "/data/SHARED.normal.bam",
      "normal_reads_index": "/data/SHARED.normal.bam.bai",
      "tumor_sample_name": "SAMPLE001",
      "normal_sample_name": "SHARED-PBMC"
    }
  ],
  "refFasta": "/data/ref/reference_genome.fasta",
  "refFastaIdx": "/data/ref/reference_genome.fasta.fai",
  "refFastaDict": "/data/ref/reference_genome.dict"
}
```

```bash
nextflow run main.nf -params-file params.json -profile local
```

## Before this runs for real

- Build/push the combined Docker image (see the earlier `docker_single/`
  scaffold) and point `--container__mutect1_pipeline` at it.
- Confirm `preprocess.py`'s `NORMAL_NAME_MARKERS = ["PBMC"]` heuristic
  actually matches how normals are labeled in whatever dataset you're
  running this on — it's copied unchanged from the Mutect2 script.
- This is hand-reviewed but execution-untested — no `nextflow` binary or
  registry network access was available while writing it. I did catch and
  fix two real bugs on manual re-read (a `val("${...}")` that needed to be
  `env(...)` instead, and unguarded `params.x = default` assignments that
  would have clobbered CLI overrides) — there could be others a live run
  would surface that static review can't.