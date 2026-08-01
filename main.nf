#!/usr/bin/env nextflow
/*
 * Combined MuTect1 pipeline -- merges NeoDisc's Mutect_v1_calling.sh (bash,
 * on-prem/HPC, fixed per-chromosome scatter) with the getzlab
 * MuTect1_Scatter_Gather.wdl (Cromwell, cloud-native, dynamic N-way scatter
 * + ContEst contamination estimation), packaged to match the conventions of
 * an existing production Cirro Mutect2 pipeline (see preprocess.py):
 *
 *   - Sample discovery happens OUTSIDE Nextflow, via a Cirro dataset
 *     preprocessing hook (preprocess.py) that injects params.mutect1_runs --
 *     not a Nextflow-side sample-sheet CSV. This is the same
 *     cirro.helpers.preprocess_dataset.PreprocessDataset pattern the
 *     Mutect2 pipeline already uses in production.
 *   - Shared-normal model: one normal BAM, shared across every tumor sample
 *     in the dataset (confirmed choice -- matches the Mutect2 pipeline's
 *     active behavior, not its single-sample variant). This removes the
 *     need for any tumor-only/paired branching: every run has a normal by
 *     construction, so it's a broadcast value like ref_fasta, not something
 *     requiring a per-row join.
 *   - Sample names are derived from BAM header SM tags at runtime (samtools
 *     + awk, same idiom as the Mutect2 wrapper's tumor/normal SM
 *     extraction) rather than threaded through as separate params.
 *   - Scatter strategy: getzlab's dynamic N-way BAM-derived split (not
 *     NeoDisc's fixed 24-way per-chromosome BAITBYCHR split, and not the
 *     Mutect2 pipeline's physical per-shard BAM subsetting via samtools --
 *     MuTect1 shards are given the full BAM restricted with -L, matching
 *     getzlab's original approach; subsetting can be added later if I/O
 *     becomes a bottleneck at scale).
 *   - Filtering: MuTect1's own PASS/REJECT call is the same underlying
 *     filter in both source pipelines (NeoDisc: "drop REJECT"; getzlab:
 *     "keep PASS"). Both outputs are produced under each pipeline's own
 *     filename convention, plus getzlab's TiN-risk call_stats subset.
 *   - NOT included: NeoDisc's multi-caller ensemble MinOverlap filter --
 *     that needs HaplotypeCaller/Mutect2/Varscan2 output too and isn't part
 *     of "the MuTect1 workflow" in either source pipeline.
 *
 * Structurally untested -- no `nextflow` binary or registry network access
 * was available while writing this. Hand-reviewed against Nextflow DSL2
 * conventions and against the working Mutect2 pipeline's own idioms, not
 * validated with `nextflow run`.
 */
nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// PARAMS
// ---------------------------------------------------------------------------
// params.mutect1_runs is injected by preprocess.py at the Cirro dataset
// level -- NOT read from a file here, and deliberately has NO default
// assigned below (it's required; missing it is a hard error in the
// workflow block). Shape (one map per tumor sample):
//   [output_prefix, tumor_reads, tumor_reads_index,
//    normal_reads, normal_reads_index, tumor_sample_name, normal_sample_name]
//
// Every other default below is guarded with containsKey() -- matching the
// Mutect2 pipeline's own idiom -- because an unconditional `params.x = ...`
// would silently clobber a value the user supplied via -params-file or
// --x on the CLI every time the script loads, not just provide a fallback.
if (!params.containsKey('outdir'))                     params.outdir = 'results'

if (!params.containsKey('ref'))                        params.ref = null  // hg19 or hg38
if (!params.containsKey('ref_fasta'))                   params.ref_fasta = null
if (!params.containsKey('ref_fai'))                params.ref_fai = null
if (!params.containsKey('ref_dict'))               params.ref_dict = null

if (!params.containsKey('scatter_count'))              params.scatter_count = 200
if (!params.containsKey('target_list'))                params.target_list = null   // optional custom BED

// ContEst-only inputs -- if not all three are set, every tumor falls back to
// params.fracContam instead (see PLAN.md from the earlier WDL conversion re:
// these being "hidden" nested-call inputs in getzlab's original WDL).
if (!params.containsKey('contest_target_intervals'))   params.contest_target_intervals = null
if (!params.containsKey('snp6_bed'))                   params.snp6_bed = null
if (!params.containsKey('hapmap_vcf'))                 params.hapmap_vcf = null
if (!params.containsKey('fracContam'))                 params.fracContam = 0.01

// Optional MuTect1 knobs (all optional in both source pipelines except
// dbsnp/cosmic, which NeoDisc treats as mandatory -- kept optional here to
// match getzlab's more flexible calling config)
if (!params.containsKey('dbsnp'))                      params.dbsnp = null
if (!params.containsKey('dbsnp_idx'))                  params.dbsnp_idx = null
if (!params.containsKey('cosmic'))                     params.cosmic = null
if (!params.containsKey('cosmic_idx'))                 params.cosmic_idx = null
if (!params.containsKey('read_group_blacklist'))       params.read_group_blacklist = null
if (!params.containsKey('normal_panel'))               params.normal_panel = null
if (!params.containsKey('downsample'))                 params.downsample = 99999
if (!params.containsKey('max_mismatch_baseq_sum'))     params.max_mismatch_baseq_sum = 100
if (!params.containsKey('force_calling'))              params.force_calling = false
if (!params.containsKey('exclude_chimeric'))           params.exclude_chimeric = false

// Single combined image (see docker_single/ from the prior Docker
// discussion). `?:`-fallback style matches the Mutect2 pipeline's own
// `params.gatk_docker ?: 'broadinstitute/gatk:4.5.0.0'` convention.
//if (!params.containsKey('container__mutect1_pipeline')) params.container__mutect1_pipeline = null

def assetsDir = new File("${workflow.projectDir}/assets")
assetsDir.mkdirs()

def makeNoFile = { String name ->
    def f = file("${workflow.projectDir}/assets/NO_${name}")
    if (!f.exists()) {
        f.text = ''
    }
    return f
}

def NO_TARGET_LIST  = makeNoFile('TARGET_LIST')
def NO_DBSNP        = makeNoFile('DBSNP')
def NO_COSMIC       = makeNoFile('COSMIC')
def NO_RG_BLACKLIST = makeNoFile('RG_BLACKLIST')
def NO_NORMAL_PANEL = makeNoFile('NORMAL_PANEL')

// ---------------------------------------------------------------------------
// PROCESSES
// ---------------------------------------------------------------------------

// Port of getzlab's ContEst task, run once per tumor sample against the
// shared normal.
process CONTEST {
    tag "${pairName}"
    label 'process_high'
    container "ghcr.io/jchen1095/contest:1.0.0"
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(pairName), path(t_bam), path(t_bai)
    path normal_bam
    path normal_bai
    path targetIntervals
    path ref_fasta
    path ref_fai
    path ref_dict
    path snp6Bed
    path hapmapVcf

    output:
    tuple val(pairName), path('fraction_contamination.txt'), emit: frac

    shell:
    '''
    set -euxo pipefail

    if [ ! -f !{t_bam}.bai ]; then
        ln -s !{t_bai} !{t_bam}.bai
    fi
    if [ ! -f !{normal_bam}.bai ]; then
        ln -s !{normal_bai} !{normal_bam}.bai
    fi

    java_mem_mb=!{task.memory.toMega() - 1024}

    /usr/local/jre1.7.0_71/bin/java -Xmx${java_mem_mb}m -Djava.io.tmpdir=$PWD -jar /usr/local/bin/GenomeAnalysisTK.jar \
        -T ContEst \
        -I:eval !{t_bam} \
        -I:genotype !{normal_bam} \
        -L !{targetIntervals} \
        -L !{snp6Bed} \
        -isr INTERSECTION \
        -R !{ref_fasta} \
        -l INFO \
        -pf !{hapmapVcf} \
        -o contamination.af.txt \
        --trim_fraction 0.03 \
        --beta_threshold 0.05 \
        -br contamination.base_report.txt \
        -mbc 100 \
        --min_genotype_depth 30 \
        --min_genotype_ratio 0.8

    extract_contamination.py \
        contamination.af.txt \
        fraction_contamination.txt \
        contamination_validation.array_free.txt \
        !{pairName}

    test -s fraction_contamination.txt
    '''
}

// Port of getzlab's split_intervals task -- dynamic N-way scatter, run once
// per tumor sample.
process SPLIT_INTERVALS {
    tag "${pairName}"
    label 'process_low'
    container "ghcr.io/jchen1095/split_intervals:v48"
    errorStrategy 'retry'
    maxRetries 4

    input:
    tuple val(pairName), path(t_bam), path(t_bai)
    val   ref
    path  target_list

    output:
    tuple val(pairName), path('picard/*.interval_list'), emit: interval_files

    shell:
    '''
    set -euxo pipefail

    if [[ "!{ref}" == "hg38" ]]; then
        selected_chrs=$(printf "chr%s," {1..22} X Y | sed 's/,$//')
    elif [[ "!{ref}" == "hg19" ]]; then
        selected_chrs=$(printf "%s," {1..22} X Y | sed 's/,$//')
    else
        echo "unrecognized reference: !{ref}" >&2
        exit 1
    fi

    target_arg=""
    if [ "!{target_list}" != "NO_TARGET_LIST" ]; then
        target_arg="-target_list !{target_list}"
    fi

    split_intervals.py \
        -bam !{t_bam} \
        -bai !{t_bai} \
        -interval_type picard \
        -N !{params.scatter_count} \
        $target_arg \
        -chrs $selected_chrs

    for interval in picard/*.interval_list; do
        test -s "$interval"
    done
    '''
}

// Port of getzlab's MuTect1 task, invoked once per (tumor sample, interval
// shard), always against the shared normal. Derives tumor/normal sample
// names from BAM header SM tags (same samtools+awk idiom as the Mutect2
// wrapper) instead of taking them as params. Adds NeoDisc's per-shard
// log-grepping safety check on top of getzlab's version, plus `test -s`
// output assertions matching the Mutect2 pipeline's own defensive style --
// GATK3-era tools are known to sometimes not propagate certain internal
// failures via exit code alone.
process MUTECT1 {
    tag "${pairName}:${task.index}"
    label 'process_medium'
    container "ghcr.io/jchen1095/mutect_1_getzlab:v36"
    stageInMode 'symlink'
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(pairName), path(intervals), val(fracContam), path(t_bam), path(t_bai)
    path normal_bam
    path normal_bai
    path ref_fasta
    path ref_fai
    path ref_dict
    path dbsnp
    path cosmic
    path readGroupBlacklist
    path normalPanel

    output:
    // env() captures the named bash variable's value as left at the end of
    // the shell block -- this is the correct way to surface a
    // runtime-derived scalar (tumor_sample/normal_sample come from samtools,
    // not from a param), NOT `val("${tumor_sample}")`, which would try to
    // interpolate a Groovy variable that doesn't exist at pipeline-parse
    // time and silently produce garbage.
    tuple val(pairName), env(tumor_sample), env(normal_sample),
          path("${pairName}.MuTect1.call_stats.${task.index}.txt"),
          path("${pairName}.Mutect1.${task.index}.vcf"), emit: shard

    shell:
    '''
    set -euxo pipefail

    if [ ! -f !{t_bam}.bai ]; then
        ln -s !{t_bai} !{t_bam}.bai
    fi
    if [ ! -f !{normal_bam}.bai ]; then
        ln -s !{normal_bai} !{normal_bam}.bai
    fi

    tumor_sample=$(samtools view -H !{t_bam} \
        | awk -F'\t' '/^@RG/ {
            for (i=1; i<=NF; i++) {
              if ($i ~ /^SM:/) { sub(/^SM:/, "", $i); print $i }
            }
          }' | sort -u)

    normal_sample=$(samtools view -H !{normal_bam} \
        | awk -F'\t' '/^@RG/ {
            for (i=1; i<=NF; i++) {
              if ($i ~ /^SM:/) { sub(/^SM:/, "", $i); print $i }
            }
          }' | sort -u)

    [[ -n "$tumor_sample" ]] || { echo "ERROR: no SM tag in tumor BAM header" >&2; exit 1; }
    [[ -n "$normal_sample" ]] || { echo "ERROR: no SM tag in normal BAM header" >&2; exit 1; }
    [[ $(printf '%s\n' "$tumor_sample" | wc -l) -eq 1 ]] || { echo "ERROR: multiple tumor SM values: $tumor_sample" >&2; exit 1; }
    [[ $(printf '%s\n' "$normal_sample" | wc -l) -eq 1 ]] || { echo "ERROR: multiple normal SM values: $normal_sample" >&2; exit 1; }
    [[ "$tumor_sample" != "$normal_sample" ]] || { echo "ERROR: tumor and normal share the same SM value: $tumor_sample" >&2; exit 1; }

    java_mem_mb=!{task.memory.toMega() - 1024}

    extra_args=""
    if [ "!{dbsnp}" != "NO_DBSNP" ]; then
        extra_args="$extra_args --dbsnp !{dbsnp}"
    fi
    if [ "!{cosmic}" != "NO_COSMIC" ]; then
        extra_args="$extra_args --cosmic !{cosmic}"
    fi
    if [ "!{readGroupBlacklist}" != "NO_RG_BLACKLIST" ]; then
        extra_args="$extra_args --read_group_black_list !{readGroupBlacklist}"
    fi
    if [ "!{normalPanel}" != "NO_NORMAL_PANEL" ]; then
        extra_args="$extra_args --normal_panel !{normalPanel}"
    fi
    if [ "!{params.force_calling}" == "true" ]; then
        extra_args="$extra_args --force_output"
    fi
    if [ "!{params.exclude_chimeric}" == "true" ]; then
        extra_args="$extra_args --exclude_chimeric_reads"
    fi

    /usr/local/java/jdk1.7.0_80/bin/java -Xmx${java_mem_mb}m -jar /app/mutect.jar \
        --analysis_type MuTect \
        --tumor_sample_name "$tumor_sample" \
        -I:tumor !{t_bam} \
        -I:normal !{normal_bam} \
        --normal_sample_name "$normal_sample" \
        --reference_sequence !{ref_fasta} \
        --fraction_contamination !{fracContam} \
        --out !{pairName}.MuTect1.call_stats.!{task.index}.txt \
        --downsample_to_coverage !{params.downsample} \
        --max_read_mismatch_quality_score_sum !{params.max_mismatch_baseq_sum} \
        --vcf !{pairName}.Mutect1.!{task.index}.vcf \
        -L !{intervals} \
        $extra_args \
        > mutect1_shard.log 2>&1

    cat mutect1_shard.log

    if grep -Piq '(error)|(killed)|(java\\.lang\\.[a-zA-Z]*(exception|error):)' mutect1_shard.log; then
        echo "ERROR: MuTect1 shard !{task.index} for !{pairName} failed -- see log above" >&2
        exit 1
    fi

    test -s !{pairName}.MuTect1.call_stats.!{task.index}.txt
    test -s !{pairName}.Mutect1.!{task.index}.vcf
    '''
}

// Port of getzlab's Gather_MuTect1 task. Merges every shard for a tumor
// sample, then applies BOTH source pipelines' filtering:
//   - NeoDisc: drop REJECT calls -> <pairName>.mutectv1.final.vcf
//   - getzlab: PASS-only, coordinate-sorted VCF (same filter, kept under its
//     own name too) + a TiN-risk call_stats subset (getzlab-only, prep for a
//     downstream DeTiN step not included here)
process GATHER_AND_FILTER {
    tag "${pairName}"
    label 'process_medium'
    container "ghcr.io/jchen1095/mutect_1_getzlab:v36"
    publishDir "${params.outdir}/${pairName}", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(pairName), val(tumorSampleName), val(normalSampleName), path(call_stats_files), path(vcf_files)

    output:
    path "${pairName}.MuTect1.call_stats.txt",       emit: call_stats
    path "${pairName}.PASS+TiNrisk.call_stats.txt",  emit: tin_call_stats
    path "${pairName}.MuTect1.vcf",                  emit: vcf
    path "${pairName}.MuTect1.PASS.vcf",             emit: pass_vcf_getzlab
    path "${pairName}.mutectv1.final.vcf",           emit: final_vcf_neodisc

    shell:
    '''
    set -euxo pipefail

    printf '%s\n' !{call_stats_files} > call_stats.list
    printf '%s\n' !{vcf_files} > vcf.list

    merge_callstats.py --FILE_ARRAY_PATH call_stats.list !{pairName}.MuTect1.call_stats.txt
    merge_callstats.py --FILE_ARRAY_PATH vcf.list !{pairName}.MuTect1.orig.vcf

    python3 - <<PYCODE
vcf_file_in = "!{pairName}.MuTect1.orig.vcf"
vcf_file_out = "!{pairName}.MuTect1.vcf"
with open(vcf_file_in) as f_in, open(vcf_file_out, 'w') as f_out:
    for line in f_in:
        if line.startswith('#CHROM'):
            f_out.write('##normal_sample=!{normalSampleName}\n')
            f_out.write('##tumor_sample=!{tumorSampleName}\n')
        f_out.write(line)
PYCODE

    # getzlab-style: coordinate-sorted PASS-only VCF
    cat \
        <(sed -n '/^#/p' !{pairName}.MuTect1.vcf) \
        <(awk -F '\t' '$7 == "PASS" && $1 ~ /[0-9]+/' !{pairName}.MuTect1.vcf | sort -k1,1n -k2,2n) \
        <(awk -F '\t' '$7 == "PASS" && $1 !~ /[0-9]+/' !{pairName}.MuTect1.vcf | sort -k1,1 -k2,2n) \
        > !{pairName}.MuTect1.PASS.vcf

    # NeoDisc-style: same underlying filter (MuTect1 only ever tags FILTER
    # as PASS or REJECT), phrased as "drop REJECT", unsorted -- kept under
    # NeoDisc's own filename convention for parity with that pipeline.
    awk '!(/REJECT/)' !{pairName}.MuTect1.vcf > !{pairName}.mutectv1.final.vcf

    # getzlab-style: TiN-risk call_stats subset, prep for a downstream
    # DeTiN step (not included in this pipeline)
    awk -F '\t' 'NR <= 2 || $51 == "KEEP" || $50 == "alt_allele_in_normal" || $50 == "normal_lod" || $50 == "normal_lod,alt_allele_in_normal" || $50 == "alt_allele_in_normal,normal_lod"' !{pairName}.MuTect1.call_stats.txt > !{pairName}.PASS+TiNrisk.call_stats.txt

    test -s !{pairName}.MuTect1.call_stats.txt
    test -s !{pairName}.MuTect1.vcf
    test -s !{pairName}.MuTect1.PASS.vcf
    test -s !{pairName}.mutectv1.final.vcf
    '''
}

// ---------------------------------------------------------------------------
// WORKFLOW
// ---------------------------------------------------------------------------
workflow {

    if (!params.mutect1_runs) {
        error "params.mutect1_runs is empty -- did the Cirro preprocess.py hook run? (see preprocess.py)"
    }
    if (!params.ref_fasta) {
        error "Missing required param: ref_fasta"
    }

    // Same defensive check as the Mutect2 pipeline: every run must declare
    // a normal, and (shared-normal model) it must be the SAME normal for
    // every run in this dataset.
    params.mutect1_runs.each { run ->
        if (!run.normal_reads || !run.normal_reads_index) {
            error "Run ${run.output_prefix} is missing normal_reads/normal_reads_index -- this pipeline requires a matched normal."
        }
    }
    shared_normal_paths       = params.mutect1_runs.collect { it.normal_reads }.unique()
    shared_normal_index_paths = params.mutect1_runs.collect { it.normal_reads_index }.unique()
    if (shared_normal_paths.size() != 1) {
        error "Expected exactly one shared normal BAM across params.mutect1_runs, found: ${shared_normal_paths}"
    }
    if (shared_normal_index_paths.size() != 1) {
        error "Expected exactly one shared normal BAM index across params.mutect1_runs, found: ${shared_normal_index_paths}"
    }

    normal_bam = file(shared_normal_paths[0], checkIfExists: true)
    normal_bai = file(shared_normal_index_paths[0], checkIfExists: true)

    ref_fasta     = file(params.ref_fasta, checkIfExists: true)
    ref_fai  = file(params.ref_fai, checkIfExists: true)
    ref_dict = file(params.ref_dict, checkIfExists: true)
    target_list = params.target_list
        ? file(params.target_list, checkIfExists: true)
        : NO_TARGET_LIST

    dbsnp = params.dbsnp
        ? file(params.dbsnp, checkIfExists: true)
        : NO_DBSNP

    cosmic = params.cosmic
        ? file(params.cosmic, checkIfExists: true)
        : NO_COSMIC

    rgBlacklist = params.read_group_blacklist
        ? file(params.read_group_blacklist, checkIfExists: true)
        : NO_RG_BLACKLIST

    normalPanel = params.normal_panel
        ? file(params.normal_panel, checkIfExists: true)
        : NO_NORMAL_PANEL

    runs_ch = Channel.fromList(params.mutect1_runs).map { run ->
        tuple(
            run.output_prefix,
            file(run.tumor_reads, checkIfExists: true),
            file(run.tumor_reads_index, checkIfExists: true)
        )
    }

    // ---- contamination estimate: ContEst per tumor against the shared normal ----
    if (params.contest_target_intervals && params.snp6_bed && params.hapmap_vcf) {
        CONTEST(
            runs_ch, normal_bam, normal_bai,
            file(params.contest_target_intervals, checkIfExists: true),
            ref_fasta, ref_fai, ref_dict,
            file(params.snp6_bed, checkIfExists: true),
            file(params.hapmap_vcf, checkIfExists: true)
        )
        frac_ch = CONTEST.out.frac.map { pairName, fracFile -> tuple(pairName, fracFile.text.trim() as Float) }
    } else {
        frac_ch = runs_ch.map { pairName, t_bam, t_bai -> tuple(pairName, params.fracContam) }
    }

    runs_with_frac = runs_ch.join(frac_ch)
    // tuple(pairName, t_bam, t_bai, fracContam)

    // ---- scatter: split intervals per tumor, then fan out MuTect1 per shard ----
    SPLIT_INTERVALS(runs_ch, params.ref, target_list)

    shards_per_pair = SPLIT_INTERVALS.out.interval_files.transpose()
    // tuple(pairName, single_interval_file) -- one row per shard

    mutect1_inputs = shards_per_pair
        .combine(runs_with_frac, by: 0)
        .map { pairName, intervalFile, t_bam, t_bai, fracContam ->
            tuple(pairName, intervalFile, fracContam, t_bam, t_bai)
        }

    MUTECT1(
        mutect1_inputs, normal_bam, normal_bai,
        ref_fasta, ref_fai, ref_dict,
        dbsnp, cosmic,
        rgBlacklist, normalPanel
    )

    // ---- gather: group shards back by tumor sample, merge + filter ----------
    gathered = MUTECT1.out.shard.groupTuple(by: [0, 1, 2])
    // tuple(pairName, tumorSampleName, normalSampleName, [call_stats...], [vcf...])

    GATHER_AND_FILTER(gathered)
}