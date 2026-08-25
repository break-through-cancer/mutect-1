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
 * --- Observability fix (this revision) -----------------------------------
 * Two silent data-loss gaps were found and fixed:
 *   1. CONTEST had no publishDir, and only ever emitted the single derived
 *      fraction_contamination.txt -- the raw ContEst table
 *      (contamination.af.txt), the base report, and the array-free
 *      validation file were computed every run and then discarded. This is
 *      exactly the evidence you'd need to catch a bad contamination
 *      estimate (e.g. extract_contamination.py grabbing the wrong row/
 *      population panel out of ContEst's multi-row output). Fixed: CONTEST
 *      now has a publishDir and emits all four files.
 *   2. MUTECT1's per-shard mutect1_shard.log was cat'd to stdout and then
 *      discarded -- never declared as an output, never published. Fixed:
 *      each shard now emits its own log file, and GATHER_AND_FILTER
 *      concatenates every shard's log into one published per-sample log --
 *      the same pattern NeoDisc's own Mutect_v1_calling.sh already used
 *      (`cat ${c}_mutect.log >> ${tumor}_mutect.log`) before this got ported
 *      to Nextflow.
 *   Also added a publishDir to SPLIT_INTERVALS (cheap, and useful for
 *   auditing exactly how a sample's genome was scattered -- e.g. confirming
 *   shard count and boundaries match what you expect).
 *
 * Why this class of bug is Nextflow/Cirro-specific and didn't show up in
 * the original WDL/Terra pipeline: Cromwell + Terra expose a task's full
 * execution directory by default -- every file that ends up in a call's
 * working directory (stdout, stderr, anything written, not just files
 * named in `output {}`) is browsable via Terra's Job Manager / can be
 * pulled straight from the call's GCS path, with no extra configuration.
 * Nextflow's model is the opposite: only files explicitly listed in a
 * process's `output:` block are tracked at all, and even those are only
 * copied somewhere permanent/user-visible if that process has a
 * `publishDir`. Everything else lives in an ephemeral, hash-named work/
 * directory that isn't part of the normal Cirro results view and can be
 * garbage-collected. Porting the WDL to Nextflow didn't carry over Terra's
 * "everything is visible for free" behavior -- each process needs to opt in
 * explicitly, and CONTEST and MUTECT1 never did.
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

if (!params.containsKey('ref'))                        params.ref = null
if (!params.containsKey('ref_fasta'))                  params.ref_fasta = null
if (!params.containsKey('ref_fai'))                    params.ref_fai = null
if (!params.containsKey('ref_dict'))                   params.ref_dict = null

if (!params.containsKey('scatter_count'))              params.scatter_count = 200
if (!params.containsKey('target_list'))                params.target_list = null

// ContEst-only inputs
if (!params.containsKey('contest_target_intervals'))   params.contest_target_intervals = null
if (!params.containsKey('snp6_bed'))                   params.snp6_bed = null
if (!params.containsKey('hapmap_vcf'))                 params.hapmap_vcf = null
if (!params.containsKey('fracContam'))                 params.fracContam = 0.01

// Optional MuTect1 knobs
if (!params.containsKey('dbsnp'))                      params.dbsnp = null
if (!params.containsKey('dbsnp_idx'))                  params.dbsnp_idx = null
if (!params.containsKey('cosmic'))                     params.cosmic = null
if (!params.containsKey('cosmic_idx'))                 params.cosmic_idx = null
if (!params.containsKey('read_group_blacklist'))       params.read_group_blacklist = null
if (!params.containsKey('normal_panel'))               params.normal_panel = null
if (!params.containsKey('normal_panel_idx'))           params.normal_panel_idx = null
if (!params.containsKey('downsample'))                 params.downsample = 99999
if (!params.containsKey('max_mismatch_baseq_sum'))     params.max_mismatch_baseq_sum = 100
if (!params.containsKey('force_calling'))              params.force_calling = false
if (!params.containsKey('exclude_chimeric'))           params.exclude_chimeric = false

def assetsDir = new File("${workflow.projectDir}/assets")
assetsDir.mkdirs()

def makeNoFile = { String name ->
    def f = file("${workflow.projectDir}/assets/NO_${name}")
    if (!f.exists()) {
        f.text = ''
    }
    return f
}

def NO_TARGET_LIST      = makeNoFile('TARGET_LIST')
def NO_SNP6_BED         = makeNoFile('SNP6_BED')
def NO_DBSNP            = makeNoFile('DBSNP')
def NO_DBSNP_IDX        = makeNoFile('DBSNP_IDX')
def NO_COSMIC           = makeNoFile('COSMIC')
def NO_COSMIC_IDX       = makeNoFile('COSMIC_IDX')
def NO_RG_BLACKLIST     = makeNoFile('RG_BLACKLIST')
def NO_NORMAL_PANEL     = makeNoFile('NORMAL_PANEL')
def NO_NORMAL_PANEL_IDX = makeNoFile('NORMAL_PANEL_IDX')

// ---------------------------------------------------------------------------
// PROCESSES
// ---------------------------------------------------------------------------

// Index the ContEst population-frequency VCF once per workflow run.
process INDEX_CONTEST_VCF {
    tag "ContEst germline resource"
    label 'process_low'
    container "ghcr.io/jchen1095/contest:1.0.0"
    errorStrategy 'retry'
    maxRetries 2

    input:
    path hapmapVcf

    output:
    tuple path(hapmapVcf), path("${hapmapVcf}.tbi"), emit: indexed_vcf

    shell:
    '''
    set -euxo pipefail

    command -v tabix >/dev/null 2>&1 || {
        echo "ERROR: tabix is not installed in the ContEst container" >&2
        exit 1
    }

    tabix -f -p vcf !{hapmapVcf}
    test -s !{hapmapVcf}.tbi
    '''
}


// Port of getzlab's ContEst task, run once per tumor sample against the
// shared normal.
//
// FIX: added publishDir + emit the raw ContEst outputs (contamination.af.txt,
// contamination.base_report.txt, contamination_validation.array_free.txt) in
// addition to the single derived fraction_contamination.txt. Previously none
// of this was published anywhere -- the derived fraction was the only thing
// that (barely) survived past the ephemeral work directory, and the raw
// table needed to audit *how* that fraction was derived was discarded every
// run.
process CONTEST {
    tag "${pairName}"
    label 'process_high'
    container "ghcr.io/jchen1095/contest:1.0.0"
    publishDir "${params.outdir}/${pairName}/contest", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(pairName), path(t_bam), path(t_bai)
    path normal_bam
    path normal_bai
    path targetIntervals, name: 'contest_targets.interval_list'
    path ref_fasta
    path ref_fai
    path ref_dict
    path snp6Bed, name: 'snp6_targets.interval_list'
    tuple path(hapmapVcf), path(hapmapVcfIndex)

    output:
    tuple val(pairName), path('fraction_contamination.txt'), emit: frac
    path 'contamination.af.txt',                             emit: raw_af
    path 'contamination.base_report.txt',                    emit: base_report
    path 'contamination_validation.array_free.txt',           emit: validation

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

    # snp6Bed is optional -- it exists to restrict ContEst to positions also
    # present on an Affymetrix SNP6 array, which only matters if you're
    # cross-checking against real SNP6 array genotyping data (this pipeline
    # doesn't do that reconciliation step). Intersecting with it unconditionally
    # is what starved a real run down to 8345bp / 194 informative sites and
    # produced a statistically meaningless contamination estimate -- so this
    # is only applied if a real snp6_bed was actually supplied.
    snp6_args=""
    if [ "!{snp6Bed}" != "NO_SNP6_BED" ]; then
        snp6_args="-L !{snp6Bed} -isr INTERSECTION"
    fi

    /usr/local/jre1.7.0_71/bin/java \
        -Xmx${java_mem_mb}m \
        -Djava.io.tmpdir=$PWD \
        -jar /usr/local/bin/GenomeAnalysisTK.jar \
        -T ContEst \
        -I:eval !{t_bam} \
        -I:genotype !{normal_bam} \
        -L !{targetIntervals} \
        $snp6_args \
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

    python /usr/local/bin/extract_contamination.py \
        contamination.af.txt \
        fraction_contamination.txt \
        contamination_validation.array_free.txt \
        !{pairName}

    test -s fraction_contamination.txt
    test -s contamination.af.txt
    test -s contamination.base_report.txt
    '''
}


// ---------------------------------------------------------------------------
// SPLIT INTERVALS
// ---------------------------------------------------------------------------
//
// FIX: added publishDir. Cheap, and useful for auditing exactly how a
// sample's genome was scattered -- e.g. confirming the shard count and total
// covered bases match what you expect (this is exactly the kind of thing
// that would have let you catch a WES-vs-WGS scope mismatch immediately
// instead of inferring it after the fact from shard counts in a log).
process SPLIT_INTERVALS {
    tag "${pairName}"
    label 'process_low'
    container "ghcr.io/jchen1095/split_intervals:v48"
    publishDir "${params.outdir}/${pairName}/intervals", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 4

    input:
    tuple val(pairName), path(t_bam), path(t_bai)
    val ref
    path target_list

    output:
    tuple val(pairName), path('picard/*.interval_list'), emit: interval_files

    shell:
    '''
    set -euxo pipefail

    # Hardcoded for Ensembl/numeric style interval lists (Twist_HCEP_V1)

    selected_chrs=$(printf "chr%s," {1..22} X Y | sed 's/,$//')

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

// ---------------------------------------------------------------------------
// MUTECT1
// ---------------------------------------------------------------------------
//
// FIX: now emits its per-shard log (renamed to include pairName + task.index
// so it doesn't collide with other shards/samples once gathered) instead of
// just cat-ing it to stdout and discarding it. Previously this log -- which
// would show GATK3's own internal warnings, retries, and effective
// parameter echoes for that specific shard -- only ever existed in
// Nextflow's ephemeral work directory.
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
    path dbsnpIdx

    path cosmic
    path cosmicIdx

    path readGroupBlacklist
    path normalPanel
    path normalPanelIdx

    output:
    tuple val(pairName),
          env(tumor_sample),
          env(normal_sample),
          path("${pairName}.MuTect1.call_stats.${task.index}.txt"),
          path("${pairName}.Mutect1.${task.index}.vcf"),
          path("${pairName}.MuTect1.shard.${task.index}.log"),
          emit: shard

    shell:
    '''
    set -euxo pipefail

    # -----------------------------------------------------------------------
    # BAM indexes
    # -----------------------------------------------------------------------

    if [ ! -f !{t_bam}.bai ]; then
        ln -s !{t_bai} !{t_bam}.bai
    fi

    if [ ! -f !{normal_bam}.bai ]; then
        ln -s !{normal_bai} !{normal_bam}.bai
    fi


    # -----------------------------------------------------------------------
    # Determine sample names from BAM headers
    # -----------------------------------------------------------------------

    tumor_sample=$(samtools view -H !{t_bam} \
        | awk -F'\t' '/^@RG/ {
            for (i=1; i<=NF; i++) {
                if ($i ~ /^SM:/) {
                    sub(/^SM:/, "", $i)
                    print $i
                }
            }
        }' \
        | sort -u)

    normal_sample=$(samtools view -H !{normal_bam} \
        | awk -F'\t' '/^@RG/ {
            for (i=1; i<=NF; i++) {
                if ($i ~ /^SM:/) {
                    sub(/^SM:/, "", $i)
                    print $i
                }
            }
        }' \
        | sort -u)

    [[ -n "$tumor_sample" ]] || {
        echo "ERROR: no SM tag in tumor BAM header" >&2
        exit 1
    }

    [[ -n "$normal_sample" ]] || {
        echo "ERROR: no SM tag in normal BAM header" >&2
        exit 1
    }

    [[ $(printf '%s\n' "$tumor_sample" | wc -l) -eq 1 ]] || {
        echo "ERROR: multiple tumor SM values: $tumor_sample" >&2
        exit 1
    }

    [[ $(printf '%s\n' "$normal_sample" | wc -l) -eq 1 ]] || {
        echo "ERROR: multiple normal SM values: $normal_sample" >&2
        exit 1
    }

    [[ "$tumor_sample" != "$normal_sample" ]] || {
        echo "ERROR: tumor and normal share the same SM value: $tumor_sample" >&2
        exit 1
    }


    # -----------------------------------------------------------------------
    # Java memory
    # -----------------------------------------------------------------------

    java_mem_mb=!{task.memory.toMega() - 1024}


    # -----------------------------------------------------------------------
    # Optional MuTect1 inputs
    # -----------------------------------------------------------------------

    extra_args=""


    # dbSNP
    #
    # Old MuTect/GATK expects the Tribble index beside the VCF using:
    #
    #     dbsnp.vcf
    #     dbsnp.vcf.idx
    #
    # Nextflow stages dbsnpIdx separately, so explicitly create the expected
    # adjacent filename.
    if [ "!{dbsnp}" != "NO_DBSNP" ]; then

        if [ "!{dbsnpIdx}" == "NO_DBSNP_IDX" ]; then
            echo "ERROR: dbSNP VCF was supplied without dbsnp_idx" >&2
            exit 1
        fi

        if [ ! -e "!{dbsnp}.idx" ]; then
            ln -s "!{dbsnpIdx}" "!{dbsnp}.idx"
        fi

        test -e "!{dbsnp}.idx"

        extra_args="$extra_args --dbsnp !{dbsnp}"
    fi


    # COSMIC
    #
    # Same requirement as dbSNP: MuTect1 discovers the Tribble index from
    # <VCF>.idx, so make sure the staged index has the expected adjacent name.
    if [ "!{cosmic}" != "NO_COSMIC" ]; then

        if [ "!{cosmicIdx}" == "NO_COSMIC_IDX" ]; then
            echo "ERROR: COSMIC VCF was supplied without cosmic_idx" >&2
            exit 1
        fi

        if [ ! -e "!{cosmic}.idx" ]; then
            ln -s "!{cosmicIdx}" "!{cosmic}.idx"
        fi

        test -e "!{cosmic}.idx"

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


    # -----------------------------------------------------------------------
    # Run MuTect1
    # -----------------------------------------------------------------------

    shard_log="!{pairName}.MuTect1.shard.!{task.index}.log"

    /usr/local/java/jdk1.7.0_80/bin/java \
        -Xmx${java_mem_mb}m \
        -jar /app/mutect.jar \
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
        > "$shard_log" 2>&1

    cat "$shard_log"


    # -----------------------------------------------------------------------
    # Defensive error checks
    # -----------------------------------------------------------------------

    if grep -Piq '(error)|(killed)|(java\\.lang\\.[a-zA-Z]*(exception|error):)' \
        "$shard_log"
    then
        echo "ERROR: MuTect1 shard !{task.index} for !{pairName} failed -- see log above" >&2
        exit 1
    fi

    test -s !{pairName}.MuTect1.call_stats.!{task.index}.txt
    test -s !{pairName}.Mutect1.!{task.index}.vcf
    test -s "$shard_log"
    '''
}


// ---------------------------------------------------------------------------
// GATHER + FILTER
// ---------------------------------------------------------------------------
//
// FIX: now also receives every shard's log file and concatenates them into
// one published per-sample log (`${pairName}.MuTect1.log`) -- the same
// pattern NeoDisc's own Mutect_v1_calling.sh used for its per-chromosome
// logs before this got ported to Nextflow, just restored.
process GATHER_AND_FILTER {
    tag "${pairName}"
    label 'process_medium'
    container "ghcr.io/jchen1095/mutect_1_getzlab:v36"

    publishDir "${params.outdir}/${pairName}", mode: 'copy'

    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(pairName),
          val(tumorSampleName),
          val(normalSampleName),
          path(call_stats_files),
          path(vcf_files),
          path(log_files)

    output:
    path "${pairName}.MuTect1.call_stats.txt",         emit: call_stats
    path "${pairName}.PASS+TiNrisk.call_stats.txt",    emit: tin_call_stats
    path "${pairName}.MuTect1.vcf",                    emit: vcf
    path "${pairName}.MuTect1.PASS.vcf",               emit: pass_vcf_getzlab
    path "${pairName}.mutectv1.final.vcf",             emit: final_vcf_neodisc
    path "${pairName}.MuTect1.log",                    emit: shard_logs

    shell:
    '''
    set -euxo pipefail

    printf '%s\n' !{call_stats_files} > call_stats.list
    printf '%s\n' !{vcf_files} > vcf.list
    printf '%s\n' !{log_files} > log.list

    merge_callstats.py \
        --FILE_ARRAY_PATH call_stats.list \
        !{pairName}.MuTect1.call_stats.txt

    merge_callstats.py \
        --FILE_ARRAY_PATH vcf.list \
        !{pairName}.MuTect1.orig.vcf


    # -----------------------------------------------------------------------
    # Concatenate every shard's MuTect1 log into one per-sample log
    # (NeoDisc's own pattern: `cat ${c}_mutect.log >> ${tumor}_mutect.log`)
    # -----------------------------------------------------------------------

    : > !{pairName}.MuTect1.log
    while IFS= read -r logfile; do
        echo "===== $logfile =====" >> !{pairName}.MuTect1.log
        cat "$logfile" >> !{pairName}.MuTect1.log
    done < log.list


    # -----------------------------------------------------------------------
    # Add tumor/normal metadata
    # -----------------------------------------------------------------------

    python3 - <<PYCODE
vcf_file_in = "!{pairName}.MuTect1.orig.vcf"
vcf_file_out = "!{pairName}.MuTect1.vcf"

with open(vcf_file_in) as f_in, open(vcf_file_out, 'w') as f_out:
    for line in f_in:
        if line.startswith('#CHROM'):
            f_out.write('##normal_sample=!{normalSampleName}\\n')
            f_out.write('##tumor_sample=!{tumorSampleName}\\n')
        f_out.write(line)
PYCODE


    # -----------------------------------------------------------------------
    # getzlab-style PASS VCF
    # -----------------------------------------------------------------------

    cat \
        <(sed -n '/^#/p' !{pairName}.MuTect1.vcf) \
        <(
            awk -F '\t' \
                '$7 == "PASS" && $1 ~ /[0-9]+/' \
                !{pairName}.MuTect1.vcf \
                | sort -k1,1n -k2,2n
        ) \
        <(
            awk -F '\t' \
                '$7 == "PASS" && $1 !~ /[0-9]+/' \
                !{pairName}.MuTect1.vcf \
                | sort -k1,1 -k2,2n
        ) \
        > !{pairName}.MuTect1.PASS.vcf


    # -----------------------------------------------------------------------
    # NeoDisc-style "drop REJECT"
    # -----------------------------------------------------------------------

    awk '!(/REJECT/)' \
        !{pairName}.MuTect1.vcf \
        > !{pairName}.mutectv1.final.vcf


    # -----------------------------------------------------------------------
    # TiN-risk call_stats subset
    # -----------------------------------------------------------------------

    awk -F '\t' \
        'NR <= 2 ||
         $51 == "KEEP" ||
         $50 == "alt_allele_in_normal" ||
         $50 == "normal_lod" ||
         $50 == "normal_lod,alt_allele_in_normal" ||
         $50 == "alt_allele_in_normal,normal_lod"' \
        !{pairName}.MuTect1.call_stats.txt \
        > !{pairName}.PASS+TiNrisk.call_stats.txt


    # -----------------------------------------------------------------------
    # Output assertions
    # -----------------------------------------------------------------------

    test -s !{pairName}.MuTect1.call_stats.txt
    test -s !{pairName}.MuTect1.vcf
    test -s !{pairName}.MuTect1.PASS.vcf
    test -s !{pairName}.mutectv1.final.vcf
    test -s !{pairName}.MuTect1.log
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


    // -----------------------------------------------------------------------
    // Shared matched normal
    // -----------------------------------------------------------------------

    params.mutect1_runs.each { run ->

        if (!run.normal_reads || !run.normal_reads_index) {
            error "Run ${run.output_prefix} is missing normal_reads/normal_reads_index -- this pipeline requires a matched normal."
        }
    }

    shared_normal_paths =
        params.mutect1_runs
            .collect { it.normal_reads }
            .unique()

    shared_normal_index_paths =
        params.mutect1_runs
            .collect { it.normal_reads_index }
            .unique()

    if (shared_normal_paths.size() != 1) {
        error "Expected exactly one shared normal BAM across params.mutect1_runs, found: ${shared_normal_paths}"
    }

    if (shared_normal_index_paths.size() != 1) {
        error "Expected exactly one shared normal BAM index across params.mutect1_runs, found: ${shared_normal_index_paths}"
    }


    normal_bam =
        file(
            shared_normal_paths[0],
            checkIfExists: true
        )

    normal_bai =
        file(
            shared_normal_index_paths[0],
            checkIfExists: true
        )


    // -----------------------------------------------------------------------
    // Reference
    // -----------------------------------------------------------------------

    ref_fasta =
        file(
            params.ref_fasta,
            checkIfExists: true
        )

    ref_fai =
        file(
            params.ref_fai,
            checkIfExists: true
        )

    ref_dict =
        file(
            params.ref_dict,
            checkIfExists: true
        )


    // -----------------------------------------------------------------------
    // Optional target list
    // -----------------------------------------------------------------------

    target_list =
        params.target_list
            ? file(params.target_list, checkIfExists: true)
            : NO_TARGET_LIST


    // -----------------------------------------------------------------------
    // dbSNP + index
    // -----------------------------------------------------------------------

    dbsnp =
        params.dbsnp
            ? file(params.dbsnp, checkIfExists: true)
            : NO_DBSNP

    dbsnpIdx =
        params.dbsnp_idx
            ? file(params.dbsnp_idx, checkIfExists: true)
            : NO_DBSNP_IDX


    // -----------------------------------------------------------------------
    // COSMIC + index
    // -----------------------------------------------------------------------

    cosmic =
        params.cosmic
            ? file(params.cosmic, checkIfExists: true)
            : NO_COSMIC

    cosmicIdx =
        params.cosmic_idx
            ? file(params.cosmic_idx, checkIfExists: true)
            : NO_COSMIC_IDX


    // -----------------------------------------------------------------------
    // Other optional resources
    // -----------------------------------------------------------------------

    rgBlacklist =
        params.read_group_blacklist
            ? file(
                params.read_group_blacklist,
                checkIfExists: true
            )
            : NO_RG_BLACKLIST

    normalPanel =
        params.normal_panel
            ? file(
                params.normal_panel,
                checkIfExists: true
            )
            : NO_NORMAL_PANEL

    normalPanelIdx =
        params.normal_panel_idx
            ? file(
                params.normal_panel_idx,
                checkIfExists: true
            )
            : NO_NORMAL_PANEL_IDX


    // -----------------------------------------------------------------------
    // Tumor runs
    // -----------------------------------------------------------------------

    runs_ch =
        Channel
            .fromList(params.mutect1_runs)
            .map { run ->

                tuple(
                    run.output_prefix,
                    file(
                        run.tumor_reads,
                        checkIfExists: true
                    ),
                    file(
                        run.tumor_reads_index,
                        checkIfExists: true
                    )
                )
            }


    // -----------------------------------------------------------------------
    // ContEst
    // -----------------------------------------------------------------------

    // snp6_bed is genuinely optional (see CONTEST's shell block) -- only
    // contest_target_intervals + hapmap_vcf actually gate whether ContEst
    // runs at all. Unconditionally intersecting with a SNP6 bed shrank a
    // real run down to 194 informative sites -- see the FIX note above
    // CONTEST -- so don't require it here either.
    if (
        params.contest_target_intervals &&
        params.hapmap_vcf
    ) {

        contest_hapmap_vcf =
            file(
                params.hapmap_vcf,
                checkIfExists: true
            )

        INDEX_CONTEST_VCF(
            contest_hapmap_vcf
        )

        snp6Bed =
            params.snp6_bed
                ? file(params.snp6_bed, checkIfExists: true)
                : NO_SNP6_BED

        CONTEST(
            runs_ch,
            normal_bam,
            normal_bai,
            file(
                params.contest_target_intervals,
                checkIfExists: true
            ),
            ref_fasta,
            ref_fai,
            ref_dict,
            snp6Bed,
            INDEX_CONTEST_VCF.out.indexed_vcf
        )

        frac_ch =
            CONTEST
                .out
                .frac
                .map {
                    pairName,
                    fracFile ->

                    tuple(
                        pairName,
                        fracFile.text.trim() as Float
                    )
                }

    } else {

        frac_ch =
            runs_ch
                .map {
                    pairName,
                    t_bam,
                    t_bai ->

                    tuple(
                        pairName,
                        params.fracContam
                    )
                }
    }


    runs_with_frac =
        runs_ch.join(frac_ch)

    // tuple(
    //     pairName,
    //     t_bam,
    //     t_bai,
    //     fracContam
    // )


    // -----------------------------------------------------------------------
    // Scatter
    // -----------------------------------------------------------------------

    SPLIT_INTERVALS(
        runs_ch,
        params.ref,
        target_list
    )

    shards_per_pair =
        SPLIT_INTERVALS
            .out
            .interval_files
            .transpose()

    // tuple(
    //     pairName,
    //     single_interval_file
    // )


    // -----------------------------------------------------------------------
    // Join intervals with tumor + contamination data
    // -----------------------------------------------------------------------

    mutect1_inputs =
        shards_per_pair
            .combine(
                runs_with_frac,
                by: 0
            )
            .map {
                pairName,
                intervalFile,
                t_bam,
                t_bai,
                fracContam ->

                tuple(
                    pairName,
                    intervalFile,
                    fracContam,
                    t_bam,
                    t_bai
                )
            }


    // -----------------------------------------------------------------------
    // MuTect1
    // -----------------------------------------------------------------------

    MUTECT1(
        mutect1_inputs,

        normal_bam,
        normal_bai,

        ref_fasta,
        ref_fai,
        ref_dict,

        dbsnp,
        dbsnpIdx,

        cosmic,
        cosmicIdx,

        rgBlacklist,
        normalPanel,
        normalPanelIdx
    )


    // -----------------------------------------------------------------------
    // Gather
    // -----------------------------------------------------------------------

    gathered =
        MUTECT1
            .out
            .shard
            .groupTuple(
                by: [0, 1, 2]
            )

    // tuple(
    //     pairName,
    //     tumorSampleName,
    //     normalSampleName,
    //     [call_stats...],
    //     [vcf...],
    //     [log...]
    // )

    GATHER_AND_FILTER(
        gathered
    )
}

// #!/usr/bin/env nextflow
// /*
//  * Combined MuTect1 pipeline -- merges NeoDisc's Mutect_v1_calling.sh (bash,
//  * on-prem/HPC, fixed per-chromosome scatter) with the getzlab
//  * MuTect1_Scatter_Gather.wdl (Cromwell, cloud-native, dynamic N-way scatter
//  * + ContEst contamination estimation), packaged to match the conventions of
//  * an existing production Cirro Mutect2 pipeline (see preprocess.py):
//  *
//  *   - Sample discovery happens OUTSIDE Nextflow, via a Cirro dataset
//  *     preprocessing hook (preprocess.py) that injects params.mutect1_runs --
//  *     not a Nextflow-side sample-sheet CSV. This is the same
//  *     cirro.helpers.preprocess_dataset.PreprocessDataset pattern the
//  *     Mutect2 pipeline already uses in production.
//  *   - Shared-normal model: one normal BAM, shared across every tumor sample
//  *     in the dataset (confirmed choice -- matches the Mutect2 pipeline's
//  *     active behavior, not its single-sample variant). This removes the
//  *     need for any tumor-only/paired branching: every run has a normal by
//  *     construction, so it's a broadcast value like ref_fasta, not something
//  *     requiring a per-row join.
//  *   - Sample names are derived from BAM header SM tags at runtime (samtools
//  *     + awk, same idiom as the Mutect2 wrapper's tumor/normal SM
//  *     extraction) rather than threaded through as separate params.
//  *   - Scatter strategy: getzlab's dynamic N-way BAM-derived split (not
//  *     NeoDisc's fixed 24-way per-chromosome BAITBYCHR split, and not the
//  *     Mutect2 pipeline's physical per-shard BAM subsetting via samtools --
//  *     MuTect1 shards are given the full BAM restricted with -L, matching
//  *     getzlab's original approach; subsetting can be added later if I/O
//  *     becomes a bottleneck at scale).
//  *   - Filtering: MuTect1's own PASS/REJECT call is the same underlying
//  *     filter in both source pipelines (NeoDisc: "drop REJECT"; getzlab:
//  *     "keep PASS"). Both outputs are produced under each pipeline's own
//  *     filename convention, plus getzlab's TiN-risk call_stats subset.
//  *   - NOT included: NeoDisc's multi-caller ensemble MinOverlap filter --
//  *     that needs HaplotypeCaller/Mutect2/Varscan2 output too and isn't part
//  *     of "the MuTect1 workflow" in either source pipeline.
//  *
//  * Structurally untested -- no `nextflow` binary or registry network access
//  * was available while writing this. Hand-reviewed against Nextflow DSL2
//  * conventions and against the working Mutect2 pipeline's own idioms, not
//  * validated with `nextflow run`.
//  */
// nextflow.enable.dsl = 2

// // ---------------------------------------------------------------------------
// // PARAMS
// // ---------------------------------------------------------------------------
// // params.mutect1_runs is injected by preprocess.py at the Cirro dataset
// // level -- NOT read from a file here, and deliberately has NO default
// // assigned below (it's required; missing it is a hard error in the
// // workflow block). Shape (one map per tumor sample):
// //   [output_prefix, tumor_reads, tumor_reads_index,
// //    normal_reads, normal_reads_index, tumor_sample_name, normal_sample_name]
// //
// // Every other default below is guarded with containsKey() -- matching the
// // Mutect2 pipeline's own idiom -- because an unconditional `params.x = ...`
// // would silently clobber a value the user supplied via -params-file or
// // --x on the CLI every time the script loads, not just provide a fallback.
// if (!params.containsKey('outdir'))                     params.outdir = 'results'

// if (!params.containsKey('ref'))                        params.ref = null
// if (!params.containsKey('ref_fasta'))                  params.ref_fasta = null
// if (!params.containsKey('ref_fai'))                    params.ref_fai = null
// if (!params.containsKey('ref_dict'))                   params.ref_dict = null

// if (!params.containsKey('scatter_count'))              params.scatter_count = 200
// if (!params.containsKey('target_list'))                params.target_list = null

// // ContEst-only inputs
// if (!params.containsKey('contest_target_intervals'))   params.contest_target_intervals = null
// if (!params.containsKey('snp6_bed'))                   params.snp6_bed = null
// if (!params.containsKey('hapmap_vcf'))                 params.hapmap_vcf = null
// if (!params.containsKey('fracContam'))                 params.fracContam = 0.01

// // Optional MuTect1 knobs
// if (!params.containsKey('dbsnp'))                      params.dbsnp = null
// if (!params.containsKey('dbsnp_idx'))                  params.dbsnp_idx = null
// if (!params.containsKey('cosmic'))                     params.cosmic = null
// if (!params.containsKey('cosmic_idx'))                 params.cosmic_idx = null
// if (!params.containsKey('read_group_blacklist'))       params.read_group_blacklist = null
// if (!params.containsKey('normal_panel'))               params.normal_panel = null
// if (!params.containsKey('normal_panel_idx'))           params.normal_panel_idx = null
// if (!params.containsKey('downsample'))                 params.downsample = 99999
// if (!params.containsKey('max_mismatch_baseq_sum'))     params.max_mismatch_baseq_sum = 100
// if (!params.containsKey('force_calling'))              params.force_calling = false
// if (!params.containsKey('exclude_chimeric'))           params.exclude_chimeric = false

// def assetsDir = new File("${workflow.projectDir}/assets")
// assetsDir.mkdirs()

// def makeNoFile = { String name ->
//     def f = file("${workflow.projectDir}/assets/NO_${name}")
//     if (!f.exists()) {
//         f.text = ''
//     }
//     return f
// }

// def NO_TARGET_LIST      = makeNoFile('TARGET_LIST')
// def NO_DBSNP            = makeNoFile('DBSNP')
// def NO_DBSNP_IDX        = makeNoFile('DBSNP_IDX')
// def NO_COSMIC           = makeNoFile('COSMIC')
// def NO_COSMIC_IDX       = makeNoFile('COSMIC_IDX')
// def NO_RG_BLACKLIST     = makeNoFile('RG_BLACKLIST')
// def NO_NORMAL_PANEL     = makeNoFile('NORMAL_PANEL')
// def NO_NORMAL_PANEL_IDX = makeNoFile('NORMAL_PANEL_IDX')

// // ---------------------------------------------------------------------------
// // PROCESSES
// // ---------------------------------------------------------------------------

// // Index the ContEst population-frequency VCF once per workflow run.
// process INDEX_CONTEST_VCF {
//     tag "ContEst germline resource"
//     label 'process_low'
//     container "ghcr.io/jchen1095/contest:1.0.0"
//     errorStrategy 'retry'
//     maxRetries 2

//     input:
//     path hapmapVcf

//     output:
//     tuple path(hapmapVcf), path("${hapmapVcf}.tbi"), emit: indexed_vcf

//     shell:
//     '''
//     set -euxo pipefail

//     command -v tabix >/dev/null 2>&1 || {
//         echo "ERROR: tabix is not installed in the ContEst container" >&2
//         exit 1
//     }

//     tabix -f -p vcf !{hapmapVcf}
//     test -s !{hapmapVcf}.tbi
//     '''
// }


// // Port of getzlab's ContEst task, run once per tumor sample against the
// // shared normal.
// process CONTEST {
//     tag "${pairName}"
//     label 'process_high'
//     container "ghcr.io/jchen1095/contest:1.0.0"
//     errorStrategy 'retry'
//     maxRetries 3

//     input:
//     tuple val(pairName), path(t_bam), path(t_bai)
//     path normal_bam
//     path normal_bai
//     path targetIntervals, name: 'contest_targets.interval_list'
//     path ref_fasta
//     path ref_fai
//     path ref_dict
//     path snp6Bed, name: 'snp6_targets.interval_list'
//     tuple path(hapmapVcf), path(hapmapVcfIndex)

//     output:
//     tuple val(pairName), path('fraction_contamination.txt'), emit: frac

//     shell:
//     '''
//     set -euxo pipefail

//     if [ ! -f !{t_bam}.bai ]; then
//         ln -s !{t_bai} !{t_bam}.bai
//     fi

//     if [ ! -f !{normal_bam}.bai ]; then
//         ln -s !{normal_bai} !{normal_bam}.bai
//     fi

//     java_mem_mb=!{task.memory.toMega() - 1024}

//     /usr/local/jre1.7.0_71/bin/java \
//         -Xmx${java_mem_mb}m \
//         -Djava.io.tmpdir=$PWD \
//         -jar /usr/local/bin/GenomeAnalysisTK.jar \
//         -T ContEst \
//         -I:eval !{t_bam} \
//         -I:genotype !{normal_bam} \
//         -L !{targetIntervals} \
//         -L !{snp6Bed} \
//         -isr INTERSECTION \
//         -R !{ref_fasta} \
//         -l INFO \
//         -pf !{hapmapVcf} \
//         -o contamination.af.txt \
//         --trim_fraction 0.03 \
//         --beta_threshold 0.05 \
//         -br contamination.base_report.txt \
//         -mbc 100 \
//         --min_genotype_depth 30 \
//         --min_genotype_ratio 0.8

//     python /usr/local/bin/extract_contamination.py \
//         contamination.af.txt \
//         fraction_contamination.txt \
//         contamination_validation.array_free.txt \
//         !{pairName}

//     test -s fraction_contamination.txt
//     '''
// }


// // ---------------------------------------------------------------------------
// // SPLIT INTERVALS
// // ---------------------------------------------------------------------------

// process SPLIT_INTERVALS {
//     tag "${pairName}"
//     label 'process_low'
//     container "ghcr.io/jchen1095/split_intervals:v48"
//     errorStrategy 'retry'
//     maxRetries 4

//     input:
//     tuple val(pairName), path(t_bam), path(t_bai)
//     val ref
//     path target_list

//     output:
//     tuple val(pairName), path('picard/*.interval_list'), emit: interval_files

//     shell:
//     '''
//     set -euxo pipefail

//     if [[ "!{ref}" == "hg38" ]]; then
//         selected_chrs=$(printf "chr%s," {1..22} X Y | sed 's/,$//')
//     elif [[ "!{ref}" == "hg19" ]]; then
//         selected_chrs=$(printf "%s," {1..22} X Y | sed 's/,$//')
//     else
//         echo "unrecognized reference: !{ref}" >&2
//         exit 1
//     fi

//     target_arg=""

//     if [ "!{target_list}" != "NO_TARGET_LIST" ]; then
//         target_arg="-target_list !{target_list}"
//     fi

//     split_intervals.py \
//         -bam !{t_bam} \
//         -bai !{t_bai} \
//         -interval_type picard \
//         -N !{params.scatter_count} \
//         $target_arg \
//         -chrs $selected_chrs

//     for interval in picard/*.interval_list; do
//         test -s "$interval"
//     done
//     '''
// }


// // ---------------------------------------------------------------------------
// // MUTECT1
// // ---------------------------------------------------------------------------

// process MUTECT1 {
//     tag "${pairName}:${task.index}"
//     label 'process_medium'
//     container "ghcr.io/jchen1095/mutect_1_getzlab:v36"
//     stageInMode 'symlink'
//     errorStrategy 'retry'
//     maxRetries 3

//     input:
//     tuple val(pairName), path(intervals), val(fracContam), path(t_bam), path(t_bai)

//     path normal_bam
//     path normal_bai

//     path ref_fasta
//     path ref_fai
//     path ref_dict

//     path dbsnp
//     path dbsnpIdx

//     path cosmic
//     path cosmicIdx

//     path readGroupBlacklist
//     path normalPanel
//     path normalPanelIdx

//     output:
//     tuple val(pairName),
//           env(tumor_sample),
//           env(normal_sample),
//           path("${pairName}.MuTect1.call_stats.${task.index}.txt"),
//           path("${pairName}.Mutect1.${task.index}.vcf"),
//           emit: shard

//     shell:
//     '''
//     set -euxo pipefail

//     # -----------------------------------------------------------------------
//     # BAM indexes
//     # -----------------------------------------------------------------------

//     if [ ! -f !{t_bam}.bai ]; then
//         ln -s !{t_bai} !{t_bam}.bai
//     fi

//     if [ ! -f !{normal_bam}.bai ]; then
//         ln -s !{normal_bai} !{normal_bam}.bai
//     fi


//     # -----------------------------------------------------------------------
//     # Determine sample names from BAM headers
//     # -----------------------------------------------------------------------

//     tumor_sample=$(samtools view -H !{t_bam} \
//         | awk -F'\t' '/^@RG/ {
//             for (i=1; i<=NF; i++) {
//                 if ($i ~ /^SM:/) {
//                     sub(/^SM:/, "", $i)
//                     print $i
//                 }
//             }
//         }' \
//         | sort -u)

//     normal_sample=$(samtools view -H !{normal_bam} \
//         | awk -F'\t' '/^@RG/ {
//             for (i=1; i<=NF; i++) {
//                 if ($i ~ /^SM:/) {
//                     sub(/^SM:/, "", $i)
//                     print $i
//                 }
//             }
//         }' \
//         | sort -u)

//     [[ -n "$tumor_sample" ]] || {
//         echo "ERROR: no SM tag in tumor BAM header" >&2
//         exit 1
//     }

//     [[ -n "$normal_sample" ]] || {
//         echo "ERROR: no SM tag in normal BAM header" >&2
//         exit 1
//     }

//     [[ $(printf '%s\n' "$tumor_sample" | wc -l) -eq 1 ]] || {
//         echo "ERROR: multiple tumor SM values: $tumor_sample" >&2
//         exit 1
//     }

//     [[ $(printf '%s\n' "$normal_sample" | wc -l) -eq 1 ]] || {
//         echo "ERROR: multiple normal SM values: $normal_sample" >&2
//         exit 1
//     }

//     [[ "$tumor_sample" != "$normal_sample" ]] || {
//         echo "ERROR: tumor and normal share the same SM value: $tumor_sample" >&2
//         exit 1
//     }


//     # -----------------------------------------------------------------------
//     # Java memory
//     # -----------------------------------------------------------------------

//     java_mem_mb=!{task.memory.toMega() - 1024}


//     # -----------------------------------------------------------------------
//     # Optional MuTect1 inputs
//     # -----------------------------------------------------------------------

//     extra_args=""


//     # dbSNP
//     #
//     # Old MuTect/GATK expects the Tribble index beside the VCF using:
//     #
//     #     dbsnp.vcf
//     #     dbsnp.vcf.idx
//     #
//     # Nextflow stages dbsnpIdx separately, so explicitly create the expected
//     # adjacent filename.
//     if [ "!{dbsnp}" != "NO_DBSNP" ]; then

//         if [ "!{dbsnpIdx}" == "NO_DBSNP_IDX" ]; then
//             echo "ERROR: dbSNP VCF was supplied without dbsnp_idx" >&2
//             exit 1
//         fi

//         if [ ! -e "!{dbsnp}.idx" ]; then
//             ln -s "!{dbsnpIdx}" "!{dbsnp}.idx"
//         fi

//         test -e "!{dbsnp}.idx"

//         extra_args="$extra_args --dbsnp !{dbsnp}"
//     fi


//     # COSMIC
//     #
//     # Same requirement as dbSNP: MuTect1 discovers the Tribble index from
//     # <VCF>.idx, so make sure the staged index has the expected adjacent name.
//     if [ "!{cosmic}" != "NO_COSMIC" ]; then

//         if [ "!{cosmicIdx}" == "NO_COSMIC_IDX" ]; then
//             echo "ERROR: COSMIC VCF was supplied without cosmic_idx" >&2
//             exit 1
//         fi

//         if [ ! -e "!{cosmic}.idx" ]; then
//             ln -s "!{cosmicIdx}" "!{cosmic}.idx"
//         fi

//         test -e "!{cosmic}.idx"

//         extra_args="$extra_args --cosmic !{cosmic}"
//     fi


//     if [ "!{readGroupBlacklist}" != "NO_RG_BLACKLIST" ]; then
//         extra_args="$extra_args --read_group_black_list !{readGroupBlacklist}"
//     fi


//     if [ "!{normalPanel}" != "NO_NORMAL_PANEL" ]; then
//         extra_args="$extra_args --normal_panel !{normalPanel}"
//     fi


//     if [ "!{params.force_calling}" == "true" ]; then
//         extra_args="$extra_args --force_output"
//     fi


//     if [ "!{params.exclude_chimeric}" == "true" ]; then
//         extra_args="$extra_args --exclude_chimeric_reads"
//     fi


//     # -----------------------------------------------------------------------
//     # Run MuTect1
//     # -----------------------------------------------------------------------

//     /usr/local/java/jdk1.7.0_80/bin/java \
//         -Xmx${java_mem_mb}m \
//         -jar /app/mutect.jar \
//         --analysis_type MuTect \
//         --tumor_sample_name "$tumor_sample" \
//         -I:tumor !{t_bam} \
//         -I:normal !{normal_bam} \
//         --normal_sample_name "$normal_sample" \
//         --reference_sequence !{ref_fasta} \
//         --fraction_contamination !{fracContam} \
//         --out !{pairName}.MuTect1.call_stats.!{task.index}.txt \
//         --downsample_to_coverage !{params.downsample} \
//         --max_read_mismatch_quality_score_sum !{params.max_mismatch_baseq_sum} \
//         --vcf !{pairName}.Mutect1.!{task.index}.vcf \
//         -L !{intervals} \
//         $extra_args \
//         > mutect1_shard.log 2>&1

//     cat mutect1_shard.log


//     # -----------------------------------------------------------------------
//     # Defensive error checks
//     # -----------------------------------------------------------------------

//     if grep -Piq '(error)|(killed)|(java\\.lang\\.[a-zA-Z]*(exception|error):)' \
//         mutect1_shard.log
//     then
//         echo "ERROR: MuTect1 shard !{task.index} for !{pairName} failed -- see log above" >&2
//         exit 1
//     fi

//     test -s !{pairName}.MuTect1.call_stats.!{task.index}.txt
//     test -s !{pairName}.Mutect1.!{task.index}.vcf
//     '''
// }


// // ---------------------------------------------------------------------------
// // GATHER + FILTER
// // ---------------------------------------------------------------------------

// process GATHER_AND_FILTER {
//     tag "${pairName}"
//     label 'process_medium'
//     container "ghcr.io/jchen1095/mutect_1_getzlab:v36"

//     publishDir "${params.outdir}/${pairName}", mode: 'copy'

//     errorStrategy 'retry'
//     maxRetries 2

//     input:
//     tuple val(pairName),
//           val(tumorSampleName),
//           val(normalSampleName),
//           path(call_stats_files),
//           path(vcf_files)

//     output:
//     path "${pairName}.MuTect1.call_stats.txt",         emit: call_stats
//     path "${pairName}.PASS+TiNrisk.call_stats.txt",    emit: tin_call_stats
//     path "${pairName}.MuTect1.vcf",                    emit: vcf
//     path "${pairName}.MuTect1.PASS.vcf",               emit: pass_vcf_getzlab
//     path "${pairName}.mutectv1.final.vcf",             emit: final_vcf_neodisc

//     shell:
//     '''
//     set -euxo pipefail

//     printf '%s\n' !{call_stats_files} > call_stats.list
//     printf '%s\n' !{vcf_files} > vcf.list

//     merge_callstats.py \
//         --FILE_ARRAY_PATH call_stats.list \
//         !{pairName}.MuTect1.call_stats.txt

//     merge_callstats.py \
//         --FILE_ARRAY_PATH vcf.list \
//         !{pairName}.MuTect1.orig.vcf


//     # -----------------------------------------------------------------------
//     # Add tumor/normal metadata
//     # -----------------------------------------------------------------------

//     python3 - <<PYCODE
// vcf_file_in = "!{pairName}.MuTect1.orig.vcf"
// vcf_file_out = "!{pairName}.MuTect1.vcf"

// with open(vcf_file_in) as f_in, open(vcf_file_out, 'w') as f_out:
//     for line in f_in:
//         if line.startswith('#CHROM'):
//             f_out.write('##normal_sample=!{normalSampleName}\\n')
//             f_out.write('##tumor_sample=!{tumorSampleName}\\n')
//         f_out.write(line)
// PYCODE


//     # -----------------------------------------------------------------------
//     # getzlab-style PASS VCF
//     # -----------------------------------------------------------------------

//     cat \
//         <(sed -n '/^#/p' !{pairName}.MuTect1.vcf) \
//         <(
//             awk -F '\t' \
//                 '$7 == "PASS" && $1 ~ /[0-9]+/' \
//                 !{pairName}.MuTect1.vcf \
//                 | sort -k1,1n -k2,2n
//         ) \
//         <(
//             awk -F '\t' \
//                 '$7 == "PASS" && $1 !~ /[0-9]+/' \
//                 !{pairName}.MuTect1.vcf \
//                 | sort -k1,1 -k2,2n
//         ) \
//         > !{pairName}.MuTect1.PASS.vcf


//     # -----------------------------------------------------------------------
//     # NeoDisc-style "drop REJECT"
//     # -----------------------------------------------------------------------

//     awk '!(/REJECT/)' \
//         !{pairName}.MuTect1.vcf \
//         > !{pairName}.mutectv1.final.vcf


//     # -----------------------------------------------------------------------
//     # TiN-risk call_stats subset
//     # -----------------------------------------------------------------------

//     awk -F '\t' \
//         'NR <= 2 ||
//          $51 == "KEEP" ||
//          $50 == "alt_allele_in_normal" ||
//          $50 == "normal_lod" ||
//          $50 == "normal_lod,alt_allele_in_normal" ||
//          $50 == "alt_allele_in_normal,normal_lod"' \
//         !{pairName}.MuTect1.call_stats.txt \
//         > !{pairName}.PASS+TiNrisk.call_stats.txt


//     # -----------------------------------------------------------------------
//     # Output assertions
//     # -----------------------------------------------------------------------

//     test -s !{pairName}.MuTect1.call_stats.txt
//     test -s !{pairName}.MuTect1.vcf
//     test -s !{pairName}.MuTect1.PASS.vcf
//     test -s !{pairName}.mutectv1.final.vcf
//     '''
// }


// // ---------------------------------------------------------------------------
// // WORKFLOW
// // ---------------------------------------------------------------------------

// workflow {

//     if (!params.mutect1_runs) {
//         error "params.mutect1_runs is empty -- did the Cirro preprocess.py hook run? (see preprocess.py)"
//     }

//     if (!params.ref_fasta) {
//         error "Missing required param: ref_fasta"
//     }


//     // -----------------------------------------------------------------------
//     // Shared matched normal
//     // -----------------------------------------------------------------------

//     params.mutect1_runs.each { run ->

//         if (!run.normal_reads || !run.normal_reads_index) {
//             error "Run ${run.output_prefix} is missing normal_reads/normal_reads_index -- this pipeline requires a matched normal."
//         }
//     }

//     shared_normal_paths =
//         params.mutect1_runs
//             .collect { it.normal_reads }
//             .unique()

//     shared_normal_index_paths =
//         params.mutect1_runs
//             .collect { it.normal_reads_index }
//             .unique()

//     if (shared_normal_paths.size() != 1) {
//         error "Expected exactly one shared normal BAM across params.mutect1_runs, found: ${shared_normal_paths}"
//     }

//     if (shared_normal_index_paths.size() != 1) {
//         error "Expected exactly one shared normal BAM index across params.mutect1_runs, found: ${shared_normal_index_paths}"
//     }


//     normal_bam =
//         file(
//             shared_normal_paths[0],
//             checkIfExists: true
//         )

//     normal_bai =
//         file(
//             shared_normal_index_paths[0],
//             checkIfExists: true
//         )


//     // -----------------------------------------------------------------------
//     // Reference
//     // -----------------------------------------------------------------------

//     ref_fasta =
//         file(
//             params.ref_fasta,
//             checkIfExists: true
//         )

//     ref_fai =
//         file(
//             params.ref_fai,
//             checkIfExists: true
//         )

//     ref_dict =
//         file(
//             params.ref_dict,
//             checkIfExists: true
//         )


//     // -----------------------------------------------------------------------
//     // Optional target list
//     // -----------------------------------------------------------------------

//     target_list =
//         params.target_list
//             ? file(params.target_list, checkIfExists: true)
//             : NO_TARGET_LIST


//     // -----------------------------------------------------------------------
//     // dbSNP + index
//     // -----------------------------------------------------------------------

//     dbsnp =
//         params.dbsnp
//             ? file(params.dbsnp, checkIfExists: true)
//             : NO_DBSNP

//     dbsnpIdx =
//         params.dbsnp_idx
//             ? file(params.dbsnp_idx, checkIfExists: true)
//             : NO_DBSNP_IDX


//     // -----------------------------------------------------------------------
//     // COSMIC + index
//     // -----------------------------------------------------------------------

//     cosmic =
//         params.cosmic
//             ? file(params.cosmic, checkIfExists: true)
//             : NO_COSMIC

//     cosmicIdx =
//         params.cosmic_idx
//             ? file(params.cosmic_idx, checkIfExists: true)
//             : NO_COSMIC_IDX


//     // -----------------------------------------------------------------------
//     // Other optional resources
//     // -----------------------------------------------------------------------

//     rgBlacklist =
//         params.read_group_blacklist
//             ? file(
//                 params.read_group_blacklist,
//                 checkIfExists: true
//             )
//             : NO_RG_BLACKLIST

//     normalPanel =
//         params.normal_panel
//             ? file(
//                 params.normal_panel,
//                 checkIfExists: true
//             )
//             : NO_NORMAL_PANEL

//     normalPanelIdx =
//         params.normal_panel_idx
//             ? file(
//                 params.normal_panel_idx,
//                 checkIfExists: true
//             )
//             : NO_NORMAL_PANEL_IDX


//     // -----------------------------------------------------------------------
//     // Tumor runs
//     // -----------------------------------------------------------------------

//     runs_ch =
//         Channel
//             .fromList(params.mutect1_runs)
//             .map { run ->

//                 tuple(
//                     run.output_prefix,
//                     file(
//                         run.tumor_reads,
//                         checkIfExists: true
//                     ),
//                     file(
//                         run.tumor_reads_index,
//                         checkIfExists: true
//                     )
//                 )
//             }


//     // -----------------------------------------------------------------------
//     // ContEst
//     // -----------------------------------------------------------------------

//     if (
//         params.contest_target_intervals &&
//         params.snp6_bed &&
//         params.hapmap_vcf
//     ) {

//         contest_hapmap_vcf =
//             file(
//                 params.hapmap_vcf,
//                 checkIfExists: true
//             )

//         INDEX_CONTEST_VCF(
//             contest_hapmap_vcf
//         )

//         CONTEST(
//             runs_ch,
//             normal_bam,
//             normal_bai,
//             file(
//                 params.contest_target_intervals,
//                 checkIfExists: true
//             ),
//             ref_fasta,
//             ref_fai,
//             ref_dict,
//             file(
//                 params.snp6_bed,
//                 checkIfExists: true
//             ),
//             INDEX_CONTEST_VCF.out.indexed_vcf
//         )

//         frac_ch =
//             CONTEST
//                 .out
//                 .frac
//                 .map {
//                     pairName,
//                     fracFile ->

//                     tuple(
//                         pairName,
//                         fracFile.text.trim() as Float
//                     )
//                 }

//     } else {

//         frac_ch =
//             runs_ch
//                 .map {
//                     pairName,
//                     t_bam,
//                     t_bai ->

//                     tuple(
//                         pairName,
//                         params.fracContam
//                     )
//                 }
//     }


//     runs_with_frac =
//         runs_ch.join(frac_ch)

//     // tuple(
//     //     pairName,
//     //     t_bam,
//     //     t_bai,
//     //     fracContam
//     // )


//     // -----------------------------------------------------------------------
//     // Scatter
//     // -----------------------------------------------------------------------

//     SPLIT_INTERVALS(
//         runs_ch,
//         params.ref,
//         target_list
//     )

//     shards_per_pair =
//         SPLIT_INTERVALS
//             .out
//             .interval_files
//             .transpose()

//     // tuple(
//     //     pairName,
//     //     single_interval_file
//     // )


//     // -----------------------------------------------------------------------
//     // Join intervals with tumor + contamination data
//     // -----------------------------------------------------------------------

//     mutect1_inputs =
//         shards_per_pair
//             .combine(
//                 runs_with_frac,
//                 by: 0
//             )
//             .map {
//                 pairName,
//                 intervalFile,
//                 t_bam,
//                 t_bai,
//                 fracContam ->

//                 tuple(
//                     pairName,
//                     intervalFile,
//                     fracContam,
//                     t_bam,
//                     t_bai
//                 )
//             }


//     // -----------------------------------------------------------------------
//     // MuTect1
//     // -----------------------------------------------------------------------

//     MUTECT1(
//         mutect1_inputs,

//         normal_bam,
//         normal_bai,

//         ref_fasta,
//         ref_fai,
//         ref_dict,

//         dbsnp,
//         dbsnpIdx,

//         cosmic,
//         cosmicIdx,

//         rgBlacklist,
//         normalPanel,
//         normalPanelIdx
//     )


//     // -----------------------------------------------------------------------
//     // Gather
//     // -----------------------------------------------------------------------

//     gathered =
//         MUTECT1
//             .out
//             .shard
//             .groupTuple(
//                 by: [0, 1, 2]
//             )

//     // tuple(
//     //     pairName,
//     //     tumorSampleName,
//     //     normalSampleName,
//     //     [call_stats...],
//     //     [vcf...]
//     // )

//     GATHER_AND_FILTER(
//         gathered
//     )
// }
