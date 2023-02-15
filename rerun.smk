import os
import sys
import pandas

# parse tsv file. Only process complete samples
def parse_samples(samples_tsv):
    return pandas.read_csv(samples_tsv, sep='\t').set_index("id", drop=False).dropna()

# get reads based on sample id
def get_sample_id(sample_df, wildcards, col):
    return sample_df.loc[wildcards.sample, [col]][0]

# load sample file
_samples = parse_samples("sample_rerun.txt")


SAMPLES = list(_samples["id"])
print (SAMPLES)

rule all:
    input:
        expand("{result_dir}/{sample}.mp3.profile", result_dir = config["assay"]["profile"], sample = SAMPLES),
        expand("{result_dir}/{sample}", result_dir = config["assembly"], sample = SAMPLES),
        expand("{result_dir}/{sample}_coverage.txt", result_dir = config["contigs"], sample = SAMPLES),
        expand("{result_dir}/{sample}/{sample}_metabat_contigs_list", result_dir = config["bins"], sample = SAMPLES),
        expand("{result_dir}/{sample}/{sample}_maxbin_contigs_list", result_dir = config["bins"], sample = SAMPLES),
        expand("{result_dir}/{sample}/{sample}_DASTool_summary.tsv", result_dir = config["bins"], sample = SAMPLES),
        expand("{result_dir}/{sample}/{sample}_checkM.txt", result_dir = config["bins"], sample = SAMPLES),
        expand("{result_dir}/{sample}/{sample}_bins_taxonomy", result_dir = config["bins"], sample = SAMPLES),
        expand("{result_dir}/{sample}/{sample}_genefamilies.tsv", result_dir = config["function"], sample = SAMPLES),
        expand("{result_dir}/{sample}/squeezeMeta_{sample}", result_dir = config["function"], sample = SAMPLES)


localrules: merge_fq,contigsMod, coverage

rule metaphlan:
    input:
        r1 = lambda wildcards: get_sample_id(_samples, wildcards, "r1"),
        r2 = lambda wildcards: get_sample_id(_samples, wildcards, "r2")
    output:
        mpa = os.path.join(config["assay"]["profile"], "{sample}.mp3.profile"),
        vir = os.path.join(config["assay"]["profile"], "{sample}.mp3_vir.profile")
    threads:
        config["params"]["metaphlan3"]["threads"]
    conda:
        "envs/metaphlan.yaml"
    params:
        bowtie2db = config["params"]["metaphlan3"]["bowtie2db"],
        index = config["params"]["metaphlan3"]["index"],
        bw2 = os.path.join(config["assay"]["profile"], "{sample}.mp3.bw2.bz2"),
    shell:
        """
        metaphlan {input.r1},{input.r2} --bowtie2db {params.bowtie2db} --index {params.index} --nproc {threads} --input_type fastq --bowtie2out {params.bw2} -t rel_ab_w_read_stats > {output.mpa}
        metaphlan {params.bw2} --bowtie2db {params.bowtie2db} --index {params.index} --nproc {threads} --input_type bowtie2out --add_viruses -t rel_ab_w_read_stats > {output.vir}
                                                                                                                             """

rule assembly:
    input:
        read1 = lambda wildcards: get_sample_id(_samples, wildcards, "r1"),
        read2 = lambda wildcards: get_sample_id(_samples, wildcards, "r2")
    output:
        directory(os.path.join(config["assembly"], "{sample}"))
    threads: 80
    conda:
        "envs/megahit.yaml"
    shell:
        """
        ## for local use low memory
        megahit -1 {input.read1} -2 {input.read2} -o {output} --k-min 25 --k-max 131 --k-step 10 -t {threads} -m 0.2 --mem-flag 1
        """

rule contigsMod:
    input:
        os.path.join(config["assembly"], "{sample}")
    output:
        os.path.join(config["contigs"], "{sample}_contigs.fa")
    shell:
        """
        sed 's/k131/{wildcards.sample}/g' {input}/final.contigs.fa > {output}
        """

rule contigIndex:
    input:
        os.path.join(config["contigs"], "{sample}_contigs.fa")
    conda:
        "envs/bwa.yaml"
    output:
        temp(multiext(os.path.join(config["contigs"],"{sample}_contigs.fa"),".amb",".ann",".bwt",".pac",".sa"))
    shell:
        """
        bwa index {input}
        """

rule mapRaws:
    input:
        read1 = lambda wildcards: get_sample_id(_samples, wildcards, "r1"),
        read2 = lambda wildcards: get_sample_id(_samples, wildcards, "r2"),
        idx = multiext(os.path.join(config["contigs"],"{sample}_contigs.fa"),".amb",".ann",".bwt",".pac",".sa")
    output:
        temp(os.path.join(config["contigs"],"{sample}.sam")),
    conda:
        "envs/bwa.yaml"
    params:
        os.path.join(config["contigs"], "{sample}_contigs.fa")
    threads: 24
    shell:
        """
        bwa mem -t {threads} {params} {input.read1} {input.read2} > {output}
        """

rule samBam:
    input:
        os.path.join(config["contigs"],"{sample}.sam")
    output:
        bam = temp(os.path.join(config["contigs"],"{sample}.bam")),
        bai = temp(os.path.join(config["contigs"],"{sample}.bam.bai"))
    conda:
        "envs/samtools.yaml"
    params:
        bam = temp(os.path.join(config["contigs"],"{sample}.bam")),
    threads: 24
    shell:
        """
        samtools view -@ {threads} -hb -F4 {input} | samtools sort -o {params.bam}
        samtools index -@ {threads} {params.bam}
        """

rule coverage:
    input:
        os.path.join(config["contigs"],"{sample}.bam")
    output:
        os.path.join(config["contigs"],"{sample}_coverage.txt")
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools depth {input} | awk '{{ sum+=$3 }} END {{print "Average coverage:","\t",sum/NR}}' > {output}
        """

rule metabat2:
    input:
        bam = os.path.join(config["contigs"],"{sample}.bam"),
        bai = os.path.join(config["contigs"],"{sample}.bam.bai"),
        fa = os.path.join(config["contigs"], "{sample}_contigs.fa")
    output:
        list = os.path.join(config["bins"],"{sample}","{sample}_metabat_contigs_list"),
        depth = os.path.join(config["bins"],"{sample}","{sample}_contigs_metabat_depth.txt"),
    params:
        bins = os.path.join(config["bins"],"{sample}","{sample}_metabat","{sample}_metabat"),
        dir = os.path.join(config["bins"],"{sample}","{sample}_metabat")
    threads: 24
    conda:
        "envs/metabat.yaml"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}
        metabat2 -i {input.fa} -a {output.depth} -o {params.bins} -t {threads} -v
        awk '/>/{{sub(">","&"FILENAME"@");sub(/\.fa/,x)}}1'  {params.dir}/*.fa | grep ">" | awk -F"@" '{{print $2,$1}}' | sed 's/ .*\//\t/g' > {output.list}
        """

rule maxbin2:
    input:
        fa = os.path.join(config["contigs"], "{sample}_contigs.fa"),
        depth = os.path.join(config["bins"],"{sample}","{sample}_contigs_metabat_depth.txt"),
    output:
        os.path.join(config["bins"],"{sample}","{sample}_maxbin_contigs_list"),
    params:
        bins = os.path.join(config["bins"],"{sample}","{sample}_maxbin", "{sample}.maxbin"),
        dir = os.path.join(config["bins"],"{sample}","{sample}_maxbin"),
        maxdepth = os.path.join(config["bins"],"{sample}","{sample}_contigs_maxbin_depth.txt")
    threads: 24
    conda:
        "envs/maxbin.yaml"
    shell:
        """
        awk -F"\t" '{{print $1"\t"$4}}' {input.depth} | sed '1d' > {params.maxdepth}
        mkdir -p {params.dir}
        run_MaxBin.pl -thread {threads} -contig {input.fa} -out {params.bins} -abund {params.maxdepth}
        awk '/>/{{sub(">","&"FILENAME"@");sub(/\.fasta/,x)}}1'  {params.dir}/*.fasta | grep ">" | awk -F"@" '{{print $2,$1}}' | sed 's/ .*\//\t/g' > {output}
        """

rule das_tools:
    input:
        fa = os.path.join(config["contigs"], "{sample}_contigs.fa"),
        maxbinBins = os.path.join(config["bins"],"{sample}","{sample}_maxbin_contigs_list"),
        metabatBins = os.path.join(config["bins"],"{sample}","{sample}_metabat_contigs_list")
    output:
        os.path.join(config["bins"],"{sample}","{sample}_DASTool_summary.tsv"),
    params:
        dir = os.path.join(config["bins"],"{sample}","{sample}"),
    threads: 24
    conda:
        "envs/dastools.yaml"
    shell:
        """
        DAS_Tool -i {input.maxbinBins},{input.metabatBins} -l Maxbin,Metabat -c {input.fa} -o {params.dir} --write_bins -t {threads}
        """

rule checkMdas:
    input:
        os.path.join(config["bins"],"{sample}","{sample}_DASTool_summary.tsv"),
    output:
        das = os.path.join(config["bins"],"{sample}","{sample}_checkM.txt"),
        tmpdas = temp(directory(os.path.join(config["bins"],"{sample}","{sample}_dastool_tmp"))),
    params:
        das = os.path.join(config["bins"],"{sample}","{sample}_DASTool_bins"),
        dasComplete = os.path.join(config["bins"],"{sample}","{sample}_dastool_checkm"),
        dasMarker = os.path.join(config["bins"],"{sample}","{sample}_dastool.marker"),
    threads: 24
    conda:
        "envs/checkm.yaml"
    shell:
        """
        mkdir -p {output.tmpdas}

        checkm tree -t {threads} -x fa {params.das} {params.dasComplete}
        checkm lineage_set {params.dasComplete} {params.dasMarker}
        checkm analyze -t {threads} -x fa {params.dasMarker} {params.das} {params.dasComplete} --tmpdir {output.tmpdas}
        checkm qa -t {threads} {params.dasMarker} {params.dasComplete} -f {output.das}
        """


rule binTaxo:
    input:
        os.path.join(config["bins"],"{sample}","{sample}_DASTool_summary.tsv"),
    output:
        directory(os.path.join(config["bins"],"{sample}","{sample}_bins_taxonomy")),
    params:
        export = config["params"]["binTaxo"],
        dir = os.path.join(config["bins"],"{sample}","{sample}_DASTool_bins"),
    threads: 8
    conda:
        "envs/gtdbkt.yaml"
    shell:
        """
        export GTDBTK_DATA_PATH={params.export}
        gtdbtk classify_wf --genome_dir {params.dir} --out_dir {output} --cpus {threads} --extension fa
        """

rule merge_fq:
    input:
        read1 = lambda wildcards: get_sample_id(_samples, wildcards, "r1"),
        read2 = lambda wildcards: get_sample_id(_samples, wildcards, "r2"),
    output:
        temp(os.path.join(config["function"],"{sample}.fq.gz"))
    shell:
        """
        cat {input} > {output}
        """

rule humann4:
    input:
        os.path.join(config["function"],"{sample}.fq.gz")
    output:
        genefam = os.path.join(config["function"],"{sample}","{sample}_genefamilies.tsv"),
        pathabund = os.path.join(config["function"],"{sample}","{sample}_pathabundance.tsv"),
        pathcov = os.path.join(config["function"],"{sample}","{sample}_pathcoverage.tsv"),
        tmpdir = temp(directory(os.path.join(config["function"],"{sample}","{sample}_humann_temp"))),
    conda:
        "envs/humann4.yaml"
    threads: 16
    params:
        outdir = os.path.join(config["function"],"{sample}"),
        nu_db = config["params"]["humann4"]["nu_db"],
        aa_db = config["params"]["humann4"]["aa_db"],
        metaphlan_command = "'--index mpa_vJan21_CHOCOPhlAnSGB_202103 --bowtie2db /ifs/data/mmbresearch/metagenomic_databases/metaphlan4_db'"
    shell:
        """
        humann -i {input} --nucleotide-database {params.nu_db} --protein-database {params.aa_db} -o {params.outdir} --metaphlan-options {params.metaphlan_command} --threads {threads}
        """

rule squeezeIn:
    input:
        read1 = lambda wildcards: get_sample_id(_samples, wildcards, "r1"),
        read2 = lambda wildcards: get_sample_id(_samples, wildcards, "r2"),
    output:
        os.path.join(config["function"],"{sample}","{sample}_squeezeSamples.txt")
    shell:
        """
        echo "{wildcards.sample} $(basename {input.read1}) pair1" | sed 's/ /\t/g' > {output}
        echo "{wildcards.sample} $(basename {input.read2}) pair2" | sed 's/ /\t/g' >> {output}
        """

rule squeezemeta:
    input:
        fa = os.path.join(config["contigs"], "{sample}_contigs.fa"),
        squeezeIn = os.path.join(config["function"],"{sample}","{sample}_squeezeSamples.txt")
    output:
        squeeze = directory(os.path.join(config["function"],"{sample}","squeezeMeta_{sample}")),
    params:
        dir = directory(os.path.join(config["function"],"{sample}","squeezeMeta_{sample}")),
        squeezeRaw = directory(os.path.join(config["function"],"{sample}","squeezeMeta_{sample}/data")),
        db = directory(config["params"]["squeezemeta"]["db"]),
        pathRaw = config["assay"]["rmhost"]
    conda:
        "envs/squeeze.yaml"
    threads: 16
    shell:
        """
        configure_nodb_alt.pl {params.db}
        SqueezeMeta.pl -m coassembly -p {params.dir} -t {threads} -s {input.squeezeIn} -f {params.pathRaw} -extassembly {input.fa}
        rm -r {params.squeezeRaw}
        """
