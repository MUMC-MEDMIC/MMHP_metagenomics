##############################################################
### giang.le@mumc.nl
### MMHP metagenomics analysis
### ref: https://github.com/ohmeta/metapi/blob/master/metapi/Snakefile
##############################################################

import os
import sys
import pandas

shell.executable("bash")
configfile: "config.yaml"

def parse_samples(samples_tsv):
    return pandas.read_csv(samples_tsv, sep='\t').set_index("id", drop=False)

def get_sample_id(sample_df, wildcards, col):
    return sample_df.loc[wildcards.sample, [col]].dropna()[0]

_samples = parse_samples("sample.txt")

localrules: filter_summary

SAMPLES = list(pandas.read_csv("sample.txt", sep='\t').set_index("id", drop=False)["id"])
print(SAMPLES)

rule all:
    input:
#        expand("{result_dir}/filter_summary.txt", result_dir = config["results"]),
        expand("{result_dir}/{sample}_coverage.txt", result_dir = config["contigs"], sample = SAMPLES),
        expand("{result_dir}/{sample}/das_tool_taxonomy", result_dir =  config["bins"], sample = SAMPLES ),
        expand("{result_dir}/{sample}/{sample}_genefamilies.tsv", result_dir = config["function"], sample = SAMPLES)

localrules: filter_summary, contigsMod, merge_fq


### trimming & remove host reads
rule filter:
    input:
        r1 = lambda wildcards: get_sample_id(_samples, wildcards, "fq1"),
        r2 = lambda wildcards: get_sample_id(_samples, wildcards, "fq2")
    output:
        html = os.path.join(config["assay"]["trimming"], "{sample}.fastp.html"),
        json = os.path.join(config["assay"]["trimming"], "{sample}.fastp.json"),
        rmhost_r1 = protected(os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.1.fq.gz")),
        rmhost_r2 = protected(os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.2.fq.gz"))
    params:
        trim_r1 = temp(os.path.join(config["assay"]["trimming"], "{sample}.trimmed.1.fq.gz")),
        trim_r2 = temp(os.path.join(config["assay"]["trimming"], "{sample}.trimmed.2.fq.gz")),
        sam = temp(os.path.join(config["assay"]["trimming"], "{sample}.sam")),
        min_len = config["params"]["fastp"]["min_len"],
        ad_r1 = config["params"]["fastp"]["adapter_r1"],
        ad_r2 = config["params"]["fastp"]["adapter_r2"],
        n_lim = config["params"]["fastp"]["n_base_limit"],
        index = config["params"]["rmhost"]["bowtie2_index"],
        maxins = config["params"]["rmhost"]["bowtie2_maxins"]
    threads:
        config["params"]["rmhost"]["threads"]
    conda:
        "envs/filter.yaml"
    log:
        fastp_log = os.path.join(config["logs"]["trimming"], "{sample}.fastp.log"),
        bowtie2_log = os.path.join(config["logs"]["rmhost"], "{sample}.rmhost.log")
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {params.trim_r1} -O {params.trim_r2} -w {threads} --n_base_limit {params.n_lim} --cut_front --cut_tail --length_required {params.min_len} --adapter_sequence={params.ad_r1} --adapter_sequence_r2={params.ad_r2} -j {output.json} -h {output.html} 2> {log.fastp_log}


        bowtie2 --end-to-end --very-sensitive -p {threads} -I 0 -X {params.maxins} -x {params.index} --mm -1 {params.trim_r1} -2 {params.trim_r2} > {params.sam} 2> {log.bowtie2_log}

        rm {params.trim_r1} {params.trim_r2} 

        samtools fastq -N -c 5 -f 12 -F 256 -1 {output.rmhost_r1} -2 {output.rmhost_r2} {params.sam}

        rm {params.sam}
        """

### summary of filtered reads
rule seqkit_stat:
    input:
        expand("{rmhost_log_dir}/{{sample}}.rmhost.{reads}.fq.gz", rmhost_log_dir = config["assay"]["rmhost"], reads = ["1","2"])
    output:
        os.path.join(config["logs"]["rmhost"], "{sample}.rmhost.reads.summary")
    conda:
        "envs/filter.yaml"
    shell:
        "seqkit stat {input} > {output}"

rule filter_summary:
    input:
        trim = expand("{trim_res}/{sample}.fastp.json", trim_res = config["assay"]["trimming"], sample = _samples.index),
        rmhost = expand("{rmhost_res}/{sample}.rmhost.reads.summary", rmhost_res = config["logs"]["rmhost"], sample = _samples.index)
    output:
        #protected(os.path.join(config["results"], "filter_summary.txt"))
        os.path.join(config["results"], "filter_summary.txt")
    params:
        trim_summary = temp(os.path.join(config["results"], "trim_summary.txt")),
        rmhost_summary = temp(os.path.join(config["results"], "rmhost_summary.txt"))
    shell:
        """
        python rules/filter_summary.py -t {input.trim} > {params.trim_summary}
        python rules/filter_summary.py -r {input.rmhost} > {params.rmhost_summary}
        python rules/merge_summary.py {params.trim_summary} {params.rmhost_summary} {output}
        rm {params.trim_summary} {params.rmhost_summary}
        """

### metaphlan3 taxonomic using markers
rule metaphlan3:
    input:
        r1 = os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.1.fq.gz"),
        r2 = os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.2.fq.gz")
    output:
        #mpa = protected(os.path.join(config["assay"]["profile"], "{sample}.mp3.profile")),
        #vir = protected(os.path.join(config["assay"]["profile"], "{sample}.mp3_vir.profile"))
        mpa = os.path.join(config["assay"]["profile"], "{sample}.mp3.profile"),
        vir = os.path.join(config["assay"]["profile"], "{sample}.mp3_vir.profile")
    threads:
        config["params"]["metaphlan3"]["threads"]
    conda:
        "envs/metaphlan.yaml"
    params:
        bowtie2db = config["params"]["metaphlan3"]["bowtie2db"],
        index = config["params"]["metaphlan3"]["index"],
        #bw2 = protected(os.path.join(config["assay"]["profile"], "{sample}.mp3.bw2.bz2")),
        bw2 = os.path.join(config["assay"]["profile"], "{sample}.mp3.bw2.bz2"),
    shell:
        """
        metaphlan {input.r1},{input.r2} --bowtie2db {params.bowtie2db} --index {params.index} --nproc {threads} --input_type fastq --bowtie2out {params.bw2} -t rel_ab_w_read_stats > {output.mpa}
        metaphlan {params.bw2} --bowtie2db {params.bowtie2db} --index {params.index} --nproc {threads} --input_type bowtie2out --add_viruses -t rel_ab_w_read_stats > {output.vir}
        """

rule assembly:
    input:
        os.path.join(config["assay"]["profile"], "{sample}.mp3.profile"),
        read1 = os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.1.fq.gz"),
        read2 = os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.2.fq.gz")
    output:
        directory(os.path.join(config["assembly"], "{sample}"))
    threads: 90
    conda:
        "envs/megahit.yaml"
    shell:
        """
        megahit -1 {input.read1} -2 {input.read2} -o {output} --k-min 25 --k-max 131 --k-step 10 -t {threads}
        """

rule contigsMod:
    input:
        ancient(os.path.join(config["assembly"], "{sample}"))
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

#https://ndombrowski.github.io/Binning_tutorial/#do-the-read-mapping
rule mapRaws:
    input:
        read1 = os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.1.fq.gz"),
        read2 = os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.2.fq.gz"),
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
        samtools index {params.bam}
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
        batdone = os.path.join(config["bins"],"{sample}","{sample}_metabat_contigs_list"),
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
        das = os.path.join(config["bins"],"{sample}","{sample}_dastool_complete.txt"),
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
        directory(os.path.join(config["bins"],"{sample}","das_tool_taxonomy")),
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
        read1 = os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.1.fq.gz"),
        read2 = os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.2.fq.gz")
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


### step4 profile_summary
rule merge_profile:
    input:
        mpa = expand("{profile_dir}/{sample}.mp3.profile", profile_dir = config["assay"]["profile"], sample = _samples.index),
        vir = expand("{profile_dir}/{sample}.mp3_vir.profile", profile_dir = config["assay"]["profile"], sample = _samples.index)
    output:
        mpa_merge = protected(os.path.join(config['results'], "metaphlan3.profile.merge.txt")),
        mpa_merge_vir = protected(os.path.join(config['results'], "metaphlan3_vir.profile.merge.txt"))
    shell:
        '''
        python rules/merge_metaphlan_tables.py {input.mpa} > {output.mpa_merge}
        python rules/merge_metaphlan_tables.py {input.vir} > {output.mpa_merge_vir}
        '''

onerror:
    shell("mail -s 'Metagenomic analysis failed' giang.le@mumc.nl < {log}")
