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
#        expand("{result_dir}/metaphlan3.profile.merge.txt", result_dir = config["results"])
        expand("{result_dir}/{sample}", result_dir = config["assembly"], sample = SAMPLES)


### trimming & remove host reads
rule filter:
    input:
        r1 = lambda wildcards: get_sample_id(_samples, wildcards, "fq1"),
        r2 = lambda wildcards: get_sample_id(_samples, wildcards, "fq2")
    output:
        trim_r1 = temp(os.path.join(config["assay"]["trimming"], "{sample}.trimmed.1.fq.gz")),
        trim_r2 = temp(os.path.join(config["assay"]["trimming"], "{sample}.trimmed.2.fq.gz")),
        html = os.path.join(config["assay"]["trimming"], "{sample}.fastp.html"),
        json = os.path.join(config["assay"]["trimming"], "{sample}.fastp.json"),
        sam = temp(os.path.join(config["assay"]["trimming"], "{sample}.sam")),
        rmhost_r1 = protected(os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.1.fq.gz")),
        rmhost_r2 = protected(os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.2.fq.gz"))
    params:
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
        fastp -i {input.r1} -I {input.r2} -o {output.trim_r1} -O {output.trim_r2} -w {threads} --n_base_limit {params.n_lim} --cut_front --cut_tail --length_required {params.min_len} --adapter_sequence={params.ad_r1} --adapter_sequence_r2={params.ad_r2} -j {output.json} -h {output.html} 2> {log.fastp_log}

        bowtie2 --end-to-end --very-sensitive -p {threads} -I 0 -X {params.maxins} -x {params.index} --mm -1 {output.trim_r1} -2 {output.trim_r2} > {output.sam} 2> {log.bowtie2_log}

        samtools fastq -N -c 5 -f 12 -F 256 -1 {output.rmhost_r1} -2 {output.rmhost_r2} {output.sam}
        """

### summary of filtered reads
rule seqkit_stat:
    input:
        expand("{rmhost_log_dir}/{{sample}}.rmhost.{reads}.fq.gz", rmhost_log_dir = config["assay"]["rmhost"], reads = ["1","2"])
    output:
        os.path.join(config["logs"]["rmhost"], "{sample}.rmhost.reads.summary")
    conda:
        "envs/seqkit.yaml"
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
        read1 = os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.1.fq.gz"),
        read2 = os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.2.fq.gz")
    conda:
        "../envs/megahit.yaml"
    output:
        directory(os.path.join(config["assembly"], "{sample}"))
    shell:
        """
        megahit -1 {input.read1} -2 {input.read1} -o {output} --k-min 25 --k-max 131 --k-step 10 
        """

uule contigIndex:
    input:
        os.path.join(config["assembly"], "{sample}")
    conda:
        "envs/maxbin.yaml"
    params:
        "{sample}.bowtie"
    output:
        temp(multiext("{sample}.bowtie",".1.bt2",".2.bt2",".3.bt2",".4.bt2",".rev.1.bt2",".rev.2.bt2"))
    threads: 12
    shell:
        """
        bowtie2-build --quiet {input}/final.contigs.fa {params}
        """

rule mapRaws:
    input:
        read1 = lambda wildcards: SAMPLES[wildcards.sample]['forward'],
        read2 = lambda wildcards: SAMPLES[wildcards.sample]['reverse'],
        idx = multiext("{sample}.bowtie",".1.bt2",".2.bt2",".3.bt2",".4.bt2",".rev.1.bt2",".rev.2.bt2")
    output:
        temp("{sample}.sam")
    conda:
        "envs/maxbin.yaml"
    params:
        "{sample}.bowtie"
    threads: 12
    shell:
        """
        bowtie2 -x {params} -1 {input.read1} -2 {input.read2} --quiet -p {threads} -S {output}
        """

rule abundContigs:
    input:
        "{sample}.sam"
    output:
        "{sample}_contigAbundance.cov"
    conda:
        "envs/maxbin.yaml"
    params:
        temp("{sample}.cov")
    shell:
        """
        pileup.sh in={input} out={params}
        awk '{{print $1"\t"$5}}' {params} | grep -v '^#' > {output}
        """
uule contigIndex:
    input:
        "{sample}_assembly"
    conda:
        "envs/maxbin.yaml"
    params:
        "{sample}.bowtie"
    output:
        temp(multiext("{sample}.bowtie",".1.bt2",".2.bt2",".3.bt2",".4.bt2",".rev.1.bt2",".rev.2.bt2"))
    threads: 12
    shell:
        """
        bowtie2-build --quiet {input}/final.contigs.fa {params}
        """

rule mapRaws:
    input:
        read1 = lambda wildcards: SAMPLES[wildcards.sample]['forward'],
        read2 = lambda wildcards: SAMPLES[wildcards.sample]['reverse'],
        idx = multiext("{sample}.bowtie",".1.bt2",".2.bt2",".3.bt2",".4.bt2",".rev.1.bt2",".rev.2.bt2")
    output:
        temp("{sample}.sam")
    conda:
        "envs/maxbin.yaml"
    params:
        "{sample}.bowtie"
    threads: 12
    shell:
        """
        bowtie2 -x {params} -1 {input.read1} -2 {input.read2} --quiet -p {threads} -S {output}
        """

rule abundContigs:
    input:
        "{sample}.sam"
    output:
        "{sample}_contigAbundance.cov"
    conda:
        "envs/maxbin.yaml"
    params:
        temp("{sample}.cov")
    shell:
        """
        pileup.sh in={input} out={params}
        awk '{{print $1"\t"$5}}' {params} | grep -v '^#' > {output}
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
