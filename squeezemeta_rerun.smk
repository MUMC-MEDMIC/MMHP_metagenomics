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
        # SqueezeMeta functional gene annotaiton from contigs
        expand("{result_dir}/{sample}/squeezeMeta_{sample}", result_dir = config["function"], sample = SAMPLES)

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
