
rule all:
    input:
        cnv="CNV/calls.tsv",


include: "src/Snakemake/workflow/CNV_WP3_workflow.smk"
