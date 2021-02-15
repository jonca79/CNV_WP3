
rule all:
    input:
        "CNV/calls.tsv",
        "CNV/Run.cov",
        "CNV/samples.txt",


include: "src/Snakemake/workflow/CNV_WP3_workflow.smk"
