
rule Call_CNV:
    input:
        samples="CNV/samples.txt",
        cov="CNV/Run.cov",
    output:
        cnv="CNV/calls.tsv",
    params:
        plotsdir="CNV/CNV_plots/"
    singularity:
        "/projects/wp2/nobackup/Twist_Myeloid/Containers/python3.6.0-pysam-xlsxwriter.simg"
    shell:
        "src/scripts/python/CNV_calling_Jonas_WP3.py {input.samples} {input.cov} {output.cnv} {params.plotsdir}", shell=True)
