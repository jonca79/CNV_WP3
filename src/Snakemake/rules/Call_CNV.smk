
rule Call_CNV:
    input:
        samples="CNV/samples.txt",
        cov="CNV/Run.cov",
        normal_samples="DATA/samples_TE6_37.txt",
        normal_cov="DATA/TE6_37_all_panels.cov",
        normal_cnv1="DATA/WW_25m_CNV_ChAS3.0.aed",
    output:
        cnv="CNV/calls.tsv",
    params:
        plotsdir="CNV/CNV_plots/",
    singularity:
        "/projects/wp2/nobackup/Twist_Myeloid/Containers/python3.6.0-pysam-xlsxwriter.simg"
    shell:
        "python3 src/scripts/python/CNV_calling_Jonas_WP3.py "
        "{input.samples} {input.cov} {input.normal_samples} {input.normal_cov} {input.normal_cnv1} "
        "{output.cnv} {params.plotsdir}"
