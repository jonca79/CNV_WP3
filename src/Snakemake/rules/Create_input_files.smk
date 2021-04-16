

rule Make_coverage_bed:
    input:
        bed="DATA/All_panels.bed",
        #bed="DATA/Twist_Exome_Target_hg19.annotated.bed",
    output:
        cov="CNV/Run.cov",
    singularity:
        "/projects/wp4/nobackup/workspace/somatic_dev/singularity/bedtools2.29.2_samtools1.9.0_fgbio1.3.0.simg"
    shell:
        "bedtools multicov -bams $PWD/../BAM/*.bam -bed {input.bed} > {output.cov}"


rule Make_sample_file:
    output:
        samples="CNV/samples.txt",
    shell:
        "ls -d $PWD/../BAM/*.bam > {output.samples}"
