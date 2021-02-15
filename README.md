# CNV_WP3

#Instructions on how to run the pipeline
cd /projects/wp3/nobackup/TWIST/OUTBOX/TEXX_XXX/CNV_WP3
git clone https://github.com/jonca79/CNV_WP3.git
cd CNV_WP3
module add snakemake
module add slurm-drmaa
module add singularity
snakemake -p -j 64 --drmaa "-A wp1 -p core -n {cluster.n} -t {cluster.time}"  -s ./CNV_WP3.smk --use-singularity --singularity-args "--bind /data --bind /beegfs-storage --bind /scratch " --cluster-config Config/Slurm/cluster.json
