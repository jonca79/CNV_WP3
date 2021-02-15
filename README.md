# CNV_WP3

#Instructions on how to run the pipeline <br>
cd /projects/wp3/nobackup/TWIST/OUTBOX/TEXX_XXX/CNV_WP3 <br>
git clone https://github.com/jonca79/CNV_WP3.git <br>
cd CNV_WP3 <br>
module add snakemake <br>
module add slurm-drmaa <br>
module add singularity <br>
snakemake -p -j 64 --drmaa "-A wp1 -p core -n {cluster.n} -t {cluster.time}"  -s ./CNV_WP3.smk --use-singularity --singularity-args "--bind /data --bind /beegfs-storage --bind /scratch " --cluster-config Config/Slurm/cluster.json <br>
