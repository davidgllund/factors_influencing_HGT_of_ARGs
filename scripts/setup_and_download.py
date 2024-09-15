import subprocess

def open_compressed():
    subprocess.run('tar -xf auxiliary_files.tar.gz', shell=True)

def download_genomes():
    subprocess.run('mkdir genomes', shell=True)
    subprocess.run('while read line; do wget $line; done<auxiliary_files/genomes_to_download.txt', shell=True)
    subprocess.run('mv *genomic.fna.gz genomes/', shell=True)
    subprocess.run('ls -d -1 genomes/* > auxiliary_files/genome_filepaths.txt', shell=True)

def download_metagenomes():
    subprocess.run('wget ftp://ftp.microbio.me/emp/release1/mapping_files/emp_qiime_mapping_qc_filtered.tsv', shell=True)
    subprocess.run('mv emp_qiime_mapping_qc_filtered.tsv auxiliary_files/', shell=True)

    subprocess.run('wget ftp://ftp.microbio.me/emp/release1/otu_tables/closed_ref_greengenes/emp_cr_gg_13_8.qc_filtered.biom', shell=True)
    subprocess.run('biom convert -i emp_cr_gg_13_8.qc_filtered.biom -o otus_gg_13_8.txt --to-tsv', shell=True)
    subprocess.run('mv otus_gg_13_8.txt auxiliary_files/', shell=True)
    subprocess.run('rm emp_cr_gg_13_8.qc_filtered.biom', shell=True)
    
    subprocess.run('wget gwmc.ou.edu/files/gwmc_16s_otutables/uparse/beforeresample/gwmc_16S_otutab.zip', shell=True)
    subprocess.run('unzip gwmc_16S_otutab.zip', shell=True)
    subprocess.run('mv GWMC_16S_otutab.txt auxiliary_files/', shell=True)
    subprocess.run('rm gwmc_16S_otutab.zip', shell=True)

def main():
    open_compressed()
    download_genomes()
    download_metagenomes()

if __name__ == '__main__':
    main()
