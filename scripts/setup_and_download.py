import subprocess

def open_compressed():
    subprocess.run('tar -xf auxiliary_files.tar.gz')

def download_genomes():
    subprocess.run('mkdir genomes')
    subprocess.run('while read line; do wget $line; done<auxiliary_files/genomes_to_download.txt')
    subprocess.run('mv *genomic.fna.gz genomes/')
    subprocess.run('ls -d -1 genomes/ > auxiliary_files/genome_filepaths.txt')

def download_metagenomes():
    subprocess.run('wget ftp://ftp.microbio.me/emp/release1/mapping_files/emp_qiime_mapping_qc_filtered.tsv')
    subprocess.run('mv emp_qiime_mapping_qc_filtered.tsv auxiliary_files/')

    subprocess.run('wget ftp://ftp.microbio.me/emp/release1/otu_tables/closed_ref_greengenes/emp_cr_gg_13_8.qc_filtered.biom')
    subprocess.run('biom convert -i emp_cr_gg_13_8.qc_filtered.biom -o otus_gg_13_8.txt --to-tsv')
    subprocess.run('mv otus_gg_13_8.txt auxiliary_files/')
    subprocess.run('rm emp_cr_gg_13_8.qc_filtered.biom')
    
    subprocess.run('wget gwmc.ou.edu/files/gwmc_16s_otutables/uparse/beforeresample/gwmc_16S_otutab.zip')
    subprocess.run('unzip gwmc_16S_otutab.zip')
    subprocess.run('mv GWMC_16S_otutab.txt auxiliary_files/')
    subprocess.run('rm gwmc_16S_otutab.zip')

def main():
    open_compressed()
    download_genomes()
    download_metagenomes()

if __name__ == '__main__':
    main()
