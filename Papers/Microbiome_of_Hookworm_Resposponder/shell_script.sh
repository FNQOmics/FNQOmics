#loads of external stuff needed
1) rebuild kraken2 (first on copyq for downloads / then in big compute)
kraken2-build --standard --db krakenStandard (or run kraken2-build --download-library LIB --db krakenStandard for human, plasmid, etc)
kraken2-build --standard --threads 32 --db krakenStandard

#Checkm
2) conda create -n checkm python=3.9
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
checkm data setRoot /g/data/u86/software/checkm

#nt blast
3) Run not on copyq as downloads corrupted for some reason
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/v4/nt_v4.*.tar.gz"

#BMTagger
4) module load blast/2.11.0
bmtool -d hg38.fa -o hg38.bitmask
srprism mkindex -i hg38.fa -o hg38.srprism -M 100000

5) Set config paths (/g/data/u86/software/metaWRAP/bin/config-metawrap)


Data
1) Standardize
for f in `cat sample_list`; do cd $f; ln -s ${f}_v2*_L001_R1_001.fastq.gz ${f}_v2_R1.fastq.gz;cd ..; done
for f in `cat sample_list`; do cd $f; ln -s ${f}_v2*_L001_R2_001.fastq.gz ${f}_v2_R2.fastq.gz;cd ..; done
for f in `cat sample_list`; do cd $f; ln -s ${f}_v5*_L001_R1_001.fastq.gz ${f}_v5_R1.fastq.gz;cd ..; done
for f in `cat sample_list`; do cd $f; ln -s ${f}_v5*_L001_R2_001.fastq.gz ${f}_v5_R2.fastq.gz;cd ..; done


Run metawrap 

export PATH="/g/data/u86/software/metaWRAP/bin/:${PATH}"
startconda
conda activate metawrap-env

#Read qc
metawrap read_qc -1 LAC018_v2_S7_L001_R1_001.fastq.gz -2 LAC018_v2_S7_L001_R2_001.fastq.gz -o v2/ -t 48

#Bulk
for f in `cat sample_list`; do cat metawrap_qc_template.qsub | sed -e "s/SAMPLE/${f}/g" > qsub/${f}_qc.qsub; done

#Assembly try metaspades
metawrap assembly -1 v2/LAC018_v2_S7_L001_R1_001_val_1.fq.gz -2 v2/LAC018_v2_S7_L001_R2_001_val_2.fq.gz -m 192 -t 48 --metaspades -o v2/

To submit jobs (in qsub):
source /g/data/u86/mxf221/miniconda3/bin/activate
conda activate metawrap-env
export PATH="/g/data/u86/software/metaWRAP/bin/:${PATH}"

#Generate biom from kraken2
#Rerun with --report
for f in `cat sample_list_v5`; do echo $f"_v5"; ~/g/software/kraken2/kraken2 --db ~/g/software/kraken2/krakenStandard/ --paired --report /g/data/u86/mxf221/hookworm/shotgun_MB/${f}/v5/KRAKEN/report/${f}_v5.report --output /g/data/u86/mxf221/hookworm/shotgun_MB/${f}/v5/KRAKEN/report/${f}_v5.krak2 --threads 24 /g/data/u86/mxf221/hookworm/shotgun_MB/${f}/v5/${f}_v5_R1_val_1.fastq /g/data/u86/mxf221/hookworm/shotgun_MB/${f}/v5/${f}_v5_R1_val_2.fastq; done

python2 ~/software/microbiome/kraken-biom/kraken_biom.py -m /drive3/work/microbiome/shotgun_T2D/meta.tsv  ALM023_v2.report ALM023_v5.report BEW029_v2.report BEW029_v5.report CEB013_v2.report CEB013_v5.report C_S221_v2.report C_S221_v5.report DAT016_v2.report DAT016_v5.report DLC028_v2.report DLC028_v5.report EKS030_v2.report EKS030_v5.report HAC056_v2.report HAC056_v5.report H_M015_v2.report H_M015_v5.report JSB049_v2.report JSB049_v5.report J_W010_v2.report J_W010_v5.report KMM057_v2.report KMM057_v5.report LAC018_v2.report LAC018_v5.report LAM205_v2.report LAM205_v5.report LCF208_v2.report LCF208_v5.report LCJ212_v2.report LCJ212_v5.report LRL008_v2.report LRL008_v5.report MEL019_v2.report MEL019_v5.report MJD223_v2.report MJD223_v5.report MJM017_v2.report MJM017_v5.report MRP214_v2.report MRP214_v5.report MRS215_v5.report PSC005_v2.report PSC005_v5.report RLC048_v2.report RLC048_v5.report SAP201_v5.report SBM218_v2.report SBM218_v5.report SLB203_v2.report SLB203_v5.report TTW209_v2.report TTW209_v5.report VLM026_v2.report

mv table.biom reads.biom (and to reassembly.biom)

Try reads/assembly/reassembly
#Create reassembly
for f in `cat sample_list_v5`; do echo $f; cat ${f}/v5/BIN_REASSEMBLY/reassembled_bins/*fa >> ${f}/v5/BIN_REASSEMBLY/reassembled_bins/${f}_v5.fa; ~/g/software/kraken2/kraken2 --use-names --db /g/data/u86/software/kraken2/krakenStandard --threads 24 --output ${f}/v5/BIN_REASSEMBLY/reassembled_bins/${f}_v5.krak2 --report ${f}/v5/BIN_REASSEMBLY/reassembled_bins/${f}_v5.report ${f}/v5/BIN_REASSEMBLY/reassembled_bins/${f}_v5.fa; done

#Create tsv bioms
biom convert -i reassembly.biom --to-tsv -o reassembly.tsv --table-type="OTU table" --header-key taxonomy

module load python3/3.12.1

module load intel-mkl/2023.2.0
