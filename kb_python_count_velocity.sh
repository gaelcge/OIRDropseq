#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=15:00:00
#SBATCH --mem=250G
#SBATCH --job-name=kbpython
#SBATCH --account=def-jsjoyal
#SBATCH --mail-type=END
#SBATCH --mail-user=gael.cagnone.1@gmail.com
#SBATCH -o job.report
#SBATCH --no-requeue


source ~/ENV_python.3.7.4_scanpy/bin/activate

cd /home/gaelcge/projects/def-jsjoyal/gaelcge/Kb_python_Velocity

##BUIld matrices


pathtoindex="/home/gaelcge/projects/def-jsjoyal/gaelcge/References/kb_python_velo/mouse"

cd /home/gaelcge/projects/def-jsjoyal/gaelcge/Kb_python_Velocity/

rm dir_list.txt

#for d in /home/gaelcge/projects/def-jsjoyal/gaelcge/Sequencing/DropSeq_Fastq/bioinfo.iric.ca/seq/*/*/fastq/*;do [[ -d "$d" ]] && echo "$d" >> dir_list.txt; done

#for d in /home/gaelcge/projects/def-jsjoyal/gaelcge/Sequencing/DropSeq_Fastq/bioinfo.iric.ca/seq/e89b397f509746458fa6de9b8eb04155/*/fastq/Sample_Gael_180830_R5_wt_2400st_N706;do [[ -d "$d" ]] && echo "$d" >> dir_list.txt; done

for d in /home/gaelcge/projects/def-jsjoyal/gaelcge/Sequencing/DropSeq_Fastq/bioinfo.iric.ca/seq/Retina_Rytvela/H3YG7BGX5_H533FBGX5/fastq/*;do [[ -d "$d" ]] && echo "$d" >> dir_list.txt; done

cat dir_list.txt | while read fastq 
do
samplenames=$(basename "$fastq")
#gzip $fastq/*fastq
mkdir /home/gaelcge/projects/def-jsjoyal/gaelcge/Kb_python_Velocity/$samplenames
output="/home/gaelcge/projects/def-jsjoyal/gaelcge/Kb_python_Velocity/"$samplenames
kb count --h5ad \
-i $pathtoindex/index.idx \
-g $pathtoindex/t2g.txt -x DROPSEQ -o $output -c1 $pathtoindex/spliced_t2c.txt -c2 \
$pathtoindex/unspliced_t2c.txt --workflow lamanno --filter bustools -t 14 \
$fastq/*.gz
done



##























