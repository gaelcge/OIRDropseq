

source ~/ENV_python.3.7.4_scanpy/bin/activate

##Build index for RNA velocity
cd /home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_EC_kbpython/

PATHTOREFERENCE="/home/gaelcge/projects/def-jsjoyal/gaelcge/References/GRCm38"

pathtoindex="/home/gaelcge/projects/def-jsjoyal/gaelcge/References/kb_python_velo/mouse"

kb ref -i $pathtoindex/index.idx -g $pathtoindex/t2g.txt -f1 $pathtoindex/cdna.fa -f2 $pathtoindex/intron.fa -c1 $pathtoindex/spliced_t2c.txt -c2 $pathtoindex/unspliced_t2c.txt --workflow lamanno \
$PATHTOREFERENCE/Mus_musculus.GRCm38.dna.primary_assembly.fa \
$PATHTOREFERENCE/Mus_musculus.GRCm38.98.gtf















