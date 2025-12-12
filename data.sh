# copy data from pipeline over
p1=/home/projects/rth/co2capture/continuation-ses/subprojects/OCyRS/development/OCyRS-pipeline_20250626-u5uvOz/data
p2=/home/projects/rth/co2capture/continuation-ses/subprojects/OCyRS/production/cyanobacteria-structures_20250505-Ft4MIO
p3=/home/projects/rth/co2capture/subprojects/OCyRS/OCyRS-pipeline/data
p4=/home/projects/rth/co2capture/subprojects/OCyRS/OCyRS-companion/Public-RNAseq/

mkdir -p data
# pathways and taxonomy associations
cp $p1/K_ko-path.tsv data/
cat $p1/K2_motifs.tsv | \
 awk '{sub(".fna.motif",""); print $0}' > data/K2_motifs.tsv
cat $p1/K_motif-path.tsv | \
 awk '{sub(".fna.motif",""); print $0}' > data/K_motif-path.tsv
cat $p1/K_motif-tax-pos.tsv | \
 awk '{sub(".fna.motif",""); print $0}' > data/K_motif-tax-pos.tsv
cat $p1/K_motif-tax.tsv | \
 awk '{sub(".fna.motif",""); print $0}' > data/K_motif-tax.tsv
# scores
cat $p1/J_novel/potentially-novel-motifs.tsv | \
 awk '{sub(".fna.motif",""); print $0}' > data/potentially-novel-motifs.tsv
cat $p1/I_fdr.tsv | \
 awk '{sub(".fna.motif",""); print $0}' > data/I_fdr.tsv
# conservation OGs
cp $p3/B_OGs.tsv data/

# CMsearch results
cp $p3/G_rfam-cmsearch.tsv.gz data/
cat $p1/J_novel/references_inside-of_intergenic_regions.tsv.gz | \
 awk '{sub(".fna.motif",""); print $0}' > data/references_inside-of_intergenic_regions.tsv.gz
cp $p3/A_representatives/taxonomy.tsv data/
cp $p3/G2_terminators.tsv.gz data/

# Results from the expression analysis
cat $p4/5-expression-ratios.tsv | \
 awk '{sub(".fna.motif",""); print $0}' > data/5-expression-ratios.tsv
cat $p4/5-maybe-interest.tsv | \
 awk '{sub(".fna.motif",""); print $0}' > data/5-maybe-interest.tsv

# iterate over motifs
tail -n +2 $p1/K_motif-tax.tsv | cut -f 1 | sort | uniq | \
 awk '{motif=$1;sub(".fna.motif","",motif); print $1"\t"motif}' > data/motifs.txt
# to load R2R figures
mkdir -p data/R2R

# get R2R images
while IFS=$'\t' read -r motif crs; do
  echo "${motif}"
  cp $p2/results/postalign-"$motif"-petfold+cmbuild+cmalign-canbp65-gapcol95-lbp0/align_qc/3/rscape-v2_0_0_j_twoset/fig/"$motif".realign-canonicalbp65_1.R2R.sto.svg data/R2R/$crs.svg
  # Remove the annoying data_1 label
  grep -v data_1 data/R2R/$crs.svg  > data/R2R/$crs.svg.new
  mv data/R2R/$crs.svg.new data/R2R/$crs.svg
done < "data/motifs.txt"

# get post-processed alignment
mkdir -p data/motifs
while IFS=$'\t' read -r motif crs; do
  echo "${motif}"
  cp $p2/results/postalign-"$motif"-petfold+cmbuild+cmalign-canbp65-gapcol95-lbp0/"$motif".realign-canonicalbp65.sto data/motifs/$crs.sto
done < "data/motifs.txt"
# build JalView images
mkdir -p data/jalview data/motifs-fasta
jalview=$(which jalview)
for i in data/motifs/*.sto ; do
    x=$(basename $i .sto)
    if [ ! -f "data/jalview/$x.svg" ] ; then
        $jalview --open="data/motifs/$x.sto" --colour=nucleotide --image="data/jalview/$x.svg" \
            --output="data/motifs-fasta/$x.fasta" \
            --props="my_jalview.properties"
    fi
done
echo "Check if this is the expected number of motifs:"
ls -1 data/jalview | wc -l

# get filtered alignment
mkdir -p data/motifs-gapfree
while IFS=$'\t' read -r motif crs; do
  echo "${motif}"
  cp $p2/results/postalign-"$motif"-petfold+cmbuild+cmalign-canbp65-gapcol95-lbp0/"$motif".realign-canonicalbp65-gapcolumn95-sorted.sto data/motifs-gapfree/$crs.sto
done < "data/motifs.txt"
# build JalView images
mkdir -p data/jalview-gapfree data/motifs-fasta-gapfree
jalview=$(which jalview)
for i in data/motifs-gapfree/*.sto ; do
    x=$(basename $i .sto)
    if [ ! -f "data/jalview-gapfree/$x.svg" ] ; then
        $jalview --open="data/motifs-gapfree/$x.sto" --colour=nucleotide --image="data/jalview-gapfree/$x.svg" \
            --output="data/motifs-fasta-gapfree/$x.fasta" \
            --props="my_jalview.properties"
    fi
done
echo "Check if this is the expected number of motifs:"
ls -1 data/jalview-gapfree | wc -l

# get covariance models
mkdir -p data/cm
while IFS=$'\t' read -r motif crs; do
  cp $p2/results/postalign-"$motif"-petfold+cmbuild+cmalign-canbp65-gapcol95-lbp0/"$motif".realign.cm data/cm/$crs.cm
done < "data/motifs.txt"

# get KO phylogenetic trees
mkdir -p data/ko_phylo
find /home/projects/rth/co2capture/subprojects/OCyRS/OCyRS-pipeline/data/C_phylo/ -type d -name "K*" | \
 xargs -I {} basename {} | \
 xargs -I {} cp /home/projects/rth/co2capture/subprojects/OCyRS/OCyRS-pipeline/data/C_phylo/{}/bootstrap-consensus.tree data/ko_phylo/{}.bootstrap-consensus.tree

# get supplementary files
cd data
tar cvfz crs-alignment.tar.gz motifs/*.sto
tar cvfz crs-covariance_models.tar.gz cm/*.cm
cd ..
