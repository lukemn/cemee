#!/bin/bash

usage()
{
    echo -e "\n\Sample GWAS null distribution using MultiTrans ( http://genetics.cs.ucla.edu/multiTrans/)."
    echo -e "\tInstallation directory is set to ~/bin/MultiTrans"
    echo -e "\tWriting to `pwd`/mult"
    echo -e "\tRequires plink for LD filtering."
    echo -e "\tRuns parallel across traits, requires >=96Gb RAM for 10 traits, 100k SNPs."
    echo -e "\n\tusage: `basename $0` <PHENO> <BFILE> <K> <PREFIX>\n"
    echo -e "\n\tPHENO\tN(+1) lines X D+1 traits (header ignored, line ID in first column)"
    echo -e "\tFILE\tpath_to_plink/binary_file_stem"
    echo -e "\tK\tN X N(+1) GRM (header ignored)"
    echo -e "\tPREFIX\toutput prefix"    
    exit 0;
}

##########################################################################################
if [[ $# -eq 0 || ! $# -eq 4 ]]; then echo "check args"; usage; exit 1; fi
PHENO=$1
BFILE=$2
K=$3
PREF=$4
MT=~/bin/MultiTrans/
MTP=~/bin/MultiTrans/Pylmm_MultiTrans
WD=`pwd`
OUTD=$WD/mult; if [ ! -s $OUTD ]; then mkdir $OUTD; fi
for i in $PHENO $SNP $K $MT; do 
    if [ ! -s $WD/$i ]; then echo -e "can't find $i in $WD\n"; exit 1; fi
done
for i in bim bed fam; do 
    if [ ! -s $BFILE.$i ]; then echo -e "can't find $BFILE.$i\n"; exit 1; fi
done
# load plink (NYU prince cluster)
module load plink/1.90b4.4
# LD pruning threshold
R2=0.99
# MultiTrans windowSize (number of markers within which local LD is considered)
W=1000
# and number of null samples to take
NSAMP=200000
##########################################################################################

# data munging
PHENO=$WD/$PHENO
K=$WD/$K
# remove GRM/pheno header if present and make space delim if not already
NR=$(cat $K | wc -l)
NC=$(head -1 $K | wc -w)
if [[ $NR -gt $NC ]]; then tail -n+2 $K | tr -s "\t" " " > $OUTD/$K; K=$OUTD/$K; fi
NR=$(cut -f1 -d" " $K | wc -w)
if [[ $NR -gt $NC ]]; then cat $K | tr -s "\t" " " > tmp; mv tmp $OUTD/$K; K=$OUTD/$K; fi
N=$(cat $K | wc -l)
NR=$(cat $PHE | wc -l)
if [[ $NR -gt $N ]]; then tail -n+2 $PHE | tr -s "\t" " " > $OUTD/pheno; PHE=$OUTD/pheno; fi
NR=$(cut -f1 -d" " $PHE | wc -w)
if [[ $NR -gt $N ]]; then cat $PHE | tr -s "\t" " " > tmp; mv tmp $OUTD/pheno; PHE=$OUTD/pheno; fi

# subset BFILE using lines in PHENO and run 2-pass LD pruning
cut -f1 $PHE > lines
grep -w -f lines $BFILE.fam | cut -f1,2 -d" " | sort -k2 > keep.id
plink --file $BFILE --keep keep.id --make-bed --indiv-sort f keep.id --out $PREF
plink --bfile $PREF --out $PREF --recode
plink --geno 0.01 --indep-pairwise 1000 500 $R2 --file $PREF
plink --file $PREF --extract plink.prune.in --keep keep.id --make-bed --indiv-sort f keep.id --out $PREF
plink --bfile $PREF --out $PREF --recode
# 2nd pass
plink --indep-pairwise 1000 500 $R --file $PREF
plink --file $PREF --extract plink.prune.in --keep keep.id --make-bed --indiv-sort f keep.id --out $PREF
plink --bfile $PREF --out $PREF --recode
rm -f plink* *nosex $PREF.log
# split by chromosome
for i in {1..6}; do 
    plink --bfile $PREF --chr $i --out $PREF.chr$i --make-bed
    plink --bfile $PREF.chr$i --out $PREF.chr$i --recode
done
# extract genotypes
for i in {1..6}; do pedExtract.R $PREF.chr$i; done

# 1. Estimate sigma_g^2, sigma_e^2 for all traits
$MTP/pylmmGWAS_multiPhHeri.py -v --bfile $PREF --emmaPHENO $PHE --kfile $K --REML $OUTD/varcomps

# 2. Estimate correlation in the rotated space, and generate correlation band matrix c
# split by chromosome
function multByChrom () {
    local PHEIX=$1
    local O=$OUTD/phe$PHEIX
    for CHROM in {1..6}; do 
        MX=$O/chr$CHROM.maxstat
        rm -f $MX
        if [ ! -s $MX ]; then 
            if [ ! -s $O/chr$CHROM ]; then mkdir -p $O/chr$CHROM; fi 
            head -n $PHEIX $OUTD/varcomps | tail -n1 > $O/$PHEIX.var
            R CMD BATCH --args -Xpath=$PREF.chr$CHROM.txt \
                -Kpath=$K -VCpath=$O/$PHEIX.var \
                -outputPath=$O/chr$CHROM $MT/generateR.R generateR.$PHEIX.$CHROM.log
            java -jar $MT/generateC/generateC.jar $W $O/chr$CHROM/r.txt $O/chr$CHROM/c.txt
            slide_1prep -C $O/chr$CHROM/c.txt $W $O/chr$CHROM/prep
            slide_2run $O/chr$CHROM/prep $O/chr$CHROM.maxstat $NSAMP 123
            rm -f $O/chr$CHROM/r.txt $O/chr$CHROM/c.txt
        fi
    done
    
    # run GWAS (approx. PPPD - variance components estimated once up front)    
    $MTP/pylmmGWAS.py --bfile $BFILE --emmaPHENO $PHE -p $((PHEIX-1)) \
        --kfile $K --REML $O/pvals
    # correct stats, and output adjusted thresholds    
    slide_3sort $O/sorted $O/*maxstat
    cut -f5 $O/pvals | tail -n+2 > tmp
    slide_4correct -p $O/sorted tmp $O/MultiTrans.output
    slide_4correct -t $O/sorted $OUTD/thresholds threshout_${PHEIX}    
}

# run for each trait in the background
echo -e "0.1\n0.05\n0.01" > $OUTD/thresholds
TIX=$(( $(head -1 $PHE | wc -w) -1 ))
for j in $(seq 1 $TIX); do multByChrom $j & done
wait $(jobs -p)    

# print mean threshold, effective number of tests
for i in $(cat $OUTD/thresholds); do 
    echo $i
    echo P
    awk -v P=$i '($1 == P)' $OUTD/threshout_* | grep -v "\--" | cut -f3 | sump.py | grep -A4 mean
    echo M_eff
    awk -v P=$i '($1 == P)' $OUTD/threshout_* | grep -v "\--" | cut -f5 | sump.py | grep -A4 mean
    echo
done    
    

