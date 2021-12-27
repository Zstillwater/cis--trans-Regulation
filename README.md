#  cis--trans-Regulation
Cis and trans regulatory evolution of cotton flowering time control
## Experimental Design
RNA-seq libraries were generated for LD and SD leaf from wild and domesticated Gossypium hirsutum accessions and their reciprocal F1 hybrids each with three replicates under long-day and short-day conditions:

**·** 2 condition: LD (7am & 9pm)& SD (7am &5 pm)
 
**·** 4 cotton accessions: Maxxa, TX2094, reciprocal F1 hybrids MxT and TxM

**·** 3-6 biological replicates
## RNA-seq mapping and high-quality SNP identification between parents
Filtered RNA-seq reads were mapped against *G. hirsutum* (TM-1) to estimate gene expressions in different accessions using Hisat2.0 and Kallisato. In order to distinguish Maxxa(Ma) and TX2094(Tx), we used GATK(4.0.0.0) to detect allelic SNPs and selected high-quality SNPs on chromosomes based three criteria:

1) At least 3 of a single parent material are of non-NN type;
2) No more than one-third of the abnormal genotypes in a single parent material;
3) All heterozygous genotypes are filtered.

### 1.SNP calling
In order to obtain high-confidence allelic SNPs between parents, a total of 41 samples (including 16 fiber samples (Bao et al. 2019) and 25 leaves samples) were used.

```
 while read line
do
 reads1=fasta/${line}_1.clean.fq.gz
 reads2=fasta/${line}_2.clean.fq.gz
 index=~/UTX/Ghir_tran
 genome=~/UTX/Ghirsutum_527_v2.0.fa

 hisat2 -p 8 --dta --min-intronlen 50 --max-intronlen 50000 -x $index -1 $reads1 -2 $reads2 -S ${line}.sam
 samtools sort -@ 8 ${line}.sam | samtools view -bhF 12 -q 30 > ${line}.unique.bam
# Remove duplicates
 sambamba markdup -r ${line}.unique.bam --overflow-list-size 600000 --tmpdir='./' ${line}.marked.bam
# Splicing site recognition
 gatk SplitNCigarReads -R $genome -I ${line}.marked.bam  -O ${line}.splited.bam
# Add sample information label 
 gatk AddOrReplaceReadGroups -LB line -PL ILLUMINA -PU line -SM line -I ${line}.splited.bam -O ${line}.splited.added.bam
# Build index
 samtools index ${line}.splited.added.bam
# SNP identification
 gatk HaplotypeCaller -R $genome -I ${line}.splited.added.bam -dont-use-soft-clipped-bases -stand-call-conf 20.0 -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation -O ${line}.g.vcf
 rm ${line}.sam
 rm ${line}.bam
done < sample.txt
# Merge and convert 
 gatk CombineGVCFs -R $genome --variant gvcf.list -O all.g.vcf
 gatk GenotypeGVCFs -R $genome -V all.g.vcf -O all.vcf
# Filter multiple alleles
 vcftools --vcf all.vcf --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out gatk --min-meanDP 3 --maf 0.05 --max-missing 0.5
 # Extract SNP information
 gatk SelectVariants --select-type-to-include SNP -R $genome -V gatk.recode.vcf -O all.snp.vcf
 # Base quality filtering
 gatk VariantFiltration -R $genome -V all.snp.vcf -window 35 -cluster 3 \-filter-name FS -filter "FS > 30.0" -filter-name QD -filter "QD < 2.0" -O all.snp.filtered.vcf
 # Filter through PASS and biased allelic SNP  
 bcftools view --threads 4 -m1 -M2 -f PASS all.snp.filtered.vcf |sed 's/|/\//g' >all.snp.filtered.pass.vcf
 # Generate snp.list file, about 100,254 SNPs on chromosome
 python ExtractVCF.py SNP all.snp.filtered.pass.vcf 
 # Generate high-quality allelic SNPs 
 python pop-specificSNPs.py -g group.txt -Q1 Q1 -Q2 Q2 -hom snp.list 

```




