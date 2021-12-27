#  cis--trans-Regulation
Cis and trans regulatory evolution of cotton flowering time control
## Experimental Design
RNA-seq libraries were generated for LD and SD leaf from wild and domesticated Gossypium hirsutum accessions and their reciprocal F1 hybrids each with three replicates under long-day and short-day conditions:

**·** 2 condition: LD (7am & 9pm)& SD (7am &5 pm)
 
**·** 4 cotton accessions: Maxxa, TX2094, reciprocal F1 hybrids MxT and TxM

**·** 3-6 biological replicates
## RNA-seq mapping and high-quality SNP identification between parents
Filtered RNA-seq reads were mapped against *G. hirsutum* (TM-1) to estimate gene expressions in different accessions using Kallisato,  which psuedo-transcripts comes from hisat2.0 and strintie. In order to distinguish Maxxa(Ma) and TX2094(Tx), we used GATK(4.0.0.0) to detect allelic SNPs and selected high-quality SNPs on chromosomes based three criteria:

1) At least 3 of a single parent material are of non-NN type;
2) No more than one-third of the abnormal genotypes in a single parent material;
3) All heterozygous genotypes are filtered.

### 1.SNP calling
In order to obtain high-confidence allelic SNPs between parents, a total of 41 samples (including 16 fiber samples (Bao et al. 2019) and 25 leaves samples) were used.

```
 hisat2-build  Ghirsutum_527_v2.0.fa Ghir_tran
ls *fq.gz |while read line
 reads1=fasta/${line}_1.clean.fq.gz
 reads2=fasta/${line}_2.clean.fq.gz
 index=Ghir_tran
 genome=Ghirsutum_527_v2.0.fa

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
done 
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
### 2.Run stringtie and identify new transcripts
```
 while read line
do
 stringtie -p 8 -f 0.1 -at -j 10 -G Ghirsutum_527_v2.1.gene_exons.gtf -o ${line}.gtf -i ${line}.unique.bam
 samtools flagstat ${line}.unique.bam >${line}.txt
done < sample.txt
stringtie --merge -p 8 -f 0.1 -F -T -i -G Ghirsutum_527_v2.1.gene_exons.gtf -o stringtie_merged.gtf mergelist.txt # mergelist file

 while read line 
do
 stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/${line}/${line}.gtf ${line}.unique.bam
done < sample.txt
# We merged the gtf files of a single sample i the ballgown folder for subsequent identification of new transcripts.
 stringtie --merge -p 8 -f 0.1 -c 0.05 -F 1 -T -i -G Ghirsutum_527_v2.1.gene_exons.gtf -o merged.gtf mergelist.txt  # mergelist file

#Comparison of assembly transcripts and referwnce annotation transcripts. The identification of new transcripts must meet: 1)coding protein; 2)covering 70% of the coding  sequence of the reference gene.
 gffcompare -r Ghirsutum_527_v2.1.gene_exons.gff3 -G merged1.gtf
 awk -F '\t' '{if($4=="=") print $0}' gffcmp.tracking >reference_transcript.txt  #107,002 transcripts 
 awk -F '\t' '{if($4=="c"|| $4=="k"|| $4=="m"|| $4=="n"|| $4=="j" ) print $0}' gffcmp.tracking >new_transcript.txt  #41,559 new transcripts
 cut -f 3 new_transcript.txt >new_transcript_id.txt
 awk -F '\t' 'NR==FNR{a[$0]}NR>FNR{if($1 in a)print $0 }' new_transcript_id.txt merged.gtf >new_merged.gtf
 gffread new_merged.gtf -o new_transcript.gff3 
 gffread -w transcripts.fa -g Ghirsutum_527_v2.0.fa new_transcript.gff3
```
### 3.Psuedo-genome preparing and new transcript screening
Based on the identified high-quality SNPs information on chromosome, we create a pseudo-genome by replacing reference genome bases to extract transcript sequence. 
```
 python changeGenomeVcf.py pseudo-genome.fa Tx2094.txt Tx2094
 python changeGenomeVcf.py pseudo-genome.fa Tx2094.txt Tx2094

#Identification of protein coding ability using CPC2 software (version 2.0 http://cpc2.gao-lab.org/download.php). 
 python CPC2.py -i new.transcripts.fa -o output.txt    #33,996 transcripts coding
```
### 4. Run Kallisato mapping
We used 140,998 transcripts (including 107,002 reference transcripts and 33,996 new trancsripts) to make the mapping database.
```
 gffread -w transcripts.fa -g $genome new_transcript1.gff3
 sed 's/>/>Ma_/g' Maxxa.fa >Maxxa.tran.fa
 sed 's/>/>Tx_/g' Tx2094.fa >Tx2094.tran.fa
 cat Maxxa.tran.fa Tx2094.tran.fa >MT.tran.fa
#build index
 kallisto index -i MT.idx MT.tran.fa
 ls *clean.fa.gz| while read line
do
 index=MT.idx
 gtf=new_transcript1.gff3
 Chr=chrom.txt
 output=${line}
 mkdir ${line}
 kallisto quant $reads1 $reads2 -i $index -o $output --bias -t 4 -b 100 --pseudobam --genomebam -g $gtf -c $Chr
  
done
```
