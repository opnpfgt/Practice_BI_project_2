# Lablog

In 2014 the vaccine includes:
- an A/California/7/2009 (H1N1) pandemic 2009-like strain (similar to swine flu)
- an A/Texas/50/2012 (H3N2)-like strain
- a B/Massachusetts/2/2012-like strain.

_A/Hong Kong/4801/2014 (H3N2) was found in our neighbour_

0. Install all the required packages through conda within my env



 1. Download the roommate's sequence
 ```bash
    wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/001/SRR1705851/SRR1705851.fastq.gz
    gunzip SRR1705851.fastq.gz
    mv SRR1705851.fastq roommate.fastq
 ```

 2. Download reference genome manually

 Send to -> File -> FASTA format

 3. Check reads content with seqkit and count reads
 
 ```bash
seqkit stats roommate.fastq

file            format  type  num_seqs     sum_len  min_len  avg_len  max_len
roommate.fastq  FASTQ   DNA    358,265  52,717,864       35    147.1      151

seqkit stats reference.fasta

file             format  type  num_seqs  sum_len  min_len  avg_len  max_len
reference.fasta  FASTA   DNA          1    1,665    1,665    1,665    1,665

cat roommate.fastq | wc -l # получилось 1433060 строк, 1433060/4 = 358265 ридов
 ```
 3.1. Check reads quality 

 ```bash
 fastqc -o. roommate.fastq
 ```

 That's fine.

 4. Reference indexing 
 ```bash
 bwa index reference.fasta
 ```

 5. Allignment, zipping, sorting, indexing of sorted file
 ```bash
 bwa mem reference.fasta roommate.fastq > alignment.sam
 samtools view -S -b alignment.sam > alignment.bam
 samtools sort alignment.bam -o roommate.bam
 samtools index roommate.bam
 ```
 ```bash
samtools flagstatroommate.bam

361349 + 0 in total (QC-passed reads + QC-failed reads)
358265 + 0 primary
0 + 0 secondary
3084 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
361116 + 0 mapped (99.94% : N/A)
358032 + 0 primary mapped (99.93% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

6. Since mpileup uses treshold = 8000 max reads per base, it's not going to be enough for our alignment.


Calculations of coverage

In order to cathc all our reads we will use max_length of reads.

mpileup_treshold = 358_265 * 147 / 1665 = 32_492 ~ 32500

7. Piling up

```bash
samtools mpileup -d 32500 -f reference.fasta roommate.bam > roommate.mpileup
```

8. SNP detection

```bash
varscan mpileup2snp roommate.mpileup --min-var-freq 0.95 --variants --output-vcf 1 > VarScan_roommate.vcf

#Общий синтез awk – это awk '/search_pattern/ {actiontotakeonmatches; another_action;}' file_to_awk. NR означает количество прочитанных строк, поэтому в команде написано «только если номер строки >24, делать то, что дальше». Awk использует знак доллара для обозначения полей (столбцов), поэтому дальше «напечатать столбцы 1, 2, 4 и 5».

cat VarScan_roommate.vcf | awk 'NR>24 {print $1, $2, $4, $5}'

KF848938.1 72 A G
KF848938.1 117 C T
KF848938.1 774 T C
KF848938.1 999 C T
KF848938.1 1260 A C
```

Found SNPs do not make sense in terms of protein sythesis.
Top reading frame was used as compatable with aminoacids in protein data from NCBI.
IGV often swiths this rows...

ACA – ACG

![SNPs visualisation with IGV](https://github.com/opnpfgt/Practice_BI_project_2/blob/main/Practice_BI_2_frequent_snps/A_to_G_snp_pos74.png)


GCC-GCT

![SNPs visualisation with IGV](https://github.com/opnpfgt/Practice_BI_project_2/blob/main/Practice_BI_2_frequent_snps/C_to_T_snp_pos117.png)


TTT-TTC

![SNPs visualisation with IGV](https://github.com/opnpfgt/Practice_BI_project_2/blob/main/Practice_BI_2_frequent_snps/TtoC_snp_pos_774.png)


GGC-GGT

![SNPs visualisation with IGV](https://github.com/opnpfgt/Practice_BI_project_2/blob/main/Practice_BI_2_frequent_snps/C_to_T_snp_pos_999.png)


CTA-CTC

![SNPs visualisation with IGV](https://github.com/opnpfgt/Practice_BI_project_2/blob/main/Practice_BI_2_frequent_snps/A_to_C_snp_pos1260.png)
