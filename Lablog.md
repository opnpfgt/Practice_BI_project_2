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

In order to catch all our reads we will use max_length of reads.
coverage = (read count * read length ) / total genome size
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
IGV often switchs this rows...

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

8.1. SNP detection
Now let's take a look at rare variants

```bash
varscan mpileup2snp roommate.mpileup --min-var-freq 0.001 --variants --output-vcf 1 > VarScan_roommate_rare.vcf

Only SNPs will be reported
Warning: No p-value threshold provided, so p-values will not be calculated
Min coverage:   8
Min reads2:     2
Min var freq:   0.001
Min avg qual:   15
P-value thresh: 0.01
Reading input from roommate.mpileup
1665 bases in pileup file
20 variant positions (18 SNP, 2 indel)
0 were failed by the strand-filter
18 variant positions reported (18 SNP, 0 indel)
```

```bash
cat VarScan_roommate_rare.vcf | awk 'NR>24 {print $1, $2, $4, $5, $10}'

KF848938.1 72 A G 1/1:255:16832:16794:6:16787:99.96%:0E0:35:36:4:2:10898:5889
KF848938.1 117 C T 1/1:255:20768:20663:36:20625:99.82%:0E0:35:37:27:9:13462:7163
KF848938.1 254 A G 0/1:25:31444:31311:31248:58:0.19%:2.7559E-3:36:36:21205:10043:36:22
KF848938.1 307 C T 0/1:255:29088:28985:28700:280:0.97%:5.4371E-54:36:35:17875:10825:148:132
KF848938.1 340 T C 0/1:21:29906:29758:29702:52:0.17%:6.9675E-3:37:35:18978:10724:34:18
KF848938.1 389 T C 0/1:37:27391:27231:27167:61:0.22%:1.8501E-4:37:36:14167:13000:41:20
KF848938.1 722 A G 0/1:41:30510:30484:30411:68:0.22%:7.715E-5:37:36:19958:10453:39:29
KF848938.1 744 A G 0/1:24:30912:30870:30805:56:0.18%:3.3182E-3:37:32:19818:10987:33:23
KF848938.1 774 T C 1/1:255:31145:31006:5:30997:99.97%:0E0:34:37:5:0:18170:12827
KF848938.1 802 A G 0/1:50:32182:32089:32008:77:0.24%:9.5375E-6:37:35:18435:13573:28:49
KF848938.1 915 T C 0/1:32:31553:31455:31384:63:0.2%:6.2548E-4:35:34:17111:14273:36:27
KF848938.1 999 C T 1/1:255:29708:29277:36:29236:99.86%:0E0:35:35:21:15:16126:13110
KF848938.1 1043 A G 0/1:26:28946:28903:28845:55:0.19%:2.0066E-3:35:33:16344:12501:18:37
KF848938.1 1086 A G 0/1:31:23993:23992:23939:51:0.21%:7.5186E-4:36:35:12540:11399:21:30
KF848938.1 1213 A G 0/1:34:25177:25093:25035:56:0.22%:3.7244E-4:37:36:8964:16071:24:32
KF848938.1 1260 A C 1/1:255:23067:23033:2:23019:99.94%:0E0:32:37:0:2:9824:13195
KF848938.1 1280 T C 0/1:20:23487:23462:23418:43:0.18%:9.2875E-3:37:35:11147:12271:24:19
KF848938.1 1458 T C 0/1:255:26333:26200:25983:214:0.82%:2.1459E-38:37:35:6834:19149:80:134
```

n=18 + 2 indels

9. Biological reps
Since we can't say whether our SNP is real or just a mistake of PCR, sequencing we need make sure that this is real sheep not the wolves. 
We can perform it using additional reference sequencing in plasmid (biologocal reps).

9.1. Download the files

```bash
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/008/SRR1705858/SRR1705858.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/009/SRR1705859/SRR1705859.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/000/SRR1705860/SRR1705860.fastq.gz
```

9.2. Check statistics

```bash
zcat SRR1705858.fastq.gz | seqkit stats
file  format  type  num_seqs     sum_len  min_len  avg_len  max_len
-     FASTQ   DNA    256,586  38,118,730       35    148.6      151

seqkit stats SRR1705859.fastq.gz
file                 format  type  num_seqs     sum_len  min_len  avg_len  max_len
SRR1705859.fastq.gz  FASTQ   DNA    233,327  34,636,567       35    148.4      151

seqkit stats SRR1705860.fastq.gz
file                 format  type  num_seqs     sum_len  min_len  avg_len  max_len
SRR1705860.fastq.gz  FASTQ   DNA    249,964  37,170,486       35    148.7      151
```

9.3. Evaluate coverage again, it's going to be the same as we have the same max length.


9.4. Align new seqs
```bash
bwa mem reference.fasta SRR1705858.fastq.gz | samtools view -S -b - | samtools sort - -o SRR1705858.bam

bwa mem reference.fasta SRR1705859.fastq.gz | samtools view -S -b - | samtools sort - -o SRR1705859.bam

bwa mem reference.fasta SRR1705860.fastq.gz | samtools view -S -b - | samtools sort - -o SRR1705860.bam
```

9.5. Index these files
```bash
samtools index SRR1705858.bam
samtools index SRR1705859.bam
samtools index SRR1705860.bam
```

9.6. Create mpileups

```bash
samtools mpileup -d 23500 -f reference.fasta index SRR1705858.bam > c58.mpileup
samtools mpileup -d 23500 -f reference.fasta index SRR1705859.bam > c59.mpileup
samtools mpileup -d 23500 -f reference.fasta index SRR17058560.bam > c60.mpileup
```

9.7. Check how many SNPs we obtained

```bash
varscan mpileup2snp c58.mpileup --min-var-freq 0.001 --variants --output-vcf 
1 > VarScan_58.vcf

1665 bases in pileup file
55 variant positions (55 SNP, 0 indel)
1 were failed by the strand-filter
54 variant positions reported (54 SNP, 0 indel)


varscan mpileup2snp c58.mpileup --min-var-freq 0.001 --variants --output-vcf 1 > VarScan_58.vcf

1665 bases in pileup file
50 variant positions (50 SNP, 0 indel)
2 were failed by the strand-filter
48 variant positions reported (48 SNP, 0 indel)


varscan mpileup2snp c58.mpileup --min-var-freq 0.001 --variants --output-vcf 1 > VarScan_58.vcf

1665 bases in pileup file
57 variant positions (57 SNP, 0 indel)
0 were failed by the strand-filter
57 variant positions reported (57 SNP, 0 indel)
```

We obtained the following picture:
58 – 54 snps
59 – 49 snps
60 – 57 snps

P.S. if we calculate it using Average reads length, we obtain 54, 49 and 57 SNP positions respectively.

9.8. Extract frequences
This thing is masked from the naked eye and we need to look for the percent subcolumn inside the 10th column.

```bash
cat VarScan_c58_rare.vcf | awk 'NR>24 {print $2, $4, $5}' #position, ref nucleotide, alt nucleotide
cat VarScan_c58_rare.vcf | awk 'NR>23, FS=":" {print $20}' #frequency
```

Next we can copy this columns to a txt file and import to Excel or R or etc or just ctrl C + ctrl V it and add to the google sheets. It will suggest you to separate columns and we need to delete percent symbols.

9.9. Count frequncy and std for each control column(it's gonna be 3 samples)
After that we get the average of these measurements and calculate average + 3std as lower bound for frequency.

10. Compare frequencies of SNPs in SOSEDs genome with ave+3std.
The ones that were higher this bound, can be estimated as real snips.

|  ID  |Ref|Alt| Freq  |
|------|---|---|-------|
|  307 | C | T | 0.95% |
| 1458 | T | C | 0.83% |

Answers:
CCG-TCG
TAT-TAC
