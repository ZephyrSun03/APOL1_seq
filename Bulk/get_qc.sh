
ls /sc/arion/projects/GOCAR/Sun/5.DRWork/6.APOL1/BulkNew/0.Data/*.fastq.gz | awk '{if(NR%2==1){split($1, tt, "/"); gsub("_R1_001.fastq.gz", "", tt[11]);printf("%s\t%s", tt[11], $1)}else{printf("\t%s\n", $1)}}' > Sam

mkdir -p QC

rm -f qc.sh

mkdir -p Clean

DIR=/sc/arion/projects/GOCAR/Sun/5.DRWork/6.APOL1/BulkNew/1.QC

awk -v DIR="$DIR" '{

	printf("module load fastqc/0.11.8");

	printf(" && fastqc -o %s/QC/ %s %s", DIR, $2, $3);

	printf(" && module load python/3.7.3 && cutadapt -m 50 -j 5 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o %s/Clean/%s.clean.R1.fq.gz -p %s/Clean/%s.clean.R2.fq.gz %s %s", DIR, $1, DIR, $1, $2, $3);

	printf(" && module load fastqc/0.11.8 && fastqc -o %s/QC/ %s/Clean/%s.clean.R1.fq.gz %s/Clean/%s.clean.R2.fq.gz", DIR, DIR, $1, DIR, $1);

	printf("\n");

}' Sam > qc.sh

nohup /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/common/common/qsub-sge.pl --queue premium --pro_code acc_zhangw09a --reqsub --jobprefix qc --resource 2000 --time 360 --convert no --verbose qc.sh &

