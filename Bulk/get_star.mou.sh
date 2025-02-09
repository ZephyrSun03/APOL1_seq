

rm -f run_star.sh

DIREC=/sc/arion/scratch/sunz04/APOL1

mkdir -p $DIREC

awk -v DIREC="$DIREC" '{

	printf("module load python/3.7.3");

	printf(" && rm -rf %s/%s_star && mkdir %s/%s_star && module load star/2.7.5b && STAR --runThreadN 5 --genomeDir /sc/arion/projects/zhangw09a/PANDA/db_ZS/Ensembl/mm39/mm39_star_149/ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within KeepPairs --readFilesIn %s %s --readFilesCommand zcat --outFileNamePrefix %s/%s_star/ ", DIREC, $1, DIREC, $1, $2, $3, DIREC, $1);

	printf(" && /sc/arion/projects/zhangw09a/PANDA/ext/samtools/samtools index %s/%s_star/Aligned.sortedByCoord.out.bam", DIREC, $1);

	printf(" && htseq-count -f bam -r pos -s no -m union --nonunique none %s/%s_star/Aligned.sortedByCoord.out.bam /sc/arion/projects/zhangw09a/PANDA/db_ZS/Ensembl/mm39/Mus_musculus.GRCm39.105.gff3.corr.gff3 > %s/%s_star/%s_count.txt", DIREC, $1, DIREC, $1, $1);

        printf(" && htseq-count -f bam -r pos -s no -m union --nonunique all %s/%s_star/Aligned.sortedByCoord.out.bam /sc/arion/projects/zhangw09a/PANDA/db_ZS/Ensembl/mm39/Mus_musculus.GRCm39.105.gff3.corr.gff3 > %s/%s_star/%s_count.multi.txt", DIREC, $1, DIREC, $1, $1);

	printf("\n");

}' Sam > run_star.sh

nohup /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/common/common/qsub-sge.pl --queue premium --convert no --pro_code acc_zhangw09a --jobprefix Star --resource 40000 --time 1440 --verbose run_star.sh &

