
ls /sc/arion/scratch/sunz04/APOL1_clean/*star/*count.txt | awk '{split($0, tt, "/"); gsub("_star", "", tt[7]); printf("\t%s", tt[7])} END {printf("\n")}' > head

ls /sc/arion/scratch/sunz04/APOL1_clean/*star/*count.txt | head -1 | awk '{system("cut -f 1 "$1)}' > Expr.xls

for expr in /sc/arion/scratch/sunz04/APOL1_clean/*star/*count.txt
do

	cut -f 2 $expr > tt

	paste Expr.xls tt > ttt

	mv ttt Expr.xls

done

cp head Expr.xls.stat.xls

cat head Expr.xls | grep '__' >> Expr.xls.stat.xls

cat head Expr.xls | grep -v '__' > ttt

mv ttt Expr.xls


Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/RNA_seq/3.Analysis/idMap.r 

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/RNA_seq/3.Analysis/pca.r Expr.xls 1 "Batch;Condition" All Factor.xls

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/variancePar.r Expr.xls Factor.xls 1 Condition "Batch;Condition"


Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/limmaVoom.r Expr.xls Factor.xls unpaired Condition "G1;G0" 0 0.05 Ensembl G10 Batch voom 1 /sc/arion/projects/zhangw09a/PANDA/db_ZS/Ensembl/mm39/Mus_musculus.GRCm39.105.gff3.corr.gff3.map.ann 

awk -F "\t" '($5 <= 0.05){printf("%s\t%s\n", $8, $2)}' G10_diff.ann.xls | sort | uniq > AllGene

awk -F "\t" '($5 <= 0.05 && $2 > 0){printf("%s\t%s\n", $8, $2)}' G10_diff.ann.xls | sort | uniq > PosGene

awk -F "\t" '($5 <= 0.05 && $2 < 0){printf("%s\t%s\n", $8, $2)}' G10_diff.ann.xls | sort | uniq > NegGene

ml R/4.0.3

for gg in *Gene

do

        awk 'BEGIN{FS="\t"; OFS = "\t"} (NR==FNR){tt[$4]=$2; next} {if(tt[$1]!=""){printf("%s\t%s\n", tt[$1], $2)}}' /sc/arion/projects/zhangw09a/PANDA/db_ZS/MGI/HOM/homo_human_mouse.rpt $gg > tt && mv tt $gg

	Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/enrichR.r $gg

done

