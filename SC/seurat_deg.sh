
ml R/4.0.3

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/seurat4_integrate_deg.r ../cluster.lst.ann_Obj_post.rds ../Sample.info.HTO "G1;WT" "Treg_1;CD4_naive;CD8_naive;CD8_early_act;Prolif;TE;TEM;Tfh;DN"


for dd in *DEG.xls
do

        awk '($2<=0.01){printf("%s\t%s\n", $1, $3)}' $dd > ${dd}.all.lst

done


for ll in *50.lst
do

        awk 'BEGIN{FS="\t"; OFS = "\t"} (NR==FNR){tt[$4]=$2; next} {if(tt[$1]!=""){printf("%s\t%s\n", tt[$1], $2)}}' /sc/arion/projects/zhangw09a/PANDA/db_ZS/MGI/HOM/homo_human_mouse.rpt $ll > tt && mv tt $ll

        Num=$(wc -l $ll | awk '{print $1}')

        if [ "$Num" -gt 20 ]
        then

                Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/enrichR.r $ll

        fi

done

