
set -e

module load R/4.0.3

ls -d /sc/arion/projects/GOCAR/Sun/5.DRWork/6.APOL1/SCNew/0.Data/Ambient/*.rds | awk '{split($1, tt, "/"); gsub("_cleaned.rds", "", tt[12]); printf("%s\t%s\tRDS\n", $1, tt[12])}' > List

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/seurat4_integrate_Raw.r List RDS /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/Markers.pbmc.mou 0.8 30 30 3 200 800000 30 All 2000 &> seurat.log

sh /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/seurat4_stat.sh

awk 'BEGIN{FS="\t"; OFS="\t"} {print $0, $2}' cluster.lst > cluster.lst.ori

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/seurat4_integrate_post.r Obj.rds No cluster.lst.ori Sample.info Null

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/seurat4_integrate_post.r Obj.rds Yes cluster.lst.ori Sample.info Null


Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/seurat4_checkMarkers.r cluster.lst.ori_Obj_post.rds Markers.imm cluster.lst.ori /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/Markers.pbmc.mou RNA Null No

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/seurat4_checkMarkers.r cluster.lst.ori_Obj_post.rds Markers.pbmc cluster.lst.ori Markers.pbmc RNA Null No

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/seurat4_checkMarkers.r cluster.lst.ori_Obj_post.rds Markers.paperNC cluster.lst.ori Markers.paperNC RNA Null No


Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/doubletFinder.r cluster.lst.ori_Obj_post.rds


Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/seurat_hash.r "/sc/arion/projects/GOCAR/Sun/5.DRWork/6.APOL1/SCNew/0.Data/SN0309466/KW12284_jchoi44/cellranger-8.0.0/outs/filtered_feature_bc_matrix" cluster.lst.ori_Obj_post.rds 0.9 "Singlet;Doublet;Negative"

awk 'BEGIN{FS = "\t"; OFS = "\t";} (NR==FNR){if($3!="RM"){tt[$1] = $3} next;} {if(tt[$2] != ""){print $0, tt[$2]}}' anno.lst cluster.lst > cluster.lst.ann

awk '(NR==FNR){if($14 == "Doublet"){tt[$1] = 1} next}{if(tt[$1] == ""){print $0}}' DBTab.xls cluster.lst.ann > tt && mv tt cluster.lst.ann

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/seurat4_integrate_post.r cluster.lst.ori_Obj_post.rds.hto.rds.hto.rds No cluster.lst.ann Sample.info.HTO Null

