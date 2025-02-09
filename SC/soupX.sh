
#ls -d /sc/arion/scratch/sunz04/PP2ASC/countAll.sh.318629.qsub/*/outs | awk '{split($0, tt, "/"); printf("%s\t%s\n", $1, tt[8])}' > Files

ml R/4.0.3

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/soupX.r Files

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/soupX.stat.r Files

