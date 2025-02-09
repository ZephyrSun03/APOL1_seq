
bsub -P acc_zhangw09a -q premium -L /bin/bash -o %J.stdout -eo %J.stderr -n 1 -R rusage[mem=30000] -R span[ptile=1] -W 1440 -J soupX <soupX.sh

