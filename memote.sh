#!/bin/sh
#SBATCH --job-name=memote
#SBATCH --time=02:00:00
#SBATCH --mem=20Gb
#SBATCH --cpus-per-task=1
#SBATCH --partition=workq
#SBATCH -e memote_ecoli_iter.%J.err
#SBATCH -o memote_ecoli_iter.%J.out


# memote report snapshot --filename "report_weighted.html" /ibex/user/niuk0a/funcarve/cobra/iAF987_gapfilled_weighted.xml
# memote report snapshot --filename "report_unweighted.html" /ibex/user/niuk0a/funcarve/cobra/iAF987_gapfilled.xml
# memote report snapshot --filename "report_ori.html" /ibex/user/niuk0a/funcarve/cobra/iAF987.xml

# /ibex/user/niuk0a/funcarve/cobra/iAF987_555_0.05_rm.xml
# memote report snapshot --filename "report_05_orirm.html" /ibex/user/niuk0a/funcarve/cobra/iAF987_555_0.05_rm.xml
# memote report snapshot --filename "report_05gfw.html" /ibex/user/niuk0a/funcarve/cobra/iAF987_0.05_gfW.xml
# memote report snapshot --filename "report_05gf.html" /ibex/user/niuk0a/funcarve/cobra/iAF987_0.05_gf.xml
# memote report snapshot --filename "report_v0_05gfw_filterU.html" /ibex/user/niuk0a/funcarve/cobra/iAF987_0.1_gfW_v0_filterU.xml
# memote report snapshot --filename "report_v0_05gfori_filterU.html" /ibex/user/niuk0a/funcarve/cobra/iAF987_0.1_gf_v0_filterU.xml
# for i in {1..9}; do
for i in {1,8};do
    echo 111111
    # memote report snapshot --filename "ecoli_inter_report_$i.html" /ibex/user/niuk0a/funcarve/cobra/NC_000913.3_ecoli_maxsepecoliINTERiter_$i.sbml
    memote run -a "--tb long" /ibex/user/niuk0a/funcarve/cobra/NC_000913.3_ecoli_maxsepecoliINTERiter_$i.sbml
done
