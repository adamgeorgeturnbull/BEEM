#!/bin/bash

subjectList="BEEM_s504 BEEM_s507 BEEM_s510 BEEM_s512 BEEM_s513 BEEM_s515 BEEM_s518 BEEM_s522 BEEM_s524 BEEM_s540 BEEM_s546 BEEM_s549 BEEM_s552 BEEM_s566 BEEM_s570 BEEM_s577 BEEM_s582 BEEM_s584 BEEM_s589 BEEM_s591 BEEM_s592 BEEM_s599 BEEM_s601 BEEM_s602 BEEM_s607 BEEM_s608 BEEM_s613 BEEM_s628 BEEM_s629 BEEM_s633 BEEM_s635 BEEM_s637 BEEM_s638 BEEM_s658 BEEM_s659 BEEM_s660 BEEM_s664 BEEM_s670 BEEM_s672 BEEM_s698"

timePoints="T1 T2 T3"

for subjectID in $subjectList

do

for timepoint in $timePoints

do
cd /scratch/tbaran2_lab/BEEM/${timepoint}/nii/${subjectID}_${timepoint}
sed -e s:subjectID:"${subjectID}":g -e s:timepoint:${timepoint}:g </scratch/tbaran2_lab/BEEM/scripts/generate_field_maps.sh >fieldmaps_${subjectID}_${timepoint}.sh

#qsub  -cwd -N pre${subjectID} preprocessing_${subjectID}.sh

done
done
