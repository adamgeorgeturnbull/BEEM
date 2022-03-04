#!/bin/bash

subjectList="BEEM_s504_T3 BEEM_s507_T3 BEEM_s510_T3 BEEM_s512_T3 BEEM_s513_T3 BEEM_s515_T3 BEEM_s518_T3 BEEM_s522_T3 BEEM_s524_T3 BEEM_s540_T3 BEEM_s546_T3 BEEM_s549_T3 BEEM_s552_T3 BEEM_s566_T3 BEEM_s570_T3 BEEM_s577_T3 BEEM_s582_T3 BEEM_s584_T3 BEEM_s589_T3 BEEM_s591_T3 BEEM_s592_T3 BEEM_s599_T3 BEEM_s601_T3 BEEM_s602_T3 BEEM_s607_T3 BEEM_s608_T3 BEEM_s613_T3 BEEM_s628_T3 BEEM_s629_T3 BEEM_s633_T3 BEEM_s635_T3 BEEM_s637_T3 BEEM_s638_T3 BEEM_s658_T3 BEEM_s659_T3 BEEM_s660_T3 BEEM_s664_T3 BEEM_s670_T3 BEEM_s672_T3 BEEM_s698_T3"

for subjectID in $subjectList
do

cd /scratch/tbaran2_lab/BEEM/T3/nii/${subjectID}
sed -e s:subjectID:"${subjectID}":g </scratch/tbaran2_lab/BEEM/RestingStateMultiband-main/preprocessing_singlesubj_3dTproject_T3.sh >preprocessing_${subjectID}.sh

#qsub  -cwd -N pre${subjectID} preprocessing_${subjectID}.sh

done
