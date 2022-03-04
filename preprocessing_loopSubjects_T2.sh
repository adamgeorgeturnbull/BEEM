#!/bin/bash

subjectList="BEEM_s504_T2 BEEM_s507_T2 BEEM_s510_T2 BEEM_s512_T2 BEEM_s513_T2 BEEM_s515_T2 BEEM_s518_T2 BEEM_s522_T2 BEEM_s524_T2 BEEM_s540_T2 BEEM_s546_T2 BEEM_s549_T2 BEEM_s552_T2 BEEM_s566_T2 BEEM_s570_T2 BEEM_s577_T2 BEEM_s582_T2 BEEM_s584_T2 BEEM_s589_T2 BEEM_s591_T2 BEEM_s592_T2 BEEM_s599_T2 BEEM_s601_T2 BEEM_s602_T2 BEEM_s607_T2 BEEM_s608_T2 BEEM_s613_T2 BEEM_s628_T2 BEEM_s629_T2 BEEM_s633_T2 BEEM_s635_T2 BEEM_s637_T2 BEEM_s638_T2 BEEM_s658_T2 BEEM_s659_T2 BEEM_s660_T2 BEEM_s664_T2 BEEM_s670_T2 BEEM_s672_T2 BEEM_s698_T2"

for subjectID in $subjectList
do

cd /scratch/tbaran2_lab/BEEM/T2/nii/${subjectID}
sed -e s:subjectID:"${subjectID}":g </scratch/tbaran2_lab/BEEM/RestingStateMultiband-main/preprocessing_singlesubj_3dTproject_T2.sh >preprocessing_${subjectID}.sh

#qsub  -cwd -N pre${subjectID} preprocessing_${subjectID}.sh

done
