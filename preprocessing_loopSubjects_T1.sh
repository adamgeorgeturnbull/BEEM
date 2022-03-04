#!/bin/bash

subjectList="BEEM_s504_T1 BEEM_s507_T1 BEEM_s510_T1 BEEM_s512_T1 BEEM_s513_T1 BEEM_s515_T1 BEEM_s518_T1 BEEM_s522_T1 BEEM_s524_T1 BEEM_s540_T1 BEEM_s546_T1 BEEM_s549_T1 BEEM_s552_T1 BEEM_s566_T1 BEEM_s570_T1 BEEM_s577_T1 BEEM_s582_T1 BEEM_s584_T1 BEEM_s589_T1 BEEM_s591_T1 BEEM_s592_T1 BEEM_s599_T1 BEEM_s601_T1 BEEM_s602_T1 BEEM_s607_T1 BEEM_s608_T1 BEEM_s613_T1 BEEM_s628_T1 BEEM_s629_T1 BEEM_s633_T1 BEEM_s635_T1 BEEM_s637_T1 BEEM_s638_T1 BEEM_s658_T1 BEEM_s659_T1 BEEM_s660_T1 BEEM_s664_T1 BEEM_s670_T1 BEEM_s672_T1 BEEM_s698_T1"

for subjectID in $subjectList
do

cd /scratch/tbaran2_lab/BEEM/T1/nii/${subjectID}
sed -e s:subjectID:"${subjectID}":g </scratch/tbaran2_lab/BEEM/RestingStateMultiband-main/preprocessing_singlesubj_3dTproject.sh >preprocessing_${subjectID}.sh

#qsub  -cwd -N pre${subjectID} preprocessing_${subjectID}.sh

done
