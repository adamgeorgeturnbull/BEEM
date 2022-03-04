#!/bin/bash

subjectList="504 507 510 512 513 515 518 522 524 540 546 549 552 566 570 577 582 584 589 591 592 599 601 602 607 608 613 628 629 633 635 637 638 658 659 660 664 670 672 698"

for subjectID in $subjectList
do

cd /scratch/tbaran2_lab/BEEM/T1/nii/BEEM_s${subjectID}_T1/bold
sed -e s:subjectID:"${subjectID}":g </scratch/tbaran2_lab/BEEM/T1/nii/first_level_template.fsf >${subjectID}_first_level.fsf

done
