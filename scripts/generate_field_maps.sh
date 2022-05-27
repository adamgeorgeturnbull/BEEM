#!/bin/bash

#SBATCH -J beem_fieldmaps
#SBATCH -c 8
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=24G
#SBATCH -t 2-
#SBATCH -o fieldmaps.log
#SBATCH -e fieldmaps.err
#SBATCH --mail-type=END
#SBATCH --mail-user=adam_turnbull@urmc.rochester.edu

subject=subjectID
tp=timepoint

module load fsl/6.0.5.1
module load afni/21.1.07
fsldir=/software/fsl/6.0.5.1/bin
afnidir=/software/afni/21.1.07

fmPA=ep2d_bold_2mm_PA
fmAP=ep2d_bold_2mm_AP
task=task_ep2d_bold_MB8_2mm

cwd=/scratch/tbaran2_lab/BEEM/${tp}/nii/${subject}_${tp}

cd $cwd

mkdir task_preproc

echo reorienting $subject $tp
fslreorient2std ${task}.nii.gz task_preproc/${task}_std.nii.gz
fslreorient2std ${fmPA}.nii.gz task_preproc/${fmPA}_std.nii.gz
fslreorient2std ${fmAP}.nii.gz task_preproc/${fmAP}_std.nii.gz

cd task_preproc

echo generating reference $subject $tp
$afnidir/3dcalc -a ${task}_std.nii.gz[11] -expr 'a' -prefix ${task}_SBRef.nii.gz -overwrite

echo motion correction $subject $tp
$fsldir/mcflirt -in ${task}_std.nii.gz -out ${task}_mcf.nii.gz -refvol ${task}_SBRef.nii.gz -plots -rmsrel

$fsldir/flirt -in ${fmPA}_std.nii.gz -ref ${task}_SBRef.nii.gz -out ${fmPA}_mc.nii.gz -dof 6 -interp spline

$fsldir/flirt -in ${fmAP}_std.nii.gz -ref ${task}_SBRef.nii.gz -out ${fmAP}_mc.nii.gz -dof 6 -interp spline

echo topup $subject $tp
$fsldir/fslmerge -t  fieldmap_taskREF.nii.gz ${fmAP}_mc.nii.gz ${fmPA}_mc.nii.gz

$fsldir/topup --imain=fieldmap_taskREF.nii.gz  --datain=../../datain.txt --config=b02b0.cnf --out=my_output

	## topup defaults to spline
$fsldir/applytopup --imain=${task}_mcf.nii.gz --inindex=1 --datain=../../datain.txt --topup=my_output --method=jac --out=${task}_tu.nii.gz
