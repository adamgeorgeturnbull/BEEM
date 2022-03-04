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

module load fsl/6.0.5.1
module load afni/21.1.07
fsldir=/software/fsl/6.0.5.1/bin
afnidir=/software/afni/21.1.07

fmPAdir=b0/003
fmAPdir=b0/004
fmPA=fieldmapPA
fmAP=fieldmapAP
taskdir=bold/006
task=fun_task

for dir in BEEM*
do
	echo reorienting $dir
	fslreorient2std ${dir}/${taskdir}/${task}.nii.gz ${dir}/${taskdir}/${task}_std.nii.gz
	fslreorient2std ${dir}/${fmPAdir}/${fmPA}.nii.gz ${dir}/${fmPAdir}/${fmPA}_std.nii.gz
	fslreorient2std ${dir}/${fmAPdir}/${fmAP}.nii.gz ${dir}/${fmAPdir}/${fmAP}_std.nii.gz

	echo generating reference $dir
	$afnidir/3dcalc -a ${dir}/${taskdir}/${task}_std.nii.gz[11] -expr 'a' -prefix ${dir}/${taskdir}/${task}_SBRef.nii.gz -overwrite

	echo motion correction $dir
	$fsldir/mcflirt -in ${dir}/${taskdir}/${task}_std.nii.gz -out ${dir}/${taskdir}/${task}_mcf.nii.gz -refvol ${dir}/${taskdir}/${task}_SBRef.nii.gz -plots -rmsrel

	$fsldir/flirt -in ${dir}/${fmPAdir}/${fmPA}_std.nii.gz -ref ${dir}/${taskdir}/${task}_SBRef.nii.gz -out ${dir}/${taskdir}/${fmPA}_mc.nii.gz -dof 6 -interp spline

	$fsldir/flirt -in ${dir}/${fmAPdir}/${fmAP}_std.nii.gz -ref ${dir}/${taskdir}/${task}_SBRef.nii.gz -out ${dir}/${taskdir}/${fmAP}_mc.nii.gz -dof 6 -interp spline

	echo topup $dir
	$fsldir/fslmerge -t  ${dir}/${taskdir}/fieldmap_taskREF.nii.gz ${dir}/${taskdir}/${fmAP}_mc.nii.gz ${dir}/${taskdir}/${fmPA}_mc.nii.gz

	$fsldir/topup --imain=${dir}/${taskdir}/fieldmap_taskREF.nii.gz  --datain=datain.txt --config=b02b0.cnf --out=${dir}/${taskdir}/my_output

	## topup defaults to spline
	$fsldir/applytopup --imain=${dir}/${taskdir}/${task}_mcf.nii.gz --inindex=1 --datain=datain.txt --topup=${dir}/${taskdir}/my_output --method=jac --out=${dir}/${taskdir}/${task}_tu.nii.gz

done

