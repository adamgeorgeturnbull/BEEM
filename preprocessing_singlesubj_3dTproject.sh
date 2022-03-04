#!/bin/bash 

#SBATCH -J preprocessing_subjectID
#SBATCH -c 8
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=24G
#SBATCH -t 2-
#SBATCH -o rest_preproc.log
#SBATCH -e rest_preproc.err
#SBATCH --mail-type=END
#SBATCH --mail-user=adam_turnbull@urmc.rochester.edu

## Controls which modules are being run:
# MODULE 1:
# motion correction, topup, and slice timing:
run_motion=1

# MODULE 2:
# anatomical registration: FNIRT to MNI
run_anatomical=1

# MODULE 3:
# Functional registration to MNI
run_functional_registration=1

# MODULE 4:
# Nuisance regression
run_nuisance_regression=1

## Modify for each subject:
subject=subjectID
FWHM=6
hp=0.009
lp=0.08

## Runs four piplines: temporal filtering (yes / no) by spatial smoothing (yes / no)

cwd=/scratch/tbaran2_lab/BEEM/T1/nii/${subject}

standard=/software/fsl/6.0.5.1/data/standard/MNI152_T1_2mm_brain.nii.gz
standard_nomask=/software/fsl/6.0.5.1/data/standard/MNI152_T1_2mm.nii.gz
standard_mask=/software/fsl/6.0.5.1/data/standard/MNI152_T1_2mm_brain_mask_dil.nii.gz
# NOTE: Need to edit the paths in the T1_2_MNI152_2mm.cnf as well

prior_dir=/scratch/tbaran2_lab/BEEM/RestingStateMultiband-main/tissuepriors
# fsl priors thresholded >0.5

anatdir=t1/002
anat=t1

restdir=rest/005
rest=fun_rest

fmPAdir=b0/003
fmAPdir=b0/004
fmPA=fieldmapPA
fmAP=fieldmapAP

# afnidir=/Users/ben/abin
afnidir=/software/afni/21.1.07

# fsldir=/Users/brisk/Applications2/fsl/bin/
fsldir=/software/fsl/6.0.5.1/bin

n_volslist=(253)
TRlist=(1.01)

cd $cwd

module load afni/21.1.07
module load fsl/6.0.5.1

mkdir preprocessed
## ######################
## 1. MOTION CORRECTION, TOP-UP, AND SLICE TIMING:
## ######################

## switch all orientations to standard

fslreorient2std ${anatdir}/${anat}.nii.gz ${anatdir}/${anat}_std.nii.gz
fslreorient2std ${restdir}/${rest}.nii.gz ${restdir}/${rest}_std.nii.gz
fslreorient2std ${fmPAdir}/${fmPA}.nii.gz ${fmPAdir}/${fmPA}_std.nii.gz
fslreorient2std ${fmAPdir}/${fmAP}.nii.gz ${fmAPdir}/${fmAP}_std.nii.gz

if [ $run_motion -eq 1 ]
then
	echo ------------------------------
	echo ---- RUNNING FUNCTIONAL PREPROCESSING ----
	echo ------------------------------

	## drop first 4 TRs

	n_vols=`echo ${n_volslist[0]}-1 | bc`

	echo "Dropping first TRs"
	$afnidir/3dcalc -a ${restdir}/${rest}_std.nii.gz[4..${n_vols}] -expr 'a' -prefix ./preprocessed/${rest}_dr.nii.gz -overwrite


	## Extract twelfth volume; this will serve as the "SBRef" for registration

	$afnidir/3dcalc -a ${restdir}/${rest}_std.nii.gz[11] -expr 'a' -prefix ${restdir}/${rest}_SBRef.nii.gz -overwrite

	## Motion correction
	## mcflirt: use -plots options to retain 6 parameters
	
	# Defaults to 6 dof, i.e., rigid body
	# use -rmsrel to obtain framewise displacement.
	$fsldir/mcflirt -in ./preprocessed/${rest}_dr.nii.gz -out ./preprocessed/${rest}_mcf.nii.gz -refvol ${restdir}/${rest}_SBRef.nii.gz -plots -rmsrel

	## Estimate distortion field
	## Topup correction
	## Note: Alternatively, we could do a single motion correction to the SBREF3 image.
	## However, this would result in inflated estimates of subject motion using conventional measures
	echo -----------------------------------------------------
	echo -------- Applying TOPUP
	echo -----------------------------------------------------


	## motion correct PA to the AP reference:
	$fsldir/flirt -in ${fmPAdir}/${fmPA}_std.nii.gz -ref ${restdir}/${rest}_SBRef.nii.gz -out ${fmPAdir}/${fmPA}_mc.nii.gz -dof 6 -interp spline

	$fsldir/flirt -in ${fmAPdir}/${fmAP}_std.nii.gz -ref ${restdir}/${rest}_SBRef.nii.gz -out ${fmAPdir}/${fmAP}_mc.nii.gz -dof 6 -interp spline

	## fslmerge: concatenate the AP and PA single band reference images
	$fsldir/fslmerge -t  b0/fieldmap_restREF.nii.gz ${fmAPdir}/${fmAP}_mc.nii.gz ${fmPAdir}/${fmPA}_mc.nii.gz

	$fsldir/topup --imain=b0/fieldmap_restREF.nii.gz  --datain=../datain.txt --config=b02b0.cnf --out=my_output

	## topup defaults to spline
	$fsldir/applytopup --imain=./preprocessed/${rest}_mcf.nii.gz --inindex=1 --datain=../datain.txt --topup=my_output --method=jac --out=./preprocessed/${rest}_tu.nii.gz

	## Slice timing correction
	## NOTE: --repetition flag does not do anything when specifying custom slice timing acquisition
	echo Slice time correction ${restdir}/${rest}.nii.gz 

	$fsldir/slicetimer -i ./preprocessed/${rest}_tu.nii.gz -o ./preprocessed/${rest}_st.nii.gz --tcustom=../slicetiming.txt

fi

## ######################
## 2. Structural registration
## ######################

if [ $run_anatomical -eq 1 ]
then
	echo ------------------------------
	echo ---- RUNNING ANATOMICAL ----
	echo ------------------------------

	## Skull-strip anatomical:
	cd $cwd

	$afnidir/3dSkullStrip -input ${anatdir}/${anat}_std.nii.gz -o_ply ${anatdir}/${anat}_surf.nii.gz -overwrite
	$afnidir/3dcalc -a ${anatdir}/${anat}_std.nii.gz -b ${anatdir}/${anat}_surf.nii.gz -expr 'a*step(b)' -prefix ${anatdir}/${anat}_std_brain.nii.gz -overwrite
	
fi

## ######################
## 3. Functional Registration
## ######################
if [ $run_functional_registration -eq 1 ]
then
	echo ------------------------------
	echo !!!! RUNNING REGISTRATION !!!!
	echo ------------------------------


	## Register FUNCTION to T1. 
	reg_dir=${cwd}/preprocessed/reg
	mkdir ${reg_dir}
	cd $cwd

	## create example_func for compatability with pipeline:
	## NOTE: data have been aligned to MB3_SBRef
	## Thus, use MB3_SBRef for all 4D datasets:
	## cp ${rest}_SBRef.nii.gz ./${rest}/example_func.nii.gz
	cp ${restdir}/${rest}_SBRef.nii.gz ${reg_dir}/example_func.nii.gz


	## 1. Copy required images into reg directory
	### copy anatomical
	cp ${anatdir}/${anat}_std_brain.nii.gz ${reg_dir}/highres.nii.gz
	cp ${anatdir}/${anat}_std.nii.gz ${reg_dir}/highres_head.nii.gz
	### copy standard
	cp ${standard} ${reg_dir}/standard.nii.gz
	cp ${standard_nomask} ${reg_dir}/standard_head.nii.gz
	cp ${standard_mask} ${reg_dir}/standard_mask.nii.gz
	#cp ${cwd}/highres2standard.mat ${reg_dir}/highres2standard.mat



	## 2. cd into reg directory
	cd ${reg_dir}

	## 3. register functional to anatomical, anatomical to standard, and functional to standard

	$fsldir/epi_reg --epi=example_func --t1=highres_head --t1brain=highres --out=example_func2highres

	$fsldir/convert_xfm -inverse -omat highres2example_func.mat example_func2highres.mat

	$fsldir/slicer example_func2highres highres -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png ; $fsldir/pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png example_func2highres1.png ; $fsldir/slicer highres example_func2highres -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png ; $fsldir/pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png example_func2highres2.png ; $fsldir/pngappend example_func2highres1.png - example_func2highres2.png example_func2highres.png; /bin/rm -f sl?.png example_func2highres2.png

	$fsldir/flirt -in highres -ref standard -out highres2standard -omat highres2standard.mat -cost corratio -dof 12 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp trilinear

	$fsldir/fnirt --iout=highres2standard_head --in=highres_head --aff=highres2standard.mat --cout=highres2standard_warp --iout=highres2standard --jout=highres2highres_jac --config=T1_2_MNI152_2mm --ref=standard_head --refmask=standard_mask --warpres=10,10,10

	$fsldir/applywarp -i highres -r standard -o highres2standard -w highres2standard_warp

	$fsldir/convert_xfm -inverse -omat standard2highres.mat highres2standard.mat

	$fsldir/slicer example_func2highres highres -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png ; $fsldir/pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png example_func2highres1.png ; $fsldir/slicer highres example_func2highres -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png ; $fsldir/pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png example_func2highres2.png ; $fsldir/pngappend example_func2highres1.png - example_func2highres2.png example_func2highres.png; /bin/rm -f sl?.png example_func2highres2.png

	/bin/rm example_func2highres1.png

	$fsldir/flirt -in highres -ref standard -out highres2standard -omat highres2standard.mat -cost corratio -dof 12 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp trilinear 

	$fsldir/fnirt --iout=highres2standard_head --in=highres_head --aff=highres2standard.mat --cout=highres2standard_warp --iout=highres2standard --jout=highres2highres_jac --config=T1_2_MNI152_2mm --ref=standard_head --refmask=standard_mask --warpres=10,10,10

	$fsldir/applywarp -i highres -r standard -o highres2standard -w highres2standard_warp

	$fsldir/convert_xfm -inverse -omat standard2highres.mat highres2standard.mat

	$fsldir/slicer highres2standard standard -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png ; $fsldir/pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png highres2standard1.png ; $fsldir/slicer standard highres2standard -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png ; $fsldir/pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png highres2standard2.png ; $fsldir/pngappend highres2standard1.png - highres2standard2.png highres2standard.png; /bin/rm -f sl?.png highres2standard2.png

	/bin/rm highres2standard1.png

	$fsldir/convert_xfm -omat example_func2standard.mat -concat highres2standard.mat example_func2highres.mat

	$fsldir/convertwarp --ref=standard --premat=example_func2highres.mat --warp1=highres2standard_warp --out=example_func2standard_warp

	$fsldir/applywarp --ref=standard --in=example_func --out=example_func2standard --warp=example_func2standard_warp

	$fsldir/convert_xfm -inverse -omat standard2example_func.mat example_func2standard.mat

	$fsldir/slicer example_func2standard standard -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png ; $fsldir/pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png example_func2standard1.png ; $fsldir/slicer standard example_func2standard -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png ; $fsldir/pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png example_func2standard2.png ; $fsldir/pngappend example_func2standard1.png - example_func2standard2.png example_func2standard.png; /bin/rm -f sl?.png example_func2standard2.png

	$fsldir/applywarp --ref=standard --in=../${rest}_st.nii.gz --out=../${rest}_res2standard_fnirt --warp=example_func2standard_warp
	$afnidir/3dTstat -mean -prefix ../${rest}_res2standard_fnirt_mean.nii.gz ../${rest}_res2standard_fnirt.nii.gz  -overwrite

	cd $cwd

	############################
	## Basic preprocessing: 
	## Skull strip
	############################

	## Remove skull/edge detect
	echo "Skull stripping ${subject} ${rest}"
	$afnidir/3dAutomask -prefix ./preprocessed/${rest}_mask.nii.gz -dilate 1 ./preprocessed/${rest}_res2standard_fnirt.nii.gz -overwrite
	$afnidir/3dcalc -a ./preprocessed/${rest}_res2standard_fnirt.nii.gz -b ./preprocessed/${rest}_mask.nii.gz -expr 'a*b' -prefix ./preprocessed/${rest}_ss.nii.gz -overwrite

	## Grandmean scaling
	echo "Grand-mean scaling ${subject} ${rest}"
	$fsldir/fslmaths ./preprocessed/${rest}_ss.nii.gz -ing 10000 ./preprocessed/${rest}_pp.nii.gz -odt float

	# Generate binary mask, this is slightly bigger than 3dAutomask for an MB 12 that I checked
	echo "Generating mask of preprocessed data for ${rest}"
	$fsldir/fslmaths ./preprocessed/${rest}_pp.nii.gz -Tmin -bin ./preprocessed/${rest}_pp_mask.nii.gz -odt char

	## BRisk: Clean-up
	## rm ${rest}_st.nii.gz ${rest}_ss.nii.gz ${rest}_sm.nii.gz ${rest}_gms.nii.gz ${rest}_filt.nii.gz ${rest}_filt_mean.nii.gz ${rest}_dt.nii.gz ${rest}_mcf_tomb3.nii.gz ${rest}_mcf.nii.gz ${rest}_mask.nii.gz
fi

##########################################################################################################################
##---NUISANCE SIGNAL REGRESSION ----------------------------------------------------------------------------------------------------##
##########################################################################################################################

## BRisk: The F1000 pipelines uses segmentation and then intersects with FSL templates.
## The fsl templates have narrowly defined wm and csf using the mask (_bin) >=0.51. The Segmentation in fast captures a 
## lot of csf between the gray matter and skull, and is a little tricky to use due to partial volume effects. 
## Since the templates are conservative, and since we are using FNIRT, I do the averaging in 
## MNI space aligned volumes using the FSL (harvard) tissue priors

if [ $run_nuisance_regression -eq 1 ]
then


	echo --------------------------------------------
	echo !!!! RUNNING NUISANCE SIGNAL REGRESSION !!!!
	echo --------------------------------------------

	segment_dir=${cwd}/preprocessed/segment
	reg_dir=${cwd}/preprocessed/reg

	cd $cwd
	# First, intersect tissue priors with global mask from the MB3 acquisition
	echo ${prior_dir}

	PRIOR_WHITE=${prior_dir}/avg152T1_white_bin.nii.gz
	PRIOR_CSF=${prior_dir}/avg152T1_csf_bin.nii.gz

	mkdir ${segment_dir}

	## 4. Copy functional mask from FSLpreproc step 5 - this is the global signal mask
	echo "Creating global mask"
	cp ./preprocessed/${rest}_pp_mask.nii.gz ${segment_dir}/global_mask.nii.gz

	## Mask CSF template by subject's global mask and copy to segment directory:
	$fsldir/fslmaths ${PRIOR_CSF} -mas ${segment_dir}/global_mask ${segment_dir}/csf_mask

	## White matter mask
	$fsldir/fslmaths ${PRIOR_WHITE} -mas ${segment_dir}/global_mask ${segment_dir}/wm_mask

	nuisance_dir=./preprocessed/nuisance

	## 1. make nuisance directory
	mkdir -p ${nuisance_dir}

	# Extract signal for global, csf, and wm
	## 2. Global
	echo "Extracting global signal for ${subject}"
	$afnidir/3dmaskave -mask ${segment_dir}/global_mask.nii.gz -quiet ./preprocessed/${rest}_pp.nii.gz > ${nuisance_dir}/global.1D -overwrite

	## 3. csf
	echo "Extracting signal from csf for ${subject}"
	$afnidir/3dmaskave -mask ${segment_dir}/csf_mask.nii.gz -quiet ./preprocessed/${rest}_pp.nii.gz > ${nuisance_dir}/csf.1D -overwrite

	## 4. wm
	echo "Extracting signal from white matter for ${subject}"
	$afnidir/3dmaskave -mask ${segment_dir}/wm_mask.nii.gz -quiet ./preprocessed/${rest}_pp.nii.gz > ${nuisance_dir}/wm.1D -overwrite

	## 5. Perform linear and quadratic detrending, nuisance regression, temporal filtering
	echo "Performing nuisance regression, spatial smoothing, and temporal filtering"

## No spatial smoothing:
	$afnidir/3dTproject -input ./preprocessed/${rest}_pp.nii.gz -prefix ./preprocessed/${rest}_pp_9p.nii.gz -polort 2 -ort ${nuisance_dir}/global.1D ${nuisance_dir}/wm.1D ${nuisance_dir}/csf.1D ./preprocessed/${rest}_mcf.nii.gz.par -automask -overwrite -verb
	$afnidir/3dTproject -input ./preprocessed/${rest}_pp.nii.gz -prefix ./preprocessed/${rest}_pp_9p_tf.nii.gz -polort 2 -ort ${nuisance_dir}/global.1D ${nuisance_dir}/wm.1D ${nuisance_dir}/csf.1D ./preprocessed/${rest}_mcf.nii.gz.par -passband $hp $lp -automask -overwrite -verb

## sdev maps
	$afnidir/3dTstat -stdev -prefix ./preprocessed/${rest}_pp_stdev.nii.gz ./preprocessed/${rest}_pp.nii.gz  -overwrite

	$afnidir/3dcalc -a ./preprocessed/${rest}_pp_9p_stdev.nii.gz -b  $base_data -expr 'b/a' -prefix ./preprocessed/${rest}_pp_9p_gfactor.nii.gz -overwrite
	$afnidir/3dcalc -a ./preprocessed/${rest}_pp_9p_tf_stdev.nii.gz -b  $base_data_tf -expr 'b/a' -prefix ./preprocessed/${rest}_pp_9p_tf_gfactor.nii.gz -overwrite


## Spatial smoothing:
	$afnidir/3dTproject -input ./preprocessed/${rest}_pp.nii.gz -prefix ./preprocessed/${rest}_pp_9p_sm.nii.gz -polort 2 -ort ${nuisance_dir}/global.1D ${nuisance_dir}/wm.1D ${nuisance_dir}/csf.1D ./preprocessed/${rest}_mcf.nii.gz.par -blur $FWHM -automask -overwrite -verb
	$afnidir/3dTproject -input ./preprocessed/${rest}_pp.nii.gz -prefix ./preprocessed/${rest}_pp_9p_tf_sm.nii.gz -polort 2 -ort ${nuisance_dir}/global.1D ${nuisance_dir}/wm.1D ${nuisance_dir}/csf.1D ./preprocessed/${rest}_mcf.nii.gz.par -passband $hp $lp -blur $FWHM -automask -overwrite -verb

## sdev maps
	$afnidir/3dTstat -stdev -prefix ./preprocessed/${rest}_pp_9p_sm_stdev.nii.gz ./preprocessed/${rest}_pp_9p_sm.nii.gz  -overwrite
	$afnidir/3dTstat -stdev -prefix ./preprocessed/${rest}_pp_9p_tf_sm_stdev.nii.gz ./preprocessed/${rest}_pp_9p_tf_sm.nii.gz -overwrite

fi
