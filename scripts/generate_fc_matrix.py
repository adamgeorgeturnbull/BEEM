#!/usr/bin/env python
# coding: utf-8

from nilearn import image
import numpy as np
import os

pp='pp_9p_tf_sm'
path='/scratch/tbaran2_lab/BEEM/timepoint/nii/subjectID_timepoint'

os.chdir(path)

atlas_filename = image.load_img('/scratch/tbaran2_lab/BEEM/AAL3/AAL3v1.nii.gz')
fmri_filename = image.load_img('rest_preproc/rest_ep2d_bold_MB8_2mm_%s.nii.gz' % pp)
myfile = open("/scratch/tbaran2_lab/BEEM/AAL3/AAL3v1.nii.txt", 'r')
labels = [line.rstrip() for line in myfile]
myfile.close()

from nilearn.input_data import NiftiLabelsMasker

masker = NiftiLabelsMasker(labels_img=atlas_filename, labels=labels, standardize=True, memory='nilearn_cache', verbose=5)

time_series = masker.fit_transform(fmri_filename)

from nilearn.connectome import ConnectivityMeasure

correlation_measure = ConnectivityMeasure(kind='correlation')
correlation_matrix = correlation_measure.fit_transform([time_series])[0]

from nilearn import plotting

np.fill_diagonal(correlation_matrix,0)
fc_matrix = plotting.plot_matrix(correlation_matrix, figure=(10,8), labels=labels, vmax=0.8, vmin=-0.8, reorder=True)
fc_matrix.write_png('subjectID_timepoint_AAL3_fc_matrix_%s.png' % pp)
np.savetxt("subjectID_timepoint_AAL3_FC_matrix_%s.csv" % pp, correlation_matrix, delimiter=",")




