# BEEM

preprocessing_singlesubj_3dTproject.sh

Primary script used for resting state preprocessing. datain.txt and slicetiming.txt required in parent directory. 

preprocessing_loopSubjects.sh

Generates a version of preprocessing_singlesubj_3dTproject.sh for each participant. 

generate_field_maps.sh

Generates field maps and performs distortion correction for task-based fMRI analysis. Rest of task-based analysis is done using FSL (.fsf design file provided).

fieldmaps_loopSubjects.sh

Generates a verrsion of generate_field_maps.sh for each participant. 

generate_fc_matrix.py

Generates fc matrix from data preprocessed using preprocessing_singlesubj_3dTproject.sh. Requires a virtual environment set up for nilearn. 

extract_fc_values.py

Generates a version of generate_fc_matrix.py for each participant.

extract_fc_values.py

Gets specific functional connectivity values from fc matrices generated by generate_fc_matrix.py

merge_fc_tp.py

Merges time points of functional connectivity values extracted using extract_fc_values.py





