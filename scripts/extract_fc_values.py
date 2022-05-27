#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 09:39:56 2022

@author: aturnbu2
"""
seed = 'Postcentral_L'
timepoint = 'T3'
dataloc = '/scratch/tbaran2_lab/BEEM/fc_results/*%s*' % timepoint

from glob import glob
import pandas as pd
import numpy as np
import re

def remove(list):
    pattern1 = '[0-9]'
    list = [re.sub(pattern1, '', i) for i in list]
    pattern2 = ' '
    list = [re.sub(pattern2, '', i) for i in list]
    return list

atlas = '/scratch/tbaran2_lab/BEEM/AAL3/AAL3v1.nii.gz'
label_file = '/scratch/tbaran2_lab/BEEM/AAL3/AAL3v1.nii.txt'

myfile = open(label_file, 'r')
labels = [line.rstrip() for line in myfile]
myfile.close()

labels = remove(labels)

cog_net = ['Frontal_Sup_Medial_L',
           'Frontal_Sup_Medial_R', 'ACC_sup_L',
           'ACC_sup_R', 'Frontal_Inf_Tri_L']
cog_net_i = [labels.index(item) for item in cog_net]

emo_net = ['OFCmed_L', 'OFCmed_R', 'Insula_L',
           'Insula_R', 'Amygdala_L', 'Amygdala_R']
emo_net_i = [labels.index(item) for item in emo_net]

NPS_net = ['Paracentral_Lobule_L', 'Frontal_Mid__L',
           'ACC_sub_L', 'ACC_pre_L', 'ACC_sup_L',
           'Caudate_L', 'Pallidum_L', 'Fusiform_L',
           'Frontal_Sup_Medial_R', 'Caudate_R',
           'Rectus_R', 'Frontal_Med_Orb_R']
NPS_net_i = [labels.index(item) for item in NPS_net]

files = sorted(glob(dataloc))

columns = cog_net + emo_net
columns.append('C3_NPS')
columns.append('within_NPS')
columns = [s + '_' + timepoint for s in columns]
columns.insert(0,'subID')
df = pd.DataFrame(columns=columns)

seed_i = labels.index(seed)

for file in files:
    data = np.arctanh(np.genfromtxt(file, delimiter=','))
    seed_con = data[seed_i,:]
    cog_vals = list(seed_con[cog_net_i])
    emo_vals = list(seed_con[emo_net_i])
    c3_nps = np.mean(seed_con[NPS_net_i])
    tmp = data[NPS_net_i,:]
    tmp2 = tmp[:,NPS_net_i]
    w_nps = np.mean(tmp2[np.triu_indices_from(tmp2, k=1)])
    subID = file[42:46]
    vals = cog_vals + emo_vals
    vals.append(c3_nps)
    vals.append(w_nps)
    vals.insert(0, subID)
    df.loc[len(df)] = vals
    
df.to_csv('/scratch/tbaran2_lab/BEEM/%s_%s_FC.csv' % (timepoint, seed), index=False)
    
