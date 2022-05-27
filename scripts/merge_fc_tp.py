#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 12:00:09 2022

@author: aturnbu2
"""
seed = 'Postcentral_L'
import pandas as pd

T1 = pd.read_csv('/scratch/tbaran2_lab/BEEM/T1_%s_FC.csv' % seed)
T2 = pd.read_csv('/scratch/tbaran2_lab/BEEM/T2_%s_FC.csv' % seed)
T3 = pd.read_csv('/scratch/tbaran2_lab/BEEM/T3_%s_FC.csv' % seed)

dfs = [T1, T2, T3]
df_final = reduce(lambda left,right:pd.merge(left,right,on='subID'),dfs)

df_final.to_csv('/scratch/tbaran2_lab/BEEM/%s_FC.csv' % seed, index=False)