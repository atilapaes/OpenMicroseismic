#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 17:42:39 2017

@author: atilapaes
"""

import pandas


catalog_file_PE='output_toc2me_3c_d305.csv'
catalog_PE=pandas.read_csv(catalog_file_PE)

catalog_PE['SNR']=catalog_PE.peak_height/catalog_PE.signal_mean

high_SNR=catalog_PE[catalog_PE['SNR']>10]

high_SNR=high_SNR.sort(['SNR'], ascending=False)


high_SNR.to_csv('high_SNR.csv', sep=',', line_terminator='\n', index=False)