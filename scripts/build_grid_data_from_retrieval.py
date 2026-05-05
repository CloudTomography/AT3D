#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 17:15:51 2026

@author: benting
"""

import grid_data_builder

grid_data_builder.main([
    "--input-nc", "../data/retrieval_1d/2019_0806_1839_N_Pxl25_3_3.nc",
    "--output-csv", "../data/synthetic_cloud_fields/jpl_les/retrieval_2019_0806_1839_extended.csv",
    "--dx-km", "0.16",
    "--dy-km", "0.16",
    "--z-levels-km", "0.01:0.5:20",
    "--mode-count", "2"
])
