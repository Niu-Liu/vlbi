#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: sol_analysis.py
"""
Created on Mon Apr 19 13:21:11 2021

@author: Neo(niu.liu@nju.edu.cn)


Perform a quick analysis of the solution:
    (1) check EOP, offset wrt. C04, any linear drift, residual Nutation term
    (2) check CRF, offset wrt. ICRF3 
    (3) check TRF
"""

import numpy as np
