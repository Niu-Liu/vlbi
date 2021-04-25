#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File Name: get_dir.py
"""
Created on Thu 15 Aug 2019 10:05:04 CST

@author: Neo (niu.liu@foxmail.com)
"""

import numpy as np
import os
from sys import platform as _platform

__all__ = ["get_home_dir", "get_aux_dir", "get_vlbi_sol_dir"]


# -------------------- FUNCTIONS --------------------
def get_home_dir():
    """Get the path of home directory.

    Return
    ------
    home_dir: string
        home directory
    """

    # Check the type of OS and get the home diretory path
    if _platform == "linux" or _platform == "linux2":
        # linux
        home_dir = os.getenv("HOME")
    elif _platform == "darwin":
        # MAC OS X
        home_dir = os.getenv("HOME")
    elif _platform == "win32" or _platform == "win64":
        # Windows
        home_dir = "C:/Users/Neo"
    else:
        print("Weird! What kind of OS do you use?")
        exit()

    return home_dir


def get_aux_dir():
    """Get the diretory of source name data

    Returns
    -------
    aux_dir: string
        diretory to put source name data

    """

    home_dir = get_home_dir()
    aux_dir = "{}/.ipython/my_progs/vlbi/aux_files".format(home_dir)

    return aux_dir


def get_vlbi_sol_dir():
    """Get the diretory of VLBI solutions

    Returns
    -------
    data_dir: string
        diretory to put source name data

    """

    # Check the type of OS and get the home diretory path
    if _platform == "linux" or _platform == "linux2":
        # linux
        sol_dir = "/data/vlbi_solutions"
    elif _platform == "darwin":
        # MAC OS X
        data_dir = os.getenv("HOME")
        sol_dir = "{}/Astronomy/data/vlbi".format(data_dir)
    elif _platform == "win32" or _platform == "win64":
        # Windows
        print("Not implemented yet")
        eixt()
    else:
        print("Weird! What kind of OS do you use?")
        exit()

    return sol_dir
# ----------------------- END -----------------------
