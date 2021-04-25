#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: souname_xmatch.py
"""
Created on Thu Apr 26 15:55:00 2018

@author: Neo(liuniu@smail.nju.edu.cn)

Read radio source name from data list. Obsolete!
"""

from astropy.io import fits
from astropy import table
from astropy.table import Table
import numpy as np
import os
import time

# My modules
from .get_dir import get_aux_dir


__all__ = ["get_souname"]


# -----------------------------  FUNCTIONS -----------------------------
def get_souname():
    """Read source names.

    Parameters
    ----------
    None

    Returns
    -------
    souname: table object
        different designations of radio source names
    """

    datadir = get_aux_dir()
    datafil = "{}/source.names".format(datadir)

    # empty array to store data
    souname = Table.read(datafil, format="ascii.fixed_width_no_header",
                         names=["ivs_name", "icrf_name", "iers_name", "class"],
                         col_starts=(0, 10, 22, 42),
                         col_ends=(8, 20, 30, 42))

    # Eliminate duplicate sources
    souname = table.unique(souname, keys="ivs_name")

    # Fill the empty filed of IERS name by the IVS name
    for i in souname["iers_name"].mask.nonzero()[0]:
        souname[i]["iers_name"] = souname[i]["ivs_name"]

    return souname


if __name__ == "__main__":
    get_souname()
# --------------------------------- END --------------------------------
