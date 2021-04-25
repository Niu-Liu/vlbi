#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: comp_cat.py
"""
Created on Fri Mar 22 15:58:02 2019

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from astropy.table import Table, join, Column
import astropy.units as u
from astropy.coordinates import Angle
import numpy as np
from numpy import cos, sqrt
from functools import reduce
import time

__all__ = {"calc_cat_offset"}


# -----------------------------  FUNCTIONS -----------------------------
def root_sum_square(x, y):
    """Calculate the root-sum-square."""

    return np.sqrt(x**2 + y**2)


def nor_sep_calc(dRA, dRA_err, dDC, dDC_err, C):
    '''Calculate the normalized seperation.

    Parameters
    ----------
    dRA/dDC : Right Ascension / Declination differences in micro-as
    dRA_err/dDC_err : formal uncertainty of dRA*cos(Dec)/dDC in micro-as
    C : correlation coeffient between dRA*cos(Dec) and dDC.

    Returns
    ----------
    ang_sep : angular seperation, in micro-as
    X_a / X_d : normalized coordinate differences in RA / DC, unit-less
    X : Normalized separations, unit-less.
    '''

    # Angular seperations
    ang_sep = sqrt(dRA**2 + dDC**2)

    # Normalised coordinate difference
    X_a = dRA / dRA_err
    X_d = dDC / dDC_err

    # Normalised separation - Mignard's statistics (considering covariance)
    X = np.zeros_like(X_a)

    for i, (X_ai, X_di, Ci) in enumerate(zip(X_a, X_d, C)):
        if Ci == -1.:
            Ci = -0.99
        if Ci == 1.0:
            Ci = 0.99
        # print(Ci)

        wgt = np.linalg.inv(np.mat([[1, Ci],
                                    [Ci, 1]]))
        Xmat = np.mat([X_ai, X_di])
        X[i] = sqrt(reduce(np.dot, (Xmat, wgt, Xmat.T)))

    return ang_sep, X_a, X_d, X


def calc_cat_offset(t_cat1,  t_cat2, souname="iers_name"):
    """Calculate the radio source position differences between two catalogs.

    Parameters
    ----------
    t_cat1, t_cat2 : astropy.table object
        radio source positions from two catalogs

    Return
    ------
    t_cat_oft : astropy.table object
        radio source position differences
    """

    # Copy the original tables and keep only the source position information.
    t_cat3 = Table(t_cat1)
    t_cat3.keep_columns([souname, "ra",
                         "dec", "ra_err", "dec_err", "ra_dec_corr",
                         "used_sess", "used_obs"])

    t_cat4 = Table(t_cat2)
    t_cat4.keep_columns([souname, "ra", "dec",
                         "ra_err", "dec_err", "ra_dec_corr"])

    # Cross-match between two tables
    soucom = join(t_cat3, t_cat4, keys=souname)

    print("There are %d and %d sources in two catalogs, respectively,"
          "between which %d are common."
          % (len(t_cat1), len(t_cat2), len(soucom)))

    # Calculate the offset and the uncertainties
    arcfac = cos(Angle(soucom["dec_1"]).radian)
    dra = (soucom["ra_1"] - soucom["ra_2"]) * arcfac

    ddec = soucom["dec_1"] - soucom["dec_2"]
    dra_err = root_sum_square(soucom["ra_err_1"], soucom["ra_err_2"])
    ddec_err = root_sum_square(soucom["dec_err_1"], soucom["dec_err_2"])

    dra_ddec_cov = soucom["ra_err_1"] * soucom["dec_err_1"] * soucom["ra_dec_corr_1"] + \
        soucom["ra_err_2"] * soucom["dec_err_2"] * soucom["ra_dec_corr_2"]

    dra_ddec_corr = dra_ddec_cov / \
        root_sum_square(soucom["ra_err_1"], soucom["dec_err_1"]) / \
        root_sum_square(soucom["ra_err_2"], soucom["dec_err_2"])

    # Convert the unit
    dra.convert_unit_to(u.uas)
    dra_err.convert_unit_to(u.uas)
    ddec.convert_unit_to(u.uas)
    ddec_err.convert_unit_to(u.uas)
    dra_ddec_corr.unit = None
    dra_ddec_cov.unit = u.mas * u.mas
    dra_ddec_cov.convert_unit_to(u.uas * u.uas)

    # Calculate the angular seperation
    ang_sep, X_a, X_d, X = nor_sep_calc(
        dra, dra_err, ddec, ddec_err, dra_ddec_corr)

    # sou_nb = len(soucom)
    ang_sep = Column(ang_sep, unit=u.uas)
    X_a = Column(X_a)
    X_d = Column(X_d)
    X = Column(X)

    # Add these columns to the combined table.
    t_cat_oft = Table([soucom[souname],
                       soucom["ra_1"], soucom["dec_1"],
                       dra, dra_err, ddec, ddec_err,
                       dra_ddec_cov, dra_ddec_corr,
                       ang_sep, X, X_a, X_d],
                      names=["ivs_name", "ra", "dec",
                             "dra", "dra_err", "ddec", "ddec_err",
                             "dra_ddec_cov", "dra_ddec_corr",
                             "ang_sep", "nor_sep",
                             "nor_sep_ra", "nor_sep_dec"])

    return t_cat_oft

 # --------------------------------- END --------------------------------
if __name__ == '__main__':
    print("Nothing to do!")
    pass
# --------------------------------- END --------------------------------
