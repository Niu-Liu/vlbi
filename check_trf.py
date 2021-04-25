#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: comp_crf.py
"""
Created on Sun Mar 24 16:14:45 2019

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from astropy.table import Table, join, Column
import astropy.units as u
import numpy as np
from numpy import cos, sqrt
from functools import reduce
import time

__all__ = {"calc_sta_offset", "calc_vel_offset"}


# -----------------------------  FUNCTIONS -----------------------------
def root_sum_square(x, y):
    """Calculate the root-sum-square."""

    return np.sqrt(x**2 + y**2)


def calc_sta_offset(t_sta1,  t_sta2):
    """Calculate station position differences.

    Parameters
    ----------
    t_sta1, t_sta2 : astropy.table object
        station positions from two solutions

    Return
    ------
    t_sta_oft : astropy.table object
        station position differences
    """

    # Copy the original tables and keep only the source position information.
    t_sta3 = Table(t_sta1)
    t_sta3.keep_columns(["station", "xp", "xp_err",
                         "yp", "yp_err", "zp", "zp_err",
                         "xp_yp_corr", "xp_zp_corr", "yp_zp_corr"])

    t_sta4 = Table(t_sta2)
    t_sta4.keep_columns(["station", "xp", "xp_err",
                         "yp", "yp_err", "zp", "zp_err",
                         "xp_yp_corr", "xp_zp_corr", "yp_zp_corr"])

    # Cross-match between two tables
    t_sta_com = join(t_sta3, t_sta4, keys="station")

    print("There are %d and %d stations in two sets, respectively, "
          "between which %d are common."
          % (len(t_sta1), len(t_sta2), len(t_sta_com)))

    # Calculate the offset and the uncertainties
    dxp = t_sta_com["xp_1"] - t_sta_com["xp_2"]
    dyp = t_sta_com["yp_1"] - t_sta_com["yp_2"]
    dzp = t_sta_com["zp_1"] - t_sta_com["zp_2"]
    dxp_err = root_sum_square(t_sta_com["xp_err_1"], t_sta_com["xp_err_2"])
    dyp_err = root_sum_square(t_sta_com["yp_err_1"], t_sta_com["yp_err_2"])
    dzp_err = root_sum_square(t_sta_com["zp_err_1"], t_sta_com["zp_err_2"])

    # Covariance
    dxp_dyp_cov = t_sta_com["xp_err_1"] * t_sta_com["yp_err_1"] \
        * t_sta_com["xp_yp_corr_1"] + t_sta_com["xp_err_2"] \
        * t_sta_com["yp_err_2"] * t_sta_com["xp_yp_corr_2"]
    dxp_dzp_cov = t_sta_com["xp_err_1"] * t_sta_com["zp_err_1"] \
        * t_sta_com["xp_zp_corr_1"] + t_sta_com["xp_err_2"] \
        * t_sta_com["zp_err_2"] * t_sta_com["xp_zp_corr_2"]
    dyp_dzp_cov = t_sta_com["yp_err_1"] * t_sta_com["zp_err_1"] \
        * t_sta_com["yp_zp_corr_1"] + t_sta_com["yp_err_2"] \
        * t_sta_com["zp_err_2"] * t_sta_com["yp_zp_corr_2"]

    # Add these columns to the combined table.
    t_sta_oft = Table([t_sta_com["station"], t_sta_com["xp_1"],
                       t_sta_com["yp_1"], t_sta_com["zp_1"],
                       dxp, dxp_err, dyp, dyp_err, dzp, dzp_err,
                       dxp_dyp_cov, dxp_dzp_cov, dyp_dzp_cov],
                      names=["station", "xp", "yp", "zp",
                             "dxp", "dxp_err", "dyp", "dyp_err",
                             "dzp", "dzp_err",
                             "dxp_dyp_cov", "dxp_dzp_cov", "dyp_dzp_cov"])

    return t_sta_oft


def calc_vel_offset(t_vel1,  t_vel2):
    """Calculate station position differences.

    Parameters
    ----------
    t_vel1, t_vel2 : astropy.table object
        station positions from two solutions

    Return
    ------
    t_vel_oft : astropy.table object
        station position differences
    """

    # Copy the original tables and keep only the source position information.
    t_vel3 = Table(t_vel1)
    t_vel3.keep_columns(["station", "xv", "xv_err",
                         "yv", "yv_err", "zv", "zv_err",
                         "xv_yv_corr", "xv_zv_corr", "yv_zv_corr"])

    t_vel4 = Table(t_vel2)
    t_vel4.keep_columns(["station", "xv", "xv_err",
                         "yv", "yv_err", "zv", "zv_err",
                         "xv_yv_corr", "xv_zv_corr", "yv_zv_corr"])

    # Cross-match between two tables
    t_vel_com = join(t_vel3, t_vel4, keys="station")

    print("There are %d and %d stations in two sets, respectively, "
          "between which %d are common."
          % (len(t_vel1), len(t_vel2), len(t_vel_com)))

    # Calculate the offset and the uncertainties
    dxv = t_vel_com["xv_1"] - t_vel_com["xv_2"]
    dyv = t_vel_com["yv_1"] - t_vel_com["yv_2"]
    dzv = t_vel_com["zv_1"] - t_vel_com["zv_2"]
    dxv_err = root_sum_square(t_vel_com["xv_err_1"], t_vel_com["xv_err_2"])
    dyv_err = root_sum_square(t_vel_com["yv_err_1"], t_vel_com["yv_err_2"])
    dzv_err = root_sum_square(t_vel_com["zv_err_1"], t_vel_com["zv_err_2"])

    # Covariance
    dxv_dyv_cov = t_vel_com["xv_err_1"] * t_vel_com["yv_err_1"] \
        * t_vel_com["xv_yv_corr_1"] + t_vel_com["xv_err_2"] \
        * t_vel_com["yv_err_2"] * t_vel_com["xv_yv_corr_2"]
    dxv_dzv_cov = t_vel_com["xv_err_1"] * t_vel_com["zv_err_1"] \
        * t_vel_com["xv_zv_corr_1"] + t_vel_com["xv_err_2"] \
        * t_vel_com["zv_err_2"] * t_vel_com["xv_zv_corr_2"]
    dyv_dzv_cov = t_vel_com["yv_err_1"] * t_vel_com["zv_err_1"] \
        * t_vel_com["yv_zv_corr_1"] + t_vel_com["yv_err_2"] \
        * t_vel_com["zv_err_2"] * t_vel_com["yv_zv_corr_2"]

    # Add these columns to the combined table.
    t_vel_oft = Table([t_vel_com["station"], t_vel_com["xv_1"],
                       t_vel_com["yv_1"], t_vel_com["zv_1"],
                       dxv, dxv_err, dyv, dyv_err, dzv, dzv_err,
                       dxv_dyv_cov, dxv_dzv_cov, dyv_dzv_cov],
                      names=["station", "xv", "yv", "zv",
                             "dxv", "dxv_err", "dyv", "dyv_err",
                             "dzv", "dzv_err",
                             "dxv_dyv_cov", "dxv_dzv_cov", "dyv_dzv_cov"])

    return t_vel_oft


def calc_trf_offset(t_trf1,  t_trf2):
    """Calculate station position and velocity differences.

    Parameters
    ----------
    t_trf1, t_trf2 : astropy.table object
        station positions from two solutions

    Return
    ------
    t_trf_oft : astropy.table object
        station position and velocity differences
    """

    # Copy the original tables and keep only the source position information.
    t_trf3 = Table(t_trf1)
    t_trf3.keep_columns(["station", "xp", "xp_err",
                         "yp", "yp_err", "zp", "zp_err",
                         "xp_yp_corr", "xp_zp_corr", "yp_zp_corr",
                         "xv", "xv_err", "yv", "yv_err", "zv", "zv_err",
                         "xv_yv_corr", "xv_zv_corr", "yv_zv_corr"])

    t_trf4 = Table(t_trf2)
    t_trf4.keep_columns(["station", "xp", "xp_err",
                         "yp", "yp_err", "zp", "zp_err",
                         "xp_yp_corr", "xp_zp_corr", "yp_zp_corr",
                         "xv", "xv_err", "yv", "yv_err", "zv", "zv_err",
                         "xv_yv_corr", "xv_zv_corr", "yv_zv_corr"])

    # Cross-match between two tables
    t_trf_com = join(t_trf3, t_trf4, keys="station")

    print("There are %d and %d stations in two sets, respectively, "
          "between which %d are common."
          % (len(t_trf1), len(t_trf2), len(t_trf_com)))

    # Calculate the positional offset and the uncertainties
    dxp = t_trf_com["xp_2"] - t_trf_com["xp_1"]
    dyp = t_trf_com["yp_2"] - t_trf_com["yp_1"]
    dzp = t_trf_com["zp_2"] - t_trf_com["zp_1"]
    dxp_err = root_sum_square(t_trf_com["xp_err_1"], t_trf_com["xp_err_2"])
    dyp_err = root_sum_square(t_trf_com["yp_err_1"], t_trf_com["yp_err_2"])
    dzp_err = root_sum_square(t_trf_com["zp_err_1"], t_trf_com["zp_err_2"])

    # Covariance
    dxp_dyp_cov = t_trf_com["xp_err_1"] * t_trf_com["yp_err_1"] \
        * t_trf_com["xp_yp_corr_1"] + t_trf_com["xp_err_2"] \
        * t_trf_com["yp_err_2"] * t_trf_com["xp_yp_corr_2"]
    dxp_dzp_cov = t_trf_com["xp_err_1"] * t_trf_com["zp_err_1"] \
        * t_trf_com["xp_zp_corr_1"] + t_trf_com["xp_err_2"] \
        * t_trf_com["zp_err_2"] * t_trf_com["xp_zp_corr_2"]
    dyp_dzp_cov = t_trf_com["yp_err_1"] * t_trf_com["zp_err_1"] \
        * t_trf_com["yp_zp_corr_1"] + t_trf_com["yp_err_2"] \
        * t_trf_com["zp_err_2"] * t_trf_com["yp_zp_corr_2"]

    # Calculate the velocity offset and the uncertainties
    dxv = t_trf_com["xv_2"] - t_trf_com["xv_1"]
    dyv = t_trf_com["yv_2"] - t_trf_com["yv_1"]
    dzv = t_trf_com["zv_2"] - t_trf_com["zv_1"]
    dxv_err = root_sum_square(t_trf_com["xv_err_1"], t_trf_com["xv_err_2"])
    dyv_err = root_sum_square(t_trf_com["yv_err_1"], t_trf_com["yv_err_2"])
    dzv_err = root_sum_square(t_trf_com["zv_err_1"], t_trf_com["zv_err_2"])

    # Covariance
    dxv_dyv_cov = t_trf_com["xv_err_1"] * t_trf_com["yv_err_1"] \
        * t_trf_com["xv_yv_corr_1"] + t_trf_com["xv_err_2"] \
        * t_trf_com["yv_err_2"] * t_trf_com["xv_yv_corr_2"]
    dxv_dzv_cov = t_trf_com["xv_err_1"] * t_trf_com["zv_err_1"] \
        * t_trf_com["xv_zv_corr_1"] + t_trf_com["xv_err_2"] \
        * t_trf_com["zv_err_2"] * t_trf_com["xv_zv_corr_2"]
    dyv_dzv_cov = t_trf_com["yv_err_1"] * t_trf_com["zv_err_1"] \
        * t_trf_com["yv_zv_corr_1"] + t_trf_com["yv_err_2"] \
        * t_trf_com["zv_err_2"] * t_trf_com["yv_zv_corr_2"]

    # Add these columns to the combined table.
    t_trf_oft = Table([t_trf_com["station"], t_trf_com["xp_1"],
                       t_trf_com["yp_1"], t_trf_com["zp_1"],
                       dxp, dxp_err, dyp, dyp_err, dzp, dzp_err,
                       dxp_dyp_cov, dxp_dzp_cov, dyp_dzp_cov,
                       t_trf_com["xv_1"],
                       t_trf_com["yv_1"], t_trf_com["zv_1"],
                       dxv, dxv_err, dyv, dyv_err, dzv, dzv_err,
                       dxv_dyv_cov, dxv_dzv_cov, dyv_dzv_cov],
                      names=["station",
                             "xp", "yp", "zp", "dxp", "dxp_err",
                             "dyp", "dyp_err", "dzp", "dzp_err",
                             "dxp_dyp_cov", "dxp_dzp_cov", "dyp_dzp_cov",
                             "xv", "yv", "zv", "dxv", "dxv_err",
                             "dyv", "dyv_err", "dzv", "dzv_err",
                             "dxv_dyv_cov", "dxv_dzv_cov", "dyv_dzv_cov"])

    return t_trf_oft


# --------------------------------- END --------------------------------
if __name__ == "__main__":
    print("Nothing to do!")
    pass
# --------------------------------- END --------------------------------
