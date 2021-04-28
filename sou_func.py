#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: read_sou.py
"""
Created on Sat Sep 15 14:38:46 2018

@author: Neo(liuniu@smail.nju.edu.cn)

Format for the .crf file
IERS ICRS-PC format: https://ivscc.gsfc.nasa.gov/products-data/crf-format.txt.

"""

from astropy.table import Table, join, unique, Column
import astropy.units as u
from astropy.units import cds
from astropy.coordinates import SkyCoord, Angle
import numpy as np
import os
import sys
import time


# My modules
from .convert_func import RA_conv, DC_conv, date2mjd
from .pos_err import pos_err_calc
from .sou_name import get_souname

__all__ = ["read_sou", "read_opa_crf", "read_cat"]


# -----------------------------  FUNCTIONS -----------------------------
def inflate_sou_err(err, scale=1.5, noise=0.03):
    """Inflate the VLBI formal error

    Parameters
    ----------
    err : ndarray
        VLBI formal error from global solution, unit is mas

    Return
    ------
    err1 : ndarray
        inflated formal error
    """

    err = np.sqrt((err * scale)**2 + noise**2)

    return err


def drop_sou_with_few_obs(t_sou, nobs_lim=3, verbose=True):
    """Drop sources with few observations

    Parameters
    ----------
    t_sou : Astropy.table.Table object
        catalog
    nobs_lim : int
        least number of observation, meaningful only if drop_few_obs is True.

    Return
    ------
    A new catalog
    """

    # NO. sources in the original solution
    N0 = len(t_sou)

    # Remove source with few observation used in the solution
    mask = (t_sou["used_obs"] >= nobs_lim)
    t_sou = Table(t_sou[mask], masked=False)

    # No. sources after elimination
    N1 = len(t_sou)

    if verbose:
        print("There are {:d} sources in the original catalog, "
              "{:d} ({:.0f}%) sources with #obs < {:.0f} dropped, leaving "
              "{:d} sources in the present catalog.".format(
                  N0, N0-N1, (N0-N1)/N0*100, nobs_lim, N1))

    return t_sou


def add_unit(t_sou, table_type="sou"):
    """Add unit information, work for .sou, .crf, and .cat file
    """

    if table_type in ["sou", "cat"]:
        t_sou["ra"].unit = u.deg
        t_sou["dec"].unit = u.deg

    if table_type in ["crf"]:
        t_sou["ra_err"] = t_sou["ra_err"] * 15 * 3.6e6  # second -> mas
        t_sou["dec_err"] = t_sou["dec_err"] * 3.6e6

    t_sou["ra_err"].unit = u.mas
    t_sou["dec_err"].unit = u.mas

    t_sou["beg_epoch"].unit = cds.MJD
    t_sou["end_epoch"].unit = cds.MJD
    if "mean_epoch" in t_sou.colnames:
        t_sou["mean_epoch"].unit = cds.MJD

    return t_sou


def ra_err_corr(t_sou):
    """Multipliy the formal error in R.A. by a factor of cos(Decl.)
    """

    factor = np.cos(Angle(t_sou["dec"]).radian)
    t_sou["ra_err"] = t_sou["ra_err"] * factor

    return t_sou


def read_sou(sou_file, drop_few_obs=False, nobs_lim=3, flate_err=False):
    """Read radio source positions

    Parameters
    ----------
    sou_file : string
        the full path of .sou file
    drop_few_obs : boolean
        flag to determine whether to remove sources with few observations

    flate_err : boolean
        flag to determine whether to inflate the formal error

    Return
    ------
    t_sou : astropy.table object
        |
        -- ivs_name : str
            IVS source name
        -- ra : float
            right ascension (degree)
        -- ra_err
            formal uncertainty in RA (mas)
        -- dec
            declination (degree)
        -- dec_err
            formal uncertainty in Dec. (mas)
        -- ra_dec_corr
            correlation coefficient between RA and Dec.
        -- pos_err
            ellipse semi-major axis of positional error (mas)
        -- used_obs
            number of observations used in the solution
        -- total_obs
            number of total observations of this source
        -- used_ses
            number of sessions used in the solution
        -- total_ses
            number of total sessions of this source
        -- beg_epoch
            epoch of first observation (MJD)
        -- end_epoch
            epoch of last observation (MJD)
    """

    if not os.path.isfile(sou_file):
        sys.exit()

    t_sou = Table.read(sou_file, format="ascii.fixed_width_no_header",
                       names=("ivs_name", "ra", "ra_err",
                              "dec", "dec_err", "ra_dec_corr",
                              "used_obs", "total_obs", "used_ses",
                              "total_ses", "beg_epoch", "end_epoch"),
                       col_starts=(10, 24, 45, 61, 82, 98, 117,
                                   132, 150, 164, 181, 202),
                       col_ends=(18, 41, 55, 78, 92, 104, 122,
                                 139, 155, 170, 191, 212))

    ra_dec_table = Table.read(sou_file, format="ascii.fixed_width_no_header",
                              names=["ra", "dec", "used_obs"],
                              col_starts=[20, 57, 117],
                              col_ends=[41, 78, 122])

    if drop_few_obs:
        t_sou = drop_sou_with_few_obs(t_sou, nobs_lim)
        ra_dec_table = drop_sou_with_few_obs(ra_dec_table, nobs_lim, verbose=False)

    # convert string into float for RA, Decl. and observing epoch
    Nsize = len(t_sou)
    ra = np.empty(shape=(Nsize,), dtype=float)
    dec = np.empty(shape=(Nsize,), dtype=float)
    date_beg = np.empty(shape=(Nsize,), dtype=float)
    date_end = np.empty(shape=(Nsize,), dtype=float)

    for i, ra_dec_tablei in enumerate(ra_dec_table):
        ra[i] = RA_conv(ra_dec_tablei["ra"])
        dec[i] = DC_conv(ra_dec_tablei["dec"])

    for i, t_soui in enumerate(t_sou):
        if t_soui["used_obs"]:
            date_beg[i] = date2mjd(t_soui["beg_epoch"])
            date_end[i] = date2mjd(t_soui["end_epoch"])
        else:
            date_beg[i] = 0
            date_end[i] = 0

    # replace original columns with new columns
    t_sou["ra"] = ra
    t_sou["dec"] = dec
    t_sou["beg_epoch"] = date_beg
    t_sou["end_epoch"] = date_end

    # Add unit information
    t_sou = add_unit(t_sou, "sou")

    # Multipliy the formal error in R.A. by a factor of cos(Decl.)
    t_sou = ra_err_corr(t_sou)

    # Inflate the formal error
    if flate_err:
        t_sou["ra_err"] = inflate_sou_err(t_sou["ra_err"])
        t_sou["dec_err"] = inflate_sou_err(t_sou["dec_err"])

    # Calculate the semi-major axis of error ellipse
    pos_err = pos_err_calc(t_sou["ra_err"], t_sou["dec_err"],
                           t_sou["ra_dec_corr"])

    # Add the semi-major axis of error ellipse to the table
    t_sou.add_column(pos_err, name="pos_err", index=6)

    # Add IERS and ICRF designations of source names
    t_souname = get_souname()
    t_sou = join(t_souname, t_sou, keys="ivs_name", join_type="right")

    # Fill the empty filed of IERS name by the IVS name
    for i in t_sou["iers_name"].mask.nonzero()[0]:
        t_sou[i]["iers_name"] = t_sou[i]["ivs_name"]

    return t_sou


def read_asi_crf(crf_file):
    """Read radio source positions in .crf file from ASI

    Parameters
    ----------
    crf_file : string
        the full path of .cat file

    Return
    ------
    t_crf : astropy.table object
    """

    t_crf = Table.read(crf_file,
                       format="ascii.fixed_width_no_header",
                       names=("ivs_name", "ra_err", "dec_err", "ra_dec_corr",
                              "mean_epoch", "beg_epoch", "end_epoch", "used_ses",
                              "used_obs", "num_delrate"),
                       col_starts=(1, 55, 66, 81, 93, 106, 119, 132, 137, 144),
                       col_ends=(9, 65, 75, 87, 105, 118, 131, 136, 143, 150))

    # Position
    t_radec = Table.read(crf_file, format="ascii.fixed_width_no_header",
                         names=["ra_dec"], col_starts=[19], col_ends=[54])

    return t_crf, t_radec


def read_opa_crf(crf_file):
    """Read radio source positions in .crf file from OPA/AUS

    Parameters
    ----------
    crf_file : string
        the full path of .cat file

    Return
    ------
    t_crf : astropy.table object
    """

    if not os.path.isfile(crf_file):
        sys.exit(1)

    t_crf = Table.read(crf_file, format="ascii",
                       names=["ivs_name", "iers_name",
                              "ra_h", "ra_m", "ra_s",
                              "dec_d", "dec_m", "dec_s",
                              "ra_err", "dec_err", "ra_dec_corr",
                              "mean_epoch", "beg_epoch", "end_epoch",
                              "used_ses", "used_obs", "num_delrate", "flag"],
                       exclude_names=["ra_h", "ra_m", "ra_s",
                                      "dec_d", "dec_m", "dec_s"])

    # Position
    t_radec = Table.read(crf_file, format="ascii.fixed_width_no_header",
                         names=["ra_dec"], col_starts=[21], col_ends=[59])

    return t_crf, t_radec


def read_usn_crf(crf_file):
    """Read radio source positions in .crf file from USNO

    Parameters
    ----------
    crf_file : string
        the full path of .cat file

    Return
    ------
    t_crf : astropy.table object
    """

    if not os.path.isfile(crf_file):
        sys.exit(1)

    t_crf = Table.read(crf_file, data_start=6,
                       format="ascii.fixed_width_no_header",
                       names=("ivs_name", "ra_err", "dec_err", "ra_dec_corr",
                              "mean_epoch", "beg_epoch", "end_epoch", "used_ses",
                              "used_obs", "num_delrate"),
                       col_starts=(9, 54, 65, 75, 82, 90, 98, 106, 112, 124),
                       col_ends=(17, 64, 74, 81, 89, 97, 105, 111, 118, 125))

    # Position
    t_radec = Table.read(crf_file, data_start=6,
                         format="ascii.fixed_width_no_header",
                         names=["ra_dec"], col_starts=[18], col_ends=[53])

    return t_crf, t_radec


def read_crf(crf_file, analy_cen="OPA", drop_few_obs=False, nobs_lim=3, flate_err=False):
    """Read radio source positions in .crf file

    Parameters
    ----------
    crf_file : string
        the full path of .cat file
    analy_cen : string
        abbreviation name of analysis center
        asi: Space Geodesy Centre of the Italian Space Agency
        aus: Geoscience Australia
        opa: Paris Observatory, France
        usn: Unites States Naval Observatory
    drop_few_obs : boolean
        flag to determine whether to remove sources with few observations
    nobs_lim : int
        least number of observation, meaningful only if drop_few_obs is True.
    flate_err : boolean
        flag to determine whether to inflate the formal error

    Return
    ------
    t_crf : astropy.table object
        |
        -- ivs_name : str
            IVS source name
        -- iers_name : str
            IERS source name
        -- ra : float
            right ascension (degree)
        -- ra_err
            formal uncertainty in RA*cos(decl.) (mas)
        -- dec
            declination (degree)
        -- dec_err
            formal uncertainty in Dec. (mas)
        -- ra_dec_corr
            correlation coefficient between RA and Dec.
        -- pos_err
            ellipse semi-major axis of positional error (mas)
        -- used_obs
            number of observations used in the solution
        -- total_obs
            number of total observations of this source
        -- used_ses
            number of sessions used in the solution
        -- total_ses
            number of total sessions of this source
        -- beg_epoch
            epoch of first observation (MJD)
        -- end_epoch
            epoch of last observation (MJD)
    """

    if analy_cen in ["asi", "ASI"]:
        t_crf, t_radec = read_asi_crf(crf_file)
    elif analy_cen in ["aus", "AUS", "opa", "OPA"]:
        t_crf, t_radec = read_opa_crf(crf_file)
    elif analy_cen in ["usn", "USN"]:
        t_crf, t_radec = read_usn_crf(crf_file)
    else:
        print("Undefined IVS analysis center")
        os.exit(1)

    # RA and decl.
    ra_dec = SkyCoord(t_radec["ra_dec"], unit=(u.hourangle, u.deg))
    ra_col = Column(ra_dec.ra.deg, name="ra", unit=u.deg)
    dec_col = Column(ra_dec.dec.deg, name="dec", unit=u.deg)
    t_crf.add_columns([ra_col, dec_col], indexes=[1, 2])

    # Multipliy the formal error in R.A. by a factor of cos(Decl.)
    factor = np.cos(ra_dec.dec.radian)
    t_crf["ra_err"] = t_crf["ra_err"] * factor

    # Drop sources with fewer No. obs
    if drop_few_obs:
        t_crf = drop_sou_with_few_obs(t_crf, nobs_lim)

    # Add unit information
    t_crf = add_unit(t_crf, "crf")

    # Inflate the formal error
    if flate_err:
        t_crf["ra_err"] = inflate_sou_err(t_crf["ra_err"])
        t_crf["dec_err"] = inflate_sou_err(t_crf["dec_err"])

    # Calculate the semi-major axis of error ellipse
    pos_err = pos_err_calc(t_crf["ra_err"], t_crf["dec_err"],
                           t_crf["ra_dec_corr"])

    # Add the semi-major axis of error ellipse to the table
    t_crf.add_column(pos_err, name="pos_err", index=7)
    # t_crf["pos_err"].unit = u.mas

    # Add IERS and ICRF designations of source names
    if analy_cen in ["asi", "ASI", "usn", "USN"]:
        t_souname = get_souname()
        t_crf = join(t_souname, t_crf, keys="ivs_name", join_type="right")

        # Fill the empty filed of IERS name by the IVS name
        for i in t_crf["iers_name"].mask.nonzero()[0]:
            t_crf[i]["iers_name"] = t_crf[i]["ivs_name"]

    return t_crf


def read_cat(cat_file, drop_few_obs=False, nobs_lim=3, flate_err=False):
    """Read radio source positions

    Parameters
    ----------
    cat_file : string
        the full path of .cat file
    drop_few_obs : boolean
        flag to determine whether to remove sources with few observations
    nobs_lim : int
        least number of observation, meaningful only if drop_few_obs is True.
    flate_err : boolean
        flag to determine whether to inflate the formal error

    Return
    ------
    t_cat : astropy.table object
        |
        -- ivs_name : str
            IVS source name
        -- iers_name : str
            IERS source name
        -- ra : float
            right ascension (degree)
        -- ra_err
            formal uncertainty in RA (mas)
        -- dec
            declination (degree)
        -- dec_err
            formal uncertainty in Dec. (mas)
        -- ra_dec_corr
            correlation coefficient between RA and Dec.
        -- pos_err
            ellipse semi-major axis of positional error (mas)
        -- used_obs
            number of observations used in the solution
        -- total_obs
            number of total observations of this source
        -- used_ses
            number of sessions used in the solution
        -- total_ses
            number of total sessions of this source
        -- beg_epoch
            epoch of first observation (MJD)
        -- end_epoch
            epoch of last observation (MJD)
    """

    if not os.path.isfile(cat_file):
        sys.exit()

    ivs_name, iers_name = np.genfromtxt(
        cat_file, dtype=str, usecols=(0, 1), unpack=True)
    ra, ra_err, dec, dec_err, ra_dec_corr = np.genfromtxt(
        cat_file, usecols=range(2, 7), unpack=True)
    used_ses, used_obs = np.genfromtxt(
        cat_file, dtype=int, usecols=(7, 8), unpack=True)
    mean_epoch, beg_epoch, end_epoch = np.genfromtxt(
        cat_file, usecols=(9, 10, 11), unpack=True)

    t_cat = Table([ivs_name, iers_name,
                   ra, dec, ra_err, dec_err, ra_dec_corr,
                   used_ses, used_obs, mean_epoch, beg_epoch, end_epoch],
                  names=["ivs_name", "iers_name",
                         "ra", "dec", "ra_err", "dec_err", "ra_dec_corr",
                         "used_ses", "used_obs", "mean_epoch", "beg_epoch", "end_epoch"])

    if drop_few_obs:
        t_cat = drop_sou_with_few_obs(t_cat, nobs_lim)

    # unit
    t_cat = add_unit(t_cat, "cat")

    # Inflate the formal error
    if flate_err:
        t_cat["ra_err"] = inflate_sou_err(t_cat["ra_err"])
        t_cat["dec_err"] = inflate_sou_err(t_cat["dec_err"])

    # Calculate the semi-major axis of error ellipse
    pos_err = pos_err_calc(t_cat["ra_err"], t_cat["dec_err"],
                           t_cat["ra_dec_corr"])

    # Add the semi-major axis of error ellipse to the table
    t_cat.add_column(pos_err, name="pos_err", index=7)
    t_cat["pos_err"].unit = u.mas

    return t_cat


# -------------------------------- MAIN --------------------------------
if __name__ == "__main__":
    read_sou(sys.argv[1])
# --------------------------------- END --------------------------------
