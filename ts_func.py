#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: read_ts.py
"""
Created on Wed Feb  5 13:44:22 2020

@author: Neo(liuniu@smail.nju.edu.cn)

Read the radio source coordinate time-series and return an astropy.table.Table object.
"""

from astropy.table import Table, Column
from astropy.time import Time

import astropy.units as u
import numpy as np

# My mudoles
from my_progs.catalog.read_icrf import read_icrf3
from my_progs.catalog.pos_diff import pa_calc
from .get_dir import get_vlbi_sol_dir


# -----------------------------  FUNCTIONS -----------------------------
def def_data_dir():
    """
    Return the default directory of radio source coordinate time series
    """

    sol_dir = get_vlbi_sol_dir()
    data_dir = "{:s}/opa/ts-sou".format(sol_dir)

    return data_dir


def mjd2jyear(mjd):
    """Convert MJD to JY.

    Parameter
    ---------
    mjd : 1darray-like
        mean Julian day

    Return
    ------
    jyear : 1darray-like
        mean Julian year
    """

    epoch = Time(mjd, format="mjd")
    jyear = epoch.jyear

    return jyear


def read_ts(data_file, data_dir=None):
    """Read radio source coordinate time-series data

    Parameters
    ----------
    data_file: string
        data file name
    data_dir: string
        directory to store the data file

    Returns
    -------
    coordts : astropy.table.Table object
        coordinate time series
    """

    if not data_dir:
        data_dir = def_data_dir()

    coordts = Table.read("{}/{}".format(data_dir, data_file), format="ascii",
                         names=["mjd", "ra", "dec", "ra_err", "dec_err", "ra_dec_corr",
                                "used_obs", "iers_name", "ivs_name", "db_name"])

    # Add unit
    coordts["ra"].unit = u.deg
    coordts["dec"].unit = u.deg
    coordts["ra_err"].unit = u.mas
    coordts["dec_err"].unit = u.mas

    jyear = mjd2jyear(coordts["mjd"])
    jyear = Column(jyear, name="jyear")
    coordts.add_column(jyear)

    return coordts


def get_sou_ts(sou_name, data_dir=None):
    """Read radio source coordinate time-series data

    Parameters
    ----------
    sou_name: string
        source name of IERS designation or file name with the data
    data_dir: string
        directory to store the data file

    Returns
    -------
    coordts : astropy.table.Table object
        coordinate time series
    """

    if not data_dir:
        data_dir = def_data_dir()

    data_file = "{}.dat".format(sou_name)
    coordts = read_ts(data_file, data_dir)

    return coordts


def ts_stat(coordts):
    """Simple statistics of coordinate time-series

    Parameter
    ---------
    sou_name: string
        source name of IERS designation or file name with the data
    """

    stat = {}

    # Number of used sessions
    numses = len(coordts)
    stat["used_ses"] = numses

    # Number of used observations
    numobs = np.sum(coordts["used_obs"])
    stat["used_obs"] = numobs

    #### TO-DO ####
    # MAYBE PRINT THESE STATISTICS
    #### TO-DO ####


def sou_ts_stat(sou_name):
    """Simple statistics of coordinate time-series

    Parameter
    ---------
    sou_name: string
        source name
    """

    coordts = get_sou_ts(sou_name)

    ts_stat(coordts)

    # NOT CLEAR YET ABOUT WHAT TO DO NEXT


def get_icrf3_pos(sou_name, wv="sx"):
    """Get radio source position from the ICRF3 catalog

    Parameters
    ----------
    sou_name: string
        source name
    wv: str,
        wavelength, SX, K, or XKa

    Returns
    -------
    pos : a dict object
        key: "ra", "dec"
        value: ra, dec in degree
    """

    icrf3 = read_icrf3(wv=wv)

    # Find the source position
    mask = (icrf3["iers_name"] == sou_name)
    sou = icrf3[mask]

    if len(sou):
        ra = sou["ra"][0]
        dec = sou["dec"][0]

        if sou["type"][0] == "C":
            sou_type = "C"
        else:
            sou_type = "O"

    else:
        print("There is no entry for {} in ICRF3 catalog,".format(sou_name),
              "use the mean position instead.")
        return False

    pos = {"ra": ra, "dec": dec, "type": sou_type}

    return pos


def calc_ts_oft(sou_name, coordts):
    """Calculate the time series offset wrt. ICRF3 position

    Parameters
    ----------
    sou_name: string
        source name
    coordts : astropy.table.Table object
        coordinate time series

    Return
    ------
    dpos : a dict object
        key: "dra", "ddec", "rho", "pa"
        value: dra, ddec, rho in mas
               pa in deg
    """

    pos = get_icrf3_pos(sou_name)

    if pos:
        ra0 = pos["ra"]
        dec0 = pos["dec"]
        sou_type = pos["type"]
    else:
        ra0 = np.mean(coordts["ra"])
        dec0 = np.mean(coordts["dec"])
        sou_type = "N"

    dra = coordts["ra"] - ra0
    ddec = coordts["dec"] - dec0

    dra.unit = u.deg
    ddec.unit = u.deg
    dra.convert_unit_to(u.mas)
    ddec.convert_unit_to(u.mas)

    rho = np.sqrt(dra**2 + ddec**2)
#     pax, pay = pa_calc(np.array(dra), np.array(ddec))
    pa = pa_calc(np.array(dra), np.array(ddec))

    dpos = {"dra": dra, "ddec": ddec, "rho": rho, "pa": pa, "type": sou_type}

    return dpos


def get_ts(sou_name, data_dir=None, calc_oft=False):
    """get the time series offset wrt. ICRF3 position

    Parameters
    ----------
    sou_name: string
        source name

    Return
    ------
    dpos : a dict object
        key: "dra", "ddec"
        value: dra, ddec in mas
    """

    if not data_dir:
        data_dir = def_data_dir()

    coordts = get_sou_ts(sou_name, data_dir)

    if calc_oft:
        dpos = calc_ts_oft(sou_name, coordts)
        dra = dpos["dra"]
        ddec = dpos["ddec"]
        rho = dpos["rho"]
        pa = dpos["pa"]

        # Check the source type in the ICRF3 S/X catalog
        type_list = [dpos["type"]] * len(coordts)
        sou_type = Column(type_list, name="type")

        coordts.add_columns([dra, ddec, rho, pa, sou_type], names=[
                            "dra", "ddec", "rho", "pa", "type"])

    return coordts


# --------------------------------- END --------------------------------
