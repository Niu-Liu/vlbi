#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: eob_comp.py
import time
"""
Created on Wed Feb 27 14:59:41 2019

@author: Neo(liuniu@smail.nju.edu.cn)

Compare the EOP series and calculate the offsets (or differences).

"""

from astropy.table import Table, join, Column
import astropy.units as u
from astropy.units import cds
import numpy as np
import sys
import os
import time


__all__ = {"calc_eop_offset", "save_eop_offset"}


# -----------------------------  FUNCTIONS -----------------------------
def root_sum_square(x, y):
    """Calculate the root-sum-square."""

    return np.sqrt(x**2 + y**2)


def save_eop_offset(eopoft, oftfile):
    """Save the eop offset series into a text file.

    Parameters
    ----------
    eopoft : Table object
        eop offset series
    oftfile : string
        name of file to store the data
    """

    # Header
    eopoft.meta["comments"] = [
        " EOP Offset series",
        " Columns  Units   Meaning",
        "    1     day     Time Tag for polar motion and UT1 (MJD)",
        "    2     day     Time Tag for Nutation (MJD)",
        "    3     uas     offset of X pole coordinate",
        "    4     uas     offset of Y pole coordinate",
        "    5     musec   offset of UT1",
        "    6     musec   offset of LOD",
        "    7     uas     offset of dX of Nutation offsets",
        "    8     uas     offset of dY of Nutation offsets",
        "    9     uas     formal uncertainty for offset of X pole coordinate",
        "   10     uas     formal uncertainty for offset of Y pole coordinate",
        "   11     musec   formal uncertainty for offset of UT1",
        "   12     musec   formal uncertainty for offset of LOD",
        "   13     uas     formal uncertainty for offset of Nutation dX",
        "   14     uas     formal uncertainty for offset of Nutation dY",
        " Created date: %s." % time.strftime("%d/%m/%Y", time.localtime())]

    eopoft.write(oftfile, format="ascii.fixed_width_no_header",
                 exclude_names=["db_name"],
                 formats={"epoch_pmr": "%14.6f", "epoch_nut": "%14.6f",
                          "dxp": "%+8.1f", "dyp": "%+8.1f",
                          "dut": "%+8.1f", "dlod": "%+8.3f",
                          "ddX": "%+8.1f", "ddY": "%+8.1f",
                          "dxp_err": "%8.1f", "dyp_err": "%8.1f",
                          "dut_err": "%8.1f", "dlod_err": "%8.3f",
                          "ddX_err": "%8.1f", "ddY_err": "%8.1f"},
                 delimiter="", overwrite=True)


def read_eop_offset(oftfile):
    """Read EOP offset data.

    Parameters
    ----------
    oftfile : string
        EOP offset file

    Returns
    ----------
    eopoft : astropy.table object
    """

    if not os.path.isfile(oftfile):
        print("Couldn't find the file", oftfile)
        sys.exit()

    eopoft = Table.read(oftfile, format="ascii",
                        names=["epoch_pmr", "epoch_nut",
                               "dxp", "dyp", "dut", "dlod", "ddX", "ddY",
                               "dxp_err", "dyp_err", "dut_err", "dlod_err",
                               "ddX_err", "ddY_err"])

    # Add unit information
    eopoft["epoch_pmr"].unit = cds.MJD
    eopoft["epoch_nut"].unit = cds.MJD

    eopoft["dxp"].unit = u.uas
    eopoft["dyp"].unit = u.uas
    eopoft["dxp_err"].unit = u.uas
    eopoft["dyp_err"].unit = u.uas

    eopoft["dut"].unit = u.second / 1e6
    eopoft["dut_err"].unit = u.second / 1e6
    eopoft["dut"].unit = u.second / 1e6
    eopoft["dut_err"].unit = u.second / 1e6

    eopoft["ddX"].unit = u.uas
    eopoft["ddY"].unit = u.uas
    eopoft["ddX_err"].unit = u.uas
    eopoft["ddY_err"].unit = u.uas

    return eopoft


def calc_eop_offset(eop1,  eop2, oftfile=None):
    """Calculate the EOP difference series between two solutions.

    Parameters
    ----------
    eop1, eop2 : astropy.table object
        EOP series from two solutions

    Return
    ------
    eopoft : astropy.table object
        EOP difference series
    """

    # Copy the original tables and keep only the EOP information
    eop3 = Table(eop1)
    eop3.keep_columns(["db_name", "epoch_pmr", "epoch_nut",
                       "xp", "yp", "ut1_tai", "dX", "dY", "lod",
                       "xp_err", "yp_err", "ut1_err", "dX_err", "dY_err",
                       "lod_err"])

    eop4 = Table(eop2)
    eop4.keep_columns(["db_name",
                       "xp", "yp", "ut1_tai", "dX", "dY", "lod",
                       "xp_err", "yp_err", "ut1_err", "dX_err", "dY_err",
                       "lod_err"])

    # Cross-match between two tables
    eopcom = join(eop3, eop4, keys="db_name")

    print("There are %d and %d points in series 1 and series 2, respectively,"
          "between which %d are common."
          % (len(eop1), len(eop2), len(eopcom)))

    # Calculate the offset and the uncertainties
    dxp = eopcom["xp_1"] - eopcom["xp_2"]
    dyp = eopcom["yp_1"] - eopcom["yp_2"]
    dut = eopcom["ut1_tai_1"] - eopcom["ut1_tai_2"]
    ddX = eopcom["dX_1"] - eopcom["dX_2"]
    ddY = eopcom["dY_1"] - eopcom["dY_2"]
    dlod = eopcom["lod_1"] - eopcom["lod_2"]

    dxperr = root_sum_square(eopcom["xp_err_1"], eopcom["xp_err_2"])
    dyperr = root_sum_square(eopcom["yp_err_1"], eopcom["yp_err_2"])
    duterr = root_sum_square(eopcom["ut1_err_1"], eopcom["ut1_err_2"])
    ddXerr = root_sum_square(eopcom["dX_err_1"], eopcom["dX_err_2"])
    ddYerr = root_sum_square(eopcom["dY_err_1"], eopcom["dY_err_2"])
    dloderr = root_sum_square(eopcom["lod_err_1"], eopcom["lod_err_2"])

    # Convert the unit
    # Time tag
    # from astropy.time import Time
    # t_pmr_mjd = Time(eopcom["epoch_pmr"], format="mjd")
    # t_pmr = Column(t_pmr_mjd.jyear, unit=u.year)
    #
    # t_nut_mjd = Time(eopcom["epoch_nut"], format="mjd")
    # t_nut = Column(t_nut_mjd.jyear, unit=u.year)

    # Polar motion (as -> uas)
    dxp.convert_unit_to(u.uas)
    dxperr.convert_unit_to(u.uas)
    dyp.convert_unit_to(u.uas)
    dyperr.convert_unit_to(u.uas)

    # UT1-UTC and lod (s -> us)
    dut.convert_unit_to(u.second / 1e6)
    duterr.convert_unit_to(u.second / 1e6)
    dlod.convert_unit_to(u.second / 1e6)
    dloderr.convert_unit_to(u.second / 1e6)

    # Nutation offset (as -> uas)
    ddX.convert_unit_to(u.uas)
    ddXerr.convert_unit_to(u.uas)
    ddY.convert_unit_to(u.uas)
    ddYerr.convert_unit_to(u.uas)

    # Add these columns to the combined table.
    eopoft = Table([eopcom["db_name"], eopcom["epoch_pmr"],
                    eopcom["epoch_nut"],
                    dxp, dyp, dut, dlod, ddX, ddY,
                    dxperr, dyperr, duterr, dloderr, ddXerr, ddYerr],
                   names=["db_name", "epoch_pmr", "epoch_nut",
                          "dxp", "dyp", "dut", "dlod", "ddX", "ddY",
                          "dxp_err", "dyp_err", "dut_err", "dlod_err",
                          "ddX_err", "ddY_err"])

    # Save the EOP offset series
    if oftfile is not None:
        print("Save the EOP offset series in", oftfile)
        save_eop_offset(eopoft, oftfile)

    return eopoft


if __name__ == '__main__':
    print('Nothing to do!')
# --------------------------------- END --------------------------------
