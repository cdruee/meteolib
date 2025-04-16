#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Funtions representing the 1975 ISO Standard Atmosphere [ISO2533]_,
wich is identical to the ICAO standard atmosphere [ICAO7488]_
'''
import numpy as np
import pandas as pd

# from ._utils import _check, _only
from .constants import R, r_earth
from .pressure import _to_Pa
from .temperature import _to_K, _to_C
from .thermodyn import gas_rho

# Table 4 from [ISO2533]_ defining the temperature profile
# with respect to geopotential height
_tab4 = pd.DataFrame.from_records([
    [-2000, 301.15, -6.5E-3],
    [0, 288.15, -6.5E-3],
    [11000, 216.65,  0.0E0],
    [20000, 216.65,  1.0E-3],
    [32000, 228.65,  2.8E-3],
    [47000, 270.65,  0.0E0],
    [51000, 270.65, -2.8E-3],
    [71000, 214.65, -2.0E-3],
    [80000, 196.65, 0.0E0]],
    columns=["H_b", "T_b", "beta"])

# Table 5 from [ISO2533]_ giving tabbed values
# with respect to geopotential height
_tab5 = pd.DataFrame.from_records([
    [-2000, 127783., 1.47816E-0],
    [0, 101325., 1.22500E-0],
    [11000, 22632.0, 3.63918E-1],
    [20000, 5474.87, 8.80345E-2],
    [32000, 868.014, 1.35551E-2],
    [47000, 110.906, 1.42752E-3],
    [51000, 66.9384, 8.61600E-4],
    [71000, 3.95639, 6.42105E-5],
    [80000, 0.886272, 1.57004E-5]],
    columns=["H", "p", "rho"])

pzero_iso = 101325.  # Pa
R_iso = 287.05287    # K

# ---------------------------------------------------------------------


def beta_iso(H):
    if H < _tab4["H_b"][0]:
        res = None
    else:
        for i in _tab4.index:
            if H <= _tab4["H_b"][i]:
                res = _tab4["beta"][i]
                break
    return res


# ----------------------------------------------------------------------
def gphi(lat=None, h=0, deg=True):
    r"""
    Lambert's function for acceleration of gravity as function of latitude.
    Standard ("normal") gravity conforms with latitude
    :math:`\phi = 45^\circ32'33"`.
    For height dependence, centrifugal acceleration is
    formally neglected and using only Newton’s gravitation law.
    Defined in ICAO/ANSI/ISO standard atmosphere [ISO2533]_

    :param lat: (float, optional) latitude in radians
        (or degrees, if `deg` is True),
        positive to the north. If missing, standard gravity is returned.
    :param deg: (float, altitude) is geometric height.
    :param deg: (bool, optional) if True `lat` is in degrees,
        if False `lat` is in radians. Defaults to True (i.e. degrees).

    :returns: gravity im :math:`m/s^{-2}` (float)
    :rtype: float
    """
    if lat is None:
        lat = 45. + 32./60. + 33./3600.
        deg = True
    if deg:
        phi = np.deg2rad(lat)
    else:
        phi = lat
    # ISO 2533 Sect. 2.1 (unnumbered eqn)
    res = 9.80616 * (1 -
                     0.0026373 * np.cos(2 * phi) +
                     0.0000059 * np.square(np.cos(2 * phi)))

    res = res * (r_earth/(r_earth + h)) ** 2
    return res

# ---------------------------------------------------------------------


def T_iso(h, gpm=False, Kelvin=True):
    """
    Temperature in the standard atmosphere at given height

    :param h: (float) height above sea level (altitude) in m
    :param gpm: (bool, optional) if True, h is given in gepotential meters,
        if False h is given in gemetric height. Default to False.
    :param Kelvin: (bool, optional) if True temperature is returned in Kelvin,
        if False, temperature is returned in degree Celsius.
        Defaults to True.
    :returns: Temperature in Kelvin or degrees Celsius,
        depending on `Kelvin`
    :rtype: float
    """
    if not gpm:
        alt = h_geopot(h)
    else:
        alt = h
    for i in reversed(_tab4.index):
        if alt >= _tab4["H_b"][i]:
            T = (_tab4["T_b"][i] +
                 (alt - _tab4["H_b"][i]) * _tab4["beta"][i])
            break
    if Kelvin is True:
        return T
    else:
        return _to_C(T, Kelvin=True)

# ---------------------------------------------------------------------


def h_geopot(h_geomet):
    """
    converts geometric altitude to geopotential altitude
    :param h_geomet: (float) geometric altitude in m

    :returns: geopotential altitude in gpm
    :rtype: float
    """
    # ISO 2533 eqn 8
    return (r_earth * h_geomet) / (r_earth + h_geomet)

# ---------------------------------------------------------------------


def h_geomet(h_geopot):
    """
    converts geometric altitude to geopotential altitude
    :param h_geopot: (float) geopotential altitude in gpm

    :returns: geometric altitude in m
    :rtype: float
    """
    # ISO 2533 eqn 9
    return (r_earth * h_geopot) / (r_earth - h_geopot)

# ---------------------------------------------------------------------


def p_barom(h, p_0=pzero_iso, T_0=_tab4["T_b"][1], h_0=0,
            beta=_tab4["beta"][1],
            gpm=False, Kelvin=None, hPa=None):
    r"""
    barometric height formula

    :param h: (float) height above sea level (altitude) (:math:`m`)
    :param p_0: (float, optional) air pressure at base level
      (in hPa or Pa, depending on `hPa`).
      Defaults to `pzero_iso`
    :param T_0: (float, optional) air temperature at base level
      in :math:`K` or :math:`^\circ C` , depending on `Kelvin` .
      Defaults to 288.15 K.
    :param h_0: (float, optional) height of the base level
      above sea level (altitude of the base level) (:math:`m`).
      Defaults to 0 m.
    :param beta: (float, optional) vertical temperature gradient
      (:math:`K m^{-1}`). Defaults to 0.0065 :math:`K m^{-1}`.
    :param gpm: (bool, optional) if True, h is given in gepotential meters,
      if False h is given in gemetric height. Default to False.
    :param Kelvin: (optional, optional)
      if ``False``, all temperatures are assumed to be Kelvin.
      If ``False``, all temperatures are assumed to be Celsius.
      If missing of ``None``, unit  temperatures are autodetected.
      Defaults to ``None``.
    :param hPa: (optional, optional)
      if ``True``,  `p`, `e`, or `ew` must be supplied in hPa.
      If ``False``, `p`, `e`, or `ew` must be supplied in Pa.
      Defaults to ``False``.

    :return: air pressure (in hPa or Pa, depending on `hPa`).
    :rtype: float

    """
    Tnull = _to_K(T_0, Kelvin)  # K
    pnull = _to_Pa(p_0, hPa)    # Pa
    if not gpm:
        hnull = h_geopot(h_0)   # gpm
        alt = h_geopot(h)       # gpm
    else:
        hnull = h_0             # gpm
        alt = h                 # gpm
    gamma = beta                # K/m

    if abs(gamma) > 1.E-6:  # i.e. is not zero
        # ISO 2533 eqn 12
        pwr = - gphi() / (R_iso * gamma)
        pp = p_0 * (1. + (gamma * (alt - hnull) / Tnull)) ** pwr  # Pa
    else:
        # ISO 2533 eqn 13
        pp = p_0 * np.exp(- (gphi() / (R * Tnull)) * (alt - hnull))  # Pa

    if hPa is None:
        # convert to the same unit as supplied p_0 (Pa if p_0 not given)
        p = pp * p_0 / pnull
    elif hPa is True:
        p = pp / 100.
    else:
        p = pp
    return p

# ---------------------------------------------------------------------


def p_iso(h, p_0=None, gpm=False, hPa=False):
    """
    air pressure in the standard atmosphere at given height

    :param h: (float) height above sea level (altitude) (:math:`m`)
    :param p_0: (float, optional) air pressure at sea level
        (in hPa or Pa, depending on `hPa`).
        Defaults to `pzero_iso`
    :param gpm: (bool, optional) if True, h is given in gepotential meters,
        if False h is given in gemetric height. Default to False.
    :param hPa: (optional, optional)
      if ``True``,  `p` must be supplied in hPa.
      If ``False``, `p` must be supplied in Pa.
      Defaults to ``False``.

    :return: air pressure (in hPa or Pa, depending on `hPa`).
    :rtype: float
    """
    if p_0 is None:
        pp_0 = pzero_iso             # Pa
    else:
        pp_0 = _to_Pa(p_0, hPa=hPa)  # Pa
    # geopotential height
    if not gpm:
        alt = h_geopot(h)  # gpm
    else:
        alt = h            # gpm
    if alt < _tab4["H_b"][0]:
        integr = np.nan
    elif _tab4["H_b"][0] < alt < 0.:
        integr = p_barom(alt, p_0=pp_0, gpm=True, hPa=False)  # Pa
    else:
        pfact = pp_0 / pzero_iso
        for i in reversed(_tab4.index[1:]):
            if alt >= _tab4["H_b"][i]:
                integr = p_barom(h=alt,
                                 p_0=_tab5["p"][i] * pfact,
                                 h_0=_tab4["H_b"][i],
                                 T_0=_tab4["T_b"][i],
                                 beta=_tab4["beta"][i],
                                 gpm=True,
                                 Kelvin=True,
                                 hPa=False,
                                 )  # Pa
                break
        else:
            return np.nan
    if hPa is True:
        p = integr / 100.
    else:
        p = integr
    return p

# ---------------------------------------------------------------------


def rho_iso(h, p_0=pzero_iso, gpm=False):
    """
    density of dry air in the standard atmosphere at given height

    :param h: altitude (float) in m or gpm, depending on `gpm`
    :param p_0: (float, optional) air pressure at base level
        (in hPa or Pa, depending on `hPa`).
        Defaults to `pzero_iso`
    :param gpm: (bool, optional)
        if ``True``, h is given in gepotential meters,
        if ``False`` h is given in gemetric height. Default to ``False``.

    :return: density in :math:`kg m^{-3}`
    """
    return gas_rho(p_iso(h, p_0=p_0, gpm=gpm, hPa=False),
                   T_iso(h, gpm=gpm, Kelvin=True),
                   q=0.,
                   Kelvin=True,
                   hPa=False)

# ---------------------------------------------------------------------


def altitude(p, p_0=None, gpm=False, hPa=False):
    """
    returns altitude (height above se level) as funtion of pressure

    :param p: (float) air pressure (in hPa or Pa, depending on `hPa`).
    :param p_0: (float, optional) air pressure at base level
        (in hPa or Pa, depending on `hPa`).
        Defaults to `pzero_iso`
    :param gpm: (bool, optional) if ``True``,
        h is returned in gepotential meters,
        if ``False`` h is returned in gemetric height. Default to False.
    :param hPa: (optional, optional)
      if ``True``,  `p` must be supplied in hPa.
      If ``False``, `p` must be supplied in Pa.
      Defaults to ``False``.
    :return: altitude (float) in m or gpm, depending on `gpm`
    """
    if p_0 is None:
        pp_0 = pzero_iso             # Pa
    else:
        pp_0 = _to_Pa(p_0, hPa=hPa)  # Pa
    pp = _to_Pa(p, hPa=hPa)  # Pa
    if pp > _tab5["p"][0]:
        H = np.nan
    else:
        for i in reversed(_tab5.index):
            pi = _tab5["p"][i] * pp_0 / pzero_iso  # Pa
            if pp <= pi:
                beta = _tab4["beta"][i]  # K/m
                p_0 = pi      # Pa
                H_0 = _tab4["H_b"][i]    # gpm
                T_0 = _tab4["T_b"][i]    # K
                if abs(beta) > 1.E-6:
                    pwr = - gphi() / (R_iso * beta)
                    H = (((pp / p_0) ** (1/pwr) - 1) * (T_0 / beta) +
                         H_0)  # gpm
                else:
                    H = (-np.log(pp / p_0) * R_iso * T_0 / gphi() +
                         H_0)  # gpm
                break
    if gpm:
        res = H  # gpm
    else:
        res = h_geomet(H)  # m
    return res


if __name__ == "__main__":
    print(p_iso(0))
