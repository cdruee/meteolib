#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Evaporation

"""
import logging

import numpy as np

import pandas as pd

from .constants import Tzero, sigma
from .humidity import esat_w
from .pressure import _to_hPa
from .temperature import _to_C


def penmanmonteith_mono(T, rh, vh, p, Kdn, time=None, df=None,
                        lat=45., ele=0.,
                        Kelvin=None, hPa=False, percent=False):
    r'''
    Calculate reference evaporation after Penman-Monteith.

    Monolithic code following the recipe of [Zot2018]_
    to implement the method as described in [FAO56]_ .

    :param T:  2-m temperature (series) in :math:`K` or :math:`^\circ C`,
      depending on `Kelvin` .
    :param rh: realtive humidity (series) in percent or fractions of unity,
      depending on `Kelvin` .
    :param vh: horizontal wind speed in :math:`m\,s^{-1}`.
    :param p: atmospheric pressure in hPa or Pa, depending on `hPa`.
    :param Kdn: short-wave downward radiation in :math:`W\,m^{-2}`.
    :param lat: Latitude (in degrees) of the location for which
      evaportaion is calculated. Defaults to 45°(N).
    :param ele: Elevation (above sea level) of the location for which
      evaportaion is calculated. Defaults to 0m.
    :param Kelvin: (optional) if ``False``, all temperatures are assumed to
      be Kelvin. If ``False``, all temperatures are assumed to be Celsius.
      If missing of ``None``, unit  temperatures are autodetected. Defaults to
      ``None``.
    :param Kelvin: (optional) if ``False``, all temperatures are assumed to
      be Kelvin. If ``False``, all temperatures are assumed to be Celsius.
      If missing of ``None``, unit  temperatures are autodetected. Defaults to
      ``None``.
    :param hPa: (optional)
      if ``True``,  `p`, `e`, or `ew` must be supplied in hPa.
      If ``False``, `p`, `e`, or `ew` must be supplied in Pa.
      Defaults to ``False``.
    :param percent: (optional)
      if ``True``, relative humidity is returned in percent.
      If ``False``, relative humidity is returned unitless (fraction of unity).
      Defaults to ``False``.

    Parameters `T`, `rh`, `vh`, `p`, snd `Kdn`  may me supplied
    all as column names of a pandas.DataFrame (requires `df`),
    all as ``pandas.Series`` (indices must be identical), or
    all as ``list`` (requires `time` as ``list``).

    :return: reference evaporation :math:`ET0` (in mm/day).
    :rtype: pandas.Series

    '''
    if all(isinstance(x, str) for x in [T, rh, vh, p, Kdn]):
        if isinstance(df, pd.DataFrame):
            # get columns from df
            o = df[[T, rh, vh, p, Kdn]]
            # rename columns
            o.columns = ['T', 'rh', 'vh', 'p', 'Kdn']
        else:
            raise ValueError(
                'string arguments require a pandas.DataFrame as df')
    elif all(isinstance(x, pd.Series) for x in [T, rh, vh, p, Kdn]):
        o = pd.DataFrame({'T': T,
                         'rh': rh,
                          'vh': vh,
                          'p': p,
                          'Kdn': Kdn
                          })
    elif all(isinstance(x, list) for x in [T, rh, vh, p, Kdn, time]):
        if time is None:
            raise ValueError(
                'list arguments require list of observation times as time')
        o = pd.DataFrame({'T': T,
                         'rh': rh,
                          'vh': vh,
                          'p': p,
                          'Kdn': Kdn,
                          },
                         index=pd.DatetimeIndex(time)
                         )
    else:
        raise ValueError('illegal or mixed types of arguments')
    #
    # convert units if needed
    #
    o['T'] = o['T'].apply(_to_C, Kelvin=Kelvin)  # -> C
    o['p'] = o['p'].apply(_to_hPa, hPa=hPa)/10.  # -> kPa
    if percent is False:
        o['rh'] = o['rh'] * 100.0  # 1 -> %
    elif percent is True:
        pass                              # % -> %
    else:
        raise ValueError('percent must be either True or False')

    # resample daily min/max values
    cur = pd.DataFrame({
        'Tmax': o['T'].resample('1D').max(),
        'Tmin': o['T'].resample('1D').min(),
        'rhmax': o['rh'].resample('1D').max(),
        'rhmin': o['rh'].resample('1D').min(),
    })
    #
    logging.debug('calculate')
    # 1.Mean daily temperature
    cur['Tmean'] = (cur['Tmax'] + cur['Tmin'])/2.
    # 2. Mean daily solar radiation
    cur['Rs'] = (
        o['Kdn'].resample('1D').mean() *  # W/m² -> J m^-2 s^-1
        86400. / 1.E6  # * 86400 / 10^6 -> MJ m^-2 day^-1
    )
    # 3. Mean daily solar radiation
    cur['umean'] = o['vh'].resample('1D').mean()
    # 4. Slope of saturation vapor pressure curve
    cur['slope'] = (4098 *
                    (0.6108 *
                     np.exp((17.27*cur['Tmean']) / (237.3+cur['Tmean']))) /
                    ((237.3+cur['Tmean'])**2)
                    )
    # 5. Atmospheric Pressure (P)
    cur['pmean'] = o['p'].resample('1D').mean()/10.
    # 6. Psychrometric constant
    cur['gamma'] = 0.000665*cur['pmean']
    # 7. Delta Term (DT) (auxiliary calculation for Radiation Term)
    cur['delta'] = cur['slope'] / \
        (cur['slope']+cur['gamma']*(1+0.34*cur['umean']))
    # 8. Psi Term (PT) (auxiliary calculation for Wind Term)
    cur['psi'] = cur['gamma'] / (cur['slope'] +
                                 cur['gamma'] * (1 + 0.34 * cur['umean']))
    # 9. Temperature Term (TT) (auxiliary calculation for Wind Term)
    cur['TT'] = (900/(273+cur['Tmean']) * cur['umean'])
    # 10. Mean saturation vapor pressure derived from air temperature
    cur['e_Tmax'] = 0.6108 * np.exp((17.27 * cur['Tmax']) /
                                    (237.3 + cur['Tmax']))
    cur['e_Tmin'] = 0.6108 * np.exp((17.27 * cur['Tmin']) /
                                    (237.3 + cur['Tmin']))
    cur['es'] = (cur['e_Tmax']+cur['e_Tmin'])/2.
    # 11. Actual vapor pressure (ea) derived from relativehumidity
    cur['ea'] = (cur['e_Tmax'] * cur['rhmax'] / 100.0 +
                 cur['e_Tmin'] * cur['rhmin'] / 100.0) / 2.
    # 12. The inverse relative distance Earth-Sun (dr) and solar declination
    cur['jul'] = [int(x) for x in cur.index.strftime('%j')]
    cur['dr'] = 1.0 + 0.033 * np.cos(cur['jul'] * 2. * np.pi / 365.0)
    cur['del'] = 0.409 * np.sin(cur['jul'] * 2. * np.pi / 365.0 - 1.39)
    # 13. Conversion of latitude in degrees to radians
    cur['phi'] = [lat * np.pi / 180.] * len(cur.index)
    # 14. Sunset hour angle
    cur['omega_s'] = np.arccos(-np.tan(cur['phi']) * np.tan(cur['del']))
    # 15. Extraterrestrial radiation (Ra)
    Gs = 0.0820
    cur['Ra'] = (24.*60./np.pi) * Gs * cur['dr'] * (
        cur['omega_s'] * np.sin(cur['phi']) * np.sin(cur['del']) +
        np.cos(cur['phi']) * np.cos(cur['del']) * np.sin(cur['omega_s'])
    )
    # 16. Clear sky solar radiation (Rso)
    cur['Rso'] = (0.75 + 2.0E-5*ele)*cur['Ra']
    # 17. Net solar or net shortwave radiation (Rns)
    albedo = 0.23  # for the hypothetical grass reference crop
    cur['Rns'] = (1.0 - albedo) * cur['Rs']
    # 18.Net outgoing long wave solar radiation (Rnl)
    # Stefan-Boltzmann constant
    sigma = 4.903E-9  # MJ K^-4 m^-2 day^-1
    cur['Rnl'] = sigma * (
        ((cur['Tmax'] + 273.16)**4 + (cur['Tmin'] + 273.16)**4) / 2.
    ) * (
        0.34 - 0.14*np.sqrt(cur['ea'])
    ) * (
        1.35*(cur['Rs']/cur['Rso']) - 0.35
    )
    # 19. Net radiation (Rn)
    cur['Rn'] = cur['Rns'] - cur['Rnl']
    # net radiation (Rn) in equivalent of evaporation (mm) (R ng)
    cur['Rng'] = 0.408 * cur['Rn']
    # FS1: Radiation term (ETrad)
    cur['ETrad'] = cur['delta'] * cur['Rng']
    # FS2: Wind term (ETwind)
    cur['ETwind'] = cur['delta'] * cur['TT'] * (cur['es']-cur['ea'])
    # Final Reference Evapotranspiration Value (ETo)
    cur['ET0'] = cur['ETwind'] + cur['ETrad']
    # in mm day^-1

    return cur['ET0']


def penmanmonteith(T, rh, vh, p, Kdn, time=None, df=None,
                   lat=45., ele=0.,
                   Kelvin=None, hPa=False, percent=False):
    r'''
    calculate reference evaporation after Penman-Monteith

    Code using otherlibrary functions, based on the steps
    described by [Zot2018]_  to implement the method as
    described in [FAO56]_ .

    :param T:  2-m temperature (series) in :math:`K` or :math:`^\circ C`,
      depending on `Kelvin` .
    :param rh: realtive humidity (series) in percent or fractions of unity,
      depending on `Kelvin` .
    :param vh: horizontal wind speed in :math:`m\,s^{-1}`.
    :param p: atmospheric pressure in hPa or Pa, depending on `hPa`.
    :param Kdn: short-wave downward radiation in :math:`W\,m^{-2}`.
    :param Kelvin: (optional) if ``False``, all temperatures are assumed to
      be Kelvin. If ``False``, all temperatures are assumed to be Celsius.
      If missing of ``None``, unit  temperatures are autodetected. Defaults to
      ``None``.
    :param Kelvin: (optional) if ``False``, all temperatures are assumed to
      be Kelvin. If ``False``, all temperatures are assumed to be Celsius.
      If missing of ``None``, unit  temperatures are autodetected. Defaults to
      ``None``.
    :param hPa: (optional)
      if ``True``,  `p`, `e`, or `ew` must be supplied in hPa.
      If ``False``, `p`, `e`, or `ew` must be supplied in Pa.
      Defaults to ``False``.
    :param percent: (optional)
      if ``True``, relative humidity is returned in percent.
      If ``False``, relative humidity is returned unitless (fraction of unity).
      Defaults to ``False``.

    Parameters `T`, `rh`, `vh`, `p`, snd `Kdn`  may me supplied
    all as column names of a pandas.DataFrame (requires `df`),
    all as ``pandas.Series`` (indices must be identical), or
    all as ``list`` (requires `time` as ``list``).

    :return: reference evaporation :math:`ET0` (in mm/day).
    :rtype: pandas.Series

    '''
    if all(isinstance(x, str) for x in [T, rh, vh, p, Kdn]):
        if isinstance(df, pd.DataFrame):
            # get columns from df
            o = df[[T, rh, vh, p, Kdn]]
            # rename columns
            o.columns = ['T', 'rh', 'vh', 'p', 'Kdn']
        else:
            raise ValueError(
                'list arguments require list of observation times as time')
    elif all(isinstance(x, pd.Series) for x in [T, rh, vh, p, Kdn]):
        o = pd.DataFrame({'T': T,
                         'rh': rh,
                          'vh': vh,
                          'p': p,
                          'Kdn': Kdn
                          })
    elif all(isinstance(x, list) for x in [T, rh, vh, p, Kdn, time]):
        if time is None:
            raise ValueError('list arguments require times as index')
        o = pd.DataFrame({'T': T,
                         'rh': rh,
                          'vh': vh,
                          'p': p,
                          'Kdn': Kdn,
                          },
                         index=pd.DatetimeIndex(time)
                         )
    else:
        raise ValueError('illegal or mixed types of arguments')
    #
    # convert units if needed
    #
    o.loc[:, 'T'] = o['T'].apply(_to_C, Kelvin=Kelvin)  # -> C
    o.loc[:, 'p'] = o['p'].apply(_to_hPa, hPa=hPa)/10.  # -> kPa
    if percent is False:
        o.loc[:, 'rh'] = o['rh'] * 100.0  # 1 -> %
    elif percent is True:
        pass                              # % -> %
    else:
        raise ValueError('percent must be either True or False')

    # resample daily min/max values
    cur = pd.DataFrame({
        'Tmax': o['T'].resample('1D').max(),
        'Tmin': o['T'].resample('1D').min(),
        'rhmax': o['rh'].resample('1D').max(),
        'rhmin': o['rh'].resample('1D').min(),
    })
    #
    logging.debug('calculate')
    # 1.Mean daily temperature
    cur['Tmean'] = (cur['Tmax'] + cur['Tmin'])/2.
    # 2. Mean daily solar radiation
    cur['Rs'] = (
        o['Kdn'].resample('1D').mean() *  # W/m² -> J m^-2 s^-1
        86400. / 1.E6  # * 86400 / 10^6 -> MJ m^-2 day^-1
    )
    # 3. Mean horizontal wind at 2m
    cur['umean'] = o['vh'].resample('1D').mean()
    # 4. Slope of saturation vapor pressure curve
    cur['slope'] = [((esat_w(x+0.1)-esat_w(x-0.1))/0.2) /
                    1000. for x in cur['Tmean']]  # Pa/K -> kPa/K
    # 5. Atmospheric Pressure (P)
    cur['pmean'] = o['p'].resample('1D').mean()/10.
    # 6. Psychrometric constant
    cur['gamma'] = 0.000665*cur['pmean']
    # 7. Delta Term (DT) (auxiliary calculation for Radiation Term)
    cur['delta'] = cur['slope'] / \
        (cur['slope']+cur['gamma']*(1+0.34*cur['umean']))
    # 8. Psi Term (PT) (auxiliary calculation for Wind Term)
    cur['psi'] = cur['gamma']/(cur['slope']+cur['gamma']*(1+0.34*cur['umean']))
    # 9. Temperature Term (TT) (auxiliary calculation for Wind Term)
    cur['TT'] = (900/(Tzero+cur['Tmean']) * cur['umean'])
    # 10. Mean saturation vapor pressure derived from air temperature
    cur['e_Tmax'] = [esat_w(x) / 1000. for x in cur['Tmax']]  # Pa -> kPa
    cur['e_Tmin'] = [esat_w(x) / 1000. for x in cur['Tmin']]  # Pa -> kPa
    cur['es'] = (cur['e_Tmax']+cur['e_Tmin'])/2.
    # 11. Actual vapor pressure (ea) derived from relativehumidity
    cur['ea'] = (cur['e_Tmax']*cur['rhmax']/100.0 +
                 cur['e_Tmin']*cur['rhmin']/100.0) / 2.
    # 12. The inverse relative distance Earth-Sun (dr) and solar declination
    cur['jul'] = [int(x) for x in cur.index.strftime('%j')]
    cur['dr'] = 1.0 + 0.033 * np.cos(cur['jul']*2.*np.pi/365.0)
    cur['del'] = 0.409 * np.sin(cur['jul']*2.*np.pi/365.0 - 1.39)
    # 13. Conversion of latitude in degrees to radians
    cur['phi'] = [lat * np.pi / 180.]*len(cur.index)
    # 14. Sunset hour angle
    cur['omega_s'] = np.arccos(-np.tan(cur['phi'])*np.tan(cur['del']))
    # 15. Extraterrestrial radiation (Ra)
    Gs = 0.0820
    cur['Ra'] = (24.*60./np.pi) * Gs * cur['dr'] * (
        cur['omega_s'] * np.sin(cur['phi']) * np.sin(cur['del']) +
        np.cos(cur['phi']) * np.cos(cur['del']) * np.sin(cur['omega_s'])
    )
    # 16. Clear sky solar radiation (Rso)
    cur['Rso'] = (0.75 + 2.0E-5*ele)*cur['Ra']
    # 17. Net solar or net shortwave radiation (Rns)
    albedo = 0.23  # for the hypothetical grass reference crop
    cur['Rns'] = (1.0 - albedo) * cur['Rs']
    # 18.Net outgoing long wave solar radiation (Rnl)
    # Stefan-Boltzmann constant
    sigmad = sigma * 86400 / 1000000  # W m-2 K-4 ->  MJ K^-4 m^-2 day^-1
    cur['Rnl'] = sigmad * (
        ((cur['Tmax'] + Tzero)**4 + (cur['Tmin'] + Tzero)**4) / 2.
    ) * (
        0.34 - 0.14*np.sqrt(cur['ea'])
    ) * (
        1.35*(cur['Rs']/cur['Rso']) - 0.35
    )
    # 19. Net radiation (Rn)
    cur['Rn'] = cur['Rns'] - cur['Rnl']
    # net radiation (Rn) in equivalent of evaporation (mm) (R ng)
    cur['Rng'] = 0.408 * cur['Rn']
    # FS1: Radiation term (ETrad)
    cur['ETrad'] = cur['delta'] * cur['Rng']
    # FS2: Wind term (ETwind)
    cur['ETwind'] = cur['delta'] * cur['TT'] * (cur['es']-cur['ea'])
    # Final Reference Evapotranspiration Value (ETo)
    cur['ET0'] = cur['ETwind'] + cur['ETrad']
    # in mm day^-1

    return cur['ET0']
