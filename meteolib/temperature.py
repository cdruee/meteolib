#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Temperature calculations an conversions
'''

import numpy as np

from ._utils import _check
from .constants import R, Tzero, cp
from .pressure import _to_Pa

"""Reference pressure for calculation of potential temperatures"""
POTENTIAL_REFERENCE_PRESSURE = 100000.  # Pa

# ---------------------------------------------------------------------


def _to_K(t, Kelvin=None):
    """
    Returns the given Temperature in Kelvin

    :param t: Temperature
    :param Kelvin: (optional) if ``False``, unit of t is assumed to
      be Kelvin. If ``False``, unit of t is assumed to be Celsius.
      If missing of ``None``, unit of t is autodetected. Defaults to
      ``None``.
    :return: temperature in K. If t is not a valid number, ``np.nan``
      is returned.
    :rtype: float
    :raises ValueError: if Kelvin is not ``True``, ``False``, or ``None``

    """
    try:
        _check('t', t, 'float', nan=True)
    except ValueError:
        return np.nan
    if Kelvin is None:
        # autodetect
        if np.all(np.less(t, 150.)):
            t = t + Tzero
    elif Kelvin is False:
        t = t + Tzero
    elif Kelvin is True:
        pass
    else:
        raise ValueError('invalid value for "Kelvin": %s' % format(Kelvin))
    return t

# ---------------------------------------------------------------------


def _to_C(t, Kelvin=None):
    '''
    Returns the given Temperature in Celsius

    :param t: Temperature
    :param Kelvin: (optional) if ``False``, unit of t is assumed to
      be Kelvin. If ``False``, unit of t is assumed to be Celsius.
      If missing of ``None``, unit of t is autodetected. Defaults to
      ``None``.
    :return: temperature in °C. If t is not a valid number, ``np.nan``
      is returned.
    :rtype: float
    :raises ValueError: if Kelvin is not ``True``, ``False``, or ``None``

    '''
    try:
        _check('t', t, 'float', nan=True)
    except ValueError:
        return np.nan
    if Kelvin is None:
        # autodetect
        if t >= 150.:
            t = t - Tzero
    elif Kelvin is False:
        pass
    elif Kelvin is True:
        t = t - Tzero
    else:
        raise ValueError('invalid value for "Kelvin": %s' % format(Kelvin))
    return t

# ---------------------------------------------------------------------


def KtoC(t):
    '''
    Converts Temperature from Kelvin to Celsius

    :param t: Temperature K
    :return: temperature in °C. If t is not a valid number, ``np.nan``
      is returned.
    :rtype: float

    '''
    return _to_C(t, Kelvin=True)

# ---------------------------------------------------------------------


def CtoK(t):
    '''
    Converts Temperature from Celsius to Kelvin

    :param t: Temperature in °C
    :return: temperature in K. If t is not a valid number, ``np.nan``
      is returned.
    :rtype: float

    '''
    return _to_K(t, Kelvin=False)

# ---------------------------------------------------------------------


def FtoC(t):
    '''
    Converts Temperature from Fahrenheit to Celsius

    :param t: Temperature in °F (float)
    :return: temperature in °C.
    :rtype: float

    Conversion:

    .. math:: T_C = (T_F + 32.0) \\cdot \\frac{5}{9}

    '''
    try:
        _check('t', t, 'float')
    except ValueError:
        return np.nan
    return (t - 32.0) * 5./9.


# ---------------------------------------------------------------------
def CtoF(t):
    '''
    Converts Temperature from Celsius to Fahrenheit

    :param t: Temperature in °C (float)
    :return: temperature in °F.
    :rtype: float

    Conversion:

    .. math:: T_F = T_C \\cdot 1.80 + 32

    '''
    try:
        _check('t', t, 'float')
    except ValueError:
        return np.nan
    return t*1.80 + 32


# ---------------------------------------------------------------------
def Tpot(t, p, Kelvin=None, hPa=False, pref=POTENTIAL_REFERENCE_PRESSURE):
    '''
    Converts Temperature from Celsius to Kelvin

    :param t: Temperature in °C
    :param p: (optional) air pressure (float).
    :param Kelvin: (optional) if ``False``, unit of t is assumed to
      be Kelvin. If ``False``, unit of t is assumed to be Celsius.
      If missing of ``None``, unit of t is autodetected. Defaults to
      ``None``.
    :param hPa: (optional)
      if ``False``, `p` must be supplied in Pa .
      If ``True``, `p` must be supplied in hPa.
      Defaults to ``False``
    :param pref: (optional)
      Defaults to 100 000 Pa.
    :return: temperature in K. If t is not a valid number, ``np.nan``
      is returned.
    :rtype: float

    '''
    T = _to_K(_check('t', t, 'float', nan=True), Kelvin=Kelvin)
    pp = _to_Pa(_check('p', p, 'float', gt=0., nan=True), hPa=hPa)

    Theta = T * (pref / pp)**(R/cp)

    if Kelvin is False:
        return _to_C(Theta, Kelvin=True)
    return _to_K(Theta, Kelvin=True)

# ---------------------------------------------------------------------


def inv_Tpot(t, p, Kelvin=None, hPa=False, pref=POTENTIAL_REFERENCE_PRESSURE):
    '''
    Converts Temperature from Celsius to Kelvin

    :param t: Temperature in °C
    :param p: (optional) air pressure (float).
    :param Kelvin: (optional) if ``False``, unit of t is assumed to
      be Kelvin. If ``False``, unit of t is assumed to be Celsius.
      If missing of ``None``, unit of t is autodetected. Defaults to
      ``None``.
    :param hPa: (optional)
      if ``False``, `p` must be supplied in Pa .
      If ``True``, `p` must be supplied in hPa.
      Defaults to ``False``
    :param pref: (optional)
      Defaults to 100 000 Pa.
    :return: temperature in K. If t is not a valid number, ``np.nan``
      is returned.
    :rtype: float

    '''
    Theta = _to_K(_check('t', t, 'float', nan=True), Kelvin=Kelvin)
    pp = _to_Pa(_check('p', p, 'float', gt=0., nan=True), hPa=hPa)

    T = Theta / (pref / pp)**(R/cp)

    if Kelvin is False:
        return _to_C(T, Kelvin=True)
    return _to_K(T, Kelvin=True)
