#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Temperature calculations an conversions
'''
import numpy as np

import pandas as pd


# ---------------------------------------------------------------------
def _to_Pa(p, hPa=None):
    '''
    Returns the given air pressure in Pascal

    :param p: air pressure
    :param nPa: (optional) if ``False``, unit of p is assumed to
      be hPa. If ``False``, unit of p is assumed to be Pa.
      If missing or ``None``, unit of p is autodetected.
      Defaults to ``None``.
    :return: pressure in Pa. If p is not a valid number, ``np.nan``
      is returned. If p is ``None``, ``None`` is returned.
    :raises ValueError: if hPa is not ``True``, ``False``, or ``None``


    Usage::

      >>> _to_Pa(1010.)
      101000.

    '''
    if p is None:
        return None
    if pd.isnull(p):
        return np.nan
    if hPa is None:
        # autodetect
        if p < 1500.:
            p = p * 100.
    elif hPa is False:
        pass
    elif hPa is True:
        p = p * 100.
    else:
        raise ValueError('invalid value passed for "hPa": {}'.format(hPa))
    return p

# ---------------------------------------------------------------------


def _to_hPa(p, hPa=None):
    '''
    Returns the given air pressure in Hectopascal

    :param p: air pressure
    :param nPa: (optional) if ``False``, unit of p is assumed to
      be hPa. If ``False``, unit of p is assumed to be Pa.
      If missing or ``None``, unit of p is autodetected.
      Defaults to ``None``.
    :return: pressure in Pa. If p is not a valid number, ``np.nan``
      is returned. If p is ``None``, ``None`` is returned.
    :raises ValueError: if hPa is not ``True``, ``False``, or ``None``


    Usage::

      >>> _to_hPa(10100.)
      101.

    '''
    if p is None:
        return None
    if pd.isnull(p):
        return np.nan
    if hPa is None:
        # autodetect
        if p >= 1500.:
            p = p / 100.
    elif hPa is False:
        p = p / 100.
    elif hPa is True:
        pass
    else:
        raise ValueError('invalid value passed for "hPa": {}'.format(hPa))
    return p

# ----------------------------------------------------------------------


def pa2mmhg(pa: float) -> float:
    """
    Converts air pressure in Pascal to mmHg
    :param pa: pressure in Pascal
    :type pa: float
    :return: air pressure in mmHg
    :rtype: float
    """
    return pa / 133.322

# ----------------------------------------------------------------------


def mmhg2pa(mmhg: float) -> float:
    """
    Converts air pressure in mmHg to Pascal
    :param mmhg: air pressure in mmHg
    :type mmhg: float
    :return: air pressure in Pascal
    :rtype: float
    """
    return mmhg * 133.322
