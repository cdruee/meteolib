#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Humidity calculations and conversions
'''
import numpy as np

from ._utils import _check, _only
from .constants import Lv, Ttriple, cp, pzero
from .pressure import _to_Pa, _to_hPa
from .temperature import Tpot, _to_C, _to_K
from .thermodyn import gas_p, gas_rho

_esat = 'magnus'

# ---------------------------------------------------------------------


def _f(p):
    '''
    Moist air correction for water vapor saturation pressure
    :param p: air pressure in hPa (float)
    :return: correction factor. If `p` is ``None``, the value 1. is returned.

    This form of the equation was adopted by Forty-second session of the
    Executive Council (EC-XLII) in 1990 [WMO8]_:
    .. math::
      f (p) = 1.001 6 + 3.15 · 10 –6 p – 0.074 p –1

    Usage::

      >>> _f(1013.)
      1.001

    '''
    if p is None:
        f = 1.
    else:
        f = 1.0016 + 3.15E-6 * p - 0.074 / p
    return (f)

# ---------------------------------------------------------------------


def magnus_w(t, Kelvin=None, p=None, hPa=False):
    '''
    Magnus Formula, returns the saturation pressure for water vapor
    over an even liquid water surface in hPa

    :param t: Temperature
    :param Kelvin: (optional)
      if ``True``, unit of `t` is assumed to be Kelvin.
      If ``False``, unit of `t` is assumed to be Celsius.
      If missing of ``None``, unit of `t` is autodetected. Defaults to
      ``None``.
    :param p: (optional) air pressure (float). If `p` is ``None`` or missing,
      the correction for moist air is omittet, i.e. :math:`f(p) = 1`.
      Defaults to ``None``.
    :param hPa: (optional)
      if ``False``, `p` must be supplied and `magnus_w is returned in Pa .
      If ``True``, `p` must be supplied and `magnus_w is returned in hPa.
      If missing or ``None`` the unit of `p` is autodetected
      and `magnus_w is returned in Pa. Defaults to ``None``.
    :return: saturation pressure in hPa. If `t` is not a valid number,
      ``np.nan`` is returned.
    :rtype: float
    :raises ValueError: if `Kelvin` is not ``True``, ``False``, ``None``,
      or missing
    :raises ValueError: if `hPa` is not ``True``, ``False``, ``None``,
      or missing

    This form of the equation was adopted by Forty-second session of the
    Executive Council (EC-XLII) in 1990 [WMO8]_:

    .. math::
      e = f(p) * 6.112 * np.exp((17.62 * t) / (243.12 + t))
      :label: magnus

    Results are valid in the range water -45°C to +60°C for pure phase water.

    '''
    t = _to_C(t, Kelvin)  # K
    if np.isnan(t):
        e = np.nan
    else:
        e = 6.112 * np.exp((17.62 * t) / (243.12 + t))  # hPa
    if p is not None:
        e = _f(_to_hPa(p, hPa)) * e
    if hPa is True:
        return e  # hPa
    else:
        return 100. * e  # Pa

# ---------------------------------------------------------------------


def magnus_i(t, Kelvin=None, p=None, hPa=False):
    '''
    Magnus Formula, returns the saturation pressure for water vapor
    over an even liquid water surface in hPa

    :param t: Temperature
    :param Kelvin: (optional)
      if ``True``, unit of  `t` is assumed to be Kelvin.
      If ``False``, unit of  `t` is assumed to be Celsius.
      If missing or ``None``, unit of  `t` is autodetected. Defaults to
      ``None``.
    :param p: (optional) air pressure (float). If `p` is ``None`` or missing,
      the correction for moist air is omittet, i.e. :math:`f(p) = 1`.
      Defaults to ``None``.
    :param hPa: (optional)
      if ``True``, `p` must be supplied and `magnus_i` is returned in hPa.
      If ``False``, `p` must be supplied and `magnus_i` is returned in Pa .
      If missing or ``None`` the unit of `p` is autodetected
      and `magnus_i` is returned in Pa. Defaults to ``None``.
    :return: saturation pressure in hPa. If `t` is not a valid number,
      ``np.nan`` is returned.
    :rtype: float
    :raises ValueError: if `Kelvin` is not ``True``, ``False``, ``None``,
      or missing
    :raises ValueError: if `hPa` is not ``True``, ``False``, ``None``,
      or missing

    This form of the equation was adopted by Forty-second session of the
    Executive Council (EC-XLII) in 1990 [WMO8]_:

    .. math::
      e = 6.112 * np.exp((22.46 * t) / (272.62 + t))
      :label: magnice

    Results are valid in the range water -65°C to 0°C for pure phase water.

    '''
    t = _to_C(t, Kelvin)  # K
    if np.isnan(t):
        e = np.nan
    else:
        e = 6.112 * np.exp((22.46 * t) / (272.62 + t))  # hPa
    if p is not None:
        e = _f(_to_hPa(p, hPa)) * e
    if hPa is True:
        return e  # hPa
    else:
        return 100. * e  # Pa

# ---------------------------------------------------------------------


def tetens_w(t, Kelvin=None, p=None, hPa=False):
    '''
    Tetens formula, identical to :py:meth:`meltib.humidity.magnus_w`, except
    Coefficients provided by [Tet1930]_, cited from [Mon2008]_:

    .. math::


      e = 6.1078 * np.exp((17.27 * t) / (237.3 + t))

    Results are within 1 Pa of exact values in the range -0°C to +35°C
    [Mon2008]_.

    '''
    t = _to_C(t, Kelvin)  # K
    if np.isnan(t):
        e = np.nan
    else:
        e = 6.1078 * np.exp((17.27 * t) / (237.3 + t))  # hPa
    if p is not None:
        e = _f(_to_hPa(p, hPa)) * e
    if hPa is True:
        return e  # hPa
    else:
        return 100. * e  # Pa

# ---------------------------------------------------------------------


def tetens_i(t, Kelvin=None, p=None, hPa=False):
    '''
    Tetens formula, identical to :py:meth:`meltib.humidity.magnus_i`, except
    Coefficients provided by [Tet1930]_, cited from [Mon2008]_:

    .. math::


      e = 6.1078 * np.exp((21.875 * t) / (265.5 + t))

    Results are within 1 Pa of exact values in the range -20°C to +0°C
    [Mon2008]_.

    '''
    t = _to_C(t, Kelvin)  # K
    if np.isnan(t):
        e = np.nan
    else:
        e = 6.1078 * np.exp((21.875 * t) / (265.5 + t))  # hPa
    if p is not None:
        e = _f(_to_hPa(p, hPa)) * e
    if hPa is True:
        return e  # hPa
    else:
        return 100. * e  # Pa

# ---------------------------------------------------------------------


def goff_gratch_w(t, Kelvin=None, p=None, hPa=False):
    '''
    Goff-Gratch formula, identical to :py:meth:`meltib.humidity.magnus_w`,
    except formula provided by [Gof1957]_, cited from [WMO49]_:

    .. math::

      log10(ew) = & 10.79574 * (1 – Ttriple/T) \\\\
                  &– 5.02800 * np.log10 (T/Ttriple) \\\\
                  &+ 1.50475E–4 * (1 – 10**(–8.2969*(T/Ttriple–1)) ) \\\\
                  &+ 0.42873E–3 * ( 10**(4.76955*(1-Ttriple/T)) – 1) \\\\
                  &+ 0.78614

    '''
    T = _to_K(t, Kelvin)  # K
    if np.isnan(T):
        e = np.nan
    else:
        loge = (10.79574 * (1 - Ttriple/T)
                - 5.02800 * np.log10(T/Ttriple)
                + 1.50475E-4 * (1 - 10**(-8.2969*(T/Ttriple-1)))
                + 0.42873E-3 * (10**(4.76955*(1-Ttriple/T)) - 1)
                + 0.78614
                )
        e = 10**loge
    if p is not None:
        e = _f(_to_hPa(p, hPa)) * e
    if hPa is True:
        return e  # hPa
    else:
        return 100. * e  # Pa

# ---------------------------------------------------------------------


def goff_gratch_i(t, Kelvin=None, p=None, hPa=False):
    '''
    Goff-Gratch formula, identical to :py:meth:`meltib.humidity.magnus_i`,
    except formula provided by [Gof1957]_, cited from [WMO49]_:

    .. math::

      log10(ei) = & -9.09718 * (Ttriple/T - 1) \\\\
                  &- 3.56654 * np.log10(Ttriple/ T) \\\\
                  &+ 0.876793 * (1 - T/Ttriple) \\\\
                  &+ np.log10(6.1071) \\\\
                  )

    '''
    T = _to_K(t, Kelvin)  # K
    if np.isnan(t):
        e = np.nan
    else:
        loge = (-9.09718 * (Ttriple/T - 1)
                - 3.56654 * np.log10(Ttriple / T)
                + 0.876793 * (1 - T/Ttriple)
                + np.log10(6.1071)
                )
        e = 10**loge
    if p is not None:
        e = _f(_to_hPa(p, hPa)) * e
    if hPa is True:
        return e  # hPa
    else:
        return 100. * e  # Pa

# ---------------------------------------------------------------------


def hyland_wexler_w(t, Kelvin=None, p=None, hPa=False):
    '''
    Formula by Hyland and Wexler,
    identical to :py:meth:`meltib.humidity.magnus_w`,
    except formula provided by [HyW1983]_:

    .. math::

      ln(ew)    =  & - 0.58002206E4 / T \\\\
                   & + 0.13914993E1 \\\\
                   & - 0.48640239E-1 * T \\\\
                   & + 0.41764768E-4 * T^2 \\\\
                   & - 0.14452093E-7 * T^3 \\\\
                   & + 0.65459673E1 * ln(T)
    '''
    T = _to_K(t, Kelvin)  # K
    if np.isnan(T):
        e = np.nan
    else:
        lne = (- 0.58002206E4 / T
               + 0.13914993E1
               - 0.48640239E-1 * T
               + 0.41764768E-4 * T**2
               - 0.14452093E-7 * T**3
               + 0.65459673E1 * np.log(T)
               )
        e = np.exp(lne)  # Pa
    if p is not None:
        e = _f(_to_hPa(p, hPa)) * e
    if hPa is True:
        return e / 100.  # hPa
    else:
        return e  # Pa

# ---------------------------------------------------------------------


def hyland_wexler_i(t, Kelvin=None, p=None, hPa=False):
    '''
    Formula by Hyland and Wexler,
    identical to :py:meth:`meltib.humidity.magnus_i`,
    except formula provided by [HyW1983]_:

    .. math::

      ln(ei)    = &- 0.56745359E4 / T \\\\
                  &+ 0.63925247E1 \\\\
                  &- 0.96778430E-2 * T \\\\
                  &+ 0.62215701E-6 * T^2 \\\\
                  &+ 0.20747825E-8 * T^3 \\\\
                  &- 0.94840240E-12 * T^4 \\\\
                  &+ 0.41635019E1 * ln(T) \\\\

    '''
    T = _to_K(t, Kelvin)  # K
    if np.isnan(t):
        e = np.nan
    else:
        lne = (- 0.56745359E4 / T
               + 0.63925247E1
               - 0.96778430E-2 * T
               + 0.62215701E-6 * T**2
               + 0.20747825E-8 * T**3
               - 0.94840240E-12 * T**4
               + 0.41635019E1 * np.log(T)
               )
        e = np.exp(lne)  # Pa
    if p is not None:
        e = _f(_to_hPa(p, hPa)) * e
    if hPa is True:
        return e / 100.  # hPa
    else:
        return e  # Pa

# ---------------------------------------------------------------------


def sonntag_w(t, Kelvin=None, p=None, hPa=False):
    '''
    Formula by Sonntag, identical to :py:meth:`meltib.humidity.magnus_w`,
    except formula provided by [Son1990]_
    (and [Son1994]_,but with fewer digits):

    .. math::

      ln(ew)    =  & - 6096.9385 / T \\\\
                    + 16.635794 \\\\
                    - 2.711193E-2 * T \\\\
                    + 1.673952E-5 * T**2 \\\\
                    + 2.433502 * ln(T) \\\\
    '''
    T = _to_K(t, Kelvin)  # K
    if np.isnan(T):
        e = np.nan
    else:
        lne = (- 6096.9385 / T
               + 16.635794
               - 2.711193E-2 * T
               + 1.673952E-5 * T**2
               + 2.433502 * np.log(T)
               )
        e = np.exp(lne)  # hPa
    if p is not None:
        e = _f(_to_hPa(p, hPa)) * e
    if hPa is True:
        return e  # hPa
    else:
        return 100. * e  # Pa

# ---------------------------------------------------------------------


def sonntag_i(t, Kelvin=None, p=None, hPa=False):
    '''
    Formula by Sonntag, identical to :py:meth:`meltib.humidity.magnus_i`,
    except formula provided by [Son1990]_
    (and [Son1994]_, but with fewer digits):

    .. math::

      ln(ei)    = &- 6024.5282 / T \\\\
                  &+ 24.721994 \\\\
                  &+ 1.0613868E-2 * T \\\\
                  &- 1.3198825E-5 * T^2 \\\\
                  &- 0.49382577 * ln(T) \\\\

    '''
    T = _to_K(t, Kelvin)  # K
    if np.isnan(t):
        e = np.nan
    else:
        lne = (- 6024.5282 / T
               + 24.721994
               + 1.0613868E-2 * T
               - 1.3198825E-5 * T**2
               - 0.49382577 * np.log(T)
               )
        e = np.exp(lne)  # hPa
    if p is not None:
        e = _f(_to_hPa(p, hPa)) * e
    if hPa is True:
        return e  # hPa
    else:
        return 100. * e  # Pa


# ---------------------------------------------------------------------
def iapws_w(t, Kelvin=None, p=None, hPa=False):
    '''
    International Association for the Properties of Water and Steam (IAPWS)
    1995 Formula for water vapor,
    published by Wagner and Pruss [WaP2002]_ eqn. 2,5a:

    .. math::

      Log (ew/22.064E6) = 647.096/T * ((-7.85951783 v
          + 1.84408259 * v**1.5
          - 11.7866497 * v**3
          + 22.6807411 * v**3.5
          - 15.9618719 * v**4
          + 1.80122502 * v**7.5))


    with T in [K] and ew in [Pa] and v = 1 - T/647.096
    '''
    T = _to_K(t, Kelvin)  # K
    if np.isnan(T):
        e = np.nan
    else:
        v = 1. - T / 647.096
        lne = 647.096/T * (
            - 7.85951783 * v
            + 1.84408259 * v**1.5
            - 11.7866497 * v**3
            + 22.6807411 * v**3.5
            - 15.9618719 * v**4
            + 1.80122502 * v**7.5
        )
        e = np.exp(lne) * 22.064E4  # hPa
    if p is not None:
        e = _f(_to_hPa(p, hPa)) * e
    if hPa is True:
        return e  # hPa
    else:
        return 100. * e  # Pa

# ---------------------------------------------------------------------


def iapws_i(t, Kelvin=None, p=None, hPa=False):
    '''
    Sublimation pressure equation by Wagner, as published by Wagner and Pruss
    [WaP2002]_ eqn. 2.21, matching the
    International Association for the Properties of Water and Steam (IAPWS)
    1995 Formula

    .. math::


      Log (ew/22.064E6) = 647.096/T * ((-7.85951783 v
          + 1.84408259 * v**1.5
          - 11.7866497 * v**3
          + 22.6807411 * v**3.5
          - 15.9618719 * v**4
          + 1.80122502 * v**7.5))

      with T in [K] and ew in [Pa] and v = 1 - T/647.096
    '''
    T = _to_K(t, Kelvin)  # K
    if np.isnan(T):
        e = np.nan
    else:
        th = T / 273.16  # 1
        lne = (
            - 13-928168 * (1 - th**-1.5)
            + 34.7078238 * (1 - th**-1.25)
        )
        empa = np.exp(lne) * 0.000611657  # MPa
        e = empa * 10000.  # hPa
    if p is not None:
        e = _f(_to_hPa(p, hPa)) * e
    if hPa is True:
        return e  # hPa
    else:
        return 100. * e  # Pa


# ---------------------------------------------------------------------
def esat_w(*args, esat=_esat, **kwargs):
    '''
    Water vapor saturation pressure over water.
    Convenience function that allows selecting actual
    formula by name

    :param esat: formula to use, possible values `None` or
        the names given under :py:meth:`meltib.humidity.set_esat`
        In case of None, the standard function is chosen.
        Defaults to `None`.
    '''
    if esat is None:
        esat = _esat
    if esat == 'goff_gratch':
        return goff_gratch_w(*args, **kwargs)
    elif esat == 'hyland_wexler':
        return hyland_wexler_w(*args, **kwargs)
    elif esat == 'iapws':
        return iapws_w(*args, **kwargs)
    elif esat == 'magnus':
        return magnus_w(*args, **kwargs)
    elif esat == 'sonntag':
        return sonntag_w(*args, **kwargs)
    elif esat == 'tetens':
        return tetens_w(*args, **kwargs)
    else:
        raise RuntimeError('esat = {} not known'.format(esat))

# ---------------------------------------------------------------------


def esat_i(*args, esat=_esat, **kwargs):
    '''
    Water vapor saturation pressure over ice.
    Convenience function that allows selecting actual
    formula by name

    :param esat: formula to use, possible values `None` or
        the names given under :py:meth:`meltib.humidity.set_esat`
        In case of None, the standard function is chosen.
        Defaults to `None`.
    '''
    if esat == 'goff_gratch':
        return goff_gratch_i(*args, **kwargs)
    elif esat == 'hyland_wexler':
        return hyland_wexler_i(*args, **kwargs)
    elif esat == 'iapws':
        return iapws_i(*args, **kwargs)
    elif esat == 'magnus':
        return magnus_i(*args, **kwargs)
    elif esat == 'sonntag':
        return sonntag_i(*args, **kwargs)
    elif esat == 'tetens':
        return tetens_i(*args, **kwargs)
    else:
        raise RuntimeError('esat = {} not known'.format(esat))

# ---------------------------------------------------------------------


def set_esat(name):
    '''
    Set standard functions `esat_w` and `esat_i`
    for the saturation pressure of water vapor.
    :param name: name of the function set (character).
    Defaults to **magnus**.

    Possible names are:
    - **goff_gratch** for `goff_gratch_w` and `goff_gratch_i`
    - **hyland_wexler** for `hyland_wexler_w` and `hyland_wexler_i`
    - **iapws** for `iapws_w` and `iapws_i`
    - **magnus** for `magnus_w` and `magnus_i`
    - **sonntag** for `sonntag_w` and `sonntag_i`
    - **tetens** for `tetens_w` and `tetens_i`
    '''
    global _esat
    if name in ['goff_gratch', 'hyland_wexler', 'magnus', 'sonntag', 'tetens']:
        _esat = name
    else:
        raise ValueError('function name {} not known'.format(name))

# ---------------------------------------------------------------------


def get_esat():
    '''
    Returns current setting of
    standard functions `esat_w` and `esat_i`
    for the saturation pressure of water vapor.

    :return: name of the function set.
    :rtype: character
    '''
    return _esat


# ---------------------------------------------------------------------
def tdew(e, p=None, hPa=False, Kelvin=True):
    '''
    Returns dewpoint, with resprect to water.

    :param e: water vaport pressuer (float)
    :param p: (optional) air pressure (float). If `p` is ``None`` or missing,
      the correction for moist air is omittet, i.e. :math:`f(p) = 1`.
      Defaults to ``None``.
    :param hPa: (optional)
      if ``True``, `p` must be supplied in hPa.
      If ``False``, `p` must be supplied in Pa.
      Defaults to ``False``.
    :param Kelvin: (optional)
      if ``True``, result is given in Kelvin.
      If ``False``, result is given in degrees Celsius.
      Defaults to ``True``.
    :return:  dewpoint
    :rtype: float
    :raises ValueError: if `e` is less equal zero.
    :raises ValueError: if `p` is less equal zero.
    :raises ValueError: if `Kelvin` is not ``True``, ``False``, or missing
    :raises ValueError: if `hPa` is not ``True``, ``False``, or missing

    Results are valid in the range water -45°C to +60°C for pure phase water:

    .. math::

      t = (243.12 * np.log(e / (6.112 * _f(p)))) /
          (17.62 - np.log(e / (6.112 * _f(p))))

    '''
    e = _check('e',  e, 'float', ge=0.)
    p = _check('p',  p, 'float', gt=0., none=True)
    hPa = _check('hPa', hPa, 'bool')
    Kelvin = _check('Kelvin', Kelvin, 'bool', none=True)

    e = _to_hPa(e, hPa)  # hPa
    p = _to_hPa(p, hPa)  # hPa

    if np.isnan(e):
        t = np.nan
    else:
        t = (243.12 * np.log(e / (6.112 * _f(p)))) / \
            (17.62 - np.log(e / (6.112*_f(p))))  # C

    if Kelvin is False:
        return _to_C(t, Kelvin=False)  # °C
    else:
        return _to_K(t, Kelvin=False)  # K

# ---------------------------------------------------------------------


def tfrost(e, p=None, hPa=False, Kelvin=True):
    '''
    Returns frostpoint, with resprect to ice.

    :param e: water vaport pressuer (float)
    :param p: (optional) air pressure (float). If `p` is ``None`` or missing,
      the correction for moist air is omittet, i.e. :math:`f(p) = 1`.
      Defaults to ``None``.
    :param hPa: (optional)
      if ``True``, `e` and `p` must be supplied in hPa.
      If ``False``, `e` and `p` must be supplied in Pa.
      Defaults to ``False``.
    :param Kelvin: (optional)
      if ``True``, result is given in Kelvin.
      If ``False``, result is given in degrees Celsius.
      Defaults to ``True``.
    :return:  dewpoint
    :rtype: float
    :raises ValueError: if `e` is less equal zero.
    :raises ValueError: if `p` is less equal zero.
    :raises ValueError: if `Kelvin` is not ``True``, ``False``, or missing
    :raises ValueError: if `hPa` is not ``True``, ``False``, or missing

    Results are valid in the range water -45°C to +60°C for pure phase water:

    .. math::

      t = (272.62 * np.log(e / (6.112 * f(p)))) /
          (22.46 - np.log(e / (6.112 * f(p))))

    '''
    e = _check('e',  e, 'float', ge=0.)
    p = _check('p',  p, 'float', gt=0., none=True)
    hPa = _check('hPa', hPa, 'bool')
    Kelvin = _check('Kelvin', Kelvin, 'bool', none=True)

    e = _to_hPa(e, hPa)  # Pa
    p = _to_hPa(p, hPa)  # hPa

    if np.isnan(e):
        t = np.nan
    else:
        t = (272.62 * np.log(e / (6.112 * _f(p)))) / \
            (22.46 - np.log(e / (6.112 * _f(p))))  # C

    if Kelvin is False:
        return _to_C(t, Kelvin=False)  # °C
    else:
        return _to_K(t, Kelvin=False)  # K

# ---------------------------------------------------------------------


def psychro(t, tw, Kelvin=None, p=None, hPa=False):
    '''
    Psychrometric formula for the Assmann psychrometer (liquid water).

    :param t: air or dry-bulb temperature (float)
    :param tw: wet-bulb temperature (float)
    :param Kelvin: (optional)
      if ``True``, `t` ant `tw` must be given in in Kelvin.
      Ff ``False``, `t` ant `tw` must be given in degrees Celsius.
      If missing or ``None``, unit of `t` and `tw` is autodetected.
      Defaults to ``None``.
    :param p: (optional) air pressure (float).
      if ``None`` or missing,
      the correction for moist air is omittet, i.e. :math:`f(p) = 1`.
      Defaults to ``None``.
    :param hPa: (optional)
      if ``True``, `p` must be given in Hectopascal.
      If ``False``, `p` must be given in Pascal.
      Defaults to ``True``.
    :return: water vapor pressure in Pa or hPa, depending on `hPa`
    :rtype: float
    :raises ValueError: if `p` is less equal zero.
    :raises ValueError: if `Kelvin` is not ``True``, ``False``, ``None``,
      or missing
    :raises ValueError: if `hPa` is not ``True``, ``False``, ``None``,
      or missing

    This form of the equation was adopted by Forty-second session of the
    Executive Council (EC-XLII) in 1990 [WMO8]_:

    .. math::

      e = esat_w(t,p) - 6.53E-4 * (1 + 0.000944 * ti) * p * (t - ti)

    '''
    t = _check('t',  t, 'float', ge=-273.15)
    tw = _check('t', tw, 'float', ge=-273.15)
    Kelvin = _check('Kelvin', Kelvin, 'bool', none=True)
    p = _check('p',  p, 'float', gt=0., none=True)
    hPa = _check('hPa', hPa, 'bool')

    t = _to_C(t,  Kelvin)  # °C
    tw = _to_C(tw, Kelvin)  # °C
    if p is None:
        p1 = _to_hPa(pzero, False)  # hPa
    else:
        p1 = _to_hPa(p, hPa)  # hPa

    e = esat_w(tw, Kelvin=False, p=p, hPa=True) - 6.53E-4 * \
        (1 + 0.000944 * tw) * p1 * (t - tw)  # hPa

    if hPa is True:
        return e
    else:
        return e * 100.

# ---------------------------------------------------------------------


def psychro_ice(t, ti, Kelvin=None, p=None, hPa=False):
    '''
    Psychrometric formula for the Assmann psychrometer (ice).

    :param t: air or dry-bulb temperature (float)
    :param ti: wet-bulb temperature (float)
    :param Kelvin: (optional)
      if ``True``, `t` ant `ti` must be given in in Kelvin.
      If ``False``, `t` ant `ti` must be given in degrees Celsius.
      If missing or ``None``, unit of `t` and `ti` is autodetected.
      Defaults to ``None``.
    :param p: (optional) air pressure (float).
      if ``None`` or missing,
      the correction for moist air is omittet, i.e. :math:`f(p) = 1`.
      Defaults to ``None``.
    :param hPa: (optional)
      if ``True``, `p` must be given in Hectopascal.
      If ``False``, `p` must be given in Pascal.
      Defaults to ``True``.
    :return: water vapor pressure (in Pa or hPa, depending on `hPa`
    :rtype: float
    :raises ValueError: if `p` is less equal zero.
    :raises ValueError: if `Kelvin` is not ``True``, ``False``, ``None``,
      or missing
    :raises ValueError: if `hPa` is not ``True``, ``False``, ``None``,
      or missing

    This form of the equation was adopted by Forty-second session of the
    Executive Council (EC-XLII) in 1990 [WMO8]_:

    .. math::

      e = esat_i(t,p) - 5.75E-4 * p * (t - tw)

    '''
    t = _check('t',  t, 'float', ge=-273.15)
    ti = _check('t', ti, 'float', ge=-273.15)
    Kelvin = _check('Kelvin', Kelvin, 'bool', none=True)
    p = _check('p',  p, 'float', gt=0., none=True)
    hPa = _check('hPa', hPa, 'bool')

    t = _to_C(t,  Kelvin)  # °C
    ti = _to_C(ti, Kelvin)  # °C
    if p is None:
        p1 = _to_hPa(pzero, False)  # hPa
    else:
        p1 = _to_hPa(p, hPa)  # hPa

    e = esat_i(ti, Kelvin=False, p=p, hPa=True) - \
        5.75E-4 * p1 * (t - ti)  # hPa

    if hPa is True:
        return e  # hPa
    else:
        return 100. * e  # Pa


# ---------------------------------------------------------------------
def relhum(e, ew, percent=False):
    '''
    returns relative humidity.

    :param e: water vapor pressure (float).
      must be supplied in the same units as `ew`.
    :param ew: saturation water vapor pressure (float)
      must be supplied in the same units as `e`.
    :param percent: (optional)
      if ``True``, relative humidity is returned in percent.
      If ``False``, relative humidity is returned unitless (fraction of 1).
      Defaults to ``False``.
    :return: relative humidity (in % or unitless, depending on `percent`
    :rtype: float
    :raises ValueError: if `ew` is less equal zero.
    :raises ValueError: if `e` is less than zero.
    :raises ValueError: if `percent` is not ``True``, ``False``, or missing

    .. math::

      relhum = e / ew

    '''
    e = _check('e', e, 'float', ge=0.)
    ew = _check('ew', ew, 'float', gt=0.)
    relhum = e / ew
    if percent is True:
        return (relhum * 100.)  # %
    else:
        return (relhum)  # 1

# ---------------------------------------------------------------------


def mixr(e, p, gkg=False):
    '''
    returns approximative mixing ratio

    :param e: water vapor pressure (float).
      must be supplied in the same units as `ew`.
    :param ew: saturation water vapor pressure (float)
      must be supplied in the same units as `e`.
    :param gkg: (optional)
      if ``True``, mixr is returned in g/kg.
      If ``False``, mixr is returned in kg/kg.
      Defaults to ``False``.

    :return: relative humidity (in % or unitless, depending on `percent`
    :rtype: float
    :raises ValueError: if `e` is less equal zero.
    :raises ValueError: if `p` is less than zero.

    .. math::

      mixr = 0.62198 * e / p


    '''
    e = _check('e', e, 'float', ge=0.)
    p = _check('p', p, 'float', ge=0.)

    m = 0.62198 * e / p

    if gkg is True:
        return m * 1000.  # g/kg
    else:
        return m  # 1


# ---------------------------------------------------------------------
def inv_mixr(m, p, hPa=False, gkg=False):
    '''
    returns water vapor pressure calculated from mixing ratio

    :param m: mixing ratio (float).
    :param p: air pressure (float)
    :param hPa: (optional)
      if ``True``,  `p`, `e`, or `ew` must be supplied in hPa.
      If ``False``, `p`, `e`, or `ew` must be supplied in Pa.
      Defaults to ``False``.
    :param gkg: (optional)
      if ``True``, `q`, or `m` must be supplied in g/kg.
      If ``False``, `q`, or `m` must be supplied unitless, i.e. in kg/kg.
      Defaults to ``False``.

    :return: water vapor pressure in Pa or hPa, depending on `hPa`.
    :rtype: float
    :raises ValueError: if `m` is less equal zero.
    :raises ValueError: if `q` is less than zero.

    .. math::

      e   =   m * q /  0.62198

    '''
    m = _check('m', m, 'float', ge=0.)
    p = _check('p', p, 'float', ge=0.)

    if gkg is True:
        mm = m / 1000.
    else:
        mm = m
    pp = _to_Pa(p, hPa=hPa)

    e = mm * pp / 0.62198

    if hPa is True:
        return e / 100.  # hPa
    else:
        return e  # Pa

# ---------------------------------------------------------------------


def spech(e, p, gkg=False):
    '''
    returns specific hmidity


    :param e: water vapor pressure (float).
      must be supplied in the same units as `ew`.
    :param ew: saturation water vapor pressure (float)
      must be supplied in the same units as `e`.
    :param gkg: (optional)
      if ``True``, spech is returned in g/kg.
      If ``False``, spech is returned in kg/kg.
      Defaults to ``False``.

    :return: relative humidity (in % or unitless, depending on `percent`
    :rtype: float
    :raises ValueError: if `e` is less equal zero.
    :raises ValueError: if `p` is less than zero.

    .. math::

      spech = 0.62198 * e / (e *0.622 + p)


    '''
    e = _check('e', e, 'float', ge=0.)
    p = _check('p', p, 'float', ge=0.)

    q = 0.62198 * e / (e * 0.622 + p)  # 1

    if gkg is True:
        return q * 1000.  # g/kg
    else:
        return q  # 1

# ---------------------------------------------------------------------


def inv_spech(q, p, hPa=False, gkg=False):
    '''
    returns water vapor pressure calculated from specific hmidity

    :param q: specific humidity (float).
    :param ew: saturation water vapor pressure (float)
    :param hPa: (optional)
      if ``True``,  `p`, `e`, or `ew` must be supplied in hPa.
      If ``False``, `p`, `e`, or `ew` must be supplied in Pa.
      Defaults to ``False``.
    :param gkg: (optional)
      if ``True``, `q`, or `m` must be supplied in g/kg.
      If ``False``, `q`, or `m` must be supplied unitless, i.e. in kg/kg.
      Defaults to ``False``.

    :return: water vapor pressure in Pa or hPa, depending on `hPa`.
    :rtype: float
    :raises ValueError: if `e` is less equal zero.
    :raises ValueError: if `q` is less than zero.

    .. math::

      e   =   p * q / ( 0.62198 * (1 -  q))

    '''
    q = _check('q', q, 'float', ge=0.)
    p = _check('p', p, 'float', ge=0.)

    if gkg is True:
        qq = q / 1000.
    else:
        qq = q
    pp = _to_Pa(p, hPa=hPa)

    e = qq * pp / (0.62198 * (1 - qq))

    if hPa is True:
        return e / 100.  # hPa
    else:
        return e  # Pa

# ---------------------------------------------------------------------


class Humidity(object):
    '''
    The class Humidity represents a state of humid air that
    my be described and requested in any possible way

    :param t: (optional) air temperature / dry-bulb temperature (float).
    :param ew: (optional) saturation water vapor pressure (float).
     One of `t` and `ew` is required

    :param e: (optional) water vapor pressure (float).
    :param rh: (optional) relative humidity (float).
    :param td: (optional) dew point temperature (float).
    :param tf: (optional) frost point temperature (float).
    :param tw: (optional) wet-bulb temperature (float).
    :param ti: (optional) ice / frozen wet-bulb temperature (float).
    :param q: (optional) specific humidity (float).
    :param m: (optional) mixing ration (float).
    :param rhow: (optional) water vapor density in :math:`kg/m^3` (float).
    :param p: (manatory with `q` or `m`, else optional) air pressure (float).
      Defaults to 101320 Pa.
    :param Kelvin: (optional) if ``False``, all temperatures are assumed to
      be Kelvin. If ``False``, all temperatures are assumed to be Celsius.
      If missing of ``None``, unit  temperatures are autodetected. Defaults to
      ``None``.
    :param hPa: (optional)
      if ``True``,  `p`, `e`, or `ew` must be supplied in hPa.
      If ``False``, `p`, `e`, or `ew` must be supplied in Pa.
      Defaults to ``False``.
    :param gkg: (optional)
      if ``True``, `q`, or `m` must be supplied in g/kg.
      If ``False``, `q`, or `m` must be supplied unitless, i.e. in kg/kg.
      Defaults to ``False``.
    :param percent: (optional)
      if ``True``, relative humidity is returned in percent.
      If ``False``, relative humidity is returned unitless (fraction of 1).
      Defaults to ``False``.
    :raises ValueError: if `Kelvin` is not ``True``, ``False``, None,
      or missing
    :raises ValueError: if `hPa` is not ``True``, ``False``, or missing
    :raises ValueError: if `gkg` is not ``True``, ``False``, or missing
    :raises ValueError: if `percent` is not ``True``, ``False``, or missing
    :raises ValueError:  if an invalid combination of input parameters
      is given, i.e. other than e, tf, tw and t, tw and t, q and p, m and p
    '''

    def __init__(self, t=None, p=None,
                 e=None, ew=None, rh=None, td=None, tf=None,
                 tw=None, ti=None, q=None, m=None, rhow=None,
                 Kelvin=None, hPa=False, gkg=False, percent=False):

        self.Kelvin = _check('Kelvin',  Kelvin, 'bool', none=True)
        self.hPa = _check('hPa',        hPa, 'bool')
        self.gkg = _check('gkg',        gkg, 'bool')
        self.percent = _check('percent', percent, 'bool')

        if p is None:
            self.p = pzero
        else:
            self.p = _check('p', p, 'float', gt=0.)

        loc = locals()

        tem = {x: loc[x] for x in ['t', 'ew']}
        if _only(tem, ['t']):
            if self.Kelvin:
                self.t = _check('t', t, 'float', ge=0)
            else:
                self.t = _check('t',   t, 'float', gt=-273.15)
        elif _only(tem, ['ew']):
            self.t = _check('ew', ew, 'float', gt=0.)
        else:
            raise ValueError('one of t or ew is required')

        num = {x: loc[x] for x in ['e', 'rh', 'td', 'tf',
                                   't', 'tw', 'ti', 'q', 'm', 'p', 'rhow']}
        if _only(num, ['rh', 't'], ['p']):
            if percent is True:
                rrh = _check('rh', rh, 'float', ge=0., le=100.) / 100.
            else:
                rrh = _check('rh', rh, 'float', ge=0., le=1.)
            self.e = rrh * esat_w(t, Kelvin=Kelvin, p=p, hPa=hPa)
        elif _only(num, ['e'], ['p', 't']):
            self.e = e
        elif _only(num, ['td'], ['p', 't']):
            self.e = esat_w(td, Kelvin=Kelvin, p=p, hPa=hPa)
        elif _only(num, ['tf'], ['p', 't']):
            self.e = esat_i(tf, Kelvin=Kelvin, p=p, hPa=hPa)
        elif _only(num, ['tw', 't'], ['p']):
            self.e = psychro(t, tw, Kelvin=Kelvin, p=p, hPa=hPa)
        elif _only(num, ['ti', 't'], ['p']):
            self.e = psychro_ice(t, ti, Kelvin=Kelvin, p=p, hPa=hPa)
        elif _only(num, ['q', 'p'], ['t']):
            self.e = inv_spech(q=q, p=p, hPa=hPa, gkg=gkg)
        elif _only(num, ['m', 'p'], ['t']):
            self.e = inv_mixr(m=m, p=p, hPa=hPa, gkg=gkg)
        elif _only(num, ['m', 'p'], ['t']):
            self.e = inv_mixr(m=m, p=p, hPa=hPa, gkg=gkg)
        elif _only(num, ['rhow', 't'], ['p']):
            self.e = gas_p(rhow, t, Kelvin=Kelvin, hPa=hPa, gas='water')
        else:
            raise ValueError('invalid combination of parameters: {}'.format(
                ' '.join(k for k, v in loc.items() if v is not None)))

    # ---------------------------------------------------------------------

    def rh(self, percent=None):
        '''
        returns relative humidity of the Humidity object.

        :param percent: (optional)
          if ``True``, relative humidity is returned in percent.
          If ``False``, relative humidity is returned unitless (fraction of 1).
          Defaults to object attribute `percent`.
        :return: relative humidity (in % or unitless, depending on `percent`).
        :rtype: float
        '''
        if percent is None:
            percent = self.percent
        else:
            percent = _check('percent', percent, 'bool')

        e_w = esat_w(self.t, Kelvin=self.Kelvin, p=self.p, hPa=self.hPa)
        rh = relhum(self.e, e_w)

        if percent is True:
            return rh * 100.
        else:
            return rh

    # ---------------------------------------------------------------------
    def m(self, gkg=None):
        '''
        returns mixing ratio of the Humidity object.

        :param gkg: (optional)
          if ``True``, `q`, or `m` must be supplied in g/kg.
          If ``False``, `q`, or `m` must be supplied unitless, i.e. in kg/kg.
          Defaults to object attribute `gkg`.
        :return: mixing ratio (in g/kg or unitless, depending on `gkg`).
        :rtype: float
        '''
        if gkg is None:
            gkg = self.gkg
        else:
            gkg = _check('gkg', gkg, 'bool')

        m = mixr(self.e, self.p, gkg=gkg)

        return m

    # ---------------------------------------------------------------------
    def q(self, gkg=None):
        '''
        returns specific humidity of the Humidity object.

        :param gkg: (optional)
          if ``True``, `q`, or `m` must be supplied in g/kg.
          If ``False``, `q`, or `m` must be supplied unitless, i.e. in kg/kg.
          Defaults to object attribute `gkg`.
        :return: specific humidity (in g/kg or unitless, depending on `gkg`).
        :rtype: float
        '''
        if gkg is None:
            gkg = self.gkg
        else:
            gkg = _check('gkg', gkg, 'bool')

        q = spech(self.e, self.p, gkg=gkg)

        return q

    # ---------------------------------------------------------------------
    def td(self, Kelvin=None):
        '''
        returns dew point of the Humidity object.

        :param Kelvin: (optional) if ``False``, all temperatures are assumed to
          be Kelvin. If ``False``, all temperatures are assumed to be Celsius.
          If missing of ``None``, unit  temperatures are autodetected.
          Defaults to ``None``.
        :return: dew point (in K or °C, depending on `Kelvin`).
        :rtype: float
        :raises ValueError: if `Kelvin` is not ``True``, ``False``, None,
          or missing
        '''
        if Kelvin is None:
            Kelvin = self.Kelvin
        else:
            Kelvin = _check('Kelvin', Kelvin, 'bool')

        td = tdew(self.e, p=self.p, hPa=self.hPa, Kelvin=Kelvin)

        if Kelvin is None or Kelvin is True:
            return _to_K(td, Kelvin=Kelvin)
        else:
            return _to_C(td, Kelvin=Kelvin)

    # ---------------------------------------------------------------------
    def tf(self, Kelvin=None):
        '''
        returns frost point of the Humidity object.

        :param Kelvin: (optional)
          if ``False``, all temperatures are assumed to be Kelvin.
          If ``False``, all temperatures are assumed to be Celsius.
          If missing of ``None``, unit  temperatures are autodetected.
          Defaults to ``None``.
        :return: frost point (in K or °C, depending on `Kelvin`).
        :rtype: float
        :raises ValueError: if `Kelvin` is not ``True``, ``False``, None,
          or missing
        '''
        if Kelvin is None:
            Kelvin = self.Kelvin
        else:
            Kelvin = _check('Kelvin', Kelvin, 'bool')

        tf = tfrost(self.e, p=self.p, hPa=self.hPa, Kelvin=Kelvin)

        if Kelvin is None or Kelvin is True:
            return _to_K(tf, Kelvin=Kelvin)
        else:
            return _to_C(tf, Kelvin=Kelvin)

    # ---------------------------------------------------------------------
    def tw(self, Kelvin=None):
        '''
        returns wet-bulb temperature of the Humidity object.

        :param Kelvin: (optional)
          if ``False``, all temperatures are assumed to be Kelvin.
          If ``False``, all temperatures are assumed to be Celsius.
          If missing of ``None``, unit  temperatures are autodetected.
          Defaults to ``None``.
        :return: wet-bulb temperatur (in K or °C, depending on `Kelvin`).
        :rtype: float
        :raises ValueError: if `Kelvin` is not ``True``, ``False``, None,
          or missing

        Since there is no closed analytical solution, we use the
        iterative method of [MAR2009]_,
        adapted for use with the Magnus formula :eq:`magnus` for
        saturation water vapor from [WMO8]_:

        1) set inital value

        .. math::
          tw_0 = t - 5

        2) calculate next value:

        .. math::
          tw' &= tl - \\frac{ ew(tw_i) - e}
          { 6.53E-4 * (1 + 0.000944 * tw_i) * p} \\\\
          tw_{i+1} &= tw_i + 0.5 * (tw' - tw_i)

        3) iterate until:

        .. math::
          \\left| tw_{i+1} - tw_i \\right| < 0.01K

        If the iteration does not converge within 100 iterations,
        ``np.nan`` is returned

        '''
        if Kelvin is None:
            Kelvin = self.Kelvin
        else:
            Kelvin = _check('Kelvin', Kelvin, 'bool')

        tl = _to_C(self.t, Kelvin=self.Kelvin)  # °C
        pp = _to_hPa(self.p, hPa=self.hPa)  # hPa
        ee = _to_hPa(self.e, hPa=self.hPa)  # hPa

        tw = tl - 5  # °C
        for i in range(100):
            twl = tw
            twn = tl - (esat_w(tw, Kelvin=False, hPa=True) - ee) / \
                (6.53E-4 * (1 + 0.000944 * tw) * pp)  # °C
            tw = tw + 0.5 * (twn - tw)
            if abs(tw - twl) <= 0.01:
                break
        else:
            tw = np.nan

        if Kelvin is False:
            return tw
        else:
            return _to_K(tw, Kelvin=False)

    # ---------------------------------------------------------------------
    def ti(self, Kelvin=None):
        '''
        returns ice / frozen wet-bulb temperature of the Humidity object.

        :param Kelvin: (optional) if ``False``, all temperatures are assumed to
          be Kelvin. If ``False``, all temperatures are assumed to be Celsius.
          If missing of ``None``, unit  temperatures are autodetected.
          Defaults to ``None``.
        :return: ice / frozen wet-bulb temperatur (in K or °C,
          depending on `Kelvin`).
        :rtype: float
        :raises ValueError: if `Kelvin` is not ``True``, ``False``, None,
          or missing

        Since there is no closed analytical solution, we use the
        iterative method of [MAR2009]_,
        adapted for use with the Magnus formula :eq:`magnice` for
        saturation water vapor from [WMO8]_:

        1) set inital value

        .. math::
          ti_0 = t - 5

        2) calculate next value:

        .. math::
          ti' &= tl - \\frac{ ew(ti_i) - e}{ 5.75E-4 * p } \\\\
          ti_{i+1} &= ti_i + 0.5 * (ti' - ti_i)

        3) iterate until:

        .. math::
          \\left| ti_{i+1} - ti_i \\right| < 0.01K

        If the iteration does not converge within 100 iterations,
        ``np.nan`` is returned
        '''
        if Kelvin is None:
            Kelvin = self.Kelvin
        else:
            Kelvin = _check('Kelvin', Kelvin, 'bool')

        tl = _to_C(self.t, Kelvin=self.Kelvin)  # °C
        pp = _to_hPa(self.p, hPa=self.hPa)  # hPa
        ee = _to_hPa(self.e, hPa=self.hPa)  # hPa

        ti = tl - 5  # °C
        for i in range(100):
            til = ti
            tin = tl - (esat_i(ti, Kelvin=False, hPa=True) - ee) / \
                (5.75E-4 * pp)  # °C
            ti = ti + 0.5 * (tin - ti)
            if abs(ti - til) <= 0.01:
                break
        else:
            ti = np.nan

        if Kelvin is False:
            return ti
        else:
            return _to_K(ti, Kelvin=False)

    # ---------------------------------------------------------------------
    def virt_inc(self):
        '''
        returns virtual temperature increment of the Humidity object.

        :return: virtual temperature increment (in K).
        :rtype: float

        '''
        T = _to_K(self.t, Kelvin=self.Kelvin)
        qq = spech(self.e, self.p, gkg=False)

        vi = T * 0.61 * qq

        return vi

    # ---------------------------------------------------------------------
    def tvirt(self, Kelvin=None):
        '''
        calculates virtual temperature of the Humidity object.

        :param Kelvin: (optional) if ``False``, all temperatures are assumed to
          be Kelvin. If ``False``, all temperatures are assumed to be Celsius.
          If missing of ``None``, unit  temperatures are autodetected.
          Defaults to ``None``.
        :return: virtual temperature (in K or °C, depending on `Kelvin`).
        :rtype: float
        :raises ValueError: if `Kelvin` is not ``True``, ``False``, None,
          or missing
        '''
        if Kelvin is None:
            Kelvin = self.Kelvin
        else:
            Kelvin = _check('Kelvin', Kelvin, 'bool')

        vi = self.virt_inc()

        T = _to_K(self.t, Kelvin=self.Kelvin)  # K
        Tv = T + vi  # K

        if Kelvin is False:
            return _to_C(Tv, Kelvin=True)
        else:
            return _to_K(Tv, Kelvin=True)

    # ---------------------------------------------------------------------
    def tequi(self, Kelvin=None):
        '''
        calculates equivalent temperature of the Humidity object.

        :param Kelvin: (optional)
          if ``False``, all temperatures are assumed to be Kelvin.
          If ``False``, all temperatures are assumed to be Celsius.
          If missing of ``None``, unit  temperatures are autodetected.
          Defaults to ``None``.
        :return: equivalent temperature (in K or °C, depending on `Kelvin`).
        :rtype: float
        :raises ValueError: if `Kelvin` is not ``True``, ``False``, None,
          or missing

        equivalent temperature is the temperature of an air parcel from which
        all the water vapor has been removed by an adiabatic process:

        .. math:: T_{e} \\approx T + \\frac{L_v}{c_p} m

        '''
        mm = self.m(gkg=False)

        T = _to_K(self.t, Kelvin=self.Kelvin)

        Te = T + mm * Lv / cp

        if Kelvin is False:
            return _to_C(Te, Kelvin=True)
        else:
            return _to_K(Te, Kelvin=True)

    # ---------------------------------------------------------------------
    def tequipot(self, Kelvin=None):
        '''
        calculates equipotential temperature of the Humidity object.

        :param Kelvin: (optional) if ``False``, all temperatures are assumed to
          be Kelvin. If ``False``, all temperatures are assumed to be Celsius.
          If missing of ``None``, unit  temperatures are autodetected.
          Defaults to ``None``.
        :return: equipotential temperature (in K or °C, depending on `Kelvin`).
        :rtype: float
        :raises ValueError: if `Kelvin` is not ``True``, ``False``, None,
          or missing

        equipotential temperature is the temperature of an air parcel
        from which all the water vapor has been removed
        by an adiabatic process [Stu1988]_:

        .. math:: \\Theta_{e} = T_e
                  \\left( \\frac{p_0}{p} \\right) ^{\\frac{R}{c_p}}
                  \\approx \\left( T + \\frac{L_v}{c_p} m \\right)
                  \\left( \\frac{p_0}{p} \\right)^{\\frac{R}{c_p}}

        '''
        te = self.tequi()

        tep = Tpot(te, p=self.p, Kelvin=self.Kelvin, hPa=self.hPa)

        return tep

    # ---------------------------------------------------------------------
    def rhow(self):
        '''
        calculates water vapor density of the Humidity object.

        :return: density (in :math:`kg/m^3` ).
        :rtype: float

        '''
        rw = gas_rho(self.e, self.t, gas='water',
                     Kelvin=self.Kelvin, hPa=self.hPa)

        return rw

# ---------------------------------------------------------------------


def __init__():
    print('tralala')
    set_esat('magnus')
