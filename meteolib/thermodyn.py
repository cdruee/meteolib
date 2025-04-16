#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''Thermodynamics laws and equations

'''
import numpy as np

from ._utils import _check
from .constants import R, Rco2, Rw
from .pressure import _to_hPa
from .temperature import _to_C, _to_K


def gas_p(rho, T, q=None, Kelvin=None, hPa=False, gkg=False, gas='air'):
    r'''
    ideal gas equation, calculate gas pressure

    :param rho: density of gas in :math:`kg/m^3` .
    :param T:  gas temperature in :math:`K` or :math:`^\circ C`,
      depending on `Kelvin` .
    :param rho: density of gas (in :math:`kg/m^3` .
    :param gas: kind of gas.
      One of `dry` (dry air), `moist` (humid air), `water` (water vapor),
      or `co2` (carbon dioxide).
      Defaults to `dry`.
    :param Kelvin: (optional) if ``False``, all temperatures are assumed to
      be Kelvin. If ``False``, all temperatures are assumed to be Celsius.
      If missing of ``None``, unit  temperatures are autodetected.
      Defaults to ``None``.
    :param hPa: (optional)
      if ``True``,  `p`, `e`, or `ew` must be supplied in hPa.
      If ``False``, `p`, `e`, or `ew` must be supplied in Pa.
      Defaults to ``False``.
    :param gkg: (optional)
      if ``True``, `q`, or `m` must be supplied in g/kg.
      If ``False``, `q`, or `m` must be supplied unitless, i.e. in kg/kg.
      Defaults to ``False``.

    :return: air pressure (in hPa or Pa, depending on `hPa`).
    :rtype: float

    '''
    if gas in ['dry', 'moist']:
        Rgas = R
    elif gas in ['water']:
        Rgas = Rw
    elif gas in ['co2']:
        Rgas = Rco2
    else:
        raise ValueError('unknow value for "gas": {}'.format(gas))
    T = _to_K(T, Kelvin)
    if gas in ['dry', 'water', 'co2']:
        Tgas = T
    elif gas == 'moist':
        #
        # virtual Temperature
        # hardcoded bevcause cannot be impoprtet mutually
        if gkg is True:
            Tgas = T * (1 + 0.00061*q)
        else:
            Tgas = T * (1 + 0.61*q)

    p = rho * Rgas * Tgas

    if hPa is True:
        return _to_hPa(p, hPa=False)
    else:
        return p


def gas_rho(p, T, q=None, Kelvin=None, hPa=False, gkg=False, gas='air'):
    r'''
    ideal gas equation, calculate gas density

    :param p: air pressure (in hPa or Pa, depending on `hPa`).
    :param T:  gas temperature in :math:`K` or :math:`^\circ C` ,
      depending on `Kelvin` .
    :param Kelvin: (optional) if ``False``, all temperatures are assumed to
      be Kelvin. If ``False``, all temperatures are assumed to be Celsius.
      If missing of ``None``, unit  temperatures are autodetected. Defaults to
      ``None``.
    :param hPa: (optional)
      if ``True``,  `p`, `e`, or `ew` must be supplied in hPa.
      If ``False``, `p`, `e`, or `ew` must be supplied in Pa.
      Defaults to ``False``.
    :param gas: kind of gas.
      One of `dry` (dry air), `moist` (humid air), `water` (water vapor),
      or `co2` (carbon dioxide).
      Defaults to `dry`.

    :return: density of gas in :math:`kg/m^3`.
    :rtype: float

    '''
    if gas in ['air', 'dry', 'moist']:
        Rgas = R
    elif gas in ['water']:
        Rgas = Rw
    elif gas in ['co2']:
        Rgas = Rco2
    else:
        raise ValueError('unknow value for "gas": {}'.format(gas))
    T = _to_K(T, Kelvin)
    if gas in ['air', 'dry', 'water', 'co2']:
        Tgas = T
    elif gas == 'moist':
        #
        # virtual Temperature
        # hardcoded bevcause cannot be impoprtet mutually
        if gkg is True:
            Tgas = T * (1 + 0.00061*q)
        else:
            Tgas = T * (1 + 0.61*q)
    if hPa is True:
        pgas = p*100.
    else:
        pgas = p

    rho = pgas / (Rgas * Tgas)

    return rho


def gas_t(p, rho, q=None, Kelvin=None, hPa=False, gas='air'):
    r'''
    ideal gas equation, calculate gas temperature

    :param p: air pressure (in hPa or Pa, depending on `hPa`).
    :param rho: density of gas in :math:`kg/m^3` .
    :param Kelvin: (optional) if ``False``, all temperatures are assumed to
      be Kelvin. If ``False``, all temperatures are assumed to be Celsius.
      If missing of ``None``, unit  temperatures are autodetected. Defaults to
      ``None``.
    :param hPa: (optional)
      if ``True``,  `p`, `e`, or `ew` must be supplied in hPa.
      If ``False``, `p`, `e`, or `ew` must be supplied in Pa.
      Defaults to ``False``.
    :param gas: kind of gas.
      One of `dry` (dry air), `moist` (humid air), `water` (water vapor),
      or `co2` (carbon dioxide).
      Defaults to `dry`.

    :return: gas temperature in :math:`K` or :math:`^\circ C`,
      depending on `Kelvin` .
    :rtype: float

    '''
    if gas in ['air', 'dry', 'moist']:
        Rgas = R
    elif gas in ['water']:
        Rgas = Rw
    elif gas in ['co2']:
        Rgas = Rco2
    else:
        raise ValueError('unknow value for "gas": {}'.format(gas))

    Tgas = p / (Rgas * rho)

    if gas in ['air', 'dry', 'water', 'co2']:
        T = Tgas
    elif gas == 'moist':
        #    Tgas = tvirt_inv(Tgas,q,Kelvin=True,gkg=gkg)
        raise RuntimeError('not yet implemented')

    if Kelvin is True:
        return T
    else:
        return _to_C(T, Kelvin=True)

# ----------------------------------------------------------------------
# not so constant constants
# ----------------------------------------------------------------------


def cpt(t, Kelvin=None):
    r'''
    Specific heat capacity of dry air, as function of temperature

    :param t:  air temperature in :math:`K` or :math:`^\circ C`,
      depending on `Kelvin` .
    :param Kelvin: (optional) if ``False``, all temperatures are assumed to
      be Kelvin. If ``False``, all temperatures are assumed to be Celsius.
      If missing of ``None``, unit  temperatures are autodetected. Defaults to
      ``None``.

    :return: specific heat capacity of **dry air** at constant pressure
      (in :math:`J kg^-1 K^-1`).
    :rtype: float

    .. math::

      c_{p} ~=~ 1005.60 + 0.017211 T + 0.000392 T^2

    valid for temperature t between –40°C and 40°C and
    for barometric pressures near one atmosphere [And2006]_.
    Minimum value of 1005.41:math:`J kg^-1 K^-1` at 21.95°C.

    '''
    t = _check('t', t, 'float', ge=0.)
    Kelvin = _check('Kelvin', Kelvin, 'bool', none=True)
    #
    th = _to_C(t, Kelvin=Kelvin)
    cp = 1005.60 + 0.017211 * th + 0.000392 * (th**2)
    # equation (A20) in Garrat, J.R., 1992. 'The Atmospheric Boundary
    #     layer', Cambridge University Press:
    # cp = 1005.16 + 0.013763 * th + 2.97265E-4 * (th**2)
    return cp

# ----------------------------------------------------------------------


def cpvt(t, Kelvin=None):
    r'''
    Specific heat capacity of water vapor, as function of temperature

    :param t:  air temperature in :math:`K` or :math:`^\circ C`,
      depending on `Kelvin`.
    :param Kelvin: (optional) if ``False``, all temperatures are assumed to
      be Kelvin. If ``False``, all temperatures are assumed to be Celsius.
      If missing of ``None``, unit  temperatures are autodetected. Defaults to
      ``None``.

    :return: specific heat capacity of **water vapor** at constant pressure
      (in :math:`J kg^-1 K^-1`).
    :rtype: float

    .. math::

      c_{pv} ~=~ 1858 + 3.820 \times 10^{-1} T + 4.220
      \times 10^{-4} T^2 - 1.996 \times 10^{-7} T^3

    "should be accurate for all near-surface atmospheric temperatures"
    [And2006]_.

    '''
    t = _check('t', t, 'float', ge=0.)
    Kelvin = _check('Kelvin', Kelvin, 'bool', none=True)
    #
    th = _to_C(t, Kelvin=Kelvin)
    cp = 1858 + 3.820E-1*th + 4.220E-4 * (th**2) - 1.996E-7*(th**3)
    return cp

# ----------------------------------------------------------------------


def cwt(t, Kelvin=None):
    r'''
    Specific heat capacity of liquid water, as function of temperature

    :param t:  air temperature in :math:`K` or :math:`^\circ C`,
      depending on `Kelvin`.
    :param Kelvin: (optional) if ``False``, all temperatures are assumed to
      be Kelvin. If ``False``, all temperatures are assumed to be Celsius.
      If missing of ``None``, unit  temperatures are autodetected. Defaults to
      ``None``.

    :return: specific heat capacity of **liquid water** at constant pressure
      (in :math:`J kg^-1 K^-1`).
    :rtype: float

    .. math::

      c_{w} ~=~ 4217.4 - 3.720283 T^2 - 2.654387
      \times 10^{-3} T^3 + 2.093236 \times 10^{-5} T^4

    "should be accurate for all near-surface atmospheric temperatures"
    [And2006]_.

    '''
    t = _check('t', t, 'float', ge=0.)
    Kelvin = _check('Kelvin', Kelvin, 'bool', none=True)
    #
    th = _to_C(t, Kelvin=Kelvin)
    cw = 4217.4 - 3.720283 * th**2 - 2.654387E-3 * th**3 + 2.093236E-5 * th**4
    return cw

# ----------------------------------------------------------------------


def cit(t, Kelvin=None):
    r'''
    Specific heat capacity of ice, as function of temperature

    :param t:  air temperature in :math:`K` or :math:`^\circ C`,
      depending on `Kelvin`.
    :param Kelvin: (optional) if ``False``, all temperatures are assumed to
      be Kelvin. If ``False``, all temperatures are assumed to be Celsius.
      If missing of ``None``, unit  temperatures are autodetected. Defaults to
      ``None``.

    :return: specific heat capacity of **ice** at constant pressure
      (in :math:`J kg^-1 K^-1`).
    :rtype: float

    .. math::

      c_{w} ~=~ -114.19 + 8.1288 T + 3.421 T~
      \exp \left[ -(T/125.1)^2 \right]

    "should be accurate for all near-surface atmospheric temperatures"
    [And2006]_.

    '''
    t = _check('t', t, 'float', ge=0.)
    Kelvin = _check('Kelvin', Kelvin, 'bool', none=True)
    #
    th = _to_K(t, Kelvin=Kelvin)
    # here : th in K !
    ci = -114.19 + 8.1288 * th + 3.421 * th * np.exp(-(th/125.1)**2)
    return ci

# ----------------------------------------------------------------------


def Lvt(t, Kelvin=None, source="andreas"):
    r'''
    Latent heat of vaporization, as function of temperature

    :param t:  air temperature in :math:`K` or :math:`^\circ C`,
      depending on `Kelvin`.
    :param Kelvin: (optional) if ``False``, all temperatures are assumed to
      be Kelvin. If ``False``, all temperatures are assumed to be Celsius.
      If missing of ``None``, unit  temperatures are autodetected. Defaults to
      ``None``.
    :param source: (optional) select source of equation, see below. Defaults to
      `andreas`.

    :return: latent heat of vaporization (in :math:`J kg^-1`).
    :rtype: float

    Sources:

    andreas
      (**default**)
      [And2006]_ gives for Use in Marine Meteorology and states
      "the Lv values are within 0.3% of the Smithsonian
      values for temperatures from 0° to 60°C".

      :math:`L_{v} ~=~ `
      :math:`\left( 25.00 - 0.02274 ~ \vartheta \right) \times 10^{5}`

    henderson-sellers
      A different for of fit is given by [HeS1984]_:

      :math:`L_{v} ~=~ 1.91846E6 \left(T/\left(T-33.91\right)\right)^2`

    ecpack
      ecpack by Arain van Dijk and Arnold Moene uses:

      :math:`L_{v} ~=~ \left(2501 - 2.375 \vartheta\right) 10^3`

    dake
      a linear fit to the Smithonian tables is presented by [Dak1972]_:

      :math:`L_{v} ~=~ \left(2501 - 2.357*\vartheta\right) 10^3`

    hyland-wexler
      The American Society of Heating, Refrigerating and
      Air-Conditioning Engineers
      in their "ASHRAE Handbook : Fundamentals" presents values
      calculated using the formulation of [HyW1983]_ which evaluates to:

      .. math::
        :nowrap:

        \begin{eqnarray*}
        L_{v} &=& 3139.817121 \\
              &-& 2.390755077 * T \\
              &+& 4.3309357163E-4 * T^2 \\
              &-& 1.573317502E-6 * T^3 \\
              &+& 2.94377462234E-9 * T^4 \\
              &-& 1.7508262E-12 * T^5 \\
              &+& 6.059E-7 * (T - 403.128)^3 \\
        \end{eqnarray*}

    '''
    t = _check('t', t, 'float', ge=0.)
    Kelvin = _check('Kelvin', Kelvin, 'bool', none=True)
    source = _check('source', source, 'str')
    # -----
    # Andreas (2006)
    # Physical Constants and Functions For Use in Marine Meteorology
    # Lv = ( 25.00 - 0.02274 * th ) * 1.0E5
    if source == "andreas":
        th = _to_C(t, Kelvin=Kelvin)
        Lv = (25.00 - 0.02274 * th) * 1.0E5
    # -----
    # Henderson‐Sellers (1984) DOI: 10.1002/qj.49711046626
    # T in K !
    # lv = 1.91846E6*(T/(T-33.91))**2
    elif source == "henderson-sellers":
        T = _to_K(t, Kelvin=Kelvin)
        Lv = 1.91846E6*(T/(T-33.91))**2
    # -----
    # Arnold Moene (EC-Pack). Source?
    # Lv = (2501 - 2.375*th)*1.0E3
    elif source == "ecpack":
        th = _to_C(t, Kelvin=Kelvin)
        Lv = (2501 - 2.375 * th) * 1.0E3
    # identical to Dake (1972) except typo
    # (there coefficient 0.5631 cal/K² -> 2.357 J/K²)
    # citation given by Drake is Simithonian tables (List 1951)
    # but there are only tabbed values "calculated after Goff_Gratch (1945)"
    elif source == "dake":
        th = _to_C(t, Kelvin=Kelvin)
        Lv = (2501 - 2.357 * th) * 1.0E3
    # -----
    # Goff_Gratch (1945) is not available
    # -----
    # (ASHRAE 2013) presents tabbed values calculated after
    # (Hyland and Wexler 2013), (ASHRAE 2013) are identical to
    # (List 1951) values, converted using 4.1839 J/cal
    elif source == "hyland-wexler":
        T = _to_K(t, Kelvin=Kelvin)
        Lv = (
            3139.817121 +
            -2.390755077 * T +
            4.3309357163E-4 * T**2 +
            -1.573317502E-6 * T**3 +
            2.94377462234E-9 * T**4 +
            -1.7508262E-12 * T**5 +
            6.059E-7 * (T - 403.128)**3
        ) * 1.0E3
    else:
        raise ValueError('evaporation heat source {} unknown'.format(source))

    return Lv

# ----------------------------------------------------------------------


def Lst(t, Kelvin=None, source="andreas"):
    r'''
    Latent heat of sublimation, as function of temperature

    :param t:  air temperature in :math:`K` or :math:`^\circ C`,
      depending on `Kelvin`.
    :param Kelvin: (optional) if ``False``, all temperatures are assumed to
      be Kelvin. If ``False``, all temperatures are assumed to be Celsius.
      If missing of ``None``, unit  temperatures are autodetected. Defaults to
      ``None``.
    :param source: (optional) select source of equation, see below. Defaults to
      `andreas`.

    :return: latent heat of sublimation (in :math:`J kg^-1`).
    :rtype: float

    Sources:

    andreas
      (**default**)
      [And2006]_ gives for Use in Marine Meteorology and states
      "the Ls values are within 0.2% of the Smithsonian
      values for temperatures from -50° to 0°C".
      :math:`L_{s} ~=~ \left( 28.34 - 0.00149 T \right) \times 10^{5}`

    hyland-wexler
      The American Society of Heating, Refrigerating and
      Air-Conditioning Engineers
      in their "ASHRAE Handbook : Fundamentals" presents values
      calculated using the formulation of [HyW1983]_ which evaluates to:

      .. math::
        :nowrap:

         \begin{eqnarray*}
         L_{v} &=&  2645.475 \\
           &+& 1.5292786 * T \\
           &-& 0.00254657837 * T^2 \\
           &-& 2.5511992E-6 * T^3 \\
           &+& 2.8726608E-9 * T^4 \\
           &-& 1.7508262E-12 * T^5
         \end{eqnarray*}


    '''
    t = _check('t', t, 'float', ge=0.)
    Kelvin = _check('Kelvin', Kelvin, 'bool', none=True)
    source = _check('source', source, 'str')
    # -----
    # Andreas (2006)
    # Physical Constants and Functions For Use in Marine Meteorology
    # Lv = ( 25.00 - 0.02274 * th ) * 1.0E5
    if source == "andreas":
        th = _to_C(t, Kelvin=Kelvin)
        Ls = (28.34 - 0.00149 * th) * 1.0E5
    # -----
    # (ASHRAE 2013) presents tabbed values calculated after
    # (Hyland and Wexler 2013), (ASHRAE 2013) are identical to
    # (List 1951) values, converted using 4.1839 J/cal
    elif source == "hyland-wexler":
        T = _to_K(t, Kelvin=Kelvin)
        Ls = (
            2645.475 +
            1.5292786 * T +
            -0.00254657837 * T**2 +
            -2.5511992E-6 * T**3 +
            2.8726608E-9 * T**4 +
            -1.7508262E-12 * T**5
        ) * 1.0E3
    else:
        raise ValueError('sublimation heat source {} unknown'.format(source))

    return Ls
