#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
physical constants

'''

#: Stefan-Boltzmann constant (:math:`W m^{-2} K^{-4}`)
#: of the identically named physical law that relatet
#: electromagnetig radion intensity to black-body temperature.
#: :math:`P = \sigma T^4`
#: Recommended value Adopted from [CODATA]_.
sigma = 5.570367E-8

#: Absolute temperature (:math:`K`)
#: of the normal ice point
#: Offset between Kelvin and Celsius Tmeperature scales.
#: Adopted by WMO in [WMO8]_.
Tzero = 273.15

#: Absolute temperature (:math:`K`)
#: of the triple point of water
#: by definition of [ITS90]_.
#: Adopted by WMO in [WMO8]_.
Ttriple = 273.16


#:  Standard acceleration of gravity (:math:`m s^{-2}`).
#:  Adopted by WMO in [WMO8]_.
gn = 9.806

#: Density of mercury (:math:`kg m^{–3}`) at 0 °C.
#: Adopted by WMO in [WMO8]_.
rhoHg = 1.35951E4

#: Reference pressure (:math:`Pa`).
#: One of the typical references levels --
#: Mean Sea Level (MSL), station altitude or the 1013.2 hPa plane.
#: -- adopted by WMO in [WMO8]_.
pzero = 101320.

# ----------------------------------------------------------------------
# WMO 180, Table 4.1: THERMODYNAMIC CONSTANTS AND FUNCTIONS
# (except constants above)

# 1. Basic constants

#: Apparent molecular weight of dry air (:math:`kg mol^{-1}`).
#: Adopted by WMO in [WMO188]_.
M = 28.9644E-3

#: Gas constant for 1 mole of ideal gas (:math:`J mol^{-1}`)
#: Adopted by WMO in [WMO188]_.
Rstar = 8.31432

#: Gas constant for dry air (:math:`J kg^-1 K^-1`).
#: Adopted by WMO in [WMO188]_.
R = 287.05

#: Molecular weight of water vapour (:math:`kg mol^{-1}`).
#: Adopted by WMO in [WMO188]_.
Mw = 18.0153E-3

#: Gas constant for water vapour (:math:`J kg^-1 K^-1`).
#: Adopted by WMO in [WMO188]_.
Rw = 461.51

# 2. Specific heat capacities

#: Specific heat capacity of **dry air**
#: at constant pressure (:math:`J kg^-1 K^-1`).
#: Recommended by WMO in [WMO188]_.
#: For temperature-dependent formulation see :py:meth:`meteolib.thermodyn.cpt`
cp = 1005.

#: Specific heat capacity of **dry air**
#: at constant volume (:math:`J kg^-1 K^-1`).
#: Recommended by WMO in [WMO188]_.
cv = 718.

#: Specific heat capacity of **liquid water** (:math:`J kg^-1 K^-1`).
#: Recommended by WMO in [WMO188]_.
#: For temperature-dependent formulation see :py:meth:`meteolib.thermodyn.cwt`
cw = 4179.

#: Specific heat capacity of frozen water / **ice**  (:math:`J kg^-1 K^-1`).
#: Recommended by WMO in [WMO188]_.
#: For temperature-dependent formulation see :py:meth:`meteolib.thermodyn.cit`
ci = 2090.

#: Specific heat capacity of **water vapor**
#: at constant pressure (:math:`J kg^-1 K^-1`).
#: Recommended by WMO in [WMO188]_.
#: For temperature-dependent formulation see :py:meth:`meteolib.thermodyn.cpvt`
cpv = 1850.

#: Specific heat capacity of **water vapor**
#: at constant volume (:math:`J kg^-1 K^-1`).
#: Recommended by WMO in [WMO188]_.
cvv = 1390.


# 3. Heats of transformation of phase of water

#: Specific heat of sublimation (:math:`J kg^-1`).
#: Recommended by WMO in [WMO188]_.
#: For temperature-dependent formulation see :py:meth:`meteolib.thermodyn.Lst`
Ls = 2.835E6

#: Specific heat of vaporization / evaportaion (:math:`J kg^-1`).
#: Mean value named by WMO in [WMO188]_.
#: For temperature-dependent formulation see :py:meth:`meteolib.thermodyn.Lvt`
Lv = 2.501E6

# ----------------------------------------------------------------------
# other constants

#: Apparent molecular weight of carbon dioxide (:math:`kg mol^{-1}`).
#: Adopted by NIST [NIST]_.
Mco2 = 44.0095E-3

#: Gas constant for carbon dioxide (:math:`J kg^-1 K^-1`).
#: calculated from `Rstar` and `Mco2`.
Rco2 = Rstar/Mco2

#: von-Karman constant. Re-evaluated value by [Hog1985]_
kappa = 0.40

#: Solar constant (:math:`W m^{-2}`)
#: a synonym for mean extra-terrestrial solar irradiance.
#: said to be defined by WMO and used in WMO documents
#: and many research papers, no immediate source.
solar = 1367.

#: Eart's nominal radius (m)
#: Defined in [ISO2533]_
r_earth = 6356766
