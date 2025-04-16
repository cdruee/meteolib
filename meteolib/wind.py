#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Wind- and turbulence-related claculations and conversions
"""
import logging

import numpy as np

from ._utils import _check, _only
from .constants import gn, kappa

logger = logging.getLogger(__name__)

DISPLACEMENT_FACTOR = 6.5


class LogWind(object):
    '''
    The class LogWind represents a logarithmic wind profile that
    my be described and requested in any possible compination of
    parameters.

    :param ust: friction velocity in m/s (float).
    :param z0:  surface rougness length in m (float).
    :param d:   displacement height in m (float). Defaults to 0 m.
    :param u:   wind speed at height `z` in m/s  (float).
    :param z:   height above ground, where `u` is measured, in m (float).
    :param u2:   wind speed at height `z2` in m/s  (float).
    :param z2:   height above ground, where `u2` is measured, in m (float).

    Possible parameter combinations for initialization are:
      - `ust` and `z0`, optionally `d`
      - `u`, `z` and `z0`, optionally `d`
      - `u`, `z` and `ust`, optionally `d`
      - `u`, `z`, `u2`, and `z2`, optionally `d`

    :ivar ust: friction velocity in m/s
    :ivar z0: surface rougness length in m
    :ivar ~.d: displacement height in m


    '''

    def __init__(self, ust=None, z0=None, d=None, u=None, z=None,
                 u2=None, z2=None):
        loc = locals()
        par = {x: loc[x] for x in ['ust', 'z0', 'd', 'u', 'z', 'u2', 'z2']}

        if _only(par, ['ust', 'z0'], ['d']):
            ust = _check('ust', ust, 'float', ge=0.)
            self.ust = ust
            z0 = _check('z0',  z0, 'float',  gt=0.)
            self.z0 = z0
            d = _check('d',    d, 'float',  ge=0., none=True)
            if d is None:
                self.d = DISPLACEMENT_FACTOR * self.z0
            else:
                self.d = d  # : displacement height in m

        elif _only(par, ['u', 'z', 'z0'], ['d']):
            z0 = _check('z0',  z0, 'float',  gt=0.)
            self.z0 = z0
            d = _check('d',   d, 'float',  ge=0., none=True)
            if d is None:
                self.d = DISPLACEMENT_FACTOR * self.z0
            else:
                self.d = d
            u = _check('u',   u, 'float',   ge=0.)
            z = _check('z',   z, 'float',   ge=self.d)
            self.ust = (kappa * u) / np.log((z-self.d)/self.z0)

        elif _only(par, ['u', 'z', 'ust'], ['d']):
            ust = _check('ust', ust, 'float', ge=0.)
            self.ust = ust
            u = _check('u',     u, 'float', gt=0.)
            d = _check('d',     d, 'float', ge=0., none=True)
            if d is None:
                z = _check('z',   z, 'float', gt=0.)
                self.z0 = z / (DISPLACEMENT_FACTOR + np.exp(kappa * u / ust))
                self.d = DISPLACEMENT_FACTOR * z0
            else:
                z = _check('z',   z, 'float', gt=self.d)
                self.z0 = (z-d) / np.exp(kappa * u / ust)
                self.d = d

        elif _only(par, ['u', 'z', 'u2', 'z2'], ['d']):
            d = _check('d',     d, 'float', ge=0., none=True)
            if d is None:
                u = _check('u',   u, 'float', ge=0.)
                z = _check('z',   z, 'float', ge=0.)
                u2 = _check('u2', u2, 'float', gt=u)
                z2 = _check('z2', z2, 'float', gt=z)
                z0l = np.min([0.5 * z, 0.01])
                if z <= DISPLACEMENT_FACTOR * z0l:
                    raise ValueError('z too low')
                for i in range(100):
                    logger.debug(format((z0l, DISPLACEMENT_FACTOR*z0l, z)))
                    zz = np.exp((u*np.log(z2-DISPLACEMENT_FACTOR*z0l) -
                                u2*np.log(z-DISPLACEMENT_FACTOR*z0l)) /
                                (u - u2))
                    z0 = z0l * (zz/z0l)**0.25
                    if abs(np.log(z0/z0l)) < 0.001:
                        break
                    z0l = z0
                else:
                    raise ValueError('did not converge')
                self.z0 = z0
                self.d = DISPLACEMENT_FACTOR*z0
                self.ust = u / np.log((z-self.d) / self.z0)
            else:
                self.d = d
                u = _check('u',   u,  'float', ge=0.)
                z = _check('z',   z,  'float', ge=self.d)
                u2 = _check('u2',  u2, 'float', gt=u)
                z2 = _check('z2',  z2, 'float', gt=z)
                if z <= d + d / DISPLACEMENT_FACTOR:
                    raise ValueError('z below displacement')
                self.z0 = np.exp(
                    (u*np.log(z2-self.d) - u2*np.log(z-self.d)) / (u - u2))
                self.ust = u / np.log((z-self.d) / self.z0)

    def u(self, z):
        '''
        returns wind speed at certain level

        :param z: height above ground in m (float)
        :return: wind speed in m/s
        :rtype: float
        '''
        z = _check('z', z, 'float', ge=0.)
        u = (self.ust / kappa) * np.log((z-self.d)/self.z0)
        return u

    def gradu(self, z):
        '''
        returns wind speed vertical gradient at certain level

        :param z: height above ground in m (float)
        :return: wind speed gradient  in 1/s
        :rtype: float
        '''
        z = _check('z', z, 'float', ge=0.)
        u = self.ust / (kappa * (z - self.d))
        return u

# -------------------------------------------------------------------


class DiabaticWind(object):
    '''
    The class LogWind represents a logarithmic wind profile that
    my be described and requested in any possible compination of
    parameters.

    :param ust: friction velocity in m/s (float).
    :param z0:  surface rougness length in m (float).
    :param d:   displacement height in m (float). Defaults to 0 m.
    :param u:   wind speed at height `z` in m/s  (float).
    :param z:   height above ground, where `u` is measured, in m (float).
    :param zoL: Monin-Boukhov stability parameter z/L, unitless (float).

    Possible parameter combinations for initialization are:
      - `ust` and `z0`, 'z', 'L', optionally `d`
      - `u`, `z` and `z0`, optionally `d`
      - `u`, `z` and `ust`, optionally `d`
      - `u`, `z`, `u2`, and `z2`, optionally `d`

    :ivar ust: friction velocity in m/s
    :ivar z0: surface rougness length in m
    :ivar ~.d: displacement height in m
    :ivar LOb: Obukhov length in m


    '''

    def __init__(self, ust=None, z0=None, d=0., u=None, z=None,
                 LOb=None, zoL=None):
        loc = locals()
        par = {x: loc[x] for x in ['ust', 'z0', 'd', 'u', 'z', 'LOb', 'zoL']}

        if _only(par, ['ust', 'z0', 'LOb'], ['d']):
            ust = _check('ust', ust, 'float', ge=0.)
            self.ust = ust
            z0 = _check('z0',  z0, 'float',  gt=0.)
            self.z0 = z0
            LOb = _check('LOb', LOb, 'float')
            self.LOb = LOb
            d = _check('d',    d, 'float',  ge=0., none=True)
            if d is None:
                self.d = DISPLACEMENT_FACTOR * self.z0
            else:
                self.d = d  # : displacement height in m

        elif _only(par, ['u', 'z', 'z0', 'LOb'], ['d']):
            z0 = _check('z0',  z0, 'float',  gt=0.)
            self.z0 = z0
            d = _check('d',   d, 'float',  ge=0., none=True)
            if d is None:
                self.d = DISPLACEMENT_FACTOR * self.z0
            else:
                self.d = d
            u = _check('u',   u, 'float',   ge=0.)
            z = _check('z',   z, 'float',   ge=self.d)
            self.ust = (kappa * u) / np.log((z-self.d)/self.z0)
            LOb = _check('LOb', LOb, 'float')
            self.LOb = LOb
            zoL = z / LOb
            self.ust = ((kappa * u) /
                        (np.log((z-self.d)/self.z0) - psi_m(zoL)))

        elif _only(par, ['u', 'z', 'z0', 'zoL'], ['d']):
            z0 = _check('z0',  z0, 'float',  gt=0.)
            self.z0 = z0
            d = _check('d',   d, 'float',  ge=0., none=True)
            if d is None:
                self.d = DISPLACEMENT_FACTOR * self.z0
            else:
                self.d = d
            u = _check('u',   u, 'float',   ge=0.)
            z = _check('z',   z, 'float',   ge=self.d)
            zoL = _check('zoL', zoL, 'float')
            self.LOb = z/zoL
            self.ust = ((kappa * u) /
                        (np.log((z-self.d)/self.z0) - psi_m(zoL)))

#    elif _only(par,['u','z','ust',],['d']):
#      ust = _check('ust', ust, 'float', ge=0.)
#      self.ust = ust
#      u   = _check('u',     u, 'float', gt=0.)
#      d   = _check('d',     d, 'float', ge=0., none=True)
#      if d is None:
#        z   = _check('z',   z, 'float', gt=0.)
#        self.z0 = z  / ( DISPLACEMENT_FACTOR + np.exp( kappa * u / ust) )
#        self.d = DISPLACEMENT_FACTOR * z0
#      else:
#        z   = _check('z',   z, 'float', gt=self.d)
#        self.z0 = (z-d) / np.exp( kappa * u / ust)
#        self.d = d

    def u(self, z):
        '''
        returns wind speed at certain level

        :param z: height above ground in m (float)
        :return: wind speed in m/s
        :rtype: float
        '''
        z = _check('z', z, 'float', ge=0.)
        zoL = z / self.LOb
        u = (self.ust / kappa) * (np.log((z-self.d) /
                                         self.z0) - psi_m(zoL))
        return u

    def gradu(self, z):
        '''
        returns wind speed vertical gradient at certain level

        :param z: height above ground in m (float)
        :return: wind speed gradient  in 1/s
        :rtype: float
        '''
        z = _check('z', z, 'float', ge=0.)
        zoL = z / self.LOb
        gu = psi_m(zoL) * self.ust / (kappa * (z - self.d))
        return gu

# -------------------------------------------------------------------


def phi_m(zoL):
    '''
    Universal functions for momentum according to
    [Bus1971]_, as recalculated by [Hog1988]_.

    :param zoL: Monin-Boukhov stability parameter z/L, unitless (float).
    :return: univesal fuction value (unitless)
    :rtype: float
    '''
    zoL = _check('zoL', zoL, 'float')
    phi_1s = 6.0
    phi_1u = 19.3
    if zoL > 0.:
        phi = 1. + phi_1s * zoL
    elif zoL < 0.:
        phi = (1. - phi_1u * zoL) ^ (-1./4.)
    else:  # if zoL == 0.
        phi = 1.
    return phi

# -------------------------------------------------------------------


def phi_H(zoL):
    '''
    Universal functions for turbulent fluxes of scalars according to
    [Bus1971]_, using the numerical "modified Kansas expression"
    valuesas recalculated by [Hog1988]_.

    :param zoL: Monin-Boukhov stability parameter z/L, unitless (float).
    :return: univesal fuction value (unitless)
    :rtype: float
    '''
    zoL = _check('zoL', zoL, 'float')
    phi_1s = 7.8
    phi_1u = 11.6
    phi_0 = 0.95
    if zoL > 0.:
        phi = phi_0 + phi_1s * zoL
    elif zoL < 0.:
        phi = phi_0 * (1. - phi_1u * zoL) ^ (-1./2.)
    else:  # if zoL == 0.
        phi = phi_0
    return phi

# -------------------------------------------------------------------


def psi_m(zoL):
    '''
    Integrated universal functions for momentum according to
    [Pau1970]_ (unstable) and [HoB1988]_ (stable),
    using the numerical "modified Kansas expression"
    values recalculated by [Hog1988]_
    :param zoL: Monin-Boukhov stability parameter z/L, unitless (float).
    :return: universal fuction value (unitless)
    :rtype: float
    '''
    zoL = _check('zoL', zoL, 'float')
    phi_1s = 6.0
    phi_1u = 19.3
    if zoL > 0.:
        psi = - phi_1s * zoL
    else:
        x = (1 - phi_1u * zoL) ** 0.25
        psi = (2. * np.log((1. + x) / 2.)
               + np.log((1. + x ** 2) / 2.)
               - 2. * np.arctan(x)
               + np.pi / 2.)
    return psi

# -------------------------------------------------------------------


def psi_H(zoL):
    '''
    Integrated universal functions for momentum according to
    [Pau1970]_ (unstable) and [HoB1988]_ (stable),
    using the numericcal "modified Kansas expression"
    values by [Hog1988]_
    :param zoL: Monin-Boukhov stability parameter z/L, unitless (float).
    :return: univesal fuction value (unitless)
    :rtype: float
    '''
    zoL = _check('zoL', zoL, 'float')
    phi_2s = 7.8
    phi_2u = 11.6
    phi_20 = 0.95
    if zoL > 0.:
        psi = - phi_2s * zoL
    else:
        x = (1 - phi_2u * zoL) ** 0.25
        psi = phi_20 * (2. * np.log((1. + x ** 2) / 2.))
    return psi


# -------------------------------------------------------------------
def charnock(ust):
    '''
    Charnock's relation: water surface roughness
    :param ust: friction velocity in m/s (float).
    :return: water surface roughness length in m
    :rtype: float
    '''
    ust = _check('ust', ust, 'float', ge=0.)
    alpha_c = 0.015  # Charnock parameter
    z0 = alpha_c * (ust**2) / gn
    return z0

# -------------------------------------------------------------------


def transfer(wind, z0, d=0.):
    '''
    returns wind profile for a different place with the same
    synoptic conditions but different roughness (and displacement) length(s)
    :param wind: logarithmic wind profile at **old** location (LogWind).
    :param z0: surface roughness length at **new** location in m (float).
    :param z0: (optional) displacement height at **new** location in m (float).
    :return: logarithmic wind profile at **new** location
    :rtype: LogWind
    '''
    if not isinstance(wind, LogWind):
        raise TypeError('agrmument wind must be of type LogWind')
    Refheight = 250.  # m
    uref = wind.u(Refheight)
    new = LogWind(u=uref, z=Refheight, z0=z0, d=d)
    return new

# -------------------------------------------------------------------


def vectormean(ff=None, dd=None, pairs=None):
    '''
    returns vector mean of wind vectors given by means of speed an direction
    or lists of speed and direction

    :param ff: wind speed in any unit,
      requires `dd` (array-like, same lenght as `dd`)
    :param dd: wind direction in degrees,
      requires `ff` (array-like, same lenght as `ff`)
    :param pairs: pairs of wind speed (any unit) and direction (in degrees)
      (list of tuples)
    :return: mean wind speed and wind direction in degrees
    :rtype: float, float
    '''
    if pairs is not None:
        ff = [x[0] for x in pairs]
        dd = [x[1] for x in pairs]
    elif ff is not None and dd is not None:
        if len(ff) != len(dd):
            raise ValueError('ff and dd mut be of same length')
    else:
        raise ValueError('einter pairs or ff and dd must be given')

    u = [-f*np.sin(np.deg2rad(d)) for f, d in zip(ff, dd)]
    v = [-f*np.cos(np.deg2rad(d)) for f, d in zip(ff, dd)]
    umean = np.nanmean(u)
    vmean = np.nanmean(v)

    fm = np.sqrt(umean**2 + vmean**2)
    dm = np.rad2deg(np.arctan2(-umean, -vmean))

    return fm, dm

# -------------------------------------------------------------------


def uv2dir(u, v, pos=True):
    '''
    converts (horizontal) wind vector components to speed and direction

    :param u: eastward wind component (float or array-like)
    :param v: northward wind component (float or array-like)
    :param pos: (optional) return positive values:
        If ``True`` returned wind directions are in the range 0..360,
        if ``False`` values are in the range -180..180.
        Defaults to ``True``.
    :return: wind speed, same unit as input, and wind direction in degrees
    :rtype: float, float or np.array, np.array
    '''
    # test if shapes are identical
    if np.shape(u) != np.shape(v):
        raise ValueError('u and v must be of the same shape')
    # convert scalars into array
    if np.isscalar(u) and np.isscalar(v):
        u = np.array([u])
        v = np.array([v])
    else:
        u = np.array(u)
        v = np.array(v)
    # calculate
    speed = np.sqrt(u**2 + v**2)
    wdir = np.rad2deg(np.arctan2(-u, -v))
    # make values (range -180..180) positive (0..360)
    if pos is True:
        wdir = wdir % 360
    # convert single values back to scalars
    if speed.shape == (1,):
        speed = speed[0]
        wdir = wdir[0]
    return speed, wdir

# -------------------------------------------------------------------


def dir2uv(ff, dd):
    '''
    converts (horizontal) wind vector from speed and direction
    to vector components

    :param ff: wind speed in any unit (float or array-like)
    :param dd: wind direction in degrees (float or array-like)
    :return: eastward wind component and northward wind component
        in the same unit `speed` is given in
    :rtype: float, float or np.array, np.array
    '''
    # test if shapes are identical
    if np.shape(ff) != np.shape(dd):
        raise ValueError('ff and dd must be of the same shape')
    # convert scalars into array
    if np.isscalar(ff) and np.isscalar(dd):
        ff = np.array([ff])
        dd = np.array([dd])
    else:
        ff = np.array(ff)
        dd = np.array(dd)
    # check values:
    if any(ff < 0.):
        raise ValueError('ff must be positive')
    # calculate
    u = - ff * np.sin(np.deg2rad(dd))
    v = - ff * np.cos(np.deg2rad(dd))
    # convert single values back to scalars
    if u.shape == (1,):
        u = u[0]
        v = v[0]
    return u, v
