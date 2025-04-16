#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculations on atmospheric radiation and astronomy

solar position functions (spa_*) are based on [ReA2008]_.
earth's clock error values are taken from [Esp2006]_.



"""

import datetime as dt
import logging

import numpy
import numpy as np

import pandas as pd

import pytz

from ._utils import _check, _expand_to_series_like
from .pressure import pa2mmhg
from .constants import sigma


# ---------------------------------------------------------------------

def _clock_error_formula(year: int):
    '''
    Calculate historical and future values of the earth's clock error
    This is the internal funtion accepting on single scalar
    '''

    if year < -500:
        u = (year - 1820.) / 100.
        delta_T = -20. + 32. * u**2

    elif year <= 500:
        u = year / 100.
        delta_T = (10583.6 - 1014.41 * u + 33.78311 * u**2
                   - 5.952053 * u**3 - 0.1798452 * u**4
                   + 0.022174192 * u**5 + 0.0090316521 * u**6)

    elif year <= 1600.:
        u = (year - 1000.) / 100.
        delta_T = (1574.2 - 556.01 * u + 71.23472 * u**2
                   + 0.319781 * u**3 - 0.8503463 * u**4
                   - 0.005050998 * u**5 + 0.0083572073 * u**6)

    elif year <= 1700.:
        t = year - 1600.
        delta_T = 120. - 0.9808 * t - 0.01532 * t**2 + t**3 / 7129.

    elif year <= 1800.:
        t = year - 1700.
        delta_T = (8.83 + 0.1603 * t - 0.0059285 * t**2
                   + 0.00013336 * t**3 - t**4 / 1174000)

    elif year <= 1860.:
        t = year - 1800.
        delta_T = (13.72 - 0.332447 * t + 0.0068612 * t**2
                   + 0.0041116 * t**3 - 0.00037436 * t**4
                   + 0.0000121272 * t**5 - 0.0000001699 * t**6
                   + 0.000000000875 * t**7)

    elif year <= 1900.:
        t = year - 1860.
        delta_T = (7.62 + 0.5737 * t - 0.251754 * t**2
                   + 0.01680668 * t**3 - 0.0004473624 * t**4
                   + t**5 / 233174.)

    elif year <= 1920.:
        t = year - 1900.
        delta_T = (-2.79 + 1.494119 * t - 0.0598939 * t**2
                   + 0.0061966 * t**3 - 0.000197 * t**4)

    elif year <= 1941.:
        t = year - 1920.
        delta_T = 21.20 + 0.84493 * t - 0.076100 * t**2 + 0.0020936 * t**3

    elif year <= 1961.:
        t = year - 1950.
        delta_T = 29.07 + 0.407 * t - t**2 / 233 + t**3 / 2547

    elif year <= 1986.:
        t = year - 1975.
        delta_T = 45.45 + 1.067 * t - t**2 / 260 - t**3 / 718

    elif year <= 2005.:
        t = year - 2000.
        delta_T = (63.86 + 0.3345 * t - 0.060374 * t**2
                   + 0.0017275 * t**3 + 0.000651814 * t**4
                   + 0.00002373599 * t**5)

    elif year <= 2050.:
        t = year - 2000.
        delta_T = 62.92 + 0.32217 * t + 0.005589 * t**2

    elif year <= 2150.:
        delta_T = -20. + 32. * ((year - 1820.) / 100.)**2 - \
            0.5628 * (2150 - year)

    else:  # after 2150
        u = (year - 1820.) / 100.
        delta_T = -20. + 32. * u**2

    return delta_T


# ======== start of solar position functions (internal) ================
#
# apparent radius of the sun (in degrees)
#
_SUN_RADIUS = 0.26667
#
# delta UT1 is a fraction of a second, positive or negative value,
# that is added to the UTC to adjust for the Earth irregular rotational rate.
# It is derived from observation, but predicted values are transmitted in
# code in some time signals, e.g. weekly by the U.S. Naval Observatory (USNO)
#
# we neglect it here
_DELTA_UT1 = 0
#
# Atmospheric refraction at sunrise and sunset (0.5667 deg is typical)
# valid range: -5   to   5 degrees
_ATMOS_REFRAC = 0.5667
#
# ---------------------------------------------------------------------
#
#  Earth Periodic Terms
#
_EARTH_PERIODIC_TERMS_DATA_L = [
    # L0
    [
        [175347046.0, 0, 0],
        [3341656.0, 4.6692568, 6283.07585],
        [34894.0, 4.6261, 12566.1517],
        [3497.0, 2.7441, 5753.3849],
        [3418.0, 2.8289, 3.5231],
        [3136.0, 3.6277, 77713.7715],
        [2676.0, 4.4181, 7860.4194],
        [2343.0, 6.1352, 3930.2097],
        [1324.0, 0.7425, 11506.7698],
        [1273.0, 2.0371, 529.691],
        [1199.0, 1.1096, 1577.3435],
        [990, 5.233, 5884.927],
        [902, 2.045, 26.298],
        [857, 3.508, 398.149],
        [780, 1.179, 5223.694],
        [753, 2.533, 5507.553],
        [505, 4.583, 18849.228],
        [492, 4.205, 775.523],
        [357, 2.92, 0.067],
        [317, 5.849, 11790.629],
        [284, 1.899, 796.298],
        [271, 0.315, 10977.079],
        [243, 0.345, 5486.778],
        [206, 4.806, 2544.314],
        [205, 1.869, 5573.143],
        [202, 2.458, 6069.777],
        [156, 0.833, 213.299],
        [132, 3.411, 2942.463],
        [126, 1.083, 20.775],
        [115, 0.645, 0.98],
        [103, 0.636, 4694.003],
        [102, 0.976, 15720.839],
        [102, 4.267, 7.114],
        [99, 6.21, 2146.17],
        [98, 0.68, 155.42],
        [86, 5.98, 161000.69],
        [85, 1.3, 6275.96],
        [85, 3.67, 71430.7],
        [80, 1.81, 17260.15],
        [79, 3.04, 12036.46],
        [75, 1.76, 5088.63],
        [74, 3.5, 3154.69],
        [74, 4.68, 801.82],
        [70, 0.83, 9437.76],
        [62, 3.98, 8827.39],
        [61, 1.82, 7084.9],
        [57, 2.78, 6286.6],
        [56, 4.39, 14143.5],
        [56, 3.47, 6279.55],
        [52, 0.19, 12139.55],
        [52, 1.33, 1748.02],
        [51, 0.28, 5856.48],
        [49, 0.49, 1194.45],
        [41, 5.37, 8429.24],
        [41, 2.4, 19651.05],
        [39, 6.17, 10447.39],
        [37, 6.04, 10213.29],
        [37, 2.57, 1059.38],
        [36, 1.71, 2352.87],
        [36, 1.78, 6812.77],
        [33, 0.59, 17789.85],
        [30, 0.44, 83996.85],
        [30, 2.74, 1349.87],
        [25, 3.16, 4690.48],
    ],
    # L1
    [
        [628331966747.0, 0, 0],
        [206059.0, 2.678235, 6283.07585],
        [4303.0, 2.6351, 12566.1517],
        [425.0, 1.59, 3.523],
        [119.0, 5.796, 26.298],
        [109.0, 2.966, 1577.344],
        [93, 2.59, 18849.23],
        [72, 1.14, 529.69],
        [68, 1.87, 398.15],
        [67, 4.41, 5507.55],
        [59, 2.89, 5223.69],
        [56, 2.17, 155.42],
        [45, 0.4, 796.3],
        [36, 0.47, 775.52],
        [29, 2.65, 7.11],
        [21, 5.34, 0.98],
        [19, 1.85, 5486.78],
        [19, 4.97, 213.3],
        [17, 2.99, 6275.96],
        [16, 0.03, 2544.31],
        [16, 1.43, 2146.17],
        [15, 1.21, 10977.08],
        [12, 2.83, 1748.02],
        [12, 3.26, 5088.63],
        [12, 5.27, 1194.45],
        [12, 2.08, 4694],
        [11, 0.77, 553.57],
        [10, 1.3, 6286.6],
        [10, 4.24, 1349.87],
        [9, 2.7, 242.73],
        [9, 5.64, 951.72],
        [8, 5.3, 2352.87],
        [6, 2.65, 9437.76],
        [6, 4.67, 4690.48],
    ],
    # L2
    [
        [52919.0, 0, 0],
        [8720.0, 1.0721, 6283.0758],
        [309.0, 0.867, 12566.152],
        [27, 0.05, 3.52],
        [16, 5.19, 26.3],
        [16, 3.68, 155.42],
        [10, 0.76, 18849.23],
        [9, 2.06, 77713.77],
        [7, 0.83, 775.52],
        [5, 4.66, 1577.34],
        [4, 1.03, 7.11],
        [4, 3.44, 5573.14],
        [3, 5.14, 796.3],
        [3, 6.05, 5507.55],
        [3, 1.19, 242.73],
        [3, 6.12, 529.69],
        [3, 0.31, 398.15],
        [3, 2.28, 553.57],
        [2, 4.38, 5223.69],
        [2, 3.75, 0.98]
    ],
    # L3
    [
        [289.0, 5.844, 6283.076],
        [35, 0, 0],
        [17, 5.49, 12566.15],
        [3, 5.2, 155.42],
        [1, 4.72, 3.52],
        [1, 5.3, 18849.23],
        [1, 5.97, 242.73]
    ],
    # L4
    [
        [114.0, 3.142, 0],
        [8, 4.13, 6283.08],
        [1, 3.84, 12566.15]
    ],
    # L5
    [
        [1, 3.14, 0]
    ]
]
_EARTH_PERIODIC_TERMS_DATA_B = [
    # B0
    [
        [280.0, 3.199, 84334.662],
        [102.0, 5.422, 5507.553],
        [80, 3.88, 5223.69],
        [44, 3.7, 2352.87],
        [32, 4, 1577.34]
    ],
    # B1
    [
        [9, 3.9, 5507.55],
        [6, 1.73, 5223.69]
    ]
]
_EARTH_PERIODIC_TERMS_DATA_R = [
    # R0
    [
        [100013989.0, 0, 0],
        [1670700.0, 3.0984635, 6283.07585],
        [13956.0, 3.05525, 12566.1517],
        [3084.0, 5.1985, 77713.7715],
        [1628.0, 1.1739, 5753.3849],
        [1576.0, 2.8469, 7860.4194],
        [925.0, 5.453, 11506.77],
        [542.0, 4.564, 3930.21],
        [472.0, 3.661, 5884.927],
        [346.0, 0.964, 5507.553],
        [329.0, 5.9, 5223.694],
        [307.0, 0.299, 5573.143],
        [243.0, 4.273, 11790.629],
        [212.0, 5.847, 1577.344],
        [186.0, 5.022, 10977.079],
        [175.0, 3.012, 18849.228],
        [110.0, 5.055, 5486.778],
        [98, 0.89, 6069.78],
        [86, 5.69, 15720.84],
        [86, 1.27, 161000.69],
        [65, 0.27, 17260.15],
        [63, 0.92, 529.69],
        [57, 2.01, 83996.85],
        [56, 5.24, 71430.7],
        [49, 3.25, 2544.31],
        [47, 2.58, 775.52],
        [45, 5.54, 9437.76],
        [43, 6.01, 6275.96],
        [39, 5.36, 4694],
        [38, 2.39, 8827.39],
        [37, 0.83, 19651.05],
        [37, 4.9, 12139.55],
        [36, 1.67, 12036.46],
        [35, 1.84, 2942.46],
        [33, 0.24, 7084.9],
        [32, 0.18, 5088.63],
        [32, 1.78, 398.15],
        [28, 1.21, 6286.6],
        [28, 1.9, 6279.55],
        [26, 4.59, 10447.39]
    ],
    # R1
    [
        [103019.0, 1.10749, 6283.07585],
        [1721.0, 1.0644, 12566.1517],
        [702.0, 3.142, 0],
        [32, 1.02, 18849.23],
        [31, 2.84, 5507.55],
        [25, 1.32, 5223.69],
        [18, 1.42, 1577.34],
        [10, 5.91, 10977.08],
        [9, 1.42, 6275.96],
        [9, 0.27, 5486.78]
    ],
    # R2
    [
        [4359.0, 5.7846, 6283.0758],
        [124.0, 5.579, 12566.152],
        [12, 3.14, 0],
        [9, 3.63, 77713.77],
        [6, 1.87, 5573.14],
        [3, 5.47, 18849.23]
    ],
    # R3
    [
        [145.0, 4.273, 6283.076],
        [7, 3.92, 12566.15]
    ],
    # R4
    [
        [4, 2.56, 6283.08]
    ]
]

_EARTH_PERIODIC_TERMS_L = [pd.DataFrame.from_records(x,
                                                     columns=['A', 'B', 'C'])
                           for x in _EARTH_PERIODIC_TERMS_DATA_L]
_EARTH_PERIODIC_TERMS_B = [pd.DataFrame.from_records(x,
                                                     columns=['A', 'B', 'C'])
                           for x in _EARTH_PERIODIC_TERMS_DATA_B]
_EARTH_PERIODIC_TERMS_R = [pd.DataFrame.from_records(x,
                                                     columns=['A', 'B', 'C'])
                           for x in _EARTH_PERIODIC_TERMS_DATA_R]

# ---------------------------------------------------------------------
#
# Periodic Terms for the nutation in longitude and obliquity
#
_NUTATION_PERIODIC_TERMS_DATA_Y = [
    [0, 0, 0, 0, 1],
    [-2, 0, 0, 2, 2],
    [0, 0, 0, 2, 2],
    [0, 0, 0, 0, 2],
    [0, 1, 0, 0, 0],
    [0, 0, 1, 0, 0],
    [-2, 1, 0, 2, 2],
    [0, 0, 0, 2, 1],
    [0, 0, 1, 2, 2],
    [-2, -1, 0, 2, 2],
    [-2, 0, 1, 0, 0],
    [-2, 0, 0, 2, 1],
    [0, 0, -1, 2, 2],
    [2, 0, 0, 0, 0],
    [0, 0, 1, 0, 1],
    [2, 0, -1, 2, 2],
    [0, 0, -1, 0, 1],
    [0, 0, 1, 2, 1],
    [-2, 0, 2, 0, 0],
    [0, 0, -2, 2, 1],
    [2, 0, 0, 2, 2],
    [0, 0, 2, 2, 2],
    [0, 0, 2, 0, 0],
    [-2, 0, 1, 2, 2],
    [0, 0, 0, 2, 0],
    [-2, 0, 0, 2, 0],
    [0, 0, -1, 2, 1],
    [0, 2, 0, 0, 0],
    [2, 0, -1, 0, 1],
    [-2, 2, 0, 2, 2],
    [0, 1, 0, 0, 1],
    [-2, 0, 1, 0, 1],
    [0, -1, 0, 0, 1],
    [0, 0, 2, -2, 0],
    [2, 0, -1, 2, 1],
    [2, 0, 1, 2, 2],
    [0, 1, 0, 2, 2],
    [-2, 1, 1, 0, 0],
    [0, -1, 0, 2, 2],
    [2, 0, 0, 2, 1],
    [2, 0, 1, 0, 0],
    [-2, 0, 2, 2, 2],
    [-2, 0, 1, 2, 1],
    [2, 0, -2, 0, 1],
    [2, 0, 0, 0, 1],
    [0, -1, 1, 0, 0],
    [-2, -1, 0, 2, 1],
    [-2, 0, 0, 0, 1],
    [0, 0, 2, 2, 1],
    [-2, 0, 2, 0, 1],
    [-2, 1, 0, 2, 1],
    [0, 0, 1, -2, 0],
    [-1, 0, 1, 0, 0],
    [-2, 1, 0, 0, 0],
    [1, 0, 0, 0, 0],
    [0, 0, 1, 2, 0],
    [0, 0, -2, 2, 2],
    [-1, -1, 1, 0, 0],
    [0, 1, 1, 0, 0],
    [0, -1, 1, 2, 2],
    [2, -1, -1, 2, 2],
    [0, 0, 3, 2, 2],
    [2, -1, 0, 2, 2],
]
_NUTATION_PERIODIC_TERMS_DATA_PE = [
    [-171996, -174.2, 92025, 8.9],
    [-13187, -1.6, 5736, -3.1],
    [-2274, -0.2, 977, -0.5],
    [2062, 0.2, -895, 0.5],
    [1426, -3.4, 54, -0.1],
    [712, 0.1, -7, 0],
    [-517, 1.2, 224, -0.6],
    [-386, -0.4, 200, 0],
    [-301, 0, 129, -0.1],
    [217, -0.5, -95, 0.3],
    [-158, 0, 0, 0],
    [129, 0.1, -70, 0],
    [123, 0, -53, 0],
    [63, 0, 0, 0],
    [63, 0.1, -33, 0],
    [-59, 0, 26, 0],
    [-58, -0.1, 32, 0],
    [-51, 0, 27, 0],
    [48, 0, 0, 0],
    [46, 0, -24, 0],
    [-38, 0, 16, 0],
    [-31, 0, 13, 0],
    [29, 0, 0, 0],
    [29, 0, -12, 0],
    [26, 0, 0, 0],
    [-22, 0, 0, 0],
    [21, 0, -10, 0],
    [17, -0.1, 0, 0],
    [16, 0, -8, 0],
    [-16, 0.1, 7, 0],
    [-15, 0, 9, 0],
    [-13, 0, 7, 0],
    [-12, 0, 6, 0],
    [11, 0, 0, 0],
    [-10, 0, 5, 0],
    [-8, 0, 3, 0],
    [7, 0, -3, 0],
    [-7, 0, 0, 0],
    [-7, 0, 3, 0],
    [-7, 0, 3, 0],
    [6, 0, 0, 0],
    [6, 0, -3, 0],
    [6, 0, -3, 0],
    [-6, 0, 3, 0],
    [-6, 0, 3, 0],
    [5, 0, 0, 0],
    [-5, 0, 3, 0],
    [-5, 0, 3, 0],
    [-5, 0, 3, 0],
    [4, 0, 0, 0],
    [4, 0, 0, 0],
    [4, 0, 0, 0],
    [-4, 0, 0, 0],
    [-4, 0, 0, 0],
    [-4, 0, 0, 0],
    [3, 0, 0, 0],
    [-3, 0, 0, 0],
    [-3, 0, 0, 0],
    [-3, 0, 0, 0],
    [-3, 0, 0, 0],
    [-3, 0, 0, 0],
    [-3, 0, 0, 0],
    [-3, 0, 0, 0],
]

_NUTATION_PERIODIC_TERMS_Y = pd.DataFrame.from_records(
    _NUTATION_PERIODIC_TERMS_DATA_Y,
    columns=['Y0', 'Y1', 'Y2', 'Y3', 'Y4']
)
_NUTATION_PERIODIC_TERMS_PE = pd.DataFrame.from_records(
    _NUTATION_PERIODIC_TERMS_DATA_PE,
    columns=['a', 'b', 'c', 'd']
)

# ---------------------------------------------------------------------
#: from [ASH1997]_
_ASHRAE1997_ch29_table8 = pd.DataFrame.from_records(
    columns=["soy", "I_r", "eot_min", "declination_deg", "A", "B", "C"],
    data=[
        (1728000.,  1416, -11.2, -20.0, 1230, 0.142, 0.058),
        (4406400.,  1401, -13.9, -10.8, 1215, 0.144, 0.060),
        (6825600.,  1381, -7.5, 0.0, 1186, 0.156, 0.071),
        (9504000.,  1356, 1.1, 11.6, 1136, 0.180, 0.097),
        (12096000., 1336, 3.3, 20.0, 1104, 0.196, 0.121),
        (14774400., 1336, -1.4, 23.45, 1088, 0.205, 0.134),
        (17366400., 1336, -6.2, 20.6, 1085, 0.207, 0.136),
        (20044800., 1338, -2.4, 12.3, 1107, 0.201, 0.122),
        (22723200., 1359, 7.5, 0.0, 1151, 0.177, 0.092),
        (25315200., 1380, 15.4, -10.5, 1192, 0.160, 0.073),
        (27993600., 1405, 13.8, -19.8, 1221, 0.149, 0.063),
        (30585600, 1417, 1.6, -23.45, 1233, 0.142, 0.057),
    ],
)

# ---------------------------------------------------------------------


def _clear_sky_ashrae1997(time, lat, lon):
    r"""
    Return direct normal irradiance for a given time an position.

    :param time: time for which to calculate clear sky irradiance
    :type time: datetime64
    :param lat: position latitude
    :type lat: float
    :param lon: position longitude
    :type lon: float
    :return: direct normal irradiance in W/m²
    :rtype: float


    Calculates direct normal irradiance, or solar irradiance,
    :math:`E_{\{mathrm{DN}}}`
    after 1997 ASHRAE Handbook [ASH97]_ chapter 29, eqn (15):

    :math:`E_{\{mathrm{DN}}} ~=~ \frac{A}{exp \left( B / sin \beta\right)}`

    where
        :math:`A` = apparent solar irradiation at air mass m = 0 (Table 8)
        :math:`B` = atmospheric extinction coefficient (Table 8)

    """
    # sun elevation, azimuth in degrees
    ele, azi = fast_sun_position(time, lat, lon)

    # estimate seconds of year to interpolate coefficients
    seconds = (
            time - pd.Timestamp(year=time.year, month=1, day=1,
                                tz=time.tzinfo)
    ).total_seconds()
    coeffs = {x: np.interp(seconds,
                           _ASHRAE1997_ch29_table8['soy'],
                           _ASHRAE1997_ch29_table8[x])
              for x in _ASHRAE1997_ch29_table8.columns
              }

    # irradiance
    #
    if ele >= 0.:
        # sun above horizon
        beta = numpy.deg2rad(ele)
        e_dn = coeffs['A'] * np.exp(- coeffs['B'] / np.sin(beta))
    else:
        # sun below horizon
        e_dn = 0.

    return e_dn
# ---------------------------------------------------------------------


def _sfc_irrad_ashrae1997(time, lat, lon, heading, slant, albedo=None,
                          _debug_angles=False):
    r"""
    Return direct and diffuse irradiance to a surface
    on a clear day

    :param time: time for which to calculate clear sky irradiance
    :type time: datetime64
    :param lat: position latitude
    :type lat: float
    :param lon: position longitude
    :type lon: float
    :param heading: surface heading angle in degrees clockwise from north
    :type heading: float
    :param slant: surface slant angle in degrees upwards from horizontal
    :type slant: float
    :param albedo: soil surface albedo in the area seen by the surface
      in 1. Defaults to 0.15 (Arable land, grazing, mixed forests,
      according to [HeS1983]_.
    :type albedo: float
    :return: direct normal irradiance in W/m²
    :rtype: float


    Calculates direct normal irradiance, or solar irradiance,
    :math:`E_{\{mathrm{DN}}}`
    after 1997 ASHRAE Handbook [ASH97]_ chapter 29, eqn (15):

    :math:`E_{\{mathrm{DN}}} ~=~ \frac{A}{exp \left( B / sin \beta\right)}`

    where
        :math:`A` = apparent solar irradiation at air mass m = 0 (Table 8)
        :math:`B` = atmospheric extinction coefficient (Table 8)

    """
    if albedo is None:
        albedo = 0.15   # typical value for Central Europe

    # sun elevation, azimuth in degrees
    ele, azi = fast_sun_position(time, lat, lon)
    #  surface-solar azimuth in radians
    gamma = numpy.deg2rad(azi - heading)
    # sun elevation in radians
    beta = numpy.deg2rad(ele)
    # surface tilt in radians (0= sfc facing up, pi/2= sfc facing horizon)
    sigma = numpy.deg2rad(90 - slant)

    # estimate seconds of year to interpolate coefficients
    e_dn = _clear_sky_ashrae1997(time, lat, lon)

    # estimate seconds of year to interpolate coefficients
    seconds = (
            time - pd.Timestamp(year=time.year, month=1, day=1,
                                tz=time.tzinfo)
    ).total_seconds()
    coeffs = {x: np.interp(seconds,
                           _ASHRAE1997_ch29_table8['soy'],
                           _ASHRAE1997_ch29_table8[x])
              for x in _ASHRAE1997_ch29_table8.columns
              }

    # surface incident angle (cosine thereof)
    cos_theta = (np.cos(beta) * np.cos(gamma) * np.sin(sigma) +
                 np.sin(beta) * np.cos(sigma))

    if cos_theta > 0:
        # sun is visible for surface
        e_dir = cos_theta * e_dn
    else:
        e_dir = 0.

    # ratio of vertical/horizontal sky diffuse
    if cos_theta > -0.2:
        yps = 0.55 + 0.437 * cos_theta + 0.313 * cos_theta**2
    else:
        yps = 0.45
    # diffuse sky irradiance
    e_ds = coeffs['C'] * yps * e_dn
    # diffuse ground reflected irradiance
    e_dg = e_dn * (coeffs['C'] +
                   np.sin(ele) * albedo * (1 - np.cos(sigma))/2.
                   )
    # diffuse irradiance
    e_diff = e_ds + e_dg

    if _debug_angles:
        return (e_dir, e_diff,
                np.rad2deg(np.arccos(cos_theta)), np.rad2deg(gamma))

    return e_dir, e_diff


# ---------------------------------------------------------------------


class _spa_location():
    def __init__(self, time: dt.datetime, lat,
                 lon, ele=0., pp=None, tk=None):

        if time.tzinfo is None:
            # assume UTC time
            time = time.replace(tzinfo=pytz.utc)
        if time.tzinfo == pytz.utc:
            self.utc = time
        else:
            self.utc = time - time.utcoffset()
        self.timezone = (
            (time.utcoffset().seconds / 3600. + 12) %
            24) - 12
        self.delta_T = _clock_error_formula(time.year)
        self.ut = self.utc + pd.Timedelta(self.delta_T, unit='seconds')
        self.jd = _spa_julian_day(self.utc)
        self.jde = _spa_julian_ephemeris_day(self)
        self.jc = _spa_julian_century(self)
        self.jce = _spa_julian_ephemeris_century(self)
        self.jme = _spa_julian_ephemeris_millennium(self)
        self.lat = lat
        self.lon = lon
        self.ele = ele
        if pp is None:
            self.pp = 1013.25
        else:
            self.pp = pp
        if tk is None:
            self.tk = 288.15
        else:
            self.tk = tk

        self._H = np.nan
        self._R = np.nan
        self._alpha = np.nan
        self._delta_alpha = np.nan
        self._delta = np.nan
        self._epsilon = np.nan
        self._delta_eps = np.nan
        self._lamda = np.nan
        self._nu = np.nan
        self._delta_psi = np.nan

    def H(self):
        if np.isnan(self._H):
            self._H = _spa_observer_local_hour_angle(self)
        return self._H

    def R(self):
        if np.isnan(self._R):
            self._R = _spa_earth_heliocentric_radius(self)
        return self._R

    def alpha(self):
        if np.isnan(self._alpha):
            self._alpha = _spa_geocentric_right_ascension(self)
        return self._alpha

    def delta_alpha(self):
        if np.isnan(self._delta_alpha):
            self._delta_alpha = _spa_topocentric_sun_right_ascension_parallax(
                self)
        return self._delta_alpha

    def delta(self):
        if np.isnan(self._delta):
            self._delta = _spa_geocentric_declination(self)
        return self._delta

    def epsilon(self):
        if np.isnan(self._epsilon):
            self._epsilon = _spa_ecliptic_true_obliquity(self)
        return self._epsilon

    def delta_eps(self):
        if np.isnan(self._delta_eps):
            self._delta_psi, self._delta_eps = _spa_nutation_lon_and_obliquity(
                self)
        return self._delta_eps

    def lamda(self):
        if np.isnan(self._lamda):
            self._lamda = _spa_apparent_sun_longitude(self)
        return self._lamda

    def nu(self):
        if np.isnan(self._nu):
            self._nu = _spa_greenwich_apparent_sidereal_time(self)
        return self._nu

    def delta_psi(self):
        if np.isnan(self._delta_psi):
            self._delta_psi, self._delta_eps = _spa_nutation_lon_and_obliquity(
                self)
        return self._delta_psi


# ---------------------------------------------------------------------

def _spa_julian_day(time: dt.datetime):
    '''
    Calculate the Julian Day (JD)

    The Julian date starts on January 1, in the year - 4712 at 12:00:00 UT.
    The Julian Day (JD) is calculated using UT and the Julian Ephemeris Day
    (JDE) is calculated using TT. In the following steps, note that there is
    a 10-day gap between the Julian and Gregorian calendar where the Julian
    calendar ends on October 4, 1582 (JD = 2299160), and after 10-days the
    Gregorian calendar starts on October 15, 1582
    '''
    Y = time.year
    M = time.month
    D = time.day + (time.hour + (time.minute +
                                 time.second / 60.) / 60.) / 24.

    if M < 3:
        M = M + 12
        Y = Y - 1

    JD = (np.int64(365.25 * (Y + 4716.))
          + np.int64(30.6001 * (M + 1))
          + D - 1524.5
          )

    if JD > 2299160.0:
        a = np.int64(Y / 100)
        JD = JD + (2 - a + np.int64(a / 4.))

    logging.debug('Julian Day (JD): %f' % JD)
    return JD

# ----------------------------------------------------


def _spa_julian_ephemeris_day(loc):
    '''
    Calculate the Julian Ephemeris Day (JDE)
    '''
    jde = loc.jd + loc.delta_T / 86400.0
    logging.debug('Julian Ephemeris Day (JDE): %f' % jde)
    return jde

# ----------------------------------------------------


def _spa_julian_century(loc):
    '''
    Calculate the Julian century (JC) for the 2000 standard epoch
    '''
    jc = (loc.jd - 2451545.) / 36525.
    logging.debug('Julian century (JC): %f' % jc)
    return jc

# ----------------------------------------------------


def _spa_julian_ephemeris_century(loc):
    '''
    Calculate the Julian Ephemeris Century (JCE) for the 2000 standard epoch
    '''
    jce = (loc.jde - 2451545.) / 36525.
    logging.debug('Julian Ephemeris Century (JCE): %f' % jce)
    return jce

# ----------------------------------------------------


def _spa_julian_ephemeris_millennium(loc):
    '''
    Calculate the Julian Ephemeris Millennium (JME) for the 2000 standard epoch
    '''
    jme = (loc.jce / 10.)
    logging.debug('Julian Ephemeris Millennium (JME): %f' % jme)
    return jme

# ----------------------------------------------------


def __spa_earth_periodic_term(eptx, jme):
    '''
    Calculate one of the Earth periodic terms from table 4.2 in _[ReA2008]
    '''
    xi = eptx['A'] * np.cos(
        eptx['B'] + eptx['C'] * jme)
    x = xi.sum()
    return x

# ----------------------------------------------------


def _spa_earth_heliocentric_longitude(loc):
    '''
    Calculate the Earth heliocentric longitude, L (in degrees),
    limited to 0°..360°
    '''
    Li = []
    for i, eptl in enumerate(_EARTH_PERIODIC_TERMS_L):
        Li.append(__spa_earth_periodic_term(eptl, loc.jme))
    Lp = np.polynomial.Polynomial(Li)
    L_rad = Lp(loc.jme) / 1.E8

    L = np.rad2deg(L_rad) % 360.
    logging.debug('Earth heliocentric longitude, L: %f' % L)
    return L

# ----------------------------------------------------


def _spa_earth_heliocentric_latitude(loc):
    '''
    Calculate the Earth heliocentric Latitude, B (in degrees),
    limited to 0°..360°
    '''
    Bi = []
    for i, eptb in enumerate(_EARTH_PERIODIC_TERMS_B):
        Bi.append(__spa_earth_periodic_term(eptb, loc.jme))
    Bp = np.polynomial.Polynomial(Bi)
    B_rad = Bp(loc.jme) / 1.E8
    B = np.rad2deg(B_rad) % 360.
    logging.debug('Earth heliocentric Latitude, B: %f' % B)
    return B

# ----------------------------------------------------


def _spa_earth_heliocentric_radius(loc):
    '''
    Calculate the Earth radius vector, R (in Astronomical Units, AU)
    '''
    Ri = []
    for i, eptr in enumerate(_EARTH_PERIODIC_TERMS_R):
        Ri.append(__spa_earth_periodic_term(eptr, loc.jme))
    Rp = np.polynomial.Polynomial(Ri)
    R = Rp(loc.jme) / 1.E8
    logging.debug('Earth radius vector, R: %f' % R)
    return R

# ----------------------------------------------------


def _spa_geocentric_longitude(loc):
    '''
    Calculate the geocentric longitude, theta (in degrees)
    '''
    L = _spa_earth_heliocentric_longitude(loc)
    theta = L + 180.
    if theta >= 360.:
        theta = theta - 360.
    logging.debug('geocentric longitude, theta: %f' % theta)
    return theta

# ----------------------------------------------------


def _spa_geocentric_latitude(loc):
    '''
    Calculate the geocentric latitude, beta (in degrees)
    '''
    B = _spa_earth_heliocentric_latitude(loc)
    beta = -B
    logging.debug('geocentric latitude, beta: %f' % beta)
    return beta

# ----------------------------------------------------


def _spa_mean_elongation_moon_sun(loc):
    '''
    Calculate the mean elongation of the moon from the sun, X0 (in degrees)
    '''
    X0p = np.polynomial.Polynomial(
        [297.85036, 445267.11148, -0.0019142, 1.0 / 189474.0])
    return X0p(loc.jce)

# ----------------------------------------------------


def _spa_mean_anomaly_sun(loc):
    '''
    Calculate the mean anomaly of the sun (Earth), X1 (in degrees)
    '''
    X1p = np.polynomial.Polynomial(
        [357.52772, 35999.05034, -0.0001603, -1.0 / 300000.0])
    X1 = X1p(loc.jce)
    logging.debug('mean anomaly of the sun, X1: %f' % X1)
    return X1

# ----------------------------------------------------


def _spa_mean_anomaly_moon(loc):
    '''
    Calculate the mean anomaly of the moon, X2 (in degrees)
    '''
    X2p = np.polynomial.Polynomial(
        [134.96298, 477198.867398, 0.0086972, 1.0 / 56250.0])
    X2 = X2p(loc.jce)
    logging.debug('mean anomaly of the moon, X2: %f' % X2)
    return X2

# ----------------------------------------------------


def _spa_argument_latitude_moon(loc):
    '''
    Calculate the moon’s argument of latitude, X3 (in degrees)
    '''
    X3p = np.polynomial.Polynomial(
        [93.27191, 483202.017538, -0.0036825, 1.0 / 327270.0])
    X3 = X3p(loc.jce)
    logging.debug('moon’s argument of latitude, X3: %f' % X3)
    return X3

# ----------------------------------------------------


def _spa_ascending_longitude_moon(loc):
    '''
    Calculate the longitude of the ascending node of the
    moon’s mean orbit on the ecliptic, measured from the
    mean equinox of the date, X4 (in degrees)
    '''
    X4p = np.polynomial.Polynomial(
        [125.04452, -1934.136261, 0.0020708, 1.0 / 450000.0])
    X4 = X4p(loc.jce)
    logging.debug('ascending longitude moon, X4: %f' % X4)
    return X4

# ----------------------------------------------------


def _spa_argument_vector_X(loc):
    X = []
    X.append(_spa_mean_elongation_moon_sun(loc))
    X.append(_spa_mean_anomaly_sun(loc))
    X.append(_spa_mean_anomaly_moon(loc))
    X.append(_spa_argument_latitude_moon(loc))
    X.append(_spa_ascending_longitude_moon(loc))
    return np.array(X)

# ----------------------------------------------------


def __spa_delta_PE(i, X, loc):
    '''
    For each row i in Table A4.3,
    calculate the terms delta psi and delta epsilon
    '''

    ai = _NUTATION_PERIODIC_TERMS_PE['a'][i]
    bi = _NUTATION_PERIODIC_TERMS_PE['b'][i]
    ci = _NUTATION_PERIODIC_TERMS_PE['c'][i]
    di = _NUTATION_PERIODIC_TERMS_PE['d'][i]
    xy = X * _NUTATION_PERIODIC_TERMS_Y.iloc[i]

    delta_psi_i = (ai + bi * loc.jce) * np.sin(xy.sum())
    delta_eps_i = (ci + di * loc.jce) * np.cos(xy.sum())

    return delta_psi_i, delta_eps_i

# ----------------------------------------------------

# def _spa_nutation_longitude(loc):
#    '''
#    Calculate the nutation in longitude, delta psi (in degrees),
#    '''
#    sum_psi = 0
#    X = _spa_argument_vector_X(loc)
#    for i in _NUTATION_PERIODIC_TERMS_Y.index:
#        delta_psi_i, _ = __spa_delta_PE(i, X, loc)
#        sum_psi = sum_psi + delta_psi_i
#    delta_psi = sum_psi / 36000000.
#    logging.debug('Nutation in longitude, delta psi: %f' % delta_psi)
#    return delta_psi
#
# ----------------------------------------------------
#
# def _spa_nutation_obliquity(loc):
#    '''
#    Calculate the nutation in obliquity, delta_eps (in degrees)
#    '''
#    sum_eps = 0
#    X = _spa_argument_vector_X(loc)
#    for i in _NUTATION_PERIODIC_TERMS_Y.index:
#        _, delta_eps_i = __spa_delta_PE(i, X, loc)
#        sum_eps = sum_eps + delta_eps_i
#
#    delta_eps = sum_eps / 36000000.
#    logging.debug('Nutation in obliquity, delta_eps: %f' % delta_eps)
#    return delta_eps


def _spa_nutation_lon_and_obliquity(loc):
    '''
    Calculate the nutation in longitude, delta psi (in degrees),
    Calculate the nutation in obliquity, delta_eps (in degrees)
    '''
    sum_psi = 0
    sum_eps = 0
    X = _spa_argument_vector_X(loc)
    for i in _NUTATION_PERIODIC_TERMS_Y.index:
        delta_psi_i, delta_eps_i = __spa_delta_PE(i, X, loc)
        sum_psi = sum_psi + delta_psi_i
        sum_eps = sum_eps + delta_eps_i
    delta_psi = sum_psi / 36000000.
    delta_eps = sum_eps / 36000000.
    logging.debug('Nutation in longitude, delta psi: %f' % delta_psi)
    logging.debug('Nutation in obliquity, delta_eps: %f' % delta_eps)
    return delta_psi, delta_eps

# ----------------------------------------------------


def _spa_ecliptic_true_obliquity(loc):
    '''
    Calculate the true obliquity of the ecliptic, epsilon (in degrees)
    '''
    delta_eps = loc.delta_eps()
    # Calculate the mean obliquity of the ecliptic, g0 (in arc seconds)
    u = loc.jme / 10.
    eps0 = np.polynomial.Polynomial(
        [84381.448, -4680.93, -1.55, 1999.25, -51.38, -249.67,
         -39.05, 7.12, 27.87, 5.79, 2.45])
    # Calculate the true obliquity of the ecliptic, g (in degrees)
    epsilon = eps0(u) / 3600. + delta_eps
    logging.debug('Ecliptic true obliquity: %f' % epsilon)
    return epsilon

# ----------------------------------------------------


def _spa_aberration_correction(loc):
    '''
    Calculate the aberration correction, delta tau (in degrees)
    '''
    R = loc.R()
    delta_tau = -20.4898 / (3600.0 * R)
    logging.debug('Aberration correction, delta tau: %f' % delta_tau)
    return delta_tau

# ----------------------------------------------------


def _spa_apparent_sun_longitude(loc):
    '''
    Calculate the apparent sun longitude, lamda (in degrees)
    '''
    theta = _spa_geocentric_longitude(loc)
    delta_psi = loc.delta_psi()
    delta_tau = _spa_aberration_correction(loc)
    lamda = theta + delta_psi + delta_tau
    logging.debug('Apparent sun longitude, lamda: %f' % lamda)
    return lamda

# ----------------------------------------------------


def _spa_greenwich_apparent_sidereal_time(loc):
    '''
    Calculate the apparent sidereal time at Greenwich
    at any given time, nu (in degrees)
    '''
    epsilon = loc.epsilon()
    delta_psi = loc.delta_psi()
    # Calculate the mean sidereal time at Greenwich, <0 (in degrees)
    nu0 = (280.46061837 + 360.98564736629 * (loc.jd - 2451545.) +
           loc.jc * loc.jc * (0.000387933 - loc.jc / 38710000.))
    nu = (nu0 % 360.) + delta_psi * np.cos(np.deg2rad(epsilon))
    logging.debug('Greenwich apparent sidereal time: %f' % nu)
    return nu

# ----------------------------------------------------


def _spa_geocentric_right_ascension(loc):
    '''
    Calculate the geocentric sun right ascension, alpha (in degrees)
    '''
    lamda_rad = np.deg2rad(loc.lamda())
    epsilon_rad = np.deg2rad(loc.epsilon())
    beta_rad = np.deg2rad(_spa_geocentric_latitude(loc))
    alpha_rad = np.arctan2(np.sin(lamda_rad) * np.cos(epsilon_rad) -
                           np.tan(beta_rad) * np.sin(epsilon_rad),
                           np.cos(lamda_rad))
    alpha = np.rad2deg(alpha_rad) % 360.
    logging.debug('Geocentric right ascension, alpha: %f' % alpha)
    return alpha

# ----------------------------------------------------


def _spa_geocentric_declination(loc):
    '''
    Calculate the geocentric sun declination, delta (in degrees)
    '''
    beta_rad = np.deg2rad(_spa_geocentric_latitude(loc))
    epsilon_rad = np.deg2rad(loc.epsilon())
    lamda_rad = np.deg2rad(loc.lamda())
    delta_rad = np.arcsin(np.sin(beta_rad) * np.cos(epsilon_rad) +
                          np.cos(beta_rad) * np.sin(epsilon_rad) *
                          np.sin(lamda_rad))
    delta = np.rad2deg(delta_rad)
    logging.debug('Geocentric sun declination, delta: %f' % delta)
    return delta

# ----------------------------------------------------


def _spa_observer_local_hour_angle(loc):
    '''
    Calculate the observer local hour angle, H (in degrees)
    '''
    nu = loc.nu()
    alpha = loc.alpha()
    H = (nu + loc.lon - alpha) % 360.
    logging.debug('observer local hour angle, H: %f' % H)
    return H

# ----------------------------------------------------


def _spa_sun_equatorial_horizontal_parallax(loc):
    lat_rad = np.deg2rad(loc.lat)
    # Calculate the equatorial horizontal parallax of the sun, xi (in
    # degrees)
    R = loc.R()
    xi = 8.794 / (3600.0 * R)

    # Calculate the term u (in radians),
    u = np.arctan(0.99664719 * np.tan(lat_rad))

    # Calculate the term x
    x = np.cos(u) + loc.ele * np.cos(lat_rad) / 6378140.

    # Calculate the term y
    y = 0.99664719 * np.sin(u) + loc.ele * np.sin(lat_rad) / 6378140.

    return xi, x, y

# ----------------------------------------------------


def _spa_topocentric_sun_right_ascension_parallax(loc):
    '''
    Calculate the parallax in the sun right ascension, )" (in degrees)
    '''
    # Calculate the equatorial horizontal parallax of the sun, xi (in
    # degrees)
    xi, x, _ = _spa_sun_equatorial_horizontal_parallax(loc)
    xi_rad = np.deg2rad(xi)

    # Calculate the parallax in the sun right ascension, delta_alpha (in
    # degrees)
    h_rad = np.deg2rad(loc.H())
    delta_rad = np.deg2rad(loc.delta())
    delta_alpha_rad = np.arctan2((-x * np.sin(xi_rad) * np.sin(h_rad)),
                                 (np.cos(delta_rad) - x * np.sin(xi_rad) *
                                  np.cos(h_rad)))
    delta_alpha = np.rad2deg(delta_alpha_rad)
    logging.debug(
        'parallax in sun right ascension, delta_alpha: %f' % delta_alpha)
    return delta_alpha

# ----------------------------------------------------


def _spa_topocentric_sun_right_ascension(loc):
    '''
    Calculate the topocentric sun right ascension alpha_prime (in degrees)
    Calculate the parallax in the sun right ascension, )" (in degrees)
    '''
    alpha = loc.alpha()
    delta_alpha = loc.delta_alpha()
    alpha_prime = alpha + delta_alpha
    logging.debug(
        'topocentric sun right ascension, alpha_prime: %f' % alpha_prime)
    return alpha_prime

# ----------------------------------------------------


def _spa_topocentric_sun_declination(loc):
    '''
    Calculate the topocentric sun declination, delta_prime (in degrees)
    '''
    delta_rad = np.deg2rad(loc.delta())
    xi, x, y = _spa_sun_equatorial_horizontal_parallax(loc)
    xi_rad = np.deg2rad(xi)
    h_rad = np.deg2rad(loc.H())
    delta_alpha_rad = np.deg2rad(loc.delta_alpha())
    delta_prime_rad = np.arctan2(
        (np.sin(delta_rad) - y * np.sin(xi_rad)) * np.cos(delta_alpha_rad),
        np.cos(delta_rad) - x * np.sin(xi_rad) * np.cos(h_rad)
    )
    delta_prime = np.rad2deg(delta_prime_rad)
    logging.debug(
        'topocentric sun declination, delta_prime: %f' %
        delta_prime)
    return delta_prime

# ----------------------------------------------------


def _spa_topocentric_local_hour_angle(loc):
    '''
    Calculate the topocentric local hour angle, H’ (in degrees)
    '''
    H = loc.H()
    delta_alpha = loc.delta_alpha()
    H_prime = H - delta_alpha
    return H_prime

# ----------------------------------------------------


def _spa_topocentric_elevation_angle(loc):
    '''
    Calculate the topocentric elevation angle, e0 (in degrees):
    tk in K
    pp in hPa
    '''
    lat_rad = np.deg2rad(loc.lat)
    delta_prime_rad = np.deg2rad(_spa_topocentric_sun_declination(loc))
    h_prime_rad = np.deg2rad(_spa_topocentric_local_hour_angle(loc))
    # Calculate the topocentric elevation angle without atmospheric refraction
    # correction, e0 (in degrees),
    e0 = np.rad2deg(np.arcsin(np.sin(lat_rad) * np.sin(delta_prime_rad) +
                              np.cos(lat_rad) * np.cos(delta_prime_rad) *
                              np.cos(h_prime_rad)))
    logging.debug('topocentric_elevation_angle e0: %f' % e0)
    return e0

# ----------------------------------------------------


def __spa_atmospheric_refraction_correction(e0, pp, tk):
    '''
    Calculate the atmospheric refraction correction, delta_e (in degrees)
    e0 in deg
    pp in hPa
    tk in K
    '''
    if pp is None or pd.isna(pp) or tk is None or pd.isna(tk):
        del_e = 0.
    # Calculate the atmospheric refraction correction, delta e (in
    # degrees),
    else:
        if e0 >= -1 * _SUN_RADIUS:
            # sun (partly) above horizon
            del_e = ((pp / 1010.0) * (283.0 / tk) *
                     1.02 /
                     (60.0 * np.tan(np.deg2rad(e0 + 10.3 / (e0 + 5.11))))
                     )
        else:
            del_e = 0.
    return del_e

# ----------------------------------------------------


def _spa_topocentric_elevation_angle_corrected(loc):
    e0 = _spa_topocentric_elevation_angle(loc)
    delta_e = __spa_atmospheric_refraction_correction(
        e0, loc.pp, loc.tk)
    e = e0 + delta_e
    logging.debug('topocentric_elevation_angle_corrected, e: %f' % e)
    return e

# ----------------------------------------------------


def _spa_topocentric_zenith_angle(loc):
    e = _spa_topocentric_elevation_angle_corrected(loc)
    return 90.0 - e

# ----------------------------------------------------


def _spa_topocentric_azimuth_angle_astro(loc):
    '''
    Calculate the topocentric astronomers azimuth angle, Gamma (in degrees)
    '''
    h_prime_rad = np.deg2rad(_spa_topocentric_local_hour_angle(loc))
    lat_rad = np.deg2rad(loc.lat)
    delta_prime_rad = np.deg2rad(
        _spa_topocentric_sun_declination(loc))

    Gamma = np.rad2deg(
        np.arctan2(
            np.sin(h_prime_rad),
            np.cos(h_prime_rad) * np.sin(lat_rad)
            - np.tan(delta_prime_rad) * np.cos(lat_rad)
        )) % 360.
    return Gamma

# ----------------------------------------------------


def _spa_topocentric_azimuth_angle(loc):
    Phi = (_spa_topocentric_azimuth_angle_astro(loc) + 180.) % 360.
    return Phi

# ---------------------------------------------------------------------


def _spa_sun_mean_longitude(loc):
    '''
    Calculate the sun’s mean longitude (in degrees)
    '''
    M = (280.4664567 + 360007.6982779 * loc.jme + 0.03032028 * loc.jme**2
         + 1. / 49931. * loc.jme**3 - 1 / 15300. * loc.jme**4
         - 1 / 2000000. * loc.jme**5)
    return M % 360.

# ----------------------------------------------------


def _spa_equation_of_time(loc):
    '''
    The Equation of Time, E, is the difference between solar apparent
    and mean time (in minutes of time)
    '''
    alpha = loc.alpha()
    delta_psi = loc.delta_psi()
    epsilon = loc.epsilon()
    M = _spa_sun_mean_longitude(loc)
    E_deg = (
        M -
        0.0057183 -
        alpha +
        delta_psi *
        np.cos(
            np.deg2rad(epsilon)))
    E = 4. * (E_deg % 360.)
    return E

# ----------------------------------------------------


def _approx_sun_transit_time(alpha_0, lon, nu):
    '''
    Calculate the approximate sun transit time, m 0 , in fraction of day
    '''
    return (alpha_0 - lon - nu) / 360.0

# ----------------------------------------------------


def __spa_sun_hour_angle_at_rise_set(lat, delta_0, h0_prime):
    lat_rad = np.deg2rad(lat)
    delta_0_rad = np.deg2rad(delta_0)
    h0_rad = np.deg2rad(h0_prime)
    frac = ((np.sin(h0_rad) - np.sin(lat_rad) * np.sin(delta_0_rad))
            / (np.cos(lat_rad) * np.cos(delta_0_rad)))
    if np.fabs(frac) <= 1:

        h0 = np.rad2deg(np.arccos(frac)) % 180
    else:  # sun is always above or below the horizon for that day
        h0 = None
    return h0

# ----------------------------------------------------


def __spa_approx_sun_rise_and_set(m0, h0):
    m_rise = (m0 - h0 / 360.) % 1.
    m_trans = m0 % 1.
    m_set = (m0 + h0 / 360.) % 1.

    return [m_rise, m_trans, m_set]

# ----------------------------------------------------


def __spa_alpha_delta_prime(ad, n):
    a = ad[1] - ad[0]
    b = ad[2] - ad[1]
    if np.fabs(a) >= 2.:
        a = a % 1.
    if np.fabs(b) >= 2.:
        b = b % 1.
    return ad[1] + n * (a + b + (b - a) * n) / 2.

# ----------------------------------------------------


def _rts_sun_altitude(lat, delta_prime, h_prime):
    lat_rad = np.deg2rad(lat)
    delta_prime_rad = np.deg2rad(delta_prime)
    h_prime_rad = np.deg2rad(h_prime)
    alt_rad = np.arcsin(np.sin(lat_rad) * np.sin(delta_prime_rad) +
                        np.cos(lat_rad) * np.cos(delta_prime_rad)
                        * np.cos(h_prime_rad))
    return np.rad2deg(alt_rad)

# ----------------------------------------------------


def __spa_rise_or_set(m, h, delta_prime, lat, h_prime, h0_prime):
    when = m + (h - h0_prime) / (360. * np.cos(np.deg2rad(delta_prime))
                                 * np.cos(np.deg2rad(lat))
                                 * np.sin(np.deg2rad(h_prime)))
    return when

# ----------------------------------------------------


def _spa_pm180(angle):
    angle = angle % 360
    if angle > 180:
        angle = angle - 360
    return angle


def _spa_day_to_hr(dayfrac, timezone):
    hr = 24.0 * ((dayfrac + timezone / 24.0) % 1.)
    return hr


def _spa_rise_transit_set(time: dt.datetime, lat, lon,
                          ele, pp, tk, approx=False):
    '''
    Calculate Sunrise, Sun Transit, and Sunset as hours of day
    This is the internal function accepting scalars only
    '''
    #
    h0_prime = -1. * (_SUN_RADIUS + _ATMOS_REFRAC)
    loc = _spa_location(time, lat, lon, ele, pp, tk)

    # location but with time 00 UT
    utc0 = dt.datetime(loc.utc.year, loc.utc.month, loc.utc.day,
                       0, 0, 0, 0, pytz.utc)
    loc0 = _spa_location(utc0, lat, lon, ele, pp, tk)

    # calculate_geocentric_sun_right_ascension_and_declination( &
    # sun_rts)
    nu = loc0.nu()

    utc_rts = [loc0.utc + pd.Timedelta(x, unit='days')
               for x in [-1, 0, 1]]
    loc_rts = [_spa_location(x, lat, lon, ele, pp, tk) for x in utc_rts]
    alpha_rts = [x.alpha() for x in loc_rts]
    delta_rts = [x.delta() for x in loc_rts]
    logging.debug('alpha_rts: %s' % repr(alpha_rts))
    logging.debug('delta_rts: %s' % repr(alpha_rts))

    m0 = _approx_sun_transit_time(alpha_rts[1], loc.lon, nu)
    logging.debug('approx Sun transit time, m0: %f' % m0)
    h0 = __spa_sun_hour_angle_at_rise_set(
        loc0.lat, delta_rts[1], h0_prime)
    logging.debug('Sunset hour angle, h0: %f' % h0)

    if h0 >= 0:
        m_rts = __spa_approx_sun_rise_and_set(m0, h0)
        logging.debug('m_rts    : %s' % repr(m_rts))

        if approx:
            sunrise = _spa_day_to_hr(m_rts[0], loc.timezone)
            transit = _spa_day_to_hr(m_rts[1], loc.timezone)
            sunset = _spa_day_to_hr(m_rts[2], loc.timezone)
        else:
            nu_rts = [nu + 360.985647 * x for x in m_rts]
            n_rts = [x + loc0.delta_T / 86400.0 for x in m_rts]
            logging.debug('nu_rts   : %s' % repr(nu_rts))
            logging.debug('n_rts    : %s' % repr(n_rts))
            alpha_prime_rts = [__spa_alpha_delta_prime(
                alpha_rts, x) for x in n_rts]
            delta_prime_rts = [__spa_alpha_delta_prime(
                delta_rts, x) for x in n_rts]
            logging.debug('alpha_prime_rts: %s' % repr(alpha_prime_rts))
            logging.debug('delta_prime_rts: %s' % repr(delta_prime_rts))
            h_prime_rts = [_spa_pm180(x + loc.lon - y)
                           for x, y in zip(nu_rts, alpha_prime_rts)]
            h_rts = [_rts_sun_altitude(loc0.lat, x, y)
                     for x, y in zip(delta_prime_rts, h_prime_rts)]
            logging.debug('h_prime_rts: %s' % repr(h_prime_rts))
            logging.debug('h_rts      : %s' % repr(h_rts))

            transit = _spa_day_to_hr(
                m_rts[1] - h_prime_rts[1] / 360., loc.timezone)
            sunrise = _spa_day_to_hr(__spa_rise_or_set(
                m_rts[0], h_rts[0], delta_prime_rts[0],
                loc0.lat, h_prime_rts[0], h0_prime), loc.timezone)
            sunset = _spa_day_to_hr(__spa_rise_or_set(
                m_rts[2], h_rts[2], delta_prime_rts[2],
                loc0.lat, h_prime_rts[2], h0_prime), loc.timezone)
    else:
        sunrise = None
        transit = None
        sunset = None
    return sunrise, transit, sunset

# ----------------------------------------------------


def _spa_position(time, lat, lon, ele=0, pp=None, tk=None):
    '''
    Calculate apparent position of the sun at time, postion lat,lon
    This is the internal function accepting scalars only
    '''
    loc = _spa_location(time, lat, lon, ele, pp, tk)
    zenith = _spa_topocentric_zenith_angle(loc)
    azimuth_astro = _spa_topocentric_azimuth_angle_astro(loc)
    azimuth = _spa_topocentric_azimuth_angle(loc)

    return zenith, azimuth, azimuth_astro


# ======== end of solar position functions (internal) ==================

def clock_error(time):
    '''
    Calculate historical and future values of the earth's clock error
    delta Tt using Polynomial expressions  for years 2000 BCE to 3000 CE,
    taken from Eclipse Predictions by Fred Espenak from NASA's GSFC
    [Esp2006]_.

    :param time: friction velocity in m/s (datetime).

    :return: clock error in seconds
    :rtype: float
    '''
    time_i = pd.to_datetime(time)
    if pd.api.types.is_scalar(time_i):
        res = _clock_error_formula(time_i.year)
    else:
        res = pd.Series([_clock_error_formula(x.year) for x in time_i])
    return res

# ----------------------------------------------------


def clear_sky_direct_normal(time, lat, lon, model="ashrae1997"):
    """
    Return direct normal irradiance on a clear day

    :param time: time for which to calculate clear sky irradiance
    :type time: datetime64
    :param lat: position latitude
    :type lat: float
    :param lon: position longitude
    :type lon: float
    :param model: name of the irradiance model to use.
      Currentlity implemented:
      - ``ashrae1997`` from [ASH1997]_

      Defaults to ``ashrae1997``.
    :type model: str
    :return: direct normal irradiance in W/m²
    :rtype: float
    """

    if model == "ashrae1997":
        return _clear_sky_ashrae1997(time, lat, lon)
    else:
        raise ValueError("unknown model: %s" % model)

# ----------------------------------------------------


def shortwave_incoming(time, lat, lon, heading, slant,
                       albedo=None, model="ashrae1997"):
    """
    Return direct and diffuse irradiance to a surface that is
    oriented in the direction described by heading and slant.

    :param time: time for which to calculate clear sky irradiance
    :type time: datetime64
    :param lat: position latitude
    :type lat: float
    :param lon: position longitude
    :type lon: float
    :param model: name of the irradiance model to use.
      Currentlity implemented:
      - ``ashrae1997`` from [ASH1997]_

      Defaults to ``ashrae1997``.
    :type model: str
    :return: direct normal irradiance in W/m²
    :rtype: float
    """

    if model == "ashrae1997":
        return _sfc_irrad_ashrae1997(time, lat, lon,
                                     heading, slant, albedo)
    else:
        raise ValueError("unknown model: %s" % model)

# ----------------------------------------------------


def longwave_incoming(t_k: float, e=0.,
                      model="angstrom"):
    r"""
    Calculate the clear-sky logwave downwelling radiation
    (counter radiation) using a formula form a choice of authors

    :param t_k: air temperature in K
    :type t_k: float
    :param e: water wapor partial pressure in Pa
    :type e: float
    :param model: name of the function to use
    :type model: str
    :return: longwave downwelling radiation in W/m
    :rtype: float

    angstrom:
        Formula by [Ang1916]_ depending on air temperature and humidity:

        :math:`L_{\mathrm{down}} = `
        :math:`0.434 - 0.158 \times 10^{-0.071 \rho} \frac{T_k^4}{293^4}`

        where :math:`\rho` is water vapor pressure in mmHg
        and :math:`T_k` is temperature in K

    swinbank:
        Formula by [Swi1963]_ depending on air temperature alone:

        :math:`L_{\mathrm{down}} = 5.31 \times 10^{-14} \, T^6`

        where :math:`T_k` is temperature in K
        and where :math:`L_{\mathrm{down}}` is in mW/cm²
        (The function returns the result is in W/m², though)

    brutsaert
        Formula by [Bru1975]_ depending on air temperature and humidity:

        :math:`L_{\mathrm{down}} = `
        :math:`\sigma T_k^4 ~ 1.24 {\frac{e}{T_k}}^{\frac{1}{7}}`

        where :math:`e` is water vapor pressure in Pa
        and :math:`T_k` is temperature in K

    aubinet
        Formula by [Aub1994]_ depending on humidity:

        :math:`L_{\mathrm{down}} =  \sigma T_s^4` with
        :math:`T_s = 147 + 18.2 log(e)`

        where :math:`e` is water vapor pressure in Pa

    """
    if model == "angstrom":
        rho = pa2mmhg(e)
        # ANGSTRÖM, A. 1916: Über die Gegenstrahlung der Atmosphäre
        # (On the counter-radiation of the atmosphere). - Meteorol. Z. 33,
        # 529-538 (translated and edited by VOLKEN, E., S. BRÖNNIMANN,
        # R. PHILIPONA). – Meteorol. Z. 22 (2013),
        # 761–769 (published online January 2014).
        l_down = ((0.434 - 0.158 * 10 ** (- 0.071 * rho)) *
                  (t_k ** 4) / (293. ** 4))
    elif model == "swinbank":
        # W. C. Swinbank, “Long-wave radiation from clear skies,”
        # Quarterly Journal of the Royal Meteorological Society,
        # vol. 89, no. 381, pp. 339–348, 1963, doi: 10.1002/qj.49708938105.
        # there: 5.31E-14 mW/cm² * T^6
        l_down = sigma * 9.36E-6 * t_k ** 6
    elif model == "brutsaert":
        # W. Brutsaert, “On a derivable formula for long‐wave radiation
        # from clear skies,” Water Resources Research,
        # vol. 11, no. 5, pp. 742–744, Oct. 1975,
        # doi: 10.1029/wr011i005p00742.
        e_hPa = e / 100.
        l_down = (sigma * t_k ** 4) * 1.24 * (e_hPa / t_k) ** (1. / 7.)
    elif model == "aubinet":
        #  M. Aubinet, “Longwave sky radiation parametrizations,”
        #  Solar Energy, vol. 53, no. 2, pp. 147–154, Aug. 1994,
        #  doi: 10.1016/0038-092x(94)90475-8.
        #  one variable model that hast the same funtional
        #  form as the recommended as three-variable model
        ts = 147. + 18.2 * np.log(e)
        l_down = (sigma * ts ** 4)
    else:
        raise ValueError("unknown method: %s" % model)
    return l_down

# ----------------------------------------------------


def spa_rise_transit_set(time, lat, lon, ele=0, pp=None, tk=None):
    '''
    Calculate Sunrise, Sun Transit, and Sunset as hours of day
    according to appendix A.2 of [Esp2006]_.

    :param time: friction velocity in m/s (datetime).
    :param lat: position latitude in degrees (float).
    :param lon: position longtude in degrees (float).
    :param ele: (optional) position elevation in m (float).
        Defaults to 0 m, if missing
    :param pp: (optional) mean atmospheric pressure at position in hPa (float).
      None causes atmospheric aberration to be neglected. Defaults to None.
    :param pp: (optional) mean air temperature at position in K (float).
      None causes atmospheric aberration to be neglected. Defaults to None.

    :return: sunrise, transit, sunset in hours of day (with fraction)
    :rtype: float
    '''
    time_i = pd.to_datetime(time)
    lat = _check('lat', lat, "float", ge=-90, le=90)
    lon = _check('lon', lon, "float", ge=-180, le=360)
    pp = _check('pp', pp, "float", ge=20, le=1100, none=True)
    tk = _check('tk', tk, "float", ge=150, le=360, none=True)
    if pd.api.types.is_scalar(time_i):
        res = _spa_rise_transit_set(time_i, lat, lon, ele, pp, tk)
    else:
        # expand scalar values to series, if provided as scalars
        lat_i = _expand_to_series_like(lat, time_i)
        lon_i = _expand_to_series_like(lon, time_i)
        ele_i = _expand_to_series_like(ele, time_i)
        pp_i = _expand_to_series_like(pp, time_i)
        tk_i = _expand_to_series_like(tk, time_i)
        # make empty output series
        sunrise = _expand_to_series_like(np.nan, time_i)
        transit = _expand_to_series_like(np.nan, time_i)
        sunset = _expand_to_series_like(np.nan, time_i)
        # calculate
        for i, arg in enumerate(
                zip(time_i, lat_i, lon_i, ele_i, pp_i, tk_i)):
            logging.info(
                'sun_rise_transit_set %i/%i' %
                (i, len(time_i)))
            sunrise[i], transit[i], sunset[i] = _spa_rise_transit_set(
                *arg)
        res = (sunrise, transit, sunset)
    return res

# ----------------------------------------------------


def fast_rise_transit_set(time, lat, lon):
    '''
    Calculate Sunrise, Sun Transit, and Sunset as hours of day
    according to [Mee1998]_,
    as described under ``https://gml.noaa.gov/grad/solcalc/calcdetails.html``.
    Caluculation is only valid for dates years 1901 and 2099
    and is less accurate than `sun_rise_transit_set`,
    but runs much (400 times) faster.

    :param time: friction velocity in m/s (datetime).
    :param lat: position latitude in degrees (float).
    :param lon: position longtude in degrees (float).
    :return: sunrise, transit, sunset in hours of day (with fraction)
    :rtype: float
    '''
    time = _check("time", time, "datetime",
                  ge=pd.to_datetime("1901-01-01", utc=True),
                  lt=pd.to_datetime("2100-01-01", utc=True))
    lat = _check('lat', lat, "float", ge=-90, le=90)
    lon = _check('lon', lon, "float", ge=-180, le=360)
    if pd.api.types.is_scalar(time):
        ts = pd.Series(pd.to_datetime(time))
        scalar = True
    else:
        ts = pd.to_datetime(time)
        scalar = False
    sr = pd.Series(np.nan, index=ts.index)
    st = pd.Series(np.nan, index=ts.index)
    ss = pd.Series(np.nan, index=ts.index)
    for i, t in enumerate(ts):
        if t.tzinfo is None:
            tz = 0
        else:
            tz = t.utcoffset().seconds / 3600.
        # Julian Day
        F2 = t.to_julian_date()
        # Julian Century
        G2 = (F2 - 2451545) / 36525
        # Geom Mean Long Sun (deg)
        I2 = (280.46646 + G2 * (36000.76983 + G2 * 0.0003032)) % 360.
        # Geom Mean Anom Sun (deg)
        J2 = 357.52911 + G2 * (35999.05029 - 0.0001537 * G2)
        # Eccent Earth Orbit
        K2 = 0.016708634 - G2 * (0.000042037 + 0.0000001267 * G2)
        # Sun Eq of Ctr
        L2 = (np.sin(np.deg2rad(J2)) * (1.914602 -
                                        G2 * (0.004817 + 0.000014 * G2)) +
              np.sin(np.deg2rad(2 * J2)) * (0.019993 - 0.000101 * G2) +
              np.sin(np.deg2rad(3 * J2)) * 0.000289
              )
        # Sun True Long (deg)
        M2 = I2 + L2
        # Sun True Anom (deg)
        # N2 = J2+L2
        # Sun Rad Vector (AUs)
        # O2 = (1.000001018*(1-K2*K2))/(1+K2*np.cos(np.deg2rad(N2)))
        # Sun App Long (deg)
        P2 = M2 - 0.00569 - 0.00478 * \
            np.sin(np.deg2rad(125.04 - 1934.136 * G2))
        # Mean Obliq Ecliptic (deg)
        Q2 = 23 + (26 + ((21.448 -
                          G2 * (46.815 +
                                G2 * (0.00059 -
                                      G2 * 0.001813)))) / 60) / 60
        # Obliq Corr (deg)
        R2 = Q2 + 0.00256 * np.cos(np.deg2rad(125.04 - 1934.136 * G2))
        # Sun Rt Ascen (deg)
        # S2 = np.rad2deg(np.arctan2(np.cos(np.deg2rad(P2)),
        #      np.cos(np.deg2rad(R2))*np.sin(np.deg2rad(P2))))
        # Sun Declin (deg)
        T2 = np.rad2deg(
            np.arcsin(np.sin(np.deg2rad(R2)) * np.sin(np.deg2rad(P2))))
        # var y
        U2 = np.tan(np.deg2rad(R2 / 2)) * np.tan(np.deg2rad(R2 / 2))
        # Eq of Time (minutes)
        V2 = 4 * np.rad2deg(U2 * np.sin(2 * np.deg2rad(I2)) -
                            2 * K2 * np.sin(np.deg2rad(J2)) +
                            4 * K2 * U2 * np.sin(
                                np.deg2rad(J2)) * np.cos(2 * np.deg2rad(I2)) -
                            0.5 * U2 * U2 * np.sin(4 * np.deg2rad(I2)) -
                            1.25 * K2 * K2 * np.sin(2 * np.deg2rad(J2)))
        # HA Sunrise (deg)
        W2 = np.rad2deg(np.arccos(
            np.cos(np.deg2rad(90.833)) /
            (np.cos(np.deg2rad(lat)) * np.cos(np.deg2rad(T2))) -
            np.tan(np.deg2rad(lat)) * np.tan(np.deg2rad(T2))))
        # Solar Noon (LST)
        X2 = (720 - 4 * lon - V2 + tz * 60) / 1440
        # Sunrise Time (LST)
        Y2 = X2 - W2 * 4 / 1440
        # Sunset Time (LST)
        Z2 = X2 + W2 * 4 / 1440
        # day -> hour of day
        sr[i] = Y2 * 24.
        st[i] = X2 * 24.
        ss[i] = Z2 * 24.
    if scalar:
        sr = sr[0]
        st = st[0]
        ss = ss[0]
    return sr, st, ss

# ----------------------------------------------------


def spa_sun_position(time, lat, lon, ele=0, pp=None, tk=None):
    '''
    Calculate apparent position of the sun at time, at postion lat,lon
    according [Esp2006]_.

    :param time: friction velocity in m/s (datetime).
    :param lat: position latitude in degrees (float).
    :param lon: position longtude in degrees (float).
    :param ele: (optional) position elevation in m (float).
        Defaults to 0 m, if missing
    :param pp: (optional) mean atmospheric pressure at position in hPa (float).
      None causes atmospheric aberration to be neglected. Defaults to None.
    :param pp: (optional) mean air temperature at position in K (float).
      None causes atmospheric aberration to be neglected. Defaults to None.

    :return: elevation angle, azimuth angle (counting eastward from north)
    :rtype: float
    '''
    #        azimuth_astro angle (counting westward from south)
    #
    time_i = pd.to_datetime(time)
    lat = _check('lat', lat, "float", ge=-90, le=90)
    lon = _check('lon', lon, "float", ge=-180, le=360)
    pp = _check('pp', pp, "float", ge=20, le=1100, none=True)
    tk = _check('tk', tk, "float", ge=150, le=360, none=True)
    if pd.api.types.is_scalar(time_i):
        zenith, azimuth, _ = _spa_position(time_i, lat, lon, ele, pp, tk)
    else:
        # expand scalar values to series, if provided as scalars
        lat_i = _expand_to_series_like(lat, time_i)
        lon_i = _expand_to_series_like(lon, time_i)
        ele_i = _expand_to_series_like(ele, time_i)
        pp_i = _expand_to_series_like(pp, time_i)
        tk_i = _expand_to_series_like(tk, time_i)
        # make empty output series
        zenith = _expand_to_series_like(np.nan, time_i)
        azimuth = _expand_to_series_like(np.nan, time_i)
        # calculate
        for i, arg in enumerate(
                zip(time_i, lat_i, lon_i, ele_i, pp_i, tk_i)):
            logging.debug(arg[0])
            zenith[i], azimuth[i], _ = _spa_position(*arg)

    res = (90 - zenith, azimuth)
    return res

# ----------------------------------------------------


def fast_sun_position(time, lat, lon):
    '''
    Calculate Sunrise, Sun Transit, and Sunset as hours of day
    according to [Mee1998]_,
    as described under ``https://gml.noaa.gov/grad/solcalc/calcdetails.html``.
    Caluculation is only valid for dates years 1901 and 2099
    and is less accurate than `sun_rise_transit_set`,
    but runs much (180 times) faster.

    :param time: friction velocity in m/s (datetime).
    :param lat: position latitude in degrees (float).
    :param lon: position longtude in degrees (float).
    :return: elevation angle, azimuth angle (counting eastward from north)
    :rtype: float
    '''
    time = _check("time", time, "datetime",
                  ge=pd.to_datetime("1901-01-01", utc=True),
                  lt=pd.to_datetime("2100-01-01", utc=True))
    lat = _check('lat', lat, "float", ge=-90, le=90)
    lon = _check('lon', lon, "float", ge=-180, le=360)
    if pd.api.types.is_scalar(time):
        ts = pd.Series(pd.to_datetime(time))
        scalar = True
    else:
        ts = pd.to_datetime(time)
        scalar = False
    ele = pd.Series(np.nan, index=ts.index)
    azi = pd.Series(np.nan, index=ts.index)
    for i, t in enumerate(ts):
        if t.tzinfo is None:
            tz = 0
        else:
            tz = t.utcoffset().seconds / 3600.
        # Time (past local midnight) in days
        E2 = t.hour / 24 + t.minute / (24 * 60) + t.second / (24 * 3600)
        # Julian Day
        F2 = t.to_julian_date()
        # Julian Century
        G2 = (F2 - 2451545) / 36525
        # Geom Mean Long Sun (deg)
        I2 = (280.46646 + G2 * (36000.76983 + G2 * 0.0003032)) % 360.
        # Geom Mean Anom Sun (deg)
        J2 = 357.52911 + G2 * (35999.05029 - 0.0001537 * G2)
        # Eccent Earth Orbit
        K2 = 0.016708634 - G2 * (0.000042037 + 0.0000001267 * G2)
        # Sun Eq of Ctr
        L2 = (np.sin(np.deg2rad(J2)) * (1.914602 -
                                        G2 * (0.004817 + 0.000014 * G2)) +
              np.sin(np.deg2rad(2 * J2)) * (0.019993 - 0.000101 * G2) +
              np.sin(np.deg2rad(3 * J2)) * 0.000289
              )
        # Sun True Long (deg)
        M2 = I2 + L2
        # Sun True Anom (deg)
        # N2 = J2+L2
        # Sun Rad Vector (AUs)
        # O2 = (1.000001018*(1-K2*K2))/(1+K2*np.cos(np.deg2rad(N2)))
        # Sun App Long (deg)
        P2 = M2 - 0.00569 - 0.00478 * \
            np.sin(np.deg2rad(125.04 - 1934.136 * G2))
        # Mean Obliq Ecliptic (deg)
        Q2 = 23 + (26 + ((21.448 -
                          G2 * (46.815 +
                                G2 * (0.00059 -
                                      G2 * 0.001813)))) / 60) / 60
        # Obliq Corr (deg)
        R2 = Q2 + 0.00256 * np.cos(np.deg2rad(125.04 - 1934.136 * G2))
        # Sun Rt Ascen (deg)
        # S2 = np.rad2deg(np.arctan2(np.cos(np.deg2rad(P2)),
        #      np.cos(np.deg2rad(R2))*np.sin(np.deg2rad(P2))))
        # Sun Declin (deg)
        T2 = np.rad2deg(
            np.arcsin(np.sin(np.deg2rad(R2)) * np.sin(np.deg2rad(P2))))
        # var y
        U2 = np.tan(np.deg2rad(R2 / 2)) * np.tan(np.deg2rad(R2 / 2))
        # Eq of Time (minutes)
        V2 = 4 * np.rad2deg(U2 * np.sin(2 * np.deg2rad(I2))
                            - 2 * K2 * np.sin(np.deg2rad(J2))
                            + 4 * K2 * U2 *
                            np.sin(np.deg2rad(J2)) * np.cos(2 * np.deg2rad(I2))
                            - 0.5 * U2 * U2 * np.sin(4 * np.deg2rad(I2))
                            - 1.25 * K2 * K2 * np.sin(2 * np.deg2rad(J2)))
        # HA Sunrise (deg)
        # W2 = np.rad2deg(np.arccos(np.cos(np.deg2rad(90.833)) /
        #                           (np.cos(np.deg2rad(lat)) *
        #                            np.cos(np.deg2rad(T2))) -
        #                           np.tan(np.deg2rad(lat)) *
        #                           np.tan(np.deg2rad(T2))))
        # Solar Noon (LST)
        # X2 = (720-4*lon-V2+tz*60)/1440
        # Sunrise Time (LST)
        # Y2 = X2-W2*4/1440
        # Sunset Time (LST)
        # Z2 = X2+W2*4/1440
        # Sunlight Duration (minutes)
        # AA2 = 8*W2
        # True Solar Time (min)
        AB2 = (E2 * 1440 + V2 + 4 * lon - 60 * tz) % 1440.
        # Hour Angle (deg)
        if AB2 / 4 < 0:
            AC2 = AB2 / 4 + 180
        else:
            AC2 = AB2 / 4 - 180
        # Solar Zenith Angle (deg)
        AD2 = np.rad2deg(np.arccos(np.sin(np.deg2rad(lat)) *
                                   np.sin(np.deg2rad(T2)) +
                                   np.cos(np.deg2rad(lat)) *
                                   np.cos(np.deg2rad(T2)) *
                                   np.cos(np.deg2rad(AC2))))
        # Solar Elevation Angle (deg)
        AE2 = 90. - AD2
        # Approx Atmospheric Refraction (deg)
        if AE2 > 85:
            AF2 = 0
        elif AE2 > 5:
            AF2 = (58.1 / np.tan(np.deg2rad(AE2))
                   - 0.07 / np.tan(np.deg2rad(AE2))**3
                   + 0.000086 / np.tan(np.deg2rad(AE2))**5
                   ) / 3600.
        elif AE2 > -0.575:
            AF2 = (1735.
                   - 518.2 * AE2
                   + 103.4 * AE2**2
                   - 12.79 * AE2**3
                   + 0.711 * AE2**4
                   ) / 3600.
        else:
            AF2 = (-20.772 / np.tan(np.deg2rad(AE2))) / 3600.
        # Solar Elevation corrected for atm refraction (deg)
        AG2 = AE2 + AF2
        # Solar Azimuth Angle (deg cw from N)
        if AC2 > 0:
            AH2 = (np.rad2deg(np.arccos(((np.sin(np.deg2rad(lat)) *
                                          np.cos(np.deg2rad(AD2))) -
                   np.sin(np.deg2rad(T2))) / (np.cos(np.deg2rad(lat)) *
                                              np.sin(np.deg2rad(AD2))
                                              ))) + 180) % 360.
        else:
            AH2 = (540 - np.rad2deg(np.arccos(((np.sin(np.deg2rad(lat)) *
                                                np.cos(np.deg2rad(AD2))) -
                   np.sin(np.deg2rad(T2))) / (np.cos(np.deg2rad(lat)) *
                                              np.sin(np.deg2rad(AD2))
                                              )))) % 360.
        # day -> hour of day
        ele[i] = AG2
        azi[i] = AH2
    if scalar:
        ele = ele[0]
        azi = azi[0]
    return ele, azi
