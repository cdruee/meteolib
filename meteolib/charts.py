#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
provides meteorogical charts
"""
import logging
import os

import meteolib.humidity
import meteolib.standard
import numpy as np
import pandas as pd

try:
    import matplotlib
    if os.name == 'posix' and "DISPLAY" not in os.environ:
        matplotlib.use('Agg')
        have_display = False
    else:
        have_display = True
    from matplotlib import pyplot as plt
    have_matplotlib = True
except ImportError:
    matplotlib = None
    plt = None
    have_matplotlib = False
    have_display = False

from .constants import R, Lv, cp
from .humidity import Humidity
from .standard import altitude, p_iso, T_iso
from .temperature import KtoC, CtoK, inv_Tpot
from .temperature import POTENTIAL_REFERENCE_PRESSURE as PREF

logger = logging.getLogger()


# ----------------------------------------------------------------------
# adiabatic process lines
# ----------------------------------------------------------------------

def dryad(t0, pmax, pmin):
    """
    Calulates the temperature of an air parcel lifted
    adiabaticially from pmax to pmin,
    i.e. a (dry) adiabate.
    The parcel is assumed to attain the (potential)
    temperature t0 at reference level (generally 1000 hPa).

    :param t0: label temperature in C
    :type t0: float
    :param pmax: upper pressure limit in hPa (lifting start level)
    :type pmax: float
    :param pmin: lower pressure limit in hPa (lifting end level)
    :type pmin: float
    :return: pressure (in hPa) and temperature (in C) of air parcel
    :rtype: list[float], list[float]
    """
    p = np.arange(pmax, pmin, - 1)
    t = [inv_Tpot(t0, x, hPa=True, Kelvin=False) for x in p]
    return p, t


def satad(th, pmax, pmin, label=None):
    r"""
    Calulates the temperature of an air parcel lifted
    (pseudo)-adiabaticially from pmax to pmin,
    assuming all humidity condensating is immedately removed,
    i.e. a moist adiabate.
    The parcel is assumed to either attain the equivalent potential
    temperature th at great height, if label is `the` or to
    attain wet-bulb potential temperature th if label is `thw`.

    To calbulate wet-bulb potential temperature (:math:`\theta_w`) from
    equivalent potential temperature (:math:`\theta_e`) equation (3.8)
    from [DavJo2008]_.
    The vertical integration follows the process described by
    equation (8) in [BaSt2913]_.

    :param th: label temperature in C
    :type th: float
    :param pmax: upper pressure limit in hPa (lifting start level)
    :type pmax: float
    :param pmin: lower pressure limit in hPa (lifting end level)
    :type pmin: float
    :param label: type of label of air parcel. Possible are
        `the` for equivalent potential temperature (:math:`\theta_e`) or
        `thw` for wet-bulb potential temperature (:math:`\theta_w`).
    :return: pressure (in hPa) and temperature (in C) of air parcel
    :rtype: list[float], list[float]
    """
    the = None
    thw = None

    if label is None:
        logging.warning('No label for satad. Assuming "the".')
        label = "the"
    if label == "thw":
        thw = th
    elif label == "the":
        the = th

    if thw is None:
        c = 273.15
        thek = CtoK(the)
        x = thek/c
        if thek >= 173.15:  # K
            a0 = 7.101574
            a1 = -20.68208
            a2 = 16.11182
            a3 = 2.574631
            a4 = -5.205688
            b1 = -3.552497
            b2 = 3.781782
            b3 = -0.6899655
            b4 = -0.5929340
            thw = (thek - c -
                   np.exp((a0 + a1 * x + a2 * x * x +
                           a3 * x ** 3 + a4 * x ** 4) /
                          (1 + b1 * x + b2 * x * x +
                           b3 * x ** 3 + b4 * x ** 4))
                   )
        else:
            thw = thek - c  # C
        thw = CtoK(thw)  # K
    else:
        thw = CtoK(thw)  # K

    p0 = PREF / 100.  # hPa
    dp = -1  # hPa
    pp = np.arange(p0, pmin, dp)  # hPa

    pout = []
    tout = []
    tk = thw
    for p in pp:  # hPa
        epsilon = 0.622
        b = 1.
        m = Humidity(t=tk, Kelvin=True, p=p, hPa=True,
                     rh=1., percent=False).m(gkg=False)  # 1
        dtdp = (b / (p * 100.)) * \
               (R * tk + Lv * m) / (
                       cp + (Lv * Lv * m * epsilon * b) / (R * tk * tk)
               )  # K/Pa
        tk = tk + dtdp * (dp * 100.)  # K (dp: Pa -> hPa)
        tc = KtoC(tk)  # C
        if pmax >= p >= pmin:  # hPa
            pout.append(p)  # hPa
            tout.append(tc)  # C
    return pout, tout


def tdm(m, pmax, pmin):
    """
    Calulates the dew point of an air parcel lifted from
    reference level (generally 1000 hPa) maintaining
    a constant mixing ratio `m`

    :param m: mixing ratio (unitless, i.e. 1. = saturation)
    :type m: float
    :param pmax: upper pressure limit in hPa (lifting start level)
    :type pmax: float
    :param pmin: lower pressure limit in hPa (lifting end level)
    :type pmin: float
    :return: pressure (in hPa) and temperature (in C) of air parcel
    :rtype: list[float], list[float]
    """

    p = np.arange(pmax, pmin, -1)
    pout = []
    tout = []
    for pp in p:
        tt = T_iso(altitude(pp, hPa=True), Kelvin=False)
        td = meteolib.humidity.Humidity(m=m, t=tt, gkg=True, p=pp,
                                        Kelvin=False, hPa=True).td()
        if pmax >= pp >= pmin:
            pout.append(pp)
            tout.append(td)
    return pout, tout


# ----------------------------------------------------------------------
# Stüve Diagram
# ----------------------------------------------------------------------


def stueve(p=None, t=None, td=None, style=None, title=None,
           fname=None, fmt=None, dpi=300):
    """
    Produces a Stüve plot in the given style.
    Optionally, a vertical sounding may be supplied by giving
    pressure `p`, temperature `t` and/or dew point `td` as
    lists of equal length.

    :param p: optional, sounding pressure levels (in hPa)
    :type p: list[float]
    :param t: optional, sounding temperature values (in C)
    :type t: list[float]
    :param td: optional, sounding dew point temperatures (in C)
    :type td: list[float]
    :param style: optional, plot style. Currently available are
      `wyoming`, `dwd_a0` and `dwd_a4`. Defaults to `dwd_a0`.
    :type style: str
    :param title: optional, plot title.
    :type style: str
    :param fname: optional, filename (optionally including path)
      to save the plot in.
      If not present, the plot is shown on the screen.
    :type fname: str or path-like object
    :param fmt: optional,
        the file format, e.g. 'png', 'pdf', 'svg',
        If missing, the file format is nfered from the filename extension
        of `fname`.
    :param dpi: optional, plot resolution in dots per inch.
        Defaults to 300.
    """
    if not have_matplotlib:
        raise EnvironmentError('matplotlib not available')
    if fname is None and not have_display:
        raise EnvironmentError('no display, cannot show chart')

    def p2z(p: float) -> float:
        r"""
        convert pressure to Stüve vertical coordinate :math:`p^\kappa`
        :param p: pressure in hPa
        :type p: float
        :return: vertical coordinate (unitless)
        :rtype:  float
        """
        return p ** 0.2859

    def h2z(h):
        r"""
        convert altitude to Stüve vertical coordinate :math:`p^\kappa`
        :param h: pressure in hPa
        :type h: float
        :return: vertical coordinate (unitless)
        :rtype:  float
        """
        return p2z(p_iso(h, hPa=True))

    #
    # check arguments
    #

    if style is None:
        style = 'dwd_a0'
    if style == 'dwd_a0':
        sty = {
            "figsize": (32, 40),
            "trange": (-80, 50),
            "prange": (1050, 1),
            "meters": [
                {"at": np.arange(0., 100000., 1000.),
                 "unit": "km",
                 }
            ],
            "isotherms": [
                {"at": list(np.arange(-80, 50, 10)),
                 "color": "green",
                 "width": "1.5",
                 "style": "-"
                 },
                {"at": list(np.arange(-80, 50, 1)),
                 "color": "green",
                 "width": "0.5",
                 "style": "-"
                 }

            ],
            "isobars": [
                {"at": (list(np.arange(1050, 250, -50)) +
                        list(np.arange(250, 100, -25)) +
                        list(np.arange(100, 20, -10)) +
                        list(np.arange(20, 10, -5)) +
                        list(np.arange(10, 1.1, -2.5))
                        ),
                 "color": "green",
                 "width": "0.5",
                 "style": "-"
                 },
                {"at": (list(np.arange(1050, 250, -10)) +
                        list(np.arange(250, 100, -5)) +
                        list(np.arange(100, 20, -2)) +
                        list(np.arange(20, 10, -1)) +
                        list(np.arange(10, 1.1, -0.5))
                        ),
                 "color": "green",
                 "width": "1.5",
                 "style": "-"
                 },
            ],
            "dryadiabates": [
                {"at": list(np.arange(-70, 120, 10)),
                 "color": "green",
                 "width": "1",
                 "style": "-"
                 }
            ],
            "moistdiabates": [
                {"at": list(np.arange(-70, 120, 10)),
                 "label": "the",
                 "color": "red",
                 "width": "1",
                 "style": "-"
                 }
            ],
            "mixingratios": [
                {"at": (list(np.arange(0.2, 1.1, 0.1)) +
                        [1.2, 1.6] +
                        list(np.arange(2.0, 5.0, 0.5)) +
                        list(np.arange(5., 20., 1.)) +
                        list(np.arange(20., 30., 4.)) +
                        list(np.arange(30., 45., 5.))),  # g/kg
                 "color": "red",
                 "width": "1",
                 "style": "-"
                 }
            ],
            "temperature": {
                "color": "black",
                "width": "2",
                "style": "-"
            },
            "dewpoint": {
                "color": "blue",
                "width": "2",
                "style": "-"
            },
        }
    elif style == 'dwd_a4':
        sty = {
            "figsize": (8, 11),
            "trange": (-60, 40),
            "prange": (1050, 200),
            "meters": [
                {"at": np.arange(0., 12000., 1000.),
                 "unit": "km",
                 }
            ],
            "isotherms": [
                {"at": list(np.arange(-50, 40, 10)),
                 "color": "green",
                 "width": "1.",
                 "style": "-"
                 },
                {"at": list(np.arange(-55, 45, 1)),
                 "color": "green",
                 "width": "0.33",
                 "style": "-"
                 }

            ],
            "isobars": [
                {"at": (list(np.arange(1050, 250, -50)) +
                        list(np.arange(250, 100, -25)) +
                        list(np.arange(100, 20, -10)) +
                        list(np.arange(20, 10, -5)) +
                        list(np.arange(10, 1.1, -2.5))
                        ),
                 "color": "green",
                 "width": "1.",
                 "style": "-"
                 },
                {"at": (list(np.arange(1050, 250, -10)) +
                        list(np.arange(250, 100, -5)) +
                        list(np.arange(100, 20, -2)) +
                        list(np.arange(20, 10, -1)) +
                        list(np.arange(10, 1.1, -0.5))
                        ),
                 "color": "green",
                 "width": "0.33",
                 "style": "-"
                 },
            ],
            "dryadiabates": [
                {"at": list(np.arange(-50, 120, 10)),
                 "color": "green",
                 "width": "0.5",
                 "style": "-"
                 }
            ],
            "moistdiabates": [
                {"at": list(np.arange(-70, 120, 10)),
                 "label": "the",
                 "color": "red",
                 "width": "0.5",
                 "style": "-"
                 }
            ],
            "mixingratios": [
                {"at": (list(np.arange(0.2, 1.1, 0.1)) +
                        [1.2, 1.6] +
                        list(np.arange(2.0, 5.0, 0.5)) +
                        list(np.arange(5., 20., 1.)) +
                        list(np.arange(20., 30., 4.)) +
                        list(np.arange(30., 45., 5.))),  # g/kg
                 "color": "red",
                 "width": "0.5",
                 "style": "--"
                 }
            ],
            "temperature": {
                "color": "black",
                "width": "1.5",
                "style": "-"
            },
            "dewpoint": {
                "color": "blue",
                "width": "1.5",
                "style": "-"
            },
        }
    elif style == 'wyoming':
        sty = {
            "figsize": (6, 6),
            "trange": (-80, 45),
            "prange": (1050, 100),
            "meters": [
                {"at": [altitude(p) for p in np.arange(1000., 100., -100.)],
                 "unit": "km",
                 }
            ],
            "isotherms": [
                {"at": list(np.arange(-80, 41, 10)),
                 "color": "blue",
                 "width": "0.5",
                 "style": "-"
                 },
                {"at": list(np.arange(-75, 46, 10)),
                 "color": "blue",
                 "width": "0.5",
                 "style": "-"
                 },
            ],
            "isobars": [
                {"at": list(np.arange(1000, 99, -100)),
                 "color": "blue",
                 "width": "0.5",
                 "style": "-"
                 },
                {"at": list(np.arange(1050, 149, -100)),
                 "color": "blue",
                 "width": "0.5",
                 "style": "-"
                 },
            ],
            "dryadiabates": [
                {"at": list(np.arange(-60, 181, 20)),
                 "color": "green",
                 "width": "0.5",
                 "style": "-"
                 }
            ],
            "moistdiabates": [
                {"at": list(np.arange(-50, 131, 20)),
                 "label": "the",
                 "color": "blue",
                 "width": "0.5",
                 "style": "-"
                 }
            ],
            "mixingratios": [
                {"at": [0.1, 0.4, 1, 2, 4, 7, 10, 16, 24, 32, 40],  # g/kg
                 "color": "purple",
                 "width": "0.5",
                 "style": "-"
                 }
            ],
            "temperature": {
                "color": "black",
                "width": "1.5",
                "style": "-"
            },
            "dewpoint": {
                "color": "black",
                "width": "1.5",
                "style": "-"
            },
        }
    else:
        raise ValueError(f'Stueve style unknown: {style}')

    # plot boundaries
    tmin, tmax = sty["trange"]
    pmax, pmin = sty["prange"]
    zmin = p2z(pmax)
    zmax = p2z(pmin)

    # start plot
    fig, ax1 = plt.subplots(figsize=sty["figsize"],
                            dpi=dpi)
    plt.subplots_adjust(left=0.15)
    ax1.set_xlim(tmin, tmax)
    ax1.set_ylim(zmin, zmax)
    ax1.set_xticks([x for x in sty["isotherms"][0]["at"]
                    if tmin <= x <= tmax])
    ax1.set_yticks([p2z(x) for x in sty["isobars"][0]["at"]
                    if pmin <= x <= pmax],
                   ["%4.0f" % x for x in sty["isobars"][0]["at"]
                    if pmin <= x <= pmax])
    meter = [x for x in sty["meters"][0]["at"]
             if pmin <= p_iso(x, hPa=True) <= pmax]
    munit = sty["meters"][0]["unit"]
    if munit == "m":
        mstr = "m"
        mfac = 1.
    elif munit == "km":
        mstr = "km"
        mfac = 1. / 1000.
    elif munit == "kft":
        mstr = "kft"
        mfac = 3.28084 / 1000.
    else:
        mstr = ""
        mfac = 1.
    ax2 = ax1.twinx()
    ax2.set_ylim(zmin, zmax)
    ax2.set_yticks([p2z(p_iso(x, hPa=True)) for x in meter],
                   ["%3.0f%s" % (x * mfac, mstr) for x in meter])

    # add isolines
    for isotherms in sty["isotherms"]:
        for x in isotherms["at"]:
            ax1.axvline(x=x,
                        color=isotherms["color"],
                        linestyle=isotherms["style"],
                        linewidth=isotherms["width"])
    for isobars in sty["isobars"]:
        for x in isobars["at"]:
            ax1.axhline(y=p2z(x),
                        color=isobars["color"],
                        linestyle=isobars["style"],
                        linewidth=isobars["width"])
    for dryads in sty["dryadiabates"]:
        for x in dryads["at"]:
            pp, tt = dryad(x, pmax, pmin)
            ax1.plot(tt, [p2z(x) for x in pp],
                     color=dryads["color"],
                     linestyle=dryads["style"],
                     linewidth=dryads["width"])
    for satads in sty["moistdiabates"]:
        for x in satads["at"]:
            pp, tt = satad(x, pmax, pmin,
                           label=satads["label"])
            ax1.plot(tt, [p2z(x) for x in pp],
                     color=satads["color"],
                     linestyle=satads["style"],
                     linewidth=satads["width"])
    for ms in sty["mixingratios"]:
        for x in ms["at"]:
            pp, tt = tdm(x, pmax,  pmin)
            ax1.plot(tt, [p2z(x) for x in pp],
                     color=ms["color"],
                     linestyle=ms["style"],
                     linewidth=ms["width"])

    # add sounding if supplied
    if p is not None and t is not None:
        ax1.plot(t, [p2z(x) for x in p],
                 color=sty["temperature"]["color"],
                 linestyle=sty["temperature"]["style"],
                 linewidth=sty["temperature"]["width"])
    if p is not None and td is not None:
        ax1.plot(td, [p2z(x) for x in p],
                 color=sty["dewpoint"]["color"],
                 linestyle=sty["dewpoint"]["style"],
                 linewidth=sty["dewpoint"]["width"])

    if title is not None:
        plt.title(label=title)
    if fname:
        plt.savefig(fname, format=fmt, dpi=dpi)
    else:
        plt.show()


# ----------------------------------------------------------------------
# plot archived radiosonde
# ----------------------------------------------------------------------


def plot_radisonde(data, diagram=None, **kwargs):
    """
    Plot a radiosonde sounding contained in `data`

    :param data: Sounding data. Either a pandas DataFrame or a dict
      of array-like objects. If data is a dict, all array-like objects
      must be one-dmensional and of the same lenght.
      `data` must contain columns / elements for
      pressure, temperature and humidity.
      Humidity can be specified as dew point, relative humidity or
      dewpoint depression.
    :type data: `pandas.DataFrame` or `dict`
    :param diagram: type of diagram to draw.
      Currently only 'stueve' is available.
    :type diagram: str
    :param kwargs: arguments passed to the diagram function


      | Accepted names for pressure are:
        'p', 'pp', 'P', 'PRES', 'PRESS', 'pressure'
      | Accepted names for temperature are:
        't', 'tt', 'T', 'TEMP', 'temperature'
      | Accepted names for dewpoint are:
        'td', 'TD', 'DWPT', 'dewpoint'
      | Accepted names for dewpoint depression are:
        'dtd', 'DTD', 'DPDP', 'spread'
      | Accepted names for relative humidity are:
        'rh', 'RH', 'RELH', 'relative humidity'

    """
    def _getcol(dataobject: dict | pd.DataFrame, names):
        """

        :type dataobject: object
        """
        for col in names:
            if col in dataobject.columns:
                res = dataobject[col]
                break
        else:
            res = None
        return res

    p_names = ['p', 'pp', 'P', 'PRES', 'PRESS', 'pressure']
    t_names = ['t', 'tt', 'T', 'TEMP', 'temperature']
    td_names = ['td', 'TD', 'DWPT', 'dewpoint']
    sp_names = ['dtd', 'DTD', 'DPDP', 'spread']
    rh_names = ['rh', 'RH', 'RELH', 'relative humidity']

    style = kwargs.pop('style', {})
    if diagram is None:
        diagram = 'stueve'
        style = 'dwd_a4'

    if isinstance(data, (pd.DataFrame, dict)):
        p = _getcol(data, p_names)
        t = _getcol(data, t_names)
        td = _getcol(data, td_names)
        sp = _getcol(data, sp_names)
        rh = _getcol(data, rh_names)
    else:
        raise ValueError('data is fo unsopported type')

    if td is None:
        if sp is not None:
            td = [xt - xsp for xt, xsp in zip(t, sp)]
        elif rh is not None:
            for xt, xp, xrh in zip(t, p, rh):
                td = Humidity(t=xt, p=xp, rh=xrh).td(Kelvin=False)
        else:
            logger.warning('no humidty info')
            td = [np.nan for x in t]

    if np.nanmax(p) > 1500.:
        t = [x/100. for x in p]    # -> hPa
    if np.nanmax(t) > 100.:
        t = [KtoC(x) for x in t]   # -> C
    if np.nanmax(td) > 100.:
        td = [KtoC(x) for x in td]  # -> C

    if diagram == 'stueve':
        stueve(p, t, td, style=style, **kwargs)
