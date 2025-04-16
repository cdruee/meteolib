#!/usr/bin/env python3
# -*-  coding: utf-8 -*-
"""
Funtions to retrieve radiosonde ascents from the
Integrated Global Radiosonde Archive (IGRA)
at `<www.ncei.noaa.gov>`_
"""

import difflib
import io
import logging
import re
import textwrap
import urllib.parse
import zipfile

import numpy as np
import pandas as pd
import requests

try:
    import requests_cache

    have_cache = True
except ImportError:
    requests_cache = None
    have_cache = False

# ------------------------------------------------------------------------

BASE_URL = "https://www.ncei.noaa.gov/" + \
           "data/integrated-global-radiosonde-archive/"
STATIONDOC = "doc/igra2-station-list.txt"
FORMATDOC = "doc/igra2-data-format.txt"
DATAFMT = "access/data-por/%11s-data.txt.zip"
VERSION = "v2.2"
STATIONLIST = None
intialized = False

FORMATDESC = {
    "v2.2": {
        "station-list": """
            ------------------------------
            Variable   Columns   Type
            ------------------------------
            ID            1-11   Character
            LATITUDE     13-20   Real
            LONGITUDE    22-30   Real
            ELEVATION    32-37   Real
            STATE        39-40   Character
            NAME         42-71   Character
            FSTYEAR      73-76   Integer
            LSTYEAR      78-81   Integer
            NOBS         83-88   Integer
            ------------------------------
            """,
        "header-record": """
            ---------------------------------
            Variable   Columns  Type
            ---------------------------------
            HEADREC       1-  1  Character
            ID            2- 12  Character
            YEAR         14- 17  Integer
            MONTH        19- 20  Integer
            DAY          22- 23  Integer
            HOUR         25- 26  Integer
            RELTIME      28- 31  Integer
            NUMLEV       33- 36  Integer
            P_SRC        38- 45  Character
            NP_SRC       47- 54  Character
            LAT          56- 62  Integer
            LON          64- 71  Integer
            ---------------------------------
            """,
        "data-record": """
            -------------------------------
            Variable        Columns Type
            -------------------------------
            LVLTYP1         1-  1   Integer
            LVLTYP2         2-  2   Integer
            ETIME           4-  8   Integer
            PRESS          10- 15   Integer
            PFLAG          16- 16   Character
            GPH            17- 21   Integer
            ZFLAG          22- 22   Character
            TEMP           23- 27   Integer
            TFLAG          28- 28   Character
            RH             29- 33   Integer
            DPDP           35- 39   Integer
            WDIR           41- 45   Integer
            WSPD           47- 51   Integer
            -------------------------------
            """,
        "units": {
            "LAT": {
                "factor": 0.0001, "unit": 'degrees', "nan": [-988888]},
            "LON": {
                "factor": 0.0001, "unit": 'degrees', "nan": [-9988888]},
            "PRESS": {
                "factor": 0.01, "unit": 'hPa', "nan": [-9999]},
            "GPH": {
                "factor": 1., "unit": 'gpm', "nan": [-8888, -9999]},
            "TEMP": {
                "factor": 0.1, "unit": 'C', "nan": [-8888, -9999]},
            "RH": {
                "factor": 0.1, "unit": '%', "nan": [-8888, -9999]},
            "DPDP": {
                "factor": 0.1, "unit": 'C', "nan": [-8888, -9999]},
            "WDIR": {
                "factor": 1., "unit": 'degrees', "nan": [-8888, -9999]},
            "WSPD": {
                "factor": 0.1, "unit": 'm/s', "nan": [-8888, -9999]},
        }
    }
}

# ------------------------------------------------------------------------

logger = logging.getLogger()

if have_cache:
    logger.debug("web request are cached")
    requests_cache.install_cache('IGRA_cache', expire_after=86400)
else:
    logger.debug("web request are NOT cached")


# ------------------------------------------------------------------------

def _parse_formatdesc(ver, item, what):
    if not intialized:
        _initialize()
    if ver not in FORMATDESC.keys():
        raise ValueError('unknown format version: %s' % ver)
    if item not in FORMATDESC[ver].keys():
        raise ValueError('unknwon format item: %s' % item)
    desc = textwrap.dedent(FORMATDESC[ver][item])
    res = {}
    for line in io.StringIO(desc).readlines():
        if "---" in line:
            continue
        if "Columns" in line:
            continue
        shrunk = re.sub(r"(\s+|-\s*)", " ", line.strip())
        if shrunk == '':
            continue
        variable, col_from, col_to, coltype = shrunk.split()
        if what == "cols":
            res[variable] = (int(col_from) - 1, int(col_to))
        elif what == "dtype":
            res[variable] = coltype
        else:
            raise ValueError('dont know what to parse: %s' % what)
    return res


def _parse_format_cols(ver, item):
    return _parse_formatdesc(ver, item, "cols")


def _parse_format_dtype(ver, item):
    return _parse_formatdesc(ver, item, "dtype")


def _get_patiently(url, timeout=5, retry=5):
    for i in range(retry):
        res = None
        logging.debug("... try to get #%i" % i)
        try:
            res = requests.get(url, timeout=timeout)
        except (requests.exceptions.ReadTimeout,
                requests.exceptions.ConnectionError):
            pass
        if res is not None and res.status_code == 200:
            break
    else:
        raise IOError('no answer from %s' %
                      urllib.parse.urlparse(url).netloc)
    return res


def _initialize():
    logger.info('trying to reach %s' % urllib.parse.urlparse(BASE_URL).netloc)
    url = urllib.parse.urljoin(BASE_URL, FORMATDOC)
    r = _get_patiently(url)
    txt = r.text
    logger.debug('downloadad format description OK')
    line = None
    for ll in io.StringIO(txt).readlines():
        if ll.strip() != "":
            line = ll
            break
    if line is None:
        raise ValueError('cannot determine data version number')
    formatstr = re.search("[vV][0-9.]+[^ ]*", line).group()

    global VERSION
    VERSION = formatstr.lower()


def _distance(lat1, lon1, lat2, lon2):
    """
    Calculate the great circle distance between two points
    (specified in decimal degrees) on a spheric earth.
    Reference:
    https://stackoverflow.com/a/29546836/7657658

    :param lat1: Position 1 latitude in degrees
    :type: float
    :param lon1: Position 1 longitude in degrees
    :type: float
    :param lat2: Position 2 latitude in degrees
    :type: float
    :param lon2: Position 2 longitude in degrees
    :type: float
    :returns: Great circle distance in km
    :rtype: float
    """
    rlat1 = np.radians(lat1)  # deg -> rad
    rlon1 = np.radians(lon1)  # deg -> rad
    rlat2 = np.radians(lat2)  # deg -> rad
    rlon2 = np.radians(lon2)  # deg -> rad

    dlon = rlon2 - rlon1  # rad
    dlat = rlat2 - rlat1  # rad
    a = (np.sin(dlat / 2.0) ** 2 +
         np.cos(rlat1) * np.cos(rlat2) * np.sin(dlon / 2.0) ** 2)
    c = 2 * np.arcsin(np.sqrt(a))
    km = 6371 * c  # km

    return km


def get_stationlist():
    """
    Get the current list of all stations.
    :return: List of all stations.
    :rtype: `pandas.DataFrame`

    Columns (as of IGRA Version 2.2) are:
      - ID (alphanumeric station ID),
      - LATITUDE, LONGITUDE (position ind degrees),
      - ELEVATION (in m),
      - STATE (two-letter US-state abbreviation, empty outside the US),
      - NAME (truncated to 30 characrters),
      - FSTYEAR, LSTYEAR (first and last year of operation),
      - NOBS (number of soundings in archive)
    """
    if not intialized:
        _initialize()
    logger.info('trying to reach %s' % urllib.parse.urlparse(BASE_URL).netloc)
    url = urllib.parse.urljoin(BASE_URL, STATIONDOC)
    r = _get_patiently(url, timeout=5)
    txt = r.text
    logger.debug('downloadad station list OK')

    global STATIONLIST
    cols = _parse_format_cols(VERSION, "station-list")
    res = pd.read_fwf(io.StringIO(txt),
                      colspecs=list(cols.values()),
                      names=cols.keys())
    return res


def get_station_id(info):
    """
    Determine the station ID from a station name,
    position or WMO-number.

    :param info: Info on the desired station (see below).
    :type info: str, int or tuple
    :return: station ID
    :rtype: str

    `info` can be given as name, number, or position.

      In case `info` is Name, it is compared (not case-sensitive) against all
      names on the list. If no name on the list matches exactly, the most
      similar name is chosen.

      In case `info` is a number (numeric or a sting containing digits),
      the numer is left padded with zeros to ten digits and then matched
      against the last ten characters of the station IDs (that contain the
      WMO-number if a station has one).

      In case `info` is a position, isst mut be given as length-2
      tuple of latitude and longitude in degrees. Then the closest
      station to this position is returned (the distance is approximated
      by the distance on a spherical earth).

    """
    if isinstance(info, tuple) and len(info) == 2:
        info_type = "latlon"
        lat, lon = info
    elif ((isinstance(info, str) and info.isnumeric()) or
          isinstance(info, int)):
        info_type = "number"
        number = "%10i" % int(info)
    elif re.match("[A-Z]{3}[0-9]{10}", info):
        info_type = "id"
        stnid = info
    else:
        info_type = "name"
        name = info

    global STATIONLIST
    if STATIONLIST is None:
        STATIONLIST = get_stationlist()

    match = None
    if info_type in ["number", "id"]:
        for k, v in STATIONLIST.to_dict(orient="index").items():
            if number in v['ID'] or stnid == v['ID']:
                match = k
                break
    elif info_type == "latlon":
        dist = STATIONLIST.apply(lambda x: _distance(
            lat, lon, STATIONLIST.LATITUDE, STATIONLIST.LONGITUDE))
        match = dist.idxmin()
    elif info_type == "name":
        for k, v in STATIONLIST.to_dict(orient="index").items():
            if name.upper() == v['NAME'].upper():
                match = k
                break
        else:
            cm = difflib.get_close_matches(
                name, STATIONLIST["NAME"], n=1, cutoff=0.5
            )
            logging.warning("not found, diy ou mean %s" % str(cm))
            match = STATIONLIST.index[STATIONLIST['NAME'] == cm]
    return STATIONLIST.loc[match, 'ID']


def get_station_file(station):
    """
    Download the entire dataset for station `id`

    :param station: Station ID
    :type station: str
    :return: the uncompressed station data file
    :rtype: str
    """
    if not intialized:
        _initialize()
    logger.info('trying to reach %s' % urllib.parse.urlparse(BASE_URL).netloc)
    filepath = DATAFMT % station
    url = urllib.parse.urljoin(BASE_URL, filepath)
    r = _get_patiently(url, timeout=5)
    zipped = r.content
    z = zipfile.ZipFile(io.BytesIO(zipped))
    txt = z.read(z.filelist[0].filename).decode()
    logger.debug('downloadad station data OK')
    return txt


def _parse_line(line, item):
    cols = _parse_format_cols(VERSION, item)
    dtype = _parse_format_dtype(VERSION, item)
    res = {}
    for k, v in cols.items():
        i, j = v
        cell = line[i:j]
        if dtype[k].lower() == "integer":
            res[k] = int(cell)
        else:
            res[k] = cell
    return res


def _parse_header_line(line):
    return _parse_line(line, "header-record")


def _parse_data_record(line):
    return _parse_line(line, "data-record")


def _file_headers(txt):
    logging.debug('scanning station data')
    res = []
    with io.StringIO(txt) as f:
        for i, l in enumerate(f.readlines()):
            if l.startswith('#'):
                h = _parse_header_line(l)
                h['line'] = i
                h['time'] = "%04i-%02i-%02iT%02i:00:00" % (
                    h["YEAR"], h["MONTH"], h["DAY"], h["HOUR"])
                res.append(h)
    logging.debug('... soundings found: %i' % len(res))
    return res


def _file_sounding(txt, time, headers=None):
    if headers is None:
        headers = _file_headers(txt)
    logging.debug('getting sounding: %s' % time)
    idx = [x['time'] for x in headers].index(time)
    lines = txt.splitlines()
    itop = headers[idx]['line'] + 1
    if idx + 1 > len(headers):
        iend = len(lines)
    else:
        iend = headers[idx + 1]['line']
    cols = _parse_format_cols(VERSION, "data-record")
    df = pd.read_fwf(io.StringIO('\n'.join(lines[itop: iend])),
                     colspecs=list(cols.values()),
                     names=cols.keys())
    return df


def get_sounding(station, time):
    """
    Return one sounding of station `station` at time `time`.

    :param station: station ID
    :type station: str
    :param time: time of sounding
    :type time: `str`, `datetime.datetime`, or `pandas.Timestamp`
    :return: The sounding data.
    :rtype: `pandas.DataFrame`

    Columns (as of IGRA Version 2.2) are:
      - LVLTYP1, LVLTYP2 (major and minor level indicator)
      - ETIME (time elapsed since launch),
      - PRESS (pressur in Pa),
      - PLFAG (pressure processing flag),
      - GPH (geopotential height ASL in m),
      - ZFLAG (geopotential height processing flag),
      - TEMP (temperature in 0.1 C),
      - TFLAG (temperature processing flag),
      - DPDP (dewpoint depression in 0.1 C),
      - WDIR (wind direction in degrees),
      - WSPD (wind speed in 0.1 m/s).

    """
    if not intialized:
        _initialize()
    time64 = pd.to_datetime(time)
    ts = time64.strftime("%Y-%m-%dT%H:%M:%S")
    station = get_station_id(station)
    txt = get_station_file(station)
    data = _file_sounding(txt, ts)
    for c in data.columns:
        if c in FORMATDESC[VERSION]['units']:
            unt = FORMATDESC[VERSION]['units'][c]
            data[c] = data[c].mask(data[c].isin(unt['nan']))
            data[c] = data[c] * unt['factor']
    return data
