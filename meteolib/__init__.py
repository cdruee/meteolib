#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
This module conatins standard equations, constants and conversions
adopted or recommended by the World meteorological Organization (WMO)
for general use in meteorology.

TBD:

* mol / density

'''

from ._version import __title__, __description__, __url__, __version__
from ._version import __author__, __author_email__
from ._version import __license__, __copyright__

from . import _utils
from . import constants
from . import evapo
from . import humidity
from . import pressure
from . import radiation
from . import standard
from . import temperature
from . import thermodyn
from . import wind
from . import charts

__all__ = ['_utils', 'constants',
           'evapo', 'humidity', 'pressure', 'radiation',
           'temperature', 'standard', 'thermodyn', 'wind', 'charts',
           '__title__', '__description__', '__url__', '__version__',
           '__author__', '__author_email__',
           '__license__', '__copyright__',
           ]
