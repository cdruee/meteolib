#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
This module conatins standard equations, constants and conversions
adopted or recommended by the World meteorological Organization (WMO)
for general use in meteorology.

TBD:

* mol / density

'''

import numpy as np
import pandas as pd

from .__version__ import __title__, __description__, __version__
from .__version__ import __url__, __author__, __author_email__
from .__version__ import __license__, __copyright__

from ._utils import *
from . import _utils

from .constants import *
from . import constants

from .evapo import *
from . import evapo

from .humidity import *
from . import humidity

from .pressure import *
from . import pressure

from .radiation import *
from . import radiation

from .temperature import *
from . import temperature

from .thermodyn import *
from . import thermodyn

from .wind import *
from . import wind




