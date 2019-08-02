#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd

from .__version__ import __title__, __description__, __version__
from .__version__ import __url__, __author__, __author_email__
from .__version__ import __license__, __copyright__

#---------------------------------------------------------------------
'''Constants

'''
Kelvin = 273.15
kappa = 0.4
cp = 1004.
gn = 9.806
h_beam = 1.2
z0 = 0.0001


#---------------------------------------------------------------------

def magnus(t):
  '''Magnus Formula
  
    blah blah
    
  '''
  if pd.isnull(t):
    return(np.nan)
  # kelvin/Celsius autodetect  
  if t > 100.: 
    t=t-273.15
  e = 6.11 * np.exp((17.1*t)/(235+t)) # hPa
  return(e)
