#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
utilities for internal use
"""
import numpy as np

import pandas as pd

KIND_DTYPES = {
    "float": "float64",
    "int": "int64",
    "str": "string",
    "datetime": "datetime64[ns, UTC]",
    "bool": "boolean",
}


def _check_scalar(name, par, kind, lt=None, le=None, ne=None,
                  ge=None, gt=None, nan=True, none=False):

    # None: is acceptable
    if par is None:
        if none is False:
            raise ValueError('parameter {} is None'.format(name))
        else:
            return None

    # check if value is of correct type or may be converted without error
    kind_ok = None
    if kind == 'float':
        try:
            val = float(par)
            kind_ok = True
        except ValueError:
            kind_ok = False
    elif kind == 'int':
        try:
            val = int(par)
            kind_ok = True
        except ValueError:
            kind_ok = False
    elif kind == 'bool':
        if isinstance(par, bool) is True:
            val = par
            kind_ok = True
        else:
            try:
                val = bool(int(par))
                kind_ok = True
            except ValueError:
                kind_ok = False
    elif kind == 'str':
        if isinstance(par, str) is True:
            val = par
            kind_ok = True
        else:
            kind_ok = False
    elif kind == 'datetime':
        try:
            val = pd.to_datetime(par)
            kind_ok = True
        except (TypeError, pd.errors.ParserError):
            kind_ok = False
    else:
        raise ValueError('unknown kind: {}'.format(kind))
    if kind_ok is False:
        raise TypeError('parameter {} has wrong type: {}'.format(name, kind))

    if kind in ['int', 'float', 'datetime']:
        if lt is not None and val >= lt:
            raise ValueError('parameter {} is not < {}'.format(name, lt))
        if le is not None and val > le:
            raise ValueError('parameter {} is not <= {}'.format(name, le))
        if ne is not None and val == ne:
            raise ValueError('parameter {} is = {}'.format(name, ne))
        if ge is not None and val < ge:
            raise ValueError('parameter {} is not >= {}'.format(name, ge))
        if gt is not None and val <= gt:
            raise ValueError('parameter {} is not > {}'.format(name, gt))

        if nan is False and np.isfinite(val) is False:
            raise ValueError('parameter {} is nan'.format(name))
    else:
        if any(x is not None for x in [lt, le, ne, ge, gt]):
            raise ValueError('parameters lt, le, ne, ge, gt are only' +
                             ' allowed with types int and float')

    return val


def _check(name, par, kind, lt=None, le=None, ne=None,
           ge=None, gt=None, nan=True, none=False):
    if par is None or pd.api.types.is_scalar(par):
        return _check_scalar(name, par, kind, lt, le,
                             ne, ge, gt, nan, none)
    else:
        if hasattr(kind, 'dtype'):
            res_type = par.dtype
        else:
            res_type = KIND_DTYPES[kind]
        res = pd.Series(index=range(len(par)), dtype=res_type)
        for i, p in enumerate(par):
            res.iloc[i] = _check_scalar(name, p, kind, lt, le,
                                        ne, ge, gt, nan, none)
        return res


def _only(pars, sel, ign=[]):
    if isinstance(sel, str):
        sel = [sel]
    elif isinstance(sel, list):
        pass
    else:
        raise TypeError('parameter "sel" must be str or list')
    if isinstance(ign, str):
        ign = [ign]
    elif isinstance(ign, list):
        pass
    else:
        raise TypeError('parameter "ign" must be str or list')
    if not isinstance(pars, dict) or not isinstance(sel, list):
        raise TypeError('parameter "pars" must be dict')
    for x in sel:
        if pars[x] is None:
            return False
    for k, v in pars.items():
        if k not in sel and k not in ign and v is not None:
            return False
    return True


def _expand_to_series_like(val, like):
    if pd.api.types.is_scalar(val):
        res = pd.Series(val,  index=range(len(like)))
    else:
        if len(val) != len(like):
            raise ValueError('lengts of series are not equal')
        res = val
    return res
