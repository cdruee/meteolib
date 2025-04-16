#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 18:39:22 2019

@author: druee
"""

import unittest

import meteolib as m

import pandas as pd

iso_tab = pd.DataFrame.from_records(
    [[-1000, -1000, 294.651, 113931., 1.34702],
     [1000,  1000, 281.651, 89876.3, 1.11166],
     [2000,  1999, 275.150, 79501.4, 1.00655],
     [5000,  4996, 255.676, 54048.3, 0.736429],
     [10000,  9984, 223.252, 26499.9, 0.413510],
     [15000, 14965, 216.650, 12111.8, 0.194755],
     [25000, 24902, 221.552, 2549.21, 4.00837E-2],
     [40000, 39750, 250.350, 287.143, 3.99566E-3],
     [55000, 54528, 260.771, 42.5249, 5.36680E-4],
     [65000, 64342, 233.292, 10.9297, 1.63209E-4],
     [75000, 74125, 208.399, 2.38814, 3.99210E-5],
     [80000, 79006, 198.639, 1.05247, 1.84580E-5]],
    columns=["h", "H", "T", "p", "rho"]
)


class Test_pressure(unittest.TestCase):
    def test_iso_zero(self):
        res = m.standard.p_iso(0)
        ref = m.standard.pzero_iso
        self.assertAlmostEqual(res, ref, delta=0.1)
        ref = m.constants.pzero
        self.assertAlmostEqual(res, ref, delta=5)

    def test_iso_table(self):
        for it in iso_tab.index:
            res = m.standard.p_iso(iso_tab["H"][it], gpm=True, hPa=False)
            ref = iso_tab["p"][it]
            self.assertAlmostEqual(res, ref, delta=ref*6.5E-5)

    def test_table_continuity(self):
        for h in m.standard._tab4["H_b"][1:-1]:
            res1 = m.standard.p_iso(h - 0.1,
                                    gpm=True, hPa=False)
            res2 = m.standard.p_iso(h + 0.1,
                                    gpm=True, hPa=False)
            self.assertAlmostEqual(res1, res2, delta=res1*3.5E-5)


class Test_temperature(unittest.TestCase):
    def test_iso_zero(self):
        res = m.standard.T_iso(0)
        ref = m.standard._tab4["T_b"][1]
        self.assertAlmostEqual(res, ref, delta=0.01)

    def test_iso_table(self):
        for it in iso_tab.index:
            res = m.standard.T_iso(iso_tab["H"][it], gpm=True, Kelvin=True)
            ref = iso_tab["T"][it]
            self.assertAlmostEqual(res, ref, delta=0.007)

    def test_table_continuity(self):
        for h in m.standard._tab4["H_b"][1:-1]:
            res1 = m.standard.T_iso(h - 0.1, gpm=True, Kelvin=False)
            res2 = m.standard.T_iso(h + 0.1, gpm=True, Kelvin=False)
            self.assertAlmostEqual(res1, res2, delta=0.002)


class Test_density(unittest.TestCase):

    def test_iso_table(self):
        for it in iso_tab.index:
            res = m.standard.rho_iso(iso_tab["H"][it], gpm=True)
            ref = iso_tab["rho"][it]
            self.assertAlmostEqual(res, ref, delta=0.000075)

    def test_table_continuity(self):
        for h in m.standard._tab4["H_b"][1:-1]:
            res1 = m.standard.rho_iso(h - 0.1, gpm=True)
            res2 = m.standard.rho_iso(h + 0.1, gpm=True)
            self.assertAlmostEqual(res1, res2, delta=0.000075)


class Test_altitude(unittest.TestCase):

    def test_iso_table(self):
        for it in iso_tab.index:
            res = m.standard.altitude(iso_tab["p"][it], gpm=True)
            ref = iso_tab["H"][it]
            self.assertAlmostEqual(res, ref, delta=max(0.5, ref*1E-5))

    def test_table_continuity(self):
        for i in m.standard._tab4.index[1:-1]:
            res1 = m.standard.altitude(m.standard._tab5["p"][i] * 0.99999,
                                       gpm=True)
            res2 = m.standard.altitude(m.standard._tab5["p"][i] * 1.00001,
                                       gpm=True)
            self.assertAlmostEqual(res1, res2, delta=max(0.5, res1*1E-5))
