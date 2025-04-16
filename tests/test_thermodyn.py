#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 02:13:26 2019

@author: clemens
"""

import unittest

import meteolib as m


class Test_gas_p(unittest.TestCase):
    '''
    def gas_p(rho, T, q=None, Kelvin=None, hPa=False, gas='air')
    '''

    def test_gas_p_dry(self):
        self.assertAlmostEqual(m.thermodyn.gas_p(
            1.25, 273.15, gas='dry'), 98009.6, delta=1.)

    def test_gas_p_moist(self):
        self.assertAlmostEqual(m.thermodyn.gas_p(
            1.25, 273.15, q=0, gas='moist'), 98009.6, delta=1.)
        self.assertAlmostEqual(m.thermodyn.gas_p(
            1.25, 273.15, q=2.E-3, gas='moist'), 98129.2, delta=1.)

    def test_gas_p_water(self):
        self.assertAlmostEqual(m.thermodyn.gas_p(
            1.25, 273.15, gas='water'), 157576.8, delta=1.)

    def test_gas_p_co2(self):
        self.assertAlmostEqual(m.thermodyn.gas_p(
            1.25, 273.15, gas='co2'), 64504.4, delta=1.)

    def test_gas_p_Celsius(self):
        self.assertAlmostEqual(m.thermodyn.gas_p(
            1.25, 0., Kelvin=False, gas='dry'), 98009.6, delta=1.)

    def test_gas_p_hPa(self):
        self.assertAlmostEqual(m.thermodyn.gas_p(
            1.25, 273.15, hPa=True, gas='dry'), 980.096, delta=0.01)

    def test_gas_p_gkg(self):
        self.assertAlmostEqual(m.thermodyn.gas_p(
            1.25, 273.15, q=2., gkg=True, gas='moist'), 98129.2, delta=1.)


class Test_gas_rho(unittest.TestCase):
    '''
    def gas_rho(p, T, q=None, Kelvin=None, hPa=False, gas='air')
    '''

    def test_gas_rho_dry(self):
        self.assertAlmostEqual(m.thermodyn.gas_rho(
            98009.6, 273.15, gas='dry'), 1.25, delta=0.01)

    def test_gas_rho_moist(self):
        self.assertAlmostEqual(m.thermodyn.gas_rho(
            98009.6, 273.15, q=0, gas='moist'), 1.25, delta=0.01)
        self.assertAlmostEqual(m.thermodyn.gas_rho(
            98129.2, 273.15, q=2.E-3, gas='moist'), 1.25, delta=0.01)

    def test_gas_rho_water(self):
        self.assertAlmostEqual(m.thermodyn.gas_rho(
            157576.8, 273.15, gas='water'), 1.25, delta=0.01)

    def test_gas_rho_co2(self):
        self.assertAlmostEqual(m.thermodyn.gas_rho(
            64504.4, 273.15, gas='co2'), 1.25, delta=0.01)

    def test_gas_rho_Celsius(self):
        self.assertAlmostEqual(m.thermodyn.gas_rho(
            98009.6, 0., Kelvin=False, gas='dry'), 1.25, delta=0.01)

    def test_gas_rho_hPa(self):
        self.assertAlmostEqual(m.thermodyn.gas_rho(
            980.096, 273.15, hPa=True, gas='dry'), 1.25, delta=0.01)

    def test_gas_rho_gkg(self):
        self.assertAlmostEqual(m.thermodyn.gas_rho(
            98129.2, 273.15, q=2., gkg=True, gas='moist'), 1.25, delta=0.01)


class Test_cpt(unittest.TestCase):
    '''
    def cpt(t, Kelvin=None):
    '''

    def test_cpt_auto(self):
        self.assertAlmostEqual(m.thermodyn.cpt(0.), m.constants.cp, delta=1.)


class Test_cpvt(unittest.TestCase):
    '''
    def cpvt(t, Kelvin=None):
    '''

    def test_cpvt_auto(self):
        self.assertAlmostEqual(m.thermodyn.cpvt(
            0.), m.constants.cpv, delta=10.)


class Test_cwt(unittest.TestCase):
    '''
    def cwt(t, Kelvin=None):
    '''

    def test_cwt_auto(self):
        self.assertAlmostEqual(m.thermodyn.cwt(
            0.), m.constants.cw, delta=0.01*m.constants.cw)


class Test_cit(unittest.TestCase):
    '''
    def cit(t, Kelvin=None):
    '''

    def test_cit_auto(self):
        self.assertAlmostEqual(m.thermodyn.cit(
            0.), m.constants.ci, delta=0.015*m.constants.ci)


class Test_Lvt(unittest.TestCase):
    '''
    def Lvt(t, Kelvin=None):
    '''
    # Smithonian tabels (1951), Table 92 (in cal/gK; J = 4.1839 cal)
    t92 = {
        0: 597.31,
        10: 591.7,
        20: 586.0,
        30: 580.4,
        40: 574.7,
        50: 569.0,
    }
    j92 = 4.1839E3

    def test_Lvt_auto(self):
        for t, c in self.t92.items():
            with self.subTest(t=t):
                self.assertAlmostEqual(m.thermodyn.Lvt(
                    t), c*self.j92, delta=0.005*c*self.j92)

    def test_Lvt_sources(self):
        for s in ['andreas', 'henderson-sellers',
                  'ecpack', 'dake', 'hyland-wexler']:
            for t, c in self.t92.items():
                with self.subTest(s=s, t=t):
                    self.assertAlmostEqual(
                        m.thermodyn.Lvt(t, source=s),
                        c*self.j92, delta=0.0075*c*self.j92)


class Test_Lst(unittest.TestCase):
    '''
    def Lst(t, Kelvin=None):
    '''

    def test_Lst_sources(self):
        for s in ['andreas', 'hyland-wexler']:
            t = 0
            with self.subTest(s=s, t=t):
                self.assertAlmostEqual(m.thermodyn.Lst(
                    t, source=s), m.constants.Ls, delta=0.005*m.constants.Ls)
