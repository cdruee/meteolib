#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 18:03:56 2019

@author: druee
"""

import unittest

import meteolib as m


class Test_to_K(unittest.TestCase):
    def test_OC_auto(self):
        """
        0°C -> 273.15 K
        """
        result = m.temperature._to_K(0.)
        self.assertEqual(result, m.constants.Tzero)

    def test_273C_auto(self):
        """
        273.15 K -> 273.15 K
        """
        result = m.temperature._to_K(m.constants.Tzero)
        self.assertEqual(result, m.constants.Tzero)

    def test_OC_C(self):
        """
        0°C -> 273.15 K
        """
        result = m.temperature._to_K(0., False)
        self.assertEqual(result, m.constants.Tzero)

    def test_273K_K(self):
        """
        273.15 K -> 273.15 K
        """
        result = m.temperature._to_K(m.constants.Tzero, True)
        self.assertEqual(result, m.constants.Tzero)


class Test_to_C(unittest.TestCase):
    def test_OC_auto(self):
        """
        0°C -> 0°C
        """
        result = m.temperature._to_C(0.)
        self.assertEqual(result, 0.)

    def test_273K_auto(self):
        """
        273.15 K -> 0°C
        """
        result = m.temperature._to_C(m.constants.Tzero)
        self.assertEqual(result, 0.)

    def test_OC_C(self):
        """
        0°C -> 0°C
        """
        result = m.temperature._to_C(0., False)
        self.assertEqual(result, 0.)

    def test_273K_K(self):
        """
        273.15 K -> 273.15 K
        """
        result = m.temperature._to_C(m.constants.Tzero, True)
        self.assertEqual(result, 0.)


class Test_KtoC(unittest.TestCase):
    def test_273K(self):
        """
        298.15 K -> 25°C
        """
        result = m.temperature.KtoC(298.15)
        self.assertEqual(result, 25.)


class Test_CtoK(unittest.TestCase):
    def test_273K(self):
        """
        273.15 K -> 0°C
        """
        result = m.temperature.CtoK(25.)
        self.assertEqual(result, 298.15)


class Test_FtoC(unittest.TestCase):
    def test_0F(self):
        """
        0 F = -17.77C
        """
        result = m.temperature.FtoC(0)
        self.assertAlmostEqual(result, -17.77, delta=0.01)

    def test_0C(self):
        """
        0 C = 32.0F
        """
        result = m.temperature.FtoC(32.0)
        self.assertAlmostEqual(result, 0., delta=0.01)


class Test_Tpot(unittest.TestCase):
    def test_pref_auto(self):
        """
        273.15 K, 100000 Pa -> 273.15 K
        """
        self.assertAlmostEqual(m.temperature.Tpot(
            273.15, 100000.), 273.15, delta=0.01)
        self.assertAlmostEqual(m.temperature.Tpot(
            0., 100000.), 273.15, delta=0.01)

    def test_pref_K(self):
        """
        273.15 K, 100000 Pa -> 273.15 K
        """
        self.assertAlmostEqual(m.temperature.Tpot(
            273.15, 100000., Kelvin=True), 273.15, delta=0.01)

    def test_pref_C(self):
        """
        0°C, 100000 Pa -> 0 °C
        """
        self.assertAlmostEqual(m.temperature.Tpot(
            0., 100000., Kelvin=True), 0., delta=0.01)

    def test_pref_hPa(self):
        """
        273.15 K, 100000 Pa -> 273.15 K
        """
        self.assertAlmostEqual(m.temperature.Tpot(
            273.15, 1000., hPa=True), 273.15, delta=0.01)

    def test_750_auto(self):
        """
        273.15 K, 100000 Pa -> 273.15 K
        """
        self.assertAlmostEqual(m.temperature.Tpot(
            293.15, 75000.), 318.28, delta=0.03)
        self.assertAlmostEqual(m.temperature.Tpot(
            20., 75000.), 318.28, delta=0.03)


class Test_inv_Tpot(unittest.TestCase):
    def test_pref_auto(self):
        """
        273.15 K, 100000 Pa -> 273.15 K
        """
        self.assertAlmostEqual(m.temperature.inv_Tpot(
            273.15, 100000.), 273.15, delta=0.01)
        self.assertAlmostEqual(m.temperature.inv_Tpot(
            0., 100000.), 273.15, delta=0.01)

    def test_pref_K(self):
        """
        273.15 K, 100000 Pa -> 273.15 K
        """
        self.assertAlmostEqual(m.temperature.inv_Tpot(
            273.15, 100000., Kelvin=True), 273.15, delta=0.01)

    def test_pref_C(self):
        """
        0°C, 100000 Pa -> 0 °C
        """
        self.assertAlmostEqual(m.temperature.inv_Tpot(
            0., 100000., Kelvin=True), 0., delta=0.01)

    def test_pref_hPa(self):
        """
        273.15 K, 100000 Pa -> 273.15 K
        """
        self.assertAlmostEqual(m.temperature.inv_Tpot(
            273.15, 1000., hPa=True), 273.15, delta=0.01)

    def test_750_auto(self):
        """
        273.15 K, 100000 Pa -> 273.15 K
        """
        self.assertAlmostEqual(m.temperature.inv_Tpot(
            318.28, 75000.), 293.15, delta=0.03)
        self.assertAlmostEqual(m.temperature.inv_Tpot(
            45.13, 75000.), 293.15, delta=0.03)


if __name__ == '__main__':
    unittest.main()
