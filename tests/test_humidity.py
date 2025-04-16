#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 18:03:56 2019

@author: druee
"""

import unittest

import meteolib as m


class Test_magnus_w(unittest.TestCase):
    '''
    magnus_w(t, Kelvin=None, p=None, hPa=False)
    '''

    def test_2p_auto(self):
        """
        0°C -> 611.2 Pa
        20°C > 2336.1 Pa
        """
        self.assertAlmostEqual(m.humidity.magnus_w(0.), 611.2, places=1)
        self.assertAlmostEqual(m.humidity.magnus_w(20.), 2332.6, places=1)

    def test_2p_C(self):
        """
        0°C -> 611.2 Pa
        20°C > 2336.1 Pa
        """
        self.assertAlmostEqual(
            m.humidity.magnus_w(0., False), 611.2, places=1)
        self.assertAlmostEqual(m.humidity.magnus_w(
            20., False), 2332.6, places=1)

    def test_2p_K(self):
        """
        0°C -> 611.2 Pa
        20°C > 2336.1 Pa
        """
        self.assertAlmostEqual(m.humidity.magnus_w(
            273.15, True), 611.2, places=1)
        self.assertAlmostEqual(m.humidity.magnus_w(
            293.15, True), 2332.6, places=1)

    def test_2p_C_hPa(self):
        """
        0°C -> 6.112 hPa
        20°C > 23.361 hPa
        """
        self.assertAlmostEqual(m.humidity.magnus_w(
            0., False, hPa=True), 6.112, places=3)
        self.assertAlmostEqual(m.humidity.magnus_w(
            20., False, hPa=True), 23.326, places=3)

    def test_2p_C_Pa(self):
        """
        0°C -> 611.2 Pa
        20°C > 2336.1 Pa
        """
        self.assertAlmostEqual(m.humidity.magnus_w(
            0., False, hPa=False), 611.2, places=1)
        self.assertAlmostEqual(m.humidity.magnus_w(
            20., False, hPa=False), 2332.6, places=1)

    def test_2p_C_1013_hPa(self):
        """
        0°C -> 6.141 Pa *
        20°C > 23.436 Pa *
        """
        self.assertAlmostEqual(m.humidity.magnus_w(
            0., False, p=1013, hPa=True), 6.141, places=2)
        self.assertAlmostEqual(m.humidity.magnus_w(
            20., False, p=1013, hPa=True), 23.436, places=1)


class Test_magnus_i(unittest.TestCase):
    '''
    magnus_i(t, Kelvin=None, p=None, hPa=False)
    '''

    def test_2p_auto(self):
        """
        0°C-> 611.2 Pa
        -20°C > 103.3 Pa
        """
        self.assertAlmostEqual(m.humidity.magnus_i(0.), 611.2, places=1)
        self.assertAlmostEqual(m.humidity.magnus_i(-20.), 103.3, places=1)

    def test_2p_C(self):
        """
        0°C-> 611.2 Pa
        -20°C > 103.3 Pa
        """
        self.assertAlmostEqual(
            m.humidity.magnus_i(0., False), 611.2, places=1)
        self.assertAlmostEqual(
            m.humidity.magnus_i(-20., False), 103.3, places=1)

    def test_2p_K(self):
        """
        0°C-> 611.2 Pa
        -20°C > 103.3 Pa
        """
        self.assertAlmostEqual(m.humidity.magnus_i(
            273.15, True), 611.2, places=1)
        self.assertAlmostEqual(m.humidity.magnus_i(
            253.15, True), 103.3, places=1)

    def test_2p_C_hPa(self):
        """
        0°C-> 6.112 hPa
        -20°C > 1.033 hPa
        """
        self.assertAlmostEqual(m.humidity.magnus_i(
            0., False, hPa=True), 6.112, places=3)
        self.assertAlmostEqual(m.humidity.magnus_i(
            -20., False, hPa=True), 1.0326, places=3)

    def test_2p_C_Pa(self):
        """
        0°C-> 611.2 Pa
        -20°C > 103.3 Pa
        """
        self.assertAlmostEqual(m.humidity.magnus_i(
            0., False, hPa=False), 611.2, places=1)
        self.assertAlmostEqual(m.humidity.magnus_i(
            -20., False, hPa=False), 103.3, places=1)

    def test_2p_C_1013_hPa(self):
        """
        0°C-> 611.2 Pa *
        -20°C > 1.032 Pa *
        """
        self.assertAlmostEqual(m.humidity.magnus_i(
            0., False, p=1013, hPa=True), 6.112, places=1)
        self.assertAlmostEqual(m.humidity.magnus_i(
            -20., False, p=1013, hPa=True), 1.032, places=1)


class Test_tetens_w(unittest.TestCase):
    '''
    tetens_w(t, Kelvin=None, p=None, hPa=False)assertErrorEqual
    '''

    def test_2p_auto(self):
        self.assertAlmostEqual(m.humidity.tetens_w(0.), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.tetens_w(20.), 2338., delta=1.)

    def test_2p_C(self):
        self.assertAlmostEqual(
            m.humidity.tetens_w(0., False), 610.78, delta=1.)
        self.assertAlmostEqual(
            m.humidity.tetens_w(20., False), 2338., delta=1.)

    def test_2p_K(self):
        self.assertAlmostEqual(m.humidity.tetens_w(
            273.15, True), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.tetens_w(
            293.15, True), 2338., delta=1.)

    def test_2p_C_hPa(self):
        self.assertAlmostEqual(m.humidity.tetens_w(
            0., False, hPa=True), 6.1078, delta=0.01)
        self.assertAlmostEqual(m.humidity.tetens_w(
            20., False, hPa=True), 23.38, delta=0.01)

    def test_2p_C_Pa(self):
        self.assertAlmostEqual(m.humidity.tetens_w(
            0., False, hPa=False), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.tetens_w(
            20., False, hPa=False), 2338., delta=1.)

    def test_2p_C_1013_hPa(self):
        self.assertAlmostEqual(m.humidity.tetens_w(
            0., False, p=1013, hPa=True), 6.1180, delta=0.02)
        self.assertAlmostEqual(m.humidity.tetens_w(
            20., False, p=1013, hPa=True), 23.49, delta=0.02)


class Test_tetens_i(unittest.TestCase):
    '''
    tetens_i(t, Kelvin=None, p=None, hPa=False)
    '''

    def test_2p_auto(self):
        self.assertAlmostEqual(m.humidity.tetens_i(
            0.), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.tetens_i(
            -20.), 103., delta=1.)

    def test_2p_C(self):
        self.assertAlmostEqual(m.humidity.tetens_i(
            0., False), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.tetens_i(
            -20., False), 103., delta=1.)

    def test_2p_K(self):
        self.assertAlmostEqual(m.humidity.tetens_i(
            273.15, True), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.tetens_i(
            253.15, True), 103., delta=1.)

    def test_2p_C_hPa(self):
        self.assertAlmostEqual(m.humidity.tetens_i(
            0., False, hPa=True), 6.1078, delta=0.01)
        self.assertAlmostEqual(m.humidity.tetens_i(
            -20., False, hPa=True), 1.03, delta=0.01)

    def test_2p_C_Pa(self):
        self.assertAlmostEqual(m.humidity.tetens_i(
            0., False, hPa=False), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.tetens_i(
            -20., False, hPa=False), 103., delta=1.)

    def test_2p_C_1013_hPa(self):
        self.assertAlmostEqual(m.humidity.tetens_i(
            0., False, p=1013, hPa=True), 6.1180, delta=0.02)
        self.assertAlmostEqual(m.humidity.tetens_i(
            -20., False, p=1013, hPa=True), 1.03, delta=0.02)


class Test_goff_gratch_w(unittest.TestCase):
    '''
    goff_gratch_w(t, Kelvin=None, p=None, hPa=False)assertErrorEqual
    '''

    def test_2p_auto(self):
        self.assertAlmostEqual(m.humidity.goff_gratch_w(
            0.), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.goff_gratch_w(
            20.), 2338., delta=1.)

    def test_2p_C(self):
        self.assertAlmostEqual(m.humidity.goff_gratch_w(
            0., False), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.goff_gratch_w(
            20., False), 2338., delta=1.)

    def test_2p_K(self):
        self.assertAlmostEqual(m.humidity.goff_gratch_w(
            273.15, True), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.goff_gratch_w(
            293.15, True), 2338., delta=1.)

    def test_2p_C_hPa(self):
        self.assertAlmostEqual(m.humidity.goff_gratch_w(
            0., False, hPa=True), 6.1078, delta=0.01)
        self.assertAlmostEqual(m.humidity.goff_gratch_w(
            20., False, hPa=True), 23.38, delta=0.01)

    def test_2p_C_Pa(self):
        self.assertAlmostEqual(m.humidity.goff_gratch_w(
            0., False, hPa=False), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.goff_gratch_w(
            20., False, hPa=False), 2338., delta=1.)

    def test_2p_C_1013_hPa(self):
        self.assertAlmostEqual(m.humidity.goff_gratch_w(
            0., False, p=1013, hPa=True), 6.1180, delta=0.02)
        self.assertAlmostEqual(m.humidity.goff_gratch_w(
            20., False, p=1013, hPa=True), 23.49, delta=0.02)


class Test_goff_gratch_i(unittest.TestCase):
    '''
    goff_gratch_i(t, Kelvin=None, p=None, hPa=False)
    '''

    def test_2p_auto(self):
        self.assertAlmostEqual(m.humidity.goff_gratch_i(
            0.), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.goff_gratch_i(
            -20.), 103., delta=1.)

    def test_2p_C(self):
        self.assertAlmostEqual(m.humidity.goff_gratch_i(
            0., False), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.goff_gratch_i(
            -20., False), 103., delta=1.)

    def test_2p_K(self):
        self.assertAlmostEqual(m.humidity.goff_gratch_i(
            273.15, True), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.goff_gratch_i(
            253.15, True), 103., delta=1.)

    def test_2p_C_hPa(self):
        self.assertAlmostEqual(m.humidity.goff_gratch_i(
            0., False, hPa=True), 6.1078, delta=0.01)
        self.assertAlmostEqual(m.humidity.goff_gratch_i(
            -20., False, hPa=True), 1.03, delta=0.01)

    def test_2p_C_Pa(self):
        self.assertAlmostEqual(m.humidity.goff_gratch_i(
            0., False, hPa=False), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.goff_gratch_i(
            -20., False, hPa=False), 103., delta=1.)

    def test_2p_C_1013_hPa(self):
        self.assertAlmostEqual(m.humidity.goff_gratch_i(
            0., False, p=1013, hPa=True), 6.1180, delta=0.02)
        self.assertAlmostEqual(m.humidity.goff_gratch_i(
            -20., False, p=1013, hPa=True), 1.03, delta=0.02)


class Test_hyland_wexler_w(unittest.TestCase):
    '''
    hyland_wexler_w(t, Kelvin=None, p=None, hPa=False)assertErrorEqual
    '''

    def test_2p_auto(self):
        self.assertAlmostEqual(m.humidity.hyland_wexler_w(
            0.), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.hyland_wexler_w(
            20.), 2338., delta=1.)

    def test_2p_C(self):
        self.assertAlmostEqual(m.humidity.hyland_wexler_w(
            0., False), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.hyland_wexler_w(
            20., False), 2338., delta=1.)

    def test_2p_K(self):
        self.assertAlmostEqual(m.humidity.hyland_wexler_w(
            273.15, True), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.hyland_wexler_w(
            293.15, True), 2338., delta=1.)

    def test_2p_C_hPa(self):
        self.assertAlmostEqual(m.humidity.hyland_wexler_w(
            0., False, hPa=True), 6.1078, delta=0.01)
        self.assertAlmostEqual(m.humidity.hyland_wexler_w(
            20., False, hPa=True), 23.38, delta=0.01)

    def test_2p_C_Pa(self):
        self.assertAlmostEqual(m.humidity.hyland_wexler_w(
            0., False, hPa=False), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.hyland_wexler_w(
            20., False, hPa=False), 2338., delta=1.)

    def test_2p_C_1013_hPa(self):
        self.assertAlmostEqual(m.humidity.hyland_wexler_w(
            0., False, p=1013, hPa=True), 6.1180, delta=0.025)
        self.assertAlmostEqual(m.humidity.hyland_wexler_w(
            20., False, p=1013, hPa=True), 23.49, delta=0.025)


class Test_hyland_wexler_i(unittest.TestCase):
    '''
    hyland_wexler_i(t, Kelvin=None, p=None, hPa=False)
    '''

    def test_2p_auto(self):
        self.assertAlmostEqual(m.humidity.hyland_wexler_i(
            0.), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.hyland_wexler_i(
            -20.), 103., delta=1.)

    def test_2p_C(self):
        self.assertAlmostEqual(m.humidity.hyland_wexler_i(
            0., False), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.hyland_wexler_i(
            -20., False), 103., delta=1.)

    def test_2p_K(self):
        self.assertAlmostEqual(m.humidity.hyland_wexler_i(
            273.15, True), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.hyland_wexler_i(
            253.15, True), 103., delta=1.)

    def test_2p_C_hPa(self):
        self.assertAlmostEqual(m.humidity.hyland_wexler_i(
            0., False, hPa=True), 6.1078, delta=0.01)
        self.assertAlmostEqual(m.humidity.hyland_wexler_i(
            -20., False, hPa=True), 1.03, delta=0.01)

    def test_2p_C_Pa(self):
        self.assertAlmostEqual(m.humidity.hyland_wexler_i(
            0., False, hPa=False), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.hyland_wexler_i(
            -20., False, hPa=False), 103., delta=1.)

    def test_2p_C_1013_hPa(self):
        self.assertAlmostEqual(m.humidity.hyland_wexler_i(
            0., False, p=1013, hPa=True), 6.1180, delta=0.025)
        self.assertAlmostEqual(m.humidity.hyland_wexler_i(
            -20., False, p=1013, hPa=True), 1.03, delta=0.025)


class Test_sonntag_w(unittest.TestCase):
    '''
    sonntag_w(t, Kelvin=None, p=None, hPa=False)assertErrorEqual
    '''

    def test_2p_auto(self):
        self.assertAlmostEqual(m.humidity.sonntag_w(
            0.), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.sonntag_w(
            20.), 2338., delta=1.5)

    def test_2p_C(self):
        self.assertAlmostEqual(m.humidity.sonntag_w(
            0., False), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.sonntag_w(
            20., False), 2338., delta=1.5)

    def test_2p_K(self):
        self.assertAlmostEqual(m.humidity.sonntag_w(
            273.15, True), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.sonntag_w(
            293.15, True), 2338., delta=1.5)

    def test_2p_C_hPa(self):
        self.assertAlmostEqual(m.humidity.sonntag_w(
            0., False, hPa=True), 6.1078, delta=0.01)
        self.assertAlmostEqual(m.humidity.sonntag_w(
            20., False, hPa=True), 23.38, delta=0.015)

    def test_2p_C_Pa(self):
        self.assertAlmostEqual(m.humidity.sonntag_w(
            0., False, hPa=False), 610.78, delta=1.5)
        self.assertAlmostEqual(m.humidity.sonntag_w(
            20., False, hPa=False), 2338., delta=1.5)

    def test_2p_C_1013_hPa(self):
        self.assertAlmostEqual(m.humidity.sonntag_w(
            0., False, p=1013, hPa=True), 6.1180, delta=0.025)
        self.assertAlmostEqual(m.humidity.sonntag_w(
            20., False, p=1013, hPa=True), 23.49, delta=0.025)


class Test_sonntag_i(unittest.TestCase):
    '''
    sonntag_i(t, Kelvin=None, p=None, hPa=False)
    '''

    def test_2p_auto(self):
        self.assertAlmostEqual(m.humidity.sonntag_i(
            0.), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.sonntag_i(
            -20.), 103., delta=1.)

    def test_2p_C(self):
        self.assertAlmostEqual(
            m.humidity.sonntag_i(0., False), 610.78, delta=1.)
        self.assertAlmostEqual(
            m.humidity.sonntag_i(-20., False), 103., delta=1.)

    def test_2p_K(self):
        self.assertAlmostEqual(m.humidity.sonntag_i(
            273.15, True), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.sonntag_i(
            253.15, True), 103., delta=1.)

    def test_2p_C_hPa(self):
        self.assertAlmostEqual(m.humidity.sonntag_i(
            0., False, hPa=True), 6.1078, delta=0.01)
        self.assertAlmostEqual(m.humidity.sonntag_i(
            -20., False, hPa=True), 1.03, delta=0.01)

    def test_2p_C_Pa(self):
        self.assertAlmostEqual(m.humidity.sonntag_i(
            0., False, hPa=False), 610.78, delta=1.)
        self.assertAlmostEqual(m.humidity.sonntag_i(
            -20., False, hPa=False), 103., delta=1.)

    def test_2p_C_1013_hPa(self):
        self.assertAlmostEqual(m.humidity.sonntag_i(
            0., False, p=1013, hPa=True), 6.1180, delta=0.025)
        self.assertAlmostEqual(m.humidity.sonntag_i(
            -20., False, p=1013, hPa=True), 1.03, delta=0.025)


class Test_tdew(unittest.TestCase):
    '''
    tdew(e, p=None, hPa=False, Kelvin=True)
    '''

    def test_2p_auto(self):
        self.assertAlmostEqual(m.humidity.tdew(
            611.2), 273.15, places=1)
        self.assertAlmostEqual(m.humidity.tdew(
            2332.6), 293.15, places=1)

    def test_2p_C(self):
        self.assertAlmostEqual(m.humidity.tdew(
            611.2, Kelvin=False), 0.0, places=1)
        self.assertAlmostEqual(m.humidity.tdew(
            2332.6, Kelvin=False), 20.0, places=1)

    def test_2p_K(self):
        self.assertAlmostEqual(m.humidity.tdew(
            611.2, Kelvin=True), 273.15, places=1)
        self.assertAlmostEqual(m.humidity.tdew(
            2332.6, Kelvin=True), 293.15, places=1)

    def test_2p_C_hPa(self):
        self.assertAlmostEqual(m.humidity.tdew(
            6.1078, hPa=True, Kelvin=True), 273.15, places=1)
        self.assertAlmostEqual(m.humidity.tdew(
            23.326, hPa=True, Kelvin=True), 293.15, places=1)

    def test_2p_C_Pa(self):
        self.assertAlmostEqual(m.humidity.tdew(
            611.2, Kelvin=False, hPa=False), 0., places=1)
        self.assertAlmostEqual(m.humidity.tdew(
            2332.6, Kelvin=False, hPa=False), 20., places=1)

    def test_2p_C_1013_hPa(self):
        self.assertAlmostEqual(m.humidity.tdew(
            614.1, 101315., Kelvin=False, hPa=False), 0., places=1)
        self.assertAlmostEqual(m.humidity.tdew(
            2343.6, 101315., Kelvin=False, hPa=False), 20., places=1)


class Test_tfrost(unittest.TestCase):
    '''
    tfrost(e, p=None, hPa=False, Kelvin=True)
    '''

    def test_2p_auto(self):
        self.assertAlmostEqual(m.humidity.tfrost(
            611.2), 273.15, places=1)
        self.assertAlmostEqual(m.humidity.tfrost(
            103.3), 253.15, places=1)

    def test_2p_C(self):
        self.assertAlmostEqual(m.humidity.tfrost(
            611.2, Kelvin=False), 0.0, places=1)
        self.assertAlmostEqual(m.humidity.tfrost(
            103.3, Kelvin=False), -20.0, places=1)

    def test_2p_K(self):
        self.assertAlmostEqual(m.humidity.tfrost(
            611.2, Kelvin=True), 273.15, places=1)
        self.assertAlmostEqual(m.humidity.tfrost(
            103.3, Kelvin=True), 253.15, places=1)

    def test_2p_C_hPa(self):
        self.assertAlmostEqual(m.humidity.tfrost(
            6.1078, hPa=True, Kelvin=True), 273.15, places=1)
        self.assertAlmostEqual(m.humidity.tfrost(
            1.033, hPa=True, Kelvin=True), 253.15, places=1)

    def test_2p_C_Pa(self):
        self.assertAlmostEqual(m.humidity.tfrost(
            611.2, Kelvin=False, hPa=False), 0., places=1)
        self.assertAlmostEqual(m.humidity.tfrost(
            103.3, Kelvin=False, hPa=False), -20., places=1)

    def test_2p_C_1013_hPa(self):
        self.assertAlmostEqual(m.humidity.tfrost(
            614.1, 101315., Kelvin=False, hPa=False), 0., places=1)
        self.assertAlmostEqual(m.humidity.tfrost(
            103.3, 101315., Kelvin=False, hPa=False), -20., places=1)


class Test_psychro(unittest.TestCase):
    '''
    psychro(t, tw, Kelvin=None, p=None, hPa=False):
    '''

    def test_4p_auto(self):
        self.assertAlmostEqual(m.humidity.psychro(
            273.15, 273.15), 611.2, delta=1.)
        self.assertAlmostEqual(m.humidity.psychro(
            273.15, 272.15), 611.2*0.821, delta=1.)
        self.assertAlmostEqual(m.humidity.psychro(
            293.15, 293.15), 2332.6, delta=1.)
        self.assertAlmostEqual(m.humidity.psychro(
            293.15, 283.15), 2332.6*0.236, delta=10.)

    def test_4p_hPa(self):
        self.assertAlmostEqual(m.humidity.psychro(
            273.15, 273.15, hPa=True), 6.112, delta=0.01)
        self.assertAlmostEqual(m.humidity.psychro(
            273.15, 272.15, hPa=True), 6.112*0.821, delta=0.01)
        self.assertAlmostEqual(m.humidity.psychro(
            293.15, 293.15, hPa=True), 23.326, delta=0.01)
        self.assertAlmostEqual(m.humidity.psychro(
            293.15, 283.15, hPa=True), 23.326*0.236, delta=0.1)

    def test_4p_C(self):
        self.assertAlmostEqual(m.humidity.psychro(
            0., 0.), 611.2, delta=1.)
        self.assertAlmostEqual(m.humidity.psychro(
            0., -1.), 611.2*0.821, delta=1.)
        self.assertAlmostEqual(m.humidity.psychro(
            20., 20.), 2332.6, delta=1.)
        self.assertAlmostEqual(m.humidity.psychro(
            20., 10.), 2332.6*0.236, delta=10.)


class Test_psychro_ice(unittest.TestCase):
    '''
    psychro_ice(t, tw, Kelvin=None, p=None, hPa=False):
    '''

    def test_4p_auto(self):
        self.assertAlmostEqual(m.humidity.psychro_ice(
            273.15, 273.15), 611.2, delta=1.)
        self.assertAlmostEqual(m.humidity.psychro_ice(
            273.15, 272.15), 611.2*0.821, delta=10.)
        self.assertAlmostEqual(m.humidity.psychro_ice(
            253.15, 253.15), 103.26, delta=1.)
        self.assertAlmostEqual(m.humidity.psychro_ice(
            253.15, 252.15), 103.26*0.344, delta=10.)

    def test_4p_hPa(self):
        self.assertAlmostEqual(m.humidity.psychro_ice(
            273.15, 273.15, hPa=True), 6.112, delta=0.01)
        self.assertAlmostEqual(m.humidity.psychro_ice(
            273.15, 272.15, hPa=True), 6.112*0.821, delta=0.1)
        self.assertAlmostEqual(m.humidity.psychro_ice(
            253.15, 253.15, hPa=True), 1.0326, delta=0.01)
        self.assertAlmostEqual(m.humidity.psychro_ice(
            253.15, 252.15, hPa=True), 1.0326*0.344, delta=0.1)

    def test_4p_C(self):
        self.assertAlmostEqual(m.humidity.psychro_ice(
            0., 0.), 611.2, delta=1.)
        self.assertAlmostEqual(m.humidity.psychro_ice(
            0., -1.), 611.2*0.821, delta=10.)
        self.assertAlmostEqual(m.humidity.psychro_ice(
            -20., -20.), 103.26, delta=1.)
        self.assertAlmostEqual(m.humidity.psychro_ice(
            -20., -21.), 103.26*0.344, delta=10.)


class Test_relhum(unittest.TestCase):
    '''
    relhum(e, ew, percent=True)
    '''

    def test_3p_auto(self):
        self.assertAlmostEqual(m.humidity.relhum(
            0., 611.2), 0.0, delta=0.001)
        self.assertAlmostEqual(m.humidity.relhum(
            305.6, 611.2), 0.5, delta=0.001)
        self.assertAlmostEqual(m.humidity.relhum(
            611.2, 611.2), 1.0, delta=0.001)

    def test_3p_percent(self):
        self.assertAlmostEqual(m.humidity.relhum(
            0., 611.2, percent=True),  0.0, delta=0.1)
        self.assertAlmostEqual(m.humidity.relhum(
            305.6, 611.2, percent=True), 50.0, delta=0.1)
        self.assertAlmostEqual(m.humidity.relhum(
            611.2, 611.2, percent=True), 100.0, delta=0.1)


class Test_mixr(unittest.TestCase):
    '''
    mixr(e, p, gkg=False)
    '''

    def test_2p_auto(self):
        self.assertAlmostEqual(m.humidity.mixr(0., 1013.2), 0.0, delta=0.0001)
        self.assertAlmostEqual(m.humidity.mixr(
            6.112, 1013.2), 0.00375, delta=0.0001)

    def test_2p_hPa(self):
        self.assertAlmostEqual(m.humidity.mixr(
            611.2, 101320), 0.00375, delta=0.0001)

    def test_3p_gkg(self):
        self.assertAlmostEqual(m.humidity.mixr(
            0., 1013.2, gkg=True), 0.0, delta=0.1)
        self.assertAlmostEqual(m.humidity.mixr(
            6.112, 1013.2, gkg=True), 3.75, delta=0.1)
        self.assertAlmostEqual(m.humidity.mixr(
            611.2, 101320, gkg=True), 3.75, delta=0.1)


class Test_inv_mixr(unittest.TestCase):
    '''
    inv_mixr(m, p, hPa=False, gkg=False)
    '''

    def test_2p_auto(self):
        self.assertAlmostEqual(m.humidity.inv_mixr(
            0., 101320), 0.0, delta=1.)
        self.assertAlmostEqual(m.humidity.inv_mixr(
            0.00375, 101320), 611.2, delta=1.)

    def test_2p_hPa(self):
        self.assertAlmostEqual(m.humidity.inv_mixr(
            0., 1013.2, hPa=True), 0.0, delta=0.01)
        self.assertAlmostEqual(m.humidity.inv_mixr(
            0.00375, 1013.2, hPa=True), 6.112, delta=0.01)

    def test_3p_gkg(self):
        self.assertAlmostEqual(m.humidity.inv_mixr(
            0., 101320, gkg=True), 0.0, delta=1.)
        self.assertAlmostEqual(m.humidity.inv_mixr(
            3.75, 101320, gkg=True), 611.2, delta=1.)


class Test_spech(unittest.TestCase):
    '''
    spech(e, p, gkg=False)
    '''

    def test_2p_auto(self):
        self.assertAlmostEqual(m.humidity.spech(
            0., 1013.2), 0.0, delta=0.0001)
        self.assertAlmostEqual(m.humidity.spech(
            6.112, 1013.2), 0.00374, delta=0.0001)

    def test_2p_hPa(self):
        self.assertAlmostEqual(m.humidity.spech(
            611.2, 101320), 0.00374, delta=0.0001)

    def test_3p_gkg(self):
        self.assertAlmostEqual(m.humidity.spech(
            0., 1013.2, gkg=True), 0.0, delta=0.1)
        self.assertAlmostEqual(m.humidity.spech(
            6.112, 1013.2, gkg=True), 3.74, delta=0.1)
        self.assertAlmostEqual(m.humidity.spech(
            611.2, 101320, gkg=True), 3.74, delta=0.1)


class Test_inv_spech(unittest.TestCase):
    '''
    inv_spech(m, p, hPa=False, gkg=False)
    '''

    def test_2p_auto(self):
        self.assertAlmostEqual(
            m.humidity.inv_spech(0., 101320), 0.0, delta=1.)
        self.assertAlmostEqual(m.humidity.inv_spech(
            0.00374, 101320), 611.2, delta=1.)

    def test_2p_hPa(self):
        self.assertAlmostEqual(m.humidity.inv_spech(
            0., 1013.2, hPa=True), 0.0, delta=0.01)
        self.assertAlmostEqual(m.humidity.inv_spech(
            0.00374, 1013.2, hPa=True), 6.112, delta=0.01)

    def test_3p_gkg(self):
        self.assertAlmostEqual(m.humidity.inv_spech(
            0., 101320, gkg=True), 0.0, delta=1.)
        self.assertAlmostEqual(m.humidity.inv_spech(
            3.74, 101320, gkg=True), 611.2, delta=1.)


class Test_to_rh(unittest.TestCase):
    '''
    Humidity(t, p,
             e=None, rh=None, td=None, tf=None,
             tw=None, ti=None, q=None, m=None,
             Kelvin=None, hPa=False, gkg=False, percent=False):
    '''
    # test input variants

    def test_rh_t_auto(self):
        self.assertAlmostEqual(m.humidity.Humidity(
            rh=0.50, t=0.).e, 305.6, delta=1.)

    def test_e_t_auto(self):
        self.assertAlmostEqual(m.humidity.Humidity(
            e=611.2, t=0.).e, 611.2, delta=1.)

    def test_e_ew_auto(self):
        self.assertAlmostEqual(m.humidity.Humidity(
            e=611.2, ew=611.2).e, 611.2, delta=1.)

    def test_e_w_auto(self):
        self.assertAlmostEqual(m.humidity.Humidity(
            e=611.2, ew=611.2).e, 611.2, delta=1.)

    def test_e_t_percent(self):
        self.assertAlmostEqual(m.humidity.Humidity(
            e=611.2, t=0., percent=True).e, 611.2, delta=1.)

    def test_e_t_hPa(self):
        self.assertAlmostEqual(m.humidity.Humidity(
            e=6.112, t=0., hPa=True).e, 6.112, delta=0.01)

    def test_e_t_Kelvin(self):
        self.assertAlmostEqual(m.humidity.Humidity(
            e=611.2, t=273.15, Kelvin=True).e, 611.2, delta=1.)

    def test_td_t_auto(self):
        self.assertAlmostEqual(m.humidity.Humidity(
            td=9.327, t=20.).e, 1172., delta=1.)

    def test_tf_t_auto(self):
        self.assertAlmostEqual(m.humidity.Humidity(
            tf=8.13, t=20.).e, 1172., delta=1.)

    def test_tw_t_auto(self):
        self.assertAlmostEqual(m.humidity.Humidity(
            tw=13.87, t=20.).e, 1172., delta=1.)

    def test_ti_t_auto(self):
        self.assertAlmostEqual(m.humidity.Humidity(
            ti=-12, t=-10.).e, 100.8, delta=1.)

    def test_q_t_auto(self):
        self.assertAlmostEqual(m.humidity.Humidity(
            q=0.00714, p=101325, t=20.).e, 1172., delta=1.)

    def test_m_t_auto(self):
        self.assertAlmostEqual(m.humidity.Humidity(
            m=0.00719, p=101325, t=20.).e, 1172., delta=1.)
    # test execptions

    def test_params_notew(self):
        self.assertRaises(ValueError, m.humidity.Humidity, e=611.2)

    def test_params_tandew(self):
        self.assertRaises(ValueError, m.humidity.Humidity,
                          e=611.2, t=0., ew=611.2)

    def test_params_none(self):
        self.assertRaises(TypeError, m.humidity.Humidity,
                          e=611.2, t=0., Kelvin='blah')
        self.assertRaises(ValueError, m.humidity.Humidity,
                          e=611.2, t=0., hPa=None)
        self.assertRaises(ValueError, m.humidity.Humidity,
                          e=611.2, t=0., gkg=None)
        self.assertRaises(ValueError, m.humidity.Humidity,
                          e=611.2, t=0., percent=None)

    def test_rh_4p_auto(self):
        self.assertAlmostEqual(m.humidity.Humidity(
            rh=1., t=0., p=101320).rh(), 1., delta=0.01)
        self.assertAlmostEqual(m.humidity.Humidity(
            rh=0.5, t=0., p=101320).rh(), 0.5, delta=0.01)
        self.assertAlmostEqual(m.humidity.Humidity(
            rh=50., t=0., percent=True).rh(), 50., delta=1.)
        self.assertAlmostEqual(m.humidity.Humidity(
            rh=0.5, t=20., p=101320).rh(), 0.5, delta=0.01)

    def test_m_auto(self):
        self.assertAlmostEqual(m.humidity.Humidity(
            e=305.6, t=273.15, p=101320).m(), 0.00187, delta=0.0001)

    def test_m_gkg(self):
        self.assertAlmostEqual(m.humidity.Humidity(
            e=305.6, t=273.15, p=101320, gkg=True).m(), 1.87, delta=0.1)

    def test_q_auto(self):
        self.assertAlmostEqual(m.humidity.Humidity(
            e=305.6, t=273.15, p=101320).q(), 0.00187, delta=0.0001)

    def test_q_gkg(self):
        self.assertAlmostEqual(m.humidity.Humidity(
            e=305.6, t=273.15, p=101320, gkg=True).q(), 1.87, delta=0.1)

    def test_tw_4p_auto(self):
        self.assertAlmostEqual(m.humidity.Humidity(
            e=611.2, t=273.15).tw(), 273.15, delta=0.1)
        self.assertAlmostEqual(m.humidity.Humidity(
            e=611.2*0.821, t=273.15).tw(), 272.15, delta=0.1)
        self.assertAlmostEqual(m.humidity.Humidity(
            e=2332.6, t=293.15).tw(), 293.15, delta=0.1)
        self.assertAlmostEqual(m.humidity.Humidity(
            e=2332.6*0.236, t=293.15).tw(), 283.15, delta=0.1)

    def test_ti_4p_auto(self):
        self.assertAlmostEqual(m.humidity.Humidity(
            e=611.2, t=273.15).ti(), 273.15, delta=0.1)
        self.assertAlmostEqual(m.humidity.Humidity(
            e=611.2*0.821, t=273.15).ti(), 272.15, delta=0.1)
        self.assertAlmostEqual(m.humidity.Humidity(
            e=103.26, t=253.15).ti(), 253.15, delta=0.1)
        self.assertAlmostEqual(m.humidity.Humidity(
            e=103.26*0.344, t=253.15).ti(), 252.15, delta=0.1)

    def test_virt_inc_4p_auto(self):
        self.assertAlmostEqual(m.humidity.Humidity(
            t=0., td=0., p=101320).virt_inc(), 0.63, delta=0.01)
        self.assertAlmostEqual(m.humidity.Humidity(
            t=0., td=-1., p=101320).virt_inc(), 0.58, delta=0.01)
        self.assertAlmostEqual(m.humidity.Humidity(
            t=20., td=10., p=101320).virt_inc(), 1.35, delta=0.01)
        self.assertAlmostEqual(m.humidity.Humidity(
            t=20., td=20., p=101320).virt_inc(), 2.58, delta=0.05)

    def test_Tvirt_4p_auto(self):
        self.assertAlmostEqual(m.humidity.Humidity(
            t=0., td=0., p=101320).tvirt(), 273.15+0.63, delta=0.01)
        self.assertAlmostEqual(m.humidity.Humidity(
            t=0., td=-1., p=101320).tvirt(), 273.15+0.58, delta=0.01)
        self.assertAlmostEqual(m.humidity.Humidity(
            t=20., td=10., p=101320).tvirt(), 273.15+21.35, delta=0.01)
        self.assertAlmostEqual(m.humidity.Humidity(
            t=20., td=20., p=101320).tvirt(), 273.15+22.58, delta=0.05)

    def test_rhow_4p_auto(self):
        self.assertAlmostEqual(m.humidity.Humidity(
            rhow=0.0013, t=0.).rhow(), 0.0013, delta=0.0001)
        self.assertAlmostEqual(m.humidity.Humidity(
            rhow=0.0013, t=20.).rhow(), 0.0013, delta=0.0001)
        self.assertAlmostEqual(m.humidity.Humidity(
            rhow=0.013, t=0.).rhow(), 0.013, delta=0.0001)
        self.assertAlmostEqual(m.humidity.Humidity(
            rhow=0.013, t=20.).rhow(), 0.013, delta=0.0001)


if __name__ == '__main__':
    unittest.main()
