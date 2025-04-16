#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 18:39:22 2019

@author: druee
"""

import unittest

import meteolib as m


class Test_dir2uv(unittest.TestCase):
    """
    def dir2uv(ff,dd)
    """

    def test_angles(self):
        res, ref = m.wind.dir2uv(1.,  0.), (0., -1.)
        for s, f in zip(res, ref):
            self.assertAlmostEqual(s, f, places=4)
        res, ref = m.wind.dir2uv(1., 180.), (0.,  1.)
        for s, f in zip(res, ref):
            self.assertAlmostEqual(s, f, places=4)
        res, ref = m.wind.dir2uv(1., 270.), (1.,  0.)
        for s, f in zip(res, ref):
            self.assertAlmostEqual(s, f, places=4)
        res, ref = m.wind.dir2uv(1., 360.), (0., -1.)
        for s, f in zip(res, ref):
            self.assertAlmostEqual(s, f, places=4)

    def test_speeds(self):
        with self.assertRaises(ValueError):
            m.wind.dir2uv(-1., 0.)
        res, ref = m.wind.dir2uv(2., 270.), (2.,  0.)
        for s, f in zip(res, ref):
            self.assertAlmostEqual(s, f, places=4)
        res, ref = m.wind.dir2uv(2., 45.), (-1.4142,  -1.4142)
        for s, f in zip(res, ref):
            self.assertAlmostEqual(s, f, places=4)

    def test_shape(self):
        with self.assertRaises(ValueError):
            m.wind.dir2uv([0], [0, 0])
        res, ref = m.wind.dir2uv([2.], [45.]), (-1.4142,  -1.4142)
        for s, f in zip(res, ref):
            self.assertAlmostEqual(s, f, places=4)
        res, ref = m.wind.dir2uv(
            [2., 2.], [45., 45.]), ([-1.4142,  -1.4142], [-1.4142,  -1.4142])
        for s, f in zip(res, ref):
            for a, b in zip(s, f):
                self.assertAlmostEqual(a, b, places=4)
