#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 18:39:22 2019

@author: druee
"""

import datetime as dt
import unittest

import meteolib as m

import pytz


CET = dt.timezone(dt.timedelta(hours=1))
# MST=dt.timezone(dt.timedelta(hours=-7))

MST = pytz.timezone('MST')


class Test_sunriseset(unittest.TestCase):
    """
    sun_rise_transit_set(time: dt.datetime, lat, lon, ele=0, pp=None, tk=None)
    """

    def test_spa_rts_values(self):
        res, ref = m.radiation.spa_rise_transit_set(
            dt.datetime(2021, 1, 4, 14, 0, 0, tzinfo=CET), 49.6, 6.9), (
            8.477979163121788, 12.624062013177186, 16.77403260994837)
        for s, f in zip(res, ref):
            self.assertAlmostEqual(s, f, places=3)

    def test_fast_rts_values(self):
        res, ref = m.radiation.fast_rise_transit_set(
            dt.datetime(2021, 1, 4, 14, 0, 0, tzinfo=CET), 49.6, 6.9), (
            8.477979163121788, 12.624062013177186, 16.77403260994837)
        for s, f in zip(res, ref):
            self.assertAlmostEqual(s, f, delta=0.25)

    def test_spa_pos_values(self):
        res, ref = m.radiation.spa_sun_position(
            dt.datetime(2010, 6, 21, 0, 6, 0, tzinfo=MST), 40, -105), (
            -26.54, 1.1)
        for s, f in zip(res, ref):
            self.assertAlmostEqual(s, f, delta=0.25)

    def test_fast_pos_values(self):
        res, ref = m.radiation.fast_sun_position(
            dt.datetime(2010, 6, 21, 0, 6, 0, tzinfo=MST), 40, -105), (
            -26.54, 1.1)
        for s, f in zip(res, ref):
            self.assertAlmostEqual(s, f, delta=0.25)

    def test_spa_rts_wrongargs(self):
        with self.assertRaises(ValueError):
            m.radiation.spa_rise_transit_set(
                dt.datetime(2021, 1, 4, 14, 0, 0, tzinfo=CET),
                49.6, 6.9, 252, pp=980, tk=-1)
        with self.assertRaises(ValueError):
            m.radiation.spa_rise_transit_set(
                dt.datetime(2021, 1, 4, 14, 0, 0, tzinfo=CET),
                49.6, 6.9, 252, pp=-1, tk=288)
        with self.assertRaises(ValueError):
            m.radiation.spa_rise_transit_set(
                dt.datetime(2021, 1, 4, 14, 0, 0, tzinfo=CET),
                49.6, 6.9, 252, pp=-1, tk=288)
        with self.assertRaises(ValueError):
            m.radiation.spa_rise_transit_set(
                dt.datetime(2021, 1, 4, 14, 0, 0, tzinfo=CET),
                100, 6.9, 252, pp=980, tk=288)
        with self.assertRaises(ValueError):
            m.radiation.spa_rise_transit_set(
                dt.datetime(2021, 1, 4, 14, 0, 0, tzinfo=CET),
                -100, 5.9, 252, pp=980, tk=288)
        with self.assertRaises(ValueError):
            m.radiation.spa_rise_transit_set(
                dt.datetime(2021, 1, 4, 14, 0, 0, tzinfo=CET),
                49.6, -181, 252, pp=-1, tk=288)
        with self.assertRaises(ValueError):
            m.radiation.spa_rise_transit_set(
                dt.datetime(2021, 1, 4, 14, 0, 0, tzinfo=CET),
                49.6, 361, 252, pp=-1, tk=288)

    def test_fast_rts_wrongargs(self):
        with self.assertRaises(ValueError):
            m.radiation.fast_rise_transit_set(
                dt.datetime(2110, 1, 1, 14, 0, 0, tzinfo=CET),
                49.6, 6.9)
        with self.assertRaises(ValueError):
            m.radiation.fast_rise_transit_set(
                dt.datetime(1890, 1, 1, 14, 0, 0, tzinfo=CET),
                49.6, 6.9)
        with self.assertRaises(ValueError):
            m.radiation.fast_rise_transit_set(
                dt.datetime(1890, 1, 1, 14, 0, 0, tzinfo=CET),
                100., 6.9)
        with self.assertRaises(ValueError):
            m.radiation.fast_rise_transit_set(
                dt.datetime(2021, 1, 4, 14, 0, 0, tzinfo=CET),
                100., 6.9)
        with self.assertRaises(ValueError):
            m.radiation.fast_rise_transit_set(
                dt.datetime(2021, 1, 4, 14, 0, 0, tzinfo=CET),
                -100., 6.9)
        with self.assertRaises(ValueError):
            m.radiation.fast_rise_transit_set(
                dt.datetime(2021, 1, 4, 14, 0, 0, tzinfo=CET),
                49.6, -181)
        with self.assertRaises(ValueError):
            m.radiation.fast_rise_transit_set(
                dt.datetime(2021, 1, 4, 14, 0, 0, tzinfo=CET),
                100, 361)

    def test_spa_pos_wrongargs(self):
        with self.assertRaises(ValueError):
            m.radiation.spa_sun_position(
                dt.datetime(2021, 1, 4, 14, 0, 0, tzinfo=CET),
                49.6, 6.9, 252, pp=980, tk=-1)
        with self.assertRaises(ValueError):
            m.radiation.spa_sun_position(
                dt.datetime(2021, 1, 4, 14, 0, 0, tzinfo=CET),
                49.6, 6.9, 252, pp=-1, tk=288)
        with self.assertRaises(ValueError):
            m.radiation.spa_sun_position(
                dt.datetime(2021, 1, 4, 14, 0, 0, tzinfo=CET),
                49.6, 6.9, 252, pp=-1, tk=288)
        with self.assertRaises(ValueError):
            m.radiation.spa_sun_position(
                dt.datetime(2021, 1, 4, 14, 0, 0, tzinfo=CET),
                100, 6.9, 252, pp=980, tk=288)
        with self.assertRaises(ValueError):
            m.radiation.spa_sun_position(
                dt.datetime(2021, 1, 4, 14, 0, 0, tzinfo=CET),
                -100, 5.9, 252, pp=980, tk=288)
        with self.assertRaises(ValueError):
            m.radiation.spa_sun_position(
                dt.datetime(2021, 1, 4, 14, 0, 0, tzinfo=CET),
                49.6, -181, 252, pp=-1, tk=288)
        with self.assertRaises(ValueError):
            m.radiation.spa_sun_position(
                dt.datetime(2021, 1, 4, 14, 0, 0, tzinfo=CET),
                49.6, 361, 252, pp=-1, tk=288)

    def test_fast_pos_wrongargs(self):
        with self.assertRaises(ValueError):
            m.radiation.fast_sun_position(
                dt.datetime(2110, 1, 1, 14, 0, 0, tzinfo=CET),
                49.6, 6.9)
        with self.assertRaises(ValueError):
            m.radiation.fast_sun_position(
                dt.datetime(1890, 1, 1, 14, 0, 0, tzinfo=CET),
                49.6, 6.9)
        with self.assertRaises(ValueError):
            m.radiation.fast_sun_position(
                dt.datetime(1890, 1, 1, 14, 0, 0, tzinfo=CET),
                100., 6.9)
        with self.assertRaises(ValueError):
            m.radiation.fast_sun_position(
                dt.datetime(2021, 1, 4, 14, 0, 0, tzinfo=CET),
                100., 6.9)
        with self.assertRaises(ValueError):
            m.radiation.fast_sun_position(
                dt.datetime(2021, 1, 4, 14, 0, 0, tzinfo=CET),
                -100., 6.9)
        with self.assertRaises(ValueError):
            m.radiation.fast_sun_position(
                dt.datetime(2021, 1, 4, 14, 0, 0, tzinfo=CET),
                49.6, -181)
        with self.assertRaises(ValueError):
            m.radiation.fast_sun_position(
                dt.datetime(2021, 1, 4, 14, 0, 0, tzinfo=CET),
                100, 361)
