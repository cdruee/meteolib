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
UTC = dt.timezone(dt.timedelta(hours=0))

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
        # reference values calculated by NREL SOLPOS Calculator
        # https://midcdmz.nrel.gov/solpos/solpos.html
        fun = m.radiation.spa_sun_position
        tuples = [
            # (mon,day, hr, lat,lon, ele, azi)
            # equator Jan
            (1, 1,  9, 0, 0, 39.8611, 120.5933),
            (1, 1, 12, 0, 0, 66.9985, 177.9096),
            # equator July
            (7, 1,  9, 0, 0, 39.8611, 59.3124),
            (7, 1, 12, 0, 0, 66.9985, 2.519),
            # 50N July
            (7, 1, 9, 50, 0, 45.3667, 109.7713),
            (7, 1, 12, 50, 0, 63.0830, 178.0475),
            # 50N 90E July
            (7, 1,  9, 50, 90, 46.5148, 248.3167),
            (7, 1, 12, 50, 90, 18.0823, 284.6427),
            # pole July
            (7, 1,  9, 89, 0, 23.7953, 133.7327),
            (7, 1, 12, 89, 0, 24.0934, 179.0427),
        ]
        for mon, day, hr, lat, lon, ele, azi in tuples:
            res = fun(
                dt.datetime(2010, mon, day, hr, 0, 0, tzinfo=UTC),
                lat, lon
            )
            ref = (ele, azi)
            for s, f in zip(res, ref):
                self.assertAlmostEqual(s, f, delta=0.5)

    def test_fast_pos_values(self):
        # reference values calculated by NREL SOLPOS Calculator
        # https://midcdmz.nrel.gov/solpos/solpos.html
        fun = m.radiation.fast_sun_position
        tuples = [
            # (mon,day, hr, lat,lon, ele, azi)
            # equator Jan
            (1, 1,  9, 0, 0, 39.8611, 120.5933),
            (1, 1, 12, 0, 0, 66.9985, 177.9096),
            # equator July
            (7, 1,  9, 0, 0, 39.8611, 59.3124),
            (7, 1, 12, 0, 0, 66.9985, 2.519),
            # 50N July
            (7, 1, 9, 50, 0, 45.3667, 109.7713),
            (7, 1, 12, 50, 0, 63.0830, 178.0475),
            # 50N 90E July
            (7, 1,  9, 50, 90, 46.5148, 248.3167),
            (7, 1, 12, 50, 90, 18.0823, 284.6427),
            # pole July
            (7, 1,  9, 89, 0, 23.7953, 133.7327),
            (7, 1, 12, 89, 0, 24.0934, 179.0427),
        ]
        for mon, day, hr, lat, lon, ele, azi in tuples:
            res = fun(
                dt.datetime(2010, mon, day, hr, 0, 0, tzinfo=UTC),
                lat, lon
            )
            ref = (ele, azi)
            for s, f in zip(res, ref):
                self.assertAlmostEqual(s, f, delta=0.5)

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


class Test_irradiance(unittest.TestCase):

    def test_clear_sky_direct_normal_values(self):
        # 11:47 UTC = 12:00 local in Jan
        # 12:01 UTC = 12:00 local in Jul
        # reference values: [ASH1997]_ table 18 page 29.32
        self.assertAlmostEqual(
            m.radiation.clear_sky_direct_normal(
                dt.datetime(2020, 1, 21, 11, 47, 0, tzinfo=UTC),
                40.0, 0), 926, delta=2)
        self.assertAlmostEqual(
            m.radiation.clear_sky_direct_normal(
                dt.datetime(2020, 1, 21, 5, 47, 0, tzinfo=UTC),
                40.0, 0), 0, delta=1)
        self.assertAlmostEqual(
            m.radiation.clear_sky_direct_normal(
                dt.datetime(2020, 6, 21, 12, 1, 0, tzinfo=UTC),
                40.0, 0), 879, delta=1)
        self.assertAlmostEqual(
            m.radiation.clear_sky_direct_normal(
                dt.datetime(2020, 6, 21, 9, 1, 0, tzinfo=UTC),
                40.0, 0), 829, delta=2)
        self.assertAlmostEqual(
            m.radiation.clear_sky_direct_normal(
                dt.datetime(2020, 6, 21, 6, 1, 0, tzinfo=UTC),
                40.0, 0), 488, delta=3)

        # reference values [ASH1997]_ table 21 page 29.35
        self.assertAlmostEqual(
            m.radiation.clear_sky_direct_normal(
                dt.datetime(2020, 1, 21, 11, 47, 0, tzinfo=UTC),
                64.0, 0), 316, delta=4)
        self.assertAlmostEqual(
            m.radiation.clear_sky_direct_normal(
                dt.datetime(2020, 1, 21, 5, 47, 0, tzinfo=UTC),
                64.0, 0), 0, delta=2)
        self.assertAlmostEqual(
            m.radiation.clear_sky_direct_normal(
                dt.datetime(2020, 6, 21, 12, 1, 0, tzinfo=UTC),
                64.0, 0), 831, delta=2)
        self.assertAlmostEqual(
            m.radiation.clear_sky_direct_normal(
                dt.datetime(2020, 6, 21, 9, 1, 0, tzinfo=UTC),
                64.0, 0), 791, delta=2)
        self.assertAlmostEqual(
            m.radiation.clear_sky_direct_normal(
                dt.datetime(2020, 6, 21, 6, 1, 0, tzinfo=UTC),
                64.0, 0), 613, delta=2)
