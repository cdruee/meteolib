#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 18:39:22 2019

@author: druee
"""

import unittest

import meteolib as m


class Test__check(unittest.TestCase):
    """
    def _check(name,par,kind,lt=None,le=None,ne=None,ge=None,gt=None)
    """

    def test_1_float(self):
        result = m._utils._check('x', 1., 'float')
        self.assertEqual(result, 1.)

    def test_1_int(self):
        result = m._utils._check('x', 1, 'int')
        self.assertEqual(result, 1)

    def test_wrong_kind(self):
        self.assertRaises(ValueError, m._utils._check, 'x', 1, 'blah')

    def test_range_bad(self):
        with self.assertRaises(ValueError):
            m._utils._check('x', 1., 'float', lt=1.)
            m._utils._check('x', 1., 'float', le=0.)
            m._utils._check('x', 1., 'float', ne=1.)
            m._utils._check('x', 1., 'float', ge=1.)
            m._utils._check('x', 1., 'float', gt=2.)

    def test_range_good(self):
        self.assertEqual(m._utils._check('x', 1., 'float', lt=2.), 1.)
        self.assertEqual(m._utils._check('x', 1., 'float', le=1.), 1.)
        self.assertEqual(m._utils._check('x', 1., 'float', ge=1.), 1.)
        self.assertEqual(m._utils._check('x', 1., 'float', gt=0.), 1.)


class Test__only(unittest.TestCase):
    '''
    def _only(pars,sel):
    '''

    def test_types(self):
        self.assertRaises(TypeError, m._utils._only, {
                          'a': None, 'b': None, 'c': None, 'd': None}, 1)
        self.assertRaises(TypeError, m._utils._only, 1, ['b'])

    def test_true_list(self):
        self.assertRaises(TypeError, m._utils._only(
            {'a': None, 'b': 1, 'c': None, 'd': None}, ['b']))
        self.assertRaises(TypeError, m._utils._only(
            {'a': None, 'b': 1, 'c': 1, 'd': None}, ['b', 'c']))

    def test_true_str(self):
        self.assertRaises(TypeError, m._utils._only(
            {'a': None, 'b': 1, 'c': None, 'd': None}, 'b'))

    def test_false_list(self):
        self.assertRaises(TypeError, m._utils._only(
            {'a': None, 'b': 1, 'c': None, 'd': None}, ['c']))
        self.assertRaises(TypeError, m._utils._only(
            {'a': None, 'b': 1, 'c': 1, 'd': None}, ['a', 'd']))
        self.assertRaises(TypeError, m._utils._only(
            {'a': None, 'b': 1, 'c': 1, 'd': None}, ['a', 'b']))
        self.assertRaises(TypeError, m._utils._only(
            {'a': None, 'b': 1, 'c': 1, 'd': None}, ['b']))

    def test_false_str(self):
        self.assertRaises(TypeError, m._utils._only(
            {'a': None, 'b': 1, 'c': None, 'd': None}, 'd'))


if __name__ == '__main__':
    unittest.main()
