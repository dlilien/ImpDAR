#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
An executable to start the picker.
"""

import argparse

from impdar import pick, load


def main():
    parser = _get_args()
    args = parser.parse_args()

    pick(


def _get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('fn', type=str, help='The file to pick. At least for now you can only do one at a time.')
    return parser


if __name__ == '__main__':
    main()
