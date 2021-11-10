#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Common utilities used by loads of multiple formats.
"""


def common_start(string_list):
    """Returns the longest common substring from the substings.
    
    Made recursive from https://stackoverflow.com/questions/18715688/find-common-substring-between-two-strings
    """
    def _cs(string_a, string_b):
        def _iter():
            for char_a, char_b in zip(string_a, string_b):
                if char_a == char_b:
                    yield char_a
                else:
                    return
        return ''.join(_iter())
    if len(string_list) == 1:
        return string_list[0]
    else:
        # copy in a 2/3 safe way
        sl = string_list[:]
        while len(sl) > 1:
            sl[-2] = _cs(sl[-2], sl[-1])
            sl.pop()
        return sl[0]
