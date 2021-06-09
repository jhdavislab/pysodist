# -*- coding: utf-8 -*-
"""
@author: Joey Davis <jhdavis@mit.edu> jhdavislab.org
@version: 0.0.5
"""

from datetime import datetime as dt
import sys


def log(msg, outfile=None):
    msg = '{} --> {}'.format(dt.now().strftime('%Y-%m-%d %H:%M:%S'), msg)
    print(msg)
    sys.stdout.flush()
    if outfile is not None:
        try:
            with open(outfile, 'a') as f:
                f.write(msg + '\n')
        except Exception as e:
            log(e)
