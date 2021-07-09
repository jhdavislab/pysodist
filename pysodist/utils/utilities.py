# -*- coding: utf-8 -*-
"""
@author: Joey Davis <jhdavis@mit.edu> jhdavislab.org
@version: 0.0.5
"""

from datetime import datetime as dt
import sys


def clean_path(path):
    path = path.replace('\\', '/')
    if path[-1] != '/':
        path += '/'
    return path


def clean_config(config_data):
    config_data.loc['output_directory']['VALUE'] = clean_path(config_data.loc['output_directory']['VALUE'])
    config_data.loc['isodist_exe']['VALUE'] = clean_path(config_data.loc['isodist_exe']['VALUE'])[:-1]
    config_data.loc['atom_file']['VALUE'] = clean_path(config_data.loc['atom_file']['VALUE'])[:-1]
    config_data.loc['res_file']['VALUE'] = clean_path(config_data.loc['res_file']['VALUE'])[:-1]
    config_data.loc['guide_file']['VALUE'] = clean_path(config_data.loc['guide_file']['VALUE'])[:-1]
    config_data.loc['mzml_file']['VALUE'] = clean_path(config_data.loc['mzml_file']['VALUE'])[:-1]
    return config_data


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
