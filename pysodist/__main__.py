# -*- coding: utf-8 -*-
"""
@author: Joey Davis <jhdavis@mit.edu> jhdavislab.org
"""

'''Pysodist: a python distribution of the isodist mass spectra fitting routine'''

def main():
    import argparse, os
    parser = argparse.ArgumentParser(description=__doc__)
    import pysodist
    parser.add_argument('--version', action='version', version='pysodist '+pysodist.__version__)

    import pysodist.commands.parse_input
    import pysodist.commands.extract_spectra
    import pysodist.commands.run_isodist
    import pysodist.commands.plot_spectra
    import pysodist.commands.full_pipeline

    modules = [pysodist.commands.parse_input,
        pysodist.commands.extract_spectra,
        pysodist.commands.run_isodist,
        pysodist.commands.plot_spectra,
        pysodist.commands.full_pipeline,
        ]

    subparsers = parser.add_subparsers(title='Choose a command')
    subparsers.required = 'True'

    def get_str_name(module):
        return os.path.splitext(os.path.basename(module.__file__))[0]

    for module in modules:
        this_parser = subparsers.add_parser(get_str_name(module), description=module.__doc__)
        module.add_args(this_parser)
        this_parser.set_defaults(func=module.main)

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
