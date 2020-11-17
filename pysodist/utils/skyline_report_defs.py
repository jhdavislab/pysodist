# -*- coding: utf-8 -*-
"""
@author: Joey Davis <jhdavis@mit.edu> jhdavislab.org
"""

import pandas as pd
########################################################
#####SKYLINE REPORT FIELD DEFINITIONS###################
########################################################
EXPECTED_REPORT_VERSION = 'davis_peptide_quant_v2'
ISOTOPE_FIND_FIELD = 'Precursor Mz'
SAMPLE_FIND_FIELD = 'Min Start Time'
PROTEIN_PREFERRED_NAME_FIELD = 'Protein Preferred Name'
PEPTIDE_MOD_SEQ_FIELD = 'Peptide Modified Sequence Monoisotopic Masses'
PEPTIDE_CHARGE_FIELD = 'Precursor Charge'
DETECT_Q_VALUE_FIELD = 'Detection Q Value'
RT_START_FIELD = 'Min Start Time'
RT_END_FIELD = 'Max End Time'
MZ_FIELD = 'Precursor Mz'
START_POS_FIELD = 'Begin Pos'
END_POS_FIELD = 'End Pos'

AVERAGE_PPM_FIELD = 'Average Mass Error PPM'
LIB_DOT_FIELD = 'Library Dot Product'
ISO_DOT_FIELD = 'Isotope Dot Product'
RATIO_DOT_FIELD = 'Ratio Dot Product'
MS1_TOTAL_AREA_FIELD = 'Total Area MS1'
FRAG_TOTAL_AREA_FIELD = 'Total Area Fragment'
CALC_RATIO_MS1_FIELD = 'MS1_ratio'
CALC_RATIO_FRAG_FIELD = 'Frag_ratio'
CALC_RATIO_TOTAL_FIELD = 'Total_ratio'
ANY_RATIO_FIELD = '_ratio'
NORM_CALC_RATIO_MS1_FIELD = 'MS1_ratio_norm'
NORM_CALC_RATIO_FRAG_FIELD = 'Frag_ratio_norm'
NORM_CALC_RATIO_TOTAL_FIELD = 'Total_ratio_norm'

########################################################
#####MASS SPEC GENERAL DEFINITIONS######################
########################################################
AAWEIGHTS = {
    'A':71.03711,
    'R':156.10111,
    'N':114.04293,
    'D':115.02694,
    'C':103.00919,
    'E':129.04259,
    'Q':128.05858,
    'G':57.02146,
    'H':137.05891,
    'I':113.08406,
    'L':113.08406,
    'K':128.09496,
    'M':131.04049,
    'F':147.06841,
    'P':97.05276,
    'S':87.03203,
    'T':101.04768,
    'W':186.07931,
    'Y':163.06333,
    'V':99.06841
}
AANITROGENS = {
    'A':1,
    'R':4,
    'N':2,
    'D':1,
    'C':1,
    'c':1, # +57 mod cys (carbamidomethylated cys, iodoacetamide treatement - typically does not get labeled differently from normal cys)
    'b':1, # +58 mod cys (carboxymethylated cys, iodoacetic acid treatement - typically does not get labeled differently from normal cys)
    'E':1,
    'Q':2,
    'G':1,
    'H':3,
    'I':1,
    'L':1,
    'K':2,
    'M':1,
    'm':1, # oxidized methionine (typically does not get labeled differently from normal met)
    'F':1,
    'P':1,
    'S':1,
    'T':1,
    'W':2,
    'Y':1,
    'V':1
}

AACARBONS = {
    'A':3,
    'R':6,
    'N':4,
    'D':4,
    'C':3,
    'c':3, # +57 mod cys (carbamidomethylated cys, iodoacetamide treatement)
    'b':3, # +58 mod cys (carboxymethylated cys, iodoacetic acid treatement)
    'E':5,
    'Q':5,
    'G':2,
    'H':6,
    'I':6,
    'L':6,
    'K':6,
    'M':5,
    'm':5, # oxidized methionine (typically does not get labeled differently from normal met)
    'F':9,
    'P':5,
    'S':3,
    'T':4,
    'W':11,
    'Y':9,
    'V':5
}

ISODIST_MOD = {
    'M[+15.994915]' : 'm',
    'C[+57.021464]' : 'c',
    'C[+58.005479]' : 'b'
}

PMS_MOD = {value:key for key, value in ISODIST_MOD.items()}

N14MASS = 14.0030740048
N15MASS = 15.0001088982

C12MASS = 12.000000
C13MASS = 12.003355

reds = ['#fee5d9', '#fcbba1', '#fc9272', '#fb6a4a', '#de2d26', '#a50f15']
blues = ['#eff3ff', '#c6dbef', '#93cae1', '#6baed6', '#3182bd', '#08519c']
diverging = ['#d7191c', '#fdae61', '#ffffbf', '#a6d96a', '#1a9641']
