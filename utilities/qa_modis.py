"""
This module makes a qc mask based on specific user input. see : MODIS Collection 6 (C6) LAI/FPAR Product User’s Guide
Table 5: Values of FparLAI_QC (8−bit)
"""

import numpy as np

MODLAND_QC = 0b00000001
SENSOR = 0b00000010
DETECTOR = 0b00000100
CLOUD_STATE = 0b00011000
SCF_QC = 0b11100000


def _get_modland_bits(value):
    # 0: good quality
    # 1: other quality, backup algorithm or fill values
    return (MODLAND_QC & value) >> 0


def _get_sensor_bits(value):
    # 0: Terra
    # 1: Aqua
    return (SENSOR & value) >> 1


def _get_detector_bits(value):
    # 0: Detectors apparently fine for up to 50% of channels 1, 2
    # 1: Dead detectors caused >50% adjacent detector retrieval
    return (DETECTOR & value) >> 2


def _get_cloud_bits(value):
    # 0: 0 Significant clouds NOT present (clear)
    # 1: 1 Significant clouds WERE present
    # 2: 2 Mixed cloud present in pixel
    # 3: 3 Cloud state not defined, assumed clear
    return (CLOUD_STATE & value) >> 3


def _get_scf_qc(value):  # return the value of digit 5 - end
    # 0: 0 Main (RT) method used, best result possible (no saturation)
    # 1: 1 Main (RT) method used with saturation. Good, very usable
    # 2: 2 Main (RT) method failed due to bad geometry, empirical algorithm used
    # 3: 3 Main (RT) method failed due to problems other than geometry, empirical algorithm used
    # 4: 4 Pixel not produced at all, value couldn’t be retrieved (possible reasons: bad L1B data, unusable
    # MOD09GA data)
    return (SCF_QC & value) >> 5


def make_LAI_FPAR_QC_mask(qc_array, parameter_name='SCF_QC', value=0):
    """
    sets 0 for all values not equal to value and sets 1 for values equal to value
    :param qc_array: the 1_original QC array
    :param parameter_name: which parameter we want to work with
    :param value: which values we should keep
    :return: array of 0 and 1 where 0 means the element was not equal to value and 1 means element was equal to value
    """
    filter_mask = np.array([])
    if parameter_name == 'MODLAND_QC':
        filter_mask = _get_modland_bits(qc_array)
    elif parameter_name == 'SENSOR':
        filter_mask = _get_sensor_bits(qc_array)
    elif parameter_name == 'DETECTOR':
        filter_mask = _get_detector_bits(qc_array)
    elif parameter_name == 'CLOUD_STATE':
        filter_mask = _get_cloud_bits(qc_array)
    elif parameter_name == 'SCF_QC':
        filter_mask = _get_scf_qc(qc_array)

    filter_mask[filter_mask != value] = 1000
    filter_mask[filter_mask == value] = 1
    filter_mask[filter_mask == 1000] = 0
    return filter_mask


def filter_non_zero(array):
    """
    sets all the nonzero values to 0 and all zeros to 1. This function is useful when make_LAI_FPAR_QC_mask is
    used multiple times and the result are summed up as a mask and some elements of the filter will not be more than
    1 so this function can set them to 1
    f1 = make_LAI_FPAR_QC_mask(qc_array, parameter_name='SCF_QC', value=0)
    f2 = make_LAI_FPAR_QC_mask(qc_array, parameter_name='SCF_QC', value=2)
    f3 = make_LAI_FPAR_QC_mask(qc_array, parameter_name='CLOUD_STATE', value=0)
    f = f1 + f2 + f3
    f = filter_non_zero(f)
    :param array:
    :return:
    """
    filter_mask = array
    filter_mask[filter_mask != 0] = 1000
    filter_mask[filter_mask == 0] = 0
    filter_mask[filter_mask == 1000] = 1
    return filter_mask
