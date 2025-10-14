#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@filename: genome_variation.py
@author: Bangjun Gong
@contact: bangjungong@foxmail.com
@license: GPL v3
@created: 2025-10-13
@updated: 2025-10-13
@version: 1.0.0
@description:
"""

import numpy as np


def calculate_boxplot_parameters(data: list, outliers_type:str = "both"):
    '''
    Calculates box plot parameters for given data and outlier type.
    :param data:
    :param outliers_type: both| upper| lower| none
    :return:
    '''

    # remove nan
    data_arr = np.array(data)
    data_arr = data_arr[~np.isnan(data_arr)]

    # Sort the input data in ascending order (for quantile computation).
    sorted_data = np.sort(data_arr)

    # eigenvalues
    median = np.median(sorted_data)
    q1 = np.percentile(sorted_data, 25)
    q3 = np.percentile(sorted_data, 75)
    iqr = q3 - q1

    mean = np.mean(sorted_data)
    std = np.std(sorted_data)
    max_value = np.max(sorted_data)
    min_value = np.min(sorted_data)

    # theoretical upper/lower bound in boxplot
    theoretical_lower_bound = q1 - 1.5 * iqr
    theoretical_upper_bound = q3 + 1.5 * iqr

    # real upper/lower bound in boxplot
    whisker_min = next((x for x in sorted_data if x >= theoretical_lower_bound), None)
    whisker_max = next((x for x in reversed(sorted_data) if x <= theoretical_upper_bound), None)

    real_lower_bound = whisker_min if whisker_min is not None else theoretical_lower_bound
    real_upper_bound = whisker_max if whisker_max is not None else theoretical_upper_bound

    # get outliers
    opt = {
        "both": lambda: [float(x) for x in data if x < theoretical_lower_bound or x > theoretical_upper_bound],
        "upper": lambda: [float(x) for x in data if x > theoretical_upper_bound],
        "lower": lambda: [float(x) for x in data if x < theoretical_lower_bound],
        "none": lambda: []
    }
    outliers = opt.get(outliers_type)()

    return _BoxPlotParameters(
        q1=float(q1),
        median=float(median),
        q3=float(q3),
        iqr=float(iqr),
        real_lower_bound=float(real_lower_bound),
        real_upper_bound=float(real_upper_bound),
        outliers=outliers,
        theoretical_lower_bound=float(theoretical_lower_bound),
        theoretical_upper_bound=float(theoretical_upper_bound),
        max_value=float(max_value),
        min_value=float(min_value),
        mean=float(mean),
        std=float(std),
    )


class _BoxPlotParameters(object):
    def __init__(
            self, q1, median, q3, iqr,
            real_lower_bound,
            real_upper_bound,
            outliers,
            theoretical_lower_bound,
            theoretical_upper_bound,
            min_value,
            max_value,
            mean,
            std,
    ):
        '''
        Store the eigenvalues of the pairwise similarity matrix from the sequence alignment.
        :param q1: 25th percentile (Q1),
        :param median: median,
        :param q3: 75th percentile (Q3)
        :param iqr: q3-q1
        :param real_lower_bound: real lower bound in boxplot
        :param real_upper_bound: real upper bound in boxplot
        :param outliers: outliers
        :param theoretical_lower_bound: theoretical lower bound in boxplot
        :param theoretical_upper_bound: theoretical upper bound in boxplot
        :param min_value: minimal value
        :param max_value: maximum value
        :param mean: mean value
        :param std: standard deviation
        '''
        self.q1 = q1
        self.median = median
        self.q3 = q3
        self.iqr = iqr
        self.real_lower_bound = real_lower_bound
        self.real_upper_bound = real_upper_bound
        self.outliers = outliers
        self.theoretical_lower_bound = theoretical_lower_bound
        self.theoretical_upper_bound = theoretical_upper_bound
        self.min_value = min_value
        self.max_value = max_value
        self.mean = mean
        self.std = std


if __name__ == "__main__":
    test_data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 30]
    result = calculate_boxplot_parameters(test_data)
    print("boxplot statics data:", result.__dict__)