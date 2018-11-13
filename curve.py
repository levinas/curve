#! /usr/bin/env python

import numpy as np
import pandas as pd

# import os, sys, inspect
# current_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
# parent_dir = os.path.dirname(current_dir)
# sys.path.insert(0, parent_dir)

import uno_data as ud


def main():
    path = '/raid/fangfang/Benchmarks/Data/Pilot1/rescaled_combined_single_drug_growth'
    df0 = pd.read_table(path, engine='c', na_values=['na', '-', ''],
                        nrows=100,
                        dtype={'SOURCE': str, 'DRUG_ID': str,
                               'CELLNAME': str, 'CONCUNIT': str,
                               'LOG_CONCENTRATION': np.float32,
                               'EXPID': str, 'GROWTH': np.float32})



if __name__ == '__main__':
    main()
