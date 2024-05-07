import argparse
import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", dest="peak_dist_bed", type=str, required=True)
    parser.add_argument('-a', dest='assay_peak_file', type=str, require=True)
    parser.add_argument('-f', dest='finemap_file', type=str, required=True)
    args = parser.parse_args()

    library = args.assay_peak_file.split('/')[-1].split('.')[0]
    peaks = pd.read_table(args.peak_dist_bed, header=None,
                          names=f"chr pos _pos _chr {library}_peak_start {library}_peak_end {library}_peak_dist".split(), index_col=[0, 1])

    # TODO: combine with finemap data