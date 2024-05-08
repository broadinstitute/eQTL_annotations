import argparse
import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", dest="peak_dist_bed", type=str, required=True)
    parser.add_argument('-a', dest='assay_peak_file', type=str, required=True)
    parser.add_argument('-f', dest='finemap_result', type=str, required=True)
    parser.add_argument('-g', dest='group_name', type=str, required=True)
    args = parser.parse_args()

    library = args.assay_peak_file.split('/')[-1].split('.')[0]
    peaks = pd.read_table(args.peak_dist_bed, header=None,
                          names=f"chr pos _pos _chr {library}_peak_start {library}_peak_end {library}_peak_dist".split(), index_col=[0, 1])

    finemapped_df = pd.read_table(args.finemap_result)

    finemapped_df.dropna(inplace=True)
    df = pd.DataFrame({'chr':finemapped_df.chr, 'pos':finemapped_df.pos.astype(int)}).set_index(['chr', 'pos'])

    min_peak_pos_all = peaks.reset_index().groupby(['chr', 'pos'])[f'{library}_peak_dist'].min()
    merged_peak_idx = df.reset_index().merge(min_peak_pos_all, how='left', on=['chr', 'pos'])
    finemapped_df[f'{library}_peak_dist']  = merged_peak_idx[f'{library}_peak_dist'].values
    print('saving')
    finemapped_df.to_parquet(f'{args.group_name}_finemap_CHIP_ATAC_overlap.parquet')

if __name__ == '__main__':
    main()