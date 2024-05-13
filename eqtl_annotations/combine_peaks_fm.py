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
                          names=f"chr pos _pos variant_id _chr {library}_peak_start {library}_peak_end {library}_peak_dist".split(), index_col=[0, 1])

    finemapped_df = pd.read_table(args.finemap_result)

    #df = pd.DataFrame({'chr':finemapped_df.chr, 'pos':finemapped_df.pos.astype(int)}).set_index(['chr', 'pos'])

    min_peak_pos_all = peaks.reset_index().groupby(['chr', 'pos', 'variant_id'])[f'{library}_peak_dist'].min()
    finemapped_df = finemapped_df.merge(min_peak_pos_all, how='left', on=['chr', 'pos', 'variant_id'])
    # add this peak group distance to the finemapped results
    print('saving')
    finemapped_df.to_csv('finemapped_results.tsv', sep='\t', header=True, index=False)

    # save off summary info for all background peaks
    peaks[f'{library}_less_than_500'] = peaks[f'{library}_peak_dist'] < 500
    peaks[f'{library}_in_a_peak'] = peaks[f'{library}_peak_dist'] == 0
    peaks.reset_index()[['chr', 'pos', 'variant_id', f'{library}_peak_dist', f'{library}_less_than_500', f'{library}_in_a_peak']].to_csv(f'{library}_peak_stats.tsv',
                                                                                                                          sep='\t', header=True, index=False)

if __name__ == '__main__':
    main()