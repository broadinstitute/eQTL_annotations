import argparse
import pandas as pd
import numpy as np
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", dest="finemapped_results", nargs='+', default=[],
                        help="Array of finemap results.")
    parser.add_argument("-n", dest="names", nargs='+', default=[])
    parser.add_argument("-v", dest="variant_annotations", help='parquet with all variant annotations', type=str, required=True)
    args = parser.parse_args()

    all_var_annotations = pd.read_parquet(args.variant_annotations)
    ABC_peak_genes_all_vars = all_var_annotations[['variant_id', 'phenotype_id', 'ABC_peak_dist']]
    non_gene_spec_peaks = all_var_annotations.drop(columns=['ABC_peak_dist', 'phenotype_id'])

    for finemapped_file, group_name in zip(args.finemapped_results, args.names):
        if finemapped_file.endswith('.tsv'):
            finemapped_df = pd.read_table(finemapped_file)
        else:
            finemapped_df = pd.read_parquet(finemapped_file)

        if 'Variant_id' in finemapped_df.columns:
            finemapped_df = finemapped_df.rename(columns={'Variant_id':'variant_id'})
        if 'variant_id' not in finemapped_df.columns:
            raise ValueError('column variant_id needs to be in data.')
        finemapped_df['variant_id'] = finemapped_df.variant_id.str.replace('_', ':')

        fm_annots = finemapped_df.merge(non_gene_spec_peaks, on='variant_id', how='left')
        merged = fm_annots.merge(ABC_peak_genes_all_vars, on=['variant_id', 'phenotype_id'], how='left')
        merged.drop_duplicates(inplace=True)
        merged.to_parquet(f'{group_name}_fm_variants_annotations.parquet')

if __name__ == '__main__':
    main()
