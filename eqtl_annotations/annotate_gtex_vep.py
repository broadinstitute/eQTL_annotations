import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", dest="finemapped_results", nargs='+', default=[],
                        help="Array of finemap results with ATAC and CHIP peaks already added, across all group names.")
    parser.add_argument("-g", dest="group_names", nargs='+', default=[])
    parser.add_argument('-g', dest='gtex_vep', help='gtex vep overlap file in tsv.gz form', type=str, required=True)
    args = parser.parse_args()

    gtex_vep = pd.read_table(args.gtex_vep, index_col=0)
    gtex_vep.index = gtex_vep.index.str.replace('_b38', '').str.replace('_', ':')

    for finemapped_file, group_name in zip(args.finemapped_results, args.group_names):
        finemapped_df = pd.read_table(finemapped_file)
        gtex_overlap = gtex_vep[gtex_vep.index.isin(finemapped_df.set_index('variant_id').index)]
        gtex_overlap = gtex_overlap.reset_index().rename(columns={'SNP':'variant_id'})
        # merge cur data with gtex vep data
        finemapped_df = finemapped_df.merge(gtex_overlap, on='variant_id', how='outer')
        finemapped_df.to_csv(f'{group_name}_finemap_CHIP_ATAC_GTEx_overlap.tsv', sep='\t', header=True, index=False)

if __name__ == '__main__':
    main()