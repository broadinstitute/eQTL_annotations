import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", dest="peak_dist_bed", type=str, required=True)
    parser.add_argument("-a", dest="abc_peak_dist_bed", type=str, required=True)
    parser.add_argument('-g', dest='gtex_vep', help='gtex vep overlap file in tsv.gz form', type=str, required=True)
    args = parser.parse_args()

    print("Reading in all variant peak overlap data")
    all_peak_overlaps = pd.read_table(args.peak_dist_bed, header=None,
                                  names=f"chr pos _pos variant_id library _chr peak_start peak_end peak_dist".split(), index_col=[0, 1])


    abc_only_peak_dists = pd.read_table(args.abc_peak_dist_bed, header=None,
                                names=f"chr pos _pos variant_id _chr peak_start peak_end phenotype_id peak_dist".split(), index_col=[0, 1])
    all_abc_peaks = abc_only_peak_dists[abc_only_peak_dists.peak_dist != -1]
    all_abc_peaks[f'ABC_peak_dist'] = all_abc_peaks.peak_dist

    grp = all_peak_overlaps.groupby('library')
    dfs = [grp.get_group(g) for g in all_peak_overlaps.groupby('library').groups]

    peak_dfs = []
    for df in dfs:
        lib_name = df.library.unique()[0]
        df[f'{lib_name}_peak_dist'] = df.peak_dist
        peak_dfs.append(df[['variant_id', f'{lib_name}_peak_dist']].reset_index().set_index(['variant_id', 'chr', 'pos']))
    all_var_peak_df = pd.concat(peak_dfs, axis=1)

    print("Reading in gtex vep.")
    gtex_vep = pd.read_table(args.gtex_vep, index_col=0)
    gtex_vep.index = gtex_vep.index.str.replace('_b38', '').str.replace('_', ':')

    print("Annotating all variants with gtex data.")
    gtex_overlap = gtex_vep[gtex_vep.index.isin(all_var_peak_df.reset_index().set_index('variant_id').index)]
    gtex_overlap = gtex_overlap.reset_index().rename(columns={'SNP':'variant_id'})
    # merge cur data with gtex vep data
    all_var_peak_gtex_df = all_var_peak_df.merge(gtex_overlap, on='variant_id', how='outer')

    # merge with abc data
    all_variant_annotations = all_var_peak_gtex_df.merge(all_abc_peaks[['variant_id', 'ABC_peak_dist', 'phenotype_id']],
                            on='variant_id', how='outer')

    all_variant_annotations.to_parquet(f'all_variant_CHIP_ATAC_GTEx_overlap.parquet')

if __name__ == '__main__':
    main()