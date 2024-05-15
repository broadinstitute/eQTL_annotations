import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.pyplot import cm

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", dest="finemapped_annotations", nargs='+', default=[],
                        help="Array of finemap results with annotations already added, across all group names.")
    parser.add_argument("-g", dest="group_names", nargs='+', default=[])
    parser.add_argument("-v", dest="variant_annotations", help='parquet with all variant annotations', type=str, required=True)
    args = parser.parse_args()

    finemapped_dfs = args.finemapped_annotations # ['ips_D0_fm_variants_annotations.parquet', 'hep_D2_fm_variants_annotations.parquet']
    group_names = args.group_names #['ips_D0', 'hep_D2']
    annotation_map = {'ATAC_peak_dist':'ATAC peak dist', 'CTCF_peak_dist':'CTCF peak dist',
                    'enhancer_d':'Enhancer', 'promoter_d':'Promoter', 'CTCF_binding_site_d':'CTCF binding site', 'TF_binding_site_d':'TF binding site',
                    '3_prime_UTR_variant_d':"3' UTR", '5_prime_UTR_variant_d':"5' UTR", 'intron_variant_d':"Intron",
                    'missense_variant_d':"Nonsynonymous", 'synonymous_variant_d':"Synonymous",
                    'open_chromatin_region_d':'Open chromatin', 'promoter_flanking_region_d':'Promoter Flanking',
                    'frameshift_variant_d':'Frameshift Variant', 'stop_gained_d':'Stop Gained',
                    'non_coding_transcript_exon_variant_d':'Non-coding transcript exon variant'}
    non_annotations = ['phenotype_id', 'variant_id', 'pip', 'af', 'cs_id', 'start_distance', 'ma_samples', 'ma_count', 'pval_nominal', 'slope', 'slope_se','bins']


    fm_dict = {}
    for fm_df_str, group_name in zip(finemapped_dfs, group_names):
        fm_df = pd.read_parquet(fm_df_str)
        fm_dict[group_name] = fm_df

    bins = [0, 0.01, 0.1, 0.5, 0.9, 1]
    labels = ['PIP<0.01' ,'0.01<PIP<0.1','0.1<PIP<0.5', '0.5<PIP<0.9', '0.9<PIP']

    for group_name in group_names:
        fm_annot_df = fm_dict[group_name].copy()

        fm_annot_df['bins'] = pd.cut(fm_annot_df.pip, bins=bins, labels=labels, right=True, include_lowest=True)
        fm_annot_df = fm_annot_df.loc[:, ~fm_annot_df.columns.str.contains('splice')]
        annotations = fm_annot_df.loc[:, ~fm_annot_df.columns.isin(non_annotations)].columns
        mean_arr = pd.DataFrame(0.0, index=annotations, columns=labels)

        for i, conseq in enumerate(annotations):
            for j, label in enumerate(labels):
                annot_bin_df = fm_annot_df.query('bins==@label').groupby(['phenotype_id','variant_id'], as_index=True)[annotations].agg('any')
                mean_arr.at[conseq, label] = annot_bin_df[conseq].mean() # ignores NaN by default

        norm = mean_arr['PIP<0.01'].values
        FE = (mean_arr.T / norm).T
        annotations = FE.index.values

        FE_melt_new = FE.reset_index().melt(id_vars='index')

        fig, ax = plt.subplots(figsize=(15,4))

        palette = ['#929591', '#FFD700', '#FFA500', '#F97306', '#FE420F']
        sns.set_palette(palette)
        sns.barplot(FE_melt_new, x='index', y='value', hue='variable')
        annotation_labels = [annotation_map[x] if x in annotation_map.keys() else x.replace('_', ' ') for x in annotations]
        ax.set_xticklabels(annotation_labels, rotation = 45, ha="right");
        ax.set_ylabel('Fold Enrichment')
        ax.set_xlabel('')
        ax.legend(title='')
        ax.set_title(group_name.replace('_', ' '))
        fig.tight_layout()
        fig.savefig(f'{group_name}_annotation_by_pip.png', dpi=300)


if __name__ == '__main__':
    main()
