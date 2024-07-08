import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os.path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", dest="finemapped_annotations", nargs='+', default=[],
                        help="Array of finemap results with annotations already added, across all group names.")
    parser.add_argument("-g", dest="group_names", nargs='+', default=[])
    parser.add_argument("-v", dest="variant_annotations", help='parquet with all variant annotations', type=str, required=True)
    args = parser.parse_args()

    annotation_map = {'ATAC_peak_dist':'ATAC peak dist', 'CTCF_peak_dist':'CTCF peak dist',
                        'enhancer_d':'Enhancer', 'promoter_d':'Promoter', 'CTCF_binding_site_d':'CTCF binding site', 'TF_binding_site_d':'TF binding site',
                        '3_prime_UTR_variant_d':"3' UTR", '5_prime_UTR_variant_d':"5' UTR", 'intron_variant_d':"Intron",
                        'missense_variant_d':"Nonsynonymous", 'synonymous_variant_d':"Synonymous",
                        'open_chromatin_region_d':'Open chromatin', 'promoter_flanking_region_d':'Promoter Flanking',
                        'frameshift_variant_d':'Frameshift Variant', 'stop_gained_d':'Stop Gained',
                        'non_coding_transcript_exon_variant_d':'Non-coding transcript exon variant'}
    non_annotations = ['phenotype_id', 'variant_id', 'pip', 'af', 'cs_id', 'start_distance', 'ma_samples', 'ma_count', 'pval_nominal', 'slope', 'slope_se','bins']

    print('Reading in all variant annotations & adding in peak info.')
    all_variant_annots = pd.read_parquet(args.variant_annotations)
    all_var_500 = all_variant_annots.loc[:, all_variant_annots.columns.str.contains('peak_dist')] < 500
    all_var_in_a_peak = all_variant_annots.loc[:, all_variant_annots.columns.str.contains('peak_dist')] == 0
    all_var_in_a_peak.columns = all_var_in_a_peak.columns.str.strip('peak_dist') + '_in_a_peak'
    all_var_500.columns = all_var_500.columns.str.strip('peak_dist') + '_500bp_from_peak'
    all_variant_annots = pd.concat([all_variant_annots, all_var_in_a_peak, all_var_500], axis=1)
    # annotate singular peak variants in all variants
    only_in_one_peak = {}
    for day in "D0 D2 D4 D7".split():
        if all_variant_annots.columns.str.contains(f'ATAC_{day}').any():
            only_in_atac = (all_variant_annots.loc[:, all_variant_annots.columns.str.contains('ATAC_D')
                                                & all_variant_annots.columns.str.endswith('in_a_peak')].sum(axis=1) == 1) & (all_variant_annots[f'ATAC_{day}_in_a_peak'])
            only_in_one_peak[f'only_ATAC_{day}'] = only_in_atac
        if all_variant_annots.columns.str.contains(f'CTCF_{day}').any():
            only_in_ctcf = (all_variant_annots.loc[:, all_variant_annots.columns.str.contains('CTCF_D')
                                                & all_variant_annots.columns.str.endswith('in_a_peak')].sum(axis=1) == 1) & (all_variant_annots[f'CTCF_{day}_in_a_peak'])
            only_in_one_peak[f'only_CTCF_{day}'] = only_in_ctcf
    only_one_peak_df = pd.concat(only_in_one_peak, axis=1)
    all_variant_annots = pd.concat([all_variant_annots, only_one_peak_df], axis=1)

    finemapped_dfs = args.finemapped_annotations
    group_names = args.group_names

    print('Reading in & adding additional info for each fm file.')
    fm_dict = {}
    for group_name in group_names:
        fm_df_str = [x for x in finemapped_dfs if os.path.basename(x) == f'{group_name}_fm_variants_annotations.parquet']
        fm_annot_df = pd.read_parquet(fm_df_str)
        # dont include splice or frameshift
        fm_annot_df = fm_annot_df.loc[:, ~fm_annot_df.columns.str.contains('splice|frameshift|CTCF_binding|open_chromatin')]
        # annotate peaks categorically as well
        peaks_500 = fm_annot_df.loc[:, fm_annot_df.columns.str.contains('peak_dist')] < 500
        in_a_peak = fm_annot_df.loc[:, fm_annot_df.columns.str.contains('peak_dist')] == 0
        in_a_peak.columns = in_a_peak.columns.str.strip('peak_dist') + '_in_a_peak'
        peaks_500.columns = peaks_500.columns.str.strip('peak_dist') + '_500bp_from_peak'
        # add that info
        # drop 500 bp info for now
        fm_annot_df = pd.concat([fm_annot_df, in_a_peak], axis=1)
        # drop peak dist info - dont want it
        fm_annot_df = fm_annot_df.loc[:, ~fm_annot_df.columns.str.contains('peak_dist')]

        # annotate day specific atac/ctcf info
        if fm_annot_df.columns.str.contains('ATAC_D').any():
            for day in "D0 D2 D4 D7".split():
                if not fm_annot_df.columns.str.contains(f'ATAC_{day}').any():
                    continue
                only_in_atac = ((fm_annot_df.loc[:, fm_annot_df.columns.str.startswith('ATAC_D') &
                                                 fm_annot_df.columns.str.endswith('in_a_peak')].sum(axis=1) == 1) & (fm_annot_df[f'ATAC_{day}_in_a_peak'])).rename(f'only_ATAC_{day}')
                fm_annot_df = pd.concat([fm_annot_df, only_in_atac],
                    axis=1)

        if fm_annot_df.columns.str.contains('CTCF_D').any():
            for day in "D2 D4 D7".split():
                if not fm_annot_df.columns.str.contains(f'CTCF_{day}').any():
                    continue
                only_in_ctcf = ((fm_annot_df.loc[:, fm_annot_df.columns.str.contains('CTCF_D') &
                                                 fm_annot_df.columns.str.endswith('in_a_peak')].sum(axis=1) == 1) & (fm_annot_df[f'CTCF_{day}_in_a_peak'])).rename(f'only_CTCF_{day}')
                fm_annot_df = pd.concat([fm_annot_df, only_in_ctcf],
                    axis=1)

        fm_dict[group_name] = fm_annot_df

    # all annotations are included
    annotations = fm_annot_df.loc[:, ~fm_annot_df.columns.isin(non_annotations)].columns

    print('Getting means in fm data compared to background (all variants).')
    # now annotate mean of each day compared to backgrounds
    mean_arr = pd.DataFrame(0.0, index=annotations, columns=np.hstack((group_names, 'background_snps')))
    for group_name in group_names:
        fm_annot_df = fm_dict[group_name]
        for i, annotation in enumerate(annotations):
                mean_arr.at[annotation, group_name] = fm_annot_df[annotation].mean()
                mean_arr.at[annotation, 'background_snps'] = all_variant_annots[annotation].mean()

    mean_arr.to_csv('raw_mean_by_group_gtex_plot.tsv', sep='\t', header=True, index=True)
    log2FE_arr = mean_arr.iloc[:,:-1] / mean_arr.iloc[:,-1].values[:,None]
    log2FE_arr = np.log2(log2FE_arr)

    print('Plotting.')
    # plot enrichment and prop. variants plot
    fig, [ax, ax1] = plt.subplots(1,2,figsize=(16, len(annotations)), sharey=True)
    colors = 'k #FFD39B #BCEE68 #556B2E #FF6A6A #CD5555 #8B393A'.split()
    group_shifts = np.linspace(0.4, -0.4, len(group_names))
    annotation_labels = [annotation_map[x] if x in annotation_map.keys() else x.replace('_', ' ') for x in annotations]

    # add in the background gray bars (this is such a kerchoo way to do this oops)
    gray_bars = np.ones(len(annotations))+3
    gray_bars[1::2] = 0
    ax.barh(range(len(annotations)), gray_bars, align='center', height=1, alpha=0.3, color='gray')
    ax.barh(range(len(annotations)), -gray_bars/4, align='center', height=1, alpha=0.3, color='gray')

    # plot the fold enrichment
    for i, group_name in enumerate(group_names):
        ax.scatter(log2FE_arr.loc[:,group_name], np.arange(len(annotations))+group_shifts[i], label=group_names[i], color=colors[i])
    ax.set_yticks(range(len(annotations)))
    ax.set_yticklabels(annotation_labels)
    ax.set_xlabel('$log_{2}$(Fold Enrichment)        ', fontsize = 30)
    for loc in np.arange(-1,5):
        ax.axvline(x=loc, c='k', ls='--', lw=.35, zorder=0, alpha=0.5)
    ax.axvline(0, c='k',  ls='--', lw=.6, zorder=0)
    # plot the proportion of variants
    bar_height = 0.1
    for i, group_name in enumerate(group_names):
        legend_label = f'{group_names[i]}, n={fm_dict[group_names[i]].shape[0]}'
        ax1.barh(np.arange(len(annotations))+group_shifts[i], mean_arr.loc[:,group_name], label=legend_label, color=colors[i], height = bar_height)
    ax1.legend(bbox_to_anchor=(1,1), loc="upper left", fontsize=20)
    ax1.set_xlabel('Prop. of Variants', fontsize = 30)
    ax1.set_xlim(0, .2)
    for loc in np.arange(0,1,.1):
        ax1.axvline(x=loc, c='k', ls='--', lw=.35, zorder=0, alpha=0.5)

    fig.tight_layout()
    fig.savefig('gtex_annot_enrich.png', dpi=300)
    print('Done.')


if __name__ == '__main__':
    main()