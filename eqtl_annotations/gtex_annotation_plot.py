import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", dest="finemapped_annotations", nargs='+', default=[],
                        help="Array of finemap results with annotations already added, across all group names.")
    parser.add_argument("-g", dest="group_names", nargs='+', default=[])
    parser.add_argument("-v", dest="variant_annotations", help='parquet with all variant annotations', type=str, required=True)
    args = parser.parse_args()

    annotations = ['enhancer_d', 'promoter_d', 'CTCF_binding_site_d', 'TF_binding_site_d', '3_prime_UTR_variant_d',
               '5_prime_UTR_variant_d', 'intron_variant_d', 'missense_variant_d', 'synonymous_variant_d']

    labels = ['Enhancer', 'Promoter', 'CTCF binding site', 'TF binding site', "3' UTR", "5' UTR", "Intron", "Nonsynonymous", "Synonymous"]

    all_variant_annots = pd.read_parquet(args.variant_annotations)

    finemapped_dfs = args.finemapped_annotations
    group_names = args.group_names
    mean_arr = pd.DataFrame(0.0, index=annotations, columns=np.hstack((group_names, 'background_snps')))

    fm_dict = {}
    for fm_df_str, group_name in zip(finemapped_dfs, group_names):
        fm_df = pd.read_parquet(fm_df_str)
        fm_dict[group_name] = fm_df
        for i, annotation in enumerate(annotations):
            mean_arr.at[annotation,group_name] = fm_df[annotation].mean()
            mean_arr.at[annotation, 'background_snps'] = all_variant_annots[annotation].mean()

    log2FE_arr = mean_arr.iloc[:,:-1] / mean_arr.iloc[:,-1].values[:,None]
    log2FE_arr = np.log2(log2FE_arr)

    fig, [ax, ax1] = plt.subplots(1,2,figsize=(16,10), sharey=True)
    colors = 'k #FFD39B #BCEE68 #556B2E #FF6A6A #CD5555 #8B393A'.split()
    group_shifts = np.linspace(0.4, -0.4, len(group_names))

    # add in the background gray bars (this is such a kerchoo way to do this oops)
    gray_bars = np.ones(len(annotations))+3
    gray_bars[1::2] = 0
    ax.barh(range(len(annotations)), gray_bars, align='center', height=1, alpha=0.3, color='gray')
    ax.barh(range(len(annotations)), -gray_bars/4, align='center', height=1, alpha=0.3, color='gray')

    # plot the fold enrichment
    for i,group_name in enumerate(group_names):
        ax.scatter(log2FE_arr.loc[:,group_name], np.arange(len(annotations))+group_shifts[i], label=group_names[i], color=colors[i])
    ax.set_yticks(range(len(annotations)))
    ax.set_yticklabels(labels)
    ax.set_xlabel('$log_{2}$(Fold Enrichment)        ', fontsize = 30)
    ax.axvline(x=0, c='k', ls='--', lw=.5)

    # plot the proportion of variants
    bar_height = 0.1
    for i, group_name in enumerate(group_names):
        legend_label = f'{group_names[i]}, n={fm_dict[group_names[i]].shape[0]}'
        ax1.barh(np.arange(len(annotations))+group_shifts[i], mean_arr.loc[:,group_name], label=legend_label, color=colors[i], height = bar_height)
    ax1.legend(bbox_to_anchor=(1,1), loc="upper left", fontsize=20)
    ax1.set_xlabel('Prop. of Variants', fontsize = 30)

    fig.tight_layout()
    fig.savefig('gtex_annot_enrich.png', dpi=300)


if __name__ == '__main__':
    main()