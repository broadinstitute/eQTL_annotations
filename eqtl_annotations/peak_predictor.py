import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.pyplot import cm
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import precision_recall_curve, PrecisionRecallDisplay
from sklearn import metrics
import os.path


"""
All helper functions for running the model and plotting the results.
"""
def abs_log1p(x):
    return np.log2(abs(x) + 1)


def chrom_train_test_split(pos, neg, chroms, high_pips, low_pips, class_weight=None):
    """
    Leave one chromosome out training and testing on high/low data only.
    """
    y_real_all = []
    y_prob_all = []
    peak_gene_order = []
    np.random.seed(1)
    for chr_str in chroms:
        print(chr_str)
        train_pos = high_pips.chr != chr_str
        train_neg = low_pips.chr != chr_str
        X_train = pd.concat([pos[train_pos], neg[train_neg]])
        X_test = pd.concat([pos[~train_pos], neg[~train_neg]])
        peak_gene_order = peak_gene_order + list(high_pips[~train_pos].index.values)
        peak_gene_order = peak_gene_order + list(low_pips[~train_neg].index.values)

        y_train = np.array([1] * sum(train_pos == 1) + [0] * sum(train_neg == 1))
        y_test = np.array([1] * sum(train_pos == 0) + [0] * sum(train_neg == 0))
        regr = LogisticRegression(
            random_state=1, max_iter=1000, class_weight=class_weight
        )

        # train the predictor
        regr.fit(X_train, y_train)
        y_prob = regr.predict_proba(X_test)
        y_prob = y_prob[:, 1]

        # save the prob
        y_real_all = y_real_all + list(y_test)
        y_prob_all = y_prob_all + list(y_prob)

    return y_real_all, y_prob_all, peak_gene_order


# Full model
def run_full_lr(high_pips, low_pips, annotation_cols, var_annots_to_pred):
    """
    Leave one chr out training, but predict on full data.
    """
    high_pip_pos = high_pips.loc[:, annotation_cols].fillna(0)
    low_pips_neg = low_pips.loc[:, annotation_cols].fillna(0)
    # and turn the log1p
    pos = high_pip_pos.apply(lambda x: abs_log1p(x))
    neg = low_pips_neg.apply(lambda x: abs_log1p(x))

    prob_dfs = {}

    np.random.seed(1)
    for chr_str in var_annots_to_pred.chr.unique():
        print(chr_str)
        train_pos = high_pips.chr != chr_str
        train_neg = low_pips.chr != chr_str
        X_train = pd.concat([pos[train_pos], neg[train_neg]])

        y_train = np.array([1] * sum(train_pos == 1) + [0] * sum(train_neg == 1))
        regr = LogisticRegression(
            random_state=1, max_iter=1000, class_weight="balanced"
        )

        # train the predictor
        regr.fit(X_train, y_train)

        # predict on all chr variants/peaks
        chr_data = var_annots_to_pred[var_annots_to_pred.chr == chr_str]
        chr_annots = (chr_data.loc[:, annotation_cols].fillna(0).apply(lambda x: abs_log1p(x)))

        # save the betas from chrom20
        if chr_str == "chr20":
            betas_df = pd.DataFrame(
                data={"betas": regr.coef_[0], "feature": chr_annots.columns}
            )

        prob_df = pd.DataFrame(
            {
                "peak_name": chr_data.peak_name,
                "phenotype_id": chr_data.phenotype_id,
                "predict_prob_high": regr.predict_proba(chr_annots)[:, 1],
            }
        )
        prob_dfs[chr_str] = prob_df

        all_peak_probs = pd.concat(prob_dfs)
        all_annots_with_predictions = pd.concat([var_annots_to_pred.set_index(["peak_name", "phenotype_id"]),
                                                 all_peak_probs.set_index(["peak_name", "phenotype_id"])],
                                                 axis=1)

    return all_annots_with_predictions, betas_df, prob_dfs, regr


def make_enrichment_plot(annots_with_pred_bin, group_name):
    bins = [0, 0.01, 0.1, 0.5, 0.9, 1]
    labels = ["PIP<0.01", "0.01<PIP<0.1", "0.1<PIP<0.5", "0.5<PIP<0.9", "0.9<PIP"]
    annots_with_pred_bin["pip_bin"] = pd.cut(
        annots_with_pred_bin.pip,
        bins=bins,
        labels=labels,
        right=True,
        include_lowest=True,
    )

    # the annotations we care about. not very user-change friendly yet.
    annotations = [
        "mean_start_distance",
        "enhancer_d",
        "open_chromatin_region_d",
        "CTCF_binding_site_d",
        "TF_binding_site_d",
        "non_coding_transcript_exon_variant_d",
        "synonymous_variant_d",
        "ABC_in_a_peak",
        "ATAC_in_a_peak",
        "ATAC_D0_in_a_peak",
        "ATAC_D2_in_a_peak",
        "ATAC_D4_in_a_peak",
        "ATAC_D7_in_a_peak",
        "CTCF_D0_in_a_peak",
        "CTCF_D2_in_a_peak",
        "CTCF_D4_in_a_peak",
        "CTCF_D7_in_a_peak",
        "H3K27ac_in_a_peak",
        "H3K4me1_in_a_peak",
        "NANOG_in_a_peak",
        "OCT4_in_a_peak",
        "SOX2_in_a_peak",
        annots_with_pred_bin["q_bins"].cat.categories[0],
        annots_with_pred_bin["q_bins"].cat.categories[1],
        annots_with_pred_bin["q_bins"].cat.categories[2],
        annots_with_pred_bin["q_bins"].cat.categories[3],
    ]

    mean_arr = pd.DataFrame(0.0, index=annotations, columns=labels)
    for i, conseq in enumerate(annotations):
        for j, label in enumerate(labels):
            if conseq == "mean_start_distance":
                mean_arr.at[conseq, label] = (
                    annots_with_pred_bin.query("pip_bin==@label")
                    .groupby(["phenotype_id"], as_index=True)[conseq]
                    .agg("mean")
                    .mean()
                )
                continue
            annot_bin_df = (
                annots_with_pred_bin.query("pip_bin==@label")
                .groupby(["phenotype_id"], as_index=True)[annotations]
                .agg("any")
            )
            mean_arr.at[conseq, label] = annot_bin_df[
                conseq
            ].mean()  # ignores NaN by default

            norm = mean_arr["PIP<0.01"].values
    FE = (mean_arr.T / norm).T
    annotations = FE.index.values
    FE_melt_new = FE.reset_index().melt(id_vars="index")

    fig, ax = plt.subplots(figsize=(15, 6))
    annotation_map = {
        "ATAC_peak_dist": "ATAC peak dist",
        "CTCF_peak_dist": "CTCF peak dist",
        "enhancer_d": "Enhancer",
        "promoter_d": "Promoter",
        "CTCF_binding_site_d": "CTCF binding site",
        "TF_binding_site_d": "TF binding site",
        "3_prime_UTR_variant_d": "3' UTR",
        "5_prime_UTR_variant_d": "5' UTR",
        "intron_variant_d": "Intron",
        "missense_variant_d": "Nonsynonymous",
        "synonymous_variant_d": "Synonymous",
        "open_chromatin_region_d": "Open chromatin",
        "promoter_flanking_region_d": "Promoter Flanking",
        "frameshift_variant_d": "Frameshift Variant",
        "stop_gained_d": "Stop Gained",
        "non_coding_transcript_exon_variant_d": "Non-coding transcript exon variant",
    }

    palette = ["#929591", "#F97306", "#FE420F"]
    sns.set_palette(palette)
    sns.barplot(FE_melt_new, x="index", y="value", hue="variable")
    annotation_labels = [
        annotation_map[x] if x in annotation_map.keys() else x.replace("_", " ")
        for x in annotations
    ]
    ax.set_xticklabels(annotation_labels, rotation=45, ha="right", fontsize=10)
    ax.set_ylabel("Fold Enrichment", fontsize=25)
    ax.set_xlabel("")
    ax.legend(title="", fontsize=15)
    ax.set_title(f'{group_name} Peak Enrichment by Max PIP, \n with Model Bin Groups', fontsize=25)
    for loc in np.arange(0, 16, 2.5):
        ax.axhline(loc, c="k", ls="--", lw=0.35, zorder=0)
    # clip axis
    ax.set_ylim(0, 8)
    fig.tight_layout()
    fig.savefig(f'{group_name}_peak_pip_enrichment_by_category.png', dpi=300)

    return mean_arr, FE, fig


def add_q_bins(all_annots_with_predictions):
    all_annots_with_predictions["q_bins"] = pd.qcut(all_annots_with_predictions.predict_prob_high, 4)
    all_annots_with_predictions["q_bins"] = (all_annots_with_predictions.q_bins.cat.rename_categories(
            [
                f"{i.left} to {i.right}"
                for i in all_annots_with_predictions["q_bins"].cat.categories
            ]
        )
    )
    all_annots_with_predictions["pip"] = all_annots_with_predictions.max_pip
    res = pd.concat([all_annots_with_predictions,
                    pd.get_dummies(all_annots_with_predictions.q_bins)],
                    axis=1)

    return res


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f",dest="finemapped_annotations",
        help="parquet with all fm var annotations",
        type=str, required=True,
    )
    parser.add_argument("-g", dest="group_names", nargs='+', default=[])
    args = parser.parse_args()

    library = "ATAC"
    group_name = None
    for name in args.group_names:
        if name in os.path.basename(args.finemapped_annotations):
            group_name = name
            break

    # read in fm annotation data
    fm_annots_parq = pd.read_parquet(args.finemapped_annotations)
    var_split = fm_annots_parq.variant_id.str.split(":")
    fm_annots = pd.concat([fm_annots_parq, pd.DataFrame({"chr": var_split.str[0], "pos": var_split.str[1]})], axis=1)

    # read in recently generated peak overlaps
    peaks = pd.read_table(
        f"peak_dists.bed.gz",
        header=None,
        names=f"chr pos _pos variant_id _chr {library}_peak_start {library}_peak_end {library}_peak_dist".split(),
        index_col=[0, 1],
    )

    # This time, we save peak name
    peaks["peak_name"] = (
        peaks["_chr"]
        + "_"
        + peaks.ATAC_peak_start.astype(str)
        + "_"
        + peaks.ATAC_peak_end.astype(str)
    )
    peaks["ATAC_peak_dist_new"] = peaks.ATAC_peak_dist
    # merge with fm variant info
    tmp = fm_annots.merge(
        peaks[["peak_name", "variant_id", "ATAC_peak_dist_new"]],
        how="left",
        on="variant_id",
    )

    # drop duplicate rows
    tmp.drop_duplicates(inplace=True)
    # get in a peak instead of distance
    in_a_peak = (tmp.loc[:, tmp.columns.str.contains("peak_dist")] == 0).astype(int)
    in_a_peak.columns = in_a_peak.columns.str.strip("peak_dist") + "_in_a_peak"
    # add that info
    var_annots = pd.concat([tmp, in_a_peak], axis=1)

    # remove vars close to tss site
    var_annots_far = var_annots[var_annots.start_distance > 250]
    print(
        "number variants dropped due to start distance: ",
        var_annots.shape[0] - var_annots_far.shape[0],
    )

    vars_in_peaks = var_annots_far.query("ATAC_peak_dist_new_in_a_peak == 1")
    print(
        "number variants in an atac peak (after distance removal): ",
        vars_in_peaks.shape[0],
    )

    print(
        "number of variants in atac peak dropped due to start distance: ",
        var_annots.query("ATAC_peak_dist_new_in_a_peak == 1").shape[0]
        - vars_in_peaks.shape[0],
    )

    #Add output of annotations df with the peak names? before taking groupby? the var_annots df?
    vars_in_peaks.to_parquet(f'{group_name}_peak_names_with_finemapped.parquet')

    peak_gene_groups = vars_in_peaks.groupby(["phenotype_id", "peak_name"]).groups.keys()
    vars_in_peaks.set_index(["phenotype_id", "peak_name"], inplace=True)
    vars_in_peaks.sort_index(inplace=True)

    print("generating the info per peak-gene pair to train the model.")
    # include all annotations for now, we will drop later
    annotation_columns = vars_in_peaks.drop(columns="pos").loc[:, "enhancer_d":].columns
    # build a df
    peak_gene_df = pd.DataFrame(
        index=peak_gene_groups,
        columns=np.hstack(("max_pip", "mean_start_distance", annotation_columns)),
    )

    # fill the DF:
    # groupby(peakname, gene)
        # take "or" of all annots
        # take average (or whatever) of gene distance
        # take max of PIP
    # output = peak-maxpip-annots-distance-gene
    for group in peak_gene_groups:
        temp_df = vars_in_peaks.loc[group]
        peak_gene_df.loc[group, annotation_columns] = (temp_df.loc[:, annotation_columns].any().astype(int))
        peak_gene_df.loc[group, "mean_start_distance"] = temp_df.start_distance.mean()
        peak_gene_df.loc[group, "max_pip"] = temp_df.pip.max()
        peak_gene_df.loc[group, "chr"] = temp_df.chr[0]

    print(peak_gene_df.head())
    print("peak gene groups shape: ", peak_gene_df.shape)

    # remove any peaks with these annotations
    remove_cols = [
        "promoter_d",
        "promoter_flanking_region_d",
        "3_prime_UTR_variant_d",
        "5_prime_UTR_variant_d",
        "frameshift_variant_d",
        "missense_variant_d",
        "splice_acceptor_variant_d",
        "splice_donor_variant_d",
        "splice_region_variant_d",
        "stop_gained_d",
    ]

    peaks_to_keep = ~peak_gene_df[remove_cols].any(axis=1)
    peak_gene_df_filt = peak_gene_df.loc[peaks_to_keep]
    # NEW BOUNDARIES
    high_pips = peak_gene_df_filt.query("max_pip >= 0.1")
    low_pips = peak_gene_df_filt.query("max_pip <= 0.03")
    print("high pip shape: ", high_pips.shape, "low pip shape: ", low_pips.shape)
    print("filtered peak gene df shape: ", peak_gene_df_filt.shape)

    print("Now running the train test split to get accuracy")
    annotation_cols = [
        "mean_start_distance",
        "enhancer_d",
        "open_chromatin_region_d",
        "CTCF_binding_site_d",
        "TF_binding_site_d",
        "non_coding_transcript_exon_variant_d",
        "synonymous_variant_d",
        "ABC_in_a_peak",
        "ATAC_D0_in_a_peak",
        "ATAC_D2_in_a_peak",
        "ATAC_D4_in_a_peak",
        "ATAC_D7_in_a_peak",
        "CTCF_D0_in_a_peak",
        "CTCF_D2_in_a_peak",
        "CTCF_D4_in_a_peak",
        "CTCF_D7_in_a_peak",
        "H3K27ac_in_a_peak",
        "H3K4me1_in_a_peak",
        "NANOG_in_a_peak",
        "OCT4_in_a_peak",
        "SOX2_in_a_peak",
        "ATAC_peak_dist_new_in_a_peak",
    ]

    high_pip_pos = high_pips.loc[:, annotation_cols].fillna(0)
    low_pips_neg = low_pips.loc[:, annotation_cols].fillna(0)

    # and turn the log1p
    pos = high_pip_pos.apply(lambda x: abs_log1p(x))
    neg = low_pips_neg.apply(lambda x: abs_log1p(x))

    # run test with this data for accuracy
    y_real, y_prob, variant_order_model = chrom_train_test_split(
        pos,
        neg,
        peak_gene_df_filt.chr.unique(),
        high_pips,
        low_pips,
        class_weight="balanced",
    )
    precision, recall, thresholds = precision_recall_curve(y_real, y_prob)
    fig, ax = plt.subplots(figsize=(10,8))
    pr_display = PrecisionRecallDisplay(precision=precision, recall=recall).plot(ax=ax)
    ax.set_ylim(0)
    ax.set_title(f'{group_name} Train/Test ROC Curve, High and Low Peaks', fontsize=30)
    fpr, tpr, thresholds = metrics.roc_curve(y_real, y_prob, pos_label=1)
    ax.text(x=.8,y=.8, s=f'AUC: {round(metrics.auc(fpr, tpr), 3)}', fontsize=20)
    fig.savefig(f'{group_name}_peak_predictor_roc.png', dpi=300)
    print("AUC: ", metrics.auc(fpr, tpr))

    # set up
    peak_gene_df_filt.reset_index(inplace=True)
    peak_gene_df_filt = peak_gene_df_filt.rename(
        columns={"level_0": "phenotype_id", "level_1": "peak_name"}
    )

    print("Running full model.")
    all_annots_with_predictions, betas_df, probs_dfs, regr = run_full_lr(
        high_pips, low_pips, annotation_cols, peak_gene_df_filt
    )

    all_annots_with_predictions = add_q_bins(all_annots_with_predictions)

    all_annots_with_predictions = all_annots_with_predictions.merge(
        vars_in_peaks.groupby(["peak_name", "phenotype_id"])
        .size()
        .rename("num_vars_in_peak"),
        left_index=True,
        right_index=True,
    )

    print("Making enrichment plot")
    mean_arr, FE, fig = make_enrichment_plot(all_annots_with_predictions, group_name)
    betas_df = pd.concat([pd.DataFrame([[regr.intercept_[0], 'intercept']], columns=betas_df.columns), betas_df], ignore_index=True)

    all_annots_with_predictions.reset_index().to_parquet(f'{group_name}_peak_preds.parquet')
    betas_df.to_parquet(f'{group_name}_peaks_betas.parquet')
    print('Done.')


if __name__ == '__main__':
    main()
