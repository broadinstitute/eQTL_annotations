version 1.0

import "https://api.firecloud.org/ga4gh/v1/tools/landerlab:copy_to_google_bucket/versions/3/plain-WDL/descriptor" as copyfiles

workflow annotate_eqtl_variants {
    input {
        Array[File] finemapped_results
        Array[String] fm_group_names
        Array[File] peakfiles
        Array[String] peakfile_names
        File ABC_peakfile
        File plink_bim_path
        File gtex_vep
        File? peak_predictor_file
        String dest_dir
        String git_branch = "main"
    }

    call peak_overlaps {
        input:
            peakfiles=peakfiles,
            peakfile_names=peakfile_names,
            plink_bim_path=plink_bim_path,
            ABC_peakfile=ABC_peakfile,
            git_branch=git_branch
    }

    call gtex_vep_overlap {
        input:
            all_variant_peak_stats=peak_overlaps.all_variant_peak_stats,
            abc_variant_peak_stats=peak_overlaps.abc_variant_peak_stats,
            gtex_vep=gtex_vep,
            git_branch=git_branch
    }

    call merge_fm_annotations {
        input:
            finemapped_results=finemapped_results,
            fm_group_names=fm_group_names,
            all_variant_peaks_gtex=gtex_vep_overlap.all_variant_peaks_gtex,
            git_branch=git_branch
    }

    call make_gtex_annotation_plot {
        input:
            finemapped_annotations=merge_fm_annotations.fm_annotations,
            fm_group_names=fm_group_names,
            all_variant_peaks_gtex=gtex_vep_overlap.all_variant_peaks_gtex,
            git_branch=git_branch
    }

    call make_pip_bin_plot {
        input:
            finemapped_annotations=merge_fm_annotations.fm_annotations,
            fm_group_names=fm_group_names,
            all_variant_peaks_gtex=gtex_vep_overlap.all_variant_peaks_gtex,
            git_branch=git_branch
    }

    if (defined(peak_predictor_file)) {
        call ATAC_peak_predictor {
            input:
                peakfile=peak_predictor_file,
                finemapped_annotations=merge_fm_annotations.fm_annotations,
                fm_group_names=fm_group_names,
                git_branch=git_branch
        }
    }

    Array[File] all_var_results = [gtex_vep_overlap.all_variant_peaks_gtex]
    Array[File] peak_res_only = [peak_overlaps.all_variant_peak_stats]
    Array[File] gtex_annotation_plot = [make_gtex_annotation_plot.gtex_annotations_plot]
    Array[File] all_files = flatten([all_var_results, peak_res_only, merge_fm_annotations.fm_annotations, gtex_annotation_plot, make_pip_bin_plot.pip_bin_plots])

    call copyfiles.copyFile {
        input:
        files_2_copy = all_files,
        output_gs_dir = dest_dir
    }

    output {
        File all_variant_peaks_gtex = gtex_vep_overlap.all_variant_peaks_gtex
        File peak_overlaps_bed = peak_overlaps.all_variant_peak_stats
        Array[File] finemapped_annotations = merge_fm_annotations.fm_annotations
        File gtex_annotations_plot = make_gtex_annotation_plot.gtex_annotations_plot
        Array[File] pip_bin_plots = make_pip_bin_plot.pip_bin_plots
        Array[File]? peak_pip_plots = ATAC_peak_predictor.peak_pip_plots
        Array[File]? roc_curves = ATAC_peak_predictor.roc_curves
        Array[File]? model_beta_dfs = ATAC_peak_predictor.model_beta_dfs
        Array[File]? model_peak_predictions = ATAC_peak_predictor.model_peak_predictions
        Array[File]? peakname_finemapped_variants = ATAC_peak_predictor.peakname_finemapped_variants
    }
}

task peak_overlaps {
    input {
        Array[File] peakfiles
        Array[String] peakfile_names
        File ABC_peakfile
        File plink_bim_path
        String git_branch
    }

    command <<<
        set -ex
        (git clone https://github.com/broadinstitute/eQTL_annotations.git /app ; cd /app ; git checkout ${git_branch})
        micromamba -n tools2 install bioconda::bedtools -y

        # convert plink.bim to a bed.gz with chr, pos, pos, variant_id
        cat ~{plink_bim_path} | awk -F '\t' -v OFS='\t' '{print $1, $4, $4, $2}' | gzip -c > bim_to_bed.bed.gz
        zcat bim_to_bed.bed.gz | head
        zcat bim_to_bed.bed.gz | sort -k1,1 -k2,2n  | gzip -c > bim_to_bed.sorted.bed.gz

        for peakfile in ~{sep=" " peakfiles}
        do
            echo $peakfile
            zcat $peakfile | sort -k1,1 -k2,2n | cut -f -3 > tmp.bed
            # overwrite
            cat tmp.bed | gzip -c > $peakfile
            #micromamba run -n tools2 python3 /app/eqtl_annotations/combine_peaks_fm.py -p peak_dists.bed.gz -a $peakfile -f finemapped_results.tsv -g ${NAME_ARRAY[$i]}
            #cat finemapped_results.tsv > ${NAME_ARRAY[$i]}_ATAC_CHIP_peaks_finemapped_results.tsv
        done
        micromamba run -n tools2 bedtools closest -d -a bim_to_bed.sorted.bed.gz -b ~{sep=" " peakfiles} -names ~{sep=" " peakfile_names} -t first | gzip -c > peak_dists.bed.gz
        zcat peak_dists.bed.gz | head

        echo ~{ABC_peakfile}
        # include genes
        zcat ~{ABC_peakfile} | sort -k1,1 -k2,2n | cut -f -3,8 > tmp.bed
        cat tmp.bed | gzip -c > ABC_peakfile.bed.gz
        # run bedtools closets with -t all for abc
        micromamba run -n tools2 bedtools closest -d -a bim_to_bed.sorted.bed.gz -b ABC_peakfile.bed.gz -t all | gzip -c > ABC_peak_dists.bed.gz
        zcat ABC_peak_dists.bed.gz | head

    >>>

    output {
        File all_variant_peak_stats = "peak_dists.bed.gz"
        File abc_variant_peak_stats = "ABC_peak_dists.bed.gz"
    }

    runtime {
        docker: 'us.gcr.io/landerlab-atacseq-200218/hgrm_multiome_cluster_processing:0.6'
        cpu: 2
        memory: "32GB"
    }
}

task gtex_vep_overlap {
    input {
        File all_variant_peak_stats
        File abc_variant_peak_stats
        File gtex_vep
        String git_branch
    }

    Int disk_size = 50 + floor(size(gtex_vep, "GB"))

    command {
        set -ex
        (git clone https://github.com/broadinstitute/eQTL_annotations.git /app ; cd /app ; git checkout ${git_branch})
        micromamba run -n tools2 python3 /app/eqtl_annotations/annotate_gtex_vep.py -p ${all_variant_peak_stats} -a ${abc_variant_peak_stats} -g ${gtex_vep}
    }

    output {
        File all_variant_peaks_gtex = "all_variant_CHIP_ATAC_GTEx_overlap.parquet"
    }

    runtime {
        docker: 'us.gcr.io/landerlab-atacseq-200218/hgrm_multiome_cluster_processing:0.6'
        cpu: 4
        memory: "128GB"
        preemptible: 1
        disks: "local-disk ~{disk_size} HDD"
    }
}

task merge_fm_annotations {
    input {
        Array[File] finemapped_results
        Array[String] fm_group_names
        File all_variant_peaks_gtex
        String git_branch
    }

    command {
        set -ex
        (git clone https://github.com/broadinstitute/eQTL_annotations.git /app ; cd /app ; git checkout ${git_branch})
        micromamba run -n tools2 python3 /app/eqtl_annotations/merge_fm_annotations.py -f ${sep=' ' finemapped_results} -n ${sep=' ' fm_group_names} -v ${all_variant_peaks_gtex}
    }

    output {
        Array[File] fm_annotations = glob("*_fm_variants_annotations.parquet")
    }

    runtime {
        docker: 'us.gcr.io/landerlab-atacseq-200218/hgrm_multiome_cluster_processing:0.6'
        cpu: 4
        memory: "32GB"
    }
}


task make_gtex_annotation_plot {
    input {
        Array[File] finemapped_annotations
        Array[String] fm_group_names
        File all_variant_peaks_gtex
        String git_branch
    }

    command {
        set -ex
        (git clone https://github.com/broadinstitute/eQTL_annotations.git /app ; cd /app ; git checkout ${git_branch})
        micromamba run -n tools2 python3 /app/eqtl_annotations/gtex_annotation_plot.py -f ${sep=' ' finemapped_annotations} -g ${sep=' ' fm_group_names} -v ${all_variant_peaks_gtex}
    }

    output {
        File gtex_annotations_plot = "gtex_annot_enrich.png"
        File enrichment_means_by_group = "raw_mean_by_group_gtex_plot.tsv"
    }

    runtime {
        docker: 'us.gcr.io/landerlab-atacseq-200218/hgrm_multiome_cluster_processing:0.6'
        cpu: 4
        memory: "32GB"
    }
}

task make_pip_bin_plot {
    input {
        Array[File] finemapped_annotations
        Array[String] fm_group_names
        File all_variant_peaks_gtex
        String git_branch
    }

    command {
        set -ex
        (git clone https://github.com/broadinstitute/eQTL_annotations.git /app ; cd /app ; git checkout ${git_branch})
        micromamba run -n tools2 python3 /app/eqtl_annotations/make_pip_bin_plot.py -f ${sep=' ' finemapped_annotations}  -g ${sep=' ' fm_group_names} -v ${all_variant_peaks_gtex}
    }

    output {
        Array[File] pip_bin_plots = glob("*_annotation_by_pip.png")
        Array[File] mean_array_by_pip = glob("*_mean_array_by_pip.tsv")
    }

    runtime {
        docker: 'us.gcr.io/landerlab-atacseq-200218/hgrm_multiome_cluster_processing:0.6'
        cpu: 1
        memory: "16GB"
    }
}

task ATAC_peak_predictor {
    input {
        File? peakfile #this has to be ? bc wdl is dumb. but it has to be defined for the task to run.
        Array[File] finemapped_annotations
        Array[String] fm_group_names
        String git_branch
    }

    command <<<
        set -ex
        (git clone https://github.com/broadinstitute/eQTL_annotations.git /app ; cd /app ; git checkout ${git_branch})
        micromamba -n tools2 install bioconda::bedtools -y

        # this is assuming .bed not .bed.gz but can edit this to make it a changeable input.
        cat ~{peakfile} | sort -k1,1 -k2,2n | cut -f -3 > peakfile.bed
        for fm_file in ~{sep=" " finemapped_annotations}
        do
            echo $fm_file
            micromamba run -n tools2 python3 /app/eqtl_annotations/get_finemap_bed.py $fm_file
            zcat sample_vars.bed.gz | sort -k1,1 -k2,2n  | gzip -c > sample_vars.sorted.bed.gz
            micromamba run -n tools2 bedtools closest -t first -d -a sample_vars.sorted.bed.gz -b peakfile.bed | gzip -c > peak_dists.bed.gz
            micromamba run -n tools2 python3 /app/eqtl_annotations/peak_predictor.py -f $fm_file -g ~{sep=" " fm_group_names}
        done

    >>>

    output {
        Array[File] peak_pip_plots = glob("*_peak_pip_enrichment_by_category.png")
        Array[File] roc_curves = glob("*_peak_predictor_roc.png")
        Array[File] model_beta_dfs = glob("*_peaks_betas.parquet")
        Array[File] model_peak_predictions = glob("*_peak_preds.parquet")
        Array[File] peakname_finemapped_variants = glob("_peak_names_with_finemapped.parquet")
    }

    runtime {
        docker: 'us.gcr.io/landerlab-atacseq-200218/hgrm_multiome_cluster_processing:0.6'
        cpu: 4
        memory: "32GB"
    }

}