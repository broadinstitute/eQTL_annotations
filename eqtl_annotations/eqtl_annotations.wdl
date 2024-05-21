version 1.0

import "https://api.firecloud.org/ga4gh/v1/tools/landerlab:copy_to_google_bucket/versions/3/plain-WDL/descriptor" as copyfiles

workflow annotate_eqtl_variants {
    input {
        Array[File] finemapped_results
        Array[String] fm_group_names
        Array[File] peakfiles
        Array[String] peakfile_names
        File plink_bim_path
        File gtex_vep
        String dest_dir
        String git_branch = "main"
    }

    call peak_overlaps {
        input:
            peakfiles=peakfiles,
            peakfile_names=peakfile_names,
            plink_bim_path=plink_bim_path,
            git_branch=git_branch
    }

    call gtex_vep_overlap {
        input:
            all_variant_peak_stats=peak_overlaps.all_variant_peak_stats,
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
    }
}

task peak_overlaps {
    input {
        Array[File] peakfiles
        Array[String] peakfile_names
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
            # micromamba run -n tools2 python3 /app/eqtl_annotations/get_finemap_bed.py ${FM_ARRAY[$i]} ${NAME_ARRAY[$i]}
            zcat $peakfile | sort -k1,1 -k2,2n | cut -f -3 > tmp.bed
            # overwrite
            cat tmp.bed | gzip -c > $peakfile
            #micromamba run -n tools2 python3 /app/eqtl_annotations/combine_peaks_fm.py -p peak_dists.bed.gz -a $peakfile -f finemapped_results.tsv -g ${NAME_ARRAY[$i]}
            #cat finemapped_results.tsv > ${NAME_ARRAY[$i]}_ATAC_CHIP_peaks_finemapped_results.tsv
        done
        micromamba run -n tools2 bedtools closest -d -a bim_to_bed.sorted.bed.gz -b ~{sep=" " peakfiles} -names ~{sep=" " peakfile_names} -t first | gzip -c > peak_dists.bed.gz
        zcat peak_dists.bed.gz | head

    >>>

    output {
        File all_variant_peak_stats = "peak_dists.bed.gz"
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
        File gtex_vep
        String git_branch
    }

    Int disk_size = 50 + floor(size(gtex_vep, "GB"))

    command {
    set -ex
    (git clone https://github.com/broadinstitute/eQTL_annotations.git /app ; cd /app ; git checkout ${git_branch})
    micromamba run -n tools2 python3 /app/eqtl_annotations/annotate_gtex_vep.py -p ${all_variant_peak_stats} -g ${gtex_vep}
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
        Array[File] fm_annotations = glob("*/*_fm_variants_annotations.parquet")
    }

    runtime {
        docker: 'us.gcr.io/landerlab-atacseq-200218/hgrm_multiome_cluster_processing:0.6'
        cpu: 1
        memory: "16GB"
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
        micromamba run -n tools2 python3 /app/eqtl_annotations/gtex_annotation_plot.py -f ${sep=' ' finemapped_annotations} -v ${all_variant_peaks_gtex}
    }

    output {
        File gtex_annotations_plot = "gtex_annot_enrich.png"
        File enrichment_means_by_group = "raw_mean_by_group_gtex_plot.tsv"
    }

    runtime {
        docker: 'us.gcr.io/landerlab-atacseq-200218/hgrm_multiome_cluster_processing:0.6'
        cpu: 1
        memory: "16GB"
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
        micromamba run -n tools2 python3 /app/eqtl_annotations/make_pip_bin_plot.py -f ${sep=' ' finemapped_annotations} -v ${all_variant_peaks_gtex}
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
