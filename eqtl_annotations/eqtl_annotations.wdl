version 1.0

workflow annotate_eqtl_variants {
    input {
        Array[File] finemapped_results
        Array[String] fm_group_names
        Array[File] peakfiles
        File gtex_vep
        String git_branch = "main"
    }

    call peak_overlaps {
        input:
            finemapped_results=finemapped_results,
            fm_group_names=fm_group_names,
            peakfiles=peakfiles,
            git_branch=git_branch
    }

    call gtex_vep_overlap {
        input:
            finemapped_results=peak_overlaps.fm_with_peak_distance,
            fm_group_names=fm_group_names,
            gtex_vep=gtex_vep,
            git_branch=git_branch
    }

    output {
        Array[File] fm_annotations = gtex_vep_overlap.fm_peaks_gtex
        Array[File] all_variant_peak_stats = peak_overlaps.all_variant_peak_stats
    }
}

task peak_overlaps {
    input {
        Array[File] finemapped_results
        Array[String] fm_group_names
        Array[File] peakfiles
        String plink_bim_path
        String git_branch
    }

    Int iters = length(finemapped_results) - 1

    command <<<
        set -ex
        (git clone https://github.com/broadinstitute/eQTL_annotations.git /app ; cd /app ; git checkout ${git_branch})
        micromamba -n tools2 install bioconda::bedtools -y

        # convert plink.bim to a bed.gz with chr, pos, pos, variant_id
        cat ~{plink_bim_path} | awk -F '\t' -v OFS='\t' '{print $1, $4, $4, $2}' | gzip -c > bim_to_bed.bed.gz
        zcat bim_to_bed.bed.gz | head
        zcat bim_to_bed.bed.gz | sort -k1,1 -k2,2n  | gzip -c > bim_to_bed.sorted.bed.gz
        FM_ARRAY=(~{sep=" " finemapped_results}) # Load array into bash variable
        NAME_ARRAY=(~{sep=" " fm_group_names})

        for i in $(seq 0 ~{iters}) #finemap_result in ${finemapped_results}
        do
            echo ${FM_ARRAY[$i]}
            micromamba run -n tools2 python3 /app/eqtl_annotations/get_finemap_bed.py ${FM_ARRAY[$i]} ${NAME_ARRAY[$i]}
            for peakfile in ~{sep=" " peakfiles}
            do
                echo $peakfile
                zcat $peakfile | sort -k1,1 -k2,2n | cut -f -3 > peak_file.bed
                micromamba run -n tools2 bedtools closest -d -a bim_to_bed.sorted.bed.gz -b peak_file.bed | gzip -c > peak_dists.bed.gz
                zcat peak_dists.bed.gz | head
                micromamba run -n tools2 python3 /app/eqtl_annotations/combine_peaks_fm.py -p peak_dists.bed.gz -a $peakfile -f finemapped_results.tsv -g ${NAME_ARRAY[$i]}
                cat finemapped_results.tsv | head
            done
            cat finemapped_results.tsv > ${NAME_ARRAY[$i]}_ATAC_CHIP_peaks_finemapped_results.tsv
        done

    >>>

    output {
        Array[File] fm_with_peak_distance = glob("*_ATAC_CHIP_peaks_finemapped_results.tsv")
        Array[File] all_variant_peak_stats = glob("*_peak_stats.tsv")
    }

    runtime {
        docker: 'us.gcr.io/landerlab-atacseq-200218/hgrm_multiome_cluster_processing:0.6'
    }
}

task gtex_vep_overlap {
    input {
        Array[String] finemapped_results
        Array[String] fm_group_names
        File gtex_vep
        String git_branch
    }

    Int disk_size = 1 + floor(size(gtex_vep, "GB"))

    command {
    set -ex
    (git clone https://github.com/broadinstitute/eQTL_annotations.git /app ; cd /app ; git checkout ${git_branch})
    micromamba run -n tools2 python3 /app/eqtl_annotations/annotate_gtex_vep.py -f ${sep=' ' finemapped_results} -n ${sep=' ' fm_group_names} -g ${gtex_vep}
    }

    output {
        Array[File] fm_peaks_gtex = glob("*_finemap_CHIP_ATAC_GTEx_overlap.tsv")
    }

    runtime {
        docker: 'us.gcr.io/landerlab-atacseq-200218/hgrm_multiome_cluster_processing:0.6'
        cpu: 1
        memory: "16GB"
        preemptible: 1
        disks: "local-disk ~{disk_size} HDD"
    }
}