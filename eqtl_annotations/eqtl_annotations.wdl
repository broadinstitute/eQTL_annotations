version 1.0

workflow annotate_eqtl_variants {
    input {
        Array[File] finemapped_results
        Array[String] fm_group_names
        Array[File] peakfiles
    }

    call peak_overlaps {
        input:
            finemapped_results=finemapped_results,
            fm_group_names=fm_group_names,
            peakfiles=peakfiles
    }

    output {
        Array[File] fm_with_peak_distance = peak_overlaps.fm_with_peak_distance
    }
}

task peak_overlaps {
    input {
        Array[File] finemapped_results
        Array[String] fm_group_names
        Array[File] peakfiles
    }

    Int fm_length = length(finemapped_results)

    command <<<
        micromamba -n tools2 install bioconda::bedtools -y
        FM_ARRAY=(~{sep=" " finemapped_results}) # Load array into bash variable
        NAME_ARRAY=(~{sep=" " fm_group_names})

        for i in $(seq 0 ${fm_length}) #finemap_result in ${finemapped_results}
        do
            echo ${FM_ARRAY[$i]}
            micromamba run -n tools2 python3 get_finemap_bed.py ${FM_ARRAY[$i]} ${NAME_ARRAY[$i]}
            zcat sample_vars.bed.gz | sort -k1,1 -k2,2n  | gzip -c > sample_vars.sorted.bed.gz
            for peakfile in ~{sep=" " peakfiles}
            do
                echo $peakfile
                zcat $peakfile | sort -k1,1 -k2,2n | cut -f -3 > temp_file.bed
                micromamba run -n tools2 bedtools closest -d -a sample_vars.sorted.bed.gz -b temp_file.bed | gzip -c > peak_dists.bed.gz
                micromamba run -n tools2 python3 combine_peaks_fm.py -p peak_dists.bed.gz -a $peakfile -f finemapped_results.tsv -g ${NAME_ARRAY[$i]}
            done
            cat finemapped_results.tsv > ${NAME_ARRAY[$i]}_ATAC_CHIP_peaks_finemapped_results.tsv
        done

    >>>

    output {
        Array[File] fm_with_peak_distance = glob("*_ATAC_CHIP_peaks_finemapped_results.tsv")
    }

    runtime {
        docker: 'us.gcr.io/landerlab-atacseq-200218/hgrm_multiome_cluster_processing:0.6'
    }
}