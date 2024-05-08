version 1.0

workflow annotate_eqtl_variants {
    input {
        Array[File] finemapped_results
        Array[File] peakfiles
    }

    call peak_overlaps {
        input:
            finemapped_results=finemapped_results,
            peakfiles=peakfiles
    }
}

task peak_overlaps {
    input {
        Array[File] finemapped_results
        Array[File] peakfiles
    }

    command {
        micromamba -n tools2 install bioconda::bedtools -y
        for finemap_result in ${finemapped_results}
        do
            echo $finemap_result
            for peakfile in ${peakfiles}
            do
                echo $peakfile
                micromamba run -n tools2 python3 get_finemap_bed.py $finemap_result
                zcat sample_vars.bed.gz | sort -k1,1 -k2,2n  | gzip -c > sample_vars.sorted.bed.gz
                zcat $peakfile | sort -k1,1 -k2,2n | cut -f -3 > temp_file.bed
                micromamba run -n tools2 bedtools closest -d -a sample_vars.sorted.bed.gz -b temp_file.bed | gzip -c > peak_dists.bed.gz
            done
        done

    }

    runtime {
        docker: 'us.gcr.io/landerlab-atacseq-200218/hgrm_multiome_cluster_processing:0.6'
    }
}