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
            for peakfile in ${peakfiles}
            do
                micromamba run -n tools2 python3 <<CODE
                    import pandas as pd
                    if '${finemap_result}'.endswith('.tsv'):
                        finemapped_df = pd.read_table('${finemap_result}')
                    else:
                        finemapped_df = pd.read_parquet('${finemap_result}')

                    if 'Variant_id' in finemapped_df.columns:
                        finemapped_df = finemapped_df.rename(columns={'Variant_id':'variant_id'})
                    if 'variant_id' not in finemapped_df.columns:
                        raiseValueError('column variant_id needs to be in data.')
                    finemapped_df['variant_id'] = finemapped_df.variant_id.str.replace('_', ':')
                    var_split = finemapped_df.variant_id.str.split(':')
                    finemapped_df = pd.concat([finemapped_df, pd.DataFrame({'chr':var_split.str[0], 'pos':var_split.str[1]})], axis=1)
                    fm_file = finemapped_df[['chr', 'pos', 'pos']]
                    fm_file.dropna(inplace=True)
                    fm_file.to_csv('sample_vars.bed.gz', sep='\t', header=None, index=None)
                CODE

                zcat sample_vars.bed.gz | sort -k1,1 -k2,2n  | gzip -c > sample_vars.sorted.bed.gz
                zcat ${peakfile} | sort -k1,1 -k2,2n | cut -f -3 > temp_file.bed
                micromamba run -n tools2 bedtools closest -d -a sample_vars.sorted.bed.gz -b temp_file.bed | gzip -c > peak_dists.bed.gz
            done
        done

    }

    runtime {
        docker: 'us.gcr.io/landerlab-atacseq-200218/hgrm_multiome_cluster_processing:0.6'
    }
}