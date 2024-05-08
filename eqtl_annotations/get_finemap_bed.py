import sys
import pandas as pd

finemap_result = sys.argv[1]
group_name = sys.argv[2]

if finemap_result.endswith('.tsv'):
    finemapped_df = pd.read_table(finemap_result)
else:
    finemapped_df = pd.read_parquet(finemap_result)

if 'Variant_id' in finemapped_df.columns:
    finemapped_df = finemapped_df.rename(columns={'Variant_id':'variant_id'})
if 'variant_id' not in finemapped_df.columns:
    raise ValueError('column variant_id needs to be in data.')

finemapped_df['variant_id'] = finemapped_df.variant_id.str.replace('_', ':')
var_split = finemapped_df.variant_id.str.split(':')
finemapped_df = pd.concat([finemapped_df, pd.DataFrame({'chr':var_split.str[0], 'pos':var_split.str[1]})], axis=1)
fm_file = finemapped_df[['chr', 'pos', 'pos']]
fm_file.dropna(inplace=True)
fm_file.to_csv('sample_vars.bed.gz', sep='\t', header=None, index=None)
finemapped_df.to_csv(f'finemapped_results.tsv', sep='\t', header=True)