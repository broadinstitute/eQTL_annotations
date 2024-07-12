# eQTL_annotations
Village finemapped eQTL annotation pipeline.

## Main pipeline
Overall: Annotates variants from our finemapped sets (eQTL pipeline) with information about if they are in peaks (e.g. ATAC or CHIP, etc.) and GTEx annotated regions (e.g. promoter, missense, splice site, etc.)
* Peak overlaps - annotates all variants in the VCF with their closest peak for each peakfile
    * Also, takes in ABC specific input, for gene-specific peak pairs (peaks are associated with specific genes in this peakfile)
* Merges all peak variant annotations with finemapped variants
* Merges GTEx annotations (enhancer, promoter, etc.) with finemapped variants
* Plots enrichment for each annotation category by pip bins: PIP < 0.01, 0.01 < PIP <0.1, 0.1 < PIP < 0.5, 0.5 < PIP < 0.9, 0.9 < PIP, for each day (each fiinemapped file)
* Plot all finemapped files together, with enrichment & proportion of variants in each annotation category across all the variants

## Peak Predictions
* Variants each have a set of annotations, a distance, a gene link, and a PIP
* Label each variant-annot-gene-pip with a peakname
* For each peak-gene pair:
    * take "or" of all annots (if any variant linked to that peak-gene is an enhancer, promoter, etc.)
    * take average of distance to gene TSS
    * take max of PIP
* Output = peak-maxpip-annots-distance-gene

Model: a logistic regression model that takes in each peak-gene pair and the output from above.
Train on high and low PIP peak-gene pairs (PIP > 0.1, PIP < 0.03>) in all other chromosomes (leave-one-out chr) to get AUC.
Score all peak-gene pairs w/ leave-one-out chromosome training.

Similar to ABC eQTL models, we remove any peaks that have variants in a promoter region, 3' or 5' UTR, frameshift, missense, stop-gain and splice variants. We also remove variants < 250 bp from the TSS.

Each peak-gene gets assigned a predicted probability of being a 'regulatory region', AKA containing a high-PIP variant.