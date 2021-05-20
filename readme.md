Script to generate the meta_PKN used in the COSMOS Paper.

Run in order:
01_recon_matlab_to_df.R
01b_list_coenzymes.R
02_STITCH_to_SIF.R
03revision_join_with_ominpath_and_stitch_met_filtered.R
04revision_solve_exchange_problem_metfiltered.R

Metabolic network that the PKN relies was obtained from: https://www.vmh.life/#downloadview

Mind that we have generated a new PKN since the publication, that is lighter and better curated. We strongly advise users to use this one instead: https://github.com/saezlab/cosmos_meta_PKN_redhuman
