# this is to filter the results from genie3. For each transcription factor, only take the top 20 highest links.

# read in the results 
CD4_pos_go_markers = CD4_pos_markers_exp %>%
    arrange(desc(tf)) %>%
    top_n(n = 20, wt = tf)
    
# write out the results.
