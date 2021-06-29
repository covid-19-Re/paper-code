HOW TO RECREATE FIG S10

1) Run the Re estimation pipeline on Swiss pipeline ignoring delay variability
(temporarily remove the line 120 from the format_CHE_data.R file to do so:
readr::write_csv(final_delay_data_FOPH, path = file.path(outDir, "CHE_data_delays.csv")))
2) Run the Re estimation pipeline on Swiss pipeline ignoring symptom onset data
   (temporarily modify format_CHE_data.R to do so)
3) Run the Re estimation doing both and run the Re pipeline doing none of the two above,
each time saving the results in a different .rds file
4) Run the make_fig_S10.R file
