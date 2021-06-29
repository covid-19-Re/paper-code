HOW TO RECREATE FIG S11

1) Run the Re estimation pipeline for Switzerland,
save the result as CHE-estimates_no_imports.rds.

2) Temporarily replace the format_CHE_data.R file of the pipeline
by format_CHE_data_S11.R

3) Run the Re estimation pipeline for Switzerland,
save the result as CHE-Estimates_imports.rds

4) Run make_fig_S11.R after filling the directory location variables
   with the appropriate values for your setup.
