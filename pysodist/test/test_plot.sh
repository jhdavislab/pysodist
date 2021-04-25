set -e
set -x

pysodist plot_spectra ./test_result/raw_spec/Glu_15N_300_DIA/U_var500N_fix998N_isodist_outputs/U_var500N_fix998N_output.csv ./test_result/raw_spec/Glu_15N_300_DIA/U_var500N_fix998N_isodist_fits ./test_result/raw_spec/Glu_15N_300_DIA/U_var500N_fix998N_plots --no_pdf --logfile ./test_plot_raw.log

pysodist plot_spectra ./test_result/sum_spec/Glu_15N_300_DIA/U_var500N_fix998N_isodist_outputs/U_var500N_fix998N_output.csv ./test_result/sum_spec/Glu_15N_300_DIA/U_var500N_fix998N_isodist_fits ./test_result/sum_spec/Glu_15N_300_DIA/U_var500N_fix998N_plots --logfile ./test_plot_raw.log
