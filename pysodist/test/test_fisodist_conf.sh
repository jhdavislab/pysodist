set -e
set -x

pysodist run_isodist ./test_result/interp_spec/Glu_15N_300_DIA/pd_exported_peaks.tsv --config test.config --threads 30 --wait_time 10 --pysodist_input ./test_result/raw_spec/Glu_15N_300_DIA/pd_parsed_report.tsv --logfile test_fisodist.log
pysodist plot_spectra ./test_result/interp_spec/Glu_15N_300_DIA/U_var500N_fix998N_isodist_outputs/U_var500N_fix998N_output.csv ./test_result/interp_spec/Glu_15N_300_DIA/U_var500N_fix998N_isodist_fits ./test_result/interp_spec/Glu_15N_300_DIA/U_var500N_fix998N_plots --logfile ./test_plot.log
