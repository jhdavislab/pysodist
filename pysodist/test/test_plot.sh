set -e
set -x

pysodist plot_spectra HFX_15Nmetabolic_test.csv --output_directory ./test_result/raw_spec --protein_list ENO2 TEF1 EFT1 RPL28 PAB1 ALD6 TPI1 PGK1 --logfile test_parse.log --isotope light --sample_list Glu_15N_300_DIA