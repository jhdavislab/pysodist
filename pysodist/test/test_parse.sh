set -e
set -x

mkdir ./test_result
pysodist parse_input HFX_15Nmetabolic_test.csv --output_directory ./test_result/raw_spec --protein_list ENO2 ADH1 TEF1 EFT1 RPL28 PAB1 ALD6 TPI1 --logfile test_parse.log --isotope light --sample_list Glu_15N_300_DIA
pysodist parse_input HFX_15Nmetabolic_test.csv --output_directory ./test_result/interp_spec --logfile test_parse.log --isotope light --sample_list Glu_15N_300_DIA
pysodist parse_input HFX_15Nmetabolic_test.csv --output_directory ./test_result/sum_spec --logfile test_parse.log
