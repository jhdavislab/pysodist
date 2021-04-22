set -e
set - x

pysodist parse_input HFX_15Nmetabolic_test --output_directory ./test_result --sample_list [Glu_15N_300_DIA] --protein_list [ADH1 SOD1 EFT1] --isotope light --logfile test_parse.log
