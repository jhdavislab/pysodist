set -e
set -x

pysodist extract_spectra Glu_15N_300_DIA.mzML ./test_result/raw_spec/Glu_15N_300_DIA/pd_parsed_report.tsv --labeling N15 --logfile test_extract.log
pysodist extract_spectra Glu_15N_300_DIA.mzML ./test_result/interp_spec/Glu_15N_300_DIA/pd_parsed_report.tsv --labeling N15 --interp_only --interp_res 0.001 --logfile test_extract.log
pysodist extract_spectra Glu_15N_300_DIA.mzML ./test_result/sum_spec/Glu_15N_300_DIA/pd_parsed_report.tsv --labeling N15 --sum_only --logfile test_extract.log
