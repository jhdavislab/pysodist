set -e
set -x

pysodist_dir=${PWD%/*}
base_dir=${pysodist_dir%/*}
pysodist run_isodist ./test_result/raw_spec/Glu_15N_300_DIA/pd_exported_peaks.tsv $base_dir/fortran/isodist $base_dir/model_files/atoms.txt $base_dir/model_files/U_var500N_fix998N.txt --threads 30 --wait_time 10 --pysodist_input ./test_result/raw_spec/Glu_15N_300_DIA/pd_parsed_report.tsv --logfile test_fisodist.log
pysodist run_isodist ./test_result/sum_spec/Glu_15N_300_DIA/pd_exported_peaks.tsv $base_dir/fortran/isodist $base_dir/model_files/atoms.txt $base_dir/model_files/U_var500N_fix998N.txt --threads 2 --wait_time 10 --pysodist_input ./test_result/sum_spec/Glu_15N_300_DIA/pd_parsed_report.tsv --no_compress --logfile test_fisodist.log
