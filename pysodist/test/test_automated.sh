set -e
set -x

rm -rf ./*.log
rm -rf ./test_result
./test_parse.sh
./test_extract.sh
./test_fisodist.sh
./test_plot.sh
