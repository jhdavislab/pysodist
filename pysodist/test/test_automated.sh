set -e
set -x

./test_parse.sh
./test_extract.sh
./test_fisodist.sh
./test_plot.sh
