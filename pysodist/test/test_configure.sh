set -e
set -x

rm -rf test.config
rm -rf test_configure.log
pysodist configure test.config --logfile test_configure.log
./test_fisodist_conf.sh
