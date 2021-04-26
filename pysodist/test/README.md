--pysodist v0.0.4--

This test suite is used to ensure pysodist is correctly installed.
It uses a small set MS1 spectra derived from peptides isolated from S. cerevisiae.
The spectra include 3 species:
* unlabeled
* ~50% 15N metabolically labeled
* ~99.8% 15N  metabolically labeled

To test pysodist, you can run each test (test_parse, test_extract, test_fisodist, test_plot) individually, or you can run an automated test using test_automated.sh
Be sure to carefully read the output and note any errors/warnings. The result will be saved in test_result.

You can also test the configuration scripts as follows:
1) run test_automated.sh
2) run test_configure.sh