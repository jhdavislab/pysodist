# Pysodist: pythonic isodist

Pysodist is a python-based tool to extract and fit mass spectra. This tool is most useful in fitting 'partially' labeled isotope envelopes, which tyipcally arise from metabolic labeling of cells using a mixture of light and heavy isotope sources. It relies on the Fourier transform convolution algorithm implemented in isodist.

## Manuscript:

The algorithm is described in 1. Utilization of a related tool can be found in 2-4

1. Sperling E, Bunner AE, Sykes MT, Williamson JR. Quantitative Analysis of Isotope Distributions in Proteomic Mass Spectrometry Using Least-Squares Fourier Transform Convolution. Analytical Chemistry, 2008. 80, pp. 4906-4917.

2. Chen SS, Sperling E, Silverman JM, Davis JH, Williamson JR. Measuring the dynamics of E. coli ribosome biogenesis using pulse-labeling and quantitative mass spectrometry. Molecular Biosystems, 2012. 8(12):3325-34.

3. Jomaa A*, Jain N*, Davis JH*, Williamson JR, Britton RA, Ortega J. Functional domains of the 50S subunit mature late in the assembly process. Nucleic Acids Research, 2014. 42(5):3419-35. 

4. Davis JH*, Tan YZ*, Carragher B, Potter CS, Lyumkis D, Williamson JR. Modular assembly of the bacterial large ribosomal subunit. Cell 2016. 167(6):1610-1622.

## New in v0.0.3
**Version 0.0.3:**
* New: Installation instructions provided in the README.md
* Last updated 11/17/2020

### Previous versions
**Version 0.0.2:**
* New: An interactive jupyter notebook to look at the results

**Version 0.0.1:**
* A pre-release beta.

## Installation/dependencies:

Pysodist requires the following python tools:
numpy, scipy, pandas, qgrid, matplotlib, seaborn, jupyter_notebooks [for interactive analysis], pyteomics [to extract from mzml files]

The follow operations will clone the repository to your machine and will install the required python libraries:

    # Create conda environment
    $ conda create --name pysodist python=3.7
    $ conda activate pysodist

    # Install dependencies
    $ conda install numpy scipy pandas qgrid matplotlib seaborn jupyter pyteomics -c bioconda
    
    # Clone source code and install
    $ git clone https://github.com/jhdavislab/pysodist.git
    $ cd pysodist
    $ git checkout v0.0.2
    $ python setup.py install

    # Test your installation
    $ pysodist --help

Additionally, one must be able to generate .mzml files from your instruments raw output. For our Thermo Q-Exactive HF-X, we use mscovert.exe, which is part of the Proteowizard toolbox, and can be downloaded here:
http://proteowizard.sourceforge.net/tools.shtml
http://proteowizard.sourceforge.net/download.html

We convert thermo .RAW files as follows:
.\msconvert.exe [YOUR_FILE] -o [DESIRED_OUTPUT_DIRECTORY]--mzML --64 -v mz64 --inten32 --noindex --filter "msLevel 1" --zlib

Finally, pysodist currently requires the Fortran implementation of the isodist fitting routine developed by Michael T. Sykes and James R. Williamson [1].
This can be installed as follows:

### Windows
    # A compiled execute is provided as pysodist/fortran/source/isodist_win.exe. The easiest option is to simply copy this executable to [your_pysodist_repo]/fortran/isodist_win.exe. Alternatively, you can compile directly below.
    For info on the compiler, go to: http://mingw-w64.org/doku.php, and http://mingw-w64.org/doku.php/download/mingw-builds
    Downoad the mingw installer from here: https://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win32/Personal%20Builds/mingw-builds/installer/mingw-w64-install.exe/download
    Run mingw-w64-install.exe
    Use the following options (note the x86_64 setting)
        Version 8.1.0
        Architecture: x86_64
        Threads: posix
        Exception: seh
        Build revision: 0
    Choose this path: C:\Program Files\mingw-w64\x86_64-8.1.0-posix-seh-rt_v6-rev0

    Start anaconda powershell
    $ cd [where_you_cloned_the_pysodist_repo]/fortran/source
    $ ./combile.bat [if this has an error, edit the file and make sure all of the paths listed are correct]
    $ cp isodist_win.exe ../.

### Linux
    # First, ensure that gfortran is installed. If not, it can be installed by your sys admin or, on Ubuntu, by running sudo apt install gfortran 
    # Compile the fortran source code on your local machine
    $ cd [your_cloned_pysodist_repo_directory]
    $ cd fortran/source
    $ make
    $ make install
    $ make clean

## Quickstart: Test isolated modules

### 1. Overview
Pysodist consists of 4 modules - parse_input, extract_spectra, run_isodist, plot_spectra. Each can be run individually as described below, or the entire pipeline can be run using the command full_pipeline. It is strongly recommended to run the modules individually when you are getting started, and to only run the full_pipeline once you are comfortable with using pysodist.

### 1. Prepare inputs (pysodist parse_input)

First, we will convert a Skyline report file to a standard pysodist input format. As requested, addtional parsers will be implemented to use other reports (e.g. from EncyclopeDIA, mascot, the TPP, etc.).

    $ pysodist parse_input -h
    usage: pysodist parse_input [-h] [--output_directory OUTPUT_DIRECTORY]
                                [--sample_list [SAMPLE_LIST [SAMPLE_LIST ...]]]
                                [--protein_list PROTEIN_LIST] [--isotope ISOTOPE]
                                [--q_value Q_VALUE]
                                input

    positional arguments:
      input                 input file to parse.

    optional arguments:
      -h, --help            show this help message and exit
      --output_directory OUTPUT_DIRECTORY
                            Output files will be saved in this folder: 1 directory
                            per sample in the skyline report. Default = ./
      --sample_list [SAMPLE_LIST [SAMPLE_LIST ...]]
                            An optional list of samples to parse. By default all
                            samples in the report are analyzed.
      --protein_list PROTEIN_LIST
                            An optional list of the proteins to parse. By default,
                            all proteins in the report are analyzed.
      --isotope ISOTOPE     Be default, it is assumed that the report contains a
                            light isotope (no special labeling), if this field is
                            not present in the report, you can specify a different
                            field here (e.g. "heavy")
      --q_value Q_VALUE     Used to optionally filter the report file based on the
                            q_value. By default, no q_value filtering is used.

This input file should be a skyline report file. A new file called pd_parsed_report.tsv will be created for each sample specified.

Example usage

    $  pysodist parse_input ./raw_data/ace_15N_pysodist_v002_trunc.csv --sample_list Ace_15N_600_DIA
    working on sample: Ace_15N_600_DIA
    reading peptides from skyline dataframe Ace_15N_600_DIA
    found 21 peptides.
    returning 11 peptides

### 2. Extract spectra (pysodist extract_spectra)

Second, we will extract each of the individual m/z vs. intensity MS1 spectra from the .mzml file. Note that depending on the number of MS1 scans you have taken per minute, and the chromatographic width of your peptide peaks, you may have many MS1 spectra per peptide. These are individually extracted. Additionally, they are interpolated at a default resolution of 0.001 m/z units, and the interpolated spectra are summed, resulting in a composite '_SUM' spectra for the entire peak. Narrowing your retention time windows in skyline will provide additional specificity if you find that you are extracting spectra that are contaminated with other peaks.
    
    $ pysodist extract_spectra -h
    usage: pysodist extract_spectra [-h] [--labeling LABELING] [--interp_only]
                                    [--sum_only] [--interp_res INTERP_RES]
                                    mzml parsed_report

    Created on Tue Oct 20 15:29:11 2020 @author: joeyd
    
    positional arguments:
      mzml                  the relative path to mzml file to be analyzed.
      parsed_report         the parsed report file generated by pysodist
                            parse_input.py

    optional arguments:
      -h, --help            show this help message and exit
      --labeling LABELING   The labeling scheme used for the highest mass isotope
                            envelope you expect to fit. E.g. N15 or C13.
      --interp_only         Only save the interpolated spectra instead of the raw
                            spectra. Optional, default is N15.
      --sum_only            Only save summed (and interpolated) spectra instead of
                            all individual spectra. Results in 1 spectra per
                            peptide. Optional, default is to extract all spectra.
      --interp_res INTERP_RES
                            Set the interpolation delta m/z - typical values from
                            0.01 to 0.001. Optional, default is 0.001.

Example usage

    $ pysodist extract_spectra ./raw_data/Ace_15N_600_DIA.mzML ./Ace_15N_600_DIA/pd_parsed_report.tsv
    reading mzml file ./raw_data/Ace_15N_600_DIA.mzML...
    total spectra in mzml: 3049
    rt_range in mzml: 0.0022458304 - 120.00355 mins
    total peptides to extract: 11
    peptide rt_range: 30.09 - 58.75 mins
    extracting 20 spectra for peptide 0 : YAGEVSHDDKHIIVDGK
    extracting 16 spectra for peptide 1 : SVVIHAGQDDLGKGDTEESLK
    extracting 32 spectra for peptide 2 : EKDIVGAVLK
    extracting 12 spectra for peptide 3 : LTSLNVVAGSDLRR
    extracting 13 spectra for peptide 4 : LIGPTSVVGR
    extracting 15 spectra for peptide 5 : ccSDVFNQVVK
    extracting 16 spectra for peptide 6 : TNNPETLVALR
    extracting 28 spectra for peptide 7 : ANELLINVK
    extracting 12 spectra for peptide 8 : FEQASESEPTTVSYEIAGNSPNAER
    extracting 13 spectra for peptide 9 : VLGIDGGEGKEELFR
    extracting 18 spectra for peptide 10 : ATDGGAHGVINVSVSEAAIEASTR

### 3. Run the isodist fitting routine on each of the extracted spectra.

Here, we will use the fortran implementation of this fitting routine, but a python-based implementation is expected in upcoming releases. For each spectra, the fitting routine requires 3 primary inputs - the raw spectra to be fit, which was generated by our extract_spectra tool above; a atom_file, which defines the mass and abundance of each atom type and is provided in model_files/atoms.txt; and a model file that defines the amino acid species you would like to define, which is described in detail below.

**Model File**
We will be inspecting U_var387N_998N.txt as an example. This name refers to unlabeled:Species1_variableLabeledNitrogenSpecies2_fixedLabeledNitrogenSpecies3. 

I strongly recommend a labeling scheme similar to this [species1_species2_speciesX.txt] for all of your model files. 

* First choose your number of isotope envelopes (species) to fit - typically 2-3. (3 = nspecies)
* Second, choose your best guess for the relative abundance of the species (U 40; L 40; F 300)
* Third, set the number of atoms types appropriately (5 + number of labeled atoms)
* Fourth, set the initial guesses for the variably labeled atoms, and whether those values are 'fixed' or 'variable'
* Fifth, ensure that each amino acid is defind correctly - the rows correspond to each of your three species (U, L, F), the columns refer to your atom times (C, H, N, O, S, NX, NY)

* Note that any modified amino acids can be specified here - we have included entries for oxidized methionine and modified cysteine, but more can be added. Please contact me directly to implement additional modifications.

run_isodist inputs:

    $ pysodist run_isodist -h
    usage: pysodist run_isodist [-h] [--threads THREADS] [--wait_time WAIT_TIME]
                                [--pysodist_input PYSODIST_INPUT] [--no_cleanup]
                                input_file isodist_command atomfile resfile

    Created on Sat Oct 24 14:07:14 2020 @author: joeyd

    positional arguments:
      input_file            path to the pd_exported_peaks.tsv (produced by
                            extract_spectra.py). This file should be in folder
                            that contains a folder called spectra holding the .tsv
                            of spectra to fit.
      isodist_command       exact fortran command to execute. e.g.
                            ~/software/pysodist/fortran/isodist or
                            C:/isodist_win/isodist_win.exe
      atomfile              Specify the path to the atom definition file (e.g.
                            exp_atom_defs.txt). You will likely not need to modify
                            this file.
      resfile               Specify the path to the residue labeling file - you
                            will likely need to edit this file based on your
                            labeling scheme. Note that your output will use the
                            name of this file to ensure you know which model file
                            produced which output.

    optional arguments:
      -h, --help            show this help message and exit
      --threads THREADS     number of threads to use. typically 1 less than the
                            number of cores available. Default=2
      --wait_time WAIT_TIME
                            number of seconds to wait between each polling to test
                            if the isodist run has finished. Default=120 seconds
      --pysodist_input PYSODIST_INPUT
                            Optional path to the "pd_parsed_report.tsv" file
                            produced by parse_input. Typically a .tsv file that
                            has extracted the relevant info from a Skyline report.
                            Useful in providing additional info such as protein
                            name for each peptide.
      --no_cleanup          Optionally do not clean up the folder my moving around
                            and deleting the intermediate files.

Example usage to fit a spectra derived from cells with unlabeled, pulsed ~39%15N nitrogen, fixed, 99.8% labeled Nitrogen spike. Note that we have initially checked how many threads are available on this machine, and have used ~75% of them.

    $ lscpu
    #check fo the "CPU(s)" field. This can be done in windows by looking at performance monitor

    $ pysodist run_isodist ./Ace_15N_600_DIA/pd_exported_peaks.tsv ~/software/pysodist/fortran/isodist ./model_files/atoms.txt ./model_files/U_var387N_fix998N.txt --wait_time 60 --threads 15
    working in directory: ./Ace_15N_600_DIA
    there are 206 total spectra to fit...
    writing 15 batch files...
    running command: /home/jhdavis/software/pysodist/fortran/isodist batch_0.in in directory: ./Ace_15N_600_DIA/
    running command: /home/jhdavis/software/pysodist/fortran/isodist batch_1.in in directory: ./Ace_15N_600_DIA/
    running command: /home/jhdavis/software/pysodist/fortran/isodist batch_10.in in directory: ./Ace_15N_600_DIA/
    running command: /home/jhdavis/software/pysodist/fortran/isodist batch_11.in in directory: ./Ace_15N_600_DIA/
    running command: /home/jhdavis/software/pysodist/fortran/isodist batch_12.in in directory: ./Ace_15N_600_DIA/
    running command: /home/jhdavis/software/pysodist/fortran/isodist batch_13.in in directory: ./Ace_15N_600_DIA/
    running command: /home/jhdavis/software/pysodist/fortran/isodist batch_14.in in directory: ./Ace_15N_600_DIA/
    running command: /home/jhdavis/software/pysodist/fortran/isodist batch_2.in in directory: ./Ace_15N_600_DIA/
    running command: /home/jhdavis/software/pysodist/fortran/isodist batch_3.in in directory: ./Ace_15N_600_DIA/
    running command: /home/jhdavis/software/pysodist/fortran/isodist batch_4.in in directory: ./Ace_15N_600_DIA/
    running command: /home/jhdavis/software/pysodist/fortran/isodist batch_5.in in directory: ./Ace_15N_600_DIA/
    running command: /home/jhdavis/software/pysodist/fortran/isodist batch_6.in in directory: ./Ace_15N_600_DIA/
    running command: /home/jhdavis/software/pysodist/fortran/isodist batch_7.in in directory: ./Ace_15N_600_DIA/
    running command: /home/jhdavis/software/pysodist/fortran/isodist batch_8.in in directory: ./Ace_15N_600_DIA/
    running command: /home/jhdavis/software/pysodist/fortran/isodist batch_9.in in directory: ./Ace_15N_600_DIA/
    waiting 1 minutes to check on processes...

### 4. Plot fits and generate interactive notebook.

This is completely optional as all of the fitting is complete and one can inspect/analyze the fits directly. I suggest using the provided jupyter notebook as a template, which is in pysodist/pysodist/utils/analysis_template.ipynb. This file can be simply copied to the directory with your fits, and you can use/edit it from there. Alternatively, you can pre-plot many of the relevant plots using this command.

    $ pysodist plot_spectra -h
    usage: pysodist plot_spectra [-h] [--numerator NUMERATOR [NUMERATOR ...]]
                                 [--denominator DENOMINATOR [DENOMINATOR ...]]
                                 [--no_png] [--no_pdf]
                                 input_file fit_folder output_folder

    Created on Thu Oct 22 16:56:18 2020 @author: joeyd
    
    positional arguments:
      input_file            path to the combiled isodist .csv file bearing all of
                            the results (1 row per fit spectra).
      fit_folder            path to the folder containing all of the isodist .fit
                            files.
      output_folder         path to a folder to save all of the fits into. Will be
                            created if it does not exist.

    optional arguments:
      -h, --help            show this help message and exit
      --numerator NUMERATOR [NUMERATOR ...]
                            list of the fields to use in the numerator of the
                            abundance ratio calculation (typically AMP_U, AMP_L,
                            AMP_F, or some combination. Default is AMP_U
      --denominator DENOMINATOR [DENOMINATOR ...]
                            list of the fields to use in the denominator of the
                            abundance ratio calculation (typically AMP_U, AMP_L,
                            AMP_F, or some combination. Default is AMP_U, AMP_F
      --no_png              By default .png files for the plots will be saved.
                            This option forces these to not be saved.
      --no_pdf              By default .pdf files for the plots will be saved.
                            This option forces these to not be saved.

Example usage:

    $ pysodist plot_spectra ./Ace_15N_600_DIA/U_var387N_fix998N_isodist_outputs/U_var387N_fix998N_output.csv ./Ace_15N_600_DIA/U_var387N_fix998N_isodist_fits/ ./Ace_15N_600_DIA/final_plots --numerator AMP_U --denominator AMP_U AMP_F

This will plot all of the spectra and summary plots for the entries in the compiled isodist ouptu (U_var387N_fix998N_output.csv). The plots will be saved in /final_plots and the ratio calculated will be [U/(U+F)]. Note that the interactive notebook will also be copied into the final_fits directory.

### 5. Running the full pipeline

An additional command is provided to run the entire pipeline. This is helpful if you are iterating over many samples with the same parameters you can you trivially produce a bash script to run this command on each of your samples. Note that this full pipeline only analyzes a single sample at a time.

    $ pysodist full_pipeline -h
    usage: pysodist full_pipeline [-h] [--threads THREADS] [--wait_time WAIT_TIME]
                                  [--output_directory OUTPUT_DIRECTORY]
                                  [--protein_list PROTEIN_LIST]
                                  [--isotope ISOTOPE] [--q_value Q_VALUE]
                                  [--labeling LABELING] [--interp_only]    
                                  [--sum_only] [--interp_res INTERP_RES]
                                  [--numerator NUMERATOR [NUMERATOR ...]]
                                  [--denominator DENOMINATOR [DENOMINATOR ...]]
                                  [--no_png] [--no_pdf]
                                  input mzml sample_name isodist_command atomfile
                                  resfile

    Created on Sat Oct 24 14:04:54 2020 @author: joeyd

    positional arguments:
      input                 input file to parse. Currently only skyline report
                            files are supported
      mzml                  the relative path to mzml file to be analyzed. For
                            Thermo instruments, one should generate the .mzml file
                            from the original .raw file using msconvert as
                            follows: .\msconvert.exe ".\[YOUR_RAW_FILE].raw" -o
                            "./" --mzML --64 -v mz64 --inten32 --noindex --filter
                            "msLevel 1" --zlib
      sample_name           name of the sample within the skyline report to be
                            analyzed.
      isodist_command       exact fortran command to execute. e.g.
                            C:\isodist\isodist.exe
      atomfile              Specify the path to the atom definition file (e.g.
                            exp_atom_defs.txt). You will likely not need to modify
                            this file.
      resfile               Specify the path to the residue labeling file - you
                            will likely need to edit this file based on your
                            labeling scheme. Note that your output will use the
                            name of this file to ensure you know which model file
                            produced which output.

    optional arguments:
      -h, --help            show this help message and exit
      --threads THREADS     number of threads to use. typically 1 less than the
                            number of cores available. Default=4
      --wait_time WAIT_TIME
                            if the isodist run has finished. Default=60 seconds
      --output_directory OUTPUT_DIRECTORY
                            Output files will be saved in this folder: 1 directory
                            per sample in the skyline report. Default = ./
      --protein_list PROTEIN_LIST
                            An optional list of the proteins to parse. By default,
                            all proteins in the report are analyzed.
      --isotope ISOTOPE     Be default, it is assumed that the report contains a
                            light isotope (no special labeling), if this field is
                            not present in the report, you can specify a different
                            field here (e.g. "heavy")
      --q_value Q_VALUE     Used to optionally filter the report file based on the
                            q_value. By default, no q_value filtering is used.
      --labeling LABELING   The labeling scheme used for the highest mass isotope
                            envelope you expect to fit. E.g. N15 or C13
      --interp_only         Only save the interpolated spectra instead of the raw
                            spectra
      --sum_only            Only save summed (and interpolated) spectra instead of
                            all individual spectra. Results in 1 spectra per
                            peptide
      --interp_res INTERP_RES
                            Set the interpolation delta m/z - typical values from
                            0.01 to 0.001
      --numerator NUMERATOR [NUMERATOR ...]
                            list of the fields to use in the numerator of the
                            abundance ratio calculation (typically AMP_U, AMP_L,
                            AMP_F, or some combination. Default is AMP_U.
      --denominator DENOMINATOR [DENOMINATOR ...]
                            list of the fields to use in the denominator of the
                            abundance ratio calculation (typically AMP_U, AMP_L,
                            AMP_F, or some combination. Default is AMP_U, AMP_F
      --no_png              By default .png files for the plots will be saved.
                            This option forces these to not be saved.
      --no_pdf              By default .pdf files for the plots will be saved.
                            This option forces these to not be saved.

Example usage:
    
    $ pysodist full_pipeline ./raw_data/ace_15N_pysodist_v002.csv ./raw_data/Ace_15N_600_DIA.mzML Ace_15N_600_DIA ~/software/pysodist/fortran/isodist ./model_files/atoms.txt ./model_files/U_var387N_fix998N.txt --threads 12 --wait_time 120 --no_pdf

This will work on all of the spectra corresponding to the Ace_15N_600_DIA sample that were specified in the ace_15N_pysodist_v002.csv skyline report file. It will use 12 threads for the isodist fitting, checking the threads every 2 minutes, and will not produce any .pdf plots.


## Contact

Please contact me directly with any bugs, feature requests or general feedback at jhdavis[at]mit[dot]edu.
