# Tutorial to accompany Fink et al. (2023)
This repository contains code and toy data associated with the paper:

Fink, L., Simola, J., Tavano, A., Lange, E. B., Wallot, S., & Laeng, B. (2023, in review). From pre-processing to advanced dynamic modeling of pupil data. https://doi.org/10.31234/osf.io/wqvue

Please note that for users who do not have MATLAB, all code can be viewed in the .html file of this repository.  
Alternatively, users can run all MATLAB code on a virtual machine in a public Code Ocean capsule, available here: 


## Contents

### Main interactive tutorial script
`Fink_pupilCodeTutorial_interactive.mlx` is the main tutorial script. It contains text, MATLAB code, and interactive user interfaces for understanding some fundamentals of time series and signal-to-signal analyses.

`Fink_pupilCodeTutorial_interactive.html` is an html print out of the code and figures in the above .mlx script, for those who want to browse online, or who do not have MATLAB. Note that the interactive controls are not available in the html version.

`Fink_toydata.mat` contains the toy data set for the tutorial.

### Helper functions
The `helperFunctions` folder contains custom and third-party functions necessary for running the interactive tutorial. 

- `genPRF.m` is a function to generate a pupillary response function, called in the interactive script.  

- `getFFT.m` is a function to calculate the Fast Fourier Transform of an input signal.

- `samps2secs.m` is a function to convert samples to seconds (useful for plotting).

- `fdeconc.m` is a function to compute fast deconvolution.
[https://www.mathworks.com/matlabcentral/fileexchange/5465-fast-deconvolution](https://www.mathworks.com/matlabcentral/fileexchange/5465-fast-deconvolution)

- `movcorr.m` is a function to compute a moving window correlation.
[https://www.mathworks.com/matlabcentral/fileexchange/65342-movcorr-x-y-k-varargin-compute-windowed-correlation-coefficient](https://www.mathworks.com/matlabcentral/fileexchange/65342-movcorr-x-y-k-varargin-compute-windowed-correlation-coefficient)

- `MDRQA.m` is a function for recurrence quantification analysis 

- `mdcrqa.m` is a function for cross-recurrence quantification analysis 

- `dcrp.m` is a function for computing a cross-recurrence profile. 

The latter three functions were written by Sebastian Wallot. Please see the citation in the code for proper reference.  


## Dependencies
This code was written in MATLAB verion 9.11.0.1769968 (R2021b) and relies on the following toolboxes:
- Signal Processing Toolbox: Version 8.7
- Statistics and Machine Learning Toolbox: Version 12.2
- Symbolic Math Toolbox: Version 9.0

Development machine details: 
Operating System: macOS  Version: 11.6.6 Build: 20G624 
Java Version: Java 1.8.0_202-b08 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode


## Running this code
Download this Github repository. Open `Fink_pupilCodeTutorial_interactive.mlx`.

As long as all the file structure and names in this repository have not been altered, MATLAB should find the required data and functions and all should run without problem.

Should trouble occur, or bugs be found, please [open an issue](https://github.com/lkfink/pupilTutorial/issues).

## Citation
If using anything from this repository, cite:

Fink, L., Simola, J., Tavano, A., Lange, E. B., Wallot, S., & Laeng, B. (2023, in press). From pre-processing to advanced dynamic modeling of pupil data. https://doi.org/10.31234/osf.io/wqvue
