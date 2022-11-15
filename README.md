# Tutorial to accompany Fink et al. (2022)
This repository contains code and toy data associated with the paper:

Fink, L., Simola, J., Tavano, A., Lange, E. B., Wallot, S., & Laeng, B. (2021). From pre-processing to advanced dynamic modeling of pupil data. https://doi.org/10.31234/osf.io/wqvue


## Contents
`Fink_pupilCodeTutorial_interactive.mlx` is the main tutorial script. It contains text, MATLAB code, and interactive user interfaces for understanding some fundamentals of time series and signal-to-signal analyses.

`Fink_toydata.mat` contains the raw data we will be working with during the workshop.

All other scripts are helper functions, called in the main interactive script:

`genPRF.m` is a function to generate a pupillary response function, called in the interactive script.  

`getFFT.m` is a function to calculate the Fast Fourier Transform of an input signal.

`samps2secs.m` is a function to convert samples to seconds (useful for plotting).

`fdeconc.m` is a function to compute fast deconvolution.
[https://www.mathworks.com/matlabcentral/fileexchange/5465-fast-deconvolution](https://www.mathworks.com/matlabcentral/fileexchange/5465-fast-deconvolution)

`movcorr.m` is a function to compute a moving window correlation.
[https://www.mathworks.com/matlabcentral/fileexchange/65342-movcorr-x-y-k-varargin-compute-windowed-correlation-coefficient](https://www.mathworks.com/matlabcentral/fileexchange/65342-movcorr-x-y-k-varargin-compute-windowed-correlation-coefficient)

`DFA_fun.m` is a function to compute a detrended fluctuation analysis.
[https://www.mathworks.com/matlabcentral/fileexchange/67889-detrended-fluctuation-analysis-dfa](https://www.mathworks.com/matlabcentral/fileexchange/67889-detrended-fluctuation-analysis-dfa)

`RQA` is a folder containing a .csv of the toy data and an R script to conduct RQA and CRQA. This code is in the process of being converted to MATLAB and integrated with the main .mlx script. We apologize for the inconvenience. Expect an update by 18.11.2022

## Dependencies
This code was written in MATLAB version 2021b and relies on the following:
- Signal Processing Toolbox
- Statistics Toolbox
- Symbolic Toolbox

## Running this code
Download this Github repository. Open `Fink_pupilCodeTutorial_interactive.mlx`.

As long as all files have stayed in the same folder, MATLAB should find the required data and functions and all should run without problem.

Should trouble occur, or bugs be found, please [open an issue](https://github.com/lkfink/pupilTutorial/issues).

## Citation
If using anything from this repository, cite:

Fink, L., Simola, J., Tavano, A., Lange, E. B., Wallot, S., & Laeng, B. (2021, December 2). From pre-processing to advanced dynamic modeling of pupil data. https://doi.org/10.31234/osf.io/wqvue

See also:

Fink, L. K., Hurley, B. K., Geng, J. J., & Janata, P. (2018). A linear oscillator model predicts dynamic temporal attention and pupillary entrainment to rhythmic patterns. Journal of Eye Movement Research, 11(2). https://doi.org/10.16910/jemr.11.2.12
