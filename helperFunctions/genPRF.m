function PRF = genPRF(Fs, motor)
% Compute Hooks & Levelt PRF, with or without McCloy adaptation

% If using this function, please cite:
% Fink, L. K., Hurley, B. K., Geng, J. J. & Janata, P. (2018). 
% A linear oscillator model predicts dynamic temporal attention and pupillary
% entrainment to rhythmic patterns. Journal of Eye Movement Research, 
% 11(2):12. DOI: 10.16910/jemr.11.2.12

% Lauren Fink (lauren.fink@ae.mpg.de)
% Max Planck Institute for Empirical Aesthetics 
% Previously: University of California, Davis 

% Input:
% - Fs: sampling frequency in Hertz, e.g. 500
% - motor: boolean specifying whether to use PRF based on motor response
% (key press) or not (no key press)

% Output:
% - Pupillary Response Function

% Empirically established constants: -------------------------------------%

% Latency of pupil response maximum
% Different depending on whether motor response involved (i.e., button
% press)
if motor
    tmax = 930; % Hoeks & Levelt
    tlim = 2500;

else
    tmax = 512; % McCloy et al. 2016
    tlim = 1300; % 1300 % See McCloy Fig. 1a
end


% Shape parameter of Erlang gamma distribution
% (proposed to be the number of signaling steps in neural pathway 
% transmitting attentional pulse to pupil)
n = 10.1; 

% Time
tscale = 0:1000/Fs:tlim;

% Pupillary Response Function  -------------------------------------------%
PRF = tscale.^n .* exp(1).^(-10.1*tscale/tmax);
% NOTE would need to divide by scaling factor if want to compare motor
% and non-motor
 
% NOTE end of function approaches infinity and does not return to
% zero. That may be a problem for some operations. But setting last point
% to 0 would introduce a small bump. 
% PRF(end) = 0;

end