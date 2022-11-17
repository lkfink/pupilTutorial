% Get fft of cleaned data 
% LF 

function [FFTresult, power, fVals] = getFFT(data, fs)

L=length(data);                                         % length of input data
exp = nextpow2(L);                                      % find exp for next closest power of 2
NFFT = 2^exp;                                           % set nfft to that value
window_hann_1 = hanning(size(data,2));                  % create window
data_1 = data'.*repmat(window_hann_1,1,size(data,1));   % window data
data_1 = data_1* 2/sqrt(1.5);                           % correction for hanning
FFTresult = fft(data_1',NFFT,2)/NFFT;                   % do fft
fVals = fs/2*linspace(0,1,NFFT/2+1);                    % get frequency values  
power = 2*(abs(FFTresult));                             % get power of each positive freq component      
power = power(1:NFFT/2+1);                              % only keep relevant data

end
