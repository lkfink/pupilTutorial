function xscale = samps2secs(data, fs)

% rescale xaxis from samples to seconds.
% input:
% - data for which we want new xscale
% - sampling rate of data

% output:
% - continuous array of time points in milliseconds, corresponding to each
% sample  (useful for plotting)

% LKF 2022

xscale = 1:numel(data);
xscale = xscale*(1000/fs); % convert to msecs
xscale = xscale / 1000; 

end