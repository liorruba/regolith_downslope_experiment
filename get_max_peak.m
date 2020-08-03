function [max_peak] = get_max_peak(time, signal,fun)

dt = mean(diff(time));
sampling_frequency = 1./dt;

p = findpeaks(signal,'MinPeakWidth',(sampling_frequency / 6)./4);
max_peak = fun(p);