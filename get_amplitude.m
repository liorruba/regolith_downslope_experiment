function [mean_max_peak, min_peak] = get_max_peak(time, signal, mean_peak_width)

p = findpeaks(time,signal,'MinPeakWidth',mean_peak_width);
mean_max_peak = mean(p);

p = findpeaks(time,-signal,'MinPeakWidth',mean_peak_width);
mean_min_peak = mean(p);