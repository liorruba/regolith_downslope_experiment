function [cohesion] = find_cohesion(time, accel_data)
% Finds the cohesion, a constant for which the total acceleration vector becomes
% zero after the slope stabilizes

mean_max_peak = get_max_peak(time, accel_data(time > 50));

plot(time,accel_data );

cohesion = mean_max_peak;
end

