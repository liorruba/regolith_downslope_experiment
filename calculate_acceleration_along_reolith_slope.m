function [ax, ay, az] = calculate_acceleration_along_reolith_slope(ax, ay, az, regolith_slope, regolith_slope_time)
    for aa = 1:length(ax)
        angle_to_rot_by = regolith_slope
       [x_comp, y_comp, z_comp] = roty_deg(alpha, x_comp, y_comp, z_comp) 
    end

end