function [binframe] = make_binary_frame(curr_frame, regolith_mean_rgb, regolith_std_rgb, rg_std_scale, spots_threshold)

    % Set default value for the spot size threshold:
    if nargin == 4
       spots_threshold = 2000;
    end

    regolith_std_rgb = regolith_std_rgb * rg_std_scale;

    currframe_r = curr_frame(:,:,1);
    currframe_r(currframe_r < (regolith_mean_rgb(1) - regolith_std_rgb(1))) = 0;
    currframe_r(currframe_r > (regolith_mean_rgb(1) + regolith_std_rgb(1))) = 0;

    currframe_g = curr_frame(:,:,2);
    currframe_g(currframe_g < (regolith_mean_rgb(2) - regolith_std_rgb(2))) = 0;
    currframe_g(currframe_g > (regolith_mean_rgb(2) + regolith_std_rgb(2))) = 0;

    currframe_b = curr_frame(:,:,3);
    currframe_b(currframe_b < (regolith_mean_rgb(3) - regolith_std_rgb(3))) = 0;
    currframe_b(currframe_b > (regolith_mean_rgb(3) + regolith_std_rgb(3))) = 0;

    binframe = currframe_r & currframe_g & currframe_b;

    % Remove unecessary parts of the frame:
%     binframe(exclusion_bounding_box(1):exclusion_bounding_box(2),exclusion_bounding_box(3):exclusion_bounding_box(4)) = 0;
   
    % Remove spots smaller than spots_threshold pixels (the drum was dirty...):
    binframe = bwareaopen(binframe,spots_threshold);
end

