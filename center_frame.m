function [curr_frame] = center_frame(curr_frame, drum_center)
    frame_center = size(curr_frame) / 2;
    
    % Adjust x axis:
    x_shift = fix(drum_center(1) - frame_center(2));
    if x_shift > 0
        curr_frame(:, 1:x_shift, :) = [];
        curr_frame = padarray(curr_frame,[0 x_shift],'post');
    else
        x_shift = abs(x_shift);
        curr_frame(:, end-x_shift:end, :) = [];
        curr_frame = padarray(curr_frame,[0 x_shift],'pre');
    end
    % Adjust y axis:
    y_shift = fix(drum_center(2) - frame_center(1));
    if y_shift > 0
        curr_frame(1:y_shift, :, :) = [];
        curr_frame = padarray(curr_frame,[y_shift 0],'post');
    else
        y_shift = abs(y_shift);
        curr_frame(end-y_shift:end, :, :) = [];
        curr_frame = padarray(curr_frame,[y_shift 0],'pre');
    end
end