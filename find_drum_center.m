function [ctr_coordinates] = find_drum_center(current_frame)
    % Try to find green sticker
    idx = current_frame(:,:,1) >= 130 & current_frame(:,:,1) < 190 &...
                    current_frame(:,:,2) >= 210 & current_frame(:,:,2) < 250 &...
                    current_frame(:,:,3) >= 160 & current_frame(:,:,3) < 210;
    % Else use green color of handle            
    if (nnz(idx(:)) < 10)
        current_frame(1:100,:) = 0; % There is a white sticker near the upper edge of the drum, so remove it
        current_frame(:,1:200) = 0; % There is a white sticker near the upper edge of the drum, so remove it
        current_frame = current_frame(:,:,1) >= 250 & current_frame(:,:,2) >= 250 & current_frame(:,:,3) >= 240;
    else
        current_frame = idx;
    end
    
    current_frame = bwareaopen(current_frame, 30);
    ctr_coordinates = regionprops(current_frame).Centroid;
    
end