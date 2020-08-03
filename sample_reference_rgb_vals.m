function [regolith_mean_rgb, regolith_std_rgb] = sample_reference_rgb_vals(ref_img, theta, tilt, sample_box)
    if nargin == 3
        sample_box = [130 160 320 350];
    end

%     ref_img = imread(['ref_frames/reference_frame_',num2str(theta),'_', num2str(tilt),'.png']);
    
    ref_img_double = double(ref_img(sample_box(1):sample_box(2),sample_box(3):sample_box(4),:));

    rmean = mean(mean(ref_img_double(:,:,1)));
    gmean = mean(mean(ref_img_double(:,:,2)));
    bmean = mean(mean(ref_img_double(:,:,3)));

    rstd = std(std(ref_img_double(:,:,1)));
    gstd = std(std(ref_img_double(:,:,2)));
    bstd = std(std(ref_img_double(:,:,3)));

    regolith_mean_rgb = [rmean gmean bmean];
    regolith_std_rgb = [rstd gstd bstd];
    
    imagesc(ref_img)
    rectangle('Position', [sample_box(3) sample_box(1) sample_box(2)-sample_box(1) sample_box(4)-sample_box(3)])
end