clear

%% Load video
% Set these parameters to change the loaded video file
theta = 0; % deviation from repose
tilt = 0; % the tilt angle of the shaker
run = 3; % trial #
t0 = 6; % sec, when to start video

% Read video
file = ['./videos/theta_',num2str(theta),'_tilt_',num2str(tilt),'_',num2str(run),'.m4v'];
v = VideoReader(file);

% Linear fit parameters
debug_fit = false;
debug_image_rec = true;

% Read reference regolith RGB colors:
% Get reference frame (image) as the first interesting frame:
while v.CurrentTime < 350
    if v.CurrentTime > t0
        ref_frame = v.readFrame();
        imagesc(ref_frame);
        break;
    end
    v.readFrame();
end
drum_center = find_drum_center(ref_frame);
ref_frame = center_frame(ref_frame, drum_center);

%% For each trial, the different illumination conditions require a different
% sensitivity and spots_threshold.
% sensitivity: the factor multiplying the standard deviation of the rgb
%              levels (see function make_binary_frame.m)
% spots_threshold: the size of the smallest cluster of binary pixels 
%                  to be removed from the binary frame constrained by the
%                  correct rgb range
switch tilt
    case 0
        sensitivity = 50;
        spots_threshold = 500;
        [regolith_mean_rgb, regolith_std_rgb] = sample_reference_rgb_vals(ref_frame, theta, tilt,[170 220 270 310]);
        mask_cols = 345:375;
        
    case 5
        sensitivity = 15;
        spots_threshold = 500;
        [regolith_mean_rgb, regolith_std_rgb] = sample_reference_rgb_vals(ref_frame, theta, tilt,[170 220 280 320]);
        
        mask_cols = 345:375;
        
        if theta == 5 || theta == 10 || theta == 15
           [regolith_mean_rgb, regolith_std_rgb] = sample_reference_rgb_vals(ref_frame, theta, tilt, [150 180 320 350]);
        end
        
    case 10
        sensitivity = 22;
        spots_threshold = 500;
        [regolith_mean_rgb, regolith_std_rgb] = sample_reference_rgb_vals(ref_frame, theta, tilt, [160 210 290 330]);
        mask_cols = [160:185 335:375];
        
    case 15
        sensitivity = 40;
        spots_threshold = 500;  
        mask_cols = [220:240 345:375];
        [regolith_mean_rgb, regolith_std_rgb] = sample_reference_rgb_vals(ref_frame, theta, tilt,[150 200 290 330]);
        
        if theta == 5 || theta == 10 || theta == 15
           sample_box = [150 180 320 350]; 
           [regolith_mean_rgb, regolith_std_rgb] = sample_reference_rgb_vals(ref_frame, theta, tilt, sample_box);
        end
        
    case 20
        sensitivity = 20;
        spots_threshold = 1000;
        
        mask_cols = [230:250 330:375];
        [regolith_mean_rgb, regolith_std_rgb] = sample_reference_rgb_vals(ref_frame, theta, tilt,[140 180 300 340]);
        
    case 25
        sensitivity = 25;
        spots_threshold = 500; 
        
        mask_cols = [130:175 230:260 355:380];
        [regolith_mean_rgb, regolith_std_rgb] = sample_reference_rgb_vals(ref_frame, theta, tilt,[140 180 300 340]);
        
        if theta == 5 || theta == 10 || theta == 15
            sample_box = [148 180 320 350];
            [regolith_mean_rgb, regolith_std_rgb] = sample_reference_rgb_vals(ref_frame, theta, tilt, sample_box);
        end
        
    case 30
        mask_cols = [220:260 350:380];
        [regolith_mean_rgb, regolith_std_rgb] = sample_reference_rgb_vals(ref_frame, theta, tilt,[145 180 300 340]);
        
        sensitivity = 23;
        spots_threshold = 2000;
        
    case 35
        mask_cols = [230:275 334:380];
        [regolith_mean_rgb, regolith_std_rgb] = sample_reference_rgb_vals(ref_frame, theta, ...
            tilt,[206 243 215 258]);
        
        sensitivity = 5;
        spots_threshold = 2000;
    
    case 40
        mask_cols = [100:110 230:285 350:410];
        [regolith_mean_rgb, regolith_std_rgb] = sample_reference_rgb_vals(ref_frame, theta, ...
            tilt,[200 240 215 258]);
        
        sensitivity = 25;
        spots_threshold = 1500;

    case 45
        mask_cols = [240:295 350:410];
        [regolith_mean_rgb, regolith_std_rgb] = sample_reference_rgb_vals(ref_frame, theta, ...
            tilt,[205 250 220 265]);
        
        sensitivity = 20;
        spots_threshold = 1500;
end

%% 
ii = 1;
jj = 1;

if debug_image_rec
%     close all
%     figure('position',[600 550 1000 250])
    clf
end

% Read all frames and place in a 4-D array:
while v.CurrentTime < 350
    if mod(ii,9) == 0 && v.CurrentTime > t0
        curr_frame = v.readFrame();
        f(:,:,:,jj) = curr_frame;
        
        drum_center = find_drum_center(curr_frame);
        curr_frame = center_frame(curr_frame, drum_center);

        binframe = make_binary_frame(curr_frame, regolith_mean_rgb, regolith_std_rgb, sensitivity, spots_threshold);
        binframe(:, mask_cols) = 0;
        
        % Fill small holes:
        binframe = imfill(binframe, 4, 'holes');
%         imagesc(curr_frame)
%         drawnow 
%         continue

        % Detect regolith surface:
        % Model function is a Gaussian:
        [h,~] = size(binframe);
        x = 1:h;
%         model_function = @(beta,x) exp(-(x-beta(1)).^2 ./ beta(2));
        model_function = @(x,b) 1 - 1./(1+exp(-100*(x-b)));
        % Iterate over all pixels columns and fit the model function:
        clear top_edge;
        vec_iterate_over = find(any(binframe ~= 0, 1)); % any() is slightly faster than all()
        for kk = 1:length(vec_iterate_over)
            % Try to remove stray pixels
            xx = binframe(:,vec_iterate_over(kk));
            xx = bwareaopen(xx, 20);
            init_guess_step_loc = find(xx,1);
            
            % If removed to much, don't remove at all:
            if isempty(init_guess_step_loc)
                xx = binframe(:,vec_iterate_over(kk));
                init_guess_step_loc = find(xx,1);
            end
            try
                b = nlinfit(x,xx',model_function,init_guess_step_loc);
            catch
                b = [nan nan];
            end
            
            if debug_fit
                clf;
                plot(binframe(:,vec_iterate_over(kk)),'linewidth',2);
                axis equal
                hold on
                plot(x,model_function(b, x),'--','linewidth',2);
                axis square
                drawnow
            end
            
            top_edge(kk) = b;
        end
        
        x = vec_iterate_over; y = top_edge;
        nanidx = isnan(x) | isnan(y);
        x(nanidx) = [];
        y(nanidx) = [];
        p = polyfit(x, y, 1);
        
        top_edge_slope(jj) = p(1);
        
        % Remove outliers
        res = polyval(p,x)-y;
        x = x(abs(res) < std(res));
        y = y(abs(res) < std(res));
        p = polyfit(x, y, 1);
        y = filloutliers(y, 'center', 'movmedian',5,'ThresholdFactor',5);
        curr_time(jj) = v.CurrentTime;
        
        % Calculate the slope of both parts of the regolith:
        % Top
        y_top = y(x > 250);
        x_top = x(x > 250);
        p_top = polyfit(x_top,y_top,1);
        top_slope(jj) = p_top(1);
        
        % Bottom
        y_bottom = y(x < 250);
        x_bottom = x(x < 250);
        p_bottom = polyfit(x_bottom,y_bottom,1);
        bottom_slope(jj) = p_bottom(1);
        
        if debug_image_rec
            x = vec_iterate_over; y = top_edge;
            nanidx = isnan(x) | isnan(y);
            x(nanidx) = [];
            y(nanidx) = [];
            
            subplot(1,2,1,'replace');
            title('Binary frame');
            imagesc(binframe); hold on
            plot(x, y, 'or');
            plot(x, polyval(p_top,x), 'w','LineWidth',2);
            plot(x, polyval(p_bottom,x), 'w','LineWidth',2);
            set(gca,'XTickLabel',[],'YTickLabel',[])
            drawnow;
            
            subplot(1,2,2,'replace');
            plot(curr_time, abs(top_slope) ,'o', curr_time, abs(bottom_slope),'o',curr_time, abs(top_edge_slope),'ko');
            ylim([ 0 1])
            drawnow;
        end
        jj = jj + 1;
        
    end
    v.readFrame();
    ii = ii + 1;
end
%%
disp('Done');

save(['results/theta_',num2str(theta),'_tilt_',num2str(tilt),'_run_',num2str(run),'.mat'], ...
    'top_edge_slope', 'curr_time','top_slope','bottom_slope');
