%% Varying the shaker tilt
clear
tilt_vec = 0:5:45;
theta_vec = 0;
run_vec = 1;
clf

model_fun = @(b,x) b(1) + b(2) * cos(b(3) .* x + b(4));
cmap1 = viridis(length(tilt_vec));
cmap2 = inferno(length(tilt_vec));
co = get(gca,'colororder');

%% Regolith slope model function
regolith_slope_model_function = @(beta, c, x) ...
    ((x < beta(1)) + (x > beta(1)) .* exp(beta(2) .* (x - beta(1))) + beta(3)) .* c;

%%
for ii = 1:length(theta_vec)
    for jj = 1:length(tilt_vec)
        for kk = 1:length(run_vec)
            %% Loads slope difference data
            load(['results/theta_',num2str(theta_vec(ii)),'_tilt_',num2str(tilt_vec(jj)),'_run_',num2str(run_vec(kk)),'.mat']);
            
            %% Load acceleration data
            acc_data = csvread(['/Users/liorr/Work/Projects/regolith_downslope_experiment/data/acceleration/theta_',...
                num2str(theta_vec(ii)),'_tilt_',num2str(tilt_vec(jj)),'_',num2str(run_vec(kk)),'.txt'],50,0);
            
            %% Unpack variables:
            t = curr_time - curr_time(1); % Time measured by video camera
            regolith_slope = abs(bottom_slope); % Slope measured by algorithm
            
            A_time = acc_data(1:end-100,1)/1000; % Time measured by accelerometer (s)
            A_x = acc_data(1:end-100,2) / 1e2; % acc in x axis (accelerometer frame of reference) (m s^-2)
            A_y = acc_data(1:end-100,3) / 1e2; % acc in y axis (accelerometer frame of reference) (m s^-2)
            A_z = acc_data(1:end-100,4) / 1e2; % acc in z axis (accelerometer frame of reference) (m s^-2)
            
            %% Set the time range to sample from
            % Define At0 as the time in which acceleration increases by
            % 130%:
            A_t0 = A_time(find(A_z > 13,1)); % sec
            A_time = A_time - A_t0;
            
            %% LSQ fit to the decrease in slope
            % Remove outliers and smooth for fitting:
            regolith_slope = filloutliers(regolith_slope,'center','movmedian',3);
            smooth_slope = movmedian(regolith_slope,50);
            dslope = gradient(smooth_slope);
            
            % Initial guess for t0 as the time in which collapse starts, using the
            % derivarive:
            [~, idx] = min(dslope(t < 50));
            t0_guess = t(idx);
            
            % Static coefficient of friction
            mu_s = mean(tan(regolith_slope(t < t0_guess)));
            
            % Initial guess for the differnece between the initial and
            % final slopes:
            intial_slope(ii,jj,kk) = mean(smooth_slope(t < t0_guess));
            final_slope(ii,jj,kk) = mean(smooth_slope(end-50:end));
            
            slope_diff_init(ii,jj,kk) = intial_slope(ii,jj,kk) - final_slope(ii,jj,kk);
            slope_diff_init_deg(ii,jj,kk) = atand(slope_diff_init(ii,jj,kk));
            
            % Non linear regression
            options = optimset('FunValCheck','off','tolfun',1e-10,'MaxFunEvals',10000);
            init_guess = [t0_guess -0.1 1];
            
            
            beta_fit_nlinfit = nlinfit(t, smooth_slope, ...
                @(beta_fit, t) regolith_slope_model_function(beta_fit, slope_diff_init(ii,jj,kk), t), init_guess ,options);
            
            beta_fit_lsq = lsqcurvefit(@(beta_fit, t) regolith_slope_model_function(beta_fit, slope_diff_init(ii,jj,kk), t), init_guess, ...
                t, smooth_slope,[t0_guess/2 -1 -1],[],options);
                        
            %% Synchronize acceleration and slope vetors
            % Choose which fit parameters to use, lsq or nlinfit:
            %             beta_fit = beta_fit_nlinfit;
            beta_fit = beta_fit_lsq;
            
            % Plot regolith angle and model fit
%                                     plot(t, regolith_slope,'color',co(1,:));
%                                     hold on
%                                     plot(t,regolith_slope_model_function(beta_fit, slope_diff_init, t),'-','color',co(2,:),'linewidth',1)
%                                     hold off
%                                     continue;
            % Get t0 from the fit:
            t0_fit = beta_fit(1);
            if ~isempty(t(find(t > t0_fit, 1)))
                t = t - t(find(t > t0_fit, 1));
            else
                % do nothing
            end
            
            % Chop-off everything before t0
            regolith_slope(t < 0) = [];
            t(t < 0) = [];
            
            A_time = A_time - A_t0;
            
            % Make all timeseries the same length
            tf = min(A_time(end), t(end)); % seconds
            regolith_slope(t >= tf) = [];
            t(t >= tf) = [];
            A_x(A_time <= 0 | A_time >= tf) = [];
            A_y(A_time <= 0 | A_time >= tf) = [];
            A_z(A_time <= 0 | A_time >= tf) = [];
            A_time(A_time <= 0 | A_time >= tf) = [];
            
            %% Project the acceleration in the regolith slope RF
            % The tilt angle of the shaker, as measured by the
            % non-vibrating component of the accelerometer
            initial_tilt_angle(ii,jj,kk) = atand(mean(A_x)./mean(A_z)); % degrees
            
            g = mean([A_x A_y A_z],1);
            
            % Check g is correct by rotating it to the lab frame of reference
            [gx,gy,gz] = (roty_deg(-initial_tilt_angle(ii,jj,kk), g(1),g(2),g(3)));
            
            A_tot = [A_x, A_y, A_z];
            A_only_shaker = A_tot - g;
            
            [A_x_lab_frame, A_y_lab_frame, A_z_lab_frame] = roty_deg(-initial_tilt_angle(ii,jj,kk), A_x, A_y, A_z);
            [A_s_x_lab_frame, A_s_y_lab_frame, A_s_z_lab_frame] = roty_deg(-initial_tilt_angle(ii,jj,kk), A_only_shaker(:, 1), A_only_shaker(:, 2), A_only_shaker(:, 3));
            
            % To calculate the acceleration along the regolith slope, we
            % need to rotate by the slope angle (varies in time)
            
            % Due to bugs in the serial connection, sometimes the time
            % vector is not monotonically increasing. Force it to be:
            A_time = linspace(A_time(1), A_time(end), length(A_time));
            clear A_s_x_regolith_frame A_s_y_regolith_frame A_s_z_regolith_frame;
            clear A_tot_x_regolith_frame A_tot_y_regolith_frame A_tot_z_regolith_frame;
            
            regolith_slope_angle = interp1(t, atand(regolith_slope), A_time);

            for aa = 1:length(A_time)
                [shaker_x, shaker_y, shaker_z] = roty_deg(regolith_slope_angle(aa) - initial_tilt_angle(ii,jj,kk), A_only_shaker(aa, 1), A_only_shaker(aa, 2), A_only_shaker(aa, 3));
                A_s_x_regolith_frame(aa) = shaker_x;
                A_s_y_regolith_frame(aa) = shaker_y;
                A_s_z_regolith_frame(aa) = shaker_z;
                
                [tot_x, tot_y, tot_z] = roty_deg(regolith_slope_angle(aa) - initial_tilt_angle(ii,jj,kk), A_tot(aa, 1), A_tot(aa, 2), A_tot(aa, 3));
                A_tot_x_regolith_frame(aa) = tot_x;
                A_tot_y_regolith_frame(aa) = tot_y;
                A_tot_z_regolith_frame(aa) = tot_z;
            end
            
            %% Record slope difference from fit:
            %             slope_diff(jj, kk) = beta_fit(4);
            slope_diff(ii, jj, kk) = slope_diff_init(ii,jj,kk);
            
            %% Calculate e-folding time from fit:
            e_folding_time(ii,jj,kk) = -1./beta_fit(2);
            
            %% Slope stability
            c = prctile(abs(A_tot_x_regolith_frame(A_time > e_folding_time(ii,jj,kk)))...
                - mu_s .* A_tot_z_regolith_frame(A_time > e_folding_time(ii,jj,kk)), 98);
            if e_folding_time(ii,jj,kk) > A_time(end)
                c = 0;
            end
%             c = 2;
%             plot(abs(A_tot_x_regolith_frame) - mu_s .* A_tot_z_regolith_frame)
%             continue
            FS = (mu_s .* A_tot_z_regolith_frame + c)./abs(A_tot_x_regolith_frame);
            [xx,yy] = binData1D(A_time, FS, 0:5:A_time(end), @(x) prctile(x,2));
            
            subplot(121)
            plot(A_time, FS, '.', xx, yy, '-r','linewidth',2);
            ylim([0 3])
            subplot(122)
            plot(initial_tilt_angle(ii,jj,kk), c,'ko'); hold on
%             title([initial_tilt_angle(ii,jj,kk) slope_diff_init_deg(ii,jj,kk)])
            
            continue
            
            %% Calculate the total parallel acceleration
            a_tot1 = norm(g) .* sind(regolith_slope_angle) + A_s_x_regolith_frame ...
                - mu_s .* (norm(g) .* cosd(regolith_slope_angle) + A_s_z_regolith_frame);
            cr = 0;
            phi = atand(mu_s);
            
            cpd = cosd(regolith_slope_angle) + sind(regolith_slope_angle) .* mu_s;
            kh = A_s_x_lab_frame ./ norm(g);
            
            
            buff = phi - regolith_slope_angle;
            ky = cr * cosd(phi) ./ cosd(phi - regolith_slope_angle) + tand(phi - regolith_slope_angle);
            kdyv = ky + A_s_z_lab_frame' ./ norm(g) .* tand(phi - regolith_slope_angle);
            
            a_tot2 = cpd .* (kh' - kdyv) .* norm(g);
            
            peaks_before_collapse(ii,jj,kk) = prctile(a_tot1(A_time < t0_fit), 98);
            peaks_after_collapse(ii,jj,kk) = prctile(a_tot1(A_time > 150), 98);
%             
%             if peaks_after_collapse(ii,jj,kk) < 0
%                 peaks_after_collapse(ii,jj,kk) = 0;
%             end
            
            cohesion = peaks_after_collapse;%-cr .* norm(g);
            
%             clf
%             subplot(211)
%             plot(A_time, a_tot_cohesion); xlim([0 20]);
%             hold on;
%             %         plot(A_time, a_tot2,'o'); xlim([0 1]);
%             plot(A_time, zeros(size(A_time)),'k--')
%             xlabel('Time (s)')
%             ylabel('Total parallel acceleration (m s^{-2})')
%             title(['Tilt = ',num2str(tilt_vec(jj)), ' Slope difference =', num2str(slope_diff_init_deg)])
%             ylim([1.05*min(a_tot_cohesion) 1.1*max(a_tot_cohesion)])
%             
%             subplot(212)
%             plot(A_time, a_tot_cohesion); xlim([150 170]);
%             hold on
%             %         plot(A_time, a_tot2,'o'); xlim([200 201]);
%             plot(A_time, zeros(size(A_time)),'k--')
%             xlabel('Time (s)')
%             ylabel('Total parallel acceleration (m s^{-2})')
%             ylim([1.05*min(a_tot_cohesion) 1.1*max(a_tot_cohesion)])

        end
    end
end
return
%%
clf
cmap = inferno(length(tilt_vec));

for ii = 1:length(theta_vec)
    for jj = 1:length(tilt_vec)
        for kk = 1:length(run_vec)
            if tilt_vec(jj) == 15 
                p(jj) = scatter(initial_tilt_angle(ii,jj,kk), squeeze(cohesion(ii,jj,kk)),'filled','MarkerFaceColor',cmap(jj,:),...
                    'MarkerEdgeColor','k','MarkerEdgeAlpha',0.4,'Marker','x');           
            else                
                p(jj) = scatter(initial_tilt_angle(ii,jj,kk), squeeze(cohesion(ii,jj,kk)),'filled','MarkerFaceColor',cmap(jj,:),...
                    'MarkerEdgeColor','k','MarkerEdgeAlpha',0.4);              
            end
            hold on
        end
    end
end

box on

l = legend(p, num2str((0:5:45)'),'location','southeast','Orientation','horizontal','NumColumns',2);
leg_ttl = get(l,'Title');
set(leg_ttl,'String',['Shaker tilt angle (',char(176),')'])

set(gca,'fontsize',14)

xlabel(['Shaker tilt angle (',char(176),')'])
ylabel('Cohesion coefficient C (m s^{-2})')

%%
clf
cmap = inferno(length(tilt_vec));

for ii = 1:length(theta_vec)
    for jj = 1:length(tilt_vec)
        for kk = 1:length(run_vec)
            if tilt_vec(jj) == 15 
                p(jj) = scatter(atand(intial_slope(ii,jj,kk)) - initial_tilt_angle(ii,jj,kk), squeeze(peaks_after_collapse(ii,jj,kk)),'filled','MarkerFaceColor',cmap(jj,:),...
                    'MarkerEdgeColor','k','MarkerEdgeAlpha',0.4,'Marker','x');
                
            else

                p(jj) = scatter(atand(intial_slope(ii,jj,kk)) - initial_tilt_angle(ii,jj,kk), squeeze(peaks_after_collapse(ii,jj,kk)),'filled','MarkerFaceColor',cmap(jj,:),...
                    'MarkerEdgeColor','k','MarkerEdgeAlpha',0.4);
                
            end
            hold on
        end
    end
end

box on

l = legend(p, num2str((0:5:45)'),'location','southeast','Orientation','horizontal','NumColumns',2);
leg_ttl = get(l,'Title');
set(leg_ttl,'String',['Shaker tilt angle (',char(176),')'])

set(gca,'fontsize',14)

xlabel(['Slope difference (',char(176),')'])
ylabel('Cohesion coefficient C (m s^{-2})')

return
%%
% errorbar(tilt_vec, mean(slopediff,2),std(slopediff,[],2),'.','MarkerSize',20,'LineWidth',2)
plot(mean(slope_diff,2))
xlabel('Time spent with a_{||}-\mu a_{\perp} > 0');
ylabel('Slope difference');
set(gca,'fontsize',16)
disp('Done');
return
%% Varying the initial slope angle
figure
clear
tilt_vec = 25;
theta_vec = 0:5:15;
run_vec = 1:5;

for ii = 1:length(theta_vec)
    for jj = 1:length(tilt_vec)
        for kk = 1:length(run_vec)
            load(['results/theta_',num2str(theta_vec(ii)),'_tilt_',num2str(tilt_vec(jj)),'_run_',num2str(run_vec(kk)),'.mat']);
            
            % Unpack variables:
            t = curr_time - curr_time(1);
            regolith_slope = abs(top_edge_slope);
            
            % Normalize slope:
            regolith_slope = regolith_slope./regolith_slope(1);
            
            plot(t, regolith_slope,'.');
            hold on
            slopediff(ii, jj, kk) = mean(regolith_slope(t < 1)) - mean(regolith_slope(t < (t(end) - 10)));
        end
        title(['\theta = ',num2str(theta_vec(ii)), ' Tilt = ',num2str(tilt_vec(jj))]);
        ylim([0 1.5])
        %         clf
    end
end

p = errorbar(repmat(theta_vec,3,1)', mean(slopediff,3),std(slopediff,[],3),'.','MarkerSize',20,'LineWidth',2)
set(gca,'xticklabel',{'Repose',['Repose - 5',char(176)],['Repose - 10',char(176)],['Repose - 15',char(176)]});
xlabel('Initial regolith slope');
ylabel('Slope difference');
set(gca,'fontsize',16)
legend(fliplr(p),['Shaker tilt = 25',char(176)],['Shaker tilt = 15',char(176)],['Shaker tilt = 5',char(176)])
disp('Done');