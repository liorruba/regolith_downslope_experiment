clf
clear
co = get(gca,'ColorOrder');

regolith_length = 0.278; % m;

tilt_vec = 25;
theta_vec = 0;
run_vec = [1];

model_fun = @(b,x) b(1) * exp(b(2) .* (x)) + b(3);

for jj = 1:length(run_vec)
    for ii = 1:length(tilt_vec)
        load(['results/theta_',num2str(theta_vec(1)),'_tilt_',num2str(tilt_vec(jj)),'_run_',num2str(run_vec(jj)),'.mat']);
%         load(file);
        
        time = curr_time - curr_time(1);
        slope = abs(top_edge_slope);
        
        time(isnan(slope)) = [];
        slope(isnan(slope)) = [];
        
        [b ,R,~,CovB] = nlinfit(time, slope, model_fun,[1 -.2 .1]);
        berr = nlparci(b,R,'covar',CovB);
        
        relax_ts(ii,jj) = b(2);
        relax_ts_err(ii,jj) = b(2) - berr(2);
        angle_diff(ii,jj) = atand(abs(model_fun(b, time(1)))) - atand(model_fun(b, time(200)));
        angle_diff_95_plus(ii,jj) = atand(abs(model_fun(berr(:,2)', time(1)))) - atand(model_fun(berr(:,2)', time(200)));
        angle_diff_95_minus(ii,jj) = atand(abs(model_fun(berr(:,1)', time(1)))) - atand(model_fun(berr(:,1)', time(200)));
        if run_vec(jj) == 1 
            plot(time,slope,'.','color',co(ii,:),'markersize',10);
            hold on
            p(ii) = plot(time, model_fun(b, time),'-','LineWidth',3,'color',co(ii,:));
        end
        
    end
end

disp_volume = (regolith_length/2).^2 * sind(mean(angle_diff,2))/2 .* 0.03;
% disp_volume_plus = (regolith_length/2).^2 * sind(mean(angle_diff_95_plus,2))/2 .* 0.03
% disp_volume_minus = (regolith_length/2).^2 * sind(mean(angle_diff_95_minus,2))/2 .* 0.03

xlim([0 210]);
ylim([0.3 1]);
xlabel('Elapsed time (sec)');
ylabel('|Slope|')
set(gca,'fontsize',16)

hold on
plot([50 50],[0.3 1],'k--')
text(5,0.95,'Diffusion Stage','fontsize',14)
text(55,0.95,'Relaxation Stage','fontsize',14)

legend(p,['\alpha_s=12.21',char(176)],['\alpha_s=24.41',char(176)],['\alpha_s=33.91',char(176)])
return

%% Plot acc
clf
amax_reg = [   3.217655688397867   2.997787302327138   2.377001102512141   1.822839580149090
   8.374060509797172   8.472064154912481   9.869110708881202   9.776333831964989] * 1e2;

subplot(2,1,1);
plot(amax_reg(1,:), disp_volume,'-o')
xlabel('Amplitude of a_{\mid\mid}');
ylabel('Displaced volume (m^3)');
subplot(2,1,2);
plot(amax_reg(2,:), disp_volume,'-o')
xlabel('Amplitude of a_{\perp}');
ylabel('Displaced volume (m^3)');