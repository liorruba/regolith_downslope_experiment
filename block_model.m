clear
clf
load ~/Desktop/oded_data.mat
% load ~/Downloads/oded_data.mat

alpha = interp1(t, alph, A_time)';

% Calculate g by averaging A
g = (mean(A,1));

% Check g is correct by rotating it to the lab frame of reference
[gx,gy,gz] = (roty_deg(-beta, g(1),g(2),g(3))); 

A_only_shaker = A - g;
v_tot = zeros(size(alpha));

for aa = 1:length(alpha)
    [a_para,~,a_perp] = roty_deg(alpha(aa) - beta, A_only_shaker(aa, 1), A_only_shaker(aa, 2), A_only_shaker(aa, 3));
    A_s_para(aa) = a_para;
    A_s_perp(aa) = a_perp;
end
A_s_perp = A_s_perp';
A_s_para = A_s_para';

mu = tand(35);
cohesion = 1;
a_tot = norm(g) .* sind(alpha) - A_s_para - mu .* (norm(g) .* cosd(alpha) + A_s_perp);

plot(A_time, a_tot,'.')
hold on
plot(A_time,zeros(size(A_time)),'k--','LineWidth',1)
%% Velocity analysis
clf
mu = tand(alpha(1));
for aa = 2:length(alpha)
    [a_para,~,a_perp] = roty_deg(alpha(aa) - beta, A_only_shaker(aa, 1), A_only_shaker(aa, 2), A_only_shaker(aa, 3));
    a_tot(aa) = norm(g) .* sind(alpha(aa)) - a_para - ...
        sign(v_tot(aa-1)) .* mu * (norm(g) .* cosd(alpha(aa)) + a_perp);
    v_tot(aa) = v_tot(aa - 1) + a_tot(aa) .* (A_time(aa) - A_time(aa - 1));
%     plot(A_time(aa), a_tot(aa),'o')
%     hold on
%     drawnow
end

hold on
plot(A_time, a_tot,'.')
% 
% % The mean A_only_shaker should be zero:
% mean(A_only_shaker,1)
% 
% % Only A_only_shaker_z should have an amplitude > 0:
% 
% 
% a_tot = norm(g) .* sind(alpha) + A_only_shaker_regolith_frame(:,1) - mu*(norm(g) .* cosd(alpha) + A_only_shaker_regolith_frame(:,3)) - A_only_shaker_regolith_frame(:,1) + A_only_shaker_regolith_frame(:,3);
% plot(A_time, a_tot,'.');
% plot(A_time, cumtrapz(A_time,a_tot),'.');
% % plot(A_time, cumtrapz(A_ti% 
% % The mean A_only_shaker should be zero:
% mean(A_only_shaker,1)
% 
% % Only A_only_shaker_z should have an amplitude > 0:
% 
% 
% a_tot = norm(g) .* sind(alpha) + A_only_shaker_regolith_frame(:,1) - mu*(norm(g) .* cosd(alpha) + A_only_shaker_regolith_frame(:,3)) - A_only_shaker_regolith_frame(:,1) + A_only_shaker_regolith_frame(:,3);
% plot(A_time, a_tot,'.');
% plot(A_time, cumtrapz(A_time,a_tot),'.');
% % plot(A_time, cumtrapz(A_time, norm(g) .* sind(alpha) - A_only_shaker_regolith_frame(:,1) - tand(alpha((find(A_time > 20, 1))))*(norm(g) .* cosd(alpha) - A_only_shaker_regolith_frame(:,3))));
% me, norm(g) .* sind(alpha) - A_only_shaker_regolith_frame(:,1) - tand(alpha((find(A_time > 20, 1))))*(norm(g) .* cosd(alpha) - A_only_shaker_regolith_frame(:,3))));

%% Pseudstatic slope analysis
% The factor of safety
FS = @(k,alph) ((norm(g).*cosd(alph) - k .* norm(g) .* sind(alph)) .* mu)...
    ./ (norm(g) .* sind(alph) + k .* norm(g) .* cosd(alph));

alphap = alpha;
% alphap(isnan(alphap)) = [];

cohesion = 0.6;

k_alpha = (mu .* cosd(alpha) - sind(alpha))./(mu .* sind(alpha) + cosd(alpha));

% for aa = 1:length(alphap)
%     fz(aa) = fzero(@(k) 1-FS(k, alphap(aa)), 1);
% end

%
t = tiledlayout(1,3);
t.TileSpacing = 'none';
nexttile
plot(A_time, (A_s_para)./norm(g))
hold on
plot(A_time(1:length(alphap)), k_alpha)
xlim([0 3])
% ylim([-0.2 0.3])
ylabel('Horizontal shaker acceleration / g');
set(gca,'fontsize',14)

nexttile
plot(A_time, (A_s_para)./norm(g))
hold on
plot(A_time(1:length(alphap)), k_alpha)
xlim([50 53])
% ylim([-0.2 0.3])
xlabel('Time (s)');
set(gca,'YTickLabel',[])
set(gca,'fontsize',14)

nexttile
plot(A_time, (A_s_para)./norm(g))
hold on
plot(A_time(1:length(alphap)), k_alpha)
xlim([100 103])
% ylim([-0.2 0.3])
set(gca,'YTickLabel',[])
set(gca,'fontsize',14)

% %% Newmark analysis
% clf
% FS = cosd(alpha).*mu./sind(alpha);
% ac = (FS - 1) .* norm(g) .* sind(alpha);
% 
% plot(A_time, A_s_para)
% hold on
% plot(A_time, ac)
% 
% greater_than_ac = A_s_para > ac;
% v = cumtrapz(A_time,(A_s_para).* greater_than_ac);
% clf
% plot(A_time, v)