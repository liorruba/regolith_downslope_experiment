clf
theta = 0; % deviation from repose
tilt = 20; % the tilt angle of the shaker
run = 1; % trial #

c = csvread(['/Users/liorr/Work/Projects/regolith_downslope_experiment/data/acceleration/theta_',...
    num2str(theta),'_tilt_',num2str(tilt),'_',num2str(run),'.txt'],20,0);

range = 3000:5000;

t = c(range,1);
ax = c(range,2);
ay = c(range,3);
az = c(range,4);

plot(t,ax,t,ay,t,az)

g = 1000; 
base_angle = asind(mean(c(10:1000,2))./g);
angle_of_repose = 27.1; % deg

[ax,az,apara,aperp] = calc_para_per(az, base_angle, angle_of_repose);
plot(t/1000,apara,t/1000,aperp);
xlabel('Time (s)')
ylabel('Acceleration (mg)');
set(gca,'fontsize',20)
% ;
% %%
% clf
% c00 = csvread('~/Downloads/a=0_theta=repose_1.txt',1,0);
% c02 = csvread('~/Downloads/a=0.2_theta=repose_1.txt',1,0);
% c04 = csvread('~/Downloads/a=0.4_theta=repose_1.txt',1,0);
% c06 = csvread('~/Downloads/a=0.6_theta=repose_1.txt',1,0);
% c08 = csvread('~/Downloads/a=0.8_theta=repose_1.txt',1,0);
% 
% % Calculate the initial tilt angle of the shaker:
% g = mean(c00(10:1000,4));
% ax_02_base_angle = asind(mean(c02(10:1000,2))./g);
% ax_04_base_angle = asind(mean(c04(10:1000,2))./g);
% ax_06_base_angle = asind(mean(c06(10:1000,2))./g);
% ax_08_base_angle = asind(mean(c08(10:1000,2))./g);
% angle_of_repose = 27.1; % deg
% 
% % Extract the z component of the acceleration vector
% range = 5000:7000;
% t04 = c04(range,1);
% ax_04 = c04(range,2);
% az_04 = c04(range,4);
% 
% subplot(2,2,1)
% [ax,az,apara,aperp] = calc_para_per(az_04, ax_04_base_angle, angle_of_repose);
% plot(t04-t04(1),ax, t04-t04(1), az)
% xlim([1000 2000])
% title(['Acc tilt: ', num2str(ax_04_base_angle),char(176)])
% 
% subplot(2,2,2)
% plot(t04-t04(1),apara, t04-t04(1), aperp)
% xlim([1000 2000])
% title(['Slope angle: ',num2str(angle_of_repose),char(176)])
% 
% % Extract the z component of the acceleration vector
% range = 5000:7000;
% t06 = c06(range,1);
% ax_06 = c06(range,2);
% az_06 = c06(range,4);
% 
% t08 = c08(range,1);
% ax_08 = c08(range,2);
% az_08 = c08(range,4);
% 
% [ax,az,apara,aperp] = calc_para_per(az_06, ax_06_base_angle, angle_of_repose);
% subplot(2,2,3)
% plot(t06-t06(1),ax, t06-t06(1), az)
% xlim([1000 2000])
% title(['Acc tilt: ', num2str(ax_06_base_angle),char(176)])
% 
% subplot(2,2,4)
% plot(t06-t06(1),apara, t06-t06(1), aperp)
% title(['Slope angle: ',num2str(angle_of_repose),char(176)])
% xlim([1000 2000])


