regolith_length = 0.278; % m;

deg_5_ini = [62.43 64.52 61.76];
deg_5_fin = [68.08 68.95 70.49];

deg_10_ini = [67.59 67.20 68.25];
deg_10_fin = [70.49 70.10 68.57];

deg_15_ini = [76.20 72.15 73.82];
deg_15_fin = [76.19 73.33 74.97];

mean_diff_5 = mean(deg_5_fin - deg_5_ini);
std_diff_5 = std(deg_5_fin - deg_5_ini);
mean_diff_10 = mean(deg_10_fin - deg_10_ini);
std_diff_10 = std(deg_10_fin - deg_10_ini);
mean_diff_15 = mean(deg_15_fin - deg_15_ini);
std_diff_15 = std(deg_15_fin - deg_15_ini);

mass_tran_5 = 0.03 * (regolith_length/2).^2 * sind((mean_diff_5))/2*1070;
mass_tran_10 = 0.03 * (regolith_length/2).^2 * sind((mean_diff_10))/2*1070;
mass_tran_15 = 0.03 * (regolith_length/2).^2 * sind((mean_diff_15))/2*1070;

mass_tran_5_sterr = 0.03 * (regolith_length/2).^2 * sind((std_diff_5))/2*1070;
mass_tran_10_sterr = 0.03 * (regolith_length/2).^2 * sind((std_diff_10))/2*1070;
mass_tran_15_sterr = 0.03 * (regolith_length/2).^2 * sind((std_diff_15))/2*1070;

clf
ang_vec = fliplr(90-mean(deg_5_ini) - (0:5:15));
mass_flux_vec = fliplr([0.44e-4*1070 mass_tran_5 mass_tran_10 mass_tran_15]);
mass_flux_vec_std =  fliplr([mass_tran_5_sterr*0.86 mass_tran_5_sterr mass_tran_10_sterr mass_tran_15_sterr]);
plot(ang_vec,mass_flux_vec,'o')
hold on
model_fun = @(b,x) b(1).*x.^(b(2));
b = nlinfit(tand(ang_vec),mass_flux_vec,model_fun,[1 1]);
plot(linspace(10,30), model_fun(b,tand(linspace(10,30))));

xlim([10 30]);
set(gca,'xtick',[10 20 30],'XTickLabel',{['tan 10',char(176)], ['tan 20',char(176)], ['tan 30',char(176)]});
set(gca,'ytick',[0:0.02:0.08])
xlabel('Slope angle');
ylabel('Total mass flux (kg)');
ax1_pos = get(gca,'Position');
set(gca,'fontsize',16)
ylim([0 0.08])
ax2 = axes('Position',ax1_pos,...
    'Color','none');
errorbar(ang_vec,mass_flux_vec,mass_flux_vec_std,'o','MarkerSize',10,'MarkerFaceColor',co(1,:))
hold on
plot(ax2,linspace(10,30), model_fun(b,tand(linspace(10,30))),'linewidth',2);
set(ax2,'xtick',ang_vec,'XAxisLocation','top','XTickLabel',...
    {['Repose-',num2str(15),char(176)],['Repose-',num2str(10),char(176)],...
    ['Repose-',num2str(5),char(176)], ['Repose']},'yticklabel',[]);
box on;
ylim([0 0.08])
legend('Data','Power law fit, 0.24 \cdot s^{2.83}','Location','northwest')
set(gca,'fontsize',18)
return
%% Craters

mass_tran_5_m2 = mass_tran_5 / 0.03 / (regolith_length/2) / 1070
mass_tran_10_m2 = mass_tran_10 / 0.03 / (regolith_length/2) / 1070
mass_tran_15_m2 = mass_tran_15 / 0.03 / (regolith_length/2) / 1070

b_m2 = nlinfit(tand(ang_vec),mass_flux_vec/ 0.03 / (regolith_length/2) ,model_fun,[1 1])

dd = 0.1;
crater_rad = 400;
Z = sphericalCrater(crater_rad, dd, 801, [0 0]); Z(Z==0) = nan;
[u,v] = gradient(Z);
U = sqrt(u.^2+v.^2); U(U<0.01) = nan;
[f,hx] = histcounts(U(:),50,'Normalization','pdf');
hx = (hx(1:end-1) + hx(2:end))/2;
% plot(hx,f)
m = model_fun(b_m2,U);

nansum(m(:))./(pi .* (crater_rad.^2 + (dd.*2.*crater_rad).^2))

% depthToDiameter * 2 * craterRadius * (1 - pow(distanceFromCraterCenter/craterRadius,2));