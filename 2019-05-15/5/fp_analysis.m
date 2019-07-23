c = 1;
B0 = 1e2;
prop_name = @pipette;
num_cycles = 1e3;
num_wells = 100;
mu = 10^-3;
label = num2str(c);
S = 0;
str = sprintf('%d\n', log10(B0));
resultsfolder = strcat(label, 'pipette',...
    'N',num2str(num_wells),'S',num2str(S),'B', str,'/');
fp0 = zeros(num_cycles, 1);
varfp0 = zeros(num_cycles, 1);
g0T = 5.5;
for i = 1 : num_cycles
    load([resultsfolder 'Winnernb' num2str(i)])
    fp0(i) = sum(newb_fp{1}.*newb_b{1}.*newb_n{1})/sum(newb_b{1}.*newb_n{1});
    varfp0(i) = sum((newb_fp{1}-fp0(i)).^2.*newb_b{1}.*newb_n{1})/sum(newb_b{1}.*newb_n{1});
end
% save('fp_data','fp0','varfp0')
%%
load('fp_data')
figure(1)
histogram(log10(varfp0),'normalization','pdf')
title('PDF of intra-group variance in selected Newborn')
xlabel('Variance x10^{-5}')
ylabel('PDF')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 4 4],'ticklength',[0.04 0.04])%,'xticklabel',[])
xlim([-7 -3])
% mu = 0;
% varfp0 = (6*mu^2*0.05^4*(1-fp0).^2*log(1.4/1e-4)).^1/3;
figure(2)
rhd = varfp0(1:end-1) .* (1 - fp0(1:end-1).*(1-fp0(1:end-1))*g0T)./(fp0(1:end-1).*(1-fp0(1:end-1)))...
    *sqrt(sqrt(2)/B0)*(2*log(num_wells)-log(log(num_wells))-log(4*pi))- varfp0(1:end-1).*g0T...
    -mu*(1-fp0(1:end-1)).*(fp0(1:end-1)*0.05^2*g0T).^2;
scatter(diff(fp0), rhd,[],(1 : num_cycles-1),'filled')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 4 4],'ticklength',[0.04 0.04])%,'xticklabel',[])
title('\Delta A*(0) comparison')
xlabel('Simulation results')
ylabel('Analytical prediction')
% axis([-2e-4, 1e-4, -2e-4, 1e-4])