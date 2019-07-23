clear
c = 1;
f = 10;
B0 = 1e2;
a = -20e-4;
b = 5e-4;
prop_name = @pipette;
num_cycles = 1e3;
num_wells = 100;
mu = 10^-3;
label = num2str(c);
S = 10;
str = sprintf('%d\n', log10(B0));
resultsfolder = strcat(num2str(f), '/', label, 'pipette',...
    'N',num2str(num_wells),'S',num2str(S),'B', str,'/');
fp0 = zeros(num_cycles, max(S,1));
varfp0 = zeros(num_cycles, max(S,1));
g0T = 5.5;
for i = 1 : num_cycles
    load([resultsfolder 'Winnernb' num2str(i)])
    for j = 1 : max(S,1)
        fp0(i, j) = sum(newb_fp{j}.*newb_b{j}.*newb_n{j})/sum(newb_b{j}.*newb_n{j});
        varfp0(i, j) = sum((newb_fp{j}-fp0(i, j)).^2.*newb_b{j}.*newb_n{j})/sum(newb_b{j}.*newb_n{j});
    end
end
save([num2str(f) '/' 'fp_data'],'fp0','varfp0')
%%
% load([num2str(f) '/' 'fp_data'])
fp0 = mean(fp0 ,2);
varfp0 = mean(varfp0 ,2);
figure(1)
histogram(log10(varfp0),'normalization','pdf')
title('PDF of intra-group variance in selected Newborn')
xlabel('log10(Variance)')
ylabel('PDF')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 4 4],'ticklength',[0.04 0.04])%,'xticklabel',[])
xlim([-8 -2])
% mu = 0;
% varfp0 = (6*mu^2*0.05^4*(1-fp0).^2*log(1.4/1e-4)).^1/3;
figure(2)
H = (2*log(num_wells)-log(log(num_wells))-log(4*pi));
rhd = varfp0(2:end) .* (1 - fp0(1:end-1).*(1-fp0(1:end-1))*g0T)./(fp0(1:end-1).*(1-fp0(1:end-1)))...
    *sqrt(sqrt(2)/B0*H)*(1-(log(max(S,1))-1)/H)- varfp0(2:end).*g0T...
    -mu*(1-fp0(1:end-1)).*(fp0(1:end-1)*0.05*g0T).^2;
% rhd = - varfp0(2:end).*g0T;
scatter(diff(fp0), rhd,[],(1 : num_cycles-1),'filled')
hold on
plot([min(diff(fp0)) max(diff(fp0))], [min(diff(fp0)) max(diff(fp0))],'k:','linewidth',2)
hold off
% axis([a b a b])
set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 4 4],'ticklength',[0.04 0.04])%,'xticklabel',[])
title('\Delta A*(0) comparison')
xlabel('Simulation results')
ylabel('Analytical prediction')
% axis([-2e-4, 1e-4, -2e-4, 1e-4])