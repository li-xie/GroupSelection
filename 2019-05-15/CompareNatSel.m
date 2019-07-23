clear
num_cycles = 1e3;
num_wells = 100;
f = 5;
c = 1;
load([num2str(f) '/' 'fp_data' num2str(c)]);
B0 = 1e4;
b0 = sqrt(2);
g0T = 5.5;
a = (fp0_all(:) - fpT_all(:) - varfp0_all(:)*g0T);
figure(1)
histogram(fp0_all(:) - fpT_all(:))
hold on
histogram(varfp0_all(:)*g0T)
hold off

figure(2)
plot(fp0_all(:) - fpT_all(:), varfp0_all(:)*g0T, '.')
hold on
plot([-2e-3, 4e-3], [-2e-3, 4e-3],'linewidth',2)
hold off
axis([-5e-3 5e-3 0 3e-3])
set(gca,'LineWidth',3,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 4 3],'ticklength',[0.04 0.04])%,'xticks''xticklabel',[1e-4,1e-3,1e-2,1e-1])
xlabel('A(0)-A(T) from simulation')
ylabel('\sigma_a^2 \times g_0T')