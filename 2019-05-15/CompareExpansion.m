clear
num_cycles = 1e3;
num_wells = 100;
f = 2;
c = 1;
load([num2str(f) '/' 'fp_data' num2str(c)]);
B0 = 10;
b0 = sqrt(2);
g0T = 5.5;
fp0_mean = mean(fp0_all, 2);

pFpA = B0 * exp((1-fp0_mean)*g0T) .* (1-fp0_mean.*(1-fp0_mean)*g0T)./ (1-fp0_mean).^2;
pFpB = fp0_mean .* exp((1-fp0_mean)*g0T)./(1-fp0_mean);

sA = std(fp0_all,1,2)*sqrt(b0/B0);
sB = sqrt(B0*b0);
P_exp = mean(P_all, 2) + pFpA*ones(1, num_wells) .* (fp0_all-fp0_mean*ones(1, num_wells))...
    + pFpB*ones(1, num_wells) .* (B0_all-B0);
P_sel_exp = mean(P_all, 2) + pFpA .* (fp0_sel-fp0_mean)...
    + pFpB .* (B0_sel-B0);

plot(P_all(:), P_exp(:), '.')
hold on
plot(P_sel(:), P_sel_exp(:), '.')
a = 0;
b = 600;
plot([a b],[a b],'k:','linewidth',2)
hold off
set(gca,'LineWidth',3,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 4 3],'ticklength',[0.04 0.04])%,'xticks''xticklabel',[1e-4,1e-3,1e-2,1e-1])
xlabel('P from simulations')
ylabel('P from linearization')
axis([a b a b])
% plot(P_sel(:), P_sel_exp(:), '.')
% hold on
% plot([0 2500],[0 2500])
% hold off