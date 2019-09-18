clear
num_cycles = 1e3;
num_wells = 100;
f = 1;
c = 1;
load([num2str(f) '/' 'fp_data' num2str(c)]);
B0 = 1e2;
b0 = sqrt(2);
g0T = 6;
k = 1;
fp0_mean = mean(fp0_all, 2);
H = 2*log(num_wells)-log(log(num_wells))-log(4*pi);
pFpA = B0 * exp((1-fp0_mean)*g0T) .* (1-fp0_mean.*(1-fp0_mean)*g0T)./ (1-fp0_mean).^2;
pFpB = fp0_mean .* exp((1-fp0_mean)*g0T)./(1-fp0_mean);
sA = std(fp0_all,1,2);
% sA = 5e-4;
sB = std(B0_all,1,2); %sqrt(B0*b0);
alpha = pFpA.*sA;
beta = pFpB.*sB;
dfp_sel = fp0_sel-fp0_mean;
dfp_sel_norm = (fp0_sel-fp0_mean)./sA;
dfp_sel_cal = alpha./sqrt(alpha.^2+beta.^2)*sqrt(H).*sA;
dB_sel_norm = (B0_sel-B0)./sB;
dB_sel_cal = beta./sqrt(alpha.^2+beta.^2)*sqrt(H);

[~,I] = sort(P_all,2,'descend');
I2 = sub2ind(size(I), (1:num_cycles)', I(:,k));
fp0_k = fp0_all(I2);

% % cycle to cycle comparison
% plot(dfp_sel, dfp_sel_cal,'.')
% hold on
% plot([0 2e-3], [0 2e-3])
% hold off
% set(gca,'LineWidth',3,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 4 3],'ticklength',[0.04 0.04])%,'xticks''xticklabel',[1e-4,1e-3,1e-2,1e-1])
% xlabel('\DeltaA* from simulation')
% ylabel('\DeltaA* from calculation')
% axis([-3e-3 3e-3 0 1e-3])

% plot(sA(2:end).^2,mean(varfpT_sel(1:end-1,:)/B0*b0,2),'.')
% hold on
% plot([0 1e-7],[0 1e-7])
% hold off


% histogram(dfp_sel_norm)
% hold on
% histogram(dfp_sel_cal)
% hold off

% compare the distribution of selection with calculation
figure(1)
histogram((dfp_sel-dfp_sel_cal)./sA./beta.*sqrt(alpha.^2+beta.^2),'normalization','pdf')
hold on
plot((-3:0.1:3),normpdf((-3:0.1:3),0,1),'linewidth',3)
hold off
set(gca,'LineWidth',3,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 4 3],'ticklength',[0.04 0.04])%,'xticks''xticklabel',[1e-4,1e-3,1e-2,1e-1])
xlabel('normalized \DeltaA* from simulation')
ylabel('PDF')
axis([-3 3 0 0.45])

% %%
% figure(2)
% histogram((fp0_sel-fp0_k)./sA,'normalization','pdf')
% hold on
% plot((-5:0.1:5),normpdf((-5:0.1:5),0,sqrt(2)))
% hold off


