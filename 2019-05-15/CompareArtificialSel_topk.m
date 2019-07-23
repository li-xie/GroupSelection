clear
num_cycles = 1e3;
num_wells = 100;
f = 3;
c = 1;
load([num2str(f) '/' 'fp_data' num2str(c)]);
B0 = 100;
b0 = sqrt(2);
g0T = 5.5;
S = 10;
k = 1;
fp0_mean = mean(fp0_all, 2);
fp0_sel_mean = mean(fp0_sel, 2);
H = 2*log(num_wells)-log(log(num_wells))-log(4*pi);
pFpA = B0 * exp((1-fp0_mean)*g0T) .* (1-fp0_mean.*(1-fp0_mean)*g0T)./ (1-fp0_mean).^2;
pFpB = fp0_mean .* exp((1-fp0_mean)*g0T)./(1-fp0_mean);
% sA = std(fp0_all,1,2);
sA = zeros(num_cycles, 1);
sA(2:end) = sqrt( S*mean(varfpT_sel(1:end-1,:), 2)*b0/B0 );
% sA(2:end) = sqrt( mean(varfpT_sel(1:end-1,:), 2)*b0/B0 + var(fpT_all(1:end-1,:),1,2)*0.9);
% sA = 5e-4;
sB = std(B0_all,1,2); %sqrt(B0*b0);
alpha = pFpA.*sA;
beta = pFpB.*sB;
dfp_sel_norm = (fp0_sel_mean-fp0_mean)./sA;
dfp_sel_knorm = (fp0_sel(:,1)-fp0_mean)./sA;
dfp_sel_cal = alpha./sqrt(alpha.^2+beta.^2)*sqrt(H)*(1-(log(S)-2)/H);%.*sA;
dfp_sel_kcal = alpha./sqrt(alpha.^2+beta.^2)*sqrt(H)*(1-(log(k)-1)/H);%.*sA;

% dB_sel_norm = (B0_sel-B0)./sB;
% dB_sel_cal = beta./sqrt(alpha.^2+beta.^2)*sqrt(H);


% histogram(dfp_sel_norm)
% hold on
% histogram(dfp_sel_cal)
% hold off

% % compare the distribution of selection with calculation
% figure(1)
% histogram((dfp_sel_norm-dfp_sel_cal)./beta.*sqrt(alpha.^2+beta.^2),'normalization','pdf')
% hold on
% plot((-3:0.1:3),normpdf((-3:0.1:3),0,1/sqrt(S)))
% hold off
% 
figure(2)
histogram((dfp_sel_knorm-dfp_sel_kcal)./beta.*sqrt(alpha.^2+beta.^2),'normalization','pdf')
hold on
plot((-3:0.1:3),normpdf((-3:0.1:3),0,1))
hold off

% plot(sA(2:end).^2,mean(varfpT_sel(1:end-1,:)/B0*b0,2),'.')
% hold on
% plot([0 1e-5],[0 1e-5])
% hold off
% %%
% figure(2)
% histogram((fp0_sel-fp0_k)./sA,'normalization','pdf')
% hold on
% plot((-5:0.1:5),normpdf((-5:0.1:5),0,sqrt(2)))
% hold off


