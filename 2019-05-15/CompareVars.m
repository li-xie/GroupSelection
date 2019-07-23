clear
num_cycles = 1e3;
num_wells = 100;
f = 3;
c = 1;
B0 = 100;
b0 = sqrt(2);
g0T = 5.5;
mu = 1e-3;
mu_sig = 5e-2;
H = 2*log(num_wells)-log(log(num_wells))-log(4*pi);
S = 10;
load([num2str(f) '/' 'fp_data' num2str(c)]);


fp0_mean = mean(fp0_all, 2);
pFpA = B0 * exp((1-fp0_mean)*g0T) .* (1-fp0_mean.*(1-fp0_mean)*g0T)./ (1-fp0_mean).^2;
pFpB = fp0_mean .* exp((1-fp0_mean)*g0T)./(1-fp0_mean);
sA = std(fp0_all,1,2);
% sA = 5e-4;
sB = std(B0_all,1,2); %sqrt(B0*b0);
alpha = pFpA.*sA;
beta = pFpB.*sB;


% plot(sA(2:end).^2, varfpT_sel(1:end-1)/B0*b0, '.')
% hold on
% plot([0 1e-5],[0 1e-5])
% hold off

histogram(var(fpT_sel,1,2)-var(fpT_all,1,2))

% figure(1)
% % plot(var(fp0_all(2:end,:),1,2),var(fpT_sel(1:end-1,:),1,2),'o')
% % plot(var(fp0_all(2:end,:),1,2),var(fpT_sel(1:end-1,:),1,2)+mean(varfpT_sel(1:end-1,:),2)/B0*1.4,'o')
% plot(var(fp0_all(2:end,:),1,2),var(fp0_sel(2:end,:),1,2),'o')
% % plot(var(fp0_all(2:end,:),1,2),beta(1:end-1)./alpha(1:end-1).*mean(varfpT_sel(1:end-1,:),2)/100*1.4,'o')
% % plot(varfp0_sel(:),varfpT_sel(:),'o')
% hold on
% plot([0 1e-5],[0 1e-5])
% hold off

% figure(1)
% for j = 1:10
%     sA = std(fp0_all(:,(j-1)*10+1:j*10),1,2);
%     plot(sA(2:end).^2, varfpT_sel(1:end-1,j)/B0*b0, '.')
%     hold on
% end
% plot([0 1e-5],[0 1e-5])
% hold off

% figure(2)
% sA_all = std(fp0_all,1,2);
% for j = 1:10
%     sA = std(fp0_all(:,(j-1)*10+1:j*10),1,2);
%     plot(sA.^2*10, sA_all.^2, '.')
%     hold on
% end
% plot([0 1e-5],[0 1e-5])
% hold off

% figure(2)
% sA_all = std(fp0_all,1,2);
% for j = 1:10
%     plot(sA_all(2:end).^2, varfpT_sel(1:end-1,j)/B0*b0, '.')
%     hold on
% end
% plot([0 1e-5],[0 1e-5])
% hold off

% figure(3)
% varA_all = var(fp0_all,1,2);
% varA_sum = mean(varfpT_sel,2)/100 + var(fp0_all,1,2);
% plot(varA_all(2:end), varA_sum(1:end-1),'.')
% hold on
% plot([0 2e-5],[0 2e-5])
% hold off
%
% figure(4)
% histogram(log10(var(fp0_all,1,2)+1e-15),'normalization','pdf')
% hold on
% histogram(log10(var(fpT_sel,1,2)+1e-15),'normalization','pdf')
% % histogram(log10(mean(varfpT_sel,2)/100*1.4+1e-15),'normalization','pdf')
% hold off
% axis([-8 -4 0 2])
% set(gca,'LineWidth',3,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 4 3],'ticklength',[0.04 0.04])%,'xticks''xticklabel',[1e-4,1e-3,1e-2,1e-1])
% xlabel('log10 of the variance')
% ylabel('PDF')