clear
num_cycles = 1e3;
num_wells = 100;
f = 1;
B0 = 1e4;
b0 = sqrt(2);
g0T = 5.5;
mu = 1e-3;
mu_sig = 5e-2;
H = 2*log(num_wells)-log(log(num_wells))-log(4*pi);
clr = {'b','k','r'};
for c = 1:3
    load([num2str(f) '/' 'fp_data' num2str(c)]);
    
    fp0_sel_cal = zeros(num_cycles+1, 1);
    fp0_sel_cal(1) = 0.13;
    Va_mean = [0; mean(varfpT_sel,2)];
    Va_cal = zeros(num_cycles+1, 1);
    
    for i = 1:num_cycles
        fp_pre = fp0_sel_cal(i);
        Va_cal(i+1) = (6*mu^2*mu_sig^4*fp_pre^4*(1-fp_pre)^2*log(1.4*g0T*B0))^(1/3);
        D = mu*mu_sig^2*fp_pre^2*(1-fp_pre)/2;
%         Va_cal(i+1) = (6*mu^2*mu_sig^4*fp_pre^4*(1-fp_pre)^2*log(B0*D^(1/3)))^(1/3);

        Va = Va_cal(i);
        fp0_sel_cal(i+1) = (1-fp_pre*(1-fp_pre)*g0T)/(fp_pre*(1-fp_pre))*sqrt(b0/B0)...
            *sqrt(H)-g0T;
%         fp0_sel_cal(i+1) = -g0T;
        fp0_sel_cal(i+1) = fp0_sel_cal(i+1)*Va - mu*(1-fp_pre)*(fp_pre*mu_sig*g0T)^2  + fp_pre;
    end
    
    figure(1)
    plot(mean(fp0_sel,2),'color',clr{c})
    hold on
    plot(fp0_sel_cal(2:end),':','color',clr{c},'linewidth',2)
    
    figure(2)
%     plot(varfpT_sel,'color',clr{c})
%     hold on
%     plot(Va_cal(2:end),'color',clr{c},'linewidth',3)
    semilogy(mean(varfpT_sel,2)+1e-15,'color',clr{c})
    hold on
    semilogy(Va_cal(2:end)+1e-15,'color',clr{c},'linewidth',3)
end

figure(1)
hold off
set(gca,'LineWidth',3,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 4 3],'ticklength',[0.04 0.04])%,'xticks''xticklabel',[1e-4,1e-3,1e-2,1e-1])
xlabel('Cycle')
ylabel('selected fp(0)')

figure(2)
hold off
ylim([1e-7 1e-3])
set(gca,'LineWidth',3,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 4 3],'ticklength',[0.04 0.04])%,'xticks''xticklabel',[1e-4,1e-3,1e-2,1e-1])
xlabel('Cycle')
ylabel('selected variance of fp')
