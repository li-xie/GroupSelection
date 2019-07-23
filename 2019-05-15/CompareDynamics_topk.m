clear
num_cycles = 1e3;
num_wells = 100;
f = 4;
B0 = 1e1;
b0 = sqrt(2);
g0T = 5.5;
mu = 1e-3;
mu_sig = 5e-2;
% H = (norminv(1-1/num_wells))^2;
H = 2*log(num_wells)-log(log(num_wells))-log(4*pi);
S = 10;

for c = 1:1
    load([num2str(f) '/' 'fp_data' num2str(c)]);
    
    fp0_sel_cal = zeros(num_cycles+1, 1);
    fp0_sel_cal(1) = 0.13;
    Va_m = zeros(num_cycles+1, 1);
%     Va_m = [0; mean(varfpT_sel, 2)];
    Va_inter_mean = zeros(num_cycles+1, 1);
    clr = {'b','k','r'};
    fp0_sel_cal(100) = mean(fp0_sel(100,:));
    for i = 100:1000%num_cycles
        
        fp_pre = fp0_sel_cal(i);
                Va_m(i) = (6*mu^2*mu_sig^4*fp_pre^4*(1-fp_pre)^2*log(1.4*g0T*B0))^(1/3);
                Va_inter_mean(i+1) = Va_m(i)*S*b0/B0;
%         fp0_all_temp = fp0_all(i, :);
%         Va_inter_mean(i+1) = var(fp0_all_temp(~isnan(fp0_all_temp)),1);
                alpha = (1-fp_pre*(1-fp_pre)*g0T) * B0 * exp((1-fp_pre)*g0T)...
                    /(1-fp_pre)^2*sqrt(Va_inter_mean(i+1));
                beta = fp_pre*sqrt(B0*b0)/(1-fp_pre)*exp((1-fp_pre)*g0T);
                b = norminv(1-1/num_wells,0,sqrt(alpha^2+beta^2));
                a = norminv(1-1/num_wells/exp(1),0,sqrt(alpha^2+beta^2))-b;
%         fp0_sel_cal(i+1) = Va_inter_mean(i)*(1-fp_pre*(1-fp_pre)*g0T)*B0/...
%             sqrt(Va_inter_mean(i)*(1-fp_pre*(1-fp_pre)*g0T)^2*B0^2+fp_pre^2*(1-fp_pre)^2*B0*b0)...
%             *sqrt(H)*(1-(log(S)-1)/H)-g0T*Va_m(i);
fun1 = @(x, y) CDF_sum_k(alpha*x + beta*y, a, b, 10).*exp(-x.^2/2-y.^2/2).*x;
fun2 = @(x, y) CDF_sum_k(alpha*x + beta*y, a, b, 10).*exp(-x.^2/2-y.^2/2);
        fp0_sel_cal(i+1) = integral2(fun1,-Inf,Inf,-Inf,Inf)/integral2(fun2,-Inf,Inf,-Inf,Inf)*sqrt(Va_inter_mean(i))-g0T*Va_m(i);
%                 fp0_sel_cal(i+1) = Va_inter_mean(i)*(1-fp_pre*(1-fp_pre)*g0T)/...
%                     sqrt(Va_inter_mean(i)*(1-fp_pre*(1-fp_pre)*g0T)^2+fp_pre^2*(1-fp_pre)^2*b0/B0)...
%                     *sqrt(H)*(1-(log(S)-1)/H)-g0T*Va_m(i);
        fp0_sel_cal(i+1) = fp0_sel_cal(i+1) + fp_pre;
    end
    
    figure(1)
    plot(mean(fp0_sel,2),'color',clr{c})
    hold on
    plot(fp0_sel_cal(2:end),':','color',clr{c},'linewidth',2)
    
    figure(2)
    plot(mean(varfpT_sel,2))
    hold on
    plot(Va_m)
end
figure(1)
hold off
% axis([1 1000 0 0.2])

figure(2)
% plot(varfpT_sel)
% hold on
% plot(Va_m)
hold off
