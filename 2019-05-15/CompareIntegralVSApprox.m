clear
num_cycles = 1e3;
num_wells = 100;
alpha = 1;
beta = 70;
B0 = 1e1;
b0 = sqrt(2);
g0T = 5.5;
mu = 1e-3;
mu_sig = 5e-2;
S = 10
b = norminv(1-1/num_wells,0,sqrt(alpha^2+beta^2));
a = norminv(1-1/num_wells/exp(1),0,sqrt(alpha^2+beta^2))-b;
fun1 = @(x, y) CDF_sum_k(alpha*x + beta*y, a, b, 10).*exp(-x.^2/2-y.^2/2).*x;
fun2 = @(x, y) CDF_sum_k(alpha*x + beta*y, a, b, 10).*exp(-x.^2/2-y.^2/2);

Va_m = zeros(num_cycles+1, 1);
Va_inter_mean = 1e-4;
fp0_sel_cal1 = zeros(num_cycles+1, 1);
fp0_sel_cal1(1) = 0.13;
fp0_sel_cal2 = zeros(num_cycles+1, 1);
fp0_sel_cal2(1) = 0.13;
for i = 1:num_cycles 
    fp_pre = fp0_sel_cal1(i)
fp0_sel_cal1(i+1) = integral2(fun1,-Inf,Inf,-Inf,Inf)/integral2(fun2,-Inf,Inf,-Inf,Inf)*sqrt(Va_inter_mean(i))-g0T*Va_m(i);
Va_m(i) = (6*mu^2*mu_sig^4*fp_pre^4*(1-fp_pre)^2*log(1.4*g0T*B0))^(1/3);
                Va_inter_mean(i+1) = Va_m(i)*S*b0/B0;