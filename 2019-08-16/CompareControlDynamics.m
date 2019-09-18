clear
f = 3;
B0 = 1e1;
num_cycles = 3e3;
num_wells = 100;
mu = 10^-3;
sigQ = 5e-2;
S = 100;
T = 5.5;
A0 = .13;

N = B0*T;
Ub = mu*(1-A0)*2;
sigma = A0 * sigQ;
% v1 = sigma^2 * (log(2*N*Ub*log(sigma/Ub)))^2 /2 ...
%     /log(sigma/Ub/sqrt(pi)*log(2*N*Ub)/sqrt(log(sigma/Ub)));
v1 = sigma^2 * (log(N*Ub))^2 /2 /log(sigma/Ub);
v2 = (24 *mu^2 *(sigQ*A0)^4 *(1-A0)^2 *log(B0*T))^(1/3);
figure(2)
for c = 1:3
    load([num2str(f) '/' 'fp_data' num2str(c)])
    plot((1:num_cycles), mean(fpT_all,2))
    hold on
end
plot((1:num_cycles), A0-(1:num_cycles)*1.5e-6*T)
% plot((1:num_cycles), A0-(1:num_cycles)*v2*T,':')
hold off