clear
num_cycles = 1e3;
num_wells = 100;
c = 1;
B0 = 1e2;
b0 = sqrt(2);
g0T = 5.5;
mu = 1e-3;
mu_sig = 5e-2;
H = 2*log(num_wells)-log(log(num_wells))-log(4*pi);

siga = (6*mu^2*mu_sig^4*log(1.4*g0T*B0))^(1/3);

f = 3;
S = 10;
load([num2str(f) '/' 'fp_data' num2str(c)]);

fp0_mean = mean(fp0_sel,2);
normsiga = mean(varfp0_sel,2)./(fp0_mean).^(4/3)./(1-fp0_mean).^(2/3);
plot(normsiga)
hold on

f = 1;
S = 1;
load([num2str(f) '/' 'fp_data' num2str(c)]);

fp0_mean = mean(fp0_sel,2);
normsiga = mean(varfp0_sel,2)./(fp0_mean).^(4/3)./(1-fp0_mean).^(2/3);
plot(normsiga)
plot([1 num_cycles], [siga siga])
hold off