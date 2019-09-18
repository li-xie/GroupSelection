function y = GoodEqs(x)
B0 = 1e1;
mu = 10^-3;
sigQ = 5e-2;
T = 5.5;
A0 = .13;
N = B0;
Ub = mu*(1-A0)*2;
sigma = A0 * sigQ;
v = x(:,1);
xc = x(:,2);
y1 = Ub .*sqrt(2*pi./v) .*(1 + sigma./xc + v./(sigma*xc-v))...
    .*exp((xc-v/sigma).^2/2./v)-2;
y2 = N*Ub .*(xc.^2./v -1 +2*xc*sigma./v +2*sigma^2./v) .*exp(-(xc-v/2/sigma)/sigma)-1;
y = [y1 y2];