function dx = EvoDynamics(t,x,para)
mu = para(1);
sig = para(2);
nu = para(3);
b0 = para(4);
B0 = para(5);
g0T = para(6);
N = para(7);
S = para(8);
dt = para(9);

siga2 = (24 *mu^2 *(sig)^4 *(1-x)^2 *log(B0*g0T))^(1/3);
sigA = sqrt(S*siga2/B0*b0);
G = (1-x)*g0T;
H = 2*log(N)-log(log(N))-log(4*pi);
abratio = sigA * (1-x*G)/x/(1-x)*sqrt(B0/b0);
dx = - siga2*g0T - mu*((sig*g0T)^2 - 2*nu*g0T)*(1-x)...
    +abratio*sigA*sqrt(H)*(1-max(0, log(S)-1)/H)/sqrt(abratio^2+1);
dx = dx/dt;
