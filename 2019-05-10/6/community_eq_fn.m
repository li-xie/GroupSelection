% community_eqs
% Alex Yuan at Fred Hutch Cancer Research Center
% 10/18/2017

function [pdt,X,t] = community_eq_fn(f_p, Bio_M_init, Bio_H_init, multiplier)

%% parameters are given mono-adapted values
K_MR = 1/3 * multiplier;
K_MB = 1/3 * 100 * multiplier;
K_HR = 1/5 * multiplier;
b_Mmax = 0.7;
b_Hmax = 0.3;
d_M = 3.5 * 10^-3;
d_H = 1.5 * 10^-3;
c_RM = 10^-4;
c_RH = 10^-4;
c_BM = 1/3;
%% variables
% R == resource. NOTE: R(t=0) is set to 1 by definition.
% B == byproduct
% Bio_H == total biomass of H
% Bio_M == total biomass of M
% f_p == fraction of manufacturer energy devoted to making product

%% analytical functions whose outputs change with time
b_H = @(R) b_Hmax * R / (R + K_HR);
R_M = @(R) R / K_MR;
B_M = @(B) B / K_MB;
b_M = @(R,B) b_Mmax * (R_M(R) * B_M(B)) / (R_M(R) + B_M(B)) *...
    (1/(R_M(R) + 1) + 1/(B_M(B) + 1));
R_prime = @(Bio_H, Bio_M, R, B) -b_H(R) * c_RH * Bio_H - b_M(R,B) * c_RM * Bio_M;
B_prime = @(Bio_H, Bio_M, R, B) b_H(R) * Bio_H - b_M(R,B) * c_BM * Bio_M;
P_prime = @(Bio_M, R, B, f_p) b_M(R,B) * f_p * Bio_M;
Bio_M_prime = @(Bio_M, R, B, f_p) b_M(R,B) * (1 - f_p)  * Bio_M - d_M * Bio_M;
Bio_H_prime = @(Bio_H, R) (b_H(R) - d_H) * Bio_H;

%% initial conditions
R_init = 1 * multiplier;
B_init = 0;
P_init = 0;
% X(1) = R
% X(2) = B
% X(3) = Bio_H
% X(4) = Bio_M
% X(5) = P
X0 = [R_init; B_init; Bio_H_init; Bio_M_init; P_init];
%% do the numerical integration

Xprime = @(X,f_p) [R_prime(X(3), X(4), X(1), X(2));...
    B_prime(X(3), X(4), X(1), X(2));...
    Bio_H_prime(X(3), X(1));...
    Bio_M_prime(X(4), X(1), X(2), f_p);...
    P_prime(X(4), X(1), X(2), f_p)];
tspan = [0, 17];
options=odeset('RelTol',1e-6,'abstol',1e-10);
[t,X] = ode15s(@(t,X) Xprime(X,f_p), tspan, X0, options);

pdt = X(end,5);
end
