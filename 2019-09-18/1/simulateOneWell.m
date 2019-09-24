function [fp, b, n, B, n_genos, P, gc]...
    = simulateOneWell(gc, ic_fp, ic_b, ic_n, ic_n_genos)
%UNTITLED2 Summary of this function goes here
% gc is a structure of global constants and is not allowed to change.

%% set up data structures

fp = ic_fp; % initial condition for fp
b = ic_b; % column vector for the initial condition for individual biomass
n = ic_n; % column vector for the initial number of cells
n_genos_curr = ic_n_genos; % current number of genotypes
B = [b' * n; zeros(gc.nsteps,1)]; % biomass at every time step
P = zeros(gc.nsteps + 1,1);
%% simulate!
for tstep = 1 : gc.nsteps
    % keep track of b and n that were calculated in previous step
    b_old = b;
    n_old = n;
    b = b .* exp((1 - fp) * gc.dt);
    % find the cells whose biomass are larger than 2 
    lgt2 = b > 2;
    % find the cells that might mutate
    pot_mut_index = and(b > 2, fp > 0); 
    n(lgt2) = n(lgt2) * 2;
    b(lgt2) = b(lgt2) / 2;
    % update product
    P(tstep + 1) = P(tstep) + sum(fp ./ (1 - fp) .* (b .* n - b_old .* n_old));
    % update n by implementing death
    n = n - fastbinorv(n, gc.d_M * gc.dt);  %round(N_manu * gc.d_M * gc.dt);
    % mutation
    [fp, b, n, n_genos_curr] = mutation(fp, b, n, n_genos_curr, pot_mut_index,...
        gc.p_mut, gc.sp0, gc.sn0, gc.fp_min, gc.fp_max);
    % update biomass
    B(tstep + 1) = b' * n;
end
n_genos = n_genos_curr;
end
