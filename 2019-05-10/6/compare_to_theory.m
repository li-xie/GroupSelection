% compare to theory
resultsfolder = 'cNphiMm100w100pu108_4/';
gen = 3;
well_num = 1;

sim_result = load(strcat(resultsfolder,'gen',num2str(gen),'.mat'));
Bio_M_sim = sim_result.Bio_M{well_num};
Bio_H_sim = sim_result.Bio_H{well_num};
R_sim = sim_result.R{well_num};
B_sim = sim_result.B{well_num};
P_sim = sim_result.P{well_num};
timestruct = load(strcat(resultsfolder,'run_conditions.mat'),'time');
time_sim = timestruct.time;
clear timestruct
if gen > 1
    nb = load(strcat(resultsfolder,'nb',num2str(gen),'.mat'));
    fp0_sim = sum(nb.newb_fp_manu{well_num} .* nb.newb_L_manu{well_num}...
        .* nb.newb_N_manu{well_num}) / sum(nb.newb_L_manu{well_num} .* nb.newb_N_manu{well_num});
else
    ic_fp_struct = load(strcat(resultsfolder,'run_conditions.mat'),'ic_fp_manu');
    fp0_sim = ic_fp_struct.ic_fp_manu(1);
end
mult_struct = load(strcat(resultsfolder,'run_conditions.mat'),'multiplier');
multiplier = mult_struct.multiplier;

[~,X,time_thry] = community_eq_fn(fp0_sim,Bio_M_sim(1),Bio_H_sim(1),multiplier);
% X(1) = R
% X(2) = B
% X(3) = Bio_H
% X(4) = Bio_M
% X(5) = P
Bio_M_thry = X(:,4);
Bio_H_thry = X(:,3);
R_thry = X(:,1);
B_thry = X(:,2);
P_thry = X(:,5);

subplot(2,3,1)
plot(time_sim,Bio_M_sim,'.')
hold on
plot(time_thry,Bio_M_thry);
title(strcat('manufacturer biomass, gen [',num2str(gen),']'))
xlabel('time')
legend('simulation','diffEQs','Location','northwest')

subplot(2,3,2)
plot(time_sim,Bio_H_sim,'.')
hold on
plot(time_thry,Bio_H_thry);
title('helper biomass')
xlabel('time')

subplot(2,3,3)
plot(time_sim,R_sim,'.')
hold on
plot(time_thry,R_thry);
title('resource')
xlabel('time')

subplot(2,3,4)
plot(time_sim,B_sim,'.')
hold on
plot(time_thry,B_thry);
title('byproduct')
xlabel('time')

subplot(2,3,5)
plot(time_sim,P_sim,'.')
hold on
plot(time_thry,P_thry);
title('product')
xlabel('time')