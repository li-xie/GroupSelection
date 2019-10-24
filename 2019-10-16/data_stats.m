clear
c = 1e3;
f = 2;
N = 1e4;
S = N/2;
fp0_mean = zeros(c, 3);
varfp0_mean = zeros(c, 3);
Intervar = zeros(c, 3);
fp0_sel_mean = zeros(c, 3);
varfp0_sel_mean = zeros(c, 3);
Intervar_sel = zeros(c, 3);
for i = 1:3
    load([num2str(f) '\fp_data' num2str(i) '.mat'])
    fp0_mean(:, i) = mean(fp0_all, 2);
    varfp0_mean(:, i) = mean(varfp0_all, 2);
    Intervar(:,i) = var(fp0_all, 1, 2);
    for j = 1:c
        [~, P_ind] = sort(P_all(j, :), 'descend');
        fp0_sel_mean(j, i) = mean(fp0_all(j, P_ind(1:S)));
        varfp0_sel_mean(j, i) = mean(varfp0_all(j, P_ind(1:S)));
        Intervar_sel(j, i) = var(fp0_all(j, P_ind(1:S)), 1, 2);
    end   
end
save([num2str(f) '/fp_stat'], 'fp0_mean', 'varfp0_mean', 'fp0_sel_mean', 'varfp0_sel_mean', 'Intervar', 'Intervar_sel')
