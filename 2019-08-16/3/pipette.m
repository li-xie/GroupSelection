function [newb_fp,newb_b,newb_n,newb_n_genos] = pipette(...
    newb_fp,newb_b,newb_n,newb_n_genos,...
    winner_fp,winner_b,winner_n, winner_n_genos,...
    target_inds,nD,cw_params,len_fp_init,B0)

%PIPETTE Summary of this function goes here
%  nD: fold of dilution, calculated from uncompressed biomass


% execute data compression if desired
if cw_params(1)
    [winner_fp, winner_b, winner_n, winner_n_genos] =...
        compress(winner_fp, winner_b, winner_n,cw_params(2:3));
end

% distribute manufacturer cells from winning well
probs = ones(length(winner_n),nD) * 1/nD;
% multinomial distribution
newb_n_mat = mnrnd(winner_n,probs);
newb_n_sum = sum(newb_n_mat, 1);
ind = find(newb_n_sum>0);
newb_n_mat = newb_n_mat(:, ind);
for i = 1 : length(target_inds)
    [nb_fp_temp,nb_b_temp,nb_n_temp] = removeZeros(winner_fp, winner_b, newb_n_mat(:,i));
    temp_n_genos = length(nb_fp_temp);
    if temp_n_genos < len_fp_init
        nb_fp_temp = [nb_fp_temp;zeros(len_fp_init-temp_n_genos,1)];
        nb_b_temp = [nb_b_temp;zeros(len_fp_init-temp_n_genos,1)];
        nb_n_temp = [nb_n_temp;zeros(len_fp_init-temp_n_genos,1)];
    end
    newb_fp{target_inds(i)} = nb_fp_temp;
    newb_b{target_inds(i)} = nb_b_temp;
    newb_n{target_inds(i)} = nb_n_temp;
    newb_n_genos{target_inds(i)} = temp_n_genos;
end

% % distribute helper cells from winning well
% 
% newb_NH_vec = mnrnd(winner_N_help,ones(1,nD) * 1/nD);
% for i = 1 : length(target_inds)
%     newb_N_help{target_inds(i)} = newb_NH_vec(i);
%     newb_L_help{target_inds(i)} = winner_L_help;
% end
end

