function [newb_fp,newb_b,newb_n,newb_n_genos]...
    = cell_sort(newb_fp,newb_b,newb_n,newb_n_genos,...
    winner_fp,winner_b,winner_n, winner_n_genos,...
    target_inds,nD,cw_params,len_fp_init,B0)


if length(target_inds) > nD
    error('The length of target_inds must be less than nD')
end

% execute data compression if desired
if cw_params(1)
    [winner_fp, winner_b, winner_n, winner_n_genos] =...
        compress(winner_fp, winner_b, winner_n,cw_params(2:3));
end

%%

probs = ones(length(winner_n),nD) * 1/nD;
% do regular pipetting, using the multinomial distribution
% where helper cells are included in the array
% cellmat is a length(winner_n) by nD matrix with the number of cells for
% each genotype
cellmat = mnrnd(winner_n,probs);
donation = zeros(length(winner_n),1);
% take "taxes" or "donations" from all columns of cellmat that have excess
% cells.
for i = 1 : nD
    % draw cells from groups whose biomass are larger than B0  
    while winner_b' * cellmat(:,i) > B0
        ne_inds = cellmat(:,i)>0;
        temp1 = cellmat(ne_inds,i);
        temp2 = donation(ne_inds);
        % draw one cell from a group to donation
        [temp1,temp2] = draw_one_cell(temp1,temp2);
        cellmat(ne_inds,i) = temp1;
        donation(ne_inds) = temp2;
    end
end
% give back cells to any wells that are lacking cells, until they achieve
% biomass over B0. Then take back the last cell you gave.

% TODO: deal with unlikely edge case where you don't have enough cells to
% give and then take. I think i took care of this edge case because there
% should always be enough biomass to top every newborn off.
for i = 1 : nD
    if nnz(donation<0)>0
        error('Not enough cells to donate!')
    end
    % give
    while winner_b' * cellmat(:,i) < B0
        ne_inds = donation > 0;
        temp1 = cellmat(ne_inds,i);
        temp2 = donation(ne_inds);
        % draw one cell from donation to a group
        [temp2,temp1,indx] = draw_one_cell(temp2,temp1);
        cellmat(ne_inds,i) = temp1;
        donation(ne_inds) = temp2;
    end
    % take out the last cell when the total biomass of the group exceeds B0
    if winner_b' * cellmat(:,i) > B0
        temp1(indx) = temp1(indx) - 1;
        temp2(indx) = temp2(indx) + 1;
        cellmat(ne_inds,i) = temp1;
        donation(ne_inds) = temp2;
    end
end

%% distribute cells
for i = 1 : length(target_inds)
    [nb_fp_temp,nb_b_temp,nb_n_temp] = removeZeros(winner_fp, winner_b, cellmat(:,i));
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
end

