clear
foldername = '1pipetteN100S100B1';
cycle = 558;
Sp = 100;
% prop_name = @pipette;
% load([foldername '/Winner' num2str(cycle)]);
load([foldername '/gen' num2str(cycle)]);
load([foldername '/run_conditions'])
load([foldername '/gc-struct'])

% P_final = zeros(num_wells, 1);
% B_final = zeros(num_wells, 1);
% for i = 1 : num_wells
%     P_final(i) = wellPlate_P{i}(end);
%     B_final(i) = wellPlate_B{i}(end);
% end
% [~, P_ind] = sort(P_final,'descend');
% B_sorted = B_final(P_ind);

rng(scurr)
cw_params = [1 5 2];
len_fp_init = 100;
B0 = gc.B0;
newb_fp = cell(num_wells, 1);
newb_b = cell(num_wells, 1);
newb_n = cell(num_wells, 1);
newb_n_genos = cell(num_wells, 1);
newb_pnum = zeros(num_wells, 1);
num_wells_filled = 0;
rep_profile = ones(Sp,1)*1;
    for i = 1 : 100
        target_inds = num_wells_filled + 1 : num_wells_filled + rep_profile(i);
        num_wells_filled = num_wells_filled + rep_profile(i);
        newb_pnum(target_inds) = i;
        [newb_fp,newb_b,newb_n, newb_n_genos,] = pipette(...
            newb_fp,newb_b,newb_n,newb_n_genos,...
            fp{i},b{i},n{i},n_genos{i},target_inds,floor(B{i}(end)/B0),cw_params,len_fp_init,B0);
    end