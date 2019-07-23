biomass_init = zeros(100,1);
avgfp = zeros(100,1);
phiM = zeros(100,1);
for well_num = 1 : 100
    biomass_init(well_num) = newb_L_manu{well_num}' * newb_N_manu{well_num} + newb_L_help{well_num} * newb_N_help{well_num};
    avgfp(well_num) = sum(newb_fp_manu{well_num} .* newb_L_manu{well_num}...
        .* newb_N_manu{well_num}) / sum(newb_L_manu{well_num} .* newb_N_manu{well_num});
    phiM(well_num) = newb_L_manu{well_num}' * newb_N_manu{well_num} / biomass_init(well_num);
end