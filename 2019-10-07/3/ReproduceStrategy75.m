function offspring_num = ReproduceStrategy75(num_wells, B_sorted, B0, S)
temp_num = [ones(num_wells/2, 1); ones(num_wells/4, 1)*2];
offspring_num = [temp_num(randperm(S)); zeros(num_wells/4, 1)];



    