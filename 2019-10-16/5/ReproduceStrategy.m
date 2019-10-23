function offspring_num = ReproduceStrategy(num_wells, B_sorted, B0, S)
offspring_num = zeros(num_wells, 1);
if S == 0
    g_counter = 0;
    newb_counter = 0;
    while newb_counter < num_wells
        g_counter = g_counter+1;
        newb_num =  floor(B_sorted(g_counter)/B0);
        offspring_num (g_counter) = min(newb_num,num_wells-newb_counter);
        newb_counter = newb_counter + offspring_num (g_counter);
    end    
else
    if mod(num_wells, S) == 0
        offspring_num(1:S) = uint32(num_wells/S);
    elseif mod(num_wells, S-1) == 0
        S = S-1;
        offspring_num(1:S) = uint32(num_wells/S);
    elseif abs(S/num_wells-0.75)<1e-6
        temp_num = [ones(num_wells/2, 1); ones(num_wells/4, 1)*2];
        offspring_num = [temp_num(randperm(S)); zeros(num_wells/4, 1)];
    else
        error('Reproduce strategy cannot be determined!')
    end
    if (B_sorted(1:S) < offspring_num(1:S)*B0)
        error ('not enough cells to generate Newborns')
    end
end

    