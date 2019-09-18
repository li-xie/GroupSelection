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
    else
        offspring_num(S) = mod(num_wells, S-1);
        offspring_num(1:S-1) = uint32((num_wells - offspring_num(S)) / (S-1));
    end
    blank_num = nnz(B_sorted(1:S) < 1e-12);
    rand_pick = randi(num_wells - blank_num, blank_num, 1);
    for i = 1:blank_num
        offspring_num(rand_pick(i)) = offspring_num(rand_pick(i))+1;
    end
    offspring_num = offspring_num(1 : num_wells - blank_num);
    S = S-blank_num;
    if nnz(B_sorted(1:S) < offspring_num(1:S)*B0)>1
        error ('not enough cells to generate Newborns')
    end
end

    