function [fp, b, n, n_genos_curr] =...
    mutation(fp, b, n, n_genos_curr, pot_mut_index,...
    p_mut, sp0, sn0, fp_min, fp_max)
% generate vector of mutant cell numbers
if sum(n(pot_mut_index)) > 0
    N_mut = fastbinorv(n(pot_mut_index),p_mut);
    if sum(N_mut) > 0
        b_mut = b(pot_mut_index, 1);
        fp_mut = fp(pot_mut_index, 1);
        % remove mutants from their original genotypes
        n(pot_mut_index) = n(pot_mut_index, 1) - N_mut;
        % remove elements of N_mut, b_mut, and fp_mut with 0 cells.
        [fp_mut,b_mut,N_mut] = removeZeros(fp_mut, b_mut, N_mut);
        
%         % generate vector of null mutant cell numbers
%         N_null = fastbinorv(N_mut,frac_null);
%         N_am = N_mut - N_null;
%         % make this more space efficient.
%         [fp_null,b_null,N_null] = removeZeros(fp_mut,b_mut,N_null);
%         % update fp, b and n, and n_genos_curr. Different cells have
%         % different biomass when they become null and are counted as
%         % different.
%         fp(n_genos_curr + 1 :n_genos_curr + length(N_null)) = fp_null * 0;
%         b(n_genos_curr + 1 :n_genos_curr + length(N_null)) = b_null;
%         n(n_genos_curr + 1 :n_genos_curr + length(N_null)) = N_null;
%         n_genos_curr = n_genos_curr + length(N_null);
        
        N_am = N_mut; 
        % generate vector of active mutant (am) cell numbers
        [fp_am,b_am,N_am] = removeZeros(fp_mut,b_mut,N_am);
        % fp_back contains information for all mutations in fp
        % b_back contains the corresponding biomass
        fp_back = zeros(sum(N_am),1);
        b_back = zeros(sum(N_am),1);
        % reformat vectors to execute non-null mutation.
        mut_count = 1; % x is the counter for mutations
        for i = 1 : length(N_am)
            fp_back(mut_count:mut_count+N_am(i)) = fp_am(i);
            b_back(mut_count:mut_count+N_am(i)) = b_am(i);
            mut_count = mut_count + N_am(i) + 1;
        end
        if any(fp_back == 0)
            error('something went wrong here')
        end
        N_back = ones(length(fp_back),1);
        % update fp_manu, L_manu, N_manu, and n_genos_curr
        params = [sp0, sn0, fp_min, fp_max];
        fp(n_genos_curr + 1 : n_genos_curr + length(fp_back), 1) = mut_spec_add(params,fp_back);
        b(n_genos_curr + 1 : n_genos_curr + length(fp_back), 1) = b_back;
        n(n_genos_curr + 1 : n_genos_curr + length(fp_back), 1) = N_back;
        n_genos_curr = n_genos_curr + length(fp_back);
    end
end
end