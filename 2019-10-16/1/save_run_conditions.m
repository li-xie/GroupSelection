function save_run_conditions(cycle_duration, num_cycles, num_wells,...
    compressWinner, gc, time, ic_fp, ic_b, ic_n, ...
    ic_n_genos,resultsfolder, seeds)

runconditions.cycle_duration = cycle_duration;
runconditions.num_cycles = num_cycles;
% runconditions.multiplier = multiplier;
runconditions.num_wells = num_wells;
runconditions.compressWinner = compressWinner;
runconditions.time = time;
runconditions.ic_fp = ic_fp;
runconditions.ic_b = ic_b;
runconditions.ic_n = ic_n;
% runconditions.ic_L_help = ic_L_help;
% runconditions.ic_N_help = ic_N_help;
runconditions.ic_n_genos = ic_n_genos;
runconditions.seeds = seeds;
%runconditions.scoremode = scoremode;
save(strcat(resultsfolder,'run_conditions'),'-struct','runconditions');
save(strcat(resultsfolder,'gc','-struct'),'gc');

end
