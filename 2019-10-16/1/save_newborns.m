function save_newborns(newb_fp, newb_b, newb_n, newb_n_genos,...
    newb_pnum, seed, cycle_num, resultsfolder)
nb.newb_fp = newb_fp;
nb.newb_b = newb_b;
nb.newb_n = newb_n;
nb.newb_n_genos = newb_n_genos;
nb.cycle = cycle_num;
nb.newb_pnum = newb_pnum;
nb.seed = seed;
save(strcat(resultsfolder,'nb',num2str(cycle_num)),'-struct','nb');
end