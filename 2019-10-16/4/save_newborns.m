function save_newborns(newb_fp,newb_b,newb_n,newb_n_genos,...
    newb_pnum, gen, resultsfolder)
nb.newb_fp = newb_fp;
nb.newb_b = newb_b;
nb.newb_n = newb_n;
nb.newb_n_genos = newb_n_genos;
nb.gen = gen + 1;
nb.pnum = newb_pnum;
save(strcat(resultsfolder,'nb',num2str(gen + 1)),'-struct','nb');
end