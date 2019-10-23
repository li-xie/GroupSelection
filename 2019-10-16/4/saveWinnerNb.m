function saveWinnerNb(newb_fp,newb_b,newb_n,newb_n_genos,...
    newb_pnum, gen, resultsfolder)
Winnernb.newb_fp = newb_fp;
Winnernb.newb_b = newb_b;
Winnernb.newb_n = newb_n;
Winnernb.newb_n_genos = newb_n_genos;
Winnernb.gen = gen;
Winnernb.pnum = newb_pnum;
save(strcat(resultsfolder,'Winnernb',num2str(gen)),'-struct','Winnernb');
end