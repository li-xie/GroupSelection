function saveWinner(winner_fp, winner_b, winner_n,...
    winner_B, winner_P, winner_n_genos, gen,...
    resultsfolder, cw_params, scurr)

% cw_params = [compressWinner, digits_fp, digits_b];
% if compressWinner = 1, compress the fp and biomass to digits 

if cw_params(1)
    for i = 1 : length(winner_fp)
        [fp_temp, b_temp, n_temp, ngenos_temp] = compress(...
            winner_fp{i}, winner_b{i}, winner_n{i}, cw_params(2:3));
        winner_fp{i} = fp_temp;
        winner_b{i} = b_temp;
        winner_n{i} = n_temp;
        winner_n_genos{i} = ngenos_temp;
    end
end

winner.fp = winner_fp;
winner.b = winner_b;
winner.n = winner_n;
winner.B = winner_B;
winner.P = winner_P;
winner.n_genos = winner_n_genos;
winner.gen = gen;
winner.DT = datetime;
winner.scurr = scurr;
save(strcat(resultsfolder,'winner',num2str(gen)),'-struct','winner');

end