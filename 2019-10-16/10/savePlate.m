function savePlate(wellPlate_fp, wellPlate_b, wellPlate_n,...
    wellPlate_B, wellPlate_P, wellPlate_n_genos, gen,...
    resultsfolder, cw_params, scurr)

% cw_params = [compressWinner, digits_fp, digits_b];
% if compressWinner = 1, compress the fp and biomass to digits 

if cw_params(1)
    for i = 1 : length(wellPlate_fp)
        [fp_temp, b_temp, n_temp, ngenos_temp] = compress(...
            wellPlate_fp{i}, wellPlate_b{i}, wellPlate_n{i}, cw_params(2:3));
        wellPlate_fp{i} = fp_temp;
        wellPlate_b{i} = b_temp;
        wellPlate_n{i} = n_temp;
        wellPlate_n_genos{i} = ngenos_temp;
    end
end

wellPlate.fp = wellPlate_fp;
wellPlate.b = wellPlate_b;
wellPlate.n = wellPlate_n;
wellPlate.B = wellPlate_B;
wellPlate.P = wellPlate_P;
wellPlate.n_genos = wellPlate_n_genos;
wellPlate.gen = gen;
wellPlate.DT = datetime;
wellPlate.scurr = scurr;
save(strcat(resultsfolder,'gen',num2str(gen)),'-struct','wellPlate');

end