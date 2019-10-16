function savePlateSummary(wellPlate_fp, wellPlate_b, wellPlate_n,...
    wellPlate_B, wellPlate_P, wellPlate_n_genos, gen,...
    resultsfolder, cw_params, scurr)

% cw_params = [compressWinner, digits_fp, digits_b];
% if compressWinner = 1, compress the fp and biomass to digits 
N = length(wellPlate_fp);
fp_mean = zeros(N, 1);
fp_var = zeros(N, 1);
B = zeros(N, 2);
P = zeros(N, 1);
if cw_params(1)
    for i = 1 : N
        B(i, :) = [wellPlate_B{i}(1) wellPlate_B{i}(end)];
        P(i) = wellPlate_P{i}(end);
        fp_mean(i) = sum(wellPlate_fp{i} .*wellPlate_b{i} .*wellPlate_n{i})/wellPlate_B{i}(end);
        fp_var(i) = sum((wellPlate_fp{i}-fp_mean(i)).^2 .*wellPlate_b{i} .*wellPlate_n{i})/wellPlate_B{i}(end);
        [~, ~, ~, ngenos_temp] = compress(...
            wellPlate_fp{i}, wellPlate_b{i}, wellPlate_n{i}, cw_params(2:3));
%         wellPlate_fp{i} = fp_temp;
%         wellPlate_b{i} = b_temp;
%         wellPlate_n{i} = n_temp;
        wellPlate_n_genos{i} = ngenos_temp;
    end
end

wellPlate.fp_mean = fp_mean;
wellPlate.fp_var = fp_var;
wellPlate.B = B;
wellPlate.P = P;
wellPlate.n_genos = wellPlate_n_genos;
wellPlate.gen = gen;
wellPlate.DT = datetime;
wellPlate.scurr = scurr;
save(strcat(resultsfolder,'gen',num2str(gen)),'-struct','wellPlateSum');

end