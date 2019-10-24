clear
f = 1;
B0 = 1e2;
N = 1e3;
S = N/2;
num_cycles = 1e3;

for c = 1:3
    fp0_all = zeros(num_cycles, N);
    fpT_all = zeros(num_cycles, N);
    P_all = zeros(num_cycles, N);
    varfp0_all = zeros(num_cycles, N);
    varfpT_all = zeros(num_cycles, N);
    
    foldername = [num2str(f) '/' num2str(c) 'pipetteN' num2str(N) 'S' num2str(S) 'B2/'];
    for i = 1:num_cycles
        filename = [foldername 'nb' num2str(i)];
        load(filename)
%         for j = 1:N
%             fp0_all(i, j) = sum(newb_fp{j}.*newb_b{j}.*newb_n{j})/sum(newb_b{j}.*newb_n{j});
%             varfp0_all(i, j) = sum((newb_fp{j}-fp0_all(i, j)).^2.*newb_b{j}.*newb_n{j})/sum(newb_b{j}.*newb_n{j});
%         end
        filename = [foldername 'gen' num2str(i)];
        load(filename)
        fpT_all(i, :) = fp_mean';
        varfpT_all(i, :) = fp_var';
        P_all(i, :) = P';
    end
%     save([num2str(f) '/' 'fp_data' num2str(c)], 'fp0_all', 'fpT_all',...
%         'varfp0_all', 'varfpT_all', 'P_all')
    save([num2str(f) '/' 'fp_data' num2str(c)], 'P_all', '-append')
end
