clear
c = 1;
f = 4;
B0 = 1e2;
a = -20e-4;
b = 5e-4;
prop_name = @pipette;
num_cycles = 1e3;
num_wells = 100;
mu = 10^-3;
label = num2str(c);
S = 0;
str = sprintf('%d\n', log10(B0));
resultsfolder = strcat(num2str(f), '/', label, 'pipette',...
    'N',num2str(num_wells),'S',num2str(S),'B', str,'/');
fp0_sel = zeros(num_cycles, max(S,1));
fp0_all = zeros(num_cycles, num_wells);
B0_sel = zeros(num_cycles, max(S,1));
B0_all = zeros(num_cycles, num_wells);
varfp0_sel = zeros(num_cycles, max(S,1));
varfp0_all = zeros(num_cycles, num_wells);
g0T = 5.5;
for i = 1 : num_cycles
    load([resultsfolder 'Winnernb' num2str(i)])
    for j = 1 : max(S,1)
        B0_sel(i, j) = sum(newb_b{j}.*newb_n{j});
        fp0_sel(i, j) = sum(newb_fp{j}.*newb_b{j}.*newb_n{j})/sum(newb_b{j}.*newb_n{j});
        varfp0_sel(i, j) = sum((newb_fp{j}-fp0_sel(i, j)).^2.*newb_b{j}.*newb_n{j})/sum(newb_b{j}.*newb_n{j});
    end
    if i == 1
        continue
    else
        load([resultsfolder 'nb' num2str(i)])
        for j = 1 : num_wells
            B0_all(i, j) = sum(newb_b{j}.*newb_n{j});
            fp0_all(i, j) = sum(newb_fp{j}.*newb_b{j}.*newb_n{j})/sum(newb_b{j}.*newb_n{j});
            varfp0_all(i, j) = sum((newb_fp{j}-fp0_all(i, j)).^2.*newb_b{j}.*newb_n{j})/sum(newb_b{j}.*newb_n{j});
        end
    end
end
save([num2str(f) '/' 'fp_data'],'fp0_sel','varfp0_sel','B0_sel','fp0_all','varfp0_all','B0_all')