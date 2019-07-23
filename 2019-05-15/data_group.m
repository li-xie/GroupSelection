clear

f = 6;
B0 = 1e4;
a = -20e-4;
b = 5e-4;
prop_name = @pipette;
num_cycles = 1e3;
num_wells = 100;
mu = 10^-3;
S = 10;

for c = 1:3
    label = num2str(c);
    str = sprintf('%d\n', log10(B0));
    resultsfolder = strcat(num2str(f), '/', label, 'pipette',...
        'N',num2str(num_wells),'S',num2str(S),'B', str,'/');
    
    fp0_sel = zeros(num_cycles, max(S,1));
    B0_sel = zeros(num_cycles, max(S,1));
    varfp0_sel = zeros(num_cycles, max(S,1));
    
    fpT_sel = zeros(num_cycles, max(S,1));
    BT_sel = zeros(num_cycles, max(S,1));
    varfpT_sel = zeros(num_cycles, max(S,1));
    
    fp0_all = zeros(num_cycles, num_wells);
    B0_all = zeros(num_cycles, num_wells);
    varfp0_all = zeros(num_cycles, num_wells);
    
    fpT_all = zeros(num_cycles, num_wells);
    varfpT_all = zeros(num_cycles, num_wells);
    BT_all = zeros(num_cycles, num_wells);
    
    P_all = zeros(num_cycles, num_wells);
    P_sel = zeros(num_cycles, max(S,1));
    
    g0T = 5.5;
    for i = 1 : num_cycles
        load([resultsfolder 'Winnernb' num2str(i)])
        load([resultsfolder 'winner' num2str(i)])
        for j = 1 : max(S,1)
            B0_sel(i, j) = sum(newb_b{j}.*newb_n{j});
            fp0_sel(i, j) = sum(newb_fp{j}.*newb_b{j}.*newb_n{j})/sum(newb_b{j}.*newb_n{j});
            varfp0_sel(i, j) = sum((newb_fp{j}-fp0_sel(i, j)).^2.*newb_b{j}.*newb_n{j})/sum(newb_b{j}.*newb_n{j});
            BT_sel(i, j) = sum(b{j}.*n{j});
            fpT_sel(i, j) = sum(fp{j}.*b{j}.*n{j})/sum(b{j}.*n{j});
            varfpT_sel(i, j) = sum((fp{j}-fpT_sel(i, j)).^2.*b{j}.*n{j})/sum(b{j}.*n{j});
            P_sel(i, j) = P{j}(end);
        end
        
        load([resultsfolder 'nb' num2str(i)])
        load([resultsfolder 'gen' num2str(i)])
        for j = 1 : num_wells
            B0_all(i, j) = sum(newb_b{j}.*newb_n{j});
            fp0_all(i, j) = sum(newb_fp{j}.*newb_b{j}.*newb_n{j})/sum(newb_b{j}.*newb_n{j});
            varfp0_all(i, j) = sum((newb_fp{j}-fp0_all(i, j)).^2.*newb_b{j}.*newb_n{j})/sum(newb_b{j}.*newb_n{j});
            BT_all(i, j) = sum(b{j}.*n{j});
            fpT_all(i, j) = sum(fp{j}.*b{j}.*n{j})/sum(b{j}.*n{j});
            varfpT_all(i, j) = sum((fp{j}-fpT_all(i, j)).^2.*b{j}.*n{j})/sum(b{j}.*n{j});
            P_all(i, j) = P{j}(end);
        end
        
    end
    save([num2str(f) '/' 'fp_data' num2str(c)],'fp0_sel','varfp0_sel','B0_sel',...
        'fpT_sel','varfpT_sel','BT_sel','fp0_all','varfp0_all','B0_all',...
        'fpT_all','varfpT_all','BT_all','P_all','P_sel')
end