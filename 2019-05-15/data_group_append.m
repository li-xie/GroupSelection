clear
c = 1;
f = 1;
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
P_all = zeros(num_cycles, num_wells);
P_sel = zeros(num_cycles, max(S,1));
for i = 1 : num_cycles
    load([resultsfolder 'winner' num2str(i)])
    for j = 1 : max(S,1)
        P_sel(i, j) = P{j}(end);
    end
    load([resultsfolder 'gen' num2str(i)])
    for j = 1 : num_wells
        P_all(i, j) = P{j}(end);
    end
    
end
save([num2str(f) '/' 'fp_data'],'P_sel','P_all','-append')