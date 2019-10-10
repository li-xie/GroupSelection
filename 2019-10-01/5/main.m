clear
for c = 1:3
    B0 = 1e6;
    prop = @pipette;
    num_cycles = 100;
    num_wells = 100;
    mu = 10^-3;
    label = num2str(c);
    S = 0;
   
    simulateGroupSelectionFn(B0,prop,S,num_cycles,...
        num_wells,mu,label)
end