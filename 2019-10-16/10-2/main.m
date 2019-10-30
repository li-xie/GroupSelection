clear
for c = 1:3
    B0 = 1e2;
    prop = @pipette;
    prev_cycles = 1e3;
    num_cycles = 1e3;
    num_wells = 10000;
    mu = 10^-3;
    label = num2str(c);
    S = num_wells/2;
    
   
    simulateGroupSelectionFn(B0,prop,S,num_cycles,prev_cycles,...
        num_wells,mu,label)
end