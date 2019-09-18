clear
for c = 1:3
    B0 = 1e2;
    prop = @pipette;
    num_cycles = 1e3;
    num_wells = 100;
    mu = 10^-3;
    label = num2str(c);
    S = 10;
   
    simulateGroupSelectionFn(B0,prop,S,num_cycles,...
        num_wells,mu,label)
end