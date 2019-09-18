clear
for c = 1:3
    B0 = 10;
    prop = @pipette;
    num_cycles = 3e3;
    num_wells = 100;
    mu = 10^-3;
    label = num2str(c);
    S = 100;
   
    simulateControlSelectionFn(B0,prop,S,num_cycles,...
        num_wells,mu,label)
end