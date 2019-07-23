mult = 1;
prop = @pipette;
num_cycles = 1e3;
num_wells = 100;
mu = 10^-3;
label = 'mono1';
simulateCommunitySelectionFn(mult,prop,num_cycles,...
    num_wells,mu,label)