function simulateGroupSelectionFn(B0,prop_method,S,num_cycles,...
    num_wells,mu,label)

% clear
% B0 = 10;
% prop_method = @pipette;
% num_cycles = 4e2;
% num_wells = 100;
% mu = 10^-3;
% label = num2str(1);
% S = 100;

if isequal(prop_method,@pipette)
    prop_name = 'pipette';
else
    if isequal(prop_method,@pipette_const_N)
        prop_name = 'const_N';
    else
        if isequal(prop_method,@pipette_const_phiM)
            prop_name = 'const_phiM';
        else
            if isequal(prop_method,@pipette_cNphiM)
                prop_name = 'cNphiM';
            else
                error('you must enter a valid propogation method')
            end
        end
    end
end
% incorporate the simulation condition within the folder name
str = sprintf('%d\n', log10(B0));
resultsfolder = strcat(label, prop_name,...
    'N',num2str(num_wells),'S',num2str(S),'B', str,'/');
[~,~,~] = mkdir(resultsfolder(1:end-1));


cycle_exponent = 5.5; % the exponent at the maximal growth rate 
compressWinner = 1; % 1 to compress, 0 not
digits_fp = 5; % number of digits after the decimal point of fp
digits_b = 2; % number of digits after the decimal point of individual biomass
cw_params = [compressWinner, digits_fp, digits_b];
gc.dt = 1e-2; % time step of the simulation
gc.nsteps = round(cycle_exponent/gc.dt); % number of time steps
gc.b_Mmax = 1; % max birth rate
gc.d_M = 5e-3;
gc.p_mut = mu; % mutation rate, probability per cell per division
% gc.frac_null = 0.5;
gc.sp0 = 0.05; % parameter for the mutation spectrum
gc.sn0 = 0.05; % parameter for the mutation spectrum
gc.g = 0;
gc.fp_max = 1; % upper bound of fp
gc.B0 = B0; %B_target
if mod(cycle_exponent, gc.dt) ~= 0
    error('cycle duration must be a multiple of dt')
end
time = 0 : gc.dt : cycle_exponent;
rng('shuffle');
seeds = randi(2^32-1,[num_cycles,1],'uint32');

ic_fp = [0.13;zeros(100 - 1,1)]; % initial distribution of fp
ic_b = [1;zeros(100 - 1,1)]; % initial distribution of individual biomass 
ic_n = [B0 ;zeros(100 - 1,1)]; % initial number of cells
ic_n_genos = 1; % initial number of genotype
% B0 = gc.B0; % B_target
% save simulation conditions 
save_run_conditions(cycle_exponent, num_cycles, num_wells,...
    cw_params, gc, time, ic_fp, ic_b, ic_n,...
    ic_n_genos, resultsfolder, seeds); 

newb_fp = cell(num_wells, 1); % list of fp within each Newborn
newb_b = cell(num_wells, 1); % list of biomass within each Newborn
newb_n = cell(num_wells, 1); % corresponding to each fp and biomass, integer number of cells.
newb_n_genos = cell(num_wells, 1); % length of the list
newb_pnum = zeros(num_wells, 1); % parent number

% initialize the simulation with identical groups
for i = 1 : num_wells
    newb_fp{i} = ic_fp;
    newb_b{i} = ic_b;
    newb_n{i} = ic_n;
    newb_n_genos{i} = ic_n_genos;
end

len_fp_init = length(ic_fp);

% data for all the wells on a Plate within one cycle
wellPlate_fp = cell(num_wells, 1);
wellPlate_b = cell(num_wells, 1);
wellPlate_n = cell(num_wells, 1);
wellPlate_B = cell(num_wells, 1);
wellPlate_P = cell(num_wells, 1);
wellPlate_n_genos = cell(num_wells, 1);

save_newborns(newb_fp,newb_b,newb_n,newb_n_genos,0,0,...
    resultsfolder);
for gen = 1 : num_cycles
    rng(seeds(gen));
    tempseeds = randi(2^32-1,[num_wells,1], 'uint32');
    parfor i = 1 : num_wells
        rng(tempseeds(i));
        [wellPlate_fp{i},wellPlate_b{i},wellPlate_n{i},...
            wellPlate_B{i}, wellPlate_n_genos{i}, wellPlate_P{i}, gc_out]...
            = simulateOneWell(gc,newb_fp{i},newb_b{i},newb_n{i},...
              newb_n_genos{i});
        if ~isequal(gc_out,gc)
            error('You fool! You must never change your global constants')
        end
    end
    scurr = rng;
    % save all the Adults
    savePlate(wellPlate_fp, wellPlate_b, wellPlate_n,...
         wellPlate_B, wellPlate_P, wellPlate_n_genos,...
         gen, resultsfolder, cw_params, scurr);
    
    P_final = zeros(num_wells, 1);
    B_final = zeros(num_wells, 1);
    for i = 1 : num_wells
        P_final(i) = wellPlate_P{i}(end);
        B_final(i) = wellPlate_B{i}(end);
    end
    [~, P_ind] = sort(P_final,'descend');
    B_sorted = B_final(P_ind);
    % reproduce strategy: 
    % 0: top-dog
    % nonzero interger S: select the top m, each reproduce N/S
    rep_profile = ReproduceStrategy(num_wells, B_sorted, B0, S);
    % Sp: the actual number of Adults that reproduced.
    Sp = nnz(rep_profile);
    winner_fp = cell(Sp, 1);
    winner_b = cell(Sp, 1);
    winner_n = cell(Sp, 1);
    winner_n_genos = cell(Sp, 1);
    winner_B = cell(Sp, 1);
    winner_P = cell(Sp, 1);
    winnernb_fp = newb_fp(P_ind(1:Sp));
    winnernb_b = newb_b(P_ind(1:Sp));
    winnernb_n = newb_n(P_ind(1:Sp));
    winnernb_n_genos = newb_n_genos(P_ind(1:Sp));
    winnernb_pnum = newb_pnum(P_ind(1:Sp));
    
    newb_fp = cell(num_wells, 1);
    newb_b = cell(num_wells, 1);
    newb_n = cell(num_wells, 1);
    newb_n_genos = cell(num_wells, 1);
    newb_pnum = zeros(num_wells, 1);
    num_wells_filled = 0;
    for i = 1 : Sp
        winner_fp{i} = wellPlate_fp{P_ind(i)};
        winner_b{i} = wellPlate_b{P_ind(i)};
        winner_n{i} = wellPlate_n{P_ind(i)};
        winner_n_genos{i} = wellPlate_n_genos{P_ind(i)};
        winner_B{i} = wellPlate_B{P_ind(i)};
        winner_P{i} = wellPlate_P{P_ind(i)};
        target_inds = num_wells_filled + 1 : num_wells_filled + rep_profile(i);
        
        num_wells_filled = num_wells_filled + rep_profile(i);
        newb_pnum(target_inds) = i;
        [newb_fp,newb_b,newb_n, newb_n_genos,] = prop_method(...
            newb_fp,newb_b,newb_n,newb_n_genos,...
            wellPlate_fp{P_ind(i)},wellPlate_b{P_ind(i)},wellPlate_n{P_ind(i)},...
            wellPlate_n_genos{P_ind(i)},target_inds,floor(B_sorted(i)/B0),cw_params,len_fp_init,B0);
    end
%     % save only the winners
%     saveWinner(winner_fp, winner_b, winner_n,...
%         winner_B, winner_P, winner_n_genos,...
%         gen, resultsfolder, cw_params, scurr);
    
%     saveWinnerNb(winnernb_fp, winnernb_b, winnernb_n,...
%         winnernb_n_genos, winnernb_pnum, gen, resultsfolder);
%     
    save_newborns(newb_fp,newb_b,newb_n,...
        newb_n_genos,newb_pnum, gen, resultsfolder);
end
% end
