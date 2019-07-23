function loadCommunitySelection(resultsfolder,prop_method,init_gen)

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

load(strcat(resultsfolder,'/run_conditions.mat'));
gc_struct = load(strcat(resultsfolder,'/gc-struct.mat'));
gc = gc_struct.gc;
cw_params = compressWinner;
if mod(cycle_duration,gc.dt) ~= 0
    error('cycle duration must be a multiple of dt')
end

N0 = ic_L_manu' * ic_N_manu + ic_L_help' * ic_N_help;

newb_fp_manu = cell(1,num_wells);
newb_L_manu = cell(1,num_wells);
newb_N_manu = cell(1,num_wells);
newb_L_help = cell(1,num_wells);
newb_N_help = cell(1,num_wells);
newb_n_genos = cell(1,num_wells);

load(strcat(resultsfolder,'/nb',num2str(init_gen),'.mat'));

len_fp_manu_init = length(ic_fp_manu);

wellPlate_fp_manu = {};
wellPlate_L_manu = {};
wellPlate_N_manu = {};
wellPlate_L_help = {};
wellPlate_N_help = {};
wellPlate_Bio_M = {};
wellPlate_Bio_H = {};
wellPlate_R = {};
wellPlate_B = {};
wellPlate_P = {};
wellPlate_n_genos = {};

resultsfolder = strcat(resultsfolder,'_ctd/');
mkdir(resultsfolder(1:end-1));

for gen = init_gen : num_cycles
    rng(seeds(gen));
    tempseeds = randi(10^7,[num_wells,1]);
    parfor i = 1 : num_wells
        rng(tempseeds(i));
        [wellPlate_fp_manu{i},wellPlate_L_manu{i},wellPlate_N_manu{i},...
            wellPlate_L_help{i},wellPlate_N_help{i},wellPlate_Bio_M{i},...
            wellPlate_Bio_H{i},wellPlate_R{i},wellPlate_B{i},wellPlate_P{i},...
            wellPlate_n_genos{i},gc_out]...
            = simulateOneWell(gc,newb_fp_manu{i},newb_L_manu{i},newb_N_manu{i},...
            newb_L_help{i},newb_N_help{i},newb_n_genos{i});
        if ~isequal(gc_out,gc)
            error('You fool! You must never change your global constants')
        end
    end
    savedata(wellPlate_fp_manu, wellPlate_L_manu, wellPlate_N_manu,...
        wellPlate_L_help, wellPlate_N_help, wellPlate_Bio_M,...
        wellPlate_Bio_H, wellPlate_R, wellPlate_B, wellPlate_P, wellPlate_n_genos,...
        gen, resultsfolder, cw_params);
    P_final = zeros(num_wells,1);
    for i = 1 : num_wells
        P_final(i) = wellPlate_P{i}(end);
    end
    P_sorted = sortrows([P_final,(1:num_wells)'],1);
    
    newb_fp_manu = cell(1,num_wells);
    newb_L_manu = cell(1,num_wells);
    newb_N_manu = cell(1,num_wells);
    newb_L_help = cell(1,num_wells);
    newb_N_help = cell(1,num_wells);
    newb_n_genos = cell(1,num_wells);
    
    num_wells_filled = 0;
    finished_pipetting = false;
    for i = 1 : num_wells
        windex = round(P_sorted(end - i + 1,2));
        if cw_params(1)
            [~, lossy_L, lossy_N] = compress(wellPlate_fp_manu{windex},...
                wellPlate_L_manu{windex},wellPlate_N_manu{windex},cw_params(2:3));
            lossy_Bio_M = lossy_L' * lossy_N;
            nD = floor((lossy_Bio_M + wellPlate_Bio_H{windex}(end)) / N0);
        else
            nD = floor((wellPlate_Bio_M{windex}(end) + wellPlate_Bio_H{windex}(end)) / N0);
        end
        target_inds = num_wells_filled + 1 : num_wells_filled + nD;
        if target_inds(end) >= num_wells
            target_inds(target_inds > num_wells) = [];
            finished_pipetting = true;
        end
        num_wells_filled = num_wells_filled + nD;
        [newb_fp_manu,newb_L_manu,newb_N_manu,newb_L_help,newb_N_help,newb_n_genos] = prop_method(...
            newb_fp_manu,newb_L_manu,newb_N_manu,newb_L_help,newb_N_help,newb_n_genos,...
            wellPlate_fp_manu{windex},wellPlate_L_manu{windex},wellPlate_N_manu{windex},...
            wellPlate_L_help{windex},wellPlate_N_help{windex}, wellPlate_n_genos{windex},...
            target_inds,nD,cw_params,len_fp_manu_init,N0);
        
        if finished_pipetting
            break
        end
    end
    save_newborns(newb_fp_manu,newb_L_manu,newb_N_manu,newb_L_help,newb_N_help,...
        newb_n_genos,gen,resultsfolder);
end
end