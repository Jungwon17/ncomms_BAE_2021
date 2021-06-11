clear; clc;
close all;

tic;

load('cellTable_v4.mat');
load('tag_v5.mat');

path_data = '';
paht_coh = '';

typeName = {'PT','IT','PC'};

%% Data preparation
D = [T.mouseNm,T.cellList, T.cellList, T.hyperLocation];

for ii = 1:size(D,1)
    temp = char(D(ii,3));
    D(ii,3) = temp(end-8:end-4);
    
    temp = char(D(ii,2));
    
    idx_dv = strfind(temp,'_');
    D(ii,2) = temp(idx_dv(3)+1:idx_dv(5)-1);
end

D_PT = D(tag.wsefr & tag.LFP & T.firingRate>0.5,:);
D_IT = D(tag.wsrxfp & tag.LFP & T.firingRate>0.5,:);
D_PC = D(tag.pc & tag.LFP & T.firingRate>0.5,:);

%%
totCoh_cell = cell(4,3);
totCoh_cell_raw = cell(2,2,3);

for type_i = 1:3
    
    type_i

    if type_i == 1
        neuronData = D_PT;
    elseif type_i == 2
        neuronData = D_IT;
    elseif type_i == 3
        neuronData = D_PC;
    end
    
    type_folder = [paht_coh typeName{type_i} '\'];
    
    temp_folder = dir(type_folder);
    neuronN = length(temp_folder)-2;

    totC_correct = nan(61,121,size(neuronData,1));
    totC_error = nan(61,121,size(neuronData,1));
    
    %%%
    totC_correct_raw = nan(61,121,size(neuronData,1));
    totC_correct_cont = nan(61,121,size(neuronData,1));
    
    totC_error_raw = nan(61,121,size(neuronData,1));
    totC_error_cont = nan(61,121,size(neuronData,1));

    
    idxTrialCountDir = nan(size(neuronData,1),4);

    for f_i = 3:length(temp_folder)
        
        disp(['%%%' num2str(f_i-2) ' / ' num2str(neuronN)]);
        
        load([type_folder temp_folder(f_i).name]);
        
        checkChannel = zeros(1,8);
        for c_i = 1:8
            if ~isempty(temp_coh_channel{1,c_i})
                checkChannel(c_i) = 1;
            end
        end
        tarChannel = find(checkChannel);
        
        if isempty(tarChannel)
            disp('Jump');
            continue;
        end
        
        if ~isempty(temp_coh_channel{1,tarChannel(1)})
            time = temp_coh_channel{1,tarChannel(1)}{3,1};
            freq = temp_coh_channel{1,tarChannel(1)}{4,1};
        end
        
        cellNum = temp_folder(f_i).name;
        cellNum(1:2) = []; cellNum(end-3:end) = []; 
        cellNum = str2double(cellNum);
        
        %%
        temp = neuronData(cellNum,:);

        %% Read behav
        tpfolder = char(temp(1));
        cd([path_data '\save_' tpfolder]);

        temp_lfp_name = dir(['LFP' '*' char(temp(2)) '*.mat']);
        load(temp_lfp_name.name);
        
        idxTrialCountDir(cellNum,1) = sum(idx_LC_f);
        idxTrialCountDir(cellNum,2) = sum(idx_RC_f);
        idxTrialCountDir(cellNum,3) = sum(idx_LE_f);
        idxTrialCountDir(cellNum,4) = sum(idx_RE_f);
        
        for correct_i = 1:2
            
            temp_coh = nan(length(freq),length(time),8);
            temp_coh_raw = nan(length(freq),length(time),8);
            temp_coh_cont = nan(length(freq),length(time),8);
            
            
            for c_i = 1:8
                if c_i == channel_i; continue; end
                if checkUseMat(c_i) == 0; continue; end
                
                if isempty(temp_coh_channel{1,c_i}); continue; end %%%% 21.02.26
                
                C_ori = temp_coh_channel{1,c_i}{1,correct_i}';
                C_control = temp_coh_channel{1,c_i}{2,correct_i};
                C_control = permute(C_control,[2 1 3]);
                
                if isempty(C_ori) || isempty(C_control); continue; end
                
                temp_coh(:,:,c_i) = C_ori - nanmean(C_control,3);
                
                temp_coh_raw(:,:,c_i) = C_ori;
                temp_coh_cont(:,:,c_i) = nanmean(C_control,3);
            end
            
            if correct_i == 1
                totC_correct(:,:,cellNum) = nanmean(temp_coh,3);
                
                totC_correct_raw(:,:,cellNum) = nanmean(temp_coh_raw,3);
                totC_correct_cont(:,:,cellNum) = nanmean(temp_coh_cont,3);
                
            elseif correct_i == 2
                totC_error(:,:,cellNum) = nanmean(temp_coh,3);
                
                totC_error_raw(:,:,cellNum) = nanmean(temp_coh_raw,3);
                totC_error_cont(:,:,cellNum) = nanmean(temp_coh_cont,3);
            end
            
        end
        
    end
    
    totCoh_cell{1,type_i} = totC_correct;
    totCoh_cell{2,type_i} = totC_error;    
    
    totCoh_cell{4,type_i} = idxTrialCountDir;
    
   %%%
    totCoh_cell_raw{1,1,type_i} = totC_correct_raw;
    totCoh_cell_raw{1,2,type_i} = totC_correct_cont;
    
    totCoh_cell_raw{2,1,type_i} = totC_error_raw;
    totCoh_cell_raw{2,2,type_i} = totC_error_cont;
end

cd(paht_coh);
save('totC_210329.mat','totCoh_cell','totCoh_cell_raw','freq','time');

toc;
