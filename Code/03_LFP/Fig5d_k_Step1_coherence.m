clear; clc;
close all;

warning('off');

tic;

load('cellTable_v4.mat');
load('tag_v5.mat');

idx_randDelay = 0;

save_coh = '';
path_main = ''; %%% Preprocessed data

cd(path_main);

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

%%% 0.5Hz
D_PT1 = D(tag.wsefr & tag.LFP & T.firingRate>0.5,:);
D_IT = D(tag.wsrxfp & tag.LFP & T.firingRate>0.5,:);
D_PC = D(tag.pc & tag.LFP & T.firingRate>0.5,:); %%% 21.02.24 

list_folder = dir();

%% Coherogram parameter
fq = 2; % [0.5 100] Hz

Fs = 1000;

params.Fs=Fs; %Sampling rate
params.tapers=[3 5]; %Taper parameters  %%%%%%%%%%%%%%%%%%%%%%%%% Fixed
params.fpass=[0.5 60]; %Look at 5-60 Hz band
params.trialave = 1;

movingwin=[1 0.1]; %Window to compute spectrum

for type_i = 1:3
    type_i
    
    if type_i == 1
        neuronData = D_PT;
        minThre = 0;
    elseif type_i == 2
        neuronData = D_IT;
        minThre = 0;
    elseif type_i == 3
        neuronData = D_PC;
        minThre = 2;
    end
    
    mkdir([save_coh typeName{1,type_i}]);
    save_coh_type = [save_coh typeName{1,type_i}];
    
    idx_temp = zeros(size(neuronData,1),1);
    
    for neuron_i = 1:size(neuronData,1)
        clearvars totLFP totCell spk_tot lfp_corr spk_corr lfp_err spk_err
        
        temp_coh_channel = cell(1,8);
        
        disp([num2str(neuron_i) ' / ' num2str(size(neuronData,1))]);
        
        temp = neuronData(neuron_i,:);
        
        channel_i = char(temp(3));
        channel_i = str2double(channel_i(3));
        
        temp_folder = char(temp(1));
        
        cd([path_main '\save_' temp_folder]);
        
        temp_lfp_name = dir(['*' char(temp(2)) '*.mat']);
        if isempty(temp_lfp_name); continue; end
        
        %% CheckTrial
        tempdata = matfile(temp_lfp_name(1).name);
        idxDelay = tempdata.delayLength > 4e3;
        
        if idx_randDelay == 0
            totIdx = cat(2,tempdata.idx_LC_f,tempdata.idx_RC_f,tempdata.idx_LE_f,tempdata.idx_RE_f);
            check = prod(sum(totIdx,1)>=minThre);
        else
            totIdx = cat(2,tempdata.idx_LC_r,tempdata.idx_RC_r,tempdata.idx_LE_r,tempdata.idx_RE_r);
            totIdx = totIdx.*repmat(idxDelay,[1 size(totIdx,2)]);
            check = prod(sum(totIdx,1)>=minThre);
        end
        trialCount = sum(totIdx,1);
        
        if check == 0; continue; end
        
        if size(temp_lfp_name,1) == 2
            load(temp_lfp_name(1).name); load(temp_lfp_name(2).name);
        else
            disp('No Files');
        end
        
        %%% Added
        check_load = exist('totCell','var');
        if ~check_load
            continue;
        else
            idx_temp(neuron_i) = 1;
        end
        
        spk_tot = [];
        for cell_i = 1:size(totCell,2)
            if isequal(totCell{1,cell_i}, char(temp(3)))
                spk_tot = totCell{2,cell_i};
            end
        end
        
        if isempty(spk_tot)
            continue;
        end
        
        %% Indexing
        if idx_randDelay == 0
            idx_corr = idx_LC_f | idx_RC_f;
            idx_err = idx_LE_f | idx_RE_f;
        else
            idx_corr = (idx_LC_r | idx_RC_r) &idxDelay;
            idx_err = (idx_LE_r | idx_RE_r) &idxDelay;
        end
        
        %% Firing rates
        tList = [0 2000;2000 4000;4000 8000;8000 13000];
        tLength = [2 2 4 5];
        FRmat = zeros(2,size(tList,1));
        burstMat = zeros(2,size(tList,1));
        burstMat_whole = zeros(2,1);
        for cor_i = 1:2
            if cor_i == 1
                tempspk = spk_tot(idx_corr,:);
            else
                tempspk = spk_tot(idx_err,:);
            end
            
            %%% FR
            for t_i = 1:size(tList,1)
                idxt = tList(t_i,1)+1:tList(t_i,2);
                
                FRmat(cor_i,t_i) = mean(sum(tempspk(:,idxt),2)/tLength(t_i));
                
                %%% Burst idx
                tp = tempspk(:,idxt);
                
                tempburst = nan(1,size(tp,1));
                for kk = 1:size(tp,1)
                    tempburst(kk) = mean(diff(find(tp(kk,:)==1)) < 100);
                end
                
                burstMat(cor_i,t_i) = nanmean(tempburst);
            end
            
            tp = tempspk(:,2001:end);
            tempburst = nan(1,size(tp,1));
            for kk = 1:size(tp,1)
                tempburst(kk) = mean(diff(find(tp(kk,:)==1)) < 100);
            end
            burstMat_whole(cor_i) = nanmean(tempburst);
        end

        parfor c_i = 1:8
            if (checkUseMat(c_i) == 1) && (c_i ~= channel_i)
                
                temp_save = cell(4,2);
                
                lfp_matrix = totLFP{1,c_i}; %%% 20.02.23
                
                %%% lfp
                lfp_corr = lfp_matrix(idx_corr,:,fq);
                spk_corr = spk_tot(idx_corr,:);
                
                checknan = logical(sum(isnan(lfp_corr),2) == size(lfp_corr,2));
                lfp_corr(checknan,:) = []; spk_corr(checknan,:) = [];
                
                lfp_err = lfp_matrix(idx_err,:,fq);
                spk_err = spk_tot(idx_err,:);
                
                checknan = logical(sum(isnan(lfp_err),2) == size(lfp_err,2));
                lfp_err(checknan,:) = []; spk_err(checknan,:) = [];
                
                ns = size(lfp_corr,1);
                ne = size(lfp_err,1);
                
                %% Coherogram -Correct
                if ~isempty(lfp_corr) && ~isempty(spk_corr)
                    
                    lfp_corr = lfp_corr - repmat(nanmean(lfp_corr,2),[1 size(lfp_corr,2)]); %%% minus mean
                    
                    lfp_corr = gpuArray(lfp_corr);
                    spk_corr = gpuArray(spk_corr);
                    
                    sav_i = 1;
                    [C,phi,S12,S1,S2,time,freq] = cohgramc_HS(spk_corr',lfp_corr',movingwin,params);
                    
                    temp_save{1,sav_i} = C;
                    
                    %%% Control
                    for tr_i = 1:size(lfp_corr,1)-1
                        lfp_shift = circshift(lfp_corr, -tr_i, 1);
                        
                        lfp_shift = gpuArray(lfp_shift);
                        
                        [C_shift,phi,S12,S1,S2,time,freq] = cohgramc_HS(spk_corr',lfp_shift',movingwin,params);
                        temp_save{2,sav_i} = cat(3,temp_save{2,sav_i},C_shift);
                    end
                    
                    temp_save{3,sav_i} = time;
                    temp_save{4,sav_i} = freq;
                    temp_save{5,sav_i} = ns;
                end
                
                %% Coherogram -Error
                if ~isempty(lfp_err) && ~isempty(spk_err)
                    
                    lfp_err = lfp_err - repmat(nanmean(lfp_err,2),[1 size(lfp_err,2)]); %%% minus mean
                    
                    lfp_err = gpuArray(lfp_err);
                    spk_err = gpuArray(spk_err);
                    
                    sav_i = 2;
                    [C,phi,S12,S1,S2,time,freq] = cohgramc_HS(spk_err',lfp_err',movingwin,params);
                    
                    temp_save{1,sav_i} = C;
                    %%% Control
                    for tr_i = 1:size(lfp_err,1)-1
                        lfp_shift = circshift(lfp_err, -tr_i, 1);
                        
                        lfp_shift = gpuArray(lfp_shift);
                        
                        [C_shift,phi,S12,S1,S2,time,freq] = cohgramc_HS(spk_err',lfp_shift',movingwin,params);
                        temp_save{2,sav_i} = cat(3,temp_save{2,sav_i},C_shift);
                        temp_save{5,sav_i} = ne;
                    end
                end
                
                temp_coh_channel{1,c_i} = temp_save;
            end
        end
        
        %%% Save results
        tempName = typeName{1,type_i};
        if neuron_i < 10
            tempName = cat(2,tempName,'000');
        elseif neuron_i >= 10 && neuron_i<100
            tempName = cat(2,tempName,'00');
        elseif neuron_i >= 100 && neuron_i<1000
            tempName = cat(2,tempName,'0');
        end
        
        save([save_coh_type '\' tempName num2str(neuron_i) '.mat'],'temp_coh_channel','channel_i','checkUseMat',...,
            'FRmat','trialCount','burstMat','burstMat_whole');
        
    end
end

toc;

