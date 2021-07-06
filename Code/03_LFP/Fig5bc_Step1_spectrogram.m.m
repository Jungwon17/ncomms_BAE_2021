clear; clc;
close all;

%% Set paths
path_main = ''; % path where data saved
save_spectrogram = ''; % path where results be saved

%% Data preparation
load('cellTable_v4.mat');
load('tag_v5.mat');

idx_randDelay = 0;
cd(path_main);

%% Param
fq = 2;
Fs = 1000;

%%% Need Chronux
params.Fs=Fs; %Sampling rate
params.tapers=[3 5]; %Taper parameters
params.fpass=[0.5 100]; %Look at 5-60 Hz band
params.trialave = 0;

movingwin=[1 0.1]; %Window to compute spectrum

filenames = dir(path_main);
for f_i = 3:length(filenames)
    
    disp([num2str(f_i-2) ' / ' num2str(length(filenames)-2)]);
    
    temp_folder = filenames(f_i).name;
    cd([path_main '\' temp_folder]);
    
    temp_lfp_name = dir('LFP*.mat');
    
    totSpec = cell(2,length(temp_lfp_name));
    for f_ii = 1:length(temp_lfp_name)
        
        disp([ '%%' num2str(f_ii) ' / ' num2str(length(temp_lfp_name))]);
        load(temp_lfp_name(f_ii).name);
        if sum(checkUseMat) == 0; disp(temp_lfp_name(f_ii).name); end

        %% Indexing
        idxDelay = delayLength > 4e3;
        if idx_randDelay == 0
            idx_corr = idx_LC_f | idx_RC_f;
            idx_err = idx_LE_f | idx_RE_f;
            
            tempcount = [sum(idx_LC_f) sum(idx_RC_f) sum(idx_LE_f) sum(idx_RE_f)];
        else
            idx_corr = (idx_LC_r | idx_RC_r) &idxDelay;
            idx_err = (idx_LE_r | idx_RE_r) &idxDelay;
            
            tempcount = [sum(idx_LC_r&idxDelay) sum(idx_RC_r&idxDelay) sum(idx_LE_r&idxDelay) sum(idx_RE_r&idxDelay)];
        end
        
        %% Spectrogram
        temp_spec_channel = cell(1,8);
        tempLFPsize = cell(1,8);
        controlSpec = cell(1,8);
        
        parfor c_i = 1:8
            temp_save = cell(5,2);
            
            lfp_matrix = totLFP{1,c_i};
            
            if checkUseMat(c_i) == 1
                lfp_corr = lfp_matrix(idx_corr,:,fq);
                
                checknan_s = logical(sum(isnan(lfp_corr),2) == size(lfp_corr,2));
                lfp_corr(checknan_s,:) = [];
                lfp_err = lfp_matrix(idx_err,:,fq);
                checknan_e = logical(sum(isnan(lfp_err),2) == size(lfp_err,2));
                lfp_err(checknan_e,:) = [];
                
                if idx_randDelay == 0
                    idxtemp = nan(1,4);
                    idxtemp(1) = sum(idx_LC_f(idx_corr) & ~checknan_s);
                    idxtemp(2) = sum(idx_RC_f(idx_corr) & ~checknan_s);
                    idxtemp(3) = sum(idx_LE_f(idx_err) & ~checknan_e);
                    idxtemp(4) = sum(idx_RE_f(idx_err) & ~checknan_e);
                else
                    idxtemp = nan(1,4);
                    idxtemp(1) = sum(idx_LC_r(idx_corr) & ~checknan_s);
                    idxtemp(2) = sum(idx_RC_r(idx_corr) & ~checknan_s);
                    idxtemp(3) = sum(idx_LE_r(idx_err) & ~checknan_e);
                    idxtemp(4) = sum(idx_RE_r(idx_err) & ~checknan_e);
                end
                tempLFPsize{1,c_i} = idxtemp;
                
                %% Correct trials
                if ~isempty(lfp_corr)
                    templfp = lfp_corr;
                    templfp = templfp - repmat(nanmean(templfp,2),[1 size(templfp,2)]);
                    
                    sav_i = 1;
                    [S,time_spec_n,freq_spec_n]=mtspecgramc(templfp',movingwin,params);

                    temp_save{1,sav_i} = S;
                    temp_save{2,sav_i} = time_spec_n;
                    temp_save{3,sav_i} = freq_spec_n;
                    temp_save{4,sav_i} = checkUseMat;
                end
                
                %% Error trial
                if ~isempty(lfp_err)
                    templfp = lfp_err;
                    templfp = templfp - repmat(nanmean(templfp,2),[1 size(templfp,2)]);

                    sav_i = 2;
                    [S,time_spec_n,freq_spec_n]=mtspecgramc(templfp',movingwin,params);
                    
                    temp_save{1,sav_i} = S;
                end
    
                %% Save
                temp_spec_channel{1,c_i} = temp_save;
            end
        end
        
        totSpec{1,f_ii} = temp_spec_channel;
        totSpec{2,f_ii} = tempcount;
        totSpec{3,f_ii} = tempLFPsize;
    end
    
    save([save_spectrogram temp_folder '_Spec.mat'],'totSpec','params');
end