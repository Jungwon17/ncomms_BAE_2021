%% 
clear; clc;
close all;

tic;
%% para
n = 4; %Controls the order of the filter
fpass=[0.5 100];

% select the raw data folder
path_data = '';
path_save_ori = '';

filenames = dir(path_data);

for f_i = 3:length(filenames)
    
    clc;
    disp([num2str(f_i) ' / ' num2str(length(filenames))]);
    
    temp_folder = filenames(f_i).name;
    
    cd(path_save_ori);
    mkdir(['save_' temp_folder]);
    path_save = [path_save_ori '\save_' temp_folder];
    
    path_neuron = [path_data '\' temp_folder];
    cd(path_neuron);
    
    list_folder = dir();
    for f_ii = 3:length(list_folder)
        
        disp([ '%%' num2str(f_ii) ' / ' num2str(length(list_folder))]);
        
        cd(path_neuron);
        current_folder = list_folder(f_ii).name;
        
        %% 
        if contains(current_folder, 'e2_dv320') || contains(current_folder, 'e9_dv480')...,
                || contains(current_folder, 'r47_dv720') || contains(current_folder, 'r54_dv80')
            continue;
        end
        
        %% Sampling frequency
        Fs_raw = 32552;
        
        if contains(temp_folder, 'rxfp 56') || contains(temp_folder, 'rxfp 57')...,
                || contains(temp_folder, 'rxfp 58') || contains(temp_folder, 'rxfp 59')...,
                || contains(temp_folder, 'rxfp 60') || contains(temp_folder, 'rxfp 62')...,
                || contains(temp_folder, 'rxfp 66') || contains(temp_folder, 'rxfp 67')...,%%% Fixed
                || contains(temp_folder, 'efr 10') || contains(temp_folder, 'efr 11')...,
                || contains(temp_folder, 'efr 13') || contains(temp_folder, 'efr 17')...,%%% Fixed 200202
                || contains(temp_folder, 'efr 18') || contains(temp_folder, 'efr 20')...,%%% Fixed 200202
                || contains(temp_folder, 'efr 21')
            Fs_raw = 1017;
        end
        
        if contains(current_folder,'r57_dv0_f4000')
            Fs_raw = 32552;
        end
        
        if contains(current_folder, 'dv')
            cd([path_neuron '\' current_folder]);
            
            if exist('Events.mat','file') == 0
                continue;
            end
            
            load('Events.mat');
            [n_trial,~]=size(eventTime);
            start_time = 0;
            
            idx_RC_f = logical(trialIndex(:,1)); idx_LC_f = logical(trialIndex(:,4));
            idx_RE_f = logical(trialIndex(:,2)); idx_LE_f = logical(trialIndex(:,5));
            
            idx_corr_f = idx_RC_f | idx_LC_f; idx_error_f = idx_RE_f | idx_LE_f;
            
            idx_RC_r = logical(trialIndex(:,7)); idx_LC_r = logical(trialIndex(:,10));
            idx_RE_r = logical(trialIndex(:,8)); idx_LE_r = logical(trialIndex(:,11));
            
            idx_corr_r = idx_RC_r | idx_LC_r; idx_error_r = idx_RE_r | idx_LE_r;
            
            delayLength = eventTime(:,4)-eventTime(:,3);
            
            totLFP = cell(1,8);
            checkUseMat = zeros(1,8); %%% added 200601
            ratioSaturate = zeros(1,8);
            
            time_ms_save = cell(1,8);
            samples_ms_save = zeros(1,8);
            
            baseT = 10;
            baseLFP = cell(1,8);
            
            parfor ch_i = 1:8
                channelNum = ch_i;
                
                %% CSC -> MAT
                checkFile = exist(['CSC' num2str(channelNum) '.ncs'],'file');
                
                if checkFile ~= 0
                    [Timestamps, Samples, Header]= Nlx2MatCSC(['CSC' num2str(channelNum) '.ncs'], [1 0 0 0 1], 1, 1, []);
                    
                    time_diff = diff(Timestamps);
                    idx_stop = find( time_diff > 1e6 );
                    
                    temp_time = Timestamps;
                    temp_time(idx_stop+1:end) = [];
                    
                    temp_sample = Samples;
                    temp_sample(:,idx_stop+1:end) = [];
                    
                    samples = temp_sample(:);
                    t_temp = (0:511)'*mean(diff(temp_time))/512;
                    tt = repmat(temp_time,[512 1])+repmat(t_temp,[1 length(temp_time)]);
                    times = tt(:);
                    
                    %%% check %%% 210218
                    chekchUse = 0;
                    thre = 32767;
                    
                    if mean(abs(samples) >= thre) < (0.01) %%% 210302
                        chekchUse = 1;
                    end
                    
                    checkUseMat(1,ch_i) = chekchUse;
                    ratioSaturate(1,ch_i) = mean(abs(samples) >= thre);
                    
                    %% Denoising
                    W0 = 60/(Fs_raw/2);  BW = W0/200; %% 21.02.18 
                    [bd,ad] = iirnotch(W0,BW); %% 21.02.18 
                    
                    samples = filter(bd,ad,samples);
                    
                    %%% FFT
                    T = 1/Fs_raw;             % Sampling period
                    L = length(samples);             % Length of signal
                    t = (0:L-1)*T;        % Time vector
                    f = Fs_raw*(0:(L/2))/L;
                    
                    Y = fft(samples);
                    
                    P2 = abs(Y/L); P1 = P2(1:L/2+1);
                    P1(2:end-1) = 2*P1(2:end-1);
                    
                    %%% Check
                    idx1 = logical(f>59.8 & f<60.2);
                    check1 = mean(P1(idx1));
                    
                    idx2 = logical(f>59.6 & f<59.8) | logical(f>60.2 & f<60.4);
                    check2 = mean(P1(idx2));
                    
                    while check1 > check2 * 1.2
                        P1(~idx1) = 0;
                        [~,idxMax] = max(P1);
                        
                        newF = f(idxMax);
                        W0 = newF/(Fs_raw/2);  BW = W0/1000; %% 21.03.02 
                        
                        [bd,ad] = iirnotch(W0,BW); %% 21.02.18 
                        
                        samples = filter(bd,ad,samples);
                        
                        %%% FFT
                        Y = fft(samples);
                        
                        P2 = abs(Y/L); P1 = P2(1:L/2+1);
                        P1(2:end-1) = 2*P1(2:end-1);
                        
                        %%% Check
                        idx1 = logical(f>59.8 & f<60.2);
                        check1 = mean(P1(idx1));
                        
                        idx2 = logical(f>59.6 & f<59.8) | logical(f>60.2 & f<60.4);
                        check2 = mean(P1(idx2));
                    end
                    
                    %% Downsampling
                    samples_deci = decimate(samples,fix(Fs_raw/1017),'FIR'); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% decimate 적용 입니다
                    times_deci = decimate(times,fix(Fs_raw/1017),'FIR');
                    
                    times_shifted = round((times_deci - start_time)/1e3); %% ms unit
                    
                    t_idx = logical(diff(times_shifted));
                    
                    if t_idx(end) == 0
                        t_idx = ([t_idx;0]) & logical(times_shifted>0);
                    elseif t_idx(end) == 1
                        t_idx = ([t_idx;1]) & logical(times_shifted>0);
                    end
                    
                    times_ms = times_shifted(t_idx);
                    samples_ms = samples_deci(t_idx); %% Fixed 191010
                    
                    time_ms_save{1,ch_i} = times_ms;
                    samples_ms_save(ch_i) = length(samples_ms);
                    
                    Fs = 1000;
                    samples_ms=samples_ms.*(1000/32767); %%% scaling
                    
                    %% Filtering
                    sample_filtered = nan(length(samples_ms), 2);
                    sample_filtered(:,1) = samples_ms;
                    
                    for filt_i=1:size(sample_filtered,2)-1
                        [b,a]=butter(n,fpass(filt_i,:)*2/Fs); % bandpass
                        sample_filtered(:,filt_i+1)=filtfilt(b,a,sample_filtered(:,1));
                    end
                    
                    %% data reshape
                    totT = 13000; % 4000 + 4000 + 5000
                    lfp_matrix = nan(n_trial,totT,size(sample_filtered,2));
                    for trial_i = 1:n_trial
                        
                        trial_info = eventTime(trial_i,:);
                        
                        [~, sample_on] = min( abs(times_ms - trial_info(1)) );
                        idx_trial = (sample_on-2000+1):(sample_on+2000+4000+5000);
                        if min(idx_trial) < 1 || max(idx_trial) > size(sample_filtered,1)
                            continue;
                        end
                        
                        lfp_matrix(trial_i,:,:) = sample_filtered(idx_trial,:); %%% 4 + 4 + 5
                    end
                    
                    %%
                    tp = lfp_matrix(:,:,2);
                    
                    tpval = std(tp,0,2)
                    valThre = nanmean(tpval) + 3 * nanstd(tpval);
                    
                    idx_del = tpval >= valThre;
                    
                    lfp_matrix(idx_del,:,:) = nan;
                    
                    totLFP{1,ch_i} = lfp_matrix;
                    
                end
            end
            
            cd(path_save);
            save(['LFP_' current_folder],'totLFP','idx_RC_f','idx_LC_f','idx_RE_f','idx_LE_f',...,
                'idx_RC_r','idx_LC_r','idx_RE_r','idx_LE_r','checkUseMat','delayLength','ratioSaturate');
            
            %% Tar channel
            checkChannel = zeros(1,8);
            for c_i = 1:8
                if ~isempty(totLFP{1,c_i})
                    checkChannel(c_i) = 1;
                end
            end
            tarChannel = find(checkChannel);
            
            times_ms = time_ms_save{1,tarChannel(1)};
            
            %% Spike data
            cd([path_neuron '\' current_folder]);
            
            celldata = dir('TT*.mat');
            tname = {celldata.name};
            t=strtok(tname,'_');
            
            nt=nan(1,8);
            for c_i=1:length(nt)
                nt(c_i)=sum(contains(t,mat2str(c_i)));
            end
            
            %%% cell anal
            n_cell = sum(nt);
            
            tname_cell = cell(1,n_cell);
            for jj = 1:length(tname)
                tname_cell{1,jj} = tname{1,jj}(1:end-4);
            end
            
            totCell = cell(2,n_cell);
            
            totT = 13000; % 4000 + 4000 + 5000
            for cell_i = 1:n_cell
                load([ tname_cell{1,cell_i}  '.mat']);
                
                spk_train = nan(n_trial,totT);
                
                temp_spk = zeros(1,length(times_ms));
                for trial_i = 1:n_trial
                    temp_idx = spikeTime{trial_i,1} + eventTime(trial_i,1);
                    idx_spk = ismember(times_ms,round(temp_idx));
                    
                    temp_spk(idx_spk) = 1;
                end
                
                for trial_i = 1:n_trial
                    trial_info = eventTime(trial_i,:);
                    
                    [~, sample_on] = min( abs(times_ms - trial_info(1)) );
                    idx_trial = (sample_on-2000+1):(sample_on+2000+4000+5000);
                    
                    if min(idx_trial) < 1 || max(idx_trial) > samples_ms_save(tarChannel(1))
                        continue;
                    end
                    
                    spk_train(trial_i,:) = temp_spk(idx_trial); %%% 2 + 4 + 2
                end
                totCell{1,cell_i} = tname_cell{1,cell_i};
                totCell{2,cell_i} = spk_train;
            end
            
            cd(path_save);
            save(['Spike_' current_folder],'totCell','nt');
            
            clearvars temp_spk spk_timing spk_timing_ms totCell;
            
        end
    end
end

toc;