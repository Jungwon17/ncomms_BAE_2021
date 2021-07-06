clear; clc;
close all;

path_spec = ''; path where data saved

cd(path_spec);
filelist = dir();

%% Main
totS = cell(1,3);
idxTrialCountDir = [];

for f_i = 3:length(filelist)
    
    clc;
    disp([num2str(f_i-2) ' / ' num2str(length(filelist)-2)]);
    
    if ~contains(filelist(f_i).name,'save'); continue; end
    load(filelist(f_i).name);
    
    for session_i = 1:length(totSpec)    
        disp(['%%%' num2str(session_i) ' / ' num2str(length(totSpec))]);
        
        checkChannel = zeros(1,8);
        for c_i = 1:8
            if ~isempty(totSpec{1,session_i}{1,c_i}); checkChannel(c_i) = 1; end
        end
        tarChannel = find(checkChannel);
        if isempty(tarChannel); continue; end
        
        check = 0; cidx = 1;
        while check == 0
            temp_channel_use = totSpec{1,session_i}{1,tarChannel(cidx)}{4,1};
            if isempty(temp_channel_use)
                cidx = cidx + 1;
            else
                check = 1;
            end
        end
        
        time_spec = totSpec{1,session_i}{1,tarChannel(cidx)}{2,1};
        freq_spec = totSpec{1,session_i}{1,tarChannel(cidx)}{3,1};
        
        tempS_correct = nan(length(freq_spec),length(time_spec),8);
        tempS_error = nan(length(freq_spec),length(time_spec),8);
        
        minThre = 2;
        for c_i = 1:8
            if temp_channel_use(c_i) == 0; continue; end
            if prod(totSpec{2,session_i}>=minThre) == 0; continue; end
            
            %%% Correct
            sav_i = 1;
            channelS_correct = totSpec{1,session_i}{1,c_i}{1,sav_i};
            
            if ~isempty(channelS_correct)
                data = channelS_correct;
                data = 10*log10(data);
                data_norm = data;
                
                tempS_correct(:,:,c_i) = nanmean(data_norm,3)';
            end
            
            %%% Error
            sav_i = 2;
            channelS_error = totSpec{1,session_i}{1,c_i}{1,sav_i};
            
            if ~isempty(channelS_error)
                data = channelS_error;
                data = 10*log10(data);
                data_norm = data;
                
                tempS_error(:,:,c_i) = nanmean(data_norm,3)';
            end
        end
        
        idxTrialCountDir = cat(1,idxTrialCountDir,totSpec{2,session_i});
        
        totS{1,1} = cat(3,totS{1,1},nanmean(tempS_correct,3));
        if sum(totSpec{2,session_i}(3:4)) > 0
            totS{1,2} = cat(3,totS{1,2},nanmean(tempS_error,3));
        else
            totS{1,2} = cat(3,totS{1,2},nan(size(nanmean(tempS_correct,3))));
        end
    end
end