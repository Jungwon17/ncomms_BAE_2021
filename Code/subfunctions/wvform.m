function wvform(fdorfile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Object: Waveform을 읽어서 분석한다.
% Author: Dohoung Kim
% First written: 2014/09/19
% Last revision: 2015/01/03
% Ver 2.0 (2015/01/03)
%   Function 실행이 용이하게 변경
%   Firing rate는 psth.m으로 이동
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function Input 확인
switch nargin
    case 0
        ttfile = FindFiles('T*.t','CheckSubdirs',0);
        if isempty(ttfile)
            disp('TT file does not exist in this folder!');
            return;
        end
    case 1
        if ~iscell(fdorfile)
            disp('Input argument is wrong. It should be cell array.');
            return;
        elseif isempty(fdorfile)
            ttfile = FindFiles('T*.t','CheckSubdirs',1);
            if isempty(ttfile)
                disp('Cluster does not exist in this folder!');
                return;
            end
        else
            nfolder = length(fdorfile);
            ttfile = cell(0,1);
            for ifolder = 1:nfolder
                if exist(fdorfile{ifolder})==7
                    cd(fdorfile{ifolder});
                    ttfile = [ttfile;FindFiles('T*.t','CheckSubdirs',1)];
                elseif strcmp(fdorfile{ifolder}(end-1:end),'.t');
                    ttfile = [ttfile;fdorfile{ifolder}];
                end
            end
            if isempty(ttfile)
                disp('TT file does not exist!');
                return;
            end
        end
end
nfile = length(ttfile);

for ifile = 1:nfile
    [cellcd,cellname,~] = fileparts(ttfile{ifile});
    ttname = strsplit(cellname,'_');
    if ttname{1}(1)~='T'; ttname{1} = ttname{1}(4:6); end;
    cd(cellcd);
    
%% Input range
    nttfile = fopen([ttname{1},'.ntt']);
    
    volts = fgetl(nttfile);
    while ~strncmp(volts,'-ADBitVolts',11)
        volts = fgetl(nttfile);
    end
    volttemp = strsplit(volts,' ');
    bitvolt = zeros(1,4);
    for ich = 1:4
        bitvolt(ich) = str2num(volttemp{ich+1});
    end
    
%% Waveform
    load([ttname{1},'.clusters'],'-mat','MClust_Clusters');
    spk_idx = FindInCluster(MClust_Clusters{str2num(ttname{2})});
    [~,wv] = LoadTT_NeuralynxNT([ttname{1},'.ntt']);
    cellwv = wv(spk_idx,:,:);
    spkwv = zeros(4,32);
    for itt = 1:4
        spkwv(itt,:) = squeeze(mean(cellwv(:,itt,:)));
    end
    spkwv = (10^6)*diag(bitvolt)*spkwv; % Unit: uV
    
%% Waveform feature
    [~,maintt] = max(max(spkwv'));
    [pkamp,pkidx] = max(spkwv(maintt,:));
    [vlamp,vlidx] = min(spkwv(maintt,pkidx:end));
    vlidx = pkidx + vlidx - 1;
    spkwth = 1000*(vlidx - pkidx)/32; % unit: us
    spkpvr = pkamp/-vlamp;
    hfvl = vlamp/2;
    hfvlfst = find(spkwv(maintt,pkidx:vlidx)<=hfvl,1,'first')+pkidx-1;
    hfvllst = find(spkwv(maintt,vlidx:end)<=hfvl,1,'last')+vlidx-1;
    hfvl1 = (hfvlfst-1) + (spkwv(maintt,hfvlfst-1)-hfvl)/(spkwv(maintt,hfvlfst-1)-spkwv(maintt,hfvlfst));
    if hfvllst<32
        hfvl2 = hfvllst + (hfvl-spkwv(maintt,hfvllst))/(spkwv(maintt,hfvllst+1)-spkwv(maintt,hfvllst));
        hfvwth = 1000*(hfvl2-hfvl1)/32;    
    else
        hfvwth = 1000*(vlidx-hfvl1)*2/32;
    end
        
    save([cellname,'.mat'],...
        'spkwv','spkwth','spkpvr','hfvwth',...
        '-append');
end