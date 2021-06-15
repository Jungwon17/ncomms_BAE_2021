clearvars; clc; close all;
%%
% cross-temporal decoding with SVM method

% Code written by Bae JW, Jeong H 
% 2021 Nat. Commun. [Parallel processing of working memory and temporal
% information by distinct types of cortical projection neurons]
% last edited by Bae JW 2021-06-11

Directory='D:\Backup\code\WM_2021';
cd(Directory); % set your directory here
%% load cell data

cd([Directory,'\Data\IT'])
IT_Neurons = FindFiles('T*.mat','CheckSubdirs',1);
[lvcd,name] = cellfun(@fileparts,IT_Neurons,'UniformOutput',false);

cd([Directory,'\Data\PT'])
PT_Neurons = FindFiles('T*.mat','CheckSubdirs',1);
[lvcd2,name2] = cellfun(@fileparts,IT_Neurons,'UniformOutput',false);
%% cross- temporal decoding

typeList = {'wsIT';'wsPT'};
nT = 2;

nTrain = 9;
nTest = 1;
nIter = 100;

winSize = 1000;
winStep = 100;
win = [0 8000];

outIterNum = [];
outCellNum = [];
if ~exist('score','var')
    [score,scoreChoice] = deal(cell(nT,2));
end

for iT = 1:nT
    if ~isempty(score{iT,1})
        continue; end
    if iT == 1 % IT
        cellList = IT_Neurons;
        cd([Directory,'\Data\IT']);
        load('hyperLoc');
    else
        cellList = PT_Neurons;
        cd([Directory,'\Data\PT']);
        load('hyperLoc');
    end
    
    [spkTotal,spkTotalChoice] = deal(cell(length(cellList),2)); %cell,trial type
    out = false(length(cellList),1);
    
    for iC = 1:length(cellList)
        load(cellList{iC},'spikeTime');
        load([fileparts(cellList{iC}),'\Events.mat'],'trialIndex','trialResult');
        if trialResult(1)<nTrain+nTest | trialResult(4)<nTrain+nTest
            out(iC) = true;
            continue; end
        
        [time,spk] = spikeBin(spikeTime,win,winSize,winStep);
       
        for i = 1:2
            spkTotal{iC,i} = spk(trialIndex(:,1+(i-1)*3),:);
        end
    end
    
    nC(iT) = sum(~out);
    spkTotal(out,:) = [];
    nBin = length(time);

    for iIter = 1:nIter
        trialInd = cellfun(@(x) randsample(size(x,1),nTrain+nTest),spkTotal,'UniformOutput',false);
        accuracy = NaN((nTrain+nTest)/nTest,nBin,nBin);
        data{1} = spkTotal;
    
        for jIter = 1:(nTrain+nTest)/nTest
            fprintf('%s : %d/%d sub-iteration of %d/%d iteration\n',...
                typeList{iT},jIter,(nTrain+nTest)/nTest,iIter,nIter);
            for iB = 1:nBin
                if jIter==1
                    spkTrain = mat2cell(cell2mat(cellfun(@(x,y) x(y(nTest*jIter+1:end),iB),...
                        data{1},trialInd,'UniformOutput',false)'),[nTrain,nTrain],nC(iT));
                elseif jIter==(nTrain+nTest)/nTest
                    spkTrain = mat2cell(cell2mat(cellfun(@(x,y) x(y(1:nTest*(jIter-1)),iB),...
                        data{1},trialInd,'UniformOutput',false)'),[nTrain,nTrain],nC(iT));
                else
                    spkTrain = mat2cell(cell2mat(cellfun(@(x,y) x(y([1:nTest*(jIter-1) jIter*nTest+1:end]),iB),...
                        data{1},trialInd,'UniformOutput',false)'),[nTrain,nTrain],nC(iT));
                end
                
                outCell = std(spkTrain{1})==0 | std(spkTrain{2})==0;
                spkTrain{1}(:,outCell) = [];
                spkTrain{2}(:,outCell) = [];
                mdl = fitcsvm([spkTrain{1};spkTrain{2}],[ones(nTrain,1);ones(nTrain,1)*2]);
                
                for jB = 1:nBin
                    spkTest = mat2cell(cell2mat(cellfun(@(x,y) x(y([1:nTest]+(jIter-1)*nTest),jB),...
                        data{1},trialInd,'UniformOutput',false)'),[nTest nTest],nC(iT));
                    spkTest{1}(:,outCell) = [];
                    spkTest{2}(:,outCell) = [];
                    
                    if sum(outCell)>0
                        outIterNum = [outIterNum;[iIter,jIter]];
                        outCellNum(size(outIterNum,1)) = sum(outCell);
                    end
                    
                    
                    accuracy(jIter,iB,jB) = nanmean(predict(mdl,[spkTest{1};spkTest{2}])==[ones(nTest,1);ones(nTest,1)*2]);
                end
            end
        end
        if isempty(score{iT,1})
            [score{iT,1}] = deal(NaN(nIter,nBin,nBin));
        end
        score{iT,1}(iIter,:,:) = squeeze(nanmean(accuracy,1));
     
    end
end



%% Main figure

for i=1:131
    for j=1:131
        [h,p]=ttest(score{1,2}(:,i,j),0.5);
        pp_IT(i,j)=p;
        hh_IT(i,j)=h;
    end
end


for i=1:131
    for j=1:131
        [h,p]=ttest(score{2,2}(:,i,j),0.5);
        pp_PT(i,j)=p;
        hh_PT(i,j)=h;
    end
end

pp_IT(pp_IT>0.05)=1;
pp_PT(pp_PT>0.05)=1;

pp_IT(squeeze(mean(score{1,2},1))<0.5)=1;
pp_PT(squeeze(mean(score{2,2},1))<0.5)=1;

clr = {[0.745 0.102 0.125];[0.055 0.451 0.725]};

time2=time+7100;
time1=[time,time2];
typeName = {'IT';'PT'};
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 9 5]);
for iT = 1:2
    axes('Position',axpt(2,1,iT,1,axpt(1,1,1,1,axpt(15,15,1:14,1:12)),[0.05 0.05]))
    hold on;
    if iT==1;pp1=pp_IT ;else pp1=pp_PT;end
         imagesc(time,time,squeeze(nanmean(score{iT,2}))');
    colormap(jet); axis xy;
    contour(time,time,pp1',[0.05 0.05],'r','Linewidth',1)
    
    logscale=[log(10^-10) log(0.1)];
    pcorrscale=[0.4 0.8];
    
    set(gca,'CLim',pcorrscale,'XLim',[1000 7000],'YLim',[1000 7000],...
        'XTick',0:2000:11000,'XTickLabel',-2:2:12,'YTick',0:2000:11000,...
        'YTickLabel',-2:2:12,'Box','off','TickDir','out','FontSize',5);
    
     plot([-1000 12000 NaN -1000 12000 NaN -1000 12000],[0 0 NaN 2000 2000 NaN 6000 6000],'k:','LineWidth',0.7);
     plot([0 0 NaN 2000 2000 NaN 6000 6000],[-1000 12000 NaN -1000 12000 NaN -1000 12000],'k:','LineWidth',0.7);
     
     
    if iT==1
         ylabel({'Testing window';'Time from delay onset (s)'},'FontSize',5);
           title(typeName{iT},'FontSize',7,'Color',clr{iT});
    else
        set(gca,'YTickLabel',[]);
        title(typeName{iT},'FontSize',7,'Color',clr{iT});
    end
    xlabel({'Time from delay onset (s)';'Training window'},'FontSize',5);
end
axes('Position',axpt(15,1,15,1))
c = colorbar('east');
set(gca,'CLim',pcorrscale);
c.Position = [0.93 0.5825 0.03 0.3];
c.XTick = pcorrscale;
% c.XTicklabel = [10^-25 1];
c.FontSize = 5;
c.Box ='off';
axis off
colormap('jet');