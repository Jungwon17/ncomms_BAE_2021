clearvars; clc; close all;
%%
% population decoding in slide with SVM method
% last edited by JWBAE 2021-06-11

load('D:\Backup\code\WM_2021\Data\Untagged\Untagged_WS.mat');

%% 
typeName = {'Untagged-WS'};
nT = 1;
nEnsemble = [365];

nTrain = 9;
nTest = 1;
nIter = 100;67

winSize = 1000;
winStep = 100;
win = [-2000 12000];
%% decoding
for Cellnum=1
[accuracyCorrect, accuracyError] = deal(NaN(nIter,(nTrain+nTest)/nTest,nBin));
for iIter = 1:nIter
    trialInd = cellfun(@(x) randsample(size(x,1),nTrain+nTest),spkTotal(:,1:2),'UniformOutput',false);
    
    data = spkTotal;
    
    for jIter = 1:(nTrain+nTest)/nTest
        fprintf('%s,# %d : %d/%d sub-iteration of %d/%d iteration\n',...
            typeName{1},Cellnum,jIter,(nTrain+nTest)/nTest,iIter,nIter);
        
        parfor iB = 1:size(data{1,1},2)
            
            if jIter==1
                spkTrain = mat2cell(cell2mat(cellfun(@(x,y) x(y(nTest*jIter+1:end),iB),...
                    data(:,1:2),trialInd,'UniformOutput',false)'),[nTrain,nTrain],nEnsemble(Cellnum));
            elseif jIter==(nTrain+nTest)/nTest
                spkTrain = mat2cell(cell2mat(cellfun(@(x,y) x(y(1:nTest*(jIter-1)),iB),...
                    data(:,1:2),trialInd,'UniformOutput',false)'),[nTrain,nTrain],nEnsemble(Cellnum));
            else
                spkTrain = mat2cell(cell2mat(cellfun(@(x,y) x(y([1:nTest*(jIter-1) jIter*nTest+1:end]),iB),...
                    data(:,1:2),trialInd,'UniformOutput',false)'),[nTrain,nTrain],nEnsemble(Cellnum));
            end
            % correct
            spkTest = mat2cell(cell2mat(cellfun(@(x,y) x(y([1:nTest]+(jIter-1)*nTest),iB),...
                data(:,1:2),trialInd,'UniformOutput',false)'),[nTest nTest],nEnsemble(Cellnum));
            % error
            trialInd2 = cellfun(@(x) randsample(size(x,1),nTest),spkTotal(:,3:4),'UniformOutput',false);
            spkTesterr = mat2cell(cell2mat(cellfun(@(x,y) x(y(:),iB),...
                data(:,3:4),trialInd2,'UniformOutput',false)'),[nTest nTest],nEnsemble(Cellnum));
            
            outCell = std(spkTrain{1})==0 | std(spkTrain{2})==0;
            spkTrain{1}(:,outCell) = [];
            spkTrain{2}(:,outCell) = [];
            spkTest{1}(:,outCell) = [];
            spkTest{2}(:,outCell) = [];
            spkTesterr{1}(:,outCell) = [];
            spkTesterr{2}(:,outCell) = [];
            
            mdl = fitcsvm([spkTrain{1};spkTrain{2}],[ones(nTrain,1);ones(nTrain,1)*2]);
            
            
            accuracyCorrect(iIter,jIter,iB) = nanmean(predict(mdl,[spkTest{1};spkTest{2}])==[ones(nTest,1);ones(nTest,1)*2]);
            accuracyError(iIter,jIter,iB) = nanmean(predict(mdl,[spkTesterr{1};spkTesterr{2}])==[ones(nTest,1);ones(nTest,1)*2]);
            
        end
    end
end
scoreCorrect{Cellnum} = squeeze(nanmean(accuracyCorrect,2));
scoreError{Cellnum} = squeeze(nanmean(accuracyError,2));

end


%% Correct + Error
for ii=1:length(time)
hh(ii)=ttest(scoreCorrect{3,1}(:,ii),scoreError{1,3}(:,ii),0.05);
hh_it(ii)=ttest(scoreCorrect{3,1}(:,ii),0.5,'Alpha',0.05);
hh_pt(ii)=ttest(scoreError{1,3}(:,ii),0.5,'Alpha',0.05);
end

hh(hh==0)=NaN;hh_it(hh_it==0)=NaN;hh_pt(hh_pt==0)=NaN;

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 5 2.5]);
i = 1
h(1) = axes('Position',axpt(1,1,1,1,axpt(15,15,2.5:15,1:12),[0.025 0.05]));
hold on;
plot([0 0 NaN 2000 2000 NaN 6000 6000],[0.2 2 NaN 0.2 2 NaN 0.2 2],'k:','LineWidth',0.35);
plot([-1000 12000],[0.5 0.5],'k:','LineWidth',0.35);
for j = 3
    ind = j+(i-1)*2;
    fill([time flip(time)],[mc{ind}+sc{ind} flip(mc{ind}-sc{ind})],clr{ind},'EdgeColor','none');
    plot(time,mc{ind},'LineWidth',0.5,'Color',clr{ind});
     fill([time flip(time)],[me{ind}+se{ind} flip(me{ind}-se{ind})],clr{ind},'EdgeColor','none');
    plot(time,me{ind},'--','LineWidth',0.5,'Color',clr{ind});
end
alpha(0.2);
ylabel('P_c_o_r_r_e_c_t','FontSize',5);
plot(time,(hh_it)*1.08,'Color',clr{3}+[0.2 0.2 0.2],'LineWidth',1);
plot(time,(hh_pt)*1.06,'Color',clr{3}+[0.4 0.4 0.4],'LineWidth',1);
plot(time,(hh)*1.1,'Color',clr{3},'LineWidth',1);

 set(h,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35,'XTick',...
    0:2000:12000,'XTickLabel',-2:2:12,'YTick',[0.5 0.75 1],'yticklabel',[50 75 100],'XLim',[-800 11000],'ylim',[0.3 1.1]);
xlabel('Time from delay onset (s)','FontSize',6);


