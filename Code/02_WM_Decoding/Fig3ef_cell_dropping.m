clearvars; clc;
%%
% cell dropping with SVM method
% last edited by JWBAE 2021-06-11

Directory='D:\Backup\code\WM_2021';
cd(Directory); % set your directory here

%% load cell data

cd([Directory,'\Data\IT'])
IT_Neurons = FindFiles('T*.mat','CheckSubdirs',1);
[lvcd,name] = cellfun(@fileparts,IT_Neurons,'UniformOutput',false);

cd([Directory,'\Data\PT'])
PT_Neurons = FindFiles('T*.mat','CheckSubdirs',1);
[lvcd2,name2] = cellfun(@fileparts,IT_Neurons,'UniformOutput',false);

%% cell dropping decoding
nT = 2;
nC = NaN(nT,1);

nTrain = 9;
nTest = 1;
nIter = 100;
nEnsemble = 50;
typeList = {'IT';'PT';'PP'};

score = cell(nT,1);
for iT = 1:nT
    if iT == 1 % IT
        cellList = IT_Neurons;
        cd([Directory,'\Data\IT']);
        load('hyperLoc');
    else
        cellList = PT_Neurons;
        cd([Directory,'\Data\PT']);
        load('hyperLoc');
    end
    
    nC(iT) = length(cellList);
    
    out = false(nC(iT),1);
    spkTotal = cell(nC(iT),2);
    for iC = 1:nC(iT)
        load(cellList{iC},'spikeTime','spikeTimeChoice');
        load([fileparts(cellList{iC}),'\Events.mat'],'trialIndex','trialResult');
        
        if trialResult(1)<nTrain+nTest | trialResult(4)<nTrain+nTest
            out(iC) = true; continue;
        end
        
        [~,spk] = spikeBin(spikeTime,[2000 6000],4000,4000);
        trialInd = [trialIndex(:,1) trialIndex(:,4)];
        for iI = 1:2
            spkTotal{iC,iI} = spk(trialInd(:,iI));
        end
    end
    spkTotal(out,:) = [];
    nC(iT) = sum(~out);
    accuracy = NaN(nEnsemble,nIter,nTrain+nTest);
    for iE = 1:nEnsemble
        if iE>nC(iT)
            break;
        end
        iE
        for iIter = 1:nIter
            fprintf('%s : %d/%d iteration of %d/%d ensemble\n',...
                typeList{iT},iIter,nIter,iE,min([nC(iT),nEnsemble]));
            ensembleInd = randsample(nC(iT),iE);
            spkIter = cellfun(@(x) x(randsample(size(x,1),nTrain+nTest)),...
                spkTotal(ensembleInd,:),'UniformOutput',false);
            for jIter = 1:nTrain+nTest
                if jIter==1
                    trainInd = 2:nTrain+nTest;
                elseif jIter==nTrain+nTest
                    trainInd = 1:nTrain;
                else
                    trainInd = [1:jIter-1, jIter+1:nTrain+nTest];
                end
                
                spkTrain = cell2mat(cellfun(@(x) x(trainInd)',spkIter,...
                    'UniformOutput',false))';
                spkTest = cellfun(@(x) x(jIter),spkIter)';
                
                zerovar = [];
                for iI = 1:2
                    zerovar = [zerovar;...
                        find(var(spkTrain([1:nTrain]+(iI-1)*nTrain,:))==0)];
                end
                spkTrain(:,zerovar) = [];
                spkTest(:,zerovar) = [];
                if isempty(spkTrain)
                    continue; end
                
                mdl = fitcsvm(spkTrain,[ones(nTrain,1);ones(nTrain,1)*2]);
                accuracy(iE,iIter,jIter) = nanmean(predict(mdl,spkTest)==...
                    [ones(nTest,1);ones(nTest,1)*2]);
            end
        end
    end
    score{iT} = squeeze(nanmean(accuracy,3));
end


%% Main -figure

nEnsemble = 50;
typeList = {'IT';'PT';'PP'};
m = cellfun(@(x) nanmean(x,2),score,'UniformOutput',false);
s = cellfun(@(x) nanstd(x,[],2)/sqrt(nIter),score,'UniformOutput',false);
clr = {[0.745 0.102 0.125];[0.055 0.451 0.725];[0.2 0.2 0.2]};
nC=[17 14 50];

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 5 5]);
axes('Position',axpt(1,1,1,1,[0.3 0.26 0.65 0.70],[0.05 0.075])); hold on;

hold on;
plot([1 nEnsemble],[0.5 0.5],'--','LineWidth',0.35,'Color',[.7 .7 .7]);
for iT = 1:2
    n = min([nC(iT),nEnsemble]);
    fill([1:n, n:-1:1],[m{iT}(1:n)+s{iT}(1:n); flip(m{iT}(1:n)-s{iT}(1:n))],...
        clr{iT},'EdgeColor','none');
    plot(1:nEnsemble,m{iT},'Color',clr{iT},'LineWidth',1);
    plot([25 30],repmat(0.57-(iT-1)*0.025,1,2),'Color',clr{iT},'LineWidth',1);
    text(32,0.57-(iT-1)*0.025,typeList{iT},'FontSize',5);
end
set(gca,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35,...
    'XTick',[1 10 20 30 40 50],'YTick',0.5:0.1:0.7,'XLim',[1 35],...
    'YLim',[0.45 0.7])
alpha(0.2);
xlabel('Number of neurons','FontSize',5);
ylabel('P_c_o_r_r_e_c_t','FontSize',5);

%%
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3 3]);
hold on;
plot([1 nEnsemble],[0.5 0.5],'k:','LineWidth',0.35);
fill([1:nEnsemble flip(1:nEnsemble)],[m{3}+s{3}; flip(m{3}-s{3})],clr{3},...
    'EdgeColor','none');
fill([2:2:50 flip(2:2:50)],[m{1}(1:25)+s{1}(1:25); flip(m{1}(1:25)-s{1}(1:25))],...
    clr{1},'EdgeColor','none');
plot(1:nEnsemble,m{3},'Color',clr{3},'LineWidth',1);
plot(2:2:50,m{1}(1:25),'Color',clr{1},'LineWidth',1);
set(gca,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35,...
    'XTick',[1 10 20 30 40 50],'YTick',0.5:0.1:0.7,'XLim',[1 nEnsemble],...
    'YLim',[0.45 0.7])
alpha(0.2);
xlabel('Number of neurons','FontSize',5);
ylabel('P_c_o_r_r_e_c_t','FontSize',5);


