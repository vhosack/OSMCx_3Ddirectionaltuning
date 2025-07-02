%% Drinking directions
L = cell(1,height(meanfiring{1}));
L(:) = {'Left'};
M = cell(1,height(meanfiring{2}));
M(:) = {'Middle'};
R = cell(1,height(meanfiring{3}));
R(:) = {'Right'};
vecI = [L M R];
vecI = vecI'; %% directions

%% Run bootstrap
% run after getting mean firing rates from directionDrink.m
nBoot = 1000;
tetas={'Left'; 'Middle'; 'Right'};

Vec = {}; p = []; m = {}; mDir = {};
% count=1;
% tCount=0;
M = []; bM = []; pBoot = []; pB = [];
Chi = []; ChiBoot = [];
for u = 1:width(meanfiring{1})
    Vec{u} = [meanfiring{1}(:,u); meanfiring{2}(:,u); meanfiring{3}(:,u)]; %% trial x neuron
    n=length(Vec{u});
    randP=round(rand(n*nBoot,1)*n+0.5);
    indx=reshape(randP,[],nBoot);
    bootRates=Vec{u}(indx);
%     index=vecI(indx);
%     tmpRates=Vec{u};
%     M(tCount)=mean(tmpRates);
%     bM(tCount,:)=mean(bootRates);
%     count=count+n;   
    [p,tbl,stats] = kruskalwallis(Vec{u},vecI','off');
    Chi(u) = cell2mat(tbl(2,5));
    for i = 1:nBoot
        y = bootRates(:,i);
%         group = index(:,i);
%         group_shuffled = group(randperm(numel(group)));
        [p,tbl,stats] = kruskalwallis(y',vecI','off');
        ChiBoot(u,i) = cell2mat(tbl(2,5));
    end
    pB(u) = sum(ChiBoot(u,:) > Chi(u))/nBoot;
end

tuned = pB < 0.05;
percentTuned = mean(tuned)*100

%% Feeding data structure
% C = horzcat(ConI{:});
C = horzcat(ConI(:));
% C = horzcat(NBI(:));
directions = {'-++', '+++', '--+', '+-+','-+-','++-','---','+--'};
% directions = {'-++', '+++', '+-+','---','+--'}; %% swallows
MeanD = {};
conD = {};
for d = 1:length(directions)
        MeanD{d} = meanfiring(catind{d},:); % Control
        conD{d} = cat(catind{d});
%         MeanD{d} = meanfiring(C(d,:),:); % sample
%         conD{d} = cat(C(d,:));

%         MeanD{d} = meanfiring(catind2{d},:); % NB
%         conD{d} = cat2(catind2{d});
%         MeanD{d} = meanfiring(NBI(d,:),:); % sample
%         conD{d} = cat2(NBI(d,:));
end
condition = [conD{1}; conD{2}; conD{3}; conD{4}; conD{5}; conD{6}; conD{7}; conD{8}];

%% Run bootstrap
% run after getting mean firing rates from directionFeed.m
nBoot = 1000;

Vec = {}; p = []; m = {}; mDir = {};
M = []; bM = []; pBoot = []; pB = [];
Chi = []; ChiBoot = [];
for u = 1:width(MeanD{1})
%     Vec{u} = [MeanD{1}(:,u); MeanD{2}(:,u); MeanD{3}(:,u); MeanD{4}(:,u); MeanD{5}(:,u)];%; MeanD{6}(:,u); MeanD{7}(:,u); MeanD{8}(:,u)]; %% trial x neuron
    for d = 1:length(MeanD)
        n=80; %60; length(Vec{u});
        randP=round(rand(n*nBoot,1)*height(MeanD{d})+0.5);  % n -> height(MeanD{d})
        indx=reshape(randP,[],nBoot);
        bRates{d}=MeanD{d}(indx); % Vec{u}(indx) -> MeanD{d}(indx)
    end
    bootRates = vertcat(bRates{:});
    
    [p,tbl,stats] = kruskalwallis(meanfiring(C,u),cat(C)','off'); % Vec{u}(C) -> meanfiring(C,u); condition(C) -> cat(C)
    Chi(u) = cell2mat(tbl(2,5));
    for i = 1:nBoot
        y = bootRates(:,i);
%         group = indx(:,i);
%         group_shuffled = group(randperm(numel(group)));
        [p,tbl,stats] = kruskalwallis(y',cat(C)','off'); % condition(C) -> cat(C)
        ChiBoot(u,i) = cell2mat(tbl(2,5));
    end
    pB(u) = sum(ChiBoot(u,:) > Chi(u))/nBoot;
end

tuned = pB < 0.05;
percentTuned = mean(tuned)*100

%% Preferred Directions
%% Drinking
%% Get meanfr of only tuned 
Left = meanfiring{1}(:,tuned);
Middle = meanfiring{2}(:,tuned);
Right = meanfiring{3}(:,tuned);

% Left
meanfr = [];
ncycles = height(Left);
for n = 1:width(Left)
        meanfr(n) = sum(Left(:,n))/ncycles;
end
FR_l = meanfr;

% Middle
ncycles = height(Middle);
for n = 1:width(Middle)
        meanfr(n) = sum(Middle(:,n))/ncycles;
end
FR_m = meanfr;

% Right
ncycles = height(Right);
for n = 1:width(Right)
        meanfr(n) = sum(Right(:,n))/ncycles;
end
FR_r = meanfr;

MeanFR = [FR_l; FR_m; FR_r]';
%% Find PD
pref = zeros(height(MeanFR),1);
for unit = 1:height(MeanFR)
    [p, d] = max(MeanFR(unit,:));
    pref(unit) = d;
end

% figure;
% histogram(pref,'BinLimits',[0.5,3.5],'normalization','pdf');
% xlabel('Preferred Direction');
% ylabel('Proportion of Neurons');
% title('Histogram of Preferred Directions of Sensory Neurons (Nerve Block)');
%% Save to compare Con with NB
pref1 = pref;





%%
L = cell(1,height(meanfiring{1}));
L(:) = {'Left'};
M = cell(1,height(meanfiring{2}));
M(:) = {'Middle'};
R = cell(1,height(meanfiring{3}));
R(:) = {'Right'};
vecI = [L M R];
vecI = vecI'; %% directions

Vec = {};
p = [];
for u = 1:width(meanfiring{1})
    Vec{u} = [meanfiring{1}(:,u); meanfiring{2}(:,u); meanfiring{3}(:,u)]; %% trial x neuron
end

%% Run Drinking
nBoot = 10000;

bootRates = {};
M = []; bM = []; pB = [];
for d = 1:3
    for u = 1:width(meanfiring{d})
        tmpRates=meanfiring{d}(:,u);
        M{d}(u)=mean(tmpRates);
        bootRates = meanfiring{d}(:,u);
        bM{d}(:,u) = bootstrp(nBoot,@mean,tmpRates);
%         st{d}(:,u) = bootstrp(1000,@std,temp);
    end
end

A = abs(M{1}-M{2}); B = abs(M{2}-M{3}); C = abs(M{1}-M{3});
bA = abs(bM{1}-bM{2}); bB = abs(bM{2}-bM{3}); bC = abs(bM{1}-bM{3});
diff = (A+B+C)/3;
diffBoot= (bA+bB+bC)/3;

% sorted = sort(diff);
sortedBoot = sort(diffBoot);

for u = 1:width(meanfiring{d})
    pB(u) = length(find(sortedBoot(:,u)>diff(u)))/nBoot;
end

tuned = pB < 0.05;
percentTuned = mean(tuned)*100

%% Run Feeding
nBoot = 1000;

bootRates = {};
M = []; bM = []; pB = [];
dv = {}; dvBoot = {}; meandiff = {}; meandiffBoot = {};
for u = 1:length(FRbyUnit)
    for d = 1:length(FRbyUnit{u})
        tmpRates=FRbyUnit{u}{d};
        M{u}(d)=mean(tmpRates);
        bootRates = FRbyUnit{u}{d};
        bM{u}(:,d) = bootstrp(nBoot,@mean,tmpRates);
%         st{d}(:,u) = bootstrp(1000,@std,temp);
%         dvBoot{u} = abs(bM{u}(:,d)-bM{u}(:,d)');
    end
    dv{u} = abs(M{u}-M{u}');
    for i = 1:nBoot
        dvBoot{u} = abs(bM{u}(i,:)-bM{u}(i,:)');
    end

%     A = M{1}-M{2}; B = M{2}-M{3}; C = M{1}-M{3};
%     bA = bM{1}-bM{2}; bB = bM{2}-bM{3}; bC = bM{1}-bM{3};
%     diff{u} = (A+B+C)/3;
%     diffBoot= (bA+bB+bC)/3;
    meandiff{u} = mean(nonzeros(dv{u}));
    meandiffBoot{u} = mean(nonzeros(dvBoot{u}));
end

% sorted = sort(diff);
sortedBoot = sort(diffBoot);

for u = 1:length(FRbyUnit)
    pB(u) = length(find(sortedBoot(:,u)<diff(u)))/nBoot;
end

tuned = pB < 0.05;
percentTuned = mean(tuned)*100

%%
L = cell(1,height(m{1}));
L(:) = {'Left'};
M = cell(1,height(m{2}));
M(:) = {'Middle'};
R = cell(1,height(m{3}));
R(:) = {'Right'};
vecI = [L M R];
vecI = vecI'; %% directions
vecI_shuffled = vecI(randperm(numel(vecI)));

%%
Vec = {};
p = [];
p_shuffled = [];
for u = 1:width(meanfiring{1})
    Vec{u} = [m{1}(:,u); m{2}(:,u); m{3}(:,u)]; %% trial x neuron
    p(u) = anova1(Vec{u},vecI','off');
%     p(u) = kruskalwallis(Vec{u},vecI','off');
%     p_shuffled(u) = kruskalwallis(Vec{u},vecI_shuffled','off');
end
ps = mean(p_shuffled);

tuned = p < 0.05;
percentTuned = mean(tuned)*100

% tuned_shuffled = p_shuffled < 0.05;
% percentTuned_shuffled = mean(tuned_shuffled)*100
% 
% p_boot = length(find(p<p_shuffled))/width(meanfiring{d})

%% Find PDs
% Assign variables
tuned = [tunedRyDrM1FCon tunedRyDrM1UCon];
spikeCountCellArray = [sVRyM1FCon sVRyM1UCon];

%% Get meanfr of only tuned
onlytuned = cell(length(spikeCountCellArray));
onlytuned(tuned) = spikeCountCellArray(tuned);
Tuned = onlytuned(~cellfun(@isempty, onlytuned));
meanfr = [];
for u = 1:length(Tuned)
    Data = Tuned{u};
    for d = 1:length(spikeCountCellArray{u})
        temp = Data{d}/0.500;  % 0.100 feeding
        meanfr(d,u) = mean(temp);
    end
end
FR = meanfr';

%% Find PD Bootstrap
% Control
nBoot = 1000;

Vec = {};
pref = zeros(length(spikeCountCellArray),1);
for u = 1:length(spikeCountCellArray)
    for d = 1: length(spikeCountCellArray{u})
        n=length(spikeCountCellArray{u}{d});
        randP=round(rand(n*nBoot,1)*n+0.5);
        indx=reshape(randP,[],nBoot);
        bootRates=spikeCountCellArray{u}{d}(indx);
        meanR{u}(d) = mean(spikeCountCellArray{u}{d});
        for i = 1:nBoot
            y = bootRates(:,i);
    %         group = index(:,i);
    %         group_shuffled = group(randperm(numel(group)));
            meanBoot{u}(i,d) = mean(y);
        end
    end
    [p, d] = max(meanR{u});
    pref(u) = d;
    [pB, dB] = max(meanBoot{u},[],2);
    prefB(u) = mode(dB);
    pBoot(u) = sum(dB==d)/nBoot;
end
mean(pBoot)
pref1 = pref;

figure;
histogram(pref,'BinLimits',[0.5,3.5],'normalization','probability','FaceAlpha',0.8,'FaceColor',[0.9290 0.6940 0.1250],...
    'Normalization','probability','BinMethod','auto');
xlabel('Preferred Direction');
ylabel('Proportion of Neurons');
title('MIo');
set(gca,'FontSize',14,'YGrid','on');
