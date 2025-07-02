%% Determine if neuron is directionally tuned & get PDs of tuned (drinking)
%% Create data matrix
L = cell(1,height(meanfiring{1}));
L(:) = {'Left'};
M = cell(1,height(meanfiring{2}));
M(:) = {'Middle'};
R = cell(1,height(meanfiring{3}));
R(:) = {'Right'};
vecI = [L M R];
vecI = vecI'; %% directions

FR = [];
Vec = {};
for u = 1:width(meanfiring{1})
    FR(:,u) = [meanfiring{1}(:,u); meanfiring{2}(:,u); meanfiring{3}(:,u)]; %% trial x neuron
    Vec{u} = {meanfiring{1}(:,u)' meanfiring{2}(:,u)' meanfiring{3}(:,u)'};
end

%% Select traning and test trials (80:20 split)
ind = {};
allind = {};
notind = {};
for i = 1:100
  ind{i} = randperm(height(FR), floorDiv(height(FR),3)); % testing trial indices
  allind{i} = 1:length(FR);
  notind{i} = setdiff(allind{i},ind{i}); % training trial indices
end
% Note: keep same indices for matching MIo and SIo data

%% Get percent tuned
Vec = {};
p = [];
for i = 1:100
    for u = 1:width(meanfiring{1})
%         fr2 = FR(ind{i},:);
%         response2 = vecI(ind{i});
        fr1 = FR(notind{i},:);
        response1 = vecI(notind{i});
    
        p(u,i) = kruskalwallis(fr1,response1','off');
    end
end

tuned = p < 0.05;
percentTuned = mean(tuned)*100;
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