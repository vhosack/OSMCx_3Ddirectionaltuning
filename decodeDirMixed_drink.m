%% Combined M1 and S1 KNN Decoding (Drinking)
%% Create data matrix
L = cell(1,height(meanfiring{1}));
L(:) = {'Left'};
M = cell(1,height(meanfiring{2}));
M(:) = {'Middle'};
R = cell(1,height(meanfiring{3}));
R(:) = {'Right'};
vecI = [L M R];
vecI = vecI'; %% corresponding directions

FR = [];
for u = 1:width(meanfiring{1})
    FR(:,u) = [meanfiring{1}(:,u); meanfiring{2}(:,u); meanfiring{3}(:,u)]; %% trial x neuron
end
%% Make M1 matrix
FR_M1 = FR;
%% Make S1 matrix
FR_S1 = FR;
%% Select traning and test trials (80:20 split)
ind = {};
allind = {};
notind = {};
for i = 1:100
  ind{i} = randperm(height(FR), floorDiv(height(FR),5)); % testing trial indices
  allind{i} = 1:length(FR);
  notind{i} = setdiff(allind{i},ind{i}); % training trial indices
end
% Note: keep same indices for matching MIo and SIo data
%% Run KNN for mixed population
percentCorrect = [];
for i = 1:100
    M1ind = randperm(width(FR_M1),25); %% select 25 random M1 neurons
    FR = FR_M1;
    FR(:,M1ind) = [];

    S1ind = randperm(width(FR_S1),25); %% select 25 random S1 neurons
    FRS = FR_S1(:,S1ind);

    FRnew = [FR FRS];

  fr2 = FRnew(ind{i},:);
  response2 = vecI(ind{i});
  fr1 = FRnew(notind{i},:);
  response1 = vecI(notind{i});

  Mdl = fitcknn(fr1,response1,'NumNeighbors',7,'NSMethod','exhaustive','Distance','euclidean','Standardize',1);

  class = predict(Mdl, fr2);
  correct = strcmp(response2, class);
  percentCorrect(i) = sum(correct)/length(correct);
  clear('Mdl');
end

Results = sum(percentCorrect)/length(percentCorrect)
stdev = std(percentCorrect)