%% KNN Classifier to decode feeding direction from firing rate of each neuron
%% Create data matrix
D1 = zeros(1,90); D1(:) = 1;
D2 = zeros(1,90); D2(:) = 2;
D3 = zeros(1,90); D3(:) = 3;
% D4 = zeros(1,90); D4(:) = 4;
% D5 = zeros(1,90); D5(:) = 5;
% D6 = zeros(1,90); D6(:) = 6;
% D7 = zeros(1,90); D7(:) = 7;
% D8 = zeros(1,90); D8(:) = 8;
vecI = [D1 D2 D3];% D4];% D5 D6 D7 D8];
vecI = vecI';

% Control
Vec = [];
for u = 1:length(FRbyUnit)
    Vec(u,:) = [FRbyUnit{u}(1,:) FRbyUnit{u}(2,:) FRbyUnit{u}(3,:) FRbyUnit{u}(4,:) FRbyUnit{u}(5,:) FRbyUnit{u}(6,:)];% FRbyUnit{u}(7,:) FRbyUnit{u}(8,:)];
end
FR = Vec';

% % Nerve Block
% Vec = [];
% for u = 1:length(FRbyUnit2)
%     Vec(u,:) = [FRbyUnit2{u}(1,:) FRbyUnit2{u}(2,:) FRbyUnit2{u}(3,:) FRbyUnit2{u}(4,:) FRbyUnit2{u}(5,:) FRbyUnit2{u}(6,:)];
% end
% FR = Vec';
%% Select training and test trials (80:20 split)
ind = {};
allind = {};
notind = {};
for i = 1:100
  ind{i} = randperm(height(FR), floorDiv(height(FR),5)); % testing trial indices
  allind{i} = 1:height(FR);
  notind{i} = setdiff(allind{i},ind{i}); % training trial indices
end
% Note: keep same indices for matching MIo and SIo data
%% Run KNN for 28 neurons
percentCorrect = [];
for i = 1:100
%   neuronind = randperm(width(FR),28); %% select 28 random neurons
%   FRnew = FR(:,neuronind);
FRnew = FR; %% all neurons

  fr2 = FRnew(ind{i},:);
  response2 = vecI(ind{i});
  fr1 = FRnew(notind{i},:);
  response1 = vecI(notind{i});

  Mdl = fitcknn(fr1,response1,'NumNeighbors',7,'NSMethod','exhaustive','Distance','euclidean','Standardize',1);

  class = predict(Mdl, fr2);
  correct = response2 == class;
  percentCorrect(i) = sum(correct)/length(correct);
  clear('Mdl');
end

Results = sum(percentCorrect)/length(percentCorrect)
stdev = std(percentCorrect)