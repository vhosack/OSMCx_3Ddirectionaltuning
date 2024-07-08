%% Get proportions of neurons that gained/lost tuning and compare across groups
%% paste list of all neurons
neurons = {}; % Control
neurons1 = {}; % Nerve Block
%% index of only neurons present in both lists
membind = ismember(neurons,neurons1);
membind1 = ismember(neurons1,neurons);
%% run and get 2 tuninginds from pref_ (tuned)
% 1 = tuned; 0 = not tuned
%% 
% neurons and whether they are tuned
tuningtbl = table(neurons,tuningind);
tuningtbl1 = table(neurons1,tuningind1);

% only neurons in both control and NB
newtbl = tuningtbl(membind,:);
newtbl1 = tuningtbl1(membind1,:);
%% remove those with no change (optional)
c = newtbl.tuningind & newtbl1.tuningind1;
tbl = newtbl(~c,:);
tbl1 = newtbl1(~c,:);
% 1 -> 0 = lost; 0 -> 1 = gained

%% Compare distributions
%% paste neurons which gained or lost tuning
glmat = {};
%% matrix of no change
nmat = cell(38,1);
nmat(:) = {'n'};
mat = [glmat; nmat];
S = mat;
%% categorize into 2 groups
Mind = cell(53,1); % total # of M1 neurons
Mind(:) = {'M1'};
Sind = cell(62,1); % total # of S1 neurons
Sind(:) = {'SC'};
%% combine
y1 = [M; S];
y2 = [Mind; Sind];
%% Chi-square test
[tbl,chi2stat,pval] = crosstab(y1,y2)