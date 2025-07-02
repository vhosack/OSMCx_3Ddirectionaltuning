%% Load data
load('20190228_Kinematics.mat')
%% Get tongue tip position
ttcols = find(contains(Kinematics.ColumnNames.points,'AnteriorM_'));
tonguetip = Kinematics.Cranium.points{2}(:,ttcols);
temptimes = Kinematics.NeuralIndex{2}(:,3);
trial = Kinematics.GapeCycleInfo(find(contains(Kinematics.GapeCycleInfo.CycleType,'hewL')),:);%contains(Kinematics.GapeCycleInfo.Trialname,'005')),:);% & contains(Kinematics.GapeCycleInfo.CycleType,'hewR')),:);
%% Plot tongue tip during chews
traj = {};
times = {};
for i = 1:height(trial)
    tonguetip = Kinematics.Cranium.points{strcmp(trial.Trialname(i),Kinematics.TrialNames)}(:,ttcols);
    traj{i} = tonguetip(Kinematics.GapeCycleInfo.MaxGapeStart(i):Kinematics.GapeCycleInfo.MaxGapeEnd(i),:);
    temptimes = Kinematics.NeuralIndex{strcmp(trial.Trialname(i),Kinematics.TrialNames)}(:,3);
    times{i} = temptimes(Kinematics.GapeCycleInfo.MaxGapeStart(i):Kinematics.GapeCycleInfo.MaxGapeEnd(i),:);
end
%%
% M1F = load('20190228_M1F_sortedspikes.mat'); %% paths to neural data
M1U = load('20190228_M1U_sortedspikes.mat');
% M1F = struct2cell(M1F); 
M1U = struct2cell(M1U);
% Units = [M1F; M1U];

Units = M1U(Ind_Fd);
%%
Counts = {};
for u = 1:length(Units)
    temp = Units{u}.times;
    for t=1:length(traj)

        if isnan(times{t}(end))
            continue
        elseif isnan(times{t}(1))
            continue
        end

        % 1ms bins for FA %
        edges = times{t}(1) : 0.001 : (times{t}(end) + 0.001); %
        allCounts = zeros(1, length(edges) - 1); %
        [counts, values] = histcounts(temp, edges); %
        allCounts = allCounts + counts; %
        CountsFd{t}(u,:) = allCounts;
    end
end
%%
data = Counts;
Data = struct();
for i = 1:length(data)
    Data(i).data = data{i};
%     Data(i).condition = vecI{i};
%     {data{i}; vecI{i}};
end

%% Dimensionality reduction
% D(itrial).data : (num_neurons x num_1ms_bins)
DataHigh(Data, 'DimReduce');
% DataHigh(D);

%%
figure;
for d = 1:length(D)
%     plot(D(d).data(1,:),D(d).data(2,:),'-');
    plot3(D(d).data(1,:),D(d).data(2,:),D(d).data(3,:),'-','LineWidth',2);%,'Color','k', 'LineWidth', 1.5);
    hold on;
%     plot3(D(d).data(1,:),-D(d).data(3,:),D(d).data(2,:),'-');
    scatter3(D(d).data(1,end),D(d).data(2,end),D(d).data(3,end),75,'>','filled');%,'MarkerFaceColor','k');
    hold on;
end
% title('Monkey R MIo Con');
xlabel('Factor 1');
ylabel('Factor 2');
zlabel('Factor 3');

%% 500 ms around min gape
%% Load data
load('20190228_Kinematics.mat')
%% Get tongue tip position
ttcols = find(contains(Kinematics.ColumnNames.points,'AnteriorM_'));
tonguetip = Kinematics.Cranium.points{2}(:,ttcols);
temptimes = Kinematics.NeuralIndex{2}(:,3);
% trial = Kinematics.GapeCycleInfo(find(contains(Kinematics.GapeCycleInfo.Trialname,'005') & contains(Kinematics.GapeCycleInfo.CycleType,'hewR')),:);
trial = Kinematics.GapeCycleInfo;
%% Plot tongue tip during chews (200 ms)
traj = {};
times = {};
for i = 1:height(trial)
    tonguetip = Kinematics.Cranium.points{strcmp(trial.Trialname(i),Kinematics.TrialNames)}(:,ttcols);
    traj{i} = tonguetip(Kinematics.GapeCycleInfo.MinGape(i)-20:Kinematics.GapeCycleInfo.MinGape(i)+20,:);
    temptimes = Kinematics.NeuralIndex{strcmp(trial.Trialname(i),Kinematics.TrialNames)}(:,3);
    times{i} = temptimes(Kinematics.GapeCycleInfo.MinGape(i)-20:Kinematics.GapeCycleInfo.MinGape(i)+20,:);
end

%% Drinking
%% Load data
load('20181210_Kinematics_Before.mat')
%% Get tongue tip position
% Min tongue protrusion to spout
ntrials = length(Kinematics.Cranium.points);
Neural = Kinematics.NeuralIndex;

ttptcols = find(contains(Kinematics.ColumnNames.points,'AnteriorM_'));

ttpt = {};
neural = {};
cycles = {}; Cycles = {};
ttAP = {};
Mins = {}; mins = {};
nAP = {}; traj = {}; times = {};
for i = 1:ntrials
    ttpt{i} = Kinematics.Cranium.points{i}(:,ttptcols);
    neural{i} = Neural{i}(:,3);
    npeaks = length(ttpt{i});
    stretch = ttpt{i}(:,1);
    [~,minstretch] = findpeaks(stretch*-1,'MinPeakProminence',10);
    [~,maxstretch] = findpeaks(stretch,'MinPeakProminence',10);
    if length(minstretch) > length(maxstretch)
        minstretch = minstretch(1:end-1);
    end
%     Mins = [ttpt{i}(:,1), ttpt{i}(:,2), ttpt{i}(:,3)]; mins = Mins(minstretch,:);
%     Maxs = [ttpt{i}(:,1), ttpt{i}(:,2), ttpt{i}(:,3)]; maxs = Maxs(maxstretch,:);
    Mins{i} = ttpt{i}(minstretch,1);
    mins{i} = minstretch;
    for l = 1:length(minstretch)
        ttAP{i}{l} = ttpt{i}(minstretch(l)-20:minstretch(l)+20,1:3);
        nAP{i}{l} = neural{i}(minstretch(l)-20:minstretch(l)+20);
    end
end
traj = horzcat(ttAP{:});
times = horzcat(nAP{:});

traj = traj(~cellfun('isempty',traj));
times = times(~cellfun('isempty',times));

%% Length of Trajs
Dists = {};
totalDist = [];
for t = 1:length(traj)

    for i = 1:height(traj{t})-1
        t1 = traj{t}(i,:);
        t2 = traj{t}(i+1,:);
        Dists{t}(i) = norm(t2-t1);
    end
    totalDist(t) = sum(Dists{t});
end

mean(totalDist,'omitnan')
std(totalDist,[],'omitnan')

%% Stable Units
% load RyDataStats2
inx = scoreW(:,2)>0.99;
% sum(inx);
stableUnits = units(inx,:);
stableM1U = stableUnits(1:34,:); % 1-96

%%
% M1F = load('20181210b_M1F_sortedspikes.mat'); %% paths to neural data
M1U = load('20181210b_M1U_sortedspikes.mat');
% M1F = struct2cell(M1F); 
M1U = struct2cell(M1U);
% Units = [M1F; M1U];

Units = M1U(Ind_Dr);
%%
CountsDr = {};
for u = 1:length(Units)
    temp = Units{u}.times;
    for t=1:length(traj)

        if isnan(times{t}(end))
            continue
        elseif isnan(times{t}(1))
            continue
        end

        % 1ms bins for FA %
        edges = times{t}(1) : 0.001 : (times{t}(end) + 0.001); %
        allCounts = zeros(1, length(edges) - 1); %
        [counts, values] = histcounts(temp, edges); %
        allCounts = allCounts + counts; %
        CountsDr{t}(u,:) = allCounts;
    end
end


%% Get Feeding units (utah)
% M1U = load('20190228_M1U_sortedspikes.mat');
% M1U = struct2cell(M1U);
% Units = M1U; %% run individual area
%%
% Counts_Dr = CountsDr;
% Counts_Fd = CountsFd;

% for s = 1:length(stableUnits)
%     stableInd = find(contains('elec' num2str(stableM1U(s,1)) '_' num2str(stableM1U(s,2), )));
% end

Counts = Counts_Fd;
count = {};
avgfiringrate = {};
for t = 1:length(Counts)
    count{t} = sum(Counts{t},2);
    firingrate{t} = count{t}/0.2;
end

meanfiringrate = mean(cat(3,firingrate{:}),3);

%% Combine drinking and feeding
Counts_Fd = Counts_Fd(~cellfun('isempty',Counts_Fd)) ; 
Counts_Dr = Counts_Dr(~cellfun('isempty',Counts_Dr)) ; 

F = cell(1,length(Counts_Fd)); F(:) = {'Feed'};
D = cell(1,length(Counts_Dr)); D(:) = {'Drink'};
vecI = [F D];
condition = vecI'; %% directions

data = [Counts_Fd Counts_Dr];
Data = struct();
for i = 1:length(data)
    Data(i).data = data{i};
    Data(i).condition = condition{i};
%     {data{i}; vecI{i}};
    Data(i).type = 'traj';
    Data(i).epochStarts = 1;
    Data(i).epochColors = [0.598590500684565,0.668487146264776,0.894564090918275];
end

%% Dimensionality reduction
% D(itrial).data : (num_neurons x num_1ms_bins)
DataHigh(Data, 'DimReduce');
% DataHigh(D);

%% Subsample trials

indices_Fd = randperm(length(Counts_Fd),10); C_Fd = Counts_Fd(indices_Fd);
indices_Dr = randperm(length(Counts_Dr),10); C_Dr = Counts_Dr(indices_Dr);

F = cell(1,length(C_Fd)); F(:) = {'Feed'};
D = cell(1,length(C_Dr)); D(:) = {'Drink'};
vecI = [F D];
condition = vecI'; %% directions

data = [C_Fd C_Dr];
Data = struct();
for i = 1:length(data)
    Data(i).data = data{i};
    Data(i).condition = condition{i};
%     {data{i}; vecI{i}};
    Data(i).type = 'traj';
    Data(i).epochStarts = 1;
%     Data(i).epochColors = [0.598590500684565,0.668487146264776,0.894564090918275];
end

%% average trajectories
% Feeding
fdata = cat(3,D(1:1118).data);
avgdata = mean(fdata,3);
figure;
plot3(avgdata(1,:), avgdata(2,:), avgdata(3,:),'-','LineWidth',2);
hold on;
scatter3(avgdata(1,end), avgdata(2,end), avgdata(3,end), 75,'>','filled');

% Drinking
ddata = cat(3,D(1119:1525).data);
avgdataD = mean(ddata,3);
plot3(avgdataD(1,:), avgdataD(2,:), avgdataD(3,:),'-','LineWidth',2);
scatter3(avgdataD(1,end), avgdataD(2,end), avgdataD(3,end), 75,'>','filled');

%% Separate feed and drink
DSF = Data(1:1118);
% DSF = Data(1119:1525);
%%
DataHigh(DSF, 'DimReduce');
%%
% figure;
hold on;
plot3(D.data(1,:),D.data(2,:),D.data(3,:),'-','LineWidth',2);
% hold on; 
scatter3(D.data(1,end),D.data(2,end),D.data(3,end),75,'>','filled');
xlabel('Factor 1');
ylabel('Factor 2');
zlabel('Factor 3');

%% fd and dr inter-traj dist
Dists = [];
for tp = 1:width(Fd.data)
    pairfirst = Fd.data(:,tp);
    pairsecond = Dr.data(:,tp);
    Dists(:,tp) = norm(pairfirst-pairsecond);
end
mean(Dists)
std(Dists)

[h,p] = ttest(Dists);

%% compare fd and dr distance travelled
Vec = [interDistsF1_Fd interDistsF1_Dr];

C = cell(1,length(interDistsF1_Fd)); C(:) = {'Feeding'};
B = cell(1,length(interDistsF1_Dr)); B(:) = {'Drinking'};
vecI = [C B];
vecI = vecI';

%% Factor 1 inter-traj distance
for p = 1:length(D)
    for di = 1:length(D)
        if p == di
            continue
        end
        for tp = 1:width(D(p).data)
            pairfirst = D(p).data(1,tp);
            pairsecond = D(di).data(1,tp);
            Dists{p,di}(tp) = norm(pairfirst-pairsecond);
        end
    end
end

dists = Dists(~cellfun(@isempty,Dists));
Distances = cell2mat(dists);
Distances = unique(Distances,'rows','stable');
Avg = mean(Distances,1);
AllAvgs = mean(Avg)
dev = std(Avg)

%% one t-test (diff from 0)
p = [];
for pair = 1:height(Distances)
    [h,p(pair)] = ttest(Distances(pair,:));
end
