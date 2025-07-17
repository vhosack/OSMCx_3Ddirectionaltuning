%% Factor Analysis
%% Create data structure (Drinking)
L = cell(1,height(meanfiring{1})); L(:) = {'Left'};
M = cell(1,height(meanfiring{2})); M(:) = {'Middle'};
R = cell(1,height(meanfiring{3})); R(:) = {'Right'};
vecI = [L M R];
vecI = vecI'; % array of directions

data = [Counts{1}'; Counts{2}'; Counts{3}'];

Data = struct();
for i = 1:length(data)
    Data(i).data = data{i};
    Data(i).condition = vecI{i};
    Data(i).type = 'traj';
    Data(i).epochStarts = 1;
    Data(i).epochColors = [0.598590500684565,0.668487146264776,0.894564090918275];
end
% Data = Data';
% DS = struct('data',data,'condition',vecI);
DSF = Data;
%% Dimensionality reduction
% D(itrial).data : (num_neurons x num_1ms_bins)
DataHigh(Data, 'DimReduce');
%% Visualise only
DataHigh(D);
%% Create structure (Feeding)
directions = {'-++', '+++', '--+', '+-+','-+-','++-','---','+--'};
CountsD = {};
conD = {};
% uncomment current dataset
for d = 1:length(directions)
        CountsD{d} = Counts(catind{d}); % Control
        conD{d} = cat(catind{d});
%         CountsD{d} = Counts(ConI(d,:)); % sample
%         conD{d} = cat(ConI(d,:));

%         CountsD{d} = NBCounts(catind2{d}); % NB
%         conD{d} = cat2(catind2{d});
%         CountsD{d} = NBCounts(NBI(d,:)); % sample
%         conD{d} = cat2(NBI(d,:));
end
data = [CountsD{1}'; CountsD{2}'; CountsD{3}'; CountsD{4}'; CountsD{5}'; CountsD{6}'; CountsD{7}'; CountsD{8}'];
condition = [conD{1}; conD{2}; conD{3}; conD{4}; conD{5}; conD{6}; conD{7}; conD{8}];
DSF = struct('data',data,'type','traj','condition',condition);
%% Dimensionality reduction
DataHigh(DSF, 'DimReduce');

% Load .mat file in FA folder
%% Plot first 2/3 latent variables
% colorvec = {[0.39 0.83 0.07 0.6], [0.00 0.00 1.00 0.6], [0.93 0.69 0.13 0.7]};
figure;
for d = 1:length(D)
%     plot(D(d).data(1,:),D(d).data(2,:),'-');
    plot3(D(d).data(1,:),D(d).data(2,:),D(d).data(3,:),'-');%, 'Color',colorvec{d}, 'LineWidth',3);
%     plot(D(d).data(3,:));
    hold on;
%     plot3(D(d).data(1,:),-D(d).data(3,:),D(d).data(2,:),'-');
    scatter3(D(d).data(1,end),D(d).data(2,end),D(d).data(3,end),75,'>','filled');
%     hold on;
end
% title('Monkey R MIo Con');
xlabel('Factor 1');
ylabel('Factor 2');
zlabel('Factor 3');
%% Calculate distance b/t pairs of trajectories
Dists = {};
for p = 1:length(D)
    for di = 1:length(D)
        if p == di
            continue
        end
        for tp = 1:width(D(p).data)
            pairfirst = D(p).data(:,tp);
            pairsecond = D(di).data(:,tp);
            Dists{p,di}(:,tp) = norm(pairfirst-pairsecond);
        end
    end
end
dists = Dists(~cellfun(@isempty,Dists));
Distances = cell2mat(dists);
Distances = unique(Distances,'rows','stable');
Avg = mean(Distances,1);
AllAvgs = mean(Avg)
std(Avg)
%% anova on pairs
[p,tbl,stats] = anova1(Distances');
[c,~,~,gnames] = multcompare(stats);
%% one sided t-test (diff from 0)
p = [];
for pair = 1:height(Distances)
    [h,p(pair)] = ttest(Distances(pair,:));
end
%% Each pair quick set
ConDist = Distances;
%% NB
NBDist = Distances;
%% Compare each pair
p = [];
for n = 1:height(Distances)
    Con = ConDist(n,:);
    NB = NBDist(n,:);
    [h,p(n)] = ttest2(Con,NB);
end
sig = p < 0.05;
mean(sig)
%% Feeding: dist b/t Ant/Post pairs of directions
% Determine across which direction there is most/least variation compared
% to others

AntPost = [Dists{1,5}; Dists{2,6}; Dists{3,7}; Dists{4,8}]; % Ant-Post
% AntPost = [Dists{1,3}; Dists{2,4}; Dists{5,7}; Dists{6,8}]; % Sup-Inf
% AntPost = [Dists{1,2}; Dists{3,4}; Dists{5,6}; Dists{7,8}]; % Left-Right
avgAP = mean(AntPost,'all')
stdAP = std(AntPost,0,'all')

Others = setdiff(Distances,AntPost,'rows');
avgOthers = mean(Others,'all')
stdOthers = std(Others,0,'all')
%% Compare groups
[h,p] = ttest2(AntPost,Others);
mean(p)
%% Length of trajectories
Dists = {};
totalDist = [];
for d = 1:length(D)

    for t = 1:width(D(d).data)-1
        t1 = D(d).data(:,t);
        t2 = D(d).data(:,t+1);
        Dists{d}(t) = norm(t2-t1);
    end
    totalDist(d) = sum(Dists{d});
end
mean(totalDist)
std(totalDist)
%% Compare total dist travelled
%% Control
Contot = totalDist;
%% NB
NBtot = totalDist;
%% Compare each
[h,p] = ttest2(Contot,NBtot);
%% Plot inter-trajectory distances over time (MIo)
% t = -0.24:0.01:0.25;
t = 0.01:0.01:0.1;
mdists = mean(Distances, 1);
figure;
plot(t,mdists);
hold on;
xlabel('Time (s)');
ylabel('Normalized distance between trajectories');
%% Add SIo
sdists = mean(Distances, 1);
plot(t,sdists);
%% Normalize inter-trajectory distances
Distances=normalize(Distances,2,'range');
%% Normalize distance travelled
% M1dist = totalDist;
Mdist_Y = totalDist;
M1dist = [Mdist_R Mdist_Y];
%%
% S1dist = totalDist;
Sdist_Y = totalDist;
S1dist = [Sdist_R Sdist_Y];
%%
distNorm = M1dist-S1dist./(M1dist+S1dist);
[p,h] = signrank(distNorm)