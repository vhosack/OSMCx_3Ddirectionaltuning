%% Get subsamples of FA
%% 3D 4/8 directions
sample = [];
sampleidx = [];
n = 80; %% number of trials per direction
Con = {};
ConI = {};
for i = 1:10
    for d = 1:length(catind)
        dir = mag{d};
        diridx = catind{d};
            m = randperm(length(catind{d}),n);
            sample(d,:) = dir(m(1:n)); 
            sampleidx(d,:) = diridx(m(1:n));
    end
    
    Con{i} = sample;
    ConI{i} = sampleidx(:,:);
end

%% MIo
% Control
M1F = load('20190228_M1F_sortedspikes.mat'); %% paths to neural data
M1U = load('20190228_M1U_sortedspikes.mat');
M1F = struct2cell(M1F); 
M1U = struct2cell(M1U);
Units = {[M1F; M1U]};
Units = Units{1};
% Units = M1F; %% run individual area

spikecount = [];
timespikes = {};
meanfiring = [];
Counts = {};
for u = 1:length(Units)
    temp = Units{u}.times;
    for t=1:length(ControlTimes)

        % 1ms bins for FA %
        edges = ControlTimes : 0.001 : (Cendtimes + 0.001); %
        allCounts = zeros(1, length(edges) - 1); %
        [counts, values] = histcounts(temp, edges); %
        allCounts = allCounts + counts; %
        Counts{t}(u,:) = allCounts;

        inx=[];
        inx=find(temp >= ControlTimes(t) & temp <= Cendtimes(t));
        if inx == 0
            spikecount(t,u) = [];
        end
        spikecount(t,u) = length(inx);
        Ts = []; Ts = temp(inx);
        if ~isempty(Ts)
            Ts = Ts - temp(1);
        end
        timespikes{t} = Ts;
        meanfiring(t,u) = spikecount(t,u)/0.100; %% 100 ms
    end
end

% % Ranges of degrees
% % directions = 0:5:30;
% directions = 1:9;
% firingrate = [];
% meanfr = [];
% FRbyUnit = {};
% M1Con = {};
% for i = 1:10
%     for d = 1:(length(directions)-1)
%         dir = ConI{i}(d,:);
%         ncycles = length(dir);
%         cycles = meanfiring(dir,:);
%         firingrate = cycles;
%         for n = 1:length(Units)
%             meanfr(d,n) = sum(firingrate(:,n))/ncycles;
%             FRbyUnit{i}{n}(d,:) = firingrate(:,n);
%         end
%     end
%     M1Con{i} = meanfr;
% end

%% Create structure (Feeding)
directions = {'-++', '+++', '--+', '+-+','-+-','++-','---','+--'};
CountsD = {};
conD = {};
DSFALL = {};
newD = {};
Lat = {};
for i = 1:10
    for d = 1:length(directions)
            CountsD{d} = Counts(catind{d}); % Control
            conD{d} = cat(catind{d});
%             CountsD{i}{d} = Counts(ConI{i}(d,:)); % sample
%             conD{i}{d} = cat(ConI{i}(d,:));
    
    %         CountsD{d} = NBCounts(catind2{d}); % NB
    %         conD{d} = cat2(catind2{d});
    %         CountsD{d} = NBCounts(NBI(d,:)); % sample
    %         conD{d} = cat2(NBI(d,:));
    end
    data = [CountsD{i}{1}'; CountsD{i}{2}'; CountsD{i}{3}'; CountsD{i}{4}'; CountsD{i}{5}'; CountsD{i}{6}'; CountsD{i}{7}'; CountsD{i}{8}'];
    condition = [conD{i}{1}; conD{i}{2}; conD{i}{3}; conD{i}{4}; conD{i}{5}; conD{i}{6}; conD{i}{7}; conD{i}{8}];
    DSFALL{i} = struct('data',data,'type','traj','condition',condition,'epochStarts',1,'epochColors',[0.598590500684565,0.668487146264776,0.894564090918275]);
end

%%
DSF = struct('data',data,'type','traj','condition',condition,'epochStarts',1,'epochColors',[0.598590500684565,0.668487146264776,0.894564090918275]);
%% Run FA (remove low firing rate neurons)
alg = 3;
handles=[];
newD = {};
Lat = {};
% for i = 1:10
%     DSF = DSFALL{1};

    % Remove low firing rate neurons
    avgFR =[];
    keep = [];
    for t = 1:length(DSF)
        avgFR(:,t) = mean(DSF(t).data,2);
    end
    avgAll = mean(avgFR,2);
    keep = avgAll > 0.001;
    
    for t = 1:length(DSF)
        DSF(t).data = DSF(t).data(keep,:);
    end
    
    %%
%     % Get binned spike counts (10ms)
%     binWidth = 10; % ms
%     D = struct('data',[],'type','traj','condition',condition,'epochStarts',1,'epochColors',[0.598590500684565,0.668487146264776,0.894564090918275]);
%     for t = 1:length(DSF)
%         ndata = DSF(t).data;
%         for b = 1:(width(ndata)-1)/binWidth
%             D(t).data(:,b) = sum(ndata(:,1:binWidth*b+1),2);
%         end
%     end

%% Run FA
% candidateDims = height(DSF(1));
% [projs, mse, like] = cvreducedims_edit(DSF, alg, candidateDims, handles);
% [~, idx] = max(like); % find q that maximises likelihood of data
% dims = candidateDims(idx);
    % Dimensionality reduction
    dims = 20;
    [newD, C, Lat, explained, params] = reducedims_edit(DSF,alg, dims, handles);
% end

%% Compare cumulative variance

%% Compare sample dists
Z = cat(3,s1_dist,s2_dist,s3_dist,s4_dist,s5_dist,s6_dist,s7_dist,s8_dist,s9_dist,s10_dist);

sampMeans = mean(Z,3);

%%
meanfull = mean(Full_dist,2);
meansamps = mean(sampMeans,2);

%%
means = [meanfull meansamps];
[p,tbl,stats] = anova1(means);
%%
[c,~,~,gnames] = multcompare(stats);

%% Each bin
p = [];
for t = 1:width(Full_dist)
    means = [Full_dist(:,t) sampMeans(:,t)];
    [p(t),tbl,stats] = anova1(means,[],'off');
end

%% Create structure
DSF = struct('data',data,'type','traj','condition',condition,'epochStarts',1,'epochColors',[0.598590500684565,0.668487146264776,0.894564090918275]);

%% Sample 24 random neurons
% D = DSF;
indices = randperm(height(DSF(1).data),100);
for t = 1:length(DSF)
    DSF(t).data = DSF(t).data(indices,:);
end
% run FA above

%% Plot cumulative variance
figure;
plot(Lat);
hold on;
xlabel('Latent Factor');
ylabel('Cumulative variance');

%%
DataHigh(DSF,'DimReduce');