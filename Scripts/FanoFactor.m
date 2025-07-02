%% Fano Factor
%% Drinking
FF = {};
ind = {};
% stab = {};
% s = {};
% figure;
for d = 1:3
    for u = 1:width(spikecount{dir})
        NumCounts_vect = spikecount{d}(:,u);
        FF{d}(u) = (std(NumCounts_vect)^2)/mean(NumCounts_vect);
%         plot(mean(NumCounts_vect),var(NumCounts_vect),'o');
%         hold on;
    end
%     FF{d} = rmmissing(FF{d});
    [R,ind{d}] = rmmissing(FF{d});
%     stab{d} = FF{d} < 1;
%     s{d} = mean(stab{d})*100;
end
Ind = ind{1} | ind{2} | ind{3};
for d = 1:3
    FF{d} = FF{d}(~Ind);
end
meanFF = mean([FF{1}; FF{2}; FF{3}]);

%% Histogram
directions = {'Left'; 'Middle'; 'Right'};
dirmeans = {};
for d = 1:3
    dirmeans{d} = mean(FF{d});

    figure(d);
    histogram(FF{d},'BinWidth',0.25);
    hold on;
    title(directions{d});
    ylim([0 40]); %80 50 70 40
end

%%
FFAll = [FF{1} FF{2} FF{3}];
histogram(FFAll,'BinWidth',0.25);
title('All');
mean(FFAll)
std(FFAll)

%% Feeding
FF= {};
% stab = {};
% s = {};
% figure;
for d = 1:8
    for u = 1:width(SpikesbyDir{d})
        NumCounts_vect = SpikesbyDir{d}(:,u);
        FF{d}(u) = (std(NumCounts_vect)^2)/mean(NumCounts_vect);
%         plot(mean(NumCounts_vect),var(NumCounts_vect),'o');
%         hold on;
    end
%     FF{d} = rmmissing(FF{d});
    [R,ind{d}] = rmmissing(FF{d});
%     stab{d} = FF{d} < 1;
%     s{d} = mean(stab{d})*100;
end
Ind = ind{1} | ind{2} | ind{3} | ind{4} | ind{5} | ind{6} | ind{7} | ind{8};
for d = 1:8
    FF{d} = FF{d}(~Ind);
end
meanFF = mean([FF{1}; FF{2}; FF{3}; FF{4}; FF{5}; FF{6}; FF{7}; FF{8}]);

%% Histogram
directions = {'-++', '+++', '--+', '+-+','-+-','++-','---','+--'};
dirmeans = [];
for d = 1:8
    dirmeans(d) = mean(FF{d});

%     figure(d);
%     histogram(FF{d},'BinWidth',0.25);
%     hold on;
%     title(directions{d});
%     ylim([0 40]); %50 70 20 40
end
mean(dirmeans)
std(dirmeans)

%%
FFAll = [FF{1} FF{2} FF{3} FF{4} FF{5} FF{6} FF{7} FF{8}];
figure;
histogram(FFAll,'BinWidth',0.25);
title('All');
mean(FFAll)
std(FFAll)

%% Kruskal-wallis test (directions)
%% Feeding
FFBoth= {};
for d = 1:8
    FFBoth{d} = [FFRyM1Con{d} FFYeM1Con{d}]';
%     FFBoth{d} = [meanFF_RyS1Con{d} meanFF_YeS1Con{d}]';
end

%%
a = cell(1,height(FFBoth{1})); a(:) = {'AntSupL'};
b = cell(1,height(FFBoth{2})); b(:) = {'AntSupR'};
c = cell(1,height(FFBoth{3})); c(:) = {'AntInfL'};
d = cell(1,height(FFBoth{4})); d(:) = {'AntInfR'};
e = cell(1,height(FFBoth{5})); e(:) = {'PostSupL'};
f = cell(1,height(FFBoth{6})); f(:) = {'PostSupR'};
g = cell(1,height(FFBoth{7})); g(:) = {'PostInfL'};
h = cell(1,height(FFBoth{8})); h(:) = {'PostInfR'};

vecI = [a b c d e f g h];
vecI = vecI'; %% directions

Vec = [FFBoth{1};FFBoth{2};FFBoth{3};FFBoth{4};FFBoth{5};FFBoth{6};FFBoth{7};FFBoth{8}]; %% trial x neuron
p = kruskalwallis(Vec,vecI');


%% Drinking
FFBoth= {};
for d = 1:3
    FFBoth{d} = [FFRyS1Con{d} FFYeS1Con{d}]';
end
%%
L = cell(1,height(FFBoth{1}));
L(:) = {'Left'};
M = cell(1,height(FFBoth{2}));
M(:) = {'Middle'};
R = cell(1,height(FFBoth{3}));
R(:) = {'Right'};
vecI = [L M R];
vecI = vecI'; %% directions

Vec = [FFBoth{1}; FFBoth{2}; FFBoth{3}]; %% trial x neuron
p = kruskalwallis(Vec,vecI');

%%
FFM1 = [FFBoth{1};FFBoth{2};FFBoth{3};FFBoth{4};FFBoth{5};FFBoth{6};FFBoth{7};FFBoth{8}];
%%
FFS1 = [FFBoth{1};FFBoth{2};FFBoth{3};FFBoth{4};FFBoth{5};FFBoth{6};FFBoth{7};FFBoth{8}];
%%
M = cell(1,height(FFM1));
M(:) = {'MIo'};
S = cell(1,height(FFS1));
S(:) = {'SIo'};
vecI = [M S];
vecI = vecI'; %% directions

Vec = [FFM1;FFS1]; %% trial x neuron
% p = kruskalwallis(Vec,vecI');
p = ranksum(FFM1,FFS1);

%% Churchland's mean-matched FF
%% Drinking
% Counts = RyDrM1Con;
Counts = [YeDrS1NB{1} YeDrS1NB{2} YeDrS1NB{3}];

times = 50:10:450;
params = struct('boxWidth',20);
Data = {};
data = struct();
for d = 1%:3
    trials = Counts;%{d};
    for t = 1:length(trials)
        allneurons = trials{t};
        for n = 1:height(allneurons)
            thisneuron = logical(allneurons(n,:));
            data(n).spikes(t,:) = thisneuron;
        end
    end
    Data{d} = data;
    Result{d} = VarVsMean(Data{d}, times);

%     outParams = plotScatter(Result, t, params);
%     outParams = plotFano(Result{d});
end

mmFFYeDrS1NB = Result;

%% Compare directions
FF = [mmFFRyDrM1Con{1}.FanoFactor mmFFRyDrM1Con{2}.FanoFactor mmFFRyDrM1Con{3}.FanoFactor];
p = kruskalwallis(FF);

%% Compare regions
M1 = [mmFFRyDrM1Con{1}.FanoFactor; mmFFRyDrM1Con{2}.FanoFactor; mmFFRyDrM1Con{3}.FanoFactor];
S1 = [mmFFRyDrS1Con{1}.FanoFactor; mmFFRyDrS1Con{2}.FanoFactor; mmFFRyDrS1Con{3}.FanoFactor];
[p,h] = ranksum(M1, S1);

%%
M1 = [mmFFRyFdS1Con{1}.FanoFactor; mmFFYeFdS1Con{1}.FanoFactor];
S1 = [mmFFRyDrS1Con{1}.FanoFactor; mmFFYeDrS1Con{1}.FanoFactor];
M = cell(1,height(M1));
M(:) = {'MIo'};
S = cell(1,height(S1));
S(:) = {'SIo'};
vecI = [M S];
vecI = vecI'; %% directions

Vec = [M1;S1]; %% trial x neuron
% p = kruskalwallis(Vec,vecI');
p = ranksum(M1,S1);
%%
[p,h] = ranksum(mmFFYeDrM1NB{1}.FanoFactor, mmFFYeDrS1NB{1}.FanoFactor);

%% Feeding
Counts = YeFdS1NB;

times = 20:10:80;
params = struct('boxWidth',20);
Data = {};
data = struct();
Result = {};
for d = 1%:8
    trials = Counts;%(catind{d});
%     trials = NBCounts(catind2{d});
    for t = 1:length(trials)
        allneurons = trials{t};
        for n = 1:height(allneurons)
            thisneuron = logical(allneurons(n,:));
            data(n).spikes(t,:) = thisneuron;
        end
    end
    Data{d} = data;
    Result{d} = VarVsMean(Data{d}, times);

%     outParams = plotScatter(Result, t, params);
%     outParams = plotFano(Result{d});
end

mmFFYeFdS1NB = Result;
%% Compare regions
M1 = [fdResults{1}.FanoFactor; fdResults{2}.FanoFactor; fdResults{3}.FanoFactor; fdResults{4}.FanoFactor; fdResults{5}.FanoFactor; fdResults{6}.FanoFactor; fdResults{7}.FanoFactor; fdResults{8}.FanoFactor];
S1 = [fdResults_S1{1}.FanoFactor; fdResults_S1{2}.FanoFactor; fdResults_S1{3}.FanoFactor; fdResults_S1{4}.FanoFactor; fdResults_S1{5}.FanoFactor; fdResults_S1{6}.FanoFactor; fdResults_S1{7}.FanoFactor; fdResults_S1{8}.FanoFactor];
[p,h] = ranksum(M1, S1);
%%
[p,h] = ranksum(mmFFRyFdM1NB{1}.FanoFactor, mmFFRyDrM1NB{1}.FanoFactor);
%% alt
spikeCounts = [sVRyM1FCon sVRyM1UCon];

times = 80:10:420;
params = struct('boxWidth',50);
Data = {};
data = struct();
for n = 1:length(spikeCounts)
    neuron = spikeCounts{n};
    for d = 1:length(neuron)
        trials = neuron{d};
        for t = 1:length(trials)
            
        end
    end
end

%% plot
figure; boxplot(Vec,vecI', 'symbol', '');
