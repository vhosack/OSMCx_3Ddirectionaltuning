%% Determine if neuron is directionally tuned & get PDs of tuned (feeding)
%% Determine if neuron is directionally tuned
p = [];
for u = 1:length(FRbyUnit)
    p(u) = kruskalwallis(FRbyUnit{u}',[],'off');
end
tuned = p < 0.05;
percentTuned = mean(tuned)*100;
%%
% NB
p = [];
for u = 1:length(FRbyUnit2)
    p(u) = kruskalwallis(FRbyUnit2{u}',[],'off');
end
tuned2 = p < 0.05;
percentTuned2 = mean(tuned2)*100;
%% Get meanfr of only tuned 
% Control
onlytuned = cell(length(FRbyUnit));
onlytuned(tuned) = FRbyUnit(tuned);
Tuned = onlytuned(~cellfun(@isempty, onlytuned));
meanfr = [];
for u = 1:length(Tuned)
    temp = Tuned{u};
    for d = 1:(length(directions)-1)
        meanfr(d,u) = sum(temp(d,:))/ncycles;
    end
end
FR = meanfr';
%%
% NB
onlytuned = cell(length(FRbyUnit2));
onlytuned(tuned2) = FRbyUnit2(tuned2);
Tuned2 = onlytuned(~cellfun(@isempty, onlytuned));
meanfr = [];
for u = 1:length(Tuned2)
    temp2 = Tuned2{u};
    for d = 1:(length(directions)-1)
        meanfr(d,u) = sum(temp2(d,:))/ncycles;
    end
end
FR_NB = meanfr';
%% Find PD 
% Control
pref = zeros(height(FR),1);
for unit = 1:height(FR)
    [p, d] = max(FR(unit,:));
    pref(unit) = d;
end
pref2 = pref;

% figure;
% histogram(pref,'BinLimits',[0.5,3.5],'normalization','probability','FaceAlpha',0.8,'FaceColor',[0.9290 0.6940 0.1250],...
%     'Normalization','probability','BinMethod','auto');
% xlabel('Preferred Direction');
% ylabel('Proportion of Neurons');
% title('MIo');
% set(gca,'FontSize',14,'YGrid','on');
%%
% Nerve Block
pref = zeros(height(FR_NB),1);
for unit = 1:height(FR_NB)
    [p, d] = max(FR_NB(unit,:));
    pref(unit) = d;
end
pref2 = pref;

% figure;
% histogram(pref,'BinLimits', [0.5,3.5],'normalization','probability','FaceAlpha',0.8,'FaceColor',[0.494117647058824 0.184313725490196 0.556862745098039],...
%     'Normalization','probability','BinMethod','auto');
% xlabel('Preferred Direction');
% ylabel('Proportion of Neurons');
% title('SIo');
% set(gca,'FontSize',14,'YGrid','on');