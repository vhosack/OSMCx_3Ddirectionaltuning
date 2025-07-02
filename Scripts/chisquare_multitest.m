%% Compare percent tuned between groups (Chi-square)
% save tuned for 2 groups to compare
%% MIo vs SIo
M1 = tunedM1';
Mind = cell(1,length(tunedM1));
Mind(:) = {'M1'};

SC = tunedSC';
Sind = cell(1,length(tunedSC));
Sind(:) = {'SC'};
y1 = [M1; SC];
y2 = [Mind Sind];
y2 = y2';

% Chi-square test
[tbl,chi2stat,pval] = crosstab(y1,y2)
%% Control vs nerve block (feeding)
Con = tuned';
Cind = cell(1,length(tuned));
Cind(:) = {'Con'};

NB = tuned2';
NBind = cell(1,length(tuned2));
NBind(:) = {'NB'};
y1 = [Con; NB];
y2 = [Cind NBind];
y2 = y2';

% Chi-square test
[tbl,chi2stat,pval] = crosstab(y1,y2)
%% Control vs nerve block (drinking)
Con = tunedC';
Cind = cell(1,length(tunedC));
Cind(:) = {'Con'};

NB = tunedNB';
NBind = cell(1,length(tunedNB));
NBind(:) = {'NB'};
y1 = [Con; NB];
y2 = [Cind NBind];
y2 = y2';

% Chi-square test
[tbl,chi2stat,pval] = crosstab(y1,y2)