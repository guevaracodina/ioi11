%% Anova Kruskal-Wallis, Tukey-Kramer
addpath(genpath('C:\spm8\toolbox\ioi11'));
% 2014 Batches J, K & L
somato = [...
0.83;
0.46;
0.58;
-0.13;
0.43;
0.79;

0.43;
0.11;
-0.35;
-0.65;
-0.25;
0.02;

-0.36;
0.05;
0.59;
0.75;
-0.44;
-0.02;
0.36;
-0.49];

motor = [...
-0.35;
0.66;
1.00;
-0.35;
-0.01;
0.43;

-0.17;
-0.31;
-0.36;
-0.80;
-0.17;
-0.71;

1.00;
-0.49;
-0.87;
0.33;
-0.34;
-0.65;
-0.74;
0.99];

% LPS = 1, LPS_dilation = 2, sham = 3
idxGroup = [...
1;
1;
2;
1;
1;
2;

3;
3;
3;
3;
3;
3;

3;
3;
3;
3;
2;
1;
2;
2];

idxCell = cell(size(idxGroup));
idxCell(idxGroup==1) = {'LPS'};
idxCell(idxGroup==2) = {'LPS_dilation'};
idxCell(idxGroup==3) = {'sham'};

%% Kruskal-Wallis for somato
% [P,ANOVATAB,STATS] = kruskalwallis([somato(idxGroup==1); somato(idxGroup==2); somato(idxGroup==3)], [ones(5,1); 2*ones(5,1); 3*ones(10,1)])
% [P,ANOVATAB,STATS] = kruskalwallis([somato(strcmp(idxCell,'LPS')); somato(strcmp(idxCell,'LPS_dilation')); somato(strcmp(idxCell,'sham'))], [ones(5,1); 2*ones(5,1); 3*ones(10,1)])
[P,ANOVATAB,STATS] = kruskalwallis(fisherz(somato), idxCell)
multcompare(STATS)

%% Kruskal-Wallis for Motor
[P,ANOVATAB,STATS] = kruskalwallis(fisherz(motor), idxCell)
multcompare(STATS)

