close all;
clearvars;
clc;

set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultAxesFontSize',18);

set(0,'DefaultTextFontname','Times New Roman');
set(0,'DefaultTextFontSize',18);

%% add folders with scripts to path
path2scripts = fullfile(pwd,'..');
addpath(genpath(path2scripts))

%% load network
load brainSC.mat % load here your favourite network

N=size(A,1); %number of nodes in the network
mask = triu(true(N,N),1); %upper triangular binary mask;

figure;
%% degree: number of connections of each node
subplot(1,2,1);
degree = sum(A>0);
bar(degree); xlabel('nodes'); ylabel('degree'); xlim([1,N]);

%% strength: total weight for all connections of each node
subplot(1,2,2);
strength = sum(A);
bar(strength); xlabel('nodes'); ylabel('strength'); xlim([1,N])


%% integration - shortest paths - efficiency - search-information

[SP.distance,SP.length,SP.B] = get_shortest_path_lengths(1./A); %distance can be expressed as -log(A), as 1./A and in other ways...
SP.I = get_information_shortest_paths_wei_und(A,SP.length,SP.B,sum(A),1); %compute pairwise search-information
%compute pairwise search-information
Isym=(SP.I+SP.I')/2; % search information is not a symmetric measurement
Isym(logical(eye(size(SP.I))))=0; %set main diagonal to zero

figure,
subplot(1,3,1); imagesc(SP.length); axis square; title('SP (edges)'); colorbar;
xlabel('nodes'); ylabel('nodes');
subplot(1,3,2); imagesc(SP.distance); axis square; title('SP (distance)'); colorbar;
xlabel('nodes'); ylabel('nodes');
subplot(1,3,3); imagesc(SP.I); axis square; title ('SP (search-information)'); colorbar;
xlabel('nodes'); ylabel('nodes');

%% diffusion: mean first passage time.

MFPT =f_mfpt(f_markov_chain(A)); % mean-first-passage-times computation
Ediff = mean(1./MFPT(mask)); % diffusion efficiency.
figure, imagesc(MFPT); axis square; colorbar;
xlabel('nodes'); ylabel('nodes');

%% modularity

gamma = 1;
[S,Q] = front_genlouvain_single(A,gamma);
[aa,bb] = sort(S);
figure, 
subplot(1,2,1); spy(A);
title('network original ordering')
subplot(1,2,2); spy(A(bb,bb));
title(sprintf('network reordered, %d modules, Q=%0.2f',length(unique(S)),Q))

%% randomization

numChanges = 500; %number of effective changes
maxIter = numChanges*100; %number of attemps
seed = 1; %seed for random-numbers generator
[randA,eff,diss] = rand_xswap_wu(A,numChanges,maxIter,seed);

figure,
subplot(1,2,1); spy(A); title('original graph')
subplot(1,2,2); spy(randA); title('randomized graph')

%% randomization experiment and assessment of network features

% step 1
% create a loop based on xswap iterations
% every 100 iterations up to 10,000, evaluate Q,Erout,Ediff,Einf
% make 4 plots showing the behavior as graph-randomization evolves

% step 2, create ensembles of randomized networks for each number of
% iterations and report mean and std-error of network measures

