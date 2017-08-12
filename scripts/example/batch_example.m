
close all;
clear all;
clc;

load SC.mat; % let's assume SC is your weighted connectome (diagonal set to zeros, symmetric matrix)

N=size(SC,1); %number of ROIs
mask = triu(true(N,N),1); %upper triangular binary mask

%% MFPT
MFPT =f_mfpt(f_markov_chain(SC)); % mean-first-passage-times computation
Ediff = mean(1./MFPT(mask)); % diffusion efficiency.

%% Search-Information
[DSPL SPL B] = get_shortest_path_lengths(-log(SC)); %distance can be expressed as -log(SC), as 1./SC and in other ways...
Iraw = get_information_shortest_paths_wei_und(SC,SPL,B,sum(SC),1); %compute pairwise search-information
I=(Iraw+Iraw')/2; % search information is not a symmetric measurement
I = log(I); %log-transform of the symmetrized SI values
I(logical(eye(size(I))))=0; %set main diagonal to zero

%% Path-transitivity
MIsp = get_shortest_path_matching_index_wei_und(SC,SPL,B);
