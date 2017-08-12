% function [rug,eff,diss] = flip(ug,numChanges,maxIter,seed)
%
% randomization technique that preserves:
% - size and density
% - connectedness 
% - degree-sequence (and hence degree distribution)
%
% ensuring that {i,j} exists, look for {i,k}{j,l} and change them for {i,l}{j,k} edges
% input parameters
% - ug: undirected adjacency matrix. Can be symmetric or upper-triangular
% - numChanges: number of effective "crossings" desired to be performed
% - maxIter: maximum number of attempts to satisfy numChanges
% - seed: specify seed for random-number-generator
% output paremeters
% - rug: randomized adjacency matrix in sparse upper-triangular format
% - eff: relative frequency of success randomizing
% - diss: dissimilarity between the original and the randomized adj matrices

% Joaquin Goni - IU
%version Jan 4th 2011
function [rug,eff,diss] = flip(ug,numChanges,maxIter,seed)

if ~issparse(ug)
    ug = sparse(ug);
end

ug = triu(ug);
rug = ug;
L = nnz(rug); %number of edges of the graph

degree_ug = full(sum(ug+ug'));
degree_rug = degree_ug;
if nargin<3
    maxIter = numChanges*100;
end

if nargin==4
    stream0 = RandStream('mt19937ar','Seed',seed); % example of creating a seed: seed = sum(clock.*([1 1 1 1 1 1000]));
    RandStream.setDefaultStream(stream0);
end

iter = 0;
changes = 0;

[nodesA,nodesB] = find(rug==1); %original pairs of connected nodes (edges)

while (changes < numChanges) && (iter <= maxIter)
    
    edge = ceil(rand*L); %select an existing edge {i,j} at random
    i = nodesA(edge); %identify i
    j = nodesB(edge); %identify j
    
    if rand>0.5
        aux=i;
        i=j;
        j=aux;
    end
    
    deg_j = degree_rug(j);
    k = find(degree_rug==(deg_j-1));
    kk0 = ((k==i)+(k==j))==0;
    k=k(kk0);
    if ~isempty(k)
        k=k(ceil(rand*length(k)));
        if ~rug(min(i,k),max(i,k))
            rug(min(i,j),max(i,j)) = 0;
            rug(min(i,k),max(i,k)) = 1;
            [nodesA,nodesB] = find(rug==1); %update original pairs of connected nodes (edges)
            degree_rug = full(sum(rug+rug'));
            changes =  changes + 1;
        end
    end
    iter = iter + 1; 
end

if L~=nnz(rug)
    warning('density of the graph altered during flip randomization')
end

if  any(sort(degree_ug)-sort(degree_rug))
    warning('degree distribution altered during flip randomization')
end

diss = full(sum(abs(ug(:)-rug(:)))/(L*2)); %relative dissimilarity (value within range [0,1])
eff = changes / iter; %relative frequency of effective changes with respect to attempts


