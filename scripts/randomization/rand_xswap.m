% function [rug,eff,diss] = xswap(ug,numChanges,maxIter,seed)
%
% randomization technique that preserves:
% - size and density
% - degree-sequence (and hence degree distribution)
%
% look for {i,k}{j,l} edges and change them for {i,l}{j,k} edges
% input parameters
% - ug: undirected adjacency matrix. Can be symmetric or upper-triangular
% - numChanges: number of effective "crossings" desired to be performed
% - maxIter: maximum number of attempts to satisfy numChanges
% - seed: specify seed for random-number-generator
% output paremeters
% - rug: randomized adjacency matrix in sparse upper-triangular format
% - eff: relative frequency of success randomizing
% - diss: dissimilarity between the original and the randomized adj matrices

% Joaquin Goni, IU
%version Jan 4th 2011
function [rug,eff,diss] = xswap(ug,numChanges,maxIter,seed)

if ~issparse(ug)
    ug = sparse(ug);
end

ug = triu(ug);
rug = ug;
L = nnz(rug); %number of edges of the graph
N = size(ug,1); %number of nodes of the graph

if nargin<3
    maxIter=numChanges*100;
end

if nargin==4
    stream0 = RandStream('mt19937ar','Seed',seed); % example of creating a seed: seed = sum(clock.*([1 1 1 1 1 1000]));
    RandStream.setDefaultStream(stream0);
end

iter = 0;
changes = 0;

[nodesA,nodesB] = find(rug==1); %original pairs of connected nodes (edges)    

while (changes < numChanges) && (iter<=maxIter)
    i = ceil(rand*N); %select i at random
    j = i; 
    while i==j
        j = ceil(rand*N); %select j at random
    end
    %%% random selection of node k within i-neighborhood
    nn = [nodesA(nodesB==i);nodesB(nodesA==i)];
    nn0 = ((nn==i)+(nn==j))==0; %neighborhood of i without {i,j}
    iNeigh = nn(nn0); %candidates to be k
    
    if ~isempty(iNeigh) %if there are candidates to become k 
        k = iNeigh(ceil(rand*numel(iNeigh)));
    
        %%% random selection of node l within j-neighborhood
        nn = [nodesA(nodesB==j);nodesB(nodesA==j)];
        nn0 = ((nn==i)+(nn==j)+(nn==k))==0; %neighborhood of j without {i,j,k}
        jNeigh = nn(nn0); %candidates to be l

        if ~isempty(jNeigh) %if there are candidates to become l 
            l = jNeigh(ceil(rand*numel(jNeigh)));

            if ~rug(min(i,l),max(i,l)) && ~rug(min(k,j),max(k,j))   %if edges {i,l} and {k,j} do not exist
                rug(min(i,k),max(i,k)) = 0; %remove edge {i,k} 
                rug(min(j,l),max(j,l)) = 0; %remove edge {j,l} 
                rug(min(i,l),max(i,l)) = 1; %add edge {i,l} 
                rug(min(j,k),max(j,k)) = 1; %add edge {j,k} 
                
                [nodesA,nodesB] = find(rug==1); %update pairs of connected nodes (edges)    
                changes = changes + 1;
            end   
        end
    end
    iter = iter + 1; 
end

if L~=nnz(rug)
    warning('density of the graph altered during randomization!')
end

diss = full(sum(abs(ug(:)-rug(:)))/(L*2)); %relative dissimilarity (value within range [0,1])
eff = changes / iter; %relative frequency of effective changes with respect to attempts

