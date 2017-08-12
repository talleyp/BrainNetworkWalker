%S: vector Nx1 with the partition found for adj matrix A
%Q: modularity score of the partition found for adj matrix A
%B: modularity matrix
function [S,Q,B,twom] = front_genlouvain_single(A,gamma)

if nargin<2
    gamma=1;
end

N = size(A,1);
numGammaValues = size(gamma,1);

A = double(A.*~eye(size(A)));

if any(A-A')
    error('matrix A is not symmetric!')
end

%% initialization
Q = nan(numGammaValues,1);
S = nan(N,numGammaValues);

%% signed matrix
if any(A(:)<0)
    if size(gamma,2)==2
        gamma_pos = gamma(:,1);
        gamma_neg = gamma(:,2);
    elseif size(gamma,2)==1
        gamma_pos = gamma;
        gamma_neg = gamma;
    else
        error('gamma cannot contain more than two columns')
    end
    %% masks for signs
    mask_pos = A>0;
    mask_neg = A<0;
    %% degrees
    k_pos = sum(A.*mask_pos);
    k_neg = sum(-A.*mask_neg);
    %% densitites (term 2*m)
    twom_pos = sum(k_pos);
    twom_neg = sum(k_neg);
    twom = twom_pos + twom_neg;

    for i=1:numGammaValues
        %% construction of null model (P)
        P = (((k_pos'*k_pos)/twom_pos)*gamma_pos(i)) - (((k_neg'*k_neg)/twom_neg)*gamma_neg(i)); %null model for signed networks. SupMat, page 5, eq 7.
        B = (A.*mask_pos) - (-A.*mask_neg) - P;
        [S(:,i),Q(i,1)] = genlouvain(B);
    end
    Q = Q./(twom);

%% unsigned matrix
else
    %% degrees
    k = sum(A);
    %% density (term 2*m)
    twom = sum(k);
    for i=1:numGammaValues
        %% construction of null model (P)
        P = ((k'*k)/twom)*gamma(i);
        B = A-P;
        [S(:,i),Q(i,1)] = genlouvain(B);
    end
    Q = Q./twom;
end

