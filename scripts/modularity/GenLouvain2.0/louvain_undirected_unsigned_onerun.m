function [S,Q,order] = louvain_undirected_unsigned_onerun(A,gamma)

if nargin<2
    gamma = 1;
end

N = size(A,1); % number of nodes
k = full(sum(A)); % degree vector
twom = sum(k); % sum of degrees is 2m
if N<=500 % small network
    B = full(A - gamma*k'*k/twom);
else % large network
    B = @(i) A(:,i) - gamma*k'*k(i)/twom;
end
[S,Q] = genlouvain(B); % run Louvain
Q = Q/twom; % normalize value of Q
[~,order] = sort(S); % order for rows and columns according to modularity
