% This function, given and adjacency matrix of an undirected connected
% graph, obtains the transition probability matrix (Markov-chain)
%
% If you are using this efficiency_package for your research, plase kindly cite the paper:
%  
%  "Exploring the Morphospace of Communication Efficiency in Complex Networks"
%  Goñi J, Avena-Koenigsberger A, Velez de Mendizabal N, van den Heuvel M,
%  Betzel RF and Sporns O.
%  PLoS ONE. 2013
%
function P = f_markov_chain(adj)

N = size(adj,1);
k = full(sum(adj));
K = k(ones(N,1),:).*eye(N,N);
P = K\adj;