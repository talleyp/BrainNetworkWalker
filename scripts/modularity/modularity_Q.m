function Q = modularity_Q(twom,B,partition)

N=length(partition);

if size(B,1)~=N || size(B,2)~=N
    error('B and partition sizes do not match!')
end

mask = false(N,N);
for i=1:N-1
    module_i = partition(i);
    for j=i+1:N
        mask(i,j) = module_i == partition(j); % assign if nodes i and j are in the same partition
    end
end

mask = mask | mask';

Q = sum(B(mask))/twom;

