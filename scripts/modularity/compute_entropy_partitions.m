function H = compute_entropy_partitions(S)

numROIs = size(S,1);
numPartitions = size(S,2);
H = nan(numROIs,1);

for i=1:numROIs
    s = S(i,:);
    states = unique(s);
    h = 0;
    if length(states)>1    
        for j=1:length(states)
            pj = nnz(s==states(j))/numPartitions;
            h = h - (pj*log(pj));
        end
    end
    H(i) = h;
end

