function [fraction_arrival,lt_mfpt]= run_random_walk(P,source,target,numSteps,numRuns)

walker = nan(numRuns,numSteps);
for run=1:numRuns
    %disp(run)
    walker(run,1) = source;
    step=2;
    while (step<=numSteps) && (walker(run,step-1)~=target)
        p = P(walker(run,step-1),:);
        p = cumsum(p);
        movement = rand;
        aux = p - movement;
        walker(run,step) = find(aux>=0,1);
        step = step + 1;
    end
end

lt_mfpt = nan(numSteps,1);
fraction_arrival = nan(numSteps,1);
for time=1:numSteps
    walker_limited=walker(:,1:time);
    [rows,cols]=find(walker_limited==target);
    fraction_arrival(time) = length(rows)/numRuns;
    if ~isempty(cols)
        lt_mfpt(time) = mean(cols);
    end
end
%figure, imagesc(walker); xlabel('time'); ylabel('runs'); colorbar;
%figure, plot(lt_mfpt); xlabel('time'); ylabel('lt_mfpt');
%figure, plot(fraction_arrival); xlabel('time'); ylabel('fraction arrival');



    


