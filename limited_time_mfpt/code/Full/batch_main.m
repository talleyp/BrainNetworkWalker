% clearvars;
close all;
clc;

path.data = fullfile(pwd,'data','subject');

load(fullfile(path.data,'FC.mat'));
load(fullfile(path.data,'SCden.mat'));
load(fullfile(path.data,'SCfa.mat'));

figure, spy(SCden); title('SC (binary)'); axis square;
figure, imagesc(log(SCden)); title('SC (weights)'); axis square;
figure, imagesc(FC,[-0.6,0.6]); axis square; colorbar;

N = size(FC,1);
Pden = SCden./repmat(sum(SCden,2),1,N);
fraction_arrival_full = nan(N,N);
lt_mfpt_full = nan(N,N);

% source = 12;
% target = 57;
numSteps = 1000;
numRuns = 100;

% loop here for every pair of source and target
for source=1:N
  for target=1:N
    [fraction_arrival,lt_mfpt] = run_random_walk(Pden,source,target,
                                   numSteps,numRuns);
    %fraction_arrival_full(source, target) = fraction_arrival;
    %lt_mfpt_full(source, target) = lt_mftp;    
    end
end

figure, plot(lt_mfpt); xlabel('time'); ylabel('lt_mfpt');
figure, plot(fraction_arrival); xlabel('time'); ylabel('fraction arrival');