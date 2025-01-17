function shuffleRNG(jobNumber)
    % script to generate matlab rng seed for parallel and batch jobs
    % since <rng('shuffle')> will automatically generate seed from current time
    % it is possible to encounter the same seed when batch jobs are created on computing clusters
    % This script to manually assign a seed based on current time and job number
    %
    c = clock;
    c = fix(c);
    x = sprintf('%d',c(2:end));
    x = str2num(x);
    seed = x+fix(jobNumber);
    fprintf('RNG seed: %d\n',seed);
    rng(seed);

end