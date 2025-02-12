function createSeeds(nSeeds,varargin)
    addpath('SC');
    addpath('lattices');
    addpath('responseMatrices');

    global plotFunctionFlag
    global verbose
    plotFunctionFlag = 0;
    verbose = 1;

    % <runParallel = 1> requires parallel computing toolbox in matlab, automatically spawns parallel processes according to the number cores available
    % <runCorrection = 1> will perform closed orbit correction on randomly initialised seed, this takes ~2.5 mins per seed
    % without runCorrection, takes about 15 seconds per seed
    %
    runParallel = getoption(varargin,'parallel',0);
    runCorrection = getoption(varargin,'correction',1);

    if runParallel
        parforArg = 4;
    else
        parforArg = 0;
    end

    outdir = './seeds';

    ring = AS2v225_15BPM_girderV4;
    RM.RM1 = load('v225_RM1.mat').RM1;
    RM.RM2 = load('v225_RM2.mat').RM2;
    RM.MCO = load('idealCORM_AS2v225_15BPM_14H19VCM.mat').MCO;

    SC = SCinit(ring);
    SC = register_AS2v2(SC);

    SC.INJ.beamSize = diag([20E-6, 10E-6, 10E-6, 5E-6, 1E-4, 1E-5].^2);

    SC.SIG.randomInjectionZ = [1E-5; 1E-6; 1E-5; 1E-6; 1E-5; 1E-5]; % [m; rad; m; rad; rel.; m]
    SC.SIG.staticInjectionZ = [1E-4; 1E-5; 1E-4; 1E-5; 1E-4; 1E-4]; % [m; rad; m; rad; rel.; m]

    SC.SIG.Circumference = 0; % relative
    SC.BPM.beamLostAt    = 0.6;  % relative


    for ord=SCgetOrds(SC.RING, 'Drift')
        SC.RING{ord}.EApertures = 12.5E+3 * [1 1]; % [m] drift
    end

    for ord=SCgetOrds(SC.RING, 'CF|CD|B1|BDS|^QDS|^QMS|^SX|^SF|^SD|^SH|^OC')
        SC.RING{ord}.EApertures = 12.5E+3 * [1 1]; % [m] magnet apertures
    end

    SCsanityCheck(SC);

    % results = {};
    parfor (seed = 1:nSeeds, parforArg)
        if verbose
            fprintf('Processing seed %d...\n',seed);
        end
        % clear results
        runTime = tic;
        shuffleRNG(seed);

        try 
            % [~,results] = evalc('runSingleSeed(SC,''RM'',RM)');
            newSeed = createOneSeed(SC,'runCorrection',runCorrection,'MCO',RM.MCO);

        catch ME
            fprintf('Seed %d failed\n',seed);
            disp(ME.message);
            continue
        end

        % convert from matlab AT to pyAT
        %
        if ~isfolder(outdir)
            mkdir(outdir);
        end
        seedRing = newSeed.preCorrection;
        outfile = sprintf('%s/seed%d_preCorrection_pyAT',outdir,seed); 
        atwritepy(seedRing,'file',outfile);
        
        % Leo: I've integrated the second for-loop into this one so that we
        % don't need to save the list of generated seeds anymore. It was causing
        % a memory leak that made it infeasible to run for large nSeeds values. 

        % The transparency violation was solved by putting the save() into a 
        % separate function. credit: random stackoverflow guy
        outfile = sprintf('%s/seed%d.mat',outdir,seed); 
        saveSeed(outfile, newSeed);

        if runCorrection
            seedRing = newSeed.postCorrection;
            outfile = sprintf('%s/seed%d_postCorrection_pyAT',outdir,seed); 
            atwritepy(seedRing,'file',outfile);
        end

        runTime = toc(runTime);
        if verbose
            fprintf('Seed %d run time: %.1f mins\n',seed,runTime/60);
        end
    end

end

function saveSeed(fname, seedToSave)
    save(fname, '-struct', 'seedToSave');
end

function newSeed = createOneSeed(varargin)
    SC = varargin{1};
    runCorrection = getoption(varargin,'runCorrection',0);
    MCO = getoption(varargin,'MCO',[]);

    newSeed = struct();

    validSeed = 0;
    while ~validSeed
        % repeat seed generation until closed orbit can be found
        %
        SC_seed = SCapplyErrors(SC);
        newRing = SC_seed.RING;
        [~,T] = evalc('findorbit6(newRing)');
        validSeed = ~any(isnan(T));
    end
    newSeed.preCorrection = SC_seed.RING;

    if runCorrection
        SC = SC_seed;
        SC.INJ.trackMode;

        BPMords = SC.ORD.BPM;
        CMords = SC.ORD.CM;

        % check supplied MCO has the correct size
        %
        if ~(numel(MCO) == length(BPMords)*2 * (length(CMords{1}) + length(CMords{2})))
            if verbose
                fprintf('MCO does not match required size, calculating ideal MCO\n');
            end
            MCO = SCgetModelRM(SC,BPMords,CMords,'trackMode','ORB','useIdealRing',1);
        end

        eta = SCgetModelDispersion(SC,BPMords,SC.ORD.Cavity,'rfStep',5);
        [SC,COexists] = runOrbitCorrection_AS2(SC,MCO,eta,'etaWeight',1e3,'fracSV',0.6);

        if ~COexists
            newSeed.postCorrection = 'Correction failed';
        else
            newSeed.postCorrection = SC.RING;
        end
    end

end