function out = createSeeds(nSeeds,varargin)
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
    version = getoption(varargin,'version','225');
    fracSV = getoption(varargin,'fracSV',1);
    scanVar = getoption(varargin,'scanVar',1);
    savePy = getoption(varargin,'savepy',0);

    if runParallel
        parforArg = Inf;
    else
        parforArg = 0;
    end

    outdir = append('./', char(version), '_seeds');
    startIdx = length(dir([outdir, '/*.mat']));
    fprintf('start index %d\n', startIdx);

    eleName = {};
    eleAperture = {};

    switch version
        case '225S01'
            ring = AS2v225_S01;
            RM.RM1 = load('v225S01_RM1.mat').RM1;
            RM.RM2 = load('v225S01_RM2.mat').RM2;
            RM.MCO = load('v225S01_MCO.mat').MCO;

            ords.drift = find(atgetcells(ring,'FamName','Drift'));
            eleAperture.drift = 12.5e3 * [1 1]; % [m] drift
            ords.magnet = find(atgetcells(ring,'FamName','CF|CD|B1|BDS|^QDS|^QMS|^SX|^SF|^SD|^SH|^OC'));
            eleAperture.magnet = 12.5e3 * [1 1]; % [m] drift

            f_register = @register_AS2v2;
            buildTargetOrbit = 0;
        case '225'
            ring = AS2v225_15BPM_girderV4;
            RM.RM1 = load('v225_RM1.mat').RM1;
            RM.RM2 = load('v225_RM2.mat').RM2;
            RM.MCO = load('idealCORM_AS2v225_15BPM_14H19VCM.mat').MCO;

            ords.drift = find(atgetcells(ring,'FamName','Drift'));
            eleAperture.drift = 12.5e3 * [1 1]; % [m] drift
            ords.magnet = find(atgetcells(ring,'FamName','CF|CD|B1|BDS|^QDS|^QMS|^SX|^SF|^SD|^SH|^OC'));
            eleAperture.magnet = 12.5e3 * [1 1]; % [m] drift

            f_register = @register_AS2v2;
            buildTargetOrbit = 0;
        case 'assr'
            ring = assr4_splitbends;
            RM.RM1 = load('assr_RM1.mat').RM1;
            RM.RM2 = load('assr_RM2.mat').RM2;
            RM.MCO = load('assr_idealCORM.mat').MCO;

            ords.drift = find(atgetcells(ring,'PassMethod','DriftPass'));
            eleAperture.drift = 12.5e3 * [1 1]; % [m] drift
            ords.magnet = find(atgetcells(ring,'Class','Bend|Quadrupole|Sextupole'));
            eleAperture.magnet = 12.5e3 * [1 1]; % [m] drift

            f_register = @register_ASSR;
            buildTargetOrbit = 1;
        otherwise
            error('lattice version not supported');
    end

    SC = SCinit(ring);
    SC = f_register(SC);

    SC.INJ.beamSize = diag([20E-6, 10E-6, 10E-6, 5E-6, 1E-4, 1E-5].^2);

    SC.SIG.randomInjectionZ = [1E-5; 1E-6; 1E-5; 1E-6; 1E-5; 1E-5]; % [m; rad; m; rad; rel.; m]
    SC.SIG.staticInjectionZ = [1E-4; 1E-5; 1E-4; 1E-5; 1E-4; 1E-4]; % [m; rad; m; rad; rel.; m]

    SC.SIG.Circumference = 0; % relative
    SC.BPM.beamLostAt    = 0.6;  % relative


    for ord=ords.drift'
        SC.RING{ord}.EApertures = eleAperture.drift; % [m] drift
    end

    for ord=ords.magnet'
        SC.RING{ord}.EApertures = eleAperture.magnet; % [m] magnet apertures
    end

    SCsanityCheck(SC);

    results = {};
    parfor (seed = startIdx:nSeeds+startIdx, parforArg)
        if verbose
            fprintf('Processing seed %d...\n',seed);
        end
        % clear results
        runTime = tic;
        shuffleRNG(seed);

        try 
            % [~,results] = evalc('runSingleSeed(SC,''RM'',RM)');
            newSeed = createOneSeed(SC,'runCorrection',runCorrection,'MCO',RM.MCO,'fracSV',fracSV,'buildTargetOrbit',buildTargetOrbit);

        catch ME
            fprintf('Seed %d failed\n',seed);
            disp(ME.message);
            continue
        end

        % results{seed} = newSeed;

        % convert from matlab AT to pyAT
        %
        if ~isfolder(outdir)
            mkdir(outdir);
        end

        if savePy
            seedRing = newSeed.preCorrection;
            outfile = sprintf('%s/%s_seed%d_preCorrection_pyAT',outdir,version,seed); 
            atwritepy(seedRing,'file',outfile);
        end

        % Leo: I've integrated the second for-loop into this one so that we
        % don't need to save the list of generated seeds anymore. It was causing
        % a memory leak that made it infeasible to run for large nSeeds values. 

        % The transparency violation was solved by putting the save() into a 
        % separate function. credit: random stackoverflow guy
        outfile = sprintf('%s/seed%d.mat',outdir,seed); 
        saveSeed(outfile, newSeed);

        if runCorrection && ~isa(newSeed.postCorrection, 'char') && savePy
            seedRing = newSeed.postCorrection;
            outfile = sprintf('%s/%s_seed%d_postCorrection_pyAT',outdir,version,seed); 
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
    fracSV = getoption(varargin,'fracSV',1);
    scanVar = getoption(varargin,'scanVar',1);
    buildTargetOrbit = getoption(varargin,'buildTargetOrbit',1);

    newSeed = struct();

    validSeed = 0;
    attemptNum = 0;
    while ~validSeed
        attemptNum = attemptNum + 1;
        if attemptNum > 10
            error('Unable to find closed orbit for uncorrected ring after 10 seeds');
        end
        % repeat seed generation until closed orbit can be found
        %
        SC_seed = SCapplyErrors(SC);
        newRing = SC_seed.RING;
        [~,T] = evalc('findorbit6(newRing)');
        validSeed = ~any(isnan(T));
    end
    
    newSeed.preCorrection = SC_seed.RING;
    % for diagnostics
    %
    % newSeed.SCpreCorrection = SC_seed;

    if runCorrection
        SC = SC_seed;
        SC.INJ.trackMode = 'ORB';

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
        warning('off','MATLAB:illConditionedMatrix');
        [SC,COexists] = runOrbitCorrection_AS2(SC,MCO,eta,'etaWeight',1,'fracSV',fracSV,'buildTargetOrbit',buildTargetOrbit);
        warning('on','MATLAB:illConditionedMatrix');

        if ~COexists
            newSeed.postCorrection = 'Correction failed';
        else
            newSeed.postCorrection = SC.RING;
            newSeed.SCpostCorrection = SC;
        end
    end

end

function makeResponseMatrices(varargin)
    SC = varargin{1};
    outdir = varargin{2};
    version = varargin{3};

    RM1 = SCgetModelRM(SC,SC.ORD.BPM,SC.ORD.CM,'nTurns',1,'useIdealRing',1);
    RM2 = SCgetModelRM(SC,SC.ORD.BPM,SC.ORD.CM,'nTurns',2,'useIdealRing',1);
    MCO = SCgetModelRM(SC,SC.ORD.BPM,SC.ORD.CM,'trackMode','ORB','useIdealRing',1);

    outfile = sprintf('%s/%s_RM1.mat',outdir,version);
    save(outfile,'RM1');
    outfile = sprintf('%s/%s_RM2.mat',outdir,version);
    save(outfile,'RM2');
    outfile = sprintf('%s/%s_idealCORM.mat',outdir,version);
    save(outfile,'MCO');
end