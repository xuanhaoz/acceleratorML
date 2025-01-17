function createSeeds(nSeeds)
    addpath('SC');
    addpath('lattices');
    addpath('responseMatrices');

    global plotFunctionFlag
    global verbose
    plotFunctionFlag = 0;
    verbose = 0;

    outdir = './seeds';

    ring = AS2v225_15BPM_girderV4;
    RM.RM1 = load('v225_RM1.mat').RM1;
    RM.RM2 = load('v225_RM2.mat').RM2;

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

    for seed = 1:nSeeds
        fprintf('Processing seed %d...\n',seed);
        clear result
        runTime = tic;
        shuffleRNG(seed);

        try 
            [~,results] = evalc('runSingleSeed(SC,''RM'',RM)');
        catch ME
            fprintf('Seed %d failed\n',seed);
            disp(ME.message);
            continue
        end

        seedRing = results.RFcorrection.RING;

        % convert from matlab AT to pyAT
        %
        if ~isfolder(outdir)
            mkdir(outdir);
        end
        outfile = sprintf('%s/seed%d_pyAT',outdir,seed); 
        atwritepy(seedRing,'file',outfile);

        runTime = toc(runTime);
        fprintf('Run time: %.1f mins\n',runTime/60);
    end
end