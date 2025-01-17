function results = runSingleSeed(SC,varargin)
    % function to run all correction steps for one seed candidate
    % contained in this way to avoid any single seed error stopping the whole job
	global plotFunctionFlag

	% lattice = 'AS2v202_15BPM_14H19VCM';
	% lattice = 'AS2v204_15BPM_14H19VCM';
	% lattice = 'AS2v214_15BPM_14H19VCM';
	% lattice = 'AS2v224_15BPM_14H19VCM';
	lattice = 'AS2v225_15BPM_14H19VCM';

    p = inputParser;
    addOptional(p,'old',0);
    addOptional(p,'fracSV',0.6);
    addOptional(p,'version',2);
    addOptional(p,'RM',0);
    addOptional(p,'nTurns',2000);
    addOptional(p,'locoIterations',4);
    addOptional(p,'CORM_name',lattice);
    addOptional(p,'injection',0);
    parse(p,varargin{:});
    par = p.Results;

	% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	% CHECK THIS BEFORE FIRING OFF JOBS
	%
	BPM2QuadPairing_version = 2;
	% nBPM = length(SC.ORD.BPM)/24;	% AS2v224_15BPM_14H19VCM
	nBPM = 19;						% AS2v224_15BPM_inj
	% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	results = struct();
	results.name = lattice;
	% ------------------------------------------------------------------------------
	% apply errors
	%
	if isstruct(par.old)
		SC = duplicateErrorSeed(par.old,SC);
	else
		SC = SCapplyErrors(SC);
	end

	if par.injection
		SC = copyGirderError(SC);
	end

	results.preCorrection = SC;

	% plot overall offsets and rolls of elements plus support girders,
	% requires defining <SCregisterSupport>
	%
	% SCplotSupport(SC);

	% ------------------------------------------------------------------------------
	% setup correction chain
	% 1. turn off cavity
	% 2. switch off sextupoles
	%
	SC.RING = SCcronoff(SC.RING,'cavityoff');

	BPMords = SC.ORD.BPM;
	CMords = SC.ORD.CM;
	
	switch par.version
		case 1
			sextOrds = SCgetOrds(SC.RING,'SF|SD|SX1|SX2|SX3|SH1|SH2|SH4|SH5');
			SC = SCsetMags2SetPoints(SC,sextOrds,2,3,0,'method','abs');
		case 2
			sextOrds = SCgetOrds(SC.RING,'^SF|^SD');
			SC = SCsetMags2SetPoints(SC,sextOrds,2,3,0,'method','abs');
			octOrds = SCgetOrds(SC.RING,'^OC');
			SC = SCsetMags2SetPoints(SC,octOrds,2,4,0,'method','abs');
		otherwise
			error('Unsupported lattice version.')
	end

	SC.INJ.nParticles = 1;
	SC.INJ.nTurns     = 1;
	SC.INJ.nShots     = 1;
	SC.INJ.trackMode  = 'TBT';

	% BPM noise level
	%
	par.eps =  1E-6;

	% turn by turn orbit target at BPM
	%
	par.TbTtarget = 100E-6;

	% tracking parameters 
	%
	par.transmissionNturns = 100;
	par.transmissionNparticles = 10;

	SCgetBPMreading(SC);

	% ------------------------------------------------------------------------------
	% setup one turn matrix and pseudo inverse matrix
	%
	% RM1 = SCgetModelRM(SC,SC.ORD.BPM,SC.ORD.CM,'nTurns',1);
	% RM1 = SCgetModelRM(SC,SC.ORD.BPM,SC.ORD.CM,'nTurns',1,'useIdealRing',0);
	RM1 = par.RM.RM1;
	Minv1 = SCgetPinv(RM1,'alpha',100);

	[CUR,ERROR] = SCfeedbackFirstTurn(SC,Minv1,'verbose',1);
	if ~ERROR; SC=CUR; else; return; end

	[CUR,ERROR] = SCfeedbackRun(SC,Minv1,'target',par.TbTtarget,'maxsteps',1000,'eps',par.eps,'verbose',1);
	if ~ERROR; SC=CUR; else; return; end
	
	results.oneTurn = SC;

	% ------------------------------------------------------------------------------
	% multi turn 
	%
	SC.INJ.nTurns = 2;
	% RM2 = SCgetModelRM(SC,SC.ORD.BPM,SC.ORD.CM,'nTurns',2);
	% RM2 = SCgetModelRM(SC,SC.ORD.BPM,SC.ORD.CM,'nTurns',2,'useIdealRing',0);
	RM2 = par.RM.RM2;
	Minv2 = SCgetPinv(RM2,'alpha',100);

	for nBPMincl=[15 10 5 3]
		[CUR,ERROR] = SCfeedbackStitch(SC,Minv2,'nBPMs',3,'maxsteps',1000,'verbose',1);
		if ~ERROR; SC=CUR; else; disp(nBPMincl); return; end
	end

	[CUR,ERROR] = SCfeedbackRun(SC,Minv2,'target',par.TbTtarget,'maxsteps',1000,'eps',par.eps,'verbose',1);
	if ~ERROR; SC=CUR; else; return; end

	[CUR,ERROR] = SCfeedbackBalance(SC,Minv2,'maxsteps',1000,'eps',par.eps,'verbose',1);
	if ~ERROR; SC=CUR; else; return; end
	
	results.twoTurn = SC;

	[maxTurns,lostCount,ERROR] = SCgetBeamTransmission(SC,'nParticles',par.transmissionNparticles,'nTurns',par.transmissionNturns,'verbose',true);
	if ERROR;return;end

	% ------------------------------------------------------------------------------
	% BBA
	%
	% quadOrds = getBPM2QuadPairing_AS2(SC,repmat(BPMords,2,1),'nBPM',nBPM,'version',BPM2QuadPairing_version);
	quadOrds = getBPM2BBAPairing_AS2(SC,repmat(BPMords,2,1),'nBPM',nBPM);
	BPMords  = repmat(SC.ORD.BPM,2,1);
	SC       = SCpseudoBBA(SC,repmat(BPMords,2,1),quadOrds,10E-6,'sigma',1);

	% find two turn again
	%
	Minv2 = SCgetPinv(RM2,'alpha',10);

	[CUR,ERROR] = SCfeedbackRun(SC,Minv2,'target',par.TbTtarget,'maxsteps',1000,'eps',par.eps,'verbose',1);
	if ~ERROR; SC=CUR; else; return; end

	[CUR,ERROR] = SCfeedbackBalance(SC,Minv2,'maxsteps',100,'eps',par.eps,'verbose',1);
	if ~ERROR; SC=CUR; else; return; end

	results.BBA = SC;

	% B1rms = sqrt(mean(SCgetBPMreading(SC ,'BPMords',BPMords).^2,2));
	% fprintf('Pre BBA reading %.2fum\n Post BBA reading %.2fum\n',1E6*B0rms,1E6*B1rms);

	% % break if no beam capture
	% %
	% [maxTurns,lostCount,ERROR] = SCgetBeamTransmission(SC,'nParticles',par.transmissionNparticles,'nTurns',par.transmissionNturns,'plotFlag',0,'verbose',true);
	% if ERROR;return;end

	% % if maxTurns < 50 | any(lostCount)
	% % 	break
	% % end

	% ------------------------------------------------------------------------------
	% ramp up sextupoles and feedback is applied after each step.
	%
	nSteps=10;
	for S = linspace(1/nSteps,1,nSteps)
		switch par.version
			case 1
				SC = SCsetMags2SetPoints(SC,sextOrds,2,3,S,'method','rel');
			case 2
				SC = SCsetMags2SetPoints(SC,sextOrds,2,3,S,'method','rel');
				SC = SCsetMags2SetPoints(SC,octOrds,2,4,S,'method','rel');
		end

		[CUR,ERROR] = SCfeedbackBalance(SC,Minv2,'maxsteps',50,'eps',par.eps,'verbose',1);
		if ~ERROR; SC=CUR; end
	end

	% SCplotPhaseSpace(SC,'nParticles',10,'nTurns',100);

	[maxTurns,lostCount,ERROR] = SCgetBeamTransmission(SC,'nParticles',par.transmissionNparticles,'nTurns',100,'plotFlag',0,'verbose',true);
	if ERROR;return;end
	
	results.sextOn = SC;

	SCgetBPMreading(SC);

	% ------------------------------------------------------------------------------
	% chromaticity correction
	%
	% SF = SCgetOrds(SC.RING,'SF');
	% SD = SCgetOrds(SC.RING,'SD');
	% sOrds = {SF SD};
	% SC = SClocoLib('fitChromaticity',SC,sOrds,'verbose',1);

	% ------------------------------------------------------------------------------
	% rf correction: phase and frequency
	%
	SC = runRFCorrection_AS2(SC);

	results.RFcorrection = SC;


	% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	% finish here for post RF correction
	% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	%{
	% Plot final phase space and check if beam capture is achieved.
	%
	% SCplotPhaseSpace(SC,'nParticles',10,'nTurns',100);

	% [maxTurns,lostCount,ERROR] = SCgetBeamTransmission(SC,...
	% 	'nParticles',par.transmissionNparticles,...
	% 	'nTurns',par.transmissionNturns,...
	% 	'verbose',true);
	% if ERROR;return;end

	% ------------------------------------------------------------------------------
	% Beam capture achieved, switch to orbit mode for tracking. Calculate the orbit
	% response matrix and the dispersion. Assume a beam based alignment procedure
	% reduces the BPM offsets to 30um rms w.r.t. their neighbouring QF/QD magnets.
	%
	SC.INJ.trackMode = 'ORB';

	BPMords = SC.ORD.BPM;
	CMords = SC.ORD.CM;
	% MCO = SCgetModelRM(SC,BPMords,CMords,'trackMode','ORB','useIdealRing',1);
	% MCO = load('./data/idealCORM_AS2v204.mat');
	% save('./data/idealCORM_AS2v204_15BPM_14H19VCM.mat','MCO');
	MCO = par.RM.MCO;
	eta = SCgetModelDispersion(SC,BPMords,SC.ORD.Cavity,'rfStep',5);

	[SC,COexists] = runOrbitCorrection_AS2(SC,MCO,eta,'etaWeight',1E3,'fracSV',par.fracSV);
	results.orbitCorrection = SC;

	% ------------------------------------------------------------------------------
	% tune correction
	%
	% tuneTarget = [70.2512 20.8108];		% empty target tune will default to ideal ring tune
	QMS1 = SCgetOrds(SC.RING,'QMS1');
	QMS2 = SCgetOrds(SC.RING,'QMS2');
	QMS3 = SCgetOrds(SC.RING,'QMS3');
	QMS4 = SCgetOrds(SC.RING,'QMS4');
	qOrds = {QMS1 QMS2 QMS3 QMS4};		% quads used to perform tune correction

    nQMSfams = 4;
    dQMS = 0.001*ones(1,nQMSfams);

	% using SClocolib function, allows for any number of quads to perform tune fit, uses the fminsearch function in matlab to minimise the objective function
	%
	SC = SClocoLib('fitTune',SC,qOrds,'verbose',1,'TolFun',1E-4,'InitStepSize',dQMS,'TolX',1E-5);
	% MCO = SCgetModelRM(SC,BPMords,CMords,'trackMode','ORB','useIdealRing',0);

	eta = SCgetModelDispersion(SC,BPMords,SC.ORD.Cavity,'rfStep',5);
	[SC,COexists] = runOrbitCorrection_AS2(SC,MCO,eta,'etaWeight',1E3,'fracSV',par.fracSV);

	results.tuneCorrection = SC;

   	runTime = tic;
	[DA.area, DA.RMAXs, DA.thetas] = SCdynamicAperture(SC.RING,0,'plot',0,'verbose',1,'bounds',[0,7e-3],'nturns',par.nTurns,'useOrbit6',1,'thetas',linspace(0,pi,16));
   	fprintf('DA area: %.4f mm^2\n',1E6*DA.area);
   	fprintf('Total run time: %.1d sec\n',toc(runTime));
	results.preLocoDA = DA;
	% ------------------------------------------------------------------------------
	% Run LOCO
	%
	fatalError = 1;		% default to stop sextupole ramp and DA search before loco is complete

    [SC,locoResults,fatalError] = runLocoAS2v2(SC,'SVmethod',1E-12,'idealCORM',MCO,'fracSV',par.fracSV,'iterations',par.locoIterations);
    % [SC,locoResults,fatalError] = runLocoAS2inj(SC,'SVmethod',1E-12,'idealCORM',MCO,'fracSV',par.fracSV,'iterations',par.locoIterations);
	results.postLoco = SC;
	results.locoResults = locoResults;
	if fatalError
		fprintf('FatalError!\n');
	end

	% ------------------------------------------------------------------------------
	% post LOCO restore ring to designed parameters
	%
	if ~fatalError
		% ------------------------------------------------------------------------------
		% ramp sextupoles
		%
		% MCO = SCgetModelRM(SC,BPMords,CMords,'trackMode','ORB','useIdealRing',0);
		% eta = SCgetModelDispersion(SC,BPMords,SC.ORD.Cavity,'rfStep',5);

		% nSteps=10;
		% for S = linspace(1/nSteps,1,nSteps)
		% 	SC = SCsetMags2SetPoints(SC,sextOrds,2,3,S,'method','rel');

		% 	MinvCO = SCgetPinv([MCO 1E3*eta],'alpha',0);
		% 	[CUR,ERROR] = SCfeedbackRun(SC,MinvCO,'target',0,'maxsteps',1000,'scaleDisp',1E3,'verbose',1,'eps',1e-6);
		% 	if ~ERROR; SC=CUR; end
		% end
		
		% ------------------------------------------------------------------------------
		% find DA
		%
    	runTime = tic;
		[DA.area, DA.RMAXs, DA.thetas] = SCdynamicAperture(SC.RING,0,'plot',0,'verbose',1,'bounds',[0,7e-3],'nturns',par.nTurns,'useOrbit6',1,'thetas',linspace(0,pi,16));
    	fprintf('DA area: %.4f mm^2\n',1E6*DA.area);
    	fprintf('Total run time: %.1d sec\n',toc(runTime));
		results.postLocoDA = DA;

	end

	%}

end