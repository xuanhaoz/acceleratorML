function [SC,COexists] = runOrbitCorrection_AS2(SC,MCO,eta,varargin)

    p = inputParser;
    addOptional(p,'alphaVec',[]);
    addOptional(p,'rampSext',0);
    addOptional(p,'checkTransmissionMode',0);   % use tracking to determine improvement in correction, use when closed orbit is bad
    addOptional(p,'etaWeight',1E2);
    addOptional(p,'SV',0);
    addOptional(p,'fracSV',1);
    addOptional(p,'verbose',0);
    addOptional(p,'buildTargetOrbit',0);
    parse(p,varargin{:});
    par = p.Results;

    BPMords = SC.ORD.BPM;
    CMords = SC.ORD.CM;

    % MCO = SCgetModelRM(SC,BPMords,CMords,'trackMode','ORB','useIdealRing',1);
    % eta = SCgetModelDispersion(SC,BPMords,SC.ORD.Cavity,'rfStep',5);

	% RMstruct.RM = RM([find(ismember(BPMords,BPMords)) length(BPMords)+find(ismember(BPMords,BPMords))],[find(ismember(CMords{1},CMords{1})) length(CMords{1})+find(ismember(CMords{2},CMords{2}))]);
    % RMstruct.eta = SCgetModelDispersion(SC,BPMords,SC.ORD.Cavity,'rfStep',5);
    % RMstruct.scaleDisp = 1E6;

    % MCO = [RMstruct.RM RMstruct.scaleDisp*RMstruct.eta];


    if isempty(par.alphaVec)
        par.alphaVec = [100 50 20 10 5:-1:1 0.9:-0.1:0];
    end

    if par.rampSext == 1
        % turning off sextupoles
        %
        sextOrds = SCgetOrds(SC.RING,'SF|SD|SX1|SX2|SX3|SH1|SH2|SH4|SH5');
        SC = SCsetMags2SetPoints(SC,sextOrds,2,3,0,'method','abs');
    end

    alphaVec = par.alphaVec;
    % flag to check if the first try alpha value yields stable machine, if not then fall back to TBT threading to orbit
    %
    for alpha = alphaVec
        if find(alphaVec==alpha) == 1
            firstPass = 1;
        else
            firstPass = 0;
        end

        MinvCO = SCgetPinv([MCO par.etaWeight*eta],'alpha',alpha,'N',par.SV,'fracSV',par.fracSV);

        if par.buildTargetOrbit
            R0 = buildTargetOrbit(SC);
        else
            R0 = zeros(size(MCO,1),1);
        end

        % check if closed orbit can be found, otherwise switch to pseudo orbit mode
        %
        [SC,COexists] = checkOrbit(SC,par.verbose);
        if ~COexists
            return
        end

        [CUR,ERROR] = SCfeedbackRun(SC,MinvCO,'target',0,'maxsteps',1000,'scaleDisp',par.etaWeight,'verbose',par.verbose,'eps',1e-6,'R0',R0);

        if ERROR & firstPass
            COexists = 0;
            return
        end

        if ~ERROR
            if par.checkTransmissionMode
                maxTurns = SCgetBeamTransmission(CUR,'nParticles',100,'nTurns',100,'verbose',par.verbose);
                if maxTurns < 50
                    if par.verbose
                        fprintf('Transmission below 50 turns, abort orbit correction.\n')
                    end
                    alpha = alphaVec(max([find(alphaVec == alpha)-1,1])); break
                end
            end

            % calculate initial and final BPM reading
            %
			B0rms = sqrt(mean(SCgetBPMreading(SC ,'BPMords',BPMords).^2,2));
			Brms  = sqrt(mean(SCgetBPMreading(CUR,'BPMords',BPMords).^2,2));
            % check if orbit improved
            %
            if mean(B0rms) < mean(Brms)
                if par.verbose
                    fprintf('No further improvements with alpha = %d\n',alpha);
                end
                alpha = alphaVec(max([find(alphaVec == alpha)-1,1])); break
            else
                SC = CUR;
                if par.verbose
				    fprintf('CM improvement with alpha = %d:\n hor: %.3fum -> %.3fum\n ver: %.3fum -> %.3fum\n',alpha,1E6*B0rms(1),1E6*Brms(1),1E6*B0rms(2),1E6*Brms(2));
                end
            end
        else
            alpha = alphaVec(max([find(alphaVec == alpha)-1,1])); break
        end
    end
    if par.verbose
	    fprintf('Final closed orbit deviation for all BPMs (hor//ver): %.3fum // %.3fum   with alpha = %d.\n',1E6*sqrt(mean(SCgetBPMreading(SC).^2,2)),alpha);
    end

    if par.rampSext
	    nSteps=10;
	    for S = linspace(1/nSteps,1,nSteps)
		    SC = SCsetMags2SetPoints(SC,sextOrds,2,3,S,'method','rel');

            MinvCO = SCgetPinv([MCO 1E3*eta],'alpha',100);
            [CUR,ERROR] = SCfeedbackRun(SC,MinvCO,'target',0,'maxsteps',1000,'scaleDisp',1E3,'verbose',par.verbose,'eps',1e-6);
		    if ~ERROR; SC=CUR; end
		    % [CUR,ERROR] = SCfeedbackBalance(SC,MinvCO,'maxsteps',1000,'eps',1e-6,'verbose',1);
		    % if ~ERROR; SC=CUR; end
	    end
    end


end

% Check if closed orbit can be found
function [SC,COexists] = checkOrbit(SC,verbose)
	if  ~any(isnan(findorbit6(SC.RING))) 
        COexists = 1;
		if ~strcmp(SC.INJ.trackMode,'ORB')
			SC.INJ.trackMode = 'ORB';
            if verbose
			    fprintf('Switch to orbit mode\n')
            end
		end
	else
		% SC.INJ.trackMode  = 'pORB';
		% SC.INJ.nTurns     = 100;
		% SC.INJ.nParticles = 50;
		% fprintf('Switch to pseudo orbit mode\n')
        if verbose
            fprintf('No closed orbit found\n');
        end
        COexists = 0;
	end
end

function out = buildTargetOrbit(varargin)
    SC = varargin{1};

    ring0 = SC.IDEALRING;
    bpm = find(atgetcells(ring0,'FamName','BPM.*'));
    [ed,~] = atlinopt(ring0,0,1:length(ring0)+1);
    CO = [ed.ClosedOrbit];

    R0 = [CO(1,bpm) CO(3,bpm)];
    out = R0;

end