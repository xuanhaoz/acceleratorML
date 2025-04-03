function out = scanSVDcorrection(varargin)
    SC = varargin{1};
    version = getoption(varargin,'version','assr');

    switch version
        case 'assr'
            ring0 = assr4_splitbends;
            MCO = load('assr_idealCORM.mat').MCO;
        otherwise
            error('lattice version not supported');
    end

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

    SVrange = 0.4:0.02:0.9;
    parfor i = 1:length(SVrange)
        fracSV = SVrange(i);
        fprintf('Running SV frac %.2f (job %02d/%02d)\n',fracSV,i,length(SVrange));
        [SC1,COexists] = runOrbitCorrection_AS2(SC,MCO,eta,'etaWeight',1,'fracSV',fracSV,'buildTargetOrbit',1);
        res{i} = SC1;

        bpmVal = plotBPMreading(SC1.RING,'plot',0);
        BPMrms(i) = bpmVal.BPMrms;
    end
    out.res = res;
    out.BPMrms = BPMrms;

    % if ~COexists
    %     newSeed.postCorrection = 'Correction failed';
    % else
    %     newSeed.postCorrection = SC.RING;
    %     newSeed.SCpostCorrection = SC;
    %     figure(3171)
    %     clf
    %     plotBPMreading(newSeed.postCorrection,'verbose',1);
    % end
end

