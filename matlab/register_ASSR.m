% matlab script as part of simulated commissioning toolkit (SC) based program
% for AS storage ring v4 splitbends
% modified based on register_AS2v2.m
% 2025/03/13, F. Zhang
%

function SC = register_ASSR(SC,varargin)
    % ---------------------------------------
    % register lattice in SC
    % NOTE: default error generation is gaussian truncated at 2 sigma
    %
    p = inputParser;
    addOptional(p,'skewQuad',1);
    addOptional(p,'sextCM',1);
    addOptional(p,'quadCM',0);
    addOptional(p,'errorScale',1);
    addOptional(p,'girders',0);
    addOptional(p,'injection',0);
    parse(p,varargin{:});
    par = p.Results;

    disp('Registering lattice in SC')

    ords = findcells(SC.RING,'Frequency');
    SC = SCregisterCAVs(SC,ords,...
    	'FrequencyOffset',0,... % [Hz]
    	'VoltageOffset',0,...   % [V]
    	'TimeLagOffset',0);     % [m]
    
    % ----------------------------------
    errorScale   = par.errorScale;

    errorRoll    = errorScale*50e-6;    % baseline 100e-6
    errorDx      = errorScale*15e-6;       % baseline 30e-6
    errorDy      = errorScale*15e-6;       % baseline 30e-6
    errorFSE     = errorScale*1e-3;       % baseline 1e-3
    errorFSE_CFD = errorScale*0;      % error on combined function quadrupole strength
    errorFSE_B   = errorScale*0;        % error on bending angle

    girderDx     = errorScale*15e-6;
    girderDy     = errorScale*15e-6;
    girderRoll   = errorScale*50e-6;% 200e-6;

    errorDxSX    = errorScale*15e-6;
    errorDySX    = errorScale*15e-6;

    errorDxBPM   = errorScale*0e-6;
    errorDyBPM   = errorScale*0e-6;
    errorRollBPM = errorScale*00e-6;

    % ----------------------------------
    % injection quads
    if par.injection
        ords = SCgetOrds(SC.RING,'Q1INJ|Q2INJ|Q3INJ');
        SC = SCregisterMagnets(SC,ords,...
        	'CalErrorB',[0 errorFSE],...            % relative
        	'MagnetOffset',[errorDx errorDy 0],...  % x, y and z, [m]
        	'MagnetRoll',errorRoll * [1 0 0]);      % az, ax and ay, [rad]
    end

    % ----------------------------------
    % girders
    % girders not implemented in lattice file (2025/03/13)
    %
    if par.girders
        ords = [SCgetOrds(SC.RING,'GirderStart'); SCgetOrds(SC.RING,'GirderEnd')];
        SC = SCregisterSupport(SC,'Girder',ords,...
            'Offset', [girderDx girderDy 0],...
            'Roll', girderRoll * [1 0 0]);
    end

    % ----------------------------------
    % BPMs
    %
    ords = SCgetOrds(SC.RING,'BPM');
    SC = SCregisterBPMs(SC, ords,...
    	'CalError', 0*0.05 * [1 1],... % x and y, relative
    	'Offset',   [errorDxBPM errorDyBPM],... % x and y, [m]
    	'Noise',    0*1e-6 * [1 1],... % x and y, [m]
    	'NoiseCO',  0*1e-7 * [1 1],... % x and y, [m]
    	'Roll',     errorRollBPM);           % az, [rad]
    
    % ----------------------------------
    % CORs
    %
    ords = SCgetOrds(SC.RING,'^COR');
    SC = SCregisterMagnets(SC,ords,...
    	'HCM',Inf,...                   % [rad]
    	'VCM',Inf,...                   % [rad]
        'CalErrorB',0,...               % Relative error in HCM
        'CalErrorA',0,...               % Relative error in VCM
    	'MagnetOffset',0 * [1 1 0],...  % x, y and z, [m]
    	'MagnetRoll',0 * [1 0 0]);      % az, ax and ay, [rad]
    
    % ----------------------------------
    % BENDs
    %
    
    ords = SCgetOrds(SC.RING,'^b_left|^b_centre|^b_right');
    SC = SCregisterMagnets(SC,ords,...
        'CalErrorB', {[0 errorFSE_CFD], 2},...      % rel error in quad strength, tied to bending angle
        'BendingAngle', {errorFSE_B, 2},...
    	'MagnetOffset',[0 0 0],...  % x, y and z, [m]
    	'MagnetRoll',0 * [1 0 0]);      % az, ax and ay, [rad]
    
    % ----------------------------------
    % Quads
    %
    ords = SCgetOrds(SC.RING,'QFA|QFB|QDA');
    SC = SCregisterMagnets(SC,ords,...
    	'CalErrorB',[0 errorFSE],...            % relative
    	'MagnetOffset',[errorDx errorDy 0],...  % x, y and z, [m]
    	'MagnetRoll',errorRoll * [1 0 0]);      % az, ax and ay, [rad]
    
    % ----------------------------------
    % Sextupoles
    %
    if par.sextCM
        ords = SCgetOrds(SC.RING,'^SF');
        SC = SCregisterMagnets(SC,ords,...
        	'HCM',Inf,...                   % [rad]
            'CalErrorB',[0 0 0],...             % relative error in sextupole strength
        	'MagnetOffset',1 * [errorDxSX errorDySX 0],...      % x, y and z, [m]
        	'MagnetRoll',0 * [1 0 0]);          % az, ax and ay, [rad]
    
        ords = SCgetOrds(SC.RING,'^SD');
        if par.skewQuad
            SC = SCregisterMagnets(SC,ords,...
        	    'VCM',Inf,...                   % [rad]
                'SkewQuad',1,...                    % add skew quad component
                'CalErrorB',[0 0 0],...             % relative error in sextupole strength
        	    'MagnetOffset',1 * [errorDxSX errorDySX 0],...      % x, y and z, [m]
            	'MagnetRoll',0 * [1 0 0]);          % az, ax and ay, [rad]
        else
            SC = SCregisterMagnets(SC,ords,...
        	    'VCM',Inf,...                   % [rad]
                'CalErrorB',[0 0 0],...             % relative error in sextupole strength
        	    'MagnetOffset',1 * [errorDxSX errorDySX 0],...      % x, y and z, [m]
            	'MagnetRoll',0 * [1 0 0]);          % az, ax and ay, [rad]
        end
    else
        ords = SCgetOrds(SC.RING,'^SF');
        SC = SCregisterMagnets(SC,ords,...
            'CalErrorB',[0 0 0],...             % relative error in sextupole strength
        	'MagnetOffset',1 * [errorDxSX errorDySX 0],...      % x, y and z, [m]
        	'MagnetRoll',0 * [1 0 0]);          % az, ax and ay, [rad]
    
        ords = SCgetOrds(SC.RING,'^SD');
        if par.skewQuad
            SC = SCregisterMagnets(SC,ords,...
                'SkewQuad',1,...                    % add skew quad component
                'CalErrorB',[0 0 0],...             % relative error in sextupole strength
        	    'MagnetOffset',1 * [errorDxSX errorDySX 0],...      % x, y and z, [m]
            	'MagnetRoll',0 * [1 0 0]);          % az, ax and ay, [rad]
        else
            SC = SCregisterMagnets(SC,ords,...
                'CalErrorB',[0 0 0],...             % relative error in sextupole strength
        	    'MagnetOffset',1 * [errorDxSX errorDySX 0],...      % x, y and z, [m]
            	'MagnetRoll',0 * [1 0 0]);          % az, ax and ay, [rad]
        end
    end


end