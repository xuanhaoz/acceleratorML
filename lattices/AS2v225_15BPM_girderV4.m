function ring = AS2v225_15BPM_girderV4
    % combined function dipole in UC, optimised to reduce CF strength, added dipole edge angles for rectangular type
    % created 2024/10/02
    % modified from AS2v224_15BPM_girderV4
    % moved working point to from 62.89/19.61 to 62.89/19.69
    %
    disp('Loading AS2 version 225 lattice with 15 BPMs girder version 4');

    bdsl = 0.2086;
    bdst = 0.29530501737265;   % variation to dispersion suppressor bending angle
    rb = -0.06;     % reverse bending angle
    cdb = 0.2;
    cdbds = 0.2;
    bdstratio = 1.60230170633766;

    CF   = atsbend('CF', 'Length',0.15, 'BendingAngle',deg2rad(rb),...
                    'K',10.109326,  'NumIntSteps',100,'MaxOrder',1);
    CD   = atsbend('CD', 'Length',0.15, 'BendingAngle',deg2rad(cdb-rb),...
                    'K',-3.742584, 'NumIntSteps',100,'MaxOrder',1);
    B1   = atrbend('B1', 'Length',0.55, 'BendingAngle',deg2rad(3.0-2*cdb),...
                    'K',-2.0,      'NumIntSteps',100,'MaxOrder',1);
    CFDS = atsbend('CFDS', 'Length',0.15, 'BendingAngle',deg2rad(rb-bdst/bdstratio),...
                    'K',9.343517,  'NumIntSteps',100,'MaxOrder',1);
    CDDS = atsbend('CDDS', 'Length',0.15, 'BendingAngle',deg2rad(cdbds-rb-bdst*(bdstratio-1)/bdstratio),...
                    'K',-3.413298, 'NumIntSteps',100,'MaxOrder',1);
    BDS  = atrbend('BDS','Length',0.275+bdsl,'BendingAngle',deg2rad(1.5-cdbds+bdst),...
                    'K',-2.138737, 'NumIntSteps',100,'MaxOrder',1);

    % don't need to explicitly define if using atrbend
    %
    % B1.('EntranceAngle')  =   B1.('BendingAngle')/2;
    % B1.('ExitAngle')      =  B1.('BendingAngle')/2;
    % BDS.('EntranceAngle') =  BDS.('BendingAngle')/2;
    % BDS.('ExitAngle')     = BDS.('BendingAngle')/2;

    % 'PassMethod','BndMPoleSymplectic4E2Pass',...
    CF.('PassMethod')    = 'BndMPoleSymplectic4Pass';
    CD.('PassMethod')    = 'BndMPoleSymplectic4Pass';
    B1.('PassMethod')    = 'BndMPoleSymplectic4Pass';
    CFDS.('PassMethod')  = 'BndMPoleSymplectic4Pass';
    CDDS.('PassMethod')  = 'BndMPoleSymplectic4Pass';
    BDS.('PassMethod')   = 'BndMPoleSymplectic4Pass';

    CF.('EntranceAngle')    = 0;
    CD.('EntranceAngle')    = 0;
    CFDS.('EntranceAngle')  = 0;
    CDDS.('EntranceAngle')  = 0;

    CF.('ExitAngle')    = 0;
    CD.('ExitAngle')    = 0;
    CFDS.('ExitAngle')  = 0;
    CDDS.('ExitAngle')  = 0;

    QMS1   = atquadrupole('QMS1','Length',0.15, 'K',  10.004469,'NumIntSteps',10,'MaxOrder',1);
    QMS2   = atquadrupole('QMS2','Length',0.15, 'K',   0.350582,'NumIntSteps',10,'MaxOrder',1);
    QMS3   = atquadrupole('QMS3','Length',0.15, 'K', -10.500000,'NumIntSteps',10,'MaxOrder',1);
    QMS4   = atquadrupole('QMS4','Length',0.15, 'K',   8.338813,'NumIntSteps',10,'MaxOrder',1);

    % 'PassMethod','StrMPoleSymplectic4RadPass',...
    QMS1.('PassMethod') = 'StrMPoleSymplectic4Pass';
    QMS2.('PassMethod') = 'StrMPoleSymplectic4Pass';
    QMS3.('PassMethod') = 'StrMPoleSymplectic4Pass';
    QMS4.('PassMethod') = 'StrMPoleSymplectic4Pass';

    % another working point??
    %
    % SF1 = atsextupole('SF1','Length',0.100,'S', 638.813956,'NumIntSteps',10,'MaxOrder',2);
    % SF2 = atsextupole('SF2','Length',0.100,'S', 591.509199,'NumIntSteps',10,'MaxOrder',2);
    % SF3 = atsextupole('SF3','Length',0.100,'S', 690.213712,'NumIntSteps',10,'MaxOrder',2);
    % SD1 = atsextupole('SD1','Length',0.100,'S',-635.647919,'NumIntSteps',10,'MaxOrder',2);
    % SD2 = atsextupole('SD2','Length',0.100,'S',-658.015477,'NumIntSteps',10,'MaxOrder',2);
    % SD3 = atsextupole('SD3','Length',0.100,'S',-659.166076,'NumIntSteps',10,'MaxOrder',2);
    % SD4 = atsextupole('SD4','Length',0.100,'S',-634.327202,'NumIntSteps',10,'MaxOrder',2);
    % SD5 = atsextupole('SD5','Length',0.100,'S',-603.436885,'NumIntSteps',10,'MaxOrder',2);

    SF1 = atsextupole('SF1','Length',0.100,'S', 631.132191,'NumIntSteps',10,'MaxOrder',2);
    SF2 = atsextupole('SF2','Length',0.100,'S', 634.082107,'NumIntSteps',10,'MaxOrder',2);
    SF3 = atsextupole('SF3','Length',0.100,'S', 612.414064,'NumIntSteps',10,'MaxOrder',2);
    SD1 = atsextupole('SD1','Length',0.100,'S',-595.250356,'NumIntSteps',10,'MaxOrder',2);
    SD2 = atsextupole('SD2','Length',0.100,'S',-619.417592,'NumIntSteps',10,'MaxOrder',2);
    SD3 = atsextupole('SD3','Length',0.100,'S',-605.519986,'NumIntSteps',10,'MaxOrder',2);
    SD4 = atsextupole('SD4','Length',0.100,'S',-627.099798,'NumIntSteps',10,'MaxOrder',2);
    SD5 = atsextupole('SD5','Length',0.100,'S',-536.558178,'NumIntSteps',10,'MaxOrder',2);

    SF1.('PolynomB')(3)  = SF1.('S');
    SF2.('PolynomB')(3)  = SF2.('S');
    SF3.('PolynomB')(3)  = SF3.('S');
    SD1.('PolynomB')(3)  = SD1.('S');
    SD2.('PolynomB')(3)  = SD2.('S');
    SD3.('PolynomB')(3)  = SD3.('S');
    SD4.('PolynomB')(3)  = SD4.('S');
    SD5.('PolynomB')(3)  = SD5.('S');

    SF1.('PassMethod')  = 'StrMPoleSymplectic4Pass';
    SF2.('PassMethod')  = 'StrMPoleSymplectic4Pass';
    SF3.('PassMethod')  = 'StrMPoleSymplectic4Pass';
    SD1.('PassMethod')  = 'StrMPoleSymplectic4Pass';
    SD2.('PassMethod')  = 'StrMPoleSymplectic4Pass';
    SD3.('PassMethod')  = 'StrMPoleSymplectic4Pass';
    SD4.('PassMethod')  = 'StrMPoleSymplectic4Pass';
    SD5.('PassMethod')  = 'StrMPoleSymplectic4Pass';

    % another working point??
    %
    % OC1 = atmultipole('OC1','length',0.05,'polynoma',[0 0 0 0],'polynomb',[0 0 0   555.08055215],'passmethod','strmpolesymplectic4pass');
    % OC2 = atmultipole('OC2','length',0.05,'polynoma',[0 0 0 0],'polynomb',[0 0 0  5332.59367287],'passmethod','strmpolesymplectic4pass');
    % OC3 = atmultipole('OC3','length',0.05,'polynoma',[0 0 0 0],'polynomb',[0 0 0 -7131.66997893],'passmethod','strmpolesymplectic4pass');
    % OC4 = atmultipole('OC4','length',0.05,'polynoma',[0 0 0 0],'polynomb',[0 0 0  3850.78562137],'passmethod','strmpolesymplectic4pass');

    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % octupole strength conversion from OPA need to change from integrated strength to normalised
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    OC1 = atmultipole('OC1','Length',0.10,'PolynomA',[0 0 0 0],'PolynomB',[0 0 0   624.39444963/0.10],'PassMethod','StrMPoleSymplectic4Pass');
    OC2 = atmultipole('OC2','Length',0.10,'PolynomA',[0 0 0 0],'PolynomB',[0 0 0  4570.06894819/0.10],'PassMethod','StrMPoleSymplectic4Pass');
    OC3 = atmultipole('OC3','Length',0.10,'PolynomA',[0 0 0 0],'PolynomB',[0 0 0 -7398.84511423/0.10],'PassMethod','StrMPoleSymplectic4Pass');
    OC4 = atmultipole('OC4','Length',0.10,'PolynomA',[0 0 0 0],'PolynomB',[0 0 0  4587.57075019/0.10],'PassMethod','StrMPoleSymplectic4Pass');

    D25    = atdrift('Drift',0.0250,'PassMethod','DriftPass');
    D37    = atdrift('Drift',0.0375,'PassMethod','DriftPass');
    D75    = atdrift('Drift',0.0750,'PassMethod','DriftPass');
    D50    = atdrift('Drift',0.0500,'PassMethod','DriftPass');
    D150   = atdrift('Drift',0.1500,'PassMethod','DriftPass');
    D100   = atdrift('Drift',0.1000,'PassMethod','DriftPass');
    D250   = atdrift('Drift',0.2500,'PassMethod','DriftPass');
    DDS3   = atdrift('Drift',0.2500-bdsl,'PassMethod','DriftPass');
    D2500  = atdrift('Drift',2.5000,'PassMethod','DriftPass');

    BPM = @(name) atmonitor(name, 'IdentityPass');
    MARK = @(name) atmarker(name, 'IdentityPass');

    COR = @(name) atcorrector(name);

    RFC  = atrfcavity('RFCav',...
        'Energy',3E9);
        % 'Voltage',2e6,...
        % 'HarmNumber',758,...
        % 'Frequency',500e6,...
    	% 'PassMethod','RFCavityPass');


    % Unit Cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    UC1 = {D50 SF3 D50 CF D75 CD D37 BPM('BPM5') D37 SD4 D50 B1 D50 SD3 D37 BPM('BPM6')  D37 CD D37 COR('COR4') D37 CF};
    UC2 = {D50 SF2 D50 CF D75 CD D37 MARK('GirderEnd') MARK('GirderStart') BPM('BPM7')  D37 SD2 D50 B1 D50 SD1 D37 BPM('BPM8')  D37 CD D37 COR('COR5') D37 CF};
    UC3 = {D50 SF1 D50 CF D75 CD D37 BPM('BPM9')  D37 SD1 D50 B1 D50 MARK('GirderEnd') MARK('GirderStart') SD2 D37 BPM('BPM10') D37 CD D37 COR('COR6') D37 CF};
    UC4 = {D50 SF2 D50 CF D75 CD D37 BPM('BPM11') D37 SD3 D50 B1 D50 SD4 D37 BPM('BPM12') D37 CD D37 COR('COR7') D37 CF};
    UC5 = {D50 SF3 D50};

    arc = [UC1 UC2 UC3 UC4 UC5];

    % Dispersion Suppressor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    dsPartA = {CFDS D37 COR('COR3') D37 CDDS D37 BPM('BPM4') D37 SD5 D50 MARK('GirderStart') MARK('GirderEnd') DDS3 BPM('BPM3') D50 BDS};
    dsPartB = {CFDS D75 CDDS D75 MARK('GirderEnd') MARK('GirderStart') SD5 D50 DDS3 BPM('BPM13') D50 BDS};

    dsA = [flip(dsPartA)];
    dsB = [dsPartB];

    % Matching Straight %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    msPartA = {D100 QMS1 D75 OC1 D75 QMS2 D75 OC2 D37 COR('COR2') BPM('BPM2')  D37 QMS3 D75 OC3 D75 QMS4 D75 OC4 D37 COR('COR1') BPM('BPM1')  D37 MARK('GirderStart') D2500};
    msPartB = {D100 QMS1 D75 OC1 D75 QMS2 D75 OC2 D37 COR('COR8') BPM('BPM14') D37 QMS3 D75 OC3 D75 QMS4 D75 OC4 D37 COR('COR9') BPM('BPM15') D37 MARK('GirderEnd')   D2500};

    msA = [flip(msPartA)];
    msB = [msPartB];

    % Ring %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %

    sector = [msA dsA arc dsB msB];
    sector = sector';

    RP = atringparam('AS2',3.0e9);
    ring = [{RP}; {RFC}; repmat(sector,24,1)];
    % --------------------------------------
    % additional cavity at midway 29/04/2024
    %
    % ring = [{RP}; {RFC}; repmat(sector,12,1); {RFC}; repmat(sector,12,1)];
    % --------------------------------------
    ring = atsetcavity(ring,3e6,1,750);

    % ring = atsetcavity(ring,'frequency','nominal',...
    %         'harmnumber',758,'voltage',2e6);

    % ring = atsetringproperties(ring);
    % ring = atsetrfcavity(ring,2e6,1,758,0);

end