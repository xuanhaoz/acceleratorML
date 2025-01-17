function QuadOrds = getBPM2BBAPairing_AS2(SC,BPMords,varargin)
	% updated function to find BPM to sextupole/quadrupole pairing
	%
    p = inputParser;
    addOptional(p,'nBPM',15);
    addOptional(p,'version',1);
    parse(p,varargin{:});
    par = p.Results;
		
	if par.nBPM == 19 
		BPM( 1).QuadName = 'QMS4'; 
		BPM( 1).updown   = 'u'; 

		BPM( 2).QuadName = 'QMS3'; 
		BPM( 2).updown   = 'd'; 

		BPM( 3).QuadName = 'SD5'; 
		BPM( 3).updown   = 'u'; 

		BPM( 4).QuadName = 'SF3'; 
		BPM( 4).updown   = 'u'; 

		BPM( 5).QuadName = 'SD4'; 
		BPM( 5).updown   = 'u'; 

		BPM( 6).QuadName = 'SF2'; 
		BPM( 6).updown   = 'u'; 

		BPM( 7).QuadName = 'SD2'; 
		BPM( 7).updown   = 'u'; 

		BPM( 8).QuadName = 'SF1'; 
		BPM( 8).updown   = 'u'; 

		BPM( 9).QuadName = 'SD1'; 
		BPM( 9).updown   = 'u'; 

		BPM(10).QuadName = 'SF2'; 
		BPM(10).updown   = 'u'; 

		BPM(11).QuadName = 'SD3'; 
		BPM(11).updown   = 'u'; 

		BPM(12).QuadName = 'SF3'; 
		BPM(12).updown   = 'u'; 

		BPM(13).QuadName = 'SD5'; 
		BPM(13).updown   = 'd'; 

		BPM(14).QuadName = 'QMS3'; 
		BPM(14).updown   = 'u'; 

		BPM(15).QuadName = 'QMS4'; 
		BPM(15).updown   = 'd'; 

		BPM(16).QuadName = 'Q1INJ'; 
		BPM(16).updown   = 'u'; 

		BPM(17).QuadName = 'Q3INJ'; 
		BPM(17).updown   = 'd'; 

		BPM(18).QuadName = 'Q3INJ'; 
		BPM(18).updown   = 'u'; 

		BPM(19).QuadName = 'Q1INJ'; 
		BPM(19).updown   = 'd'; 

	elseif par.nBPM == 15 
		BPM( 1).QuadName = 'QMS4'; 
		BPM( 1).updown   = 'u'; 

		BPM( 2).QuadName = 'QMS3'; 
		BPM( 2).updown   = 'd'; 

		BPM( 3).QuadName = 'SD5'; 
		BPM( 3).updown   = 'u'; 

		BPM( 4).QuadName = 'SF3'; 
		BPM( 4).updown   = 'u'; 

		BPM( 5).QuadName = 'SD4'; 
		BPM( 5).updown   = 'u'; 

		BPM( 6).QuadName = 'SF2'; 
		BPM( 6).updown   = 'u'; 

		BPM( 7).QuadName = 'SD2'; 
		BPM( 7).updown   = 'u'; 

		BPM( 8).QuadName = 'SF1'; 
		BPM( 8).updown   = 'u'; 

		BPM( 9).QuadName = 'SD1'; 
		BPM( 9).updown   = 'u'; 

		BPM(10).QuadName = 'SF2'; 
		BPM(10).updown   = 'u'; 

		BPM(11).QuadName = 'SD3'; 
		BPM(11).updown   = 'u'; 

		BPM(12).QuadName = 'SF3'; 
		BPM(12).updown   = 'u'; 

		BPM(13).QuadName = 'SD5'; 
		BPM(13).updown   = 'd'; 

		BPM(14).QuadName = 'QMS3'; 
		BPM(14).updown   = 'u'; 

		BPM(15).QuadName = 'QMS4'; 
		BPM(15).updown   = 'd'; 

	elseif par.nBPM == 14 
		BPM( 1).QuadName = 'QMS4'; 
		BPM( 1).updown   = 'u'; 

		BPM( 2).QuadName = 'QMS3'; 
		BPM( 2).updown   = 'd'; 

		BPM( 3).QuadName = 'SD5'; 
		BPM( 3).updown   = 'u'; 

		BPM( 4).QuadName = 'SD4'; 
		BPM( 4).updown   = 'u'; 

		BPM( 5).QuadName = 'SD3'; 
		BPM( 5).updown   = 'd'; 

		BPM( 6).QuadName = 'SD2'; 
		BPM( 6).updown   = 'u'; 

		BPM( 7).QuadName = 'SD1'; 
		BPM( 7).updown   = 'd'; 

		BPM( 8).QuadName = 'SD1'; 
		BPM( 8).updown   = 'u'; 

		BPM( 9).QuadName = 'SD2'; 
		BPM( 9).updown   = 'd'; 

		BPM(10).QuadName = 'SD3'; 
		BPM(10).updown   = 'u'; 

		BPM(11).QuadName = 'SD4'; 
		BPM(11).updown   = 'd'; 

		BPM(12).QuadName = 'SD5'; 
		BPM(12).updown   = 'd'; 

		BPM(13).QuadName = 'QMS3'; 
		BPM(13).updown   = 'u'; 

		BPM(14).QuadName = 'QMS4'; 
		BPM(14).updown   = 'd'; 

	else
		error('BPM to Quad pairing not defined!')
	end
	
	for nDim=1:size(BPMords,1)
		i = 1;
		for ord=BPMords(nDim,:)
			BPMind  = str2num(SC.RING{ord}.FamName(4:end));
			tmpQuad = SCgetOrds(SC.RING,BPM(BPMind).QuadName);
			
			if isempty(tmpQuad)
				error('No quadrupole name could be found.')
			end
			
			switch BPM(BPMind).updown
				case 'u' % quad is upstream of BPM
					QuadOrds(nDim,i) = tmpQuad(find(tmpQuad>ord,1));
				case 'd' % quad is downstream of BPM
					QuadOrds(nDim,i) = tmpQuad(find(tmpQuad<ord,1,'last'));
			end
			i = i+1;
		end
	end
end
