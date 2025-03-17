function out = plotBPMreading(varargin)
    ring = varargin{1};
    verbose = getoption(varargin,'verbose',0);

    bpm = find(atgetcells(ring,'Class','Monitor'));

    % get actual trajectory
    %
    [ed,~] = atlinopt(ring,0,1:length(ring)+1);

    spos = [ed.SPos];
    CO = [ed.ClosedOrbit];

    if verbose
        fprintf('BPM rms reading: %.4f [micron]\n',1e6*rms(CO(1,bpm)));
    end

    plot(spos,CO(1,:),'black','LineWidth',1.5);
    hold on
    plot(spos(bpm),CO(1,bpm),'ro','LineWidth',1.5);

    out.CO = CO(1,:);
    out.spos = spos;

    xlim([0 spos(end)]);
    grid on
    xlabel('s [m]');
    ylabel('x [m]');
end
