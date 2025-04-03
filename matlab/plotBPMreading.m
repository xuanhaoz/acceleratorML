function out = plotBPMreading(varargin)
    ring = varargin{1};
    verbose = getoption(varargin,'verbose',0);
    makeplot = getoption(varargin,'plot',1);

    bpm = find(atgetcells(ring,'Class','Monitor'));

    % get actual trajectory
    %
    [ed,~] = atlinopt(ring,0,1:length(ring)+1);

    spos = [ed.SPos];
    CO = [ed.ClosedOrbit];

    if verbose
        fprintf('BPM rms reading: %.4f [micron]\n',1e6*rms(CO(1,bpm)));
    end

    out.CO = CO(1,:);
    out.spos = spos;
    out.BPMrms = rms(CO(1,bpm));

    if makeplot
        l1 = plot(spos,CO(1,:),'black','LineWidth',1.5);
        hold on
        plot(spos(bpm),CO(1,bpm),'ro','LineWidth',1.5);
        legend([l1],sprintf('rms: %.2f [micron]',1e6*rms(CO(1,:))));

        xlim([0 spos(end)]);
        grid on
        xlabel('s [m]');
        ylabel('x [m]');
    end
end
