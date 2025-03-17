for i = 1:length(sf)
    ele = ring{sf(i)};
    HCM(i) = ele.PolynomB(1)*ele.Length;
end

for i = 1:length(sd)
    ele = ring{sd(i)};
    VCM(i) = ele.PolynomA(1)*ele.Length;
end

figure(1)
clf
subplot(2,1,1)
plot(1:length(HCM),HCM,'ro','LineWidth',1.5);
grid on
xlabel('HCM number');ylabel('HCM kick [rad]')

subplot(2,1,2)
plot(1:length(VCM),VCM,'ro','LineWidth',1.5);
grid on
xlabel('VCM number');ylabel('VCM kick [rad]')

% res = {};
% varList = [1];
% for i = 1:length(varList)
%     ring = createSeeds(1,'version','225');
%     res{i} = ring;
% end
% out = res;

% ring = assr4_splitbends;

% diff = [];
% % for i = 1:length(varList)
% for i = 1:1
%     clf
%     subplot(5,1,1:2)
%     CO0 = plotBPMreading(ring);
%     legend('Ideal ring');

%     subplot(5,1,3:4)
%     CO1 = plotBPMreading(res{i});
%     legend(sprintf('fracSV = %.2f',i/100+0.6));

%     diff(i,:) = CO1.CO - CO0.CO;
%     subplot(5,1,5)
%     plot(CO1.spos,diff(i,:),'LineWidth',1.5);
%     grid on
%     xlim([0 CO1.spos(end)]);
%     % pause(2)
% end
