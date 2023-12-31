function plt_templates_bp(greenf,stas,lowlet,hiwlet,sps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plt_templates_bp(greenf,lowlet,hiwlet,sps)
%
% this function simply plot the bandpassed filtered templates (or green's 
% functions), either optimal component or orthogonal comp.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/09/12
% Last modified date:   2022/09/12 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
lwlet = size(greenf,1);
if size(greenf,2) == 3
  plot((1:lwlet)/sps,greenf(:,1),'r');
  plot((1:lwlet)/sps,greenf(:,2),'b');
  plot((1:lwlet)/sps,greenf(:,3),'k');
else
  for ista = 1: size(greenf,2)
    plot((1:lwlet)/sps,greenf(:,ista));
  end
end  
text(0.95,0.9,sprintf('%.1f-%.1f Hz opt.',lowlet,hiwlet),'Units','normalized','HorizontalAlignment',...
  'right');
mx=max(max(abs(greenf(:,:))));
%%%running CC using a window length of 'cclen'
mwlen=sps/2;
[ircc,rcc12] = RunningCC(greenf(:,1), greenf(:,2), mwlen);
[~,rcc13] = RunningCC(greenf(:,1), greenf(:,3), mwlen);
[~,rcc23] = RunningCC(greenf(:,2), greenf(:,3), mwlen);
rcc = (rcc12+rcc13+rcc23)/3;
%alln(alln<0)=-10^4*yma; %just so they don't plot.
% plot(samples(cclen/2+1:greenlen-cclen/2),mx*alln,'co','markersize',2); % scale with data amp.
plot(ircc/sps,mx*rcc,'co','markersize',2); % scale with data amp.
legend([stas; 'RCC  '],'Location','southeast');
greenlen = size(greenf,1);
xlim([0 greenlen]/sps)
ylim([-mx mx])
box on
xlabel('Time (s)');

