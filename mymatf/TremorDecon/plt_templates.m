function plt_templates(green,greenf,greenort,greenfort,lowlet,hiwlet,sps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plt_templates(green,greenf,greenort,greenfort,lowlet,hiwlet)
%
% this function simply plot the templates (or green's functions), in raw
% waveform or bandpass filtered; both optimal component and orthogonal
% comp. Frequently used in 'deconvbursts?????', etc.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/09/12
% Last modified date:   2022/09/12 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(2,2,1)
hold on
plot(green(:,1),'r')
plot(green(:,2),'b')
plot(green(:,3),'k')
text(0.95,0.9,'Raw opt.','Units','normalized','HorizontalAlignment',...
  'right');
mx=max(max(abs(green(:,:))));
greenlen = size(green,1);
xlim([0 greenlen])
ylim([-mx mx])
box on

subplot(2,2,2)
hold on
plot(greenf(:,1),'r')
plot(greenf(:,2),'b')
plot(greenf(:,3),'k')
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
plot(ircc,mx*rcc,'co','markersize',2); % scale with data amp.
xlim([0 greenlen])
ylim([-mx mx])
box on
xc23=xcorr(greenf(:,2),greenf(:,3),10,'coeff');
xc13=xcorr(greenf(:,1),greenf(:,3),10,'coeff');
xc12=xcorr(greenf(:,1),greenf(:,2),10,'coeff');
[ccmax23,imax23]=max(xc23);
[ccmax13,imax13]=max(xc13);
[ccmax12,imax12]=max(xc12);

subplot(2,2,3)
hold on
plot(greenort(:,1),'r')
plot(greenort(:,2),'b')
plot(greenort(:,3),'k')
text(0.95,0.9,'Raw ort.','Units','normalized','HorizontalAlignment',...
  'right');
mx=max(max(abs(green(:,:))));
xlim([0 greenlen])
ylim([-mx mx])
box on

subplot(2,2,4)
hold on
plot(greenfort(:,1),'r')
plot(greenfort(:,2),'b')
plot(greenfort(:,3),'k')
text(0.95,0.9,sprintf('%.1f-%.1f Hz ort.',lowlet,hiwlet),'Units','normalized','HorizontalAlignment',...
  'right');
mx=max(max(abs(greenf(:,:))));

%%%running CC using a window length of 'cclen'
mwlen=sps/2;
[ircc,rcc12] = RunningCC(greenfort(:,1), greenfort(:,2), mwlen);
[~,rcc13] = RunningCC(greenfort(:,1), greenfort(:,3), mwlen);
[~,rcc23] = RunningCC(greenfort(:,2), greenfort(:,3), mwlen);
rcc = (rcc12+rcc13+rcc23)/3;
%alln(alln<0)=-10^4*yma; %just so they don't plot.
% plot(samples(cclen/2+1:greenlen-cclen/2),mx*alln,'co','markersize',2); % scale with data amp.
plot(ircc,mx*rcc,'co','markersize',2); % scale with data amp.
xlim([0 greenlen])
ylim([-mx mx])
box on
xc23=xcorr(greenfort(:,2),greenfort(:,3),10,'coeff');
xc13=xcorr(greenfort(:,1),greenfort(:,3),10,'coeff');
xc12=xcorr(greenfort(:,1),greenfort(:,2),10,'coeff');
[ccmax23,imax23]=max(xc23);
[ccmax13,imax13]=max(xc13);
[ccmax12,imax12]=max(xc12);




