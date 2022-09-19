function [f1,f2,clppk,clnpk] = plt_deconpk_sigpk_comp(sigsta,zcrssrc,ppksrc,npksrc,greenf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f] = plt_deconpk_sigpk(sigsta,pkindep,indremove,fpltremv)
%
% Different from 'plt_deconpk_sigpk' which only show the location of deconvolved
% positive peaks and waveform peaks to see if they are close enough, thus 
% reasonable, this script inherits most of it, but also shows the amplitude
% of the deconvolved, to see how comparable they are to the nearest waveform
% peaks. This comparison may give us some insights on if it is necessary to bring
% the waveform amplitude ratio information into the grouping algorithm of the
% the deconvolved sources.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/09/08
% Last modified date:   2022/09/08 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1.fig = figure;
f1.fig.Renderer = 'painters';
widin = 12;  % maximum width allowed is 8.5 inches
htin = 6;   % maximum height allowed is 11 inches
[scrsz, resol] = pixelperinch(1);
set(f1.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);
nrow = size(sigsta,2); 
ncol = 1;
for isub = 1:nrow*ncol
  f1.ax(isub) = subplot(nrow,ncol,isub);
end
pltxran = [0.1 0.9]; pltyran = [0.15 0.9];
pltxsep = 0.02; pltysep = 0.05;
%get the locations for each axis
optaxpos(f1,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

color = ['r';'b';'k'];

ym = max(abs(sigsta(:)));
yran=1.2*[-ym ym];

lsig = size(sigsta,1); 

ppkhtwf = cell(size(sigsta,2),1);  %positve peak height of the waveform
ppkwf = cell(size(sigsta,2),1); %positive peak index of the waveform
npkhtwf = cell(size(sigsta,2),1);  %negative peak height of the waveform
npkwf = cell(size(sigsta,2),1); %negative peak index of the waveform
for i = 1: nrow
  ax=f1.ax(i);
  hold(ax,'on');
  plot(ax, sigsta(:,i), '-','Color',color(i,:)); 
  xlim(ax,[0,lsig]); ylim(ax,yran);
  [ppkhtwf{i}, ppkwf{i}] = findpeaks(sigsta(:,i));
  [npkhtwf{i}, npkwf{i}] = findpeaks(-sigsta(:,i));
  npkhtwf{i} = -npkhtwf{i}; %reverse back the sign
  scatter(ax,ppkwf{i},ppkhtwf{i},8,color(i,:)); %plot the negative waveform peaks only
  for j = 1: size(zcrssrc,1)
    plot(ax,[zcrssrc(j,(i-1)*2+1) zcrssrc(j,(i-1)*2+1)], ax.YLim, '--','Color',[.7 .7 .7]);
  end
%   for j = 1: size(ppkisave,1)
%     plot(ax,[ppkisave(j,(i-1)*2+1) ppkisave(j,(i-1)*2+1)], [0 ppkisave(j,i*2).*max(greenf(:,i))],...
%       '-','Color',[.5 .5 .5],'linew',1);
%   end
%   scatter(ax,ppkisave(:,(i-1)*2+1),ppkisave(:,i*2).*max(greenf(:,i)),10,[.5 .5 .5],'filled');
%   for j = 1: size(npkisave,1)
%     plot(ax,[npkisave(j,(i-1)*2+1) npkisave(j,(i-1)*2+1)], [0 npkisave(j,i*2).*min(greenf(:,i))],...
%       '-','Color',[.5 .5 .5],'linew',1);
%   end
%   scatter(ax,npkisave(:,(i-1)*2+1),npkisave(:,i*2).*min(greenf(:,i)),10,[.5 .5 .5],'filled');
  stem(ax,ppksrc(:,(i-1)*2+1),ppksrc(:,i*2).*max(greenf(:,i)),'LineStyle','-','MarkerFaceColor',...
    [.5 .5 .5],'MarkerEdgeColor',[.5 .5 .5],'Color',[.5 .5 .5],'linew',1,'MarkerSize',2);
  stem(ax,npksrc(:,(i-1)*2+1),npksrc(:,i*2).*min(greenf(:,i)),'LineStyle','-','MarkerFaceColor',...
    [.5 .5 .5],'MarkerEdgeColor',[.5 .5 .5],'Color',[.5 .5 .5],'linew',1,'MarkerSize',2); 

  ax.Box='on'; 
%   grid(ax,'on');
end
xlabel(f1.ax(3),'Samples');
ylabel(f1.ax(3),'Amplitude');

f2.fig = figure;
f2.fig.Renderer = 'painters';
widin = 12;  % maximum width allowed is 8.5 inches
htin = 6;   % maximum height allowed is 11 inches
set(f2.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*resol htin*resol]);
nrow = 4; 
ncol = 1;
for isub = 1:nrow*ncol
  f2.ax(isub) = subplot(nrow,ncol,isub);
end
%get the locations for each axis
optaxpos(f2,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);


ax=f2.ax(1);
hold(ax,'on');
%obtain the mean envelope
tmp = detrend(sigsta);
envup = envelope(tmp);
menv = median(envup,1);
for j = 1: size(sigsta,2)
  plot(ax, sigsta(:,j)./menv(j), '-','Color',color(j,:));
end
xlim(ax,[0,lsig/2]); %ylim(ax,yran);
for j = 1: size(zcrssrc,1)
  plot(ax,[zcrssrc(j,(2-1)*2+1) zcrssrc(j,(2-1)*2+1)], ax.YLim, '--','Color',[.7 .7 .7]);
end
text(ax,0.05,0.9,'each normalized by its median envelope',...
  'HorizontalAlignment','left','Units','normalized');
ax.Box='on';
%   grid(ax,'on');

ax=f2.ax(2);
hold(ax,'on');
for j = 1: size(zcrssrc,1)
  plot(ax,[zcrssrc(j,(2-1)*2+1) zcrssrc(j,(2-1)*2+1)], ax.YLim, '--','Color',[.7 .7 .7]);
end
for j = 1: size(sigsta,2)
  scatter(ax,zcrssrc(:,(j-1)*2+1),zcrssrc(:,j*2),10,color(j,:),'filled');
end
ax.Box='on';
%   grid(ax,'on');
xlim(ax,[0,lsig/2]);

ax=f2.ax(3);
hold(ax,'on');
%obtain the mean envelope
tmp = detrend(sigsta);
envup = envelope(tmp);
menv = median(envup,1);
for j = 1: size(sigsta,2)
  plot(ax, sigsta(:,j)./menv(j), '-','Color',color(j,:));
end
xlim(ax,[lsig/2,lsig]); %ylim(ax,yran);
for j = 1: size(zcrssrc,1)
  plot(ax,[zcrssrc(j,(2-1)*2+1) zcrssrc(j,(2-1)*2+1)], ax.YLim, '--','Color',[.7 .7 .7]);
end
ax.Box='on';
ylabel(ax,'Amplitude');

ax=f2.ax(4);
hold(ax,'on');
for j = 1: size(zcrssrc,1)
  plot(ax,[zcrssrc(j,(2-1)*2+1) zcrssrc(j,(2-1)*2+1)], ax.YLim, '--','Color',[.7 .7 .7]);
end
for j = 1: size(sigsta,2)
  scatter(ax,zcrssrc(:,(j-1)*2+1),zcrssrc(:,j*2),10,color(j,:),'filled');
end
ax.Box='on';
%   grid(ax,'on');
xlabel(ax,'Samples');
ylabel(ax,'Deconvolved source amplitude');
xlim(ax,[lsig/2,lsig]);


%%%beside the plot, we want to identify the nearest, ie, related pair of waveform pos&neg peaks and 
%%%the deconvolved sources shown by the scaled templates that also have a pos&neg peaks
mindisp = zeros(size(ppksrc,1), size(sigsta,2)); %min distance samples between src pos peaks and closest waveform pos peaks
clppkwf = zeros(size(ppksrc,1), size(sigsta,2));  %closest waveform pos peak index
clppkhtwf = zeros(size(ppksrc,1), size(sigsta,2));  %closest waveform pos peak height
ppkhtsrc = zeros(size(ppksrc,1), size(sigsta,2));  %*SCALED* height of src pos peaks
mindisn = zeros(size(ppksrc,1), size(sigsta,2)); %min distance samples between src neg peaks and closest waveform pos peaks
clnpkwf = zeros(size(ppksrc,1), size(sigsta,2));  %closest waveform neg peak index
clnpkhtwf = zeros(size(ppksrc,1), size(sigsta,2));  %closest waveform neg peak height
npkhtsrc = zeros(size(ppksrc,1), size(sigsta,2));  %*SCALED* height of src neg peaks

%identify closest peaks
for i = 1: size(sigsta,2) 
  for j = 1: size(ppksrc,1)
    %for pos peaks
    pki = ppkwf{i};
    pkhti = ppkhtwf{i};
    [mindisp(j,i),tempind] = min(abs(pki-ppksrc(j,(i-1)*2+1)));
    clppkwf(j,i) = pki(tempind);
    clppkhtwf(j,i) = pkhti(tempind);
    ppkhtsrc(j,i) = ppksrc(j,i*2).*max(greenf(:,i));
    %for pos peaks
    pki = npkwf{i};
    pkhti = npkhtwf{i};
    [mindisn(j,i),tempind] = min(abs(pki-npksrc(j,(i-1)*2+1)));
    clnpkwf(j,i) = pki(tempind);
    clnpkhtwf(j,i) = pkhti(tempind);
    npkhtsrc(j,i) = npksrc(j,i*2).*min(greenf(:,i));
  end
end

%make a structure array for fewer outputs
clppk.mindisp = mindisp;
clppk.clppkwf = clppkwf;
clppk.clppkhtwf = clppkhtwf;
clppk.ppkhtsrc = ppkhtsrc;
clnpk.mindisn = mindisn;
clnpk.clnpkwf = clnpkwf;
clnpk.clnpkhtwf = clnpkhtwf;
clnpk.npkhtsrc = npkhtsrc;
 
% keyboard


