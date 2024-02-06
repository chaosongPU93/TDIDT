function [ampplt,dtplt,indimpdtcut]=med_amp_incluster(nbst,imp,nsrc,m,dtcut,sps,timetype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ampplt,dtplt,indimpdtcut]=med_amp_incluster(nbst,imp,nsrc,m,dtcut,sps,timetype)
%
% This function is to compute the median amp of the 
% cluster defined by consecutive events within N & N-m source pairs.
% It also returns the diff time between source pair N & N-m, for a certain m.
% Deal with data and synthetic noise only, NOT for synthetics.
%
% --Added another col to 'ampplt' and 'dtplt' to recognize which burst is it
% from. Also added the indices of consecutive evts N,N+1,...,N+m inside the 
% cluster that is defined by evt pair N and N-m. This is returned only when a
% diff time cutoff in sec 'dtcut' is assigned. 
% --Added option for choosing the type of time of each event, 'timetype' 
% either using 'tarvl' or 'tori'. 
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/11/04
% Last modified date:   2024/01/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('dtcut',[]);
defval('sps',160);
defval('timetype','tarvl');

ftrans = 'interpchao';
[imploc0, ~] = off2space002([0 0],sps,ftrans,0);  % a ref source at 0,0

% nsrcsep = nsrc-m;
% nsrcsep(nsrcsep<0) = 0;
ampplt = []; % the median amp for all sources in between N & N-m pair
dtplt = []; %differential time, can be arrival or origin time
inddtcutbst = cell(nbst,1);
for i = 1: nbst
  if nsrc(i) == 0
    continue
  end
  ist = sum(nsrc(1:i-1))+1;
  ied = ist+nsrc(i)-1;
  impi = imp(ist:ied,:);
  tevt = impi(:,1);

  %which type of time to use
  if strcmp(timetype,'tori')
    [tevt,impi]=tarvl2tori(impi,sps,ftrans,1);  %return the origin time, sorted  
  end

  %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
  dtevt = diffcustom(tevt, m,'forward');
  if isempty(dtevt)
    continue
  end
  
  %meidan amp of the inclusive cluster
  ampcont = [];
  for j = 1: nsrc(i)-m
    impcont = impi(j: j+m,:);
    ampcont(j,1) = median(mean(impcont(:,[2 4 6]),2));
  end
  
  %if assigned any dtcut, then return the indices of evts in the clusters defined by N&N-m
  if ~isempty(dtcut)
    ind = find(dtevt <= dtcut*sps); 
    dtevt = dtevt(ind);
    ampcont = ampcont(ind); %cut correspondingly
    indimp = [];  %indices of evts in the clusters whose diff time between start and end is w/i dtcut
    for j = 1: length(ind)
      indimp{j,1} = reshape(ind(j):ind(j)+m, [], 1);
    end
    inddtcutbst{i} = indimp;  %save separately by burst
  end
  
  dtplt = [dtplt; [dtevt i*ones(size(dtevt,1),1)]]; %add 2nd col as burst #
  ampplt = [ampplt; [ampcont i*ones(size(ampcont,1),1)]]; %add 2nd col as burst #
  
end

if ~isempty(dtcut)
  indimpdtcut = cat(1,inddtcutbst{:});  %concatenate all bursts
%   indimpdtcut = inddtcutbst;
else
  indimpdtcut = [];
end

% for i = 1: nbst
%   aa(i) = size(inddtcutbst{i},1);
% end
% sum(aa)

% keyboard

