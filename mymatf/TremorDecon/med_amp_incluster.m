function [ampplt,dtarvlplt,indimpdtcut]=med_amp_incluster(nbst,imp,nsrc,m,dtcutspl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ampplt,dtarvlplt]=med_amp_incluster(nbst,imp,nsrc,mmax,m)
%
% This function is to compute the median amp of the 
% cluster defined by consecutive events within N & N-m source pairs.
% It also returns the diff time between source pair N & N-m, for a certain m.
% Deal with data and synthetic noise only, NOT for synthetics.
%
% --Added another col to 'ampplt' and 'dtarvlplt' to recognize which burst is it
% from. Also added the indices of consecutive evts N,N+1,...,N+m inside the 
% cluster that is defined by evt pair N and N-m. This is returned only when a
% diff time cutoff in samples 'dtcutspl' is assigned. 
% 
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/11/04
% Last modified date:   2023/11/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('dtcutspl',[]);

% nsrcsep = nsrc-m;
% nsrcsep(nsrcsep<0) = 0;
ampplt = []; % the median amp for all sources in between N & N-m pair
dtarvlplt = [];
inddtcutbst = cell(nbst,1);
for i = 1: nbst
  if nsrc(i) == 0
    continue
  end
  ist = sum(nsrc(1:i-1))+1;
  ied = ist+nsrc(i)-1;
  impi = imp(ist:ied,:);
  %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
  dtarvl = diffcustom(impi(:,1), m,'forward');
%   dtarvl = srcdistNtoNm(impi(:,1),impi(:,7:8),mmax);
  %   dtarvlnnsepall = [dtarvlnnsepall; dtarvl{nsep}];
  if isempty(dtarvl)
    continue
  end
  
  ampcont = [];
  for j = 1: nsrc(i)-m
    impcont = impi(j: j+m,:);
    ampcont(j,1) = median(mean(impcont(:,[2 4 6]),2));
  end
  
  %if assigned any dtcut, then return the indices of evts in the clusters defined by N&N-m
  if ~isempty(dtcutspl)
    ind = find(dtarvl <= dtcutspl); 
    dtarvl = dtarvl(ind);
    ampcont = ampcont(ind); %cut correspondingly
    indimp = [];  %indices of evts in the clusters whose diff time between start and end is w/i dtcut
    for j = 1: length(ind)
      indimp{j,1} = reshape(ind(j):ind(j)+m, [], 1);
    end
    inddtcutbst{i} = indimp;  %save separately by burst
  end
  
  dtarvlplt = [dtarvlplt; [dtarvl i*ones(size(dtarvl,1),1)]]; %add 2nd col as burst #
  ampplt = [ampplt; [ampcont i*ones(size(ampcont,1),1)]]; %add 2nd col as burst #
  
end

if ~isempty(dtcutspl)
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

