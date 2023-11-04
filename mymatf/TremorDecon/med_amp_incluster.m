function [ampplt,dtarvlplt]=med_amp_incluster(nbst,imp,nsrc,mmax,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ampplt,dtarvlplt]=med_amp_incluster(nbst,imp,nsrc,mmax,m)
%
% This function is to compute the median amp of the 
% cluster composed by consecutive events within N & N-m source pairs.
% It also returns the diff time between source pair N & N-m.x
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/11/04
% Last modified date:   2023/11/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  nsrcsep = nsrc-m;
  nsrcsep(nsrcsep<0) = 0;
  ampdt = []; % the median amp for all sources in between N & N-m pair
  dtarvlplt = [];
  for i = 1: nbst
    if nsrc(i) == 0
      continue
    end
    ist = sum(nsrc(1:i-1))+1;
    ied = ist+nsrc(i)-1;
    impi = imp(ist:ied,:);
    %between Nth and (N-1)th source; Nth and (N-2)th; Nth and (N-3)th
    dtarvl = srcdistNtoNm(impi(:,1),impi(:,7:8),mmax);
    %   dtarvlnnsepall = [dtarvlnnsepall; dtarvl{nsep}];
    if isempty(dtarvl{m})
      continue
    end
    dtarvlplt = [dtarvlplt; dtarvl{m}];
    ampcont = [];
    for j = 1: nsrc(i)-m
      impcont = impi(j: j+m,:);
      ampcont(j,1) = median(mean(impcont(:,[2 4 6]),2));
    end
    ampdt = [ampdt; ampcont];
  end
  ampplt = ampdt;
