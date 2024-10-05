function [dupinds,induni,inddup]=finddupind(mat)
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
% --Now, for a cluster, not only the time separation between N and N-m needs to
% be smaller than 'dtcut', but also the max time separation between each
% consecutive events needs to smaller than 0.25+0.125 s. 
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/11/04
% Last modified date:   2024/01/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[matuni, induni] = unique(mat,'row','stable');

inddup = setdiff((1: size(mat,1))', induni);
ndup = length(inddup);

dupinds =  cell(ndup,1);
for j = 1: ndup
  dupinds{j} = find(ismember(mat, mat(inddup(j),:), 'rows'));
end


