function [dloc,dt,dloc_spl,dloc2all,dt2all,dloc2all_spl]=...
  dloc_evtcustom(impiwin,implociwin,sps,ftrans,m,timetype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [dloc,dt,dloc_spl]=dloc_evtcustom(impiwin,implociwin,sps,ftrans,m,n,timetype)
%
% This function is to obtain the location and time difference between events
% in a custom cluster. Events and their locations are stored in 'impiwin' and
% 'implociwin'. 
% --It is similar to 'dloc_evtcluster'.
%
% 2024/02/06, mapping from samples to relation location already contains error,
% so when you try to obtain the location difference, or distance, start with 
% the difference in samples first. For example, if you ask the difference in 
% certain directions, first obtain the difference in off12 and off13, then map
% the diff to location, then you can project along any direction to get the
% loc difference along that direction, or the distance (abs of loc diff) in
% that direction, etc. In all, it is NOT recommended to use output 'dloc' 
% unless a lof of data points would average out the error. 
% 
% 2024/02/15, looks like if you need the loc difference in map-view, you have
% to first invert time offsets to map locations, then obtain the difference.
% Despite the possible location error propagation, you need it. Since if you 
% start from difference in samples, then transform them to map locations,
% you are essentially assuming one of the two sources is located at the 
% origin, which is not the case in reality, even if the actual difference in
% the two ways may be small.
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/02/22
% Last modified date:   2024/02/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('timetype','tarvl');

if strcmp(timetype,'tori')
  %return the origin time, sorted srcs, and indices
  [tevt,impiwin,indsort]=tarvl2tori(impiwin,sps,ftrans,1);
  implociwin = implociwin(indsort, :);
else
  tevt=impiwin(:,1);
end

%%%between N and N-n, n<=m
dt = cell(m, 1); %differential time, since it is sorted by time, it is always positive
dloc = cell(m, 1); %differential location map-view
dloc_spl = cell(m, 1); %differential location in samples
for i = 1: m
  dt{i} = diffcustom(tevt,i,'forward');
  dloc{i} = diffcustom(implociwin(:,1:2),i,'forward'); %relative loc in km
  dloc_spl{i} = diffcustom(impiwin(:,7:8),i,'forward'); %time offset in samples
end

%%%between each source to all others
[dt2all,dloc2all] = srcdistall(tevt,implociwin(:,1:2));
[~,dloc2all_spl] = srcdistall(tevt,impiwin(:,7:8));

% nsrc = size(impiwin,1);
% dt2all = cell(nsrc, 1); %differential time, since it is sorted by time, it is always positive
% dloc2all = cell(nsrc, 1); %differential location
% dloc2all_spl = cell(nsrc, 1); %absolute Euclidean distance
% 
% %loop for every source
% for i = 1: nsrc-1
%   othersrc = i+1:nsrc;
%   tempt = tevt(othersrc);
%   temploc = implociwin(othersrc,:);
%   temploc_spl = impiwin(othersrc,:);
% 
%   dt2all{i} = abs(tempt-tevt(i));  
%   dloc2all{i} = [temploc(:,1)-implociwin(i,1) temploc(:,2)-implociwin(i,2)];
%   dloc2all_spl{i} = [temploc_spl(:,7)-impiwin(i,7) temploc_spl(:,8)-impiwin(i,8)];
%   
% end
% 
% dt2all = cat(1,dt2all{:});
% dloc2all = cat(1,dloc2all{:});
% dloc2all_spl = cat(1,dloc2all_spl{:});
