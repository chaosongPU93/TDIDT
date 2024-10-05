function [impnew,nsrcnew,lsignew] = discard_badwins_decon_synth(imp,badwins,windows,sps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function to read the misaligned windows 'badwins', and discard the
% detections that belong to those windows from the whole set 'impi',
% according to the complete window assignment 'windows'.
%
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/07/27
% Last modified date:   2024/07/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('sps',160);

%misaligned windows
% badwinsi = badwins{insat,ireg};

%window segmentations
% windows = irccran{insat,ireg};
nwin =  size(windows,1);
goodwinsi = setdiff((1:nwin)', badwins); %normal windows
ngoodw = length(goodwinsi);

%range of good windows, useful to know the length of time
windowsgood = windows(goodwinsi,:);
lsignew= sum((windowsgood(:,2)-windowsgood(:,1)+1)/sps);

%corresponding windows that detections belong to
iwin = findwhichrange(imp(:,1), windows);

%find all detections that belong to the good wins
impnew = [];
for j = 1: ngoodw
  goodwin = goodwinsi(j);
  impiwin = imp(iwin==goodwin, :);
  impnew = [impnew; impiwin];
end

%return sources and count
impnew = sortrows(impnew, 1);
nsrcnew = size(impnew, 1);
