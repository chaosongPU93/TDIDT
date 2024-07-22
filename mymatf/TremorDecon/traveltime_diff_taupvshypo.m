% function traveltime_diff_taupvshypo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to use the customized velocity model Puget Sound P3 
% velocity model created by 'taupcreate' to compute the theoretical travel
% time difference at the Cascadia stations, given a set of source locations.
% The idea is to get a map of differential travel time relative to the
% main station at any other station from a grid of locations (which come
% from the inversion using hypoinverse of a grid of differential travel times
% between the trio stations, off12 and off13). Since we are using different
% packages for raytracing, it is necessary to compare the traveltime estimate
% from Taup to Hypoinverse.
% 
% 
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/02/17
% Last modified date:   2022/02/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clc;
clear;
close all;

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);

%% generate a grid of offsets around fam 002, then transform the offset grid to relative locations
sps = 160;
offmax = round(4*sps/40); % make it compatible to different sps
%grid of time offsets
i = 1;
for off12 = -offmax: 1: offmax
  for off13 = -offmax: 1: offmax
    offset(i,1:2) = [off12 off13];  % offset in samples at requested sampling rate.
    i = i+1;
  end
end

%transform to locations
ftrans = 'interpchao';
[evtloc, indinput] = off2space002(offset,sps,ftrans,0); 
nevt = size(evtloc,1);

%% compute the travel time difference using Taup at all stations
stas=['PGC  '
      'SSIB '
      'SILB '
      'KLNB '
      'LZB  '
      'TWKB '
      'MGCB '
      ];
nsta = size(stas,1);

fid = fopen(fullfile('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary','staloc.txt'), 'r');
datalist = textscan(fid, '%s %s %s \n');  % 3 cell arrays, station name, lat, lon 
stanm = char(length(datalist{1}), 5);  % station name array
% stalist = [];
for i = 1: length(datalist{1})
  str = datalist{1}{i};
  if length(str)==3
    str = strcat(str,{'  '});
  elseif length(str)==4
    str = strcat(str,{' '});
  end
  stanm(i,1:5) = char(str);
end

staloc = [];  % station loc, lon, lat
for ista = 1: size(stas,1)
  [~,idx] = ismember(stas(ista,:),stanm,'rows');
  staloc(ista,:) = [str2double(datalist{3}{idx}) str2double(datalist{2}{idx})];
end

pugetp3 = '/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/pugetp3.taup';  % pre-created vel model for Taup
tts = zeros(nevt, nsta);  % S-wave travel time
for i = 1: nevt
  for j = 1: nsta
    tt=tauptime('mod',pugetp3,'dep',evtloc(i,5),'ph','s,S','evt',[evtloc(i,4) evtloc(i,3)],...
      'sta',[staloc(j,2) staloc(j,1)]); 
    tts(i,j) = tt(1).time;   % choose the first S arrival
    if ~strcmp(tt(1).phase,'s')
      fprintf('At i=%d, j=%d, the take-off angle might be <90 degrees, ie., downgoing wave \n.',...
        i,j);
      break
    end
  end
end
          
dtts_taup = tts(:,1)-tts;  % differential traveltime realative to main sta PGC, in sec
dttsspl_taup = dtts_taup*sps; % differential traveltime in samples from taup

%%
PERMSTA=['PGC  '        % permanent station names
    'LZB  '];
POLSTA=['SSIB '           % polaris station names
    'SILB '
    'KLNB '
    'MGCB '
    'TWKB '];

%%%We need the 4th col of the rot params, for the traveltime difference at stations from the lfe fam
workpath = getenv('ALLAN');
datapath = strcat(workpath,'/data-no-resp');
CATA = 'new';
FLAG = 'PGC'; % detector
fam = '002';   % family number
sft2=0;     % centroid shift of station 2
sft3=0;     % centroid shift of station 3
[PERMROTS, POLROTS] = GetRotsCommon(FLAG,fam,CATA,datapath,sft2,sft3);
if ~isequal(CATA, 'fixed')
  reftime = PERMROTS(1,4);     % reftime is the reference time at the 1st station,depends on the choice of 1st station
  PERMROTS(:,4) = PERMROTS(:,4)-reftime;    % to make sure that 1st station is 0
  POLROTS(:,4) = POLROTS(:,4)-reftime;
end

constpgc = zeros(1,nsta);  % differential traveltime between each sta and main sta at the origin, ie. offset 0,0
%Assign to each station correctly
for i = 1: nsta
  stas(i,:)
  [LIA,idx] = ismember(stas(i,:),PERMSTA,'rows');
  if LIA
    constpgc(i) = PERMROTS(idx,4)*sps/40;  % in samples, they were computed at 40 sps
  else
    [LIA,idx] = ismember(stas(i,:),POLSTA,'rows');
    if LIA
      constpgc(i) = POLROTS(idx,4)*sps/40;
    else
      fprintf('Station %s does not have a match',stas(i,:));
    end
  end
end

% differential traveltime in samples from hypoinverse
dttsspl_hypo = [offset(:,1)-constpgc(2) offset(:,2)-constpgc(3)];  % similar to hypoinverse input
dtts_hypo = dttsspl_hypo/sps;   % differential traveltime in sec

%% plot the difference of prediction on differential traveltime between Hypoinverse and Taup
figure
difdtts = [dttsspl_hypo(:,1)-dttsspl_taup(:,2) dttsspl_hypo(:,2)-dttsspl_taup(:,3)];
plot(difdtts(:,1),'b-'); hold on
plot(difdtts(:,2),'r-');
ax = gca;
plot(ax.XLim,[median(difdtts(:,1)) median(difdtts(:,1))], 'b--');
text(sum(ax.XLim)*0.1, sum(ax.YLim)*0.5, sprintf('Median: %.4f sec',median(difdtts(:,1))/sps),...
  'FontSize',12);
plot(ax.XLim,[median(difdtts(:,2)) median(difdtts(:,2))], 'r--');
text(sum(ax.XLim)*0.5, sum(ax.YLim)*0.5, sprintf('%.4f sec',median(difdtts(:,2))/sps),...
  'FontSize',12);
text(sum(ax.XLim)*0.1, sum(ax.YLim)*0.4, sprintf('Std: %.4f sec',std(difdtts(:,1))/sps),...
  'FontSize',12);
text(sum(ax.XLim)*0.5, sum(ax.YLim)*0.4, sprintf('%.4f sec',std(difdtts(:,2))/sps),...
  'FontSize',12);
%below is the difference of prediction on differential traveltime from the origin, hopefully it is
%close to the median above.
ind = find(offset(:,1)==0 & offset(:,2)==0);  % index of the origin, ie. offset 0,0
difdttsori = [dttsspl_hypo(ind,1)-dttsspl_taup(ind,2) dttsspl_hypo(ind,2)-dttsspl_taup(ind,3)];
plot(ax.XLim,[difdttsori(1) difdttsori(1)], 'b:');
plot(ax.XLim,[difdttsori(2) difdttsori(2)], 'r:');
xlabel('Event index');
ylabel(sprintf('Difference in samples at %d Hz',sps));

%%%relative difference, looks like hypo predict a larger differential travel 
figure
rdifdtts = [(dttsspl_hypo(:,1)-dttsspl_taup(:,2))./dttsspl_taup(:,2) ...
            (dttsspl_hypo(:,2)-dttsspl_taup(:,3))./dttsspl_taup(:,3)];
plot(rdifdtts(:,1),'b-'); hold on
plot(rdifdtts(:,2),'r-');
ax = gca;
plot(ax.XLim,[median(rdifdtts(:,1)) median(rdifdtts(:,1))], 'b--');
text(sum(ax.XLim)*0.1, sum(ax.YLim)*0.5, sprintf('Median: %.2f %%',100*median(rdifdtts(:,1))),...
  'FontSize',12);
plot(ax.XLim,[median(rdifdtts(:,2)) median(rdifdtts(:,2))], 'r--');
text(sum(ax.XLim)*0.5, sum(ax.YLim)*0.5, sprintf('%.2f %%',100*median(rdifdtts(:,2))),...
  'FontSize',12);
text(sum(ax.XLim)*0.1, sum(ax.YLim)*0.4, sprintf('Std: %.2f %%',100*std(rdifdtts(:,1))),...
  'FontSize',12);
text(sum(ax.XLim)*0.5, sum(ax.YLim)*0.4, sprintf('%.2f %%',100*std(rdifdtts(:,2))),...
  'FontSize',12);

%below is the difference of prediction on differential traveltime from the origin, hopefully it is
%close to the median above.
rdifdttsori = [(dttsspl_hypo(ind,1)-dttsspl_taup(ind,2))/dttsspl_taup(ind,2) ...
               (dttsspl_hypo(ind,2)-dttsspl_taup(ind,3))/dttsspl_taup(ind,3)];
plot(ax.XLim,[rdifdttsori(1) rdifdttsori(1)], 'b:');
plot(ax.XLim,[rdifdttsori(2) rdifdttsori(2)], 'r:');
xlabel('Event index');
ylabel(sprintf('Relative difference at %d Hz',sps));

%% Use the difference to correct Taup prediction at other stations, and scatter the results
%Correct the Taup prediction to get the the diff traveltime under the framework of hypoinverse 
%Looks like the difference, so is the correction may be non-linear, so we may have to use the
%difference from the source at the center 0,0

%%%Use absolute correction
difori = constpgc-dttsspl_taup(ind,:);  % diff in samples, the first 3 are the same as difdttsori
dttsspl_taupc1 = dttsspl_taup;
for i = 1: nsta
  dttsspl_taupc1(:,i) = dttsspl_taupc1(:,i)+difori(i);
  offpred1(:,i)=constpgc(i)-dttsspl_taupc1(:,i);
end
%At trio stations, horizontal or vertical band is expected
figure
for i = 1:3
  subplot(2,2,i)
  scatter(offset(:,1)/sps,offset(:,2)/sps,50,offpred1(:,i)/sps,'filled');
%   oldc = colormap('kelicol');
%   newc = flipud(oldc);
%   colormap(newc);
  colormap('bluewhitered');
  c=colorbar;
  title(sprintf('Offset (s) PGC-%s', stas(i,:)));
end          
xlabel(sprintf('Offset (s) PGC-SSIB'));
ylabel(sprintf('Offset (s) PGC-SILB'));
          
figure
for i = 4:nsta
  subplot(2,2,i-3)
  scatter(offset(:,1)/sps,offset(:,2)/sps,50,offpred1(:,i)/sps,'filled');
%   oldc = colormap('kelicol');
%   newc = flipud(oldc);
%   colormap(newc);
  colormap('bluewhitered');
  c=colorbar;
  title(sprintf('Offset (s) PGC-%s', stas(i,:)));
end          
xlabel(sprintf('Offset (s) PGC-SSIB'));
ylabel(sprintf('Offset (s) PGC-SILB'));
          
%%
%%%Use relative correction
rdifori = (constpgc-dttsspl_taup(ind,:))./dttsspl_taup(ind,:);  % relative difference
rdifori(1) = 0;
%Use the average from sta1-2 and 1-3 to represente the general difference for all other stations
% rdif = (median(rdifdtts(:,1))+median(rdifdtts(:,2)))/2; 
% rdif = (rdifdttsori(1)+rdifdttsori(2))/2; 
dttsspl_taupc2 = dttsspl_taup;
for i = 1: nsta
  dttsspl_taupc2(:,i) = dttsspl_taupc2(:,i).*(1+rdifori(i));
  offpred2(:,i)=constpgc(i)-dttsspl_taupc2(:,i);
end
%At trio stations, horizontal or vertical band is expected
figure
for i = 1:3
  subplot(2,2,i)
  scatter(offset(:,1)/sps,offset(:,2)/sps,50,offpred2(:,i)/sps,'filled');
%   oldc = colormap('kelicol');
%   newc = flipud(oldc);
%   colormap(newc);
  colormap('bluewhitered');
  c=colorbar;
  title(sprintf('Offset (s) PGC-%s', stas(i,:)));
end          
xlabel(sprintf('Offset (s) PGC-SSIB'));
ylabel(sprintf('Offset (s) PGC-SILB'));
          
figure
for i = 4:nsta
  subplot(2,2,i-3)
  scatter(offset(:,1),offset(:,2),50,offpred2(:,i)/sps,'filled');
%   oldc = colormap('kelicol');
%   newc = flipud(oldc);
%   colormap(newc);
  colormap('bluewhitered');
  c=colorbar;
  title(sprintf('Offset (s) PGC-%s', stas(i,:)));
end
xlabel(sprintf('Offset (s) PGC-SSIB'));
ylabel(sprintf('Offset (s) PGC-SILB'));


          
          
          
          
          
          
          
          
          