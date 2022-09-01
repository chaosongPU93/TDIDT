% function puget_sound_p3_vel_model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to analyze the Puget Sound P3 velocity model. Visualize
% and modify it to create an eligible input for taup_create. We want to 
% generate a taup-compatible customized model based on it. Then we could
% use Taup to compute the travel time difference at stations in Cascadia given
% a set of source locations.
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/02/17
% Last modified date:   2020/02/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
clc;
clear;
close all;

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);

%% load the velocity model
data = load('/home/data2/chaosong/Seisbasics/hypoinverse/forsummary/nyvel');
vs = data(:,1);
dep = data(:,2);
vp = [5.40
      6.38
      6.59
      6.73
      6.86
      6.95
      7.80];
oldmod = [dep vp vs];

% plot it      
f.fig=figure;
ax = gca;
hold(ax,'on');
stairs(ax,vp,dep,'linew',0.5,'color','b');
stairs(ax,vs,dep,'linew',0.5,'color','r');
set(ax,'YDir','reverse');
set(ax,'XAxisLocation','top');
box on;
grid on;
text(0.83,0.91,'Puget Sound P3','FontSize',12,'unit','normalized','horizontalalignment','right',...
  'EdgeColor','k','Margin',2);
xlabel('Velocity (km/s)');
ylabel('Depth (km)');

%% write into file
fid = fopen('/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/ak135.tvel', 'r');
format = '%f %f %f %f \n';
dcell = textscan(fid,format,'headerlines',2); 
ak135(:,1) = dcell{1,1};    % dep
ak135(:,2) = dcell{1,2};    % vp
ak135(:,3) = dcell{1,3};    % vs

indak = find(ak135(:,1)>200,1,'first');

newmod = [dep vp vs;
          ak135(indak:end,1) ak135(indak:end,2) ak135(indak:end,3)];
fid = fopen('/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/pugetp3.tvel', 'w+');
fprintf(fid,'%s \n', 'pugetp3 - P');
fprintf(fid,'%s \n', 'pugetp3 - S');
fprintf(fid,'%10.4f %10.4f %10.4f\n',newmod');
fclose(fid);

f.fig=figure;
ax = gca;
hold(ax,'on');
stairs(ax,newmod(:,2),newmod(:,1),'linew',0.5,'color','b');
stairs(ax,newmod(:,3),newmod(:,1),'linew',0.5,'color','r');
set(ax,'YDir','reverse');
set(ax,'XAxisLocation','top');
box on;
grid on;
text(0.83,0.91,'Puget Sound P3','FontSize',12,'unit','normalized','horizontalalignment','right',...
  'EdgeColor','k','Margin',2);
xlabel('Velocity (km/s)');
ylabel('Depth (km)');

%% use taupcreate to make a new taup-style model for further use
% In default, taupcreate in mattaup creat the model in the path same as the mfile, the .taup
% velocity model can ONLY been called by absolute path;
% This is different from the linux version taup_create, which can be called as long as the model is
% placed into the StdModel folder. And models from the two seem UNcompatible with each other. 
% '.tvel' file can be read by linux and matlab, but you need to create '.taup' separately in
% different platforms

intvel1 = '/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/pugetp3.tvel';
outtaup1 = '/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/pugetp3.taup';
taupcreate(intvel1, outtaup1);    

%% check if it works
tic
pugetp3 = '/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/pugetp3.taup';
tt=tauptime('mod',outtaup1,'dep',36.0400,'ph','s,S','evt',[48.415167 -123.632167 ],...
            'sta',[48+39.0/60 -(123+27.05/60)]);
tts = tt(1).time   % choose the first S arrival
rayp = tt(1).rayparameter  % unit in [s/deg]toc


toc
% tt=tauptime('mod',pugetp3,'dep',30,'ph','s,S','dist',100);


% jma2001 = '/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/jma2001_sparse.taup';
% tt=tauptime('mod',jma2001,'dep',30,'ph','s,S','dist',100);
% 









