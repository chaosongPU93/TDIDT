% test envelope

%%
%%% load seismic trace file in that day.
datapath = strcat(getenv('ALLAN'),'/data-no-resp/PGCtrio');
datestr = '2005255';
yr = datestr(1:4);
day = datestr(5:end);
sps=40;
winlensec=4;
winoffsec=1;        % window offset in sec, which is the step of a moving window
winlen=winlensec*sps;      % length in smaples
hi=6.5;
lo=1.25;
loff = 2.1;
ccmin = 0.44;
fam = '002';
npo = 2;
npa = 2;
mshift = 29;
nstanew = 4;
IDENTIF = strcat(yr,'.',day,'.',fam,'.loff',num2str(loff),'.ccmin',num2str(ccmin),'.nponpa', ...
    num2str(npo),num2str(npa),'.ms',num2str(mshift));
oritrace = load(strcat(datapath, '/MAPS/seistraceoneday_',IDENTIF,'_',num2str(lo),'-', ...
    num2str(hi),'_',num2str(winlen/sps),'s',num2str(sps),'sps'));   % original trace

%%
close all
figure
plot(oritrace(100:end,1), oritrace(100:end,2), 'b-','linewidth',1); hold on

[upana,loana] = envelope(oritrace(100:end,2),500,'analytic');   % use the hilbert transform
[upana2,loana2] = envelope(oritrace(100:end,2));    % the default would use a pretty long hilbert transform window

% [uppeak,lopeak] = envelope(oritrace(100:end,2),2,'peak');     % the peak envelop simply get the local min/max usually shorter than analytic option

plot(oritrace(100:end,1),upana, 'k-','linewidth',2); hold on
plot(oritrace(100:end,1),loana, 'k-','linewidth',2); hold on
% 
plot(oritrace(100:end,1),upana2, 'r-','linewidth',2); hold on
plot(oritrace(100:end,1),loana2, 'r-','linewidth',2); hold on

xlim([10 20]);

