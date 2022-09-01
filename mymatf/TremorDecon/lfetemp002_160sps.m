% lfetemp002_160sps.m
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this script, this script just simply create the LFE templates at fam 
% 002 under 160 sps.
% 
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2022/04/14
% Last modified date:   2022/04/14
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

fam = '002';
CATA = 'new';
sps = 160;
templensec = 60;
ccmethod = 2; % only using lfes passed the cc thresholds (ccmethod=1); using all lfes (ccmethod=2)  
ccbp = [2 8];
plflag = 0;

[dstack,ccstack,dstackort,ccstackort] = mk_bbtemp_PGC(fam,CATA,sps,templensec,ccmethod,ccbp,plflag);
           
%write into files, NOTE the remade stacks contain all 7 stations
allstas=['PGC  '
  'SSIB '
  'SILB '
  'LZB  '
  'TWKB '
  'MGCB '
  'KLNB '];
allnsta = size(allstas,1);


%%% plot the 2 kinds of template together
figure;
mid = size(dstack,2)/2;
for ista = 1: allnsta
  subplot(allnsta,1,ista);
  plot(dstack(ista, mid-5*sps: mid+5*sps), 'k', 'linewidth', 1); hold on
  plot(ccstack(ista, mid-5*sps: mid+5*sps), 'r', 'linewidth', 1);
  title(strcat('Stacked template with/without CC at station', {' '}, strtrim(allstas(ista, :))));
  legend('Direct stack', 'CC stack');
  xlabel('Samples');
  ylabel('Amplitude');
  box on
  grid on
end

figure;
for ista = 1: allnsta
  subplot(allnsta,1,ista);
  plot(dstackort(ista, mid-5*sps: mid+5*sps), 'k', 'linewidth', 1); hold on
  plot(ccstackort(ista, mid-5*sps: mid+5*sps), 'r', 'linewidth', 1);
  title(strcat('Stacked template (orthogonal) with/without CC at station', {' '}, ...
    strtrim(allstas(ista, :))));
  legend('Direct stack', 'CC stack');
  xlabel('Samples');
  ylabel('Amplitude');
  box on
  grid on
end

%write into files
workpath = getenv('ALLAN');
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/PGCtrio/');
suffix = '_catnew';
for ista = 1: allnsta
  fidds = fopen(strcat(temppath, fam, '_', strtrim(allstas(ista, :)), '_', ...
    num2str(sps), 'sps_', num2str(templensec), 's_', ...
    'BBDS_', 'opt_Nof_Non_Chao',suffix), 'w+');  % bandpassed direct stack, no filter, no norm
  fidccs = fopen(strcat(temppath, fam, '_', strtrim(allstas(ista, :)), '_', ...
    num2str(sps), 'sps_', num2str(templensec), 's_', ...
    'BBCCS_', 'opt_Nof_Non_Chao',suffix), 'w+');  % bandpassed cc stack, no filter, no norm
  fprintf(fidds, '%f \n', dstack(ista, :)');
  fclose(fidds);
  fprintf(fidccs, '%f \n', ccstack(ista, :)');
  fclose(fidccs);
  
  fidds = fopen(strcat(temppath, fam, '_', strtrim(allstas(ista, :)), '_', ...
    num2str(sps), 'sps_', num2str(templensec), 's_', ...
    'BBDS_', 'ort_Nof_Non_Chao',suffix), 'w+');  % bandpassed direct stack, no filter, no norm
  fidccs = fopen(strcat(temppath, fam, '_', strtrim(allstas(ista, :)), '_', ...
    num2str(sps), 'sps_', num2str(templensec), 's_', ...
    'BBCCS_', 'ort_Nof_Non_Chao',suffix), 'w+');  % bandpassed cc stack, no filter, no norm
  fprintf(fidds, '%f \n', dstackort(ista, :)');
  fclose(fidds);
  fprintf(fidccs, '%f \n', ccstackort(ista, :)');
  fclose(fidccs);
end


