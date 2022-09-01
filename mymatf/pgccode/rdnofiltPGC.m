function [STAopt,timeperm] = ...
         rdnofiltPGC(fam,timoffrot,datapath,stas,sps,PERMSTA,PERMROTS,POLSTA,POLROTS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the function to make it easier to read one day of station-response-free
% data without filtering 
%
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2019/11/02
% Last modified date:   2019/11/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% 'scalefact' scales templates; 'scaleseisms' scales seisms.  Strategy changes with family.
if isequal(fam,'002')
    scaleseisms=[1.0 0.76 0.95];       % scaleseisms scales seismograms
elseif isequal(fam,'068')
    scaleseisms=[1.0 1.0 1.0];
end

%Which days of data to read?
year=timoffrot(1);
YEAR=int2str(year);
jday=timoffrot(2);
if jday <= 9
    JDAY=['00',int2str(jday)];
elseif jday<= 99
    JDAY=['0',int2str(jday)];
else
    JDAY=int2str(jday);
end
%     MO = 'SEP';
MO=day2month(jday,year);     % EXTERNAL function, day2month, get the month of one particular date
direc=[datapath, '/', YEAR,'/',MO,'/'];     % directory name
fprintf('%s \n', direc);
prename=[direc,YEAR,'.',JDAY,'.00.00.00.0000.CN'];    %  path plus prefix of data file,

if year==2003   % in 2003, station KLNB is named KELB, so use KELB to replace KLNB in 2003
    POLSTA(3,:)='KELB ';
else
    POLSTA(3,:)='KLNB ';  % remember to change it back
end

%Read the data; find glitches (many consecutive zeros) and flag those windows in STAnzeros.
%Get timsSTA from the permanent stations (last one over-writes):
nsta = size(stas,1);
lenofday = 24*3600;
STAopt=zeros(nsta,sps * lenofday);
%STAort=STAopt;
fileflag = 1;   % 1 means all files exist, the file status is normal
for ista=1:nsta
    found=0;
    [LIA,idx]=ismember(stas(ista,:),PERMSTA,'rows');     % ismember, to determine whether each row of stas is contained in PERMSTA, return logical value 1/0 and index
    if LIA
        found=found+LIA;
        if strcmp(PERMSTA(idx,1:3),'PGC')     % string compare
            %%% WAHT is the meanning of 'fact', similar to instrument
            %%% response
            fact=1.0e-3;
        elseif strcmp(PERMSTA(idx,1:3),'LZB')
            fact=1.8e-3;
        end
        %%% readperms is an EXTERNAL FUNCTION
        % opt: optimal seismogram after rotations
        % ort: orthogonal seismogram after rotations
        % nzeros: number of zeros in the trace
        % timsSTA: time sequence
        fname = strcat(prename,'.',PERMSTA(idx,1:3),'..BHE.D.SAC');
        if isfile(fname)    % if have the data file
            %                 % 1. this is for data without removing station response
            %                 [opt,ort,nzeros,timsSTA]=readperms(prename,PERMSTA,PERMROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
            % 2. this is for data with no response
            [opt, ort, timeperm] = readperm_nofilterv2(prename,PERMSTA,PERMROTS,idx,sps,fact);
        else
            fileflag = 0;   % change the file flag to 0, meaning abnormal
            fprintf('No data for station %s in day %s %s, this day will be omitted. \n',PERMSTA(idx,1:3), YEAR, JDAY);
            break   % break the entire station loop
        end
    end
    
    [LIA,idx]=ismember(stas(ista,:),POLSTA,'rows');
    if LIA        % if are in POLSTA
        found=found+LIA; %better be 1
        if year==2003 && jday<213        % should result from some criteria
            fact=7.5e-3;
        else
            fact=1.5e-3;
        end
        
        %%% readpols is an EXTERNAL FUNCTION
        % opt:
        % ort:
        % nzeros:
        fname = strcat(prename,'.',POLSTA(idx,1:4),'..HHE.D.SAC');
        if isfile(fname)    % if have the data file
            %                 % 1. this is for data without removing station response
            %                 [opt,ort,nzeros]=readpols(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
            % 2. this is for data with no response
            [opt, ort, timepola] = readpola_nofilterv2(prename,POLSTA,POLROTS,idx,sps,fact);
        else
            fileflag = 0;   % change the file flag to 0, meaning abnormal
            fprintf('No data for station %s in day %s / %s. \n',POLSTA(idx,1:4), YEAR, JDAY);
            break   % break the entire station loop
        end
    end
    %         found=found         % could be a benchmark
    %factr1(ista)=prctile(abs(opt),90); %Not normalized
    %factr2(ista)=factr1(ista)/factr1(1); %This is what is used; keeps 1st station unchanged but scales the others
    STAopt(ista,:)=opt/scaleseisms(ista);   % what is STAopt for??
    %STAort(ista,:)=ort;
end
if fileflag == 0    % means there are missing files
    fprintf('Day %s / %s will be omitted because of missing files. \n', YEAR, JDAY);
end




