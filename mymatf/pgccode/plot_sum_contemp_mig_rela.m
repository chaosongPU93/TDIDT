% function plot_sum_contemp_mig_rela
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to plot the summary of contemporaneous hf & lf 
% migrations using different window length, loopoffmax, ccmin in lf
% detections, while assuming hf uses the same parameters all the time. But
% technically it just reads the file specified.
% 
% Looks like figure 4b in Rubin&Armbruster (2013)
% 
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/09/19
% Last modified date:   2019/09/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
% close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

scrsz=get(0,'ScreenSize');

workpath = getenv('MHOME');
rstpath = strcat(workpath, '/Seisbasics/hypoinverse/forsummary');
datapath = strcat(getenv('ALLAN'),'/data-no-resp/PGCtrio');

% load files to find lat0, lat0
fam = '002';
winlenhf = 4;

%% plot the comparison with different lf detections
trange = [2004199,0.20e+4,0.37e+4;
          2005255,3.42e+4,3.65e+4;
          2005255,5.80e+4,5.96e+4;
          2005255,6.15e+4,6.26e+4;
          2005255,6.70e+4,6.90e+4;
          2005255,7.50e+4,7.60e+4;
          2005256,0.36e+4,0.55e+4;
          2005256,7.70e+4,8.55e+4;];           
meddif1 = zeros(length(trange), 4);   

for i = 1: length(trange)
% i=3;

    %%%% 1st combination, start
    %%% load original hf migration
    winlenlf1 = 16;
    lofflf1 = 4;
    ccminlf1 = 0.5;
    
    SUFFIXhf = strcat('hf.alldays.time.',num2str(winlenhf),'_',num2str(winlenlf1),'.',num2str(lofflf1),...
                      '.',num2str(ccminlf1));
    fname = strcat(rstpath, '/eventloc.',fam,'.',SUFFIXhf);
    hfmaplocall = load(fname);
    
    % convert absolute loc to relative loc
    [ind,~] = ind2sub(size(hfmaplocall), find(hfmaplocall(:,4)==0 & hfmaplocall(:,5)==0, 1, 'first'));
    lon0 = hfmaplocall(ind,1);
    lat0 = hfmaplocall(ind,2);
    if isempty(lon0)
        lon0 = -123.5850;
        lat0 = 48.4367;
    end
    
    hfrelalocall = hfmaplocall;
    [dx, dy] = absloc2relaloc(hfmaplocall(:,1),hfmaplocall(:,2),lon0,lat0);
    hfrelalocall(:,1) = dx;
    hfrelalocall(:,2) = dy;

    %%% load contemporal migration results and convert abs map loc to
    %%% rela loc
    fname = strcat(rstpath,'/eventloc.contemp_mig_hf.',fam,'_',num2str(trange(i,1)),...
                   '_',num2str(trange(i,2)),'-',num2str(trange(i,3)),'.',num2str(winlenhf),...
                   '_',num2str(winlenlf1),'.',num2str(lofflf1),'.',num2str(ccminlf1));
    conmigmaphf1 = load(fname);
    fname = strcat(rstpath,'/eventloc.contemp_mig_lf.',fam,'_',num2str(trange(i,1)),...
                   '_',num2str(trange(i,2)),'-',num2str(trange(i,3)),'.',num2str(winlenhf),...
                   '_',num2str(winlenlf1),'.',num2str(lofflf1),'.',num2str(ccminlf1));
    conmigmaplf1 = load(fname);
    
    conmigrelahf1 = conmigmaphf1;
    [dx, dy] = absloc2relaloc(conmigmaphf1(:,1),conmigmaphf1(:,2),lon0,lat0);
    conmigrelahf1(:,1) = dx;
    conmigrelahf1(:,2) = dy;
    conmigrelahf1(:,end) = conmigrelahf1(:,end)/3600;
    
    conmigrelalf1 = conmigmaplf1;
    [dx, dy] = absloc2relaloc(conmigmaplf1(:,1),conmigmaplf1(:,2),lon0,lat0);
    conmigrelalf1(:,1) = dx;
    conmigrelalf1(:,2) = dy;
    conmigrelalf1(:,end) = conmigrelalf1(:,end)/3600;
    
    %%% calculate the magnitude (length) of difference
    difx1 = conmigrelalf1(:,1)-conmigrelahf1(:,1);
    dify1 = conmigrelalf1(:,2)-conmigrelahf1(:,2);
    mag1 = sqrt( difx1.^2 + dify1.^2 );
    [mag1sort,ind1] = sort(mag1);
    difx1sort = difx1(ind1);
    dify1sort = dify1(ind1);
    wtdifx1 = zeros(length(mag1),1);
    wtdify1 = zeros(length(mag1),1);
    wtmag1 = zeros(length(mag1),1);
    for j = 1:length(mag1)
% j=2;
        fact = 50;
        magthr = mag1sort(j)+3/fact;
        wt1 = 1/2*erfc((mag1sort-magthr)*fact);
        wtdifx1(j) = sum(wt1.*difx1sort)/sum(wt1);
        wtdify1(j) = sum(wt1.*dify1sort)/sum(wt1);
        wtmag1(j) = sqrt(wtdifx1(j)^2+wtdify1(j)^2);
    end
    varwtdifx1 = wtdifx1(2:end)-wtdifx1(1:end-1);
    varwtdify1 = wtdify1(2:end)-wtdify1(1:end-1);
    varmag1 = sqrt(varwtdifx1.^2+ varwtdify1.^2);
    medmag1 = median(mag1);
    
    %%% remember to examine the results, compared to hf_lf_temopral
    hflfdif1 = [];
    hflfdif1(:,1:2) = conmigrelahf1(:,1:2) - conmigrelalf1(:,1:2);
    hflfdif1(:,3:4) = conmigrelahf1(:,4:5) - conmigrelalf1(:,4:5)/20*40;   % offset in lf is based on 20 sps
    meddif1(i,:) = median(hflfdif1);
    
    if(min(conmigmaplf1(:,10))< trange(i,2) || min(conmigmaphf1(:,10)) <trange(i,2) ...
            ||  max(conmigmaplf1(:,10)) >trange(i,3) || max(conmigmaphf1(:,10)) >trange(i,3))
        disp('wrong');
    else
        disp('right');
    end
    %%%% 1st combination, end
    
    %%%% 2nd combination, start
    winlenlf2 = 16;
    lofflf2 = 2;
    ccminlf2 = 0.55;
    SUFFIXhf = strcat('hf.alldays.time.',num2str(winlenhf),'_',num2str(winlenlf2),'.',num2str(lofflf2),...
                      '.',num2str(ccminlf2));
    fname = strcat(rstpath, '/eventloc.',fam,'.',SUFFIXhf);
    hfmaplocall = load(fname);    
    % convert absolute loc to relative loc
    [ind,~] = ind2sub(size(hfmaplocall), find(hfmaplocall(:,4)==0 & hfmaplocall(:,5)==0, 1, 'first'));
    lon0 = hfmaplocall(ind,1);
    lat0 = hfmaplocall(ind,2);
    if isempty(lon0)
        lon0 = -123.5850;
        lat0 = 48.4367;
    end
    hfrelalocall = hfmaplocall;
    [dx, dy] = absloc2relaloc(hfmaplocall(:,1),hfmaplocall(:,2),lon0,lat0);
    hfrelalocall(:,1) = dx;
    hfrelalocall(:,2) = dy;
    fname = strcat(rstpath,'/eventloc.contemp_mig_hf.',fam,'_',num2str(trange(i,1)),...
                   '_',num2str(trange(i,2)),'-',num2str(trange(i,3)),'.',num2str(winlenhf),...
                   '_',num2str(winlenlf2),'.',num2str(lofflf2),'.',num2str(ccminlf2));
    conmigmaphf2 = load(fname);
    fname = strcat(rstpath,'/eventloc.contemp_mig_lf.',fam,'_',num2str(trange(i,1)),...
                   '_',num2str(trange(i,2)),'-',num2str(trange(i,3)),'.',num2str(winlenhf),...
                   '_',num2str(winlenlf2),'.',num2str(lofflf2),'.',num2str(ccminlf2));
    conmigmaplf2 = load(fname);
    
    conmigrelahf2 = conmigmaphf2;
    [dx, dy] = absloc2relaloc(conmigmaphf2(:,1),conmigmaphf2(:,2),lon0,lat0);
    conmigrelahf2(:,1) = dx;
    conmigrelahf2(:,2) = dy;
    conmigrelahf2(:,end) = conmigrelahf2(:,end)/3600;
    
    conmigrelalf2 = conmigmaplf2;
    [dx, dy] = absloc2relaloc(conmigmaplf2(:,1),conmigmaplf2(:,2),lon0,lat0);
    conmigrelalf2(:,1) = dx;
    conmigrelalf2(:,2) = dy;
    conmigrelalf2(:,end) = conmigrelalf2(:,end)/3600;
    
    %%% calculate the magnitude (length) of difference
    difx2 = conmigrelalf2(:,1)-conmigrelahf2(:,1);
    dify2 = conmigrelalf2(:,2)-conmigrelahf2(:,2);
    mag2 = sqrt( difx2.^2 + dify2.^2 );
    [mag2sort,ind2] = sort(mag2);
    difx2sort = difx2(ind2);
    dify2sort = dify2(ind2);
    wtdifx2 = zeros(length(mag2),1);
    wtdify2 = zeros(length(mag2),1);
    wtmag2 = zeros(length(mag2),1);
    for j = 1:length(mag2)
        fact = 50;
        magthr = mag2sort(j)+3/fact;
        wt2 = 1/2*erfc((mag2sort-magthr)*fact);
        wtdifx2(j) = sum(wt2.*difx2sort)/sum(wt2);
        wtdify2(j) = sum(wt2.*dify2sort)/sum(wt2);
        wtmag2(j) = sqrt(wtdifx2(j)^2+wtdify2(j)^2);
    end
    varwtdifx2 = wtdifx2(2:end)-wtdifx2(1:end-1);
    varwtdify2 = wtdify2(2:end)-wtdify2(1:end-1);
    varmag2 = sqrt(varwtdifx2.^2+ varwtdify2.^2);
    medmag2 = median(mag2);
    
    %%% remember to examine the results, compared to hf_lf_temopral
    hflfdif2 = [];
    hflfdif2(:,1:2) = conmigrelahf2(:,1:2) - conmigrelalf2(:,1:2);
    hflfdif2(:,3:4) = conmigrelahf2(:,4:5) - conmigrelalf2(:,4:5)/20*40;   % offset in lf is based on 20 sps
    meddif2(i,:) = median(hflfdif2);
    
    if(min(conmigmaplf2(:,10))< trange(i,2) || min(conmigmaphf2(:,10)) <trange(i,2) ...
            ||  max(conmigmaplf2(:,10)) >trange(i,3) || max(conmigmaphf2(:,10)) >trange(i,3))
        disp('wrong');
    else
        disp('right');
    end
    %%%% 2nd combination, end
    
    %%%% 3rd combination, start
    winlenlf3 = 24;
    lofflf3 = 2;
    ccminlf3 = 0.45;
    SUFFIXhf = strcat('hf.alldays.time.',num2str(winlenhf),'_',num2str(winlenlf3),'.',num2str(lofflf3),...
                      '.',num2str(ccminlf3));
    fname = strcat(rstpath, '/eventloc.',fam,'.',SUFFIXhf);
    hfmaplocall = load(fname);
    
    % convert absolute loc to relative loc
    [ind,~] = ind2sub(size(hfmaplocall), find(hfmaplocall(:,4)==0 & hfmaplocall(:,5)==0, 1, 'first'));
    lon0 = hfmaplocall(ind,1);
    lat0 = hfmaplocall(ind,2);
    if isempty(lon0)
        lon0 = -123.5850;
        lat0 = 48.4367;
    end
    
    hfrelalocall = hfmaplocall;
    [dx, dy] = absloc2relaloc(hfmaplocall(:,1),hfmaplocall(:,2),lon0,lat0);
    hfrelalocall(:,1) = dx;
    hfrelalocall(:,2) = dy;
    
    fname = strcat(rstpath,'/eventloc.contemp_mig_hf.',fam,'_',num2str(trange(i,1)),...
                   '_',num2str(trange(i,2)),'-',num2str(trange(i,3)),'.',num2str(winlenhf),...
                   '_',num2str(winlenlf3),'.',num2str(lofflf3),'.',num2str(ccminlf3));
    conmigmaphf3 = load(fname);
    fname = strcat(rstpath,'/eventloc.contemp_mig_lf.',fam,'_',num2str(trange(i,1)),...
                   '_',num2str(trange(i,2)),'-',num2str(trange(i,3)),'.',num2str(winlenhf),...
                   '_',num2str(winlenlf3),'.',num2str(lofflf3),'.',num2str(ccminlf3));
    conmigmaplf3 = load(fname);
    
    conmigrelahf3 = conmigmaphf3;
    [dx, dy] = absloc2relaloc(conmigmaphf3(:,1),conmigmaphf3(:,2),lon0,lat0);
    conmigrelahf3(:,1) = dx;
    conmigrelahf3(:,2) = dy;
    conmigrelahf3(:,end) = conmigrelahf3(:,end)/3600;
    
    conmigrelalf3 = conmigmaplf3;
    [dx, dy] = absloc2relaloc(conmigmaplf3(:,1),conmigmaplf3(:,2),lon0,lat0);
    conmigrelalf3(:,1) = dx;
    conmigrelalf3(:,2) = dy;
    conmigrelalf3(:,end) = conmigrelalf3(:,end)/3600;
    
    %%% calculate the magnitude (length) of difference
    difx3 = conmigrelalf3(:,1)-conmigrelahf3(:,1);
    dify3 = conmigrelalf3(:,2)-conmigrelahf3(:,2);
    mag3 = sqrt( difx3.^2 + dify3.^2 );
    [mag3sort,ind3] = sort(mag3);
    difx3sort = difx3(ind3);
    dify3sort = dify3(ind3);
    wtdifx3 = zeros(length(mag3),1);
    wtdify3 = zeros(length(mag3),1);
    wtmag3 = zeros(length(mag3),1);
    for j = 1:length(mag3)
        fact = 50;
        magthr = mag3sort(j)+3/fact;
        wt3 = 1/2*erfc((mag3sort-magthr)*fact);
        wtdifx3(j) = sum(wt3.*difx3sort)/sum(wt3);
        wtdify3(j) = sum(wt3.*dify3sort)/sum(wt3);
        wtmag3(j) = sqrt(wtdifx3(j)^2+wtdify3(j)^2);
    end
    varwtdifx3 = wtdifx3(2:end)-wtdifx3(1:end-1);
    varwtdify3 = wtdify3(2:end)-wtdify3(1:end-1);
    varmag3 = sqrt(varwtdifx3.^2+ varwtdify3.^2);
    medmag3 = median(mag3);
    
    %%% remember to examine the results, compared to hf_lf_temopral
    hflfdif3 = [];
    hflfdif3(:,1:2) = conmigrelahf3(:,1:2) - conmigrelalf3(:,1:2);
    hflfdif3(:,3:4) = conmigrelahf3(:,4:5) - conmigrelalf3(:,4:5)/20*40;   % offset in lf is based on 20 sps
    meddif3(i,:) = median(hflfdif3);
    
    if(min(conmigmaplf3(:,10))< trange(i,2) || min(conmigmaphf3(:,10)) <trange(i,2) ...
            ||  max(conmigmaplf3(:,10)) >trange(i,3) || max(conmigmaphf3(:,10)) >trange(i,3))
        disp('wrong');
    else
        disp('right');
    end
    %%%% 3rd combination, end
    
    %%%% 4th combination, start
    winlenlf4 = 20;
    lofflf4 = 2;
    ccminlf4 = 0.5;
    SUFFIXhf = strcat('hf.alldays.time.',num2str(winlenhf),'_',num2str(winlenlf4),'.',num2str(lofflf4),...
                      '.',num2str(ccminlf4));
    fname = strcat(rstpath, '/eventloc.',fam,'.',SUFFIXhf);
    hfmaplocall = load(fname);
    
    % convert absolute loc to relative loc
    [ind,~] = ind2sub(size(hfmaplocall), find(hfmaplocall(:,4)==0 & hfmaplocall(:,5)==0, 1, 'first'));
    lon0 = hfmaplocall(ind,1);
    lat0 = hfmaplocall(ind,2);
    if isempty(lon0)
        lon0 = -123.5850;
        lat0 = 48.4367;
    end
    
    hfrelalocall = hfmaplocall;
    [dx, dy] = absloc2relaloc(hfmaplocall(:,1),hfmaplocall(:,2),lon0,lat0);
    hfrelalocall(:,1) = dx;
    hfrelalocall(:,2) = dy;
    
    fname = strcat(rstpath,'/eventloc.contemp_mig_hf.',fam,'_',num2str(trange(i,1)),...
                   '_',num2str(trange(i,2)),'-',num2str(trange(i,3)),'.',num2str(winlenhf),...
                   '_',num2str(winlenlf4),'.',num2str(lofflf4),'.',num2str(ccminlf4));
    conmigmaphf4 = load(fname);
    fname = strcat(rstpath,'/eventloc.contemp_mig_lf.',fam,'_',num2str(trange(i,1)),...
                   '_',num2str(trange(i,2)),'-',num2str(trange(i,3)),'.',num2str(winlenhf),...
                   '_',num2str(winlenlf4),'.',num2str(lofflf4),'.',num2str(ccminlf4));
    conmigmaplf4 = load(fname);
    
    conmigrelahf4 = conmigmaphf4;
    [dx, dy] = absloc2relaloc(conmigmaphf4(:,1),conmigmaphf4(:,2),lon0,lat0);
    conmigrelahf4(:,1) = dx;
    conmigrelahf4(:,2) = dy;
    conmigrelahf4(:,end) = conmigrelahf4(:,end)/3600;
    
    conmigrelalf4 = conmigmaplf4;
    [dx, dy] = absloc2relaloc(conmigmaplf4(:,1),conmigmaplf4(:,2),lon0,lat0);
    conmigrelalf4(:,1) = dx;
    conmigrelalf4(:,2) = dy;
    conmigrelalf4(:,end) = conmigrelalf4(:,end)/3600;
    
    %%% calculate the magnitude (length) of difference
    difx4 = conmigrelalf4(:,1)-conmigrelahf4(:,1);
    dify4 = conmigrelalf4(:,2)-conmigrelahf4(:,2);
    mag4 = sqrt( difx4.^2 + dify4.^2 );
    [mag4sort,ind4] = sort(mag4);
    difx4sort = difx4(ind4);
    dify4sort = dify4(ind4);
    wtdifx4 = zeros(length(mag4),1);
    wtdify4 = zeros(length(mag4),1);
    wtmag4 = zeros(length(mag4),1);
    for j = 1:length(mag4)
        fact = 50;
        magthr = mag4sort(j)+3/fact;
        wt4 = 1/2*erfc((mag4sort-magthr)*fact);
        wtdifx4(j) = sum(wt4.*difx4sort)/sum(wt4);
        wtdify4(j) = sum(wt4.*dify4sort)/sum(wt4);
        wtmag4(j) = sqrt(wtdifx4(j)^2+wtdify4(j)^2);
    end
    varwtdifx4 = wtdifx4(2:end)-wtdifx4(1:end-1);
    varwtdify4 = wtdify4(2:end)-wtdify4(1:end-1);
    varmag4 = sqrt(varwtdifx4.^2+ varwtdify4.^2);
    medmag4 = median(mag4);
    
    %%% remember to examine the results, compared to hf_lf_temopral
    hflfdif4 = [];
    hflfdif4(:,1:2) = conmigrelahf4(:,1:2) - conmigrelalf4(:,1:2);
    hflfdif4(:,3:4) = conmigrelahf4(:,4:5) - conmigrelalf4(:,4:5)/20*40;   % offset in lf is based on 20 sps
    meddif4(i,:) = median(hflfdif4);
    
    if(min(conmigmaplf4(:,10))< trange(i,2) || min(conmigmaphf4(:,10)) <trange(i,2) ...
            ||  max(conmigmaplf4(:,10)) >trange(i,3) || max(conmigmaphf4(:,10)) >trange(i,3))
        disp('wrong');
    else
        disp('right');
    end
    %%%% 4th combination, end
    
    %%%% 5th combination, start
    winlenlf5 = 12;
    lofflf5 = 2;
    ccminlf5 = 0.6;
    SUFFIXhf = strcat('hf.alldays.time.',num2str(winlenhf),'_',num2str(winlenlf5),'.',num2str(lofflf5),...
                      '.',num2str(ccminlf5));
    fname = strcat(rstpath, '/eventloc.',fam,'.',SUFFIXhf);
    hfmaplocall = load(fname);
    
    % convert absolute loc to relative loc
    [ind,~] = ind2sub(size(hfmaplocall), find(hfmaplocall(:,4)==0 & hfmaplocall(:,5)==0, 1, 'first'));
    lon0 = hfmaplocall(ind,1);
    lat0 = hfmaplocall(ind,2);
    if isempty(lon0)
        lon0 = -123.5850;
        lat0 = 48.4367;
    end
    
    hfrelalocall = hfmaplocall;
    [dx, dy] = absloc2relaloc(hfmaplocall(:,1),hfmaplocall(:,2),lon0,lat0);
    hfrelalocall(:,1) = dx;
    hfrelalocall(:,2) = dy;
    
    fname = strcat(rstpath,'/eventloc.contemp_mig_hf.',fam,'_',num2str(trange(i,1)),...
                   '_',num2str(trange(i,2)),'-',num2str(trange(i,3)),'.',num2str(winlenhf),...
                   '_',num2str(winlenlf5),'.',num2str(lofflf5),'.',num2str(ccminlf5));
    conmigmaphf5 = load(fname);
    fname = strcat(rstpath,'/eventloc.contemp_mig_lf.',fam,'_',num2str(trange(i,1)),...
                   '_',num2str(trange(i,2)),'-',num2str(trange(i,3)),'.',num2str(winlenhf),...
                   '_',num2str(winlenlf5),'.',num2str(lofflf5),'.',num2str(ccminlf5));
    conmigmaplf5 = load(fname);
    
    conmigrelahf5 = conmigmaphf5;
    [dx, dy] = absloc2relaloc(conmigmaphf5(:,1),conmigmaphf5(:,2),lon0,lat0);
    conmigrelahf5(:,1) = dx;
    conmigrelahf5(:,2) = dy;
    conmigrelahf5(:,end) = conmigrelahf5(:,end)/3600;
    
    conmigrelalf5 = conmigmaplf5;
    [dx, dy] = absloc2relaloc(conmigmaplf5(:,1),conmigmaplf5(:,2),lon0,lat0);
    conmigrelalf5(:,1) = dx;
    conmigrelalf5(:,2) = dy;
    conmigrelalf5(:,end) = conmigrelalf5(:,end)/3600;
    
    %%% calculate the magnitude (length) of difference
    difx5 = conmigrelalf5(:,1)-conmigrelahf5(:,1);
    dify5 = conmigrelalf5(:,2)-conmigrelahf5(:,2);
    mag5 = sqrt( difx5.^2 + dify5.^2 );
    [mag5sort,ind5] = sort(mag5);
    difx5sort = difx5(ind5);
    dify5sort = dify5(ind5);
    wtdifx5 = zeros(length(mag5),1);
    wtdify5 = zeros(length(mag5),1);
    wtmag5 = zeros(length(mag5),1);
    for j = 1:length(mag5)
        fact = 50;
        magthr = mag5sort(j)+3/fact;
        wt5 = 1/2*erfc((mag5sort-magthr)*fact);
        wtdifx5(j) = sum(wt5.*difx5sort)/sum(wt5);
        wtdify5(j) = sum(wt5.*dify5sort)/sum(wt5);
        wtmag5(j) = sqrt(wtdifx5(j)^2+wtdify5(j)^2);
    end
    varwtdifx5 = wtdifx5(2:end)-wtdifx5(1:end-1);
    varwtdify5 = wtdify5(2:end)-wtdify5(1:end-1);
    varmag5 = sqrt(varwtdifx5.^2+ varwtdify5.^2);
    medmag5 = median(mag5);
    
    %%% remember to examine the results, compared to hf_lf_temopral
    hflfdif5 = [];
    hflfdif5(:,1:2) = conmigrelahf5(:,1:2) - conmigrelalf5(:,1:2);
    hflfdif5(:,3:4) = conmigrelahf5(:,4:5) - conmigrelalf5(:,4:5)/20*40;   % offset in lf is based on 20 sps
    meddif5(i,:) = median(hflfdif5);
    
    if(min(conmigmaplf5(:,10))< trange(i,2) || min(conmigmaphf5(:,10)) <trange(i,2) ...
            ||  max(conmigmaplf5(:,10)) >trange(i,3) || max(conmigmaphf5(:,10)) >trange(i,3))
        disp('wrong');
    else
        disp('right');
    end
    %%%% 5th combination, enconmigmaplf2d
    
    %%% save all into a summary matrix
    detnum(1) = size(conmigrelahf1,1);
    detnum(2) = size(conmigrelahf2,1);
    detnum(3) = size(conmigrelahf3,1);
    detnum(4) = size(conmigrelahf4,1);
    detnum(5) = size(conmigrelahf5,1);
    
    maxcol = max([detnum(1),detnum(2),detnum(3),...
                 detnum(4),detnum(5)]);
    wtdifx = zeros(5,maxcol);
    wtdifx(1,1:detnum(1)) = wtdifx1;
    wtdifx(2,1:detnum(2)) = wtdifx2;
    wtdifx(3,1:detnum(3)) = wtdifx3;
    wtdifx(4,1:detnum(4)) = wtdifx4;
    wtdifx(5,1:detnum(5)) = wtdifx5;
    
    wtdify = zeros(5,maxcol);
    wtdify(1,1:detnum(1)) = wtdify1;
    wtdify(2,1:detnum(2)) = wtdify2;
    wtdify(3,1:detnum(3)) = wtdify3;
    wtdify(4,1:detnum(4)) = wtdify4;
    wtdify(5,1:detnum(5)) = wtdify5;    
    
    wtmag = zeros(5,maxcol);
    wtmag(1,1:detnum(1)) = wtmag1;
    wtmag(2,1:detnum(2)) = wtmag2;
    wtmag(3,1:detnum(3)) = wtmag3;
    wtmag(4,1:detnum(4)) = wtmag4;
    wtmag(5,1:detnum(5)) = wtmag5;
    
    varmag = zeros(5,maxcol-1);
    varmag(1,1:detnum(1)-1) = varmag1;
    varmag(2,1:detnum(2)-1) = varmag2;
    varmag(3,1:detnum(3)-1) = varmag3;
    varmag(4,1:detnum(4)-1) = varmag4;
    varmag(5,1:detnum(5)-1) = varmag5;    
    
    %%% define and position the figure frame and axes of each plot 
    f.fig=figure(i);
    set(f.fig,'Position',[scrsz(3)/10 scrsz(4)/10 3*scrsz(3)/4 4*scrsz(4)/5]);
    nrow = 2;
    ncol = 3;    
    for isub = 1:nrow*ncol
        f.ax(isub) = subplot(nrow,ncol,isub);
        % set some common features
        f.ax(isub).Box = 'on';
        grid(f.ax(isub), 'on');
        axis(f.ax(isub), 'equal');
        f.ax(isub).GridLineStyle = '--';
        xlabel(f.ax(isub),'E (km)','fontsize',12);
        ylabel(f.ax(isub),'N (km)','fontsize',12);
        f.ax(isub).XAxisLocation = 'top';
        xran = [-8 4];
        yran = [-4 6];
%         xran = [-10 10];
%         yran = [-10 10];
        xlen = xran(2)-xran(1);
        ylen = yran(2)-yran(1);
    end
    
    %%% reposition
    set(f.ax(1), 'position', [ (1-0.9)/ncol, 0.55, 1/ncol*0.8 0.4]);
    set(f.ax(2), 'position', [ (2-0.9)/ncol, 0.55, 1/ncol*0.8 0.4]);
    set(f.ax(3), 'position', [ (3-0.9)/ncol, 0.55, 1/ncol*0.8 0.4]);
    set(f.ax(4), 'position', [ (1-0.9)/ncol, 0.1, 1/ncol*0.8 0.4]);
    set(f.ax(5), 'position', [ (2-0.9)/ncol, 0.1, 1/ncol*0.8 0.4]);
    set(f.ax(6), 'position', [ (3-0.9)/ncol, 0.1, 1/ncol*0.8 0.4]);
    
    % subplot 1 of figure i
    hold(f.ax(1),'on');
    scatter(f.ax(1),conmigrelahf1(:,1),conmigrelahf1(:,2), 40, conmigrelahf1(:,end), 'filled','o', 'MarkerEdgeColor', 'k');   %, 'MarkerEdgeColor', 'w')
    scatter(f.ax(1),conmigrelalf1(:,1),conmigrelalf1(:,2), 60, conmigrelalf1(:,end), 'filled','s', 'MarkerEdgeColor', 'k');   %, 'MarkerEdgeColor', 'w')       
    colormap(f.ax(1),'jet');
    c=colorbar(f.ax(1),'SouthOutside');
    c.Position = [(1-0.9)/ncol+1/ncol*0.8/6, 0.55, 1/ncol*0.8*4/6, 0.02];
    caxis(f.ax(1),[trange(i,2)/3600 trange(i,3)/3600]);
    plot(f.ax(1),[-100 100],[0 0],'k--'); 
    plot(f.ax(1),[0 0],[-100 100],'k--'); 
    for j = 1: length(conmigrelahf1(:,1))
        x = [conmigrelahf1(j,1) conmigrelalf1(j,1)];
        y = [conmigrelahf1(j,2) conmigrelalf1(j,2)];        
        drawArrow(f.ax(1),x,y,xran,yran,'color',[140/255 140/255 140/255]);    
    end
    x = [0 -meddif1(i,1)];
    y = [0 -meddif1(i,2)];
    drawArrow(f.ax(1),x,y,xran,yran, 'color',[50/255 50/255 50/255],'linewidth',1.5);   % draw the median difference from hf to lf 
    text(f.ax(1),xran(1)+0.35*xlen, yran(1)+0.9*ylen,strcat(num2str(winlenhf),'s-',num2str(winlenlf1),'s-',num2str(lofflf1),'-',num2str(ccminlf1)),'fontsize',14);
    text(f.ax(1),xran(1)+0.4*xlen, yran(1)+0.8*ylen, strcat({'med vector= ('},sprintf('%.2f',-meddif1(i,1)),{', '},sprintf('%.2f',-meddif1(i,2)),')'),'fontsize',12);
    drawArrow(f.ax(1),[0 0],[0,medmag1],xran,yran,'linewidth',2);
    text(f.ax(1),xran(1)+0.4*xlen, yran(1)+0.7*ylen, strcat({'med magitude= '},sprintf('%.2f',medmag1)),'fontsize',12);
    hold(f.ax(1),'off');
    
    % subplot 2 of figure i
    hold(f.ax(2),'on');
    scatter(f.ax(2),conmigrelahf2(:,1),conmigrelahf2(:,2), 40, conmigrelahf2(:,end), 'filled','o', 'MarkerEdgeColor', 'k');   %, 'MarkerEdgeColor', 'w')
    scatter(f.ax(2),conmigrelalf2(:,1),conmigrelalf2(:,2), 60, conmigrelalf2(:,end), 'filled','s', 'MarkerEdgeColor', 'k');   %, 'MarkerEdgeColor', 'w')       
    colormap(f.ax(2),'jet');
    c=colorbar(f.ax(2),'SouthOutside');
    c.Position = [(2-0.9)/ncol+1/ncol*0.8/6, 0.55, 1/ncol*0.8*4/6, 0.02];
    caxis(f.ax(2),[trange(i,2)/3600 trange(i,3)/3600]);
    plot(f.ax(2),[-100 100],[0 0],'k--'); 
    plot(f.ax(2),[0 0],[-100 100],'k--'); 
    for j = 1: length(conmigrelahf2(:,1))
        x = [conmigrelahf2(j,1) conmigrelalf2(j,1)];
        y = [conmigrelahf2(j,2) conmigrelalf2(j,2)];        
        drawArrow(f.ax(2),x,y,xran,yran,'color',[140/255 140/255 140/255]);    
    end
    x = [0 -meddif2(i,1)];
    y = [0 -meddif2(i,2)];
    drawArrow(f.ax(2),x,y,xran,yran, 'color',[50/255 50/255 50/255],'linewidth',1.5);   % draw the median difference from hf to lf 
    text(f.ax(2),xran(1)+0.35*xlen, yran(1)+0.9*ylen,strcat(num2str(winlenhf),'s-',num2str(winlenlf2),'s-',num2str(lofflf2),'-',num2str(ccminlf2)),'fontsize',14);
    text(f.ax(2),xran(1)+0.4*xlen, yran(1)+0.8*ylen, strcat({'med vector= ('},sprintf('%.2f',-meddif2(i,1)),{', '},sprintf('%.2f',-meddif2(i,2)),')'),'fontsize',12);
    drawArrow(f.ax(2),[0 0],[0,medmag2],xran,yran,'linewidth',2);
    text(f.ax(2),xran(1)+0.4*xlen, yran(1)+0.7*ylen, strcat({'med magitude= '},sprintf('%.2f',medmag2)),'fontsize',12);
    hold(f.ax(2),'off');
    
    % subplot 3 of figure i
    hold(f.ax(3),'on');
    scatter(f.ax(3),conmigrelahf3(:,1),conmigrelahf3(:,2), 40, conmigrelahf3(:,end), 'filled','o', 'MarkerEdgeColor', 'k');   %, 'MarkerEdgeColor', 'w')
    scatter(f.ax(3),conmigrelalf3(:,1),conmigrelalf3(:,2), 60, conmigrelalf3(:,end), 'filled','s', 'MarkerEdgeColor', 'k');   %, 'MarkerEdgeColor', 'w')       
    colormap(f.ax(3),'jet');
    c=colorbar(f.ax(3),'SouthOutside');
    c.Position = [(3-0.9)/ncol+1/ncol*0.8/6, 0.55, 1/ncol*0.8*4/6, 0.02];
    caxis(f.ax(3),[trange(i,2)/3600 trange(i,3)/3600]);
    plot(f.ax(3),[-100 100],[0 0],'k--'); 
    plot(f.ax(3),[0 0],[-100 100],'k--'); 
    for j = 1: length(conmigrelahf3(:,1))
        x = [conmigrelahf3(j,1) conmigrelalf3(j,1)];
        y = [conmigrelahf3(j,2) conmigrelalf3(j,2)];        
        drawArrow(f.ax(3),x,y,xran,yran,'color',[140/255 140/255 140/255]);    
    end
    x = [0 -meddif3(i,1)];
    y = [0 -meddif3(i,2)];
    drawArrow(f.ax(3),x,y,xran,yran, 'color',[50/255 50/255 50/255],'linewidth',1.5);   % draw the median difference from hf to lf 
    text(f.ax(3),xran(1)+0.35*xlen, yran(1)+0.9*ylen,strcat(num2str(winlenhf),'s-',num2str(winlenlf3),'s-',num2str(lofflf3),'-',num2str(ccminlf3)),'fontsize',14);
    text(f.ax(3),xran(1)+0.4*xlen, yran(1)+0.8*ylen, strcat({'med vector= ('},sprintf('%.2f',-meddif3(i,1)),{', '},sprintf('%.2f',-meddif3(i,2)),')'),'fontsize',12);
    drawArrow(f.ax(3),[0 0],[0,medmag3],xran,yran,'linewidth',2);
    text(f.ax(3),xran(1)+0.4*xlen, yran(1)+0.7*ylen, strcat({'med magitude= '},sprintf('%.2f',medmag3)),'fontsize',12);
    hold(f.ax(3),'off');
    
    % subplot 4 of figure i
    hold(f.ax(4),'on');
    scatter(f.ax(4),conmigrelahf4(:,1),conmigrelahf4(:,2), 40, conmigrelahf4(:,end), 'filled','o', 'MarkerEdgeColor', 'k');   %, 'MarkerEdgeColor', 'w')
    scatter(f.ax(4),conmigrelalf4(:,1),conmigrelalf4(:,2), 60, conmigrelalf4(:,end), 'filled','s', 'MarkerEdgeColor', 'k');   %, 'MarkerEdgeColor', 'w')       
    colormap(f.ax(4),'jet');
    c=colorbar(f.ax(4),'SouthOutside');
    c.Position = [(1-0.9)/ncol+1/ncol*0.8/6, 0.1, 1/ncol*0.8*4/6, 0.02];
    caxis(f.ax(4),[trange(i,2)/3600 trange(i,3)/3600]);
    plot(f.ax(4),[-100 100],[0 0],'k--'); 
    plot(f.ax(4),[0 0],[-100 100],'k--'); 
    for j = 1: length(conmigrelahf4(:,1))
        x = [conmigrelahf4(j,1) conmigrelalf4(j,1)];
        y = [conmigrelahf4(j,2) conmigrelalf4(j,2)];        
        drawArrow(f.ax(4),x,y,xran,yran,'color',[140/255 140/255 140/255]);    
    end
    x = [0 -meddif4(i,1)];
    y = [0 -meddif4(i,2)];
    drawArrow(f.ax(4),x,y,xran,yran, 'color',[50/255 50/255 50/255],'linewidth',1.5);   % draw the median difference from hf to lf 
    text(f.ax(4),xran(1)+0.35*xlen, yran(1)+0.9*ylen,strcat(num2str(winlenhf),'s-',num2str(winlenlf4),'s-',num2str(lofflf4),'-',num2str(ccminlf4)),'fontsize',14);
    text(f.ax(4),xran(1)+0.4*xlen, yran(1)+0.8*ylen, strcat({'med vector= ('},sprintf('%.2f',-meddif4(i,1)),{', '},sprintf('%.2f',-meddif4(i,2)),')'),'fontsize',12);
    drawArrow(f.ax(4),[0 0],[0,medmag4],xran,yran,'linewidth',2);
    text(f.ax(4),xran(1)+0.4*xlen, yran(1)+0.7*ylen, strcat({'med magitude= '},sprintf('%.2f',medmag4)),'fontsize',12);
    hold(f.ax(4),'off');
    
    % subplot 5 of figure i
    hold(f.ax(5),'on');
    scatter(f.ax(5),conmigrelahf5(:,1),conmigrelahf5(:,2), 40, conmigrelahf5(:,end), 'filled','o', 'MarkerEdgeColor', 'k');   %, 'MarkerEdgeColor', 'w')
    scatter(f.ax(5),conmigrelalf5(:,1),conmigrelalf5(:,2), 60, conmigrelalf5(:,end), 'filled','s', 'MarkerEdgeColor', 'k');   %, 'MarkerEdgeColor', 'w')       
    colormap(f.ax(5),'jet');
    c=colorbar(f.ax(5),'SouthOutside');
    c.Position = [(2-0.9)/ncol+1/ncol*0.8/6, 0.1, 1/ncol*0.8*4/6, 0.02];
    caxis(f.ax(5),[trange(i,2)/3600 trange(i,3)/3600]);
    plot(f.ax(5),[-100 100],[0 0],'k--'); 
    plot(f.ax(5),[0 0],[-100 100],'k--'); 
    for j = 1: length(conmigrelahf5(:,1))
        x = [conmigrelahf5(j,1) conmigrelalf5(j,1)];
        y = [conmigrelahf5(j,2) conmigrelalf5(j,2)];        
        drawArrow(f.ax(5),x,y,xran,yran,'color',[140/255 140/255 140/255]);    
    end
    x = [0 -meddif5(i,1)];
    y = [0 -meddif5(i,2)];
    drawArrow(f.ax(5),x,y,xran,yran, 'color',[50/255 50/255 50/255],'linewidth',1.5);   % draw the median difference from hf to lf 
    text(f.ax(5),xran(1)+0.35*xlen, yran(1)+0.9*ylen,strcat(num2str(winlenhf),'s-',num2str(winlenlf5),'s-',num2str(lofflf5),'-',num2str(ccminlf5)),'fontsize',14);
    text(f.ax(5),xran(1)+0.4*xlen, yran(1)+0.8*ylen, strcat({'med vector= ('},sprintf('%.2f',-meddif5(i,1)),{', '},sprintf('%.2f',-meddif5(i,2)),')'),'fontsize',12);
    drawArrow(f.ax(5),[0 0],[0,medmag5],xran,yran,'linewidth',2);
    text(f.ax(5),xran(1)+0.4*xlen, yran(1)+0.7*ylen, strcat({'med magitude= '},sprintf('%.2f',medmag5)),'fontsize',12);
    hold(f.ax(5),'off');
    
    
     % subplot 6 of figure i 
    hold(f.ax(6),'on');
    ind = find(hfrelalocall(:,8)==trange(i,1) & hfrelalocall(:,10)>=trange(i,2) & hfrelalocall(:,10)<=trange(i,3));
    mighf = sortrows(hfrelalocall(ind,:),-10);
    scatter(f.ax(6),mighf(:,1),mighf(:,2), 40, mighf(:,end)/3600, 'filled','o');  %, 'MarkerEdgeColor', 'w')
    colormap(f.ax(6),'jet');
    c=colorbar(f.ax(6),'SouthOutside');
    c.Position = [(3-0.9)/ncol+1/ncol*0.8/6, 0.1, 1/ncol*0.8*4/6, 0.02];
    datestr = num2str(trange(i,1));
    yr = datestr(1:4);
    day = datestr(5:end);
    c.Label.String = strcat({'Hours on '},yr,{' '},day);
    c.Label.FontSize = 12;
    caxis(f.ax(6),[trange(i,2)/3600 trange(i,3)/3600])
    plot(f.ax(6),[-100 100],[0 0],'k--');
    plot(f.ax(6),[0 0],[-100 100],'k--');
    xlim(f.ax(6),xran);
    ylim(f.ax(6),yran);
    xlabel(f.ax(6),'E (km)','fontsize',12);
    ylabel(f.ax(6),'N (km)','fontsize',12);   
    hold(f.ax(6),'off');
    
    
    %%% plot the vector summation variation
    f2.fig = figure(i+length(trange));
    set(f2.fig,'Position',[scrsz(3)/10 scrsz(4)/10 3*scrsz(3)/4 3*scrsz(4)/5]);
    nrow = 3;
    ncol = 5;    
    for isub = 1:nrow*ncol
        f2.ax(isub) = subplot(nrow,ncol,isub);
        % set some common features
        f2.ax(isub).Box = 'on';
        grid(f2.ax(isub), 'on');
        f2.ax(isub).GridLineStyle = '--';
    end
    
    %%% reposition
    set(f2.ax(1), 'position', [ (1-0.9)/ncol, 0.55, 1/ncol*0.8 0.4]);
    set(f2.ax(2), 'position', [ (2-0.9)/ncol, 0.55, 1/ncol*0.8 0.4]);
    set(f2.ax(3), 'position', [ (3-0.9)/ncol, 0.55, 1/ncol*0.8 0.4]);
    set(f2.ax(4), 'position', [ (4-0.9)/ncol, 0.55, 1/ncol*0.8 0.4]);
    set(f2.ax(5), 'position', [ (5-0.9)/ncol, 0.55, 1/ncol*0.8 0.4]);
    
    set(f2.ax(6), 'position', [ (1-0.9)/ncol, 0.3, 1/ncol*0.8 0.2]);
    set(f2.ax(7), 'position', [ (2-0.9)/ncol, 0.3, 1/ncol*0.8 0.2]);
    set(f2.ax(8), 'position', [ (3-0.9)/ncol, 0.3, 1/ncol*0.8 0.2]);
    set(f2.ax(9), 'position', [ (4-0.9)/ncol, 0.3, 1/ncol*0.8 0.2]);
    set(f2.ax(10), 'position', [ (5-0.9)/ncol, 0.3, 1/ncol*0.8 0.2]);
    
    set(f2.ax(11), 'position', [ (1-0.9)/ncol, 0.05, 1/ncol*0.8 0.2]);
    set(f2.ax(12), 'position', [ (2-0.9)/ncol, 0.05, 1/ncol*0.8 0.2]);
    set(f2.ax(13), 'position', [ (3-0.9)/ncol, 0.05, 1/ncol*0.8 0.2]);
    set(f2.ax(14), 'position', [ (4-0.9)/ncol, 0.05, 1/ncol*0.8 0.2]);
    set(f2.ax(15), 'position', [ (5-0.9)/ncol, 0.05, 1/ncol*0.8 0.2]);
    


    for j = 1:5
%         j=1
        hold(f2.ax(j),'on');
        axis(f2.ax(j),'equal');
        alpha(f2.ax(j),0.5);
        xran = [-2 2];
        yran = [-2 2];
        plot(f2.ax(j),[-100 100],[0 0],'k--');
        plot(f2.ax(j),[0 0],[-100 100],'k--');
        inc = 200/detnum(j);
        for k = 1: detnum(j)
%             k=1
            x = [0 wtdifx(j,k)];
            y = [0 wtdify(j,k)];          
            drawArrow(f2.ax(j),x,y,xran,yran,'color',[inc*k/255 inc*k/255 inc*k/255]);
        end
        text(f2.ax(j),-1,-1.5,'weighted sum');
        hold(f2.ax(j),'off');
    end
    
    for j = 1:5
        hold(f2.ax(j+5),'on');
        plot(f2.ax(j+5),wtmag(j,1:detnum(j)),'k-');
        xran = [0 detnum(j)+1];
        yran = [min(wtmag(j,:)) max(wtmag(j,:))];
        xlen = xran(2)-xran(1);
        ylen = yran(2)-yran(1);
        xlim(f2.ax(j+5),xran);
        ylim(f2.ax(j+5),yran);
        text(f2.ax(j+5),xran(1)+xlen/4, yran(2)-ylen/5,'mag var');
        hold(f2.ax(j+5),'off');
    end
    
    for j = 1:5
        hold(f2.ax(j+10),'on');
        plot(f2.ax(j+10),varmag(j,1:detnum(j)-1),'k-');
        xran = [0 detnum(j)];
        yran = [min(varmag(j,:)) max(varmag(j,:))];
        xlen = xran(2)-xran(1);
        ylen = yran(2)-yran(1);
        xlim(f2.ax(j+10),xran);
        ylim(f2.ax(j+10),yran);
        text(f2.ax(j+10),xran(1)+xlen/4, yran(2)-ylen/5,'mag var of diff vector');
        hold(f2.ax(j+10),'off');
    end

end




































