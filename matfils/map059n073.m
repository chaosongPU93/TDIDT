%Sits at one spot (rotations; delays) and looks for the cross-correlation.
clear all
close all
format short e

% IDENTIF4='map2003.269.59.-73.95.90.50_2-8-4s';
% IDENTIF128='map2003.269.59.-73.95.90.50_2-8-128s';
% ALL128='all.59.-73.128s.fort.13.sor';
% OLDER=0;
% tees(1,:)=[47500 71000];
% tees(2,:)=[47500 84000];
% tees(3,:)=[60200 60400];
% tees(4,:)=[60700 61300];
% tees(5,:)=[78400 79100];
% tees(6,:)=[0 86400];
% tees(7,:)=[47500 86400];
% ntees=7;
% % IDENTIF4='map2003.270.59.-73.95.90.50_2-8-4s';
% % IDENTIF128='map2003.270.59.-73.95.90.50_2-8-128s';
% % ALL128='all.59.-73.128s.fort.13.sor';
% % OLDER128='map2003.269.59.-73.95.90.50_2-8-128s';
% % OLDER=1;
% % tees(1,:)=[1400 1800];
% % tees(2,:)=[2000 3600];
% % tees(3,:)=[10000 11400];
% % tees(4,:)=[47200 51400];
% % tees(5,:)=[47200 48400];
% % tees(6,:)=[48800 50200];
% % tees(7,:)=[50500 51200];
% % tees(8,:)=[0 86400];
% % ntees=8;
% IDENTIF4='map2003.271.59.-73.95.90.50_2-8-4s';  %right to left
% IDENTIF128='map2003.271.59.-73.95.90.50_2-8-128s';
% ALL128='all.59.-73.128s.fort.13.sor';
% OLDER128='map2003.269-270.59.-73.95.90.50_2-8-128s';
% OLDER=1;
% tees(1,:)=[15550 22000];
% tees(2,:)=[22750 27915];
% tees(3,:)=[27915 30850];
% tees(4,:)=[30870 32030];
% tees(5,:)=[30870 32330];
% tees(6,:)=[32330 32681];
% tees(7,:)=[32681 32900];
% tees(8,:)=[32900 33333];
% tees(9,:)=[32330 32900]; 
% tees(10,:)=[32330 33333];
% tees(11,:)=[33333 36805];
% tees(12,:)=[36809 37187];
% tees(13,:)=[36809 39351];
% tees(14,:)=[39351 39900];
% tees(15,:)=[39351 46670];
% tees(16,:)=[46670 48920];
% tees(17,:)=[48920 51100];
% tees(18,:)=[51100 52070];
% tees(19,:)=[52070 62000];
% tees(20,:)=[54669 54850];
% tees(21,:)=[54669 54780];
% tees(22,:)=[54780 54850];
% tees(23,:)=[62000 66020];
% tees(24,:)=[66080 66880];
% tees(25,:)=[68990 86400];
% tees(26,:)=[0 86400];
% ntees=26;
% % IDENTIF4='map2003.272.59.-73.95.90.50_2-8-4s';
% % IDENTIF128='map2003.272.59.-73.95.90.50_2-8-128s';
% % ALL128='all.59.-73.128s.fort.13.sor';
% % OLDER128='map2003.269-271.59.-73.95.90.50_2-8-128s';
% % OLDER=1;
% % tees(1,:)=[100 4400];
% % tees(2,:)=[4600 7150];
% % tees(3,:)=[8000 17000];
% % tees(4,:)=[17700 20000];
% % tees(5,:)=[20000 24000];
% % tees(6,:)=[17700 24000];
% % tees(7,:)=[24000 29000];
% % tees(8,:)=[29000 40000];
% % tees(9,:)=[43000 48000];
% % tees(10,:)=[49000 61000];
% % tees(11,:)=[0 86400];
% % ntees=11;

% % IDENTIF4='map2005.245.59.-73.95.90.50_2-8-4s';
% % IDENTIF128='map2005.245.59.-73.95.90.50_2-8-128s';
% % ALL128='all.59.-73.128s.fort.13.sor';
% % OLDER=0;
% % tees(1,:)=[0 86400];
% % tees(2,:)=[66000 86400];
% % ntees=2;
% IDENTIF4='map2005.246.59.-73.95.90.50_2-8-4s';
% IDENTIF128='map2005.246.59.-73.95.90.50_2-8-128s';
% ALL128='all.59.-73.128s.fort.13.sor';
% OLDER=0;
% tees(1,:)=[3600 4000];
% %tees(2,:)=[0 15000]; %tees(3,:)=[16500 24000]; %tees(4,:)=[29000 34000]; %tees(9,:)=[65000 86400]; %nonsense!
% tees(2,:)=[15100 16500];
% tees(3,:)=[24400 29000];
% tees(4,:)=[34700 37600];
% tees(5,:)=[40000 40320];
% tees(6,:)=[44000 47000];
% tees(7,:)=[46000 47000];
% tees(8,:)=[62200 64500];
% tees(9,:)=[0 86400];
% ntees=9;
% % IDENTIF4='map2005.247.59.-73.95.90.50_2-8-4s';  %right to left
% % IDENTIF128='map2005.247.59.-73.95.90.50_2-8-128s';
% % ALL128='all.59.-73.128s.fort.13.sor';
% % OLDER128='map2005.246.59.-73.95.90.50_2-8-128s';
% % OLDER=1;
% % tees(1,:)=[16500 19200];
% % tees(2,:)=[23040 32000];
% % tees(3,:)=[32000 33000];
% % tees(4,:)=[34500 34900];
% % tees(5,:)=[35600 37885];
% % tees(6,:)=[37885 38230];
% % tees(7,:)=[38660 41210];
% % tees(8,:)=[38660 44000];
% % tees(9,:)=[23040 44000];
% % tees(10,:)=[44000 46500];
% % tees(11,:)=[46500 47440];
% % tees(12,:)=[47440 50170];
% % tees(13,:)=[50170 51060];
% % tees(14,:)=[51200 53000];
% % tees(15,:)=[53000 57000];
% % tees(16,:)=[55800 56500];
% % tees(17,:)=[51200 57000];
% % tees(18,:)=[57000 61000];
% % tees(19,:)=[59200 59950];
% % tees(20,:)=[61000 62500];
% % tees(21,:)=[61900 62350];
% % tees(22,:)=[62600 67704];
% % tees(23,:)=[68220 69120];
% % tees(24,:)=[69300 71800];
% % tees(25,:)=[72000 76100];
% % tees(26,:)=[69300 76100];
% % tees(27,:)=[76500 78210];
% % tees(28,:)=[78900 80100];
% % tees(29,:)=[80760 81200];
% % tees(30,:)=[84050 85430];
% % tees(31,:)=[85430 86100];
% % tees(32,:)=[84000 86100];
% % tees(33,:)=[0 86100];
% % %tees(13,:)=[51200 86400];
% % ntees=33;
% IDENTIF4='map2005.248.59.-73.95.90.50_2-8-4s';
% IDENTIF128='map2005.248.59.-73.95.90.50_2-8-128s';
% ALL128='all.59.-73.128s.fort.13.sor';
% OLDER128='map2005.246-247.59.-73.95.90.50_2-8-128s';
% OLDER=1;
% tees(1,:)=[0 4000];
% tees(2,:)=[4800 6140];
% tees(3,:)=[7020 8700];
% tees(4,:)=[9290 11050];
% tees(5,:)=[11290 12700];
% tees(6,:)=[13000 15440];
% tees(7,:)=[15440 16310];
% tees(8,:)=[16570 16995];
% tees(9,:)=[16995 17550];
% tees(10,:)=[17550 18570];
% % tees(10,:)=[17550 17820];
% % tees(11,:)=[17820 18660];
% tees(11,:)=[18570 19500];
% tees(12,:)=[19500 20040];
% tees(13,:)=[20150 21300];
% tees(14,:)=[21300 22500];
% tees(15,:)=[20150 22500];
% tees(16,:)=[22500 28000];
% %bad data gap
% tees(17,:)=[32420 34600];
% tees(18,:)=[34600 36300];
% tees(19,:)=[36300 38860];
% tees(20,:)=[38860 45700];
% tees(21,:)=[45700 46600];
% tees(22,:)=[46600 48400];
% tees(23,:)=[48400 49820];
% tees(24,:)=[49200 49820];
% tees(25,:)=[50050 54530];
% tees(26,:)=[49200 54530];
% tees(27,:)=[54890 57140];
% %tees(28,:)=[57720 58750];
% tees(28,:)=[57720 59400];
% tees(29,:)=[59400 72400];
% tees(30,:)=[72400 72920];
% tees(31,:)=[72920 74620];
% tees(32,:)=[74620 75380];
% tees(33,:)=[72920 75380];
% tees(34,:)=[75380 86000];
% tees(35,:)=[0 86000];
% ntees=35;

% IDENTIF4='map2006.119.59.-73.95.90.50_2-8-4s';  %Down from "NW"
% IDENTIF128='map2006.119.59.-73.95.90.50_2-8-128s';
% ALL128='all.59.-73.128s.fort.13.sor';
% OLDER=0;
% tees(1,:)=[0 9000];
% tees(2,:)=[9000 21800];
% tees(3,:)=[21800 28000];
% tees(4,:)=[28000 30500];
% tees(5,:)=[30500 34600];
% tees(6,:)=[34600 36000];
% tees(7,:)=[36000 36610];
% tees(8,:)=[36610 40200];
% tees(9,:)=[40200 42500];
% tees(10,:)=[42500 44000];
% tees(11,:)=[44000 44200];
% tees(12,:)=[44200 45250];
% tees(13,:)=[44000 45250];
% tees(14,:)=[45250 47500];
% tees(15,:)=[47500 49500];
% tees(16,:)=[49500 49750];
% tees(17,:)=[49750 51425];
% tees(18,:)=[51425 51931];
% tees(19,:)=[51931 52550];
% tees(20,:)=[52550 62000];
% tees(21,:)=[62000 64750];
% tees(22,:)=[64000 64550];
% tees(23,:)=[64550 71900];
% tees(24,:)=[71900 72300];
% tees(25,:)=[72300 75800];
% tees(26,:)=[75800 78500]; %77900
% tees(27,:)=[78500 86400];
% tees(28,:)=[0 86400];
% ntees=28;
% % IDENTIF4='map2006.120.59.-73.95.90.50_2-8-4s';
% % IDENTIF128='map2006.120.59.-73.95.90.50_2-8-128s';
% % ALL128='all.59.-73.128s.fort.13.sor';
% % OLDER128='map2006.119.59.-73.95.90.50_2-8-128s';
% % OLDER=1;
% % tees(1,:)=[0 4000];
% % tees(2,:)=[4000 10000];
% % tees(3,:)=[10000 11600];
% % tees(4,:)=[11600 13100];
% % tees(5,:)=[10000 13100];
% % tees(6,:)=[13100 16000];
% % tees(7,:)=[16000 30000];
% % tees(8,:)=[30000 31200];
% % tees(9,:)=[31200 43000];
% % tees(10,:)=[43000 47000];
% % tees(11,:)=[47000 47900];
% % tees(12,:)=[47000 55000];
% % tees(13,:)=[43000 55000];
% % tees(14,:)=[55000 86400];
% % tees(15,:)=[0 86400];
% % ntees=15;
% IDENTIF4='map2006.121.59.-73.95.90.50_2-8-4s';
% IDENTIF128='map2006.121.59.-73.95.90.50_2-8-128s';
% ALL128='all.59.-73.128s.fort.13.sor';
% OLDER128='map2006.119-120.59.-73.95.90.50_2-8-128s';
% OLDER=1;
% tees(1,:)=[0 6100];
% tees(2,:)=[6100 13000];
% tees(3,:)=[13000 20800];
% tees(4,:)=[22000 23100];
% tees(5,:)=[23100 25000];
% tees(6,:)=[25000 39000];
% tees(7,:)=[39000 54000];
% tees(8,:)=[54000 71500];
% tees(9,:)=[71500 73500];
% tees(10,:)=[73500 86400];
% tees(11,:)=[0 86400];
% ntees=11;

  %%Different direction(???), but amplitude too low to see everything that happens.
  %%Good one for slow and rapid progression through lower right region (day 72).
  %%Can see amp increase w/time in that region (day 72; 73?).
% IDENTIF4='map2005.072.59.-73.95.90.50_2-8-4s';
% IDENTIF128='map2005.072.59.-73.95.90.50_2-8-128s';
% ALL128='all.59.-73.128s.fort.13.sor';
% OLDER=0;
% tees(1,:)=[0 11000]; 
% tees(2,:)=[11000 20000]; 
% tees(3,:)=[20000 34300]; 
% tees(4,:)=[34300 36000]; 
% tees(5,:)=[36000 38000]; 
% tees(6,:)=[38000 40000]; 
% tees(7,:)=[40000 41000]; 
% tees(8,:)=[41000 44000]; 
% tees(9,:)=[45600 46400]; 
% tees(10,:)=[50200 51500]; 
% tees(11,:)=[51500 55500]; 
% tees(12,:)=[55500 58200]; 
% tees(13,:)=[58200 63000]; 
% tees(14,:)=[63000 64500]; 
% tees(15,:)=[71000 76000]; 
% tees(16,:)=[0 77000]; %77-82 garbage
% ntees=16;
% % IDENTIF4='map2005.073.59.-73.95.90.50_2-8-4s';
% % IDENTIF128='map2005.073.59.-73.95.90.50_2-8-128s';
% % ALL128='all.59.-73.128s.fort.13.sor';
% % OLDER128='map2005.072.59.-73.95.90.50_2-8-128s';
% % OLDER=1;
% % tees(1,:)=[7000 19000]; 
% % tees(2,:)=[19000 24100]; 
% % tees(3,:)=[19000 24100]; 
% % tees(4,:)=[30000 41000]; 
% % tees(5,:)=[41000 84000]; 
% % tees(6,:)=[0 86400]; 
% % ntees=6;
% IDENTIF4='map2005.075.59.-73.95.90.50_2-8-4s';
% IDENTIF128='map2005.075.59.-73.95.90.50_2-8-128s';
% ALL128='all.59.-73.128s.fort.13.sor';
% OLDER128='map2005.072-073.59.-73.95.90.50_2-8-128s';
% OLDER=1;
% tees(1,:)=[0 9000]; 
% tees(2,:)=[9000 15000]; 
% tees(3,:)=[15000 19000]; 
% tees(4,:)=[23000 24300]; 
% tees(5,:)=[25000 28500]; 
% tees(6,:)=[28500 31200]; 
% tees(7,:)=[31200 37000]; 
% tees(8,:)=[37000 42000]; 
% tees(9,:)=[42000 46100]; 
% tees(10,:)=[47000 50000]; 
% tees(11,:)=[52000 86000]; 
% tees(12,:)=[0 86400]; 
% ntees=12;
% % IDENTIF4='map2005.077.59.-73.95.90.50_2-8-4s';
% % IDENTIF128='map2005.077.59.-73.95.90.50_2-8-128s';
% % ALL128='all.59.-73.128s.fort.13.sor';
% % OLDER128='map2005.072-073-075.59.-73.95.90.50_2-8-128s';
% % OLDER=1;
% % tees(1,:)=[0 22000]; 
% % tees(2,:)=[22000 25000]; 
% % tees(3,:)=[25000 34000]; 
% % tees(4,:)=[29500 34000]; 
% % tees(5,:)=[29500 44500]; 
% % tees(6,:)=[44500 50500]; 
% % tees(7,:)=[50500 52500]; 
% % tees(8,:)=[53500 60000]; 
% % tees(9,:)=[61500 65500]; 
% % tees(10,:)=[0 86400]; 
% % ntees=10;
% IDENTIF4='map2005.078.59.-73.95.90.50_2-8-4s';  %Pretty much nada, except evidence of cycle-skipping
% IDENTIF128='map2005.078.59.-73.95.90.50_2-8-128s';
% ALL128='all.59.-73.128s.fort.13.sor';
% OLDER128='map2005.072-073-075-077.59.-73.95.90.50_2-8-128s';
% OLDER=1;
% tees(1,:)=[0 86400]; %77-82
% ntees=1;

IDENTIF4='map2004.191.59.-73.95.90.50_2-8-4s';
IDENTIF128='map2004.191.59.-73.95.90.50_2-8-128s';
ALL128='all.59.-73.128s.fort.13.sor';
OLDER=0;
tees(1,:)=[66000 77000];
tees(2,:)=[77000 82200];
tees(3,:)=[82200 85000];
tees(4,:)=[60000 86400];
ntees=4;
% % IDENTIF4='map2004.192.59.-73.95.90.50_2-8-4s';
% % IDENTIF128='map2004.192.59.-73.95.90.50_2-8-128s';
% % ALL128='all.59.-73.128s.fort.13.sor';
% % OLDER128='map2004.191.59.-73.95.90.50_2-8-128s';
% % OLDER=1;
% % tees(1,:)=[0 2000];
% % tees(2,:)=[2000 3000];
% % tees(3,:)=[3000 5500];
% % tees(4,:)=[6000 7000];
% % tees(5,:)=[7000 10000];
% % tees(6,:)=[11800 13000];
% % tees(7,:)=[16900 18400];
% % tees(8,:)=[18400 36000];
% % tees(9,:)=[36000 38000];
% % tees(10,:)=[38000 54000];
% % tees(11,:)=[54000 59000];
% % tees(12,:)=[59000 64000];
% % tees(13,:)=[64000 66000];
% % tees(14,:)=[66000 78000];
% % tees(15,:)=[78000 86400];
% % tees(16,:)=[0 86400];
% % ntees=16;
% IDENTIF4='map2004.193.59.-73.95.90.50_2-8-4s';
% IDENTIF128='map2004.193.59.-73.95.90.50_2-8-128s';
% ALL128='all.59.-73.128s.fort.13.sor';
% OLDER128='map2004.191-192.59.-73.95.90.50_2-8-128s';
% OLDER=1;
% tees(1,:)=[1000 4000];
% tees(2,:)=[5000 8000];
% tees(3,:)=[8000 14000];
% tees(4,:)=[14000 22000];
% tees(5,:)=[22000 22600];
% tees(6,:)=[22600 31000];
% tees(7,:)=[31000 32000];
% tees(8,:)=[32000 35000];
% tees(9,:)=[35000 37000];
% tees(10,:)=[37000 41200];
% tees(11,:)=[41200 42000];
% tees(12,:)=[42000 44500];
% tees(13,:)=[44500 50800];
% tees(14,:)=[50800 52500];
% tees(15,:)=[52500 55000];
% tees(16,:)=[55000 58300];
% tees(17,:)=[58300 73000];
% tees(18,:)=[73000 86400];
% tees(19,:)=[0 86400];
% ntees=19;
% % IDENTIF4='map2004.194.59.-73.95.90.50_2-8-4s';  %Nada, really
% % IDENTIF128='map2004.194.59.-73.95.90.50_2-8-128s';
% % ALL128='all.59.-73.128s.fort.13.sor';
% % OLDER128='map2004.191-193.59.-73.95.90.50_2-8-128s';
% % OLDER=1;
% % tees(1,:)=[1000 9000];
% % tees(2,:)=[9000 15000];
% % tees(3,:)=[15000 27000]; %Data gap next
% % tees(4,:)=[41000 51000];
% % tees(5,:)=[51000 64000];
% % tees(6,:)=[64000 72000];
% % tees(7,:)=[72000 86400];
% % tees(8,:)=[0 86400];
% % ntees=8;
% IDENTIF4='map2004.195.59.-73.95.90.50_2-8-4s';
% IDENTIF128='map2004.195.59.-73.95.90.50_2-8-128s';
% ALL128='all.59.-73.128s.fort.13.sor';
% OLDER128='map2004.191-194.59.-73.95.90.50_2-8-128s';
% OLDER=1;
% tees(1,:)=[20500 22500]; %low amplitude; window in activity elsewhere
% tees(2,:)=[0 42000];
% ntees=2;

  %%Looks like some cycle skipping here.  
  %%Comes down from the top.  Growth of amp pretty obvious.
% IDENTIF4='map2004.117.59.-73.95.90.50_2-8-4s';
% IDENTIF128='map2004.117.59.-73.95.90.50_2-8-128s';
% ALL128='all.59.-73.128s.fort.13.sor';
% OLDER=0;
% tees(1,:)=[12000 44000];
% tees(2,:)=[44000 60000];
% tees(3,:)=[60000 70500];
% tees(4,:)=[70900 72100];
% tees(5,:)=[73000 86400];
% tees(6,:)=[44000 86400];
% ntees=6;
% % IDENTIF4='map2004.118.59.-73.95.90.50_2-8-4s';
% % IDENTIF128='map2004.118.59.-73.95.90.50_2-8-128s';
% % ALL128='all.59.-73.128s.fort.13.sor';
% % OLDER128='map2004.117.59.-73.95.90.50_2-8-128s';
% % OLDER=1;
% % tees(1,:)=[1000 11500];
% % tees(2,:)=[11900 14000];
% % tees(3,:)=[16000 18000];
% % tees(4,:)=[19000 21500];
% % tees(5,:)=[21500 25000];
% % tees(6,:)=[25000 27000];
% % tees(7,:)=[27000 30000];
% % tees(8,:)=[31000 34000];
% % tees(9,:)=[34000 36000];
% % tees(10,:)=[36000 44000];
% % tees(11,:)=[47000 58000];
% % tees(12,:)=[58000 61000];
% % tees(13,:)=[61000 85000];
% % tees(14,:)=[0 86400];
% % ntees=14;

  %%Bad data in Day 1; looks lie it might start in NW.
  %%Maybe trend in amp on day 1, but harder to be confident we have front.
  %%Not obviously present on day 2.
% IDENTIF4='map2004.362.59.-73.95.90.50_2-8-4s';
% IDENTIF128='map2004.362.59.-73.95.90.50_2-8-128s';
% ALL128='all.59.-73.128s.fort.13.sor';
% OLDER=0;
% tees(1,:)=[65000 85000];
% tees(2,:)=[85000 86400];
% tees(3,:)=[65000 86400];
% ntees=3;
% % IDENTIF4='map2004.363.59.-73.95.90.50_2-8-4s';
% % IDENTIF128='map2004.363.59.-73.95.90.50_2-8-128s';
% % ALL128='all.59.-73.128s.fort.13.sor';
% % OLDER128='map2004.362.59.-73.95.90.50_2-8-128s';
% % OLDER=1;
% % tees(1,:)=[0 2000];
% % tees(2,:)=[2000 3600];
% % tees(3,:)=[4300 5400];
% % tees(4,:)=[5400 17000];
% % tees(5,:)=[17000 19500];
% % tees(6,:)=[19500 21500];
% % tees(7,:)=[21500 34500];
% % tees(8,:)=[34500 37000];
% % tees(9,:)=[40000 58000];
% % tees(10,:)=[58000 62500];
% % tees(11,:)=[62500 66000];
% % tees(12,:)=[71000 86400];
% % tees(13,:)=[0 86400];
% % ntees=13;
% IDENTIF4='map2004.364.59.-73.95.90.50_2-8-4s';
% IDENTIF128='map2004.364.59.-73.95.90.50_2-8-128s';
% ALL128='all.59.-73.128s.fort.13.sor';
% OLDER128='map2004.362-363.59.-73.95.90.50_2-8-128s';
% OLDER=1;
% tees(1,:)=[6500 11800];
% tees(2,:)=[11800 18000];
% tees(3,:)=[18000 24800];
% tees(4,:)=[0 24500];
% ntees=4;

locs4=load(IDENTIF4);
locs128=load(IDENTIF128);
locsAll=load(ALL128);
if OLDER==1
    olderlocs=load(OLDER128);
end
nin4=length(locs4);
nin128=length(locs128);

scrsz=get(0,'ScreenSize');
for k=ntees:-1:1
    if locs4(1,1) >= tees(k,1)-1
        istart4=1;
    end
    for i=2:nin4-1
        if locs4(i-1,1) < tees(k,1)-1 && locs4(i,1) >= tees(k,1)-1
            istart4=i;
        end
        if locs4(i+1,1) > tees(k,2)+1 && locs4(i,1) <= tees(k,2)+1
            iend4=i;
        end
    end
    if locs4(nin4,1) <= tees(k,2)-1
        iend4=nin4;
    end
    if locs128(1,1) >= tees(k,1)-1
        istart128=1;
    end
    for i=2:nin128-1
        if locs128(i-1,1) < tees(k,1)-1 && locs128(i,1) >= tees(k,1)-1
            istart128=i;
        end
        if locs128(i+1,1) > tees(k,2)+1 && locs128(i,1) <= tees(k,2)+1
            iend128=i;
        end
    end
    if locs128(nin128,1) <= tees(k,2)-1
        iend128=nin128;
    end

    h=figure('Position',[scrsz(3)/3.1 1 scrsz(3)/3.1 scrsz(4)]);
    subplot(2,1,1,'align')
    colormap(jet) %PGSI=x, PGSS=y
    plot(locsAll(:,1),locsAll(:,2),'ko','MarkerSize',2,'linewidth',0.2)
    hold on
    plot(locs128(1:istart128,2),locs128(1:istart128,3),'ro','MarkerSize',2,'linewidth',0.2)
    if OLDER==1
        plot(olderlocs(:,2),olderlocs(:,3),'ro','MarkerSize',2,'linewidth',0.2)
    end
    scatter(locs128(istart128:iend128,2),locs128(istart128:iend128,3),90,locs128(istart128:iend128,1),'linewidth',1)
    % First on top
    scatter3(locs4(istart4:iend4,2),locs4(istart4:iend4,3),86400-locs4(istart4:iend4,1),25,locs4(istart4:iend4,1),'filled')
    axis equal
    axis([-16 16 -16 16])
    xlabel('PGSI')
    ylabel('PGSS')
    caxis([tees(k,1) tees(k,2)])
    colorbar
    title(IDENTIF4)
    box on

    subplot(2,1,2,'align')
    colormap(jet) %PGSI=x, PGSS=y
    plot(locsAll(:,1),locsAll(:,2),'ko','MarkerSize',2,'linewidth',0.2)
    hold on
    plot(locs128(1:istart128,2),locs128(1:istart128,3),'ro','MarkerSize',2,'linewidth',0.2)
    if OLDER==1
        plot(olderlocs(:,2),olderlocs(:,3),'ro','MarkerSize',2,'linewidth',0.2)
    end
    scatter3(locs4(istart4:iend4,2),locs4(istart4:iend4,3),86400-locs4(istart4:iend4,1),40,log(locs4(istart4:iend4,7)),'filled')
    axis equal
    axis([-16 16 -16 16])
    xlabel('PGSI')
    ylabel('PGSS')
    caxis([-5 1])
    colorbar
    box on
    
    set(h,'PaperPosition',[0.25 0.25 8 10.5])
    print(h,'-depsc',['tmp',int2str(k),'.eps'])
    
end

%plot whole day (or last window)
k=ntees;
if locs4(1,1) >= tees(k,1)-1
    istart4=1;
end
for i=2:nin4-1
    if locs4(i-1,1) < tees(k,1)-1 && locs4(i,1) >= tees(k,1)-1
        istart4=i;
    end
    if locs4(i+1,1) > tees(k,2)+1 && locs4(i,1) <= tees(k,2)+1
        iend4=i;
    end
end
if locs4(nin4,1) <= tees(k,2)-1
    iend4=nin4;
end
if locs128(1,1) >= tees(k,1)-1
    istart128=1;
end
for i=2:nin128-1
    if locs128(i-1,1) < tees(k,1)-1 && locs128(i,1) >= tees(k,1)-1
        istart128=i;
    end
    if locs128(i+1,1) > tees(k,2)+1 && locs128(i,1) <= tees(k,2)+1
        iend128=i;
    end
end
if locs128(nin128,1) <= tees(k,2)-1
    iend128=nin128;
end
if locs128(1,1) >= tees(k,2)-1
    iend128=istart128;
end
if locs128(nin128,1) <= tees(k,1)-1
    istart128=iend128;
end

h=figure('Position',[1 1 scrsz(3)/3.1 scrsz(4)]);
subplot(2,1,1,'align')
colormap(jet) 
plot(locsAll(:,1),locsAll(:,2),'ko','MarkerSize',2,'linewidth',0.2)
hold on
plot(locs128(1:istart128,2),locs128(1:istart128,3),'ro','MarkerSize',2,'linewidth',0.2)
if OLDER==1
    plot(olderlocs(:,2),olderlocs(:,3),'ro','MarkerSize',2,'linewidth',0.2)
end
% Last on top
scatter3(locs4(istart4:iend4,2),locs4(istart4:iend4,3),86400-locs4(istart4:iend4,1),25,locs4(istart4:iend4,1),'filled')
axis equal
axis([-16 16 -16 16])
xlabel('PGSI')
ylabel('PGSS')
caxis([tees(k,1) tees(k,2)])
colorbar
title(IDENTIF4)
box on
subplot(2,1,2,'align')
colormap(jet) 
plot(locsAll(:,1),locsAll(:,2),'ko','MarkerSize',2,'linewidth',0.2)
hold on
plot(locs128(1:istart128,2),locs128(1:istart128,3),'ro','MarkerSize',2,'linewidth',0.2)
if OLDER==1
    plot(olderlocs(:,2),olderlocs(:,3),'ro','MarkerSize',2,'linewidth',0.2)
end
scatter3(locs4(istart4:iend4,2),locs4(istart4:iend4,3),86400-locs4(istart4:iend4,1),25,log(locs4(istart4:iend4,7)),'filled')
axis equal
axis([-16 16 -16 16])
xlabel('PGSI')
ylabel('PGSS')
caxis([-5 1])
colorbar
box on
set(h,'PaperPosition',[0.25 0.25 8 10.5])
print(gcf,'-depsc',['tmp',int2str(k+1),'.eps'],'-zbuffer')

%Reverse order of same
h=figure('Position',[2*scrsz(3)/3 1 scrsz(3)/3 scrsz(4)]);
subplot(2,1,1,'align')
colormap(jet) 
plot(locsAll(:,1),locsAll(:,2),'ko','MarkerSize',2,'linewidth',0.2)
hold on
plot(locs128(1:istart128,2),locs128(1:istart128,3),'ro','MarkerSize',2,'linewidth',0.2)
if OLDER==1
    plot(olderlocs(:,2),olderlocs(:,3),'ro','MarkerSize',2,'linewidth',0.2)
end
% Last on top
scatter3(locs4(istart4:iend4,2),locs4(istart4:iend4,3),locs4(istart4:iend4,1),25,locs4(istart4:iend4,1),'filled')
axis equal
axis([-16 16 -16 16])
xlabel('PGSI')
ylabel('PGSS')
caxis([tees(k,1) tees(k,2)])
colorbar
title(IDENTIF4)
box on
subplot(2,1,2,'align')
colormap(jet) 
plot(locsAll(:,1),locsAll(:,2),'ko','MarkerSize',2,'linewidth',0.2)
hold on
plot(locs128(1:istart128,2),locs128(1:istart128,3),'ro','MarkerSize',2,'linewidth',0.2)
if OLDER==1
    plot(olderlocs(:,2),olderlocs(:,3),'ro','MarkerSize',2,'linewidth',0.2)
end
scatter3(locs4(istart4:iend4,2),locs4(istart4:iend4,3),locs4(istart4:iend4,1),25,log(locs4(istart4:iend4,7)),'filled')
axis equal
axis([-16 16 -16 16])
xlabel('PGSI')
ylabel('PGSS')
caxis([-5 1])
colorbar
box on
set(h,'PaperPosition',[0.25 0.25 8 10.5])
print(gcf,'-depsc',['tmp',int2str(k+2),'.eps'],'-zbuffer')
%print(h,'-depsc', ['tmp',int2str(k+1),'.eps'])
%print(h,'-djpeg','-r600',['tmp',int2str(k+1),'.jpg'])
%print(h,'-dpdf', ['tmp',int2str(k+1),'.pdf'])
     
