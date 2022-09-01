%%% try to understand if the frequency band used in calculating the rotation para will bias the
%%% result, and see if it is necessary to use different rotation parameters in different occasions,
%%% say LF and HF cases.

%%
[rothf] = calc_rots_noresp('002',0.5,6.5,40,25,35,12,8,4);
% rothf = [0	0	14	-10;
%          6	75	18	76;
%          0	0	39	10;
%          0	0	53	0;
%          4	70	36	-14;
%          5	65	40	-36;
%          2	55	15	-15];

%%
[rotlf1] = calc_rots_noresp('002',0.5,1.25,20,25,35,12,8,4);
% rotlf1 = [0	0	17	-4;
%          4	75	18	38;
%          0	0	33	6;
%          2	105	56	0;
%          2	90	41	-6;
%          5	135	75	-14;
%          0	0	21	-7];
     
%%
[rotlf2] = calc_rots_noresp('002',0.5,1.25,40,50,70,20,15,4);
% rotlf2 = [0	0	17	-9;
%           9	80	21	77;
%           0	0	33	13;
%           4	105	56	0;
%           3	85	41	-11;
%           9	130	72	-28;
%           2	70	21	-14];
