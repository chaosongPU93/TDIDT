%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% it seems that the LFE detection catalog in Michael Bostock's newest catalog
% has several fams that are close in space, but we need to make sure if they
% are indeed the same fam or not. One way is to check the number and exact 
% components of each fam.
%
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2020/09/01
% Last modified date:   2020/09/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
close all
clear

% [timoffrot158,nLFE158,lfe158] = testnlfe('158');
% 
% [timoffrot234,nLFE234,lfe234] = testnlfe('234');
% 
% % it seems that they don't share any same LFE
% [lfesame,I158,I234] = intersect(lfe158(2:4,:),lfe234(2:4,:),'rows','stable');

catapath= ('/home/data2/chaosong/matlab/allan/BOSTOCK/');
cataname = ('total_mag_detect_0000_cull_NEW.txt');
catalog = load(strcat(catapath, cataname));
%%% find the families that have LFE detections
famall = unique(catalog(:, 1));

%%
nmin = 10000;
famobj = '002';
for i = 1: length(famall)
    if famall(i) < 10
        fam = strcat('00',num2str(famall(i)));
    elseif famall(i) >= 10 && famall(i) < 100
        fam = strcat('0',num2str(famall(i)));
    else
        fam = num2str(famall(i));
    end
    
    [~,nLFE(i),~] = testnlfe(fam);
    
    if nLFE(i) <= nmin
        nmin = nLFE(i);
        famobj = fam;
    end
    
end

%% choose the above found fam with least LFE detections 269
famobj = '269';
[~,nLFEobj,lfeobj] = testnlfe(famobj);
lfeclose = zeros(size(lfeobj));
for i = 1: nLFEobj
    tmp = catalog(catalog(:,2)==lfeobj(i,2) & catalog(:,3)==lfeobj(i,3) & catalog(:,1)~=lfeobj(i,1),...
                  :);
    tmp1 = setdiff(tmp, lfeobj(i,:), 'rows', 'stable');
    lfeclose(i,:) = tmp1(abs(tmp1(:,4)-lfeobj(i,4)) == min(abs(tmp1(:,4)-lfeobj(i,4))), :);
    
end

[minwid, ind] = min(abs(lfeclose(:,4)- lfeobj(:,4))); 
timediff = sort(abs(lfeclose(:,4)- lfeobj(:,4)));
unitimed = unique(timediff);
counts = zeros(length(unitimed),1);
for i = 1: length(unitimed)
    counts(i) = sum(timediff == unitimed(i));
end
[maxct,ind1] = max(counts);
maxtdif = unitimed(ind1);

figure
plot(abs(lfeclose(:,4)- lfeobj(:,4)));

figure
h=histogram(timediff,'binwidth',0.025);
xlim([0 40]);


%% choose fam 158 or fam 234
famobj = '158';
[~,nLFEobj,lfeobj] = testnlfe(famobj);
lfeclose = zeros(size(lfeobj));
for i = 1: nLFEobj
    tmp = catalog(catalog(:,2)==lfeobj(i,2) & catalog(:,3)==lfeobj(i,3) & catalog(:,1)~=lfeobj(i,1),...
                  :);
    tmp1 = setdiff(tmp, lfeobj(i,:), 'rows', 'stable');
    tmp2 = tmp1(abs(tmp1(:,4)-lfeobj(i,4)) == min(abs(tmp1(:,4)-lfeobj(i,4))), :);
    lfeclose(i,:) = tmp2(1,:);      % in case there are more than 1
end

[minwid, ind] = min(abs(lfeclose(:,4)- lfeobj(:,4))) 
timediff = sort(abs(lfeclose(:,4)- lfeobj(:,4)));
unitimed = unique(timediff);
counts = zeros(length(unitimed),1);
for i = 1: length(unitimed)
    counts(i) = sum(timediff == unitimed(i));
end
[maxct,ind] = max(counts);
maxtdif = unitimed(ind);

lfemaxtdif = lfeclose(abs(lfeclose(:,4)- lfeobj(:,4)) == maxtdif, :);
% n234 = sum(lfemaxtdif(:,1)==234);
n158 = sum(lfemaxtdif(:,1)==158);

figure
plot(abs(lfeclose(:,4)- lfeobj(:,4)));

figure
h=histogram(timediff,'binwidth',0.025);
xlim([0 40]);


%% choose a random fam 
% rng('default');
indsel = randi(length(famall));
if famall(indsel) < 10
    famobj = strcat('00',num2str(famall(indsel)));
elseif famall(indsel) >= 10 && famall(indsel) < 100
    famobj = strcat('0',num2str(famall(indsel)));
else
    famobj = num2str(famall(indsel));
end

[~,nLFEobj,lfeobj] = testnlfe(famobj);
lfeclose = zeros(size(lfeobj));
for i = 1: nLFEobj
    tmp = catalog(catalog(:,2)==lfeobj(i,2) & catalog(:,3)==lfeobj(i,3) & catalog(:,1)~=lfeobj(i,1),...
                  :);
    tmp1 = setdiff(tmp, lfeobj(i,:), 'rows', 'stable');
    tmp2 = tmp1(abs(tmp1(:,4)-lfeobj(i,4)) == min(abs(tmp1(:,4)-lfeobj(i,4))), :);
    lfeclose(i,:) = tmp2(1,:);      % in case there are more than 1
end

[minwid, ind] = min(abs(lfeclose(:,4)- lfeobj(:,4))) 
timediff = sort(abs(lfeclose(:,4)- lfeobj(:,4)));
unitimed = unique(timediff);
counts = zeros(length(unitimed),1);
for i = 1: length(unitimed)
    counts(i) = sum(timediff == unitimed(i));
end
[maxct,ind] = max(counts);
maxtdif = unitimed(ind);


% figure
% plot(abs(lfeclose(:,4)- lfeobj(:,4)));
% 
% figure
% h=histogram(timediff,'binwidth',0.025);
% xlim([0 40]);








