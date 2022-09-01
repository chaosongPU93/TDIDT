set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);

workpath = getenv('MHOME');
rstpath = strcat(workpath, '/Seisbasics/hypoinverse/lzbrst');

%%
%%% assume the real case is that HF and LF should locate at the same place, and set a number of
%%% points to get the uncorrected HF and LF points, vector from HF to LF would represent the
%%% filterring effect, and i want to see the random vector in time domain which can be put anywhere
%%% would also be transformed into space domain freely as well, i.e., the orientation of arrows
%%% should be same
a = -3:2:5;
b = -3:2:5;
[aa,bb] = meshgrid(a,b);
aaa = reshape(aa,[],1);
bbb = reshape(bb,[],1);

ccc = aaa+(1);
ddd = bbb+(-7);

sps = 40;
hfarr = [aaa/sps bbb/sps];
lfarr = [ccc/sps ddd/sps];

fid = fopen(fullfile(rstpath,'offset_filtercortest_hf'),'w+');
fprintf(fid,'%.4f %.4f \n',hfarr');
fclose(fid);    

fid = fopen(fullfile(rstpath,'offset_filtercortest_lf'),'w+');
fprintf(fid,'%.4f %.4f \n',lfarr');
fclose(fid);    

fam = '144';
locfile = strcat('eventloc.',fam,'.filtercortest_hf');
loccont = load(fullfile(rstpath,locfile));

%%% convert lfe location to km relative to 043
hf = loccont;
[dx, dy] = absloc2relaloc(loccont(:,1),loccont(:,2),-123.8933, 48.5365);
hf(:,1) = dx;
hf(:,2) = dy;

locfile = strcat('eventloc.',fam,'.filtercortest_lf');
loccont = load(fullfile(rstpath,locfile));

%%% convert lfe location to km relative to 043
lf = loccont;
[dx, dy] = absloc2relaloc(loccont(:,1),loccont(:,2),-123.8933, 48.5365);
lf(:,1) = dx;
lf(:,2) = dy;


f.fig=figure;
widin = 8;  % maximum width allowed is 8.5 inches
htin = 4;   % maximum height allowed is 11 inches
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);

nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f.ax(1), 'position', [ 0.1, 0.1, 0.35, 0.75]);
set(f.ax(2), 'position', [ 0.55, 0.1, 0.35, 0.75]);


ax = f.ax(1);
hold(ax,'on');

plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
% plot(ax,offcont(:,1)*sps,offcont(:,2)*sps,'b-','linew',1);
scatter(ax,aaa,bbb, 30, [0.5 0.5 0.5], 'filled','o');
scatter(ax,ccc,ddd, 10, [0.1 0.1 0.1], 'filled','o');
xran = [-10 10];
yran = [-10 10];
for i = 1: length(aaa)
    xvect = [aaa(i) ccc(i)];
    yvect = [bbb(i) ddd(i)];
    drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1,'color','b');
end
% strlbl = ['A';'B';'C';'D';'E';'F';'G';'H';'I']; 
% text(ax,offcont(:,1)*sps,offcont(:,2)*sps+2,strlbl(:,:),'fontsize',12);
text(ax,0.05,0.05,strcat({'fam: '},fam),'FontSize',11,'unit','normalized');

ax.Box='on';
grid(ax,'on');
ax.GridLineStyle = '--';
axis(ax,'equal');
axis(ax,[-10 10 -10 10]);
xlabel(ax,'off12 (sample) 40 sps');
ylabel(ax,'off13 (sample) 40 sps');
hold(ax,'off');


ax = f.ax(2);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
% plot(ax,relacont(:,1),relacont(:,2),'b-','linew',1);
scatter(ax,hf(:,1),hf(:,2), 30, [0.5 0.5 0.5], 'filled','o');
scatter(ax,lf(:,1),lf(:,2), 10, [0.1 0.1 0.1], 'filled','o');
xran = [-10 10];
yran = [-10 10];
for i = 1: length(aaa)
    xvect = [hf(i,1) lf(i,1)];
    yvect = [hf(i,2) lf(i,2)];
    drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1,'color','b');
end
% strlbl = ['A';'B';'C';'D';'E';'F';'G';'H';'I'];
% text(ax,relacont(:,1),relacont(:,2)+1,strlbl(:,:),'fontsize',12);
text(ax,0.05,0.05,strcat({'fam: '},fam),'FontSize',11,'unit','normalized');

ax.Box='on';
grid(ax,'on');
ax.GridLineStyle = '--';
axis(ax,'equal');
axis(ax,[-10 10 -10 10]);
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
hold(ax,'off');


%%
%%% assume a number of points which are uncorrected detections in LF (or, HF, doesn't matter tho),
%%% now i want to know if using filtering correction from different fams would make a difference.
%%% Of course it will, but i want to make sure if the difference vector in space domain would differ
%%% wrt the point location in time domain.  
a = -10:1:10;
b = -10:1:10;
[aa,bb] = meshgrid(a,b);
aaa = reshape(aa,[],1);
bbb = reshape(bb,[],1);

ccc = aaa-(-3/2);    % fam 006
ddd = bbb-(-7/2);

eee = aaa-(-2/2);    % fam 099
fff = bbb-(-9/2);

sps = 40;
arr006 = [ccc/sps ddd/sps];
arr099 = [eee/sps fff/sps];

fid = fopen(fullfile(rstpath,'offset_filtercortest_006'),'w+');
fprintf(fid,'%.4f %.4f \n',arr006');
fclose(fid);    

fid = fopen(fullfile(rstpath,'offset_filtercortest_099'),'w+');
fprintf(fid,'%.4f %.4f \n',arr099');
fclose(fid);    

fam = '006';
locfile = strcat('eventloc.',fam,'.filtercortest_006');
loccont = load(fullfile(rstpath,locfile));

%%% convert lfe location to km relative to 043
lf006 = loccont;
[dx, dy] = absloc2relaloc(loccont(:,1),loccont(:,2),-123.908000, 48.494167);
lf006(:,1) = dx;
lf006(:,2) = dy;

locfile = strcat('eventloc.',fam,'.filtercortest_099');
loccont = load(fullfile(rstpath,locfile));

%%% convert lfe location to km relative to 043
lf099 = loccont;
[dx, dy] = absloc2relaloc(loccont(:,1),loccont(:,2),-123.908000, 48.494167);
lf099(:,1) = dx;
lf099(:,2) = dy;


f.fig=figure;
widin = 12;  % maximum width allowed is 8.5 inches
htin = 7;   % maximum height allowed is 11 inches
set(f.fig,'Position',[1*scrsz(3)/20 scrsz(4)/10 widin*res htin*res]);

nrow = 1;
ncol = 2;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

%%% reposition
set(f.ax(1), 'position', [ 0.1, 0.1, 0.35, 0.75]);
set(f.ax(2), 'position', [ 0.55, 0.1, 0.35, 0.75]);


ax = f.ax(1);
hold(ax,'on');

plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
% plot(ax,offcont(:,1)*sps,offcont(:,2)*sps,'b-','linew',1);
scatter(ax,eee,fff, 30, [0.5 0.5 0.5], 'filled','o');   % 099
scatter(ax,ccc,ddd, 10, [0.1 0.1 0.1], 'filled','o');   % 006
xran = [-15 15];
yran = [-10 20];
for i = 1: length(aaa)
    xvect = [eee(i) ccc(i)];
    yvect = [fff(i) ddd(i)];
    drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1,'color','b');
end
% strlbl = ['A';'B';'C';'D';'E';'F';'G';'H';'I']; 
% text(ax,offcont(:,1)*sps,offcont(:,2)*sps+2,strlbl(:,:),'fontsize',12);
text(ax,0.05,0.05,strcat({'fam: '},fam),'FontSize',11,'unit','normalized');

ax.Box='on';
grid(ax,'on');
ax.GridLineStyle = '--';
axis(ax,'equal');
axis(ax,[-15 15 -10 20]);
xlabel(ax,'off12 (sample) 40 sps');
ylabel(ax,'off13 (sample) 40 sps');
hold(ax,'off');


ax = f.ax(2);
hold(ax,'on');
plot(ax,[-100 100],[0 0],'k--');
plot(ax,[0 0],[-100 100],'k--');
% plot(ax,relacont(:,1),relacont(:,2),'b-','linew',1);
scatter(ax,lf099(:,1),lf099(:,2), 30, [0.5 0.5 0.5], 'filled','o');
scatter(ax,lf006(:,1),lf006(:,2), 10, [0.1 0.1 0.1], 'filled','o');
xran = [-15 21];
yran = [-23 13];
for i = 1: length(aaa)
    xvect = [lf099(i,1) lf006(i,1)];
    yvect = [lf099(i,2) lf006(i,2)];
    drawArrow(ax,xvect,yvect,xran,yran,'linewidth',1,'color','b');
end
% strlbl = ['A';'B';'C';'D';'E';'F';'G';'H';'I'];
% text(ax,relacont(:,1),relacont(:,2)+1,strlbl(:,:),'fontsize',12);
text(ax,0.05,0.05,strcat({'fam: '},fam),'FontSize',11,'unit','normalized');

ax.Box='on';
grid(ax,'on');
ax.GridLineStyle = '--';
axis(ax,'equal');
axis(ax,[-15 21 -23 13]);
xlabel(ax,'E (km)');
ylabel(ax,'N (km)');
hold(ax,'off');


















