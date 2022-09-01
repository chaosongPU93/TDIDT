%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is the script to plot the shear wave splitting parameter calculated
% by Chao (anyone)
%
%   
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/09/02
% Last modified date:   2019/09/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all
set(0,'DefaultFigureVisible','on');
scrsz=get(0,'ScreenSize');

% set path
workpath = getenv('ALLAN');
temppath = strcat(workpath, '/templates/');
if ~exist(temppath, 'dir')
    mkdir(temppath);
end
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
split = load(strcat(datapath, '/split_chao/lzbtrio/','rot_para_11fam0.5-6.5_25_35LZB_40sps7stas'));
% split = load(strcat(datapath, '/split_chao/','rot_para_11fam0.5-6.5_20_30LZB_40sps7stas'));


dtdir = ('/home/data2/chaosong/matlab/allan/2004/JUL/');
flist= ('datalist');
outf = ('staloc.txt');
stainfo = getstainfo(dtdir, flist, outf);
ind=[3,11,4];
stainfo = stainfo(ind,:);

bosdir = ('/home/data2/chaosong/matlab/allan/BOSTOCK');
lfefnm = ('newlfeloc');
lfeloc = load(fullfile(bosdir, lfefnm));    
%%% lfeloc is loctions of all fams, but not necessarily each fam has detections

nfampool = ['002';
            '043';
            '141';
            '047';
            '010';
            '144';
            '099';
            '068';
            '125';
            '147';
            '017'];   
nfam = length(nfampool);

[~,~,ind] = intersect(str2num(nfampool(:,:)), lfeloc(:, 1),'stable');

[scrsz, res] = pixelperinch(1);
f.fig=figure;
f.fig.Renderer='Painters';
widin = 8;  % maximum width allowed is 8.5 inches
htin = 2.5;   % maximum height allowed is 11 inches
set(f.fig,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);

nrow = 3;
ncol = 1;
for isub = 1:nrow*ncol
    f.ax(isub) = subplot(nrow,ncol,isub);
end

% reposition
set(f.ax(1), 'position', [ 0.1, 0.1, 0.25, 0.9]);
set(f.ax(2), 'position', [ 0.4, 0.1, 0.25, 0.9]);
set(f.ax(3), 'position', [ 0.7, 0.1, 0.25, 0.9]);

staind = [4,5,6];
for ii=1:3
    hold(f.ax(ii),'on');
    %%% plot the location of new fam pool
    scatter(f.ax(ii),lfeloc(ind,3), lfeloc(ind,2), 15,'filled','MarkerEdgeColor',[0.2 0.2 0.2],...
            'MarkerFaceColor',[0.2 0.2 0.2]);
    % text(lfeloc(ind,3),  lfeloc(ind,2), num2str(lfeloc(ind,1)), 'fontsize',12,'color','m');
    refsamp = 4;   % 0.1s in 40 sps
    for i = 1: nfam
        irow = (i-1)*7+staind(ii);
        angle = split(irow,2);
        offset = split(irow,1);
        x0 = lfeloc(ind(i),3);
        y0 = lfeloc(ind(i),2);
        reflen = 0.03;
        
        a = 0+1i*reflen*offset/refsamp;  % complex on positive y axis, reference vector points to N
        b = a*exp(1i * deg2rad(angle));   % rotate counter-clockwise by angle to the fast direction (axis)
        dx = real(b);
        dy = imag(b);
        
        x = [x0-0.5*dx; x0+0.5*dx];
        y = [y0-0.5*dy; y0+0.5*dy];
        
        plot(f.ax(ii),x,y,'color',[0.6 0.6 0.6],'linewidth',1.5);
        text(f.ax(ii),x0-0.02,y0+0.02,nfampool(i,:),'fontsize',7);
    end
    plot(f.ax(ii),[-123.95 -123.95+reflen],[48.65 48.65],'color',[0.6 0.6 0.6],'linewidth',1.5);
    text(f.ax(ii),-123.965,48.63,'0.1 s','fontsize',8);
    for jj=1:3
        if jj == 1
            plot(f.ax(ii),str2num(stainfo(jj,3)),str2num(stainfo(jj,2)),'^',...
                'MarkerSize',5,'MarkerEdgeColor','k');
        else
            plot(f.ax(ii),str2num(stainfo(jj,3)),str2num(stainfo(jj,2)),'s',...
                'MarkerSize',5,'MarkerEdgeColor','k');
        end
    end
    text(f.ax(ii),str2num(stainfo(1,3))-0.02,str2num(stainfo(1,2))+0.025,stainfo(1,1),...
        'fontsize',8);
    text(f.ax(ii),str2num(stainfo(2,3))-0.03,str2num(stainfo(2,2))+0.025,stainfo(2,1),...
        'fontsize',8);
    text(f.ax(ii),str2num(stainfo(3,3))+0.01,str2num(stainfo(3,2))+0.01,stainfo(3,1),...
        'fontsize',8);
    if ii == 1
        plot(f.ax(ii),str2num(stainfo(ii,3)),str2num(stainfo(ii,2)),'^','MarkerFaceColor','k',...
            'MarkerSize',5,'MarkerEdgeColor', 'k');
    else
        plot(f.ax(ii),str2num(stainfo(ii,3)),str2num(stainfo(ii,2)),'s','MarkerFaceColor','k',...
            'MarkerSize',5,'MarkerEdgeColor', 'k');
    end
    f.ax(ii).Box='on';
    grid(f.ax(ii),'on');
    axis(f.ax(ii),'equal');
    axis(f.ax(ii),[-124 -123.55 48.3 48.7]);
    hold(f.ax(ii),'off');
end

% f.ax(2).YAxis.Visible = 'off';
% f.ax(3).YAxis.Visible = 'off';
xlabel(f.ax(1),'Longitude (^o)');
ylabel(f.ax(1),'Latitude (^o)');
xlabel(f.ax(2),'Longitude (^o)');
xlabel(f.ax(3),'Longitude (^o)');

print(f.fig,'-depsc2',strcat('/home/data2/chaosong/CurrentResearch/Song_Rubin_2020/Figures/splitpara.eps'));


% %%% how to plot a circle
% th = 0:pi/50:2*pi;
% xunit = reflen/6 * cos(th) + x0;
% yunit = reflen/6 * sin(th) + y0;
% plot(xunit, yunit, 'k','linewidth',1);
% 
% xunit = reflen/3 * cos(th) + x0;
% yunit = reflen/3 * sin(th) + y0;
% plot(xunit, yunit, 'k','linewidth',1.5);