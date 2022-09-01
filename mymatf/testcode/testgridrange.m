
grid = load('/home/data2/chaosong/Seisbasics/hypoinverse/lzbrst/slab1.0.grid');

loccont = [-123.492667 48.451500 38.1400; 
           -123.772167 48.493000 35.5900; 
           -123.863167 48.528167 35.2100;
           -123.603333 48.440167 36.7100;
           -123.800167 48.408833 34.5200;
           -123.893333 48.536500 35.0700;
           -123.864500 48.498667 34.8800;
           -123.753333 48.525667 36.2000;
           -123.703667 48.502667 36.4100;
           -123.814333 48.538667 35.7900;
           -123.838500 48.544833 35.6600];
       
[dx, dy] = absloc2relaloc(-grid(:,1),grid(:,2),loccont(2,1),loccont(2,2));
grid = [dx dy grid];     % now increases to 23 cols

relacont = loccont;
[dx, dy] = absloc2relaloc(loccont(:,1),loccont(:,2),loccont(2,1),loccont(2,2));
relacont(:,1) = dx;
relacont(:,2) = dy;

e = -30:0.1:50;
n = -30:0.1:40;
[egrd, ngrd] = meshgrid(e,n);
depgrd = griddata(grid(:,1),grid(:,2),grid(:,5),egrd,ngrd, 'cubic');

figure
hold on
scatter(grid(:,1),grid(:,2),300,grid(:,5),'filled','s');
scatter(relacont(:,1),relacont(:,2),20,'ko','LineWidth',1);
axis([-30 50 -30 40]);
colormap jet
colorbar
caxis([32 40]);

workpath = getenv('MHOME');
rstpath = strcat(workpath, '/Seisbasics/hypoinverse/lzbrst');

print('-dpdf',strcat(rstpath,'/','slabstructurescatter.pdf'));


f.fig = figure;
f.fig.Renderer = 'Painters';
hold on
contourf(egrd,ngrd,depgrd,30:2:50,'k-',...
         'showtext','on','LabelSpacing',500); %,'linec','k'
% contour(egrd,ngrd,depgrd,32:2:40,'w--',...
%         'linew',1.5);
scatter(relacont(:,1),relacont(:,2),20,'ko','LineWidth',1);
axis([-30 50 -30 40]);
colormap jet
colorbar
caxis([30 50]);
print(f.fig,'-dpdf',strcat(rstpath,'/','slabstructurecoutour.pdf'));



