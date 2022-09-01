 % This is Allan's script to generate the grid file 'xygridArmb' [off12 off13 dx dy]
 % This script interpolate John's original grid 'p085p020' that is evenly spaced in offset to a grid
 % that is evenly spaced in dx and dy, using the function 'griddata'.
 clear all
 close all
 scrsz=get(0,'ScreenSize');
 wid=scrsz(3);
 hite=scrsz(4);
 scrat=wid/hite; 
 
 rads=pi/180.;
 erad=6372.028;
 lat0=48.0+26.32/60.;
 lon0=123.0+35.07/60.;
 srad=erad*cos(lat0*rads)
 
 A=load('p085p020');
 dlat=A(:,3); %degrees
 mlat=A(:,4); %minutes
 lat=dlat+mlat/60.;
 dlon=A(:,5); %degrees
 mlon=A(:,6); %minutes
 lon=dlon+mlon/60.;
 dep=A(:,7);

 centPGSS=86 %5 %I think this 86 is correct ...
 centPGSI=20
 offPGSS=centPGSS-A(:,1); %PGSS is uniform in columns; changes across rows.
 offPGSI=centPGSI-A(:,2); %PGSI is uniform in rows; changes down columns.
% In p085p020, order is PGSS PGSI.  PGSI changes more rapidly.
% in output of wiggledaywig2, order is [timswin(n) xmaxPGSIntmp(n) xmaxPGSSntmp(n) ...
 dy=rads*(lat-lat0)*erad; %x and y coords in km
 dx=-rads*(lon-lon0)*srad; %(minus bcause 123 is -123)
%  fid = fopen('newgrid','w');
 ng=[offPGSS,offPGSI,dx,dy];
%  fprintf(fid,'%4i %4i %6.2f %6.2f\n',ng');
%  fclose(fid);

 mm=25; nn=25;  %original grid size; uniformly spaced station offsets
 dxr=reshape(dx,mm,nn); %degrees
 dyr=reshape(dy,mm,nn); %minutes
 offPGSSr=reshape(offPGSS,mm,nn);
 offPGSIr=reshape(offPGSI,mm,nn);
 xends=[-7 3];
 yends=[-4 4];
 %xends=[-5 5];
 %yends=[-5 5];
 [xq,yq]=meshgrid(xends(1):0.05:xends(2), yends(1):0.05:yends(2));
 oPGSS=griddata(dxr,dyr,offPGSSr,xq,yq);
 oPGSI=griddata(dxr,dyr,offPGSIr,xq,yq);

 h=figure('Position',[wid/6 1 2*wid/3 hite]); %center
 len=size(xq,1)*size(xq,2);
 dxvect=reshape(xq,len,1);
 dyvect=reshape(yq,len,1);
 oPGSSvect=reshape(oPGSS,len,1);
 oPGSIvect=reshape(oPGSI,len,1);
 nng=[oPGSSvect,oPGSIvect,dxvect,dyvect]; %nng (written to xygrid) is PGSS, PGSI, dx, dy
 maxPGSS=max(oPGSSvect); %These are for colorbar limits
 minPGSS=min(oPGSSvect);
 maxPGSI=max(oPGSIvect);
 minPGSI=min(oPGSIvect);
 maxSSSI=max(oPGSIvect-oPGSSvect); % SSIB-SILB = (PGC-SILB)-(PGC-SSIB) ?
 minSSSI=min(oPGSIvect-oPGSSvect);
 ranges=[maxPGSS-minPGSS maxPGSI-minPGSI maxSSSI-minSSSI];
 mids=[(maxPGSS+minPGSS)/2 (maxPGSI+minPGSI)/2 (maxSSSI+minSSSI)/2]
 [maxrange,imax]=max(ranges)
 colormins=mids-maxrange/2
 colormaxs=mids+maxrange/2
 subplot(3,2,2)
 scatter(nng(:,3),nng(:,4),10,nng(:,1),'filled')
 axis equal
 axis([xends(1) xends(2) yends(1) yends(2)])
 %caxis([minPGSS maxPGSS])
 caxis([colormins(1) colormaxs(1)])
 colorbar
 subplot(3,2,4)
 scatter(nng(:,3),nng(:,4),10,nng(:,2),'filled')
 axis equal
 axis([xends(1) xends(2) yends(1) yends(2)])
 %caxis([minPGSI maxPGSI])
 caxis([colormins(2) colormaxs(2)])
 colorbar
 subplot(3,2,6)
 scatter(nng(:,3),nng(:,4),10,nng(:,2)-nng(:,1),'filled')
 axis equal
 axis([xends(1) xends(2) yends(1) yends(2)])
 %caxis([minSSSI maxSSSI])
 caxis([colormins(3) colormaxs(3)])
 colorbar

 subplot(3,2,1)
 scatter(ng(:,3),ng(:,4),30,ng(:,1),'filled')
 axis equal
 axis([xends(1) xends(2) yends(1) yends(2)])
 %caxis([minPGSS maxPGSS])
 caxis([colormins(1) colormaxs(1)])
 colorbar
 ylabel('PGC-SSIB')
 subplot(3,2,3)
 scatter(ng(:,3),ng(:,4),30,ng(:,2),'filled')
 axis equal
 axis([xends(1) xends(2) yends(1) yends(2)])
 %caxis([minPGSI maxPGSI])
 caxis([colormins(2) colormaxs(2)])
 colorbar
 ylabel('PGC-SILB')
 subplot(3,2,5)
 scatter(ng(:,3),ng(:,4),30,ng(:,2)-ng(:,1),'filled')
 axis equal
 axis([xends(1) xends(2) yends(1) yends(2)])
 %caxis([minSSSI maxSSSI])
 caxis([colormins(3) colormaxs(3)])
 colorbar
 ylabel('SSIB-SILB')
 
%  fid = fopen('xygrid','w');
%  fprintf(fid,'%7.3f %7.3f %7.2f %6.2f\n',nng');
%  fclose(fid);

 nng(sqrt(nng(:,3).^2+nng(:,4).^2)>2,:)=[];
