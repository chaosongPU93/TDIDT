  



clear all
close all
 
 wr1(:,1)=[1;2;3];

 
 

 ntrio=1:1
 


 wr=wr1(:,ntrio);


sta=['URSH';'MASH';'MGWH'];


%nsh=120;

nsamp=100;






loc=[];


     %loc=load(['chat10s_kii_',num2str(wr(1)),num2str(wr(2)),num2str(wr(3))]);
    loc=load('chat10s_kii_123.txt');
   %    loc=load('chat10s_kii123_john.txt');
    
loc(:,15)=[];

z=0;






    for j=-200:2:200
        for k=-500:2:500
            z=z+1;
            loc(z,15:16)=[j,k];
        end
    end


load('/u/yajun/Desktop/tryout/Japan/Kii/staloc.mat')

staloc(:,3:4)=[(staloc(:,1)-34.5)/180*pi*6372.028,...
    (staloc(:,2)-136.5)/180*pi*cos(34.5/180*pi)*6372.028];


loc(:,17:18)=[(loc(:,4)+loc(:,5)/60-34.5)/180*pi*6372.028,...
    (loc(:,6)+loc(:,7)/60-136.5)/180*pi*cos(34.5/180*pi)*6372.028];


figure;scatter(loc(:,18),loc(:,17),15,loc(:,13)*nsamp,'filled')
axis equal
colorbar
caxis([0 10])
axis([-100 100 -100 100])
hold on
plot(staloc(:,4),staloc(:,3),'r^','markersize',10)


figure;plot(loc(:,18),loc(:,17),'.')
hold on
text(loc(:,18),loc(:,17),num2str(loc(:,16)))
hold on
plot(staloc(:,4),staloc(:,3),'r^','markersize',10,'markerfacecolor','r')
axis([-100 100 -100 100])
axis equal



xx=min(loc(:,15)):2:max(loc(:,15));
yy=min(loc(:,16)):2:max(loc(:,16));

loc1=zeros(size(xx,2),size(yy,2));
loc2=zeros(size(xx,2),size(yy,2));

%loc3=zeros(size(xx,2),size(yy,2));

for i=1:size(loc,1)
    if loc1((loc(i,15)-min(loc(:,15)))/2+1,(loc(i,16)-min(loc(:,16)))/2+1)==0 && loc(i,13)*nsamp<=1
        loc1((loc(i,15)-min(loc(:,15)))/2+1,(loc(i,16)-min(loc(:,16)))/2+1)=loc(i,4)+loc(i,5)/60;
        loc2((loc(i,15)-min(loc(:,15)))/2+1,(loc(i,16)-min(loc(:,16)))/2+1)=loc(i,6)+loc(i,7)/60;
       % loc3((loc(i,15)-min(loc(:,15)))/2+1,(loc(i,16)-min(loc(:,16)))/2+1)=loc(i,8);
    end
end

z=0;locnn=[];
for i=1:size(xx,2)
    for j=1:size(yy,2)
        if loc1(i,j)~=0
        z=z+1;
        %locnn(z,:)=[xx(i),yy(j),loc1(i,j),loc2(i,j),loc3(i,j)];
        
        locnn(z,:)=[xx(i),yy(j),loc1(i,j),loc2(i,j)];
        end
    end
end
        
figure;
plot(-loc(:,6)-loc(:,7)/60,loc(:,4)+loc(:,5)/60,'r.')
hold on
plot(-locnn(:,4),locnn(:,3),'.')


% figure;
% scatter(-locnn(:,4),locnn(:,3),15,locnn(:,5),'filled')
% colorbar


[xq,yq]=meshgrid(xx,yy);

zq=griddata(locnn(:,1),locnn(:,2),locnn(:,3),xq,yq);
zq1=griddata(locnn(:,1),locnn(:,2),locnn(:,4),xq,yq);
%zq2=griddata(locnn(:,1),locnn(:,2),locnn(:,5),xq,yq);



locnx=[];z=0;
for i=1:size(xq,1)
    for j=1:size(xq,2)
        z=z+1;
        %locnx(z,:)=[xq(i,j),yq(i,j),zq(i,j),zq1(i,j),zq2(i,j)];
        
        locnx(z,:)=[xq(i,j),yq(i,j),zq(i,j),zq1(i,j)];
    end
end

figure;
plot(-loc(:,6)-loc(:,7)/60,loc(:,4)+loc(:,5)/60,'r.')
hold on
plot(-locnx(:,4),locnx(:,3),'+')










xx=min(loc(:,15)):2:max(loc(:,15));
yy=min(loc(:,16)):2:max(loc(:,16));

loc1=zeros(size(yy,2),size(xx,2));
loc2=zeros(size(yy,2),size(xx,2));
loc3=zeros(size(yy,2),size(xx,2));

for i=1:size(loc,1)
    if loc1((loc(i,16)-min(loc(:,16)))/2+1,(loc(i,15)-min(loc(:,15)))/2+1)==0 %&& loc(i,13)*40<=2
        loc1((loc(i,16)-min(loc(:,16)))/2+1,(loc(i,15)-min(loc(:,15)))/2+1)=loc(i,4)+loc(i,5)/60;
        loc2((loc(i,16)-min(loc(:,16)))/2+1,(loc(i,15)-min(loc(:,15)))/2+1)=loc(i,6)+loc(i,7)/60;
        loc3((loc(i,16)-min(loc(:,16)))/2+1,(loc(i,15)-min(loc(:,15)))/2+1)=loc(i,8);
    end
end



zqn=zeros(size(yy,2),size(xx,2));zqn1=zeros(size(yy,2),size(xx,2));zqn2=zeros(size(yy,2),size(xx,2));

zqn(zqn==0)=NaN;zqn1(zqn1==0)=NaN;zqn2(zqn2==0)=NaN;

for i=1:size(xq,1)
    for j=1:size(xq,2)
        if zq(i,j)>8 || zq(i,j)<=8 
            if loc1(i,j)~=0
            zqn(i,j)=loc1(i,j);
            zqn1(i,j)=loc2(i,j);
            zqn2(i,j)=loc3(i,j);
            end
        end
    end
end



locnx=[];z=0;
for i=1:size(xq,1)
    for j=1:size(xq,2)
        z=z+1;
        locnx(z,:)=[xq(i,j),yq(i,j),zqn(i,j),zqn1(i,j),zqn2(i,j)];
    end
end

% figure;
% plot(-loc(:,6)-loc(:,7)/60,loc(:,4)+loc(:,5)/60,'r.')
% hold on
% plot(-locnx(:,4),locnx(:,3),'+')


figure;
plot(-locnn(:,4),locnn(:,3),'r.')
hold on
plot(-locnx(:,4),locnx(:,3),'+')


figure;scatter(locnx(:,4),locnx(:,3),15,locnx(:,2),'filled')
colorbar;axis equal


hold on
plot(staloc(:,2),staloc(:,1),'r^','markersize',10)




hold on
text(-locnx(:,4),locnx(:,3),num2str(locnx(:,2)),'fontsize',2)


figure;scatter(-locnx(:,4),locnx(:,3),15,locnx(:,5),'filled')
colorbar;axis equal

hold on
text(-locnx(:,4),locnx(:,3),num2str(locnx(:,1)),'fontsize',2)



%[lat2((fcatal(j,6)-fcatal(j,8)-off22(1,1))/0.25+1,(fcatal(j,5)-fcatal(j,9)-off12(1,1))/0.25+1),...
%             lon2((fcatal(j,6)-fcatal(j,8)-off22(1,1))/0.25+1,(fcatal(j,5)-fcatal(j,9)-off12(1,1))/0.25+1)]


zq=[];zq1=[];zq2=[];
zq=zqn;zq1=zqn1;zq2=zqn2;



save(['Japan/Kii/kii_',num2str(wr(1)),num2str(wr(2)),num2str(wr(3)),'_hypoinv_depth_100sps'],'xq','yq','zq','zq1','zq2')




















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
  
 ix=1:1;



Kii_time=load('Kii_ide_time_new');


 
 wr1(:,1)=[1;2;3];

 
 wr=wr1(:,ix);



sta=['URSH';'MASH';'MGWH'];




%nnd=[datenum(2005,1,5)+8,datenum(2007,4,7)];

kiiden1=[];


for njap=1:size(Kii_time,1)
        
        nyr=round(floor(Kii_time(njap)/10000));
        nmon=round(floor(Kii_time(njap)/100)-floor(Kii_time(njap)/10000)*100);
        nda=round(Kii_time(njap)-floor(Kii_time(njap)/100)*100);
        nday=datenum(nyr,nmon,nda)-datenum(nyr-1,12,31);
        
        
         da=datevec(datenum(nyr,nmon,nda));
        
    
    
                    
if da(1)-2000 <= 9
    JYr1=['0',int2str(da(1)-2000)];

else
    JYr1=int2str(da(1)-2000);
end
                
if da(2) <= 9
    JMon1=['0',int2str(da(2))];

else
    JMon1=int2str(da(2));
end

if da(3) <= 9
    JDay1=['0',int2str(da(3))];

else
    JDay1=int2str(da(3));
end

   
    
    
    %nday=126:2:ytime(1,3)-1
  
    
 %   load(['mexico_maze/',num2str(nyr),'_',num2str(nday),'_',num2str(nday+1),...
 %       '/maze',num2str(wr(1)),num2str(wr(2)),num2str(wr(3)),'_','mazeden'])
    
 
    fna=['Japan/Kii/',JYr1,JMon1,JDay1,'/'];
fname1=[fna,'kii',num2str(wr(1)),num2str(wr(2)),num2str(wr(3)),'_','kiiden_100sps.mat'];

if exist(fname1)~=0
load(fname1)
    
    if size(kiiden,2)~=0
    
   kiiden(:,2)=kiiden(:,2)/86400+datenum(nyr,nmon,nda)-datenum(2003,12,31);
    kiiden1=[kiiden1;kiiden];
    
    end
end
    
    nday
    
end
   
ncc=0.1;

xxx=[];
xxx=kiiden1(:,6);



xxx(xxx<ncc/2)=0;
xxx(xxx>=ncc/2)=1;

xxx1=[];
xxx1=kiiden1(:,7);

xxx1(xxx1<ncc/2)=0;
xxx1(xxx1>=ncc/2)=1;

xxx2=[];
xxx2=kiiden1(:,8);

xxx2(xxx2<ncc/2)=0;
xxx2(xxx2>=ncc/2)=1;


xxx3=[];
xxx3=kiiden1(:,9);

xxx3(xxx3<ncc)=0;
xxx3(xxx3>=ncc)=1;


xxx4=[];
xxx4=kiiden1(:,9);

xxx4(xxx4<=0.9)=-1;
xxx4(xxx4>0.9)=0;

xxx4(xxx4==-1)=1;



xxx5=[];
xxx5=kiiden1(:,14);

xxx5(xxx5<0.5)=0;
xxx5(xxx5>=0.5)=1;



xxx6=[];
xxx6=kiiden1(:,16);

xxx6(xxx6<1.5 | xxx6>6)=0;
xxx6(xxx6~=0)=1;



xxx7=[];
xxx7=kiiden1(:,5);

xxx7(xxx7<6)=1;
xxx7(xxx7==6)=0;


% aa=[];
% aa=kiide1(:,22);
% aa(aa>2 | aa<1)=0;
% aa(aa<=2 & aa>=1)=1;
% 
% aa1=[];
% aa1=kiide1(:,19);
% aa1(aa1>=0.2)=0;
% aa1(aa1<0.2 & aa1~=0)=1;



xxxn=xxx.*xxx1.*xxx2.*xxx3.*xxx4.*xxx5.*xxx6.*xxx7;
xxxn=xxxn.*(1:size(kiiden1,1))';

xxxn(xxxn==0)=[];

kiide11=[];
kiide11=kiiden1(xxxn,:);

figure;
plot(kiide11(:,2),kiide11(:,4),'o','markersize',2)


figure;plot(kiide11(:,3),kiide11(:,4),'o','markersize',2);axis equal



%save('Japan/Kii/kiide11_123_cc01_soff5_100sps.mat','kiide11')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

  
 ix=1:1;



Kii_time=load('Kii_ide_time_new');


 
 wr1(:,1)=[1;2;3];

 
 wr=wr1(:,ix);



sta=['URSH';'MASH';'MGWH'];


load('Japan/Kii/kiide11_123_cc01_soff5_100sps.mat')


load(['Japan/Kii/kii_',num2str(wr(1)),num2str(wr(2)),num2str(wr(3)),'_hypoinv_depth_100sps'])

 lat2=interp2(zq,3);
 lon2=interp2(zq1,3);
 off12=interp2(xq,3);
 off22=interp2(yq,3);


kiide11(:,17:18)=NaN;

for j=1:size(kiide11,1)
%     if (kiide11(i,3)-xq(1,1))+1>=1 && (kiide11(i,3)-xq(1,1))+1<=size(xq,2) && ...
%             (kiide11(i,4)-yq(1,1))+1>=1 && (kiide11(i,4)-yq(1,1))+1<=size(yq,1)
%     kiide11(i,13:14)=[zq((kiide11(i,4)-yq(1,1))+1,(kiide11(i,3)-xq(1,1))+1), ...
%         zq1((kiide11(i,4)-yq(1,1))+1,(kiide11(i,3)-xq(1,1))+1)];
%     end
    
    if (kiide11(j,4)-off22(1,1))/0.25+1>=1 && (kiide11(j,4)-off22(1,1))/0.25+1<=size(lat2,1)
        if (kiide11(j,3)-off12(1,1))/0.25+1>=1 && (kiide11(j,3)-off12(1,1))/0.25+1<=size(lat2,2)
     kiide11(j,17:18)=[lat2((kiide11(j,4)-off22(1,1))/0.25+1,(kiide11(j,3)-off12(1,1))/0.25+1),...
             lon2((kiide11(j,4)-off22(1,1))/0.25+1,(kiide11(j,3)-off12(1,1))/0.25+1)];
        end
    end
    
    
end

kiide11(:,19:20)=[(kiide11(:,17)-34.5)/180*pi*6372.028,...
    (kiide11(:,18)-136.5)/180*pi*cos(34.5/180*pi)*6372.028];


% 
% tokai1=load('Japan/kii.mgw.urs.mas.hypo');
% 
% 
% tokai1(:,15:16)=[(tokai1(:,4)+tokai1(:,5)/60-34.5)/180*pi*6372.028,...
%     (tokai1(:,6)+tokai1(:,7)/60-136.5)/180*pi*cos(34.5/180*pi)*6372.028];




% 
% figure;plot(tokai1(:,6)+tokai1(:,7)/60,tokai1(:,4)+tokai1(:,5)/60,'.')
% hold on
% plot(kiide11(:,18),kiide11(:,17),'r+')
% axis equal
% 
% 
% figure;plot(tokai1(:,16),tokai1(:,15),'.')
% hold on
% plot(kiide11(:,20),kiide11(:,19),'r+')
% axis equal




    fc2=[];
    fc2=[cos(pi/3),sin(pi/3);-sin(pi/3),cos(pi/3)]*[kiide11(:,20),kiide11(:,19)]';

    kiide11(:,21:22)=fc2';
    
    
    
    load('/u/yajun/Desktop/tryout/Japan/Kii/staloc.mat')

staloc(:,3:4)=[(staloc(:,1)-34.5)/180*pi*6372.028,...
    (staloc(:,2)-136.5)/180*pi*cos(34.5/180*pi)*6372.028];
    
figure;plot(kiide11(:,20),kiide11(:,19),'o','markersize',2);axis equal


hold on
plot(staloc(:,4),staloc(:,3),'r^','markersize',10)
    
     %save(['Japan/Kii/','kiide11_',num2str(wr(1)),num2str(wr(2)),num2str(wr(3)),'_cc012'],'kiide11')
     
    save('Japan/Kii/kiide11_123_cc01_soff5_100sps.mat','kiide11') 
     
%%%%%%%%%%%%%%%%%%%%%%%%%%make floc%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
     
     
     
     clear
close all

    
  ix=1:1;



Kii_time=load('Kii_ide_time_new');


 
 wr1(:,1)=[1;2;3];

 
 wr=wr1(:,ix);



sta=['URSH';'MASH';'MGWH'];


load('/u/yajun/Desktop/tryout/Japan/Kii/staloc.mat')

staloc(:,3:4)=[(staloc(:,1)-34.5)/180*pi*6372.028,...
    (staloc(:,2)-136.5)/180*pi*cos(34.5/180*pi)*6372.028];



 load(['Japan/Kii/','kiide11_',num2str(wr(1)),num2str(wr(2)),num2str(wr(3)),'_cc01_soff5_100sps'])

load(['Japan/Kii/kii_',num2str(wr(1)),num2str(wr(2)),num2str(wr(3)),'_hypoinv_depth_100sps'])

 lat2=interp2(zq,3);
 lon2=interp2(zq1,3);
 off12=interp2(xq,3);
 off22=interp2(yq,3);



detemp=[];
centrloc=[];

nsamp=50;
 kiide11(:,23)=0;
 

%nmax=3000;
nmax=1200;

for xx=-210:nsamp:210 %-50:10:80
    for yy=-510:nsamp:510 %-110:10:100
        
        if (yy-off22(1,1))/0.25+1>=1 && (yy-off22(1,1))/0.25+1<=size(lat2,1)
        if (xx-off12(1,1))/0.25+1>=1 && (xx-off12(1,1))/0.25+1<=size(lat2,2)
            
     yy1=lat2((yy-off22(1,1))/0.25+1,(xx-off12(1,1))/0.25+1);
     xx1=lon2((yy-off22(1,1))/0.25+1,(xx-off12(1,1))/0.25+1);
        
    
        
        yy2=(yy1-34.5)/180*pi*6372.028;
        xx2=(xx1-136.5)/180*pi*cos(34.5/180*pi)*6372.028;
        
        if xx2>=8 || xx2<8
            % if xx1<=-99
        
        detemp1=[];z=0;
        
%         for i=1:size(arrayden1,1)
%             if abs(arrayden1(i,3)-xx)<=nsamp && abs(arrayden1(i,4)-yy)<=nsamp 
%                 z=z+1;
%                 detemp1(z,:)=[arrayrotn1(i,:),arrayoff1(i,:)];
%                 anum(z,:)=[1,arraycc1(i,:)];
%                 
%             end
%         end
        
        
        for i=1:size(kiide11,1)
            if sqrt((kiide11(i,20)-xx2)^2+(kiide11(i,19)-yy2)^2)<=6
                
                if abs(kiide11(i,3)-xx)<=nsamp && abs(kiide11(i,4)-yy)<=nsamp && ...
                abs(kiide11(i,3)-xx-(kiide11(i,4)-yy))<=nsamp
                
                z=z+1;
                detemp1(z,:)=kiide11(i,10:12);
                %kiide11(i,19)=-1;
               % anum(z,:)=[1,arraycc1(i,:)];
                end
                
            end
        end
        
        
%         if ix==3
%         
%         if xx1>=-99.38 || yy1>=17.8
%             nmax=500;
%         else
%             nmax=1800;
%         end
%         end
%         
        
        
        if size(detemp1,1)>=nmax
            
            
    %        figure;scatter(1:size(detemp1,1),detemp1(:,1),15,detemp1(:,size(detemp1,2)),'filled')
    %    caxis([0.1 0.3])
    %    colorbar
     
            
        meantemp1=[];stdtemp1=[];
       %meantemp1=mean(detemp1,1);
       % stdtemp1=std(detemp1,0,1);
       
      % meantemp1(:,1)=mean(detemp1(:,1));
      % stdtemp1(:,1)=std(detemp1(:,1));
       
       meantemp1=median(detemp1,1);
       stdtemp1=mad(detemp1,1,1);
       
        

        
        meantemp2=[];stdtemp2=[];
        
        
        for k=1:size(detemp1,2)
            
            detemp2=[];
        for j=1:size(detemp1,1)
            
               if  abs(detemp1(j,k)-meantemp1(k))<=3*stdtemp1(k) 
                   detemp2=[detemp2;detemp1(j,k)];
               end
        end
        

        
        
        meantemp2(:,k)=median(detemp2);
        stdtemp2(:,k)=mad(detemp2,1);
        
        end
        

%      
        
        
       detemp=[detemp;xx,yy,meantemp2,stdtemp2,size(detemp1,1)];
       
       centrloc=[centrloc;xx1,yy1,xx2,yy2];
       
        end
        
        
          %  end
        end
        
        
        
        end
        end
       
    end
    xx
end





xx=[2;6;7;8;14;24;38;40];detemp(xx,:)=[];centrloc(xx,:)=[];


for ii=1:3   
figure;
plot(kiide11(:,20),kiide11(:,19),'ro','markersize',3)
hold on
plot(centrloc(:,3),centrloc(:,4),'.')
hold on
quiver(centrloc(:,3),centrloc(:,4),0.3*cosd(detemp(:,ii+2)),0.3*sind(detemp(:,ii+2)))
    axis equal
    axis([-40 20 -30 40])
    hold on
    text(centrloc(:,3),centrloc(:,4),num2str(detemp(:,9)),'fontsize',5)

  hold on
plot(staloc(:,4),staloc(:,3),'k^','markersize',10)



end

figure;
plot(kiide11(:,20),kiide11(:,19),'ro','markersize',3)
hold on
plot(centrloc(:,3),centrloc(:,4),'.')
hold on
text(centrloc(:,3),centrloc(:,4),num2str((1:size(centrloc,1))'))

    axis equal
    axis([-40 20 -30 40])
    
    
    

mm=[];z=0;
for i=1:size(kiide11,1) 
    for j=1:size(detemp,1)
        %if j~=1 || j~=4
        if abs(kiide11(i,3)-detemp(j,1))<=nsamp && abs(kiide11(i,4)-detemp(j,2))<=nsamp && ...
                abs(kiide11(i,3)-detemp(j,1)-(kiide11(i,4)-detemp(j,2)))<=nsamp
            if sqrt((kiide11(i,20)-centrloc(j,3))^2+(kiide11(i,19)-centrloc(j,4))^2)<=6
                z=z+1;
                mm(z,:)=kiide11(i,:);
                break
            end
        end
       % end
    end
    
    
   i 
end

     
     
    figure;
plot(kiide11(:,20),kiide11(:,19),'bo','markersize',3)
hold on
plot(mm(:,20),mm(:,19),'r.')
%hold on
%quiver(centrloc(:,3),centrloc(:,4),0.3*cosd(detemp(:,ii+2)),0.3*sind(detemp(:,ii+2)))
    axis equal
    axis([-40 20 -30 40]) 
    
      hold on
plot(staloc(:,4),staloc(:,3),'k^','markersize',10)
hold on
plot(centrloc(:,3),centrloc(:,4),'k^','markerfacecolor','k','markersize',6)

   

figure;
plot(kiide11(:,2),kiide11(:,21),'bo','markersize',3)
hold on
plot(mm(:,2),mm(:,21),'r.')


   
floc=detemp;




save(['Japan/Kii/','floc_kii',num2str(wr(1)),num2str(wr(2)),num2str(wr(3)),'_100sps_new'],'floc')    
     
     
     
     
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % 
% % clear
% %   
% %  ix=1:1;
% % 
% % 
% % 
% % Kii_time=load('../Kii_data_time.xx');
% % 
% % 
% %  
% %  wr1(:,1)=[1;2;3];
% % 
% %  
% %  wr=wr1(:,ix);
% % 
% % 
% % 
% % sta=['URSH';'MASH';'MGWH'];
% % 
% % 
% % 
% % 
% % %nnd=[datenum(2005,1,5)+8,datenum(2007,4,7)];
% % 
% % kiiden1=[];
% % 
% % 
% % for njap=1:size(Kii_time,1)
% %         
% %         nyr=round(floor(Kii_time(njap)/10000));
% %         nmon=round(floor(Kii_time(njap)/100)-floor(Kii_time(njap)/10000)*100);
% %         nda=round(Kii_time(njap)-floor(Kii_time(njap)/100)*100);
% %         nday=datenum(nyr,nmon,nda)-datenum(nyr-1,12,31);
% %         
% %         
% %          da=datevec(datenum(nyr,nmon,nda));
% %         
% %     
% %     
% %                     
% % if da(1)-2000 <= 9
% %     JYr1=['0',int2str(da(1)-2000)];
% % 
% % else
% %     JYr1=int2str(da(1)-2000);
% % end
% %                 
% % if da(2) <= 9
% %     JMon1=['0',int2str(da(2))];
% % 
% % else
% %     JMon1=int2str(da(2));
% % end
% % 
% % if da(3) <= 9
% %     JDay1=['0',int2str(da(3))];
% % 
% % else
% %     JDay1=int2str(da(3));
% % end
% % 
% %    
% %     
% %     
% %     %nday=126:2:ytime(1,3)-1
% %   
% %     
% %  %   load(['mexico_maze/',num2str(nyr),'_',num2str(nday),'_',num2str(nday+1),...
% %  %       '/maze',num2str(wr(1)),num2str(wr(2)),num2str(wr(3)),'_','mazeden'])
% %     
% %  
% %     fna=['Japan/Kii/',JYr1,JMon1,JDay1,'/'];
% % fname1=[fna,'kii',num2str(wr(1)),num2str(wr(2)),num2str(wr(3)),'_','kiiden.mat'];
% % load(fname1)
% %     
% %     if size(kiiden,2)~=0
% %     
% %    kiiden(:,2)=kiiden(:,2)/86400+datenum(nyr,nmon,nda)-datenum(2003,12,31);
% %     kiiden1=[kiiden1;kiiden];
% %     
% %     end
% %     
% %     nday
% %     
% % end
% %    
% % ncc=0.12;
% % 
% % xxx=[];
% % xxx=kiiden1(:,6);
% % 
% % 
% % 
% % xxx(xxx<ncc/2)=0;
% % xxx(xxx>=ncc/2)=1;
% % 
% % xxx1=[];
% % xxx1=kiiden1(:,7);
% % 
% % xxx1(xxx1<ncc/2)=0;
% % xxx1(xxx1>=ncc/2)=1;
% % 
% % xxx2=[];
% % xxx2=kiiden1(:,8);
% % 
% % xxx2(xxx2<ncc/2)=0;
% % xxx2(xxx2>=ncc/2)=1;
% % 
% % 
% % xxx3=[];
% % xxx3=kiiden1(:,9);
% % 
% % xxx3(xxx3<ncc)=0;
% % xxx3(xxx3>=ncc)=1;
% % 
% % 
% % xxx4=[];
% % xxx4=kiiden1(:,9);
% % 
% % xxx4(xxx4<=0.9)=-1;
% % xxx4(xxx4>0.9)=0;
% % 
% % xxx4(xxx4==-1)=1;
% % 
% % 
% % % aa=[];
% % % aa=kiide1(:,22);
% % % aa(aa>2 | aa<1)=0;
% % % aa(aa<=2 & aa>=1)=1;
% % % 
% % % aa1=[];
% % % aa1=kiide1(:,19);
% % % aa1(aa1>=0.2)=0;
% % % aa1(aa1<0.2 & aa1~=0)=1;
% % 
% % 
% % 
% % xxxn=xxx.*xxx1.*xxx2.*xxx3.*xxx4;
% % xxxn=xxxn.*(1:size(kiiden1,1))';
% % 
% % xxxn(xxxn==0)=[];
% % 
% % kiide11=[];
% % kiide11=kiiden1(xxxn,:);
% % 
% % 
% % % figure;hist(kiide11(:,16),50)
% % % 
% % % figure;hist(kiide11(:,10:12),50)
% % % 
% % % figure;plot(kiide11(:,3),kiide11(:,4),'.');axis equal
% % % 
% % % figure;scatter(kiide11(:,2),kiide11(:,4),15,kiide11(:,14),'filled');colorbar
% % % 
% % % 
% % % figure;plot(kiiden1(:,2),kiiden1(:,4),'.')
% % % figure;plot(kiide11(:,2),kiide11(:,4),'.')
% % % figure;plot(kiide11(:,2),kiide11(:,3),'.')
% % 

     
     
     
     
     