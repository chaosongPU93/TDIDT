clear



 
 
 
Kii_time1=load('Kii_ide_time_new');


Kii_time2=load('Kii_JMA_time');

Kii_time=[Kii_time1;Kii_time2];
 
 wr1(:,1)=[1;2;3];

 



sta=['URSH';'MASH';'MGWH'];

 

 ntrio=1;
 wr=wr1(:,ntrio)


 load(['Japan/Kii/','floc_kii',num2str(wr(1)),num2str(wr(2)),num2str(wr(3)),'_100sps_new'])





 





xspara=[];
xspara=[(1:size(floc,1))',wr(1)*ones(size(floc,1),1),wr(2)*ones(size(floc,1),1),wr(3)*ones(size(floc,1),1),floc(:,1:5)];








 
 
ntar=(1:size(floc,1))';

load(['Japan/Kii/kii_',num2str(wr(1)),num2str(wr(2)),num2str(wr(3)),'_hypoinv_depth_100sps'])

 lat2=interp2(zq,3);
 lon2=interp2(zq1,3);
 off12=interp2(xq,3);
 off22=interp2(yq,3);
 
 
 for j=1:size(floc,1)
    
    if (floc(j,2)-off22(1,1))/0.25+1>=1 && (floc(j,2)-off22(1,1))/0.25+1<=size(lat2,1)
        if (floc(j,1)-off12(1,1))/0.25+1>=1 && (floc(j,1)-off12(1,1))/0.25+1<=size(lat2,2)
     centrloc44(j,:)=[lat2((floc(j,2)-off22(1,1))/0.25+1,(floc(j,1)-off12(1,1))/0.25+1),...
             lon2((floc(j,2)-off22(1,1))/0.25+1,(floc(j,1)-off12(1,1))/0.25+1)];
        end
    end
    
    
 end





fna='kii_'



fna1=['kii_',num2str(wr(1)),num2str(wr(2)),num2str(wr(3)),'_'];


nfm15x=[];


nfm15x=[centrloc44,...
    (centrloc44(:,1)-34.5)/180*pi*6372.028,...
    (centrloc44(:,2)-136.5)/180*pi*cos(34.5/180*pi)*6372.028];



nsamp=100;
hi=6;
lo=1;

winlen=nsamp*8;
cohlen=nsamp*2;

nn=size(ntar,1);



for njap=1:size(Kii_time,1) %1
        
%         nyr=round(floor(Kii_time(njap)/10000));
%         nmon=round(floor(Kii_time(njap)/100)-floor(Kii_time(njap)/10000)*100);
%         nday=round(Kii_time(njap)-floor(Kii_time(njap)/100)*100);
% 
%  da=datevec(datenum(nyr,nmon,nda));

 
 
    nyr=round(floor(Kii_time(njap)/10000));
        nmon=round(floor(Kii_time(njap)/100)-floor(Kii_time(njap)/10000)*100);
        nda=round(Kii_time(njap)-floor(Kii_time(njap)/100)*100);
        nday=datenum(nyr,nmon,nda)-datenum(nyr-1,12,31);
        
        
         da=datevec(datenum(nyr,nmon,nda));
        

nre=19;






                    
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








    
    filen=['Japan/Kii/',JYr1,JMon1,JDay1];

if exist(filen)~=0 && exist([filen,'/kii_1_123PGSSSImap.mat'])~=0    
for nr=1:nn
    
    nfami=ntar(nr);
    temp=[];fcat=[];eve=0;
    temp=xspara(nr,:);
    
%     for nxs=1:size(xspara,1)
%         if xspara(nxs,1)==nfami
%             temp=[temp;xspara(nxs,:)];
%         else
%             if xspara(nxs,1)>nfami
%                 break
%             end
%         end
%     end
    
   % if size(temp,1)~=0
        
    
        
    for i=1:size(temp,1)
        PGSSSImap=[];PGSSSIwa=[];

    filename1=[filen,'/',fna,num2str(temp(i,1)),'_',num2str(temp(i,2)),num2str(temp(i,3)),num2str(temp(i,4)),'PGSSSImap.mat'];
    filename2=[filen,'/',fna,num2str(temp(i,1)),'_',num2str(temp(i,2)),num2str(temp(i,3)),num2str(temp(i,4)),'PGSSSIwa.mat'];
%    filename3=[filen,'/',fna,num2str(temp(i,1)),'_',num2str(temp(i,2)),num2str(temp(i,3)),num2str(temp(i,4)),'PGSSSIde.mat'];
%    filename4=[filen,'/',fna,num2str(temp(i,1)),'_',num2str(temp(i,2)),num2str(temp(i,3)),num2str(temp(i,4)),'PGSSSIwa1.mat'];
     filename5=[filen,'/',fna,num2str(temp(i,1)),'_',num2str(temp(i,2)),num2str(temp(i,3)),num2str(temp(i,4)),'peak2.mat'];
   
   
   
   fcat=[];
   
    if exist(filename2)~=0 && exist(filename1)~=0
    load(filename1)
    load(filename2)
%    load(filename3)
%    load(filename4)
    load(filename5)
    
    if size(peak2,2)==6
    peak2=[zeros(size(peak2,1),1),peak2];
    end
    
    if size(PGSSSImap,2)~=0


    
    

    for j=1:size(PGSSSImap,1)

        
        if   abs(peak2(j,5))<=3 && abs(peak2(j,6))<=3 && abs(peak2(j,7))<=3 && abs(peak2(j,5)-peak2(j,6)-peak2(j,7))<=3 && ...
                 peak2(j,2)<0 && peak2(j,3)<0 && peak2(j,4)<0 && ...
                 PGSSSImap(j,14)>lo && PGSSSImap(j,14)<hi  && PGSSSImap(j,12)>=0.65 && PGSSSImap(j,11)<=0.15 &&...        %
             PGSSSImap(j,10)<=1 && PGSSSImap(j,10)>=0.4 && PGSSSImap(j,7)>=0.2 && PGSSSImap(j,8)>=0.2 && PGSSSImap(j,9)>=0.2
                
                

         cumsumPGSS=[];
        cumsumPGSI=[];
        cumsumSISS=[];
        cumsumPGSSSI=[];
        winmax=[];p=[];q=[];p1=[];q1=[];
        
        cumsumPGSS=cumsum(PGSSSIwa((j-1)*winlen+1:j*winlen,2).*PGSSSIwa((j-1)*winlen+1:j*winlen,3));
        cumsumPGSI=cumsum(PGSSSIwa((j-1)*winlen+1:j*winlen,2).*PGSSSIwa((j-1)*winlen+1:j*winlen,4));
        cumsumSISS=cumsum(PGSSSIwa((j-1)*winlen+1:j*winlen,3).*PGSSSIwa((j-1)*winlen+1:j*winlen,4));
        cumsumPGSSSI=(cumsumPGSS*PGSSSImap(j,7)/cumsumPGSS(winlen)+cumsumPGSI*PGSSSImap(j,8)/cumsumPGSI(winlen)+...
            cumsumSISS*PGSSSImap(j,9)/cumsumSISS(winlen))/3;
    
    
    
    winmax(1,:)=cumsumPGSSSI(cohlen,:);
    winmax(2:winlen-cohlen+1,:)=cumsumPGSSSI(cohlen+1:winlen,:)-cumsumPGSSSI(1:winlen-cohlen,:);
    
    [p,q]=max(winmax);
    ratio=[];
    ratio=(sum(PGSSSIwa((j-1)*winlen+q:(j-1)*winlen+q+cohlen-1,2).^2)/sum(PGSSSIwa((j-1)*winlen+q:(j-1)*winlen+q+cohlen-1,5).^2)+...
        sum(PGSSSIwa((j-1)*winlen+q:(j-1)*winlen+q+cohlen-1,3).^2)/sum(PGSSSIwa((j-1)*winlen+q:(j-1)*winlen+q+cohlen-1,6).^2)+...
        sum(PGSSSIwa((j-1)*winlen+q:(j-1)*winlen+q+cohlen-1,4).^2)/sum(PGSSSIwa((j-1)*winlen+q:(j-1)*winlen+q+cohlen-1,7).^2))/3;
    

    
    
    amp=[];
    amp=(sum(PGSSSIwa((j-1)*winlen+q:(j-1)*winlen+q+cohlen-1,2).*PGSSSIwa((j-1)*winlen+q:(j-1)*winlen+q+cohlen-1,3))+...
        sum(PGSSSIwa((j-1)*winlen+q:(j-1)*winlen+q+cohlen-1,2).*PGSSSIwa((j-1)*winlen+q:(j-1)*winlen+q+cohlen-1,4))+...
        sum(PGSSSIwa((j-1)*winlen+q:(j-1)*winlen+q+cohlen-1,3).*PGSSSIwa((j-1)*winlen+q:(j-1)*winlen+q+cohlen-1,4)))/3;
    
    [p1,q1]=max((PGSSSIwa((j-1)*winlen+q:(j-1)*winlen+q+cohlen-1,2).*PGSSSIwa((j-1)*winlen+q:(j-1)*winlen+q+cohlen-1,3))+...
        (PGSSSIwa((j-1)*winlen+q:(j-1)*winlen+q+cohlen-1,2).*PGSSSIwa((j-1)*winlen+q:(j-1)*winlen+q+cohlen-1,4))+...
        (PGSSSIwa((j-1)*winlen+q:(j-1)*winlen+q+cohlen-1,3).*PGSSSIwa((j-1)*winlen+q:(j-1)*winlen+q+cohlen-1,4)));
    
    

    
    eve=eve+1;
    fcat(eve,1:18)=[temp(i,1:6),PGSSSIwa((j-1)*winlen+1,1)+(q-1+q1-1)/nsamp,PGSSSImap(j,3:4),PGSSSImap(j,10),...
        p/cumsumPGSSSI(winlen),temp(i,7:9),ratio,amp,PGSSSImap(j,5),j];
    
    fcat(eve,27:30)=[max(peak2(j,2:4)),PGSSSImap(j,11),PGSSSImap(j,12),PGSSSImap(j,14)];
 %   fcat(eve,28)=f1(qqxx);


        
        end
         
        
    end
    end
    
    
    


   if size(fcat,2)~=0
     for jj=1:size(fcat,1) 
         
         if (fcat(jj,6)+fcat(jj,9)-off22(1,1))/0.25+1>0 && (fcat(jj,6)+fcat(jj,9)-off22(1,1))/0.25+1<=size(lat2,1) &&...
             (fcat(jj,5)+fcat(jj,8)-off12(1,1))/0.25+1>0 &&  (fcat(jj,5)+fcat(jj,8)-off12(1,1))/0.25+1<=size(lat2,2)
   
         fcat(jj,nre:nre+1)=[lat2((fcat(jj,6)+fcat(jj,9)-off22(1,1))/0.25+1,(fcat(jj,5)+fcat(jj,8)-off12(1,1))/0.25+1),...
             lon2((fcat(jj,6)+fcat(jj,9)-off22(1,1))/0.25+1,(fcat(jj,5)+fcat(jj,8)-off12(1,1))/0.25+1)];
         end
         
    
     end
     

    fcat(:,21:22)=[(fcat(:,19)-34.5)/180*pi*6372.028,...
    (fcat(:,20)-136.5)/180*pi*cos(34.5/180*pi)*6372.028];
% nfm15x=[centrloc44,...
%     (centrloc44(:,1)-34.5)/180*pi*6372.028,...
%     (centrloc44(:,2)-136.5)/180*pi*cos(34.5/180*pi)*6372.028];


    fc2=[];
    fc2=[cos(pi/3),sin(pi/3);-sin(pi/3),cos(pi/3)]*[fcat(:,22),fcat(:,21)]';
 
    fcat(:,23:24)=fc2';
    fcat(:,32)=sqrt((nfm15x(fcat(:,1),3)-fcat(:,21)).^2+(nfm15x(fcat(:,1),4)-fcat(:,22)).^2);
    

    
    
        zz=0;fcat1=[];
    for jj1=1:size(fcat,1)
        if fcat(jj1,15)>=2.5 && fcat(jj1,32)<=4 %&& log10(fcat(jj1,16))>=-13.5
        zz=zz+1;
        fcat1(zz,:)=fcat(jj1,:);
        end
    end



    
    fcat=[];fcat=fcat1;fcat1=[];
    
    
    
    
   end
    
    
    fn=[filen,'/',fna1,'fcat',num2str(nfami),'.mat'];
    
    
    save(fn,'fcat')
    
    
    end
    end
    
   % nfami
   % end
end




nnx=size(ntar,1); 


ths=1;



 
 
fcatal=[];denum=[];
 for mm=1:nnx
     m=ntar(mm);
      fcatnn=[];filename3=[];
      
      filename3=[filen,'/',fna1,'fcat',num2str(m),'.mat'];
      load(filename3)
      if size(fcat,2)~=0  
          
           



    

    
        



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      fct=[];
      fct=fcat(2:size(fcat,1),7)-fcat(1:size(fcat,1)-1,7);



nfc=1;

fcatnn=[];
z=0;

while nfc<=size(fcat,1)
    
    nfce=nfc;
    for i=nfce:size(fcat,1)-1
        if fct(i,1)<ths
            nfc=nfc+1;
        else
            break
        end
    end
    
    
    fcna_temp=[];
    fcna_temp=fcat(nfce:nfc,:);
    
    posn=size(fcna_temp,1);
    
    while posn>=1
        
       % fcna_t1=fcna_temp;
       fcm=[];fcmin=[];
    [fcm,fcmin]=max(fcna_temp(:,10));
    z=z+1;
    fcatnn(z,:)=fcna_temp(fcmin,:);
    
    tmax=[];
    tmax=fcna_temp(fcmin,7);
    
    for i=1:size(fcna_temp,1)
     
        if abs(fcna_temp(i,7)-tmax)<ths
            fcna_temp(i,7)=-1000;
            fcna_temp(i,10)=-1000;
        end
        
    end
    
         
    
    ps=[];
    ps=fcna_temp(:,7);
    ps(ps<0)=[];
    
    if size(ps,2)==0
 
    posn=size(ps,2);
    else
        posn=size(ps,1);
    end
    
    
    end
    
    nfc=nfc+1;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      
      fcatal=[fcatal;fcatnn];
      denum=[denum;size(fcatnn,1)];
      else
         denum=[denum;0]; 
      end
      
     % m
                      
                      
      
 end
  
 %save('fcatal','fcatal')
  cdenum=cumsum(denum);


 
 rads=pi/180.;
 erad=6372.028;
 lat0=48.0+26.32/60.;
 lon0=123.0+35.07/60.;
 srad=erad*cos(lat0*rads);
 

 
 cdenum=[0;cdenum];
 
  save([filen,'/',fna1,'denum.mat'],'denum')
 save([filen,'/',fna1,'cdenum.mat'],'cdenum')
 


 
 
save([filen,'/',fna1,'fcatal'],'fcatal')





end


njap  %nday

end

%figure;plot(fcatal(:,7),fcatal(:,23),'.')
% figure;scatter(fcatal(:,7),fcatal(:,22),10,fcatal(:,10),'filled');colorbar;caxis([0.55 0.7])
%  figure;scatter(fcatal(:,20),fcatal(:,19),10,fcatal(:,10),'filled');colorbar;caxis([0.55 0.8]);axis equal
%  figure;plot(fcatal(:,20),fcatal(:,19),'ro','markersize',3);axis equal


fcna1=[];




for njap=1:size(Kii_time,1) %1
        nyr=round(floor(Kii_time(njap)/10000));
        nmon=round(floor(Kii_time(njap)/100)-floor(Kii_time(njap)/10000)*100);
        nda=round(Kii_time(njap)-floor(Kii_time(njap)/100)*100);
        nday=datenum(nyr,nmon,nda)-datenum(2003,12,31);
        
        
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

       filen=['Japan/Kii/',JYr1,JMon1,JDay1];     

if exist([filen,'/',fna1,'fcatal.mat'])~=0       
load([filen,'/',fna1,'fcatal'])
if size(fcatal,2)~=0
fcatal(:,25:26)=[fcatal(:,6)+fcatal(:,9),fcatal(:,5)+fcatal(:,8)];
fcatal(:,7)=fcatal(:,7)+nday*3600*24;

fcna1=[fcna1;fcatal];
end
end



end








fcn=[];fcnin=[];

[fcn,fcnin]=sort(fcna1(:,7));
fcna1=fcna1(fcnin,:);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
            fct=[];
      fct=fcna1(2:size(fcna1,1),7)-fcna1(1:size(fcna1,1)-1,7);



nfc=1;

fcna=[];
z=0;




xx=0;fx=[];
while nfc<=size(fcna1,1)
    
    nfce=nfc;
    for i=nfce:size(fcna1,1)-1
        if fct(i,1)<ths
            nfc=nfc+1;
        else
            break
        end
    end
    
    
    fcna_temp=[];
    fcna_temp=fcna1(nfce:nfc,:);
    
    posn=size(fcna_temp,1);
      fxc=[];fxc=fcna_temp;
    while posn>=1
        
        
        fcna_temp1=[];
%         for i=1:size(fcna_temp,1)
%             if fcna_temp(i,32)>=4
%                fcna_temp1=[fcna_temp1;i];
%             end
%         end
        
        fmax=[];
       % if size(fcna_temp1,2)==0
            fcna_temp1=(1:size(fcna_temp,1))';
             fmax=fcna_temp(fcna_temp1,10);
      %  else
           %  fmax=fcna_temp(fcna_temp1,31);
      %  end
        
        
       % fcna_t1=fcna_temp;
       fcm=[];fcmin=[];
       
     % 
     % fmax=fcna_temp(fcna_temp1,10);
       
    [fcm,fcmin]=max(fmax);
    z=z+1;
    fcna(z,:)=fcna_temp(fcna_temp1(fcmin),:);
    
    tmax=[];
    tmax=fcna_temp(fcna_temp1(fcmin),7);
    
    for i=1:size(fcna_temp,1)
     
        if fcna_temp(i,7)~=-1000 && abs(fcna_temp(i,7)-tmax)<ths 
            fcna_temp(i,7)=-1000;
            fcna_temp(i,10)=-1000;
           % fcna_temp(i,27:32)=-1000;
            %fcna_temp(i,28)=-1000;
%             
%             fcna_temp1(i,7)=-1000;
%             fcna_temp1(i,10)=-1000;
%             fcna_temp1(i,27)=1000;
            
        end
        
    end
    
         
    
    ps=[];
    ps=fcna_temp(:,7);
    ps(ps<0)=[];
    
    if size(ps,2)==0
    posn=size(ps,2);
    else
        posn=size(ps,1);
    end
    
%     for xz=1:size(fxc,1)
%         
%         if abs(fxc(xz,7)-tmax)<ths
%         xx=xx+1;
%     
%     fx(xx,:)=[sqrt((fxc(xz,21)-fcna_temp(fcna_temp1(fcmin),21))^2+(fxc(xz,22)-fcna_temp(fcna_temp1(fcmin),22))^2),...
%         fxc(xz,1),fcna_temp(fcna_temp1(fcmin),1)];
%         end
%     end
%     
    
        
    end
    
    nfc=nfc+1
end
fcn=[];fcnin=[];

[fcn,fcnin]=sort(fcna(:,7));
fcna=fcna(fcnin,:);

% ff=fcna(:,10);
% ff(ff<=0.45)=NaN;
% ff(ff>0.45)=1;


 figure;plot(fcna(:,7)/86400,fcna(:,23),'.')
 
 figure;plot(fcna(:,22),fcna(:,21),'.');axis equal
 
%  aa=fcna(:,15);
%  aa(aa<2.5)=NaN;
%  aa(aa>=2.5)=1;
 
  figure;plot(fcna(:,22),fcna(:,21),'.');axis equal
   figure;plot(fcna(:,7)/86400,fcna(:,23),'.')
 
save('Japan/Kii/kii_123_fcna_6s_100sps_cc04_r25_peak0_1to6Hz_new.mat','fcna')

 figure;scatter(fcna(:,7)/86400,fcna(:,23),15,fcna(:,15),'filled');colorbar;caxis([2 6])











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear



 
 
 
%Kii_time=load('Kii_ide_time_new');

Kii_time1=load('Kii_ide_time_new');


Kii_time2=load('Kii_JMA_time');

Kii_time=[Kii_time1;Kii_time2];


 
 wr1(:,1)=[1;2;3];

 



sta=['URSH';'MASH';'MGWH'];

 

 ntrio=1;
 wr=wr1(:,ntrio)


 load(['Japan/Kii/','floc_kii',num2str(wr(1)),num2str(wr(2)),num2str(wr(3)),'_100sps_new'])





 





xspara=[];
xspara=[(1:size(floc,1))',wr(1)*ones(size(floc,1),1),wr(2)*ones(size(floc,1),1),wr(3)*ones(size(floc,1),1),floc(:,1:5)];








 
 
ntar=(1:size(floc,1))';

load(['Japan/Kii/kii_',num2str(wr(1)),num2str(wr(2)),num2str(wr(3)),'_hypoinv_depth_100sps'])

 lat2=interp2(zq,3);
 lon2=interp2(zq1,3);
 off12=interp2(xq,3);
 off22=interp2(yq,3);
 
 
 for j=1:size(floc,1)
    
    if (floc(j,2)-off22(1,1))/0.25+1>=1 && (floc(j,2)-off22(1,1))/0.25+1<=size(lat2,1)
        if (floc(j,1)-off12(1,1))/0.25+1>=1 && (floc(j,1)-off12(1,1))/0.25+1<=size(lat2,2)
     centrloc44(j,:)=[lat2((floc(j,2)-off22(1,1))/0.25+1,(floc(j,1)-off12(1,1))/0.25+1),...
             lon2((floc(j,2)-off22(1,1))/0.25+1,(floc(j,1)-off12(1,1))/0.25+1)];
        end
    end
    
    
 end





fna='kii_128s_'



fna1=['kii_128s_',num2str(wr(1)),num2str(wr(2)),num2str(wr(3)),'_'];


nfm15x=[];


nfm15x=[centrloc44,...
    (centrloc44(:,1)-34.5)/180*pi*6372.028,...
    (centrloc44(:,2)-136.5)/180*pi*cos(34.5/180*pi)*6372.028];





for njap=1:size(Kii_time,1) %1
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



nre=19;

    Fs=40;

L=128*Fs;
NFFT=L;

f1=Fs/2*linspace(0,1,NFFT/2+1);
f2=Fs/2*linspace(1,0,NFFT/2+1);
f=[f1,-f2(2:NFFT/2)];
f=f';
    





nn=size(ntar,1);





    
    %filen=['mexico_maze/',num2str(nyear),'_',num2str(nday),'_',num2str(nday+1)];
    filen=['Japan/Kii/',JYr1,JMon1,JDay1];
if exist(filen)~=0 && exist([filen,'/kii_128s_1_123PGSSSImap.mat'])~=0     
for nr=1:nn
    
    nfami=ntar(nr);
    temp=[];fcat=[];eve=0;
    
    for nxs=1:size(xspara,1)
        if xspara(nxs,1)==nfami
            temp=[temp;xspara(nxs,:)];
        else
            if xspara(nxs,1)>nfami
                break
            end
        end
    end
    
   % if size(temp,1)~=0
        
    
        
    for i=1:size(temp,1)
        PGSSSImap=[];PGSSSIwa=[];PGSSSIde=[];
       % filen=[num2str(nyear),'_',num2str(nday),'_',num2str(nday+1),'_nc'];
       % filen=['array_346_EN/',num2str(nyear),'_',num2str(nday),'_',num2str(nday+1)];

    filename1=[filen,'/',fna,num2str(temp(i,1)),'_',num2str(temp(i,2)),num2str(temp(i,3)),num2str(temp(i,4)),'PGSSSImap.mat'];
    %filename2=[filen,'/',fna,num2str(temp(i,1)),'_',num2str(temp(i,2)),num2str(temp(i,3)),num2str(temp(i,4)),'PGSSSIwa.mat'];
    %filename3=[filen,'/',fna,num2str(temp(i,1)),'_',num2str(temp(i,2)),num2str(temp(i,3)),num2str(temp(i,4)),'PGSSSIde.mat'];

   fcat=[];
   
    if exist(filename1)~=0
    load(filename1)
    %load(filename2)
    %load(filename3)
    %load(filename4)
    
    if size(PGSSSImap,2)~=0


    
    

    for j=1:size(PGSSSImap,1)
        
        
%                     pgfft=[];ssfft=[];sifft=[];psssum=[];
%             
%             pgfft=fft(realPGC(istart:iend),NFFT);
%             ssfft=fft(realSSIB(istart-imaxPGSSwr:iend-imaxPGSSwr),NFFT);
%             sifft=fft(realSILB(istart-imaxPGSIwr:iend-imaxPGSIwr),NFFT);
%    
%    psssum=1/max(abs(pgfft(1:NFFT/2+1)))*abs(pgfft(1:NFFT/2+1))+...
%        1/max(abs(ssfft(1:NFFT/2+1)))*abs(ssfft(1:NFFT/2+1))+...
%        1/max(abs(sifft(1:NFFT/2+1)))*abs(sifft(1:NFFT/2+1));
%    
%    
%    pp=[];qq=[];
%    [pp,qq]=max(psssum);  
        

        
        if PGSSSImap(j,14)>=1.5 && PGSSSImap(j,14)<=6  && PGSSSImap(j,12)>=0.65 && PGSSSImap(j,11)<=0.2 &&...        %
             PGSSSImap(j,10)<=1 && PGSSSImap(j,10)>=0.15 && PGSSSImap(j,7)>=0.075 && PGSSSImap(j,8)>=0.075 && PGSSSImap(j,9)>=0.075
            
            
            %abs(PGSSSImap(j,4))<=1   && PGSSSImap(j,13)>=1.5 && PGSSSImap(j,13)<=6 &&... %&& PGSSSImap(j,10)<0.2
             %   PGSSSImap(j,4)>=0.14 && PGSSSIde(j,6)>=0.07 && PGSSSIde(j,8)>=0.07 && PGSSSIde(j,10)>=0.07
        
           eve=eve+1;
    fcat(eve,1:10)=[temp(i,1:6),PGSSSImap(j,2)+128/2,PGSSSImap(j,3:4),PGSSSImap(j,10)];
      % fcat(eve,27:30)=PGSSSImap(j,10:13);
      
       fcat(eve,27:30)=[PGSSSImap(j,5),PGSSSImap(j,11),PGSSSImap(j,12),PGSSSImap(j,14)];
    

        
        end
         
        
    end
    end
    
    
    
    end
    end

   if size(fcat,2)~=0
       
       %fcat(:,11:32)=0;
       
     for jj=1:size(fcat,1) %cdenum(im,1)+1:cdenum(im+1,1)
         if (fcat(jj,6)+fcat(jj,9)-off22(1,1))/0.25+1>0 && (fcat(jj,6)+fcat(jj,9)-off22(1,1))/0.25+1<=size(lat2,1) &&...
             (fcat(jj,5)+fcat(jj,8)-off12(1,1))/0.25+1>0 &&  (fcat(jj,5)+fcat(jj,8)-off12(1,1))/0.25+1<=size(lat2,2)
         fcat(jj,11:12)=[lat2((fcat(jj,6)+fcat(jj,9)-off22(1,1))/0.25+1,(fcat(jj,5)+fcat(jj,8)-off12(1,1))/0.25+1),...
             lon2((fcat(jj,6)+fcat(jj,9)-off22(1,1))/0.25+1,(fcat(jj,5)+fcat(jj,8)-off12(1,1))/0.25+1)];
         
         else
             
           fcat(jj,11:12)=[NaN,NaN];  
         end
    
     end
     

    fcat(:,13:14)=[(fcat(:,11)-34.5)/180*pi*6372.028,...
    (fcat(:,12)-136.5)/180*pi*cos(34.5/180*pi)*6372.028];



    fc2=[];
    fc2=[cos(pi/3),sin(pi/3);-sin(pi/3),cos(pi/3)]*[fcat(:,14),fcat(:,13)]';
 
    fcat(:,15:16)=fc2';
    fcat(:,32)=sqrt((nfm15x(fcat(:,1),3)-fcat(:,13)).^2+(nfm15x(fcat(:,1),4)-fcat(:,14)).^2);
    
        zz=0;fcat1=[];
    for jj1=1:size(fcat,1)
        if fcat(jj1,32)<=5
        zz=zz+1;
        fcat1(zz,:)=fcat(jj1,:);
        end
    end
    
    fcat=[];fcat=fcat1;fcat1=[];
   end
    
    
    fn=[filen,'/',fna1,'fcat',num2str(nfami),'.mat'];
    
    
    save(fn,'fcat')
    
   % nfami
   % end
end




nnx=size(ntar,1); 


ths=1;



 
 
fcatal=[];denum=[];
 for mm=1:nnx
     m=ntar(mm);
     filename3=[];
      
      filename3=[filen,'/',fna1,'fcat',num2str(m),'.mat'];
      load(filename3)
      if size(fcat,2)~=0  
          
           



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      

      
      fcatal=[fcatal;fcat];
      denum=[denum;size(fcat,1)];
      else
         denum=[denum;0]; 
      end
      
     % m
                      
                      
      
 end
 
 
 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcatal1=[];
z=0;

if size(fcatal,2)~=0


 fcn=[];fcnin=[];

[fcn,fcnin]=sort(fcatal(:,7));
fcatal=fcatal(fcnin,:);

      fct=[];
      fct=fcatal(2:size(fcatal,1),7)-fcatal(1:size(fcatal,1)-1,7);



nfc=1;

while nfc<=size(fcatal,1)
    
    nfce=nfc;
    for i=nfce:size(fcatal,1)-1
        if fct(i,1)<ths
            nfc=nfc+1;
        else
            break
        end
    end
    
    
    fcna_temp=[];
    fcna_temp=fcatal(nfce:nfc,:);
    
    posn=size(fcna_temp,1);
    
    while posn>=1
        
       % fcna_t1=fcna_temp;
       fcm=[];fcmin=[];
    [fcm,fcmin]=max(fcna_temp(:,10));
    z=z+1;
    fcatal1(z,:)=fcna_temp(fcmin,:);
    
    tmax=[];
    tmax=fcna_temp(fcmin,7);
    
    for i=1:size(fcna_temp,1)
     
        if abs(fcna_temp(i,7)-tmax)<ths
            fcna_temp(i,7)=-1000;
            fcna_temp(i,10)=-1000;
        end
        
    end
    
         
    
    ps=[];
    ps=fcna_temp(:,7);
    ps(ps<0)=[];
    
    if size(ps,2)==0
 
    posn=size(ps,2);
    else
        posn=size(ps,1);
    end
    
    
    end
    
    nfc=nfc+1;
    
end
 
 
end
 
 
 
 
  

  cdenum=cumsum(denum);
 

 
 cdenum=[0;cdenum];
 
  save([filen,'/',fna1,'denum.mat'],'denum')
 save([filen,'/',fna1,'cdenum.mat'],'cdenum')
 


 
 
save([filen,'/',fna1,'fcatal1'],'fcatal1')







end
njap

end




fcna1=[];




for njap=1:size(Kii_time,1) %1
        nyr=round(floor(Kii_time(njap)/10000));
        nmon=round(floor(Kii_time(njap)/100)-floor(Kii_time(njap)/10000)*100);
        nda=round(Kii_time(njap)-floor(Kii_time(njap)/100)*100);
        nday=datenum(nyr,nmon,nda)-datenum(2003,12,31);
        
        
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

       filen=['Japan/Kii/',JYr1,JMon1,JDay1];     

if exist([filen,'/',fna1,'fcatal1.mat'])~=0       
load([filen,'/',fna1,'fcatal1'])
if size(fcatal1,2)~=0
fcatal1(:,17:18)=[fcatal1(:,6)+fcatal1(:,9),fcatal1(:,5)+fcatal1(:,8)];
fcatal1(:,7)=fcatal1(:,7)+nday*3600*24;

fcna1=[fcna1;fcatal1];
end
end



end







fcn=[];fcnin=[];

[fcn,fcnin]=sort(fcna1(:,7));
fcna1=fcna1(fcnin,:);

fcna=fcna1;

fcna_128=fcna;
%fcna_128(:,19:31)=[];


figure;plot(fcna_128(:,14),fcna_128(:,13),'o','markersize',2);axis equal

figure;plot(fcna_128(:,7)/86400,fcna_128(:,15),'o','markersize',2);

save(['Japan/Kii/','kii_123','_fcna_128s_100sps_cc015_1to6Hz_new'],'fcna_128')












