
clear

%addpath /media/Elements/Japan/Kii/Hinet


Kii_time=load('Kii_ide_time1');

%Kii_time=load('../Kii_data_time.xx');


%addpath ./
% 
%  wr1(:,1)=[3;4;6];
%  wr1(:,2)=[3;1;4];
%  wr1(:,3)=[4;5;6];
%  wr1(:,4)=[1;2;3];
%  
%   wr1(:,5)=[2;1;7];
%  wr1(:,6)=[6;4;8];
%  
 
 
 
 
 
 wr1(:,1)=[1;2;3];
 wr1(:,2)=[2;1;4];

 %wr1(:,3)=[5;6;7];
 
 
 
 



sta=['URSH';'MASH';'MGWH';'WATH'];



nsamp=100;


mshift=7*nsamp;

 %for win_length=1:2
 win_length=2;
    if win_length==1
       winlen=4*nsamp;

       winoff=nsamp/1;

       xcmaxAVEnmin=0.4%0.4 for 4s 1.5-6 Hz; 0.36 for 4s2-8 Hz ;  0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz
       fna1='/05_'
       
    else
       winlen=128*nsamp;

       winoff=nsamp*8;
       
       xcmaxAVEnmin=0.08%0.4 for 4s 1.5-6 Hz; 0.36 for 4s2-8 Hz ;  0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz
       
       loopmax=6;
       
       fna1='/05_128s_'
    end
    
 
    
    
    
    
    
   % timsPGC=0:0.025:(86400*2-0.025);
    
      timsPGC=0:1/nsamp:(86400*2-1/nsamp);


ss=0;

        
     wlength=86400*nsamp;   
     
     cxx=(3*nsamp:winoff:(wlength-mshift-3*nsamp-winlen))';
     
     cxx=cxx-mshift+15;
     



 nwr=2;
 
     wr=wr1(:,nwr)   
     
        

Fs=nsamp;
L=winlen;
NFFT=L;

f1=Fs/2*linspace(0,1,NFFT/2+1);
f2=Fs/2*linspace(1,0,NFFT/2+1);
f=[f1,-f2(2:NFFT/2)];
f=f';    

    
    
    for njap=251:251 %247:255 %233:440 %1:220 %181:181 %261:size(Kii_time,1)
        
        Kii_time(njap)
        
        nyr=round(floor(Kii_time(njap)/10000));
        nmon=round(floor(Kii_time(njap)/100)-floor(Kii_time(njap)/10000)*100);
        nday=round(Kii_time(njap)-floor(Kii_time(njap)/100)*100);
        
    da=datevec(datenum(nyr,nmon,nday));
        
    
    
                    
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




 
            
            
     
             %   prename1=['../array/rdseed/rdseed/mexico/',YEAR,'.',JDAY1,'.00.00.00.0000.TO'];

          %  prename1=['/media/Elements/Japan/Kii/Hinet/',JYr1,JMon1,JDay1,'/N'];
          
          prename1=['/data2/yajun/Japan/Kii/hinet/',JYr1,JMon1,JDay1,'/N'];
             prename2=['/data2/yajun/Japan/Kii/hinet_add/',JYr1,JMon1,JDay1,'/N'];
            
    fna=['Japan/Kii/',JYr1,JMon1,JDay1,'/'];
    
for narr1=1:3
        
    narr=wr(narr1);
    
        PGCEfdx=zeros(wlength,1);
        PGCNfdx=zeros(wlength,1);

                
            
if narr~=4
                
PGCEdat1=[prename1,'.',sta(narr,:),'.E'];
PGCNdat1=[prename1,'.',sta(narr,:),'.N'];
else
    
 PGCEdat1=[prename2,'.',sta(narr,:),'.E'];
PGCNdat1=[prename2,'.',sta(narr,:),'.N'];
end


if exist(PGCEdat1)~=0 


%if exist(PGCEdat1)~=0 && exist(PGCNdat1)~=0 && exist(PGCZdat1)~=0 
[PGCE1,HdrDataPGC1,tnuPGC1,pobjPGC1,timsPGC1]=readsac(PGCEdat1,0,'l');
[PGCN1,~,~,~,~]=readsac(PGCNdat1,0,'l');
%[PGCZ1,~,~,~,~]=readsac(PGCZdat1,0,'l');
clear timsPGC1

if narr~=4
if size(PGCE1,1)~=100*86400
    PGCE1=[PGCE1;zeros(100*86400-size(PGCE1,1),1)];
end

if size(PGCN1,1)~=100*86400
    PGCN1=[PGCN1;zeros(100*86400-size(PGCN1,1),1)];
end


else
    
if size(PGCE1,1)~=100*86400
    PGCE1=[-PGCE1;zeros(100*86400-size(PGCE1,1),1)];
end

if size(PGCN1,1)~=100*86400
    PGCN1=[-PGCN1;zeros(100*86400-size(PGCN1,1),1)];
end
   
    
end


tracelenPGC1=length(PGCE1);



  
%cosine taper before filtering:
x=(0:pi/200:pi/2-pi/200)';
% %Seems to be necessary at the start of each day for PGCE:
    PGCE1(1:100)=sin(x).*PGCE1(1:100);

    PGCN1(1:100)=sin(x).*PGCN1(1:100);

    
%     nnz=100*86400;
%     
%     PGCE1(1+nnz:100+nnz)=sin(x).*PGCE1(1+nnz:100+nnz);
% 
%     PGCN1(1+nnz:100+nnz)=sin(x).*PGCN1(1+nnz:100+nnz);
%     


x=flipud(x);
PGCE1(tracelenPGC1-99:tracelenPGC1)=sin(x).*PGCE1(tracelenPGC1-99:tracelenPGC1);
PGCN1(tracelenPGC1-99:tracelenPGC1)=sin(x).*PGCN1(tracelenPGC1-99:tracelenPGC1);

% 
% PGCE1(tracelenPGC1-99-nnz:tracelenPGC1-nnz)=sin(x).*PGCE1(tracelenPGC1-99-nnz:tracelenPGC1-nnz);
% PGCN1(tracelenPGC1-99-nnz:tracelenPGC1-nnz)=sin(x).*PGCN1(tracelenPGC1-99-nnz:tracelenPGC1-nnz);





hi=6;
lo=1.5;
%Filter data:
npo=2;
npa=1;

% [PGCEf]=bandpass(PGCE1,100,lo,hi,npo,npa,'butter');
% [PGCNf]=bandpass(PGCN1,100,lo,hi,npo,npa,'butter');
% 
% 
% clear PGCE1 PGCN1% PGCZ1
% 
%  PGCEfdx = resample(PGCEf,2,5);
%  PGCNfdx = resample(PGCNf,2,5);
%  
% 
% 
% 
%  clear PGCEf PGCNf 
 

 [PGCEfdx]=bandpass(PGCE1,nsamp,lo,hi,npo,npa,'butter');
[PGCNfdx]=bandpass(PGCN1,nsamp,lo,hi,npo,npa,'butter');


clear PGCE1 PGCN1% PGCZ1

% PGCEfdx = resample(PGCEf,2,5);
% PGCNfdx = resample(PGCNf,2,5);
 



 %clear PGCEf PGCNf 
 

 

end

        
        
if narr1==1

PGCEfd=PGCEfdx;

PGCNfd=PGCNfdx; 
elseif narr1==2

SSIBEfd=PGCEfdx;

SSIBNfd=PGCNfdx; 
elseif narr1==3
    
SILBEfd=PGCEfdx;

SILBNfd=PGCNfdx; 
end


clear PGCEfdx PGCNfdx
end

        
        
   

        


     zz=0;kiide=[];nsect=200;
     
     
if max(abs(PGCEfd))~=0 && max(abs(SSIBEfd))~=0 && max(abs(SILBEfd))~=0 
    nslen=ceil((size(cxx,1)-2)/nsect);
    
    for ii=1:nsect
        
        
        
        z1=(ii-1)*nslen+1+1;
        z2=min(ii*nslen+1,size(cxx,1)-1);
        
        PGC=[];
        PGC=[PGCEfd(cxx(z1)+mshift+1-mshift:cxx(z2)+mshift+winlen+mshift),...
                PGCNfd(cxx(z1)+mshift+1-mshift:cxx(z2)+mshift+winlen+mshift)];
            
         SSIB=[];
        SSIB=[SSIBEfd(cxx(z1)+mshift+1-mshift:cxx(z2)+mshift+winlen+mshift),...
                SSIBNfd(cxx(z1)+mshift+1-mshift:cxx(z2)+mshift+winlen+mshift)];
        
         SILB=[];
        SILB=[SILBEfd(cxx(z1)+mshift+1-mshift:cxx(z2)+mshift+winlen+mshift),...
                SILBNfd(cxx(z1)+mshift+1-mshift:cxx(z2)+mshift+winlen+mshift)];
        
            
            PGC2=cumsum(PGC(:,1).^2+PGC(:,2).^2,1);
            PGC2=[zeros(1,size(PGC2,2));PGC2];
            SSIB2=cumsum(SSIB(:,1).^2+SSIB(:,2).^2,1);
             SSIB2=[zeros(1,size(SSIB2,2));SSIB2];
             SILB2=cumsum(SILB(:,1).^2+SILB(:,2).^2,1);
              SILB2=[zeros(1,size(SILB2,2));SILB2];
              
              
              
     lenx=size(PGC,1)-2*mshift;        
PGSSx=zeros(lenx, 2*mshift+1);
PGSIx=zeros(lenx, 2*mshift+1);
SSSIx=zeros(lenx, 2*mshift+1);
for n=-mshift:mshift;
    PGSSx(:,n+mshift+1)=PGC(1+mshift:size(PGC,1)-mshift,1).* ...
        SSIB(1+mshift+n:size(PGC,1)-mshift+n,1)+PGC(1+mshift:size(PGC,1)-mshift,2).* ...
        SSIB(1+mshift+n:size(PGC,1)-mshift+n,2);
    
    
    PGSIx(:,n+mshift+1)=PGC(1+mshift:size(PGC,1)-mshift,1).* ...
        SILB(1+mshift+n:size(PGC,1)-mshift+n,1)+PGC(1+mshift:size(PGC,1)-mshift,2).* ...
        SILB(1+mshift+n:size(PGC,1)-mshift+n,2);
    
    
    SSSIx(:,n+mshift+1)=SSIB(1+mshift:size(PGC,1)-mshift,1).* ...
        SILB(1+mshift+n:size(PGC,1)-mshift+n,1)+SSIB(1+mshift:size(PGC,1)-mshift,2).* ...
        SILB(1+mshift+n:size(PGC,1)-mshift+n,2);
end

    clear PGC SSIB SILB

cumsumPGSS=cumsum(PGSSx);
cumsumPGSS=[zeros(1,size(cumsumPGSS,2));cumsumPGSS];
clear PGSSx
cumsumPGSI=cumsum(PGSIx);
cumsumPGSI=[zeros(1,size(cumsumPGSI,2));cumsumPGSI];
clear PGSIx
cumsumSSSI=cumsum(SSSIx);
cumsumSSSI=[zeros(1,size(cumsumSSSI,2));cumsumSSSI];
clear SSSIx 
            
  
nwin=z2-z1+1;

sumsPGSS=zeros(nwin,2*mshift+1);
sumsPGSI=zeros(nwin,2*mshift+1);
sumsSSSI=zeros(nwin,2*mshift+1);
sumsPGC2=zeros(nwin,2*mshift+1);
sumsSSIB2=zeros(nwin,2*mshift+1);
sumsSSIB2b=zeros(nwin,2*mshift+1);
sumsSILB2=zeros(nwin,2*mshift+1);



for n=1:nwin
    istart=1+(n-1)*winoff;
    iend=istart+winlen;
    
    sumsPGSS(n,:)=cumsumPGSS(iend,:)-cumsumPGSS(istart,:); 
    sumsPGSI(n,:)=cumsumPGSI(iend,:)-cumsumPGSI(istart,:);
    sumsSSSI(n,:)=cumsumSSSI(iend,:)-cumsumSSSI(istart,:);
    sumsPGC2(n,:)=PGC2(iend+mshift)-PGC2(istart+mshift);  %PGC2 is cumsummed. Yes, +mshift.
    sumsSSIB2b(n,:)=SSIB2(iend+mshift)-SSIB2(istart+mshift); %Similar, for the SILB-SSIB connection.
    for m=-mshift:mshift;
        sumsSSIB2(n,m+mshift+1)=SSIB2(iend+mshift+m)-SSIB2(istart+mshift+m); %+m??? (yes).
        sumsSILB2(n,m+mshift+1)=SILB2(iend+mshift+m)-SILB2(istart+mshift+m);
    end
end
clear PGC2 SSIB2 SILB2 cumsumPGSS cumsumPGSI cumsumSSSI


sumsPGSSn=sumsPGSS./sqrt(sumsPGC2.*sumsSSIB2);
sumsPGSIn=sumsPGSI./sqrt(sumsPGC2.*sumsSILB2);
sumsSSSIn=sumsSSSI./sqrt(sumsSSIB2b.*sumsSILB2);

clear sumsPGSS sumsPGSI sumsSSSI sumsPGC2 sumsSSIB2 sumsSSIB2b sumsSILB2

        
        xcmax1=[];imax1=[];
        [xcmax1,imax1]=max(sumsPGSSn,[],2);
        
        xcmax2=[];imax2=[];
        [xcmax2,imax2]=max(sumsPGSIn,[],2);
        
        xcmax3=[];imax3=[];
        [xcmax3,imax3]=max(sumsSSSIn,[],2);
        
        
        imax1=imax1-mshift-1;
        imax2=imax2-mshift-1;
        imax3=imax3-mshift-1;
        imax=abs(imax1-imax2+imax3);
        xcmax=(xcmax1+xcmax2+xcmax3)/3;
        
        
        for jj=1:size(imax,1)
           if imax(jj)<=loopmax && xcmax(jj)>=xcmaxAVEnmin
               
                xcmaxcc=-999;PGSSimax=0;PGSIimax=0;xcmaxcc123=[];
               
                for iz1=max(imax1(jj)-3,-mshift):min(imax1(jj)+3,mshift)
                    for iz2=max(imax2(jj)-3,-mshift):min(imax2(jj)+3,mshift)
                        if iz2-iz1+mshift+1<=2*mshift+1 && iz2-iz1+mshift+1>=1
                         xcmaxcc1=sumsPGSSn(jj,iz1+mshift+1);
                         xcmaxcc2=sumsPGSIn(jj,iz2+mshift+1);
                         
                         
                         xcmaxcc3=sumsSSSIn(jj,iz2-iz1+mshift+1);
                              if xcmaxcc1+xcmaxcc2+xcmaxcc3>xcmaxcc
                                  xcmaxcc=xcmaxcc1+xcmaxcc2+xcmaxcc3;
                                    PGSSimax=iz1;
                                    PGSIimax=iz2;
                                   xcmaxcc123=[xcmaxcc1,xcmaxcc2,xcmaxcc3];
        
                
                              end
                        end
                    end
                end
               
        if xcmaxcc>=xcmaxAVEnmin*3
            zz=zz+1;
        
%           pgfft=[];ssfft=[];sifft=[];psssum=[];
%             
%             pgfft=fft(PGC(1+mshift:winlen+mshift),NFFT);
%             ssfft=fft(SSIB(1+mshift+PGSSimax:winlen+mshift+PGSSimax),NFFT);
%             sifft=fft(SILB(1+mshift+PGSIimax:winlen+mshift+PGSIimax),NFFT);
%    
%    psssum=1/max(abs(pgfft(1:NFFT/2+1)))*abs(pgfft(1:NFFT/2+1))+...
%        1/max(abs(ssfft(1:NFFT/2+1)))*abs(ssfft(1:NFFT/2+1))+...
%        1/max(abs(sifft(1:NFFT/2+1)))*abs(sifft(1:NFFT/2+1));
%    
%    
%    pp=[];qq=[];
%    [pp,qq]=max(psssum);  
            
            
            kiide(zz,:)=[z1+jj-1,(cxx(z1+jj-1)+mshift+1)/nsamp,PGSSimax,PGSIimax,imax(jj),xcmaxcc123,xcmaxcc/3,ii,jj];
           
            
        end
                
            
           end
            
        end

     ii   
    end
    
    











clear PGCEfd PGCNfd SSIBEfd SSIBNfd SILBEfd SILBNfd







            
     
 
            
            
       
        
        

        

% 

% %     
% %     if exist(['Japan/Kii/',JYr1,JMon1,JDay1])==0
% %         mkdir(['Japan/Kii/',JYr1,JMon1,JDay1])
% %     end
% %         
% % 
% % 
% % 
% % 
% % 
% % fname1=[];%fname2=[];
% % 
% % fname1=[fna,'kii',num2str(wr(1)),num2str(wr(2)),num2str(wr(3)),'_','kiide_100sps.mat'];
% % save(fname1,'kiide')
% % 















        
        
        
        
     njap 
end

    end
    
    
    
    
    
    
    
    
    
% %   
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % 
% % clear
% % 
% % %addpath /media/Elements/Japan/Kii/Hinet
% % 
% % 
% % Kii_time=load('Kii_ide_time_new');
% % 
% % 
% % %addpath ./
% % % 
% % %  wr1(:,1)=[3;4;6];
% % %  wr1(:,2)=[3;1;4];
% % %  wr1(:,3)=[4;5;6];
% % %  wr1(:,4)=[1;2;3];
% % %  
% % %   wr1(:,5)=[2;1;7];
% % %  wr1(:,6)=[6;4;8];
% % %  
% %  
% %  
% %  
% %  
% %  
% %  wr1(:,1)=[1;2;3];
% %  %wr1(:,2)=[4;5;6];
% % 
% %  %wr1(:,3)=[5;6;7];
% %  
% %  
% %  
% %  
% % 
% % 
% % 
% % sta=['URSH';'MASH';'MGWH'];
% % 
% % 
% % 
% % nsamp=100;
% % 
% % 
% % mshift=7*nsamp;
% % 
% %  %for win_length=1:2
% %  win_length=2;
% %     if win_length==1
% %        winlen=4*nsamp;
% % 
% %        winoff=nsamp/1;
% % 
% %        xcmaxAVEnmin=0.4%0.4 for 4s 1.5-6 Hz; 0.36 for 4s2-8 Hz ;  0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz
% %        fna1='/05_'
% %        
% %     else
% %        winlen=128*nsamp;
% % 
% %        winoff=nsamp*8;
% %        
% %        xcmaxAVEnmin=0.08%0.4 for 4s 1.5-6 Hz; 0.36 for 4s2-8 Hz ;  0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz
% %        
% %        loopmax=6;
% %        
% %        fna1='/05_128s_'
% %     end
% %     
% %  
% %     
% %     
% %     
% %     
% %     
% %    % timsPGC=0:0.025:(86400*2-0.025);
% %     
% %       timsPGC=0:1/nsamp:(86400*2-1/nsamp);
% % 
% % 
% % ss=0;
% % 
% %         
% %      wlength=86400*nsamp;   
% %      
% %      cxx=(3*nsamp:winoff:(wlength-mshift-3*nsamp-winlen))';
% %      
% %      cxx=cxx-mshift+15;
% %      
% % 
% % 
% % 
% %  nwr=1:1 
% %      wr=wr1(:,nwr)   
% %      
% %         
% % 
% % Fs=nsamp;
% % L=winlen;
% % NFFT=L;
% % 
% % f1=Fs/2*linspace(0,1,NFFT/2+1);
% % f2=Fs/2*linspace(1,0,NFFT/2+1);
% % f=[f1,-f2(2:NFFT/2)];
% % f=f';  
% % 
% % nff=f1';
% % nff(nff>6)=[];
% % zff2=size(nff,1);
% % nff(nff>1.5)=[];
% % zff1=size(nff,1);
% % 
% % 
% %         
% % 
% % 
% %     
% %     
% %     for njap=254:345 %2:size(Kii_time,1) %1
% %         
% %         Kii_time(njap)
% %         
% %         nyr=round(floor(Kii_time(njap)/10000));
% %         nmon=round(floor(Kii_time(njap)/100)-floor(Kii_time(njap)/10000)*100);
% %         nday=round(Kii_time(njap)-floor(Kii_time(njap)/100)*100);
% %         
% %     da=datevec(datenum(nyr,nmon,nday));
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
% %  
% %             
% %             
% %      
% %              %   prename1=['../array/rdseed/rdseed/mexico/',YEAR,'.',JDAY1,'.00.00.00.0000.TO'];
% % 
% %          %   prename1=['/media/Elements/Japan/Kii/Hinet/',JYr1,JMon1,JDay1,'/N'];
% %           prename1=['/data2/yajun/Japan/Kii/hinet/',JYr1,JMon1,JDay1,'/N'];
% %             
% %     fna=['Japan/Kii/',JYr1,JMon1,JDay1,'/'];
% %     
% %     
% %       
% % fname1=[];%fname2=[];
% % 
% % fname1=[fna,'kii',num2str(wr(1)),num2str(wr(2)),num2str(wr(3)),'_','kiide_100sps.mat'];
% % 
% % if exist(fname1)~=0;
% % 
% % load(fname1)
% % 
% % kiiden=[];
% % 
% % 
% % if size(kiide,2)~=0    
% %     
% % for narr1=1:3
% %         
% %     narr=wr(narr1);
% %     
% %         PGCEfdx=zeros(wlength,1);
% %         PGCNfdx=zeros(wlength,1);
% % 
% %                 
% %             
% % 
% %                 
% % PGCEdat1=[prename1,'.',sta(narr,:),'.E'];
% % PGCNdat1=[prename1,'.',sta(narr,:),'.N'];
% % 
% % 
% % if exist(PGCEdat1)~=0 
% % 
% % 
% % %if exist(PGCEdat1)~=0 && exist(PGCNdat1)~=0 && exist(PGCZdat1)~=0 
% % [PGCE1,HdrDataPGC1,tnuPGC1,pobjPGC1,timsPGC1]=readsac(PGCEdat1,0,'l');
% % [PGCN1,~,~,~,~]=readsac(PGCNdat1,0,'l');
% % %[PGCZ1,~,~,~,~]=readsac(PGCZdat1,0,'l');
% % clear timsPGC1
% % 
% % if size(PGCE1,1)~=100*86400
% %     PGCE1=[PGCE1;zeros(100*86400-size(PGCE1,1),1)];
% % end
% % 
% % if size(PGCN1,1)~=100*86400
% %     PGCN1=[PGCN1;zeros(100*86400-size(PGCN1,1),1)];
% % end
% % 
% % 
% % tracelenPGC1=length(PGCE1);
% % 
% % 
% % 
% %   
% % %cosine taper before filtering:
% % x=(0:pi/200:pi/2-pi/200)';
% % % %Seems to be necessary at the start of each day for PGCE:
% %     PGCE1(1:100)=sin(x).*PGCE1(1:100);
% % 
% %     PGCN1(1:100)=sin(x).*PGCN1(1:100);
% % 
% %     
% % %     nnz=100*86400;
% % %     
% % %     PGCE1(1+nnz:100+nnz)=sin(x).*PGCE1(1+nnz:100+nnz);
% % % 
% % %     PGCN1(1+nnz:100+nnz)=sin(x).*PGCN1(1+nnz:100+nnz);
% % %     
% % 
% % 
% % x=flipud(x);
% % PGCE1(tracelenPGC1-99:tracelenPGC1)=sin(x).*PGCE1(tracelenPGC1-99:tracelenPGC1);
% % PGCN1(tracelenPGC1-99:tracelenPGC1)=sin(x).*PGCN1(tracelenPGC1-99:tracelenPGC1);
% % 
% % % 
% % % PGCE1(tracelenPGC1-99-nnz:tracelenPGC1-nnz)=sin(x).*PGCE1(tracelenPGC1-99-nnz:tracelenPGC1-nnz);
% % % PGCN1(tracelenPGC1-99-nnz:tracelenPGC1-nnz)=sin(x).*PGCN1(tracelenPGC1-99-nnz:tracelenPGC1-nnz);
% % 
% % 
% % 
% % 
% % 
% % hi=6;
% % lo=1.5;
% % %Filter data:
% % npo=2;
% % npa=1;
% % 
% % % [PGCEf]=bandpass(PGCE1,100,lo,hi,npo,npa,'butter');
% % % [PGCNf]=bandpass(PGCN1,100,lo,hi,npo,npa,'butter');
% % % 
% % % 
% % % clear PGCE1 PGCN1% PGCZ1
% % % 
% % %  PGCEfdx = resample(PGCEf,2,5);
% % %  PGCNfdx = resample(PGCNf,2,5);
% % %  
% % % 
% % % 
% % % 
% % %  clear PGCEf PGCNf 
% %  
% %  
% %   [PGCEfdx]=bandpass(PGCE1,nsamp,lo,hi,npo,npa,'butter');
% % [PGCNfdx]=bandpass(PGCN1,nsamp,lo,hi,npo,npa,'butter');
% % 
% % 
% % clear PGCE1 PGCN1% PGCZ1
% %  
% % 
% %  
% %  
% % 
% %  
% % 
% % end
% % 
% %         
% %         
% % if narr1==1
% % 
% % PGCEfd=PGCEfdx;
% % 
% % PGCNfd=PGCNfdx; 
% % elseif narr1==2
% % 
% % SSIBEfd=PGCEfdx;
% % 
% % SSIBNfd=PGCNfdx; 
% % elseif narr1==3
% %     
% % SILBEfd=PGCEfdx;
% % 
% % SILBNfd=PGCNfdx; 
% % end
% % 
% % 
% % clear PGCEfdx PGCNfdx
% % end
% % 
% %         
% %          
% % 
% % 
% % % Fs=nsamp;
% % % L=winlen;
% % % NFFT=L;
% % % 
% % % f1=Fs/2*linspace(0,1,NFFT/2+1);
% % % f2=Fs/2*linspace(1,0,NFFT/2+1);
% % % f=[f1,-f2(2:NFFT/2)];
% % % f=f';    
% % %           
% %    
% % 
% %         ncc=0.1;
% % 
% % 
% %      zz=0;
% %      
% %      
% % if max(abs(PGCEfd))~=0 && max(abs(SSIBEfd))~=0 && max(abs(SILBEfd))~=0          
% %     for im=1:size(kiide,1) 
% %         if kiide(im,6)>=ncc/2 && kiide(im,7)>=ncc/2 && kiide(im,8)>=ncc/2 && kiide(im,9)>=ncc
% %         
% %         ii=kiide(im,1); 
% %         
% %         
% %        % niter=0;
% %    % while niter<=3
% %    
% %    
% %    
% %         narr=[];narr=wr(1);  
% %         pg_stack=zeros((winlen+2*mshift)*2,1);
% %         
% %         pg_stack=[PGCEfd(cxx(ii)+mshift+1-mshift:cxx(ii)+mshift+winlen+mshift);...
% %                 PGCNfd(cxx(ii)+mshift+1-mshift:cxx(ii)+mshift+winlen+mshift)];
% %         
% % 
% % 
% %         
% % 
% %         
% %    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% %         
% %          narr=[];narr=wr(2);  
% %         ss_stack=zeros((winlen+2*mshift)*2,1);
% %         
% %         ss_stack=[SSIBEfd(cxx(ii)+mshift+1-mshift:cxx(ii)+mshift+winlen+mshift);...
% %                 SSIBNfd(cxx(ii)+mshift+1-mshift:cxx(ii)+mshift+winlen+mshift)];
% %         
% % 
% %         
% %         
% % 
% % 
% %         
% %         
% %    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% %          narr=[];narr=wr(3);  
% %         si_stack=zeros((winlen+2*mshift)*2,1);
% %         
% %         si_stack=[SILBEfd(cxx(ii)+mshift+1-mshift:cxx(ii)+mshift+winlen+mshift);...
% %                 SILBNfd(cxx(ii)+mshift+1-mshift:cxx(ii)+mshift+winlen+mshift)];
% %         
% % 
% %         
% % 
% %         
% %         
% %         
% %         
% %         
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %         
% %            
% %         
% %         
% %         
% %        PGSSimax=kiide(im,3);
% %        PGSIimax=kiide(im,4);
% %        
% %                     
% %            % end
% %             
% %             
% %         %if xcmaxcc>=xcmaxAVEnmin*3
% %             
% %             
% %             xcmaxccn=-999;rot123=[];
% %             zz=zz+1;
% %             
% %             for rot1=40:10:210
% %                 pstn1=[];      
% %             pstn1=real((pg_stack(mshift+1:mshift+winlen)+1i*pg_stack(mshift+1+winlen+2*mshift:mshift+winlen+winlen+2*mshift))...
% %                 *exp(-1i*rot1/180*pi));  
% %             
% %                 
% %                 for rot2=40:10:210
% %                     
% %                                   
% %             pstn2=[];
% %             pstn2=real((ss_stack(mshift+1+PGSSimax:mshift+winlen+PGSSimax)+1i*ss_stack(mshift+1+winlen+2*mshift+PGSSimax:mshift+winlen+winlen+2*mshift+PGSSimax))...
% %                 *exp(-1i*rot2/180*pi));
% %                     
% %                     for rot3=40:10:210
% %             
% %                   
% %       
% %             
% %             pstn3=[];
% %             pstn3=real((si_stack(mshift+1+PGSIimax:mshift+winlen+PGSIimax)+1i*si_stack(mshift+1+winlen+2*mshift+PGSIimax:mshift+winlen+winlen+2*mshift+PGSIimax))...
% %                 *exp(-1i*rot3/180*pi));
% %             
% %             
% % 
% %             
% % %             
% %                         xcmaxcc1=[];%xcmaxcc2=[];xcmaxcc3=[];
% %             xcmaxcc1=sum(pstn1.*pstn2)+sum(pstn1.*pstn3)+sum(pstn3.*pstn2);
% %             
% % %             xcmaxcc2=sum(pstn1.*pstn3)/sqrt(sum(pstn1.^2)*sum(pstn3.^2));
% % %             xcmaxcc3=sum(pstn3.*pstn2)/sqrt(sum(pstn3.^2)*sum(pstn2.^2));
% %             
% %             if xcmaxcc1>xcmaxccn
% %                 xcmaxccn=xcmaxcc1;
% %                
% %                % xcmaxcc123=[xcmaxcc1,xcmaxcc2,xcmaxcc3];
% %                 rot123=[rot1,rot2,rot3];
% %                 
% %                 
% %             end
% %             
% %             
% %             
% %             
% %             
% %             
% %                     end
% %                 end
% %             end
% %         
% % 
% %             
% %             
% %               pstn1=[];      
% %             pstn1=real((pg_stack(mshift+1:mshift+winlen)+1i*pg_stack(mshift+1+winlen+2*mshift:mshift+winlen+winlen+2*mshift))...
% %                 *exp(-1i*rot123(1)/180*pi));        
% %                     
% %             pstn2=[];
% %             pstn2=real((ss_stack(mshift+1+PGSSimax:mshift+winlen+PGSSimax)+1i*ss_stack(mshift+1+winlen+2*mshift+PGSSimax:mshift+winlen+winlen+2*mshift+PGSSimax))...
% %                 *exp(-1i*rot123(2)/180*pi));
% %             
% %             pstn3=[];
% %             pstn3=real((si_stack(mshift+1+PGSIimax:mshift+winlen+PGSIimax)+1i*si_stack(mshift+1+winlen+2*mshift+PGSIimax:mshift+winlen+winlen+2*mshift+PGSIimax))...
% %                 *exp(-1i*rot123(3)/180*pi));
% %        
% %       
% %             pgfft=[];ssfft=[];sifft=[];
% %             
% %             pgfft=fft(pstn1,NFFT);
% %    ssfft=fft(pstn2,NFFT);
% %    sifft=fft(pstn3,NFFT);
% %    
% %    psssum=1/max(abs(pgfft(1:NFFT/2+1)))*abs(pgfft(1:NFFT/2+1))+...
% %        1/max(abs(ssfft(1:NFFT/2+1)))*abs(ssfft(1:NFFT/2+1))+...
% %        1/max(abs(sifft(1:NFFT/2+1)))*abs(sifft(1:NFFT/2+1));
% %    
% %    
% %    pp=[];qq=[];
% %    [pp,qq]=max(psssum);
% %            
% %         
% % 
% %         
% %         
% %          kiiden(zz,:)=[kiide(im,1:9),rot123,...
% %              sum(psssum(1:zff1-1))/sum(psssum),sum(psssum(zff1:zff2))/sum(psssum),sum(psssum(zff2:size(psssum,1)))/sum(psssum),f1(qq)];
% %         
% %         
% % 
% %     
% %    
% %     %im
% %         end
% %     end
% % 
% % end
% % 
% % 
% % 
% % 
% % 
% % 
% % clear PGCEfd PGCNfd SSIBEfd SSIBNfd SILBEfd SILBNfd
% % 
% % 
% % end
% % 
% % 
% % 
% % 
% %             
% %      
% %  
% %             
% %             
% %        
% %         
% %         
% % 
% % %         
% % % 
% % % 
% % 
% %     
% %     if exist(['Japan/Kii/',JYr1,JMon1,JDay1])==0
% %         mkdir(['Japan/Kii/',JYr1,JMon1,JDay1])
% %     end
% %         
% % 
% % 
% % 
% % 
% % 
% % fname1=[];%fname2=[];
% % 
% % fname1=[fna,'kii',num2str(wr(1)),num2str(wr(2)),num2str(wr(3)),'_','kiiden_100sps.mat'];
% % save(fname1,'kiiden')
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % end
% % 
% % 
% % 
% % 
% %         
% %         
% %         
% %         
% %      njap 
% % 
% % 
% %     end  
% %     
% %     
% %     
% %     
% %     
% %     
% %     
% %     
% %     
% %     