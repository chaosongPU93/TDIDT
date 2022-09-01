
clear

%addpath /media/Elements/Japan/Kii/Hinet


%Kii_time=load('Kii_ide_time_new');

Kii_time=load('Kii_JMA_time');

%Kii_time=20120411;

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
 %wr1(:,2)=[4;5;6];

 %wr1(:,3)=[5;6;7];
 
 




sta=['URSH';'MASH';'MGWH'];



nsamp=100;


mshift=50;

for win_length=1:2
% win_length=1;
    if win_length==1
       winlen=8*nsamp;

       winoff=nsamp/1;

       xcmaxAVEnmin=0.38%0.4 for 4s 1.5-6 Hz; 0.36 for 4s2-8 Hz ;  0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz
       fna1='kii_'
       
    else
       winlen=128*nsamp;

       winoff=nsamp*8;
       
       xcmaxAVEnmin=0.12%0.4 for 4s 1.5-6 Hz; 0.36 for 4s2-8 Hz ;  0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz
       
    
       
       fna1='kii_128s_'
    end
    
 
       loopmax=3;
    
      iup=2;
    
    
   % timsPGC=0:0.025:(86400*2-0.025);
    
      timsPGC=(0:1/nsamp:(86400*2-1/nsamp))';


ss=0;

        
     wlength=86400*nsamp;   
     
     cxx=(3*nsamp:winoff:(wlength-mshift-3*nsamp-winlen))';
     
     cxx=cxx-mshift+15;
     



 nwr=1:1 
     wr=wr1(:,nwr)   
     
%load(['Japan/Kii/','floc_kii',num2str(wr(1)),num2str(wr(2)),num2str(wr(3))]) 
 load(['Japan/Kii/','floc_kii',num2str(wr(1)),num2str(wr(2)),num2str(wr(3)),'_100sps_new'])
     




    
    
    for njap=1:size(Kii_time,1)                 % %%282:size(Kii_time,1)%% 8s
        
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
            
    fna=['Japan/Kii/',JYr1,JMon1,JDay1,'/'];
    
for narr1=1:3
        
    narr=wr(narr1);
    
        PGCEfdx=zeros(wlength,1);
        PGCNfdx=zeros(wlength,1);

                
            

                
PGCEdat1=[prename1,'.',sta(narr,:),'.E'];
PGCNdat1=[prename1,'.',sta(narr,:),'.N'];


if exist(PGCEdat1)~=0 


%if exist(PGCEdat1)~=0 && exist(PGCNdat1)~=0 && exist(PGCZdat1)~=0 
[PGCE1,HdrDataPGC1,tnuPGC1,pobjPGC1,timsPGC1]=readsac(PGCEdat1,0,'l');
[PGCN1,~,~,~,~]=readsac(PGCNdat1,0,'l');
%[PGCZ1,~,~,~,~]=readsac(PGCZdat1,0,'l');
clear timsPGC1

if size(PGCE1,1)~=100*86400
    PGCE1=[PGCE1;zeros(100*86400-size(PGCE1,1),1)];
end

if size(PGCN1,1)~=100*86400
    PGCN1=[PGCN1;zeros(100*86400-size(PGCN1,1),1)];
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
lo=1;
%Filter data:
npo=2;
npa=1;




%if nsamp==40

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
%end
 

 [PGCEfdx]=bandpass(PGCE1,nsamp,lo,hi,npo,npa,'butter');
[PGCNfdx]=bandpass(PGCN1,nsamp,lo,hi,npo,npa,'butter');


clear PGCE1 PGCN1% PGCZ1


 

 

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



Fs=nsamp;
L=winlen;
NFFT=L;

f1=Fs/2*linspace(0,1,NFFT/2+1);
f2=Fs/2*linspace(1,0,NFFT/2+1);
f=[f1,-f2(2:NFFT/2)];
f=f';  

nff=f1';
nff(nff>hi)=[];
zff2=size(nff,1);
nff(nff>lo)=[];
zff1=size(nff,1);     
        
        
for nfami=1:size(floc,1)
  
    
    
    
             rotPGC=floc(nfami,3);
             rotSSIB=floc(nfami,4);
             rotSILB=floc(nfami,5);
          
              
            
                PGCsoff=0;
                SSIBsoff=floc(nfami,1);
                SILBsoff=floc(nfami,2);

        
 
realPGC=real((PGCEfd+1i*PGCNfd)*exp(-1i*rotPGC/180*pi));

imgPGC=imag((PGCEfd+1i*PGCNfd)*exp(-1i*rotPGC/180*pi)); 


realSSIB=real((SSIBEfd+1i*SSIBNfd)*exp(-1i*rotSSIB/180*pi));

 imgSSIB=imag((SSIBEfd+1i*SSIBNfd)*exp(-1i*rotSSIB/180*pi));

    
    realSILB=real((SILBEfd+1i*SILBNfd)*exp(-1i*rotSILB/180*pi));

 imgSILB=imag((SILBEfd+1i*SILBNfd)*exp(-1i*rotSILB/180*pi));       

    
  %  
     
     
tracelenPGC=length(realPGC);
tracelenSSIB=length(realSSIB);
tracelenSILB=length(realSILB);


if PGCsoff > 0
    realPGC(1:tracelenPGC-PGCsoff)=realPGC(PGCsoff+1:tracelenPGC);
    realPGC(tracelenPGC-PGCsoff+1:tracelenPGC)=0;
    
    imgPGC(1:tracelenPGC-PGCsoff)=imgPGC(PGCsoff+1:tracelenPGC);
    imgPGC(tracelenPGC-PGCsoff+1:tracelenPGC)=0;
   
else
    realPGC(-PGCsoff+1:tracelenPGC)=realPGC(1:tracelenPGC+PGCsoff);
    realPGC(1:-PGCsoff)=0;
    
 imgPGC(-PGCsoff+1:tracelenPGC)=imgPGC(1:tracelenPGC+PGCsoff);
    imgPGC(1:-PGCsoff)=0;
end




if SSIBsoff > 0
    realSSIB(1:tracelenSSIB-SSIBsoff)=realSSIB(SSIBsoff+1:tracelenSSIB);
    realSSIB(tracelenSSIB-SSIBsoff+1:tracelenSSIB)=0;
    
    imgSSIB(1:tracelenSSIB-SSIBsoff)=imgSSIB(SSIBsoff+1:tracelenSSIB);
    imgSSIB(tracelenSSIB-SSIBsoff+1:tracelenSSIB)=0;
   
else
    realSSIB(-SSIBsoff+1:tracelenSSIB)=realSSIB(1:tracelenSSIB+SSIBsoff);
    realSSIB(1:-SSIBsoff)=0;
    
    imgSSIB(-SSIBsoff+1:tracelenSSIB)=imgSSIB(1:tracelenSSIB+SSIBsoff);
    imgSSIB(1:-SSIBsoff)=0;
   
end





if SILBsoff > 0
    realSILB(1:tracelenSILB-SILBsoff)=realSILB(SILBsoff+1:tracelenSILB);
    realSILB(tracelenSILB-SILBsoff+1:tracelenSILB)=0;
    
      imgSILB(1:tracelenSILB-SILBsoff)=imgSILB(SILBsoff+1:tracelenSILB);
    imgSILB(tracelenSILB-SILBsoff+1:tracelenSILB)=0;
   
else
    realSILB(-SILBsoff+1:tracelenSILB)=realSILB(1:tracelenSILB+SILBsoff);
    realSILB(1:-SILBsoff)=0;
    
   imgSILB(-SILBsoff+1:tracelenSILB)=imgSILB(1:tracelenSILB+SILBsoff);
    imgSILB(1:-SILBsoff)=0;
   
end


    






        


     zz=0;PGSSSImap=[];PGSSSIwa=[];nsect=4;
     
      
     
     
     
     
if max(abs(realPGC))~=0 && max(abs(realSSIB))~=0 && max(abs(realSILB))~=0 
    nslen=ceil((size(cxx,1)-2)/nsect);
    
    for ii=1:nsect
        
        
        
        z1=(ii-1)*nslen+1+1;
        z2=min(ii*nslen+1,size(cxx,1)-1);
        
        PGC=[];
        PGC=realPGC(cxx(z1)+mshift+1-mshift:cxx(z2)+mshift+winlen+mshift);
            
         SSIB=[];
        SSIB=realSSIB(cxx(z1)+mshift+1-mshift:cxx(z2)+mshift+winlen+mshift);
        
         SILB=[];
        SILB=realSILB(cxx(z1)+mshift+1-mshift:cxx(z2)+mshift+winlen+mshift);
        
            
            PGC2=cumsum(PGC(:,1).^2,1);
            PGC2=[zeros(1,size(PGC2,2));PGC2];
            SSIB2=cumsum(SSIB(:,1).^2,1);
             SSIB2=[zeros(1,size(SSIB2,2));SSIB2];
             SILB2=cumsum(SILB(:,1).^2,1);
              SILB2=[zeros(1,size(SILB2,2));SILB2];
              
              
              
     lenx=size(PGC,1)-2*mshift;        
PGSSx=zeros(lenx, 2*mshift+1);
PGSIx=zeros(lenx, 2*mshift+1);
SSSIx=zeros(lenx, 2*mshift+1);
for n=-mshift:mshift;
    PGSSx(:,n+mshift+1)=PGC(1+mshift:size(PGC,1)-mshift,1).* ...
        SSIB(1+mshift+n:size(PGC,1)-mshift+n,1);
    
    
    PGSIx(:,n+mshift+1)=PGC(1+mshift:size(PGC,1)-mshift,1).* ...
        SILB(1+mshift+n:size(PGC,1)-mshift+n,1);
    
    
    SSSIx(:,n+mshift+1)=SSIB(1+mshift:size(PGC,1)-mshift,1).* ...
        SILB(1+mshift+n:size(PGC,1)-mshift+n,1);
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
        
        
        
[xmaxPGSSn,ymaxPGSSn,aPGSS]=parabol(nwin,mshift,sumsPGSSn,imax1+mshift+1); %Interpolated max cross-correlation
[xmaxPGSIn,ymaxPGSIn,aPGSI]=parabol(nwin,mshift,sumsPGSIn,imax2+mshift+1);
[xmaxSSSIn,ymaxSSSIn,aSSSI]=parabol(nwin,mshift,sumsSSSIn,imax3+mshift+1);
    

xmaxPGSSn=xmaxPGSSn-mshift-1;
xmaxPGSIn=xmaxPGSIn-mshift-1;
xmaxSSSIn=xmaxSSSIn-mshift-1;
 loopoff=xmaxPGSSn-xmaxPGSIn+xmaxSSSIn;  
 
 
      
        
        for jj=1:size(loopoff,1) %1:size(imax,1)
           if abs(loopoff(jj))<=loopmax && xcmax(jj)>=xcmaxAVEnmin ...
                   && abs(imax1(jj))~=mshift && abs(imax2(jj))~=mshift && abs(imax3(jj))~=mshift
               
               
        interpPGSSn=interp(sumsPGSSn(jj,:),iup,3);
        interpPGSIn=interp(sumsPGSIn(jj,:),iup,3);
        interpSSSIn=interp(sumsSSSIn(jj,:),iup,3);
        
        
                leninterp=length(interpPGSSn);
        [xcmaxinterpPGSSn,imaxinterpPGSS]=max(interpPGSSn(1:leninterp-(iup-1)));
        [xcmaxinterpPGSIn,imaxinterpPGSI]=max(interpPGSIn(1:leninterp-(iup-1)));
        [xcmaxinterpSSSIn,imaxinterpSSSI]=max(interpSSSIn(1:leninterp-(iup-1)));
               
                xcmaxcc=-999;PGSSimax=0;PGSIimax=0;xcmaxcc123=[];
               
                for iz1=max(1,imaxinterpPGSS-3*iup):min(imaxinterpPGSS+3*iup,iup*(2*mshift+1)-(iup-1)) %max(imax1(jj)-3,-mshift):min(imax1(jj)+3,mshift)
                    for iz2=max(1,imaxinterpPGSI-3*iup):min(imaxinterpPGSI+3*iup,iup*(2*mshift+1)-(iup-1))  %max(imax2(jj)-3,-mshift):min(imax2(jj)+3,mshift)
                        
                         iz3 = (iup*mshift+1)+iz2-iz1;
                        
                        if iz3<=iup*(2*mshift+1) && iz3>=1 %iz2-iz1+mshift+1<=2*mshift+1 && iz2-iz1+mshift+1>=1
                         xcmaxcc1=interpPGSSn(iz1);
                         xcmaxcc2=interpPGSIn(iz2);
                         
                         
                         xcmaxcc3=interpSSSIn(iz3);
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
           
        
            PGSSimax=(PGSSimax-(iup*mshift+1))/iup;
            PGSIimax=(PGSIimax-(iup*mshift+1))/iup;
            
            
            pgfft=[];ssfft=[];sifft=[];
            
            pgfft=fft(realPGC(cxx(z1+jj-1)+mshift+1:cxx(z1+jj-1)+mshift+winlen),NFFT);
   ssfft=fft(realSSIB(cxx(z1+jj-1)+mshift+1+round(PGSSimax):cxx(z1+jj-1)+mshift+winlen+round(PGSSimax)),NFFT);
   sifft=fft(realSILB(cxx(z1+jj-1)+mshift+1+round(PGSIimax):cxx(z1+jj-1)+mshift+winlen+round(PGSIimax)),NFFT);
   
   psssum=1/max(abs(pgfft(1:NFFT/2+1)))*abs(pgfft(1:NFFT/2+1))+...
       1/max(abs(ssfft(1:NFFT/2+1)))*abs(ssfft(1:NFFT/2+1))+...
       1/max(abs(sifft(1:NFFT/2+1)))*abs(sifft(1:NFFT/2+1));
   
   
   pp=[];qq=[];
   [pp,qq]=max(psssum);
           
        
            
    if f1(qq)>lo && f1(qq)<hi && sum(psssum(1:zff1-1))/sum(psssum)<0.2 && sum(psssum(zff1:zff2))/sum(psssum)>0.6
        
         zz=zz+1;
            PGSSSImap(zz,:)=[z1+jj-1,(cxx(z1+jj-1)+mshift+1)/nsamp-1/nsamp,PGSSimax,PGSIimax,loopoff(jj),imax(jj),xcmaxcc123,xcmaxcc/3,...
                sum(psssum(1:zff1-1))/sum(psssum),sum(psssum(zff1:zff2))/sum(psssum),sum(psssum(zff2:size(psssum,1)))/sum(psssum),f1(qq),...
                ii,jj];
           
            if win_length==1
            
              PGSSSIwa((zz-1)*winlen+1:zz*winlen,:)=[timsPGC(cxx(z1+jj-1)+mshift+1:cxx(z1+jj-1)+mshift+winlen),...
                  realPGC(cxx(z1+jj-1)+mshift+1:cxx(z1+jj-1)+mshift+winlen),...
                  realSSIB(cxx(z1+jj-1)+mshift+1+round(PGSSimax):cxx(z1+jj-1)+mshift+winlen+round(PGSSimax)),...
                  realSILB(cxx(z1+jj-1)+mshift+1+round(PGSIimax):cxx(z1+jj-1)+mshift+winlen+round(PGSIimax)),...
                  imgPGC(cxx(z1+jj-1)+mshift+1:cxx(z1+jj-1)+mshift+winlen),...
                  imgSSIB(cxx(z1+jj-1)+mshift+1+round(PGSSimax):cxx(z1+jj-1)+mshift+winlen+round(PGSSimax)),...
                  imgSILB(cxx(z1+jj-1)+mshift+1+round(PGSIimax):cxx(z1+jj-1)+mshift+winlen+round(PGSIimax))];
            end
                 
    end
            
        end
                
            
           end
            
        end

    % ii   
    end
    
    
    if exist(['Japan/Kii/',JYr1,JMon1,JDay1])==0
        mkdir(['Japan/Kii/',JYr1,JMon1,JDay1])
    end
        

fname1=[];
fname1=[fna,fna1,num2str(nfami),'_',num2str(wr(1)),num2str(wr(2)),num2str(wr(3)),'PGSSSImap.mat'];
 save(fname1,'PGSSSImap')

 if win_length==1
 
fname2=[];
fname2=[fna,fna1,num2str(nfami),'_',num2str(wr(1)),num2str(wr(2)),num2str(wr(3)),'PGSSSIwa.mat'];
 save(fname2,'PGSSSIwa')
 end



            








        
        
        
        
     
end

nfami
end


 clear PGCEfd PGCNfd SSIBEfd SSIBNfd SILBEfd SILBNfd
 
njap 
    end
end
    
%     fname1=[];
% fname1=[fna,fna1,num2str(nfami),'_',num2str(wr(1)),num2str(wr(2)),num2str(wr(3)),'PGSSSImap_1to6Hz.mat'];
%  save(fname1,'PGSSSImap')
    
 %figure;plot(mapfile(:,2),mapfile(:,4),'.')   
    
    