function STAstretch=stretchtemps(STAt,minstretch,maxstretch,winlen,sps) %sps added 11/2/2017
STAstretch=zeros(maxstretch-minstretch+1,winlen);
for istretch=minstretch:maxstretch
    STAdum=resample(STAt,istretch,sps); %Downsample to minstretch...maxstretch sps
    [~,imaxSTA]=max(STAdum); 
    [~,iminSTA]=min(STAdum); 
    for i=iminSTA:imaxSTA-1
        if STAdum(i)*STAdum(i+1)<0
            if abs(STAdum(i))<abs(STAdum(i+1)) %zero-crossing
                izero=i;
            else
                izero=i+1;
            end
        end
    end
    if izero-winlen/2 >= 1
        STAdum=STAdum(izero-winlen/2:izero+winlen/2-1); %STAdum now centered on the zero crossing.
    else
        STAdum(1:winlen/2-izero+1)=0;
        STAdum(winlen-(winlen/2-izero):winlen)=0;
        STAdum(winlen/2-izero+2:izero+winlen/2-1)=STAdum(1:2*izero-2); %STAdum now centered on the zero crossing.
    end
    [~,imax]=max(STAdum); 
    [~,imin]=min(STAdum); 
    i=imin;
    while STAdum(i)<0 %assumes first excursion is negative (broadband, 2-pass filter)
        i=i-1;
    end
    ikeepmin=i+1;
    i=imax;
    while STAdum(i)>0 %assumes 2nd excursion is positive (broadband, 2-pass filter)
        i=i+1;
    end
    ikeepmax=i-1;
    STAstretch(istretch-minstretch+1,ikeepmin:ikeepmax)=STAdum(ikeepmin:ikeepmax); %Still centered on the zero crossing, but samples outside of dipole are zero.
end
