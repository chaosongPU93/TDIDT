function nzeros=glitches2(trace,nwin,winlen,winoff,igstart,spsfactor)
% In function glitches2, we're counting CONSECUTIVE zeros.

nzeros=zeros(nwin,1);
igstart=round((igstart-1)*spsfactor)+1;
winoff=round(winoff*spsfactor);
winlen=round(winlen*spsfactor);
% Could expand window, in case filtered glitch extends into this window from adjacent zeros.
for n=1:nwin
    istart=igstart+(n-1)*winoff;
    iend=istart+winlen;
    nzeros(n)=0;
    ii=istart-1;
    while ii<iend-1
        ii=ii+1;
        nz=0;
        while (abs(trace(ii)) <= 1.e-07 && ii<iend-1)
            nz=nz+1;
            ii=ii+1;
        end
        nzeros(n)=max(nzeros(n),nz);
%         if nzeros(n) >= 5
%             [n nzeros(n)]
%         end
    end
end
