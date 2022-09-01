function nzeros=glitches(trace,nwin,winlen,winoff,igstart,spsfactor)
nzeros=zeros(nwin,1);
igstart=round((igstart-1)*spsfactor)+1;
winoff=round(winoff*spsfactor);
winlen=round(winlen*spsfactor);
for n=1:nwin
    istart=igstart+(n-1)*winoff;
    iend=istart+winlen;
    for i=istart:iend
        if abs(trace(i)) <= 1.e-08
            nzeros(n)=nzeros(n)+1;
        end
    end
end
