%Reads a Bostock file, converts mags to moment using ... , sums moments and
%writes it out including date.
%246  30226  6 3518.450 1.387  9
%  2  30226 11 1151.450 1.404 12
%  2  30226 22 3204.000 1.290  4
%  2  30228  6 2473.625 2.007  7
magfile=load('BOSTOCK/NEW/002-246_culled.mags');
afile(:,1)=magfile(:,2);
afile(:,2)=magfile(:,5);
afile(:,3)=10.^(1.5*(6+magfile(:,5)));
afile(:,4)=cumsum(afile(:,3));
fid=fopen('BFIGS/cummom','w');
fprintf(fid,'%6i %6.3f %12.5e %12.5e\n',afile');
fclose(fid);
