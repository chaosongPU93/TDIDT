function [gauss]=gaussian(dt,nf,width)

%function to create a FREQ-domain Gaussian filter with sampling intervall (dt),
%number of points in frequency domain (nf) and width of the filter (width).
%Output is a Gaussian filter of form G=exp(-w^2/(4*width^2)).
% Note also that this should be the same as 
% '/mymatf/TremorDecon/gauss_zerophase.m'.


df=1/(dt*nf);
nft=0.5*nf+1;

f=df*(0:1:nft-1);
w=2*pi*f;

gauss=zeros(1,nf);
gauss(1:nf/2+1)=exp(-w.^2/(4.*width^2))/dt;
gauss(nf/2+2:end)=fliplr(gauss(2:nf/2));
gauss(nf/2+1)=0.;

end
