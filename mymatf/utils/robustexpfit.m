function fitstruct=robustexpfit(xdata,ydata,logoption,startpt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fitstruct=robustlnfit(xdata,ydata)
%
% This function is to apply a bi-square robust linear fit to data 'x' and 'y'.
% Return the result in a structure that contains the 'fitobj','gof', 'output',
% and 'stats' that summarizes some statistics of the fit. 
% Related functions are 'fit', 'fittype', and 'statsofrobustlnfit'.
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/01/18
% Last modified date:   2024/01/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('logoption','log');
defval('startpt',[.1 .1]);


xdata=reshape(xdata,[],1);
ydata=reshape(ydata,[],1);

if strcmp(logoption,'log10')
  fttpfree = fittype( @(a,b,x) a*10.^(b*x) );
elseif strcmp(logoption,'log')
  fttpfree = fittype( @(a,b,x) a*exp(b*x) );
end
[fitobj,gof,output] = fit(xdata,ydata,fttpfree,'Robust','Bisquare',...
  'StartPoint',startpt);

fitstruct.fitobj = fitobj;
fitstruct.gof = gof;
fitstruct.output = output;
