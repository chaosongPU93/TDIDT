function color = kelicmap(n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% color = kelicmap(n)
%
% Generate n colors following the colormap of 'kelicol'.
% 
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2023/09/19
% Last modified date:   2023/09/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colorraw = kelicol; %first generate denser colors

ncolor = size(colorraw,1);

if n >= ncolor
  color = colorraw;
else
  ind = round(linspace(1,ncolor,n));
  color = colorraw(ind,:);
end

