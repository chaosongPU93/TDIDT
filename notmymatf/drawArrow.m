function [ h ] = drawArrow(ax,x,y,xlimits,ylimits,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to plot arrows from x(1),y(1) to x(2),y(2)
% 
% from https://stackoverflow.com/questions/25729784/how-to-draw-an-arrow-in-matlab
% 
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2019/08/22
% Last modified date:   2019/08/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlim(ax,xlimits)
ylim(ax,ylimits)

h = annotation('arrow');
set(h,'parent', ax, ...
    'position', [x(1),y(1),x(2)-x(1),y(2)-y(1)], ...
    'HeadLength', 8, 'HeadWidth', 5, 'HeadStyle', 'cback1', ...
    varargin{:} );

end