function [map,num] = tab20(N)
% Qualitative colormap from MatPlotLib, for plots using the line ColorOrder.
% In MatPlotLib 2 it is named VEGA20, for MatPlotLib 3 was renamed TAB20.
%
% Copyright (c) 2017-2024 Stephen Cobeldick
%
%%% Syntax:
%  map = tab20
%  map = tab20(N)
%
% For MatPlotLib 2.0 improved colormaps were created for plot lines of
% categorical data. The new colormaps are introduced here:
% <http://matplotlib.org/2.0.0rc2/users/dflt_style_changes.html>
% VEGA10/TAB10 is the default Line Color Order for MatPlotLib 2 and 3.
%
% MATLAB axes ColorOrder (note that this is NOT the axes COLORMAP):
% <https://www.mathworks.com/help/matlab/creating_plots/defining-the-color-of-lines-for-plotting.html>
%
%% Examples %%
%
%%% PLOT using matrices:
% N = 20;
% axes('ColorOrder',tab20(N),'NextPlot','replacechildren')
% X = linspace(0,pi*3,1000);
% Y = bsxfun(@(x,n)sqrt(n)*sin(x+2*n*pi/N), X(:), 1:N);
% plot(X,Y, 'linewidth',4)
%
%%% PLOT in a loop:
% N = 20;
% set(0,'DefaultAxesColorOrder',tab20(N))
% X = linspace(0,pi*3,1000);
% Y = bsxfun(@(x,n)sqrt(n)*sin(x+2*n*pi/N), X(:), 1:N);
% for n = 1:N
%     plot(X(:),Y(:,n), 'linewidth',4);
%     hold all
% end
%
%%% LINE using matrices:
% N = 20;
% set(0,'DefaultAxesColorOrder',tab20(N))
% X = linspace(0,pi*3,1000);
% Y = bsxfun(@(x,n)sqrt(n)*cos(x+2*n*pi/N), X(:), 1:N);
% line(X(:),Y)
%
%% Input and Output Arguments %%
%
%%% Inputs (**=default):
% N = NumericScalar, N>=0, an integer to define the colormap length.
%   = []**, same length as MATLAB's inbuilt colormap functions.
%
%%% Outputs:
% map = NumericMatrix, size Nx3, a colormap of RGB values between 0 and 1.
%
% See also TAB10 TAB20B TAB20C
% CIVIDIS INFERNO MAGMA PLASMA TWILIGHT TWLIGHT_SHIFTED VIRIDIS
% MAXDISTCOLOR BREWERMAP CUBEHELIX CMOCEAN LBMAP LINES PARULA COLORMAP SET

if nargin<1 || isnumeric(N)&&isequal(N,[])
	N = cmDefaultN();
else
	assert(isscalar(N)&&isfinite(N)&&isreal(N),...
		'SC:tab20:N:NotRealFiniteScalarNumeric',...
		'First argument <N> must be a real finite numeric scalar.')
end
%
hex = ['#1f77b4';'#aec7e8';'#ff7f0e';'#ffbb78';'#2ca02c';'#98df8a';'#d62728';'#ff9896';'#9467bd';'#c5b0d5';'#8c564b';'#c49c94';'#e377c2';'#f7b6d2';'#7f7f7f';'#c7c7c7';'#bcbd22';'#dbdb8d';'#17becf';'#9edae5'];
raw = sscanf(hex.','#%2x%2x%2x',[3,Inf]).';
num = size(raw,1);
%
map = raw(1+mod(0:N-1,num),:) / 255;
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%tab20
function N = cmDefaultN()
% Get the colormap size from the current figure or default colormap.
try
	F = get(groot,'CurrentFigure');
catch %#ok<CTCH> pre HG2
	N = size(get(gcf,'colormap'),1);
	return
end
if isempty(F)
	N = size(get(groot,'DefaultFigureColormap'),1);
else
	N = size(F.Colormap,1);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%cmDefaultN