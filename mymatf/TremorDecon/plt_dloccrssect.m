function [ax,mu,sigma,mdistcrssec]=plt_dloccrssect(ax,F,crssecx,crssecy,angle,...
  color,xran,normalizer,crssecy1,crssecy2,gaussnum)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ax,mu,sigma,mdistcrssec]=plt_dloccrssect(ax,crssecx,crssecy,angle,...
%   color,normalizer,crssecy1,crssecy2)
%
% This function is to plot the distribution of value on a cross section,
% a custom line defined by 'crssecx, crssecy' and its angle. Value on 
% the line is interpolated using a scattered Interpolant 'F'. 
% Distribution value can be normalized. Color of the line can be changed.
% able to plot other lines as well
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/03/05
% Last modified date:   2024/03/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('normalizer',1);
defval('crssecy1',[]);
defval('crssecy2',[]);
defval('gaussnum',1);

hold(ax,'on'); ax.Box='on'; grid(ax,'on');
count = F(crssecx,crssecy);
dprojx = customprojection([crssecx crssecy],angle);
dprojx = dprojx(~isnan(count));
count = count(~isnan(count));
plot(ax,dprojx,count / normalizer,'-','Color',color(1,:),'linew',2);
%gaussian fit
if gaussnum == 1
  fitobj = fit(dprojx,count,'gauss1','StartPoint',[10 0 1]);
  Y1=feval(fitobj,dprojx);
  coef = coeffvalues(fitobj);
  mu = coef(2);
  sigma = coef(3);
elseif gaussnum == 2
  fitobj = fit(dprojx,count,'gauss2','StartPoint',[10 -0.2 1 10 0.2 1]);
  coef = coeffvalues(fitobj);
end
plot(ax,dprojx,Y1 / normalizer,'-','Color',color(2,:),'linew',2);
[absdprojx, sumcount] = absxsumy(dprojx,count);
mdistcrssec = median_of_wt_data_pdf(absdprojx, sumcount);
if ~isempty(crssecy1)
  countopt1 = F(crssecx,crssecy1);
  dprojx = customprojection([crssecx crssecy1],angle);
  plot(ax,dprojx,countopt1 / normalizer,'--','Color',color(1,:),'linew',1);
end
if ~isempty(crssecy2)
  countopt2 = F(crssecx,crssecy2);
  dprojx = customprojection([crssecx crssecy2],angle);
  plot(ax,dprojx,countopt2 / normalizer,'-.','Color',color(1,:),'linew',1);
end
text(ax,0.02,0.95,sprintf('\\mu=%.2f, \\sigma=%.2f, med(|x|)=%.2f',...
  mu,sigma,mdistcrssec),...
  'Units','normalized','HorizontalAlignment','left','color',color(2,:));
xlim(ax,xran);
xlabel(ax,'Location difference (km)');
ylabel(ax,'Normalized count');  
title(ax,sprintf('%d ^o',round(angle)));
longticks(ax,2);

