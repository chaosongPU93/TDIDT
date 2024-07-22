function [ax,avg1,std1,mu1,sigma1,mdistcrssec,pks]=...
  plt_loccrssect(ax,F,crssecx,crssecy,angle,...
  color,linestyle,xran,which2plt,normalizer,gaussnum)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ax,mu,sigma,mdistcrssec]=plt_dloccrssect(ax,crssecx,crssecy,angle,...
%   color,normalizer,crssecy1,crssecy2)
%
% Very similar to to 'plt_dloccrssect.m', this function is to plot the 
% distribution of value on a cross section,
% a custom line defined by 'crssecx, crssecy' and its angle. Value on 
% the line is interpolated using a scattered Interpolant 'F'. 
% Distribution value can be normalized. Color of the line can be changed.
% Different from 'plt_dloccrssect.m' which was called in 'plt_srcdloc.m'
% and 'lfedloc_incluster.m', this function would plot EITHER the interpolated
% value, OR the Gaussian fit, determined by 'which2plt'. 
% 
% 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/04/21
% Last modified date:   2024/04/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('normalizer',1);
defval('gaussnum',1);

count = F(crssecx,crssecy);
dprojx = customprojection([crssecx crssecy],angle);
dprojx = dprojx(~isnan(count));
count = count(~isnan(count));
% countn = count/sum(count);  %essentially, normalized to 'probability'
countn = count;
%gaussian fit
if gaussnum == 1
  fit3par = fit(dprojx,countn,'gauss1','StartPoint',[1 0 1]);
  Y1=feval(fit3par,dprojx);
  coef1 = coeffvalues(fit3par);
  mu1 = coef1(2);
  sigma1 = coef1(3);
  avg1=wt_mean(dprojx,countn);
  std1=sqrt(wt_var(dprojx,countn));
  
%   fttpnorm = fittype( @(mu,sigma,x) 1/sigma/sqrt(2*pi)*exp(-0.5*(x-mu).^2/sigma^2) );
%   fit2par = fit(dprojx,count,fttpnorm,'StartPoint',[0 1]);
%   Y2=feval(fit2par,dprojx)*sum(count)*(dprojx(2)-dprojx(1));
%   coef2 = coeffvalues(fit2par);
%   mu2 = coef2(1);
%   sigma2 = coef2(2);
  
%   data=[];
%   for k =1:length(dprojx)
%     data = [data; dprojx(k)*ones(round(count(k)),1)];
%   end
%   [MUHAT1,SIGMAHAT1] = normfit(data);
%   Y2=normpdf(dprojx,MUHAT1,SIGMAHAT1)*size(data,1)*(dprojx(2)-dprojx(1));
elseif gaussnum == 2
  fitobj = fit(dprojx,countn,'gauss2','StartPoint',[10 -0.2 1 10 0.2 1]);
  coef1 = coeffvalues(fitobj);
end
[absdprojx, sumcount] = absxsumy(dprojx,countn);
% mdistcrssec = median_of_wt_data_pdf(absdprojx, sumcount);
mdistcrssec = wt_median(absdprojx, sumcount);
% keyboard

[pkhgt,locs] = findpeaks(countn,dprojx);
pks=[reshape(locs,[],1) reshape(pkhgt,[],1)];

if ~isempty(ax)
  hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  if strcmp(which2plt,'data')
    plot(ax,dprojx,countn / normalizer,linestyle,'Color',color,'linew',2);
  elseif strcmp(which2plt,'gsfit')
    plot(ax,dprojx,Y1 / normalizer,linestyle,'Color',color,'linew',2);
  end
  % text(ax,0.99,0.85,sprintf('\\mu=%.2f',mu1),'Units','normalized',...
  %   'HorizontalAlignment','right','FontSize',10);%,'color',color(2,:)
  % text(ax,0.99,0.78,sprintf('\\sigma=%.2f',sigma1),'Units','normalized',...
  %   'HorizontalAlignment','right','FontSize',10);%,'color',color(2,:)
  % text(ax,0.99,0.71,sprintf('med(|x|)=%.2f',mdistcrssec),'Units','normalized',...
  %   'HorizontalAlignment','right','FontSize',10);
  xlim(ax,xran);
  xlabel(ax,'Projected location (km)');
  if isequal(normalizer,1) && sum(countn)-1<1e-6
    ylabel(ax,'Probability');  
  else
    ylabel(ax,'Normalized count');  
  end
%   title(ax,sprintf('%d\\circ',round(angle)));
%   text(ax,0.99,0.95,strcat(num2str(round(angle)),'$^{\,\circ}$'),'interpreter',...
%     'latex','Units','normalized','FontSize',12,'HorizontalAlignment','right');
  longticks(ax,2);
%   xticks(ax,xran(1):1:xran(2));
%   ax.XMinorGrid = 'on';
%   ax.XMinorTick = 'on';
end

% keyboard