function [ax,mu1,sigma1,mdistcrssec,dprojx,countn,Y1,pks]=plt_dloccrssect(ax,F,crssecx,crssecy,angle,...
  linetype,color,xran,which2plt,normalizer,gaussnum)
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
defval('which2plt','both');
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
  fit3par = fit(dprojx,countn,'gauss1','StartPoint',[10 0 1]);
  Y1=feval(fit3par,dprojx);
  coef1 = coeffvalues(fit3par);
  mu1 = coef1(2);
  sigma1 = coef1(3);
  
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

[pkhgt,locs] = findpeaks(countn,dprojx);
pks=[reshape(locs,[],1) reshape(pkhgt,[],1)];

if ~isempty(ax)
  hold(ax,'on'); ax.Box='on'; grid(ax,'on');
  if strcmp(which2plt,'data')
    plot(ax,dprojx,countn / normalizer,linetype(1,:),'Color',color(1,:),'linew',2);
  elseif strcmp(which2plt,'gsfit')
    plot(ax,dprojx,Y1 / normalizer,linetype(2,:),'Color',color(2,:),'linew',2);
  elseif strcmp(which2plt,'both')
    plot(ax,dprojx,countn / normalizer,linetype(1,:),'Color',color(1,:),'linew',2);
    plot(ax,dprojx,Y1 / normalizer,linetype(2,:),'Color',color(2,:),'linew',2);
  end
%   plot(ax,dprojx,Y2 / normalizer,'-','Color',color(2,:),'linew',1);
%   scatter(ax,pks(:,1),pks(:,2) / normalizer,8,'ko','filled');  
%   text(ax,0.99,0.85,sprintf('\\mu=%.2f',mu1),'Units','normalized',...
%     'HorizontalAlignment','right','FontSize',10);%,'color',color(2,:)
%   text(ax,0.99,0.78,sprintf('\\sigma=%.2f',sigma1),'Units','normalized',...
%     'HorizontalAlignment','right','FontSize',10);%,'color',color(2,:)
%   text(ax,0.99,0.71,sprintf('med(|x|)=%.2f',mdistcrssec),'Units','normalized',...
%     'HorizontalAlignment','right','FontSize',10);
  xlim(ax,xran);
  xlabel(ax,'Location difference (km)');
  % keyboard
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