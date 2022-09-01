% testfitting

% Create a baseline sinusoidal signal:
xdata = (0:0.1:2*pi)'; 
y0 = sin(xdata);

% Add noise to the signal with nonconstant variance.
% Response-dependent Gaussian noise
gnoise = y0.*randn(size(y0));

% Salt-and-pepper noise
spnoise = zeros(size(y0)); 
p = randperm(length(y0));
sppoints = p(1:round(length(p)/5));
spnoise(sppoints) = 5*sign(y0(sppoints));

ydata = y0 + gnoise + spnoise;

% Fit the noisy data with a baseline sinusoidal model, and specify 3 output arguments to get fitting
% information including residuals.
f = fittype('a*sin(b*x)'); 
[fit1,gof,fitinfo] = fit(xdata,ydata,f,'StartPoint',[1 1]);

% Examine the information in the fitinfo structure.

fitinfo

% Get the residuals from the fitinfo structure.

residuals = fitinfo.residuals;

% Identify "outliers" as points at an arbitrary distance greater than 1.5 standard deviations from
% the baseline model, and refit the data with the outliers excluded.

I = abs( residuals) > 1.5 * std( residuals );
outliers = excludedata(xdata,ydata,'indices',I);

fit2 = fit(xdata,ydata,f,'StartPoint',[1 1],...
           'Exclude',outliers);

% Compare the effect of excluding the outliers with the effect of giving them lower bisquare weight
% in a robust fit.

fit3 = fit(xdata,ydata,f,'StartPoint',[1 1],'Robust','on');

% Plot the data, the outliers, and the results of the fits. Specify an informative legend.
figure 
plot(fit1,'r-',xdata,ydata,'k.',outliers,'m*') 
hold on
plot(fit2,'c--')
plot(fit3,'b:')
xlim([0 2*pi])
legend( 'Data', 'Data excluded from second fit', 'Original fit',...
    'Fit with points excluded', 'Robust fit' )
hold off

%Plot the residuals for the two fits considering outliers:
figure 
plot(fit2,xdata,ydata,'co','residuals') 
hold on
plot(fit3,xdata,ydata,'bx','residuals')
hold off


    %%% test section
%     fttpfree = fittype( @(a,b,x) a*x+b,'Robust','Bisquare');
    fttpfree = fittype( @(a,b,x) a*x+b);
    [fittest1,goftest1,outptest1] = fit(mighfdum(:,end)/3600, mighfdum(:,1),fttpfree,'Robust','Bisquare');
    testval1 = feval(fittest1,mighfdum(:,end)/3600);
    plot(f.ax(3),mighfdum(:,end)/3600,testval1,'k--','linewidth',2);
    
    %     fttpfix = fittype( @(a,b,x) a*x+b,'Robust','Bisquare','problem','b');
    a=slopeprophf;
    fttpfix = fittype( @(b,x) a*x+b);
    [fittest2,goftest2,outptest2] = fit(miglfdum(:,end)/3600, miglfdum(:,1),fttpfix,'Robust','Bisquare');
    testval2 = feval(fittest2,miglfdum(:,end)/3600);
    plot(f.ax(3),miglfdum(:,end)/3600,testval2,'-','linewidth',2,'color',[153/255 255/255 255/255]);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    