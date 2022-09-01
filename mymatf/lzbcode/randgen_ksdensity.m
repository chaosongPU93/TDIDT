function [randnum,meanrand,pdfmisfit] = randgen_ksdensity(seed,N,data,wt,pdfxloc,pdfval)
% [randnum,meanrand,pdfmisfit] = randgen_ksdensity(seed,N,data,wt,pdfxloc,pdfval)
%
% Using kernel smoothing function estimate for univariate and bivariate data, 
% this function is to first generate random numbers from a Gaussian fit to 
% the data 'data' with non-equal weights 'wt', with the same dimension
% size of the summed weights, by N times (realizations), calling a built-in  
% function 'ksdensity'. 
% It takes in random seed state for reproduction.
% Then it evaluate the how well the PDF of generated random numbers matches
% that of data.
% Output the median of randoms and the misfit between randoms PDF and data
% PDF in each round of generation. 
%
% See also 'randgen_pdf', 'randgen_norm', they do the same thing but utilizing
% different methods.
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/10/27
% Last modified date:   2022/04/06

pdfmisfit = zeros(N, 1);    % misfits of PDFs between data and generated random numbers
dimsize = [ceil(sum(wt)),1];  % size of desired randum numbers per sampling
randnum = zeros(ceil(sum(wt)), N);
for i = 1: N
% i = randi(1000,1);

    %use 'ksdensity', kernel smoothing function estimate for univariate and bivariate data
    rng(seed{i});
    %return random numbers
    randnum(:,i) = ksdensity(data, rand(dimsize), 'function', 'icdf', 'weight', wt);
    
    %%% output the misfits of PDFs as well but no figures
    for ii = 1: size(randnum(:,i),2)
        binw = pdfxloc(2)-pdfxloc(1);
        refval = wt_mean(data, wt);
        pdfmisfit(i) = check_randgen(randnum(:,i),pdfxloc,pdfval,refval,binw);
    end
    
end

meanrand = mean(randnum,1);   % mean of the generated random numbers
meanrand = meanrand';

% keyboard


