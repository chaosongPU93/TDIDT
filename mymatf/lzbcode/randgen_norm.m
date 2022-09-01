function [randnum,meanrand,pdfmisfit] = randgen_norm(seed,N,data,wt,pdfxloc,pdfval)
% [randnum,meanrand,pdfmisfit] = randgen_norm(seed,N,data,wt,pdfxloc,pdfval)
%
% Assuming the underlying PDF of data is Gaussian (normal)-like, this
% function is to first generate random numbers from a Gaussian fit to the 
% data 'data' with non-equal weights 'wt', with the same dimension
% size of the summed weights, by N times (realizations), calling function 
% 'utils/randn_wt'. 
% It takes in random seed state for reproduction.
% Then it evaluate the how well the PDF of generated random numbers matches
% that of data.
% Output the median of randoms and the misfit between randoms PDF and data
% PDF in each round of generation. 
%
% See also 'randgen_pdf', 'randgen_ksdensity', they do the same thing but utilizing
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

    %assume pdf is normal-like, and try to fit it with a gaussian 
    rng(seed{i});
    %return random numbers, and estimates of normal distribution parameters, mu and sigma
    [randnum(:,i),muHat,sigmaHat,~,~] = randn_wt(data, wt, dimsize);
%     some sense of the misfit
%     [pdfmisfit1, pdfmisfit2, f66] = compare_distribution2(randnum,muHat,sigmaHat,vdist,wt,pdfxloc,pdfval);    
    
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