function [randnum,meanrand,pdfmisfit,exppdf] = randgen_pdf(seed,N,data,wt,pdfxloc,pdfval)
% [randnum,meanrand,pdfmisfit,exppdf] = randgen_pdf(seed,N,data,wt,pdfxloc,pdfval)
%
% Taking the underlying PDF of data as custom, this function is to 
% first generate random numbers to the data 'data' with non-equal
% weights 'wt', with the same dimension size of the summed weights,
% by N times (realizations), calling function 'notmymatf/randpdf'. 
% It takes in random seed state for reproduction.
% Then it evaluate the how well the PDF of generated random numbers matches
% that of data.
% Output the median of randoms and the misfit between randoms PDF and data
% PDF in each round of generation. 
%
% See also 'randgen_norm', 'randgen_ksdensity', they do the same thing but utilizing
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
    
    % 1. use randpdf
    %%% NOTE:
    %%% the shifting arises bc. the obtained pdf requires certain binning, so if you want to shift
    %%% the data by amount 'shift' by directly add it to the data, it is fine; but if the bin width
    %%% and edge remains unchanged, the histogram (pdf) usually would NOT have a same shape as the
    %%% unshifted histogram (pdf), so preserve the shape, you MUST add the shift to the pdf x
    %%% location, rather than do binning again to the shifted new data
    rng(seed{i});
    %return random numbers and expectation of mean of the input pdf
    [randnum(:,i), exppdf] = randpdf(pdfxloc, pdfval, dimsize);
    
    %%% output the misfits of PDFs as well but no figures
    for ii = 1: size(randnum(:,i),2)
        binw = pdfxloc(2)-pdfxloc(1);
        refval = exppdf;
        pdfmisfit(i) = check_randgen(randnum(:,i),pdfxloc,pdfval,refval,binw);
    end
    
end

meanrand = mean(randnum,1);   % mean of the generated random numbers
meanrand = meanrand';

% keyboard
