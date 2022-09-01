%test
npass=0;
lo = 1;
hi = 5;
for istack= 1: nstack
windata = zeros(1, templen+ 2*extra);
windata(:) = datamat(1, istack, :);
template = stackex(1, :)/ nstack;

fwindata = Bandpass(windata, sps, lo, hi, npo, npa, 'butter');
ftemplate = Bandpass(template, sps, lo, hi, npo, npa, 'butter');
d1 = fwindata(templen/2+extra-4*sps+1:templen/2+extra+4*sps);
t1 = ftemplate(templen/2+extra-4*sps+1:templen/2+extra+4*sps);
[coef, lag] = xcorr(d1, t1, offsetmax, 'coeff');
[maxcoef1, idx] = max(coef);
lagsamp1 = lag(idx);

[coef, lag] = xcorr(fwindata, ftemplate, offsetmax, 'coeff');
[maxcoef2, idx] = max(coef);
lagsamp2 = lag(idx);
difflag(istack) = lagsamp1-lagsamp2;
if abs(lagsamp1-lagsamp2)<=1
    npass = npass+1;

end

end
npass
figure
scatter(1: nstack, difflag); hold on
plot([0, 2500], [-1, -1], 'r-');
plot([0, 2500], [1, 1], 'r-');
