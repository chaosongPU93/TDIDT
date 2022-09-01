% this code to do a simple test to determine the timeshift between the fast
% & the slow shear wave.
% 
% chaosong@princeton.edu
% Last modified? 2019/02/14

close all
data=load("/home/data2/chaosong/matlab/allan/split_correction/MGCB_splt");
pgc = data.MGCB_splt;
slowb=pgc(:, 1);    % slow component,before shift due to correction
fastb=pgc(:, 2);    % after shift
slowa=pgc(:, 3);    % fast component
fasta=pgc(:, 4);

sps = 40;
winsec = 4;
nfam = floor(length(slowb)/sps/winsec);
i = 1;
inds = (i-1)*sps*winsec+1;
inde = i*sps*winsec;
sb1 = slowb(inds: inde);
fb1 = fastb(inds: inde);
sa1 = slowa(inds: inde);
fa1 = fasta(inds: inde);

figure
plot(1:sps*winsec, sb1, 'b'); hold on
plot(1:sps*winsec, fb1, 'r');
% plot(sa1, 'k');
% plot(fa1, 'c');

[coeff,lag] = xcorr(sb1, fb1, 'coeff');
[mcoe, I] = max(abs(coeff));
shift = lag(I)

figure
plot(1:sps*winsec, sb1, 'b'); hold on
plot(1+shift:sps*winsec+shift, fb1, 'r');



