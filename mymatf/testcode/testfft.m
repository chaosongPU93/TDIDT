% this is for testing fft

close all
clear
%% test spectral analysis, self-made 
aaa=[0 0 0 0 0 1 10 500];
bbb=[aaa fliplr(aaa(1:end-1))];
figure
plot(bbb);
ns = 1;
n = 1:ns:length(bbb);
N = length(n);
k = -(N-1)/2:(N-1)/2;
f = k/(N*ns);

y = fftshift(fft(bbb));
subplot(211),stem(n,bbb),title('???');
subplot(212),stem(f,abs(y)/N),title('?????');

%%%%%%%%%%%%%%%%%%%%
%% an example from matlab library
fs = 100;                                % sample frequency (Hz)
t = 0:1/fs:10-1/fs;                      % 10 second span time vector
x = (1.3)*sin(2*pi*15*t) ...             % 15 Hz component
  + (1.7)*sin(2*pi*40*(t-2)) ...         % 40 Hz component
  ;

y = fft(x);
n = length(x);          % number of samples
f = (0:n-1)/n*fs;     % frequency range
power = abs(y).^2/n;    % power of the DFT

figure
plot(f,power)
xlabel('Frequency')
ylabel('Power')

y0 = fftshift(y);         % shift y values
f0 = (-n/2:n/2-1)*(fs/n); % 0-centered frequency range
power0 = abs(y0).^2/n;    % 0-centered power

figure
plot(f0,power0)
xlabel('Frequency')
ylabel('Power')

%% compare with peridogram
[xf,f,amp,pha,power,psd] = fftspectrum(x, n, fs,'onesided');
figure
plot(f,amp)
xlabel('Frequency')
ylabel('Amplitude')

[pxx,fp] = periodogram(x,rectwin(n),n,fs,'power','onesided');
figure
plot(fp,pxx); hold on
plot(f,power);
xlabel('Frequency')
ylabel('Power')

[pxx,fp] = periodogram(x,rectwin(n),n,fs,'onesided');
figure
plot(fp,pxx); hold on
plot(f,psd);
xlabel('Frequency')
ylabel('PSD')

%% compare with peridogram
[xf,f,amp,pha,power,psd] = fftspectrum(x, n, fs,'twosided');
figure
plot(f,amp)
xlabel('Frequency')
ylabel('Amplitude')

[pxx,fp] = periodogram(x,rectwin(n),n,fs,'power','twosided');
figure
plot(fp,pxx); hold on
plot(f,power);
xlabel('Frequency')
ylabel('Power')

[pxx,fp] = periodogram(x,rectwin(n),n,fs,'twosided');
figure
plot(fp,pxx); hold on
plot(f,psd);
xlabel('Frequency')
ylabel('PSD')

%%%%%%%%%%%%%%%%%%%
%% an example from matlab library
whaleFile = 'bluewhale.au';
[x,fs] = audioread(whaleFile);

figure
plot(x)
xlabel('Sample Number')
ylabel('Amplitude')
moan = x(2.45e4:3.10e4);
t = 10*(0:1/fs:(length(moan)-1)/fs);

figure
plot(t,moan)
xlabel('Time (seconds)')
ylabel('Amplitude')
xlim([0 t(end)])
m = length(moan);       % original sample length
n = pow2(nextpow2(m));  % transform length
y = fft(moan,n);        % DFT of signal
f = (0:n-1)/n*fs/10;
power = abs(y).^2/n;      

figure
plot(f(1:floor(n/2)),power(1:floor(n/2)))
xlabel('Frequency')
ylabel('Power')