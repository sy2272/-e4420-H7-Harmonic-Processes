clear all; clc; close all;
[y,Fs]=wavread('speech.wav');

K=10; % harmonic order K
P=15; % LPC order P

y=y(15000:15255);
T=length(y);
t=linspace(0,16,T);
%% harmonic model
f=[100:200];
for i=1:length(f)
    phase=[0:T-1]'*[1:K]*2*pi*f(i)/Fs;
    S=[sin(phase) cos(phase)]';
    b=S'\y;
    s(i)=std(y-S'*b);
end

[m,i]=min(s); w=f(i);
figure(1); subplot(2,1,1),plot(f,s),
title(['Error of fit = ' num2str(m)]), xlabel(['pitch frequency in Hz = ' num2str(w)]), ylabel('error');
% estimate for pitch frequency is 165 Hz
phase=[0:T-1]'*[1:K]*2*pi*w/Fs;
S=[sin(phase) cos(phase)]';
b=S'\y;
h=S'*b;
subplot(2,1,2),plot(t,y,t,h),title(['harmonic fit with K = ' num2str(K)]);
xlabel('time in ms'),ylabel('speech amplitude'),legend('original','harmonic fit');
saveas(gcf,'harmonic_fit.jpg')
% K=4 seems appropriate for signal

%% check orthogonal vectors
n=y-h; c=dot(n,h); 
% c=10^-16; n and h are orthogonal vectors 

%% noise model using lpc
% Akaike Information Criterion (AIC)
N=length(n);
for P=1:20
    a = lpc(n,P);
    est_n = filter(-1*[0 a(2:end)],1,n);
    err = est_n - n;
    Perror = (1/N) * (abs(fft(err,N))).^2;
    errsig = mean(Perror);
    AIC(P)=N*log(errsig)+2*P;
end
[mAIC,P]=min(AIC);
figure(2),plot(AIC)
title(['best noise model order P = ' num2str(P) ', AIC = ' num2str(mAIC)]);
xlabel('model order P'),ylabel('Akaike Value');
saveas(gcf,'findPAkaike.jpg')
% minimum AIC is AIC = -2784,  P=15

a=lpc(n,P);
est_n = filter(-1*[0 a(2:end)],1,n);
%% Power spectrum of combined model
ydft = fft(y);
psdy = (1/(Fs*length(y))).*abs(ydft(1:length(ydft)/2+1)).^2;
psdy(2:end-1) = 2*psdy(2:end-1);
Py=10*log10(psdy);

errordft = fft(h);
psdh = (1/(Fs*length(h))).*abs(errordft(1:length(errordft)/2+1)).^2;
psdh(2:end-1) = 2*psdh(2:end-1);    
Ph=10*log10(psdh);

ndft = fft(est_n);
psdn = (1/(Fs*length(est_n))).*abs(ndft(1:length(ndft)/2+1)).^2;
psdn(2:end-1) = 2*psdn(2:end-1);    

Pnh=10*log10(psdn+psdh);

f = 0:Fs/length(y):Fs/2;
figure(3),plot(f,Py,f,Ph,f,Pnh);
title(['combined model with K = ' num2str(K) ', P = ' num2str(P)]);
xlabel('Frequency in Hz');
ylabel('Power spectrum in dB');
legend('|X(e^iw)|^2','|H(e^(iw)|^2','|N(e^iw)^2+|H(e^iw)|^2');
saveas(gcf,'powerspectrum_combined.jpg');

%% Power spectrum of model error
err=est_n+h-y;
figure(4), subplot(2,1,1), plot(err);
title(['estimation error with K = ' num2str(K) ', P = ' num2str(P)]);
ylabel('amplitude'),xlabel('time in ms');
errordft = fft(err);
psderr = (1/(Fs*length(err))).*abs(errordft(1:length(errordft)/2+1)).^2;
psderr(2:end-1) = 2*psderr(2:end-1);    
Perr=10*log10(psderr);
subplot(2,1,2), plot(f,Perr)
ylabel('Power spectrum in dB'),xlabel('frequency in Hz');
saveas(gcf,'powerspectrum_error.jpg');

%% SNR of combined model
modelSNR=sum(est_n+h.^2)/sum(err.^2);
% SNR for K=6, P=15 is 7

%% SNR in dB as function of model order P
SNR=[]; dBSNR=[];

for i=1:20
    a=lpc(y,i);
    est_y = filter(-1*[0 a(2:end)],1,y);
    err=y-est_y;
    SNR(i)=sum(est_y.^2)/sum(err.^2);
    dBSNR(i)=db(SNR(i));  
end
PSNR=dBSNR(P);
figure(5),plot(1:20,dBSNR)
title(['SNR in dB for P = ' num2str(P) ', dBSNR = ' num2str(PSNR)]);
ylabel('SNR in dB'); xlabel('model order P');
saveas(gcf,'lpc_dbSNR.jpg');
