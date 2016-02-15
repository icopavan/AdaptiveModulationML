clc;
close all;
clear all;

load 'channel_instantiation.mat'; 

N=256;
iter = 1000; % number of iterations (i.e., number of OFDM symbols simulated)

EbNodB = 0:20; % Eb/No in dB
EbNo = 10.^(EbNodB/10); % convert to linear
sigma = sqrt(1./(2*EbNo)); % noise standard deviation

errors = zeros(length(EbNo),1);

n_taps =60; % plot of abs(h) shows that after roughly 60 samples, abs(h) is insignificant. 
Cp=60; %CP should at least be equal to delay spread.


h(1:n_taps) = h(1:n_taps)/sqrt(sum(abs(h(1:n_taps)).^2)); %Normalizing
tau = t(1:n_taps)*1e-6;
pdb = 10*log10(abs(h(1:n_taps)).^2);

% chan = rayleighchan(1/20e6, 0, tau, pdb);
% chan.StoreHistory = 0;
% chan.StorePathGains = 1;
% chan.ResetBeforeFiltering = 1;
% chan.NormalizePathGains = 1;
% 
% for k=1:length(EbNo)
%     k
%     for i=1:iter
%         b = sign(randn(1,N)); % create data
%         x = sqrt(N)*ifft(b);% create tx signal
%         x_cp = [x(:,end-Cp+1:end) x]; %Add cyclic prefix
%         hx = filter(chan,x_cp);
%         n = sigma(k)*randn(1,length(hx))+1i*sigma(k)*randn(1,length(hx)); % create noise
%         r_cp = hx + n;
%         r = r_cp(:,Cp+1:end); %Remove Cyclic Prefix
%         R = fft(r)/sqrt(N);
%         
%         Hx = chan.PathGains;
%         Hx = fft(Hx, N)/sqrt(N); 
%         
%         x_hat = R./Hx;
%         bhat = sign(real(x_hat)); % create data estimates
%         errors(k) = errors(k) + sum(abs(b-bhat))/2; % count errors
%     end
% end

%% Using Convolution 
%%% (gives a different output each time I run qn2.1)

h = h(:,1:n_taps);
h = h/sqrt(sum(abs(h).^2));
H = fftshift(fft(h,N));

for k=1:length(EbNo)
    k
    for i=1:iter
        b = rand(1,N) > 0.5; % create data
        b_bpsk = 2*b - 1;
        x = sqrt(N)*ifft(fftshift(b_bpsk));% create tx signal
        x_cp = [x(:,end-Cp+1:end) x];
        hx = conv(h, x_cp);
        n = sigma(k)*randn(1,length(hx))+1i*sigma(k)*randn(1,length(hx)); % create noise
        r_cp = hx + n;
        r = r_cp(:,Cp+1:Cp+N); %Remove Cyclic Prefix
        R = fftshift(fft(r))/sqrt(N);
        x_hat = R./H;
        bhat =  2*floor(real(x_hat/2)) + 1; % create data estimates
        bhat(find(bhat>1)) = +1;
        bhat(find(bhat<-1)) = -1;
        bhat = (bhat+1)/2;
        errors(k) = errors(k) + biterr(b,bhat); % count errors
    end
end

ber = (errors/iter/N)'; % calculate BER
thy = 0.5*(1-sqrt(EbNo./(EbNo+1))); %Theoretical BER

figure(2);
a=semilogy(EbNodB,ber,'-r*','LineWidth',2);
hold on
grid on
a1=semilogy(EbNodB,thy,'-ks','LineWidth',2);
legend([a a1],'Simulated BER','Theoretical BER')

