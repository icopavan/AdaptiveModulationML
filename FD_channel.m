clear all;
clc;
%close all;
% clf;

Pt = 1;

%decay rate per microsec
gamma = 0.03;

N = 256;

% sampling period (in Frequency domain) in MHz
D_f = 20/N;

K = 0;
sigma = sqrt(Pt/N);
f = 0:D_f:(N-1)*D_f;
t = [0:N-1]/(N*D_f); %(micro seconds)

% create random frequency channel values
W = sigma*randn(1,N);

% transform to the time domain
w = sqrt(N)*ifft(W);

g = (1/(N*D_f))*Pt*gamma/(K+1)*exp(-gamma*abs(t-(N/2)/(N*D_f)));

h_tmp = fftshift(w).*sqrt(g);
H_tmp = 1/sqrt(N)*fft(h_tmp);

H_tmp2 = hilbert(real(H_tmp));

H = H_tmp + j*imag(H_tmp2);

Ht = sqrt(1/(K+1))*H/sqrt(sum(abs(H).^2));

H = Ht + sqrt(K/(K+1))*1/sqrt(N);

tau_rms = (1/gamma)/(N*D_f)*sqrt(2*K+1)/(K+1)

% figure(1)
% plot(f*10^6, abs(H))
% xlabel('Frequency (Hz)')
% ylabel('|H(f)|')
% figure(2)
% plot(f*10^6, 10*log10(abs(H)))
% xlabel('Frequency (Hz)')
% ylabel('|H(f)|')
% figure(3)
% plot(f, unwrap(angle(H)))
% xlabel('Frequency (Hz)')
% ylabel('angle(H(f))')

h = ifft(H);
h = fftshift(h);
h = h/sqrt(sum(abs(h).^2));
%figure(4)
% clf
% plot(t, abs(h).^2)
% hold on
% plot(t-128/(N*D_f), 1/(K+1)*abs(g),'r')
% axis([0 1/2/D_f 0 0.25])
% xlabel('Delay (\mu sec)')
% ylabel('PDP')
% legend('Channel Instantiation','\phi_h(\tau)')

% compute and save values for problem 2.2
no_of_MPcomp = 60;
pdb = 10*log10(abs(h(1:no_of_MPcomp)).^2);
tau = t(1:no_of_MPcomp)*1e-6;
save('channel_instantiation.mat', 'pdb', 'tau');