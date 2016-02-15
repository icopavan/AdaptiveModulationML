clear all;
clc;
clf;

N = 256; % FFT length or the number of active subcarriers
iter = 1000; % number of iterations (i.e., number of OFDM symbols simulated)
M = 16;

EbNodB_all = 0:12;
ber = zeros(1, length(EbNodB_all));
thy1 = zeros(1, length(EbNodB_all));

for itr = 1:length(EbNodB_all)
    EbNodB = EbNodB_all(itr); % Eb/No in dB
    EsNo = log2(M)*10.^(EbNodB/10); % convert to linear
    sigma = sqrt(10)*sqrt(1/(2*EsNo)); % noise standard deviation
    errors = 0;
    for i = 1:iter
        data_int = randi([0 M-1], 1, N); % create data
        data_bin = reshape(de2bi(data_int, 'left-msb', sqrt(M))', 1, length(data_int)*sqrt(M));
        s = qammod(data_int, M, 0, 'gray');
        x = sqrt(N)*ifft(s);% create tx signal
        n = sigma*randn(1,N)+1i*sigma*randn(1,N); % create noise
        r = x + n; % create received signal
        Z = fft(r)/sqrt(N); % create decision variables
        bhat_int = qamdemod(Z, M, 0, 'gray'); % create data estimates
        bhat_bin = reshape(de2bi(bhat_int, 'left-msb', sqrt(M))', 1, length(bhat_int)*sqrt(M));
        errors = errors + biterr(data_bin, bhat_bin); % count errors
    end
    ber(itr) = errors/(iter*length(data_bin)); % calculate BER
    thy1(itr) = (4/sqrt(M))*(1-1/sqrt(M))*(1/2)*erfc(sqrt(3*sqrt(M)*10^(EbNodB/10)/(M-1))/sqrt(2)); % theoretical BER
end

semilogy(EbNodB_all, ber, '-ro', 'Linewidth', 2, 'Markersize', 10);
hold on;
semilogy(EbNodB_all, thy1, '--b*', 'Linewidth', 2);
xlabel('Eb/No (dB)')
ylabel('BER')
legend('Simulated', 'Theoretical')
grid on