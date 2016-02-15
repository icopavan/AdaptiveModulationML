clear all;
clc;
% clf;

N = 256; % FFT length or the number of active subcarriers
iter = 1; % number of iterations (i.e., number of OFDM symbols simulated)
Ncp = 67;
EbNodB_all = 10;

load 'channel_instantiation.mat';
Ts = 1/20e6; % sampling period of the channel
Fd = 0; % maximum Doppler frequency
h = rayleighchan(Ts, Fd, tau, pdb);
h.StoreHistory = 0;
h.StorePathGains = 1;
h.ResetBeforeFiltering = 1;
h.NormalizePathGains = 1;

M_all = [4 16 64 256];

for M_itr = 1:length(M_all)
    M = M_all(M_itr);
    ber = zeros(1, length(EbNodB_all));
    ser = zeros(1, length(EbNodB_all));
    thy = zeros(1, length(EbNodB_all));

    hMod = comm.RectangularQAMModulator(M);
    const = step(hMod, (0:M-1)');
    scale = modnorm(const, 'avpow',1);
    for itr = 1:length(EbNodB_all)
        EbNodB = EbNodB_all(itr); % Eb/No in dB
        EbNo = 10.^(EbNodB/10); % convert to linear
        sigma = sqrt(1/2/EbNo); % noise standard deviation

        errors_bin = 0;
        errors_int = 0;
        for i = 1:iter
            %% BPSK
%             data_bin = sign(randn(1,N)); % create data
%             x = sqrt(N)*ifft(data_bin);% create tx signal
%             x_wCP = [x((N-Ncp + 1):256) x]; %Add CP
%             hx = filter(h, x_wCP); % Pass through Rayleigh channel
%             g = h.PathGains; % Channel coefficients
%             G = fft(g, N)/sqrt(N); % DFT of channel coefficients
%             n = sigma*randn(size(hx))+1i*sigma*randn(size(hx)); % create noise
%             r = hx + n; % create received signal
%             r_woCP = r(Ncp + 1:end); %remove CP
%             r_fft = (fft(r_woCP)/sqrt(N)); % convert to frequency domain
%             Z = r_fft./G;  % create decision variables by doing frequency equalization
%             bhat = sign(real(Z)); % create data estimates
%             errors_bin = errors_bin + sum(abs(data_bin-bhat))/2; % count errors

            %% M-QAM
            data_int = randi([0 M-1], 1, N); % create data
            data_bin = reshape(de2bi(data_int, 'left-msb', log2(M))', 1, length(data_int)*log2(M));
            s = scale*qammod(data_int, M, 0, 'gray');
            
            x = ifft(s);% create tx signal
            x_wCP = [x((N-Ncp + 1):N) x]; %Add CP
            
            ch_pow = 10*rand(1,1);
            H = ch_pow*(randn(size(x)) + 1i*randn(size(x)));
            h = ifft(H);
            %hx=ifft(fft(x)*H);
            HX=fft(x).*H;
            hx=ifft(HX);
            %hx = conv(h, x_wCP);
            
            SNR = 10*log10(real(1/sigma^2*ctranspose(H)*H));
            [min(sort(diag(SNR))) max(sort(diag(SNR)))]
            %g = h.PathGains; % Channel coefficients
            %H = fft(g, N); % DFT of channel coefficients
            n = sigma*randn(size(hx)) + 1i*sigma*randn(size(hx)); % create noise
            r = hx + n; % create received signal
            %r_woCP = r(Ncp + 1:end); %remove CP
            r_fft = fft(r); % convert to frequency domain
            Z = r_fft./H;  % create decision variables by doing frequency equalization

            bhat_int = qamdemod(Z/scale, M, 0, 'gray'); % create data estimates
            bhat_bin = reshape(de2bi(bhat_int, 'left-msb', log2(M))', 1, length(bhat_int)*log2(M));
            errors_bin = errors_bin + biterr(data_bin, bhat_bin); % count errors

            errors_int = errors_int + length(find(data_int ~= bhat_int));
        end
        ber(itr) = errors_bin/length(data_bin)/iter; % calculate BER
        ser(itr) = errors_int/length(data_int)/iter; % calculate BER
    end

    [M ber(1) ser(1)]
end