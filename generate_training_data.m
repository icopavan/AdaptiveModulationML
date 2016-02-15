clear all;
%clc;

N = 256;  % FFT length or the number of active subcarriers
iter = 100; % number of iterations (i.e., number of OFDM symbols simulated)
Ncp = 67;
EbNodB_all = randi([0 35],1,1);

M_all = [2 4 16 64 256];

ch_pow = 25*rand(1,1);
H = ch_pow*(randn(1,N) + 1i*randn(1, N));

ber_th = 0.05;

for M_itr = 1:length(M_all)
    if M_itr == 1
        M = M_all(M_itr);
        %% BPSK
        for itr = 1:length(EbNodB_all)
            EbNodB = EbNodB_all(itr); % Eb/No in dB
            EbNo = 10.^(EbNodB/10); % convert to linear
            sigma = sqrt(1/2/EbNo); % noise standard deviation

            SNR = sort(diag(10*log10(real(1/sigma^2*ctranspose(H)*H))));
            
            errors_bin = 0;
            errors_int = 0;

            for i = 1:iter
                data_int = sign(randn(1,N)); % create data
                data_bin = data_int;
                x = ifft(data_bin);% create tx signal

                h = ifft(H);
                HX = fft(x).*H;
                hx = ifft(HX);

                n = sigma*randn(size(hx))+1i*sigma*randn(size(hx)); % create noise
                r = hx + n; % create received signal

                r_fft = fft(r); % convert to frequency domain
                Z = r_fft./H;  % create decision variables by doing frequency equalization

                bhat = sign(real(Z)); % create data estimates
                bhat_int = bhat;
                errors_bin = errors_bin + sum(abs(data_bin-bhat))/2; % count errors
                errors_int = errors_int + length(find(data_int ~= bhat_int));
            end
            ber(itr) = errors_bin/length(data_bin)/iter; % calculate BER
            ser(itr) = errors_int/length(data_int)/iter; % calculate BER
        end
    else
    %% M-QAM
        M = M_all(M_itr);

        hMod = comm.RectangularQAMModulator(M);
        const = step(hMod, (0:M-1)');
        scale = modnorm(const, 'avpow',1);
        for itr = 1:length(EbNodB_all)
            EbNodB = EbNodB_all(itr); % Eb/No in dB
            EbNo = 10.^(EbNodB/10); % convert to linear
            sigma = sqrt(1/2/EbNo); % noise standard deviation

            SNR = sort(diag(10*log10(real(1/sigma^2*ctranspose(H)*H))));
            
            errors_bin = 0;
            errors_int = 0;

            for i = 1:iter
                %% M-QAM
                data_int = randi([0 M-1], 1, N); % create data
                data_bin = reshape(de2bi(data_int, 'left-msb', log2(M))', 1, length(data_int)*log2(M));
                s = scale*qammod(data_int, M, 0, 'gray');

                x = ifft(s);% create tx signal

                h = ifft(H);
                HX=fft(x).*H;
                hx=ifft(HX);

                n = sigma*randn(size(hx)) + 1i*sigma*randn(size(hx)); % create noise
                r = hx + n; % create received signal

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
    end
    ber_all(M_itr) = ber;
    ser_all(M_itr) = ser;
end

MCS_ind = max(find(ber_all <= ber_th));
if isempty(MCS_ind)
    MCS_ind = 1;
end

ber_all;
ber_th;
MCS_ind;

dlmwrite('true_data_w_BER.csv',[SNR' MCS_ind ber_all], '-append');