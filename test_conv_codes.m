clear all;
clc;

SNR = 5;

msg = randi([0 1], 1, 1000000);

% trel = poly2trellis(2,[2 3]);   % Rate = 1/2
%trel = poly2trellis([5 4],[23 35 0; 0 5 13]);  % Rate = 2/3
trel = poly2trellis(7, [133 171]);   % Rate = 3/4

%%
coded = convenc(msg, trel);

coded(coded == 0) = -1;

rx_coded = awgn(coded, SNR, 'measured');

rx_coded(rx_coded <= 0) = 0;
rx_coded(rx_coded > 0) = 1;

tblen = 50;
dec_rx_coded = vitdec(rx_coded, trel, tblen,'term','hard');

[length(msg)/length(coded) biterr(msg, dec_rx_coded)/length(dec_rx_coded)]

%%
uncoded = msg;

uncoded(uncoded == 0) = -1;

rx_uncoded = awgn(uncoded, SNR, 'measured');

rx_uncoded(rx_uncoded <= 0) = 0;
rx_uncoded(rx_uncoded > 0) = 1;

[length(msg)/length(uncoded) biterr(msg, rx_uncoded)/length(rx_uncoded)]