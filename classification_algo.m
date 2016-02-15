clear all;
clc;

ber_th = 0.05;

tot_OFDM_sym = 1000;
bits_per_sym = [1 2 4 6 8];

data = dlmread('true_data_w_BER.csv');

tr_set = floor(size(data,1)*0.6);   %Percentage of data for training
tt_set = size(data,1) - tr_set;     %Percentage of data for testing
feat = 256;            %Number of features

tr_feat = data(1:tr_set, 1:feat);
tr_label = data(1:tr_set, feat+1);

tt_feat = data(tr_set+1:end, 1:feat);
tt_label = data(tr_set+1:end, feat+1);

tt_ber_all = data(tr_set+1:end, feat+2:end);

for k = 1:256
    pred_label = knnclassify(tt_feat, tr_feat, tr_label, k);
    accuracy(k) = length(find(tt_label == pred_label))/tt_set;
end

%%
[~, K] = max(accuracy);
pred_label = knnclassify(tt_feat, tr_feat, tr_label, K);

for i = 1:length(tt_ber_all)
    tt_ber(i) = tt_ber_all(i, tt_label(i));
    pred_ber(i) = tt_ber_all(i, pred_label(i));
    
    [tt_label(i) bits_per_sym(tt_label(i))]
    tt_throughput(i) = tot_OFDM_sym*bits_per_sym(tt_label(i))*(1 - tt_ber(i));
    pred_throughput(i) = tot_OFDM_sym*bits_per_sym(pred_label(i))*(1 - pred_ber(i));
end


figure(1); clf;
plot(accuracy, 'Linewidth', 2);
grid on;
xlabel('k');
ylabel('Classification accuracy on test set');

[xdata_val, xdata_ind] = sort(mean(tt_feat, 2));


figure(2); clf;
plot(xdata_val, pred_ber(xdata_ind), 'r-', 'Linewidth', 2);
hold on
grid on
plot(xdata_val, tt_ber(xdata_ind), 'b.', 'Linewidth', 2);
plot(xdata_val, ber_th*ones(size(xdata_val)), 'k--', 'Linewidth', 2.5);
xlabel('Mean SNR across subcarriers (dB)')
ylabel('BER')
legend('Proposed', 'Optimal', 'Target', 'Location', 'Best')


figure(3); clf;
tt_throughput = tt_throughput/max(tt_throughput);
pred_throughput = pred_throughput/max(pred_throughput);

plot(xdata_val, tt_throughput(xdata_ind), 'r-', 'Linewidth', 2);
hold on
grid on
plot(xdata_val, pred_throughput(xdata_ind), 'b.', 'Linewidth', 2);
xlabel('Mean SNR across subcarriers (dB)')
ylabel('Normalized throughput')
legend('Proposed', 'Optimal', 'Location', 'Best')


figure(4); clf;
stem(xdata_val, abs(pred_label-tt_label), 'ro', 'Linewidth', 2);
xlabel('Mean SNR across subcarriers (dB)')
ylabel('Difference in modulation order (optimal - proposed)')
axis([min(xdata_val) max(xdata_val) 0 length(bits_per_sym)])
grid on