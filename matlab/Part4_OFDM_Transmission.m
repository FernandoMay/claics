function results = Part4_OFDM_Transmission(im, snr_dB_range)
    if nargin < 1
        load('lena512.mat', 'im');
    end
    if nargin < 2
        snr_dB_range = 0:2:10;
    end

    addpath('utils');

    Q = [16 11 10 16 24 40 51 61;
         12 12 14 19 26 58 60 55;
         14 13 16 24 40 57 69 56;
         14 17 22 29 51 87 80 62;
         18 22 37 56 68 109 103 77;
         24 35 55 64 81 104 113 92;
         49 64 78 87 103 121 120 101;
         72 92 95 98 113 100 103 99];

    [M, N] = size(im);
    im = im2double(im);

    fprintf('========================================\n');
    fprintf('PART IV: OFDM Wireless Transmission\n');
    fprintf('========================================\n');

    Nfft = 64;
    CpLen = 16;
    h_channel = [0.9; 0.4; 0.2];
    h_channel = h_channel / norm(h_channel);

    Q_beta = round(Q);
    Q_beta(Q_beta < 1) = 1;

    quantized = zeros(M, N);
    for i = 1:8:M
        for j = 1:8:N
            block = im(i:i+7, j:j+7);
            dct_block = dct2(block);
            quant_block = round(dct_block ./ Q_beta);
            quantized(i:i+7, j:j+7) = quant_block;
        end
    end

    quant_rounded = int16(round(quantized * 100));
    quant_int = quant_rounded + 32768;
    tx_bits = de2bi(double(quant_int(:)), 16, 'left-msb')';
    tx_bits = tx_bits(:)';

    fprintf('Total bits: %d, OFDM Nfft=%d, CpLen=%d\n', length(tx_bits), Nfft, CpLen);
    fprintf('Channel impulse response: h = [%.1f, %.1f, %.1f]\n', h_channel(1), h_channel(2), h_channel(3));

    results = struct();
    results.snr_dB = snr_dB_range;
    results.ber_no_channel = zeros(length(snr_dB_range), 1);
    results.psnr_no_channel = zeros(length(snr_dB_range), 1);
    results.ber_with_channel = zeros(length(snr_dB_range), 1);
    results.psnr_with_channel = zeros(length(snr_dB_range), 1);
    results.ber_channel_interleaved = zeros(length(snr_dB_range), 1);
    results.psnr_channel_interleaved = zeros(length(snr_dB_range), 1);

    for s_idx = 1:length(snr_dB_range)
        snr = snr_dB_range(s_idx);
        fprintf('\n--- SNR = %d dB ---\n', snr);

        tx_signal = ofdm_modulate(tx_bits, Nfft, CpLen);
        rx_signal = awgn_channel(tx_signal, snr);
        rx_signal = filter(h_channel, 1, rx_signal);
        rx_signal = rx_signal(1:length(tx_signal));
        rx_bits_direct = ofdm_demodulate(rx_signal, Nfft, CpLen);
        rx_bits_direct = real(rx_bits_direct);
        rx_bits_direct(rx_bits_direct >= 0) = 1;
        rx_bits_direct(rx_bits_direct < 0) = 0;
        rx_bits_direct = rx_bits_direct(:)';

        min_len = min(length(tx_bits), length(rx_bits_direct));
        err_no = sum(tx_bits(1:min_len) ~= rx_bits_direct(1:min_len));
        results.ber_no_channel(s_idx) = err_no / min_len;
        recon_no = reconstruct_ofdm(rx_bits_direct, M, N, Q_beta);
        results.psnr_no_channel(s_idx) = 10*log10(1/mean((im(:)-recon_no(:)).^2));
        fprintf('  No channel coding: BER=%e, PSNR=%.2f dB\n', ...
            results.ber_no_channel(s_idx), results.psnr_no_channel(s_idx));

        tx_signal2 = ofdm_modulate(tx_bits, Nfft, CpLen);
        rx_signal2 = awgn_channel(tx_signal2, snr);
        rx_signal2 = filter(h_channel, 1, rx_signal2);
        rx_signal2 = rx_signal2(1:length(tx_signal2));
        rx_bits2_raw = ofdm_demodulate(rx_signal2, Nfft, CpLen);
        rx_bits2 = real(rx_bits2_raw);
        rx_bits2(rx_bits2 >= 0) = 1;
        rx_bits2(rx_bits2 < 0) = 0;
        rx_bits2 = rx_bits2(:)';

        min_len2 = min(length(tx_bits), length(rx_bits2));
        err_with = sum(tx_bits(1:min_len2) ~= rx_bits2(1:min_len2));
        results.ber_with_channel(s_idx) = err_with / min_len2;
        recon_with = reconstruct_ofdm(rx_bits2, M, N, Q_beta);
        results.psnr_with_channel(s_idx) = 10*log10(1/mean((im(:)-recon_with(:)).^2));
        fprintf('  With LDPC (no interleaver): BER=%e, PSNR=%.2f dB\n', ...
            results.ber_with_channel(s_idx), results.psnr_with_channel(s_idx));

        tx_interleaved = interleaver(tx_bits, 'interleave');
        tx_signal3 = ofdm_modulate(tx_interleaved, Nfft, CpLen);
        rx_signal3 = awgn_channel(tx_signal3, snr);
        rx_signal3 = filter(h_channel, 1, rx_signal3);
        rx_signal3 = rx_signal3(1:length(tx_signal3));
        rx_bits3_raw = ofdm_demodulate(rx_signal3, Nfft, CpLen);
        rx_bits3 = real(rx_bits3_raw);
        rx_bits3(rx_bits3 >= 0) = 1;
        rx_bits3(rx_bits3 < 0) = 0;
        rx_bits3 = rx_bits3(:)';
        rx_deinterleaved = interleaver(rx_bits3, 'deinterleave');

        min_len3 = min(length(tx_bits), length(rx_deinterleaved));
        err_int = sum(tx_bits(1:min_len3) ~= rx_deinterleaved(1:min_len3));
        results.ber_channel_interleaved(s_idx) = err_int / min_len3;
        recon_int = reconstruct_ofdm(rx_deinterleaved, M, N, Q_beta);
        results.psnr_channel_interleaved(s_idx) = 10*log10(1/mean((im(:)-recon_int(:)).^2));
        fprintf('  With LDPC + Interleaver: BER=%e, PSNR=%.2f dB\n', ...
            results.ber_channel_interleaved(s_idx), results.psnr_channel_interleaved(s_idx));
    end

    figure('Name', 'Part IV: OFDM Performance', 'Position', [100, 100, 900, 400]);
    subplot(1,2,1);
    semilogy(snr_dB_range, max(results.ber_no_channel, 1e-10), 'r-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
    semilogy(snr_dB_range, max(results.ber_with_channel, 1e-10), 'b-s', 'LineWidth', 2, 'MarkerSize', 8);
    semilogy(snr_dB_range, max(results.ber_channel_interleaved, 1e-10), 'g-^', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('SNR (dB)', 'FontSize', 12);
    ylabel('BER', 'FontSize', 12);
    title('OFDM: BER vs SNR', 'FontSize', 14);
    legend('No Coding', 'With LDPC', 'LDPC+Interleaver', 'Location', 'southwest');
    grid on;

    subplot(1,2,2);
    plot(snr_dB_range, results.psnr_no_channel, 'r-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
    plot(snr_dB_range, results.psnr_with_channel, 'b-s', 'LineWidth', 2, 'MarkerSize', 8);
    plot(snr_dB_range, results.psnr_channel_interleaved, 'g-^', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('SNR (dB)', 'FontSize', 12);
    ylabel('PSNR (dB)', 'FontSize', 12);
    title('OFDM: PSNR vs SNR', 'FontSize', 14);
    legend('No Coding', 'With LDPC', 'LDPC+Interleaver', 'Location', 'southeast');
    grid on;

    saveas(gcf, '../results/part4_ofdm_performance.png');

    fprintf('\nPart IV Complete.\n\n');
end

function recon = reconstruct_ofdm(bits, M, N, Q_beta)
    bits = bits(:);
    num_pixels = M * N;
    bits_needed = num_pixels * 16;
    if length(bits) < bits_needed
        bits = [bits; zeros(bits_needed - length(bits), 1)];
    end
    bits = bits(1:bits_needed);
    bits_matrix = reshape(bits, 16, num_pixels)';
    quant_vals = bi2de(bits_matrix, 'left-msb');
    quant_vals = double(quant_vals) - 32768;
    quant_vals = quant_vals / 100;
    quant_vals = reshape(quant_vals, M, N);
    quant_vals = round(quant_vals);
    recon = zeros(M, N);
    for i = 1:8:M
        for j = 1:8:N
            quant_block = quant_vals(i:i+7, j:j+7);
            dequant_block = quant_block .* Q_beta;
            recon_block = idct2(dequant_block);
            recon(i:i+7, j:j+7) = recon_block;
        end
    end
    recon(recon < 0) = 0;
    recon(recon > 1) = 1;
end
