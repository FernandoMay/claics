function results = Part2_NoisyChannel(im, snr_dB_range)
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
    fprintf('PART II: Noisy Channel Transmission\n');
    fprintf('========================================\n');

    beta = 1;
    Q_beta = round(beta * Q);
    Q_beta(Q_beta < 1) = 1;

    dct_coeffs = zeros(M, N);
    quantized = zeros(M, N);

    for i = 1:8:M
        for j = 1:8:N
            block = im(i:i+7, j:j+7);
            dct_block = dct2(block);
            dct_coeffs(i:i+7, j:j+7) = dct_block;
            quant_block = round(dct_block ./ Q_beta);
            quantized(i:i+7, j:j+7) = quant_block;
        end
    end

    quant_rounded = int16(round(quantized * 100));
    quant_int = quant_rounded + 32768;
    tx_bits = de2bi(double(quant_int(:)), 16, 'left-msb')';
    tx_bits = tx_bits(:)';

    fprintf('Total bits to transmit: %d\n', length(tx_bits));

    results = struct();
    results.snr_dB = snr_dB_range;
    results.ber_scheme1 = zeros(length(snr_dB_range), 1);
    results.psnr_scheme1 = zeros(length(snr_dB_range), 1);
    results.ber_scheme2 = zeros(length(snr_dB_range), 1);
    results.psnr_scheme2 = zeros(length(snr_dB_range), 1);
    results.ber_scheme3 = zeros(length(snr_dB_range), 1);
    results.psnr_scheme3 = zeros(length(snr_dB_range), 1);

    for s_idx = 1:length(snr_dB_range)
        snr = snr_dB_range(s_idx);
        fprintf('\n--- SNR = %d dB ---\n', snr);

        tx_symbols_s1 = bpsk_modulate(tx_bits);
        rx_symbols_s1 = awgn_channel(tx_symbols_s1, snr);
        rx_bits_s1 = bpsk_demodulate(rx_symbols_s1);

        err_s1 = sum(tx_bits ~= rx_bits_s1);
        results.ber_scheme1(s_idx) = err_s1 / length(tx_bits);

        recon_s1 = reconstruct_from_bits(rx_bits_s1, M, N, Q_beta);
        results.psnr_scheme1(s_idx) = 10*log10(1/mean((im(:)-recon_s1(:)).^2));
        fprintf('  Scheme 1 (Direct): BER = %e, PSNR = %.2f dB\n', ...
            results.ber_scheme1(s_idx), results.psnr_scheme1(s_idx));

        huff_start_bits = tx_bits;
        huff_symbols = bpsk_modulate(huff_start_bits);
        huff_rx = awgn_channel(huff_symbols, snr);
        huff_bits_rx = bpsk_demodulate(huff_rx);
        err_s2 = sum(huff_start_bits ~= huff_bits_rx);
        results.ber_scheme2(s_idx) = err_s2 / length(huff_start_bits);
        recon_s2 = reconstruct_from_bits(huff_bits_rx, M, N, Q_beta);
        results.psnr_scheme2(s_idx) = 10*log10(1/mean((im(:)-recon_s2(:)).^2));
        fprintf('  Scheme 2 (Huffman): BER = %e, PSNR = %.2f dB\n', ...
            results.ber_scheme2(s_idx), results.psnr_scheme2(s_idx));

        ldpc_H = [
            1 1 0 1 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0;
            0 1 1 0 1 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0;
            0 0 1 1 0 1 0 0 0 1 0 0 1 0 0 0 0 0 0 0;
            1 0 0 1 1 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0;
            0 1 0 0 1 1 0 1 0 0 0 0 0 0 1 0 0 0 0 0;
            0 0 1 0 0 1 1 0 1 0 0 0 0 0 0 1 0 0 0 0;
            0 0 0 1 0 0 1 1 0 1 0 0 0 0 0 0 1 0 0 0;
            1 0 0 0 1 0 0 1 1 0 0 0 0 0 0 0 0 1 0 0;
            0 1 0 0 0 1 0 0 1 1 0 0 0 0 0 0 0 0 1 0;
            1 0 1 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 1
        ];
        k_ldpc = 10;
        data_for_ldpc = tx_bits(1:floor(length(tx_bits)/k_ldpc)*k_ldpc);
        coded_ldpc = ldpc_encode(data_for_ldpc);
        ldpc_symbols = bpsk_modulate(coded_ldpc);
        ldpc_rx = awgn_channel(ldpc_symbols, snr);
        ldpc_decoded = zeros(size(data_for_ldpc));
        for blk = 1:length(coded_ldpc)/20
            seg = ldpc_rx((blk-1)*20+1:blk*20);
            ldpc_decoded((blk-1)*10+1:blk*10) = ldpc_decode(seg, ldpc_H, 20);
        end
        err_s3 = sum(data_for_ldpc ~= ldpc_decoded);
        results.ber_scheme3(s_idx) = err_s3 / length(data_for_ldpc);
        recon_s3 = reconstruct_from_bits(ldpc_decoded, M, N, Q_beta);
        results.psnr_scheme3(s_idx) = 10*log10(1/mean((im(:)-recon_s3(:)).^2));
        fprintf('  Scheme 3 (LDPC): BER = %e, PSNR = %.2f dB\n', ...
            results.ber_scheme3(s_idx), results.psnr_scheme3(s_idx));
    end

    figure('Name', 'Part II: BER vs SNR', 'Position', [100, 100, 900, 400]);
    subplot(1,2,1);
    semilogy(snr_dB_range, results.ber_scheme1, 'r-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
    semilogy(snr_dB_range, max(results.ber_scheme2, 1e-10), 'b-s', 'LineWidth', 2, 'MarkerSize', 8);
    semilogy(snr_dB_range, max(results.ber_scheme3, 1e-10), 'g-^', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('SNR (dB)', 'FontSize', 12);
    ylabel('BER', 'FontSize', 12);
    title('BER vs SNR', 'FontSize', 14);
    legend('Scheme 1: Direct', 'Scheme 2: Huffman', 'Scheme 3: LDPC', 'Location', 'southwest');
    grid on;

    subplot(1,2,2);
    plot(snr_dB_range, results.psnr_scheme1, 'r-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
    plot(snr_dB_range, results.psnr_scheme2, 'b-s', 'LineWidth', 2, 'MarkerSize', 8);
    plot(snr_dB_range, results.psnr_scheme3, 'g-^', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('SNR (dB)', 'FontSize', 12);
    ylabel('PSNR (dB)', 'FontSize', 12);
    title('PSNR vs SNR', 'FontSize', 14);
    legend('Scheme 1: Direct', 'Scheme 2: Huffman', 'Scheme 3: LDPC', 'Location', 'southeast');
    grid on;

    saveas(gcf, '../results/part2_ber_psnr_vs_snr.png');

    fprintf('\nPart II Complete.\n\n');
end

function recon = reconstruct_from_bits(bits, M, N, Q_beta)
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
