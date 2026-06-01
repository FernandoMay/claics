function results = Part3_ScalableCoding(im, snr_weak_dB, snr_diff_dB)
    if nargin < 1
        load('lena512.mat', 'im');
    end
    if nargin < 2
        snr_weak_dB = 5;
    end
    if nargin < 3
        snr_diff_dB = 10;
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
    im = double(im);

    fprintf('========================================\n');
    fprintf('PART III: Scalable Image Coding (FGS)\n');
    fprintf('========================================\n');

    Q_beta = round(1 * Q);
    Q_beta(Q_beta < 1) = 1;

    quantized_all = zeros(M, N);
    for i = 1:8:M
        for j = 1:8:N
            block = im(i:i+7, j:j+7);
            dct_block = dct2(block);
            quant_block = round(dct_block ./ Q_beta);
            quantized_all(i:i+7, j:j+7) = quant_block;
        end
    end

    min_val = min(quantized_all(:));
    max_val = max(quantized_all(:));
    quant_norm = quantized_all - min_val;
    n_bits_total = ceil(log2(max_val - min_val + 1));
    fprintf('Quantized range: [%d, %d], bits needed: %d\n', min_val, max_val, n_bits_total);

    fprintf('Using %d bits for base layer, remaining for enhancement\n', ...
        min(4, max(1, ceil(n_bits_total/2))));

    n_base_bits = min(4, max(1, ceil(n_bits_total/2)));
    n_enhance_bits = n_bits_total - n_base_bits;

    base_layer = floor(quant_norm / (2^n_enhance_bits));
    enhance_layer = mod(quant_norm, 2^n_enhance_bits);

    base_layer_bits = de2bi(double(base_layer(:)), n_base_bits, 'left-msb')';
    base_layer_bits = base_layer_bits(:)';

    snr_strong_dB = snr_weak_dB + snr_diff_dB;

    fprintf('Weak user SNR: %d dB\n', snr_weak_dB);
    fprintf('Strong user SNR: %d dB\n', snr_strong_dB);

    C_weak = log2(1 + 10^(snr_weak_dB/10));
    C_strong = log2(1 + 10^(snr_strong_dB/10));
    fprintf('Weak user channel capacity: %.2f bits/symbol\n', C_weak);
    fprintf('Strong user channel capacity: %.2f bits/symbol\n', C_strong);

    tx_base_symbols = bpsk_modulate(base_layer_bits);
    rx_base_symbols = awgn_channel(tx_base_symbols, snr_weak_dB);
    rx_base_bits = bpsk_demodulate(rx_base_symbols);

    rx_base_bits = rx_base_bits(:)';
    base_layer_rx = reshape(rx_base_bits(1:M*N*n_base_bits), n_base_bits, M*N)';
    base_layer_rx = bi2de(base_layer_rx, 'left-msb');

    recon_weak = base_layer_rx * (2^n_enhance_bits);
    recon_weak = reshape(recon_weak, M, N) + min_val;
    recon_weak(recon_weak < min_val) = min_val;
    recon_weak(recon_weak > max_val) = max_val;

    recon_weak_img = zeros(M, N);
    for i = 1:8:M
        for j = 1:8:N
            quant_block = recon_weak(i:i+7, j:j+7);
            dequant_block = quant_block .* Q_beta;
            recon_block = idct2(dequant_block);
            recon_weak_img(i:i+7, j:j+7) = recon_block;
        end
    end
    recon_weak_img(recon_weak_img < 0) = 0;
    recon_weak_img(recon_weak_img > 255) = 255;

    enhance_layer_bits = de2bi(double(enhance_layer(:)), n_enhance_bits, 'left-msb')';
    enhance_layer_bits = enhance_layer_bits(:)';
    tx_enhance_symbols = bpsk_modulate(enhance_layer_bits);
    rx_enhance_symbols = awgn_channel(tx_enhance_symbols, snr_strong_dB);
    rx_enhance_bits = bpsk_demodulate(rx_enhance_symbols);

    rx_enhance_bits = rx_enhance_bits(:)';
    enhance_len = M*N*n_enhance_bits;
    if length(rx_enhance_bits) >= enhance_len
        enhance_layer_rx = reshape(rx_enhance_bits(1:enhance_len), n_enhance_bits, M*N)';
        enhance_layer_rx = bi2de(enhance_layer_rx, 'left-msb');
    else
        enhance_layer_rx = zeros(M*N, 1);
    end

    recon_strong = (base_layer_rx * (2^n_enhance_bits)) + enhance_layer_rx;
    recon_strong = reshape(recon_strong, M, N) + min_val;
    recon_strong(recon_strong < min_val) = min_val;
    recon_strong(recon_strong > max_val) = max_val;

    recon_strong_img = zeros(M, N);
    for i = 1:8:M
        for j = 1:8:N
            quant_block = recon_strong(i:i+7, j:j+7);
            dequant_block = quant_block .* Q_beta;
            recon_block = idct2(dequant_block);
            recon_strong_img(i:i+7, j:j+7) = recon_block;
        end
    end
    recon_strong_img(recon_strong_img < 0) = 0;
    recon_strong_img(recon_strong_img > 255) = 255;

    psnr_weak = 10*log10(255^2/mean((im(:)-recon_weak_img(:)).^2));
    psnr_strong = 10*log10(255^2/mean((im(:)-recon_strong_img(:)).^2));

    fprintf('\nWeak user PSNR: %.2f dB\n', psnr_weak);
    fprintf('Strong user PSNR: %.2f dB\n', psnr_strong);

    results = struct();
    results.snr_weak_dB = snr_weak_dB;
    results.snr_strong_dB = snr_strong_dB;
    results.C_weak = C_weak;
    results.C_strong = C_strong;
    results.psnr_weak = psnr_weak;
    results.psnr_strong = psnr_strong;
    results.recon_weak = recon_weak_img;
    results.recon_strong = recon_strong_img;

    figure('Name', 'Part III: Scalable Coding', 'Position', [100, 100, 1400, 400]);
    subplot(1,4,1);
    imshow(uint8(im));
    title('Original', 'FontSize', 12);

    base_recon_img = zeros(M, N);
    for i = 1:8:M
        for j = 1:8:N
            qblock = base_layer(i:i+7, j:j+7) * (2^n_enhance_bits);
            qblock = qblock .* Q_beta;
            recon_block = idct2(qblock);
            base_recon_img(i:i+7, j:j+7) = recon_block;
        end
    end
    base_recon_img(base_recon_img < 0) = 0;
    base_recon_img(base_recon_img > 255) = 255;
    psnr_base_only = 10*log10(255^2/mean((im(:)-base_recon_img(:)).^2));
    subplot(1,4,2);
    imshow(uint8(base_recon_img));
    title(sprintf('Base Layer Only\nPSNR=%.2f dB', psnr_base_only), 'FontSize', 11);

    subplot(1,4,3);
    imshow(uint8(recon_weak_img));
    title(sprintf('Weak User\nPSNR=%.2f dB', psnr_weak), 'FontSize', 11);

    subplot(1,4,4);
    imshow(uint8(recon_strong_img));
    title(sprintf('Strong User\nPSNR=%.2f dB', psnr_strong), 'FontSize', 11);

    print(gcf, '../results/part3_scalable_coding.png', '-dpng');

    fprintf('\nPart III Complete.\n\n');
end
