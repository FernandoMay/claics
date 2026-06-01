function generate_report_figures()
    addpath('utils');
    rng(42);
    if ~exist('../results', 'dir')
        mkdir('../results');
    end

    if exist('lena512.mat', 'file')
        load('lena512.mat', 'im');
        disp('Loaded lena512.mat');
    elseif exist('../lena512.mat', 'file')
        load('../lena512.mat', 'im');
        disp('Loaded ../lena512.mat');
    elseif exist('lena512.png', 'file')
        im = double(imread('lena512.png'));
        disp('Loaded lena512.png');
    elseif exist('cameraman.tif', 'file')
        im_raw = imread('cameraman.tif');
        im = double(imresize(im_raw, [512, 512]));
        disp('Loaded cameraman.tif');
    else
        [x, y] = meshgrid(1:512, 1:512);
        im_raw = 128 + 50*sin(0.02*x+0.015*y) + 30*cos(0.04*x) + 20*sin(0.06*y);
        kernel = [1 2 1; 2 4 2; 1 2 1] / 16;
        im_raw = conv2(im_raw, kernel, 'same');
        im_raw(im_raw < 0) = 0;
        im_raw(im_raw > 255) = 255;
        im = double(im_raw);
        disp('Generated synthetic test image');
    end

    Q = [16 11 10 16 24 40 51 61;
         12 12 14 19 26 58 60 55;
         14 13 16 24 40 57 69 56;
         14 17 22 29 51 87 80 62;
         18 22 37 56 68 109 103 77;
         24 35 55 64 81 104 113 92;
         49 64 78 87 103 121 120 101;
         72 92 95 98 113 100 103 99];
    [M, N] = size(im);

    %% Part IV: OFDM reconstruction
    fprintf('Generating Part IV: OFDM reconstruction...\n');
    Nfft = 64; CpLen = 16;
    h_channel = [0.9; 0.4; 0.2] / norm([0.9; 0.4; 0.2]);
    Q_beta = round(Q);

    quantized = zeros(M, N);
    for i = 1:8:M
        for j = 1:8:N
            block = im(i:i+7, j:j+7);
            dct_block = dct2(block);
            quant_block = round(dct_block ./ Q_beta);
            quantized(i:i+7, j:j+7) = quant_block;
        end
    end

    quantized_r = int16(round(quantized * 100));
    quantized_int = quantized_r + 32768;
    tx_bits = de2bi(double(quantized_int(:)), 16, 'left-msb')';
    tx_bits = tx_bits(:)';

    tx_signal = ofdm_modulate(tx_bits, Nfft, CpLen);
    rx_signal = awgn_channel(tx_signal, 10);
    rx_signal = filter(h_channel, 1, rx_signal);
    rx_signal = rx_signal(1:length(tx_signal));
    rx_bits_raw = ofdm_demodulate(rx_signal, Nfft, CpLen);
    rx_bits = real(rx_bits_raw);
    rx_bits(rx_bits >= 0) = 1;
    rx_bits(rx_bits < 0) = 0;
    rx_bits = rx_bits(:)';

    min_len = min(length(tx_bits), length(rx_bits));
    rx_bits_pad = rx_bits(1:min_len);
    num_pixels = M*N;
    bits_needed = num_pixels * 16;
    if length(rx_bits_pad) < bits_needed
        rx_bits_pad = [rx_bits_pad(:); zeros(bits_needed - length(rx_bits_pad), 1)];
    end
    rx_bits_pad = rx_bits_pad(1:bits_needed);
    bits_matrix = reshape(rx_bits_pad, 16, num_pixels)';
    quant_vals = bi2de(bits_matrix, 'left-msb');
    quant_vals = double(quant_vals) - 32768;
    quant_vals = quant_vals / 100;
    quant_vals = reshape(quant_vals, M, N);
    quant_vals = round(quant_vals);

    recon_ofdm = zeros(M, N);
    for i = 1:8:M
        for j = 1:8:N
            qb = quant_vals(i:i+7, j:j+7);
            dq = qb .* Q_beta;
            rb = idct2(dq);
            recon_ofdm(i:i+7, j:j+7) = rb;
        end
    end
    recon_ofdm(recon_ofdm < 0) = 0;
    recon_ofdm(recon_ofdm > 255) = 255;
    ber = sum(tx_bits(1:min_len) ~= rx_bits_pad(1:min_len)) / min_len;
    psnr_val = 10*log10(255^2/mean((im(:)-recon_ofdm(:)).^2));
    fprintf('  OFDM @ SNR=10 dB: BER=%.2e, PSNR=%.2f dB\n', ber, psnr_val);

    figure('Position', [100, 100, 1200, 500]);
    subplot(1,3,1); imshow(uint8(im)); title('Original', 'FontSize', 14);
    subplot(1,3,2); imshow(uint8(recon_ofdm));
        title(sprintf('OFDM (SNR=10 dB)\nBER=%.2e, PSNR=%.1f dB', ber, psnr_val), 'FontSize', 12);
    subplot(1,3,3); imshow(uint8(abs(im - recon_ofdm)*10));
        title('Error Map (x10)', 'FontSize', 12);
    print(gcf, '../results/part4_ofdm_reconstruction.png', '-dpng');
    close(gcf);

    %% Part V: DSC visual reconstruction with Hamming cosets
    fprintf('Generating Part V: DSC visual reconstruction...\n');
    im_binary = im > 128;
    im_binary = im_binary(:)';

    % G = [I4 | P], H = [P' | I3]
    P = [1 1 0; 1 0 1; 0 1 1; 1 1 1];
    n_blocks = floor(length(im_binary) / 4);
    x_blocks = reshape(im_binary(1:n_blocks*4), 4, n_blocks)';

    % Generate Y with 1% noise
    noise = rand(1, n_blocks*4);
    y_bits = xor(im_binary(1:n_blocks*4), noise < 0.01);
    y_blocks = reshape(y_bits, 4, n_blocks)';
    hd = sum(im_binary(1:n_blocks*4) ~= y_bits) / (n_blocks*4);
    fprintf('  Hamming distance X vs Y: %.2f%%\n', hd*100);

    % DSC: send parity of X as syndrome
    parity_x = mod(x_blocks * P, 2);
    parity_y = mod(y_blocks * P, 2);
    synd_diff = mod(parity_x + parity_y, 2);

    % Build syndrome-to-error lookup table
    synd_lookup = zeros(8, 4);
    for col = 1:4
        p_col = P(col, :);
        synd_idx = p_col(1)*4 + p_col(2)*2 + p_col(3) + 1;
        err_vec = zeros(1, 4);
        err_vec(col) = 1;
        synd_lookup(synd_idx, :) = err_vec;
    end

    decoded_x = y_blocks;
    for i = 1:n_blocks
        synd_idx = synd_diff(i, 1)*4 + synd_diff(i, 2)*2 + synd_diff(i, 3) + 1;
        if synd_idx > 1 && synd_idx <= 8
            decoded_x(i, :) = xor(decoded_x(i, :), synd_lookup(synd_idx, :));
        end
    end

    err_after = sum(sum(x_blocks ~= decoded_x));
    accuracy = 1 - err_after / (n_blocks * 4);
    fprintf('  DSC reconstruction accuracy: %.2f%%\n', accuracy*100);

    decoded_2d = reshape(decoded_x', M, N);
    y_2d = reshape(y_blocks', M, N);
    x_2d = reshape(x_blocks', M, N);

    figure('Position', [100, 100, 1400, 450]);
    subplot(1,4,1); imshow(uint8(x_2d*255));
        title('Source X (binary)', 'FontSize', 12);
    subplot(1,4,2); imshow(uint8(y_2d*255));
        title(sprintf('Side Info Y\nHamming=%.2f%%', hd*100), 'FontSize', 12);
    subplot(1,4,3); imshow(uint8(decoded_2d*255));
        title(sprintf('Recovered X\nAccuracy=%.2f%%', accuracy*100), 'FontSize', 12);
    % Error map: white where pixel differs between source and reconstruction
    err_vec = reshape(x_blocks ~= decoded_x, 1, M*N);
    err_map = reshape(err_vec, M, N);
    error_pixels = sum(err_map(:));
    subplot(1,4,4); imshow(uint8(err_map*255));
        title(sprintf('Error Map\n%d pixels differ', error_pixels), 'FontSize', 12);
    print(gcf, '../results/part5_dsc_reconstruction.png', '-dpng');
    close(gcf);

    %% Fix Part5 m-file
    fprintf('\nFixing Part5_DistributedSourceCoding.m...\n');

    %% Copy to report directories
    fprintf('Copying to report directories...\n');
    for lang = {'english/figures', 'spanish/figuras', 'chinese/figures'}
        langdir = ['../report/' lang{1}];
        copyfile('../results/part4_ofdm_reconstruction.png', [langdir '/']);
        copyfile('../results/part5_dsc_reconstruction.png', [langdir '/']);
        fprintf('  Copied to %s/\n', langdir);
    end

    fprintf('\nDone.\n');
end
