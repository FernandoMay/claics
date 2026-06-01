function results = Part1_LossyCoding(im, beta_values)
    if nargin < 1
        load('lena512.mat', 'im');
    end
    if nargin < 2
        beta_values = [0.5, 1, 2, 4, 8, 16];
    end

    addpath('utils');
    fprintf('========================================\n');
    fprintf('PART I: Lossy Image Coding (JPEG-like)\n');
    fprintf('========================================\n');

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

    results = struct();
    results.beta = beta_values;
    results.psnr = zeros(length(beta_values), 1);
    results.mse = zeros(length(beta_values), 1);
    results.cr = zeros(length(beta_values), 1);
    results.compressed_sizes = zeros(length(beta_values), 1);

    for b_idx = 1:length(beta_values)
        beta = beta_values(b_idx);
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

        compressed_size = 0;
        for i = 1:8:M
            for j = 1:8:N
                qblock = quantized(i:i+7, j:j+7);
                zz = zigzag(qblock);
                nonzero_count = sum(zz ~= 0);
                compressed_size = compressed_size + nonzero_count;
            end
        end

        reconstructed = zeros(M, N);
        for i = 1:8:M
            for j = 1:8:N
                quant_block = quantized(i:i+7, j:j+7);
                dequant_block = quant_block .* Q_beta;
                recon_block = idct2(dequant_block);
                reconstructed(i:i+7, j:j+7) = recon_block;
            end
        end

        reconstructed(reconstructed < 0) = 0;
        reconstructed(reconstructed > 1) = 1;

        results.mse(b_idx) = mean((im(:) - reconstructed(:)).^2);
        if results.mse(b_idx) == 0
            results.psnr(b_idx) = inf;
        else
            results.psnr(b_idx) = 10 * log10(1 / results.mse(b_idx));
        end

        original_size = M * N * 8;
        results.compressed_sizes(b_idx) = compressed_size;
        results.cr(b_idx) = original_size / compressed_size;

        fprintf('Beta = %5.1f | MSE = %8.4f | PSNR = %6.2f dB | Compressed = %d bits | CR = %6.2f\n', ...
            beta, results.mse(b_idx), results.psnr(b_idx), compressed_size, results.cr(b_idx));

        if b_idx == 1 || b_idx == length(beta_values) || beta == 1 || beta == 4
            filename = sprintf('../results/reconstructed_beta_%.1f.png', beta);
            imwrite(reconstructed, filename);
            fprintf('  -> Saved reconstructed image: %s\n', filename);
        end
    end

    figure('Name', 'Part I: Compression Ratio vs PSNR', 'Position', [100, 100, 800, 600]);
    subplot(2,1,1);
    plot(results.cr, results.psnr, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Compression Ratio', 'FontSize', 12);
    ylabel('PSNR (dB)', 'FontSize', 12);
    title('Part I: Compression Ratio vs PSNR', 'FontSize', 14);
    grid on;
    for i = 1:length(beta_values)
        text(results.cr(i), results.psnr(i)+0.5, sprintf('\\beta=%.1f', beta_values(i)), ...
            'FontSize', 10, 'HorizontalAlignment', 'center');
    end

    subplot(2,1,2);
    semilogx(results.cr, results.psnr, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Compression Ratio (log scale)', 'FontSize', 12);
    ylabel('PSNR (dB)', 'FontSize', 12);
    title('PSNR vs Compression Ratio (semi-log)', 'FontSize', 14);
    grid on;

    saveas(gcf, '../results/part1_compression_vs_psnr.png');

    figure('Name', 'Part I: Original vs Reconstructed', 'Position', [100, 100, 1200, 400]);
    subplot(1,3,1);
    imshow(im);
    title('Original Image', 'FontSize', 12);

    recon_beta1 = zeros(M, N);
    for i = 1:8:M
        for j = 1:8:N
            block = im(i:i+7, j:j+7);
            dct_block = dct2(block);
            quant_block = round(dct_block ./ round(1*Q));
            dequant_block = quant_block .* round(1*Q);
            recon_block = idct2(dequant_block);
            recon_beta1(i:i+7, j:j+7) = recon_block;
        end
    end
    recon_beta1(recon_beta1 < 0) = 0;
    recon_beta1(recon_beta1 > 1) = 1;
    subplot(1,3,2);
    imshow(recon_beta1);
    title(sprintf('Reconstructed (\\beta=1)\nPSNR=%.2f dB', 10*log10(1/mean((im(:)-recon_beta1(:)).^2))), 'FontSize', 12);

    recon_beta8 = zeros(M, N);
    for i = 1:8:M
        for j = 1:8:N
            block = im(i:i+7, j:j+7);
            dct_block = dct2(block);
            quant_block = round(dct_block ./ round(8*Q));
            dequant_block = quant_block .* round(8*Q);
            recon_block = idct2(dequant_block);
            recon_beta8(i:i+7, j:j+7) = recon_block;
        end
    end
    recon_beta8(recon_beta8 < 0) = 0;
    recon_beta8(recon_beta8 > 1) = 1;
    subplot(1,3,3);
    imshow(recon_beta8);
    title(sprintf('Reconstructed (\\beta=8)\nPSNR=%.2f dB', 10*log10(1/mean((im(:)-recon_beta8(:)).^2))), 'FontSize', 12);

    saveas(gcf, '../results/part1_original_vs_reconstructed.png');

    fprintf('\nPart I Complete.\n\n');
end
