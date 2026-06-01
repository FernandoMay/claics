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

    im = double(im);

    Q = [16 11 10 16 24 40 51 61;
         12 12 14 19 26 58 60 55;
         14 13 16 24 40 57 69 56;
         14 17 22 29 51 87 80 62;
         18 22 37 56 68 109 103 77;
         24 35 55 64 81 104 113 92;
         49 64 78 87 103 121 120 101;
         72 92 95 98 113 100 103 99];

    [M, N] = size(im);

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

        quantized = zeros(M, N);
        for i = 1:8:M
            for j = 1:8:N
                block = im(i:i+7, j:j+7);
                dct_block = dct2(block);
                quant_block = round(dct_block ./ Q_beta);
                quantized(i:i+7, j:j+7) = quant_block;
            end
        end

        compressed_size = 0;
        for i = 1:8:M
            for j = 1:8:N
                qblock = quantized(i:i+7, j:j+7);
                zz = zigzag(qblock);
                compressed_size = compressed_size + nnz(zz);
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
        reconstructed(reconstructed > 255) = 255;

        results.mse(b_idx) = mean((im(:) - reconstructed(:)).^2);
        if results.mse(b_idx) == 0
            results.psnr(b_idx) = inf;
        else
            results.psnr(b_idx) = 10 * log10(255^2 / results.mse(b_idx));
        end

        original_size = M * N * 8;
        results.compressed_sizes(b_idx) = compressed_size;
        results.cr(b_idx) = original_size / max(compressed_size, 1);

        fprintf('Beta = %5.1f | MSE = %8.4f | PSNR = %6.2f dB | Compressed = %d bits | CR = %6.2f\n', ...
            beta, results.mse(b_idx), results.psnr(b_idx), compressed_size, results.cr(b_idx));

        if any(b_idx == [1, 2, 4, 6])
            filename = sprintf('../results/reconstructed_beta_%.1f.png', beta);
            imwrite(uint8(reconstructed), filename);
            fprintf('  -> Saved: %s\n', filename);
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
    print(gcf, '../results/part1_compression_vs_psnr.png', '-dpng');

    figure('Name', 'Part I: Original vs Reconstructed', 'Position', [100, 100, 1200, 400]);
    subplot(1,3,1);
    imshow(uint8(im));
    title('Original Image', 'FontSize', 12);

    Q_beta = round(1 * Q);
    Q_beta(Q_beta < 1) = 1;
    recon1 = zeros(M, N);
    for i = 1:8:M
        for j = 1:8:N
            block = im(i:i+7, j:j+7);
            dct_block = dct2(block);
            quant_block = round(dct_block ./ Q_beta);
            dequant_block = quant_block .* Q_beta;
            recon_block = idct2(dequant_block);
            recon1(i:i+7, j:j+7) = recon_block;
        end
    end
    recon1(recon1 < 0) = 0;
    recon1(recon1 > 255) = 255;
    subplot(1,3,2);
    imshow(uint8(recon1));
    mse1 = mean((im(:)-recon1(:)).^2);
    title(sprintf('Reconstructed (\\beta=1)\nPSNR=%.2f dB', 10*log10(255^2/mse1)), 'FontSize', 12);

    Q_beta8 = round(8 * Q);
    Q_beta8(Q_beta8 < 1) = 1;
    recon8 = zeros(M, N);
    for i = 1:8:M
        for j = 1:8:N
            block = im(i:i+7, j:j+7);
            dct_block = dct2(block);
            quant_block = round(dct_block ./ Q_beta8);
            dequant_block = quant_block .* Q_beta8;
            recon_block = idct2(dequant_block);
            recon8(i:i+7, j:j+7) = recon_block;
        end
    end
    recon8(recon8 < 0) = 0;
    recon8(recon8 > 255) = 255;
    subplot(1,3,3);
    imshow(uint8(recon8));
    mse8 = mean((im(:)-recon8(:)).^2);
    title(sprintf('Reconstructed (\\beta=8)\nPSNR=%.2f dB', 10*log10(255^2/mse8)), 'FontSize', 12);

    print(gcf, '../results/part1_original_vs_reconstructed.png', '-dpng');
    close(gcf);

    fprintf('\nPart I Complete.\n\n');
end
