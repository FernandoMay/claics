function results = Part5_DistributedSourceCoding(im, correlation_levels)
    if nargin < 1
        load('lena512.mat', 'im');
    end
    if nargin < 2
        correlation_levels = [0.01, 0.05, 0.15];
    end

    addpath('utils');

    fprintf('========================================\n');
    fprintf('PART V: Distributed Source Coding\n');
    fprintf('       (Slepian-Wolf Coding)\n');
    fprintf('========================================\n');

    fprintf('\n--- Slepian-Wolf Theorem ---\n');
    fprintf('For two correlated sources X and Y:\n');
    fprintf('  R_X >= H(X|Y)\n');
    fprintf('  R_Y >= H(Y|X)\n');
    fprintf('  R_X + R_Y >= H(X,Y)\n');
    fprintf('Separate encoding, joint decoding.\n');

    fprintf('\n--- Hamming Distance ---\n');
    fprintf('d_H(x,y) = sum_i (x_i XOR y_i)\n');
    fprintf('Measures correlation between two binary sources.\n');

    im = double(im);
    im_binary = im > 128;
    im_binary = im_binary(:)';

    % Hamming(7,4) parity submatrix P (4x3)
    % G = [I4 | P],  H = [P' | I3]
    P = [1 1 0;
         1 0 1;
         0 1 1;
         1 1 1];

    results = struct();
    results.correlation_levels = correlation_levels;
    results.hamming_distance = zeros(length(correlation_levels), 1);
    results.correlation_rate = zeros(length(correlation_levels), 1);
    results.reconstruction_quality = zeros(length(correlation_levels), 1);
    results.slepian_wolf_rate = zeros(length(correlation_levels), 1);

    % Build syndrome-to-error lookup for Hamming(7,4) DSC
    % Each column of P maps to a syndrome that identifies a 1-bit error
    synd_lookup = zeros(8, 4);
    for col = 1:4
        p_col = P(col, :);
        synd_idx = p_col(1)*4 + p_col(2)*2 + p_col(3) + 1;
        err_vec = zeros(1, 4);
        err_vec(col) = 1;
        synd_lookup(synd_idx, :) = err_vec;
    end

    for c_idx = 1:length(correlation_levels)
        corr_level = correlation_levels(c_idx);
        fprintf('\n--- Correlation Level: %.3f ---\n', corr_level);

        noise = rand(size(im_binary));
        y_binary = xor(im_binary, noise < corr_level);
        hamming_dist = sum(im_binary ~= y_binary);
        results.hamming_distance(c_idx) = hamming_dist / length(im_binary);
        fprintf('Hamming distance: %d (%.4f%%)\n', hamming_dist, 100*results.hamming_distance(c_idx));

        n_blocks = floor(length(im_binary) / 4);
        x_blocks = reshape(im_binary(1:n_blocks*4), 4, n_blocks)';
        y_blocks = reshape(y_binary(1:n_blocks*4), 4, n_blocks)';

        % DSC: send only parity bits of X as syndrome (3 bits per 4 data bits)
        parity_x = mod(x_blocks * P, 2);
        parity_y = mod(y_blocks * P, 2);
        synd_diff = mod(parity_x + parity_y, 2);

        % Joint decoder: use syndrome difference to correct Y
        decoded_x = y_blocks;
        for i = 1:n_blocks
            synd_idx = synd_diff(i, 1)*4 + synd_diff(i, 2)*2 + synd_diff(i, 3) + 1;
            if synd_idx > 1 && synd_idx <= 8
                decoded_x(i, :) = xor(decoded_x(i, :), synd_lookup(synd_idx, :));
            end
        end

        sw_rate = 3 / 4;
        results.slepian_wolf_rate(c_idx) = sw_rate;
        fprintf('Slepian-Wolf rate (syndrome/X): %.3f (parity bits / data bits)\n', sw_rate);

        errors_after = sum(sum(x_blocks ~= decoded_x));
        results.reconstruction_quality(c_idx) = 1 - (errors_after / (n_blocks * 4));
        fprintf('Reconstruction accuracy: %.4f (errors corrected: %d -> %d)\n', ...
            results.reconstruction_quality(c_idx), hamming_dist, errors_after);

        results.correlation_rate(c_idx) = corr_level;
    end

    figure('Name', 'Part V: Distributed Source Coding', 'Position', [100, 100, 1000, 600]);

    subplot(2,2,1);
    plot(correlation_levels, results.hamming_distance * 100, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Noise Level (Correlation)', 'FontSize', 12);
    ylabel('Hamming Distance (%)', 'FontSize', 12);
    title('Hamming Distance vs Correlation', 'FontSize', 14);
    grid on;

    subplot(2,2,2);
    plot(correlation_levels, results.reconstruction_quality * 100, 'g-^', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Noise Level (Correlation)', 'FontSize', 12);
    ylabel('Reconstruction Accuracy (%)', 'FontSize', 12);
    title('DSC: Reconstruction Quality', 'FontSize', 14);
    grid on;

    subplot(2,2,3);
    bar(correlation_levels, [results.hamming_distance' * 100; (1-results.reconstruction_quality)' * 100]');
    xlabel('Noise Level', 'FontSize', 12);
    ylabel('Error Rate (%)', 'FontSize', 12);
    legend('Before correction', 'After correction', 'Location', 'northwest');
    title('DSC Error Correction Performance', 'FontSize', 14);
    grid on;

    subplot(2,2,4);
    text(0.1, 0.9, 'Slepian-Wolf Coding with Hamming(7,4):', 'FontSize', 12, 'FontWeight', 'bold');
    text(0.1, 0.75, 'G = [I4 | P],  syndrome = parity bits', 'FontSize', 11);
    text(0.1, 0.60, 'TX sends: 3 syndrome bits / 4 data bits', 'FontSize', 11);
    text(0.1, 0.45, 'Rate: 3/4 = 0.75 (compression)', 'FontSize', 11);
    text(0.1, 0.30, 'Decoder: Y + syndrome -> corrects 1-bit errors', 'FontSize', 11);
    text(0.1, 0.15, 'Higher correlation -> better reconstruction', 'FontSize', 11);
    axis off;

    print(gcf, '../results/part5_distributed_source_coding.png', '-dpng');
    close(gcf);

    fprintf('\n--- Comparison of Channel Codes for DSC ---\n');
    fprintf('Hamming (7,4):  Rate=4/7=0.571,  Corrects 1 error/block\n');
    fprintf('LDPC (20,10):   Rate=10/20=0.50, Corrects multiple errors\n');
    fprintf('Turbo codes:    Rate~1/3-1/2,    Near Shannon limit\n');
    fprintf('Polar codes:    Rate=1/2,        Capacity-achieving\n');
    fprintf('Low correlation -> need stronger channel codes\n');
    fprintf('High correlation -> simple Hamming may suffice\n');

    fprintf('\nPart V Complete.\n\n');
end
