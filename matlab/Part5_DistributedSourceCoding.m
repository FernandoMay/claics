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

    results = struct();
    results.correlation_levels = correlation_levels;
    results.hamming_distance = zeros(length(correlation_levels), 1);
    results.correlation_rate = zeros(length(correlation_levels), 1);
    results.reconstruction_quality = zeros(length(correlation_levels), 1);
    results.slepian_wolf_rate = zeros(length(correlation_levels), 1);

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

        G = [1 0 0 0 1 1 0;
             0 1 0 0 1 0 1;
             0 0 1 0 0 1 1;
             0 0 0 1 1 1 1];
        H = [1 1 0 1 1 0 0;
             1 0 1 1 0 1 0;
             0 1 1 1 0 0 1];

        coded_x = mod(x_blocks * G, 2);
        syndrome_x = mod(coded_x * H', 2);
        coded_y = mod(y_blocks * G, 2);
        syndrome_y = mod(coded_y * H', 2);
        syndrome_diff = mod(syndrome_x + syndrome_y, 2);

        synd_lookup = zeros(8, 7);
        for err_pos = 0:7
            err_vec = zeros(1, 7);
            if err_pos > 0
                err_vec(err_pos) = 1;
            end
            test_synd = mod(err_vec * H', 2);
            synd_idx = test_synd(1)*4 + test_synd(2)*2 + test_synd(3) + 1;
            synd_lookup(synd_idx, :) = err_vec;
        end

        decoded_y_coded = coded_y;
        for i = 1:n_blocks
            synd_idx = syndrome_diff(i, 1)*4 + syndrome_diff(i, 2)*2 + syndrome_diff(i, 3) + 1;
            err_pattern = synd_lookup(synd_idx, :);
            decoded_y_coded(i, :) = xor(coded_y(i, :), err_pattern);
        end
        Ginv = G(1:4, 1:4);
        decoded_x = mod(decoded_y_coded(:, 1:4) * inv(Ginv), 2);
        decoded_x = round(decoded_x);
        decoded_x = mod(decoded_x, 2);

        sw_rate = size(syndrome_x, 2) / size(x_blocks, 2);
        results.slepian_wolf_rate(c_idx) = sw_rate;
        fprintf('Slepian-Wolf rate (syndrome/X): %.3f\n', sw_rate);

        errors_after = sum(sum(x_blocks ~= decoded_x));
        results.reconstruction_quality(c_idx) = 1 - (errors_after / (n_blocks * 4));
        fprintf('Reconstruction accuracy: %.4f\n', results.reconstruction_quality(c_idx));

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
    plot(correlation_levels, results.slepian_wolf_rate, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Noise Level (Correlation)', 'FontSize', 12);
    ylabel('Slepian-Wolf Rate (syndrome bits / data bits)', 'FontSize', 12);
    title('Slepian-Wolf Coding Rate', 'FontSize', 14);
    grid on;

    subplot(2,2,3);
    plot(correlation_levels, results.reconstruction_quality * 100, 'g-^', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Noise Level (Correlation)', 'FontSize', 12);
    ylabel('Reconstruction Accuracy (%)', 'FontSize', 12);
    title('Reconstruction Quality', 'FontSize', 14);
    grid on;

    subplot(2,2,4);
    text(0.1, 0.9, 'Slepian-Wolf Coding Principle:', 'FontSize', 12, 'FontWeight', 'bold');
    text(0.1, 0.75, 'R_X >= H(X|Y)  and  R_Y >= H(Y|X)', 'FontSize', 11);
    text(0.1, 0.60, 'Separate encoders, joint decoder', 'FontSize', 11);
    text(0.1, 0.45, 'X = Lena, Y = Lena + noise', 'FontSize', 11);
    text(0.1, 0.30, 'Hamming codes used for syndrome', 'FontSize', 11);
    text(0.1, 0.15, 'Decoder uses correlation to reconstruct', 'FontSize', 11);
    axis off;

    print(gcf, '../results/part5_distributed_source_coding.png', '-dpng');

    fprintf('\n--- Comparison of Channel Codes for DSC ---\n');
    fprintf('Hamming (7,4):  Rate=4/7=0.571,  Corrects 1 error\n');
    fprintf('LDPC (20,10):   Rate=10/20=0.50, Corrects multiple errors\n');
    fprintf('Turbo codes:    Rate~1/3-1/2,    Near Shannon limit\n');
    fprintf('Polar codes:    Rate=1/2,        Capacity-achieving\n');
    fprintf('Low correlation -> need stronger channel codes\n');
    fprintf('High correlation -> simple Hamming may suffice\n');

    fprintf('\nPart V Complete.\n\n');
end
