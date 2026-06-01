function coded_bits = ldpc_encode(data_bits)
    n = 20; k = 10;
    H = [
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
    num_blocks = floor(length(data_bits) / k);
    data_bits = data_bits(1:num_blocks*k);
    data_blocks = reshape(data_bits, k, num_blocks)';
    H1 = H(:, 1:k);
    H2 = H(:, k+1:end);
    P = (H2 \ H1)';
    P = mod(round(P), 2);
    G = [eye(k), P];
    coded_blocks = mod(data_blocks * G, 2);
    coded_bits = coded_blocks(:)';
end
