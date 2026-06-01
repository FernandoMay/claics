function coded_bits = hamming_encode(data_bits)
    n = 7; k = 4;
    num_blocks = floor(length(data_bits) / k);
    data_bits = data_bits(1:num_blocks*k);
    data_blocks = reshape(data_bits, k, num_blocks)';
    G = [1 0 0 0 1 1 0;
         0 1 0 0 1 0 1;
         0 0 1 0 0 1 1;
         0 0 0 1 1 1 1];
    coded_blocks = mod(data_blocks * G, 2);
    coded_bits = coded_blocks(:)';
end
