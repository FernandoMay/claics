function decoded_bits = hamming_decode(received_bits)
    n = 7; k = 4;
    num_blocks = floor(length(received_bits) / n);
    received_bits = received_bits(1:num_blocks*n);
    recv_blocks = reshape(received_bits, n, num_blocks)';
    H = [1 1 0 1 1 0 0;
         1 0 1 1 0 1 0;
         0 1 1 1 0 0 1];
    syndrome = mod(recv_blocks * H', 2);
    error_table = [
        0 0 0 0 0 0 0;
        0 0 0 0 0 0 1;
        0 0 0 0 0 1 0;
        0 0 0 1 0 0 0;
        0 0 0 0 1 0 0;
        0 0 1 0 0 0 0;
        0 1 0 0 0 0 0;
        1 0 0 0 0 0 0
    ];
    for i = 1:num_blocks
        synd_vec = bi2de(syndrome(i,:), 'left-msb');
        if synd_vec > 0 && synd_vec <= 7
            recv_blocks(i,:) = xor(recv_blocks(i,:), error_table(synd_vec+1,:));
        end
    end
    decoded_bits = recv_blocks(:, 1:k)';
    decoded_bits = decoded_bits(:)';
end
