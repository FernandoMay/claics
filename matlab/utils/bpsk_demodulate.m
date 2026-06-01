function bits = bpsk_demodulate(symbols)
    bits = zeros(size(symbols));
    bits(symbols >= 0) = 1;
    bits(symbols < 0) = 0;
end
