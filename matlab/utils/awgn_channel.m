function output = awgn_channel(signal, snr_dB)
    signal_power = mean(abs(signal).^2);
    snr_linear = 10^(snr_dB/10);
    noise_power = signal_power / snr_linear;
    noise = sqrt(noise_power/2) * (randn(size(signal)) + 1j*randn(size(signal)));
    if isreal(signal)
        noise = real(noise);
    end
    output = signal + noise;
end
