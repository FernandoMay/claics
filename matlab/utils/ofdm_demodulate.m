function rx_symbols = ofdm_demodulate(rx_signal, Nfft, CpLen)
    if nargin < 2
        Nfft = 64;
    end
    if nargin < 3
        CpLen = 16;
    end
    Ndata = Nfft - 2;
    total_symbols = floor(length(rx_signal) / (Nfft + CpLen));
    rx_frames = reshape(rx_signal(1:total_symbols*(Nfft+CpLen)), Nfft+CpLen, total_symbols);
    rx_frames = rx_frames(CpLen+1:end, :);
    freq_frames = fftshift(fft(rx_frames, Nfft, 1), 1);
    data_frames = freq_frames(2:Ndata+1, :);
    rx_symbols = data_frames(:).';
end
