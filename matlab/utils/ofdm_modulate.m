function tx_signal = ofdm_modulate(data_symbols, Nfft, CpLen)
    if nargin < 2
        Nfft = 64;
    end
    if nargin < 3
        CpLen = 16;
    end
    Ndata = Nfft - 2;
    num_symbols = ceil(length(data_symbols) / Ndata);
    data_padded = [data_symbols(:); zeros(Ndata * num_symbols - length(data_symbols), 1)];
    data_frames = reshape(data_padded, Ndata, num_symbols);
    ofdm_frames = zeros(Nfft, num_symbols);
    ofdm_frames(2:Ndata+1, :) = data_frames;
    time_frames = ifft(ifftshift(ofdm_frames, 1), Nfft, 1);
    time_frames_cp = [time_frames(end-CpLen+1:end, :); time_frames];
    tx_signal = time_frames_cp(:).';
end
