function interleaved = interleaver(data, mode)
    if nargin < 2
        mode = 'interleave';
    end
    if strcmp(mode, 'interleave')
        n = length(data);
        cols = ceil(sqrt(n));
        rows = ceil(n / cols);
        data_padded = [data, zeros(1, rows*cols - n)];
        matrix = reshape(data_padded, rows, cols);
        interleaved = matrix.';
        interleaved = interleaved(:)';
    else
        n = length(data);
        cols = ceil(sqrt(n));
        rows = ceil(n / cols);
        data_padded = [data, zeros(1, rows*cols - n)];
        matrix = reshape(data_padded, cols, rows);
        interleaved = matrix.';
        interleaved = interleaved(:)';
    end
end
