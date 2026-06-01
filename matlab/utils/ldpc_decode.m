function decoded_bits = ldpc_decode(received_signal, H, max_iter)
    if nargin < 3
        max_iter = 50;
    end
    [m, n] = size(H);
    LLR = 2 * received_signal / mean(abs(received_signal))^2;
    L_q = zeros(m, n);
    L_r = zeros(m, n);
    L_total = LLR;
    for iter = 1:max_iter
        for i = 1:m
            cols = find(H(i, :));
            for j_idx = 1:length(cols)
                j = cols(j_idx);
                temp = L_total(j);
                for k_idx = 1:length(cols)
                    if cols(k_idx) ~= j
                        temp = temp + L_q(i, cols(k_idx));
                    end
                end
                L_r(i, j) = temp;
            end
        end
        for j = 1:n
            rows = find(H(:, j));
            for i_idx = 1:length(rows)
                i = rows(i_idx);
                temp = 0;
                for k_idx = 1:length(rows)
                    if rows(k_idx) ~= i
                        temp = temp + L_r(rows(k_idx), j);
                    end
                end
                L_q(i, j) = temp;
            end
        end
        L_total = LLR + sum(L_r(:, :), 1);
        decisions = zeros(1, n);
        decisions(L_total > 0) = 1;
        if mod(H * decisions', 2) == 0
            break;
        end
    end
    decoded_bits = decisions;
end
