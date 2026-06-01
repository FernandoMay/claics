function decoded_bits = ldpc_decode(received, H, max_iter)
    if nargin < 3
        max_iter = 100;
    end
    [m, n] = size(H);
    L_ch = 2 * received(:)';
    L_q = zeros(m, n);
    L_r = zeros(m, n);
    L_Q = L_ch;
    for iter = 1:max_iter
        for i = 1:m
            idx = find(H(i, :));
            for j_idx = 1:length(idx)
                j = idx(j_idx);
                temp = L_Q(j);
                for k_idx = 1:length(idx)
                    if idx(k_idx) ~= j
                        temp = temp - L_r(i, idx(k_idx));
                    end
                    L_q(i, idx(k_idx)) = temp;
                end
                prod_alpha = 1;
                sum_phi = 0;
                for k_idx = 1:length(idx)
                    alpha_ik = sign(L_q(i, idx(k_idx)));
                    phi_ik = tanh(abs(L_q(i, idx(k_idx)))/2);
                    if idx(k_idx) ~= j
                        prod_alpha = prod_alpha * alpha_ik;
                        sum_phi = sum_phi + phi_ik;
                    end
                end
                sum_phi = min(sum_phi, 1-eps);
                L_r(i, j) = prod_alpha * log((1 + sum_phi) / (1 - sum_phi));
            end
        end
        L_Q = L_ch + sum(L_r, 1);
        hard_dec = zeros(1, n);
        hard_dec(L_Q > 0) = 1;
        parity = mod(H * hard_dec', 2);
        if sum(parity) == 0
            break;
        end
    end
    decoded_bits = hard_dec;
end
