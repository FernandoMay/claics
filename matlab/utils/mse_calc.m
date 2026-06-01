function mse_val = mse_calc(original, reconstructed)
    mse_val = mean((double(original(:)) - double(reconstructed(:))).^2);
end
