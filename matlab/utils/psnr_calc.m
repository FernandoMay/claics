function psnr_val = psnr_calc(original, reconstructed)
    mse_val = mean((double(original(:)) - double(reconstructed(:))).^2);
    if mse_val == 0
        psnr_val = inf;
    else
        psnr_val = 10 * log10(255^2 / mse_val);
    end
end
