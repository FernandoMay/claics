%% ========================================================================
%  Cross-Layer Analysis of Image Communication Systems
%  From JPEG-Based Compression to Wireless Transmission
%  and Distributed Source Coding
%  ========================================================================
%  This script executes all five parts of the project:
%   Part I:   Lossy Image Coding (JPEG-like)
%   Part II:  Noisy Channel Transmission (BPSK + AWGN)
%   Part III: Scalable Image Coding (FGS)
%   Part IV:  OFDM Wireless Transmission
%   Part V:   Distributed Source Coding (Slepian-Wolf)
%  ========================================================================

clear; close all; clc;
fprintf('============================================================\n');
fprintf('Cross-Layer Analysis of Image Communication Systems\n');
fprintf('============================================================\n\n');

addpath('utils');
if ~exist('../results', 'dir')
    mkdir('../results');
end

%% ------------------------------------------------------------------------
%  Load Lena image (512x512)
%  ------------------------------------------------------------------------
%  If you have lena512.mat, load it.
%  Otherwise, generate a synthetic test image.

fprintf('Loading Lena image...\n');
if exist('lena512.mat', 'file')
    load('lena512.mat', 'im');
    fprintf('Loaded lena512.mat (512x512)\n');
elseif exist('../lena512.mat', 'file')
    load('../lena512.mat', 'im');
    fprintf('Loaded ../lena512.mat (512x512)\n');
else
    fprintf('lena512.mat not found. Generating synthetic image...\n');
    im = imread('cameraman.tif');
    im = imresize(im, [512, 512]);
    im = im2double(im);
    fprintf('Generated 512x512 test image\n');
end

[M, N] = size(im);
fprintf('Image size: %d x %d\n\n', M, N);

%% ========================================================================
%  PART I: Lossy Image Coding
%  ========================================================================
fprintf('============================================================\n');
fprintf('PART I: Lossy Image Coding (JPEG-like)\n');
fprintf('============================================================\n\n');

beta_values = [0.5, 1, 2, 4, 8, 16];
results1 = Part1_LossyCoding(im, beta_values);

%% ========================================================================
%  PART II: Noisy Channel Transmission
%  ========================================================================
fprintf('============================================================\n');
fprintf('PART II: Noisy Channel Transmission\n');
fprintf('============================================================\n\n');

snr_range_2 = 0:2:10;
results2 = Part2_NoisyChannel(im, snr_range_2);

%% ========================================================================
%  PART III: Scalable Image Coding
%  ========================================================================
fprintf('============================================================\n');
fprintf('PART III: Scalable Image Coding (FGS)\n');
fprintf('============================================================\n\n');

snr_weak = 5;
snr_diff = 10;
results3 = Part3_ScalableCoding(im, snr_weak, snr_diff);

%% ========================================================================
%  PART IV: OFDM Wireless Transmission
%  ========================================================================
fprintf('============================================================\n');
fprintf('PART IV: OFDM Wireless Transmission\n');
fprintf('============================================================\n\n');

snr_range_4 = 0:2:10;
results4 = Part4_OFDM_Transmission(im, snr_range_4);

%% ========================================================================
%  PART V: Distributed Source Coding (Slepian-Wolf - Option A)
%  ========================================================================
fprintf('============================================================\n');
fprintf('PART V: Distributed Source Coding\n');
fprintf('============================================================\n\n');

correlation_levels = [0.01, 0.05, 0.15];
results5 = Part5_DistributedSourceCoding(im, correlation_levels);

%% ========================================================================
%  Summary
%  ========================================================================
fprintf('============================================================\n');
fprintf('SIMULATION COMPLETE\n');
fprintf('============================================================\n');
fprintf('\nResults Summary:\n');
fprintf('  Part I:   Compression ratios from %.1f to %.1f\n', ...
    min(results1.cr), max(results1.cr));
fprintf('            PSNR from %.2f to %.2f dB\n', ...
    min(results1.psnr), max(results1.psnr));
fprintf('  Part II:  BER range: %.2e to %.2e (Scheme 1)\n', ...
    min(results2.ber_scheme1), max(results2.ber_scheme1));
fprintf('            BER range: %.2e to %.2e (Scheme 3)\n', ...
    min(results2.ber_scheme3), max(results2.ber_scheme3));
fprintf('  Part III: Weak user PSNR = %.2f dB\n', results3.psnr_weak);
fprintf('            Strong user PSNR = %.2f dB\n', results3.psnr_strong);
fprintf('  Part IV:  OFDM BER with interleaver: %.2e to %.2e\n', ...
    min(results4.ber_channel_interleaved), max(results4.ber_channel_interleaved));
fprintf('  Part V:   DSC accuracy from %.1f%% to %.1f%%\n', ...
    min(results5.reconstruction_quality)*100, max(results5.reconstruction_quality)*100);
fprintf('\nAll results saved to ../results/\n');
fprintf('============================================================\n');
