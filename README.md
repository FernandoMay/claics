# Cross-Layer Analysis of Image Communication Systems

## From JPEG-Based Compression to Wireless Transmission and Distributed Source Coding

### Project Overview

This project implements a comprehensive simulation framework for modern image communication systems, covering five major areas:

1. **Part I**: Lossy Image Coding (JPEG-like compression using DCT and quantization)
2. **Part II**: Noisy Channel Transmission (BPSK + AWGN with/without channel coding)
3. **Part III**: Scalable Image Coding (Fine Granularity Scalability via bit-plane coding)
4. **Part IV**: OFDM Wireless Transmission (frequency-selective fading channel)
5. **Part V (Option A)**: Distributed Source Coding (Slepian-Wolf coding using Hamming codes)

### Repository Structure

```
claics/
├── README.md                          # This file
├── AGENTS.md                          # Agent configuration
├── matlab/
│   ├── main.m                         # Main entry point - runs all parts
│   ├── Part1_LossyCoding.m            # Part I: Lossy JPEG-like coding
│   ├── Part2_NoisyChannel.m           # Part II: Noisy channel transmission
│   ├── Part3_ScalableCoding.m         # Part III: Scalable (FGS) coding
│   ├── Part4_OFDM_Transmission.m      # Part IV: OFDM wireless transmission
│   ├── Part5_DistributedSourceCoding.m # Part V: Slepian-Wolf DSC
│   └── utils/
│       ├── awgn_channel.m             # AWGN channel model
│       ├── bpsk_modulate.m            # BPSK modulation
│       ├── bpsk_demodulate.m          # BPSK demodulation (hard decision)
│       ├── compression_ratio.m        # Compression ratio calculator
│       ├── hamming_encode.m           # Hamming (7,4) encoder
│       ├── hamming_decode.m           # Hamming (7,4) decoder
│       ├── interleaver.m              # Block interleaver/deinterleaver
│       ├── izigzag.m                  # Inverse zigzag scan
│       ├── ldpc_encode.m              # LDPC (20,10) encoder
│       ├── ldpc_decode.m              # LDPC decoder (BP algorithm)
│       ├── mse_calc.m                 # Mean Square Error calculator
│       ├── ofdm_demodulate.m          # OFDM demodulator
│       ├── ofdm_modulate.m            # OFDM modulator
│       ├── psnr_calc.m               # PSNR calculator
│       └── zigzag.m                   # Zigzag scan for 8x8 blocks
├── report/
│   ├── english/
│   │   └── report_en.tex             # IEEE report in English
│   ├── spanish/
│   │   └── report_es.tex             # IEEE report in Spanish
│   └── chinese/
│       └── report_zh.tex             # IEEE report in Chinese
└── results/                           # Generated figures and data
```

### Requirements

- MATLAB R2018b or newer
- Image Processing Toolbox (for `dct2`, `idct2`, `imread`, `imresize`)
- Communications Toolbox (optional, for enhanced LDPC)
- Lena image (512x512) - if not available, the script generates a test image

### Quick Start

1. Open MATLAB and navigate to the `matlab/` directory:
   ```matlab
   cd /path/to/claics/matlab
   ```

2. Run the main simulation:
   ```matlab
   main
   ```

3. Individual parts can be run separately:
   ```matlab
   Part1_LossyCoding(im)
   Part2_NoisyChannel(im, 0:2:10)
   Part3_ScalableCoding(im, 5, 10)
   Part4_OFDM_Transmission(im, 0:2:10)
   Part5_DistributedSourceCoding(im, [0.01, 0.05, 0.15])
   ```

### Output

All generated figures are saved to the `results/` directory:
- `part1_compression_vs_psnr.png` - Compression ratio vs PSNR curves
- `part1_original_vs_reconstructed.png` - Original vs reconstructed images
- `part2_ber_psnr_vs_snr.png` - BER and PSNR vs SNR for 3 schemes
- `part3_scalable_coding.png` - Scalable coding results for weak/strong users
- `part4_ofdm_performance.png` - OFDM BER/PSNR with/without coding/interleaver
- `part5_distributed_source_coding.png` - DSC performance analysis

### LaTeX Reports

IEEE-format reports are available in three languages:
- `report/english/report_en.tex`
- `report/spanish/report_es.tex`
- `report/chinese/report_zh.tex`

Compile with:
```bash
pdflatex report_en.tex
bibtex report_en
pdflatex report_en.tex
pdflatex report_en.tex
```

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| Beta (β) | 0.5, 1, 2, 4, 8, 16 | Quantization multiplier |
| SNR range | 0-10 dB | Signal-to-noise ratio |
| Nfft | 64 | OFDM FFT size |
| CpLen | 16 | Cyclic prefix length |
| Channel h | [0.9, 0.4, 0.2] | Multipath channel taps |

### Academic Context

This project was developed for a graduate-level course in Modern Image Communication. It covers:

- **Signal Processing**: DCT, quantization, image reconstruction
- **Information Theory**: Entropy, channel capacity, Slepian-Wolf theorem
- **Channel Coding**: Hamming codes, LDPC codes
- **Digital Communications**: BPSK, AWGN, OFDM
- **Source Coding**: Huffman coding, run-length coding, bit-plane coding
- **Wireless Communications**: Frequency-selective fading, interleaving

### References

1. G. K. Wallace, "The JPEG still picture compression standard," IEEE Trans. Consumer Electronics, 1992.
2. T. Cover and J. Thomas, *Elements of Information Theory*, 2nd ed., Wiley, 2006.
3. J. Proakis and M. Salehi, *Digital Communications*, 5th ed., McGraw-Hill, 2008.
4. T. Richardson and R. Urbanke, *Modern Coding Theory*, Cambridge Univ. Press, 2008.
5. D. Slepian and J. Wolf, "Noiseless coding of correlated information sources," IEEE Trans. IT, 1973.

### License

Academic project - for educational purposes only.
