bc7enc - Fast BC1-7 GPU texture encoders with optional Rate Distortion Optimization (RDO)

This repo contains fast texture encoders for BC1-7. All formats support a simple post-processing transform on the encoded texture data designed to trade off quality for smaller compressed file sizes using LZ compression. Significant (10-50%) size reductions are possible. The BC7 encoder also supports a "reduced entropy" mode using the -e option which causes the output to be biased/weighted in various ways which minimally impact quality, which results in 5-10% smaller file sizes with no slowdowns in encoding time.

Currently, the entropy reduction transform is tuned for Deflate, LZHAM, or LZMA. The method used to control the rate distortion tradeoff is the Lagrangian multiplier method. Rate is approximated. The post-processing transform applied to the encoded texture data tries to introduce the longest match it can into every encoded output block. It also tries to continue matches between blocks and (specifically for codecs like LZHAM/LZMA/Zstd) it tries to utilize REP0 (repeat) matches.

You can see examples of the RDO BC7 encoder's current output [here](https://richg42.blogspot.com/2021/02/more-rdo-bc7-encoding.html). Some examples on how to use the command line tool for various BC7 modes are on my blog, [here](https://richg42.blogspot.com/2021/02/how-to-use-bc7encrdo.html).

Note the BC7 encoder in bc7enc.cpp only supports modes 1/5/6/7, but the RDO post-processor function supports all BC7 modes. More modes are coming. My very high quality/fast ISPC encoder, [bc7e](https://github.com/BinomialLLC/bc7e), supports all the BC7 modes. I will be adding mode/p-bit/partition weighting to bc7e soon.

The next major focus will be improving the default smooth block handling, and adding in bc7e to get all the BC7 modes.

### Compiling

To compile (tested MSVC 2019 x64 and clang 6.0.0):

```
cmake .
make
```

Note the MSVC build enables OpenMP for faster compression.

### Examples

To encode to non-RDO BC7 using entropy reduced or quantized/weighted BC7 (super fast, slightly reduced quality, but 5-10% better LZ compression):

```
./bc7enc blah.png -e
```

To encode to RDO BC7 using the entropy reduction transform combined with reduced entropy BC7 encoding, with a slightly larger window size than the default which is 128 bytes:

```
./bc7enc -zc256 blah.png -e -z1.0
```

Same, except disable ultra-smooth block handling:

```
./bc7enc -zc256 blah.png -e -z1.0 -zu
```

To encode to RDO BC7 using the entropy reduction transform at lower quality, combined with reduced entropy BC7 encoding, with a slightly larger window size than the default which is 128 bytes:

```
./bc7enc -zc256 blah.png -e -z2.0
```

To encode to RDO BC7 using the entropy reduction transform at higher effectivenes using a larger window size, without using reduced entropy BC7 encoding:

```
./bc7enc -zc1024 blah.png -z1.0
```

To encode to RDO BC7 using the entropy reduction transform at higher effectivenes using a larger window size, with a manually specified max smooth block max error scale:

```
./bc7enc -zc1024 blah.png -z2.0 -zb30.0
```

To encode to RDO BC7 using the entropy reduction transform at higher effectivenes using a larger window size, using only mode 6 (more block artifacts, but better rate-distortion performance as measured by PSNR):

```
./bc7enc -zc1024 blah.png -6 -z1.0 -e
```

To encode to BC1:
```
./bc7enc -1 blah.png
```

To encode to BC1 with Rate Distortion Optimization (RDO) at lambda=1.0:
```
./bc7enc -1 -z1.0 blah.png
```

The -z option controls lambda, or the rate vs. distortion tradeoff. 0 = maximum quality, higher values=lower bitrates but lower quality. Try values [.25-8].

To encode to BC1 with RDO, with RDO debug output, to monitor the percentage of blocks impacted:
```
./bc7enc -1 -z1.0 -zd blah.png
```

To encode to BC1 with RDO with a higher then default smooth block scale factor (which is 10.0):
```
./bc7enc -1 -z1.0 -zb20.0 blah.png
```

Use -zb1.0 to disable smooth block error scaling completely, which increases RDO performance but can result in noticeable artifacts on smooth/flat blocks at higher lambdas.

Use -zc# to control the RDO window size in bytes. Good values are 128-8192. 
Use -zt to disable RDO multithreading. Currently, OpenMP is enabled in MSVC builds, but not in other builds.

To encode to BC1 with RDO at the highest achievable quality/effectiveness (this is extremely slower):

```
./bc7enc -1 -z1.0 -zc32768 blah.png
```

This sets the window size to 32KB (the highest setting that makes sense for Deflate). Window sizes of 2KB (the default) to 8KB are way faster and in practice are almost as effective. The maximum window size setting supported by the command line tool is 64KB, but this would be very slow.

RDO mode is supported for all the BC formats. The command line tool uses the same value of lambda for all blocks.

### Dependencies
There are no 3rd party code or library dependencies. utils.cpp/.h is only needed by the example command line tool. It uses C++11 (although it's mostly C++03).

For RDO post-processing of any block-based format: ert.cpp/.h. You'll need to supply a block decoder function for your format as a callback. It must return false if the passed in block data is invalid.

For BC1-5 encoding/decoding: rgbcx.cpp/.h

For BC7 encoding: bc7enc.cpp/.h

For BC7 decoding: bc7decomp.cpp/.h

### Features:
- Rate Distortion Optimization (RDO) using a simple transform applied to the encoded texture data.

- BC1/3 encoder (in [rgbcx.h](https://github.com/richgel999/bc7enc/blob/master/rgbcx.h)) uses a new algorithm (which we've named "prioritized cluster fit") which is 3-4x faster than traditional cluster fit (as implemented in [libsquish](https://github.com/svn2github/libsquish) with SSE2) at the same or slightly higher average quality using scalar CPU instructions. This algorithm is suitable for GPU encoder implementations.

The BC1/BC3 encoder also implements [Castano's optimal endpoint rounding improvement](https://gist.github.com/castano/c92c7626f288f9e99e158520b14a61cf).

rgbcx's BC1 encoder is faster than both AMD Compressonator and libsquish at the same average quality.

- BC7 encoder (in bc7enc.c/.h) has perceptual colorspace metric support, and is very fast compared to ispc_texcomp (see below) for RGB textures. Important: The BC7 encoder included in this repo is still a work in progress. I took bc7enc16 and added more modes for better alpha support, but it needs more testing and development.

- Full decoders for BC1-5/7. BC7 decoder is in bc7decomp.cpp/.h, BC1-5 decoders in rgbcx.h.

bc7enc currently only supports modes 1 and 6 for RGB, and modes 1, 5, 6, and 7 for alpha. The plan is to add all the modes. See the [bc7enc16](https://github.com/richgel999/bc7enc16) project for the previous version (which only supports modes 1 and 6). Note this readme still refers to "bc7enc16", but bc7enc is the same encoder but with more alpha modes.

This codec supports a perceptual mode when encoding BC7, where it computes colorspace error in
weighted YCbCr space (like etc2comp), and it also supports weighted RGBA
metrics. It's particular strong in perceptual mode, beating the current state of
the art CPU encoder (Intel's ispc_texcomp) by a wide margin when measured by
Luma PSNR, even though it only supports 2 modes and isn't vectorized.

Why only modes 1 and 6 for opaque BC7?
Because with these two modes you have a complete encoder that supports both
opaque and transparent textures in a small amount (~1400 lines) of
understandable plain C code. Mode 6 excels on smooth blocks, and mode 1 is
strong with complex blocks, and a strong encoder that combines both modes can be
quite high quality. Fast mode 6-only encoders will have noticeable block
artifacts which this codec avoids by fully supporting mode 1.

Modes 1 and 6 are typically the most used modes on many textures using other
encoders. Mode 1 has two subsets, 64 possible partitions, and 3-bit indices,
while mode 6 has large 4-bit indices and high precision 7777.1 endpoints. This
codec produces output that is far higher quality than any BC1 encoder, and
approaches (or in perceptual mode exceeds!) the quality of other full BC7
encoders.

Why is bc7enc16 so fast in perceptual mode?
Computing error in YCbCr space is more expensive than in RGB space, yet bc7enc16
in perceptual mode is stronger than ispc_texcomp (see the benchmark below) -
even without SSE/AVX vectorization and with only 2 modes to work with!

Most BC7 encoders only support linear RGB colorspace metrics, which is a
fundamental weakness. Some support weighted RGB metrics, which is better. With
linear RGB metrics, encoding error is roughly balanced between each channel, and
encoders have to work *very* hard (examining large amounts of RGB search space)
to get overall quality up. With perceptual colorspace metrics, RGB error tends
to become a bit unbalanced, with green quality favored more highly than red and
blue, and blue quality favored the least. A perceptual encoder is tuned to
prefer exploring solutions along the luma axis, where it's much less work to find
solutions with less luma error. bc7enc16 is, as far as I know, the first BC7
codec to support computing error in weighted YCbCr colorspace.

Note: Most of the timings here (except for the ispc_texcomp "fast" mode timings at the very bottom)
are for the *original* release, before I added several more optimizations. The latest version of 
bc7enc16.c is around 8-27% faster than the initial release at same quality (when mode 1 is enabled - 
there's no change with just mode 6).

Some benchmarks across 31 images (kodim corpus+others):

Perceptual (average REC709 Luma PSNR - higher is better quality):
```
iscp_texcomp slow vs. bc7enc16 uber4/max_partitions 64
iscp_texcomp:   355.4 secs 48.6 dB
bc7enc16:       122.6 secs 50.0 dB

iscp_texcomp slow vs. bc7enc16 uber0/max_partitions 64
iscp_texcomp:   355.4 secs 48.6 dB
bc7enc16:       38.3 secs 49.6 dB

iscp_texcomp basic vs. bc7enc16 uber0/max_partitions 16
ispc_texcomp:   100.2 secs 48.3 dB
bc7enc16:       20.8 secs 49.3 dB 

iscp_texcomp fast vs. bc7enc16 uber0/max_partitions 16
iscp_texcomp:   41.5 secs 48.0 dB 
bc7enc16:       20.8 secs 49.3 dB

iscp_texcomp ultrafast vs. bc7enc16 uber0/max_partitions 0
iscp_texcomp:   1.9 secs 46.2 dB
bc7enc16:       8.9 secs 48.4 dB 

Non-perceptual (average RGB PSNR):

iscp_texcomp slow vs. bc7enc16 uber4/max_partitions 64
iscp_texcomp:   355.4 secs 46.8 dB 
bc7enc16:       51 secs 46.1 dB

iscp_texcomp slow vs. bc7enc16 uber0/max_partitions 64
iscp_texcomp:   355.4 secs 46.8 dB
bc7enc16:       29.3 secs 45.8 dB

iscp_texcomp basic vs. bc7enc16 uber4/max_partitions 64
iscp_texcomp:   99.9 secs 46.5 dB
bc7enc16:       51 secs 46.1 dB

iscp_texcomp fast vs. bc7enc16 uber1/max_partitions 16
ispc_texcomp:   41.5 secs 46.1 dB
bc7enc16:       19.8 secs 45.5 dB

iscp_texcomp fast vs. bc7enc16 uber0/max_partitions 8
ispc_texcomp:   41.5 secs 46.1 dB
bc7enc16:       10.46 secs 44.4 dB

iscp_texcomp ultrafast vs. bc7enc16 uber0/max_partitions 0
ispc_texcomp:   1.9 secs 42.7 dB 
bc7enc16:       3.8 secs 42.7 dB

DirectXTex CPU in "mode 6 only" mode vs. bc7enc16 uber1/max_partions 0 (mode 6 only), non-perceptual:

DirectXTex:     466.4 secs 41.9 dB 
bc7enc16:       6.7 secs 42.8 dB

DirectXTex CPU in (default - no 3 subset modes) vs. bc7enc16 uber1/max_partions 64, non-perceptual:

DirectXTex:     9485.1 secs 45.6 dB 
bc7enc16:       36 secs 46.0 dB
```
(Note this version of DirectXTex has a key pbit bugfix which I've submitted but
is still waiting to be accepted. Non-bugfixed versions will be slightly lower
quality.)

UPDATE: To illustrate how strong the mode 1+6 implementation is in bc7enc16, let's compare ispc_texcomp 
fast vs. the latest version of bc7enc16 uber4/max_partitions 64:

Without filterbank optimizations:
```
                Time       RGB PSNR   Y PSNR
ispc_texcomp:   41.45 secs 46.09 dB   48.0 dB
bc7enc16:       41.42 secs 46.03 dB   48.2 dB

With filterbank optimizations enabled:
bc7enc16:       38.78 secs 45.94 dB   48.12 dB
```
They both have virtually the same average RGB PSNR with these settings (.06 dB is basically noise), but 
bc7enc16 is just as fast as ispc_texcomp fast, even though it's not vectorized. Interestingly, our Y PSNR is better, 
although bc7enc16 wasn't using perceptual metrics in these benchmarks. 

This was a multithreaded benchmark (using OpenMP) on a dual Xeon workstation.
ispc_texcomp was called with 64-blocks at a time and used AVX instructions.
Timings are for encoding only.
