// test.cpp - Command line example/test app
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#define BC7ENC_VERSION "1.07"

#define LZHAM_STATS (0)
#define DECODE_BC4_TO_GRAYSCALE (0)
#define COMPUTE_SSIM (0)

#ifndef SUPPORT_BC7E
#define SUPPORT_BC7E (0)
#endif

#if _OPENMP
#include <omp.h>
#endif

#include "utils.h"
#include "ert.h"

#include "bc7enc.h"
#include "bc7decomp.h"

#define RGBCX_IMPLEMENTATION
#include "rgbcx.h"

#include "miniz.h"

#if SUPPORT_BC7E
#include "bc7e_ispc.h"
#endif

using namespace utils;

static int print_usage()
{
	fprintf(stderr, "Reads PNG files (with or without alpha channels) and packs them to BC1-5 or BC7/BPTC (default) using\nmodes 1, 6 (opaque blocks) or modes 1, 5, 6, and 7 (alpha blocks).\n");
	fprintf(stderr, "Supports optional reduced entropy BC7 encoding (using -e) and Rate Distortion Optimization (RDO) for BC1-7 (using -z# where # is lambda).\n");
	fprintf(stderr, "By default, this tool compresses to BC7. A DX10 DDS file and a unpacked PNG file will be written to the source\ndirectory with the .dds/_unpacked.png/_unpacked_alpha.png suffixes.\n");
	fprintf(stderr, "This tool does not yet support generating mipmaps (yet).\n");
	fprintf(stderr, "\nUsage: bc7enc [-apng_filename] [options] input_filename.png [compressed_output.dds] [unpacked_output.png]\n\n");
	fprintf(stderr, "-apng_filename Load G channel of PNG file into alpha channel of source image\n");
	fprintf(stderr, "-g Don't write unpacked output PNG files (this disables PSNR metrics too).\n");
	fprintf(stderr, "-y Flip source image along Y axis before packing\n");
	fprintf(stderr, "-o Write output files to the source file's directory\n");
	fprintf(stderr, "-1 Encode to BC1. -u[0,5] controls quality vs. perf. tradeoff for RGB.\n");
	fprintf(stderr, "-3 Encode to BC3. -u[0,5] controls quality vs. perf. tradeoff for RGB.\n");
	fprintf(stderr, "-4 Encode to BC4\n");
	fprintf(stderr, "-5 Encode to BC5\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "-X# BC4/5: Set first color channel (defaults to 0 or red)\n");
	fprintf(stderr, "-Y# BC4/5: Set second color channel (defaults to 1 or green)\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "-s BC7: Use perceptual colorspace metrics instead of linear. The default for all formats is to use linear RGB/RGBA metrics. BC7 RDO mode is always linear.\n");
	fprintf(stderr, "-uX BC7: Higher BC7 quality levels, X ranges from [0,4] for BC7. Default is 4.\n");
	fprintf(stderr, "-pX BC7: Scan X partitions in mode 1, X ranges from [0,64], use 0 to disable mode 1 entirely (faster)\n");
	fprintf(stderr, "-LX BC1: Set encoding level, where 0=fastest and 18=slowest but highest quality. Default is 18.\n");
	fprintf(stderr, "\nBC3-5 alpha block encoding options:\n");
	fprintf(stderr, "-hl BC3-5: Use lower quality BC4 block encoder (much faster, but lower quality, only uses 8 value mode)\n");
	fprintf(stderr, "-h6 BC3-5: Use 6 value mode only for BC4 blocks\n");
	fprintf(stderr, "-h8 BC3-5: Use 8 value mode only for BC4 blocks\n");
	fprintf(stderr, "-hr# BC3-5: Set search radius, default is 5, larger=higher quality but slower compression\n");
	fprintf(stderr, "\nRDO encoding options:\n");
	fprintf(stderr, "-e BC7: Quantize/weight BC7 output for lower entropy (no slowdown but only 5-10%% gains, can be combined with -z# for more gains)\n");
	fprintf(stderr, "-z# BC1-7: Set RDO lambda factor (quality), lower=higher quality/larger LZ compressed files, try .1-4, combine with -e for BC7 for more gains\n");
	fprintf(stderr, "-zb# BC1-7: Manually set smooth block scale factor, higher values = less distortion on smooth blocks, try 5-70\n");
	fprintf(stderr, "-zc# BC1: Set RDO lookback window size in bytes (higher=more effective but slower, default=128, try 64-16384)\n");
	fprintf(stderr, "-zn BC1-7: Inject up to 2 matches into each block vs. 1 (a little slower, but noticeably higher compression)\n");
	fprintf(stderr, "-zm BC1-7: Allow byte sequences to be moved inside blocks (much slower)\n");
	fprintf(stderr, "-zu BC1/3/7: Disable RGB ultrasmooth block detection/handling\n");
	fprintf(stderr, "RDO debugging/development:\n");
	fprintf(stderr, "-zd BC1-7: Enable debug output\n");
	fprintf(stderr, "-zt BC1-7: Disable RDO multithreading\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "-b BC1: Don't use 3-color mode transparent texels on blocks containing black or very dark pixels. By default this mode is now enabled.\n");
	fprintf(stderr, "-c BC1: Disable 3-color mode\n");
	fprintf(stderr, "-n BC1: Encode/decode for NVidia GPU's\n");
	fprintf(stderr, "-m BC1: Encode/decode for AMD GPU's\n");
	fprintf(stderr, "-r BC1: Encode/decode using ideal BC1 formulas with rounding for 4-color block colors 2,3 (same as AMD Compressonator)\n");
	fprintf(stderr, "-f Force writing DX10-style DDS files (otherwise for BC1-5 it uses DX9-style DDS files)\n");
	fprintf(stderr, "\nBy default, this tool encodes to BC1 *without rounding* 4-color block colors 2,3, which may not match the output of some software decoders.\n");
	fprintf(stderr, "\nFor BC4 and BC5: Not all tools support reading DX9-style BC4/BC5 format files (or BC4/5 files at all). AMD Compressonator does.\n");
	fprintf(stderr, "\nFor BC1, the engine/shader must ignore decoded texture alpha because the encoder utilizes transparent texel to get black/dark texels. Use -b to disable.\n");
	fprintf(stderr, "\nReduced entropy/RDO encoding examples:\n");
	fprintf(stderr, "\n\"bc7enc -e blah.png\" - Reduced entropy BC7 encoding (fast, but only 5-10%% gains)\n");
	fprintf(stderr, "\"bc7enc -z1.0 -zc256 blah.png\" - RDO BC7 with lambda 1.0, window size 256 bytes (default window is only 128)\n");
	fprintf(stderr, "\"bc7enc -z1.0 -e -zc1024 blah.png\" - RDO BC7 with lambda 1.0, window size 1024 bytes for more gains (but slower), combined with reduced entropy BC7\n");
	fprintf(stderr, "\"bc7enc -1 -z1.0 blah.png\" - RDO BC1 with lambda 1.0\n");
			
	return EXIT_FAILURE;
}

static std::vector<float> compute_block_mse_scales(const image_u8& source_image, uint32_t blocks_x, uint32_t blocks_y, uint32_t total_blocks, bool rdo_debug_output)
{
	const float ULTRASMOOTH_BLOCK_STD_DEV_THRESHOLD = 2.9f;
	const float DARK_THRESHOLD = 13.0f;
	const float BRIGHT_THRESHOLD = 222.0f;
	const float ULTRAMOOTH_BLOCK_MSE_SCALE = 120.0f;
	const uint32_t ULTRASMOOTH_REGION_TOO_SMALL_THRESHOLD = 64;

	image_u8 ultrasmooth_blocks_vis(blocks_x, blocks_y);

	for (uint32_t by = 0; by < blocks_y; by++)
	{
		for (uint32_t bx = 0; bx < blocks_x; bx++)
		{
			color_quad_u8 block_pixels[16];
			source_image.get_block(bx, by, 4, 4, block_pixels);

			tracked_stat y_stats;
			for (uint32_t y = 0; y < 4; y++)
				for (uint32_t x = 0; x < 4; x++)
				{
					int l = block_pixels[x + y * 4].get_luma();
					y_stats.update(l);
				}

			float max_std_dev = compute_block_max_std_dev((color_quad_u8*)block_pixels, 4, 4, 3);

			float yl = max_std_dev / ULTRASMOOTH_BLOCK_STD_DEV_THRESHOLD;

			yl = clamp(yl, 0.0f, 1.0f);
			yl *= yl;

			float y_avg = y_stats.get_mean();

			if ((y_avg < DARK_THRESHOLD) || (y_avg >= BRIGHT_THRESHOLD))
				yl = 1.0f;

			int k = std::min<int>((int)(yl * 255.0f + .5f), 255);

			ultrasmooth_blocks_vis.set_rect_clipped(bx, by, 1, 1, color_quad_u8((uint8_t)k, 255));
		}
	}

	for (int pass = 0; pass < 1; pass++)
	{
		image_u8 next_vis(ultrasmooth_blocks_vis);

		for (int y = 0; y < (int)blocks_y; y++)
		{
			for (int x = 0; x < (int)blocks_x; x++)
			{
				int m = 0;

				for (int dy = -1; dy <= 1; dy++)
					for (int dx = -1; dx <= 1; dx++)
					{
						if (ultrasmooth_blocks_vis.get_clamped(x + dx, y + dy).r() == 255)
							m = std::max<int>(m, ultrasmooth_blocks_vis.get_clamped(x + dx, y + dy).r());
					}

				next_vis(x, y).set((uint8_t)m, 255);
			}
		}

		ultrasmooth_blocks_vis.swap(next_vis);
	}

	for (uint32_t pass = 0; pass < 32; pass++)
	{
		image_u8 next_vis(ultrasmooth_blocks_vis);
		for (int y = 0; y < (int)blocks_y; y++)
		{
			for (int x = 0; x < (int)blocks_x; x++)
			{
				if (ultrasmooth_blocks_vis.get_clamped(x, y).r() < 255)
				{
					int m = 0;

					for (int dy = -1; dy <= 1; dy++)
						for (int dx = -1; dx <= 1; dx++)
							if (ultrasmooth_blocks_vis.get_clamped(x + dx, y + dy).r() == 255)
								m++;

					if (m >= 5)
						next_vis.set_clipped(x, y, color_quad_u8(255, 255, 255, 255));
				}
			}
		}
		ultrasmooth_blocks_vis.swap(next_vis);
	}

	image_u8 orig_ultrasmooth_blocks_vis(ultrasmooth_blocks_vis);

	if (rdo_debug_output)
	{
		save_png("ultrasmooth_block_mask_pre_filter.png", ultrasmooth_blocks_vis, false);
	}

	for (uint32_t by = 0; by < blocks_y; by++)
	{
		for (uint32_t bx = 0; bx < blocks_x; bx++)
		{
			const bool is_ultrasmooth = ultrasmooth_blocks_vis(bx, by).r() == 0;
			if (!is_ultrasmooth)
				continue;

			std::vector<image_u8::pixel_coord> filled_pixels;
			filled_pixels.reserve(256);

			uint32_t total_set_pixels = ultrasmooth_blocks_vis.flood_fill(bx, by, color_quad_u8(255, 255, 255, 255), color_quad_u8(0, 0, 0, 255), &filled_pixels);
						
			if (total_set_pixels < ULTRASMOOTH_REGION_TOO_SMALL_THRESHOLD)
			{
				for (uint32_t i = 0; i < filled_pixels.size(); i++)
					orig_ultrasmooth_blocks_vis(filled_pixels[i].m_x, filled_pixels[i].m_y) = color_quad_u8(255, 255, 255, 255);
			}
		
		} // bx
	} // by

	ultrasmooth_blocks_vis = orig_ultrasmooth_blocks_vis;
	
	if (rdo_debug_output)
	{
		save_png("ultrasmooth_block_mask.png", ultrasmooth_blocks_vis, false);
	}
		
	std::vector<float> block_mse_scales(total_blocks);
	
	uint32_t total_ultrasmooth_blocks = 0;
	for (uint32_t by = 0; by < blocks_y; by++)
	{
		for (uint32_t bx = 0; bx < blocks_x; bx++)
		{
			const bool is_ultrasmooth = ultrasmooth_blocks_vis(bx, by).r() == 0;

			block_mse_scales[bx + by * blocks_x] = is_ultrasmooth ? ULTRAMOOTH_BLOCK_MSE_SCALE : -1.0f;
			
			total_ultrasmooth_blocks += is_ultrasmooth;
		}
	}

	printf("Total ultrasmooth blocks: %3.2f%%\n", total_ultrasmooth_blocks * 100.0f / total_blocks);

	return block_mse_scales;
}

int main(int argc, char* argv[])
{
	int max_threads = 1;
#if _OPENMP
	max_threads = std::min(std::max(1, omp_get_max_threads()), 128);
#endif

	printf("bc7enc v%s - RDO BC1-7 Texture Compressor\n", BC7ENC_VERSION);

	if (argc < 2)
		return print_usage();
		
	printf("Max threads: %u\n", max_threads);

	std::string src_filename, src_alpha_filename, dds_output_filename, png_output_filename, png_alpha_output_filename;

	bool no_output_png = false;
	bool out_cur_dir = true;

	int bc7_uber_level = BC7ENC_MAX_UBER_LEVEL;
	int max_partitions_to_scan = BC7ENC_MAX_PARTITIONS;
	bool perceptual = false;
	bool y_flip = false;
	uint32_t bc45_channel0 = 0;
	uint32_t bc45_channel1 = 1;
	
	rgbcx::bc1_approx_mode bc1_mode = rgbcx::bc1_approx_mode::cBC1Ideal;
	bool use_bc1_3color_mode = true;
	
	// We're just turning this on by default now, like NVDXT.EXE used to do back in the old original Xbox days.
	bool use_bc1_3color_mode_for_black = true; // false; 

	int bc1_quality_level = rgbcx::MAX_LEVEL;

	DXGI_FORMAT dxgi_format = DXGI_FORMAT_BC7_UNORM;
	uint32_t pixel_format_bpp = 8;
	bool force_dx10_dds = false;
		
	float rdo_lambda = 0.0f;
	bool rdo_debug_output = false;
	float rdo_smooth_block_error_scale = 15.0f;
	bool custom_rdo_smooth_block_error_scale = false;
	uint32_t m_lookback_window_size = 128;
	bool custom_lookback_window_size = false;
	bool rdo_bc7_quant_mode6_endpoints = true;
	bool rdo_bc7_weight_modes = true;
	bool rdo_bc7_weight_low_frequency_partitions = true;
	bool rdo_bc7_pbit1_weighting = true;
	float rdo_max_smooth_block_std_dev = 18.0f;
	bool rdo_allow_relative_movement = false;
	bool rdo_try_2_matches = false;
	bool rdo_ultrasmooth_block_handling = true;
	
	bool use_hq_bc345 = true;
	int bc345_search_rad = 5;
	uint32_t bc345_mode_mask = rgbcx::BC4_USE_ALL_MODES;

	bool bc7_mode6_only = false;
	bool rdo_multithreading = true;

	bool bc7_reduce_entropy = false;

	bool use_bc7e = false;
			
	FILE* pCSV_file = nullptr;
			
	for (int i = 1; i < argc; i++)
	{
		const char *pArg = argv[i];
		if (pArg[0] == '-')
		{
			switch (pArg[1])
			{
				case 'U':
				{
					use_bc7e = true;
					break;
				}
				case 'e':
				{
					bc7_reduce_entropy = true;
					break;
				}
				case 'h':
				{
					if (strcmp(pArg, "-hl") == 0)
						use_hq_bc345 = false;
					else if (strcmp(pArg, "-h6") == 0)
						bc345_mode_mask = rgbcx::BC4_USE_MODE6_FLAG;
					else if (strcmp(pArg, "-h8") == 0)
						bc345_mode_mask = rgbcx::BC4_USE_MODE8_FLAG;
					else if (strncmp(pArg, "-hr", 3) == 0)
					{
						bc345_search_rad = atoi(pArg + 3);
						bc345_search_rad = std::max(0, std::min(32, bc345_search_rad));
					}
											
					break;
				}
				case '6':
				{
					bc7_mode6_only = true;
					break;
				}
				case '1':
				{
					dxgi_format = DXGI_FORMAT_BC1_UNORM;
					pixel_format_bpp = 4;
					printf("Compressing to BC1\n");
					break;
				}
				case '3':
				{
					dxgi_format = DXGI_FORMAT_BC3_UNORM;
					pixel_format_bpp = 8;
					printf("Compressing to BC3\n");
					break;
				}
				case '4':
				{
					dxgi_format = DXGI_FORMAT_BC4_UNORM;
					pixel_format_bpp = 4;
					printf("Compressing to BC4\n");
					break;
				}
				case '5':
				{
					dxgi_format = DXGI_FORMAT_BC5_UNORM;
					pixel_format_bpp = 8;
					printf("Compressing to BC5\n");
					break;
				}
				case 'y':
				{
					y_flip = true;
					break;
				}
				case 'a':
				{
					src_alpha_filename = pArg + 2;
					break;
				}
				case 'X':
				{
					bc45_channel0 = atoi(pArg + 2);
					if ((bc45_channel0 < 0) || (bc45_channel0 > 3))
					{
						fprintf(stderr, "Invalid argument: %s\n", pArg);
						return EXIT_FAILURE;
					}
					break;
				}
				case 'Y':
				{
					bc45_channel1 = atoi(pArg + 2);
					if ((bc45_channel1 < 0) || (bc45_channel1 > 3))
					{
						fprintf(stderr, "Invalid argument: %s\n", pArg);
						return EXIT_FAILURE;
					}
					break;
				}
				case 'f':
				{
					force_dx10_dds = true;
					break;
				}
				case 'u':
				{
					bc7_uber_level = atoi(pArg + 2);
					if ((bc7_uber_level < 0) || (bc7_uber_level > 6)) //BC7ENC_MAX_UBER_LEVEL))
					{
						fprintf(stderr, "Invalid argument: %s\n", pArg);
						return EXIT_FAILURE;
					}
					break;

				}
				case 'L':
				{
					bc1_quality_level = atoi(pArg + 2);
					if (((int)bc1_quality_level < (int)rgbcx::MIN_LEVEL) || ((int)bc1_quality_level > (int)(rgbcx::MAX_LEVEL + 1)))
					{
						fprintf(stderr, "Invalid argument: %s\n", pArg);
						return EXIT_FAILURE;
					}
					break;

				}
				case 'g':
				{
					no_output_png = true;
					break;
				}
				case 's':
				{
					perceptual = true;
					break;
				}
				case 'p':
				{
					max_partitions_to_scan = atoi(pArg + 2);
					if ((max_partitions_to_scan < 0) || (max_partitions_to_scan > BC7ENC_MAX_PARTITIONS))
					{
						fprintf(stderr, "Invalid argument: %s\n", pArg);
						return EXIT_FAILURE;
					}
					break;
				}
				case 'n':
				{
					bc1_mode = rgbcx::bc1_approx_mode::cBC1NVidia;
					break;
				}
				case 'm':
				{
					bc1_mode = rgbcx::bc1_approx_mode::cBC1AMD;
					break;
				}
				case 'r':
				{
					bc1_mode = rgbcx::bc1_approx_mode::cBC1IdealRound4;
					break;
				}
				case 'z':
				{
					if (strncmp(pArg, "-zt", 3) == 0)
					{
						rdo_multithreading = false;
					}
					else if (strncmp(pArg, "-zd", 3) == 0)
					{
						rdo_debug_output = true;
					}
					else if (strncmp(pArg, "-zq", 3) == 0)
					{
						rdo_bc7_quant_mode6_endpoints = false;
					}
					else if (strncmp(pArg, "-zw", 3) == 0)
					{
						rdo_bc7_weight_modes = false;
					}
					else if (strncmp(pArg, "-zp", 3) == 0)
					{
						rdo_bc7_weight_low_frequency_partitions = false;
					}
					else if (strncmp(pArg, "-zo", 3) == 0)
					{
						rdo_bc7_pbit1_weighting = false;
					}
					else if (strncmp(pArg, "-zm", 3) == 0)
					{
						rdo_allow_relative_movement = true;
					}
					else if (strncmp(pArg, "-zn", 3) == 0)
					{
						rdo_try_2_matches = true;
					}
					else if (strncmp(pArg, "-zu", 3) == 0)
					{
						rdo_ultrasmooth_block_handling = false;
					}
					else if (strncmp(pArg, "-zb", 3) == 0)
					{
						rdo_smooth_block_error_scale = (float)atof(pArg + 3);
						rdo_smooth_block_error_scale = std::min<float>(std::max<float>(rdo_smooth_block_error_scale, 1.0f), 500.0f);
						custom_rdo_smooth_block_error_scale = true;
					}
					else if (strncmp(pArg, "-zc", 3) == 0)
					{
						m_lookback_window_size = atoi(pArg + 3);
						m_lookback_window_size = std::min<int>(std::max<int>(m_lookback_window_size, 8), 65536*2);
						custom_lookback_window_size = true;
					}
					else if (strncmp(pArg, "-zv", 3) == 0)
					{
						rdo_max_smooth_block_std_dev = (float)atof(pArg + 3);
						rdo_max_smooth_block_std_dev = std::min<float>(std::max<float>(rdo_max_smooth_block_std_dev, .000125f), 256.0f);
					}
					else
					{
						rdo_lambda = (float)atof(pArg + 2);
						rdo_lambda = std::min<float>(std::max<float>(rdo_lambda, 0.0f), 500.0f);
					}
					break;
				}
				case 'o':
				{
					out_cur_dir = false;
					break;
				}
				case 'b':
				{
					use_bc1_3color_mode_for_black = false;
					break;
				}
				case 'c':
				{
					use_bc1_3color_mode = false;
					break;
				}
				case 'v':
				{
					if (pCSV_file)
						fclose(pCSV_file);
					
					pCSV_file = fopen(pArg + 2, "a");
					if (!pCSV_file)
					{
						fprintf(stderr, "Failed opening file %s\n", pArg + 2);
						return EXIT_FAILURE;
					}
					break;
				}
				default:
				{
					fprintf(stderr, "Invalid argument: %s\n", pArg);
					return EXIT_FAILURE;
				}
			}
		}
		else
		{
			if (!src_filename.size())
				src_filename = pArg;
			else if (!dds_output_filename.size())
				dds_output_filename = pArg;
			else if (!png_output_filename.size())
				png_output_filename = pArg;
			else
			{
				fprintf(stderr, "Invalid argument: %s\n", pArg);
				return EXIT_FAILURE;
			}
		}
	}

	const uint32_t bytes_per_block = (16 * pixel_format_bpp) / 8;
	assert(bytes_per_block == 8 || bytes_per_block == 16);

	if (!src_filename.size())
	{
		fprintf(stderr, "No source filename specified!\n");
		return EXIT_FAILURE;
	}

	if (!dds_output_filename.size())
	{
		dds_output_filename = src_filename;
		strip_extension(dds_output_filename);
		if (out_cur_dir)
			strip_path(dds_output_filename);
		dds_output_filename += ".dds";
	}

	if (!png_output_filename.size())
	{
		png_output_filename = src_filename;
		strip_extension(png_output_filename);
		if (out_cur_dir)
			strip_path(png_output_filename);
		png_output_filename += "_unpacked.png";
	}

	png_alpha_output_filename = png_output_filename;
	strip_extension(png_alpha_output_filename);
	png_alpha_output_filename += "_alpha.png";
		
	image_u8 source_image;
	if (!load_png(src_filename.c_str(), source_image))
		return EXIT_FAILURE;

	printf("Source image: %s %ux%u\n", src_filename.c_str(), source_image.width(), source_image.height());

	if (src_alpha_filename.size())
	{
		image_u8 source_alpha_image;
		if (!load_png(src_alpha_filename.c_str(), source_alpha_image))
			return EXIT_FAILURE;

		printf("Source alpha image: %s %ux%u\n", src_alpha_filename.c_str(), source_alpha_image.width(), source_alpha_image.height());

		const uint32_t w = std::min(source_alpha_image.width(), source_image.width());
		const uint32_t h = std::min(source_alpha_image.height(), source_image.height());
		
		for (uint32_t y = 0; y < h; y++)
			for (uint32_t x = 0; x < w; x++)
				source_image(x, y)[3] = source_alpha_image(x, y)[1];
	}

#if 0
	// HACK HACK
	for (uint32_t y = 0; y < source_image.height(); y++)
		for (uint32_t x = 0; x < source_image.width(); x++)
			source_image(x, y)[3] = 254;
#endif

	const uint32_t orig_width = source_image.width();
	const uint32_t orig_height = source_image.height();

	if (y_flip)
	{
		image_u8 temp;
		temp.init(orig_width, orig_height);

		for (uint32_t y = 0; y < orig_height; y++)
			for (uint32_t x = 0; x < orig_width; x++)
				temp(x, (orig_height - 1) - y) = source_image(x, y);

		temp.swap(source_image);
	}

	source_image.crop_dup_borders((source_image.width() + 3) & ~3, (source_image.height() + 3) & ~3);

	const uint32_t blocks_x = source_image.width() / 4;
	const uint32_t blocks_y = source_image.height() / 4;
	const uint32_t total_blocks = blocks_x * blocks_y;
	const uint32_t total_texels = total_blocks * 16;

	bool has_alpha = false;
	for (int by = 0; by < ((int)blocks_y) && !has_alpha; by++)
	{
		for (uint32_t bx = 0; bx < blocks_x; bx++)
		{
			color_quad_u8 pixels[16];
			source_image.get_block(bx, by, 4, 4, pixels);
			for (uint32_t i = 0; i < 16; i++)
			{
				if (pixels[i].m_c[3] < 255)
				{
					has_alpha = true;
					break;
				}
			}
		}
	}
		
	if (has_alpha)
		printf("Source image has an alpha channel.\n");
	else
		printf("Source image is opaque.\n");

	block16_vec packed_image16(total_blocks);
	block8_vec packed_image8(total_blocks);

	bc7enc_compress_block_params bc7_pack_params;
	bc7enc_compress_block_params_init(&bc7_pack_params);
	if (!perceptual)
		bc7enc_compress_block_params_init_linear_weights(&bc7_pack_params);
	bc7_pack_params.m_max_partitions = max_partitions_to_scan;
	bc7_pack_params.m_uber_level = std::min(BC7ENC_MAX_UBER_LEVEL, bc7_uber_level);
				
	if (bc7_mode6_only)
		bc7_pack_params.m_mode_mask = 1 << 6;

	if ((dxgi_format == DXGI_FORMAT_BC7_UNORM) && (rdo_lambda > 0.0f))
	{
		// Slam off perceptual in RDO mode - we don't support it (too slow).
		perceptual = false;
		bc7_pack_params.m_perceptual = false;
		bc7enc_compress_block_params_init_linear_weights(&bc7_pack_params);
	}
		
	if ((dxgi_format == DXGI_FORMAT_BC7_UNORM) && (bc7_reduce_entropy))
	{
		// Configure the BC7 encoder with some decent parameters for later RDO post-processing.
		// Textures with alpha are harder for BC7 to handle, so we use more conservative defaults.
				
		bc7_pack_params.m_mode17_partition_estimation_filterbank = false;
								
		if (rdo_bc7_weight_modes)
		{
			// Weight modes 5 and especially 6 more highly than the other modes.
			if (has_alpha)
			{
				bc7_pack_params.m_mode5_error_weight = .7f;
				bc7_pack_params.m_mode6_error_weight = .6f;
			}
			else
			{
				bc7_pack_params.m_mode6_error_weight = .4f;
			}
		}

		if (rdo_bc7_weight_low_frequency_partitions)
		{
			// Slightly prefer the lower frequency partition patterns.
			bc7_pack_params.m_low_frequency_partition_weight = .9999f;
		}

		if (rdo_bc7_quant_mode6_endpoints)
		{
			// As a good default, don't quantize mode 6 endpoints if the texture has alpha. This isn't required, but helps mask textures.
			//if (!has_alpha)
				bc7_pack_params.m_quant_mode6_endpoints = true;
		}

		if (rdo_bc7_pbit1_weighting)
		{
			// Favor p-bit 0 vs. 1, to slightly lower the entropy of output blocks with p-bits
			bc7_pack_params.m_pbit1_weight = 1.3f;
		}
	} 

#if SUPPORT_BC7E
	ispc::bc7e_compress_block_init();

	// Now initialize the BC7 compressor's parameters.
	ispc::bc7e_compress_block_params pack_params;
	memset(&pack_params, 0, sizeof(pack_params));
	switch (bc7_uber_level)
	{
	case 0:
		ispc::bc7e_compress_block_params_init_ultrafast(&pack_params, perceptual);
		break;
	case 1:
		ispc::bc7e_compress_block_params_init_veryfast(&pack_params, perceptual);
		break;
	case 2:
		ispc::bc7e_compress_block_params_init_fast(&pack_params, perceptual);
		break;
	case 3:
		ispc::bc7e_compress_block_params_init_basic(&pack_params, perceptual);
		break;
	case 4:
		ispc::bc7e_compress_block_params_init_slow(&pack_params, perceptual);
		break;
	case 5:
		ispc::bc7e_compress_block_params_init_veryslow(&pack_params, perceptual);
		break;
	case 6:
	default:
		ispc::bc7e_compress_block_params_init_slowest(&pack_params, perceptual);
		break;
	}
#endif

	if (dxgi_format == DXGI_FORMAT_BC7_UNORM)
	{
		if ((SUPPORT_BC7E) && (use_bc7e))
			printf("BC7E uber level: %u, perceptual: %u\n", bc7_uber_level, perceptual);
		else
		{
			printf("\nbc7enc parameters:\n");
			bc7_pack_params.print();
		}
	}
	else
	{
		printf("BC1 level: %u, use 3-color mode: %u, use 3-color mode for black: %u, bc1_mode: %u\nrdo_q: %f, lookback_window_size: %u, rdo_smooth_block_error_scale: %f\n", 
			bc1_quality_level, use_bc1_3color_mode, use_bc1_3color_mode_for_black, (int)bc1_mode, rdo_lambda, m_lookback_window_size, rdo_smooth_block_error_scale);
	}

	if ((dxgi_format == DXGI_FORMAT_BC3_UNORM) || (dxgi_format == DXGI_FORMAT_BC4_UNORM) || (dxgi_format == DXGI_FORMAT_BC5_UNORM))
		printf("Use high quality BC4 block encoder: %u, BC4 block radius: %u, use 6 value mode: %u, use 8 value mode: %u\n",
			use_hq_bc345, bc345_search_rad, (bc345_mode_mask & 2) != 0, (bc345_mode_mask & 1) != 0);

	// Compress all the blocks to BC1-7
	bc7enc_compress_block_init();
	rgbcx::init(bc1_mode);

	clock_t start_t = clock();
	uint32_t bc7_mode_hist[8];
	memset(bc7_mode_hist, 0, sizeof(bc7_mode_hist));

#if SUPPORT_BC7E
	if ((dxgi_format == DXGI_FORMAT_BC7_UNORM) && (use_bc7e))
	{
		printf("Using bc7e: ");

#pragma omp parallel for
		for (int32_t by = 0; by < static_cast<int32_t>(blocks_y); by++)
		{
			// Process 64 blocks at a time, for efficient SIMD processing.
			// Ideally, N >= 8 (or more) and (N % 8) == 0.
			const int N = 64;

			for (uint32_t bx = 0; bx < blocks_x; bx += N)
			{
				const uint32_t num_blocks_to_process = std::min<uint32_t>(blocks_x - bx, N);

				color_quad_u8 pixels[16 * N];

				// Extract num_blocks_to_process 4x4 pixel blocks from the source image and put them into the pixels[] array.
				for (uint32_t b = 0; b < num_blocks_to_process; b++)
					source_image.get_block(bx + b, by, 4, 4, pixels + b * 16);

				// Compress the blocks to BC7.
				// Note: If you've used Intel's ispc_texcomp, the input pixels are different. BC7E requires a pointer to an array of 16 pixels for each block.
				block16* pBlock = &packed_image16[bx + by * blocks_x];
				ispc::bc7e_compress_blocks(num_blocks_to_process, reinterpret_cast<uint64_t*>(pBlock), reinterpret_cast<const uint32_t*>(pixels), &pack_params);
			}

			if ((by & 63) == 0)
				printf(".");
		}

		for (int by = 0; by < (int)blocks_y; by++)
		{
			for (uint32_t bx = 0; bx < blocks_x; bx++)
			{
				block16* pBlock = &packed_image16[bx + by * blocks_x];

				uint32_t mode = ((uint8_t*)pBlock)[0];
				for (uint32_t m = 0; m <= 7; m++)
				{
					if (mode & (1 << m))
					{
						bc7_mode_hist[m]++;
						break;
					}
				}
			}
		}
	}
	else
#endif
	{
#pragma omp parallel for
		for (int by = 0; by < (int)blocks_y; by++)
		{
			for (uint32_t bx = 0; bx < blocks_x; bx++)
			{
				color_quad_u8 pixels[16];

				source_image.get_block(bx, by, 4, 4, pixels);

				switch (dxgi_format)
				{
				case DXGI_FORMAT_BC1_UNORM:
				{
					block8* pBlock = &packed_image8[bx + by * blocks_x];

					rgbcx::encode_bc1(bc1_quality_level, pBlock, &pixels[0].m_c[0], use_bc1_3color_mode, use_bc1_3color_mode_for_black);
					break;
				}
				case DXGI_FORMAT_BC3_UNORM:
				{
					block16* pBlock = &packed_image16[bx + by * blocks_x];

					if (use_hq_bc345)
						rgbcx::encode_bc3_hq(bc1_quality_level, pBlock, &pixels[0].m_c[0], bc345_search_rad, bc345_mode_mask);
					else
						rgbcx::encode_bc3(bc1_quality_level, pBlock, &pixels[0].m_c[0]);
					break;
				}
				case DXGI_FORMAT_BC4_UNORM:
				{
					block8* pBlock = &packed_image8[bx + by * blocks_x];

					if (use_hq_bc345)
						rgbcx::encode_bc4_hq(pBlock, &pixels[0].m_c[bc45_channel0], 4, bc345_search_rad, bc345_mode_mask);
					else
						rgbcx::encode_bc4(pBlock, &pixels[0].m_c[bc45_channel0], 4);
					break;
				}
				case DXGI_FORMAT_BC5_UNORM:
				{
					block16* pBlock = &packed_image16[bx + by * blocks_x];

					if (use_hq_bc345)
						rgbcx::encode_bc5_hq(pBlock, &pixels[0].m_c[0], bc45_channel0, bc45_channel1, 4, bc345_search_rad, bc345_mode_mask);
					else
						rgbcx::encode_bc5(pBlock, &pixels[0].m_c[0], bc45_channel0, bc45_channel1, 4);
					break;
				}
				case DXGI_FORMAT_BC7_UNORM:
				{
					block16* pBlock = &packed_image16[bx + by * blocks_x];

					bc7enc_compress_block(pBlock, pixels, &bc7_pack_params);

#pragma omp critical
					{
						uint32_t mode = ((uint8_t*)pBlock)[0];
						for (uint32_t m = 0; m <= 7; m++)
						{
							if (mode & (1 << m))
							{
								bc7_mode_hist[m]++;
								break;
							}
						}
					}

					break;
				}
				default:
				{
					assert(0);
					break;
				}
				}
			}

			if ((by & 127) == 0)
				printf(".");
		}
	}

	clock_t end_t = clock();
	
	printf("\nTotal encoding time: %f secs\n", (double)(end_t - start_t) / CLOCKS_PER_SEC);
		
	if (dxgi_format == DXGI_FORMAT_BC7_UNORM)
	{
		printf("BC7 mode histogram:\n");
		for (uint32_t i = 0; i < 8; i++)
			printf("%u: %u\n", i, bc7_mode_hist[i]);
	}

	// Compress the output data losslessly using Deflate
	const void* pOutput_data = (bytes_per_block == 16) ? (void*)&packed_image16[0] : (void*)&packed_image8[0];
	const uint32_t output_data_size = total_blocks * bytes_per_block;

	size_t pre_rdo_comp_size = 0;
	void* pPre_RDO_Comp_data = tdefl_compress_mem_to_heap(pOutput_data, output_data_size, &pre_rdo_comp_size, TDEFL_MAX_PROBES_MASK);// TDEFL_DEFAULT_MAX_PROBES);
	mz_free(pPre_RDO_Comp_data);

	float pre_rdo_lz_bits_per_texel = pre_rdo_comp_size * 8.0f / total_texels;

	printf("Pre-RDO output data size: %u, LZ (Deflate) compressed file size: %u, %3.2f bits/texel\n",
		output_data_size,
		(uint32_t)pre_rdo_comp_size,
		pre_rdo_lz_bits_per_texel);
	
	//save_dds("before_rdo.dds", orig_width, orig_height, (bytes_per_block == 16) ? (void*)&packed_image16[0] : (void*)&packed_image8[0], pixel_format_bpp, dxgi_format, perceptual, force_dx10_dds);
	
	// Post-process the data with Rate Distortion Optimization
	if (rdo_lambda > 0.0f)
	{
		const uint32_t MIN_RDO_MULTITHREADING_BLOCKS = 4096;
		const int rdo_total_threads = (rdo_multithreading && (max_threads > 1) && (total_blocks >= MIN_RDO_MULTITHREADING_BLOCKS)) ? max_threads : 1;
		
		printf("rdo_total_threads: %u\n", rdo_total_threads);
		
		int blocks_remaining = total_blocks, cur_block_index = 0;
		std::vector<int> blocks_to_do(rdo_total_threads), first_block_index(rdo_total_threads);
		for (int p = 0; p < rdo_total_threads; p++)
		{
			const int num_blocks = (p == (rdo_total_threads - 1)) ? blocks_remaining : (total_blocks / rdo_total_threads);

			blocks_to_do[p] = num_blocks;
			first_block_index[p] = cur_block_index;

			cur_block_index += num_blocks;
			blocks_remaining -= num_blocks;
		}

		assert(!blocks_remaining && cur_block_index == (int)total_blocks);

		ert::reduce_entropy_params ert_p;

		ert_p.m_lambda = rdo_lambda;
		ert_p.m_lookback_window_size = m_lookback_window_size;
		ert_p.m_smooth_block_max_mse_scale = rdo_smooth_block_error_scale;
		ert_p.m_max_smooth_block_std_dev = rdo_max_smooth_block_std_dev;
		ert_p.m_debug_output = rdo_debug_output;
		ert_p.m_try_two_matches = rdo_try_2_matches;
		ert_p.m_allow_relative_movement = rdo_allow_relative_movement;
		ert_p.m_skip_zero_mse_blocks = false;

		std::vector<float> block_rgb_mse_scales(compute_block_mse_scales(source_image, blocks_x, blocks_y, total_blocks, rdo_debug_output));
						
		std::vector<rgbcx::color32> block_pixels(total_blocks * 16);

		for (uint32_t by = 0; by < blocks_y; by++)
			for (uint32_t bx = 0; bx < blocks_x; bx++)
				source_image.get_block(bx, by, 4, 4, (color_quad_u8*)&block_pixels[(bx + by * blocks_x) * 16]);

		struct unpacker_funcs
		{
			rgbcx::bc1_approx_mode m_mode;
			bool m_allow_3color_mode;
			bool m_use_bc1_3color_mode_for_black;

			static bool unpack_bc1_block(const void* pBlock, ert::color_rgba* pPixels, uint32_t block_index, void* pUser_data)
			{
				(void)block_index;

				const unpacker_funcs* pState = (const unpacker_funcs*)pUser_data;
								
				bool used_3color_mode = rgbcx::unpack_bc1(pBlock, pPixels, true, pState->m_mode);

				if (used_3color_mode)
				{
					if (!pState->m_allow_3color_mode)
						return false;

					if (!pState->m_use_bc1_3color_mode_for_black)
					{
						rgbcx::bc1_block* pBC1_block = (rgbcx::bc1_block*)pBlock;

						for (uint32_t y = 0; y < 4; y++)
						{
							for (uint32_t x = 0; x < 4; x++)
							{
								if (pBC1_block->get_selector(x, y) == 3)
									return false;
							} // x
						} // y
					}
				}

				return true;
			}

			static bool unpack_bc4_block(const void* pBlock, ert::color_rgba* pPixels, uint32_t block_index, void* pUser_data)
			{
				(void)block_index;
				(void)pUser_data;
				memset(pPixels, 0, sizeof(ert::color_rgba) * 16);
				rgbcx::unpack_bc4(pBlock, (uint8_t*)pPixels, 4);
				return true;
			}

			static bool unpack_bc7_block(const void* pBlock, ert::color_rgba* pPixels, uint32_t block_index, void* pUser_data)
			{
				(void)block_index;
				(void)pUser_data;
				return bc7decomp::unpack_bc7(pBlock, (bc7decomp::color_rgba*)pPixels);
			}
		};

		unpacker_funcs block_unpackers;
		block_unpackers.m_allow_3color_mode = use_bc1_3color_mode;
		block_unpackers.m_use_bc1_3color_mode_for_black = use_bc1_3color_mode_for_black;
		block_unpackers.m_mode = bc1_mode;

		if (dxgi_format == DXGI_FORMAT_BC7_UNORM)
		{
			ert_p.m_lookback_window_size = std::max(16U, m_lookback_window_size);

			// BC7 RDO
			const uint32_t NUM_COMPONENTS = 4;

			if (!custom_rdo_smooth_block_error_scale)
			{
				// Attempt to compute a decent conservative smooth block MSE max scaling factor.
				// No single smooth block scale setting can work for all textures (unless it's ridiuclously large, killing efficiency).
				ert_p.m_smooth_block_max_mse_scale = lerp(15.0f, 50.0f, std::min(1.0f, ert_p.m_lambda / 4.0f));

				printf("Using an automatically computed smooth block error scale of %f (use -zb# to override)\n", ert_p.m_smooth_block_max_mse_scale);
			}

			for (uint32_t by = 0; by < blocks_y; by++)
				for (uint32_t bx = 0; bx < blocks_x; bx++)
				{
					float& s = block_rgb_mse_scales[bx + by * blocks_x];
					if (s > 0.0f)
						s = std::max(ert_p.m_smooth_block_max_mse_scale, s * std::min(ert_p.m_lambda, 3.0f));
				}
									
			printf("\nERT parameters:\n");
			ert_p.print();
			printf("\n");
			
			uint32_t total_modified = 0;

			clock_t rdo_start_t = clock();
												
#pragma omp parallel for
			for (int p = 0; p < rdo_total_threads; p++)
			{
				const int first_block_to_encode = first_block_index[p];
				const int num_blocks_to_encode = blocks_to_do[p];
				if (!num_blocks_to_encode)
					continue;

				uint32_t total_modified_local = 0;

				std::vector<float> local_block_rgb_mse_scales(num_blocks_to_encode);
				for (int i = 0; i < num_blocks_to_encode; i++)
					local_block_rgb_mse_scales[i] = block_rgb_mse_scales[first_block_to_encode + i];
										
				ert::reduce_entropy(&packed_image16[first_block_to_encode], num_blocks_to_encode,
					16, 16, 4, 4, NUM_COMPONENTS,
					(ert::color_rgba*)&block_pixels[16 * first_block_to_encode], ert_p, total_modified_local,
					unpacker_funcs::unpack_bc7_block, &block_unpackers, 
					rdo_ultrasmooth_block_handling ? &local_block_rgb_mse_scales : nullptr);

#pragma omp critical
				{
					total_modified += total_modified_local;
				}
			} // p
			
			clock_t rdo_end_t = clock();

			printf("Total RDO time: %f secs\n", (double)(rdo_end_t - rdo_start_t) / CLOCKS_PER_SEC);

			printf("Total blocks modified: %u %3.2f%%\n", total_modified, total_modified * 100.0f / total_blocks);
			
			memset(bc7_mode_hist, 0, sizeof(bc7_mode_hist));

			for (int by = 0; by < (int)blocks_y; by++)
			{
				for (uint32_t bx = 0; bx < blocks_x; bx++)
				{
					block16* pBlock = &packed_image16[bx + by * blocks_x];

					uint32_t mode = ((uint8_t*)pBlock)[0];
					for (uint32_t m = 0; m <= 7; m++)
					{
						if (mode & (1 << m))
						{
							bc7_mode_hist[m]++;
							break;
						}
					}
				}
			}

			printf("BC7 mode histogram:\n");
			for (uint32_t i = 0; i < 8; i++)
				printf("%u: %u\n", i, bc7_mode_hist[i]);
		}
		else if (dxgi_format == DXGI_FORMAT_BC5_UNORM)
		{
			// BC5 RDO - One BC4 block for R followed by one BC4 block for G

			ert_p.m_lookback_window_size = std::max(16U, m_lookback_window_size);
			
			std::vector<rgbcx::color32> block_pixels_r(total_blocks * 16), block_pixels_g(total_blocks * 16);

			for (uint32_t by = 0; by < blocks_y; by++)
			{
				for (uint32_t bx = 0; bx < blocks_x; bx++)
				{
					color_quad_u8 orig_block[16];
					source_image.get_block(bx, by, 4, 4, orig_block);

					color_quad_u8* pDst_block_r = (color_quad_u8*)&block_pixels_r[(bx + by * blocks_x) * 16];
					color_quad_u8* pDst_block_g = (color_quad_u8*)&block_pixels_g[(bx + by * blocks_x) * 16];

					for (uint32_t i = 0; i < 16; i++)
					{
						pDst_block_r[i].set(orig_block[i].r(), 0, 0, 0);
						pDst_block_g[i].set(orig_block[i].g(), 0, 0, 0);
					}
				}
			}

			const uint32_t NUM_COMPONENTS = 1;

			ert_p.m_color_weights[1] = 0;
			ert_p.m_color_weights[2] = 0;
			ert_p.m_color_weights[3] = 0;

			if (!custom_rdo_smooth_block_error_scale)
			{
				// Attempt to compute a decent conservative smooth block MSE max scaling factor.
				// No single smooth block scale setting can work for all textures (unless it's ridiuclously large, killing efficiency).
				ert_p.m_smooth_block_max_mse_scale = lerp(10.0f, 30.0f, std::min(1.0f, ert_p.m_lambda / 4.0f));
								
				printf("Using an automatically computed smooth block error scale of %f (use -zb# to override)\n", ert_p.m_smooth_block_max_mse_scale);
			}

			printf("\nERT parameters:\n");
			ert_p.print();
			printf("\n");

			uint32_t total_modified_r = 0, total_modified_g = 0;
						
			clock_t rdo_start_t = clock();

#pragma omp parallel for
			for (int p = 0; p < rdo_total_threads; p++)
			{
				const int first_block_to_encode = first_block_index[p];
				const int num_blocks_to_encode = blocks_to_do[p];
				if (!num_blocks_to_encode)
					continue;

				uint32_t total_modified_local_r = 0, total_modified_local_g = 0;

				ert::reduce_entropy(&packed_image16[first_block_to_encode], num_blocks_to_encode,
					2 * sizeof(rgbcx::bc4_block), sizeof(rgbcx::bc4_block), 4, 4, NUM_COMPONENTS,
					(ert::color_rgba*)&block_pixels_r[16 * first_block_to_encode], ert_p, total_modified_local_r,
					unpacker_funcs::unpack_bc4_block, &block_unpackers);

				ert::reduce_entropy((uint8_t*)&packed_image16[first_block_to_encode] + sizeof(rgbcx::bc4_block), num_blocks_to_encode,
					2 * sizeof(rgbcx::bc4_block), sizeof(rgbcx::bc4_block), 4, 4, NUM_COMPONENTS,
					(ert::color_rgba*)&block_pixels_g[16 * first_block_to_encode], ert_p, total_modified_local_g,
					unpacker_funcs::unpack_bc4_block, &block_unpackers);

#pragma omp critical
				{
					total_modified_r += total_modified_local_r;
					total_modified_g += total_modified_local_g;
				}
			} // p

			clock_t rdo_end_t = clock();

			printf("Total RDO time: %f secs\n", (double)(rdo_end_t - rdo_start_t) / CLOCKS_PER_SEC);

			printf("Total blocks modified R: %u %3.2f%%\n", total_modified_r, total_modified_r * 100.0f / total_blocks);
			printf("Total blocks modified G: %u %3.2f%%\n", total_modified_g, total_modified_g * 100.0f / total_blocks);
		}
		else if (dxgi_format == DXGI_FORMAT_BC4_UNORM) 
		{
			// BC4 RDO - One BC4 block for R

			const uint32_t NUM_COMPONENTS = 1;

			ert_p.m_color_weights[1] = 0;
			ert_p.m_color_weights[2] = 0;
			ert_p.m_color_weights[3] = 0;

			if (!custom_rdo_smooth_block_error_scale)
			{
				// Attempt to compute a decent conservative smooth block MSE max scaling factor.
				// No single smooth block scale setting can work for all textures (unless it's ridiuclously large, killing efficiency).
				ert_p.m_smooth_block_max_mse_scale = lerp(10.0f, 30.0f, std::min(1.0f, ert_p.m_lambda / 4.0f));

				printf("Using an automatically computed smooth block error scale of %f (use -zb# to override)\n", ert_p.m_smooth_block_max_mse_scale);
			}

			printf("\nERT parameters:\n");
			ert_p.print();
			printf("\n");

			uint32_t total_modified = 0;
						
			clock_t rdo_start_t = clock();
						
#pragma omp parallel for
			for (int p = 0; p < rdo_total_threads; p++)
			{
				const int first_block_to_encode = first_block_index[p];
				const int num_blocks_to_encode = blocks_to_do[p];
				if (!num_blocks_to_encode)
					continue;

				uint32_t total_modified_local = 0;

				ert::reduce_entropy(&packed_image8[first_block_to_encode], num_blocks_to_encode,
					sizeof(rgbcx::bc4_block), sizeof(rgbcx::bc4_block), 4, 4, NUM_COMPONENTS,
					(ert::color_rgba*)&block_pixels[16 * first_block_to_encode], ert_p, total_modified_local,
					unpacker_funcs::unpack_bc4_block, &block_unpackers);

#pragma omp critical
				{
					total_modified += total_modified_local;
				}
			} // p

			clock_t rdo_end_t = clock();

			printf("Total RDO time: %f secs\n", (double)(rdo_end_t - rdo_start_t) / CLOCKS_PER_SEC);

			printf("Total blocks modified: %u %3.2f%%\n", total_modified, total_modified * 100.0f / total_blocks);
		}
		else if (dxgi_format == DXGI_FORMAT_BC1_UNORM)
		{
			// BC1 RDO - One BC1 block
			const uint32_t NUM_COMPONENTS = 3;

			ert_p.m_color_weights[3] = 0;
						
			if (!custom_rdo_smooth_block_error_scale)
			{
				// This is just a hack - no single setting can work for all textures.
				ert_p.m_smooth_block_max_mse_scale = lerp(15.0f, 50.0f, std::min(1.0f, ert_p.m_lambda / 8.0f));

				printf("Using an automatically computed smooth block error scale of %f (use -zb# to override)\n", ert_p.m_smooth_block_max_mse_scale);
			}

			for (uint32_t by = 0; by < blocks_y; by++)
				for (uint32_t bx = 0; bx < blocks_x; bx++)
				{
					float& s = block_rgb_mse_scales[bx + by * blocks_x];
					if (s > 0.0f)
						s = std::max(ert_p.m_smooth_block_max_mse_scale, s * std::min(ert_p.m_lambda, 3.0f));
				}

			printf("\nERT parameters:\n");
			ert_p.print();
			printf("\n");

			uint32_t total_modified = 0;
						
			clock_t rdo_start_t = clock();
						
#pragma omp parallel for
			for (int p = 0; p < rdo_total_threads; p++)
			{
				const int first_block_to_encode = first_block_index[p];
				const int num_blocks_to_encode = blocks_to_do[p];
				if (!num_blocks_to_encode)
					continue;

				uint32_t total_modified_local = 0;

				std::vector<float> local_block_rgb_mse_scales(num_blocks_to_encode);
				for (int i = 0; i < num_blocks_to_encode; i++)
					local_block_rgb_mse_scales[i] = block_rgb_mse_scales[first_block_to_encode + i];
					
				ert::reduce_entropy(&packed_image8[first_block_to_encode], num_blocks_to_encode,
					sizeof(rgbcx::bc1_block), sizeof(rgbcx::bc1_block), 4, 4, NUM_COMPONENTS,
					(ert::color_rgba*)&block_pixels[16 * first_block_to_encode], ert_p, total_modified_local,
					unpacker_funcs::unpack_bc1_block, &block_unpackers, 
					rdo_ultrasmooth_block_handling ? &local_block_rgb_mse_scales : nullptr);
					
#pragma omp critical
				{
					total_modified += total_modified_local;
				}
			} // p

			clock_t rdo_end_t = clock();

			printf("Total RDO time: %f secs\n", (double)(rdo_end_t - rdo_start_t) / CLOCKS_PER_SEC);

			printf("Total blocks modified: %u %3.2f%%\n",
				total_modified, total_modified * 100.0f / total_blocks);
		}
		else if (dxgi_format == DXGI_FORMAT_BC3_UNORM)
		{
			// BC3 RDO - One BC4 block followed by one BC1 block
			
			ert_p.m_lookback_window_size = std::max(16U, m_lookback_window_size);

			std::vector<rgbcx::color32> block_pixels_a(total_blocks * 16);

			for (uint32_t by = 0; by < blocks_y; by++)
			{
				for (uint32_t bx = 0; bx < blocks_x; bx++)
				{
					color_quad_u8 orig_block[16];
					source_image.get_block(bx, by, 4, 4, orig_block);

					color_quad_u8* pDst_block_a = (color_quad_u8*)&block_pixels_a[(bx + by * blocks_x) * 16];
					for (uint32_t i = 0; i < 16; i++)
						pDst_block_a[i].set(orig_block[i].a(), 0, 0, 0);
				}
			}

			ert_p.m_color_weights[3] = 0;

			ert::reduce_entropy_params ert_alpha_p(ert_p);
			ert_alpha_p.m_color_weights[1] = 0;
			ert_alpha_p.m_color_weights[2] = 0;
			ert_alpha_p.m_color_weights[3] = 0;
						
			if (!custom_rdo_smooth_block_error_scale)
			{
				// This is just a hack - no single setting can work for all textures.
				ert_p.m_smooth_block_max_mse_scale = lerp(15.0f, 50.0f, std::min(1.0f, ert_p.m_lambda / 8.0f));

				printf("Using an automatically computed smooth block error scale of %f (use -zb# to override) for RGB\n", ert_p.m_smooth_block_max_mse_scale);
				
				ert_alpha_p.m_smooth_block_max_mse_scale = lerp(10.0f, 30.0f, std::min(1.0f, ert_alpha_p.m_lambda / 4.0f));
				printf("Using an automatically computed smooth block error scale of %f for Alpha\n", ert_alpha_p.m_smooth_block_max_mse_scale);
			}

			for (uint32_t by = 0; by < blocks_y; by++)
				for (uint32_t bx = 0; bx < blocks_x; bx++)
				{
					float& s = block_rgb_mse_scales[bx + by * blocks_x];
					if (s > 0.0f)
						s = std::max(ert_p.m_smooth_block_max_mse_scale, s * std::min(ert_p.m_lambda, 3.0f));
				}

			printf("\nERT RGB parameters:\n");
			ert_p.print();

			printf("\nERT Alpha parameters:\n");
			ert_alpha_p.print();
			printf("\n");

			uint32_t total_modified_rgb = 0, total_modified_alpha = 0;

			block_unpackers.m_allow_3color_mode = false;
			block_unpackers.m_use_bc1_3color_mode_for_black = false;
						
			clock_t rdo_start_t = clock();
						
#pragma omp parallel for
			for (int p = 0; p < rdo_total_threads; p++)
			{
				const int first_block_to_encode = first_block_index[p];
				const int num_blocks_to_encode = blocks_to_do[p];
				if (!num_blocks_to_encode)
					continue;

				uint32_t total_modified_local_rgb = 0, total_modified_local_alpha = 0;
					
				ert::reduce_entropy((uint8_t*)&packed_image16[first_block_to_encode], num_blocks_to_encode,
					sizeof(rgbcx::bc1_block) * 2, sizeof(rgbcx::bc4_block), 4, 4, 1,
					(ert::color_rgba*)&block_pixels_a[16 * first_block_to_encode], ert_alpha_p, total_modified_local_alpha,
					unpacker_funcs::unpack_bc4_block, &block_unpackers);

				std::vector<float> local_block_rgb_mse_scales(num_blocks_to_encode);
				for (int i = 0; i < num_blocks_to_encode; i++)
					local_block_rgb_mse_scales[i] = block_rgb_mse_scales[first_block_to_encode + i];

				ert::reduce_entropy((uint8_t *)&packed_image16[first_block_to_encode] + sizeof(rgbcx::bc1_block), num_blocks_to_encode,
					sizeof(rgbcx::bc1_block) * 2, sizeof(rgbcx::bc1_block), 4, 4, 3,
					(ert::color_rgba*)&block_pixels[16 * first_block_to_encode], ert_p, total_modified_local_rgb,
					unpacker_funcs::unpack_bc1_block, &block_unpackers, 
					rdo_ultrasmooth_block_handling ? &local_block_rgb_mse_scales : nullptr);
					
#pragma omp critical
				{
					total_modified_rgb += total_modified_local_rgb;
					total_modified_alpha += total_modified_local_alpha;
				}
			} // p
			
			clock_t rdo_end_t = clock();

			printf("Total RDO time: %f secs\n", (double)(rdo_end_t - rdo_start_t) / CLOCKS_PER_SEC);

			printf("Total RGB blocks modified: %u %3.2f%%\n", total_modified_rgb, total_modified_rgb * 100.0f / total_blocks);
			printf("Total Alpha blocks modified: %u %3.2f%%\n", total_modified_alpha, total_modified_alpha * 100.0f / total_blocks);
		}

	} // if (rdo_lambda > 0.0f)

	clock_t overall_end_t = clock();

	printf("Total processing time: %f secs\n", (double)(overall_end_t - start_t) / CLOCKS_PER_SEC);
	
	size_t comp_size = 0;
	void* pComp_data = tdefl_compress_mem_to_heap(pOutput_data, output_data_size, &comp_size, TDEFL_MAX_PROBES_MASK);// TDEFL_DEFAULT_MAX_PROBES);
	size_t decomp_size = 0;
	void* pDecomp_data = tinfl_decompress_mem_to_heap(pComp_data, comp_size, &decomp_size, 0);
	if ((decomp_size != output_data_size) || (memcmp(pDecomp_data, pOutput_data, decomp_size) != 0))
	{
		fprintf(stderr, "miniz compression or decompression failed!\n");
		return EXIT_FAILURE;
	}
	mz_free(pComp_data);
	mz_free(pDecomp_data);
		
	float lz_bits_per_texel = comp_size * 8.0f / total_texels;

	printf("Output data size: %u, LZ (Deflate) compressed file size: %u, %3.2f bits/texel\n",
		output_data_size,
		(uint32_t)comp_size,
		lz_bits_per_texel);
			
	bool failed = false;
	if (!save_dds(dds_output_filename.c_str(), orig_width, orig_height, (bytes_per_block == 16) ? (void*)&packed_image16[0] : (void*)&packed_image8[0], pixel_format_bpp, dxgi_format, perceptual, force_dx10_dds))
		failed = true;
	else
		printf("Wrote DDS file %s\n", dds_output_filename.c_str());

	float csv_psnr = 0.0f, csv_ssim = 0.0f;
	(void)csv_ssim;

	if ((!no_output_png) && (png_output_filename.size()))
	{
		image_u8 unpacked_image(source_image.width(), source_image.height());

		bool bc1_punchthrough_flag = false;
		bool used_bc1_transparent_texels_for_black = false;

#pragma omp parallel for
		for (int by = 0; by < (int)blocks_y; by++)
		{
			for (uint32_t bx = 0; bx < blocks_x; bx++)
			{
				void* pBlock = (bytes_per_block == 16) ? (void *)&packed_image16[bx + by * blocks_x] : (void*)&packed_image8[bx + by * blocks_x];

				color_quad_u8 unpacked_pixels[16];
				for (uint32_t i = 0; i < 16; i++)
					unpacked_pixels[i].set(0, 0, 0, 255);

				switch (dxgi_format)
				{
				case DXGI_FORMAT_BC1_UNORM:
				{
					const bool used_punchthrough = rgbcx::unpack_bc1(pBlock, unpacked_pixels, true, bc1_mode);

					if (used_punchthrough)
					{
						bc1_punchthrough_flag = true;

						const rgbcx::bc1_block* pBC1_block = (const rgbcx::bc1_block*)pBlock;

						for (uint32_t y = 0; y < 4; y++)
							for (uint32_t x = 0; x < 4; x++)
								if (pBC1_block->get_selector(x, y) == 3)
									used_bc1_transparent_texels_for_black = true;
					}

					break;
				}
				case DXGI_FORMAT_BC3_UNORM:
					if (!rgbcx::unpack_bc3(pBlock, unpacked_pixels, bc1_mode))
						bc1_punchthrough_flag = true;
					break;
				case DXGI_FORMAT_BC4_UNORM:
					rgbcx::unpack_bc4(pBlock, &unpacked_pixels[0][0], 4);
#if DECODE_BC4_TO_GRAYSCALE
					for (uint32_t i = 0; i < 16; i++)
					{
						unpacked_pixels[i][1] = unpacked_pixels[i][0];
						unpacked_pixels[i][2] = unpacked_pixels[i][0];
					}
#endif
					break;
				case DXGI_FORMAT_BC5_UNORM:
					rgbcx::unpack_bc5(pBlock, &unpacked_pixels[0][0], 0, 1, 4);
					break;
				case DXGI_FORMAT_BC7_UNORM:
					if (!bc7decomp::unpack_bc7((const uint8_t*)pBlock, (bc7decomp::color_rgba*)unpacked_pixels))
						printf("bc7decomp::unpack_bc7() failed!\n");

					break;
				default:
					assert(0);
					break;
				}

				unpacked_image.set_block(bx, by, 4, 4, unpacked_pixels);
			} // bx
		} // by

		// Sanity check the BC1/BC3 output
		if (dxgi_format == DXGI_FORMAT_BC3_UNORM)
		{
			if (bc1_punchthrough_flag)
				fprintf(stderr, "WARNING: BC3 mode selected, but rgbcx::unpack_bc3() returned one or more blocks using 3-color mode!\n");
		}
		else if (dxgi_format == DXGI_FORMAT_BC1_UNORM)
		{
			if ((bc1_punchthrough_flag) && (!use_bc1_3color_mode))
				fprintf(stderr, "WARNING: BC1 output used 3-color mode, when this was disabled!\n");

			if ((used_bc1_transparent_texels_for_black) && (!used_bc1_transparent_texels_for_black))
				fprintf(stderr, "WARNING: BC1 output used the transparent selector for black, when this was disabled!\n");
		}

		if ((dxgi_format == DXGI_FORMAT_BC1_UNORM) || (dxgi_format == DXGI_FORMAT_BC3_UNORM))
			printf("Output used 3-color mode: %u, output used transparent texels for black: %u\n", bc1_punchthrough_flag, used_bc1_transparent_texels_for_black);

		if ((dxgi_format != DXGI_FORMAT_BC4_UNORM) && (dxgi_format != DXGI_FORMAT_BC5_UNORM))
		{
			image_metrics y_metrics;
			y_metrics.compute(source_image, unpacked_image, 0, 0);
			printf("Luma  Max error: %3.0f RMSE: %f PSNR %03.3f dB, PSNR per bits/texel: %f\n", y_metrics.m_max, y_metrics.m_root_mean_squared, y_metrics.m_peak_snr, y_metrics.m_peak_snr / lz_bits_per_texel * 10000.0f);

			image_metrics rgb_metrics;
			rgb_metrics.compute(source_image, unpacked_image, 0, 3);
			printf("RGB   Max error: %3.0f RMSE: %f PSNR %03.3f dB, PSNR per bits/texel: %f\n", rgb_metrics.m_max, rgb_metrics.m_root_mean_squared, rgb_metrics.m_peak_snr, rgb_metrics.m_peak_snr / lz_bits_per_texel * 10000.0f);

			csv_psnr = (float)rgb_metrics.m_peak_snr;

			image_metrics rgba_metrics;
			rgba_metrics.compute(source_image, unpacked_image, 0, 4);
			printf("RGBA  Max error: %3.0f RMSE: %f PSNR %03.3f dB, PSNR per bits/texel: %f\n", rgba_metrics.m_max, rgba_metrics.m_root_mean_squared, rgba_metrics.m_peak_snr, rgba_metrics.m_peak_snr / lz_bits_per_texel * 10000.0f);

#if COMPUTE_SSIM
			vec4F ssim_y(compute_ssim(source_image, unpacked_image, true));
			vec4F ssim_rgba(compute_ssim(source_image, unpacked_image, false));

			printf("R       SSIM: %f\n", ssim_rgba[0]);
			printf("G       SSIM: %f\n", ssim_rgba[1]);
			printf("B       SSIM: %f\n", ssim_rgba[2]);
			printf("RGB Avg SSIM: %f\n", (ssim_rgba[0] + ssim_rgba[1] + ssim_rgba[2]) / 3.0f);
			printf("A       SSIM: %f\n", ssim_rgba[3]);
						
			printf("Luma    SSIM: %f\n", ssim_y[0]);

			csv_ssim = (ssim_rgba[0] + ssim_rgba[1] + ssim_rgba[2]) / 3.0f;
#endif
		}
						
		for (uint32_t chan = 0; chan < 4; chan++)
		{
			if (dxgi_format == DXGI_FORMAT_BC4_UNORM)
			{
				if (chan != bc45_channel0)
					continue;
			}
			else if (dxgi_format == DXGI_FORMAT_BC5_UNORM)
			{
				if ((chan != bc45_channel0) && (chan != bc45_channel1))
					continue;
			}

			image_metrics c_metrics;
			c_metrics.compute(source_image, unpacked_image, chan, 1);
			static const char *s_chan_names[4] = { "Red  ", "Green", "Blue ", "Alpha" };
			printf("%s Max error: %3.0f RMSE: %f PSNR %03.3f dB\n", s_chan_names[chan], c_metrics.m_max, c_metrics.m_root_mean_squared, c_metrics.m_peak_snr);

			if (dxgi_format == DXGI_FORMAT_BC4_UNORM)
				csv_psnr = (float)c_metrics.m_peak_snr;
		}

		if (dxgi_format == DXGI_FORMAT_BC5_UNORM)
		{
			image_metrics c_metrics;
			c_metrics.compute(source_image, unpacked_image, 0, 2);
			
			printf("RG Max error: %3.0f RMSE: %f PSNR %03.3f dB\n", c_metrics.m_max, c_metrics.m_root_mean_squared, c_metrics.m_peak_snr);

			csv_psnr = (float)c_metrics.m_peak_snr;
		}

		if (bc1_mode != rgbcx::bc1_approx_mode::cBC1Ideal)
			printf("Note: BC1/BC3 RGB decoding was done with the specified vendor's BC1 approximations.\n");

		image_u8 unpacked_image_cropped(unpacked_image);
		unpacked_image_cropped.crop(orig_width, orig_height);
		if (!save_png(png_output_filename.c_str(), unpacked_image_cropped, false))
			failed = true;
		else
			printf("Wrote PNG file %s\n", png_output_filename.c_str());
				
		if (png_alpha_output_filename.size())
		{
			image_u8 unpacked_image_alpha(unpacked_image);
			for (uint32_t y = 0; y < unpacked_image_alpha.height(); y++)
				for (uint32_t x = 0; x < unpacked_image_alpha.width(); x++)
					unpacked_image_alpha(x, y).set(unpacked_image_alpha(x, y)[3], 255);
			unpacked_image_alpha.crop(orig_width, orig_height);

			if (!save_png(png_alpha_output_filename.c_str(), unpacked_image_alpha, false))
				failed = true;
			else
				printf("Wrote PNG file %s\n", png_alpha_output_filename.c_str());
		}
	}

	if (pCSV_file)
	{
		fprintf(pCSV_file, "%f,%f,%f\n", rdo_lambda, lz_bits_per_texel, csv_psnr);
		fclose(pCSV_file);
		pCSV_file = nullptr;
	}

#if LZHAM_STATS
	remove("__lzham_compressed");
	char buf[1024];
	sprintf(buf, "lzhamtest_x64 -x8 -o -e -h4 c \"%s\" __lzham_compressed", dds_output_filename.c_str());
	system(buf);

	{
		FILE* p = fopen("__lzham_compressed", "rb");
		if (p)
		{
			fseek(p, 0, SEEK_END);
			uint32_t lzham_size = (uint32_t)ftell(p);
			fclose(p);

			float lzham_bits_per_texel = lzham_size * 8.0f / total_texels;

			printf("LZHAM compressed file size: %u, %3.2f bits/texel\n",
				(uint32_t)lzham_size,
				lzham_bits_per_texel);
		}
	}
#endif
		
	return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
