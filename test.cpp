// test.cpp - Command line example/test app
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#if _OPENMP
#include <omp.h>
#endif

#include "utils.h"

#include "bc7enc.h"
#include "bc7decomp.h"

#define RGBCX_IMPLEMENTATION
#include "rgbcx.h"

#include "miniz.h"

#define COMPUTE_SSIM (0)

const int MAX_UBER_LEVEL = 5;

static int print_usage()
{
	fprintf(stderr, "bc7enc\n");
	fprintf(stderr, "Reads PNG files (with or without alpha channels) and packs them to BC1-5 or BC7/BPTC (default) using\nmodes 1, 6 (opaque blocks) or modes 1, 5, 6, and 7 (alpha blocks).\n");
	fprintf(stderr, "Supports optional reduced entropy BC7 encoding (using -e) and Rate Distortion Optimization (RDO) for BC1-7 (using -z# where # is lambda).\n");
	fprintf(stderr, "By default, this tool compresses to BC7. A DX10 DDS file and a unpacked PNG file will be written to the source\ndirectory with the .dds/_unpacked.png/_unpacked_alpha.png suffixes.\n\n");
	fprintf(stderr, "Usage: bc7enc [-apng_filename] [options] input_filename.png [compressed_output.dds] [unpacked_output.png]\n\n");
	fprintf(stderr, "-apng_filename Load G channel of PNG file into alpha channel of source image\n");
	fprintf(stderr, "-g Don't write unpacked output PNG files (this disables PSNR metrics too).\n");
	fprintf(stderr, "-y Flip source image along Y axis before packing\n");
	fprintf(stderr, "-o Write output files to the current directory\n");
	fprintf(stderr, "-1 Encode to BC1. -u[0,5] controls quality vs. perf. tradeoff for RGB.\n");
	fprintf(stderr, "-3 Encode to BC3. -u[0,5] controls quality vs. perf. tradeoff for RGB.\n");
	fprintf(stderr, "-4 Encode to BC4\n");
	fprintf(stderr, "-5 Encode to BC5\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "-X# BC4/5: Set first color channel (defaults to 0 or red)\n");
	fprintf(stderr, "-Y# BC4/5: Set second color channel (defaults to 1 or green)\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "-l BC7: Use linear colorspace metrics instead of perceptual (the default is perceptual or sRGB!). BC7 RDO mode is always linear.\n");
	fprintf(stderr, "-uX BC7: Higher quality levels, X ranges from [0,4] for BC7\n");
	fprintf(stderr, "-pX BC7: Scan X partitions in mode 1, X ranges from [0,64], use 0 to disable mode 1 entirely (faster)\n");
	fprintf(stderr, "-LX BC1: Set encoding level, where 0=fastest and 18=slowest but highest quality\n");
	fprintf(stderr, "\nRDO mode options:\n");
	fprintf(stderr, "-z# BC1-7: Set RDO lambda factor (quality), lower=higher quality/larger LZ compressed files, try .1-4\n");
	fprintf(stderr, "-zb# BC1-7: Manually set smooth block scale factor, higher values = less distortion on smooth blocks, try 5-70\n");
	fprintf(stderr, "-zc# BC1: Set RDO lookback window size in bytes (higher=more effective but slower, default=256 for BC7 and 1024 for BC1-5, try 128-8192)\n");
	fprintf(stderr, "-e BC7: Quantize/weight BC7 output for lower entropy (no slowdown but only 5-10%% gains, can be combined with -z# for more gains)\n");
	fprintf(stderr, "RDO debugging/development:\n");
	fprintf(stderr, "-zd BC1-7: Enable debug output\n");
	fprintf(stderr, "-zr BC1: Disable RDO endpoint/selector refinement stages (lowers quality)\n");
	fprintf(stderr, "-zs BC1: Disable selector RDO (lowers avg quality per output bit)\n");
	fprintf(stderr, "-ze BC1: Disable endpoint RDO (lowers avg quality per output bit)\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "-b BC1: Enable 3-color mode for blocks containing black or very dark pixels. (Important: engine/shader MUST ignore decoded texture alpha if this flag is enabled!)\n");
	fprintf(stderr, "-c BC1: Disable 3-color mode for solid color blocks\n");
	fprintf(stderr, "-n BC1: Encode/decode for NVidia GPU's\n");
	fprintf(stderr, "-m BC1: Encode/decode for AMD GPU's\n");
	fprintf(stderr, "-r BC1: Encode/decode using ideal BC1 formulas with rounding for 4-color block colors 2,3 (same as AMD Compressonator)\n");
	fprintf(stderr, "-f Force writing DX10-style DDS files (otherwise for BC1-5 it uses DX9-style DDS files)\n");
	fprintf(stderr, "\nBy default, this tool encodes to BC1 *without rounding* 4-color block colors 2,3, which may not match the output of some software decoders.\n");
	fprintf(stderr, "\nFor BC4 and BC5: Not all tools support reading DX9-style BC4/BC5 format files (or BC4/5 files at all). AMD Compressonator does.\n");
	fprintf(stderr, "\nReduced entropy/RDO examples:\n");
	fprintf(stderr, "\n\"bc7enc -o -u4 -e blah.png\" - Reduced entropy BC7 encoding (fast, but only 5-10%% gains)\n");
	fprintf(stderr, "\"bc7enc -o -u4 -z1.0 -zc256 blah.png\" - RDO BC7 with lambda 1.0, window size 256 bytes\n");
	fprintf(stderr, "\"bc7enc -o -u4 -z1.0 -e -zc256 blah.png\" - RDO BC7 with lambda 1.0, window size 256 bytes, combined with reduced entropy BC7\n");
	fprintf(stderr, "\"bc7enc -o -1 -L18 -z1.0 blah.png\" - RDO BC1 with lambda 1.0\n");
			
	return EXIT_FAILURE;
}

int main(int argc, char *argv[])
{
	if (argc < 2)
		return print_usage();

	int max_threads = 1;
#if _OPENMP
	max_threads = std::max(1, omp_get_max_threads());
#endif
	printf("Max threads: %u\n", max_threads);

	std::string src_filename;
	std::string src_alpha_filename;
	std::string dds_output_filename;
	std::string png_output_filename;
	std::string png_alpha_output_filename;

	bool no_output_png = false;
	bool out_cur_dir = false;

	int uber_level = 0;
	int max_partitions_to_scan = BC7ENC_MAX_PARTITIONS;
	bool perceptual = true;
	bool y_flip = false;
	uint32_t bc45_channel0 = 0;
	uint32_t bc45_channel1 = 1;
	
	rgbcx::bc1_approx_mode bc1_mode = rgbcx::bc1_approx_mode::cBC1Ideal;
	bool use_bc1_3color_mode = true;
	bool use_bc1_3color_mode_for_black = false;
	int bc1_quality_level = 10;

	DXGI_FORMAT dxgi_format = DXGI_FORMAT_BC7_UNORM;
	uint32_t pixel_format_bpp = 8;
	bool force_dx10_dds = false;

	float rdo_q = 0.0f, rdo_q_alpha = 0.0f;
	bool rdo_refinement = true;
	bool selector_rdo = true;
	bool endpoint_rdo = true;
	bool rdo_debug_output = false;
	float rdo_smooth_block_error_scale = rgbcx::BC1_RDO_DEFAULT_SMOOTH_BLOCK_ERROR_SCALE;
	float rdo_alpha_smooth_block_error_scale = rgbcx::BC4_RDO_DEFAULT_SMOOTH_BLOCK_ERROR_SCALE;
	bool custom_rdo_smooth_block_error_scale = false;
	uint32_t m_lookback_window_size = 256;
	bool custom_lookback_window_size = false;
	bool rdo_bc7_quant_mode6_endpoints = true;
	bool rdo_bc7_weight_modes = true;
	bool rdo_bc7_weight_low_frequency_partitions = true;
	bool rdo_bc7_pbit1_weighting = true;
	float rdo_max_smooth_block_std_dev = 18.0f;
	
	bool use_hq_bc345 = false;
	int bc345_search_rad = 5;
	uint32_t bc345_mode_mask = rgbcx::BC4_USE_ALL_MODES;

	bool bc7_mode6_only = false;
	bool rdo_multithreading = true;

	bool bc7_reduce_entropy = false;
			
	FILE* pCSV_file = nullptr;
			
	for (int i = 1; i < argc; i++)
	{
		const char *pArg = argv[i];
		if (pArg[0] == '-')
		{
			switch (pArg[1])
			{
				case 'e':
				{
					bc7_reduce_entropy = true;
					break;
				}
				case 'h':
				{
					if (strcmp(pArg, "-h6") == 0)
						bc345_mode_mask = rgbcx::BC4_USE_MODE6_FLAG;
					else if (strcmp(pArg, "-h8") == 0)
						bc345_mode_mask = rgbcx::BC4_USE_MODE8_FLAG;
					else if (strncmp(pArg, "-hr", 3) == 0)
					{
						bc345_search_rad = atoi(pArg + 3);
						bc345_search_rad = std::max(0, std::min(32, bc345_search_rad));
					}
					else
						use_hq_bc345 = true;
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
					uber_level = atoi(pArg + 2);
					if ((uber_level < 0) || (uber_level > MAX_UBER_LEVEL))
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
				case 'l':
				{
					perceptual = false;
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
					else if (strncmp(pArg, "-zr", 3) == 0)
					{
						rdo_refinement = false;
					}
					else if (strncmp(pArg, "-zs", 3) == 0)
					{
						selector_rdo = false;
					}
					else if (strncmp(pArg, "-ze", 3) == 0)
					{
						endpoint_rdo = false;
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
					else if (strncmp(pArg, "-zb", 3) == 0)
					{
						rdo_smooth_block_error_scale = (float)atof(pArg + 3);
						rdo_smooth_block_error_scale = std::min<float>(std::max<float>(rdo_smooth_block_error_scale, 1.0f), 500.0f);
						custom_rdo_smooth_block_error_scale = true;
					}
					else if (strncmp(pArg, "-zba", 4) == 0)
					{
						rdo_alpha_smooth_block_error_scale = (float)atof(pArg + 4);
						rdo_alpha_smooth_block_error_scale = std::min<float>(std::max<float>(rdo_alpha_smooth_block_error_scale, 1.0f), 500.0f);
					}
					else if (strncmp(pArg, "-za", 3) == 0)
					{
						rdo_q_alpha = (float)atof(pArg + 3);
						rdo_q_alpha = std::min<float>(std::max<float>(rdo_q_alpha, 32), 65536 * 2);
					}
					else if (strncmp(pArg, "-zc", 3) == 0)
					{
						m_lookback_window_size = atoi(pArg + 3);
						m_lookback_window_size = std::min<int>(std::max<int>(m_lookback_window_size, 32), 65536*2);
						custom_lookback_window_size = true;
					}
					else if (strncmp(pArg, "-zv", 3) == 0)
					{
						rdo_max_smooth_block_std_dev = (float)atof(pArg + 3);
						rdo_max_smooth_block_std_dev = std::min<float>(std::max<float>(rdo_max_smooth_block_std_dev, .000125f), 256.0f);
					}
					else
					{
						rdo_q = (float)atof(pArg + 2);
						rdo_q = std::min<float>(std::max<float>(rdo_q, 0.0f), 100.0f);
					}
					break;
				}
				case 'o':
				{
					out_cur_dir = true;
					break;
				}
				case 'b':
				{
					use_bc1_3color_mode_for_black = true;
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
	bc7_pack_params.m_uber_level = std::min(BC7ENC_MAX_UBER_LEVEL, uber_level);
				
	if (bc7_mode6_only)
		bc7_pack_params.m_mode_mask = 1 << 6;

	if ((dxgi_format == DXGI_FORMAT_BC7_UNORM) && (rdo_q > 0.0f))
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
			
	if (dxgi_format == DXGI_FORMAT_BC7_UNORM)
	{
		printf("\nbc7enc parameters:\n");
		bc7_pack_params.print();
	}
	else
	{
		if (!custom_lookback_window_size)
		{
			// Use a better default for BC1. 
			m_lookback_window_size = 1024;
		}

		printf("BC1 level: %u, use 3-color mode: %u, use 3-color mode for black: %u, bc1_mode: %u\nrdo_q: %f, rdo_q_alpha: %f, rdo_refinement: %u, selector_rdo: %u, endpoint_rdo: %u\nm_lookback_window_size: %u, rdo_smooth_block_error_scale: %f, rdo_alpha_smooth_block_error_scale: %f\n", 
			bc1_quality_level, use_bc1_3color_mode, use_bc1_3color_mode_for_black, (int)bc1_mode, rdo_q, rdo_q_alpha, rdo_refinement, selector_rdo, endpoint_rdo, m_lookback_window_size, rdo_smooth_block_error_scale, rdo_alpha_smooth_block_error_scale);
	}

	// Compress all the blocks to BC1-7
	bc7enc_compress_block_init();
	rgbcx::init(bc1_mode);

	clock_t start_t = clock();
	uint32_t bc7_mode_hist[8];
	memset(bc7_mode_hist, 0, sizeof(bc7_mode_hist));

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
	if (rdo_q > 0.0f)
	{
		const uint32_t MIN_RDO_MULTITHREADING_BLOCKS = 4096;
		const int rdo_total_threads = rdo_multithreading ? max_threads : 1;
		
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

		if (dxgi_format == DXGI_FORMAT_BC7_UNORM)
		{
			bc7enc_rdo_params bc7_rdo_p;

			bc7_rdo_p.m_lambda = rdo_q;
			bc7_rdo_p.m_lookback_window_size = m_lookback_window_size;
			bc7_rdo_p.m_smooth_block_max_mse_scale = rdo_smooth_block_error_scale;
			bc7_rdo_p.m_max_smooth_block_std_dev = rdo_max_smooth_block_std_dev;
			bc7_rdo_p.m_debug_output = rdo_debug_output;

			if (!custom_rdo_smooth_block_error_scale)
			{
				// Attempt to compute a decent conservative smooth block MSE max scaling factor.
				// No single smooth block scale setting can work for all textures (unless it's ridiuclously large, killing efficiency).
				bc7_rdo_p.m_smooth_block_max_mse_scale = lerp(15.0f, 50.0f, std::min(1.0f, bc7_rdo_p.m_lambda / 4.0f));

				printf("Using an automatically computed smooth block error scale of %f\n", bc7_rdo_p.m_smooth_block_max_mse_scale);
			}

			std::vector<rgbcx::color32> block_pixels(total_blocks * 16);

			for (uint32_t by = 0; by < blocks_y; by++)
				for (uint32_t bx = 0; bx < blocks_x; bx++)
					source_image.get_block(bx, by, 4, 4, (color_quad_u8*)&block_pixels[(bx + by * blocks_x) * 16]);

			bc7_rdo_p.print();

			clock_t rdo_start_t = clock();
						
			uint32_t total_modified = 0;

			if ((rdo_multithreading) && (rdo_total_threads > 1) && (total_blocks >= MIN_RDO_MULTITHREADING_BLOCKS))
			{
#pragma omp parallel for
				for (int p = 0; p < rdo_total_threads; p++)
				{
					const int first_block_to_encode = first_block_index[p];
					const int blocks_to_encode = blocks_to_do[p];
					if (!blocks_to_encode)
						continue;

					uint32_t total_modified_local = 0;
					bool status = bc7enc_reduce_entropy(&packed_image16[first_block_to_encode], blocks_to_encode, (color_rgba*)&block_pixels[16 * first_block_to_encode], bc7_rdo_p, total_modified_local);

#pragma omp critical
					{
						total_modified += total_modified_local;
					}

					if (!status)
					{
						fprintf(stderr, "bc7enc_reduce_entropy() failed!\n");
						exit(EXIT_FAILURE);
					}
				} // p
			}
			else
			{
				bool status = bc7enc_reduce_entropy(&packed_image16[0], total_blocks, (color_rgba*)&block_pixels[0], bc7_rdo_p, total_modified);
				if (!status)
				{
					fprintf(stderr, "bc7enc_reduce_entropy() failed!\n");
					exit(EXIT_FAILURE);
				}
			}

			printf("Total modified blocks: %u %3.2f%%\n", total_modified, total_modified * 100.0f / total_blocks);

			clock_t rdo_end_t = clock();

			printf("Total RDO postprocess time: %f secs\n", (double)(rdo_end_t - rdo_start_t) / CLOCKS_PER_SEC);
		}
		else if (dxgi_format == DXGI_FORMAT_BC4_UNORM)
		{
			if (rdo_q_alpha > 0.0f)
			{
				rgbcx::bc4_rdo_params bc4_rdo_p;
				bc4_rdo_p.m_langrangian_multiplier = rdo_q_alpha;
				bc4_rdo_p.m_refine_endpoints = rdo_refinement;
				bc4_rdo_p.m_refine_selectors = rdo_refinement;
				bc4_rdo_p.m_debug_output = rdo_debug_output;
				bc4_rdo_p.m_lz_dict_size = m_lookback_window_size;
				bc4_rdo_p.m_smooth_block_error_scale = rdo_alpha_smooth_block_error_scale;
				bc4_rdo_p.m_block_max_std_dev_rdo_quality_scaler = rdo_max_smooth_block_std_dev;

				std::vector<rgbcx::color32> block_pixels(total_blocks * 16);

				for (uint32_t by = 0; by < blocks_y; by++)
					for (uint32_t bx = 0; bx < blocks_x; bx++)
						source_image.get_block(bx, by, 4, 4, (color_quad_u8*)&block_pixels[(bx + by * blocks_x) * 16]);

				uint32_t total_skipped = 0, total_endpoints_refined = 0, total_selectors_refined = 0, total_merged = 0;

				bc4_rdo_postprocess((rgbcx::bc4_block*)&packed_image8[0], 1, total_blocks, &block_pixels[0], bc45_channel0, bc4_rdo_p, total_skipped, total_endpoints_refined, total_selectors_refined, total_merged);

				printf("Total blocks skipped: %u %3.2f%%, Total blocks merged: %u %3.2f%%, Total block endpoint refined: %u %3.2f%%, Total block selectors refined: %u %3.2f%%\n",
					total_skipped, total_skipped * 100.0f / total_blocks,
					total_merged, total_merged * 100.0f / total_blocks,
					total_endpoints_refined, total_endpoints_refined * 100.0f / total_blocks,
					total_selectors_refined, total_selectors_refined * 100.0f / total_blocks);

				printf("\n");
			}
		}
		else if (dxgi_format == DXGI_FORMAT_BC1_UNORM)
		{
			rgbcx::bc1_rdo_params bc1_rdo_p;
			bc1_rdo_p.m_langrangian_multiplier = rdo_q;
			bc1_rdo_p.m_use_transparent_texels_for_black = use_bc1_3color_mode_for_black;
			bc1_rdo_p.m_mode = bc1_mode;
			bc1_rdo_p.m_level = bc1_quality_level;
			bc1_rdo_p.m_refine_endpoints = rdo_refinement;
			bc1_rdo_p.m_refine_selectors = rdo_refinement;
			bc1_rdo_p.m_debug_output = rdo_debug_output;
			bc1_rdo_p.m_lz_dict_size = m_lookback_window_size;
			bc1_rdo_p.m_smooth_block_error_scale = rdo_smooth_block_error_scale;
			bc1_rdo_p.m_block_max_std_dev_rdo_quality_scaler = rdo_max_smooth_block_std_dev;

			if (!custom_rdo_smooth_block_error_scale)
			{
				// This is just a hack - no single setting can work for all textures.
				bc1_rdo_p.m_smooth_block_error_scale = lerp(8.0f, 30.0f, std::min(1.0f, bc1_rdo_p.m_langrangian_multiplier / 10.0f));
				printf("Using an automatically computed smooth block error scale of %f (use -zb to override)\n", bc1_rdo_p.m_smooth_block_error_scale);
			}

			std::vector<rgbcx::color32> block_pixels(total_blocks * 16);

			for (uint32_t by = 0; by < blocks_y; by++)
				for (uint32_t bx = 0; bx < blocks_x; bx++)
					source_image.get_block(bx, by, 4, 4, (color_quad_u8*)&block_pixels[(bx + by * blocks_x) * 16]);

			if (selector_rdo)
			{
				printf("--- Begin selector RDO post-process:\n");

				uint32_t total_skipped = 0, total_endpoints_refined = 0, total_selectors_refined = 0, total_merged = 0;

				clock_t rdo_start_t = clock();

				bc1_rdo_postprocess((rgbcx::bc1_block*)&packed_image8[0], 1, total_blocks, &block_pixels[0], bc1_rdo_p, total_skipped, total_endpoints_refined, total_selectors_refined, total_merged, true);

				clock_t rdo_end_t = clock();

				printf("Total time: %f secs\n", (double)(rdo_end_t - rdo_start_t) / CLOCKS_PER_SEC);

				printf("Total blocks skipped: %u %3.2f%%, Total blocks merged: %u %3.2f%%, Total block endpoints refined: %u %3.2f%%, Total block selectors refined: %u %3.2f%%\n",
					total_skipped, total_skipped * 100.0f / total_blocks,
					total_merged, total_merged * 100.0f / total_blocks,
					total_endpoints_refined, total_endpoints_refined * 100.0f / total_blocks,
					total_selectors_refined, total_selectors_refined * 100.0f / total_blocks);
			}

			if (endpoint_rdo)
			{
				printf("--- Begin endpoint RDO post-process:\n");

				uint32_t total_skipped = 0, total_endpoints_refined = 0, total_selectors_refined = 0, total_merged = 0;

				clock_t rdo_start_t = clock();

				bc1_rdo_postprocess((rgbcx::bc1_block*)&packed_image8[0], 1, total_blocks, &block_pixels[0], bc1_rdo_p, total_skipped, total_endpoints_refined, total_selectors_refined, total_merged, false);

				clock_t rdo_end_t = clock();

				printf("Total time: %f secs\n", (double)(rdo_end_t - rdo_start_t) / CLOCKS_PER_SEC);

				printf("Total blocks skipped: %u %3.2f%%, Total blocks merged: %u %3.2f%%, Total block endpoints refined: %u %3.2f%%, Total block selectors refined: %u %3.2f%%\n",
					total_skipped, total_skipped * 100.0f / total_blocks,
					total_merged, total_merged * 100.0f / total_blocks,
					total_endpoints_refined, total_endpoints_refined * 100.0f / total_blocks,
					total_selectors_refined, total_selectors_refined * 100.0f / total_blocks);
			}

			printf("\n");
		}

	} // if (rdo_q > 0.0f)

	clock_t overall_end_t = clock();

	printf("Processing time: %f secs\n", (double)(overall_end_t - start_t) / CLOCKS_PER_SEC);
	
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

	if ((!no_output_png) && (png_output_filename.size()))
	{
		image_u8 unpacked_image(source_image.width(), source_image.height());

		bool punchthrough_flag = false;

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
					rgbcx::unpack_bc1(pBlock, unpacked_pixels, true, bc1_mode);
					break;
				case DXGI_FORMAT_BC3_UNORM:
					if (!rgbcx::unpack_bc3(pBlock, unpacked_pixels, bc1_mode))
						punchthrough_flag = true;
					break;
				case DXGI_FORMAT_BC4_UNORM:
					rgbcx::unpack_bc4(pBlock, &unpacked_pixels[0][0], 4);
					for (uint32_t i = 0; i < 16; i++)
					{
						unpacked_pixels[i][1] = unpacked_pixels[i][0];
						unpacked_pixels[i][2] = unpacked_pixels[i][0];
					}
					break;
				case DXGI_FORMAT_BC5_UNORM:
					rgbcx::unpack_bc5(pBlock, &unpacked_pixels[0][0], 0, 1, 4);
					break;
				case DXGI_FORMAT_BC7_UNORM:
					bc7decomp::unpack_bc7((const uint8_t*)pBlock, (bc7decomp::color_rgba*)unpacked_pixels);
					break;
				default:
					assert(0);
					break;
				}

				unpacked_image.set_block(bx, by, 4, 4, unpacked_pixels);
			} // bx
		} // by

		if ((punchthrough_flag) && (dxgi_format == DXGI_FORMAT_BC3_UNORM))
			fprintf(stderr, "Warning: BC3 mode selected, but rgbcx::unpack_bc3() returned one or more blocks using 3-color mode!\n");

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
		fprintf(pCSV_file, "%f,%f,%f\n", rdo_q, lz_bits_per_texel, csv_psnr);
		fclose(pCSV_file);
		pCSV_file = nullptr;
	}
		
	return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
