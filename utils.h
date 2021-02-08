// File: utils.h
#pragma once
#ifdef _MSC_VER
#pragma warning (push)
#pragma warning (disable:4127) // conditional expression is constant
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <time.h>
#include <vector>
#include <string>
#include "dds_defs.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

inline int iabs(int i) { if (i < 0) i = -i; return i; }
inline uint8_t clamp255(int32_t i) { return (uint8_t)((i & 0xFFFFFF00U) ? (~(i >> 31)) : i); }
template <typename S> inline S clamp(S value, S low, S high) { return (value < low) ? low : ((value > high) ? high : value); }
template<typename F> inline F lerp(F a, F b, F s) { return a + (b - a) * s; }

struct color_quad_u8
{
	uint8_t m_c[4];

	inline color_quad_u8(uint8_t r, uint8_t g, uint8_t b, uint8_t a)
	{
		set(r, g, b, a);
	}

	inline color_quad_u8(uint8_t y = 0, uint8_t a = 255)
	{
		set(y, a);
	}

	inline color_quad_u8& set(uint8_t y, uint8_t a = 255)
	{
		m_c[0] = y;
		m_c[1] = y;
		m_c[2] = y;
		m_c[3] = a;
		return *this;
	}

	inline color_quad_u8& set(uint8_t r, uint8_t g, uint8_t b, uint8_t a)
	{
		m_c[0] = r;
		m_c[1] = g;
		m_c[2] = b;
		m_c[3] = a;
		return *this;
	}

	inline uint8_t& operator[] (uint32_t i) { assert(i < 4);  return m_c[i]; }
	inline uint8_t operator[] (uint32_t i) const { assert(i < 4); return m_c[i]; }

	inline uint8_t r() const { return m_c[0]; }
	inline uint8_t g() const { return m_c[1]; }
	inline uint8_t b() const { return m_c[2]; }
	inline uint8_t a() const { return m_c[3]; }

	inline int get_luma() const { return (13938U * m_c[0] + 46869U * m_c[1] + 4729U * m_c[2] + 32768U) >> 16U; } // REC709 weightings
};
typedef std::vector<color_quad_u8> color_quad_u8_vec;

class image_u8
{
public:
	image_u8() :
		m_width(0), m_height(0)
	{
	}

	image_u8(uint32_t width, uint32_t height) :
		m_width(width), m_height(height)
	{
		m_pixels.resize(width * height);
	}

	inline const color_quad_u8_vec& get_pixels() const { return m_pixels; }
	inline color_quad_u8_vec& get_pixels() { return m_pixels; }

	inline uint32_t width() const { return m_width; }
	inline uint32_t height() const { return m_height; }
	inline uint32_t total_pixels() const { return m_width * m_height; }

	inline color_quad_u8& operator()(uint32_t x, uint32_t y) { assert(x < m_width&& y < m_height);  return m_pixels[x + m_width * y]; }
	inline const color_quad_u8& operator()(uint32_t x, uint32_t y) const { assert(x < m_width&& y < m_height);  return m_pixels[x + m_width * y]; }

	image_u8& clear()
	{
		m_width = m_height = 0;
		m_pixels.clear();
		return *this;
	}

	image_u8& init(uint32_t width, uint32_t height)
	{
		clear();

		m_width = width;
		m_height = height;
		m_pixels.resize(width * height);
		return *this;
	}

	image_u8& set_all(const color_quad_u8& p)
	{
		for (uint32_t i = 0; i < m_pixels.size(); i++)
			m_pixels[i] = p;
		return *this;
	}

	inline const color_quad_u8& get_clamped(int x, int y) const { return (*this)(clamp<int>(x, 0, m_width - 1), clamp<int>(y, 0, m_height - 1)); }
	inline color_quad_u8& get_clamped(int x, int y) { return (*this)(clamp<int>(x, 0, m_width - 1), clamp<int>(y, 0, m_height - 1)); }

	inline image_u8& set_clipped(int x, int y, const color_quad_u8& c)
	{
		if ((static_cast<uint32_t>(x) < m_width) && (static_cast<uint32_t>(y) < m_height))
			(*this)(x, y) = c;
		return *this;
	}

	image_u8& crop_dup_borders(uint32_t w, uint32_t h)
	{
		const uint32_t orig_w = m_width, orig_h = m_height;

		crop(w, h);

		if (orig_w && orig_h)
		{
			if (m_width > orig_w)
			{
				for (uint32_t x = orig_w; x < m_width; x++)
					for (uint32_t y = 0; y < m_height; y++)
						set_clipped(x, y, get_clamped(std::min(x, orig_w - 1U), std::min(y, orig_h - 1U)));
			}

			if (m_height > orig_h)
			{
				for (uint32_t y = orig_h; y < m_height; y++)
					for (uint32_t x = 0; x < m_width; x++)
						set_clipped(x, y, get_clamped(std::min(x, orig_w - 1U), std::min(y, orig_h - 1U)));
			}
		}
		return *this;
	}

	image_u8& crop(uint32_t new_width, uint32_t new_height)
	{
		if ((m_width == new_width) && (m_height == new_height))
			return *this;

		image_u8 new_image(new_width, new_height);

		const uint32_t w = std::min(m_width, new_width);
		const uint32_t h = std::min(m_height, new_height);

		for (uint32_t y = 0; y < h; y++)
			for (uint32_t x = 0; x < w; x++)
				new_image(x, y) = (*this)(x, y);

		return swap(new_image);
	}

	image_u8& swap(image_u8& other)
	{
		std::swap(m_width, other.m_width);
		std::swap(m_height, other.m_height);
		std::swap(m_pixels, other.m_pixels);
		return *this;
	}

	inline void get_block(uint32_t bx, uint32_t by, uint32_t width, uint32_t height, color_quad_u8* pPixels)
	{
		assert((bx * width + width) <= m_width);
		assert((by * height + height) <= m_height);

		for (uint32_t y = 0; y < height; y++)
			memcpy(pPixels + y * width, &(*this)(bx * width, by * height + y), width * sizeof(color_quad_u8));
	}

	inline void set_block(uint32_t bx, uint32_t by, uint32_t width, uint32_t height, const color_quad_u8* pPixels)
	{
		assert((bx * width + width) <= m_width);
		assert((by * height + height) <= m_height);

		for (uint32_t y = 0; y < height; y++)
			memcpy(&(*this)(bx * width, by * height + y), pPixels + y * width, width * sizeof(color_quad_u8));
	}

	image_u8& swizzle(uint32_t r, uint32_t g, uint32_t b, uint32_t a)
	{
		assert((r | g | b | a) <= 3);
		for (uint32_t y = 0; y < m_height; y++)
		{
			for (uint32_t x = 0; x < m_width; x++)
			{
				color_quad_u8 tmp((*this)(x, y));
				(*this)(x, y).set(tmp[r], tmp[g], tmp[b], tmp[a]);
			}
		}

		return *this;
	}

private:
	color_quad_u8_vec m_pixels;
	uint32_t m_width, m_height;
};

bool load_png(const char* pFilename, image_u8& img);

bool save_png(const char* pFilename, const image_u8& img, bool save_alpha);

class image_metrics
{
public:
	double m_max, m_mean, m_mean_squared, m_root_mean_squared, m_peak_snr;

	image_metrics()
	{
		clear();
	}

	void clear()
	{
		memset(this, 0, sizeof(*this));
	}

	void compute(const image_u8& a, const image_u8& b, uint32_t first_channel, uint32_t num_channels)
	{
		const bool average_component_error = true;

		const uint32_t width = std::min(a.width(), b.width());
		const uint32_t height = std::min(a.height(), b.height());

		assert((first_channel < 4U) && (first_channel + num_channels <= 4U));

		// Histogram approach originally due to Charles Bloom.
		double hist[256];
		memset(hist, 0, sizeof(hist));

		for (uint32_t y = 0; y < height; y++)
		{
			for (uint32_t x = 0; x < width; x++)
			{
				const color_quad_u8& ca = a(x, y);
				const color_quad_u8& cb = b(x, y);

				if (!num_channels)
					hist[iabs(ca.get_luma() - cb.get_luma())]++;
				else
				{
					for (uint32_t c = 0; c < num_channels; c++)
						hist[iabs(ca[first_channel + c] - cb[first_channel + c])]++;
				}
			}
		}

		m_max = 0;
		double sum = 0.0f, sum2 = 0.0f;
		for (uint32_t i = 0; i < 256; i++)
		{
			if (!hist[i])
				continue;

			m_max = std::max<double>(m_max, i);

			double x = i * hist[i];

			sum += x;
			sum2 += i * x;
		}

		// See http://richg42.blogspot.com/2016/09/how-to-compute-psnr-from-old-berkeley.html
		double total_values = width * height;

		if (average_component_error)
			total_values *= clamp<uint32_t>(num_channels, 1, 4);

		m_mean = clamp<double>(sum / total_values, 0.0f, 255.0f);
		m_mean_squared = clamp<double>(sum2 / total_values, 0.0f, 255.0f * 255.0f);

		m_root_mean_squared = sqrt(m_mean_squared);

		if (!m_root_mean_squared)
			m_peak_snr = 100.0f;
		else
			m_peak_snr = clamp<double>(log10(255.0f / m_root_mean_squared) * 20.0f, 0.0f, 100.0f);
	}
};

enum eZero { cZero };

template <uint32_t N, typename T>
class vec
{
protected:
	T m_v[N];

public:
	enum { num_elements = N };

	inline vec() { }

	inline vec(eZero) { set_zero(); }

	explicit inline vec(T val) { set(val); }
	inline vec(T v0, T v1) { set(v0, v1); }
	inline vec(T v0, T v1, T v2) { set(v0, v1, v2); }
	inline vec(T v0, T v1, T v2, T v3) { set(v0, v1, v2, v3); }
	inline vec(const vec& other) { for (uint32_t i = 0; i < N; i++) m_v[i] = other.m_v[i]; }
	template <uint32_t OtherN, typename OtherT> inline vec(const vec<OtherN, OtherT>& other) { set(other); }

	inline T operator[](uint32_t i) const { assert(i < N); return m_v[i]; }
	inline T& operator[](uint32_t i) { assert(i < N); return m_v[i]; }

	inline T getX() const { return m_v[0]; }
	inline T getY() const { static_assert(N >= 2, "N too small"); return m_v[1]; }
	inline T getZ() const { static_assert(N >= 3, "N too small"); return m_v[2]; }
	inline T getW() const { static_assert(N >= 4, "N too small"); return m_v[3]; }

	inline bool operator==(const vec& rhs) const { for (uint32_t i = 0; i < N; i++) if (m_v[i] != rhs.m_v[i]) return false;	return true; }
	inline bool operator<(const vec& rhs) const { for (uint32_t i = 0; i < N; i++) { if (m_v[i] < rhs.m_v[i]) return true; else if (m_v[i] != rhs.m_v[i]) return false; } return false; }

	inline void set_zero() { for (uint32_t i = 0; i < N; i++) m_v[i] = 0; }

	template <uint32_t OtherN, typename OtherT>
	inline vec& set(const vec<OtherN, OtherT>& other)
	{
		uint32_t i;
		if ((const void*)(&other) == (const void*)(this))
			return *this;
		const uint32_t m = std::min(OtherN, N);
		for (i = 0; i < m; i++)
			m_v[i] = static_cast<T>(other[i]);
		for (; i < N; i++)
			m_v[i] = 0;
		return *this;
	}

	inline vec& set_component(uint32_t index, T val) { assert(index < N); m_v[index] = val; return *this; }
	inline vec& set(T val) { for (uint32_t i = 0; i < N; i++) m_v[i] = val; return *this; }
	inline void clear_elements(uint32_t s, uint32_t e) { assert(e <= N); for (uint32_t i = s; i < e; i++) m_v[i] = 0; }

	inline vec& set(T v0, T v1)
	{
		m_v[0] = v0;
		if (N >= 2)
		{
			m_v[1] = v1;
			clear_elements(2, N);
		}
		return *this;
	}

	inline vec& set(T v0, T v1, T v2)
	{
		m_v[0] = v0;
		if (N >= 2)
		{
			m_v[1] = v1;
			if (N >= 3)
			{
				m_v[2] = v2;
				clear_elements(3, N);
			}
		}
		return *this;
	}

	inline vec& set(T v0, T v1, T v2, T v3)
	{
		m_v[0] = v0;
		if (N >= 2)
		{
			m_v[1] = v1;
			if (N >= 3)
			{
				m_v[2] = v2;

				if (N >= 4)
				{
					m_v[3] = v3;
					clear_elements(5, N);
				}
			}
		}
		return *this;
	}

	inline vec& operator=(const vec& rhs) { if (this != &rhs) for (uint32_t i = 0; i < N; i++) m_v[i] = rhs.m_v[i]; return *this; }
	template <uint32_t OtherN, typename OtherT> inline vec& operator=(const vec<OtherN, OtherT>& rhs) { set(rhs); return *this; }

	inline const T* get_ptr() const { return reinterpret_cast<const T*>(&m_v[0]); }
	inline T* get_ptr() { return reinterpret_cast<T*>(&m_v[0]); }

	inline vec operator- () const { vec res; for (uint32_t i = 0; i < N; i++) res.m_v[i] = -m_v[i]; return res; }
	inline vec operator+ () const { return *this; }
	inline vec& operator+= (const vec& other) { for (uint32_t i = 0; i < N; i++) m_v[i] += other.m_v[i]; return *this; }
	inline vec& operator-= (const vec& other) { for (uint32_t i = 0; i < N; i++) m_v[i] -= other.m_v[i]; return *this; }
	inline vec& operator/= (const vec& other) { for (uint32_t i = 0; i < N; i++) m_v[i] /= other.m_v[i]; return *this; }
	inline vec& operator*=(const vec& other) { for (uint32_t i = 0; i < N; i++) m_v[i] *= other.m_v[i]; return *this; }
	inline vec& operator/= (T s) { for (uint32_t i = 0; i < N; i++) m_v[i] /= s; return *this; }
	inline vec& operator*= (T s) { for (uint32_t i = 0; i < N; i++) m_v[i] *= s; return *this; }

	friend inline vec operator+(const vec& lhs, const vec& rhs) { vec res; for (uint32_t i = 0; i < N; i++) res.m_v[i] = lhs.m_v[i] + rhs.m_v[i]; return res; }
	friend inline vec operator-(const vec& lhs, const vec& rhs) { vec res; for (uint32_t i = 0; i < N; i++) res.m_v[i] = lhs.m_v[i] - rhs.m_v[i]; return res; }
	friend inline vec operator*(const vec& lhs, T val) { vec res; for (uint32_t i = 0; i < N; i++) res.m_v[i] = lhs.m_v[i] * val; return res; }
	friend inline vec operator*(T val, const vec& rhs) { vec res; for (uint32_t i = 0; i < N; i++) res.m_v[i] = val * rhs.m_v[i]; return res; }
	friend inline vec operator/(const vec& lhs, T val) { vec res; for (uint32_t i = 0; i < N; i++) res.m_v[i] = lhs.m_v[i] / val; return res; }
	friend inline vec operator/(const vec& lhs, const vec& rhs) { vec res; for (uint32_t i = 0; i < N; i++) res.m_v[i] = lhs.m_v[i] / rhs.m_v[i]; return res; }

	static inline T dot_product(const vec& lhs, const vec& rhs) { T res = lhs.m_v[0] * rhs.m_v[0]; for (uint32_t i = 1; i < N; i++) res += lhs.m_v[i] * rhs.m_v[i]; return res; }

	inline T dot(const vec& rhs) const { return dot_product(*this, rhs); }

	inline T norm() const { return dot_product(*this, *this); }
	inline T length() const { return sqrt(norm()); }

	inline T squared_distance(const vec& other) const { T d2 = 0; for (uint32_t i = 0; i < N; i++) { T d = m_v[i] - other.m_v[i]; d2 += d * d; } return d2; }
	inline double squared_distance_d(const vec& other) const { double d2 = 0; for (uint32_t i = 0; i < N; i++) { double d = (double)m_v[i] - (double)other.m_v[i]; d2 += d * d; } return d2; }

	inline T distance(const vec& other) const { return static_cast<T>(sqrt(squared_distance(other))); }
	inline double distance_d(const vec& other) const { return sqrt(squared_distance_d(other)); }

	inline vec& normalize_in_place() { T len = length(); if (len != 0.0f) *this *= (1.0f / len);	return *this; }

	inline vec& clamp(T l, T h)
	{
		for (uint32_t i = 0; i < N; i++)
			m_v[i] = clamp(m_v[i], l, h);
		return *this;
	}

	static vec component_min(const vec& a, const vec& b)
	{
		vec res;
		for (uint32_t i = 0; i < N; i++)
			res[i] = std::min(a[i], b[i]);
		return res;
	}

	static vec component_max(const vec& a, const vec& b)
	{
		vec res;
		for (uint32_t i = 0; i < N; i++)
			res[i] = std::max(a[i], b[i]);
		return res;
	}
};

typedef vec<4, float> vec4F;
typedef std::vector<vec4F> vec4F_vec;

inline int posmod(int x, int y)
{
	if (x >= 0)
		return (x < y) ? x : (x % y);
	int m = (-x) % y;
	return (m != 0) ? (y - m) : m;
}

class imagef
{
public:
	imagef() :
		m_width(0), m_height(0), m_pitch(0)
	{
	}

	imagef(uint32_t w, uint32_t h, uint32_t p = UINT32_MAX) :
		m_width(0), m_height(0), m_pitch(0)
	{
		resize(w, h, p);
	}

	imagef(const imagef& other) :
		m_width(0), m_height(0), m_pitch(0)
	{
		*this = other;
	}

	imagef& swap(imagef& other)
	{
		std::swap(m_width, other.m_width);
		std::swap(m_height, other.m_height);
		std::swap(m_pitch, other.m_pitch);
		m_pixels.swap(other.m_pixels);
		return *this;
	}

	imagef& operator= (const imagef& rhs)
	{
		if (this != &rhs)
		{
			m_width = rhs.m_width;
			m_height = rhs.m_height;
			m_pitch = rhs.m_pitch;
			m_pixels = rhs.m_pixels;
		}
		return *this;
	}

	imagef& clear()
	{
		m_width = 0;
		m_height = 0;
		m_pitch = 0;
		m_pixels.resize(0);
		return *this;
	}

	imagef& set(const image_u8& src, const vec4F& scale = vec4F(1), const vec4F& bias = vec4F(0))
	{
		const uint32_t width = src.width();
		const uint32_t height = src.height();

		resize(width, height);

		for (int y = 0; y < (int)height; y++)
		{
			for (uint32_t x = 0; x < width; x++)
			{
				const color_quad_u8& src_pixel = src(x, y);
				(*this)(x, y).set((float)src_pixel.r() * scale[0] + bias[0], (float)src_pixel.g() * scale[1] + bias[1], (float)src_pixel.b() * scale[2] + bias[2], (float)src_pixel.a() * scale[3] + bias[3]);
			}
		}

		return *this;
	}

	imagef& resize(const imagef& other, uint32_t p = UINT32_MAX, const vec4F& background = vec4F(0, 0, 0, 1))
	{
		return resize(other.get_width(), other.get_height(), p, background);
	}

	imagef& resize(uint32_t w, uint32_t h, uint32_t p = UINT32_MAX, const vec4F& background = vec4F(0, 0, 0, 1))
	{
		return crop(w, h, p, background);
	}

	imagef& set_all(const vec4F& c)
	{
		for (uint32_t i = 0; i < m_pixels.size(); i++)
			m_pixels[i] = c;
		return *this;
	}

	imagef& fill_box(uint32_t x, uint32_t y, uint32_t w, uint32_t h, const vec4F& c)
	{
		for (uint32_t iy = 0; iy < h; iy++)
			for (uint32_t ix = 0; ix < w; ix++)
				set_clipped(x + ix, y + iy, c);
		return *this;
	}

	imagef& crop(uint32_t w, uint32_t h, uint32_t p = UINT32_MAX, const vec4F& background = vec4F(0, 0, 0, 1))
	{
		if (p == UINT32_MAX)
			p = w;

		if ((w == m_width) && (m_height == h) && (m_pitch == p))
			return *this;

		if ((!w) || (!h) || (!p))
		{
			clear();
			return *this;
		}

		vec4F_vec cur_state;
		cur_state.swap(m_pixels);

		m_pixels.resize(p * h);

		for (uint32_t y = 0; y < h; y++)
		{
			for (uint32_t x = 0; x < w; x++)
			{
				if ((x < m_width) && (y < m_height))
					m_pixels[x + y * p] = cur_state[x + y * m_pitch];
				else
					m_pixels[x + y * p] = background;
			}
		}

		m_width = w;
		m_height = h;
		m_pitch = p;

		return *this;
	}

	inline const vec4F& operator() (uint32_t x, uint32_t y) const { assert(x < m_width&& y < m_height); return m_pixels[x + y * m_pitch]; }
	inline vec4F& operator() (uint32_t x, uint32_t y) { assert(x < m_width&& y < m_height); return m_pixels[x + y * m_pitch]; }

	inline const vec4F& get_clamped(int x, int y) const { return (*this)(clamp<int>(x, 0, m_width - 1), clamp<int>(y, 0, m_height - 1)); }
	inline vec4F& get_clamped(int x, int y) { return (*this)(clamp<int>(x, 0, m_width - 1), clamp<int>(y, 0, m_height - 1)); }

	inline const vec4F& get_clamped_or_wrapped(int x, int y, bool wrap_u, bool wrap_v) const
	{
		x = wrap_u ? posmod(x, m_width) : clamp<int>(x, 0, m_width - 1);
		y = wrap_v ? posmod(y, m_height) : clamp<int>(y, 0, m_height - 1);
		return m_pixels[x + y * m_pitch];
	}

	inline vec4F& get_clamped_or_wrapped(int x, int y, bool wrap_u, bool wrap_v)
	{
		x = wrap_u ? posmod(x, m_width) : clamp<int>(x, 0, m_width - 1);
		y = wrap_v ? posmod(y, m_height) : clamp<int>(y, 0, m_height - 1);
		return m_pixels[x + y * m_pitch];
	}

	inline imagef& set_clipped(int x, int y, const vec4F& c)
	{
		if ((static_cast<uint32_t>(x) < m_width) && (static_cast<uint32_t>(y) < m_height))
			(*this)(x, y) = c;
		return *this;
	}

	// Very straightforward blit with full clipping. Not fast, but it works.
	imagef& blit(const imagef& src, int src_x, int src_y, int src_w, int src_h, int dst_x, int dst_y)
	{
		for (int y = 0; y < src_h; y++)
		{
			const int sy = src_y + y;
			if (sy < 0)
				continue;
			else if (sy >= (int)src.get_height())
				break;

			for (int x = 0; x < src_w; x++)
			{
				const int sx = src_x + x;
				if (sx < 0)
					continue;
				else if (sx >= (int)src.get_height())
					break;

				set_clipped(dst_x + x, dst_y + y, src(sx, sy));
			}
		}

		return *this;
	}

	const imagef& extract_block_clamped(vec4F* pDst, uint32_t src_x, uint32_t src_y, uint32_t w, uint32_t h) const
	{
		for (uint32_t y = 0; y < h; y++)
			for (uint32_t x = 0; x < w; x++)
				*pDst++ = get_clamped(src_x + x, src_y + y);
		return *this;
	}

	imagef& set_block_clipped(const vec4F* pSrc, uint32_t dst_x, uint32_t dst_y, uint32_t w, uint32_t h)
	{
		for (uint32_t y = 0; y < h; y++)
			for (uint32_t x = 0; x < w; x++)
				set_clipped(dst_x + x, dst_y + y, *pSrc++);
		return *this;
	}

	inline uint32_t get_width() const { return m_width; }
	inline uint32_t get_height() const { return m_height; }
	inline uint32_t get_pitch() const { return m_pitch; }
	inline uint32_t get_total_pixels() const { return m_width * m_height; }

	inline uint32_t get_block_width(uint32_t w) const { return (m_width + (w - 1)) / w; }
	inline uint32_t get_block_height(uint32_t h) const { return (m_height + (h - 1)) / h; }
	inline uint32_t get_total_blocks(uint32_t w, uint32_t h) const { return get_block_width(w) * get_block_height(h); }

	inline const vec4F_vec& get_pixels() const { return m_pixels; }
	inline vec4F_vec& get_pixels() { return m_pixels; }

	inline const vec4F* get_ptr() const { return &m_pixels[0]; }
	inline vec4F* get_ptr() { return &m_pixels[0]; }

private:
	uint32_t m_width, m_height, m_pitch;  // all in pixels
	vec4F_vec m_pixels;
};

enum
{
	cComputeGaussianFlagNormalize = 1,
	cComputeGaussianFlagPrint = 2,
	cComputeGaussianFlagNormalizeCenterToOne = 4
};

// size_x/y should be odd
void compute_gaussian_kernel(float* pDst, int size_x, int size_y, float sigma_sqr, uint32_t flags);

void gaussian_filter(imagef& dst, const imagef& orig_img, uint32_t odd_filter_width, float sigma_sqr, bool wrapping = false, uint32_t width_divisor = 1, uint32_t height_divisor = 1);

vec4F compute_ssim(const imagef& a, const imagef& b);

vec4F compute_ssim(const image_u8& a, const image_u8& b, bool luma);

struct block8
{
	uint64_t m_vals[1];
};

typedef std::vector<block8> block8_vec;

struct block16
{
	uint64_t m_vals[2];
};

typedef std::vector<block16> block16_vec;

bool save_dds(const char* pFilename, uint32_t width, uint32_t height, const void* pBlocks, uint32_t pixel_format_bpp, DXGI_FORMAT dxgi_format, bool srgb, bool force_dx10_header);

void strip_extension(std::string& s);
void strip_path(std::string& s);

uint32_t hash_hsieh(const uint8_t* pBuf, size_t len);

#ifdef _MSC_VER
#pragma warning (pop)
#endif