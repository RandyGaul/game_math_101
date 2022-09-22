#include <math.h>
#include <float.h>
#include <stdint.h>
#include <string.h>

float min(float a, float b) { return a < b ? a : b; }
float max(float a, float b) { return a > b ? a : b; }

struct v2
{
	v2() { }
	v2(float x, float y) { this->x = x; this->y = y; }
	float x;
	float y;
};

v2 operator-(v2 a) { return v2(-a.x, -a.y); }
v2 operator+(v2 a, v2 b) { return v2(a.x + b.x, a.y + b.y); }
v2 operator+=(v2& a, v2 b) { a = v2(a.x + b.x, a.y + b.y); return a; }
v2 operator-(v2 a, v2 b) { return v2(a.x - b.x, a.y - b.y); }
v2 operator-=(v2& a, v2 b) { a = v2(a.x - b.x, a.y - b.y); return a; }
v2 operator*(v2 a, float b) { return v2(a.x * b, a.y * b); }
v2 operator*=(v2& a, float b) { a = v2(a.x * b, a.y * b); return a; }
v2 operator/(v2 a, float b) { return v2(a.x / b, a.y / b); }
v2 operator/=(v2& a, float b) { a = v2(a.x / b, a.y / b); return a; }
float dot(v2 a, v2 b) { return a.x * b.x + a.y * b.y; }
float det2(v2 a, v2 b) { return a.x * b.y - a.y * b.x; }
float len(v2 v) { return sqrtf(dot(v, v)); }
float len_squared(v2 v) { return dot(v, v); }
v2 norm(v2 v) { float l = len(v); return v * (1.0f / l); }
v2 skew(v2 v) { return v2(-v.y, v.x); }

float shortest_arc(v2 a, v2 b)
{
	a = norm(a);
	b = norm(b);
	float c = dot(a, b);
	float s = det2(a, b);
	float theta = acosf(c);
	if (s > 0) {
		return theta;
	} else {
		return -theta;
	}
}

struct rotation
{
	rotation() { }
	rotation(float angle) { s = sinf(angle); c = sinf(angle); }
	rotation(float s, float c) { this->s = s; this->c = c; }
	float s;
	float c;
};

rotation sincos(float a) { rotation r; r.c = cosf(a); r.s = sinf(a); return r; }
float atan2_360(float y, float x) { return atan2f(-y, -x) + 3.14159265f; }
float atan2_360(rotation r) { return atan2_360(r.s, r.c); }
float atan2_360(v2 v) { return atan2f(-v.y, -v.x) + 3.14159265f; }
v2 mul(rotation a, v2 b) { return v2(a.c * b.x - a.s * b.y, a.s * b.x + a.c * b.y); }

struct halfspace
{
	halfspace() { }
	halfspace(v2 n, v2 p) { this->n = n; c = dot(n, p); }
	halfspace(v2 n, float c) { this->n = n; this->c = c; }
	v2 n;
	float c;
};

float distance(halfspace h, v2 p) { return dot(h.n, p) - h.c; }

v2 bezier(v2 a, v2 b, v2 c, float t)
{
	float u = 1.0f - t;
	float ut = u * t;
	v2 auu = a * u * u;
	v2 but2 = b * ut * 2.0f;
	v2 ctt = c * t * t;
	return auu + but2 + ctt;
}

v2 bezier(v2 a, v2 b, v2 c, v2 d, float t)
{
	float u = 1 - t;
	float tt = t * t;
	float uu = u * u;
	v2 auuu = a * uu * u;
	v2 buut3 = b * uu * t * 3.0f;
	v2 cutt3 = c * u * tt * 3.0f;
	v2 dttt = d * tt * t;
	return auuu + buut3 + cutt3 + dttt;
}

struct m2
{
	v2 x;
	v2 y;
};

m2 m2_identity()
{
	m2 m;
	m.x = v2(1, 0);
	m.y = v2(0, 1);
	return m;
}

m2 m2_rotation(float angle)
{
	float c = cosf(angle);
	float s = sinf(angle);
	m2 m;
	m.x = v2(c, -s);
	m.y = v2(s, c);
	return m;
}
m2 m2_scale(float x_scale, float y_scale)
{
	m2 m = m2_identity();
	m.x *= x_scale;
	m.y *= y_scale;
	return m;
}

v2 mul(m2 m, v2 v) { return v2(m.x.x * v.x + m.y.x * v.y, m.x.y * v.x + m.y.y * v.y); }
m2 mul(m2 a, m2 b) { m2 c; c.x = mul(a, b.x); c.y = mul(a, b.y); return c; }

struct aabb
{
	aabb() { }
	aabb(v2 min, v2 max) { this->min = min; this->max = max; }
	v2 min;
	v2 max;
};

float width(aabb box) { return box.max.x - box.min.x; }
float height(aabb box) { return box.max.y - box.min.y; }
v2 center(aabb box) { return (box.min + box.max) * 0.5f; }

struct circle
{
	circle() { }
	circle(v2 p, float r) { this->p = p; this->r = r; }
	circle(float x, float y, float r) { this->p = v2(x, y); this->r = r; }
	v2 p;
	float r;
};

#define POLYGON_MAX_VERTS 8

struct polygon
{
	int count = 0;
	v2 verts[POLYGON_MAX_VERTS];
	v2 norms[POLYGON_MAX_VERTS];
	void compute_norms()
	{
		for (int i = 0; i < count; ++i) {
			int j = i + 1 < count ? i + 1 : 0;
			norms[i] = norm(skew(verts[i] - verts[j]));
		}
	}
};

struct raycast_output
{
	float t; // Time of impact.
	v2 n;    // Normal of the surface at impact (unit length).
};

struct ray
{
	ray() { }
	ray(v2 p, v2 d, float t) { this->p = p; this->d = d; this->t = t; }
	v2 p;    // Start position.
	v2 d;    // Direction of the ray (normalized)
	float t; // Distance along d the ray travels.
	v2 endpoint() { return p + d * t; }
	v2 impact(raycast_output hit_data) { return p + d * hit_data.t; }
};

// Reflect vector d across vector n. See: http://paulbourke.net/geometry/reflected/
v2 reflect(v2 d, v2 n) { return d - n * 2 * dot(d, n); }

bool raycast(ray r, circle circ, raycast_output* out)
{
	v2 e = r.p - circ.p;
	float rr = circ.r * circ.r;
	float b = dot(e, r.d);
	float c = dot(e, e) - rr;
	float discriminant = b * b - c;
	bool missed_circle = discriminant < 0;
	if (missed_circle) return false;
	float t = -b - sqrtf(discriminant);
	if (t >= 0 && t <= r.t) {
		out->t = t;
		v2 impact = r.p + r.d * t;
		out->n = norm(impact - circ.p);
		return true;
	} else {
		return false;
	}
}

bool raycast(ray r, polygon poly, raycast_output* out)
{
	float lo = 0;
	float hi = r.t;
	int index = ~0;

	for (int i = 0; i < poly.count; ++i) {
		// Calculate distance of point to plane.
		// This is a slight variation of dot(n, p) - c, where instead of pre-computing c as scalar,
		// we form a vector pointing from the ray's start point to a point on the plane. This is
		// functionally equivalent to dot(n, p) - c, but we don't have to store c inside of our poly.
		float distance = dot(poly.norms[i], poly.verts[i] - r.p);
		float denominator = dot(poly.norms[i], r.d);
		if (denominator == 0) {
			// Ray direction is parallel to this poly's face plane.
			// If the ray's start direction is outside the plane we know there's no intersection.
			if (distance > 0) return false;
		} else {
			float t = distance / denominator;
			if (denominator < 0) {
				// Ray is entering the plane.
				lo = max(lo, t);
				index = i;
			}
			else hi = min(hi, t); // Ray is exiting the plane.
			bool ray_clipped_away = lo > hi;
			if (ray_clipped_away) return false;
		}
	}

	if (index != ~0) {
		out->t = lo;
		out->n = poly.norms[index];
		return true;
	} else {
		return false;
	}
}

// Based on Andrew's Algorithm from Ericson's Real-Time Collision Detection book.
int convex_hull(v2* verts, int count)
{
	count = min(count, POLYGON_MAX_VERTS);
	if (count < 3) {
		return 0;
	}

	// Sort lexicographically (on x-axis, then y-axis).
	for (int i = 0; i < count; ++i) {
		int lo = i;
		for (int j = i+1; j < count; ++j) {
			if (verts[j].x < verts[lo].x) {
				lo = j;
			} else if (verts[j].x == verts[lo].x && verts[j].y < verts[lo].y) {
				lo = j;
			}
		}
		v2 swap = verts[i];
		verts[i] = verts[lo];
		verts[lo] = swap;
	}

	int j = 2;
	int hull[POLYGON_MAX_VERTS + 1];
	hull[0] = 0;
	hull[1] = 1;

	// Find lower-half of hull.
	for (int i = 2; i < count; ++i) {
		while (j >= 2) {
			v2 e0 = verts[hull[j-1]] - verts[hull[j-2]];
			v2 e1 = verts[i] - verts[hull[j-2]];
			if (det2(e0, e1) <= 0) --j;
			else break;
		}
		hull[j++] = i;
	}

	// Find top-half of hull.
	for (int i = count-2, k = j+1; i >= 0; --i) {
		while (j >= k) {
			v2 e0 = verts[hull[j-1]] - verts[hull[j-2]];
			v2 e1 = verts[i] - verts[hull[j-2]];
			if (det2(e0, e1) <= 0) --j;
			else break;
		}
		hull[j++] = i;
	}

	--j; // Pop the last vert off as it's a duplicate.
	if (j < 3) return 0;
	v2 hull_verts[POLYGON_MAX_VERTS];
	for (int i = 0; i < j; ++i) hull_verts[i] = verts[hull[i]];
	memcpy(verts, hull_verts, sizeof(v2) * j);
	return j;
}

/*
 * A random number generator of the type LFSR (linear feedback shift registers). This specific
 * implementation uses the XorShift+ variation, and returns 64-bit random numbers. More information
 * can be found on Wikipedia.
 * https://en.wikipedia.org/wiki/Xorshift
 *
 * This implementation comes from Mattias Gustavsson's single-file header collection.
 * https://github.com/mattiasgustavsson/libs/blob/main/rnd.h
 */
struct rnd_t
{
	uint64_t state[2];
};

uint64_t internal_rnd_murmur3_avalanche64(uint64_t h)
{
	h ^= h >> 33;
	h *= 0xff51afd7ed558ccd;
	h ^= h >> 33;
	h *= 0xc4ceb9fe1a85ec53;
	h ^= h >> 33;
	return h;
}

rnd_t rnd_seed(uint64_t seed)
{
	rnd_t rnd;
	uint64_t value = internal_rnd_murmur3_avalanche64((seed << 1ULL) | 1ULL);
	rnd.state[0] = value;
	rnd.state[1] = internal_rnd_murmur3_avalanche64(value);
	return rnd;
}

uint64_t rnd_next(rnd_t* rnd)
{
	uint64_t x = rnd->state[0];
	uint64_t y = rnd->state[1];
	rnd->state[0] = y;
	x ^= x << 23;
	x ^= x >> 17;
	x ^= y ^ (y >> 26);
	rnd->state[1] = x;
	return x + y;
}

float rnd_next_float(rnd_t* rnd)
{
	uint32_t value = (uint32_t)(rnd_next(rnd) >> 32);

	// Convert a randomized uint32_t value to a float value x in the range 0.0f <= x < 1.0f.
	// Contributed by Jonatan Hedborg.
	uint32_t exponent = 127;
	uint32_t mantissa = value >> 9;
	uint32_t result = (exponent << 23) | mantissa;
	return *(float*)&result - 1.0f;
}

double rnd_next_double(rnd_t* rnd)
{
	uint64_t value = rnd_next(rnd);
	uint64_t exponent = 1023;
	uint64_t mantissa = value >> 12;
	uint64_t result = (exponent << 52) | mantissa;
	return *(double*)&result - 1.0;
}

int rnd_next_range_int(rnd_t* rnd, int min, int max)
{
	int range = (max - min) + 1;
	int value = (int)(rnd_next(rnd) % range);
	return min + value;
}

uint64_t rnd_next_range_uint64(rnd_t* rnd, uint64_t min, uint64_t max)
{
	uint64_t range = (max - min) + 1;
	uint64_t value = rnd_next(rnd) % range;
	return min + value;
}

float rnd_next_range_float(rnd_t* rnd, float min, float max)
{
	float range = max - min;
	float value = rnd_next_float(rnd) * range;
	return min + value;
}

double rnd_next_range_double(rnd_t* rnd, double min, double max)
{
	double range = max - min;
	double value = rnd_next_float(rnd) * range;
	return min + value;
}

uint64_t rnd_next(rnd_t& rnd) { return rnd_next(&rnd); }
float    rnd_next_float(rnd_t& rnd) { return rnd_next_float(&rnd); }
double   rnd_next_double(rnd_t& rnd) { return rnd_next_double(&rnd); }
int      rnd_next_range(rnd_t& rnd, int min, int max) { return rnd_next_range_int(&rnd, min, max); }
uint64_t rnd_next_range(rnd_t& rnd, uint64_t min, uint64_t max) { return rnd_next_range_uint64(&rnd, min, max); }
float    rnd_next_range(rnd_t& rnd, float min, float max) { return rnd_next_range_float(&rnd, min, max); }
double   rnd_next_range(rnd_t& rnd, double min, double max) { return rnd_next_range_double(&rnd, min, max); }