#include <math.h>
#include <float.h>
#include <stdint.h>
#include <string.h>

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
v2 safe_norm(v2 v) { float l = len(v); return l == 0 ? v2(0,0) : v * (1.0f / l); }
v2 skew(v2 v) { return v2(-v.y, v.x); }

int min(int a, int b) { return a < b ? a : b; }
int max(int a, int b) { return a > b ? a : b; }
float min(float a, float b) { return a < b ? a : b; }
float max(float a, float b) { return a > b ? a : b; }
float abs(float a) { return a < 0 ? -a : a; }
v2 abs(v2 a) { return v2(abs(a.x), abs(a.y)); }
float approach(float t, float target, float delta) { return t < target ? min(t + delta, target) : max(t - delta, target); }
float map(float t, float lo, float hi, float old_lo = 0, float old_hi = 1) { return lo + ((t - old_lo) / (old_hi - old_lo)) * (hi - lo); }
float smoothstep(float x) { return x * x * (3.0f - 2.0f * x); }

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
v2 intersect(v2 a, v2 b, float da, float db) { return a + (b - a) * (da / (da - db)); }
v2 intersect(halfspace h, v2 a, v2 b) { return intersect(a, b, distance(h, a), distance(h, b)); }

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
	aabb(v2 p, float w, float h) { min = p - v2(w, h); max = p + v2(w, h); }
	v2 min;
	v2 max;
};

float width(aabb box) { return box.max.x - box.min.x; }
float height(aabb box) { return box.max.y - box.min.y; }
v2 center(aabb box) { return (box.min + box.max) * 0.5f; }
v2 top(aabb box) { return v2((box.min.x + box.max.x) * 0.5f, box.max.y); }
v2 bottom(aabb box) { return v2((box.min.x + box.max.x) * 0.5f, box.min.y); }
v2 left(aabb box) { return v2(box.min.x, (box.min.y + box.max.y) * 0.5f); }
v2 right(aabb box) { return v2(box.max.x, (box.min.y + box.max.y) * 0.5f); }

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

bool point_in_poly(polygon poly, v2 p)
{
	for (int i = 0; i < poly.count; ++i) {
		float c = dot(poly.norms[i], poly.verts[i]);
		float dist = dot(poly.norms[i], p) - c;
		if (dist >= 0) return false;
	}
	return true;
}

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

struct collision_data
{
	bool hit = false;
	v2 hit_spot;
	float depth;
	v2 normal;
};

collision_data circle_to_circle(circle circ_a, circle circ_b)
{
	collision_data out;
	v2 d = circ_b.p - circ_a.p;
	float d2 = dot(d, d);
	float r = circ_a.r + circ_b.r;
	if (d2 < r * r) {
		float l = len(d);
		d = l == 0 ? v2(0, 1) : d * (1.0f / l);
		out.hit = true;
		out.hit_spot = circ_b.p - d * circ_b.r;
		out.depth = r - l;
		out.normal = d;
	}
	return out;
}

collision_data aabb_to_aabb(aabb a, aabb b)
{
	collision_data out;
	v2 mid_a = (a.min + a.max) * 0.5f;
	v2 mid_b = (b.min + b.max) * 0.5f;
	v2 ea = abs((a.max - a.min) * 0.5f);
	v2 eb = abs((b.max - b.min) * 0.5f);
	v2 d = mid_b - mid_a;

	// calc overlap on x and y axes
	float dx = ea.x + eb.x - abs(d.x);
	if (dx < 0) return out;
	float dy = ea.y + eb.y - abs(d.y);
	if (dy < 0) return out;

	v2 n;
	float depth;
	v2 p;

	if (dx < dy) {
		// x axis overlap is smaller
		depth = dx;
		if (d.x < 0) {
			n = v2(-1.0f, 0);
			p = mid_a - v2(ea.x, 0);
		} else {
			n = v2(1.0f, 0);
			p = mid_a + v2(ea.x, 0);
		}
	}  else {
		// y axis overlap is smaller
		depth = dy;
		if (d.y < 0) {
			n = v2(0, -1.0f);
			p = mid_a - v2(0, ea.y);
		} else {
			n = v2(0, 1.0f);
			p = mid_a + v2(0, ea.y);
		}
	}

	out.hit = true;
	out.hit_spot = p;
	out.depth = depth;
	out.normal = n;
	return out;
}

float clamp(float v, float lo, float hi) { return min(max(v, lo), hi); }
v2 clamp(v2 v, v2 lo, v2 hi) { return v2(clamp(v.x, lo.x, hi.x), clamp(v.y, lo.y, hi.y)); }

collision_data circle_to_aabb(circle a, aabb b)
{
	collision_data out;
	v2 L = clamp(a.p, b.min, b.max);
	v2 ab = a.p - L;
	float d2 = dot(ab, ab);
	float r2 = a.r * a.r;

	if (d2 < r2) {
		if (d2 != 0) {
			// shallow (center of circle not inside of AABB)
			float d = sqrtf(d2);
			v2 n = norm(ab);
			out.hit = true;
			out.depth = a.r - d;
			out.hit_spot = a.p + n * d;
			out.normal = n;
		} else {
			// deep (center of circle inside of AABB)
			// clamp circle's center to edge of AABB, then form the manifold
			v2 mid = (b.min + b.max) * 0.5f;
			v2 e = (b.max - b.min) * 0.5f;
			v2 d = a.p - mid;
			v2 abs_d = abs(d);

			float x_overlap = e.x - abs_d.x;
			float y_overlap = e.y - abs_d.y;

			float depth;
			v2 n;

			if (x_overlap < y_overlap) {
				depth = x_overlap;
				n = v2(1.0f, 0);
				n *= (d.x < 0 ? 1.0f : -1.0f);
			} else {
				depth = y_overlap;
				n = v2(0, 1.0f);
				n *= (d.y < 0 ? 1.0f : -1.0f);
			}

			out.hit = true;
			out.depth = a.r + depth;
			out.hit_spot = a.p - n * depth;
			out.normal = n;
		}
	}

	return out;
}

struct sutherland_hodgman_output
{
	polygon front;
	polygon back;
};

bool in_front(float distance, float epsilon) { return distance > epsilon; }
bool behind(float distance, float epsilon) { return distance < -epsilon; }
bool on(float distance, float epsilon) { return !in_front(distance, epsilon) && !behind(distance, epsilon); }

// See: https://gamedevelopment.tutsplus.com/tutorials/how-to-dynamically-slice-a-convex-shape--gamedev-14479
sutherland_hodgman_output sutherland_hodgman(halfspace split, polygon in, const float k_epsilon = 1.e-4f)
{
	sutherland_hodgman_output out;
	v2 a = in.verts[in.count - 1];
	float da = distance(split, a);

	for(int i = 0; i < in.count; ++i) {
		v2 b = in.verts[i];
		float db = distance(split, b);

		if(in_front(db, k_epsilon)) {
			if(behind(da, k_epsilon)) {
				v2 i = intersect(b, a, db, da);
				out.front.verts[out.front.count++] = i;
				out.back.verts[out.back.count++] = i;
			}
			out.front.verts[out.front.count++] = b;
		} else if(behind(db, k_epsilon)) {
			if(in_front(da, k_epsilon)) {
				v2 i = intersect(a, b, da, db);
				out.front.verts[out.front.count++] = i;
				out.back.verts[out.back.count++] = i;
			} else if(on(da, k_epsilon)) {
				out.back.verts[out.back.count++] = a;
			}
			out.back.verts[out.back.count++] = b;
		} else {
			out.front.verts[out.front.count++] = b;
			if(on(da, k_epsilon)) {
				out.back.verts[out.back.count++] = b;
			}
		}

		a = b;
		da = db;
	}

	return out;
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