#include "tigr.h"
#include "math_101.h"
#include "draw.h"

v2 random_v2(rnd_t& rnd)
{
	return v2(rnd_next_range(rnd, -1.0f, 1.0f), rnd_next_range(rnd, -1.0f, 1.0f));
}

int main()
{
	screen = tigrWindow(640, 480, "Math 101", 0);
	rnd_t rnd = rnd_seed(20);

	v2 verts[] = {
		v2(-100, 0),
		v2(-25, 0),
		v2(25, 0),
		v2(100, 0),
		v2(0, 50),
	};
	int count = sizeof(verts) / sizeof(*verts);

	v2 old[sizeof(verts) / sizeof(*verts)];
	memcpy(old, verts, sizeof(verts));

	count = convex_hull(verts, count);
	polygon poly;
	memcpy(poly.verts, verts, sizeof(v2) * count);
	poly.count = count;
	poly.compute_norms();

	float t = 0;
	while (!tigrClosed(screen) && !tigrKeyDown(screen, TK_ESCAPE)) {
		float dt = tigrTime();
		t += dt;
		tigrClear(screen, color_black());

		for (int i = 0; i < sizeof(old) / sizeof(*old); ++i) {
			draw_circle(circle(old[i], 5.0f), color_red());
		}

		v2 a = verts[count - 1];
		for (int i = 0; i < count; ++i) {
			v2 b = verts[i];
			draw_circle(circle(b, 5.0f), color_white());
			draw_line(a, b, color_white());
			a = b;
		}

		v2 m = mouse();
		ray r = ray(m, norm(v2(0, 55) - m), 400.0f);
		raycast_output out;
		bool hit = raycast(r, poly, &out);
		if (hit) draw_vector(r.p, r.impact(out) - r.p, color_red());
		else draw_vector(r.p, r.endpoint() - r.p, color_white());

		tigrUpdate(screen);
	}

	tigrFree(screen);

	return 0;
}
