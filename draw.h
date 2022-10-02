Tigr* screen;

v2 world_to_screen(v2 p)
{
	p.y = -p.y;
	float half_screen_width = 640.0f / 2.0f;
	float half_screen_height = 480.0f / 2.0f;
	p.x += half_screen_width;
	p.y += half_screen_height;
	return p;
}

v2 screen_to_world(int x, int y)
{
	v2 p = v2((float)x, (float)y);
	float half_screen_width = 640.0f / 2.0f;
	float half_screen_height = 480.0f / 2.0f;
	p.x -= half_screen_width;
	p.y -= half_screen_height;
	p.y = -p.y;
	return p;
}

v2 mouse()
{
	int x, y, buttons;
	tigrMouse(screen, &x, &y, &buttons);
	return screen_to_world(x, y);
}

void draw_point(v2 p, TPixel color)
{
	p = world_to_screen(p);
	tigrPlot(screen, (int)p.x, (int)p.y, color);
}

void draw_line(v2 a, v2 b, TPixel color)
{
	a = world_to_screen(a);
	b = world_to_screen(b);
	tigrLine(screen, (int)a.x, (int)a.y, (int)b.x, (int)b.y, color);
}

void draw_box(aabb box, TPixel color)
{
	float w = width(box) + 1;
	float h = height(box) + 1;
	draw_line(box.min, box.min + v2(w, 0), color);
	draw_line(box.min + v2(0,1), box.min + v2(0, h-1), color);
	draw_line(box.max, box.max - v2(w, 0), color);
	draw_line(box.max - v2(0,1), box.max - v2(0, h-1), color);
}

void draw_circle(circle c, TPixel color)
{
	c.p = world_to_screen(c.p);
	tigrCircle(screen, (int)c.p.x, (int)c.p.y, (int)c.r, color);
}

void draw_circle_fill(circle c, TPixel color)
{
	c.p = world_to_screen(c.p);
	tigrFillCircle(screen, (int)c.p.x, (int)c.p.y, (int)c.r, color);
}

void draw_vector(v2 p, v2 v, TPixel color)
{
	v2 arrow[] = {
		v2(0.0f, 0.0f),
		v2(-5.0f, 5.0f),
		v2(0.0f, 0.0f),
		v2(-5.0f, -5.0f),
	};
	v2 n = norm(v);
	rotation r = rotation(n.y, n.x);
	for (int i = 0; i < 4; ++i) {
		arrow[i] = mul(r, arrow[i]);
		arrow[i] += p + v;
	}
	draw_line(arrow[0], arrow[1], color);
	draw_line(arrow[2], arrow[3], color);
	draw_line(p, p + v, color);
}

void draw_polygon(polygon poly, TPixel color)
{
	v2 a = poly.verts[poly.count - 1];
	for (int i = 0; i < poly.count; ++i) {
		v2 b = poly.verts[i];
		draw_line(a, b, color);
		a = b;
	}
}

TPixel color_white(int alpha = 0xFF) { return tigrRGBA(0xFF, 0xFF, 0xFF, alpha); }
TPixel color_black(int alpha = 0xFF) { return tigrRGBA(0, 0, 0, alpha); }
TPixel color_red(int alpha = 0xFF) { return tigrRGBA(0xFF, 0, 0, alpha); }
TPixel color_green(int alpha = 0xFF) { return tigrRGBA(0, 0xFF, 0, alpha); }
TPixel color_blue(int alpha = 0xFF) { return tigrRGBA(0, 0, 0xFF, alpha); }