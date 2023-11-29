/// @description Insert description here
// You can write your code in this editor
draw_set_color(c_white);
draw_set_alpha(1);
for( var i = 0; i < numParticles; i++ ) {
	var _x = particlePos[i * 2];
	var _y = particlePos[i * 2 + 1];
	
	var _r = particleColor[i * 3]     * 255;
	var _g = particleColor[i * 3 + 1] * 255;
	var _b = particleColor[i * 3 + 2] * 255;
	
	//draw_set_color(make_color_rgb(_r, _g, _b));
	//draw_circle_color(_x, _y, particleRadius, c_white, c_black, false);
	draw_circle(_x, _y, particleRadius, false);
}

for( var i = 0; i < fNumCells; i++ ) {
	var _d  = particleDensity[i];
	var _ix = floor(i / fNumY);
	var _iy =      (i % fNumY);
	if(_ix == 0 || _iy == 0 || _ix == fNumX - 2 || _iy == fNumY - 2) continue;
	
	var _x = _ix * spacing;
	var _y = _iy * spacing;
	var _u = u[i];
	var _v = v[i];
	
	draw_set_color(c_red);
	draw_line(_x, _y, _x + _u * 0.25, _y + _v * 0.25);
}

draw_set_color(c_red);
draw_circle(mouse_x, mouse_y, obstacleRadius, false);

draw_set_color(c_white);
draw_set_halign(fa_right);
draw_text(room_width - 16, 16, fps_real);