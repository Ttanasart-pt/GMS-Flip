/// @description Insert description here
// You can write your code in this editor

draw_set_color(c_white);
draw_set_alpha(1);
for( var i = 0; i < numParticles; i++ ) {
	var _x = particlePos[i * 2];
	var _y = particlePos[i * 2 + 1];
	
	draw_circle(_x, _y, particleRadius * 2, false);
}

draw_set_color(c_red);
draw_circle(mouse_x, mouse_y, obstacleRadius, false);

draw_set_color(c_white);
draw_set_halign(fa_right);
draw_text(room_width - 16, 16, fps_real);