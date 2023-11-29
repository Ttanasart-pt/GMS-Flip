/// @description Insert description here
// You can write your code in this editor
FLIP_setObstacle(mouse_x, mouse_y, obstacleRadius, false);
FLIP_simulate(dt);										

FLIP_setParticleBuffer(buffer_get_address(particlePosBuff));
buffer_seek(particlePosBuff, buffer_seek_start, 0);
for(var i = 0; i < 2 * maxParticles; i++)
	particlePos[i] = buffer_read(particlePosBuff, buffer_f64);
