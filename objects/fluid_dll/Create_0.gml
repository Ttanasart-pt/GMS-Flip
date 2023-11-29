/// @description Insert description here
// You can write your code in this editor
show_debug_overlay(true);

width   = 500;
height  = 500;
spacing = 5;
density = 100;

maxParticles = 10000;
velocityDamping = 0.9;

numParticles = 5000;

dt = 0.25;
iteration = 4;

g  = 1;
flipRatio         = 0.5;
numPressureIters  = 4;
numParticleIters  = 8;
overRelaxation    = 1.8;

obstacleRadius = 32;

particlePos     = array_create(2 * maxParticles);
particlePosBuff = buffer_create(maxParticles * 2 * 8, buffer_fixed, 8);

FLIP_setParticleBuffer(buffer_get_address(particlePosBuff));
FLIP_initDomain(width, height, spacing, density, maxParticles);
FLIP_spawnRandomFluid(numParticles);
FLIP_setQuality(iteration, numPressureIters, numParticleIters);
FLIP_setGravity(g);
FLIP_setFlipRatio(flipRatio);
FLIP_setVelocityDamping(velocityDamping);
FLIP_setOverRelaxation(overRelaxation);

particleRadius = FLIP_getParticleRadius();