/// @description Insert description here
// You can write your code in this editor
//show_debug_overlay(true);

enum CELL {
	solid,
	fluid,
	air
}

width   = 500;
height  = 500;
spacing = 20;
density = 100;
fNumX   = floor(width / spacing) + 1;
fNumY   = floor(height / spacing) + 1;
fNumX2  = fNumX - 2;
fNumY2  = fNumY - 2;
h       = max(width / fNumX, height / fNumY);
fInvSpacing = 1.0 / h;
fNumCells   = fNumX * fNumY;

u     = array_create(fNumCells);
v     = array_create(fNumCells);
du    = array_create(fNumCells);
dv    = array_create(fNumCells);
prevU = array_create(fNumCells);
prevV = array_create(fNumCells);
p     = array_create(fNumCells);
s     = array_create(fNumCells);
cellType  = array_create(fNumCells);
cellColor = array_create(3 * fNumCells);

maxParticles = 10000;

particlePos   = array_create(2 * maxParticles);
particleColor = array_create(3 * maxParticles);
for (var i = 0; i < maxParticles; i++)
	particleColor[3 * i + 2] = 1.0;

particleVel         = array_create(2 * maxParticles);
particleDensity     = array_create(fNumCells);
particleRestDensity = 0;

particleRadius = 0.3 * h;
pInvSpacing    = 1.0 / (2.2 * particleRadius);
pNumX          = floor(width * pInvSpacing) + 1;
pNumY          = floor(height * pInvSpacing) + 1;
pNumCells      = pNumX * pNumY;

velocityDamping = 0.9;

numCellParticles  = array_create(pNumCells);
firstCellParticle = array_create(pNumCells + 1);
cellParticleIds   = array_create(maxParticles);

numParticles = 500;
for (var i = 0; i < numParticles; i++) {
	particlePos[i * 2]     = random(width);
	particlePos[i * 2 + 1] = random(height);
}

dt = 0.5;
iteration = 4;

g  = 50;
flipRatio         = 0.5;
numPressureIters  = 4;
numParticleIters  = 8;
overRelaxation    = 1.8;

compensateDrift   = true;
separateParticles = true;

obstacleX      = 0;
obstacleY      = 0;
obstacleVelX   = 0;
obstacleVelY   = 0;
obstacleRadius = 32;

calculateColor = false;

function integrateParticles(dt, g) { #region
	for (var i = 0; i < numParticles; i++) {
		particleVel[2 * i + 1] += dt * g;
		
		particlePos[2 * i]     += particleVel[2 * i] * dt;
		particlePos[2 * i + 1] += particleVel[2 * i + 1] * dt;
	}
} #endregion

function pushParticlesApart(numIters) { #region
	var colorDiffusionCoeff = 0.001;

	// count particles per cell
	for( var i = 0; i < pNumCells; i++ ) 
		numCellParticles[i] = 0;

	for (var i = 0; i < numParticles; i++) {
		var _x = particlePos[2 * i];
		var _y = particlePos[2 * i + 1];
		
		var xi = clamp(floor(_x * pInvSpacing), 0, pNumX - 1);
		var yi = clamp(floor(_y * pInvSpacing), 0, pNumY - 1);
		var cellNr = xi * pNumY + yi;
		numCellParticles[cellNr]++;
	}

	// partial sums
	
	var first = 0;

	for (var i = 0; i < pNumCells; i++) {
		first += numCellParticles[i];
		firstCellParticle[i] = first;
	}
	firstCellParticle[pNumCells] = first;		// guard

	// fill particles into cells

	for (var i = 0; i < numParticles; i++) {
		var _x = particlePos[2 * i];
		var _y = particlePos[2 * i + 1];

		var xi = clamp(floor(_x * pInvSpacing), 0, pNumX - 1);
		var yi = clamp(floor(_y * pInvSpacing), 0, pNumY - 1);
		var cellNr = xi * pNumY + yi;
		firstCellParticle[cellNr]--;
		cellParticleIds[firstCellParticle[cellNr]] = i;
	}

	// push particles apart

	var minDist  = 2.0 * particleRadius;
	var minDist2 = minDist * minDist;

	repeat(numIters) {
		for (var i = 0; i < numParticles; i++) {
			var px = particlePos[2 * i];
			var py = particlePos[2 * i + 1];

			var pxi = floor(px * pInvSpacing);
			var pyi = floor(py * pInvSpacing);
			var x0  = max(pxi - 1, 0);
			var y0  = max(pyi - 1, 0);
			var x1  = min(pxi + 1, pNumX - 1);
			var y1  = min(pyi + 1, pNumY - 1);

			for (var xi = x0; xi <= x1; xi++)
			for (var yi = y0; yi <= y1; yi++) {
				var cellNr = xi * pNumY + yi;
				var first  = firstCellParticle[cellNr];
				var last   = firstCellParticle[cellNr + 1];
				
				for (var j = first; j < last; j++) {
					var _id = cellParticleIds[j];
					if (_id == i) continue;
					
					var qx = particlePos[2 * _id];
					var qy = particlePos[2 * _id + 1];

					var dx = qx - px;
					var dy = qy - py;
					var d2 = dx * dx + dy * dy;
					if (d2 > minDist2 || d2 == 0) continue;
					
					var d = sqrt(d2);
					var s = 0.5 * (minDist - d) / d;
					dx *= s;
					dy *= s;
					particlePos[2 * i]       -= dx;
					particlePos[2 * i + 1]   -= dy;
					particlePos[2 * _id]     += dx;
					particlePos[2 * _id + 1] += dy;
					
					// diffuse colors

					//for (var k = 0; k < 3; k++) {
					//	var color0 = particleColor[3 * i + k];
					//	var color1 = particleColor[3 * _id + k];
					//	var color = (color0 + color1) * 0.5;
						
					//	particleColor[3 * i + k]   = color0 + (color - color0) * colorDiffusionCoeff;
					//	particleColor[3 * _id + k] = color1 + (color - color1) * colorDiffusionCoeff;
					//}
				}
			}
		}
	}
} #endregion

function handleParticleCollisions(obstacleX, obstacleY, obstacleRadius) { #region
	var h   = 1.0 / fInvSpacing;
	var r   = particleRadius;
	var _or = obstacleRadius;
	var or2 = _or * _or;
	var minDist  = obstacleRadius + r;
	var minDist2 = minDist * minDist;
	
	var minX = h + r;
	var maxX = (fNumX - 1) * h - r;
	var minY = h + r;
	var maxY = (fNumY - 1) * h - r;
	
	for (var i = 0; i < numParticles; i++) {
		var _x = particlePos[2 * i];
		var _y = particlePos[2 * i + 1];

		var dx = _x - obstacleX;
		var dy = _y - obstacleY;
		var d2 = dx * dx + dy * dy;

		// obstacle collision

		if (d2 < minDist2) {
			particleVel[2 * i]     = obstacleVelX;
			particleVel[2 * i + 1] = obstacleVelY;
		}

		// wall collisions

		if (_x < minX) {
			_x = minX;
			particleVel[2 * i] = 0;

		}
		if (_x > maxX) {
			_x = maxX;
			particleVel[2 * i] = 0;
		}
		if (_y < minY) {
			_y = minY;
			particleVel[2 * i + 1] = 0;
		}
		if (_y > maxY) {
			_y = maxY;
			particleVel[2 * i + 1] = 0;
		}
		particlePos[2 * i] = _x;
		particlePos[2 * i + 1] = _y;
	}
} #endregion

function updateParticleDensity() { #region
	var n  = fNumY;
	var h1 = fInvSpacing;
	var h2 = 0.5 * h;

	var d = particleDensity;
	for( var i = 0; i < fNumCells; i++ ) 
		d[i] = 0;

	for (var i = 0; i < numParticles; i++) {
		var _x = particlePos[2 * i];
		var _y = particlePos[2 * i + 1];

		_x = clamp(_x, h, (fNumX - 1) * h);
		_y = clamp(_y, h, (fNumY - 1) * h);

		var x0 = floor((_x - h2) * h1);
		var tx = ((_x - h2) - x0 * h) * h1;
		var x1 = min(x0 + 1, fNumX2);
			
		var y0 = floor((_y - h2) * h1);
		var ty = ((_y - h2) - y0 * h) * h1;
		var y1 = min(y0 + 1, fNumY2);
		
		var sx = 1.0 - tx;
		var sy = 1.0 - ty;
		
		if (x0 < fNumX && y0 < fNumY) d[x0 * n + y0] += sx * sy;
		if (x1 < fNumX && y0 < fNumY) d[x1 * n + y0] += tx * sy;
		if (x1 < fNumX && y1 < fNumY) d[x1 * n + y1] += tx * ty;
		if (x0 < fNumX && y1 < fNumY) d[x0 * n + y1] += sx * ty;
	}
	
	if (particleRestDensity == 0) {
		var sum = 0;
		var numFluidCells = 0;
		
		for (var i = 0; i < fNumCells; i++) {
			if (cellType[i] == CELL.fluid) {
				sum += d[i];
				numFluidCells++;
			}
		}
		
		if (numFluidCells > 0)
			particleRestDensity = sum / numFluidCells;
	}
} #endregion

function transferVelocities(toGrid, flipRatio) { #region
	var n  = fNumY;
	var h1 = fInvSpacing;
	var h2 = 0.5 * h;

	if (toGrid) {
		for( var i = 0; i < fNumCells; i++ ) {
			prevU[i] = u[i] * velocityDamping;
			prevV[i] = v[i] * velocityDamping;
			
			du[i] = 0;
			dv[i] = 0;
			u[i]  = 0;
			v[i]  = 0;
		}
		
		for (var i = 0; i < fNumCells; i++) 
			cellType[i] = s[i] == 0 ? CELL.solid : CELL.air;

		for (var i = 0; i < numParticles; i++) {
			var _x = particlePos[2 * i];
			var _y = particlePos[2 * i + 1];
			var xi = clamp(floor(_x * h1), 0, fNumX - 1);
			var yi = clamp(floor(_y * h1), 0, fNumY - 1);
			var cellNr = xi * n + yi;
			if (cellType[cellNr] == CELL.air)
				cellType[cellNr] = CELL.fluid;
		}
	}

	var _fNxH = (fNumX - 1) * h;
	var _fNyH = (fNumY - 1) * h;

	for (var _com = 0; _com < 2; _com++) {
		var dx = _com == 0? 0 : h2;
		var dy = _com == 0? h2 : 0;

		var f     = _com == 0? u : v;
		var prevF = _com == 0? prevU : prevV;
		var d     = _com == 0? du : dv;
		
		for (var i = 0; i < numParticles; i++) {
			var _x = particlePos[2 * i];
			var _y = particlePos[2 * i + 1];

			_x = clamp(_x, h, _fNxH);
			_y = clamp(_y, h, _fNyH);
			
			var x0 = min(floor((_x - dx) * h1), fNumX2);
			var tx = ((_x - dx) - x0 * h) * h1;
			var x1 = min(x0 + 1, fNumX2);
				
			var y0 = min(floor((_y-dy) * h1), fNumY2);
			var ty = ((_y - dy) - y0 * h) * h1;
			var y1 = min(y0 + 1, fNumY2);

			var sx = 1.0 - tx;
			var sy = 1.0 - ty;
			
			var d0 = sx * sy;
			var d1 = tx * sy;
			var d2 = tx * ty;
			var d3 = sx * ty;

			var nr0 = x0 * n + y0;
			var nr1 = x1 * n + y0;
			var nr2 = x1 * n + y1;
			var nr3 = x0 * n + y1;

			if (toGrid) {
				var pv = particleVel[2 * i + _com];
				f[nr0] += pv * d0; d[nr0] += d0;
				f[nr1] += pv * d1; d[nr1] += d1;
				f[nr2] += pv * d2; d[nr2] += d2;
				f[nr3] += pv * d3; d[nr3] += d3;
			} else {
				var offset = _com == 0 ? n : 1;
				var valid0 = cellType[nr0] != CELL.air || cellType[nr0 - offset] != CELL.air ? 1.0 : 0;
				var valid1 = cellType[nr1] != CELL.air || cellType[nr1 - offset] != CELL.air ? 1.0 : 0;
				var valid2 = cellType[nr2] != CELL.air || cellType[nr2 - offset] != CELL.air ? 1.0 : 0;
				var valid3 = cellType[nr3] != CELL.air || cellType[nr3 - offset] != CELL.air ? 1.0 : 0;

				var _v = particleVel[2 * i + _com];
				var _d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

				if (_d > 0) {
					var picV = (valid0 * d0 * f[nr0] + 
					            valid1 * d1 * f[nr1] + 
								valid2 * d2 * f[nr2] + 
								valid3 * d3 * f[nr3]) / _d;
								
					var corr = (valid0 * d0 * (f[nr0] - prevF[nr0]) + 
					            valid1 * d1 * (f[nr1] - prevF[nr1]) + 
								valid2 * d2 * (f[nr2] - prevF[nr2]) + 
								valid3 * d3 * (f[nr3] - prevF[nr3])) / _d;
					var flipV = _v + corr;
					
					particleVel[2 * i + _com] = (1.0 - flipRatio) * picV + flipRatio * flipV;
				}
			}
		}
		
		if (toGrid) {
			for (var i = 0; i < fNumCells; i++)
				if (d[i] > 0) f[i] /= d[i];
			
			for (var i = 0; i < fNumX; i++) 
			for (var j = 0; j < fNumY; j++) {
				var _solid = cellType[i * n + j] == CELL.solid;
				if (_solid || (i > 0 && cellType[(i - 1) * n + j] == CELL.solid))
					u[i * n + j] = prevU[i * n + j];
				if (_solid || (j > 0 && cellType[i * n + j - 1] == CELL.solid))
					v[i * n + j] = prevV[i * n + j];
			}
		}
	}
} #endregion

function solveIncompressibility(numIters, dt, overRelaxation, compensateDrift = true) { #region
	for( var i = 0; i < fNumCells; i++ ) {
		prevU[i] = u[i];
		prevV[i] = v[i];
		p[i] = 0;
	}

	var n  = fNumY;
	var cp = density * h / dt;
	
	for (var iter = 0; iter < numIters; iter++) {
		for (var i = 1; i < fNumX - 1; i++) {
			for (var j = 1; j < fNumY - 1; j++) {
				if (cellType[i * n + j] != CELL.fluid)
					continue;

				var center = i * n + j;
				var left   = (i - 1) * n + j;
				var right  = (i + 1) * n + j;
				var bottom = i * n + j - 1;
				var top    = i * n + j + 1;
				
				var sx0 = s[left];
				var sx1 = s[right];
				var sy0 = s[bottom];
				var sy1 = s[top];
				var _s  = sx0 + sx1 + sy0 + sy1;
				if (_s == 0) continue;
				
				var _div = u[right] - u[center] + v[top] - v[center];
				
				if (particleRestDensity > 0 && compensateDrift) {
					var k = 1.0;
					var compression = particleDensity[i*n + j] - particleRestDensity;
					if (compression > 0)
						_div = _div - k * compression;
				}

				var _p = -_div / _s;
				_p *= overRelaxation;
				p[center] += cp * _p;

				u[center] -= sx0 * _p;
				u[right]  += sx1 * _p;
				v[center] -= sy0 * _p;
				v[top]    += sy1 * _p;
			}
		}
	}
} #endregion

function updateParticleColors() { #region
	var h1 = fInvSpacing;

	for (var i = 0; i < numParticles; i++) {
		var s = 0.01;

		particleColor[3 * i]     = clamp(particleColor[3 * i] - s, 0, 1.0);
		particleColor[3 * i + 1] = clamp(particleColor[3 * i + 1] - s, 0, 1.0);
		particleColor[3 * i + 2] = clamp(particleColor[3 * i + 2] + s, 0, 1.0);

		var _x = particlePos[2 * i];
		var _y = particlePos[2 * i + 1];
		var xi = clamp(floor(_x * h1), 1, fNumX - 1);
		var yi = clamp(floor(_y * h1), 1, fNumY - 1);
		var cellNr = xi * fNumY + yi;

		var d0 = particleRestDensity;

		if (d0 > 0) {
			var relDensity = particleDensity[cellNr] / d0;
			if (relDensity < 0.7) {
				var s = 0.8;
				particleColor[3 * i] = s;
				particleColor[3 * i + 1] = s;
				particleColor[3 * i + 2] = 1.0;
			}
		}
	}
} #endregion

function setSciColor(cellNr, val, minVal, maxVal) { #region
	val     = min(max(val, minVal), maxVal- 0.0001);
	var d   = maxVal - minVal;
	val     = d == 0 ? 0.5 : (val - minVal) / d;
	var m   = 0.25;
	var num = floor(val / m);
	var s   = (val - num * m) / m;
	var r, g, b;

	switch (num) {
		case 0 : r = 0; g = s;       b = 1.0;     break;
		case 1 : r = 0; g = 1.0;     b = 1.0 - s; break;
		case 2 : r = s;   g = 1.0;     b = 0;     break;
		case 3 : r = 1.0; g = 1.0 - s; b = 0;     break;
	}

	cellColor[3 * cellNr]     = r;
	cellColor[3 * cellNr + 1] = g;
	cellColor[3 * cellNr + 2] = b;
} #endregion

function updateCellColors() { #region
	for (var i = 0; i < fNumCells; i++) 
		cellColor[i] = 0;

	for (var i = 0; i < fNumCells; i++) {

		if (cellType[i] == CELL.solid) {
			cellColor[3 * i]     = 0.5;
			cellColor[3 * i + 1] = 0.5;
			cellColor[3 * i + 2] = 0.5;
		} else if (cellType[i] == CELL.fluid) {
			var d = particleDensity[i];
			if (particleRestDensity > 0)
				d /= particleRestDensity;
			setSciColor(i, d, 0, 2.0);
		}
	}
} #endregion

function simulate(dt, g, flipRatio, numPressureIters, numParticleIters, overRelaxation, compensateDrift, separateParticles, obstacleX, abstacleY, obstacleRadius) { #region
	var sdt = dt / iteration;
	
	for (var step = 0; step < iteration; step++) {
		integrateParticles(sdt, g);
		if (separateParticles) pushParticlesApart(numParticleIters); 
		handleParticleCollisions(obstacleX, abstacleY, obstacleRadius)
		transferVelocities(true);
		updateParticleDensity();
		solveIncompressibility(numPressureIters, sdt, overRelaxation, compensateDrift);
		transferVelocities(false, flipRatio);
	}

	//updateParticleColors();
	//updateCellColors();
} #endregion
	
function setObstacle(_x, _y, reset) { #region
	var vx = 0;
	var vy = 0;

	if (!reset) {
		vx = (_x - obstacleX) * 5 * dt;
		vy = (_y - obstacleY) * 5 * dt;
	}
	
	obstacleX = _x;
	obstacleY = _y;
	var r  = obstacleRadius;
	var n  = fNumY;
	var cd = sqrt(2) * h;

	for (var i = 1; i < fNumX2; i++)
	for (var j = 1; j < fNumY2; j++) {
		s[i * n + j] = 1.0;
			
		var dx = (i + 0.5) * h - _x;
		var dy = (j + 0.5) * h - _y;
		
		if (dx * dx + dy * dy < r * r) {
			s[i * n + j]       = 0;
			u[i * n + j]       = vx;
			u[(i + 1) * n + j] = vx;
			v[i * n + j]       = vy;
			v[i * n + j + 1]   = vy;
		}
	}
	
	obstacleVelX = vx;
	obstacleVelY = vy;
} #endregion