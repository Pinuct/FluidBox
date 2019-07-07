#include <stdio.h>
#include "stdlib.h"
#include <unistd.h>
#include <math.h>
#include "CollisionList.h"
#include <windows.h> 

// Implementation of sph as described in
//    Muller, et al, Particle-based fluid simulation
// Adapted for 2D in
//    http://www.cs.cornell.edu/~bindel/class/cs5220-f11/code/sph.pdf

#define h 					(0.016)          	// Particle radius
#define h2  				h * h          	// Radius squared
#define h4  				h2 * h2        	// Radius to the 4
#define h8  				h4 * h4        	// Radius to the 8
#define N 					50 			// Number of particles

#define height 					20 			// alto
#define width 					16 			// ancho

#define alto 					20 			// alto
#define ancho 					16 			// ancho

// Resistance to compression
// speed of sound = sqrt(k / rho0)
#define k 					30             	// Bulk modulus (1000)
#define mu 					3             	// Viscosity (0.1)
#define rho0 				1000          	// Reference density
#define rho02 				rho0 * 2           
#define dt2 				dt / 2       	// Half time step in seconds
#define restitution 		0.2         	// Coefficient of restituion for boundary



#define edge1 				h * 0.5
#define edge2 				1 - edge1
#define edge3 				height / width - edge1

float mass;
float i, j;

float Cp 			=		15 * k;
float Cv 			=		-40 * mu;
float dt 			=		18e-4;       	// Time step in seconds
float gravity[2] 	=		{9.8,-9.8};      	// Gravity

float C0, C1, C2;

float x[N], y[N], vx[N], vy[N], vhx[N], vhy[N], ax[N], ay[N], rho[N];

float random(float min, float max);
void gotoxy(int x, int y);
void clearScreen();

void particlesInMesh(float y1, float y2) {
    float xp = h * 0.5 + 0.01;
    float yp = y1;
    float r = h;
    
	int i;
    for (i = 0; i < N; i++) {
        // Initialize particle positions
        x[i] = xp;
        y[i] = yp;
        yp += r;
        
        if (yp > y2) {
            yp = y1;
            xp += r;
        }
        
        // Initialize particle velocities
        vx[i] = random(-0.02, 0.02);
        vy[i] = random(-0.02, 0.02);
    }
}

void computeDensities() {

	int i,j;

    // Find new densities
    float dx, dy, r2, z, rho_ij;
    float C1 = 4 * mass / (3.1416 * h2);
    float C2 = 4 * mass / (3.1416 * h8);

    // Initialise densities
    for (i = N; i--;) {
        rho[i] = C1;
    }

    for (i = 0; i < N; i++) {
        for (j = i + 1; j < N; j++) {
            dx = x[i] - x[j];
            dy = y[i] - y[j];
            r2 = dx * dx + dy * dy;
            z  = h2 - r2;
            
            if (z > 0) {
                rho_ij = C2 * z * z * z;
                rho[i] += rho_ij;
                rho[j] += rho_ij;
            }
        }
    }

}

void computeAccelerations() {

	int i,j;

    computeDensities();
    
    // Start with gravity and surface forces
    for (i = N; i--;) {
        ax[i] = gravity[0];
        ay[i] = -gravity[1];
    }
    
    // Find new densities
    float dx, dy, r2, rhoi, rhoj, q, u, w0, wp, wv, dvx, dvy;
    
    for (i = N; i--;) {
        rhoi = rho[i];
        
        for (j = i; j--;) {
            dx = x[i] - x[j];
            dy = y[i] - y[j];
            r2 = dx * dx + dy * dy;
            
            if (r2 < h2) {
                rhoj = rho[j];
                q = sqrt(r2) / h;
                u = 1 - q;
                w0 = C0 * u / (rhoi * rhoj);
                wp = w0 * Cp * (rhoi + rhoj - rho02) * u / q;
                wv = w0 * Cv;
                
                dvx = vx[i] - vx[j];
                dvy = vy[i] - vy[j];
                
                ax[i] += wp * dx + wv * dvx;
                ay[i] += wp * dy + wv * dvy;
                ax[j] -= wp * dx + wv * dvx;
                ay[j] -= wp * dy + wv * dvy;
            }
        }
    }
    
}

void updateParticles() {

	int i,j;

    CollisionList_t *collisions = CollisionList_create();
    float dx, dy, r2;
    
    // Reset properties and find collisions
    for (i = N; i--;) {
        // Reset density
        rho[i] = C1;
        
        // Reset acceleration
        ax[i] = gravity[0];
        ay[i] = -gravity[1];
        
        // Calculate which particles overlap
        for (j = i; j--;) {
            dx = x[i] - x[j];
            dy = y[i] - y[j];
            r2 = dx * dx + dy * dy;
            if (r2 < h2) {
                
                CollisionList_addElement(collisions, i, j, dx, dy, r2);
            }
        }
    }
    
    // Calculate densities
    float rho_ij, z;
    for (i = collisions->listSize; i--;) {
        Collision_t *c = CollisionList_getElement(collisions, i);
        z = h2 - c->r2;
        rho_ij = C2 * z * z * z;
        rho[c->element_i] += rho_ij;
        rho[c->element_j] += rho_ij;
    }
    
    // TODO: Find max density
    
    // Calculate accelerations

    float q, u, w0, wp, wv, dvx, dvy;
    int pi, pj;
    for (i = collisions->listSize; i--;) {
        Collision_t *c = CollisionList_getElement(collisions, i);
        pi = c->element_i;
        pj = c->element_j;

        q = sqrt(c->r2) / h;
        u = 1 - q;
        w0 = C0 * u / (rho[pi] * rho[pj]);
        wp = w0 * Cp * (rho[pi] + rho[pj] - rho02) * u / q;
        wv = w0 * Cv;
        
        dvx = vx[pi] - vx[pj];
        dvy = vy[pi] - vy[pj];
        
        ax[pi] += wp * c->distance_x + wv * dvx;
        ay[pi] += wp * c->distance_y + wv * dvy;
        ax[pj] -= wp * c->distance_x + wv * dvx;
        ay[pj] -= wp * c->distance_y + wv * dvy;
    }

}

void normalizeMass() {

	int i,j;

    mass = 1;
    computeDensities();
    
    float rho2s = 0;
    float rhos  = 0;
    for (i = N; i--;) {
        rho2s += rho[i] * rho[i];
        rhos  += rho[i];
    }
    
    mass = rho0 * rhos / rho2s;
    // Constants for interaction term
    C0 = mass / (3.1416 * h4);
    C1 = 4 * mass / (3.1416 * h2);
    C2 = 4 * mass / (3.1416 * h8);
}

void leapfrogInit() {
	int i;

    for (i = N; i--;) {
        // Update half step velocity
        vhx[i] = vx[i] + ax[i] * dt2;
        vhy[i] = vy[i] + ay[i] * dt2;
        
        // Update velocity
        vx[i] = ax[i] * dt2;
        vy[i] = ay[i] * dt2;
        
        // Update position
        x[i] += vhx[i] * dt;
        y[i] += vhy[i] * dt;
    }
}

void leapfrogStep() {
	int i;

    for (i = N; i--;) {
        // Update half step velocity
        vhx[i] += ax[i] * dt;
        vhy[i] += ay[i] * dt;
        
        // Update velocity
        vx[i] = vhx[i] + ax[i] * dt2;
        vy[i] = vhy[i] + ay[i] * dt2;
        
        // Update position
        x[i] += vhx[i] * dt;
        y[i] += vhy[i] * dt;
        
        // Handle boundaries
        if (x[i] < edge1) {
            x[i] = edge1;// + random(0.0001, 0.0005);
            vx[i] *= -restitution;
            vhx[i] *= -restitution;
        } else if (x[i] > edge2) {
            x[i] = edge2;// - random(0.0001, 0.0005);
            vx[i] *= -restitution;
            vhx[i] *= -restitution;
        }
        
        /*if (y[i] < edge1) {
            y[i] = edge1 + random(0.0001, 0.0005);
            vy[i] *= -restitution;
            vhy[i] *= -restitution;
        } else */if (y[i] > edge3) {
            y[i] = edge3 ;//- random(0.0001, 0.0005);
            vy[i] *= -restitution;
            vhy[i] *= -restitution;
        }
    }
}

void update() {
    //computeAccelerations();
    updateParticles();
    leapfrogStep();
}

void initialiseSystem() {
    //initialiseArrays();
    
    particlesInMesh(0.05, height / width - 0.01);

    normalizeMass();
    computeAccelerations();
    leapfrogInit();
}

float random(float min, float max) {
    return min + (float)((float)(max - min) * ((float)rand() / ((float)RAND_MAX + 1.0)));
}

void gotoxy(int x, int y) {
    COORD coord;
  coord.X = x;
  coord.Y = y;
  SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), coord);
}

#include <unistd.h>

void clearScreen()
{
  system("cls");
}

char lastCoord[ancho][alto] = {' '};
int main()
{

    printf("Inicializando sistema");
    printf("\r\n");
    initialiseSystem();

    Cv = -40 * 3;
    dt = min(0.003, 0.05 / sqrt(-Cp * Cv));

    float v = 101 - 60;
    Cp =  15 * v;
    dt = min(0.003, 0.05 / sqrt(-Cp * Cv));

    printf("Iniciando programa");
    printf("\r\n");

    while(1){

        int count = 0;
        int MAX_COUNT = 30;
        while (count < MAX_COUNT) {
            update();
            count++;
        }

        //imprimir valores
        //printf("Particula[0]: x = %f y =%f", x[0], y[0]);
        //printf("\r\n");
        
        char coord[ancho][alto] = {' '};
        int i;
        for (i = N; i--;) {
            coord[(int)(x[i]*ancho)][(int)(y[i]*alto)] = 'O'; 
        }
        for (i = ancho; i--;) {
            for (j = alto; j--;) {
                if(coord[(int)i][(int)j] != lastCoord[(int)i][(int)j]){
                    gotoxy(i,j);
                    printf("%c", coord[(int)i][(int)j]);
                }
                
            }
        }
        Sleep(3);

        memcpy(lastCoord, coord, sizeof(lastCoord));

        //clearScreen();
//while(1){
    }

	return 0;
}

/**************************************
 *      Set-up system
***************************************/

//initialiseSystem();

/**************************************
 *      Main loop
***************************************/

// draw = function(){
  
//     var m = millis();
  
//     // Find maxRho
//     var maxRho = max.apply(null, rho);
  
//     // Draw particles
//     strokeWeight(diameter);
//     for (var i = N; i--;) {
//         stroke(lerpColor(minCol, maxCol, (rho[i] - rho0) / maxRho));
//         point(_scale * x[i], _scale * y[i]);
//     }
    
//     var count = 0;
//     var MAX_COUNT = 30;
//     while (count < MAX_COUNT && millis() - m < 40) {
//         update();
//         count++;
//     }
//     //println(count);
    
//     toolbar.draw();
    
// };

