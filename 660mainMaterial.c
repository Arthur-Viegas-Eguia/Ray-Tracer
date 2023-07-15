/* On macOS, compile with...
    clang 640mainSpheres.c 040pixel.o -lglfw -framework OpenGL -framework Cocoa -framework IOKit
On Ubuntu, compile with...
    cc 640mainSpheres.c 040pixel.o -lglfw -lGL -lm -ldl
*/
#include <stdio.h>
#include <math.h>
#include <GLFW/glfw3.h>
#include "040pixel.h"

#include "650vector.c"
#include "150texture.c"
#include "280matrix.c"
#include "300isometry.c"
#include "300camera.c"
#include "660ray.c"

#define SCREENWIDTH 512
#define SCREENHEIGHT 512



/*** SPHERES ******************************************************************/

/* Given a sphere of radius r, centered at the origin in its local coordinates. 
Given a modeling isometry that places that sphere in the scene. Given a ray 
x(t) = p + t d in world coordinates. Outputs a rayIntersection, whose t member 
is the least t, in the interval [rayEPSILON, bound], when the ray intersects the 
sphere. If there is no such t, then the t member is instead rayNONE. */

void getIntersection(
        double r, const isoIsometry *isom, const double p[3], const double d[3], 
        double bound, rayIntersection* inter) {
    double c[3];
    vecCopy(3, isom->translation, c);
    double pMinusC[3], dPMinusC, dD, rSq, disc, t;
    vecSubtract(3, p, c, pMinusC);
    dPMinusC = vecDot(3, d, pMinusC);
    dD = vecDot(3, d, d);
    rSq = r * r;
    disc = dPMinusC * dPMinusC - dD * (vecDot(3, pMinusC, pMinusC) - rSq);
    if (disc <= 0) {
    inter->t = rayNONE;
    return;
    }
    disc = sqrt(disc);
    t = (-dPMinusC - disc) / dD;
    if (rayEPSILON <= t && t <= bound) {
    inter->t = t;
    return;
    }
    t = (-dPMinusC + disc) / dD;
    if (rayEPSILON <= t && t <= bound) {
    inter->t = t;
    return;
    }
    inter->t = rayNONE;
}


/* Given the sphere that just produced the given rayIntersection. Outputs the 
sphere's texture coordinates at the intersection point. Also outputs the 
sphere's unit outward-pointing normal vector there, in world coordinates. */
void getTexCoordsAndNormal(double r, const isoIsometry *isom, const double p[3], const double d[3], const rayIntersection* inter, double texCoords[2], double normal[3]) {
    double unitOutwardNormal[3], x[3], xLocal[3], rho, phi, theta, s, t, isoLocal[3], dLocal[3];
    vecScale(3, inter->t, d, dLocal);
    vecAdd(3, p, dLocal, x);
    isoUntransformPoint(isom, x, xLocal);
    vec3Rectangular(xLocal, &rho, &phi, &theta);
    s = theta/ (2 * M_PI);
    t = 1 - (phi/M_PI);
    texCoords[0] = s;
    texCoords[1] = t;
    isoUntransformPoint(isom, isom->translation, isoLocal);
    vecSubtract(3,xLocal, isoLocal,normal);
    vecUnit(3, unitOutwardNormal, unitOutwardNormal);
}


/*** ARTWORK ******************************************************************/

camCamera camera;
double cameraTarget[3] = {0.0, 0.0, 0.0};
double cameraRho = 10.0, cameraPhi = M_PI / 3.0, cameraTheta = M_PI / 3.0;

/* Four Spheres */
#define BODYNUM 4
#define unifNUM 1
#define texNUM 1
isoIsometry isoms[BODYNUM];
double radii[BODYNUM] = {1.0, 0.5, 0.5, 0.5};
double unif[1] = {1};
double cAmbient[3] = {0.75, 0.75, 0.75};
texTexture texture;
const texTexture *textures[texNUM] = {&texture};
rayMaterial material;

/* Based on the uniforms, textures, rayIntersection, and texture coordinates, 
outputs a material. */
void getMaterial(
        int unifDim, const double unif[], int texNum, const texTexture *tex[], 
        const rayIntersection *inter, const double texCoords[2], 
        rayMaterial *material){
    material->hasAmbient = 1;
    material->hasDiffuse = 0;
    material->hasSpecular = 0;
    material->hasMirror = 0;
    material->hasTransmission = 0;
    texSample(textures[0], texCoords[0], texCoords[1], material->cDiffuse);
}

int initializeArtwork(void) {
    if (texInitializeFile(&texture, "meme.jpg") != 0) {
        return 1;
	}
    camSetProjectionType(&camera, camPERSPECTIVE);
    camSetFrustum(
        &camera, M_PI / 6.0, cameraRho, 10.0, SCREENWIDTH, SCREENHEIGHT);
    camLookAt(&camera, cameraTarget, cameraRho, cameraPhi, cameraTheta);
    double rot[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    for (int k = 0; k < BODYNUM; k += 1)
        isoSetRotation(&(isoms[k]), rot);
    double transl[3] = {0.0, 0.0, 0.0};
    isoSetTranslation(&(isoms[0]), transl);
    vec3Set(1.0, 0.0, 0.0, transl);
    isoSetTranslation(&(isoms[1]), transl);
    vec3Set(0.0, 1.0, 0.0, transl);
    isoSetTranslation(&(isoms[2]), transl);
    vec3Set(0.0, 0.0, 1.0, transl);
    isoSetTranslation(&(isoms[3]), transl);
    return 0;
}

void finalizeArtwork(void) {
    texFinalize(&texture);
    return;
}



/*** RENDERING ****************************************************************/
/* Given a ray x(t) = p + t d. Finds the color where that ray hits the scene (or 
the background) and loads the color into the rgb parameter. */


void getSceneColor(const double p[3], const double d[3], double rgb[3]) {
    double bound = rayINFINITY, texCoords[2], normal[3], sample[3];
    rayIntersection intersection;
    int winningI = -1;
    rayIntersection bestIntersection;
    bestIntersection.t = rayNONE;
    for(int i = 0; i < BODYNUM; i++){
        getIntersection(radii[i], &isoms[i], p, d, bound, &intersection);
        if(intersection.t != rayNONE){
            bound = intersection.t;
            winningI = i;
            bestIntersection = intersection;
        } 
    }
    if(bestIntersection.t != rayNONE){
        getTexCoordsAndNormal(radii[winningI], &isoms[winningI], p, d, &bestIntersection, texCoords, normal);
        getMaterial(unifNUM, unif, texNUM, textures, &bestIntersection, texCoords, &material);
        rgb[0] = 0.0;
        rgb[1] = 0.0;
        rgb[2] = 0.0;
        if(material.hasAmbient){
            double modulatedAmbient[3];
            vecModulate(3, cAmbient, material.cDiffuse, modulatedAmbient);
            vecAdd(3, modulatedAmbient, rgb, rgb);
        }
    }
    else{
        rgb[0] = 1.0;
        rgb[1] = 1.0;
        rgb[2] = 0.0;
    }
}
void render(void) {
    /* Build a 4x4 matrix that (along with homogeneous division) takes screen 
    coordinates (x0, x1, 0, 1) to the corresponding world coordinates. */
    double inverseViewport[4][4], inverseProjection[4][4], pv[4][4], homogeneous[4][4], mat[4][4];
    mat44InverseViewport(SCREENWIDTH, SCREENHEIGHT, inverseViewport);
    if(camera.projectionType == camPERSPECTIVE){
        camGetInversePerspective(&camera, inverseProjection);
    }
    else{
        camGetInverseOrthographic(&camera, inverseProjection);
    }
    mat444Multiply(inverseProjection, inverseViewport, pv);
    isoGetHomogeneous(&camera.isometry, homogeneous);
    mat444Multiply(homogeneous, pv, mat);
    /* Declare p and maybe compute d. */
    double p[4], d[3];
    if(camera.projectionType == camORTHOGRAPHIC){
        double point[3] = {0.0, 0.0, -1.0};
        isoRotateDirection(&camera.isometry, point, d);
    }
    /* Each screen point is chosen to be on the near plane. */
    double screen[4] = {0.0, 0.0, 0.0, 1.0};
    for (int i = 0; i < SCREENWIDTH; i += 1) {
        screen[0] = i;
        for (int j = 0; j < SCREENHEIGHT; j += 1) {
            screen[1] = j;
            /* Compute p and maybe also d. */
            mat441Multiply(mat, screen, p);
            vecScale(4, 1/p[3], p, p);
            if(camera.projectionType == camPERSPECTIVE){
                vecSubtract(3,p, camera.isometry.translation, d);
            }
            /* Set the pixel to the color of that ray. */
            double rgb[3];
            getSceneColor(p, d, rgb);
            pixSetRGB(i, j, rgb[0], rgb[1], rgb[2]);
        }
    }
}

/*** USER INTERFACE ***********************************************************/

void handleKey(
        int key, int shiftIsDown, int controlIsDown, int altOptionIsDown, 
        int superCommandIsDown) {
    if (key == GLFW_KEY_I)
        cameraPhi -= 0.1;
    else if (key == GLFW_KEY_K)
        cameraPhi += 0.1;
    else if (key == GLFW_KEY_J)
        cameraTheta -= 0.1;
    else if (key == GLFW_KEY_L)
        cameraTheta += 0.1;
    else if (key == GLFW_KEY_U)
        cameraRho *= 1.1;
    else if (key == GLFW_KEY_O)
        cameraRho *= 0.9;
    else if (key == GLFW_KEY_P) {
        if (camera.projectionType == camORTHOGRAPHIC)
            camSetProjectionType(&camera, camPERSPECTIVE);
        else
            camSetProjectionType(&camera, camORTHOGRAPHIC);
    }
    camSetFrustum(
        &camera, M_PI / 6.0, cameraRho, 10.0, SCREENWIDTH, SCREENHEIGHT);
    camLookAt(&camera, cameraTarget, cameraRho, cameraPhi, cameraTheta);
}

void handleTimeStep(double oldTime, double newTime) {
    double rotAxis[3] = {1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};
    double rotMatrix[3][3];
    mat33AngleAxisRotation(newTime, rotAxis, rotMatrix);
    for (int k = 0; k < BODYNUM; k += 1)
        isoSetRotation(&(isoms[k]), rotMatrix);
    render();
}

int main(void) {
    if (pixInitialize(SCREENWIDTH, SCREENHEIGHT, "640mainSpheres") != 0)
        return 1;
    if (initializeArtwork() != 0) {
        pixFinalize();
        return 2;
    }
    pixSetKeyDownHandler(handleKey);
    pixSetKeyRepeatHandler(handleKey);
    pixSetTimeStepHandler(handleTimeStep);
    pixRun();
    finalizeArtwork();
    pixFinalize();
    return 0;
}


