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
#include "670body.c"
#include "670sphere.c"

#define SCREENWIDTH 512
#define SCREENHEIGHT 512




/*** ARTWORK ******************************************************************/

camCamera camera;
double cameraTarget[3] = {0.0, 0.0, 0.0};
double cameraRho = 10.0, cameraPhi = M_PI / 3.0, cameraTheta = M_PI / 3.0;

/* Four Spheres */
#define BODYNUM 4
#define unifNUM 1
#define texNUM 4
bodyBody bodies[BODYNUM];
double unif[1] = {1};
double cAmbient[3] = {0.75, 0.75, 0.75};
texTexture texture1, texture2, texture3, texture4;
texTexture *textures[texNUM] = {&texture1, &texture2, &texture3, &texture4};
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
    texSample(tex[0], texCoords[0], texCoords[1], material->cDiffuse);
}

int initializeArtwork(void) {
    if ((texInitializeFile(&texture1, "meme.jpg") != 0) || (texInitializeFile(&texture2, "test2.jpg") != 0) || (texInitializeFile(&texture3, "test3.png") != 0) || (texInitializeFile(&texture4, "orange.jpg") != 0)) {
        return 1;
	}
    if(bodyInitialize(&bodies[0], sphUNIFDIM, 0, 1, sphGetIntersection, sphGetTexCoordsAndNormal, getMaterial) != 0 
    || bodyInitialize(&bodies[1], sphUNIFDIM, 0, 1, sphGetIntersection, sphGetTexCoordsAndNormal, getMaterial) != 0 
    || bodyInitialize(&bodies[2], sphUNIFDIM, 0, 1, sphGetIntersection, sphGetTexCoordsAndNormal, getMaterial) != 0 
    || bodyInitialize(&bodies[3], sphUNIFDIM, 0, 1, sphGetIntersection, sphGetTexCoordsAndNormal, getMaterial) != 0){
        texFinalize(&texture1);
        texFinalize(&texture2);
        texFinalize(&texture3);
        texFinalize(&texture4);
        return 2;
    }

    //Setting the textures of the spheres
    for(int i = 0; i < BODYNUM; i++){
        bodySetTexture(&bodies[i], 0, textures[i]);
    } 

    //Set the radius of the bodies
    double radius[1] = {1.0};
    bodySetGeometryUniforms(&bodies[0], 0, radius,1);
    radius[0] = 0.5;
    bodySetGeometryUniforms(&bodies[1], 0, radius,1);
    bodySetGeometryUniforms(&bodies[2], 0, radius,1);
    bodySetGeometryUniforms(&bodies[3], 0, radius,1);

    //Configuring the camera
    camSetProjectionType(&camera, camPERSPECTIVE);
    camSetFrustum(
        &camera, M_PI / 6.0, cameraRho, 10.0, SCREENWIDTH, SCREENHEIGHT);
    camLookAt(&camera, cameraTarget, cameraRho, cameraPhi, cameraTheta);

    //Setting the translations of the spheres
    double transl[3] = {0.0, 0.0, 0.0};
    isoSetTranslation(&(bodies[0].isometry), transl);
    vec3Set(1.0, 0.0, 0.0, transl);
    isoSetTranslation(&(bodies[1].isometry), transl);
    vec3Set(0.0, 1.0, 0.0, transl);
    isoSetTranslation(&(bodies[2].isometry), transl);
    vec3Set(0.0, 0.0, 1.0, transl);
    isoSetTranslation(&(bodies[3].isometry), transl);

    return 0;
}

void finalizeArtwork(void) {
    texFinalize(&texture1);
    texFinalize(&texture2);
    texFinalize(&texture3);
    texFinalize(&texture4);
    for(int i = 0; i < BODYNUM; i++){
        bodyFinalize(&bodies[i]);
    }
    return;
}



/*** RENDERING ****************************************************************/
/* Given a ray x(t) = p + t d. Finds the color where that ray hits the scene (or 
the background) and loads the color into the rgb parameter. */


void getSceneColor(int bodyNum, const bodyBody bodies[], const double p[3], const double d[3], double rgb[3]) {
    double bound = rayINFINITY, texCoords[2], normal[3], sample[3];
    rayIntersection intersection;
    int winningI = -1;
    rayIntersection bestIntersection;
    bestIntersection.t = rayNONE;
    for(int i = 0; i < BODYNUM; i++){
        bodyGetIntersection(&bodies[i], p, d, bound, &intersection);
        if(intersection.t != rayNONE){
            bound = intersection.t;
            winningI = i;
            bestIntersection = intersection;
        } 
    }
    if(bestIntersection.t != rayNONE){
        bodyGetTexCoordsAndNormal(&bodies[winningI], p, d, &bestIntersection, texCoords, normal);
        bodyGetMaterial(&bodies[winningI], &bestIntersection, texCoords, &material);
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
            getSceneColor(BODYNUM, bodies, p, d, rgb);
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
        isoSetRotation(&(bodies[k].isometry), rotMatrix);
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


