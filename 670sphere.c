


/* A sphere has exactly one geometry uniform: its radius. */
#define sphUNIFDIM 1

/* An implementation of getIntersection for bodies that are spheres. */
void sphGetIntersection(
        int unifDim, const double unif[], const isoIsometry *isom, 
        const double p[3], const double d[3], double bound, 
        rayIntersection* inter) {
        double c[3];
    vecCopy(3, isom->translation, c);
    double pMinusC[3], dPMinusC, dD, rSq, disc, t;
    vecSubtract(3, p, c, pMinusC);
    dPMinusC = vecDot(3, d, pMinusC);
    dD = vecDot(3, d, d);
    rSq = unif[0] * unif[0];
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

/* An implementation of getTexCoordsAndNormal for bodies that are spheres. */
void sphGetTexCoordsAndNormal(
        int unifDim, const double unif[], const isoIsometry *isom, 
        const double p[3], const double d[3], const rayIntersection *inter, 
        double texCoords[2], double normal[3]) {
    double unitOutwardNormal[3], x[3], xLocal[3], rho, phi, theta, s, t, dLocal[3];
    vecScale(3, inter->t, d, dLocal);
    vecAdd(3, p, dLocal, x);
    vecSubtract(3, x, isom->translation, normal);
    vecUnit(3, normal, normal);
    isoUntransformPoint(isom, x, xLocal);
    vec3Rectangular(xLocal, &rho, &phi, &theta);
    s = theta/ (2 * M_PI);
    t = 1 - (phi/M_PI);
    texCoords[0] = s;
    texCoords[1] = t;
}


