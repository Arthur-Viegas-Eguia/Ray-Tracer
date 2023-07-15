


/* A plane has no geometry uniforms. */
#define plaUNIFDIM 0

/* An implementation of getIntersection for bodies that are planes. */
void plaGetIntersection(
        int unifDim, const double unif[], const void *data, const isoIsometry *isom, 
        const double p[3], const double d[3], double bound, 
        rayIntersection* inter) {
        double pLocal[3], dLocal[3], t;
        isoUntransformPoint(isom, p, pLocal);
        isoUnrotateDirection(isom, d, dLocal);
        if(dLocal[2] == 0){
                inter->t = rayNONE;
                return;
        }
        t = (-pLocal[2])/dLocal[2];
        if(t >= rayEPSILON && t <= bound){
                inter->t = t;
                return;
        }
        inter->t = rayNONE;
}

/* An implementation of getTexCoordsAndNormal for bodies that are planes. */
void plaGetTexCoordsAndNormal(
        int unifDim, const double unif[], const void *data, const isoIsometry *isom, 
        const double p[3], const double d[3], const rayIntersection *inter, 
        double texCoords[2], double normal[3]) {
        double pLocal[3], dLocal[3], localNormal[3] = {0.0, 0.0, 1.0}, td[3];
        isoRotateDirection(isom, localNormal, normal);
        isoUntransformPoint(isom, p, pLocal);
        isoUnrotateDirection(isom, d, dLocal);
        vecScale(3, inter->t, dLocal, td);
        texCoords[0] = pLocal[0] + td[0];
        texCoords[1] = pLocal[1] + td[1];
    
}


