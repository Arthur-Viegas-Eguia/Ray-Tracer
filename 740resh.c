


/* A resh is a ray-tracing mesh. It has no geometry uniforms outside the 
attached meshMesh. */
#define reshUNIFDIM 0

/* Given vectors a, b - a, and c - a describing a triangle, with the first three 
entries being XYZ. Given point x, with the first three entries being XYZ, such 
that (up to numerical precision) x - a = p (b - a) + q (c - a). Computes p and 
q. Returns 1 on success or 0 on failure. */
int reshGetPQ(
        const double a[], const double bMinusA[], const double cMinusA[], 
        const double x[], double pq[2]) {
    /* For the 3x2 matrix A with columns b - a, c - a, compute A^T A. */
    double aTA[2][2];
    aTA[0][0] = vecDot(3, bMinusA, bMinusA);
    aTA[0][1] = vecDot(3, bMinusA, cMinusA);
    aTA[1][0] = aTA[0][1];
    aTA[1][1] = vecDot(3, cMinusA, cMinusA);
    /* Compute the 2x2 matrix (A^T A)^-1 if possible. */
    double aTAInv[2][2];
    if (mat22Invert(aTA, aTAInv) == 0.0)
        return 0;
    /* Compute the 2x3 matrix (A^T A)^-1 A^T. */
    double aTAInvAT[2][3];
    for (int i = 0; i < 2; i += 1)
        for (int j = 0; j < 3; j += 1)
            aTAInvAT[i][j] = 
                aTAInv[i][0] * bMinusA[j] + aTAInv[i][1] * cMinusA[j];
    /* Then pq = (A^T A)^-1 A^T (x - a). */
    double xMinusA[3];
    vecSubtract(3, x, a, xMinusA);
    pq[0] = vecDot(3, aTAInvAT[0], xMinusA);
    pq[1] = vecDot(3, aTAInvAT[1], xMinusA);
    return 1;
}

/* An implementation of getIntersection for bodies that are reshes. Assumes that 
the data parameter points to an underlying meshMesh with attribute structure 
XYZSTNOP. */
void reshGetIntersection(
        int unifDim, const double unif[], const void *data, 
        const isoIsometry *isom, const double p[3], const double d[3], 
        double bound, rayIntersection* inter) {
    meshMesh *mesh = (meshMesh *)data;
    int *triangle;
	double *a, *b, *c, bMinusA[3], cMinusA[3], norm[3];
    double aMinusP[3], nDotAP, nDotD, pLocal[3], dLocal[3], td[3], x[3], pq[2], t;
    double tTemp;
    int tempIndex;
    isoUntransformPoint(isom, p, pLocal);
    isoUnrotateDirection(isom, d, dLocal);
    inter->t = rayNONE;
    inter->index = -1;
	for(int i = 0; i < mesh->triNum; i++){ 
		triangle = meshGetTrianglePointer(mesh, i);
		a = meshGetVertexPointer(mesh, triangle[0]);
		b = meshGetVertexPointer(mesh, triangle[1]);
		c = meshGetVertexPointer(mesh, triangle[2]);
        vecSubtract(3, b, a, bMinusA);
        vecSubtract(3, c, a, cMinusA);
        vec3Cross(bMinusA, cMinusA, norm);
        nDotD = vecDot(3, norm, dLocal);
        if(nDotD != 0){
            vecSubtract(3, a, pLocal, aMinusP);
            nDotAP = vecDot(3, norm, aMinusP);
            t= nDotAP/nDotD;
            vecScale(3, t, dLocal, td);
            vecAdd(3, pLocal, td, x);
            if(t >= rayEPSILON && t <= bound){
                reshGetPQ(a, bMinusA, cMinusA, x, pq);
                if((pq[0] >= 0 && pq[1] >= 0 && (pq[0] + pq[1]) <= 1)){
                    inter->t = t;
                    inter->index = i;
                    bound = t;
                }
            }
        }
	}
}

/* An implementation of getTexCoordsAndNormal for bodies that are reshes. 
Assumes that the data parameter points to an underlying meshMesh with attribute 
structure XYZSTNOP. */
void reshGetTexCoordsAndNormal(
        int unifDim, const double unif[], const void *data, 
        const isoIsometry *isom, const double p[3], const double d[3], 
        const rayIntersection *inter, double texCoords[2], double normal[3]) {
    meshMesh *mesh = (meshMesh *)data;
    int *triangle;
	double *a, *b, *c, bMinusA[mesh->attrDim], cMinusA[mesh->attrDim];
    double pBMinusA[mesh->attrDim], qCMinusA[mesh->attrDim], scaledSum[mesh->attrDim], result[mesh->attrDim];
    double td[3], x[3], pLocal[3], dLocal[3], pq[2], uNorm[3];
    isoUntransformPoint(isom, p, pLocal);
    isoUnrotateDirection(isom, d, dLocal);
    vecScale(3, inter->t, dLocal, td);
    vecAdd(3, pLocal, td, x);
    triangle = meshGetTrianglePointer(mesh, inter->index);
	a = meshGetVertexPointer(mesh, triangle[0]);
	b = meshGetVertexPointer(mesh, triangle[1]);
	c = meshGetVertexPointer(mesh, triangle[2]);
    vecSubtract(mesh->attrDim, b, a, bMinusA);
    vecSubtract(mesh->attrDim, c, a, cMinusA);
    reshGetPQ(a, bMinusA, cMinusA, x, pq);
    vecScale(mesh->attrDim, pq[0], bMinusA, pBMinusA);
    vecScale(mesh->attrDim, pq[1], cMinusA, qCMinusA);
    vecAdd(mesh->attrDim, pBMinusA, qCMinusA, scaledSum);
    vecAdd(mesh->attrDim, scaledSum, a, result);
    texCoords[0] = result[3];
    texCoords[1] = result[4];
    vecUnit(3, (&result[5]), uNorm);
    isoRotateDirection(isom, uNorm, normal);
}


