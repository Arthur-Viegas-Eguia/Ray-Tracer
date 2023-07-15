


#define rayNONE -1.0
#define rayEPSILON 0.000001
#define rayINFINITY 1000000.0

/* Feel free to read and write this data structure's members. */
typedef struct rayIntersection rayIntersection;
struct rayIntersection {
    double t;   // rayNONE or the first intersection time in [rayEPSILON, tEnd]
};


