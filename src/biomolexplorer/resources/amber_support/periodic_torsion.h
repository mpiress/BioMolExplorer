/* BioMolExplorer: exact periodic AMBER torsions with forward-mode derivatives.
 * E = k [1 + cos(n phi - phase)]. Original DOCK sources are not redistributed.
 * The signed dihedral convention matches NAB eff.c's improper-angle convention.
 */
#ifndef BIOMOL_PERIODIC_TORSION_H
#define BIOMOL_PERIODIC_TORSION_H
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define BME_N 12
typedef struct { double v; double d[BME_N]; } bme_dual;
static bme_dual bme_constant(double v) {
    bme_dual r; int i; r.v=v; for(i=0;i<BME_N;i++) r.d[i]=0.; return r;
}
static bme_dual bme_add(bme_dual a,bme_dual b) {
    bme_dual r; int i; r.v=a.v+b.v; for(i=0;i<BME_N;i++)r.d[i]=a.d[i]+b.d[i]; return r;
}
static bme_dual bme_sub(bme_dual a,bme_dual b) {
    bme_dual r; int i; r.v=a.v-b.v; for(i=0;i<BME_N;i++)r.d[i]=a.d[i]-b.d[i]; return r;
}
static bme_dual bme_mul(bme_dual a,bme_dual b) {
    bme_dual r; int i; r.v=a.v*b.v; for(i=0;i<BME_N;i++)r.d[i]=a.d[i]*b.v+a.v*b.d[i]; return r;
}
static bme_dual bme_div(bme_dual a,bme_dual b) {
    bme_dual r; int i; r.v=a.v/b.v; for(i=0;i<BME_N;i++)r.d[i]=(a.d[i]-r.v*b.d[i])/b.v; return r;
}
static bme_dual bme_sqrt(bme_dual a) {
    bme_dual r; int i; r.v=sqrt(a.v); for(i=0;i<BME_N;i++)r.d[i]=a.d[i]/(2.*r.v); return r;
}
static bme_dual bme_atan2(bme_dual y,bme_dual x) {
    bme_dual r; int i; double den=x.v*x.v+y.v*y.v; r.v=atan2(y.v,x.v);
    for(i=0;i<BME_N;i++)r.d[i]=(x.v*y.d[i]-y.v*x.d[i])/den; return r;
}
static bme_dual bme_cos(bme_dual a) {
    bme_dual r; int i; r.v=cos(a.v); for(i=0;i<BME_N;i++)r.d[i]=-sin(a.v)*a.d[i]; return r;
}
static bme_dual bme_dot(bme_dual a[3],bme_dual b[3]) {
    return bme_add(bme_add(bme_mul(a[0],b[0]),bme_mul(a[1],b[1])),bme_mul(a[2],b[2]));
}
static void bme_cross(bme_dual a[3],bme_dual b[3],bme_dual r[3]) {
    r[0]=bme_sub(bme_mul(a[1],b[2]),bme_mul(a[2],b[1]));
    r[1]=bme_sub(bme_mul(a[2],b[0]),bme_mul(a[0],b[2]));
    r[2]=bme_sub(bme_mul(a[0],b[1]),bme_mul(a[1],b[0]));
}
static int bme_extended_torsion(double n,double phase) {
    /* Native LEaP prints pi as 3.141594; retain its native fast path. */
    double p=remainder(phase,6.2831853071795864769); int k=(int)fabs(n);
    return fabs(fabs(n)-k)>1e-8 || (k!=1 && k!=2 && k!=3 && k!=4 && k!=6) ||
        (fabs(p)>1e-5 && fabs(fabs(p)-3.14159265358979323846)>1e-5);
}
static double bme_periodic_torsion(double coords[12],double n,double k,double phase,double gradient[12]) {
    bme_dual p[4][3],a[3],b[3],c[3],u[3],v[3],w[3],bb,uu,vv,phi,e; int i,j;
    n=fabs(n);
    if(!isfinite(n) || n<1. || fabs(n-round(n))>1e-8 || !isfinite(k) || !isfinite(phase)) {
        fprintf(stderr,"Invalid periodic torsion parameter\n"); exit(1);
    }
    for(i=0;i<4;i++)for(j=0;j<3;j++) { p[i][j]=bme_constant(coords[3*i+j]); p[i][j].d[3*i+j]=1.; }
    for(j=0;j<3;j++) { a[j]=bme_sub(p[1][j],p[0][j]); b[j]=bme_sub(p[2][j],p[1][j]); c[j]=bme_sub(p[3][j],p[2][j]); }
    bme_cross(a,b,u); bme_cross(b,c,v); bme_cross(u,v,w);
    bb=bme_dot(b,b); uu=bme_dot(u,u); vv=bme_dot(v,v);
    if(bb.v<1e-20 || uu.v<1e-20 || vv.v<1e-20) { fprintf(stderr,"Degenerate periodic torsion geometry\n"); exit(1); }
    phi=bme_atan2(bme_div(bme_dot(w,b),bme_sqrt(bb)),bme_dot(u,v));
    e=bme_mul(bme_constant(k),bme_add(bme_constant(1.),bme_cos(bme_sub(bme_mul(bme_constant(n),phi),bme_constant(phase)))));
    for(i=0;i<12;i++)gradient[i]=e.d[i];
    return e.v;
}
#endif
