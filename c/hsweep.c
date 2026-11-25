// H(lam): number of slanted primitive units on the upper chain of conv S(lam),
// S(lam) = {(a,b) in Z>=0^2 : (a+1)^2+(b+1)^2 <= lam}.
// r(lam) = #{(m,k) in Z>=1^2 : m^2+k^2 = lam}.
// Sweep: check ceil(H/2) >= r-1  (equivalently H >= 2r-3) at every lam <= X with r >= 2.
// Usage: ./hsweep X            (full sweep)
//        ./hsweep -l lam ...   (single levels; r by trial division)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
typedef long long ll;
static ll isq(ll n){ if(n<0) return -1; ll s=(ll)sqrtl((long double)n); while(s*s>n) s--; while((s+1)*(s+1)<=n) s++; return s; }
static ll gcdl(ll a, ll b){ if(a<0)a=-a; if(b<0)b=-b; while(b){ll t=a%b;a=b;b=t;} return a; }
// returns H; uses caller buffers
static ll H_of(ll lam, ll *hx, ll *hy){
  ll W=-1; int n=0;
  for(ll x=0;;x++){ ll rem=lam-(x+1)*(x+1); if(rem<1) break; ll c=isq(rem)-1; if(c<0) break; W=x;
    // monotone chain upper hull, points in increasing x
    while(n>=2){ ll ax=hx[n-1]-hx[n-2], ay=hy[n-1]-hy[n-2], bx=x-hx[n-2], by=c-hy[n-2];
      if(ax*by-ay*bx>=0) n--; else break; }   // remove if not strictly right turn
    hx[n]=x; hy[n]=c; n++; }
  if(W<0) return 0;
  ll H=0; for(int i=1;i<n;i++){ ll dx=hx[i]-hx[i-1], dy=hy[i]-hy[i-1]; if(dy!=0) H+=gcdl(dx,dy); }
  return H;
}
static int r_trial(ll lam){ int c=0; for(ll m=1;m*m<lam;m++){ ll k2=lam-m*m; ll k=isq(k2); if(k>=1&&k*k==k2) c++; } return c; }
int main(int argc,char**argv){
  if(argc>2 && !strcmp(argv[1],"-l")){
    ll *hx=malloc(sizeof(ll)*4000000), *hy=malloc(sizeof(ll)*4000000);
    for(int i=2;i<argc;i++){ ll lam=atoll(argv[i]); ll H=H_of(lam,hx,hy); int r=r_trial(lam);
      printf("%lld H=%lld r=%d ceilH2=%lld need=%d slack=%lld H/lam^(1/3)=%.4f\n",lam,H,r,(H+1)/2,r-1,(H+1)/2-(r-1),H/cbrt((double)lam)); }
    return 0; }
  ll X=atoll(argv[1]);
  uint8_t *cnt=calloc(X+1,1);
  for(ll m=1;m*m<X;m++) for(ll k=1;m*m+k*k<=X;k++) cnt[m*m+k*k]++;
  // per-decade stats
  int D=0; while(pow(10,D)<X) D++;
  double minratio[32]; ll minratio_lam[32]; ll minslack[32]; ll minslack_lam[32]; ll nlev[32]; ll fails=0;
  for(int d=0;d<32;d++){minratio[d]=1e9;minslack[d]=1LL<<60;nlev[d]=0;}
  #pragma omp parallel
  {
    ll *hx=malloc(sizeof(ll)*(isq(X)+10)), *hy=malloc(sizeof(ll)*(isq(X)+10));
    double mr[32]; ll mrl[32], ms[32], msl[32], nl[32], f=0;
    for(int d=0;d<32;d++){mr[d]=1e9;ms[d]=1LL<<60;nl[d]=0;}
    #pragma omp for schedule(dynamic,4096)
    for(ll lam=2;lam<=X;lam++){
      int r=cnt[lam]; if(r==0) continue;
      ll H=H_of(lam,hx,hy);
      int d=(int)floor(log10((double)lam)+1e-12);
      double ratio=H/cbrt((double)lam);
      if(ratio<mr[d]){mr[d]=ratio;mrl[d]=lam;}
      ll s=(H+1)/2-(r-1);
      nl[d]++;
      if(s<ms[d]||(s==ms[d]&&lam<msl[d])){ms[d]=s;msl[d]=lam;}
      if(s<0){ f++;
        #pragma omp critical
        printf("FAIL lam=%lld r=%d H=%lld\n",lam,r,H); }
    }
    #pragma omp critical
    { for(int d=0;d<32;d++){ if(mr[d]<minratio[d]){minratio[d]=mr[d];minratio_lam[d]=mrl[d];}
        if(ms[d]<minslack[d]||(ms[d]==minslack[d]&&msl[d]<minslack_lam[d])){minslack[d]=ms[d];minslack_lam[d]=msl[d];} nlev[d]+=nl[d]; }
      fails+=f; }
  }
  printf("X=%lld fails=%lld\n",X,fails);
  for(int d=0;d<=D;d++) if(nlev[d]) printf("[1e%d,1e%d): levels=%lld  min H/lam^(1/3)=%.4f at %lld   min slack ceil(H/2)-(r-1)=%lld at %lld\n",d,d+1,nlev[d],minratio[d],minratio_lam[d],minslack[d],minslack_lam[d]);
  return 0;
}
