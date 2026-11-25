// LB(lambda) from Proposition 7.4 and the dynamic programme in Appendix A of the paper. Written from the paper only.
// S(l) = {(a,b) in Z>=0^2 : (a+1)^2+(b+1)^2 <= l}
// LB = min over down-closed lattice polygons D in conv S of miss(D) + ceil((n(D)+phi(o(D)))/2)
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

static ll isqrt_ll(ll v){ if(v<0) return -1; ll r=(ll)sqrtl((long double)v); while(r*r>v) r--; while((r+1)*(r+1)<=v) r++; return r; }
static const int PHI[4]={0,1,2,1};

struct Shape { ll L; int W; vector<int> Hc; ll N; };
Shape shape(ll L){
  Shape s; s.L=L; s.N=0;
  for(int x=0;;x++){ ll rem=L-(ll)(x+1)*(x+1); if(rem<1) break; int h=(int)isqrt_ll(rem)-1; if(h<0) break; s.Hc.push_back(h); s.N+=h+1; }
  s.W=(int)s.Hc.size()-1; return s;
}

// ---------- Method 1: dynamic programme ----------
ll lb_dp(ll L){
  Shape s=shape(L); if(s.W<0) return 0;
  int W=s.W, H=s.Hc[0];
  vector<ll> pref(W+2,0); // pref[i] = sum_{x<i} H_x
  for(int x=0;x<=W;x++) pref[x+1]=pref[x]+s.Hc[x];
  vector<ll> tailPts(W+2,0); // sum_{x>=i}(H_x+1)
  for(int x=W;x>=0;x--) tailPts[x]=tailPts[x+1]+s.Hc[x]+1;
  const ll INF=LLONG_MAX/4;
  int stride=H+1;
  vector<ll> dp((size_t)(W+1)*stride*8, INF);
  auto id=[&](int x,int y,int p){ return ((size_t)x*stride+y)*8+p; };
  for(int x=0;x<=W;x++) for(int y=0;y<=s.Hc[x];y++) dp[id(x,y,0)] = 2*(tailPts[x+1] + s.Hc[x]-y);
  // directions (-a,b), a>=1, b>=0 primitive, ordered by decreasing slope b/a
  vector<pair<int,int>> dirs;
  for(int a=1;a<=W;a++) for(int b=0;b<=H;b++) if(__gcd(a,b)==1) dirs.push_back({a,b});
  sort(dirs.begin(),dirs.end(),[](const pair<int,int>&p,const pair<int,int>&q){ return (ll)p.second*q.first > (ll)q.second*p.first; });
  for(auto [a,b]:dirs){
    int cls = (b==0)?-1 : ((a&1) | ((b&1)<<1)); // (1,0)->1,(0,1)->2,(1,1)->3
    int bit = (cls<0)?0 : (cls==1?1: cls==2?2:4);
    ll inc = (cls<0)?0:1;
    ll fl = (ll)(a-1)*(b-1)/2 + b; // sum_{j=1..a} floor(b j / a)
    for(int cx=W; cx>=a; cx--){
      int nx=cx-a;
      for(int cy=0; cy<=s.Hc[cx]; cy++){
        int ny=cy+b; if(ny>s.Hc[nx]) break;
        ll missStep = (pref[cx]-pref[nx]) - (ll)a*cy - fl; // columns nx..cx-1
        ll add = 2*missStep + inc;
        for(int p=0;p<8;p++){ ll v=dp[id(cx,cy,p)]; if(v>=INF) continue; ll &t=dp[id(nx,ny,p^bit)]; if(v+add<t) t=v+add; }
      }
    }
  }
  ll best=INF;
  for(int y=0;y<=s.Hc[0];y++) for(int p=0;p<8;p++){ ll v=dp[id(0,y,p)]; if(v>=INF) continue; int o=__builtin_popcount(p); ll val=(v+PHI[o]+1)/2; best=min(best,val);}
  return best;
}

// ---------- Method 2: brute force over profiles ----------
struct BF { Shape s; ll best; vector<int> c; ll count=0; };
static ll cross(ll ox,ll oy,ll ax,ll ay,ll bx,ll by){ return (ax-ox)*(by-oy)-(ay-oy)*(bx-ox); }
void bf_eval(BF& B){
  int w=(int)B.c.size()-1;
  // points: (x,c_x) for x=0..w plus (w,0); upper hull from left to right
  vector<pair<ll,ll>> pts; for(int x=0;x<=w;x++) pts.push_back({x,B.c[x]}); pts.push_back({w,0});
  // sort by x then y desc for upper hull? build monotone chain upper hull
  sort(pts.begin(),pts.end());
  vector<pair<ll,ll>> h;
  for(auto&p:pts){ while(h.size()>=2 && cross(h[h.size()-2].first,h[h.size()-2].second,h.back().first,h.back().second,p.first,p.second)>=0) h.pop_back(); h.push_back(p);}
  // h: from (0,c0) to (w,0) upper hull; includes vertical drop at x=w if c_w>0 ((w,0) after (w,c_w))
  ll pts_in=0;
  for(int x=0;x<=w;x++){
    // height of hull at column x: find segment
    ll best=-1;
    for(size_t i=0;i+1<h.size();i++){ ll x1=h[i].first,y1=h[i].second,x2=h[i+1].first,y2=h[i+1].second; if(x1==x2){ if(x1==x) best=max(best,max(y1,y2)); continue;} if(x<x1||x>x2) continue; // y = y1 + (y2-y1)(x-x1)/(x2-x1), floor
      ll num=y1*(x2-x1)+(y2-y1)*(x-x1), den=x2-x1; ll f = num>=0? num/den : -((-num+den-1)/den); best=max(best,f);}
    if(x>=(int)B.s.Hc.size() || best>B.s.Hc[x]) { return; } // outside S (cannot happen for valid profile)
    pts_in+=best+1;
  }
  ll miss=B.s.N-pts_in;
  ll n=0; int par=0;
  for(size_t i=0;i+1<h.size();i++){ ll dx=h[i+1].first-h[i].first, dy=h[i].second-h[i+1].second; // edge traversed right-to-left: direction (-dx, dy)
    if(dx==0||dy==0) continue; ll g=__gcd(dx,dy); ll a=dx/g,b=dy/g; n+=g; if(g&1){ int cls=(a&1)|((b&1)<<1); int bit=cls==1?1:cls==2?2:4; par^=bit; } }
  int o=__builtin_popcount(par);
  ll val = miss + (n+PHI[o]+1)/2;
  B.count++;
  if(val<B.best) B.best=val;
}
void bf_rec(BF& B,int x,int prev){
  // choose c_x <= prev with (x,c_x) in S, then either stop or continue
  if(x>B.s.W) return;
  for(int cx=min(prev,B.s.Hc[x]); cx>=0; cx--){ B.c.push_back(cx); bf_eval(B); bf_rec(B,x+1,cx); B.c.pop_back(); }
}
ll lb_bf(ll L, ll* cnt=nullptr){ BF B; B.s=shape(L); B.best=LLONG_MAX; if(B.s.W<0) return 0; bf_rec(B,0,B.s.Hc[0]); if(cnt)*cnt=B.count; return B.best; }

int r_of(ll L){ int r=0; for(ll m=1;m*m<L;m++){ ll k2=L-m*m; ll k=isqrt_ll(k2); if(k>=1&&k*k==k2) r++; } return r; }

int main(int argc,char**argv){
  string mode=argc>1?argv[1]:"table";
  if(mode=="table"){ ll maxL=argc>2?atoll(argv[2]):500;
    for(ll L=2;L<=maxL;L++){ int r=r_of(L); if(!r) continue; ll lb=lb_dp(L); Shape s=shape(L); printf("%lld N=%lld r=%d LB=%lld slack=%lld\n",L,s.N,r,lb,lb-(r-1)); }
  } else if(mode=="bf"){ ll maxL=atoll(argv[2]);
    int bad=0, tot=0; for(ll L=2;L<=maxL;L++){ if(!r_of(L)) continue; ll c; ll a=lb_dp(L), b=lb_bf(L,&c); tot++; if(a!=b){bad++; printf("MISMATCH %lld dp=%lld bf=%lld\n",L,a,b);} }
    printf("bf check: %d eigenvalues, %d mismatches\n",tot,bad);
  } else { for(int i=1;i<argc;i++){ ll L=atoll(argv[i]); auto t=chrono::steady_clock::now(); ll lb=lb_dp(L); double sec=chrono::duration<double>(chrono::steady_clock::now()-t).count(); printf("%lld N=%lld r=%d LB=%lld slack=%lld (%.1fs)\n",L,shape(L).N,r_of(L),lb,lb-(r_of(L)-1),sec);} }
}
