// U(lambda) from Definitions 2.1 (cap) and 2.2 (U) of the paper. Written from the paper only.
#include <bits/stdc++.h>
using namespace std;
typedef long long ll; typedef pair<ll,ll> P;
static ll isq(ll v){ ll r=(ll)sqrtl((long double)v); while(r*r>v) r--; while((r+1)*(r+1)<=v) r++; return r; }
ll L; vector<int> Hc; int W,H;
bool inS(ll x,ll y){ return x>=0&&y>=0&&(x+1)*(x+1)+(y+1)*(y+1)<=L; }
static ll cr(P o,P a,P b){ return (a.first-o.first)*(b.second-o.second)-(a.second-o.second)*(b.first-o.first); }
vector<P> hull(vector<P> p){ sort(p.begin(),p.end()); p.erase(unique(p.begin(),p.end()),p.end()); if(p.size()<=2) return p;
  vector<P> h(2*p.size()); int k=0; for(auto&q:p){ while(k>=2&&cr(h[k-2],h[k-1],q)<=0) k--; h[k++]=q;} for(int i=p.size()-2,t=k+1;i>=0;i--){ while(k>=t&&cr(h[k-2],h[k-1],p[i])<=0) k--; h[k++]=p[i];} h.resize(k-1); return h; }
ll area2(const vector<P>& h){ if(h.size()<3) return 0; ll s=0; for(size_t i=0;i<h.size();i++){ auto a=h[i],b=h[(i+1)%h.size()]; s+=a.first*b.second-a.second*b.first;} return llabs(s); }
ll bdry(const vector<P>& h){ if(h.size()==1) return 0; if(h.size()==2) return 2*__gcd(llabs(h[0].first-h[1].first),llabs(h[0].second-h[1].second)); ll s=0; for(size_t i=0;i<h.size();i++){ auto a=h[i],b=h[(i+1)%h.size()]; s+=__gcd(llabs(a.first-b.first),llabs(a.second-b.second));} return s; }
vector<P> msum(const vector<P>& a,const vector<P>& b){ vector<P> p; for(auto&u:a) for(auto&v:b) p.push_back({u.first+v.first,u.second+v.second}); return hull(p); }
bool inside(const vector<P>& h,P q){ // q in convex polygon (ccw) or segment/point
  if(h.size()==1) return q==h[0];
  if(h.size()==2){ if(cr(h[0],h[1],q)!=0) return false; return min(h[0].first,h[1].first)<=q.first&&q.first<=max(h[0].first,h[1].first)&&min(h[0].second,h[1].second)<=q.second&&q.second<=max(h[0].second,h[1].second);}
  for(size_t i=0;i<h.size();i++) if(cr(h[i],h[(i+1)%h.size()],q)<0) return false; return true; }
int beta(int d,int u,int v,int s){ if(d==0) return 0; int sd=(d&1)?-s:s; return d-(v!=s)-(u!=sd); }
int cap(const vector<P>& A){ int dx=0,dy=0; for(auto&p:A){dx=max<ll>(dx,p.first);dy=max<ll>(dy,p.second);} bool top=true,rc=true;
  for(auto&p:A){ if(p.second==dy && (p.first&1)) top=false; if(p.first==dx && (p.second&1)) rc=false; }
  int best=INT_MIN; for(int sg=0;sg<16;sg++){ auto S=[&](int x,int y){ int i=((x>0)?1:0)*2+((y>0)?1:0); return (sg>>i&1)?1:-1; };
    for(int ls=0;ls<16;ls++){ int R=(ls&1)?1:-1,Lf=(ls&2)?1:-1,T=(ls&4)?1:-1,Bt=(ls&8)?1:-1; if(top&&Lf!=R) continue; if(rc&&Bt!=T) continue;
      int v=beta(dy,S(1,-1),S(1,1),R)+beta(dy,S(-1,-1),S(-1,1),Lf)+beta(dx,S(-1,1),S(1,1),T)+beta(dx,S(-1,-1),S(1,-1),Bt); best=max(best,v);} }
  return best; }
struct Fac{ vector<P> h; int I,eta; ll a2; };
vector<Fac> F; vector<P> Delta; int halfCapS;
map<pair<int,string>,int> memo;
int G(vector<int> e){ sort(e.rbegin(),e.rend()); int best=0,s=0; for(size_t j=0;j<e.size();j++){ s+=e[j]; best=max(best,min(halfCapS,s)-(int)(j+1)); } return best; }
bool fits(const vector<P>& h){ for(auto&p:h) if(!inS(p.first,p.second)) return false; return true; }
int dfs(int idx,const vector<P>& sum,vector<int>& etas){
  string key; for(auto&p:sum) key+=to_string(p.first)+","+to_string(p.second)+";"; key+="|"; {auto e=etas; sort(e.begin(),e.end()); for(int x:e) key+=to_string(x)+",";}
  auto it=memo.find({idx,key}); if(it!=memo.end()) return it->second;
  int best=G(etas);
  for(int j=idx;j<(int)F.size();j++){ auto ns=msum(sum,F[j].h); if(!fits(ns)) continue; ll mv2=area2(ns)-area2(sum)-F[j].a2; // doubled MV(sum,D_j)
    etas.push_back(F[j].eta); int v=F[j].I+1+(int)(mv2/2)+dfs(j,ns,etas); etas.pop_back(); best=max(best,v); }
  memo[{idx,key}]=best; return best; }
int main(int argc,char**argv){ int maxL=argc>1?atoi(argv[1]):50;
 for(L=2;L<=maxL;L++){ int r=0; for(ll m=1;m*m<L;m++){ ll k=isq(L-m*m); if(k>=1&&k*k==L-m*m) r++; } if(!r) continue;
  Hc.clear(); for(int x=0;(ll)(x+1)*(x+1)+1<=L;x++) Hc.push_back((int)isq(L-(ll)(x+1)*(x+1))-1); W=Hc.size()-1; H=Hc[0];
  vector<P> Spts; for(int x=0;x<=W;x++) for(int y=0;y<=Hc[x];y++) Spts.push_back({x,y}); ll N=Spts.size();
  ll kappa=1; for(ll m=1;m*m<L;m++) for(ll k=1;m*m+k*k<L;k++) kappa++;
  halfCapS=cap(Spts)/2;
  // enumerate down-closed polygons via profiles
  set<vector<P>> seen; F.clear();
  function<void(vector<int>&)> rec=[&](vector<int>& c){ int x=c.size(); if(x>0){ vector<P> p; for(int i=0;i<x;i++){p.push_back({i,c[i]});p.push_back({i,0});} auto h=hull(p);
      if(!seen.count(h)){ seen.insert(h); bool isPoint=h.size()==1; bool axisSeg=false; if(h.size()==2){ ll len=max(llabs(h[0].first-h[1].first),llabs(h[0].second-h[1].second)); axisSeg=(len>=2);}
        if(!isPoint&&!axisSeg){ Fac f; f.h=h; f.a2=area2(h); ll B=bdry(h); f.I=(int)((f.a2-B+2)/2); if(h.size()==2) f.I=0;
          vector<P> A; int wx=0,wy=0; for(auto&q:Spts) if(inside(h,q)){A.push_back(q); wx=max<ll>(wx,q.first); wy=max<ll>(wy,q.second);}
          int hh; auto eqv=[&](vector<P> g){ return hull(g)==h; };
          if(eqv({{0,0},{1,0},{0,1}})) hh=1; else if(wx==1&&wy>=2&&eqv({{0,0},{1,0},{0,wy}})) hh=wy; else if(wy==1&&wx>=2&&eqv({{0,0},{wx,0},{0,1}})) hh=wx; else hh=wx+wy;
          f.eta=min(hh,cap(A)/2); F.push_back(f);} } }
    if(x>W) return; int prev=x?c.back():H; for(int v=min(prev,Hc[x]);v>=0;v--){ c.push_back(v); rec(c); c.pop_back(); } };
  vector<int> c; rec(c);
  memo.clear(); vector<int> et; int U=1+dfs(0,{{0,0}},et);
  printf("%lld N=%lld kappa=%lld U=%d N-U=%lld |D|=%zu %s\n",L,N,kappa,U,N-U,F.size(), U<=kappa?(U==kappa?"U=kappa":"ok"):"VIOLATION");
 }
}
