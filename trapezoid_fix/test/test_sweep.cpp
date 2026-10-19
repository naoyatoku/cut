#define protected public
#define private public
#include "../trapezoid.h"
#undef protected
#undef private
#include <vector>
// tick-level reconstruction: returns total distance moved using the profile semantic (speed per tick, raw/accy)
static void show(const char*n,trapezoid&t,long long accy){
  printf("%s: vt=%lld(raw) ",n,t._vt.raw());
  for(int i=0;i<3;i++){auto&p=t._prof[i]; printf("| %s n=%d l=%ld am=%ld ",i==0?"ACC":i==1?"CONST":"DEC",p.n,p.l,p.l_amari);}
  long s=0; for(int i=0;i<3;i++) s+=t._prof[i].l_all();
  printf("| sum=%ld\n",s);
}
int main(){
  long long accy=1000;
  acc a; accurized<> mx(1000LL*accy, accy); a.set(mx,100,100);   // 1000 pulse/tick max, 100 ticks to max
  long long VT=1000LL*accy;
  show("trapezoid L=300000     ",*new trapezoid(300000,0,VT,a,accy),accy);
  show("boundary  L=100000     ",*new trapezoid(100000,0,VT,a,accy),accy);
  show("triangle  L=50000      ",*new trapezoid(50000,0,VT,a,accy),accy);
  show("tiny      L=10         ",*new trapezoid(10,0,VT,a,accy),accy);
  show("vc>0      L=300000 vc=200",*new trapezoid(300000,200*accy,VT,a,accy),accy);
  show("vc>0 short L=1000 vc=500 ",*new trapezoid(1000,500*accy,VT,a,accy),accy);
  show("vt<vc     L=300000 vc=800 vt=300",*new trapezoid(300000,800*accy,300*accy,a,accy),accy);
  // sweep: ACC-end speed vs vt, DEC-start speed vs vt; check speed-gap vs accel and amari bound
  int bad=0,tot=0,exc=0; long maxam=0;
  for(long long L=1; L<=400000; L=L*13/10+1) for(long long vc: {0LL,100LL*accy,600LL*accy}) for(long long vt: {50LL*accy,333LL*accy+7,1000LL*accy}){
    tot++;
    try{ trapezoid t(L,vc,vt,a,accy);
      long s=0; for(int i=0;i<3;i++) s+=t._prof[i].l_all();
      long am=t._prof[1].l_amari; if(am>maxam)maxam=am;
      long vtick=t._vt.raw()/accy;
      bool decOnly=(t._prof[0].n==0&&t._prof[1].n==0&&t._prof[2].n>0&&s>L); bool b=!decOnly&&((s!=L)||am<0||am>vtick+1||t._vt.raw()>std::max(vt,vc));
      if(b){ if(bad<10) printf("BAD L=%lld vc=%lld vt=%lld sum=%ld am=%ld vt_real=%lld\n",L,vc/accy,vt/accy,s,am,t._vt.raw()/accy); bad++; }
    }catch(std::exception&e){ if(exc++<5) printf("EXC L=%lld vc=%lld vt=%lld: %s\n",L,vc/accy,vt/accy,e.what()); }
  }
  printf("sweep tot=%d bad=%d exc=%d maxamari=%ld\n",tot,bad,exc,maxam);
}
