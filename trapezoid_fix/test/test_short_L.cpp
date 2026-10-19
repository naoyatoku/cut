#define protected public
#define private public
#include "../trapezoid.h"
#undef protected
#undef private
int main(){
  for(long long accy: {1000LL,1LL}) for(int ms:{100,50}){
    acc a; accurized<> mx(1000LL*accy, accy); a.set(mx,ms,ms);
    long long firstOK=-1, lastBad=-1;
    for(long long L=1; L<=3000; L++){
      try{ trapezoid t(L,0,1000,a,accy); if(firstOK<0) firstOK=L; }catch(std::exception&e){ lastBad=L; }
    }
    printf("accy=%lld ms=%d a1u=%lld(raw/tick) : vc=0 vt=1000, L=1..3000 : last throwing L=%lld, first OK L=%lld\n",accy,ms,a.a1u.raw(),lastBad,firstOK);
  }
  // realistic: 1000 pulse/tick max, 100ms, vc=0 -> vt=1000, L=300000 
  acc a; accurized<> mx(1000LL*1000, 1000); a.set(mx,100,100);
  trapezoid t(300000,0,1000,a,1000);
  for(int i=0;i<3;i++){auto&p=t._prof[i]; printf("phase%d n=%d l=%ld amari=%ld vc=%lld vt=%lld\n",i,p.n,p.l,p.l_amari,p.vc.raw(),p.vt.raw());}
  // triangle: L=50000 (decel+accel = 2*(1000*~100/2)=100000 >50000)
  trapezoid u(50000,0,1000,a,1000);
  printf("tri vt=%lld/1000\n",u._vt.raw());
  for(int i=0;i<3;i++){auto&p=u._prof[i]; printf("phase%d n=%d l=%ld amari=%ld vc=%lld vt=%lld\n",i,p.n,p.l,p.l_amari,p.vc.raw(),p.vt.raw());}
}
