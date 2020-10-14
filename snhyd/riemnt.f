c
cn    name:     r i e m n t
c
c---- solve the riemann problem by newton's method
      subroutine riemnt( n, ihyd, gl, gr, g1l, g1r, pl, pr, taul, taur,
     *     vell, velr, pn)
c$ use omp_lib
      include 'inclm1.f'
      include 'omp_lib.h'

      integer i, it, maxrmn
      real*8 epsrmn, undflw, gmin, gmax
      parameter ( maxrmn = 100,
     &     epsrmn = 1d-8, undflw = 1d-15,
     &     gmin = 1.d0, gmax = 1.66667d0 )
      integer n
      real*8 ga, psprev(mn), pa, va,
     *     wl(mn), wr(mn), wlp(mn), wrp(mn)
     &     ,gb(mn), g1b(mn), gsl(mn), gsr(mn)
      real*8 gl(*), gr(*), g1l(*),g1r(*), pl(*), pr(*), taul(*), taur(*)
     *     ,vell(*), velr(*), ps(mn), us(mn)

      common /riem/ ps, us

c---- initial guess for resolved pressure
      do 10 i = 1, n
         ga = 0.5*(gl(i)+gr(i))
         pa = 0.5d0*(pl(i)+pr(i))
         va = 0.5d0*(taul(i)+taur(i))
         ps(i) = max(pa,pa+(vell(i)-velr(i))*sqrt(ga*pa/va))
         gb(i) = 0.5d0*(gl(i)+gr(i))
         g1b(i) = 0.5d0*(g1l(i)+g1r(i))
 10   continue
      ps(n) = pn

!      call omp_set_num_threads(12)

c$OMP PARALLEL
c$OMP& PRIVATE(i,it,gslc,gsrc,psl,psr,alp,psla,arp,psra,ggl,
c$OMP&                 pslg,fgm,wlg,ggr,psrg,wrg,fdel,delps,
c$OMP&                 error,ga,pa,va)
c$OMP DO
      do 20 i = 2, n
!         if(ihyd.le.1000)then
!         write(*,*)'THREAD=',omp_get_thread_num()
!         end if
         do 30 it = 1, maxrmn
            gslc = gl(i)+2.d0*(1.d0-gb(i)/g1b(i))*(gb(i)-1.d0)
     &           *(ps(i)-pl(i))/(ps(i)+pl(i))
            gsl(i) = max(gmin, min(gslc, gmax))
            gsrc = gr(i)+2.d0*(1.d0-gb(i)/g1b(i))*(gb(i)-1.d0)
     &           *(ps(i)-pr(i))/(ps(i)+pr(i))
            gsr(i) = max(gmin, min(gsrc, gmax))
            psl = ps(i)/pl(i)
            psr = ps(i)/pr(i)
            if(psl.eq.1.d0)psl=psl + 2.d-15
            if(psr.eq.1.d0)psr=psr + 2.d-15
            alp = 0.5d0*(g1l(i)-1.d0)/g1l(i)
            psla = psl**alp
            arp = 0.5d0*(g1r(i)-1.d0)/g1r(i)
            psra = psr**arp
            if(psla.ge.1.d0)then
cexpl -------- for shock wave
               ggl = -sign(0.5d0,-abs(gsl(i)/gl(i)-1.d0))+0.5d0
               pslg = psl*ggl
               fgm = sqrt(abs((1.d0-pslg)/((gsl(i)-1.d0)/(gl(i)-1.d0)
     &              -pslg+undflw)))
               if(i.eq.499.and.fgm.eq.0.and.it.eq.1)write(*,*)"e1",pslg
               wlg = sqrt(abs(0.5*pl(i)/taul(i)*((gsl(i)+1.0)*psl
     &              +gsl(i)-1.d0)))
               wl(i) = wlg*fgm
               wlp(i) = 0.25*(gsl(i)+1.0)/pl(i)*fgm*sqrt(abs(
     &              0.5*pl(i)/taul(i)/((gsl(i)+1.0)*psl+gsl(i)-1.0)))
     &              +0.5d0*ggl*wlg/(fgm+undflw)*(gsl(i)-gl(i)+undflw)/
     &              (pl(i)*((1.d0-psl)**2+undflw)*(gl(i)-1.d0))
            else
cexpl -------- for rarefaction wave
               wl(i) = sqrt(pl(i)/taul(i)*g1l(i))*alp*
     &              (1.d0-psl+undflw)/(1.d0-psla+undflw)
               wlp(i) = alp*sqrt(g1l(i)/taul(i)*pl(i))/pl(i)
     &              *((1.d0-alp)*psla+alp*psla/psl-1.d0)
     &              /(1.d0-psla+undflw)**2
            endif
            if(psra.ge.1.d0)then
cexpl -------- for shock wave
               ggr = -sign(0.5d0,-abs(gsr(i)/gr(i)-1.d0))+0.5d0
               psrg = psr*ggr
               fgm = sqrt(abs((1.d0-psrg)/((gsr(i)-1.d0)/(gr(i)-1.d0)
     &              -psrg+undflw)))
               if(i.eq.499.and.fgm.eq.0.and.it.eq.1)write(*,*)"e3",psrg
               wrg = sqrt(abs(0.5*pr(i)/taur(i)*((gsr(i)+1.0)*psr
     &              +gsr(i)-1.0)))
               wr(i) = wrg*fgm
               wrp(i) = 0.25*(gsr(i)+1.0)/pr(i)*fgm*sqrt(abs(
     &              0.5*pr(i)/taur(i)/((gsr(i)+1.0)*psr+gsr(i)-1.0)))
     &              +0.5d0*wrg*ggr/fgm*(gsr(i)-gr(i)+undflw)
     &              /(pr(i)*((1.d0-psr)**2+undflw)*(gr(i)-1.d0))
            else
cexpl -------- for rarefaction wave
               wr(i) = sqrt(pr(i)/taur(i)*g1r(i))*arp*
     &              (1.d0-psr+undflw)/(1.d0-psra+undflw)
               wrp(i) = arp*sqrt(g1r(i)/taur(i)*pr(i))/pr(i)
     &              *((1.d0-arp)*psra+arp*psra/psr-1.d0)
     &              /(1.d0-psra+undflw)**2
            endif
            psprev(i) = ps(i)
            fdel = 1.d0-wr(i)*(pr(i)-pl(i)-wr(i)*(velr(i)-vell(i)))/
     *           (wl(i)+wr(i))**2*wlp(i)
     &           -wl(i)*(pl(i)-pr(i)-wl(i)*(velr(i)-vell(i)))/
     *           (wr(i)+wl(i))**2*wrp(i)
            delps = ((wl(i)*wr(i)*(vell(i) - velr(i))
     *           + wr(i)*pl(i) + wl(i)*pr(i))
     *           / (wl(i) + wr(i))-psprev(i))/fdel

            ps(i) = psprev(i)+max(-0.9d0*ps(i),delps)

            if(i.eq.499.and.it.eq.1)then
              write(*,*)"error ps",wl(i),wr(i),vell(i),velr(i),
     &                                 pl(i),pr(i)
            end if

            error = abs(delps/ps(i))
            if(i.eq.n)then
               ps(n) = pn
               error = 0.d0
            end if
            if(error.lt.epsrmn) then
               us(i) = (pl(i)-pr(i)+wl(i)*vell(i)+wr(i)*velr(i))
     $              /(wl(i)+wr(i))
               goto 20
            endif
 30      continue
         ga = 0.5*(gl(i)+gr(i))
         pa = 0.5d0*(pl(i)+pr(i))
         va = 0.5d0*(taul(i)+taur(i))
         if(error.gt.1d-1)then
            write(*,'((a),i10,i5,1p,4e10.3)') 'riemnt: ',ihyd,i,error,ga,pa
     $           ,va
            ps(i) = max(pa,pa+(vell(i)-velr(i))*sqrt(ga*pa/va))
            us(i)=0.5d0*(vell(i)+velr(i))
         end if
 20   continue
c$OMP END DO
c$OMP END PARALLEL

      if(ps(99).ne.ps(99).or.ps(100).ne.ps(100))then
        write(*,*)"ps error in riemnt"
      end if

      ps(2)=0.5d0*(pl(2)+pr(2))
      us(2) = 0.d0
      ps(1) = ps(3)
      us(1) = -us(3)

      return
      end
