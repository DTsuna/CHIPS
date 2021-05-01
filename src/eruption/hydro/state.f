      subroutine state(n,dt,psl,psr,usl,usr)

      include 'inclm1.f'
      real*8 encm(mn), dmass(mn)
      real*8 psl(*), psr(*), usl(*), usr(*)
      real*8 pint(mn), uint(mn), tauint(mn), arint(mn),
     &       gint(mn), g1int(mn)
      real*8 taul(mn),pl(mn),ul(mn),arl(mn),gl(mn),g1l(mn)
      real*8 taur(mn),pr(mn),ur(mn),arr(mn),gr(mn),g1r(mn)
      real*8 taum(mn), pm(mn), um(mn), arm(mn), gm(mn), g1m(mn)
      real*8 taup(mn), pp(mn), up(mn), arp(mn), gp(mn), g1p(mn)
      real*8 dltau(mn), dlp(mn), dlu(mn), dlar(mn), dlg(mn), dlg1(mn)
      common /delta/ dltau, dlp, dlu, dlar, dlg, dlg1
      common /mass/encm, dmass
      common /intf/pint, uint, tauint, arint, gint, g1int
      common /avarg/ taum, pm, um, arm, gm, g1m, taup, pp,
     &               up, arp, gp, g1p
      common /intp/ taul, pl, ul, arl, gl, g1l
     &            , taur, pr, ur, arr, gr, g1r
      call intfce(n, dmass)
      call alar(n)
      call avrge(n,dt,dmass)

      do j = 2, n
         psl(j) = pp(j) 
     $        + dt*sqrt(g1p(j)*pp(j)/taup(j))*grv(j)
         psr(j) = pm(j) 
     $        - dt*sqrt(g1m(j)*pm(j)/taum(j))*grv(j+1)
c$$$         usl(j) = up(j)
c$$$         usr(j) = um(j)
         usl(j) = ur(j) + 2.0*(arp(j)*up(j)-ar(j)*ur(j))/(arp(j)+ar(j))
         usr(j) = ul(j+1) + 2.0*(arm(j)*um(j)-ar(j)*ul(j+1))
     $           /(arm(j)+ar(j))
      enddo
      psl(1) = psr(3)
      psr(1) = psl(3)
      usl(1) = -usr(3)
      usr(1) = -usl(3)
c$$$      usl(2) = 0
c$$$      usr(2) = 0
      return
      end
