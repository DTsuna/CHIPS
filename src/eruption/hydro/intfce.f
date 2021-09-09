      subroutine intfce(n, dxi)

      include 'inclm1.f'

      real*8 pi
      parameter ( pi = 3.1415926d0 )
cex    calculate the values of physical quantities a(j+1/2) at the
cex    boundary between the jth and the j+1st zones.
cex    xi(j) <==> xi(j+1/2)
      real*8 dxi(*)
      real*8 dltau(mn), dlp(mn), dlu(mn), dlar(mn), dlg(mn), dlg1(mn)
      real*8 pint(mn), uint(mn), tauint(mn), arint(mn),
     &       gint(mn), g1int(mn)

      common /intf/pint, uint, tauint, arint, gint, g1int
      common /delta/ dltau, dlp, dlu, dlar, dlg, dlg1
      call deltaa(n,dxi)
      sigma = 0.d0
      do 20 i = 2, n-2
         sigma =  dxi(i-1)+dxi(i)+dxi(i+1)+dxi(i+2)
         dxi1 = dxi(i+1) + dxi(i)
         dxi2 = dxi(i+2) + dxi(i+1)
         dxim = dxi(i-1) + dxi(i)
         dxi20 = 2.d0*dxi(i) + dxi(i+1)
         dxi21 = 2.d0*dxi(i+1) + dxi(i)
         dpj = p(i+1) - p(i)
         pint(i) = p(i)+dxi(i)/dxi1*dpj+1.d0/sigma*
     &             (2.d0*dxi(i+1)/dxi1*dpj*dxi(i)*(dxim/dxi20
     &              -dxi2/dxi21)
     &              -dxi(i)/dxi20*dxim*dlp(i+1)
     &              +dxi(i+1)/dxi21*dxi2*dlp(i))
         duj = u(i+1) - u(i)
         uint(i) = u(i)+dxi(i)/dxi1*duj+1.d0/sigma*
     &             (2.d0*dxi(i+1)/dxi1*duj*dxi(i)*(dxim/dxi20
     &              -dxi2/dxi21)
     &              -dxi(i)/dxi20*dxim*dlu(i+1)
     &              +dxi(i+1)/dxi21*dxi2*dlu(i))
         dtauj = tau(i+1) - tau(i)
         tauint(i) = tau(i)+dxi(i)/dxi1*dtauj+1.d0/sigma*
     &             (2.d0*dxi(i+1)/dxi1*dtauj*dxi(i)*(dxim/dxi20
     &              -dxi2/dxi21)
     &              -dxi(i)/dxi20*dxim*dltau(i+1)
     &              +dxi(i+1)/dxi21*dxi2*dltau(i))
         dgj = g(i+1) - g(i)
         gint(i) = g(i)+dxi(i)/dxi1*dgj+1.d0/sigma*
     &             (2.d0*dxi(i+1)/dxi1*dgj*dxi(i)*(dxim/dxi20
     &              -dxi2/dxi21)
     &              -dxi(i)/dxi20*dxim*dlg(i+1)
     &              +dxi(i+1)/dxi21*dxi2*dlg(i))
         dg1j = g1(i+1) - g1(i)
         g1int(i) = g1(i)+dxi(i)/dxi1*dg1j+1.d0/sigma*
     &             (2.d0*dxi(i+1)/dxi1*dg1j*dxi(i)*(dxim/dxi20
     &              -dxi2/dxi21)
     &              -dxi(i)/dxi20*dxim*dlg1(i+1)
     &              +dxi(i+1)/dxi21*dxi2*dlg1(i))
         arint(i) = 4*pi*rad(i)*rad(i)
 20   continue
      dxi2 = dxi(n)+dxi(n-1)
      pint(n-1) = (dxi(n)*p(n-1)+dxi(n-1)*p(n))/dxi2
      uint(n-1) = (dxi(n)*u(n-1)+dxi(n-1)*u(n))/dxi2
      tauint(n-1) = (dxi(n)*tau(n-1)+dxi(n-1)*tau(n))/dxi2
      arint(n-1) = (dxi(n)*ar(n-1)+dxi(n-1)*ar(n))/dxi2
      gint(n-1) = (dxi(n)*g(n-1)+dxi(n-1)*g(n))/dxi2
      g1int(n-1) = (dxi(n)*g1(n-1)+dxi(n-1)*g1(n))/dxi2
      dxi2n = 2.d0*dxi(n)+dxi(n-1)
      pint(n) = (dxi(n)*p(n-1)+dxi2n*p(n))/dxi2
      uint(n) = (dxi(n)*u(n-1)+dxi2n*u(n))/dxi2
      arint(n) = (dxi(n)*ar(n-1)+dxi2n*ar(n))/dxi2
      tauint(n) = (dxi(n)*tau(n-1)+dxi2n*tau(n))/dxi2
      gint(n) = (dxi(n)*g(n-1)+dxi2n*g(n))/dxi2
      g1int(n) = (dxi(n)*g1(n-1)+dxi2n*g1(n))/dxi2
      do 30 i = 1, 1
         pint(i) = pint(4-i)
         uint(i) = -uint(4-i)
         arint(i) = -arint(4-i)
         tauint(i) = tauint(4-i)
         gint(i) = gint(4-i)
         g1int(i) = g1int(4-i)
 30   continue
      return
      end
