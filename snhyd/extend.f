      subroutine extend(n,nadd)

      include 'inclm1.f'

      integer n, nadd

      real*8 encm(mn),dmass(mn)

      common /mass/encm, dmass

      if(nadd.eq.0)return

c     dmass(n) = dmass(n)/10.d0
      taun10 = log10(tau(n))+10.d0
c     taun10 = -20.d0
      pn10 = 2.d0
      taun11 = log10(tau(n))+3.d0
c     taun11 = -20.d0
      pn11 = log10(temp(n))
      n = n + nadd
      do 14 j = n-nadd+1, n
         dmass(j) = dmass(j-1)*1.03
         tau(j) = 10.d0**(((j-n+nadd-1)*taun10+(n-j)*taun11)
     $        /float(nadd-1))
         encm(j) = encm(j-1)+dmass(j)
         temp(j) = 10.d0**(((j-n+nadd-1)*pn10+(n-j)*pn11)
     $        /float(nadd-1))
         do k = 1, nelem
            x(j,k) = x(n-nadd,k)
         enddo
         p(j) = 8d7*temp(j)/tau(j) + 7.56d-15*temp(j)**4
         u(j) = max(0.d0,u(n-nadd)*1.1d0)
         e(j) = 1.5d0*8d7*temp(j) + 3.d0*7.56d-15*temp(j)**4*tau(j)
     $        +0.5d0*u(j)**2+grv(j)*rad(j)
c      write(*,'(i4,1p4e12.4)')j,p(j),tau(j),temp(j),e(j)
 14   continue

      return
      end
