      subroutine opac(n,kap,iphoto)
      include 'inclm1.f'
      include 'inclcnst.f'
      real*8 kap(*), muinv, ni, ne, z2
      real*8 tauo(mn)
      common/opdept/tauo
      tauo(n) = 0.d0
      do j = 3, n
         muinv = x(j,1)+0.25*x(j,2)
         ni = muinv/(tau(j)*h)
         ne = ye(j)/(tau(j)*h)
         z2 = x(j,1) + x(j,2)
         kap(j) = 0.4d0*ye(j) + x(j,2)*9d-3+ 1.7d-25*z2*ne*ni
     $        /sqrt(temp(j))**7*tau(j) 
     $        + 1d-2*(1.d0-x(j,1)-x(j,2))
c     $        + 2d-1*(1.d0-x(j,1)-x(j,2))
      enddo
      do j = n-1, 3, -1
         tauo(j) = tauo(j+1) + kap(j+1)/tau(j+1)*(rad(j+1)-rad(j))
      enddo

      do j = n-1, 3, -1
         iphoto = j - 1
         if(tauo(j).gt.0.66666667d0)return
      enddo
      return
      end
