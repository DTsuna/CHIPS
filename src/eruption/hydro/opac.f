      subroutine opac(n,kap,iphoto)
      include 'inclm1.f'
      include 'inclcnst.f'
      real*8 kap(*), muinv, ni, ne, z2
      real*8 k_m, k_h, k_e, k_k, k_grain
      real*8 tauo(mn)
      common/opdept/tauo
      tauo(n) = 0.d0


      do  j = 3,n

         kap(j) = 0.d0
      end do




!for https://www.astro.princeton.edu/~gk/A403/opac.pdf
      do j = 3, n
         k_m = 1.d-1*(1.d0 - x(j,1) - x(j,2) - x(j,3))
         k_h = 1.1d-25*sqrt(1.d0 - x(j,1) - x(j,2) - x(j,3))*
     $         sqrt(1.d0/tau(j))*
     $              (temp(j)**7.7d0)
         k_e = 2.d-1*(1.d0+x(j,1))/
     $     ((1.d0 + 2.7d11/(tau(j)*temp(j)*temp(j)))*
     $     (1.d0 + (temp(j)/4.5d8)**0.86))
         k_k = 4.d25*(1.d0+x(j,1))*(1.d0-x(j,1)-x(j,2)-x(j,3)+0.001)/
     $       (tau(j)*(temp(j)**3.5))


         kap(j) = k_m + 1.d0/(1.d0/k_h + 1.d0/(k_e + k_k))
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
