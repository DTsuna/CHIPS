      subroutine atmos(n,natmos)
      include 'inclm1.f'

      real*8 encm(mn),  rho(mn)
     $     , dmass(mn), M, R, scaleheight, tempn, muinv
      common /mass/encm, dmass

      muinv = x(n,1)+x(n,2)*0.25
      M= encm(n)
      R= rad(n)
      tempn = temp(n)
      scaleheight=1.565857d16*tempn*r**4*muinv/m
      rho(n) = 1.d0/tau(n)
      i = 0
      do j = n+1, n+natmos
         dmass(j) = dmass(j-1)*0.8
         encm(j) = encm(j-1) + dmass(j-1)
         rho(j) = rho(j-1)-dmass(j)/scaleheight
         tau(j) = 1.d0/rho(j)
         temp(j) = tempn
         p(j) = rho(j)*temp(j)*muinv*8.31d7
         rad(j) = (rad(j-1)**3+dmass(j)/(12.56637*rho(j)))**(1.d0/3.d0)
         if(rho(j).lt.0.d0)then
            i = j - n - 1
            go to 98
         end if
      enddo

 98   if(i.gt.0)natmos = i
      n = n + natmos

      return
      end

