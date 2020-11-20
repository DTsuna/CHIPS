c+
c
cn    name:     e o s
c
      subroutine eos(n,mode,dedlt,kap)
      include 'inclm1.f'
      include 'inclcnst.f'
                              
      integer           n, mode
c
c mode: 0--For given Temperature, 1--For given internal energy
c       2--For given Pressure.
c       
      real*8 dedlt(*), am(nelem), muinv(mn), kap(*)
      integer   i
      data am/1.d0, 3.d0, 4.d0, 14.0, 12.d0, 16.d0, 20.d0, 24.d0, 28.d0,
     $     32.d0, 36.d0, 40.d0, 44.d0, 48.d0, 56.d0, 52.d0, 54.d0,
     $     56.d0, 56.d0/
                 
c       ------- evaluate eos

      data dkh,a/8.3147515e+07,7.56464e-15/
      data eps/1e-06/
      do i = 1, n			!from 1 to 4
         ye(i) = 0.5d0
         eu(i) = e(i)-0.5d0*u(i)**2-grv(i)*rad(i)
         temp(i) = (eu(i)/(tau(i)*a))**0.25d0
         muinv(i) = 0.d0
         do k = 1, nelem
            muinv(i) = muinv(i) + x(i,k)/am(k)
         enddo
      enddo

      eu(1)  = eu(4)
      eu(2)  = eu(3)
      temp(1)  = temp(4)
      temp(2)  = temp(3)
      if(mode.eq.1)then
         dtt = 1d30
         do it = 1,20
            if(abs(dtt).le.eps)goto 25
            dtt = -1.0
            do i = 3, n			!from 3 to 4
               muie = muinv(i) + ye(i)
               erad = a*temp(i)**4*tau(i)
               dt = (eu(i)-erad-1.5*dkh*muie*temp(i))/(4.0*erad/temp(i)
     $              +1.5*dkh*muie)
               temp(i) = temp(i)+dt
               dtt = max(dtt,abs(dt/temp(i)))
            enddo
         enddo
 25      if(abs(dtt).gt.eps)write(6,*)'iteration is not converged.',dtt
      end if
      if(mode.eq.2)then
         dtt = 1d30
         do it = 1,20
            if(abs(dtt).le.eps)goto 35
            dtt = -1.0
            do i = 3, n			!from 3 to 4
               muie = muinv(i) + ye(i)
               prad = a*temp(i)**4/3.d0
               dt = (p(i)-prad-dkh*muie*temp(i)/tau(i))/
     $              (4.0*prad/temp(i)+dkh*muie/tau(i))
               temp(i) = temp(i)+dt
               dtt = max(dtt,abs(dt/temp(i)))
            enddo
         enddo
 35      if(abs(dtt).gt.eps)write(6,*)'iteration is not converged.',dtt
      end if
      do i = 3, n			!from 3 to 4
         muie = muinv(i) + ye(i)
         p(i) = a*temp(i)**4/3.0+dkh*muie/tau(i)*temp(i)
         eu(i) = a*temp(i)**4*tau(i)+1.50*dkh*muie*temp(i)
         beta = dkh*muie*temp(i)/(p(i)*tau(i))
         g1(i) = (8.d0-3.d0*beta)/(6.d0-3.d0*beta)
         g(i) = p(i)*tau(i)/eu(i)+1.d0
         cs(i) = sqrt( g1(i)*p(i)*tau(i) )
         e(i) = eu(i)+0.5d0*u(i)**2+grv(i)*rad(i)
         dedlt(i) = 1.5d0*dkh*muie*temp(i)+4*temp(i)**4*tau(i)
!         write(*,'(i4,1p3e12.4)')i,e(i),temp(i),tau(i)  一旦消しておく
!         pause  なんのためにpauseしているのか?
      enddo

      e(1) = e(4)
      e(2) = e(3)
      p(1) = p(4)
      p(2) = p(3)
      g(1) = g(4)
      g(2) = g(3)
      g1(1) = g1(4)
      g1(2) = g1(3)
      cs(1) = cs(4)
      cs(2) = cs(3)
      return
            
      end
