      subroutine eos(n,mode,dedlt,kap)
      
      include 'inclm1.f'
      include 'inclcnst.f'
      include 'inclion.f'

      integer mode
      real*8 tauo(mn), lum(mn), prad(mn)

      real*8 e0(mn), yet(mn), yer(mn), xiont(mn,3)
     $     , xionr(mn,3), dedlt(mn), kap(mn), dltemp(mn), ek(mn)
     $     , a11, a12, a21, a22, b1, b2, dpdlt(mn), p0(mn)
     $     , yion(mn,3), yiont(mn,3), yionr(mn,3), x1, x3
     $     , eion, muinv, error, dlntde, dedlnr

      real*8 encm(mn), dmass(mn)

      common /mass/encm, dmass
      common / lumi / lum
      common/opdept/tauo
c
c mode: 0--For given Temperature, 1--For given internal energy
c       2--For given Pressure.
c       
      do j = 4, n			!form 3 to 4
         rho(j) = 1.d0/tau(j)
         ek(j) = 0.5*u(j)*u(j)+grv(j)*rad(j)
         eu(j) = max(1d-5*(ek(j)-e(j)),e(j)-ek(j))
      enddo
      
      if(mode.eq.1.and.temp(3).eq.0.d0)then
         do j = 4, n			!from 3 to 4
            xheavy = 1.d0-x(j,1)-x(j,2)
            eion = eu(j) - Na*(x(j,1)*ih
     $           + 0.25*x(j,2)*(ihe1+ihe2))
            if(eion.le.0.d0)then
               temp(j) = eu(j)/(1.5*R*(x(j,1)+0.25*x(j,2)+xheavy/0.16))
            else
               temp(j) = max(sqrt(sqrt(eion/(arad*tau(j)))),
     $              eion/(1.5*R*(x(j,1)+0.25*x(j,2)+xheavy/0.16)))
            end if
            temp(j)  = max(t_min,temp(j))
         enddo
      else if(mode.eq.2)then
         do j = 4, n			!from 3 to 4
            temp(j) = min(sqrt(sqrt(3.d0*p(j)/arad)),p(j)/(R*rho(j)
     $           *(x(j,1) + 0.25*x(j,2)+xheavy/0.16)))
         enddo
      end if
      
      if(ye(n).eq.0.d0)then
         do j = 4, n			!from 3 to 4
c            ye(j) = max(1d-1,ye(j))
            ye(j) = 0.5d0
            xion(j,1) = x(j,1)
            xion(j,3) = 0.25*x(j,2)
         enddo

         do j = 4, n			!from 3 to 4
            if(temp(j).gt.4d3.and.temp(j).lt.7d3)then
               ye(j) = x(j,1)
            else if(temp(j).lt.4d3)then
               ye(j) = 1d-3*x(j,1)
            else if(temp(j).ge.7d3.and.temp(j).lt.1.8d4)then
               ye(j) = x(j,1) + 0.25*x(j,2)
            else if(temp(j).ge.1.8d4)then
               ye(j) = x(j,1) + 0.5*x(j,2)
            end if
         enddo
      end if
      
      do j = 4, n			!from 3 to 4
         ne(j) = ye(j)*Na*rho(j)
c         ne(j) =0
      enddo
      
      call saha(n)

      if(mode.eq.0)then

         do j = 4, n			!from 3 to 4
            x1 = x(j,1)-xion(j,1)
            x3 = 0.25*x(j,2)-xion(j,2)-xion(j,3)
            a11 = 2.d0*x3+xion(j,2)
            a12 = xion(j,3)-x3
            a21 = ye(j)-xion(j,2)+x1*xion(j,1)/x(j,1)
            a22 = -(ye(j)+2.d0*xion(j,3)+x1*xion(j,1)/x(j,1))
            b1 = (ihe1-ihe2)*x3/(kb*temp(j))
            b2 = (x1*xion(j,1)*(ih-ihe2)/x(j,1)
     $           -ye(j)*(ihe2+1.5*kb*temp(j)))/(kb*temp(j))
            det = a11*a22-a12*a21
            xiont(j,2) = (a22*b1-a12*b2)/det
            xiont(j,3) = (a11*b2-a21*b1)/det
            xiont(j,1) = x1*(ye(j)*(ih/(kb*temp(J))+1.5)
     $           -xion(j,2)*xiont(j,2)-2.d0*xion(j,3)*xiont(j,3))/
     $           (x1*ye(j)+ye(j)*xion(j,1)+x1*xion(j,1))
            yet(j) = xiont(j,1)*xion(j,1)+xiont(j,2)*xion(j,2)
     $           +2.d0*xiont(j,3)*xion(j,3)
c     write(*,'(1p6e12.4)')a11,a12,b1, a21, a22, b2
c$$$  write(*,'(1p5e12.4)')temp(j),det,xion(j,3),xiont(j,3),
c$$$  $        (xion(j+1,3)-xion(j,3))*(temp(j)+temp(j+1))
c$$$  $        /((temp(j+1)-temp(j))*(xion(j,3)+xion(j+1,3)))
         enddo
         
         b1 = 0.d0
         do j = 4, n			!from 3 to 4		
            x1 = x(j,1)-xion(j,1)
            x3 = 0.25*x(j,2)-xion(j,2)-xion(j,3)
            a11 = 2.d0*x3+xion(j,2)
            a12 = xion(j,3)-x3
            a21 = ye(j)-xion(j,2)+x1*xion(j,1)/x(j,1)
            a22 = -(ye(j)+2.d0*xion(j,3)+x1*xion(j,1)/x(j,1))
            b2 = ye(j)
            det = a11*a22-a12*a21
            xionr(j,2) = (a22*b1-a12*b2)/det
            xionr(j,3) = (a11*b2-a21*b1)/det
            yer(j) = xion(j,2)*xionr(j,2)+2*xion(j,3)*xionr(j,3)
            xionr(j,1) = -x1*yer(j)/(x(j,1)*ye(j))
c$$$  write(*,'(1p5e12.4)')rho(j),det,xion(j,1),xionr(j,1),
c$$$  $           (xion(j+1,1)-xion(j,1))*(rho(j)+rho(j+1))
c$$$  $           /((rho(j+1)-rho(j))*(xion(j,1)+xion(j+1,1)))
         enddo
         
      else 
c            write(*,'(1p6e12.4)')(TEMP(j),(xion(j,k),k=1,3),ye(j),oye(j)
c     $        ,j=3,n)
         do it = 1, 30
            
            error = 0.d0
            
            do j = 4, n			!from 3 to 4
               x1 = x(j,1)-xion(j,1)
               x3 = 0.25*x(j,2)-xion(j,2)-xion(j,3)
               a11 = 2.d0*x3+xion(j,2)
               a12 = xion(j,3)-x3
               a21 = ye(j)-xion(j,2)+x1*xion(j,1)/x(j,1)
               a22 = -(ye(j)+2.d0*xion(j,3)+x1*xion(j,1)/x(j,1))
               b1 = (ihe1-ihe2)*x3/(kb*temp(j))
               b2 = (x1*xion(j,1)*(ih-ihe2)/x(j,1)
     $              -ye(j)*(ihe2+1.5*kb*temp(j)))/(kb*temp(j))
               det = a11*a22-a12*a21
               xiont(j,2) = (a22*b1-a12*b2)/det
               xiont(j,3) = (a11*b2-a21*b1)/det
               xiont(j,1) = x1*(ye(j)*(ih/(kb*temp(J))+1.5)
     $              -xion(j,2)*xiont(j,2)-2.d0*xion(j,3)*xiont(j,3))/
     $              (x1*ye(j)+ye(j)*xion(j,1)+x1*xion(j,1))
               yet(j) = xiont(j,1)*xion(j,1)+xiont(j,2)*xion(j,2)
     $              +2.d0*xiont(j,3)*xion(j,3)
c               write(*,'(1p6e12.4)')a11,a12,b1, a21, a22, b2
            enddo

            do k = 1, 3
               do j = 4, n		!from 3 to 4
                  yion(j,k) = Na*xion(j,k)
               enddo
            enddo
            do k = 1, 3
               do j = 3, n
                  yiont(j,k) = xiont(j,k)*yion(j,k)
               enddo
            enddo
            
            if(mode.eq.1)then
               do j = 3, n
                  muinv = x(j,1)+0.25*x(j,2)+(1.d0-x(j,1)-x(j,2))/16.d0
     $                 +ye(j)
                  e0(j) = 1.5*muinv*R*temp(j)
     $                 + arad*temp(j)**4*tau(j)
     $                 + ih*yion(j,1) + ihe1*yion(j,2)
     $                 + (ihe1+ihe2)*yion(j,3)
               enddo
               
               do j = 3, n
                  dedlt(j) = (1.5*R*temp(j)*(ye(j)
     $                 +x(j,1)+0.25*x(j,2)+yet(j))+4*arad*temp(j)**4
     $                 *tau(j)+ih*yiont(j,1)
     $                 +ihe1*yiont(j,2)+(ihe1+ihe2)*yiont(j,3))
               enddo
               
               do j = 3, n
                  dltemp(j) = max(-0.2d0,(eu(j)-e0(j))/dedlt(j))
c                  yet(j) = max(-0.2d0*ye(j),yet(j)*dltemp(j))
               enddo
               
            else if(mode.eq.2)then
               do j = 3, n
                  muinv = x(j,1)+0.25*x(j,2)+(1.d0-x(j,1)-x(j,2))/16.d0
     $                 +ye(j)
                  p0(j) = rho(j)*(muinv*R*temp(j)
     $                 + ih*yion(j,1) + ihe1*yion(j,2)
     $                 + (ihe1+ihe2)*yion(j,3))
     $                 + arad*temp(j)**4/3.d0
c$$$                  write(*,'(i4,1p7e10.2)')j,p0(j),rho(j),ye(j),temp(j)
c$$$     $                 ,(xion(j,k),k=1,3)
               enddo
                  
               do j = 3, n
                  dpdlt(j) = rho(j)*(R*temp(j)*(ye(j)
     $                 +x(j,1)+0.25*x(j,2)+yet(j))
     $                 +ih*yiont(j,1)
     $                 +ihe1*yiont(j,2)+(ihe1+ihe2)*yiont(j,3))
     $                 +arad*temp(j)**4/0.75d0
               enddo
               
               do j = 3, n
                  dltemp(j) = max(-0.2d0,(p(j)-p0(j))/dpdlt(j))
c                  yet(j) = min(0.2d0*ye(j),yet(j)*dltemp(j))
               enddo
            end if
c$$$  write(*,'(1p4e12.4)')(temp(j),e0(j),dedlt(j),
c$$$  $        0.5*(temp(1)+temp(2))*(e0(2)-e0(1))/
c$$$  $        (temp(2)-temp(1)),j=1,n)
            
            do j = 3, n
               temp(j) = temp(j)*(1.d0+min(0.2d0,dltemp(j)))
            enddo
  
            do j = 3, n
               temp(j) = max(t_min,temp(j))
            enddo
c            do j = 3, n
c               write(*,'(i4,1p4e12.4)')j,temp(j),dltemp(j),ye(j),yet(j)
c               ye(j) = oye(j) + yet(j)
c               ne(j) = ye(j)*rho(j)*Na
c            enddo
            
            call saha(n)

c$$$            write(*,'(1p6e12.4)')(eu(j),e0(j),arad*temp(j)**4*tau(j)
c$$$     $           ,(xiont(j,k)*dltemp(j),k=1,3),j=3,n)

            do j = 3, n
               if(temp(j).eq.t_min)dltemp(j) = 0.d0
            enddo
            
            do j = 3, n
               error = max(error,abs(dltemp(j)))
            enddo
            
            if(error.lt.1d-10)go to 99
            
         enddo
         
 99      if(error.ge.1d-5)then
c$$$            print *,'no convergence in eos ',error,mode
            do j = 3, n
               if(abs(dltemp(j)).gt.1d-3)write(*,'(i4,1p7e12.4)')
     $              j,temp(j),dltemp(j),ye(j),eu(j),e0(j),p0(j),p(j)
            enddo
c     if(error.gt.0.1d0)stop
         endif
      end if
      if(mode.eq.0)then
         do k = 1, 3
            do j = 3, n
               yion(j,k) = Na*xion(j,k)
            enddo
         enddo
         do k = 1, 3
            do j = 3, n
               yiont(j,k) = xiont(j,k)*yion(j,k)
               yionr(j,k) = xionr(j,k)*yion(j,k)
            enddo
         enddo
      end if

c      if(mode.eq.2)then
         do j = 3, n
            temp(j)=max(1d2,temp(j))
            prad(j) = arad*temp(j)**4/3.d0
         enddo
c$$$      else
c$$$         prad(n) = 0.5d0*kap(n)*lum(n)*dmass(n)
c$$$     $        /(16.d0*pi**2*rad(n)**4*cl)
c$$$         do j = n-1, 3, -1
c$$$            kapj = 0.5d0*(kap(j)+kap(j+1))
c$$$            mgr = 0.5d0*(dmass(j)+dmass(j+1))
c$$$            prad(j) = prad(j+1) + kapj*lum(j)/(16.d0*pi**2*rad(j)**4*cl)
c$$$     $           *mgr
c$$$         enddo
c$$$      end if

      do j = 3, n
         muinv = ye(j)+x(j,1)+0.25*x(j,2)+(1.d0-x(j,1)-x(j,2))/16.d0
         p(j) = R*rho(j)*muinv*temp(j)
     $        +prad(j)
         eu(j) = 1.5*muinv*R*temp(j) + arad*temp(j)**4*tau(j)
     $        + ih*yion(j,1) + ihe1*yion(j,2) + (ihe1+ihe2)*yion(j,3)
         e(j) = ek(j) + eu(j)
         dedlt(j) = (1.5*R*temp(j)*(ye(j)
     $        +x(j,1)+0.25*x(j,2)+yet(j))+4*arad*temp(j)**4
     $        *tau(j)+ih*yiont(j,1)
     $        +ihe1*yiont(j,2)+(ihe1+ihe2)*yiont(j,3))
         dlntde = 1.d0/dedlt(j)
         dedlnr = ih*yionr(j,1) + ihe1*yionr(j,2)
     $        + (ihe1+ihe2)*yionr(j,3) - arad*temp(j)**4*tau(j)
c$$$         g(j) = 1+(tau(j)*p(j)-dedlnr)*dlntde
c$$$         g1(j) = R*rho(j)*temp(j)*(muinv+yer(j))/p(j) 
c$$$     $        + (tau(j)*p(j)-dedlnr)*((R*rho(j)*temp(j)
c$$$     $        *(muinv+yet(j))+4.d0*prad(j))/p(j))*dlntde
         g(j) = 4.d0/3.d0
         g1(j) = 4.d0/3.d0
         cs(j) = sqrt(g1(j)*p(j)/rho(j))
c$$$         print *,j,p(j),eu(j),e(j)
c$$$         pause
      enddo

      do j = 1, 2
         p(j) = p(5-j)
         e(j) = e(5-j)
         g(j) = g(5-j)
         g1(j) = g1(5-j)
         cs(j) = cs(5-j)
         temp(j) = temp(5-j)
c$$$         p(j+n-2) = p(n-1-j)
c$$$         e(j+n-2) = e(n-1-j)
c$$$         g(j+n-2) = g(n-1-j)
c$$$         g1(j+n-2) = g1(n-1-j)
c$$$         cs(j+n-2) = cs(n-1-j)
c$$$         temp(j+n-2) = temp(n-1-j)
      enddo
      
      do j = 3, n
         xion(j,1) = xion(j,1)/(x(j,1))
         xion(j,2) = xion(j,2)*4/(x(j,2))
         xion(j,3) = xion(j,3)*4/(x(j,2))
      enddo
      
      return
      end
