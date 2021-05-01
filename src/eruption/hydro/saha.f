      subroutine saha(n)
      
      include 'inclm1.f'
      include 'inclcnst.f'
      include 'inclion.f'

      real*8 t32, ex1, ex2, ex3

         
      do j = 3, n
         do it = 1, 100
            oye(j) = ye(j)
c$$$  if(ne(j).gt.1d26)then
c$$$  xion(j,1) = x(j,1)
c$$$  xion(j,2) = 0.d0
c$$$  xion(j,3) = 0.25*x(j,2)
c$$$  else
            t32 = temp(j)*sqrt(temp(J))
            ex1 = t32*exp(max(-100.d0,-ih/(kb*temp(j))))/ci
            ex2 = t32*exp(max(-100.d0,-ihe1/(kb*temp(j))))/ci
            ex3 = t32*exp(max(-100.d0,-ihe2/(kb*temp(j))))/ci
            xion(j,1) = x(j,1)*ex1/(ex1+2.d0*ne(j))
            xion(j,3) = 0.25*x(j,2)*ex2*ex3/(4.d0*ne(j)**2+
     $           2.d0*ne(j)*ex2 + ex2*ex3)
            xion(j,2) = 2.d0*xion(j,3)*ne(j)/ex3
c     end if
            ye(j) = sqrt(oye(j)
     $        *(xion(j,1) + xion(j,2) + 2.d0*xion(j,3)))
            ne(j) = ye(j)*rho(j)*Na 
            if(abs((1.d0-ye(j))/(1.d0-oye(j))-1.d0).le.1d-7)go to 99
         enddo
 99      continue
      enddo
      
c$$$      do j = 3, n
c$$$         t32 = temp(j)*sqrt(temp(J))
c$$$         ex1 = t32*exp(-ih/(kb*temp(j)))
c$$$         ex2 = t32*exp(-ihe1/(kb*temp(j)))
c$$$         ex3 = t32*exp(-ihe2/(kb*temp(j)))
c$$$         write(*,'(1p4e12.4)')temp(j),x(j,1)-xion(j,1)
c$$$     $        , 2.d0*ci*xion(j,1)*ne(j)/ex1,ye(j)
c$$$         write(*,'(1p3e12.4)')temp(j),0.25*x(j,2)-xion(j,2)-xion(j,3)
c$$$     $        , 2.d0*ci*xion(j,2)*ne(j)/ex2
c$$$         write(*,'(1p3e12.4)')temp(j),xion(j,2), 2.d0*ci*xion(j,3)*ne(j)
c$$$     $        /ex3
c$$$      enddo
c$$$      stop

      return
      end
