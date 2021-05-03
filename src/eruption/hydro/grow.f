      subroutine grow(n, finish, dt, time, encm)
      include 'inclm1.f'
      integer n, j
      real*8 dt, time, encm(*), amp(mn), ps(mn), us(mn), rho(mn) 
      logical initial, finish
      data initial/.true./
      common /riem/ ps, us

      if(initial)then
         open(25,file='EruptionFiles/grow.d',
     $            form='formatted',status='unknown')
         do j = 3, n
            amp(j) = 0.d0
         enddo
         initial = .false.
      end if

      do j = 3, n-1
         rho(j) = 1.d0/tau(j)
      enddo

      do j = 3, n-1
         if((rho(j-1)-rho(j+1))*(ps(j)-ps(j-1)).gt.0.d0)then
            gr = sqrt((rho(j-1)-rho(j+1))*(ps(j)-ps(j-1)))
     $           /(rho(j)*(rad(j)-rad(j-1)))
            amp(j) = amp(j) + gr*dt
         end if
      enddo

      write(25,'(1p3e15.7)')time,exp(amp(59)),exp(amp(168))
      if(finish)then
         open(45,file='EruptionFiles/amp_RT.d',
     $          status='unknown',form='formatted')
         write(45,'(i5,1p2e15.7)')
     $        (j,encm(j)/1.989d33,exp(amp(j)),j=3,n-1)
         close(45)
         close(25)
      end if

      return
      end
