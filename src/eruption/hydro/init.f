      subroutine init(n,hyd,alpha,cut,istart,time,encmg,eje,nadd,
     $                    dynamicalTime)

      include 'inclm1.f'

      character*20 hyd

      integer nadd, na

      real*8  encmg(*),  rho(mn), kdeg, solm, cv(mn), kap(mn), x1
      real*8  dynamicalTime

      parameter ( tau0 = 1d-2, solm = 1.989d33, pi = 3.1415926d0,
     $     as = 7.564d-15, r = 8.3147515d7, kdeg = 0.d0, x1=0.2 )

      real*8 encm(mn),dmass(mn)
      integer dummyInt

      common /mass/encm, dmass
      common /neutrn/ tmns
      alpha = 2.d0
      al1 = alpha+1.d0
      dummyInt = 3

      if(istart.le.1)then
         time = 0.d0
         open(15,file='EruptionFiles/InitForHydro.txt', status='old')
         read(15,*)n
         read(15,*)
         print *,n
         if(n.ge.mn)write(*,*)'you should define larger matrices'
         i = 3
         do j = 1, n
            read(15,*)kk,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,
     $     r12,r13,r14,r15,r16,r17,r18,r19,r20,r21,r22,r23,
     $     r24,r25,r26
            if(j.eq.2)tmns = r2
            if(j.ge.3)then
               dmass(i) = r3
               rho(i) = r5
               p(i) = r6
               temp(i) =r7
               rad(i) = r4
               x(i,1)=max(1d-40,r8)
               x(i,2)=r9
               x(i,3)=r10
               x(i,4)=r11
               x(i,5)=r12
               x(i,6)=r13
               x(i,7)=r14
               x(i,8)=r15
               x(i,9)=r16
               x(i,10)=r17
               x(i,11)=r18
               x(i,12)=r19
               x(i,13)=r20
               x(i,14)=r21
               x(i,15)=r22
               x(i,16)=r23
               x(i,17)=r24
               x(i,18)=r25
               x(i,19)=r26

               i = i + 1
            end if
         enddo

         temp(2) = temp(3)


 99      n = i-1
         print *,"n=",n
         do j = 3, n
            encm(j) = encm(j-1) + dmass(j)
         enddo
         write(*,*)"encm(3)=",encm(3)
         cut=tmns

         dmass(1) = dmass(4)
         dmass(2) = dmass(3)
         encm(2) = 0
         rho(1) = rho(4)
         rho(2) = rho(3)
         p(1) = p(4)
         p(2) = p(3)
         rad(2) = 1d-15

         do j = 1, n
            tau(j) = 1.d0/rho(j)
         enddo

         do j = 1, n
            encmg(j+1) = 0.5d0*(encm(j)+encm(j+1))
         enddo

         dynamicalTime = (3.141593*3.141593*
     $           rad(n)*rad(n)*rad(n)/4.d0/6.67e-8/encm(n))**(0.5d0)

         write(*,*)"encm(3)=",encm(3)
         call grav(n,encmg)

         call eoshelm_p(n,cv,temp,e,tau,p,x,grv,rad,eu,g,g1,cs,u,mn,
     $           nelem)

         e(1) = e(4)
         e(2) = e(3)



         rad(1) = -rad(3)

         do 20 j = 1, n
            j1 = max(1,j-1)
            if(j.ge.3)then
               ar(j) = 4*pi*(rad(j)*rad(j)+rad(j1)*rad(j)
     $              +rad(j1)*rad(j1))/al1
            endif
 20      continue
         call grav(n,encmg)

         call eoshelm(n,cv,temp,e,tau,p,x,grv,rad,eu,g,g1,cs,u,
     $        mn,nelem,time,dummyInt)

         print *,temp(3)

         e(1) = e(4)
         e(2) = e(3)

         u(1) = -u(4)
         u(2) = -u(3)
         ar(1) = -ar(4)
         ar(2) = -ar(3)



      else
         tmns = cut*solm
         write(*,*)' ns mass ',tmns
         open(15,file=hyd,status='old')
         do 40 it = 1, 100
            read(15,'(i5,5x,1pe12.4,7x,1pe12.4,12x,i6)', end=50)
     $           n, time, dt, istep
            write(*,'(i5,5x,1pe12.4,7x,1pe12.4,12x,i6)')
     $           n, time, dt, istep
            read(15,*)
            read(15,*)(k,rad(j),encm(j),dmass(j),rho(j),u(j),p(j),e(j)
     $           ,temp(j),dlr,ye(j),j=1,n)
            if(istep.eq.istart)go to 50
 40      continue
 50      continue
         close(15)

         open(65,file='abundance.d',status='old',form='formatted')
         read(65,*)na
         read(65,*)
         read(65,*)
         read(65,*)(idum,dum,x(j,3),dum,x(j,2),x(j,1),j=3,na)
         close(65)


         do k = 1, nelem
            do j = na+1, n
               x(j,k) = x(na,k)
            enddo
         enddo

         do j = 3, n
            x(j,1) = max(1d-10,x(j,1))
            x(j,2) = max(1d-10,x(j,2))

         enddo

         do j = 1, n
            tau(j) = 1.d0/rho(j)
         enddo
         do 60 j = 1, n
            j1 = max(j-1,1)
            if(j.ge.2)then
               ar(j) = 4*pi*(rad(j)*rad(j)+rad(j1)*rad(j)
     $              +rad(j1)*rad(j1))/al1
               rad(j+1) = (al1*dmass(j+1)*tau(j+1)/(4.d0*pi)+rad(j)
     $              **al1)**(1.d0/al1)
            end if
            encmg(j) = 0.5d0*(encm(max(2,j-1))+encm(j))
 60      continue
         call grav(n,encmg)
         encm(n+1) = encm(n)+dmass(n)
         ar(1) = -ar(4)
         ar(2) = -ar(3)
         pause
         !call eos(n,1,cv,kap)

      end if
      write(*,'(12htotal mass =,0pf9.4,4hmsol)')encm(n)/solm
      return
      end
