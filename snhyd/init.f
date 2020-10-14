      subroutine init(n,hyd,alpha,cut,istart,time,encmg,eje,nadd)

      include 'inclm1.f'

      character*20 hyd

      integer nadd, na 

      real*8  encmg(*),  rho(mn), kdeg, solm, cv(mn), kap(mn), x1
c      real*8 r1, r2, r3, r4, r5, r6
      parameter ( tau0 = 1d-2, solm = 1.989d33, pi = 3.1415926d0,
     $     as = 7.564d-15, r = 8.3147515d7, kdeg = 0.d0, x1=0.2 )
c     &            as = 7.564d-15, r = 8.3147515d7, kdeg = 2.384d14 )
      real*8 encm(mn),dmass(mn)

      common /mass/encm, dmass
      common /neutrn/ tmns
      alpha = 2.d0
      al1 = alpha+1.d0

      if(istart.le.1)then
         time = 0.d0
         open(15,file='InitForHydro.txt', status='old')		!名前変更
c$$$         open(15,file='Heger/hyd.rezone',status='old')
         read(15,*)n
         read(15,*)							!読み込み
         print *,n
c$$$         read(15,*)
         if(n.ge.mn)write(*,*)'you should define larger matrices'
         i = 3
c$$$         j = 1
c$$$c         read(15,*)kk,r3,r1,r2,r4,r6,r5,r7,r8
c$$$         read(15,*)r1,r2,r3,r4,r5,r6
c$$$         if(r1.gt.1.d0)then
c$$$            tmns = cut*solm
c$$$         else
c$$$            tmns = cut
c$$$         end if
         do j = 1, n
c$$$            read(15,*)kk,r1,r2,r3,r4,r5,r6,r7
            read(15,*)kk,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,
     $     r12,r13,r14,r15,r16,r17,r18,r19,r20,r21,r22,r23,
     $     r24,r25,r26							!読み込み
c$$$            print *,"grid#=",kk
c$$$ These are for 15 Msun model by Heger              
c$$$            if(j.eq.226)tmns = r2
c$$$            if(j.ge.227)then
c$$$ These are for 20 Msun model by Heger              
            if(j.eq.2)tmns = r2
!            if(j.eq.3)temp(3)=r7					!ここだけは温度も読み込み
            if(j.ge.3)then
c$$$ These are for 25 Msun model by Heger              
c$$$            if(j.eq.172)tmns = r2
c$$$            if(j.ge.173)then
c$$$            if(j.eq.1)tmns = r2
               dmass(i) = r3
               rho(i) = r5
               p(i) = r6						!この辺も読み込みのために変更
               temp(i) =r7                  !温度も読み込んでいいはずだよね
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

c$$$               print *,i,dmass(i),rho(i),p(i),rad(i),x(i,1),
c$$$     $              x(i,2)
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
c$$$         read(15,'(1P6E16.9)')
c$$$     $        (encm(j),dmass(j),rad(j),rho(j),p(j),temp(j),j=3,n+2)
c$$$         close(15)
c$$$         n = n + 2
c         if(tmns.lt.10)tmns = tmns*solm
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

c$$$         open(65,file='abundance.d',status='old',form='formatted')
c$$$         read(65,*)na
c$$$         read(65,*)
c$$$         read(65,*)
c$$$         read(65,*)(idum,dum,x(j,3),dum,x(j,2),x(j,1),dum,j=3,na)
c$$$         close(65)

c$$$         open(65,file='../sndata/pre16mcmp.d',status='old',
c$$$     $        form='formatted')
c$$$         read(65,*)ni,nf
c$$$         read(65,*)
c$$$         read(65,*)(kk,(x(j-ni+3,k),k=1,5),j=ni,nf)
c$$$         read(65,*)
c$$$         read(65,*)(kk,(x(j-ni+3,k),k=6,9),j=ni,nf)
c$$$         do k = 2, nelem
c$$$            do j = 1, n
c$$$               x(j,1)=1.0
c$$$               x(j,k)=0.0
c$$$c$$$               x(j,k) = x(nf-ni+3,k)
c$$$            enddo
c$$$         enddo

         do j = 1, n
            encmg(j+1) = 0.5d0*(encm(j)+encm(j+1))
         enddo
 


         write(*,*)"encm(3)=",encm(3)
         call grav(n,encmg)

         call eoshelm_p(n,cv,temp,e,tau,p,x,grv,rad,eu,g,g1,cs,u,mn,
     $           nelem)

!         do j = 3, n
!            e(j) = p(j)*3.d0*tau(j)+grv(j)*rad(j)
!         enddo
         e(1) = e(4)
         e(2) = e(3)
!         call eos(n,2,cv,kap) !実はこれいらないのでは?
c$$$         call extend(n,nadd)

!	print *,temp(3) !コメント
!	print *,rad(4) !コメント
!	print *,rho(3)	!コメント
!	print *,(al1*dmass(3)*tau(3)/(4.d0*pi)+rad(2)**al1)
!    $           **(1.d0/al1)!コメント


!	print *,(al1*dmass(4)*tau(4)/(4.d0*pi)+rad(3)**al1)
!     $           **(1.d0/al1)!コメント
!	print *,rho(3) !コメント
         
         rad(1) = -rad(3)
c     rad(1) = (al1*dmass(1)/(4*pi*rho(1)))**(1/al1)
         
         do 20 j = 1, n
            j1 = max(1,j-1)
            if(j.ge.3)then
               ar(j) = 4*pi*(rad(j)*rad(j)+rad(j1)*rad(j)
     $              +rad(j1)*rad(j1))/al1
            endif
 20      continue
         call grav(n,encmg)
!         call eos(n,1,cv,kap)
         call eoshelm(n,cv,temp,e,tau,p,x,grv,rad,eu,g,g1,cs,u,
     $        mn,nelem,time)

         print *,temp(3) !コメント

!         do 30 j = 3, 8
!            e(j) = e(j)+eje/(6*dmass(j))
c            u(j) = sqrt(eje/(6*dmass(j)))*rad(j-1)/rad(9)
! 30      continue
         e(1) = e(4)
         e(2) = e(3)
c     ar(1) = 4*pi*rad(1)*rad(1)/3
         u(1) = -u(4)
         u(2) = -u(3)
         ar(1) = -ar(4)
         ar(2) = -ar(3)
!         call grav(n,encmg)
!         call eos(n,1,cv,kap)
!         call eoshelm(n,cv,temp,e,tau,p,x,grv,rad,eu,g,g1,cs,u,mn,nelem)


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
c     $           ,temp(j),j=1,n)
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

c$$$         do j = 144, na
c$$$            x(j,2) = max(x(j,2),0.8d0)
c$$$            x(j,1) = 1.d0-x(j,2)-0.02d0
c$$$         enddo
c$$$         do j = 3, na
c$$$            x(j,1) = 0.5d0*x(j,1)
c$$$            x(j,2) = x(j,2) + x(j,1)
c$$$            if(x(j,1).gt.1d-5)write(*,'(1p3e12.5)')encm(j),
c$$$     $           (x(j,k),k=1,2)
c$$$         enddo

         do k = 1, nelem
            do j = na+1, n
               x(j,k) = x(na,k)
            enddo
         enddo

         do j = 3, n
            x(j,1) = max(1d-10,x(j,1))
            x(j,2) = max(1d-10,x(j,2))
c$$$            if(encm(j)-encm(3).lt.0.37*1.989d33)x(j,3) = 1.d0
         enddo

         do j = 1, n
c$$$            if(encm(j).gt.3.3d0*solm-tmns)then
c$$$c            if(encm(j).gt.4.d0*solm-tmns)then
c$$$c$$$               x(j,1) = x1*(encm(n)-4.d0*solm+tmns)
c$$$c$$$     $              /(encm(n)-4.d0*solm+tmns)
c$$$               x(j,1) = x1*(encm(n)-3.3d0*solm+tmns)
c$$$     $              /(encm(n)-3.3d0*solm+tmns)
c$$$               x(j,2) = 1.d0-x(j,1)-0.02d0
c$$$            else if(encm(j).gt.1.8d0*solm-tmns)then
c$$$c            else if(encm(j).gt.2.d0*solm-tmns)then
c$$$               x(j,1) = 1d-10
c$$$               x(j,2) = 1.d0-1.5d0*solm/(3.3d0*solm-tmns)-x(j,1)
c$$$c               x(j,2) = 1.d0-2.d0*solm/(4.d0*solm-tmns)-x(j,1)
c$$$            else
c$$$               x(j,1) = 1d-10
c$$$               x(j,2) = 1.d0-1.5d0*solm/(3.3d0*solm-tmns)-x(j,1)
c$$$c               x(j,2) = 1.d0-2.d0*solm/(4.d0*solm-tmns)-x(j,1)
c$$$            end if
            tau(j) = 1.d0/rho(j)
         enddo
c$$$         natmos = nadd/2
c$$$         call atmos(n,natmos)
c$$$         nadd = nadd - natmos
c$$$         call extend(n,nadd)
         do 60 j = 1, n
            j1 = max(j-1,1)
c            encm(j) = encm(j)*solm
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
c     ar(1) = 4*pi*rad(1)*rad(1)/3
         ar(1) = -ar(4)
         ar(2) = -ar(3)
         call eos(n,1,cv,kap)
c$$$         do j = 3, 4
c$$$            e(j) = e(j)+eje/(2*dmass(j))
c$$$            temp(j) = sqrt(sqrt(e(j)/(as*tau(j))))
c$$$c            u(j) = sqrt(eje/(6*dmass(j)))*rad(j-1)/rad(9)
c$$$         enddo
c$$$         e(1) = e(4)
c$$$         e(2) = e(3)
      end if
      write(*,'(12htotal mass =,0pf9.4,4hmsol)')encm(n)/solm
      return
      end

