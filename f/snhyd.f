      program snhyd
cexpl   one dimensional lagrangian hydrodynamic scheme.
c---  initial data is required.

      include 'inclm1.f'

      integer idev, ni, ipass, jj, kk
      common/nickel/ni
      common / neutrn / tmns
      real*8 msol, kh, eje, AU
      character*20 hyd, lightc
      character*128 filename 
      parameter ( msol = 1.989e33, kh = 8.3147515e7, AU=1.495978707e13 )
      real*8 dmass(mn), encm(mn), encmg(mn), tp(100)
      real*8 taum(mn), pm(mn), um(mn), arm(mn), gm(mn), g1m(mn)
     &     ,taup(mn), pp(mn), up(mn), arp(mn), gp(mn), g1p(mn)
      real*8 taul(mn), pl(mn), ul(mn), arl(mn), gl(mn), g1l(mn)
     &     ,taur(mn), pr(mn), ur(mn), arr(mn), gr(mn), g1r(mn)
      real*8 psl(mn), psr(mn), usl(mn), usr(mn)
      real*8 kap(mn), lum(mn), cv(mn)
      real*8 ps(mn), us(mn), pn
      real*8 tauo(mn)

      real*8 u_eos(mn)
      real*8 old_e(mn)
      real*8 old_eu(mn),integrate
      real*8 e_without
      real*8 boundr
      integer output_do
      real*8 when_out(99)
      integer output_init, dummyInt


      logical finish

      common/opdept/tauo
      common /riem/ ps, us
      common /mass/encm, dmass
      common /avarg/ taum, pm, um, arm, gm, g1m, taup, pp,
     &              up, arp, gp, g1p
      common /lumi/ lum
      common /intp/ taul, pl, ul, arl, gl, g1l
     &            , taur, pr, ur, arr, gr, g1r
      common /massn/ am(14)
      data iarrv, finish/0, .false./

      real*8  scaleDeposition
      logical scaleDepositionFlag
      real*8 scalingRate

      real*8 e_charge_tot, injection_time, time_to_cc, dynamicalTime
      integer ejectaCut
      integer fixedCell, innerCell
      real*8 fixedRad
      integer year
      ! Calculate only ejecta after dynamcal time if true.
      logical EjectaOnly 
      EjectaOnly = .true.
      dummyInt = 3
      year = 1

      output_do = 1
      ejectaCut = 0
      fixedCell = 0
      fixedRad = 1.d13

      pn = 0.d0
      am(1) = 1.d0
      am(2) = 4.d0
      do 5 k = 3, 14
         am(k) = 4.d0*k
 5       continue
      alpha = 0.d0
      open(17,file='snhydOutput/start.time',status='unknown')
      write(17,*)'just starting'
      close(17)
      open(18,file='f/para1.d',status='old')
      read(18,*)
      read(18,*)istart,ihydm, jw, iout, idev
      read(18,*)
      read(18,*)cut, dtcfac, eje, nadd
      read(18,*)
      read(18,*)ntp,(tp(k),k=1,ntp)
      read(18,*)
      read(18,*)hyd
      read(18,*)lightc
      write(*,*)' iout = ',iout
      close(18)

        
      open(21,file='f/eruptPara.d',status='old')
      read(21,*)
      read(21,*)time_to_cc, e_charge_tot, injection_time,
     $     scaleDeposition, scalingRate
      close(21)
      if(scaleDeposition.eq.0)scaleDepositionFlag = .false.
      if(scaleDeposition.eq.1)scaleDepositionFlag = .true.
      write(*,*)time_to_cc, e_charge_tot, injection_time,
     $          scaleDepositionFlag, scalingRate
      open(66,file='snhydOutput/passage@0.1AU.txt',form='formatted')
      write(66,*)' no. time radius mass density velocity pressure'
     $     ,' temperature'
c$$$      do k = 1, ntp
c$$$         tp(k) = tp(k)*8.64d4
c$$$      enddo
cexpl  construct the initial model
c      time = 0.d0
      call init(n, hyd, alpha, cut, istart, time, encmg, eje, nadd,
     $                    dynamicalTime)
      do jj = 1, 90
        when_out(jj) = (dynamicalTime*2/90.d0)*(jj-1)
      end do
      do jj = 1, 9
        when_out(jj+90) = (dynamicalTime*2/90.d0)*89.d0+
     $  ((time_to_cc -(dynamicalTime*2/90.d0)*89.d0)/9.d0)*jj
      end do 
      write(*,*)"********* OUTPUT TIME ***********"
      do jj = 1, 99
        write(*,*),jj,when_out(jj)
      end do

!      dynamicalTime = 1.d3


      boundr = (rad(3)+rad(4))/2.d0

      nna = n-nadd
      ipass = n
      write(*,*)' number of meshes ',n,' mass cut ',cut/1.989e33
      pn = 0

      call grav(n,encmg)

      l = max(1,1-istart)
      print *,e(3)
      call tote(n,nadd,e,dmass,rad,grv,u,te,tet)



      write(*,'(''total energy ='',1pe12.4,''erg, etherm =''
     $     ,e12.4,'' erg'')')te, tet

      if(scaleDepositionFlag)e_charge_tot = -1.d0*te*scalingRate
      write(*,*)time_to_cc, e_charge_tot, injection_time
c$$$      if(idev.ne.0)call view(nna,idev,time,rad,tau,p,u,ye,lum,temp)

!      call eos(n,1,cv,kap)
      call eoshelm(n,cv,temp,e,tau,p,x,
     %        grv,rad,eu,g,g1,cs,u,mn,nelem,time,dummyInt)

      call tote(n,nadd,e,dmass,rad,grv,u,te,tet)

      write(*,'(''total energy ='',1pe12.4,''erg, etherm =''
     $     ,e12.4,'' erg'')')te, tet

      write(*,66)(j,encm(j)/msol,rad(j),1.d0/tau(j),p(j),u(j),e(j)
     &            ,temp(j),j=max(1,jw-10),max(10,min(jw+10,n)))
 66    format(' no.',4x,'mr',4x,'rad',8x,'rho',6x,'p',5x,'vel'
     &       ,8x,'e',8x,'temp',/,(i4,1p7e9.2))

      call state(n,dt,psl,psr,usl,usr)


      call riemnt( n,ihyd,gl,gr,g1l,g1r,psl,psr,taup,taum,usl,
     *     usr,pn)
      open (11,file='snhydOutput/hyd.d',
     $               status='unknown',form='formatted')
      WRITe(11,'(''total energy ='',1pe12.4,''erg, total mass =''
     $     ,e12.4,'' Msun'')')te, encm(n)/1.989e33
      close(11)
      open (12,file='snhydOutput/lightc.d',
     $               status='unknown',form='formatted')
      dt = 0.d0
      kp = 1
      do 9 kpp = 1, ntp
         if(tp(kpp).gt.time)then
            kp = kpp
            go to 95
         end if
 9    continue
cexpl  start the hydrodynamical calculation
 95   write(*,*)'calculation starts here'
      ihyd = istart
c$$$      do 10 ihyd = istart, ihydm
      do while(time.le.tp(ntp))
      print *,time,ihyd,dt
      write(*,*)"fixedCell=",fixedCell


      if(EjectaOnly)then
        if(ejectaCut.eq.0)then
          if(time.gt.dynamicalTime*2.d0)then
            fixedCell = 5
            do jj = 6, n
              if(rad(jj).le.fixedRad)then
                fixedCell = jj
              end if
            end do
            boundr = rad(fixedCell)
            ejectaCut = 1
            write(*,*)"fixedcell=",fixedCell,"at",rad(fixedCell)
!            pause
          end if
        end if
      end if


      if(output_do.le.99)then
        if(time.gt.when_out(output_do))then
           write (filename, '("snhydOutput/result", i2.2, ".txt")')
     $         output_do
           open(91, file=filename,status='unknown',form='formatted')

           if(ejectaCut.eq.0)then
             write(91,93)n,time,te,ihyd,(j,rad(j),encm(j),dmass(j),
     $       1./tau(j), u(j), p(j), e(j), temp(j), lum(j)*1d-40, ye(j),
     $       j= 3, n)
           end if
           if(ejectaCut.eq.1)then
             write(91,93)n,time,te,ihyd,(j,rad(j),encm(j),dmass(j),
     $       1./tau(j), u(j), p(j), e(j), temp(j), lum(j)*1d-40, ye(j),
     $       j= fixedCell, n)
           end if
 93     format(i5,' time', 1pe12.4,' sec',' te',e12.4,' erg  istep ',i8
     &  ,/,' no.',5x,'rad',10x,'encm',13x,'dm',12x,'rho',14x,'v'
     &         ,14x,'p',14x,'e',14x,'t',/,(i5,1p10e15.7))
           close(91)
           output_do = output_do + 1
        end if
      end if

      innerCell = 3
        if(ejectaCut.eq.1)then
        innerCell = fixedCell
      end if 
      call cournt( n, dtcfac, time, dtc ,innerCell)
      dt = min(tp(kp)-time,dtc)
      if(dt.lt.dtc)kp = kp+1
      if(dt.le.0.d0)then
         write(*,*)dtc,dt,time,kp
         stop' due to negative time step'
      end if

      if(dt.gt.1000.d0)dt = 1000.d0

      call state(n,dt,psl,psr,usl,usr)


      call riemnt( n,ihyd,gl,gr,g1l,g1r,psl,psr,taup,taum,usl,
     *   usr,pn)

      write(*,*)"before advanc"
      call tote(n,nadd,e,dmass,rad,grv,u,te,tet)
      write(*,*)"e_tot eu_tot",te,tet

      innerCell = 3
      if(ejectaCut.eq.1)then
        innerCell = fixedCell
      end if
      call advanc(n,alpha,nadd,dt,dmass,encmg,time,boundr,
     $                e_charge_tot,injection_time,innerCell)

c      call grow(n, finish, dt, time, encmg)
      time = time + dt
      write(*,*)"timetocc=",time_to_cc


      do jj = 3,n
        u_eos(jj)= (us(jj - 1)+us(jj))/2.d0
      end do

      call tote(n,nadd,e,dmass,rad,grv,u,te,tet)


      innerCell = 3
      if(ejectaCut.eq.1)then
        innerCell = fixedCell
      end if      
      call eoshelm(n,cv,temp,e,tau,p,x,grv,rad,eu,g,g1,cs,u,mn,
     $   nelem,time,innerCell)
      call tote(n,nadd,e,dmass,rad,grv,u,te,tet)



      if(ejectaCut.eq.0)then
        call opac(n, kap,iphoto)
        call radtra(n,ihyd,dt,time,cv,kap)
      end if


      call tote(n,nadd,e,dmass,rad,grv,u,te,tet)

      if(time-dt.gt.year*86400.d0*365.25d0)then
        write (filename, '("snhydOutput/intermediate", i2.2, "yr.txt")')
     $         year
        year = year + 1
        output_init = 3
        if(ejectaCut.eq.1)then
          output_init = fixedCell
        end if
        open(98, file=filename,status='unknown',form='formatted')
        write(98,*)"j EnclosedM[g] Rad[cm] Vel[cm/s] Den[g/cc] X_H X_H
     $P[erg/cc] t=",time
        do jj = output_init, n
           write(98,'(i0, e18.10, e18.10, e18.10,
     $                   e18.10, e18.10, e18.10, e18.10)'),jj,
     $encm(jj)-encm(output_init-1),rad(jj),u(jj),
     $1.d0/tau(jj),x(jj,1),x(jj,3),p(jj)
        end do
        close(98)
      end if

      if(time-dt.gt.time_to_cc)then
        output_init = 3
        if(ejectaCut.eq.1)then
          output_init = fixedCell
        end if
        open(98,file='snhydOutput/atCCSN.txt',status='unknown'
     $               ,form='formatted')
        write(98,*)"j EnclosedM[g] Rad[cm] Vel[cm/s] Den[g/cc] X_H X_He
     $ P[erg/cc]"
        do jj = output_init, n 
           write(98,'(i0, e18.10, e18.10, e18.10,
     $                   e18.10, e18.10, e18.10, e18.10)'),jj,
     $encm(jj)-encm(output_init-1),rad(jj),u(jj),
     $1.d0/tau(jj),x(jj,1),x(jj,3),p(jj)
        end do
        close(98)
        finish = .true.
        go to 99
      end if

      if(ihyd.ge.0)then
!      if(ihyd/iout*iout.eq.ihyd)then
c$$$         jw = iphoto
         call tote(n,nadd,e,dmass,rad,grv,u,te,tet)
         write(*,'(''total energy ='',1pe12.4,''erg, etherm =''
     $        ,e12.4,'' erg'')')te, tet
         write(*,
     $        '(8h time = ,1pe12.4,4h sec,6h dt = ,e12.4,4h sec,
     $        6h ihyd ,i8)')
     $        time, dt, ihyd
         if(ejectaCut.eq.0)then
           write(*,'(5h no.  ,8h  mr    ,9hradius   9hdensity  ,
     $        9hpressure ,9hvelocity ,9h     e   ,9h temp    ,
     $        9h    u   ,/,(i5,1p8e9.2))')
     $        (j,encm(j)/msol,rad(j),1.d0/tau(j),ps(j),us(j),e(j)
     $        ,temp(j),u(j),j=max(1,jw-10),max(10,min(jw+10,n)))
         end if
         if(ejectaCut.eq.1)then
           write(*,'(5h no.  ,8h  mr    ,9hradius   9hdensity  ,
     $        9hpressure ,9hvelocity ,9h     e   ,9h temp    ,
     $        9h    u   ,/,(i5,1p8e9.2))')
     $        (j,encm(j)/msol,rad(j),1.d0/tau(j),ps(j),us(j),e(j)
     $        ,temp(j),u(j),
     $        j=max(fixedCell,jw-10),max(20+fixedCell,min(jw+10,n)))
         end if

c$$$         if(idev.ne.0)call view(nna,idev,time,rad,tau,p,u,ye,lum,temp)
      endif

      kp1 = max(kp-1,1)
      if(time.eq.tp(kp1)) then
         open (11,file='snhydOutput/hyd.d',
     $           access='append',form='formatted')
!         call output(n, alpha, ihyd, time, dt)
         close(11)
      end if
!      if(u(n).gt.2.d9.and.iarrv.eq.0)then
      if(u(n).gt.1.d9.and.iarrv.eq.0)then
         open (11,file='snhydOutput/hyd.d',
     $        access='append',form='formatted')
!         call output(n, alpha, ihyd, time, dt)
         close(11)
         iarrv = 1
         print *,"shock breaks out at t=",time
      end if
      if(iarrv.eq.1)then
         jn=n
         do j = 3, n
            if(u(j).gt.3.d9.or.1.d0/tau(j).lt.1e-20)then
               jn=j-2
               print *,"jn=",jn
               exit
            end if
         enddo
         if(jn.lt.n)then
!               jw = jw+(jn-n)
            n = jn
            print *,"n changes to",n," rho(j)=",1.d0/tau(n+1),
     $              " rho(j-1)=",1.d0/tau(n)
         end if
      endif
      ihyd = ihyd+1
      enddo
c$$$ 10   continue 
      finish = .true.
      call output(n, alpha, ihyd, time, dt)
 99   close(12)
      close(11)
      write(*,*)"at 99"
      open(19,file='snhydOutput/finish.time',status='unknown')
      write(19,*)'just finished'
      close(19)
      stop' normal end.'
      end
