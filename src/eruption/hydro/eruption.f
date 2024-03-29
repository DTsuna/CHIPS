      program eruption
cexpl   one dimensional lagrangian hydrodynamic scheme.
c---  initial data is required.

      include 'inclm1.f'

      integer idev, ni, ipass, jj, kk
      common/nickel/ni
      common / neutrn / tmns
      real*8 msol, kh, eje, AU
      character*128 filename
      character*20 hyd
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
      real*8 relt
      integer output_do
      real*8 when_out(99)
      integer output_init, dummyInt
      logical flag
      integer zero

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

      integer  scaleDeposition, continueTransfer
      integer  useOpacityTable
      integer  discriminant
      character*128 OpacityTable
      logical scaleDepositionFlag
      real*8 scalingRate
      real*8 optical_depth

      real*8 e_charge_tot, injection_time, time_to_cc, dynamicalTime
      integer ejectaCut
      integer fixedCell, innerCell, icell
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
      open(18,file='src/eruption/hydro/para1.d',status='old')
      read(18,*)
      read(18,*)istart,ihydm, jw, iout, idev
      read(18,*)
      read(18,*)cut, dtcfac, eje, nadd
      read(18,*)
      read(18,*)ntp,(tp(k),k=1,ntp)
      read(18,*)
      read(18,*)hyd
      write(*,*)' iout = ',iout
      close(18)


      open(21,file='src/eruption/hydro/eruptPara.d',status='old')
      read(21,*)
      read(21,*)time_to_cc, e_charge_tot, injection_time,
     $     scaleDeposition, scalingRate, continueTransfer,
     $     useOpacityTable, OpacityTable, discriminant
      close(21)
      if(scaleDeposition.eq.0)scaleDepositionFlag = .false.
      if(scaleDeposition.eq.1)scaleDepositionFlag = .true.
      write(*,*)time_to_cc, e_charge_tot, injection_time,
     $          scaleDepositionFlag, scalingRate, continueTransfer,
     $          useOpacityTable, OpacityTable, discriminant

cexpl  construct the initial model
      call init(n, hyd, alpha, cut, istart, time, encmg, eje, nadd,
     $                    dynamicalTime)

      if(discriminant.eq.0)then
        relt = 5.d1
        do jj = 1, 90
          when_out(jj) = (dynamicalTime*2/90.d0)*(jj-1)
        end do
        do jj = 1, 9
          when_out(jj+90) = (dynamicalTime*2/90.d0)*89.d0+
     $    ((time_to_cc -(dynamicalTime*2/90.d0)*89.d0)/9.d0)*jj
        end do
      else
        relt = 5.d0
        do jj = 1, 99
          when_out(jj) = (dynamicalTime*2/99.d0)*(jj-1)
        enddo
      endif
      write(*,*)"********* OUTPUT TIME ***********"
      do jj = 1, 99
        write(*,*),jj,when_out(jj)
      end do

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


      open(97, file='EruptionFiles/parameter.txt'
     $       ,status='unknown',form='formatted')
      write(97,*)"time_to_cc injected_energy inject_duration"
      write(97,*)time_to_cc, e_charge_tot, injection_time
      close(97)


      open(63, file='EruptionFiles/photosphere.txt'
     $       ,status='unknown',form='formatted')
      write(63,*)"time, L_ph, L_edge, r_ph, v_ph, rho_ph, T_ph"

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
          end if
        end if
      end if

      
      do jj=3,n
        call checknan(e(jj),zero,flag)
        if(flag.eqv..true.)then
          print *, 'detected nan. stop.'
          goto 99
          stop
        endif
      enddo


      if(output_do.le.99)then
        if(time.gt.when_out(output_do))then
           write (filename, '("EruptionFiles/result", i2.2, ".txt")')
     $         output_do
           open(91, file=filename,status='unknown',form='formatted')

           do j=3,n
             if(abs(lum(j)).lt.1.d-20)then
               lum(j)=0.d0
             endif
           enddo
           if(ejectaCut.eq.0)then
             write(91,93)n,time,te,ihyd,(j,rad(j),encm(j)+tmns,dmass(j),
     $       1./tau(j), u(j), p(j), e(j), temp(j), lum(j)*1d-40,
     $       kap(j), j= 3, n)
           end if
           if(ejectaCut.eq.1)then
             write(91,93)n,time,te,ihyd,(j,rad(j),encm(j)+tmns,dmass(j),
     $       1./tau(j), u(j), p(j), e(j), temp(j), lum(j)*1d-40,
     $       kap(j), j= fixedCell, n)
           end if
 93     format(i5,' time', 1pe12.4,' sec',' te',e12.4,' erg  istep ',i8
     &  ,/,' no.',5x,'rad',10x,'encm',13x,'dm',12x,'rho',14x,'v'
     &         ,14x,'p',14x,'e',14x,'t',/,(i5,1p10e15.7))
           close(91)
           output_do = output_do + 1
        end if
      else if(discriminant.ne.0)then
        goto 99
        stop
      endif

      innerCell = 3
        if(ejectaCut.eq.1)then
        innerCell = fixedCell
      end if
      call cournt( n, dtcfac, time, dtc ,innerCell)
      dt = min(tp(kp)-time,dtc)
      if(dt.lt.dtc)kp = kp+1
c     stop when timestep becomes too small (currently conservative)
      if(dt.le.1.d-7)then
         write(*,*)dtc,dt,time,kp
         stop' due to too low time step'
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
      call advanc(n,alpha,nadd,dt,dmass,encmg,time,relt,boundr,
     $                e_charge_tot,injection_time,innerCell)

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
        icell = 1
          call opac(n,kap,iphoto,ihyd,useOpacityTable, OpacityTable)
          call radtra(n,ihyd,dt,time,cv,kap,icell)
      end if
      if(ejectaCut.eq.1)then
        icell = innerCell
        if(continueTransfer.eq.1)then
          call opac(n,kap,iphoto,ihyd,useOpacityTable, OpacityTable)
          call radtra(n,ihyd,dt,time,cv,kap,icell)
        end if
      end if

      call tote(n,nadd,e,dmass,rad,grv,u,te,tet)

      if(time-dt.gt.year*86400.d0*365.25d0)then
        write (filename, '("EruptionFiles/intermediate",i2.2,"yr.txt")')
     $         year
        year = year + 1
        output_init = 3
        if(ejectaCut.eq.1)then
          output_init = fixedCell
        end if
        open(98, file=filename,status='unknown',form='formatted')
        write(98,*)"j EnclosedM[g] Rad[cm] Vel[cm/s] Den[g/cc] X_H Y_He
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
99      open(98,file='EruptionFiles/atCCSN.txt',status='unknown'
     $               ,form='formatted')
        write(98,*)"j EnclosedM[g] Rad[cm] Vel[cm/s] Den[g/cc] X_H Y_He
     $ P[erg/cc]"
        if(discriminant.ne.0)then
          output_init = 3
        endif
        do jj = output_init, n
           write(98,'(i0, e18.10, e18.10, e18.10,
     $                   e18.10, e18.10, e18.10, e18.10)'),jj,
     $encm(jj)-encm(output_init-1),rad(jj),u(jj),
     $1.d0/tau(jj),x(jj,1),x(jj,3),p(jj)
        end do
        close(98)
        finish = .true.
        stop
      end if

      if(ihyd.ge.0)then

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
     $        9h    L   ,/,(i5,1p8e9.2))')
     $        (j,encm(j)/msol,rad(j),1.d0/tau(j),ps(j),us(j),e(j)
     $        ,temp(j),lum(j),j=max(1,jw-10),max(10,min(jw+10,n)))
         end if
         if(ejectaCut.eq.1)then
           write(*,'(5h no.  ,8h  mr    ,9hradius   9hdensity  ,
     $        9hpressure ,9hvelocity ,9h     e   ,9h temp    ,
     $        9h    u   ,/,(i5,1p8e9.2))')
     $        (j,encm(j)/msol,rad(j),1.d0/tau(j),ps(j),us(j),e(j)
     $        ,temp(j),u(j),
     $        j=max(fixedCell,jw-10),max(20+fixedCell,min(jw+10,n)))
         end if

      endif

      kp1 = max(kp-1,1)
      if(time.eq.tp(kp1)) then
      end if
      if(u(n).gt.1.d9.and.iarrv.eq.0)then
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
            n = jn
            print *,"n changes to",n," rho(j)=",1.d0/tau(n+1),
     $              " rho(j-1)=",1.d0/tau(n)
         end if
      endif

c     record photosphere information
      optical_depth = 0.0d0
      do j = n, 3, -1
         optical_depth = optical_depth+kap(j)*(rad(j)-rad(j-1))/tau(j)
         if(optical_depth > 2.d0/3.d0) then
           write(63,*)time,lum(j),lum(n),rad(j),u(j),1.d0/tau(j),temp(j)
           exit
         end if
      end do
      ihyd = ihyd+1
      enddo
      finish = .true.
      call output(n, alpha, ihyd, time, dt)
      close(63)
      write(*,*)"at 99"
      stop' normal end.'
      end
