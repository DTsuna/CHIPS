      subroutine radtra(nn, istep, dt, time, dedlt, kap, icell)

      include 'inclm1.f'
      include 'inclcnst.f'
      include 'inclold.f'

      integer  nn, n, ni, icell
      common /nickel/ni
c     physical quantities
      real*8 encm(mn), dmass(mn), dedlt(mn), oeu(mn)
      common /mass/encm, dmass
c     interface values
      real*8 us(mn), ps(mn)
      common /riem/ps,us
c     variations
      real*8   dltemp(mn), dll(mn)
c     luminosity
      real*8 difl(mn), lum(mn), srce(mn), srceni(mn)
      common / lumi / lum
      real*8 tauo(mn)
      common/opdept/tauo
c     step
      integer  istep
c     coefficients
      real*8   gam1(mn), gam2(mn), lgam1(mn), lgam2(mn)
      real*8   a11, a12, a21, a22, b1, b2
      real*8   alpha, alpha1, beta, beta1, c1, c2, detm, eps, test
      real*8   r, df(mn), ddf, dr, dr1, mni
c    mass mesh at grid point
      real*8   mgr
c     constant
      real*8 as, kh
c     loop variables
      integer  i, j, it , iter
c     time step and time
      real*8 dt, time
c     source function
      real*8  st
c     opacity for gamma ray
      real*8  kapg, extau, taug(mn)
c     life time of ni and co
      real*8 tni, tco
c     opacities
      real*8  kap(*), kapj
c     normalization factor
      real*8 radn, rnm(mn)

      integer test_n, tesl_n, which, erorr_n
      real*8 test_temp, tesl_temp, which_temp

      parameter ( tni = 7.6032e5, tco = 9.61632e6,
     &   sigma = 5.67032e-5, as = 7.56464e-15, kh = 8.3147515d7 )
      data irad / 1 /
      data kapg / 0.03d0 /


      n = nn
      eps = 1e-6
      epsl = 3e-5
      if(irad.eq.1)then
         write(6,*)' radtra starts now '
         mni = 0.d0
         do i = icell+2, n
            mni = mni + dmass(i)*x(i,3)
            if(x(i,3).eq.0.d0)then
               ni = i
               go to 93
            end if
         enddo
 93      irad=0
         print *,' ni ',ni, mni
      endif
      radn = rad(n)
      radn2 = radn*radn
      taug(n) = 0.0
      st = 3.9e+10*exp(-time/tni)
     &     +7.03e+09*(exp(-time/tco)-exp(-time/tni))

      Do j = icell, n
         rnm(j) = rad(j)/radn
         oeu(j) = eu(j)
         gam1(j) = 0
         gam2(j) = 0
         lgam1(j) = 0
         lgam2(j) = 0
      enddo
      do i = icell, n-2
         rn2 = 0.25*(rnm(n-i+1)+rnm(n-i))**2*radn2
         taug(n-i) = taug(n-i+1) + kapg/(rn2*pi)*dmass(n-i+1)*0.25
      enddo
      taug(icell) = taug(icell+2)
      if(difl(n).eq.0.0)then
         do i = icell+2, n
            rad2 = rnm(i)*rnm(i)
            t2 = temp(i)*temp(i)
            t2m = temp(i-1)*temp(i-1)
            t4 = t2*t2
            t4m = t2m*t2m
            mgr = 0.5*(dmass(i)+dmass(i-1))
            kapj = 0.5*(kap(i)+kap(i-1))
            alpha = 64.0*pi*pi*sigma*rad2*rad2/(mgr*kapj)
            beta = 4.0*pi*rad2*radn2/(kapj*mgr)
            r = min(abs(beta*(t4-t4m)
     $           /(temp(i)+temp(i-1))**4),1d18)
            df(i) = (2.0+r)/(6.0+3.0*r+r*r)
            if(tauo(i).gt.1.d0)then
               lum(i) = alpha*df(i)*(t4m-t4)*radn2*radn2
            else
               lum(i) = 16.d0*rad(i)**2*pi*sigma*temp(i)**4/(tauo(i)
     $              +2.d0/3.d0)
            end if
c           lum(i) = st*dmass(i)*1d-2
         enddo
      else
         do i = icell, n
            lum(i) = difl(i)
         enddo
      endif
      lum(n) = radn2*8.d0*pi*temp(n)**4*sigma
      if(lum(n-1).eq.0.0)lum(n-1)=lum(n)

      do j = icell, n
         srce(j) = 0.d0
         srceni(j) = 0.d0
      enddo

      do j = icell+2, n
         j1 = j - 1
         r2 = rad(j)*rad(j)
         or2 = orad(j)*orad(j)
         r21 = rad(j1)*rad(j1)
         or21 = orad(j1)*orad(j1)
         ab = (r2+or2+rad(j)*orad(j))/3.d0
         ab1 = (r21+or21+rad(j1)*orad(j1))/3.d0
         srce(j) = srceni(j)
     $           -4.d0*pi*(ab*ps(j)*us(j)-ab1*ps(j1)*us(j1))
     $           /dmass(j)-(0.5d0*(u(j)*u(j)-ou(j)*ou(j))
     $           +(rad(j)*grv(j)- orad(j)*ogrv(j)))/dt
      enddo
      srce(icell+2) = srce(icell+2)
      srce(icell)=srce(icell+3)
      srce(icell+1)=srce(icell+2)
      lum(icell+1) = 0.0
c     iteration starts here
      iter = 0


      do it = 1, 100
         rad32 = rnm(icell+2)*rnm(icell+2)
         rad34 = rad32*rad32
         btl = temp(icell+1)
         btl2 = btl*btl
         btl4 = btl2*btl2
         t2 = temp(icell+2)*temp(icell+2)
         t4 = t2*t2
         tl3= t2*temp(icell+2)
         tb2 = t2+btl2+2*temp(icell+2)*btl
         tb4 = 0.0625*tb2*tb2
         bluml = lum(icell+1)
         alpha1 = 64.0*pi*pi*sigma/kap(icell+2)*rad34/dmass(icell+2)
         beta1 = 4.0*pi*rad32*radn2/kap(icell+2)/dmass(icell+2)
         r = abs(beta1*(t4-btl4)/tb4)
         if(r.lt.1e+18)then
            df(icell+2) = (3.0+r)/(6.0+3.0*r+r*r)
            ddf = -r*(4.0+r)/(6.0+3.0*r+r*r)**2
         else
            df(icell+2) = 1.0/r
            ddf = -1.0/(r*r)
         endif
         dr1 = -4.0*beta1*temp(icell+2)*btl*
     $          (t2-temp(icell+2)*btl+btl2)/tb4
         a11 = dt/dmass(icell+2)
         a12 = dedlt(icell+2)
         a21 = 1.0/radn2/radn2
         a22 = -alpha1*(4.0*df(icell+2)*btl4+dr1*ddf*(t4-btl4))
         b1 = -dt/dmass(icell+2)
         b2 = alpha1*(4.0*t4*df(icell+2)+(t4-btl4)*ddf*dr1)
         c1 = (lum(icell+2)-bluml)/dmass(icell+2)*dt +
     $      (eu(icell+2)-oeu(icell+2))
     $        -srce(icell+2)*dt
         c1 = (lum(icell+2)-bluml)/dmass(icell+2)*dt +
     $           (eu(icell+2)-oeu(icell+2))
         c2 = bluml/radn2/radn2
     &        + alpha1*df(icell+2)*(t4-btl4)
         detmi = 1.d0/(b1*b2-a12*a21)
         gam1(icell+2) = b2*a11*detmi
         gam2(icell+2) = -a21*a11*detmi
         lgam1(icell+2) = (b2*c1-a12*c2)*detmi
         lgam2(icell+2) = (b1*c2-a21*c1)*detmi

         gam1(icell+2) = 0.d0
         lgam1(icell+2) = 0.d0
         gam2(icell+2) = a11/a12
         lgam2(icell+2) = c1/a12

         do j = icell+3, n
            radj2 = rnm(j-1)*rnm(j-1)
            radj4 = radj2*radj2
            t2 = temp(j)*temp(j)
            t2m = temp(j-1)*temp(j-1)
            t4 = t2*t2
            t4m = t2m*t2m
            tt2 = t2+t2m+2*temp(j)*temp(j-1)
            tt4 = 0.0625*tt2*tt2
            mgr = 0.5*(dmass(j)+dmass(j-1))
            kapj = 0.5*(kap(j)+kap(j-1))
            alpha = 64.d0*pi*pi*sigma*radj4/(mgr*kapj)
            beta = 4.d0*pi*radj2*radn2/(mgr*kapj)
            r = abs(beta*(t4-t4m)/tt4)
            df(j) = (2.0+r)/(6.0+3.0*r+r*r)
            ddf = -r*(4.0+r)/(6.0+3.0*r+r*r)**2
            dr = 4.0*beta*temp(j)*temp(j-1)*(t2-temp(j)
     $           *temp(j-1)+t2m)/tt4
            if(t4m.le.t4)then
              dr = -4.0*beta*temp(j)*temp(j-1)*(t2-temp(j)
     $             *temp(j-1)+t2m)/tt4
            end if

            a11 = 1.0/dmass(j)*dt
            deu=eu(j)-oeu(j)
            a12 = dedlt(j)
            a21 = 1.0/radn2/radn2
            a22 = -alpha*(4.0*df(j)*t4m+ddf*dr*(t4-t4m))
            b1 = -1.0/dmass(j)*dt
            b2 = alpha*(4.0*t4*df(j)+(t4-t4m)*ddf*dr)
            c1 = (lum(j) - lum(j-1))/dmass(j)*dt + deu
     $           -srce(j)*dt
            c1 = (lum(j) - lum(j-1))/dmass(j)*dt + deu
            c2 = lum(j-1)/radn2/radn2
     &          + alpha*df(j)*(t4-t4m)
            detm = b1*b2-a12*(a21-gam2(j-1)*a22)
            detmi = 1.d0/detm
            gam1(j) = a11*b2*detmi
            gam2(j) = a11*(a22*gam2(j-1)-a21)*detmi
            lgam1(j) = (b2*c1+a12*(a22*lgam2(j-1)-c2))*detmi
            lgam2(j) = ((a22*gam2(j-1)-a21)*c1+b1*(c2-a22*lgam2(j-1)))
     $           *detmi
         enddo
         t2n = temp(n)*temp(n)
         t4n = t2n*t2n

         dll(n) = (radn2*16.d0*pi*sigma*temp(n)**4)
     $    /3.d0/((kap(n)*dmass(n)/4.d0/pi/radn2)+2.d0/3.d0) - lum(n)
         dltemp(n) = -gam2(n)*dll(n)-lgam2(n)
         do j = n-1, icell+2, -1
            j1 = j + 1
            dll(j) = -gam1(j1)*dll(j1) - lgam1(j1)
            dltemp(j) = -gam2(j)*dll(j) - lgam2(j)
         enddo
         test = 0.d0
         tesl = 0.d0
         do j = icell+2, n
            temp(j) = max(t_min,temp(j)*(1.0+max(-0.1d0,dltemp(j))))

            ts = sqrt(temp(j))
            t2 = temp(j)*temp(j)
            lum(j) = lum(j)+dll(j)
            difl(j)=lum(j)
            if(temp(j).eq.t_min)dltemp(j) = 0.d0

            if(abs(dltemp(j)).gt.test)then
              test_n = j
              test_temp = temp(j)
            end if
            if(abs(dll(j)/(lum(j)+1.d0)).gt.tesl)then
              tesl_n = j
              tesl_temp = temp(j)
            end if

            test = max(test,abs(dltemp(j)))
            tesl = max(tesl,abs(dll(j)/(lum(j)+1.d0)))
         enddo
         temp(icell) = temp(icell+3)
         temp(icell+1) = temp(icell+2)
         call eoshelm_e(n,dedlt,temp,e,tau,p,x,grv,rad,eu,g,g1,cs,u,mn,
     $           nelem,time)
         iter = iter + 1
         tesm = max(test,tesl)
         if(test.gt.tesl)then
           which = 1
           which_temp = test_temp
           error_n = test_n
         end if
         if(tesl.ge.test)then
           which = 2
           which_temp = tesl_temp
           error_n = tesl_n
         end if


         if(tesm.lt.eps)goto 400
      enddo
 400  continue
      if(istep.eq.0)then
         open(84,file='EruptionFiles/rtReport.d', status='unknown')
         write(84,*)"*****rtReport*****"
         close(84)
      end if

      if(istep.eq.0)then
         open(82,file='EruptionFiles/rtDebuger.d', status='unknown')
         write(82,*)"*****rtDebuger*****"
         close(82)
      end if

      if(test.gt.epsl)then
         open(84,file='EruptionFiles/rtReport.d', access='append')
         write(84,*),test,istep,time,which,error_n,which_temp
         close(84)
      end if
      if(test.gt.10000.d0)then!epsl to 10000
         write(6,'(a,1p3e12.4,i6)')'test = ',test, tesl, time, istep
         do j = icell+2, n
            if(abs(dltemp(j)).gt.epsl.and.tauo(j).gt.1d-2)then
               write(*,'(i4,1p6e12.4)')j,temp(j),lum(j),dltemp(j),dll(j)
     $              ,eu(j),oeu(j)
            end if
         enddo
         stop
      endif
498   format(' no',7x,'rho',4x,'dltemp',9x,' t',7x,' tau',5x,'lum'
     &       ,7x,'  dll',/,(i5,1p6e11.3))

      return
      end
