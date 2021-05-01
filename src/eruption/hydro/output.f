      subroutine output(n, alpha, ihyd, time, dt)

      include 'inclm1.f'

      real*8 encm(mn), dmass(mn), ps(mn), us(mn), tauo(mn), lum(mn)

      common /mass/encm, dmass
      common /riem/ ps, us
      common /lumi/ lum
      common /opdept/tauo
      write(*,*)'restore at ',ihyd,' time is ',time,' sec'
      write(11,68)n,time,dt,ihyd,(j,rad(j),encm(j),dmass(j),
     $     1./tau(j), u(j), p(j), e(j), temp(j), lum(j)*1d-40, ps(j),
     $     j= 1, n)
c      write(12,69)time,ihyd,(j,(x(j,i),i=1,14),j=1,n)
 68   format(i5,' time', 1pe12.4,' sec',' dt',e12.4,' sec  istep ',i8
     &,/,' no.',5x,'rad',10x,'encm',13x,'dm',12x,'rho',14x,'v'
     &       ,14x,'p',14x,'e',14x,'t',/,(i5,1p10e15.7))
c69   format(1pe12.4,i5,/,(i5,1p13e9.2))
c68    format(i5,' time',1pe10.3,' sec',' dt',e8.2,' sec  istep ',i6
c    &,/,' no.',5x,'rad',6x,'encm', 8x,'dm', 8x,'rho', 7x,'v',
c    &       ,8x,'p',8x,'e',/,(i5,1p7e10.3))
69    format(1pe8.2,i5,/,(i5,1p13e5.3))
      return
      end
