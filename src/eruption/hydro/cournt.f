      subroutine cournt( n, dtcfac, time, dtc, innerCell )

      include 'inclm1.f'
      integer meshc, meshr, innercell
      real*8 dtcc_old, dtcr_old

      real*8  ps(mn), us(mn)
      real*8  dtcc, dtcr, dtccr, dtcfac, dtc, time
      common/riem/ps, us
      dtcc = 0
      dtcr = 1d-15
      do i = innerCell+1, n				!nを3から4に変更
         dtcc_old = dtcc
         dtcr_old = dtcr
         dradi = 1.d0/(rad(i)-rad(i-1))
         dtcc = max(dtcc,cs(i)*abs(dradi)
     $        +max(0.d0,(us(i-1)-us(i))*dradi))
         dtcr = max(dtcr,100.d0*abs(us(i)/rad(i)))
         if(dtcc.ne.dtcc_old)meshc = i
         if(dtcr.ne.dtcr_old)meshr = i

      enddo
!      write(*,*)'Time step is decided by mesh',mesh

      dtc = min(dtcfac/dtcc,1/dtcr)
      if(dtcfac/dtcc.lt.1/dtcr)then
        write(*,*)'Time step is decided by mesh',meshc
      end if
      if(dtcfac/dtcc.gt.1/dtcr)then
        write(*,*)'Time step is decided by mesh',meshr
      end if
      return
      end
