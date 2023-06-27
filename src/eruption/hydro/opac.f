      subroutine opac(n,kap,iphoto,ihyd,useOpacityTable, OpacityTable)
      include 'inclm1.f'
      include 'inclcnst.f'
      real*8 kap(*), muinv, ni, ne, z2
      real*8 k_m, k_h, k_e, k_k, k_grain
      real*8 tauo(mn)
      real*8 logR,logT,logK, dummy
      integer jj, kk, rownum, colnum
      real*8 rowfrac, colfrac
      integer useOpacityTable
      character*8 dummy_char
      character*128 OpacityTable
      common/opdept/tauo
      tauo(n) = 0.d0


      do  j = 3,n
         kap(j) = 0.d0
      end do

      if((useOpacityTable.eq.1) .and. (ihyd.eq.0))then
         ! build opacity array
         write(*,*) 'reading opacity table...'
         open(15, file = OpacityTable, status = 'old')
         read(15,*) dummy_char, (Rarray(jj),jj=1,ncol)
         read(15,*)
         do jj = 1, nrow
            read(15,*) tarray(jj), (kappa(jj,kk),kk=1,ncol)
         end do
         close(15)
      end if

      if(useOpacityTable.eq.1)then
         ! obtain from opacity array
         do j = 3, n
            logR = log10(1.d18/tau(j)/temp(j)**3)
            logT = log10(temp(j))
            ! find which R to interp and how much
            if(logR < Rarray(1)) then
                colnum = 1
                colfrac = 0.0
            else if(logR > Rarray(ncol)) then
                colnum = ncol
                colfrac = 0.0
            else
            do jj = 2, ncol
                if(logR < Rarray(jj)) then
                  colnum = jj-1
                  colfrac=(logR-Rarray(jj-1))/(Rarray(jj)-Rarray(jj-1))
                  exit
                end if
            end do
            end if
            ! find which T to interp and how much
            if(logT < tarray(1)) then
                rownum = 1
                rowfrac = 0.0
            else if(logT > tarray(nrow)) then
                rownum = nrow
                rowfrac = 0.0
            else
            do jj = 2, nrow
                if(logT < tarray(jj)) then
                  rownum = jj-1
                  rowfrac=(logT-tarray(jj-1))/(tarray(jj)-tarray(jj-1))
                  exit
                end if
            end do
            end if
            ! do the interpolation
            if((colnum .eq. ncol) .or. (rownum .eq. nrow)) then
                logK = kappa(rownum, colnum)
            else
                logK = kappa(rownum, colnum)
     $       +rowfrac*(kappa(rownum+1, colnum)-kappa(rownum, colnum))
     $       +colfrac*(kappa(rownum, colnum+1)-kappa(rownum, colnum))
            end if
            kap(j) = 1.d1**logK
         end do
      else
         ! obtain from analytical formula (Kuriyama & Shigeyama 20)
         ! in https://www.astro.princeton.edu/~gk/A403/opac.pdf
         do j = 3, n
            k_m = 1.d-1*(1.d0 - x(j,1) - x(j,2) - x(j,3))
            k_h = 1.1d-25*sqrt(1.d0 - x(j,1) - x(j,2) - x(j,3))*
     $            sqrt(1.d0/tau(j))*
     $              (temp(j)**7.7d0)
            k_e = 2.d-1*(1.d0+x(j,1))/
     $        ((1.d0 + 2.7d11/(tau(j)*temp(j)*temp(j)))*
     $        (1.d0 + (temp(j)/4.5d8)**0.86))
            k_k = 4.d25*(1.d0+x(j,1))*(1.d0-x(j,1)-x(j,2)-x(j,3)+0.001)/
     $          (tau(j)*(temp(j)**3.5))

            kap(j) = k_m + 4.d0/(1.d0/k_h + 1.d0/(k_e + k_k))
         enddo
      end if


      do j = n-1, 3, -1
         tauo(j) = tauo(j+1) + kap(j+1)/tau(j+1)*(rad(j+1)-rad(j))
      enddo

      do j = n-1, 3, -1
         iphoto = j - 1
         if(tauo(j).gt.0.66666667d0)return
      enddo
      return
      end
