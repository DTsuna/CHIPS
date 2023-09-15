      subroutine checknan(a,zero,flag)
        implicit none
        real*8 a,zero
        logical flag
        flag = .false.
        if(a*zero/=zero)then
          flag = .true.
        endif
      end subroutine

