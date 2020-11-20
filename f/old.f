      subroutine old(n, nu)

      include 'inclm1.f'
      include 'inclold.f'

      real*8 nu(*)

      do 10 j = 1, n
         oe(j) = e(j)
         op(j) = p(j)
         ou(j) = u(j)
         otau(j) = tau(j)
         ogrv(j) = grv(j)
         orad(j) = rad(j)
         onu(j) = nu(j)
 10      continue

      return
      end
