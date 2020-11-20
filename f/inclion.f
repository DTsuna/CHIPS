      real*8 ne(mn), xion(mn,3), R, Na, ih, ihe1, ihe2, ci
     $     , oye(mn), arad, rho(mn)
      parameter ( R = 1.380622d8/1.660531, Na = 6.022169d23
     $     , ev = 1.60217733d-12, arad = 7.5646d-15
     $     , ih = 13.6*ev, ihe1 = 24.6*ev, ihe2 = 54.4*ev
     $     , ci = 2.07d-16 )

      common/ion/ne,xion,oye,rho

      integer n, j
