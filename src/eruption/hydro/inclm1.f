      implicit real*8( a-h, o-z )

      include 'inclmn.f'

      real*8 tau(mn), e(mn), p(mn), u(mn), x(mn,nelem), cs(mn), rad(mn),
     &        g(mn), g1(mn), temp(mn), ye(mn), eu(mn), grv(mn), ar(mn),
     &        kappa(nrow,ncol), tarray(nrow), Rarray(ncol)

      common /physic/ tau, e, p, u, ar, x, cs, temp, g, g1, ye, eu, grv
     $     , rad


