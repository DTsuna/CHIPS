      SUBROUTINE TOTE(N,nadd,E,DMASS,RAD,GRV,U,TENERG,THERME)
      IMPLICIT REAL*8( A-H, O-Z )
      REAL*8 E(*), DMASS(*), RAD(*), GRV(*), U(*), TENERG, THERME
      TENERG  = 0.
      THERME  = 0.
      DO j = 3, N-nadd
         TENERG  = TENERG +E(J)*DMASS(J)
         THERME  = THERME +(E(J)-0.5D0*U(J)*U(J)-RAD(J)*GRV(J))*DMASS(J)
c$$$         write(*,'(i4,1p4e12.4)')j,e(j),u(j),rad(j),grv(j)
c$$$         pause
      enddo

      RETURN
      END



















