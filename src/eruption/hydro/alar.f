      SUBROUTINE ALAR(N)

      include 'inclm1.f'

      REAL*8 PINT(MN), UINT(MN), TAUINT(MN), ARINT(MN),
     &       GINT(MN), G1INT(MN)
      REAL*8 TAUL(MN),PL(MN),UL(MN),ARL(MN),GL(MN),G1L(MN)
      REAL*8 TAUR(MN),PR(MN),UR(MN),ARR(MN),GR(MN),G1R(MN)
      REAL*8 PLD(MN), ULD(MN), GLD(MN), G1LD(MN), TAULD(MN)
      REAL*8 PRD(MN), URD(MN), GRD(MN), G1RD(MN), TAURD(MN)
      REAL*8 FLAT(MN)

      COMMON /INTF/PINT, UINT, TAUINT, ARINT, GINT, G1INT
      COMMON /INTP/ TAUL, PL, UL, ARL, GL, G1L
     &            , TAUR, PR, UR, ARR, GR, G1R
      DO 10 J = 1, N-1
         PL(J+1) = PINT(J)
         PR(J) = PINT(J)
         ARL(J+1) = ARINT(J)
         ARR(J) = ARINT(J)
         TAUL(J+1) = TAUINT(J)
         TAUR(J) = TAUINT(J)
         GL(J+1) = GINT(J)
         GR(J) = GINT(J)
         G1L(J+1) = G1INT(J)
         G1R(J) = G1INT(J)
         UL(J+1) = UINT(J)
         UR(J) = UINT(J)
         FLAT(J) = 0
 10   CONTINUE
      DO 15 J = 3, N-2
         TFLAT = ABS(ABS(P(J+1)-P(J-1))
     $        /ABS(1D-99+(P(J+2)-P(J-2))-1.D0))
c         if(j.le.5)print *,j,tflat
         IF(TFLAT.LT.1D-2)then
            FLAT(J) = 1
         end if
 15   CONTINUE
      TAUL(1) = TAUR(4)
      TAUR(N) = TAUINT(N)
      PL(1) = PR(4)
      PR(N) = PINT(N)
      GL(1) = GR(4)
      GR(N) = GINT(N)
      G1L(1) = G1R(4)
      G1R(N) = G1INT(N)
      UL(1) = -UR(4)
      UR(N) = UINT(N)
      UL(N+1) = UINT(N)
      ARL(1) = -ARR(4)
      ARL(N+1) = ARINT(N)
      ARR(N) = ARINT(N)
      DO 20 J = 1, N
         PLD(J) = PL(J)
         PRD(J) = PR(J)
         TAULD(J) = TAUL(J)
         TAURD(J) = TAUR(J)
         GLD(J) = GL(J)
         GRD(J) = GR(J)
         G1LD(J) = G1L(J)
         G1RD(J) = G1R(J)
         ULD(J) = UL(J)
         URD(J) = UR(J)
 20   CONTINUE
      DO 30 J = 1, N
         PRJ = PRD(J)-P(J)
         PJL = P(J)-PLD(J)
         PRL = PRD(J)-PLD(J)
         PJLR2 = P(J)-0.5*(PLD(J)+PRD(J))
         IF(PRJ*PJL.LT.0.D0)THEN
            PL(J) = P(J)
            PR(J) = P(J)
         ELSE
            IF(PRL*PJLR2.GT.PRL**2/6.D0)PL(J) = 3.D0*P(J)-2.D0*PRD(J)
            IF(PRL*PJLR2+PRL**2/6.D0.LT.0)PR(J) = 3.D0*P(J)-2.D0*PLD(J)
         ENDIF
         TAURJ = TAURD(J)-TAU(J)
         TAUJL = TAU(J)-TAULD(J)
         TAURL = TAURD(J)-TAULD(J)
         TAUJLR  = TAU(J)-0.5*(TAULD(J)+TAURD(J))
         IF(TAURJ*TAUJL.LT.0.D0)THEN
            TAUL(J) = TAU(J)
            TAUR(J) = TAU(J)
         ELSE
            IF(TAURL*TAUJLR .GT.TAURL**2/6.D0)
     &           TAUL(J) = 3.D0*TAU(J)-2.D0*TAURD(J)
            IF(TAURL*TAUJLR +TAURL**2/6.D0.LT.0.D0)
     &           TAUR(J) = 3.D0*TAU(J)-2.D0*TAULD(J)
         ENDIF
         GRJ = GRD(J)-G(J)
         GJL = G(J)-GLD(J)
         GRL = GRD(J)-GLD(J)
         GJLR2 = G(J)-0.5*(GLD(J)+GRD(J))
         IF(GRJ*GJL.LT.0.D0)THEN
            GL(J) = G(J)
            GR(J) = G(J)
         ELSE
            IF(GRL*GJLR2.GT.GRL**2/6.D0)GL(J) = 3.D0*G(J)-2.D0*GRD(J)
            IF(GRL*GJLR2+GRL**2/6.D0.LT.0.D0)GR(J)=3*G(J)-2*GLD(J)
         ENDIF
         G1RJ = G1RD(J)-G1(J)
         G1JL = G1(J)-G1LD(J)
         G1RL = G1RD(J)-G1LD(J)
         G1JLR2 = G1(J)-0.5*(G1LD(J)+G1RD(J))
         IF(G1RJ*G1JL.LT.0.D0)THEN
            G1L(J) = G1(J)
            G1R(J) = G1(J)
         ELSE
            IF(G1RL*G1JLR2.GT.G1RL**2/6.D0)
     &           G1L(J) = 3.D0*G1(J)-2.D0*G1RD(J)
            IF(G1RL*G1JLR2+G1RL**2/6.D0.LT.0.D0)
     &           G1R(J) = 3.D0*G1(J)-2.D0*G1LD(J)
         ENDIF
         URJ = URD(J)-U(J)
         UJL = U(J)-ULD(J)
         URL = URD(J)-ULD(J)
         UJLR2 = U(J)-0.5*(ULD(J)+URD(J))
         IF(URJ*UJL.LT.0.D0)THEN
            UL(J) = U(J)
            UR(J) = U(J)
         ELSE
            IF(URL*UJLR2.GT.URL**2/6.D0)UL(J) = 3.D0*U(J)-2.D0*URD(J)
            IF(URL*UJLR2+URL**2/6.D0.LT.0.D0)UR(J) = 3*U(J)-2*ULD(J)
         ENDIF
 30   CONTINUE
      DO 40 J = 1, N
         TAUL(J) = TAU(J)*FLAT(J)+TAUL(J)*(1-FLAT(J))
         TAUR(J) = TAU(J)*FLAT(J)+TAUR(J)*(1-FLAT(J))
         PL(J) = P(J)*FLAT(J)+PL(J)*(1-FLAT(J))
         PR(J) = P(J)*FLAT(J)+PR(J)*(1-FLAT(J))
         UL(J) = U(J)*FLAT(J)+UL(J)*(1-FLAT(J))
         UR(J) = U(J)*FLAT(J)+UR(J)*(1-FLAT(J))
         GL(J) = G(J)*FLAT(J)+GL(J)*(1-FLAT(J))
         GR(J) = G(J)*FLAT(J)+GR(J)*(1-FLAT(J))
         G1L(J) = G1(J)*FLAT(J)+G1L(J)*(1-FLAT(J))
         G1R(J) = G1(J)*FLAT(J)+G1R(J)*(1-FLAT(J))
 40   CONTINUE
      RETURN
      END
