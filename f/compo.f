      SUBROUTINE COMPO(N,ENCM,X,DMASS,ICO,ye)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER JZONE, KATM, MN, ICO
      PARAMETER ( IZONE=10, KATM = 16, MN =1001 )
      REAL*8 ENCM(1),X(MN,14), DMASS(1), ye(1)
      REAL*8 MATM(IZONE,KATM),MASS(IZONE)
      INTEGER NATM(KATM)
      OPEN(22,FILE='inputs/compo.87.d',STATUS='OLD')
c      OPEN(22,FILE='inputs/cmpo1bm3.d',STATUS='OLD')
C      OPEN(22,FILE='inputs/compo.rsg900.d',STATUS='OLD')
      READ(22,*)JZONE
      DO 1 J = 1, JZONE
         READ(22,*)MASS(J),(NATM(L),L=1,6)
         READ(22,*)(MATM(J,L),L=1,6)
         READ(22,*)(NATM(L),L=7,12)
         READ(22,*)(MATM(J,L),L=7,12)
         READ(22,*)(NATM(L),L=13,16)
         READ(22,*)(MATM(J,L),L=13,16)
1        CONTINUE
        CLOSE(22)
      J=1
      DO 3 I = 1, N
        TESM = ENCM(I)/1.989E33-MASS(J)
          IF(TESM.LE.0.0.AND.J.LT.JZONE)THEN
             X(I,1) = MATM(J,1)
             X(I,2) = MATM(J,2)
             X(I,3) = MATM(J,3)
             X(I,4) = MATM(J,5)
             X(I,5) = MATM(J,6)
             X(I,6) = MATM(J,8)
             X(I,7) = MATM(J,10)
             X(I,8) = MATM(J,11)
             X(I,9) = MATM(J,12)
             X(I,10) = MATM(J,13)
             X(I,13) = MATM(J,14)
             X(I,14) = MATM(J,15)
             ENDIF
          IF(TESM.GT.0.0.AND.J.LT.JZONE)THEN
             X(I,1) = MATM(J,1)
             X(I,2) = MATM(J,2)
             X(I,3) = MATM(J,3)
             X(I,4) = MATM(J,5)
             X(I,5) = MATM(J,6)
             X(I,6) = MATM(J,8)
             X(I,7) = MATM(J,10)
             X(I,8) = MATM(J,11)
             X(I,9) = MATM(J,12)
             X(I,10) = MATM(J,13)
             X(I,13) = MATM(J,14)
             X(I,14) = MATM(J,15)
             ENDIF
          IF(J.GE.JZONE)THEN
             X(I,1) = MATM(JZONE,1)
             X(I,2) = MATM(JZONE,2)
             X(I,3) = MATM(JZONE,3)
             X(I,4) = MATM(JZONE,5)
             X(I,5) = MATM(JZONE,6)
             X(I,6) = MATM(JZONE,8)
             X(I,7) = MATM(JZONE,10)
             X(I,8) = MATM(JZONE,11)
             X(I,9) = MATM(JZONE,12)
             X(I,10) = MATM(JZONE,13)
             X(I,13) = MATM(JZONE,14)
             X(I,14) = MATM(JZONE,15)
             ENDIF
          IF(TESM.GT.0.0)J = MIN(JZONE,J+1)
 3     CONTINUE
       DO 15 K = 1, 14
          X(1,K) = X(2,K)
 15    CONTINUE
      CO56M = 0.D0

      DO 10 J = 1, N
         YE(J) = 0.5D0*(1.d0+x(j,1))
         IF(X(J,14).GT.0.D0)ICO = J
         CO56M = CO56M + X(J,14)*DMASS(J)
 10   CONTINUE
      WRITE(*,'(16HTOTAL MASS OF CO,1PE12.4)')CO56M/1.989E33
      RETURN
      END
