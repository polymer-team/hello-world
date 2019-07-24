        SUBROUTINE RMARIN(IJKL)
C
C       INITIALIZING ROUTINE FOR RANMAR. THE INPUT VALUE SHOULD
C       BE IN THE RANGE:   0 <= IJKL <= 900 000 000
C       TO GET THE STANDARD VALUES IN THE MARSAGLIA - ZAMAN PAPER
C       (I=12, J=34, K=56, L=78) PUT  IJKL = 54217137
C
        COMMON/RASET1/U(97),C,CD,CM,I97,J97,IM_RAN(97)
C
        IJ = IJKL / 30082
        KL = IJKL - IJ * 30082
        I  = MOD(IJ/177,177) + 2
        J  = MOD(IJ,177) + 2
        K  = MOD(KL/169,178) + 1
        L  = MOD(KL,169)
C        WRITE(*,*) 'RANMAR INITIALIZED: ',IJKL,I,J,K,L
        DO   2    II=1,97
         S = 0.
         T = 0.5
         DO   3    JJ=1,24
          M = MOD(MOD(I*J,179)*K,179)
          I = J
          J = K
          K = M
          L = MOD(53*L+1,169)
          IF(MOD(L*M,64).GE.32) S = S + T
3         T = 0.5 * T
2        U(II) = S
        C = 362436. / 16777216.
        CD = 7654321. / 16777216.
        CM = 16777213. / 16777216.
        I97 = 97
        J97 = 33
        DO   4    II=1,97
         IM_RAN(II) = II-1
 4      CONTINUE
        IM_RAN(1) = 97
        RETURN
        END
C
 

