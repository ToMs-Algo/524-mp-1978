C--**--CH3541--524--C:ID--24:5:2000
C--**--CH179--524--Fix--26:5:1999
C--**--CH174--524--P:CAP--25:5:1999
C--**--CH173--524--A:1--25:5:1999
C--**--CH168--524--U:D--20:5:1999
C--**--CH167--524--A:H--20:5:1999
C--**--CH166--524--C:D--20:5:1999
      SUBROUTINE MPABS(X,Y)
C SETS Y = ABS(X) FOR MP NUMBERS X AND Y
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. External Subroutines ..
      EXTERNAL MPSTR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IABS
C     ..
      CALL MPSTR(X,Y)
      Y(1) = IABS(Y(1))
      RETURN

      END
      SUBROUTINE MPADD(X,Y,Z)
C ADDS X AND Y, FORMING RESULT IN Z, WHERE X, Y AND Z ARE MP
C NUMBERS.  FOUR GUARD DIGITS ARE USED, AND THEN R*-ROUNDING.
C     .. Array Arguments ..
      INTEGER X(*),Y(*),Z(*)
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADD2
C     ..
      CALL MPADD2(X,Y,Z,Y,0)
      RETURN

      END
      SUBROUTINE MPADDI(X,IY,Z)
C ADDS MULTIPLE-PRECISION X TO INTEGER IY
C GIVING MULTIPLE-PRECISION Z.
C DIMENSION OF R IN CALLING PROGRAM MUST BE
C AT LEAST 2T+6 (BUT Z(1) MAY BE R(T+5)).
C DIMENSION R(6) BECAUSE RALPH COMPILER ON UNIVAC 1100 COMPUTERS
C OBJECTS TO DIMENSION R(1).
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER IY
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Z(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADD,MPCHK,MPCIM
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(2,6)
      CALL MPCIM(IY,R(T+5))
      CALL MPADD(X,R(T+5),Z)
      RETURN

      END
      SUBROUTINE MPADDQ(X,I,J,Y)
C ADDS THE RATIONAL NUMBER I/J TO MP NUMBER X, MP RESULT IN Y
C DIMENSION OF R MUST BE AT LEAST 2T+6
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER I,J
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADD,MPCHK,MPCQM
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(2,6)
      CALL MPCQM(I,J,R(T+5))
      CALL MPADD(X,R(T+5),Y)
      RETURN

      END
      SUBROUTINE MPADD2(X,Y,Z,Y1,TRUNC)
C CALLED BY MPADD, MPSUB ETC.
C X, Y AND Z ARE MP NUMBERS, Y1 AND TRUNC ARE INTEGERS.
C TO FORCE CALL BY REFERENCE RATHER THAN VALUE/RESULT, Y1 IS
C DECLARED AS AN ARRAY, BUT ONLY Y1(1) IS EVER USED.
C SETS Z = X + Y1(1)*ABS(Y), WHERE Y1(1) = +- Y(1).
C IF TRUNC.EQ.0 R*-ROUNDING IS USED, OTHERWISE TRUNCATION.
C R*-ROUNDING IS DEFINED IN KUKI AND CODI, COMM. ACM
C 16(1973), 223.  (SEE ALSO BRENT, IEEE TC-22(1973), 601.)
C CHECK FOR X OR Y ZERO
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER TRUNC
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*),Y1(*),Z(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER ED,J,MED,RE,RS,S
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADD3,MPCHK,MPERR,MPNZR,MPSTR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IABS
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      IF (X(1).NE.0) GO TO 20
C X = 0 OR NEGLIGIBLE, SO RESULT = +-Y
   10 CALL MPSTR(Y,Z)
      Z(1) = Y1(1)
      RETURN

   20 IF (Y1(1).NE.0) GO TO 40
C Y = 0 OR NEGLIGIBLE, SO RESULT = X
   30 CALL MPSTR(X,Z)
      RETURN
C COMPARE SIGNS
   40 S = X(1)*Y1(1)
      IF (IABS(S).LE.1) GO TO 50
      CALL MPCHK(1,4)
      WRITE (LUN,FMT=9000)
      CALL MPERR
      Z(1) = 0
      RETURN
C COMPARE EXPONENTS
   50 ED = X(2) - Y(2)
      MED = IABS(ED)
      IF (ED) 80,60,110
C EXPONENTS EQUAL SO COMPARE SIGNS, THEN FRACTIONS IF NEC.
   60 IF (S.GT.0) GO TO 90
      DO 70 J = 1,T
          IF (X(J+2)-Y(J+2)) 90,70,120
   70 CONTINUE
C RESULT IS ZERO
      Z(1) = 0
      RETURN
C HERE EXPONENT(Y) .GE. EXPONENT(X)
   80 IF (MED.GT.T) GO TO 10
   90 RS = Y1(1)
      RE = Y(2)
      CALL MPADD3(X,Y,S,MED,RE)
C NORMALIZE, ROUND OR TRUNCATE, AND RETURN
  100 CALL MPNZR(RS,RE,Z,TRUNC)
      RETURN
C ABS(X) .GT. ABS(Y)
  110 IF (MED.GT.T) GO TO 30
  120 RS = X(1)
      RE = X(2)
      CALL MPADD3(Y,X,S,MED,RE)
      GO TO 100

 9000 FORMAT (' *** SIGN NOT 0, +1 OR -1 IN CALL TO MPADD2,',' POSSIBL',
     +       'E OVERWRITING PROBLEM ***')
      END
      SUBROUTINE MPADD3(X,Y,S,MED,RE)
C CALLED BY MPADD2, DOES INNER LOOPS OF ADDITION
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER MED,RE,S
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER C,I,I2,I2P,J,TED
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      TED = T + MED
      I2 = T + 4
      I = I2
      C = 0
C CLEAR GUARD DIGITS TO RIGHT OF X DIGITS
   10 IF (I.LE.TED) GO TO 20
      R(I) = 0
      I = I - 1
      GO TO 10

   20 IF (S.LT.0) GO TO 130
C HERE DO ADDITION, EXPONENT(Y) .GE. EXPONENT(X)
      IF (I.LE.T) GO TO 40
   30 J = I - MED
      R(I) = X(J+2)
      I = I - 1
      IF (I.GT.T) GO TO 30
   40 IF (I.LE.MED) GO TO 60
      J = I - MED
      C = Y(I+2) + X(J+2) + C
      IF (C.LT.B) GO TO 50
C CARRY GENERATED HERE
      R(I) = C - B
      C = 1
      I = I - 1
      GO TO 40
C NO CARRY GENERATED HERE
   50 R(I) = C
      C = 0
      I = I - 1
      GO TO 40

   60 IF (I.LE.0) GO TO 90
      C = Y(I+2) + C
      IF (C.LT.B) GO TO 70
      R(I) = 0
      C = 1
      I = I - 1
      GO TO 60

   70 R(I) = C
      I = I - 1
C NO CARRY POSSIBLE HERE
   80 IF (I.LE.0) RETURN
      R(I) = Y(I+2)
      I = I - 1
      GO TO 80

   90 IF (C.EQ.0) RETURN
C MUST SHIFT RIGHT HERE AS CARRY OFF END
      I2P = I2 + 1
      DO 100 J = 2,I2
          I = I2P - J
          R(I+1) = R(I)
  100 CONTINUE
      R(1) = 1
      RE = RE + 1
      RETURN
C HERE DO SUBTRACTION, ABS(Y) .GT. ABS(X)
  110 J = I - MED
      R(I) = C - X(J+2)
      C = 0
      IF (R(I).GE.0) GO TO 120
C BORROW GENERATED HERE
      C = -1
      R(I) = R(I) + B
  120 I = I - 1
  130 IF (I.GT.T) GO TO 110
  140 IF (I.LE.MED) GO TO 160
      J = I - MED
      C = Y(I+2) + C - X(J+2)
      IF (C.GE.0) GO TO 150
C BORROW GENERATED HERE
      R(I) = C + B
      C = -1
      I = I - 1
      GO TO 140
C NO BORROW GENERATED HERE
  150 R(I) = C
      C = 0
      I = I - 1
      GO TO 140

  160 IF (I.LE.0) RETURN
      C = Y(I+2) + C
      IF (C.GE.0) GO TO 70
      R(I) = C + B
      C = -1
      I = I - 1
      GO TO 160

      END
      SUBROUTINE MPART1(N,Y)
C COMPUTES MP Y = ARCTAN(1/N), ASSUMING INTEGER N .GT. 1.
C USES SERIES ARCTAN(X) = X - X**3/3 + X**5/5 - ...
C DIMENSION OF R IN CALLING PROGRAM MUST BE
C AT LEAST 2T+6
C CHECK LEGALITY OF B, T, M, MXR AND LUN
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      INTEGER Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER B2,I,I2,ID,TS
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADD2,MPCHK,MPCQM,MPDIVI,MPERR,MPMULQ,MPSTR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX0,MIN0
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(2,6)
      IF (N.GT.1) GO TO 10
      WRITE (LUN,FMT=9000)
      CALL MPERR
      Y(1) = 0
      RETURN

   10 I2 = T + 5
      TS = T
C SET SUM TO X = 1/N
      CALL MPCQM(1,N,Y)
C SET ADDITIVE TERM TO X
      CALL MPSTR(Y,R(I2))
      I = 1
      ID = 0
C ASSUME AT LEAST 16-BIT WORD.
      B2 = MAX0(B,64)
      IF (N.LT.B2) ID = (7*B2*B2)/ (N*N)
C MAIN LOOP.  FIRST REDUCE T IF POSSIBLE
   20 T = TS + 2 + R(I2+1) - Y(2)
      IF (T.LT.2) GO TO 50
      T = MIN0(T,TS)
C IF (I+2)*N**2 IS NOT REPRESENTABLE AS AN INTEGER THE DIVISION
C FOLLOWING HAS TO BE PERFORMED IN SEVERAL STEPS.
      IF (I.GE.ID) GO TO 30
      CALL MPMULQ(R(I2),-I, (I+2)*N*N,R(I2))
      GO TO 40

   30 CALL MPMULQ(R(I2),-I,I+2,R(I2))
      CALL MPDIVI(R(I2),N,R(I2))
      CALL MPDIVI(R(I2),N,R(I2))
   40 I = I + 2
C RESTORE T
      T = TS
C ADD TO SUM, USING MPADD2 (FASTER THAN MPADD)
      CALL MPADD2(R(I2),Y,Y,Y,0)
      IF (R(I2).NE.0) GO TO 20
   50 T = TS
      RETURN

 9000 FORMAT (' *** N .LE. 1 IN CALL TO MPART1 ***')
      END
      SUBROUTINE MPASIN(X,Y)
C RETURNS Y = ARCSIN(X), ASSUMING ABS(X) .LE. 1,
C FOR MP NUMBERS X AND Y.
C Y IS IN THE RANGE -PI/2 TO +PI/2.
C METHOD IS TO USE MPATAN, SO TIME IS O(M(T)T).
C DIMENSION OF R MUST BE AT LEAST 5T+12
C CHECK LEGALITY OF B, T, M, MXR AND LUN
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I2,I3
C     ..
C     .. External Functions ..
      INTEGER MPCOMP
      EXTERNAL MPCOMP
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADD,MPATAN,MPCHK,MPCIM,MPDIVI,MPERR,MPMUL,MPPI,MPROOT,
     +         MPSTR,MPSUB
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(5,12)
      I3 = 4*T + 11
      IF (X(1).EQ.0) GO TO 20
      IF (X(2).LE.0) GO TO 30
C HERE ABS(X) .GE. 1.  SEE IF X = +-1
      CALL MPCIM(X(1),R(I3))
      IF (MPCOMP(X,R(I3)).NE.0) GO TO 10
C X = +-1 SO RETURN +-PI/2
      CALL MPPI(Y)
      CALL MPDIVI(Y,2*R(I3),Y)
      RETURN

   10 WRITE (LUN,FMT=9000)
      CALL MPERR
   20 Y(1) = 0
      RETURN
C HERE ABS(X) .LT. 1 SO USE ARCTAN(X/SQRT(1 - X**2))
   30 I2 = I3 - (T+2)
      CALL MPCIM(1,R(I2))
      CALL MPSTR(R(I2),R(I3))
      CALL MPSUB(R(I2),X,R(I2))
      CALL MPADD(R(I3),X,R(I3))
      CALL MPMUL(R(I2),R(I3),R(I3))
      CALL MPROOT(R(I3),-2,R(I3))
      CALL MPMUL(X,R(I3),Y)
      CALL MPATAN(Y,Y)
      RETURN

 9000 FORMAT (' *** ABS(X) .GT. 1 IN CALL TO MPASIN ***')
      END
      SUBROUTINE MPATAN(X,Y)
C RETURNS Y = ARCTAN(X) FOR MP X AND Y, USING AN O(T.M(T)) METHOD
C WHICH COULD EASILY BE MODIFIED TO AN O(SQRT(T)M(T))
C METHOD (AS IN MPEXP1). Y IS IN THE RANGE -PI/2 TO +PI/2.
C FOR AN ASYMPTOTICALLY FASTER METHOD, SEE - FAST MULTIPLE-
C PRECISION EVALUATION OF ELEMENTARY FUNCTIONS
C (BY R. P. BRENT), J. ACM 23 (1976), 242-251,
C AND THE COMMENTS IN MPPIGL.
C DIMENSION OF R IN CALLING PROGRAM MUST BE AT LEAST 5T+12
C CHECK LEGALITY OF B, T, M, MXR AND LUN
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      REAL RX,RY
      INTEGER I,I2,I3,IE,Q,TS
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADD,MPADDI,MPCHK,MPCMR,MPDIV,MPERR,MPMUL,MPMULI,MPMULQ,
     +         MPSQRT,MPSTR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,IABS,MIN0
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(5,12)
      I2 = 3*T + 9
      I3 = I2 + T + 2
      IF (X(1).NE.0) GO TO 10
      Y(1) = 0
      RETURN

   10 CALL MPSTR(X,R(I3))
      IE = IABS(X(2))
      IF (IE.LE.2) CALL MPCMR(X,RX)
      Q = 1
C REDUCE ARGUMENT IF NECESSARY BEFORE USING SERIES
   20 IF (R(I3+1).LT.0) GO TO 30
      IF ((R(I3+1).EQ.0) .AND. ((2* (R(I3+2)+1)).LE.B)) GO TO 30
      Q = 2*Q
      CALL MPMUL(R(I3),R(I3),Y)
      CALL MPADDI(Y,1,Y)
      CALL MPSQRT(Y,Y)
      CALL MPADDI(Y,1,Y)
      CALL MPDIV(R(I3),Y,R(I3))
      GO TO 20
C USE POWER SERIES NOW ARGUMENT IN (-0.5, 0.5)
   30 CALL MPSTR(R(I3),Y)
      CALL MPMUL(R(I3),R(I3),R(I2))
      I = 1
      TS = T
C SERIES LOOP.  REDUCE T IF POSSIBLE.
   40 T = TS + 2 + R(I3+1)
      IF (T.LE.2) GO TO 50
      T = MIN0(T,TS)
      CALL MPMUL(R(I3),R(I2),R(I3))
      CALL MPMULQ(R(I3),-I,I+2,R(I3))
      I = I + 2
      T = TS
      CALL MPADD(Y,R(I3),Y)
      IF (R(I3).NE.0) GO TO 40
C RESTORE T, CORRECT FOR ARGUMENT REDUCTION, AND EXIT
   50 T = TS
      CALL MPMULI(Y,Q,Y)
C CHECK THAT RELATIVE ERROR LESS THAN 0.01 UNLESS EXPONENT
C OF X IS LARGE (WHEN ATAN MIGHT NOT WORK)
      IF (IE.GT.2) RETURN
      CALL MPCMR(Y,RY)
      IF (ABS(RY-ATAN(RX)).LT. (0.01*ABS(RY))) RETURN
      WRITE (LUN,FMT=9000)
C THE FOLLOWING MESSAGE MAY INDICATE THAT B**(T-1) IS TOO SMALL.
      CALL MPERR
      RETURN

 9000 FORMAT (' *** ERROR OCCURRED IN MPATAN, RESULT INCORRECT ***')
      END
      INTEGER FUNCTION MPBASA(X)
C RETURNS THE MP BASE (FIRST WORD IN COMMON).
C X IS A DUMMY MP ARGUMENT.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      MPBASA = B
      RETURN

      END
      SUBROUTINE MPBASB(I,X)
C SETS THE MP BASE (FIRST WORD OF COMMON) TO I.
C I SHOULD BE AN INTEGER SUCH THAT I .GE. 2
C AND (8*I*I-1) IS REPRESENTABLE AS A SINGLE-PRECISION INTEGER.
C X IS A DUMMY MP ARGUMENT (AUGMENT EXPECTS ONE).
C SET BASE TO I, THEN CHECK VALIDITY
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER I
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      B = I
      CALL MPCHK(1,4)
      RETURN

      END
      SUBROUTINE MPBERN(N,P,X)
C COMPUTES THE BERNOULLI NUMBERS B2 = 1/6, B4 = -1/30,
C B6 = 1/42, B8 = -1/30, B10 = 5/66, B12 = -691/2730, ETC.,
C DEFINED BY THE GENERATING FUNCTION Y/(EXP(Y)-1).
C N AND P ARE SINGLE-PRECISION INTEGERS, WITH 2*P .GE. T+2.
C X SHOULD BE A ONE-DIMENSIONAL INTEGER ARRAY OF DIMENSION AT
C LEAST P*N.  THE BERNOULLI NUMBERS B2, B4, ... , B(2N) ARE
C RETURNED IN PACKED FORMAT IN X, WITH B(2J) IN LOCATIONS
C X((J-1)*P+1), ... , X(P*J).  THUS, TO GET B(2J) IN USUAL
C MP FORMAT IN Y, ONE SHOULD CALL MPUNPK (X(IX), Y) AFTER
C CALLING MPBERN, WHERE IX = (J-1)*P+1.
C
C ALTERNATIVELY (SIMPLER BUT NONSTANDARD) -
C X MAY BE A TWO-DIMENSIONAL INTEGER ARRAY DECLARED WITH
C DIMENSION (P, N1), WHERE N1 .GE. N AND 2*P .GE. T+2.
C THEN B2, B4, ... , B(2N) ARE RETURNED IN PACKED FORMAT IN
C X, WITH B(2J) IN X(1,J), ... , X(P,J).  THUS, TO GET
C B(2J) IN USUAL MP FORMAT IN Y ONE SHOULD
C CALL MPUNPK (X(1, J), Y) AFTER CALLING MPBERN.
C
C THE WELL-KNOWN RECURRENCE IS UNSTABLE (LOSING ABOUT 2J BITS
C OF RELATIVE ACCURACY IN THE COMPUTED B(2J)), SO WE USE A
C DIFFERENT RECURRENCE DERIVED BY EQUATING COEFFICIENTS IN
C (EXP(Y)+1)*(2Y/(EXP(2Y)-1)) = 2*(Y/(EXP(Y)-1)).  THE RELATION
C B(2J) = -2*((-1)**J)*FACTORIAL(2J)*ZETA(2J)/((2PI)**(2J))
C IS USED IF ZETA(2J) IS EQUAL TO 1 TO WORKING ACCURACY.
C A DIFFERENT METHOD IS GIVEN BY KNUTH AND BUCKHOLTZ IN
C MATH. COMP. 21 (1967), 663-688.
C THE RELATIVE ERROR IN B(2J) IS O((J**2)*(B**(1-T))).
C TIME IS O(T*(MIN(N, T)**2) + N*M(T)), SPACE = 8T+18.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER N,P
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER B2,I,I2,I3,I4,I5,IX,J,N2
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADDI,MPCHK,MPCIM,MPCQM,MPDIV,MPDIVI,MPERR,MPMUL,MPMULI,
     +         MPPACK,MPPI,MPPWR,MPSTR,MPSUB,MPUNPK
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ALOG,FLOAT,INT,MAX0,MIN0
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      IF (N.LE.0) RETURN
C CHECK LEGALITY OF B, T, M, LUN AND MXR
      CALL MPCHK(8,18)
      IF ((2*P).GE. (T+2)) GO TO 10
      WRITE (LUN,FMT=9000)
      CALL MPERR
      RETURN

   10 I2 = 4*T + 11
      I3 = I2 + T + 2
      I4 = I3 + T + 2
      I5 = I4 + T + 2
      B2 = MAX0(B/2,32)
C COMPUTE UPPER LIMIT FOR RECURRENCE RELATION METHOD.
      N2 = MIN0(N,INT(0.5E0+ALOG(FLOAT(B))*FLOAT(T)/ALOG(4E0)))
C SET ALL RESULTS TO ZERO
      DO 20 I = 1,N2
          IX = (I-1)*P + 2
          X(IX) = 0
   20 CONTINUE
      CALL MPCQM(1,8,R(I2))
      CALL MPSTR(R(I2),R(I3))
      CALL MPCIM(-1,R(I5))
C MAIN LOOP TO GENERATE SCALED BERNOULLI NUMBERS
      DO 60 J = 1,N2
          CALL MPDIVI(R(I3),2,R(I3))
          CALL MPDIVI(R(I5),4,R(I5))
          CALL MPADDI(R(I5),1,R(I4))
          CALL MPDIV(R(I3),R(I4),R(I3))
          IX = (J-1)*P + 1
          CALL MPPACK(R(I3),X(IX))
          IF (J.GE.N2) GO TO 70
          CALL MPDIVI(R(I2),4*J-2,R(I2))
          CALL MPDIVI(R(I2),4*J+4,R(I2))
          CALL MPSTR(R(I2),R(I3))
          DO 50 I = 1,J
              IX = (I-1)*P + 1
              CALL MPUNPK(X(IX),R(I4))
              IF ((J-I).GE.B2) GO TO 30
              CALL MPDIVI(R(I4),8* (2* (J-I)+1)* (J+1-I),R(I4))
              GO TO 40
C HERE SPLIT UP IN CASE WOULD GET OVERFLOW IN ONE CALL TO MPDIVI
   30         CALL MPDIVI(R(I4),4* (J+1-I),R(I4))
              CALL MPDIVI(R(I4),4* (J-I)+2,R(I4))
   40         CALL MPPACK(R(I4),X(IX))
              CALL MPSUB(R(I3),R(I4),R(I3))
   50     CONTINUE
   60 CONTINUE
C NOW UNSCALE RESULTS
   70 CALL MPCIM(1,R(I2))
      IF (N2.LE.1) GO TO 90
      I = N2
   80 CALL MPMULI(R(I2), (4* (N2-I)+4),R(I2))
      CALL MPMULI(R(I2), (4* (N2-I)+2),R(I2))
      I = I - 1
      IX = (I-1)*P + 1
      CALL MPUNPK(X(IX),R(I4))
      CALL MPMUL(R(I2),R(I4),R(I4))
      CALL MPPACK(R(I4),X(IX))
      IF (I.GT.1) GO TO 80
C NOW HAVE B(2J)/FACTORIAL(2J) IN X
      CALL MPCIM(1,R(I2))
   90 DO 100 I = 1,N2
          CALL MPMULI(R(I2),2*I-1,R(I2))
          CALL MPMULI(R(I2),2*I,R(I2))
          IX = (I-1)*P + 1
          CALL MPUNPK(X(IX),R(I4))
          CALL MPMUL(R(I2),R(I4),R(I4))
          CALL MPPACK(R(I4),X(IX))
  100 CONTINUE
C RETURN IF FINISHED
      IF (N.LE.N2) RETURN
C ELSE COMPUTE REMAINING NUMBERS
      CALL MPPI(R(I3))
      CALL MPPWR(R(I3),-2,R(I3))
      CALL MPDIVI(R(I3),-4,R(I3))
      N2 = N2 + 1
      DO 110 I = N2,N
          CALL MPMUL(R(I4),R(I3),R(I4))
          CALL MPMULI(R(I4),2*I-1,R(I4))
          CALL MPMULI(R(I4),2*I,R(I4))
          IX = (I-1)*P + 1
          CALL MPPACK(R(I4),X(IX))
  110 CONTINUE
      RETURN

 9000 FORMAT (' *** P TOO SMALL IN CALL TO MPBERN ***')
      END
      SUBROUTINE MPBESJ(X,NU,Y)
C RETURNS Y = J(NU,X), THE FIRST-KIND BESSEL FUNCTION OF ORDER NU,
C FOR SMALL INTEGER NU, MP X AND Y.  ABS(NU) MUST BE
C .LE. MAX(B, 64).  METHOD IS HANKELS ASYMPTOTIC EXPANSION IF
C ABS(X) LARGE, THE POWER SERIES IF ABS(X) SMALL, AND THE
C BACKWARD RECURRENCE METHOD OTHERWISE.
C RESULTS FOR NEGATIVE ARGUMENTS ARE DEFINED BY
C J(-NU,X) = J(NU,-X) = ((-1)**NU)*J(NU,X).
C ERROR COULD BE INDUCED BY O(B**(1-T)) PERTURBATIONS
C IN X AND Y.  TIME IS O(T.M(T)) FOR FIXED X AND NU, INCREASES
C AS X AND NU INCREASE, UNLESS X LARGE ENOUGH FOR ASYMPTOTIC
C SERIES TO BE USED.   SPACE = 14T+156
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER NU
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      REAL RX
      INTEGER B2,ERROR,I2,I3,I4,IE,K,NUA,TM,TS,TS2
C     ..
C     .. External Subroutines ..
      EXTERNAL MPABS,MPADD,MPBES2,MPCHK,MPCIM,MPCLR,MPCMR,MPDIV,MPDIVI,
     +         MPERR,MPGAMQ,MPHANK,MPMUL,MPPWR,MPSTR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ALOG,FLOAT,IABS,INT,MAX0,MIN0,MOD
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(14,156)
      TS = T
      B2 = MAX0(B,64)
      NUA = IABS(NU)
C CHECK THAT ABS(NU) IS .LE. MAX(B, 64).  THIS RESTRICTION
C ENSURES THAT 4*(NU**2) IS REPRESENTABLE AS AN INTEGER.
      IF (NUA.LE.B2) GO TO 10
      WRITE (LUN,FMT=9000)
      GO TO 100
C CHECK FOR X ZERO
   10 IF (X(1).NE.0) GO TO 20
C J(NU,0) = 0 IF NU .EQ. 0, 1 IF NU .NE. 0
      Y(1) = 0
      IF (NU.EQ.0) CALL MPCIM(1,Y)
      RETURN
C SEE IF ABS(X) SO LARGE THAT NO ACCURACY POSSIBLE
   20 IF (X(2).GE.T) GO TO 90
C X NONZERO SO TRY HANKEL ASYMPTOTIC SERIES WITH ONE GUARD DIGIT
      I2 = 11*T + 36
      CALL MPCLR(R(I2),T+1)
      CALL MPSTR(X,R(I2))
      T = T + 1
      CALL MPHANK(R(I2),NUA,R(I2),ERROR)
      T = TS
      CALL MPSTR(R(I2),Y)
C RETURN IF ASYMPTOTIC SERIES WAS ACCURATE ENOUGH
      IF (ERROR.EQ.0) GO TO 80
C ASYMPTOTIC SERIES INADEQUATE HERE, SO USE POWER SERIES
C MAY NEED TO INCREASE T LATER SO PREPARE FOR THIS
C MAX ALLOWABLE T IS APPROXIMATELY DOUBLE
      TM = 2*T + 20
      I2 = 4*TM + 11
      I3 = I2 + TM + 2
      I4 = I3 + TM + 2
C ZERO TRAILING DIGITS OF R(I2) AND R(I4)
      CALL MPCLR(R(I2),TM)
      CALL MPCLR(R(I4),TM)
      TS2 = T
C NO APPRECIABLE CANCELLATION IN POWER SERIES IF ABS(X) .LT. 1
      IF (X(2).LE.0) GO TO 30
C SHOULD BE OK TO CONVERT TO REAL HERE AS X NOT TOO LARGE OR SMALL.
      CALL MPCMR(X,RX)
C ESTIMATE NUMBER OF DIGITS REQUIRED TO COMPENSATE FOR CANCELLATION
      TS2 = MAX0(TS,T+1+INT((ABS(RX)+ (FLOAT(NUA)+
     +      0.5E0)*ALOG(0.5E0*ABS(RX)))/ALOG(FLOAT(B))))
C IF NEED MORE DIGITS THAN SPACE ALLOWS FOR POWER SERIES THEN
C USE RECURRENCE METHOD INSTEAD
      IF (TS2.GT.TM) GO TO 110
C PREPARE FOR POWER SERIES LOOP
   30 CALL MPDIVI(X,2,R(I4))
      CALL MPPWR(R(I4),NUA,R(I4))
      CALL MPGAMQ(NUA+1,1,R(I3))
      CALL MPDIV(R(I4),R(I3),R(I4))
      CALL MPMUL(X,X,R(I2))
      CALL MPDIVI(R(I2),-4,R(I2))
      T = TS2
      CALL MPSTR(R(I4),R(I3))
      IE = R(I3+1)
      K = 0
C POWER SERIES LOOP, REDUCE T IF POSSIBLE
   40 T = MIN0(TS2,TS2+2+R(I4+1)-IE)
      IF (T.LT.2) GO TO 70
      CALL MPMUL(R(I2),R(I4),R(I4))
      K = K + 1
C MAY NEED TO SPLIT UP CALL TO MPDIVI
      IF (K.GT.B2) GO TO 50
      CALL MPDIVI(R(I4),K* (K+NUA),R(I4))
      GO TO 60
C HERE IT IS SPLIT UP TO AVOID OVERFLOW
   50 CALL MPDIVI(R(I4),K,R(I4))
      CALL MPDIVI(R(I4),K+NUA,R(I4))
C RESTORE T FOR ADDITION
   60 T = TS2
      CALL MPADD(R(I3),R(I4),R(I3))
      IF ((R(I4).NE.0) .AND. (R(I4+1).GE. (R(I3+1)-TS))) GO TO 40
C RESTORE T AND MOVE FINAL RESULT
   70 T = TS
      CALL MPSTR(R(I3),Y)
C CORRECT SIGN IF NU ODD AND NEGATIVE
   80 IF ((NU.LT.0) .AND. (MOD(NUA,2).NE.0)) Y(1) = -Y(1)
      RETURN
C HERE ABS(X) SO LARGE THAT NO SIGNIFICANT DIGITS COULD BE
C GUARANTEED
   90 WRITE (LUN,FMT=9010)
  100 CALL MPERR
      T = TS
      Y(1) = 0
      RETURN
C HERE USE BACKWARD RECURRENCE METHOD WITH TWO GUARD DIGITS
  110 CALL MPABS(X,R(I4))
      T = T + 2
      CALL MPBES2(R(I4),NUA,R(I3))
C CORRECT SIGN IF NUA ODD
      IF (MOD(NUA,2).NE.0) R(I3) = X(1)*R(I3)
      GO TO 70

 9000 FORMAT (' *** ABS(NU) TOO LARGE IN CALL TO MPBESJ ***')
 9010 FORMAT (' *** ABS(X) TOO LARGE IN CALL TO MPBESJ ***')
      END
      SUBROUTINE MPBES2(X,NU,Y)
C USES THE BACKWARD RECURRENCE METHOD TO EVALUATE Y = J(NU,X),
C WHERE X AND Y ARE MP NUMBERS, NU (THE INDEX) IS AN INTEGER,
C AND J IS THE BESSEL FUNCTION OF THE FIRST KIND.  ASSUMES THAT
C 0 .LE. NU .LE. MAX(B,64) AND X .GT. 0. ALSO ASSUMED THAT X
C CAN BE CONVERTED TO REAL WITHOUT FLOATING-POINT OVERFLOW OR
C UNDERFLOW.  FOR NORMALIZATION THE IDENTITY
C J(0,X) + 2*J(2,X) + 2*J(4,X) + ... = 1  IS USED.
C CALLED BY MPBESJ AND NOT RECOMMENDED FOR INDEPENDENT USE.
C SPACE = 8T+18
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER NU
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      REAL FLNU,RT,RX,RY
      INTEGER I,I2,I3,I3S,I4,I5,I6,NU1
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADD,MPCHK,MPCIM,MPCMR,MPDIV,MPERR,MPMUL,MPMULI,MPREC,
     +         MPSTR,MPSUB
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ALOG,AMAX1,FLOAT,INT,MAX0,MOD
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(8,18)
C CHECK LEGALITY OF NU AND X
      IF ((NU.GE.0) .AND. (NU.LE.MAX0(B,64)) .AND.
     +    (X(1).EQ.1)) GO TO 10
      WRITE (LUN,FMT=9000)
      CALL MPERR
      Y(1) = 0
      RETURN
C ASSUME CONVERSION TO REAL IS POSSIBLE WITHOUT OVERFLOW (TRUE
C WHEN CALLED BY MPBESJ ELSE MPHANK OR POWER SERIES WOULD BE USED).
   10 CALL MPCMR(X,RX)
C COMPUTE STARTING POINT NU1 FOR BACKWARD RECURRENCE
      FLNU = FLOAT(MAX0(1,NU))
      RY = AMAX1(1E0,ALOG(2E0*FLNU/RX)-1E0)
C 1.35914 IS E/2 ROUNDED DOWN, 1.35915 IS E/2 ROUNDED UP
      RY = (FLNU*RY+0.5E0*FLOAT(T)*ALOG(FLOAT(B)))/ (1.35914E0*RX)
      RY = AMAX1(2E0,RY)
      RT = RY
C ITERATE AN EVEN NUMBER OF TIMES TO OVERESTIMATE NU1
      DO 20 I = 1,4
          RT = AMAX1(2E0,RY/ALOG(RT))
   20 CONTINUE
      NU1 = 2 + INT(1.35915E0*RX*RT)
      I2 = 3*T + 9
      I3 = I2 + T + 2
      I4 = I3 + T + 2
      I5 = I4 + T + 2
      I6 = I5 + T + 2
      CALL MPCIM(MOD(NU1+1,2),R(I6))
      CALL MPREC(X,R(I2))
      CALL MPMULI(R(I2),2,R(I2))
      R(I3) = 0
      CALL MPCIM(1,R(I4))
C BACKWARD RECURRENCE LOOP
   30 CALL MPMUL(R(I4),R(I2),R(I5))
      CALL MPMULI(R(I5),NU1,R(I5))
      CALL MPSUB(R(I5),R(I3),R(I5))
      NU1 = NU1 - 1
C FASTER TO INTERCHANGE POINTERS THAN MP NUMBERS
      I3S = I3
      I3 = I4
      I4 = I5
      I5 = I3S
      IF (MOD(NU1,2).NE.0) GO TO 40
C NU1 EVEN SO UPDATE NORMALIZING SUM
      IF (NU1.EQ.0) CALL MPMULI(R(I6),2,R(I6))
      CALL MPADD(R(I6),R(I4),R(I6))
C SAVE UNNORMALIZED RESULT IF NU1 .EQ. NU
   40 IF (NU1.EQ.NU) CALL MPSTR(R(I4),Y)
      IF (NU1.GT.0) GO TO 30
C NORMALIZE RESULT AND RETURN
      CALL MPDIV(Y,R(I6),Y)
      RETURN

 9000 FORMAT (' *** NU .LT. 0 OR NU TOO LARGE OR X .LE. 0 IN CALL',' T',
     +       'O MPBES2 ***')
      END
      SUBROUTINE MPCAM(A,X)
C CONVERTS THE HOLLERITH STRING A TO AN MP NUMBER X.
C A CAN BE A STRING OF DIGITS ACCEPTABLE TO ROUTINE MPIN
C AND TERMINATED BY A DOLLAR ($), E.G. 7H-5.367$,
C OR ONE OF THE FOLLOWING SPECIAL STRINGS -
C            EPS  (MP MACHINE-PRECISION, SEE MPEPS),
C            EUL  (EULERS CONSTANT 0.5772..., SEE MPEUL),
C            MAXR (LARGEST VALID MP NUMBER, SEE MPMAXR),
C            MINR (SMALLEST POSTIVE MP NUMBER, SEE MPMINR),
C            PI   (PI = 3.14..., SEE MPPI).
C ONLY THE FIRST TWO CHARACTERS OF THESE STRINGS ARE CHECKED.
C SPACE REQUIRED IS NO MORE THAN 5*T+L+14, WHERE L IS THE
C NUMBER OF CHARACTERS IN THE STRING A (EXCLUDING $).
C IF SPACE IS LESS 3*T+L+11 THE STRING A WILL EFFECTIVELY BE TRUNCATED
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER A(*),X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER ERROR,I,I2,N
C     ..
C     .. Local Arrays ..
      INTEGER C(6),D(2)
C     ..
C     .. External Subroutines ..
      EXTERNAL MPEPS,MPERR,MPEUL,MPIN,MPMAXR,MPMINR,MPPI,MPUPK
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
C     .. Data statements ..
      DATA C(1)/1HA/,C(2)/1HE/,C(3)/1HI/
      DATA C(4)/1HM/,C(5)/1HP/,C(6)/1HU/
C     ..
C UNPACK FIRST 2 CHARACTERS OF A
      CALL MPUPK(A,D,2,N)
      IF (N.NE.2) GO TO 10
C SET X TO ZERO AFTER SAVING A(1) IN CASE A AND X COINCIDE
      I = A(1)
      X(1) = 0
C CHECK FOR SPECIAL STRINGS
      IF ((D(1).EQ.C(2)) .AND. (D(2).EQ.C(5))) CALL MPEPS(X)
      IF ((D(1).EQ.C(2)) .AND. (D(2).EQ.C(6))) CALL MPEUL(X)
      IF ((D(1).EQ.C(4)) .AND. (D(2).EQ.C(1))) CALL MPMAXR(X)
      IF ((D(1).EQ.C(4)) .AND. (D(2).EQ.C(3))) CALL MPMINR(X)
      IF ((D(1).EQ.C(5)) .AND. (D(2).EQ.C(3))) CALL MPPI(X)
C RETURN IF X NONZERO (SO ONE OF ABOVE TESTS SUCCEEDED)
      IF (X(1).NE.0) RETURN
C RESTORE A(1) AND UNPACK, THEN CALL MPIN TO DECODE.
      A(1) = I
   10 I2 = 3*T + 12
      CALL MPUPK(A,R(I2),MXR+1-I2,N)
      CALL MPIN(R(I2),X,N,ERROR)
      IF (ERROR.EQ.0) RETURN
      WRITE (LUN,FMT=9000)
      CALL MPERR
      RETURN

 9000 FORMAT (' *** ERROR IN HOLLERITH CONSTANT IN CALL TO MPCAM ***')
      END
      SUBROUTINE MPCDM(DX,Z)
C CONVERTS DOUBLE-PRECISION NUMBER DX TO MULTIPLE-PRECISION Z.
C SOME NUMBERS WILL NOT CONVERT EXACTLY ON MACHINES
C WITH BASE OTHER THAN TWO, FOUR OR SIXTEEN.
C THIS ROUTINE IS NOT CALLED BY ANY OTHER ROUTINE IN MP,
C SO MAY BE OMITTED IF DOUBLE-PRECISION IS NOT AVAILABLE.
C CHECK LEGALITY OF B, T, M, MXR AND LUN
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION DX
C     ..
C     .. Array Arguments ..
      INTEGER Z(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DB,DJ
      INTEGER I,I2,IB,IE,K,RE,RS,TP
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPDIVI,MPMULI,MPNZR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE,FLOAT,IDINT,MAX0
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(1,4)
      I2 = T + 4
C CHECK SIGN
      IF (DX) 20,10,30
C IF DX = 0D0 RETURN 0
   10 Z(1) = 0
      RETURN
C DX .LT. 0D0
   20 RS = -1
      DJ = -DX
      GO TO 40
C DX .GT. 0D0
   30 RS = 1
      DJ = DX
   40 IE = 0
   50 IF (DJ.LT.1D0) GO TO 60
C INCREASE IE AND DIVIDE DJ BY 16.
      IE = IE + 1
      DJ = 0.0625D0*DJ
      GO TO 50

   60 IF (DJ.GE.0.0625D0) GO TO 70
      IE = IE - 1
      DJ = 16D0*DJ
      GO TO 60
C NOW DJ IS DY DIVIDED BY SUITABLE POWER OF 16
C SET EXPONENT TO 0
   70 RE = 0
C DB = DFLOAT(B) IS NOT ANSI STANDARD SO USE FLOAT AND DBLE
      DB = DBLE(FLOAT(B))
C CONVERSION LOOP (ASSUME DOUBLE-PRECISION OPS. EXACT)
      DO 80 I = 1,I2
          DJ = DB*DJ
          R(I) = IDINT(DJ)
          DJ = DJ - DBLE(FLOAT(R(I)))
   80 CONTINUE
C NORMALIZE RESULT
      CALL MPNZR(RS,RE,Z,0)
      IB = MAX0(7*B*B,32767)/16
      TP = 1
C NOW MULTIPLY BY 16**IE
      IF (IE) 90,130,110
   90 K = -IE
      DO 100 I = 1,K
          TP = 16*TP
          IF ((TP.LE.IB) .AND. (TP.NE.B) .AND. (I.LT.K)) GO TO 100
          CALL MPDIVI(Z,TP,Z)
          TP = 1
  100 CONTINUE
      RETURN

  110 DO 120 I = 1,IE
          TP = 16*TP
          IF ((TP.LE.IB) .AND. (TP.NE.B) .AND. (I.LT.IE)) GO TO 120
          CALL MPMULI(Z,TP,Z)
          TP = 1
  120 CONTINUE
  130 RETURN

      END
      SUBROUTINE MPCHK(I,J)
C CHECKS LEGALITY OF B, T, M, MXR AND LUN WHICH SHOULD BE SET
C IN COMMON.
C THE CONDITION ON MXR (THE DIMENSION OF R IN COMMON) IS THAT
C MXR .GE. (I*T + J)
C FIRST CHECK THAT LUN IN RANGE 1 TO 99, IF NOT PRINT ERROR
C MESSAGE ON LOGICAL UNIT 6.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER I,J
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER IB,MX
C     ..
C     .. External Subroutines ..
      EXTERNAL MPERR
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      IF ((0.LT.LUN) .AND. (LUN.LT.100)) GO TO 10
      WRITE (6,FMT=9000) LUN
      LUN = 6
      CALL MPERR
C NOW CHECK LEGALITY OF B, T AND M
   10 IF (B.GT.1) GO TO 20
      WRITE (LUN,FMT=9010) B
      CALL MPERR
   20 IF (T.GT.1) GO TO 30
      WRITE (LUN,FMT=9020) T
      CALL MPERR
   30 IF (M.GT.T) GO TO 40
      WRITE (LUN,FMT=9030)
      CALL MPERR
C 8*B*B-1 SHOULD BE REPRESENTABLE, IF NOT WILL OVERFLOW
C AND MAY BECOME NEGATIVE, SO CHECK FOR THIS
   40 IB = 4*B*B - 1
      IF ((IB.GT.0) .AND. ((2*IB+1).GT.0)) GO TO 50
      WRITE (LUN,FMT=9040)
      CALL MPERR
C CHECK THAT SPACE IN COMMON IS SUFFICIENT
   50 MX = I*T + J
      IF (MXR.GE.MX) RETURN
C HERE COMMON IS TOO SMALL, SO GIVE ERROR MESSAGE.
      WRITE (LUN,FMT=9050) I,J,MX,MXR,T
      CALL MPERR
      RETURN

 9000 FORMAT (' *** LUN =',I10,' ILLEGAL IN CALL TO MPCHK,',' PERHAPS ',
     +       'NOT SET BEFORE CALL TO AN MP ROUTINE ***')
 9010 FORMAT (' *** B =',I10,' ILLEGAL IN CALL TO MPCHK,',/' PERHAPS N',
     +       'OT SET BEFORE CALL TO AN MP ROUTINE ***')
 9020 FORMAT (' *** T =',I10,' ILLEGAL IN CALL TO MPCHK,',/' PERHAPS N',
     +       'OT SET BEFORE CALL TO AN MP ROUTINE ***')
 9030 FORMAT (' *** M .LE. T IN CALL TO MPCHK,',/' PERHAPS NOT SET BEF',
     +       'ORE CALL TO AN MP ROUTINE ***')
 9040 FORMAT (' *** B TOO LARGE IN CALL TO MPCHK ***')
 9050 FORMAT (' *** MXR TOO SMALL OR NOT SET TO DIM(R) BEFORE CALL',
     +       ' TO AN MP ROUTINE ***',/' *** MXR SHOULD BE AT LEAST',I3,
     +       '*T +',I4,' =',I6,'  ***',/' *** ACTUALLY MXR =',I10,', A',
     +       'ND T =',I10,'  ***')
      END
      SUBROUTINE MPCIM(IX,Z)
C CONVERTS INTEGER IX TO MULTIPLE-PRECISION Z.
C CHECK LEGALITY OF B, T, M, MXR AND LUN
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER IX
C     ..
C     .. Array Arguments ..
      INTEGER Z(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I,N
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPMUL2
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(1,4)
      N = IX
      IF (N) 20,10,30
   10 Z(1) = 0
      RETURN

   20 N = -N
      Z(1) = -1
      GO TO 40

   30 Z(1) = 1
C SET EXPONENT TO T
   40 Z(2) = T
C CLEAR FRACTION
      DO 50 I = 2,T
          Z(I+1) = 0
   50 CONTINUE
C INSERT N
      Z(T+2) = N
C NORMALIZE BY CALLING MPMUL2
      CALL MPMUL2(Z,1,Z,1)
      RETURN

      END
      SUBROUTINE MPCLR(X,N)
C SETS X(T+3), ... , X(N+2) TO ZERO, USEFUL
C IF PRECISION IS GOING TO BE INCREASED.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I,I2,I3
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      IF (N.LE.T) RETURN
      I2 = T + 3
      I3 = N + 2
      DO 10 I = I2,I3
          X(I) = 0
   10 CONTINUE
      RETURN

      END
      SUBROUTINE MPCMD(X,DZ)
C CONVERTS MULTIPLE-PRECISION X TO DOUBLE-PRECISION DZ.
C ASSUMES X IS IN ALLOWABLE RANGE FOR DOUBLE-PRECISION
C NUMBERS.   THERE IS SOME LOSS OF ACCURACY IF THE
C EXPONENT IS LARGE.
C CHECK LEGALITY OF B, T, M, MXR AND LUN
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION DZ
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DB,DZ2
      INTEGER I,TM
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPERR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS,DBLE,DLOG,FLOAT
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(1,4)
      DZ = 0D0
      IF (X(1).EQ.0) RETURN
C DB = DFLOAT(B) IS NOT ANSI STANDARD, SO USE FLOAT AND DBLE
      DB = DBLE(FLOAT(B))
      DO 10 I = 1,T
          DZ = DB*DZ + DBLE(FLOAT(X(I+2)))
          TM = I
C CHECK IF FULL DOUBLE-PRECISION ACCURACY ATTAINED
          DZ2 = DZ + 1D0
C TEST BELOW NOT ALWAYS EQUIVALENT TO - IF (DZ2.LE.DZ) GO TO 20,
C FOR EXAMPLE ON CYBER 76.
          IF ((DZ2-DZ).LE.0D0) GO TO 20
   10 CONTINUE
C NOW ALLOW FOR EXPONENT
   20 DZ = DZ* (DB** (X(2)-TM))
C CHECK REASONABLENESS OF RESULT.
      IF (DZ.LE.0D0) GO TO 30
C LHS SHOULD BE .LE. 0.5 BUT ALLOW FOR SOME ERROR IN DLOG
      IF (DABS(DBLE(FLOAT(X(2)))- (DLOG(DZ)/DLOG(DBLE(FLOAT(B)))+
     +    0.5D0)).GT.0.6D0) GO TO 30
      IF (X(1).LT.0) DZ = -DZ
      RETURN
C FOLLOWING MESSAGE INDICATES THAT X IS TOO LARGE OR SMALL -
C TRY USING MPCMDE INSTEAD.
   30 WRITE (LUN,FMT=9000)
      CALL MPERR
      RETURN

 9000 FORMAT (' *** FLOATING-POINT OVER/UNDER-FLOW IN MPCMD ***')
      END
      SUBROUTINE MPCMDE(X,N,DX)
C RETURNS INTEGER N AND DOUBLE-PRECISION DX SUCH THAT MP
C X = DX*10**N (APPROXIMATELY), WHERE 1 .LE. ABS(DX) .LT. 10
C UNLESS DX = 0.    SPACE = 6T+14
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION DX
      INTEGER N
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I2
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPCMD,MPCMEF
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS,DBLE,FLOAT
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      IF (X(1).NE.0) GO TO 10
      N = 0
      DX = 0D0
      RETURN
C CHECK LEGALITY OF B, T, M, LUN AND MXR
   10 CALL MPCHK(6,14)
      I2 = 5*T + 13
      CALL MPCMEF(X,N,R(I2))
      CALL MPCMD(R(I2),DX)
      IF (DABS(DX).LT.10D0) RETURN
C HERE DX WAS ROUNDED UP TO TEN
      N = N + 1
      DX = DBLE(FLOAT(R(I2)))
      RETURN

      END
      SUBROUTINE MPCMEF(X,N,Y)
C GIVEN MP X, RETURNS INTEGER N AND MP Y SUCH THAT X = (10**N)*Y
C AND 1 .LE. ABS(Y) .LT. 10 (UNLESS X .EQ. 0, WHEN N .EQ. 0 AND
C Y .EQ. 0).
C IT IS ASSUMED THAT X IS NOT SO LARGE OR SMALL THAT N OVERFLOWS.
C SPACE = 5T+12
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      REAL RLY,RY
      INTEGER I2,IY,J,N1,TEN
C     ..
C     .. External Functions ..
      INTEGER MPCMPI
      EXTERNAL MPCMPI
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPCIM,MPCMR,MPDIV,MPDIVI,MPERR,MPMUL,MPMULI,MPPWR,
     +         MPSTR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ALOG,FLOAT,IABS,INT
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
C     .. Data statements ..
C FOR OCTAL OUTPUT CHANGE 10 TO 8 BELOW, ETC.
      DATA TEN/10/
C     ..
C CHECK FOR X ZERO
      IF (X(1).NE.0) GO TO 10
      N = 0
      Y(1) = 0
      RETURN
C X NONZERO, CHECK LEGALITY OF B, T, M, LUN AND MXR
   10 CALL MPCHK(5,12)
      CALL MPSTR(X,Y)
      N = 0
      I2 = 4*T + 11
C LOOP UP TO FOUR TIMES (USUALLY ONE IS SUFFICIENT)
      DO 20 J = 1,4
          IY = Y(2)
          Y(2) = 0
          CALL MPCMR(Y,RY)
          Y(2) = IY
C ESTIMATE LOG10 (ABS(Y))
          RLY = (FLOAT(IY)*ALOG(FLOAT(B))+ALOG(ABS(RY)))/
     +          ALOG(FLOAT(TEN))
          N1 = INT(RLY)
C CHECK IF N1 OBVIOUSLY OVERFLOWED
          IF (ABS(RLY-FLOAT(N1)).GT.16E0) GO TO 30
C FOLLOWING AVOIDS POSSIBILITY OF R(I2) OVERFLOWING BELOW
          IF ((J.EQ.1) .AND. (IABS(N1).GT. (M/4))) N1 = N1/2
C LEAVE J LOOP IF N1 SMALL
          IF (IABS(N1).LE.2) GO TO 40
C DIVIDE BY TEN**N1
          N = N + N1
          CALL MPCIM(TEN,R(I2))
          CALL MPPWR(R(I2),IABS(N1),R(I2))
          IF (R(I2).EQ.0) GO TO 30
          IF (N1.GT.0) CALL MPDIV(Y,R(I2),Y)
          IF (N1.LT.0) CALL MPMUL(R(I2),Y,Y)
   20 CONTINUE
   30 WRITE (LUN,FMT=9000)
      CALL MPERR
      RETURN

   40 IF (Y(1).EQ.0) GO TO 30
C LOOP DIVIDING BY TEN UNTIL ABS(Y) .LT. 1
   50 IF (Y(2).LE.0) GO TO 70
      N = N + 1
      CALL MPDIVI(Y,TEN,Y)
      GO TO 50
C LOOP MULTIPLYING BY TEN UNTIL ABS(Y) .GE. 1
   60 IF (Y(2).GT.0) GO TO 80
   70 N = N - 1
      CALL MPMULI(Y,TEN,Y)
      GO TO 60
C CHECK FOR POSSIBILITY THAT ROUNDING UP WAS TO TEN
   80 IY = Y(1)
      Y(1) = 1
      IF (MPCMPI(Y,TEN).LT.0) GO TO 90
C IT WAS, SO SET Y TO 1 AND ADD ONE TO EXPONENT
      CALL MPCIM(1,Y)
      N = N + 1
C RESTORE SIGN OF Y AND RETURN
   90 Y(1) = IY
      RETURN

 9000 FORMAT (' *** ERROR OCCURRED IN MPCMEF, PROBABLY OVERFLOW',' CAU',
     +       'SED BY LARGE EXPONENT ***')
      END
      SUBROUTINE MPCMF(X,Y)
C FOR MP X AND Y, RETURNS FRACTIONAL PART OF X IN Y,
C I.E., Y = X - INTEGER PART OF X (TRUNCATED TOWARDS 0).
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I,IL,IP,X2,XS
C     ..
C     .. External Subroutines ..
      EXTERNAL MPNZR,MPSTR
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      IF (X(1).NE.0) GO TO 20
C RETURN 0 IF X = 0
   10 Y(1) = 0
      RETURN

   20 X2 = X(2)
C RETURN 0 IF EXPONENT SO LARGE THAT NO FRACTIONAL PART
      IF (X2.GE.T) GO TO 10
C IF EXPONENT NOT POSITIVE CAN RETURN X
      IF (X2.GT.0) GO TO 30
      CALL MPSTR(X,Y)
      RETURN
C CLEAR ACCUMULATOR
   30 DO 40 I = 1,X2
          R(I) = 0
   40 CONTINUE
      IL = X2 + 1
C MOVE FRACTIONAL PART OF X TO ACCUMULATOR
      DO 50 I = IL,T
          R(I) = X(I+2)
   50 CONTINUE
      DO 60 I = 1,4
          IP = I + T
          R(IP) = 0
   60 CONTINUE
      XS = X(1)
C NORMALIZE RESULT AND RETURN
      CALL MPNZR(XS,X2,Y,1)
      RETURN

      END
      SUBROUTINE MPCMI(X,IZ)
C CONVERTS MULTIPLE-PRECISION X TO INTEGER IZ,
C ASSUMING THAT X NOT TOO LARGE (ELSE USE MPCMIM).
C X IS TRUNCATED TOWARDS ZERO.
C IF INT(X)IS TOO LARGE TO BE REPRESENTED AS A SINGLE-
C PRECISION INTEGER, IZ IS RETURNED AS ZERO.  THE USER
C MAY CHECK FOR THIS POSSIBILITY BY TESTING IF
C ((X(1).NE.0).AND.(X(2).GT.0).AND.(IZ.EQ.0)) IS TRUE ON
C RETURN FROM MPCMI.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER IZ
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I,IZS,J,J1,K,KX,X2,XS
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      XS = X(1)
      IZ = 0
      IF (XS.EQ.0) RETURN
      IF (X(2).LE.0) RETURN
      X2 = X(2)
      DO 10 I = 1,X2
          IZS = IZ
          IZ = B*IZ
          IF (I.LE.T) IZ = IZ + X(I+2)
C CHECK FOR SIGNS OF INTEGER OVERFLOW
          IF ((IZ.LE.0) .OR. (IZ.LE.IZS)) GO TO 30
   10 CONTINUE
C CHECK THAT RESULT IS CORRECT (AN UNDETECTED OVERFLOW MAY
C HAVE OCCURRED).
      J = IZ
      DO 20 I = 1,X2
          J1 = J/B
          K = X2 + 1 - I
          KX = 0
          IF (K.LE.T) KX = X(K+2)
          IF (KX.NE. (J-B*J1)) GO TO 30
          J = J1
   20 CONTINUE
      IF (J.NE.0) GO TO 30
C RESULT CORRECT SO RESTORE SIGN AND RETURN
      IZ = XS*IZ
      RETURN
C HERE OVERFLOW OCCURRED (OR X WAS UNNORMALIZED), SO
C RETURN ZERO.
   30 IZ = 0
      RETURN

      END
      SUBROUTINE MPCMIM(X,Y)
C RETURNS Y = INTEGER PART OF X (TRUNCATED TOWARDS 0), FOR MP X AND Y.
C USE IF Y TOO LARGE TO BE REPRESENTABLE AS A SINGLE-PRECISION INTEGER.
C (ELSE COULD USE MPCMI).
C CHECK LEGALITY OF B, T, M, MXR AND LUN
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I,IL
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPSTR
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(1,4)
      CALL MPSTR(X,Y)
      IF (Y(1).EQ.0) RETURN
      IL = Y(2) + 1
C IF EXPONENT LARGE ENOUGH RETURN Y = X
      IF (IL.GT.T) RETURN
C IF EXPONENT SMALL ENOUGH RETURN ZERO
      IF (IL.GT.1) GO TO 10
      Y(1) = 0
      RETURN
C SET FRACTION TO ZERO
   10 DO 20 I = IL,T
          Y(I+2) = 0
   20 CONTINUE
      RETURN

      END
      INTEGER FUNCTION MPCMPA(X,Y)
C COMPARES ABS(X) WITH ABS(Y) FOR MP X AND Y,
C RETURNING +1 IF ABS(X) .GT. ABS(Y),
C           -1 IF ABS(X) .LT. ABS(Y),
C AND        0 IF ABS(X) .EQ. ABS(Y)
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Local Scalars ..
      INTEGER XS,YS
C     ..
C     .. External Functions ..
      INTEGER MPCOMP
      EXTERNAL MPCOMP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IABS
C     ..
      XS = X(1)
      X(1) = IABS(XS)
      YS = Y(1)
      Y(1) = IABS(YS)
      MPCMPA = MPCOMP(X,Y)
      X(1) = XS
      Y(1) = YS
      RETURN

      END
      INTEGER FUNCTION MPCMPI(X,I)
C COMPARES MP NUMBER X WITH INTEGER I, RETURNING
C     +1 IF X .GT. I,
C      0 IF X .EQ. I,
C     -1 IF X .LT. I
C DIMENSION OF R IN COMMON AT LEAST 2T+6
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER I
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. External Functions ..
      INTEGER MPCOMP
      EXTERNAL MPCOMP
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPCIM
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(2,6)
C CONVERT I TO MULTIPLE-PRECISION AND COMPARE
      CALL MPCIM(I,R(T+5))
      MPCMPI = MPCOMP(X,R(T+5))
      RETURN

      END
      INTEGER FUNCTION MPCMPR(X,RI)
C COMPARES MP NUMBER X WITH REAL NUMBER RI, RETURNING
C     +1 IF X .GT. RI,
C      0 IF X .EQ. RI,
C     -1 IF X .LT. RI
C DIMENSION OF R IN COMMON AT LEAST 2T+6
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      REAL RI
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. External Functions ..
      INTEGER MPCOMP
      EXTERNAL MPCOMP
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPCRM
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(2,6)
C CONVERT RI TO MULTIPLE-PRECISION AND COMPARE
      CALL MPCRM(RI,R(T+5))
      MPCMPR = MPCOMP(X,R(T+5))
      RETURN

      END
      SUBROUTINE MPCMR(X,RZ)
C CONVERTS MULTIPLE-PRECISION X TO SINGLE-PRECISION RZ.
C ASSUMES X IN ALLOWABLE RANGE.  THERE IS SOME LOSS OF
C ACCURACY IF THE EXPONENT IS LARGE.
C CHECK LEGALITY OF B, T, M, MXR AND LUN
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      REAL RZ
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      REAL RB,RZ2
      INTEGER I,TM
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPERR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ALOG,FLOAT
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(1,4)
      RZ = 0E0
      IF (X(1).EQ.0) RETURN
      RB = FLOAT(B)
      DO 10 I = 1,T
          RZ = RB*RZ + FLOAT(X(I+2))
          TM = I
C CHECK IF FULL SINGLE-PRECISION ACCURACY ATTAINED
          RZ2 = RZ + 1E0
          IF (RZ2.LE.RZ) GO TO 20
   10 CONTINUE
C NOW ALLOW FOR EXPONENT
   20 RZ = RZ* (RB** (X(2)-TM))
C CHECK REASONABLENESS OF RESULT
      IF (RZ.LE.0E0) GO TO 30
C LHS SHOULD BE .LE. 0.5, BUT ALLOW FOR SOME ERROR IN ALOG
      IF (ABS(FLOAT(X(2))- (ALOG(RZ)/ALOG(FLOAT(B))+0.5E0)).GT.
     +    0.6E0) GO TO 30
      IF (X(1).LT.0) RZ = -RZ
      RETURN
C FOLLOWING MESSAGE INDICATES THAT X IS TOO LARGE OR SMALL -
C TRY USING MPCMRE INSTEAD.
   30 WRITE (LUN,FMT=9000)
      CALL MPERR
      RETURN

 9000 FORMAT (' *** FLOATING-POINT OVER/UNDER-FLOW IN MPCMR ***')
      END
      SUBROUTINE MPCMRE(X,N,RX)
C RETURNS INTEGER N AND SINGLE-PRECISION RX SUCH THAT MP
C X = RX*10**N (APPROXIMATELY), WHERE 1 .LE. ABS(RX) .LT. 10
C UNLESS RX = 0.    SPACE = 6T+14
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      REAL RX
      INTEGER N
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I2
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPCMEF,MPCMR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,FLOAT
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      IF (X(1).NE.0) GO TO 10
      N = 0
      RX = 0E0
      RETURN
C CHECK LEGALITY OF B, T, M, LUN AND MXR
   10 CALL MPCHK(6,14)
      I2 = 5*T + 13
      CALL MPCMEF(X,N,R(I2))
      CALL MPCMR(R(I2),RX)
      IF (ABS(RX).LT.10E0) RETURN
C HERE RX WAS ROUNDED UP TO TEN
      N = N + 1
      RX = FLOAT(R(I2))
      RETURN

      END
      INTEGER FUNCTION MPCOMP(X,Y)
C COMPARES THE MULTIPLE-PRECISION NUMBERS X AND Y,
C RETURNING +1 IF X .GT. Y,
C           -1 IF X .LT. Y,
C AND        0 IF X .EQ. Y.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I,T2
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      IF (X(1)-Y(1)) 10,30,20
C X .LT. Y
   10 MPCOMP = -1
      RETURN
C X .GT. Y
   20 MPCOMP = 1
      RETURN
C SIGN(X) = SIGN(Y), SEE IF ZERO
   30 IF (X(1).NE.0) GO TO 40
C X = Y = 0
      MPCOMP = 0
      RETURN
C HAVE TO COMPARE EXPONENTS AND FRACTIONS
   40 T2 = T + 2
      DO 50 I = 2,T2
          IF (X(I)-Y(I)) 60,50,70
   50 CONTINUE
C NUMBERS EQUAL
      MPCOMP = 0
      RETURN
C ABS(X) .LT. ABS(Y)
   60 MPCOMP = -X(1)
      RETURN
C ABS(X) .GT. ABS(Y)
   70 MPCOMP = X(1)
      RETURN

      END
      SUBROUTINE MPCOS(X,Y)
C RETURNS Y = COS(X) FOR MP X AND Y, USING MPSIN AND MPSIN1.
C DIMENSION OF R IN COMMON AT LEAST 5T+12.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I2
C     ..
C     .. External Functions ..
      INTEGER MPCMPI
      EXTERNAL MPCMPI
C     ..
C     .. External Subroutines ..
      EXTERNAL MPABS,MPCHK,MPCIM,MPDIVI,MPPI,MPSIN,MPSIN1,MPSUB
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      IF (X(1).NE.0) GO TO 10
C COS(0) = 1
      CALL MPCIM(1,Y)
      RETURN
C CHECK LEGALITY OF B, T, M, LUN AND MXR
   10 CALL MPCHK(5,12)
      I2 = 3*T + 12
C SEE IF ABS(X) .LE. 1
      CALL MPABS(X,Y)
      IF (MPCMPI(Y,1).LE.0) GO TO 20
C HERE ABS(X) .GT. 1 SO USE COS(X) = SIN(PI/2 - ABS(X)),
C COMPUTING PI/2 WITH ONE GUARD DIGIT.
      T = T + 1
      CALL MPPI(R(I2))
      CALL MPDIVI(R(I2),2,R(I2))
      T = T - 1
      CALL MPSUB(R(I2),Y,Y)
      CALL MPSIN(Y,Y)
      RETURN
C HERE ABS(X) .LE. 1 SO USE POWER SERIES
   20 CALL MPSIN1(Y,Y,0)
      RETURN

      END
      SUBROUTINE MPCOSH(X,Y)
C RETURNS Y = COSH(X) FOR MP NUMBERS X AND Y, X NOT TOO LARGE.
C USES MPEXP, DIMENSION OF R IN COMMON AT LEAST 5T+12
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I2
C     ..
C     .. External Subroutines ..
      EXTERNAL MPABS,MPADD,MPCHK,MPCIM,MPDIVI,MPEXP,MPREC
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      IF (X(1).NE.0) GO TO 10
C COSH(0) = 1
      CALL MPCIM(1,Y)
      RETURN
C CHECK LEGALITY OF B, T, M, LUN AND MXR
   10 CALL MPCHK(5,12)
      I2 = 4*T + 11
      CALL MPABS(X,R(I2))
C IF ABS(X) TOO LARGE MPEXP WILL PRINT ERROR MESSAGE
C INCREASE M TO AVOID OVERFLOW WHEN COSH(X) REPRESENTABLE
      M = M + 2
      CALL MPEXP(R(I2),R(I2))
      CALL MPREC(R(I2),Y)
      CALL MPADD(R(I2),Y,Y)
C RESTORE M.  IF RESULT OVERFLOWS OR UNDERFLOWS, MPDIVI WILL
C ACT ACCORDINGLY.
      M = M - 2
      CALL MPDIVI(Y,2,Y)
      RETURN

      END
      SUBROUTINE MPCQM(I,J,Q)
C CONVERTS THE RATIONAL NUMBER I/J TO MULTIPLE PRECISION Q.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER I,J
C     ..
C     .. Array Arguments ..
      INTEGER Q(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I1,J1
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCIM,MPDIVI,MPERR,MPGCD
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      I1 = I
      J1 = J
      CALL MPGCD(I1,J1)
      IF (J1) 20,10,30
   10 WRITE (LUN,FMT=9000)
      CALL MPERR
      Q(1) = 0
      RETURN

   20 I1 = -I1
      J1 = -J1
   30 CALL MPCIM(I1,Q)
      IF (J1.NE.1) CALL MPDIVI(Q,J1,Q)
      RETURN

 9000 FORMAT (' *** J = 0 IN CALL TO MPCQM ***')
      END
      SUBROUTINE MPCRM(RX,Z)
C CONVERTS SINGLE-PRECISION NUMBER RX TO MULTIPLE-PRECISION Z.
C SOME NUMBERS WILL NOT CONVERT EXACTLY ON MACHINES
C WITH BASE OTHER THAN TWO, FOUR OR SIXTEEN.
C CHECK LEGALITY OF B, T, M, MXR AND LUN
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      REAL RX
C     ..
C     .. Array Arguments ..
      INTEGER Z(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      REAL RB,RJ
      INTEGER I,I2,IB,IE,K,RE,RS,TP
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPDIVI,MPMULI,MPNZR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC FLOAT,INT,MAX0
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(1,4)
      I2 = T + 4
C CHECK SIGN
      IF (RX) 20,10,30
C IF RX = 0E0 RETURN 0
   10 Z(1) = 0
      RETURN
C RX .LT. 0E0
   20 RS = -1
      RJ = -RX
      GO TO 40
C RX .GT. 0E0
   30 RS = 1
      RJ = RX
   40 IE = 0
   50 IF (RJ.LT.1E0) GO TO 60
C INCREASE IE AND DIVIDE RJ BY 16.
      IE = IE + 1
      RJ = 0.0625E0*RJ
      GO TO 50

   60 IF (RJ.GE.0.0625E0) GO TO 70
      IE = IE - 1
      RJ = 16E0*RJ
      GO TO 60
C NOW RJ IS DY DIVIDED BY SUITABLE POWER OF 16.
C SET EXPONENT TO 0
   70 RE = 0
      RB = FLOAT(B)
C CONVERSION LOOP (ASSUME SINGLE-PRECISION OPS. EXACT)
      DO 80 I = 1,I2
          RJ = RB*RJ
          R(I) = INT(RJ)
          RJ = RJ - FLOAT(R(I))
   80 CONTINUE
C NORMALIZE RESULT
      CALL MPNZR(RS,RE,Z,0)
      IB = MAX0(7*B*B,32767)/16
      TP = 1
C NOW MULTIPLY BY 16**IE
      IF (IE) 90,130,110
   90 K = -IE
      DO 100 I = 1,K
          TP = 16*TP
          IF ((TP.LE.IB) .AND. (TP.NE.B) .AND. (I.LT.K)) GO TO 100
          CALL MPDIVI(Z,TP,Z)
          TP = 1
  100 CONTINUE
      RETURN

  110 DO 120 I = 1,IE
          TP = 16*TP
          IF ((TP.LE.IB) .AND. (TP.NE.B) .AND. (I.LT.IE)) GO TO 120
          CALL MPMULI(Z,TP,Z)
          TP = 1
  120 CONTINUE
  130 RETURN

      END
      SUBROUTINE MPDAW(X,Y)
C RETURNS Y = DAWSONS INTEGRAL (X)
C           = EXP(-X**2)*(INTEGRAL FROM 0 TO X OF EXP(U**2)DU),
C FOR MP X AND Y.    SPACE = 5T+17.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER B2,I,I2,I3,I4,IER,TS,XS
C     ..
C     .. External Subroutines ..
      EXTERNAL MPABS,MPADD,MPCHK,MPCLR,MPDIVI,MPERF3,MPEXP,MPMUL,MPMULQ,
     +         MPNEG,MPSTR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX0,MIN0
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      XS = X(1)
      IF (XS.NE.0) GO TO 10
C DAW(0) = 0
      Y(1) = 0
      RETURN
C CHECK LEGALITY OF B, T, M, LUN AND MXR
   10 CALL MPCHK(5,17)
      I2 = 2*T + 9
      I3 = I2 + T + 3
      I4 = I3 + T + 3
      B2 = 2*MAX0(B,64)
C WORK WITH ABS(X)
      CALL MPABS(X,R(I4))
C TRY ASYMPTOTIC SERIES
      CALL MPERF3(R(I4),R(I4),1,IER)
      IF (IER.NE.0) GO TO 20
      CALL MPSTR(R(I4),Y)
      Y(1) = XS*Y(1)
      RETURN
C ASYMPTOTIC SERIES NOT ACCURATE ENOUGH SO USE POWER SERIES
C WITH ONE GUARD DIGIT.
   20 CALL MPCLR(R(I4),T+1)
      CALL MPSTR(X,R(I4))
      T = T + 1
      CALL MPMUL(R(I4),R(I4),R(I4))
      CALL MPNEG(R(I4),R(I4))
      CALL MPEXP(R(I4),R(I4))
      T = T - 1
      CALL MPCLR(R(I2),T+1)
      CALL MPSTR(X,R(I2))
      T = T + 1
      CALL MPMUL(R(I2),R(I4),R(I3))
      CALL MPMUL(R(I2),R(I2),R(I4))
      CALL MPSTR(R(I3),R(I2))
      I = 0
      TS = T
C POWER SERIES LOOP, REDUCE T IF POSSIBLE
   30 T = TS + 2 + R(I3+1) - R(I2+1)
      IF (T.LE.2) GO TO 60
      T = MIN0(T,TS)
      I = I + 1
      CALL MPMUL(R(I4),R(I3),R(I3))
C SEE IF NEXT CALL TO MPMULQ HAS TO BE SPLIT UP
      IF (I.GE.B2) GO TO 40
      CALL MPMULQ(R(I3),2*I-1,I* (2*I+1),R(I3))
      GO TO 50

   40 CALL MPMULQ(R(I3),2*I-1,2*I+1,R(I3))
      CALL MPDIVI(R(I3),I,R(I3))
C RESTORE T FOR ADDITION
   50 T = TS
      CALL MPADD(R(I2),R(I3),R(I2))
      IF (R(I3).NE.0) GO TO 30
C RESTORE T AND RETURN
   60 T = TS - 1
      CALL MPSTR(R(I2),Y)
      RETURN

      END
      INTEGER FUNCTION MPDGA(X,N)
C RETURNS THE N-TH DIGIT OF THE MP NUMBER X FOR 1 .LE. N .LE. T.
C RETURNS ZERO IF X IS ZERO OR N .LE. 0 OR N .GT. T.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      MPDGA = 0
      IF ((X(1).NE.0) .AND. (N.GT.0) .AND. (N.LE.T)) MPDGA = X(N+2)
      RETURN

      END
      SUBROUTINE MPDGB(I,X,N)
C SETS THE N-TH DIGIT OF THE MP NUMBER X TO I.
C N MUST BE IN THE RANGE 1 .LE. N .LE T,
C I MUST BE IN THE RANGE 0 .LE. I .LT. B
C (AND I .NE. 0 IF N .EQ. 1).
C THE SIGN AND EXPONENT OF X ARE UNCHANGED.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER I,N
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. External Subroutines ..
      EXTERNAL MPERR
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      IF ((N.GT.0) .AND. (N.LE.T)) GO TO 10
      WRITE (LUN,FMT=9000)
      GO TO 20

   10 IF ((I.GE.0) .AND. (I.LT.B) .AND. ((I+N).GT.1)) GO TO 30
      WRITE (LUN,FMT=9010)
   20 CALL MPERR
      RETURN

   30 X(N+2) = I
      RETURN

 9000 FORMAT (' *** DIGIT POSITION ILLEGAL IN CALL TO MPDGB ***')
 9010 FORMAT (' *** DIGIT VALUE ILLEGAL IN CALL TO MPDGB ***')
      END
      INTEGER FUNCTION MPDIGA(X)
C RETURNS THE NUMBER OF MP DIGITS (SECOND WORD IN COMMON).
C X IS A DUMMY MP ARGUMENT.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      MPDIGA = T
      RETURN

      END
      SUBROUTINE MPDIGB(I,X)
C SETS THE NUMBER OF MP DIGITS (SECOND WORD OF COMMON) TO I.
C I SHOULD BE AN INTEGER SUCH THAT I .GE. 2
C X IS A DUMMY MP ARGUMENT (AUGMENT EXPECTS ONE).
C WARNING *** MP NUMBERS MUST BE DECLARED AS INTEGER ARRAYS OF
C         *** DIMENSION AT LEAST I+2. MPDIGB DOES NOT CHECK THIS.
C SET DIGITS TO I, THEN CHECK VALIDITY
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER I
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      T = I
      CALL MPCHK(1,4)
      RETURN

      END
      SUBROUTINE MPDIV(X,Y,Z)
C SETS Z = X/Y, FOR MP X, Y AND Z.
C DIMENSION OF R IN CALLING PROGRAM MUST BE AT LEAST 4T+10
C (BUT Z(1) MAY BE R(3T+9)).
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*),Z(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I,I2,IE,IZ3
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPERR,MPEXT,MPMUL,MPOVFL,MPREC,MPUNFL
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(4,10)
C CHECK FOR DIVISION BY ZERO
      IF (Y(1).NE.0) GO TO 10
      WRITE (LUN,FMT=9000)
      CALL MPERR
      Z(1) = 0
      RETURN
C SPACE USED BY MPREC IS 4T+10 WORDS, BUT CAN OVERLAP SLIGHTLY.
   10 I2 = 3*T + 9
C CHECK FOR X = 0
      IF (X(1).NE.0) GO TO 20
      Z(1) = 0
      RETURN
C INCREASE M TO AVOID OVERFLOW IN MPREC
   20 M = M + 2
C FORM RECIPROCAL OF Y
      CALL MPREC(Y,R(I2))
C SET EXPONENT OF R(I2) TO ZERO TO AVOID OVERFLOW IN MPMUL
      IE = R(I2+1)
      R(I2+1) = 0
      I = R(I2+2)
C MULTIPLY BY X
      CALL MPMUL(X,R(I2),Z)
      IZ3 = Z(3)
      CALL MPEXT(I,IZ3,Z)
C RESTORE M, CORRECT EXPONENT AND RETURN
      M = M - 2
      Z(2) = Z(2) + IE
      IF (Z(2).GE. (-M)) GO TO 30
C UNDERFLOW HERE
      CALL MPUNFL(Z)
      RETURN

   30 IF (Z(2).LE.M) RETURN
C OVERFLOW HERE
      WRITE (LUN,FMT=9010)
      CALL MPOVFL(Z)
      RETURN

 9000 FORMAT (' *** ATTEMPTED DIVISION BY ZERO IN CALL TO MPDIV ***')
 9010 FORMAT (' *** OVERFLOW OCCURRED IN MPDIV ***')
      END
      SUBROUTINE MPDIVI(X,IY,Z)
C DIVIDES MP X BY THE SINGLE-PRECISION INTEGER IY GIVING MP Z.
C THIS IS MUCH FASTER THAN DIVISION BY AN MP NUMBER.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER IY
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Z(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER B2,C,C2,I,I2,IQ,IQJ,IR,J,J1,J11,J2,K,KH,R1,RE,RS
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPERR,MPNZR,MPSTR,MPUNFL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX0
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      RS = X(1)
      J = IY
      IF (J) 20,10,30
   10 WRITE (LUN,FMT=9000)
      GO TO 210

   20 J = -J
      RS = -RS
   30 RE = X(2)
C CHECK FOR ZERO DIVIDEND
      IF (RS.EQ.0) GO TO 110
C CHECK FOR DIVISION BY B
      IF (J.NE.B) GO TO 40
      CALL MPSTR(X,Z)
      IF (RE.LE. (-M)) GO TO 220
      Z(1) = RS
      Z(2) = RE - 1
      RETURN
C CHECK FOR DIVISION BY 1 OR -1
   40 IF (J.NE.1) GO TO 50
      CALL MPSTR(X,Z)
      Z(1) = RS
      RETURN

   50 C = 0
      I2 = T + 4
      I = 0
C IF J*B NOT REPRESENTABLE AS AN INTEGER HAVE TO SIMULATE
C LONG DIVISION.   ASSUME AT LEAST 16-BIT WORD.
      B2 = MAX0(8*B,32767/B)
      IF (J.GE.B2) GO TO 120
C LOOK FOR FIRST NONZERO DIGIT IN QUOTIENT
   60 I = I + 1
      C = B*C
      IF (I.LE.T) C = C + X(I+2)
      R1 = C/J
      IF (R1) 200,60,70
C ADJUST EXPONENT AND GET T+4 DIGITS IN QUOTIENT
   70 RE = RE + 1 - I
      R(1) = R1
      C = B* (C-J*R1)
      KH = 2
      IF (I.GE.T) GO TO 90
      KH = 1 + T - I
      DO 80 K = 2,KH
          I = I + 1
          C = C + X(I+2)
          R(K) = C/J
          C = B* (C-J*R(K))
   80 CONTINUE
      IF (C.LT.0) GO TO 200
      KH = KH + 1
   90 DO 100 K = KH,I2
          R(K) = C/J
          C = B* (C-J*R(K))
  100 CONTINUE
      IF (C.LT.0) GO TO 200
C NORMALIZE AND ROUND RESULT
  110 CALL MPNZR(RS,RE,Z,0)
      RETURN
C HERE NEED SIMULATED DOUBLE-PRECISION DIVISION
  120 C2 = 0
      J1 = J/B
      J2 = J - J1*B
      J11 = J1 + 1
C LOOK FOR FIRST NONZERO DIGIT
  130 I = I + 1
      C = B*C + C2
      C2 = 0
      IF (I.LE.T) C2 = X(I+2)
      IF (C-J1) 130,140,150
  140 IF (C2.LT.J2) GO TO 130
C COMPUTE T+4 QUOTIENT DIGITS
  150 RE = RE + 1 - I
      K = 1
      GO TO 170
C MAIN LOOP FOR LARGE ABS(IY) CASE
  160 K = K + 1
      IF (K.GT.I2) GO TO 110
      I = I + 1
C GET APPROXIMATE QUOTIENT FIRST
  170 IR = C/J11
C NOW REDUCE SO OVERFLOW DOES NOT OCCUR
      IQ = C - IR*J1
      IF (IQ.LT.B2) GO TO 180
C HERE IQ*B WOULD POSSIBLY OVERFLOW SO INCREASE IR
      IR = IR + 1
      IQ = IQ - J1
  180 IQ = IQ*B - IR*J2
      IF (IQ.GE.0) GO TO 190
C HERE IQ NEGATIVE SO IR WAS TOO LARGE
      IR = IR - 1
      IQ = IQ + J
  190 IF (I.LE.T) IQ = IQ + X(I+2)
      IQJ = IQ/J
C R(K) = QUOTIENT, C = REMAINDER
      R(K) = IQJ + IR
      C = IQ - J*IQJ
      IF (C.GE.0) GO TO 160
C CARRY NEGATIVE SO OVERFLOW MUST HAVE OCCURRED
  200 CALL MPCHK(1,4)
      WRITE (LUN,FMT=9010)
  210 CALL MPERR
      Z(1) = 0
      RETURN
C UNDERFLOW HERE
  220 CALL MPUNFL(Z)
      RETURN

 9000 FORMAT (' *** ATTEMPTED DIVISION BY ZERO IN CALL TO MPDIVI ***')
 9010 FORMAT (' *** INTEGER OVERFLOW IN MPDIVI, B TOO LARGE ***')
      END
      SUBROUTINE MPDUMP(X)
C DUMPS OUT THE MP NUMBER X (SIGN, EXPONENT, FRACTION DIGITS),
C USEFUL FOR DEBUGGING PURPOSES.
C EMBEDDED BLANKS SHOULD BE INTERPRETED AS ZEROS. (THEY COULD BE
C AVOIDED BY USING J INSTEAD OF I FORMAT, BUT THIS IS NONSTANDARD.)
C CHECK LEGALITY OF B, T, M, MXR AND LUN
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I,T2
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(1,4)
      IF (X(1).NE.0) GO TO 10
C IF X = 0 JUST WRITE SIGN AS REMAINDER UNDEFINED
      WRITE (LUN,FMT=9000) X(1)
      RETURN

   10 T2 = T + 2
      IF (B.LE.10) WRITE (LUN,FMT=9000) (X(I),I=1,T2)
      IF ((B.GT.10) .AND. (B.LE.100)) WRITE (LUN,FMT=9010) (X(I),I=1,
     +    T2)
      IF ((B.GT.100) .AND. (B.LE.1000)) WRITE (LUN,FMT=9020) (X(I),I=1,
     +    T2)
      IF (B.GT.1000) WRITE (LUN,FMT=9030) (X(I),I=1,T2)
C ASSUME RECORDS OF UP TO 79 CHARACTERS OK ON UNIT LUN
      RETURN

 9000 FORMAT (1X,I2,I12,4X,60I1,/ (19X,60I1))
 9010 FORMAT (1X,I2,I12,4X,30I2,/ (19X,30I2))
 9020 FORMAT (1X,I2,I12,4X,20I3,/ (19X,20I3))
 9030 FORMAT (1X,I2,I16,4X,8I7,/ (23X,8I7))
      END
      SUBROUTINE MPEI(X,Y)
C RETURNS Y = EI(X) = -E1(-X)
C           = (PRINCIPAL VALUE INTEGRAL FROM -INFINITY TO X OF
C              EXP(U)/U DU),
C FOR MP NUMBERS X AND Y,
C USING THE POWER SERIES FOR SMALL ABS(X), THE ASYMPTOTIC SERIES FOR
C LARGE ABS(X), AND THE CONTINUED FRACTION FOR INTERMEDIATE NEGATIVE
C X.  RELATIVE ERROR IN Y IS SMALL EXCEPT IF X IS VERY CLOSE TO THE
C ZERO  0.37250741078136663446... OF EI(X),
C AND THEN THE ABSOLUTE ERROR IN Y IS O(B**(1-T)).
C IN ANY CASE THE ERROR IN Y COULD BE INDUCED BY AN
C O(B**(1-T)) RELATIVE PERTURBATION IN X.
C TIME IS O(T.M(T)), SPACE = 10T+38.
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      REAL C,CP,RAX,RB,RT,RX
      INTEGER B2,I,I2,I3,I4,J,K,TD,TM,TM2,TS,TS1,TS2,XS
C     ..
C     .. External Functions ..
      INTEGER MPCMPA,MPCMPR
      EXTERNAL MPCMPA,MPCMPR
C     ..
C     .. External Subroutines ..
      EXTERNAL MPABS,MPADD,MPCHK,MPCIM,MPCLR,MPCMR,MPDIV,MPDIVI,MPEUL,
     +         MPEXP,MPLN,MPMUL,MPMULI,MPMULQ,MPOVFL,MPREC,MPSTR,MPSUB
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ALOG,FLOAT,INT,MAX0,MIN0
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(10,38)
      XS = X(1)
      IF (XS.NE.0) GO TO 10
      WRITE (LUN,FMT=9000)
C EI(0) IS UNDEFINED, TREAT AS OVERFLOW
      CALL MPOVFL(Y)
      RETURN
C SAVE T ETC
   10 TS = T
      TM2 = (11*T+19)/10
      TM = (6*T+9)/5
      B2 = 2*MAX0(B,64)
      I = 0
C ALLOW SPACE FOR MPEUL
      I2 = 5*TM2 + 19
      I3 = I2 + TM + 2
      I4 = I3 + TM + 2
C CLEAR DIGITS OF R(I3)
      CALL MPCLR(R(I3),T+1)
      CALL MPABS(X,R(I3))
      RB = FLOAT(B)
      RT = FLOAT(T)*ALOG(RB)
C SEE IF ABS(X) LARGE ENOUGH TO USE ASYMPTOTIC SERIES
      IF (MPCMPR(R(I3),RT).GT.0) GO TO 70
C SEE IF X NEGATIVE AND CONTINUED FRACTION USABLE
C THE CONSTANT 0.1 WAS DETERMINED EMPIRICALLY AND MAY BE
C DECREASED (BUT NOT INCREASED) IF DESIRED.
      IF ((XS.LT.0) .AND. (MPCMPR(R(I3),0.1*RT).GT.0)) GO TO 100
C USE POWER SERIES HERE, BUT NEED TO INCREASE T IF X NEGATIVE
C TO COMPENSATE FOR CANCELLATION.
      T = T + 1
      TS1 = T
      TS2 = T
      IF (XS.GT.0) GO TO 20
      CALL MPCMR(R(I3),RAX)
C IF X NEGATIVE RESULT ABOUT B**(-TD) AND TERMS ABOUT B**TD SO
C NEED UP TO 2*TD EXTRA DIGITS TO COMPENSATE FOR CANCELLATION
      TD = 1 + INT(RAX/ALOG(RB))
      TS2 = T + TD
      TS1 = MIN0(TS2+TD,TM)
      TS2 = MIN0(TS2,TM2)
C CLEAR TRAILING DIGITS OF R(I2) AND R(I3)
      CALL MPCLR(R(I2),TS1)
      CALL MPCLR(R(I3),TS1)
C USE TS2 DIGITS FOR LN AND EULERS CONSTANT COMPUTATION
      T = TS2
C NOW PREPARE TO SUM POWER SERIES
   20 CALL MPLN(R(I3),R(I4))
C MPEI COULD BE SPEEDED UP IF EULERS CONSTANT WERE
C PRECOMPUTED AND SAVED
      CALL MPEUL(R(I2))
      CALL MPADD(R(I2),R(I4),R(I2))
C NOW USE TS1 DIGITS FOR SUMMING POWER SERIES
      T = TS1
C RESTORE SIGN OF R(I3)
      R(I3) = XS
      CALL MPADD(R(I2),R(I3),R(I2))
      CALL MPSTR(R(I3),R(I4))
C LOOP TO SUM POWER SERIES, REDUCING T IF POSSIBLE
   30 IF (XS.GE.0) T = TS1 + 2 + R(I4+1) - R(I2+1)
      IF ((XS.LT.0) .AND. (R(I4+1).LE.0)) T = TS2 + 2 + R(I4+1)
      T = MIN0(T,TS1)
      IF (T.LE.2) GO TO 60
      CALL MPMUL(R(I3),R(I4),R(I4))
      I = I + 1
C IF I LARGE NEED TO SPLIT UP CALL TO MPMULQ
      IF (I.GE.B2) GO TO 40
      CALL MPMULQ(R(I4),I, (I+1)**2,R(I4))
      GO TO 50

   40 CALL MPMULQ(R(I4),I,I+1,R(I4))
      CALL MPDIVI(R(I4),I+1,R(I4))
C RESTORE T FOR ADDITION
   50 T = TS1
      CALL MPADD(R(I2),R(I4),R(I2))
      IF (R(I4).NE.0) GO TO 30
C RESTORE T, MOVE RESULT AND RETURN
   60 T = TS
      CALL MPSTR(R(I2),Y)
      RETURN
C HERE WE CAN USE ASYMPTOTIC SERIES, AND NO NEED TO INCREASE T
   70 CALL MPREC(X,R(I3))
C MPEXP GIVES ERROR MESSAGE IF X TOO LARGE HERE
      CALL MPEXP(X,Y)
      IF (Y(1).EQ.0) RETURN
      CALL MPMUL(Y,R(I3),Y)
      CALL MPSTR(Y,R(I2))
C LOOP TO SUM ASYMPTOTIC SERIES, REDUCING T IF POSSIBLE
   80 T = TS + 2 + R(I2+1) - Y(2)
C RETURN IF TERMS SMALL ENOUGH TO BE NEGLIGIBLE
      IF (T.LE.2) GO TO 90
      T = MIN0(T,TS)
      CALL MPSTR(R(I2),R(I4))
      CALL MPMUL(R(I2),R(I3),R(I2))
      I = I + 1
      CALL MPMULI(R(I2),I,R(I2))
C RETURN IF TERMS INCREASING
      IF (MPCMPA(R(I2),R(I4)).GE.0) GO TO 90
C RESTORE T FOR ADDITION
      T = TS
      CALL MPADD(Y,R(I2),Y)
      IF (R(I2).NE.0) GO TO 80
C RESTORE T AND RETURN
   90 T = TS
      RETURN
C HERE 0.1*T*LN(B) .LT. -X .LE T*LN(B), SO USE CONTINUED FRACTION.
  100 CALL MPCMR(X,RX)
      C = -RX
      CP = 1.0
      K = T
      J = 0
C USE FORWARD RECURRENCE WITH SINGLE-PRECISION TO FIND HOW
C MANY TERMS NEEDED FOR FULL MP ACCURACY.
  110 J = J + 1
      CP = CP + C/FLOAT(J)
      C = C - RX*CP
C SCALE TO AVOID OVERFLOW
  120 IF (CP.LT.RB) GO TO 110
      C = C/RB
      CP = CP/RB
      K = K - 2
      IF (K.GT.0) GO TO 120
C NOW USE BACKWARD RECURRENCE WITH MP ARITHMETIC
      CALL MPCIM(1,R(I2))
  130 CALL MPDIVI(R(I3),J,R(I4))
      CALL MPADD(R(I2),R(I4),R(I2))
      CALL MPMUL(X,R(I2),R(I4))
      CALL MPSUB(R(I3),R(I4),R(I3))
C SCALE TO AVOID OVERFLOW
      R(I2+1) = R(I2+1) - R(I3+1)
      R(I3+1) = 0
      J = J - 1
      IF (J.GT.0) GO TO 130
      CALL MPDIV(R(I2),R(I3),R(I2))
      CALL MPEXP(X,Y)
      CALL MPMUL(Y,R(I2),Y)
      Y(1) = -Y(1)
      RETURN

 9000 FORMAT (' *** X ZERO IN CALL TO MPEI ***')
      END
      SUBROUTINE MPEPS(X)
C SETS MP X TO THE (MULTIPLE-PRECISION) MACHINE PRECISION,
C THAT IS THE SMALLEST POSITIVE NUMBER X SUCH THAT
C THE COMPUTED VALUE OF 1 + X IS GREATER THAN 1
C CHECK LEGALITY OF B, T, M, MXR AND LUN
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN0,MOD
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(1,4)
C SET SIGN AND EXPONENT
      X(1) = 1
      X(2) = 1 - T
C SET FRACTION DIGITS TO ZERO
      DO 10 I = 2,T
          X(I+2) = 0
   10 CONTINUE
C SEE IF BASE IS EVEN OR ODD
      IF (MOD(B,2).NE.0) GO TO 20
C EVEN BASE HERE SO X = 0.5*B**(1-T)
      X(3) = B/2
      RETURN
C ODD BASE HERE, SET X SLIGHTLY LARGER (NOTE THAT
C FOUR GUARD DIGITS ARE USED IN MPADD)
   20 I = 1
   30 X(I+2) = B/2
      I = I + 1
      IF (I.LT.MIN0(4,T)) GO TO 30
      X(I+2) = B/2 + 1
      RETURN

      END
      LOGICAL FUNCTION MPEQ(X,Y)
C RETURNS LOGICAL VALUE OF (X .EQ. Y) FOR MP X AND Y.
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. External Functions ..
      INTEGER MPCOMP
      EXTERNAL MPCOMP
C     ..
      MPEQ = (MPCOMP(X,Y).EQ.0)
      RETURN

      END
      SUBROUTINE MPERF(X,Y)
C RETURNS Y = ERF(X) = SQRT(4/PI)*(INTEGRAL FROM 0 TO X OF
C EXP(-U**2) DU) FOR MP X AND Y,  SPACE = 5T+12.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      REAL RX
      INTEGER I2,IER,TS,XS
C     ..
C     .. External Functions ..
      INTEGER MPCMPR
      EXTERNAL MPCMPR
C     ..
C     .. External Subroutines ..
      EXTERNAL MPABS,MPADDI,MPCHK,MPCIM,MPCLR,MPCMR,MPERF2,MPERF3,MPEXP,
     +         MPMUL,MPMULI,MPPI,MPROOT,MPSTR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ALOG,FLOAT,INT,MAX0,MIN0,SQRT
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      XS = X(1)
      IF (XS.NE.0) GO TO 10
C ERF(0) = 0
      Y(1) = 0
      RETURN
C CHECK LEGALITY OF B, T, M, LUN AND MXR
   10 CALL MPCHK(5,12)
      I2 = 4*T + 11
      CALL MPABS(X,R(I2))
      IF (MPCMPR(R(I2),SQRT(FLOAT(T)*ALOG(FLOAT(B)))).LT.0) GO TO 20
C HERE ABS(X) SO LARGE THAT ERF(X) = +-1 TO FULL ACCURACY
      CALL MPCIM(XS,Y)
      RETURN
C HERE SAVE T AND TRY USING ASYMPTOTIC SERIES
   20 TS = T
      CALL MPCMR(X,RX)
C CAN POSSIBLY REDUCE T TEMPORARILY
      IF (B.GE.64) T = MIN0(TS,MAX0(4,T-INT(RX*RX/ALOG(FLOAT(B)))))
C TRY ASYMPTOTIC SERIES
      CALL MPERF3(R(I2),R(I2),0,IER)
      IF (IER.EQ.0) GO TO 30
C ASYMPTOTIC SERIES INSUFFICIENT, SO USE POWER SERIES
C WITH ONE GUARD DIGIT.  SPACE REQUIRED BY MPERF2 IS
C ONLY 3(T+1)+8 = 3T+11 AS ABS(X) SMALL
      T = TS
      CALL MPCLR(R(I2),T+1)
      CALL MPSTR(X,R(I2))
      T = T + 1
      CALL MPERF2(R(I2),R(I2))
C NOW RESTORE T
      T = TS
C IN BOTH CASES MULTIPLY BY SQRT(4/PI)*EXP(-X**2)
   30 CALL MPMUL(X,X,Y)
      Y(1) = -Y(1)
      CALL MPEXP(Y,Y)
      CALL MPMUL(Y,R(I2),R(I2))
      CALL MPPI(Y)
      CALL MPROOT(Y,-2,Y)
      CALL MPMUL(Y,R(I2),R(I2))
      IF (IER.EQ.0) GO TO 40
C USED POWER SERIES SO CAN RETURN
      CALL MPMULI(R(I2),2,Y)
      RETURN
C USED ASYMPTOTIC SERIES SO SUBTRACT FROM 1
   40 CALL MPMULI(R(I2),-2,R(I2))
      T = TS
      CALL MPADDI(R(I2),1,Y)
      Y(1) = XS*Y(1)
      RETURN

      END
      SUBROUTINE MPERFC(X,Y)
C RETURNS Y = ERFC(X) = 1 - ERF(X) FOR MP NUMBERS X AND Y,
C USING MPERF AND MPERF3.   SPACE = 12T+26
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      REAL RX
      INTEGER I2,IER,TS,TS2
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADDI,MPCHK,MPCLR,MPCMR,MPERF,MPERF3,MPEXP,MPMUL,MPMULI,
     +         MPPI,MPROOT,MPSTR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ALOG,FLOAT,INT,MIN0
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      IF (X(1).GT.0) GO TO 10
C FOR X .LE. 0 NO LOSS OF ACCURACY IN USING ERF(X)
      CALL MPERF(X,Y)
      Y(1) = -Y(1)
      CALL MPADDI(Y,1,Y)
      RETURN
C CHECK LEGALITY OF B, T, M, LUN AND MXR
   10 CALL MPCHK(12,26)
      TS = T
      TS2 = 2*T + 2
      I2 = 5*TS2 + 13
C TRY ASYMPTOTIC SERIES
      CALL MPERF3(X,R(I2),0,IER)
      IF (IER.NE.0) GO TO 20
C ASYMPTOTIC SERIES WORKED, SO MULTIPLY BY
C SQRT(4/PI)*EXP(-X**2) AND RETURN
      CALL MPMUL(X,X,Y)
      Y(1) = -Y(1)
      CALL MPEXP(Y,Y)
      CALL MPMUL(Y,R(I2),R(I2))
      CALL MPPI(Y)
      CALL MPROOT(Y,-2,Y)
      CALL MPMUL(Y,R(I2),R(I2))
      CALL MPMULI(R(I2),2,Y)
      RETURN
C HERE ASYMPTOTIC SERIES INACCURATE SO HAVE TO
C USE MPERF, INCREASING PRECISION TO COMPENSATE FOR
C CANCELLATION.  AN ALTERNATIVE METHOD (POSSIBLY FASTER) IS
C TO USE THE CONTINUED FRACTION FOR EXP(X**2)*ERFC(X).
   20 CALL MPCMR(X,RX)
C CLEAR DIGITS OF R(I2)
      CALL MPCLR(R(I2),TS2)
C MOVE X TO R(I2)
      CALL MPSTR(X,R(I2))
C COMPUTE NEW T FOR MPERF COMPUTATION
      T = MIN0(TS2,TS+2+INT(RX*RX/ALOG(FLOAT(B))))
      CALL MPERF(R(I2),R(I2))
      R(I2) = -R(I2)
      CALL MPADDI(R(I2),1,R(I2))
C RESTORE T AND MOVE RESULT TO Y
      T = TS
      CALL MPSTR(R(I2),Y)
      RETURN

      END
      SUBROUTINE MPERF2(X,Y)
C RETURNS Y = EXP(X**2)*(INTEGRAL FROM 0 TO X OF EXP(-U*U) DU)
C FOR MP NUMBERS X AND Y, USING THE POWER SERIES FOR SMALL X,
C AND MPEXP FOR LARGE X.  SPACE = 5T+12 (OR 3T+8 FOR
C SMALL X).   CALLED BY MPERF.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I,I2,I3,I4,TS,XS
C     ..
C     .. External Functions ..
      INTEGER MPCMPR
      EXTERNAL MPCMPR
C     ..
C     .. External Subroutines ..
      EXTERNAL MPABS,MPADD,MPCHK,MPDIVI,MPEXP,MPMUL,MPMULI,MPPI,MPSQRT,
     +         MPSTR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ALOG,FLOAT,MIN0,SQRT
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      IF (X(1).NE.0) GO TO 10
C RETURN 0 IF X .EQ. 0
      Y(1) = 0
      RETURN
C CHECK LEGALITY OF B, T, M, LUN AND MXR
   10 CALL MPCHK(3,8)
      I2 = T + 5
      I3 = I2 + T + 2
      CALL MPABS(X,R(I3))
      IF (MPCMPR(R(I3),SQRT(FLOAT(T)*ALOG(FLOAT(B)))).GT.0) GO TO 40
C USE THE POWER SERIES HERE
      CALL MPSTR(X,Y)
      CALL MPMUL(X,X,R(I2))
      CALL MPMULI(R(I2),2,R(I2))
      CALL MPSTR(X,R(I3))
      TS = T
      I = 1
C LOOP TO SUM SERIES, REDUCING T IF POSSIBLE
   20 T = TS + 2 + R(I3+1) - Y(2)
      IF (T.LE.2) GO TO 30
      T = MIN0(T,TS)
      CALL MPMUL(R(I2),R(I3),R(I3))
      I = I + 2
      CALL MPDIVI(R(I3),I,R(I3))
C RESTORE T FOR ADDITION
      T = TS
      CALL MPADD(Y,R(I3),Y)
      IF (R(I3).NE.0) GO TO 20
C RESTORE T AND RETURN
   30 T = TS
      RETURN
C HERE ABS(X) LARGE, SO INTEGRAL IS +-SQRT(PI/4)
   40 CALL MPCHK(5,12)
      I4 = 4*T + 11
      CALL MPMUL(X,X,R(I4))
C IF ABS(X) TOO LARGE MPEXP GIVES ERROR MESSAGE
      CALL MPEXP(R(I4),R(I4))
      XS = X(1)
      CALL MPPI(Y)
      CALL MPSQRT(Y,Y)
      CALL MPDIVI(Y,2*XS,Y)
      CALL MPMUL(Y,R(I4),Y)
      RETURN

      END
      SUBROUTINE MPERF3(X,Y,IND,ERROR)
C IF IND .EQ. 0, RETURNS Y = EXP(X**2)*(INTEGRAL FROM X TO
C                INFINITY OF EXP(-U**2) DU),
C IF IND .NE. 0, RETURNS Y = EXP(-X**2)*(INTEGRAL FROM 0 TO
C                X OF EXP(U**2) DU),
C IN BOTH CASES USING THE ASYMPTOTIC SERIES.
C X AND Y ARE MP NUMBERS, IND AND ERROR ARE INTEGERS.
C ERROR IS RETURNED AS 0 IF X IS LARGE ENOUGH FOR
C THE ASYMPTOTIC SERIES TO GIVE FULL ACCURACY,
C OTHERWISE ERROR IS RETURNED AS 1 AND Y AS ZERO.
C THE CONDITION ON X FOR ERROR .EQ. 0 IS APPROXIMATELY THAT
C X .GT. SQRT(T*LOG(B)).
C CALLED BY MPERF, MPERFC AND MPDAW, SPACE = 4T+10
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER ERROR,IND
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I,I2,I3,IE,TS
C     ..
C     .. External Functions ..
      INTEGER MPCMPR
      EXTERNAL MPCMPR
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADD,MPCHK,MPDIVI,MPMUL,MPMULI,MPREC,MPSTR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ALOG,FLOAT,MIN0,SQRT
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      ERROR = 0
C CHECK LEGALITY OF B, T, M, LUN AND MXR
      CALL MPCHK(4,10)
      TS = T
C CHECK THAT CAN GET AT LEAST T-2 DIGITS ACCURACY
      IF (MPCMPR(X,SQRT(FLOAT(T-2)*ALOG(FLOAT(B)))).GT.0) GO TO 30
C HERE X IS TOO SMALL FOR ASYMPTOTIC SERIES TO GIVE
C FULL ACCURACY, SO RETURN WITH ERROR .EQ. 1
   10 Y(1) = 0
      ERROR = 1
   20 T = TS
      RETURN

   30 CALL MPREC(X,Y)
      I2 = T + 5
      I3 = I2 + T + 2
      CALL MPMUL(Y,Y,R(I2))
      CALL MPDIVI(R(I2),2,R(I2))
      IF (IND.EQ.0) R(I2) = -R(I2)
      CALL MPDIVI(Y,2,Y)
      CALL MPSTR(Y,R(I3))
      I = 1
C LOOP TO SUM SERIES, REDUCING T IF POSSIBLE
   40 IE = R(I3+1)
      T = TS + 2 + IE - Y(2)
      IF (T.LE.2) GO TO 20
      T = MIN0(T,TS)
      CALL MPMUL(R(I2),R(I3),R(I3))
      CALL MPMULI(R(I3),I,R(I3))
      I = I + 2
C RESTORE T FOR ADDITION
      T = TS
C CHECK IF TERMS ARE GETTING LARGER - IF SO X IS TOO
C SMALL FOR ASYMPTOTIC SERIES TO BE ACCURATE
      IF (R(I3+1).GT.IE) GO TO 10
      CALL MPADD(Y,R(I3),Y)
      IF (R(I3).NE.0) GO TO 40
      GO TO 20

      END
      SUBROUTINE MPERR
C THIS ROUTINE IS CALLED WHEN A FATAL ERROR CONDITION IS
C ENCOUNTERED, AND AFTER A MESSAGE HAS BEEN WRITTEN ON
C LOGICAL UNIT LUN.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      WRITE (LUN,FMT=9000)
C AT PRESENT JUST STOP, BUT COULD DUMP B, T, ETC. HERE.
C ACTION COULD EASILY BE CONTROLLED BY A FLAG IN LABELLED COMMON.
C ANSI VERSION USES STOP, UNIVAC 1108 VERSION USES
C RETURN 0 IN ORDER TO GIVE A TRACE-BACK.
C FOR DEBUGGING PURPOSES IT MAY BE USEFUL SIMPLY TO
C RETURN HERE.  MOST MP ROUTINES RETURN WITH RESULT
C ZERO AFTER CALLING MPERR.
      STOP

 9000 FORMAT (' *** EXECUTION TERMINATED BY CALL TO MPERR',' IN MP VER',
     +       'SION 780802 ***')
      END
      SUBROUTINE MPEUL(G)
C RETURNS MP G = EULERS CONSTANT (GAMMA = 0.57721566...)
C TO ALMOST FULL MULTIPLE-PRECISION ACCURACY.
C THE METHOD IS BASED ON BESSEL FUNCTION IDENTITIES AND WAS
C DISCOVERED BY EDWIN MC MILLAN AND R. BRENT.  IT IS FASTER THAN THE
C METHOD OF SWEENEY (MATH. COMP. 17, 1963, 170) USED IN EARLIER
C VERSIONS OF MPEUL.  TIME O(T**2),  SPACE = 5T+18.
C CHECK LEGALITY OF B, T, M ETC
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER G(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      REAL RG
      INTEGER B2,I2,I3,I4,I5,K,N,N2,TS
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADD,MPCHK,MPCIM,MPCMR,MPDIVI,MPERR,MPLNI,MPMUL,MPMULI,
     +         MPMULQ,MPREC,MPSTR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ALOG,FLOAT,INT,MAX0,MIN0
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(5,18)
      B2 = 2*MAX0(B,64)
C USE ONE GUARD DIGIT TO GIVE ALMOST FULL-PRECISION RESULT
      TS = T
      T = T + 1
      I2 = T + 6
      I3 = I2 + T + 2
      I4 = I3 + T + 2
      I5 = I4 + T + 2
C COMPUTE N SO TRUNCATION ERROR O(EXP(-4*N)) IS O(B**(-T))
      N = INT(0.25*FLOAT(TS)*ALOG(FLOAT(B))) + 1
      IF (N.LE.B2) N2 = N*N
      CALL MPLNI(N,R(I4))
      R(I4) = -R(I4)
      CALL MPSTR(R(I4),R(I3))
      CALL MPCIM(1,R(I2))
      CALL MPSTR(R(I2),R(I5))
      K = 0
C MAIN LOOP STARTS HERE
   10 K = K + 1
C REDUCE T HERE IF POSSIBLE
      IF (K.GT.N) T = MIN0(T,T+1+R(I2+1)-R(I5+1))
C TEST FOR CONVERGENCE
      IF (T.LT.2) GO TO 40
C SPLIT UP CALLS TO MPMULQ IF NECESSARY
      IF ((N.GT.B2) .OR. (K.GT.B2)) GO TO 20
      CALL MPMULQ(R(I2),N2,K*K,R(I2))
      CALL MPMULQ(R(I3),N2,K,R(I3))
      GO TO 30
C HERE CALLS TO MPMULQ ARE SPLIT UP
   20 CALL MPMULQ(R(I2),N,K,R(I2))
      CALL MPMULQ(R(I2),N,K,R(I2))
      CALL MPMULQ(R(I3),N,K,R(I3))
      CALL MPMULI(R(I3),N,R(I3))
   30 CALL MPADD(R(I3),R(I2),R(I3))
      CALL MPDIVI(R(I3),K,R(I3))
C INCREASE T HERE
      T = TS + 1
      CALL MPADD(R(I5),R(I2),R(I5))
      CALL MPADD(R(I4),R(I3),R(I4))
C END OF MAIN LOOP
      IF (R(I2).NE.0) GO TO 10
C RESTORE T AND COMPUTE FINAL QUOTIENT
C R(I4) (EXCEPT LAST DIGIT) WILL BE OVERWRITTEN BY MPREC
   40 T = TS
      CALL MPSTR(R(I4),G)
      T = TS + 1
      CALL MPREC(R(I5),R(I5))
      T = TS
      CALL MPSTR(G,R(I4))
      T = TS + 1
      CALL MPMUL(R(I4),R(I5),R(I4))
      T = TS
      CALL MPSTR(R(I4),G)
C CHECK REASONABLENESS OF RESULT, ASSUMING B AND T LARGE
C ENOUGH TO GIVE ERROR LESS THAN 0.01
      CALL MPCMR(G,RG)
      IF (ABS(RG-0.5772).LT.0.01) RETURN
      WRITE (LUN,FMT=9000)
C THE FOLLOWING MESSAGE MAY INDICATE THAT
C B**(T-1) IS TOO SMALL.
      CALL MPERR
      RETURN

 9000 FORMAT (' *** ERROR OCCURRED IN MPEUL, RESULT INCORRECT ***')
      END
      SUBROUTINE MPEXP(X,Y)
C RETURNS Y = EXP(X) FOR MP X AND Y.
C EXP OF INTEGER AND FRACTIONAL PARTS OF X ARE COMPUTED
C SEPARATELY.  SEE ALSO COMMENTS IN MPEXP1.
C TIME IS O(SQRT(T)M(T)).
C DIMENSION OF R MUST BE AT LEAST 4T+10 IN CALLING PROGRAM
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      REAL RLB,RX,RY
      INTEGER I,I2,I3,IE,IX,TS,TSS,XS
C     ..
C     .. External Functions ..
      INTEGER MPCMPR
      EXTERNAL MPCMPR
C     ..
C     .. External Subroutines ..
      EXTERNAL MPABS,MPADD2,MPADDI,MPCHK,MPCIM,MPCMF,MPCMI,MPCMR,MPDIVI,
     +         MPERR,MPEXP1,MPMUL,MPOVFL,MPPWR,MPUNFL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ALOG,EXP,FLOAT,MIN0
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(4,10)
      I2 = 2*T + 7
      I3 = I2 + T + 2
C CHECK FOR X = 0
      IF (X(1).NE.0) GO TO 10
      CALL MPCIM(1,Y)
      RETURN
C CHECK IF ABS(X) .LT. 1
   10 IF (X(2).GT.0) GO TO 20
C USE MPEXP1 HERE
      CALL MPEXP1(X,Y)
      CALL MPADDI(Y,1,Y)
      RETURN
C SEE IF ABS(X) SO LARGE THAT EXP(X) WILL CERTAINLY OVERFLOW
C OR UNDERFLOW.  1.01 IS TO ALLOW FOR ERRORS IN ALOG.
   20 RLB = 1.01E0*ALOG(FLOAT(B))
      IF (MPCMPR(X,-FLOAT(M+1)*RLB).GE.0) GO TO 40
C UNDERFLOW SO CALL MPUNFL AND RETURN
   30 CALL MPUNFL(Y)
      RETURN

   40 IF (MPCMPR(X,FLOAT(M)*RLB).LE.0) GO TO 60
C OVERFLOW HERE
   50 WRITE (LUN,FMT=9000)
      CALL MPOVFL(Y)
      RETURN
C NOW SAFE TO CONVERT X TO REAL
   60 CALL MPCMR(X,RX)
C SAVE SIGN AND WORK WITH ABS(X)
      XS = X(1)
      CALL MPABS(X,R(I3))
C IF ABS(X) .GT. M POSSIBLE THAT INT(X) OVERFLOWS,
C SO DIVIDE BY 32.
      IF (ABS(RX).GT.FLOAT(M)) CALL MPDIVI(R(I3),32,R(I3))
C GET FRACTIONAL AND INTEGER PARTS OF ABS(X)
      CALL MPCMI(R(I3),IX)
      CALL MPCMF(R(I3),R(I3))
C ATTACH SIGN TO FRACTIONAL PART AND COMPUTE EXP OF IT
      R(I3) = XS*R(I3)
      CALL MPEXP1(R(I3),Y)
      CALL MPADDI(Y,1,Y)
C COMPUTE E-2 OR 1/E USING TWO EXTRA DIGITS IN CASE ABS(X) LARGE
C (BUT ONLY ONE EXTRA DIGIT IF T .LT. 4)
      TSS = T
      TS = T + 2
      IF (T.LT.4) TS = T + 1
      T = TS
      I2 = T + 5
      I3 = I2 + T + 2
      R(I3) = 0
      CALL MPCIM(XS,R(I2))
      I = 1
C LOOP FOR E COMPUTATION. DECREASE T IF POSSIBLE.
   70 T = MIN0(TS,TS+2+R(I2+1))
      IF (T.LE.2) GO TO 80
      I = I + 1
      CALL MPDIVI(R(I2),I*XS,R(I2))
      T = TS
      CALL MPADD2(R(I3),R(I2),R(I3),R(I2),0)
      IF (R(I2).NE.0) GO TO 70
C RAISE E OR 1/E TO POWER IX
   80 T = TS
      IF (XS.GT.0) CALL MPADDI(R(I3),2,R(I3))
      CALL MPPWR(R(I3),IX,R(I3))
C RESTORE T NOW
      T = TSS
C MULTIPLY EXPS OF INTEGER AND FRACTIONAL PARTS
      CALL MPMUL(Y,R(I3),Y)
C MUST CORRECT RESULT IF DIVIDED BY 32 ABOVE.
      IF ((ABS(RX).LE.FLOAT(M)) .OR. (Y(1).EQ.0)) GO TO 100
      DO 90 I = 1,5
C SAVE EXPONENT TO AVOID OVERFLOW IN MPMUL
          IE = Y(2)
          Y(2) = 0
          CALL MPMUL(Y,Y,Y)
          Y(2) = Y(2) + 2*IE
C CHECK FOR UNDERFLOW AND OVERFLOW
          IF (Y(2).LT. (-M)) GO TO 30
          IF (Y(2).GT.M) GO TO 50
   90 CONTINUE
C CHECK THAT RELATIVE ERROR LESS THAN 0.01 UNLESS ABS(X) LARGE
C (WHEN EXP MIGHT OVERFLOW OR UNDERFLOW)
  100 IF (ABS(RX).GT.10.0) RETURN
      CALL MPCMR(Y,RY)
      IF (ABS(RY-EXP(RX)).LT. (0.01*RY)) RETURN
      WRITE (LUN,FMT=9010)
C THE FOLLOWING MESSAGE MAY INDICATE THAT
C B**(T-1) IS TOO SMALL, OR THAT M IS TOO SMALL SO THE
C RESULT UNDERFLOWED.
      CALL MPERR
      RETURN

 9000 FORMAT (' *** OVERFLOW IN SUBROUTINE MPEXP ***')
 9010 FORMAT (' *** ERROR OCCURRED IN MPEXP, RESULT INCORRECT ***')
      END
      INTEGER FUNCTION MPEXPA(X)
C RETURNS THE EXPONENT OF THE MP NUMBER X
C (OR LARGEST NEGATIVE EXPONENT IF X IS ZERO).
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      MPEXPA = -M
C RETURN -M IF X ZERO, X(2) OTHERWISE
      IF (X(1).NE.0) MPEXPA = X(2)
      RETURN

      END
      SUBROUTINE MPEXPB(I,X)
C SETS EXPONENT OF MP NUMBER X TO I UNLESS X IS ZERO
C (WHEN EXPONENT IS UNCHANGED).
C X MUST BE A VALID MP NUMBER (EITHER ZERO OR NORMALIZED).
C RETURN IF X IS ZERO
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER I
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. External Subroutines ..
      EXTERNAL MPERR,MPOVFL,MPUNFL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IABS
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      IF (X(1).EQ.0) RETURN
C CHECK FOR VALID MP SIGN AND LEADING DIGIT
      IF ((IABS(X(1)).LE.1) .AND. (X(3).GT.0) .AND.
     +    (X(3).LT.B)) GO TO 10
      WRITE (LUN,FMT=9000)
      CALL MPERR
      X(1) = 0
      RETURN
C SET EXPONENT OF X TO I
   10 X(2) = I
C CHECK FOR OVERFLOW AND UNDERFLOW
      IF (I.GT.M) CALL MPOVFL(X)
      IF (I.LT. (-M)) CALL MPUNFL(X)
      RETURN

 9000 FORMAT (' *** X NOT VALID MP NUMBER IN CALL TO MPEXPB ***')
      END
      SUBROUTINE MPEXP1(X,Y)
C ASSUMES THAT X AND Y ARE MP NUMBERS,  -1 .LT. X .LT. 1.
C RETURNS Y = EXP(X) - 1 USING AN O(SQRT(T).M(T)) ALGORITHM
C DESCRIBED IN - R. P. BRENT, THE COMPLEXITY OF MULTIPLE-
C PRECISION ARITHMETIC (IN COMPLEXITY OF COMPUTATIONAL PROBLEM
C SOLVING, UNIV. OF QUEENSLAND PRESS, BRISBANE, 1976, 126-165).
C ASYMPTOTICALLY FASTER METHODS EXIST, BUT ARE NOT USEFUL
C UNLESS T IS VERY LARGE. SEE COMMENTS TO MPATAN AND MPPIGL.
C DIMENSION OF R IN CALLING PROGRAM MUST BE AT LEAST 3T+8
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      REAL RLB
      INTEGER I,I2,I3,IB,IC,Q,TS
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADD2,MPADDI,MPCHK,MPDIVI,MPERR,MPMUL,MPSTR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ALOG,FLOAT,INT,MIN0,SQRT
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(3,8)
      I2 = T + 5
      I3 = I2 + T + 2
C CHECK FOR X = 0
      IF (X(1).NE.0) GO TO 20
   10 Y(1) = 0
      RETURN
C CHECK THAT ABS(X) .LT. 1
   20 IF (X(2).LE.0) GO TO 30
      WRITE (LUN,FMT=9000)
      CALL MPERR
      GO TO 10

   30 CALL MPSTR(X,R(I2))
      RLB = ALOG(FLOAT(B))
C COMPUTE APPROXIMATELY OPTIMAL Q (AND DIVIDE X BY 2**Q)
      Q = INT(SQRT(0.48E0*FLOAT(T)*RLB)+1.44E0*FLOAT(X(2))*RLB)
C HALVE Q TIMES
      IF (Q.LE.0) GO TO 50
      IB = 4*B
      IC = 1
      DO 40 I = 1,Q
          IC = 2*IC
          IF ((IC.LT.IB) .AND. (IC.NE.B) .AND. (I.LT.Q)) GO TO 40
          CALL MPDIVI(R(I2),IC,R(I2))
          IC = 1
   40 CONTINUE
   50 IF (R(I2).EQ.0) GO TO 10
      CALL MPSTR(R(I2),Y)
      CALL MPSTR(R(I2),R(I3))
      I = 1
      TS = T
C SUM SERIES, REDUCING T WHERE POSSIBLE
   60 T = TS + 2 + R(I3+1) - Y(2)
      IF (T.LE.2) GO TO 70
      T = MIN0(T,TS)
      CALL MPMUL(R(I2),R(I3),R(I3))
      I = I + 1
      CALL MPDIVI(R(I3),I,R(I3))
      T = TS
      CALL MPADD2(R(I3),Y,Y,Y,0)
      IF (R(I3).NE.0) GO TO 60
   70 T = TS
      IF (Q.LE.0) RETURN
C APPLY (X+1)**2 - 1 = X(2 + X) FOR Q ITERATIONS
      DO 80 I = 1,Q
          CALL MPADDI(Y,2,R(I2))
          CALL MPMUL(R(I2),Y,Y)
   80 CONTINUE
      RETURN

 9000 FORMAT (' *** ABS(X) NOT LESS THAN 1 IN CALL TO MPEXP1 ***')
      END
      SUBROUTINE MPEXT(I,J,X)
C ROUTINE CALLED BY MPDIV AND MPSQRT TO ENSURE THAT
C RESULTS ARE REPRESENTED EXACTLY IN T-2 DIGITS IF THEY
C CAN BE.  X IS AN MP NUMBER, I AND J ARE INTEGERS.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER I,J
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER Q,S
C     ..
C     .. External Subroutines ..
      EXTERNAL MPMULI
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      IF ((X(1).EQ.0) .OR. (T.LE.2) .OR. (I.EQ.0)) RETURN
C COMPUTE MAXIMUM POSSIBLE ERROR IN THE LAST PLACE
      Q = (J+1)/I + 1
      S = B*X(T+1) + X(T+2)
      IF (S.GT.Q) GO TO 10
C SET LAST TWO DIGITS TO ZERO
      X(T+1) = 0
      X(T+2) = 0
      RETURN

   10 IF ((S+Q).LT. (B*B)) RETURN
C ROUND UP HERE
      X(T+1) = B - 1
      X(T+2) = B
C NORMALIZE X (LAST DIGIT B IS OK IN MPMULI)
      CALL MPMULI(X,1,X)
      RETURN

      END
      SUBROUTINE MPGAM(X,Y)
C COMPUTES MP Y = GAMMA(X) FOR MP ARGUMENT X, USING
C MPGAMQ IF ABS(X) .LE. 100 AND 240*X IS AN INTEGER,
C OTHERWISE USING MPLNGM.  SPACE REQUIRED IS THE SAME
C AS FOR MPLNGM (THOUGH ONLY 9T+20 IF MPGAMQ IS USED).
C TIME IS O(T**3).
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I2,I3,IX,KT,N
C     ..
C     .. External Functions ..
      INTEGER MPCMPI,MPCMPR
      EXTERNAL MPCMPI,MPCMPR
C     ..
C     .. External Subroutines ..
      EXTERNAL MPABS,MPADDI,MPADDQ,MPCHK,MPCMF,MPCMI,MPDIV,MPDIVI,MPEXP,
     +         MPGAMQ,MPLNGM,MPMUL,MPMULI,MPOVFL,MPPI,MPREC,MPSIN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ALOG,FLOAT,IABS,MOD
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(9,20)
      I2 = 7*T + 17
      I3 = I2 + T + 2
      CALL MPABS(X,R(I3))
      IF (MPCMPI(R(I3),100).GT.0) GO TO 20
C HERE ABS(X) .LE. 100, SEE IF 240*X IS ALMOST AN INTEGER
      CALL MPMULI(X,240,R(I3))
      CALL MPCMI(R(I3),IX)
C COMPARE WITH IX AND IX+1 BECAUSE R(I3) COULD BE JUST
C BELOW AN INTEGER.
      DO 10 KT = 1,2
          CALL MPADDI(R(I3),-IX,R(I2))
          IF ((R(I2).EQ.0) .OR. (((R(I3+1)-R(I2+1)).GE. (T-1)).AND.
     +        (R(I3+2).GE.R(I2+2)))) GO TO 30
          IX = IX + 1
   10 CONTINUE
C HERE X IS LARGE OR NOT SIMPLE RATIONAL,
C CHECK IF ABS(X) VERY SMALL.
      IF (X(2).LE. (-T)) GO TO 90
C NOW CHECK SIGN OF X
   20 IF (X(1).LT.0) GO TO 40
C X IS POSITIVE SO USE MPLNGM DIRECTLY
      CALL MPLNGM(X,Y)
C SEE IF MPEXP WILL GIVE OVERFLOW
      IF (MPCMPR(Y,FLOAT(M)*ALOG(FLOAT(B))).GE.0) GO TO 70
C SAFE TO CALL MPEXP HERE EXCEPT IN VERY UNLIKELY CIRCUMSTANCES
      CALL MPEXP(Y,Y)
      RETURN
C X = IX/240 SO USE MPGAMQ UNLESS X ZERO OR NEGATIVE INTEGER.
   30 IF ((IX.LE.0) .AND. (MOD(IABS(IX),240).EQ.0)) GO TO 50
      CALL MPGAMQ(IX,240,Y)
      RETURN
C HERE X IS NEGATIVE, SO USE REFLECTION FORMULA
   40 CALL MPABS(X,Y)
C SUBTRACT EVEN INTEGER TO AVOID ERRORS NEAR POLES
      CALL MPDIVI(Y,2,R(I3))
      CALL MPCMF(R(I3),R(I3))
      CALL MPMULI(R(I3),2,R(I3))
      CALL MPADDQ(R(I3),1,2,R(I2))
      CALL MPCMI(R(I2),N)
C CHECK FOR INTEGER OVERFLOW IN MPCMI
      IF ((R(I3).NE.0) .AND. (R(I3+1).GT.0) .AND. (N.EQ.0)) GO TO 70
      CALL MPADDI(R(I3),-N,R(I3))
C NOW ABS(R(I3)) .LE. 1/2 AND SIGN DETERMINED BY N
      IF (R(I3).NE.0) GO TO 60
   50 WRITE (LUN,FMT=9000)
C TREAT AS OVERFLOW
      GO TO 80

   60 CALL MPPI(R(I2))
      CALL MPMUL(R(I3),R(I2),R(I3))
      CALL MPSIN(R(I3),R(I3))
      CALL MPMUL(R(I3),Y,R(I3))
      IF (R(I3).EQ.0) GO TO 70
      CALL MPDIV(R(I2),R(I3),R(I3))
      R(I3) = - ((-1)**N)*R(I3)
C NOTE THAT MPLNGM PRESERVES R(I3), ... , R(I3+T+1)
      CALL MPLNGM(Y,Y)
      Y(1) = -Y(1)
      CALL MPEXP(Y,Y)
      CALL MPMUL(Y,R(I3),Y)
      RETURN
C HERE X WAS TOO LARGE OR TOO CLOSE TO A POLE
   70 WRITE (LUN,FMT=9010)
   80 CALL MPOVFL(Y)
      RETURN
C HERE ABS(X) IS VERY SMALL
   90 CALL MPREC(X,Y)
      RETURN

 9000 FORMAT (' *** X ZERO OR NEGATIVE INTEGER IN CALL TO MPGAM ***')
 9010 FORMAT (' *** OVERFLOW IN MPGAM ***')
      END
      SUBROUTINE MPGAMQ(I,J,X)

C RETURNS X = GAMMA (I/J) WHERE X IS MULTIPLE-PRECISION AND
C I, J ARE SMALL INTEGERS.   THE METHOD USED IS REDUCTION OF
C THE ARGUMENT TO (0, 1) AND THEN A DIRECT
C EXPANSION OF THE DEFINING INTEGRAL TRUNCATED AT A
C SUFFICIENTLY HIGH LIMIT, USING 2T DIGITS TO
C COMPENSATE FOR CANCELLATION.
C TIME IS O(T**2) IF I/J IS NOT TOO LARGE.
C IF I/J .GT. 100 (APPROXIMATELY) IT IS FASTER TO USE
C MPGAM (IF ENOUGH SPACE IS AVAILABLE).
C MPGAMQ IS VERY SLOW IF I/J IS VERY LARGE, BECAUSE
C THE RELATION GAMMA(X+1) = X*GAMMA(X) IS USED REPEATEDLY.
C MPGAMQ COULD BE SPEEDED UP BY USING THE ASYMPTOTIC SERIES OR
C CONTINUED FRACTION FOR (INTEGRAL FROM N TO INFINITY OF
C U**(I/J-1)*EXP(-U)DU).
C IF I OR J IS TOO LARGE, INTEGER OVERFLOW WILL OCCUR, AND
C THE RESULT WILL BE INCORRECT.  THIS WILL USUALLY (BUT NOT
C ALWAYS) BE DETECTED AND AN ERROR MESSAGE GIVEN.
C DIMENSION OF R IN CALLING PROGRAM AT LEAST 6T+12.
C CHECK LEGALITY OF B, T, M, MXR AND LUN
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER I,J
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I2,I3,IBT,IBTN,ID,IJ,IL,IN,IS,IS2,JS,N,TS,TS2,TS3
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADD2,MPCHK,MPCIM,MPERR,MPGCD,MPMUL,MPMULQ,MPOVFL,MPPI,
     +         MPQPWR,MPSQRT,MPSTR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ALOG,FLOAT,IABS,INT,MAX0,MIN0
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(6,12)
      IS = I
      JS = J
C LOOK AT SIGN OF J
      IF (JS) 20,10,30
   10 WRITE (LUN,FMT=9000)
      GO TO 60
C J NEGATIVE HERE SO REVERSE SIGNS OF IS AND JS
   20 IS = -IS
      JS = -JS
C NOW JS IS POSITIVE.  REDUCE TO LOWEST TERMS.
   30 CALL MPGCD(IS,JS)
      IBT = MAX0(7*B*B,32767)
      IJ = IBT/JS
C SEE IF JS = 1, 2, OR .GT. 2
      IF (JS-2) 50,40,90
C JS = 2 HERE, FOR SPEED TREAT AS SPECIAL CASE
   40 CALL MPPI(X)
      CALL MPSQRT(X,X)
      GO TO 80
C JS = 1 HERE, CHECK THAT IS IS POSITIVE
   50 IF (IS.GT.0) GO TO 70
      WRITE (LUN,FMT=9010)
C TREAT AS OVERFLOW
   60 CALL MPOVFL(X)
      RETURN
C I/J = POSITIVE INTEGER HERE
   70 CALL MPCIM(1,X)
   80 IS2 = 1
      GO TO 130
C JS .GT. 2 HERE SO REDUCE TO (0, 1)
   90 IS2 = IS - (IS/JS)*JS
      IF (IS2.LT.0) IS2 = IS2 + JS
C NOW 0 .LT. IS2 .LT. JS.   COMPUTE UPPER LIMIT OF INTEGRAL
      N = INT(FLOAT(T)*ALOG(FLOAT(B)))
      IBTN = IBT/N
      I3 = 4*T + 11
      TS = T
      TS3 = T + 2
C INCREASE T TO COMPENSATE FOR CANCELLATION
      T = 2*T
      TS2 = T
      I2 = I3 - (T+2)
      CALL MPCIM(N,R(I2))
      CALL MPSTR(R(I2),R(I3))
      IL = 0
      IN = JS - IS2
      ID = IS2
C MAIN LOOP
  100 IL = IL + 1
C IF TERMS DECREASING MAY DECREASE T
      IF (IL.GE.N) T = R(I3+1) + TS3
      T = MAX0(2,MIN0(T,TS2))
      IN = IN - JS
      ID = ID + JS
C CHECK FOR OVERFLOW HERE (ID SHOULD BE POSITIVE)
      IF (ID.LE.0) GO TO 180
C SPLIT UP CALL TO MPMULQ IF IN OR ID TOO LARGE
      IF ((IABS(IN).GT.IBTN) .OR. (ID.GT. (IBT/IL))) GO TO 110
      CALL MPMULQ(R(I3),N*IN,IL*ID,R(I3))
      GO TO 120

  110 CALL MPMULQ(R(I3),N,IL,R(I3))
      CALL MPMULQ(R(I3),IN,ID,R(I3))
  120 T = MAX0(T,TS3)
      CALL MPADD2(R(I2),R(I3),R(I2),R(I3),0)
C LOOP UNTIL EXPONENT SMALL
      IF ((R(I3).NE.0) .AND. (R(I3+1).GE. (-TS))) GO TO 100
C RESTORE T
      T = TS
      CALL MPMULQ(R(I2),JS,IS2,X)
      CALL MPQPWR(N,1,IS2-JS,JS,R(I3))
      CALL MPMUL(X,R(I3),X)
C NOW X IS GAMMA (IS2/JS), SO USE THE RECURRENCE RELATION
C REPEATEDLY TO GET GAMMA (I/J)  (SLOW IF I/J IS LARGE).
  130 IN = 1
      ID = 1
      IF (IS-IS2) 170,140,150
  140 RETURN

  150 IN = IN*IS2
      ID = ID*JS
      IS2 = IS2 + JS
      IF ((ID.LE.IJ) .AND. (IABS(IN).LE. (IBT/IABS(IS2))) .AND.
     +    (IS.NE.IS2)) GO TO 150
  160 CALL MPMULQ(X,IN,ID,X)
      GO TO 130

  170 IN = IN*JS
      ID = ID*IS
      IS = IS + JS
      IF ((IN.LE.IJ) .AND. (IABS(ID).LE. (IBT/IABS(IS))) .AND.
     +    (IS.NE.IS2)) GO TO 170
      GO TO 160
C HERE INTEGER OVERFLOW OCCURRED, J MUST HAVE BEEN TOO LARGE
  180 WRITE (LUN,FMT=9020)
      CALL MPERR
      X(1) = 0
      RETURN

 9000 FORMAT (' *** J = 0 IN CALL TO MPGAMQ ***')
 9010 FORMAT (' *** I/J = ZERO OR NEGATIVE INTEGER IN CALL',' TO MPGAM',
     +       'Q ***')
 9020 FORMAT (' *** INTEGER OVERFLOW OCCURRED,',' J TOO LARGE IN CALL ',
     +       'TO MPGAMQ ***')
      END
      SUBROUTINE MPGCD(K,L)
C RETURNS K = K/GCD AND L = L/GCD, WHERE GCD IS THE
C GREATEST COMMON DIVISOR OF K AND L.
C SAVE INPUT PARAMETERS IN LOCAL VARIABLES
C     .. Scalar Arguments ..
      INTEGER K,L
C     ..
C     .. Local Scalars ..
      INTEGER I,IS,J,JS
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IABS,MOD
C     ..
      I = K
      J = L
      IS = IABS(I)
      JS = IABS(J)
      IF (JS.EQ.0) GO TO 30
C EUCLIDEAN ALGORITHM LOOP
   10 IS = MOD(IS,JS)
      IF (IS.EQ.0) GO TO 20
      JS = MOD(JS,IS)
      IF (JS.NE.0) GO TO 10
      JS = IS
C HERE JS IS THE GCD OF I AND J
   20 K = I/JS
      L = J/JS
      RETURN
C IF J = 0 RETURN (1, 0) UNLESS I = 0, THEN (0, 0)
   30 K = 1
      IF (IS.EQ.0) K = 0
      L = 0
      RETURN

      END
      SUBROUTINE MPGCDA(X,Y,Z)
C RETURNS Z = GREATEST COMMON DIVISOR OF X AND Y.
C GCD (X, 0) = GCD (0, X) = ABS(X), GCD (X, Y) .GE. 0.
C X, Y AND Z ARE INTEGERS REPRESENTED AS MP NUMBERS,
C AND MUST SATISFY ABS(X) .LT. B**T, ABS(Y) .LT. B**T
C TIME O(T**2), SPACE = 4T+10.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*),Z(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I,I2,I3,I4,IQ,IS,J,TS
C     ..
C     .. External Functions ..
      INTEGER MPCOMP
      EXTERNAL MPCOMP
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPCIM,MPCMI,MPCMIM,MPERR,MPMULI,MPSTR,MPSUB
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IABS,MAX0,MOD
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(3,8)
      TS = T
      I2 = T + 5
      I3 = I2 + T + 2
      I4 = I3 + T + 2
      CALL MPCMIM(X,R(I2))
C CHECK THAT X EXACT INTEGER
      IF (MPCOMP(X,R(I2)).NE.0) GO TO 40
      R(I2) = IABS(R(I2))
      CALL MPCMIM(Y,R(I3))
C CHECK THAT Y EXACT INTEGER
      IF (MPCOMP(Y,R(I3)).NE.0) GO TO 40
      R(I3) = IABS(R(I3))
C CHECK FOR X OR Y ZERO
      IF (X(1).NE.0) GO TO 20
   10 T = TS
      CALL MPSTR(R(I3),Z)
      RETURN

   20 IF (Y(1).NE.0) GO TO 30
      CALL MPSTR(R(I2),Z)
      RETURN
C CHECK THAT ABS(X), ABS(Y) .LT. B**T
   30 IF ((R(I2+1).LE.T) .AND. (R(I3+1).LE.T)) GO TO 50
   40 WRITE (LUN,FMT=9000)
      CALL MPERR
      Z(1) = 0
      RETURN
C START OF MAIN EUCLIDEAN ALGORITHM LOOP
   50 IF (R(I2).EQ.0) GO TO 10
      IF (MPCOMP(R(I2),R(I3))) 60,10,70
C EXCHANGE POINTERS ONLY
   60 IS = I2
      I2 = I3
      I3 = IS
C CHECK FOR SMALL EXPONENT
   70 IF (R(I2+1).LE.2) GO TO 100
C REDUCE T (TRAILING DIGITS MUST BE ZERO)
      T = R(I2+1)
      CALL MPSTR(R(I3),R(I4))
C FORCE EXPONENTS TO BE EQUAL
      R(I4+1) = R(I2+1)
C GET FIRST TWO DIGITS
      IQ = B*R(I2+2) + R(I2+3)
      IF (MPCOMP(R(I2),R(I4)).GE.0) GO TO 80
C REDUCE EXPONENT BY ONE
      R(I4+1) = R(I4+1) - 1
C UNDERESTIMATE QUOTIENT
      IQ = IQ/ (R(I4+2)+1)
      GO TO 90
C LEHMERS METHOD WOULD SAVE SOME MP OPERATIONS BUT NOT VERY
C MANY UNLESS WE COULD USE DOUBLE-PRECISION SAFELY.
   80 IQ = MAX0(1,IQ/ (B*R(I4+2)+R(I4+3)+1))
   90 CALL MPMULI(R(I4),IQ,R(I4))
      CALL MPSUB(R(I2),R(I4),R(I2))
      GO TO 50
C HERE SAFE TO USE INTEGER ARITHMETIC
  100 CALL MPCMI(R(I2),I)
      CALL MPCMI(R(I3),J)
      T = TS
  110 I = MOD(I,J)
      IF (I.EQ.0) GO TO 120
      J = MOD(J,I)
      IF (J.NE.0) GO TO 110
      J = I
  120 CALL MPCIM(J,Z)
      RETURN

 9000 FORMAT (' *** X OR Y NON-INTEGER OR TOO LARGE',' IN CALL TO MPGC',
     +       'DA ***')
      END
      SUBROUTINE MPGCDB(X,Y)
C RETURNS (X, Y) AS (X/Z, Y/Z) WHERE Z IS THE GCD OF X AND Y.
C X AND Y ARE INTEGERS REPRESENTED AS MP NUMBERS,
C AND MUST SATISFY ABS(X) .LT. B**T, ABS(Y) .LT. B**T
C TIME O(T**2), SPACE = 5T+12.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I2,IS,IZ
C     ..
C     .. External Functions ..
      INTEGER MPCMPI,MPCOMP
      EXTERNAL MPCMPI,MPCOMP
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADDQ,MPCHK,MPCIM,MPCMI,MPCMIM,MPDIVI,MPGCDA,MPMUL,
     +         MPREC,MPSTR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX0
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(5,12)
      I2 = 4*T + 11
C FIND GCD OF X AND Y USING MPGCDA
      CALL MPGCDA(X,Y,R(I2))
C CHECK FOR X AND Y EQUAL (WHEN MAY COINCIDE)
      IF (MPCOMP(X,Y).NE.0) GO TO 10
      IS = X(1)
      CALL MPCIM(IS,X)
      CALL MPSTR(X,Y)
      RETURN
C CHECK IF GCD IS SMALL.
   10 IF (MPCMPI(R(I2),7*MAX0(B*B,4096)).GT.0) GO TO 20
      CALL MPCMI(R(I2),IZ)
      IF (IZ.EQ.1) RETURN
      CALL MPDIVI(X,IZ,X)
      CALL MPDIVI(Y,IZ,Y)
      RETURN
C HERE GCD IS LARGE
   20 CALL MPREC(R(I2),R(I2))
      CALL MPMUL(X,R(I2),X)
      CALL MPMUL(Y,R(I2),Y)
C ADD SIGN/2 AND TRUNCATE TO GET CORRECT INTEGER
      IS = X(1)
      CALL MPADDQ(X,IS,2,X)
      CALL MPCMIM(X,X)
      IS = Y(1)
      CALL MPADDQ(Y,IS,2,Y)
      CALL MPCMIM(Y,Y)
      RETURN

      END
      LOGICAL FUNCTION MPGE(X,Y)
C RETURNS LOGICAL VALUE OF (X .GE. Y) FOR MP X AND Y.
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. External Functions ..
      INTEGER MPCOMP
      EXTERNAL MPCOMP
C     ..
      MPGE = (MPCOMP(X,Y).GE.0)
      RETURN

      END
      LOGICAL FUNCTION MPGT(X,Y)
C RETURNS LOGICAL VALUE OF (X .GT. Y) FOR MP X AND Y.
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. External Functions ..
      INTEGER MPCOMP
      EXTERNAL MPCOMP
C     ..
      MPGT = (MPCOMP(X,Y).GT.0)
      RETURN

      END
      SUBROUTINE MPHANK(X,NU,Y,ERROR)
C TRIES TO COMPUTE THE BESSEL FUNCTION J(NU,X) USING HANKELS
C ASYMPTOTIC SERIES.  NU IS A NONNEGATIVE INTEGER .LE. MAX(B,64),
C ERROR IS AN INTEGER, X AND Y ARE MP NUMBERS.
C RETURNS ERROR = 0 IF SUCCESSFUL (RESULT IN Y),
C         ERROR = 1 IF UNSUCCESSFUL (Y UNCHANGED)
C ERROR COULD BE INDUCED BY O(B**(1-T)) PERTURBATIONS IN
C X AND Y.   TIME IS O(T**3).
C CALLED BY MPBESJ, SPACE = 11T+24
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER ERROR,NU
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER B2,I2,I3,I4,I5,I6,I7,IE,K
C     ..
C     .. External Functions ..
      INTEGER MPCMPR
      EXTERNAL MPCMPR
C     ..
C     .. External Subroutines ..
      EXTERNAL MPABS,MPADD,MPART1,MPCHK,MPCIM,MPCOS,MPDIV,MPDIVI,MPMUL,
     +         MPMULI,MPMULQ,MPPWR,MPROOT,MPSIN,MPSTR,MPSUB
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ALOG,FLOAT,MAX0,MOD
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(11,24)
      ERROR = 1
      B2 = MAX0(B,64)
C GIVE ERROR RETURN IF NU IS NEGATIVE OR TOO LARGE.
      IF ((NU.LT.0) .OR. (NU.GT.B2)) RETURN
      I2 = 5*T + 13
      I3 = I2 + T + 2
      I4 = I3 + T + 2
      I5 = I4 + T + 2
      I6 = I5 + T + 2
      I7 = I6 + T + 2
C WORK WITH ABS(X)
      CALL MPABS(X,R(I2))
C CHECK IF ABS(X) CLEARLY TOO SMALL FOR ASYMPTOTIC SERIES
      IF (MPCMPR(R(I2),0.5E0*FLOAT(T)*ALOG(FLOAT(B))).LE.0) RETURN
      CALL MPPWR(X,-2,R(I3))
      CALL MPDIVI(R(I3),-64,R(I3))
      CALL MPCIM(1,R(I4))
      R(I5) = 0
      CALL MPSTR(R(I4),R(I6))
      IE = 1
      K = 0
C LOOP TO SUM TWO ASYMPTOTIC SERIES
   10 K = K + 2
C ERROR RETURN IF TERMS INCREASING
      IF (R(I6+1).GT.IE) RETURN
      IE = R(I6+1)
      IF (K.GT.B2) GO TO 20
      CALL MPMULQ(R(I6), (2* (NU+K)-3)* (2* (NU-K)+3),K-1,R(I6))
      CALL MPADD(R(I5),R(I6),R(I5))
      CALL MPMULQ(R(I6), (2* (NU+K)-1)* (2* (NU-K)+1),K,R(I6))
      GO TO 30
C HERE NEED TO SPLIT UP CALLS TO MPMULQ
   20 CALL MPMULQ(R(I6),2* (NU+K)-3,K-1,R(I6))
      CALL MPMULI(R(I6),2* (NU-K)+3,R(I6))
      CALL MPADD(R(I5),R(I6),R(I5))
      CALL MPMULQ(R(I6),2* (NU+K)-1,K,R(I6))
      CALL MPMULI(R(I6),2* (NU-K)+1,R(I6))
   30 CALL MPMUL(R(I6),R(I3),R(I6))
      CALL MPADD(R(I4),R(I6),R(I4))
C LOOP IF TERMS NOT SUFFICIENTLY SMALL YET
      IF ((R(I6).NE.0) .AND. (R(I6+1).GT. (-T))) GO TO 10
C END OF ASYMPTOTIC SERIES, NOW COMPUTE RESULT
      CALL MPDIV(R(I5),R(I2),R(I5))
      CALL MPDIVI(R(I5),8,R(I5))
C COMPUTE PI/4 (SLIGHTLY MORE ACCURATE THAN CALLING
C MPPI AND DIVIDING BY FOUR)
      CALL MPART1(5,R(I6))
      CALL MPMULI(R(I6),4,R(I6))
      CALL MPART1(239,R(I3))
      CALL MPSUB(R(I6),R(I3),R(I3))
C AVOID TOO MUCH CANCELLATION IN SUBTRACTING MULTIPLE OF PI
      CALL MPMULI(R(I3),MOD(2*NU+1,8),R(I6))
      CALL MPSUB(R(I2),R(I6),R(I6))
C COULD SAVE SOME TIME BY NOT COMPUTING BOTH SIN AND COS
      CALL MPCOS(R(I6),R(I7))
      CALL MPMUL(R(I4),R(I7),R(I4))
      CALL MPSIN(R(I6),R(I7))
      CALL MPMUL(R(I5),R(I7),R(I5))
      CALL MPSUB(R(I4),R(I5),R(I4))
      CALL MPMUL(R(I3),R(I2),R(I3))
      CALL MPMULI(R(I3),2,R(I3))
      CALL MPROOT(R(I3),-2,R(I3))
      CALL MPMUL(R(I3),R(I4),R(I3))
C CORRECT SIGN OF RESULT
      IF (MOD(NU,2).NE.0) R(I3) = R(I3)*X(1)
      ERROR = 0
      CALL MPSTR(R(I3),Y)
      RETURN

      END
      SUBROUTINE MPIN(C,X,N,ERROR)
C CONVERTS THE FIXED-POINT DECIMAL NUMBER (READ UNDER NA1
C FORMAT) IN C(1) ... C(N) TO A MULTIPLE-PRECISION NUMBER
C IN X.   IF C REPRESENTS A VALID NUMBER, ERROR IS RETURNED
C AS 0.  IF C DOES NOT REPRESENT A VALID NUMBER, ERROR
C IS RETURNED AS 1 AND X AS ZERO.
C LEADING AND TRAILING BLANKS ARE ALLOWED, EMBEDDED BLANKS
C (EXCEPT BETWEEN THE NUMBER AND ITS SIGN) ARE FORBIDDEN.
C IF THERE IS NO DECIMAL POINT ONE IS ASSUMED TO LIE JUST TO
C THE RIGHT OF THE LAST DECIMAL DIGIT.
C FOR EFFICIENCY CHOOSE B A POWER OF 10.
C X IS AN MP NUMBER, C AN INTEGER ARRAY, N AND ERROR INTEGERS.
C DIMENSION OF R IN CALLING PROGRAM .GE. 3T+11.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER ERROR,N
C     ..
C     .. Array Arguments ..
      INTEGER C(*),X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER FIRST,I,I2,I3,IB,IP,J,K,S,TEN,TP
C     ..
C     .. Local Arrays ..
      INTEGER NUM(14)
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADDI,MPCHK,MPCLR,MPDIVI,MPERR,MPMULI
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX0
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
C     .. Data statements ..
C TO READ INPUT IN OCTAL CHANGE 10 TO 8 IN NEXT
C DATA STATEMENT.  SIMILARLY FOR OTHER BASES LESS THAN 10.
      DATA NUM(1),NUM(2),NUM(3)/1H0,1H1,1H2/
      DATA NUM(4),NUM(5),NUM(6)/1H3,1H4,1H5/
      DATA NUM(7),NUM(8),NUM(9)/1H6,1H7,1H8/
      DATA NUM(10),NUM(11),NUM(12)/1H9,1H ,1H./
      DATA NUM(13),NUM(14)/1H+,1H-/
      DATA TEN/10/
C     ..
C CHECK LEGALITY OF B, T, M, MXR AND LUN
      CALL MPCHK(3,11)
C CHECK FOR N .GT. 0
      IF (N.GT.0) GO TO 10
      WRITE (LUN,FMT=9000)
      CALL MPERR
      X(1) = 0
      RETURN

   10 I2 = 2*T + 9
C USE ONE GUARD DIGIT
      CALL MPCLR(R(I2),T+1)
      T = T + 1
      FIRST = 1
      R(I2) = 0
      S = 1
      ERROR = 0
      IP = 0
      K = 0
C SCAN C FROM LEFT, SKIPPING BLANKS
   20 K = K + 1
      IF (C(K).NE.NUM(11)) GO TO 50
   30 IF (K.LT.N) GO TO 20
C FIELD WAS ALL BLANK - TREAT AS ERROR CONDITION
   40 X(1) = 0
      ERROR = 1
      T = T - 1
      RETURN
C NONBLANK CHARACTER FOUND
   50 DO 60 I = 1,14
          J = I
          IF (C(K).EQ.NUM(I)) GO TO 70
   60 CONTINUE
C ILLEGAL CHARACTER, SO ERROR
      GO TO 40
C LEGAL CHARACTER, SEE IF DIGIT OR POINT
   70 IF (J.GT.10) GO TO 90
C MUST BE DIGIT, SO CONTINUE FORMING NUMBER
      CALL MPMULI(R(I2),TEN,R(I2))
      CALL MPADDI(R(I2),J-1,R(I2))
      FIRST = 0
      K = K + 1
      IF (K.LE.N) GO TO 50
C RESTORE T, ROUND RESULT AND RETURN
   80 I3 = I2 + T
      T = T - 1
      IF ((2*R(I3+1)).GT.B) R(I3) = R(I3) + 1
C MULTIPLICATION BY +-1 ALSO FIXES UP LAST DIGIT IF NECESSARY
      CALL MPMULI(R(I2),S,X)
      RETURN
C NONDIGIT FOUND, IS IT SIGN, BLANK, OR POINT
   90 IF (J.EQ.12) GO TO 100
      IF (J.EQ.11) GO TO 160
C MUST BE SIGN, ONLY LEGAL IF FIRST = 1
      IF (FIRST.EQ.0) GO TO 40
      IF (J.EQ.14) S = -1
      FIRST = 0
      GO TO 30
C POINT ENCOUNTERED
  100 IP = 0
  110 K = K + 1
      IF (K.GT.N) GO TO 140
C LOOK AT C(K)
      DO 120 I = 1,11
          J = I
          IF (C(K).EQ.NUM(I)) GO TO 130
  120 CONTINUE
C ILLEGAL CHARACTER
      GO TO 40
C IF BLANK GO TO 170
  130 IF (J.EQ.11) GO TO 160
C DIGIT (AFTER POINT)
      IP = IP + 1
      CALL MPMULI(R(I2),TEN,R(I2))
      CALL MPADDI(R(I2),J-1,R(I2))
      GO TO 110
C END OF INPUT FIELD, MULTIPLY BY TEN**(-IP)
  140 IF (IP.LE.0) GO TO 80
      IB = MAX0(7*B*B,32767)/TEN
      TP = 1
      DO 150 I = 1,IP
          TP = TEN*TP
          IF ((TP.LE.IB) .AND. (TP.NE.B) .AND. (I.LT.IP)) GO TO 150
          CALL MPDIVI(R(I2),TP,R(I2))
          TP = 1
  150 CONTINUE
      GO TO 80
C TRAILING BLANK, CHECK THAT ALL TO RIGHT ARE BLANKS
  160 DO 170 I = K,N
          IF (C(I).NE.NUM(11)) GO TO 40
  170 CONTINUE
      GO TO 140

 9000 FORMAT (' *** N NOT POSITIVE IN CALL TO MPIN ***')
      END
      SUBROUTINE MPINE(C,X,N,J,ERROR)
C SAME AS MPIN EXCEPT THAT THE RESULT (X) IS MULTIPLIED BY
C 10**J, WHERE J IS A SINGLE-PRECISION INTEGER.  FOR DETAILS
C OF THE OTHER ARGUMENTS, SEE MPIN.
C USEFUL FOR FLOATING-POINT INPUT OF MP NUMBERS.  THE USER CAN
C READ THE EXPONENT INTO J (USING ANY SUITABLE FORMAT) AND
C THE FRACTION INTO C (USING A1 FORMAT), THEN CALL MPINE TO
C CONVERT TO MULTIPLE-PRECISION.
C SPACE = 5T+12
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER ERROR,J,N
C     ..
C     .. Array Arguments ..
      INTEGER C(*),X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I,I2,IB,JA,JP,TEN
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPCIM,MPDIV,MPDIVI,MPIN,MPMUL,MPMULI,MPPWR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IABS,MAX0
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
C     .. Data statements ..
C CHANGE NEXT DATA STATEMENT IF INPUT RADIX NOT 10
      DATA TEN/10/
C     ..
C CHECK LEGALITY OF B, T, M, LUN AND MXR
      CALL MPCHK(5,12)
C CALL MPIN TO CONVERT C TO MP FORMAT
      CALL MPIN(C,X,N,ERROR)
C RETURN IF J ZERO OR X ZERO
      IF ((J.EQ.0) .OR. (X(1).EQ.0)) RETURN
C OTHERWISE MULTIPLY BY TEN**J
      JA = IABS(J)
C THE NUMBERS -500 AND 100 WERE DETERMINED EMPIRICALLY.  THE OPTIMUM
C CHOICE DEPENDS ON B AND T.
      IF ((J.GT. (-500)) .AND. (J.LT.100)) GO TO 10
C HERE EXPONENT LARGE, SO USE MPPWR TO COMPUTE TEN**ABS(J)
C LEAVE SPACE FOR MPDIV
      I2 = 4*T + 11
      CALL MPCIM(TEN,R(I2))
      CALL MPPWR(R(I2),JA,R(I2))
      IF (J.LT.0) CALL MPDIV(X,R(I2),X)
      IF (J.GE.0) CALL MPMUL(X,R(I2),X)
      RETURN
C HERE ABS(J) IS SMALL SO PROBABLY FASTER TO USE MPDIVI OR MPMULI
   10 JP = 1
      IB = MAX0(7*B*B,32767)/TEN
      DO 20 I = 1,JA
          JP = TEN*JP
          IF ((JP.LE.IB) .AND. (JP.NE.B) .AND. (I.LT.JA)) GO TO 20
          IF (J.LT.0) CALL MPDIVI(X,JP,X)
          IF (J.GE.0) CALL MPMULI(X,JP,X)
          JP = 1
   20 CONTINUE
      RETURN

      END
      SUBROUTINE MPINF(X,N,UNIT,IFORM,ERR)
C READS N WORDS FROM LOGICAL UNIT IABS(UNIT) USING FORMAT IN IFORM,
C THEN CONVERTS TO MP NUMBER X USING ROUTINE MPIN.
C IFORM SHOULD CONTAIN A FORMAT WHICH ALLOWS FOR READING N WORDS
C IN A1 FORMAT, E.G. 6H(80A1)
C ERR RETURNED AS TRUE IF MPIN COULD NOT INTERPRET INPUT AS
C AN MP NUMBER OR IF N NOT POSITIVE, OTHERWISE FALSE.
C IF ERR IS TRUE THEN X IS RETURNED AS ZERO.
C SPACE REQUIRED 3T+N+11.
C CHECK THAT ENOUGH SPACE AVAILABLE
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER N,UNIT
      LOGICAL ERR
C     ..
C     .. Array Arguments ..
      INTEGER IFORM(*),X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I2,IER
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPIN,MPIO
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IABS
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(3,N+11)
      I2 = 3*T + 12
C READ N WORDS UNDER FORMAT IFORM.
      CALL MPIO(R(I2),N, (-IABS(UNIT)),IFORM,ERR)
      X(1) = 0
C RETURN IF ERROR
      IF (ERR) RETURN
C ELSE CONVERT TO MP NUMBER.
      CALL MPIN(R(I2),X,N,IER)
C RETURN ERROR FLAG IF MPIN OBJECTED
      ERR = (IER.NE.0)
      RETURN

      END
      SUBROUTINE MPINIT(X)
C DECLARES BLANK COMMON (USED BY MP PACKAGE) AND
C CALLS MPSET TO INITIALIZE PARAMETERS
C THE AUGMENT DECLARATION
C       INITIALIZE MP
C CAUSES A CALL TO MPINIT TO BE GENERATED.
C *** ASSUMES OUTPUT UNIT 6, 43 DECIMAL PLACES,
C *** 10 MP DIGITS, SPACE 296 WORDS.  IF THE AUGMENT
C *** DESCRIPTION DECK IS CHANGED THIS ROUTINE SHOULD
C *** BE CHANGED ACCORDINGLY.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. External Subroutines ..
      EXTERNAL MPSET
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPSET(6,43,12,296)
C ARE A SPECIAL CASE OF
C     INTEGER R(MXR)
C     CALL MPSET (LUN, IDECPL, T+2, MXR)
C WHERE LUN IS THE LOGICAL UNIT FOR OUTPUT,
C IDECPL IS THE EQUIVALENT NUMBER OF DECIMAL PLACES REQUIRED,
C T IS THE NUMBER OF MP DIGITS, AND
C MXR IS THE SIZE OF THE WORKING AREA USED BY MP
C (MXR = MAX (T*T+15*T+27, 14*T+156) IS SUFFICIENT).
C TO CHANGE THE PRECISION, MODIFY THE DIMENSIONS IN THE
C DECLARE STATEMENTS IN THE AUGMENT DESCRIPTION DECK -
C THE DIMENSION FOR TYPE MULTIPLE SHOULD BE T+2 AND
C FOR TYPE MULTIPAK SHOULD BE INT ((T+3)/2).
C SEE COMMENTS IN ROUTINE MPSET FOR THE NUMBER OF MP
C DIGITS REQUIRED TO GIVE THE EQUIVALENT OF ANY DESIRED
C NUMBER OF DECIMAL PLACES.
C *** ON SOME SYSTEMS A DECLARATION OF BLANK COMMON IN THE MAIN
C *** PROGRAM MAY BE NECESSARY.  IF SO, DECLARE
C ***       COMMON MPWORK(301)
C *** OR, MORE GENERALLY,
C ***       COMMON MPWORK(MXR+5)
C *** IN THE MAIN PROGRAM.
      RETURN

      END
      SUBROUTINE MPIO(C,N,UNIT,IFORM,ERR)
C IF UNIT .GT. 0 WRITES C(1), ... , C(N) IN FORMAT IFORM
C IF UNIT .LE. 0 READS  C(1), ... , C(N) IN FORMAT IFORM
C IN BOTH CASES USES LOGICAL UNIT IABS(UNIT).
C ERR IS RETURNED AS TRUE IF N NON-POSITIVE, OTHERWISE FALSE.
C WE WOULD LIKE TO RETURN ERR AS TRUE IF READ/WRITE ERROR DETECTED,
C BUT THIS CAN NOT BE DONE WITH ANSI STANDARD FORTRAN (1966).
C *** UNIVAC ASCII FORTRAN (FTN 5R1AE) DOES NOT WORK IF IFORM
C *** IS DECLARED WITH DIMENSION 1.  MOST FORTRANS DO THOUGH.
C     .. Scalar Arguments ..
      INTEGER N,UNIT
      LOGICAL ERR
C     ..
C     .. Array Arguments ..
      INTEGER C(N),IFORM(20)
C     ..
C     .. Local Scalars ..
      INTEGER IU
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IABS
C     ..
      ERR = (N.LE.0)
      IF (ERR) RETURN
      IU = IABS(UNIT)
      IF (UNIT.GT.0) WRITE (IU,FMT=IFORM) C
      IF (UNIT.LE.0) READ (IU,FMT=IFORM) C
      RETURN

      END
      SUBROUTINE MPKSTR(X,Y)
C SETS Y = X FOR PACKED MP NUMBERS X AND Y.
C ASSUMES SAME PACKED FORMAT AS MPPACK AND MPUNPK.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I,N
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      Y(2) = X(2)
C CHECK FOR ZERO
      IF (Y(2).EQ.0) RETURN
C HERE X NONZERO SO MOVE PACKED NUMBER
      N = (T+3)/2
      DO 10 I = 1,N
          Y(I) = X(I)
   10 CONTINUE
      RETURN

      END
      LOGICAL FUNCTION MPLE(X,Y)
C RETURNS LOGICAL VALUE OF (X .LE. Y) FOR MP X AND Y.
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. External Functions ..
      INTEGER MPCOMP
      EXTERNAL MPCOMP
C     ..
      MPLE = (MPCOMP(X,Y).LE.0)
      RETURN

      END
      SUBROUTINE MPLI(X,Y)
C RETURNS Y = LI(X) = LOGARITHMIC INTEGRAL OF X
C           = (PRINCIPAL VALUE INTEGRAL FROM 0 TO X OF
C              DU/LOG(U)),
C USING MPEI.  X AND Y ARE MP NUMBERS, X .GE. 0, X .NE. 1.
C ERROR IN Y COULD BE INDUCED BY AN O(B**(1-T)) RELATIVE
C PERTURBATION IN X FOLLOWED BY SIMILAR PERTURBATION IN Y.
C THUS RELATIVE ERROR IN Y IS SMALL UNLESS X IS CLOSE TO
C 1 OR TO THE ZERO 1.45136923488338105028... OF LI(X).
C TIME IS O(T.M(T)), SPACE = 10T+38
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. External Functions ..
      INTEGER MPCMPI
      EXTERNAL MPCMPI
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPEI,MPERR,MPLN,MPOVFL
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(10,38)
      IF (X(1)) 10,20,30
C HERE X NEGATIVE, GIVE ERROR MESSAGE
   10 WRITE (LUN,FMT=9000)
      CALL MPERR
C LI(0) = 0
   20 Y(1) = 0
      RETURN
C HERE X IS POSITIVE, SEE IF EQUAL TO 1
   30 IF (MPCMPI(X,1).NE.0) GO TO 40
C HERE X EXACTLY EQUAL TO 1, GIVE ERROR MESSAGE AND
C TREAT AS MP OVERFLOW
      WRITE (LUN,FMT=9010)
      CALL MPOVFL(Y)
      RETURN
C HERE X POSITIVE AND .NE. 1, SO USE EI(LN(X))
   40 CALL MPLN(X,Y)
      CALL MPEI(Y,Y)
      RETURN

 9000 FORMAT (' *** X NEGATIVE IN CALL TO MPLI ***')
 9010 FORMAT (' *** X .EQ. 1 IN CALL TO MPLI ***')
      END
      SUBROUTINE MPLN(X,Y)
C RETURNS Y = LN(X), FOR MP X AND Y, USING MPLNS.
C RESTRICTION - INTEGER PART OF LN(X) MUST BE REPRESENTABLE
C AS A SINGLE-PRECISION INTEGER.  TIME IS O(SQRT(T).M(T)).
C FOR SMALL INTEGER X, MPLNI IS FASTER.
C ASYMPTOTICALLY FASTER METHODS EXIST (EG THE GAUSS-SALAMIN
C METHOD, SEE MPLNGS), BUT ARE NOT USEFUL UNLESS T IS LARGE.
C SEE COMMENTS TO MPATAN, MPEXP1 AND MPPIGL.
C DIMENSION OF R IN CALLING PROGRAM MUST BE AT LEAST 6T+14.
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      REAL RLX,RX
      INTEGER E,I2,I3,K
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADD,MPADDI,MPCHK,MPCMR,MPCRM,MPERR,MPEXP,MPLNS,MPMUL,
     +         MPSTR,MPSUB
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ALOG,FLOAT
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(6,14)
      I2 = 4*T + 11
      I3 = I2 + T + 2
C CHECK THAT X IS POSITIVE
      IF (X(1).GT.0) GO TO 10
      WRITE (LUN,FMT=9000)
      CALL MPERR
      Y(1) = 0
      RETURN
C MOVE X TO LOCAL STORAGE
   10 CALL MPSTR(X,R(I2))
      Y(1) = 0
      K = 0
C LOOP TO GET APPROXIMATE LN(X) USING SINGLE-PRECISION
   20 CALL MPADDI(R(I2),-1,R(I3))
C IF POSSIBLE GO TO CALL MPLNS
      IF ((R(I3).EQ.0) .OR. ((R(I3+1)+1).LE.0)) GO TO 30
C REMOVE EXPONENT TO AVOID FLOATING-POINT OVERFLOW
      E = R(I2+1)
      R(I2+1) = 0
      CALL MPCMR(R(I2),RX)
C RESTORE EXPONENT AND COMPUTE SINGLE-PRECISION LOG
      R(I2+1) = E
      RLX = ALOG(RX) + FLOAT(E)*ALOG(FLOAT(B))
      CALL MPCRM(-RLX,R(I3))
C UPDATE Y AND COMPUTE ACCURATE EXP OF APPROXIMATE LOG
      CALL MPSUB(Y,R(I3),Y)
      CALL MPEXP(R(I3),R(I3))
C COMPUTE RESIDUAL WHOSE LOG IS STILL TO BE FOUND
      CALL MPMUL(R(I2),R(I3),R(I2))
C MAKE SURE NOT LOOPING INDEFINITELY
      K = K + 1
      IF (K.LT.10) GO TO 20
      WRITE (LUN,FMT=9010)
      CALL MPERR
      RETURN
C COMPUTE FINAL CORRECTION ACCURATELY USING MPLNS
   30 CALL MPLNS(R(I3),R(I3))
      CALL MPADD(Y,R(I3),Y)
      RETURN

 9000 FORMAT (' *** X NONPOSITIVE IN CALL TO MPLN ***')
 9010 FORMAT (' *** ERROR IN MPLN, ITERATION NOT CONVERGING ***')
      END
      SUBROUTINE MPLNGM(X,Y)
C RETURNS MP Y = LN(GAMMA(X)) FOR POSITIVE MP X, USING STIRLINGS
C ASYMPTOTIC APPROXIMATION.  SLOWER THAN MPGAMQ (UNLESS X LARGE)
C AND USES MORE SPACE, SO USE MPGAMQ AND MPLN IF X IS RATIONAL AND
C NOT TOO LARGE, SAY X .LE. 100.  TIME IS O(T**3).
C SPACE REQUIRED IS 11T+24+NL*((T+3)/2), WHERE NL IS THE NUMBER
C OF TERMS USED IN THE ASYMPTOTIC EXPANSION,
C NL .LE. 2+AL*T*LN(B), WHERE AL IS GIVEN BELOW.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      REAL AB,AL,RX,XL,XLN
      INTEGER I,I2,I4,I5,I6,IP,IR2,NL,NLP,P
C     ..
C     .. External Functions ..
      INTEGER MPCMPR
      EXTERNAL MPCMPR
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADD,MPADDI,MPADDQ,MPBERN,MPCHK,MPCIM,MPCMR,MPDIVI,
     +         MPERR,MPLN,MPMUL,MPMULI,MPPI,MPPWR,MPSTR,MPSUB,MPUNPK
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ALOG,EXP,FLOAT,INT,MIN0
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
C     .. Data statements ..
C AL .LT. 1 IS CHOSEN EMPIRICALLY TO MINIMIZE TIME.
      DATA AL/0.125E0/
C     ..
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C ASSUMING THAT NL IS ZERO
      CALL MPCHK(11,24)
C MAKE PRELIMINARY ESTIMATE OF NL
      AB = ALOG(FLOAT(B))
      NLP = INT(AL*FLOAT(T)*AB) + 2
C ESTIMATE HOW LARGE X NEEDS TO BE FOR SUFFICIENT ACCURACY
      XL = FLOAT(NLP)*EXP(0.5E0/AL-1E0)/3.14159E0
C CHECK THAT X IS POSITIVE
      IF (X(1).GT.0) GO TO 10
      WRITE (LUN,FMT=9000)
      CALL MPERR
      Y(1) = 0
      RETURN
C ALLOW SPACE FOR MPBERN AND LEAVE R(8T+19) TO R(9T+20) FOR MPGAM
   10 I2 = 7*T + 17
      I4 = 9*T + 21
      I5 = I4 + T + 2
      I6 = I5 + T + 2
C MOVE X AND SET Y = 0
      CALL MPSTR(X,R(I4))
      Y(1) = 0
C SEE IF X LARGE ENOUGH TO USE ASYMPTOTIC SERIES
      IF (MPCMPR(R(I4),XL).GE.0) GO TO 30
C HERE X NOT LARGE ENOUGH, SO INCREASE USING THE
C IDENTITY GAMMA(X+1) = X*GAMMA(X) TO CORRECT RESULT
      CALL MPCIM(1,Y)
   20 CALL MPMUL(Y,R(I4),Y)
      CALL MPADDI(R(I4),1,R(I4))
      IF (MPCMPR(R(I4),XL).LT.0) GO TO 20
      CALL MPLN(Y,Y)
      Y(1) = -Y(1)
C COMPUTE FIRST TERMS IN STIRLINGS APPROXIMATION
   30 CALL MPLN(R(I4),R(I5))
      CALL MPADDQ(R(I4),-1,2,R(I2))
      CALL MPMUL(R(I2),R(I5),R(I5))
      CALL MPSUB(R(I5),R(I4),R(I5))
      CALL MPADD(Y,R(I5),Y)
      CALL MPPI(R(I5))
      CALL MPMULI(R(I5),2,R(I5))
      CALL MPLN(R(I5),R(I5))
      CALL MPDIVI(R(I5),2,R(I5))
      CALL MPADD(Y,R(I5),Y)
C IF X VERY LARGE CAN RETURN HERE
      IF (R(I4+1).GE.T) RETURN
C DEPENDING ON HOW LARGE X IS, MAY BE ABLE TO DECREASE NL HERE
      IR2 = R(I4+1)
C CONVERT TO REAL AFTER ENSURING NO OVERFLOW
      R(I4+1) = 0
      CALL MPCMR(R(I4),RX)
      R(I4+1) = IR2
      XLN = 1.0E0 + ALOG(3.14159E0*RX/FLOAT(NLP)) + AB*FLOAT(IR2)
      NL = MIN0(NLP,INT(0.5E0*FLOAT(T)*AB/XLN))
      IF (NL.LE.0) RETURN
      CALL MPPWR(R(I4),-2,R(I5))
      P = (T+3)/2
C CHECK THAT MXR LARGE ENOUGH
      CALL MPCHK(11,NL*P+24)
C COMPUTE BERNOULLI NUMBERS REQUIRED (MUCH TIME COULD BE
C SAVED IF THESE WERE PRECOMPUTED)
      CALL MPBERN(NL,P,R(I6))
C SUM ASYMPTOTIC SERIES
      DO 40 I = 1,NL
          IP = I6 + (I-1)*P
          CALL MPUNPK(R(IP),R(I2))
          CALL MPDIVI(R(I2),2*I,R(I2))
          CALL MPDIVI(R(I2),2*I-1,R(I2))
          CALL MPMUL(R(I4),R(I5),R(I4))
          CALL MPMUL(R(I4),R(I2),R(I2))
          IF ((R(I2).EQ.0) .OR. (R(I2+1).LE. (-T))) RETURN
          CALL MPADD(Y,R(I2),Y)
   40 CONTINUE
      RETURN

 9000 FORMAT (' *** X NONPOSITIVE IN CALL TO MPLNGM ***')
      END
      SUBROUTINE MPLNGS(X,Y)
C RETURNS Y = LN(X) FOR MP X AND Y, USING THE GAUSS-SALAMIN
C ALGORITHM BASED ON THE ARITHMETIC-GEOMETRIC MEAN ITERATION
C (SEE ANALYTIC COMPUTATIONAL COMPLEXITY (ED. BY J. F. TRAUB),
C ACADEMIC PRESS, 1976, 151-176) UNLESS X IS CLOSE TO 1.
C SPACE = 6T+26, TIME = O(LOG(T)M(T)) + O(T**2) IF
C ABS(X-1) .GE. 1/B AND AS FOR MPLNS OTHERWISE.
C SLOWER THAN MPLN UNLESS T IS LARGE (.GE. ABOUT 500) SO
C MAINLY USEFUL FOR TESTING PURPOSES.
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      REAL B2L
      INTEGER E,I,I2,I3,I4,IT,N,T2
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADD,MPADDI,MPCHK,MPCIM,MPCLR,MPDIV,MPDIVI,MPERR,MPLNI,
     +         MPLNS,MPMUL,MPMULI,MPPI,MPSQRT,MPSTR,MPSUB
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ALOG,FLOAT,IABS,INT
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(6,26)
      IF (X(1).GT.0) GO TO 10
      WRITE (LUN,FMT=9000)
      CALL MPERR
      Y(1) = 0
      RETURN
C ALLOW SPACE FOR MPSQRT WITH T+2 DIGITS (OVERLAP R(I2))
   10 I2 = 3*T + 15
      I3 = I2 + T + 4
      I4 = I3 + T + 4
C SEE IF X CLOSE TO 1
      CALL MPADDI(X,-1,R(I4))
      IF ((R(I4).NE.0) .AND. (R(I4+1).GE.0)) GO TO 20
C HERE ABS(X-1) .LT. 1/B SO GAUSS-SALAMIN ALGORITHM COULD BE
C INACCURATE BECAUSE OF CANCELLATION.  THE PRECISION COULD BE
C INCREASED TO COMPENSATE FOR THIS, BUT SIMPLER TO USE MPLNS.
      CALL MPLNS(R(I4),Y)
      RETURN
C PREPARE TO USE 2 GUARD DIGITS (BECAUSE SOME CANCELLATION)
   20 CALL MPCLR(R(I2),T+2)
      CALL MPSTR(X,R(I2))
      T = T + 2
      T2 = (T+1)/2
      E = X(2)
C MODIFY EXPONENT TO MAKE RI2 SUFFICIENTLY SMALL THAT
C ERROR WILL BE NEGLIGIBLE
      R(I2+1) = -T2
      CALL MPMULI(R(I2),4,R(I2))
      CALL MPCIM(1,R(I3))
C COMPUTE NUMBER OF ITERATIONS REQUIRED.  THE CONSTANT
C 2.36... IS ALMOST OPTIMAL.
      B2L = ALOG(FLOAT(B))/ALOG(2E0)
      N = INT(ALOG(FLOAT(T2+1)*B2L* (3E0+FLOAT(T)*B2L))/ALOG(2E0)-
     +    2.36E0)
C ARITHMETIC-GEOMETRIC MEAN LOOP
      DO 30 I = 1,N
          CALL MPADD(R(I2),R(I3),R(I4))
          CALL MPDIVI(R(I4),2,R(I4))
          CALL MPMUL(R(I2),R(I3),R(I3))
          CALL MPSQRT(R(I3),R(I2))
C FASTER TO EXCHANGE POINTERS THAN MP NUMBERS
          IT = I3
          I3 = I4
          I4 = IT
   30 CONTINUE
C CHECK THAT CONVERGENCE OCCURRED
      CALL MPSUB(R(I2),R(I3),R(I4))
      IF ((R(I4).EQ.0) .OR. ((R(I2+1)-R(I4+1)).GE. (T-3))) GO TO 40
      WRITE (LUN,FMT=9010)
      CALL MPERR
C COULD SAVE SOME TIME BY PRECOMPUTING PI AND LN(B)
   40 CALL MPPI(R(I4))
      CALL MPDIV(R(I4),R(I3),R(I3))
      CALL MPDIVI(R(I3),2,R(I3))
      CALL MPLNI(IABS(B),R(I4))
      CALL MPMULI(R(I4),E+T2,R(I4))
C ALLOW FOR MODIFIED EXPONENT
      CALL MPSUB(R(I4),R(I3),R(I3))
C RESTORE T AND RETURN
      T = T - 2
      CALL MPSTR(R(I3),Y)
      RETURN

 9000 FORMAT (' *** X NONPOSITIVE IN CALL TO MPLNGS ***')
 9010 FORMAT (' *** ITERATION FAILED TO CONVERGE IN MPLNGS ***')
      END
      SUBROUTINE MPLNI(N,X)
C RETURNS MULTIPLE-PRECISION X = LN(N) FOR SMALL POSITIVE
C INTEGER N, TIME IS O(T**2).
C METHOD IS TO USE A RAPIDLY CONVERGING SERIES AND MPL235.
C DIMENSION OF R IN CALLING PROGRAM AT LEAST 3T+8.
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      REAL RX
      INTEGER B2,I,I2,I3,IA,IK,IM,IP,IP2,IQ,IQ2,J,N1,N2,TS
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADD,MPADD2,MPCHK,MPCMR,MPCQM,MPDIVI,MPERR,MPGCD,MPL235,
     +         MPMULQ,MPSTR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ALOG,FLOAT,IABS,MAX0,MIN0
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(3,8)
C CHECK FOR N = 1 AND N .LT. 1
      IF (N-1) 10,20,30
   10 WRITE (LUN,FMT=9000)
      GO TO 90
C LN(1) = 0
   20 X(1) = 0
      RETURN
C HERE N .GE. 2
   30 IF (N.GT.2) GO TO 40
C N = 2 IS A SPECIAL CASE
      CALL MPL235(1,0,0,X)
      GO TO 150
C HERE N .GE. 3
   40 B2 = MAX0(B,64)
      IF (N.GT. (3*B2*B2)) GO TO 80
      J = 3
      IA = 0
      N2 = N/2
   50 IF (J.GT.N2) GO TO 60
      IA = IA + 1
      J = 2*J
      GO TO 50
C NOW J = 3*(2**IA) .LE. N .LT. 6*(2**IA)
   60 J = J/3
      IM = N
      IK = 0
      DO 70 I = 3,6
          N1 = I*J
          IF (IABS(N1-N).GT.IM) GO TO 70
          IM = IABS(N1-N)
          IK = I
   70 CONTINUE
      N1 = IK*J
C NOW N IS CLOSE TO N1 = IK*(2**IA)
C AND IK = 3, 4, 5 OR 6, SO MPL235 GIVES LN(N1).
      IF (IK.EQ.3) CALL MPL235(IA,1,0,X)
      IF (IK.EQ.4) CALL MPL235(IA+2,0,0,X)
      IF (IK.EQ.5) CALL MPL235(IA,0,1,X)
      IF (IK.EQ.6) CALL MPL235(IA+1,1,0,X)
      IF (N.EQ.N1) GO TO 150
C NOW NEED LN(N/N1).
      N2 = N
      CALL MPGCD(N2,N1)
      IP = N2 - N1
      IQ = N2 + N1
C CHECK FOR POSSIBLE INTEGER OVERFLOW
      IF (IQ.GT.14) GO TO 100
   80 WRITE (LUN,FMT=9010)
   90 CALL MPERR
      X(1) = 0
      RETURN
C REDUCE TO LOWEST TERMS
  100 CALL MPGCD(IP,IQ)
      TS = T
      I2 = T + 5
      I3 = I2 + T + 2
      CALL MPCQM(2*IP,IQ,R(I2))
      CALL MPSTR(R(I2),R(I3))
      CALL MPADD(X,R(I3),X)
      I = 1
      IF (IQ.GT.B2) GO TO 110
      IQ2 = IQ**2
      IP2 = IP**2
C LOOP TO SUM SERIES FOR LN(N2/N1)
  110 I = I + 2
      IF (R(I2).EQ.0) GO TO 140
C REDUCE T IF POSSIBLE, DONE IF CAN REDUCE BELOW 2
      T = TS + R(I2+1) + 2
      IF (T.LE.2) GO TO 140
      T = MIN0(T,TS)
C SPLIT UP CALL TO MPMULQ IF IQ TOO LARGE
      IF (IQ.GT.B2) GO TO 120
      CALL MPMULQ(R(I2),IP2,IQ2,R(I2))
      GO TO 130
C HERE IQ TOO LARGE FOR ONE CALL TO MPMULQ
  120 CALL MPMULQ(R(I2),IP,IQ,R(I2))
      CALL MPMULQ(R(I2),IP,IQ,R(I2))
  130 CALL MPDIVI(R(I2),I,R(I3))
C RESTORE T AND ACCUMULATE SUM
      T = TS
      CALL MPADD2(R(I3),X,X,X,0)
      GO TO 110

  140 T = TS
C RETURN IF RESULT ACCURATE TO RELATIVE ERROR 0.01
  150 CALL MPCMR(X,RX)
      IF (ABS(RX-ALOG(FLOAT(N))).LT. (0.01*RX)) RETURN
      WRITE (LUN,FMT=9020)
C THE FOLLOWING MESSAGE MAY INDICATE
C THAT B**(T-1) IS TOO SMALL OR THAT N IS TOO LARGE.
      CALL MPERR
      RETURN

 9000 FORMAT (' *** N NOT POSITIVE IN CALL TO MPLNI ***')
 9010 FORMAT (' *** N TOO LARGE IN CALL TO MPLNI ***')
 9020 FORMAT (' *** ERROR OCCURRED IN MPLNI, RESULT INCORRECT ***')
      END
      SUBROUTINE MPLNS(X,Y)
C RETURNS MP Y = LN(1+X) IF X IS AN MP NUMBER SATISFYING THE
C CONDITION ABS(X) .LT. 1/B, ERROR OTHERWISE.
C USES NEWTONS METHOD TO SOLVE THE EQUATION
C EXP1(-Y) = X, THEN REVERSES SIGN OF Y.
C (HERE EXP1(Y) = EXP(Y) - 1 IS COMPUTED USING MPEXP1).
C TIME IS O(SQRT(T).M(T)) AS FOR MPEXP1, SPACE = 5T+12.
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I2,I3,I4,IT0,TS,TS2,TS3
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADD,MPADDI,MPADDQ,MPCHK,MPDIVI,MPERR,MPEXP1,MPMUL,
     +         MPSTR,MPSUB
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX0
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(5,12)
      I2 = 2*T + 7
      I3 = I2 + T + 2
      I4 = I3 + T + 2
C CHECK FOR X = 0 EXACTLY
      IF (X(1).NE.0) GO TO 10
      Y(1) = 0
      RETURN
C CHECK THAT ABS(X) .LT. 1/B
   10 IF ((X(2)+1).LE.0) GO TO 20
      WRITE (LUN,FMT=9000)
      CALL MPERR
      Y(1) = 0
      RETURN
C SAVE T AND GET STARTING APPROXIMATION TO -LN(1+X)
   20 TS = T
      CALL MPSTR(X,R(I3))
      CALL MPDIVI(X,4,R(I2))
      CALL MPADDQ(R(I2),-1,3,R(I2))
      CALL MPMUL(X,R(I2),R(I2))
      CALL MPADDQ(R(I2),1,2,R(I2))
      CALL MPMUL(X,R(I2),R(I2))
      CALL MPADDI(R(I2),-1,R(I2))
      CALL MPMUL(X,R(I2),Y)
C START NEWTON ITERATION USING SMALL T, LATER INCREASE
      T = MAX0(5,13-2*B)
      IF (T.GT.TS) GO TO 60
      IT0 = (T+5)/2
   30 CALL MPEXP1(Y,R(I4))
      CALL MPMUL(R(I3),R(I4),R(I2))
      CALL MPADD(R(I4),R(I2),R(I4))
      CALL MPADD(R(I3),R(I4),R(I4))
      CALL MPSUB(Y,R(I4),Y)
      IF (T.GE.TS) GO TO 50
C FOLLOWING LOOP COMPUTES NEXT VALUE OF T TO USE.
C BECAUSE NEWTONS METHOD HAS 2ND ORDER CONVERGENCE,
C WE CAN ALMOST DOUBLE T EACH TIME.
      TS3 = T
      T = TS
   40 TS2 = T
      T = (T+IT0)/2
      IF (T.GT.TS3) GO TO 40
      T = TS2
      GO TO 30
C CHECK THAT NEWTON ITERATION WAS CONVERGING AS EXPECTED
   50 IF ((R(I4).EQ.0) .OR. ((2*R(I4+1)).LE. (IT0-T))) GO TO 60
      WRITE (LUN,FMT=9010)
      CALL MPERR
C REVERSE SIGN OF Y AND RETURN
   60 Y(1) = -Y(1)
      T = TS
      RETURN

 9000 FORMAT (' *** ABS(X) .GE. 1/B IN CALL TO MPLNS ***')
 9010 FORMAT (' *** ERROR OCCURRED IN MPLNS, NEWTON ITERATION NOT',' C',
     +       'ONVERGING PROPERLY ***')
      END
      LOGICAL FUNCTION MPLT(X,Y)
C RETURNS LOGICAL VALUE OF (X .LT. Y) FOR MP X AND Y.
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. External Functions ..
      INTEGER MPCOMP
      EXTERNAL MPCOMP
C     ..
      MPLT = (MPCOMP(X,Y).LT.0)
      RETURN

      END
      SUBROUTINE MPL235(I,J,K,X)
C RETURNS MP X = LN((2**I)*(3**J)*(5**K)), FOR INTEGER I, J AND K.
C THE METHOD REQUIRES TIME O(T**2).  LN(81/80), LN(25/24) AND
C LN(16/15) ARE CALCULATED FIRST.  MPL235 COULD BE SPEEDED
C UP IF THESE CONSTANTS WERE PRECOMPUTED AND SAVED.
C ASSUMED THAT I, J AND K NOT TOO LARGE.
C DIMENSION OF R IN CALLING PROGRAM MUST BE AT LEAST 3T+8
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER I,J,K
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER B2,D2,I2,I3,N,Q,TS
C     ..
C     .. Local Arrays ..
      INTEGER C(3),D(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADD,MPADD2,MPCHK,MPCQM,MPDIVI,MPERR,MPSTR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IABS,MAX0,MIN0
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
C     .. Data statements ..
      DATA D(1),D(2),D(3)/161,49,31/
C     ..
C CHECK LEGALITY OF B, T, M, LUN AND MXR
      CALL MPCHK(3,8)
      X(1) = 0
      IF (MAX0(IABS(I),IABS(J),IABS(K)).LT.
     +    (MAX0(B*B,4096)/9)) GO TO 10
      WRITE (LUN,FMT=9000)
      CALL MPERR
      RETURN

   10 B2 = 2*MAX0(B,64)
      I2 = T + 5
      I3 = I2 + T + 2
      C(1) = 3*I + 5*J + 7*K
      C(2) = 5*I + 8*J + 12*K
      C(3) = 7*I + 11*J + 16*K
      TS = T
      DO 60 Q = 1,3
          CALL MPCQM(2*C(Q),D(Q),R(I2))
          CALL MPSTR(R(I2),R(I3))
          CALL MPADD(X,R(I3),X)
          IF (D(Q).LE.B2) D2 = D(Q)**2
          N = 1
   20     N = N + 2
          IF (R(I2).EQ.0) GO TO 50
C REDUCE T IF POSSIBLE
          T = TS + R(I2+1) + 2 - X(2)
          IF (T.LE.2) GO TO 50
          T = MIN0(T,TS)
C IF D(Q)**2 NOT REPRESENTABLE AS AN INTEGER, THE FOLLOWING
C DIVISION MUST BE SPLIT UP
          IF (D(Q).GT.B2) GO TO 30
          CALL MPDIVI(R(I2),D2,R(I2))
          GO TO 40

   30     CALL MPDIVI(R(I2),D(Q),R(I2))
          CALL MPDIVI(R(I2),D(Q),R(I2))
   40     CALL MPDIVI(R(I2),N,R(I3))
          T = TS
          CALL MPADD2(R(I3),X,X,X,0)
          GO TO 20

   50     T = TS
   60 CONTINUE
      RETURN

 9000 FORMAT (' *** I, J OR K TOO LARGE IN CALL TO MPLNI ***')
      END
      SUBROUTINE MPMAX(X,Y,Z)
C SETS Z = MAX (X, Y) WHERE X, Y AND Z ARE MULTIPLE-PRECISION
C     .. Array Arguments ..
      INTEGER X(*),Y(*),Z(*)
C     ..
C     .. External Functions ..
      INTEGER MPCOMP
      EXTERNAL MPCOMP
C     ..
C     .. External Subroutines ..
      EXTERNAL MPSTR
C     ..
      IF (MPCOMP(X,Y).GE.0) GO TO 10
C HERE X .LT. Y
      CALL MPSTR(Y,Z)
      RETURN
C HERE X .GE. Y
   10 CALL MPSTR(X,Z)
      RETURN

      END
      SUBROUTINE MPMAXR(X)
C SETS X TO THE LARGEST POSSIBLE POSITIVE MP NUMBER
C CHECK LEGALITY OF B, T, M, MXR AND LUN
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I,IT
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(1,4)
      IT = B - 1
C SET FRACTION DIGITS TO B-1
      DO 10 I = 1,T
          X(I+2) = IT
   10 CONTINUE
C SET SIGN AND EXPONENT
      X(1) = 1
      X(2) = M
      RETURN

      END
      INTEGER FUNCTION MPMEXA(X)
C RETURNS THE MAXIMUM ALLOWABLE EXPONENT OF MP NUMBERS (THE THIRD
C WORD OF COMMON).  X IS A DUMMY MP ARGUMENT.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      MPMEXA = M
      RETURN

      END
      SUBROUTINE MPMEXB(I,X)
C SETS THE MAXIMUM ALLOWABLE EXPONENT OF MP NUMBERS (I.E. THE
C THIRD WORD OF COMMON) TO I.
C I SHOULD BE GREATER THAN T, AND 4*I SHOULD BE REPRESENTABLE
C AS A SINGLE-PRECISION INTEGER.
C X IS A DUMMY MP ARGUMENT (AUGMENT EXPECTS ONE).
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER I
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. External Subroutines ..
      EXTERNAL MPERR
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      M = I
C CHECK LEGALITY OF M.  IF TOO LARGE, 4*M MAY OVERFLOW AND TEST .LE. 0
      IF ((M.GT.T) .AND. ((4*M).GT.0)) RETURN
      WRITE (LUN,FMT=9000)
      CALL MPERR
      RETURN

 9000 FORMAT (' *** ATTEMPT TO SET ILLEGAL MAXIMUM EXPONENT',' IN CALL',
     +       ' TO MPMEXB ***')
      END
      SUBROUTINE MPMIN(X,Y,Z)
C SETS Z = MIN (X, Y) WHERE X, Y AND Z ARE MULTIPLE-PRECISION
C     .. Array Arguments ..
      INTEGER X(*),Y(*),Z(*)
C     ..
C     .. External Functions ..
      INTEGER MPCOMP
      EXTERNAL MPCOMP
C     ..
C     .. External Subroutines ..
      EXTERNAL MPSTR
C     ..
      IF (MPCOMP(X,Y).GE.0) GO TO 10
C HERE X .LT. Y
      CALL MPSTR(X,Z)
      RETURN
C HERE X .GE. Y
   10 CALL MPSTR(Y,Z)
      RETURN

      END
      SUBROUTINE MPMINR(X)
C SETS X TO THE SMALLEST POSITIVE NORMALIZED MP NUMBER
C CHECK LEGALITY OF B, T, M, MXR AND LUN
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(1,4)
C SET FRACTION DIGITS TO ZERO
      DO 10 I = 2,T
          X(I+2) = 0
   10 CONTINUE
C SET SIGN, EXPONENT AND FIRST FRACTION DIGIT
      X(1) = 1
      X(2) = -M
      X(3) = 1
      RETURN

      END
      SUBROUTINE MPMLP(U,V,W,J)
C PERFORMS INNER MULTIPLICATION LOOP FOR MPMUL
C NOTE THAT CARRIES ARE NOT PROPAGATED IN INNER LOOP,
C WHICH SAVES TIME AT THE EXPENSE OF SPACE.
C     .. Scalar Arguments ..
      INTEGER J,W
C     ..
C     .. Array Arguments ..
      INTEGER U(*),V(*)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
      DO 10 I = 1,J
          U(I) = U(I) + W*V(I)
   10 CONTINUE
      RETURN

      END
      SUBROUTINE MPMUL(X,Y,Z)
C MULTIPLIES X AND Y, RETURNING RESULT IN Z, FOR MP X, Y AND Z.
C THE SIMPLE O(T**2) ALGORITHM IS USED, WITH
C FOUR GUARD DIGITS AND R*-ROUNDING.
C ADVANTAGE IS TAKEN OF ZERO DIGITS IN X, BUT NOT IN Y.
C ASYMPTOTICALLY FASTER ALGORITHMS ARE KNOWN (SEE KNUTH,
C VOL. 2), BUT ARE DIFFICULT TO IMPLEMENT IN FORTRAN IN AN
C EFFICIENT AND MACHINE-INDEPENDENT MANNER.
C IN COMMENTS TO OTHER MP ROUTINES, M(T) IS THE TIME
C TO PERFORM T-DIGIT MP MULTIPLICATION.   THUS
C M(T) = O(T**2) WITH THE PRESENT VERSION OF MPMUL,
C BUT M(T) = O(T.LOG(T).LOG(LOG(T))) IS THEORETICALLY POSSIBLE.
C CHECK LEGALITY OF B, T, M, MXR AND LUN
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*),Z(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER C,I,I2,I2P,J,J1,RE,RI,RS,XI
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPERR,MPMLP,MPNZR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN0
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(1,4)
      I2 = T + 4
      I2P = I2 + 1
C FORM SIGN OF PRODUCT
      RS = X(1)*Y(1)
      IF (RS.NE.0) GO TO 10
C SET RESULT TO ZERO
      Z(1) = 0
      RETURN
C FORM EXPONENT OF PRODUCT
   10 RE = X(2) + Y(2)
C CLEAR ACCUMULATOR
      DO 20 I = 1,I2
          R(I) = 0
   20 CONTINUE
C PERFORM MULTIPLICATION
      C = 8
      DO 40 I = 1,T
          XI = X(I+2)
C FOR SPEED, PUT THE NUMBER WITH MANY ZEROS FIRST
          IF (XI.EQ.0) GO TO 40
          CALL MPMLP(R(I+1),Y(3),XI,MIN0(T,I2-I))
          C = C - 1
          IF (C.GT.0) GO TO 40
C CHECK FOR LEGAL BASE B DIGIT
          IF ((XI.LT.0) .OR. (XI.GE.B)) GO TO 80
C PROPAGATE CARRIES AT END AND EVERY EIGHTH TIME,
C FASTER THAN DOING IT EVERY TIME.
          DO 30 J = 1,I2
              J1 = I2P - J
              RI = R(J1) + C
              IF (RI.LT.0) GO TO 70
              C = RI/B
              R(J1) = RI - B*C
   30     CONTINUE
          IF (C.NE.0) GO TO 80
          C = 8
   40 CONTINUE
      IF (C.EQ.8) GO TO 60
      IF ((XI.LT.0) .OR. (XI.GE.B)) GO TO 80
      C = 0
      DO 50 J = 1,I2
          J1 = I2P - J
          RI = R(J1) + C
          IF (RI.LT.0) GO TO 70
          C = RI/B
          R(J1) = RI - B*C
   50 CONTINUE
      IF (C.NE.0) GO TO 80
C NORMALIZE AND ROUND RESULT
   60 CALL MPNZR(RS,RE,Z,0)
      RETURN

   70 WRITE (LUN,FMT=9000)
      GO TO 90

   80 WRITE (LUN,FMT=9010)
   90 CALL MPERR
      Z(1) = 0
      RETURN

 9000 FORMAT (' *** INTEGER OVERFLOW IN MPMUL, B TOO LARGE ***')
 9010 FORMAT (' *** ILLEGAL BASE B DIGIT IN CALL TO MPMUL,',' POSSIBLE',
     +       ' OVERWRITING PROBLEM ***')
      END
      SUBROUTINE MPMULI(X,IY,Z)
C MULTIPLIES MP X BY SINGLE-PRECISION INTEGER IY GIVING MP Z.
C THIS IS FASTER THAN USING MPMUL.  RESULT IS ROUNDED.
C MULTIPLICATION BY 1 MAY BE USED TO NORMALIZE A NUMBER
C EVEN IF THE LAST DIGIT IS B.
C     .. Scalar Arguments ..
      INTEGER IY
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Z(*)
C     ..
C     .. External Subroutines ..
      EXTERNAL MPMUL2
C     ..
      CALL MPMUL2(X,IY,Z,0)
      RETURN

      END
      SUBROUTINE MPMULQ(X,I,J,Y)
C MULTIPLIES MP X BY I/J, GIVING MP Y
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER I,J
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER IS,JS
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPDIVI,MPERR,MPGCD,MPMUL2
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IABS
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      IF (J.NE.0) GO TO 10
      CALL MPCHK(1,4)
      WRITE (LUN,FMT=9000)
      CALL MPERR
      GO TO 20

   10 IF (I.NE.0) GO TO 30
   20 Y(1) = 0
      RETURN
C REDUCE TO LOWEST TERMS
   30 IS = I
      JS = J
      CALL MPGCD(IS,JS)
      IF (IABS(IS).EQ.1) GO TO 40
      CALL MPDIVI(X,JS,Y)
      CALL MPMUL2(Y,IS,Y,0)
      RETURN
C HERE IS = +-1
   40 CALL MPDIVI(X,IS*JS,Y)
      RETURN

 9000 FORMAT (' *** ATTEMPTED DIVISION BY ZERO IN MPMULQ ***')
      END
      SUBROUTINE MPMUL2(X,IY,Z,TRUNC)
C MULTIPLIES MP X BY SINGLE-PRECISION INTEGER IY GIVING MP Z.
C MULTIPLICATION BY 1 MAY BE USED TO NORMALIZE A NUMBER
C EVEN IF SOME DIGITS ARE GREATER THAN B-1.
C RESULT IS ROUNDED IF TRUNC.EQ.0, OTHERWISE TRUNCATED.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER IY,TRUNC
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Z(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER C,C1,C2,I,IJ,IS,IX,J,J1,J2,RE,RI,RS,T1,T3,T4
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPERR,MPNZR,MPOVFL,MPSTR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX0
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      RS = X(1)
      IF (RS.EQ.0) GO TO 10
      J = IY
      IF (J) 20,10,40
C RESULT ZERO
   10 Z(1) = 0
      RETURN

   20 J = -J
      RS = -RS
C CHECK FOR MULTIPLICATION BY B
      IF (J.NE.B) GO TO 40
      IF (X(2).LT.M) GO TO 30
      CALL MPCHK(1,4)
      WRITE (LUN,FMT=9000)
      CALL MPOVFL(Z)
      RETURN

   30 CALL MPSTR(X,Z)
      Z(1) = RS
      Z(2) = X(2) + 1
      RETURN
C SET EXPONENT TO EXPONENT(X) + 4
   40 RE = X(2) + 4
C FORM PRODUCT IN ACCUMULATOR
      C = 0
      T1 = T + 1
      T3 = T + 3
      T4 = T + 4
C IF J*B NOT REPRESENTABLE AS AN INTEGER WE HAVE TO SIMULATE
C DOUBLE-PRECISION MULTIPLICATION.
      IF (J.GE.MAX0(8*B,32767/B)) GO TO 100
      DO 50 IJ = 1,T
          I = T1 - IJ
          RI = J*X(I+2) + C
          C = RI/B
          R(I+4) = RI - B*C
   50 CONTINUE
C CHECK FOR INTEGER OVERFLOW
      IF (RI.LT.0) GO TO 120
C HAVE TO TREAT FIRST FOUR WORDS OF R SEPARATELY
      DO 60 IJ = 1,4
          I = 5 - IJ
          RI = C
          C = RI/B
          R(I) = RI - B*C
   60 CONTINUE
      IF (C.EQ.0) GO TO 90
C HAVE TO SHIFT RIGHT HERE AS CARRY OFF END
   70 DO 80 IJ = 1,T3
          I = T4 - IJ
          R(I+1) = R(I)
   80 CONTINUE
      RI = C
      C = RI/B
      R(1) = RI - B*C
      RE = RE + 1
      IF (C) 120,90,70
C NORMALIZE AND ROUND OR TRUNCATE RESULT
   90 CALL MPNZR(RS,RE,Z,TRUNC)
      RETURN
C HERE J IS TOO LARGE FOR SINGLE-PRECISION MULTIPLICATION
  100 J1 = J/B
      J2 = J - J1*B
C FORM PRODUCT
      DO 110 IJ = 1,T4
          C1 = C/B
          C2 = C - B*C1
          I = T1 - IJ
          IX = 0
          IF (I.GT.0) IX = X(I+2)
          RI = J2*IX + C2
          IS = RI/B
          C = J1*IX + C1 + IS
          R(I+4) = RI - B*IS
  110 CONTINUE
      IF (C) 120,90,70
C CAN ONLY GET HERE IF INTEGER OVERFLOW OCCURRED
  120 CALL MPCHK(1,4)
      WRITE (LUN,FMT=9010)
      CALL MPERR
      GO TO 10

 9000 FORMAT (' *** OVERFLOW OCCURRED IN MPMUL2 ***')
 9010 FORMAT (' *** INTEGER OVERFLOW IN MPMUL2, B TOO LARGE ***')
      END
      LOGICAL FUNCTION MPNE(X,Y)
C RETURNS LOGICAL VALUE OF (X .NE. Y) FOR MP X AND Y.
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. External Functions ..
      INTEGER MPCOMP
      EXTERNAL MPCOMP
C     ..
      MPNE = (MPCOMP(X,Y).NE.0)
      RETURN

      END
      SUBROUTINE MPNEG(X,Y)
C SETS Y = -X FOR MP NUMBERS X AND Y
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. External Subroutines ..
      EXTERNAL MPSTR
C     ..
      CALL MPSTR(X,Y)
      Y(1) = -Y(1)
      RETURN

      END
      SUBROUTINE MPNZR(RS,RE,Z,TRUNC)
C ASSUMES LONG (I.E. (T+4)-DIGIT) FRACTION IN
C R, SIGN = RS, EXPONENT = RE.  NORMALIZES,
C AND RETURNS MP RESULT IN Z.  INTEGER ARGUMENTS RS AND RE
C ARE NOT PRESERVED. R*-ROUNDING IS USED IF TRUNC.EQ.0
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER RE,RS,TRUNC
C     ..
C     .. Array Arguments ..
      INTEGER Z(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER B2,I,I2,I2M,I2P,IS,IT,J,K
C     ..
C     .. External Subroutines ..
      EXTERNAL MPERR,MPOVFL,MPUNFL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IABS,MOD
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      I2 = T + 4
      IF (RS.NE.0) GO TO 20
C STORE ZERO IN Z
   10 Z(1) = 0
      RETURN
C CHECK THAT SIGN = +-1
   20 IF (IABS(RS).LE.1) GO TO 30
      WRITE (LUN,FMT=9000)
      CALL MPERR
      GO TO 10
C LOOK FOR FIRST NONZERO DIGIT
   30 DO 40 I = 1,I2
          IS = I - 1
          IF (R(I).GT.0) GO TO 50
   40 CONTINUE
C FRACTION ZERO
      GO TO 10

   50 IF (IS.EQ.0) GO TO 80
C NORMALIZE
      RE = RE - IS
      I2M = I2 - IS
      DO 60 J = 1,I2M
          K = J + IS
          R(J) = R(K)
   60 CONTINUE
      I2P = I2M + 1
      DO 70 J = I2P,I2
          R(J) = 0
   70 CONTINUE
C CHECK TO SEE IF TRUNCATION IS DESIRED
   80 IF (TRUNC.NE.0) GO TO 140
C SEE IF ROUNDING NECESSARY
C TREAT EVEN AND ODD BASES DIFFERENTLY
      B2 = B/2
      IF ((2*B2).NE.B) GO TO 120
C B EVEN.  ROUND IF R(T+1).GE.B2 UNLESS R(T) ODD AND ALL ZEROS
C AFTER R(T+2).
      IF (R(T+1)-B2) 140,90,100
   90 IF (MOD(R(T),2).EQ.0) GO TO 100
      IF ((R(T+2)+R(T+3)+R(T+4)).EQ.0) GO TO 140
C ROUND
  100 DO 110 J = 1,T
          I = T + 1 - J
          R(I) = R(I) + 1
          IF (R(I).LT.B) GO TO 140
          R(I) = 0
  110 CONTINUE
C EXCEPTIONAL CASE, ROUNDED UP TO .10000...
      RE = RE + 1
      R(1) = 1
      GO TO 140
C ODD BASE, ROUND IF R(T+1)... .GT. 1/2
  120 DO 130 I = 1,4
          IT = T + I
          IF (R(IT)-B2) 140,130,100
  130 CONTINUE
C CHECK FOR OVERFLOW
  140 IF (RE.LE.M) GO TO 150
      WRITE (LUN,FMT=9010)
      CALL MPOVFL(Z)
      RETURN
C CHECK FOR UNDERFLOW
  150 IF (RE.LT. (-M)) GO TO 170
C STORE RESULT IN Z
      Z(1) = RS
      Z(2) = RE
      DO 160 I = 1,T
          Z(I+2) = R(I)
  160 CONTINUE
      RETURN
C UNDERFLOW HERE
  170 CALL MPUNFL(Z)
      RETURN

 9000 FORMAT (' *** SIGN NOT 0, +1 OR -1 IN CALL TO MPNZR,',' POSSIBLE',
     +       ' OVERWRITING PROBLEM ***')
 9010 FORMAT (' *** OVERFLOW OCCURRED IN MPNZR ***')
      END
      SUBROUTINE MPOUT(X,C,P,N)
C CONVERTS MULTIPLE-PRECISION X TO FP.N FORMAT IN C,
C WHICH MAY BE PRINTED UNDER PA1 FORMAT.  NOTE THAT
C N = -1 IS ALLOWED, AND EFFECTIVELY GIVES IP FORMAT.
C DIGITS AFTER THE DECIMAL POINT ARE BLANKED OUT IF
C THEY COULD NOT BE SIGNIFICANT.
C EFFICIENCY IS HIGHER IF B IS A POWER OF 10 THAN IF NOT.
C DIMENSION OF C MUST BE AT LEAST P.
C C IS AN INTEGER ARRAY, P AND N ARE INTEGERS.
C DIMENSION OF R IN COMMON MUST BE AT LEAST 3T+11
C     .. Scalar Arguments ..
      INTEGER N,P
C     ..
C     .. Array Arguments ..
      INTEGER C(*),X(*)
C     ..
C     .. External Subroutines ..
      EXTERNAL MPOUT2
C     ..
      CALL MPOUT2(X,C,P,N,10)
      RETURN

      END
      SUBROUTINE MPOUTE(X,C,J,P)
C ASSUMES X IS AN MP NUMBER AND C AN INTEGER ARRAY OF DIMENSION AT
C LEAST P .GE. 4.  ON RETURN J IS THE EXPONENT (TO BASE TEN) OF X
C AND THE FRACTION IS IN C, READY TO BE PRINTED IN A1 FORMAT.
C FOR EXAMPLE, WE COULD PRINT J AND C IN I10, 1X, PA1 FORMAT.
C THE FRACTION HAS ONE PLACE BEFORE DECIMAL POINT AND P-3 AFTER.
C J AND P ARE INTEGERS.    SPACE = 6T+14
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER J,P
C     ..
C     .. Array Arguments ..
      INTEGER C(*),X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I2,IBL,IMIN,IR
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPCIM,MPCMEF,MPERR,MPOUT
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
C     .. Data statements ..
      DATA IBL/1H /,IMIN/1H-/
C     ..
C CHECK LEGALITY OF B, T, M, LUN AND MXR
      CALL MPCHK(6,14)
      IF (P.GE.4) GO TO 10
      WRITE (LUN,FMT=9000)
      CALL MPERR
      RETURN

   10 I2 = 5*T + 13
      CALL MPCMEF(X,J,R(I2))
      CALL MPOUT(R(I2),C,P,P-3)
C SEE IF OUTPUT OF MPOUT WAS ROUNDED UP TO TEN
      IF ((C(1).EQ.IBL) .OR. (C(1).EQ.IMIN)) RETURN
C IT WAS, SO ADD 1 TO J AND CONVERT SIGN TO MP
      J = J + 1
C AVOID POSSIBLY UNSAFE REFERENCE (SEE SOFTWARE PRACTICE
C AND EXPERIENCE, VOL. 4, 359-378).
      IR = R(I2)
      CALL MPCIM(IR,R(I2))
      CALL MPOUT(R(I2),C,P,P-3)
      RETURN

 9000 FORMAT (' *** P .LT. 4 IN CALL TO MPOUTE ***')
      END
      SUBROUTINE MPOUTF(X,P,N,IFORM,ERR)
C WRITES MP NUMBER X ON LOGICAL UNIT LUN (FOURTH WORD OF COMMON)
C IN FORMAT IFORM AFTER CONVERTING TO FP.N DECIMAL REPRESENTATION
C USING ROUTINE MPOUT. FOR FURTHER DETAILS SEE COMMENTS IN MPOUT.
C IFORM SHOULD CONTAIN A FORMAT WHICH ALLOWS FOR OUTPUT OF P
C WORDS IN A1 FORMAT, PLUS ANY DESIRED HEADINGS, SPACING ETC.
C E.G. 24H(8H1HEADING/(11X,100A1))
C ERR RETURNED AS TRUE IF P NOT POSITIVE, OTHERWISE FALSE.
C SPACE REQUIRED 3T+P+11 WORDS.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER N,P
      LOGICAL ERR
C     ..
C     .. Array Arguments ..
      INTEGER IFORM(*),X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I2
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPIO,MPOUT
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      ERR = .TRUE.
C RETURN WITH ERROR FLAG SET IF OUTPUT FIELD WIDTH P NOT POSITIVE
      IF (P.LE.0) RETURN
C CHECK THAT ENOUGH SPACE IS AVAILABLE
      CALL MPCHK(3,P+11)
      I2 = 3*T + 12
C CONVERT X TO DECIMAL FORM
      CALL MPOUT(X,R(I2),P,N)
C AND WRITE ON UNIT LUN WITH FORMAT IFORM
      CALL MPIO(R(I2),P,LUN,IFORM,ERR)
      RETURN

      END
      SUBROUTINE MPOUT2(X,C,P,N,NB)

C SAME AS MPOUT EXCEPT THAT OUTPUT REPRESENTATION IS IN
C BASE NB, WHERE 2 .LE. NB .LE. 16,
C EG NB = 8 GIVES OCTAL OUTPUT, NB = 16 GIVES HEXADECIMAL.
C OUTPUT DIGITS ARE 0123456789ABCDEF.
C X IS AN MP NUMBER, C AN INTEGER ARRAY, P, N AND NB ARE INTEGERS.
C DIMENSION OF C MUST BE AT LEAST P, DIMENSION OF R
C IN CALLING PROGRAM MUST BE AT LEAST 3T+11
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER N,NB,P
C     ..
C     .. Array Arguments ..
      INTEGER C(*),X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I,I1,I2,I3,IB,IP,IS,ISZ,ITP,IZ,J,JD,JP,NMAX,NP,TP
C     ..
C     .. Local Arrays ..
      INTEGER NUM(20)
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADD2,MPCHK,MPCLR,MPCMF,MPCMI,MPCQM,MPDIVI,MPERR,MPMUL2,
     +         MPSTR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ALOG,AMAX1,FLOAT,MAX0,MIN1
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
C     .. Data statements ..
      DATA NUM(1),NUM(2),NUM(3)/1H0,1H1,1H2/
      DATA NUM(4),NUM(5),NUM(6)/1H3,1H4,1H5/
      DATA NUM(7),NUM(8),NUM(9)/1H6,1H7,1H8/
      DATA NUM(10),NUM(11),NUM(12)/1H9,1HA,1HB/
      DATA NUM(13),NUM(14),NUM(15)/1HC,1HD,1HE/
      DATA NUM(16),NUM(17),NUM(18)/1HF,1H ,1H./
      DATA NUM(19),NUM(20)/1H*,1H-/
C     ..
C CHECK LEGALITY OF B, T, M, MXR AND LUN
      CALL MPCHK(3,11)
C CHECK LEGALITY OF P, N AND NB
      IF ((N.GE. (-1)) .AND. (P.GT.N) .AND. (P.GT.0) .AND.
     +    (NB.GT.1) .AND. (NB.LE.16)) GO TO 10
      WRITE (LUN,FMT=9000)
      CALL MPERR
      RETURN
C COMPUTE DISPLACEMENTS, MOVE X
   10 I2 = T + 6
      I3 = I2 + T + 3
      CALL MPSTR(X,R(I2))
      NP = P - N
      IP = NP - 1
C COMPUTE POWER OF NB WHICH WE CAN SAFELY MULTIPLY AND
C DIVIDE BY (THIS SAVES TIME).
      IB = MAX0(7*B*B,32767)/NB
      TP = NB
      ITP = 1
   20 IF ((TP.GT.IB) .OR. (TP.EQ.B)) GO TO 30
      TP = TP*NB
      ITP = ITP + 1
      GO TO 20
C PUT FORMATTED ZERO IN C
   30 DO 40 I = 1,P
          C(I) = NUM(17)
          IF (I.GE.IP) C(I) = NUM(1)
          IF (I.EQ.NP) C(I) = NUM(18)
   40 CONTINUE
C GET SIGN OF X, CHECK FOR ZERO
      IS = R(I2)
      IF (IS.EQ.0) RETURN
      R(I2) = 1
C COMPUTE MAXIMUM NUMBER OF NONZERO DIGITS WHICH WE CAN
C MEANINGFULLY GIVE AFTER DECIMAL POINT.
      NMAX = MIN1(FLOAT(N)+0.001E0,AMAX1(0E0,
     +       FLOAT(T-R(I2+1))*ALOG(FLOAT(B))/ALOG(FLOAT(NB))+0.001E0))
C WORK WITH ONE GUARD DIGIT
      CALL MPCLR(R(I2),T+1)
      T = T + 1
C COMPUTE ROUNDING CONSTANT
      CALL MPCQM(1,2,R(I3))
      IF (NMAX.LE.0) GO TO 60
      JP = 1
      DO 50 I = 1,NMAX
          JP = NB*JP
          IF ((JP.LE.IB) .AND. (JP.NE.B) .AND. (I.LT.NMAX)) GO TO 50
          CALL MPDIVI(R(I3),JP,R(I3))
          JP = 1
   50 CONTINUE
C ADD ROUNDING CONSTANT TO ABS(X), TRUNCATING RESULT
   60 CALL MPADD2(R(I2),R(I3),R(I2),R(I3),1)
C IP PLACES BEFORE POINT, SO DIVIDE BY NB**IP
      IF (IP.LE.0) GO TO 80
      JP = 1
      DO 70 I = 1,IP
          JP = NB*JP
          IF ((JP.LE.IB) .AND. (I.LT.IP)) GO TO 70
          CALL MPDIVI(R(I2),JP,R(I2))
          JP = 1
   70 CONTINUE
   80 IZ = 0
C CHECK THAT NUMBER IS LESS THAN ONE
      IF (R(I2+1).GT.0) GO TO 160
      IF (IP.LE.0) GO TO 130
C PUT DIGITS BEFORE POINT IN
      JD = 1
      DO 120 I = 1,IP
          IF (JD.GT.1) GO TO 110
          IF ((I+ITP).LE. (IP+1)) GO TO 90
C MULTIPLY BY NB, TRUNCATING RESULT
          CALL MPMUL2(R(I2),NB,R(I2),1)
          JD = NB
          GO TO 100
C HERE WE CAN MULTIPLY BY A POWER OF NB TO SAVE TIME
   90     CALL MPMUL2(R(I2),TP,R(I2),1)
          JD = TP
C GET INTEGER PART
  100     CALL MPCMI(R(I2),JP)
C AND FRACTIONAL PART
          CALL MPCMF(R(I2),R(I2))
  110     JD = JD/NB
C GET NEXT DECIMAL DIGIT
          J = JP/JD
          JP = JP - J*JD
          ISZ = IZ
          IF ((J.GT.0) .OR. (I.EQ.IP)) IZ = 1
          IF (IZ.GT.0) C(I) = NUM(J+1)
          IF ((IZ.EQ.ISZ) .OR. (IS.GT.0)) GO TO 120
          IF (I.EQ.1) GO TO 160
          C(I-1) = NUM(20)
  120 CONTINUE
  130 IF (NMAX.LE.0) GO TO 180
C PUT IN DIGITS AFTER DECIMAL POINT
      JD = 1
      DO 150 I = 1,NMAX
          IF (JD.GT.1) GO TO 140
          CALL MPMUL2(R(I2),TP,R(I2),1)
          CALL MPCMI(R(I2),JP)
          CALL MPCMF(R(I2),R(I2))
          JD = TP
  140     JD = JD/NB
          J = JP/JD
          JP = JP - J*JD
          I1 = NP + I
          C(I1) = NUM(J+1)
  150 CONTINUE
      GO TO 180
C ERROR OCCURRED, RETURN ASTERISKS.
  160 DO 170 I = 1,P
          C(I) = NUM(19)
  170 CONTINUE
C RESTORE T
  180 T = T - 1
C BLANK OUT ANY NONSIGNIFICANT TRAILING ZEROS
      IF ((NMAX.GE.N) .OR. (C(1).EQ.NUM(19))) RETURN
      I1 = NP + NMAX + 1
      DO 190 I = I1,P
          C(I) = NUM(17)
  190 CONTINUE
      RETURN

 9000 FORMAT (' *** PARAMETERS P, N AND/OR NB',' ILLEGAL IN CALL TO SU',
     +       'BROUTINE MPOUT2 ***')
      END
      SUBROUTINE MPOVFL(X)
C CALLED ON MULTIPLE-PRECISION OVERFLOW, IE WHEN THE
C EXPONENT OF MP NUMBER X WOULD EXCEED M.
C AT PRESENT EXECUTION IS TERMINATED WITH AN ERROR MESSAGE
C AFTER CALLING MPMAXR(X), BUT IT WOULD BE POSSIBLE TO RETURN,
C POSSIBLY UPDATING A COUNTER AND TERMINATING EXECUTION AFTER
C A PRESET NUMBER OF OVERFLOWS.  ACTION COULD EASILY BE DETERMINED
C BY A FLAG IN LABELLED COMMON.
C M MAY HAVE BEEN OVERWRITTEN, SO CHECK B, T, M ETC.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPERR,MPMAXR
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(1,4)
C SET X TO LARGEST POSSIBLE POSITIVE NUMBER
      CALL MPMAXR(X)
      WRITE (LUN,FMT=9000)
C TERMINATE EXECUTION BY CALLING MPERR
      CALL MPERR
      RETURN

 9000 FORMAT (' *** CALL TO MPOVFL, MP OVERFLOW OCCURRED ***')
      END
      SUBROUTINE MPPACK(X,Y)
C ASSUMES THAT X IS AN MP NUMBER STORED AS USUAL IN AN INTEGER
C ARRAY OF DIMENSION AT LEAST T+2, AND Y IS AN INTEGER ARRAY
C OF DIMENSION AT LEAST INT((T+3)/2).
C X IS STORED IN A COMPACT FORMAT IN Y, AND MAY BE RETRIEVED
C BY CALLING MPUNPK (Y, X).
C MPPACK AND MPUNPK ARE USEFUL IF SPACE IS CRITICAL, FOR EXAMPLE
C WHEN WORKING WITH LARGE ARRAYS OF MP NUMBERS.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I,IS,J
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      IS = X(1)
      IF (IS.NE.0) GO TO 10
C X ZERO HERE
      Y(2) = 0
      RETURN
C X NONZERO.  FIRST MOVE EXPONENT TO Y(1).
   10 Y(1) = X(2)
      J = T/2
C NOW PACK TWO DIGITS OF X IN EACH WORD OF Y.
      DO 20 I = 1,J
          Y(I+1) = B*X(2*I+1) + X(2*I+2)
   20 CONTINUE
C FIX UP LAST DIGIT IF T ODD, AND CORRECT SIGN.
      IF ((2*J).LT.T) Y(J+2) = B*X(T+2)
      Y(2) = IS*Y(2)
      RETURN

      END
      SUBROUTINE MPPI(X)
C SETS MP X = PI TO THE AVAILABLE PRECISION.
C USES PI/4 = 4.ARCTAN(1/5) - ARCTAN(1/239).
C TIME IS O(T**2).
C DIMENSION OF R MUST BE AT LEAST 3T+8
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      REAL RX
      INTEGER I2
C     ..
C     .. External Subroutines ..
      EXTERNAL MPART1,MPCHK,MPCMR,MPERR,MPMULI,MPSUB
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(3,8)
C ALLOW SPACE FOR MPART1
      I2 = 2*T + 7
      CALL MPART1(5,R(I2))
      CALL MPMULI(R(I2),4,R(I2))
      CALL MPART1(239,X)
      CALL MPSUB(R(I2),X,X)
      CALL MPMULI(X,4,X)
C RETURN IF ERROR IS LESS THAN 0.01
      CALL MPCMR(X,RX)
      IF (ABS(RX-3.1416).LT.0.01) RETURN
      WRITE (LUN,FMT=9000)
C FOLLOWING MESSAGE MAY INDICATE THAT B**(T-1) IS TOO SMALL
      CALL MPERR
      RETURN

 9000 FORMAT (' *** ERROR OCCURRED IN MPPI, RESULT INCORRECT ***')
      END
      SUBROUTINE MPPIGL(PI)
C SETS MP PI = 3.14159... TO THE AVAILABLE PRECISION.
C USES THE GAUSS-LEGENDRE ALGORITHM.
C THIS METHOD REQUIRES TIME O(LN(T)M(T)), SO IT IS SLOWER
C THAN MPPI IF M(T) = O(T**2), BUT WOULD BE FASTER FOR
C LARGE T IF A FASTER MULTIPLICATION ALGORITHM WERE USED
C (SEE COMMENTS IN MPMUL).
C FOR A DESCRIPTION OF THE METHOD, SEE - MULTIPLE-PRECISION
C ZERO-FINDING AND THE COMPLEXITY OF ELEMENTARY FUNCTION
C EVALUATION (BY R. P. BRENT), IN ANALYTIC COMPUTATIONAL
C COMPLEXITY (EDITED BY J. F. TRAUB), ACADEMIC PRESS, 1976,
C 151-176.  DIMENSION OF R MUST BE AT LEAST 6T+14
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER PI(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I2,I3,I4,IX,R1,R2
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADD,MPCHK,MPCIM,MPCQM,MPDIV,MPDIVI,MPMUL,MPMULI,MPSQRT,
     +         MPSTR,MPSUB
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(6,14)
C MPSQRT AND MPDIV REQUIRE SPACE 4T+10, BUT CAN OVERLAP
      I2 = 3*T + 9
      I3 = I2 + T + 2
      I4 = I3 + T + 2
      CALL MPCIM(1,PI(1))
      CALL MPCQM(1,2,R(I4))
      CALL MPSQRT(R(I4),R(I4))
      CALL MPCQM(1,4,R(I3))
      IX = 1
   10 CALL MPSTR(PI(1),R(I2))
      CALL MPADD(PI(1),R(I4),PI(1))
      CALL MPDIVI(PI(1),2,PI(1))
      CALL MPMUL(R(I2),R(I4),R(I4))
      CALL MPSUB(PI(1),R(I2),R(I2))
      CALL MPMUL(R(I2),R(I2),R(I2))
      CALL MPMULI(R(I2),IX,R(I2))
      CALL MPSUB(R(I3),R(I2),R(I3))
C SAVE ARRAY ELEMENTS WHICH WILL BE OVERWRITTEN BY MPSQRT
      R1 = R(I2)
      R2 = R(I2+1)
      CALL MPSQRT(R(I4),R(I4))
      IX = 2*IX
C CHECK FOR CONVERGENCE
      IF ((R1.NE.0) .AND. (R2.GE. (-T))) GO TO 10
      CALL MPMUL(PI(1),R(I4),PI(1))
      CALL MPDIV(PI(1),R(I3),PI(1))
      RETURN

      END
      SUBROUTINE MPPOLY(X,Y,IC,N)
C SETS Y = IC(1) + IC(2)*X + ... + IC(N)*X**(N-1),
C WHERE X AND Y ARE MULTIPLE-PRECISION NUMBERS AND
C IC IS AN INTEGER ARRAY OF DIMENSION AT LEAST N .GT. 0
C DIMENSION OF R IN CALLING PROGRAM MUST BE AT LEAST 3T+8
C (BUT Y(1) MAY BE R(2T+7))
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      INTEGER IC(*),X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I,I2
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADDI,MPCHK,MPCIM,MPERR,MPMUL,MPSTR
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(3,8)
      I2 = 2*T + 7
      IF (N.GT.0) GO TO 10
      WRITE (LUN,FMT=9000)
      CALL MPERR
      Y(1) = 0
      RETURN

   10 CALL MPCIM(IC(N),R(I2))
      I = N - 1
      IF (I.LE.0) GO TO 30
   20 CALL MPMUL(R(I2),X,R(I2))
      CALL MPADDI(R(I2),IC(I),R(I2))
      I = I - 1
      IF (I.GT.0) GO TO 20
   30 CALL MPSTR(R(I2),Y)
      RETURN

 9000 FORMAT (' *** N NOT POSITIVE IN CALL TO MPPOLY ***')
      END
      SUBROUTINE MPPWR(X,N,Y)
C RETURNS Y = X**N, FOR MP X AND Y, INTEGER N, WITH 0**0 = 1.
C R MUST BE DIMENSIONED AT LEAST 4T+10 IN CALLING PROGRAM
C (2T+6 IS ENOUGH IF N NONNEGATIVE).
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I2,N2,NS
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPCIM,MPERR,MPMUL,MPREC,MPSTR
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      I2 = T + 5
      N2 = N
      IF (N2) 20,10,30
C N = 0, RETURN Y = 1.
   10 CALL MPCIM(1,Y)
      RETURN
C N .LT. 0
   20 CALL MPCHK(4,10)
      N2 = -N2
      IF (X(1).NE.0) GO TO 50
      WRITE (LUN,FMT=9000)
      CALL MPERR
      GO TO 40
C N .GT. 0
   30 CALL MPCHK(2,6)
      IF (X(1).NE.0) GO TO 50
C X = 0, N .GT. 0, SO Y = 0
   40 Y(1) = 0
      RETURN
C MOVE X
   50 CALL MPSTR(X,Y)
C IF N .LT. 0 FORM RECIPROCAL
      IF (N.LT.0) CALL MPREC(Y,Y)
      CALL MPSTR(Y,R(I2))
C SET PRODUCT TERM TO ONE
      CALL MPCIM(1,Y)
C MAIN LOOP, LOOK AT BITS OF N2 FROM RIGHT
   60 NS = N2
      N2 = N2/2
      IF ((2*N2).NE.NS) CALL MPMUL(Y,R(I2),Y)
      IF (N2.LE.0) RETURN
      CALL MPMUL(R(I2),R(I2),R(I2))
      GO TO 60

 9000 FORMAT (' *** ATTEMPT TO RAISE ZERO TO NEGATIVE POWER IN',' CALL',
     +       ' TO SUBROUTINE MPPWR ***')
      END
      SUBROUTINE MPPWR2(X,Y,Z)
C RETURNS Z = X**Y FOR MP NUMBERS X, Y AND Z, WHERE X IS
C POSITIVE (X .EQ. 0 ALLOWED IF Y .GT. 0).  SLOWER THAN
C MPPWR AND MPQPWR, SO USE THEM IF POSSIBLE.
C DIMENSION OF R IN COMMON AT LEAST 7T+16
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*),Z(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I2
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPERR,MPEXP,MPLN,MPMUL
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(7,16)
      IF (X(1)) 10,20,50
   10 WRITE (LUN,FMT=9000)
      GO TO 30
C HERE X IS ZERO, RETURN ZERO IF Y POSITIVE, OTHERWISE ERROR
   20 IF (Y(1).GT.0) GO TO 40
      WRITE (LUN,FMT=9010)
   30 CALL MPERR
C RETURN ZERO HERE
   40 Z(1) = 0
      RETURN
C USUAL CASE HERE, X POSITIVE
C USE MPLN AND MPEXP TO COMPUTE POWER
   50 I2 = 6*T + 15
      CALL MPLN(X,R(I2))
      CALL MPMUL(Y,R(I2),Z)
C IF X**Y IS TOO LARGE, MPEXP WILL PRINT ERROR MESSAGE
      CALL MPEXP(Z,Z)
      RETURN

 9000 FORMAT (' *** X NEGATIVE IN CALL TO MPPWR2 ***')
 9010 FORMAT (' *** X ZERO AND Y NONPOSITIVE IN CALL TO MPPWR2 ***')
      END
      SUBROUTINE MPQPWR(I,J,K,L,X)
C SETS MULTIPLE-PRECISION X = (I/J)**(K/L) FOR INTEGERS
C I, J, K AND L.  USES MPROOT IF ABS(L) SMALL, OTHERWISE
C USES MPLNI AND MPEXP.
C SPACE (DIMENSION OF R IN COMMON) = 4T+10
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER I,J,K,L
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I2,IS,JS,KS,LS
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPCIM,MPCQM,MPERR,MPEXP,MPGCD,MPLNI,MPMULQ,MPPWR,
     +         MPROOT,MPSUB
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IABS,MAX0,MOD
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(4,10)
      IS = I
      JS = J
      KS = K
      LS = L
C FOR EFFICIENCY MAKE KS POSITIVE AND LS NEGATIVE
C (SEE COMMENTS IN MPROOT AND MPPWR).
      IF (KS) 30,10,40
C (I/J)**(0/L) = 1 IF J AND L ARE NONZERO.
   10 CALL MPCIM(1,X)
      IF ((JS.NE.0) .AND. (LS.NE.0)) RETURN
   20 WRITE (LUN,FMT=9000)
      GO TO 70

   30 KS = -KS
      LS = -LS
C NOW KS IS POSITIVE, SO LOOK AT LS
   40 IF (LS) 60,20,50
C LS POSITIVE SO INTERCHANGE IS AND JS TO MAKE NEGATIVE
   50 IS = J
      JS = I
      LS = -LS
C NOW KS POSITIVE, LS NEGATIVE
   60 IF (IS.NE.0) GO TO 80
      WRITE (LUN,FMT=9010)
   70 CALL MPERR
      X(1) = 0
      RETURN

   80 X(1) = 0
C (I/0)**(NEGATIVE) = 0 IF I NONZERO
      IF (JS.EQ.0) RETURN
C TO SAVE TIME IN MPROOT AND MPPWR, FIND GCD OF KS AND LS
      CALL MPGCD(KS,LS)
C CHECK FOR LS = -1, TREAT AS SPECIAL CASE
      IF (LS.NE. (-1)) GO TO 90
      CALL MPCQM(JS,IS,X)
      GO TO 100
C USUAL CASE HERE, LS .NE. -1
   90 CALL MPCQM(IS,JS,X)
C USE MPROOT IF ABS(LS) .LE. MAX(B,64), OTHERWISE LOG AND EXP
      IF (IABS(LS).GT.MAX0(B,64)) GO TO 110
      CALL MPROOT(X,LS,X)
  100 CALL MPPWR(X,KS,X)
      RETURN
C HERE USE LOG AND EXP (SLOWER THAN MPROOT)
  110 I2 = 3*T + 9
      CALL MPLNI(IABS(IS),R(I2))
      CALL MPLNI(IABS(JS),X)
C SOME CANCELLATION BUT NOT SERIOUS HERE
      CALL MPSUB(R(I2),X,X)
      CALL MPMULQ(X,KS,LS,X)
      CALL MPEXP(X,X)
C CORRECT SIGN IF NECESSARY
      IF (((IS.GE.0).AND. (JS.GE.0)) .OR.
     +    ((IS.LT.0).AND. (JS.LT.0))) RETURN
C HERE IS/JS IS NEGATIVE
      X(1) = -X(1)
      IF (MOD(LS,2).NE.0) RETURN
      WRITE (LUN,FMT=9020)
      GO TO 70

 9000 FORMAT (' *** J = 0 OR L = 0 IN CALL TO MPQPWR ***')
 9010 FORMAT (' *** I = 0 AND K/L NEGATIVE OR',' J = 0 AND K/L POSITIV',
     +       'E IN CALL TO MPQPWR ***')
 9020 FORMAT (' *** I/J NEGATIVE AND L EVEN IN CALL TO MPQPWR ***')
      END
      SUBROUTINE MPREC(X,Y)
C RETURNS Y = 1/X, FOR MP X AND Y.
C MPROOT (X, -1, Y) HAS THE SAME EFFECT.
C DIMENSION OF R MUST BE AT LEAST 4*T+10 IN CALLING PROGRAM
C (BUT Y(1) MAY BE R(3T+9)).
C NEWTONS METHOD IS USED, SO FINAL ONE OR TWO DIGITS MAY
C NOT BE CORRECT.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      REAL RX
      INTEGER EX,I2,I3,IT0,TS,TS2,TS3
C     ..
C     .. Local Arrays ..
      INTEGER IT(9)
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADDI,MPCHK,MPCMR,MPCRM,MPERR,MPMUL,MPOVFL,MPSTR,MPSUB
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN0
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
C     .. Data statements ..
      DATA IT(1),IT(2),IT(3),IT(4),IT(5)/0,8,6,5,4/
      DATA IT(6),IT(7),IT(8),IT(9)/4,4,4,4/
C     ..
C CHECK LEGALITY OF B, T, M, LUN AND MXR
      CALL MPCHK(4,10)
C MPADDI REQUIRES 2T+6 WORDS.
      I2 = 2*T + 7
      I3 = I2 + T + 2
      IF (X(1).NE.0) GO TO 10
      WRITE (LUN,FMT=9000)
      CALL MPERR
      Y(1) = 0
      RETURN

   10 EX = X(2)
C TEMPORARILY INCREASE M TO AVOID OVERFLOW
      M = M + 2
C SET EXPONENT TO ZERO SO RX NOT TOO LARGE OR SMALL.
      X(2) = 0
      CALL MPCMR(X,RX)
C USE SINGLE-PRECISION RECIPROCAL AS FIRST APPROXIMATION
      CALL MPCRM(1E0/RX,R(I2))
C RESTORE EXPONENT
      X(2) = EX
C CORRECT EXPONENT OF FIRST APPROXIMATION
      R(I2+1) = R(I2+1) - EX
C SAVE T (NUMBER OF DIGITS)
      TS = T
C START WITH SMALL T TO SAVE TIME. ENSURE THAT B**(T-1) .GE. 100
      T = 3
      IF (B.LT.10) T = IT(B)
      IT0 = (T+4)/2
      IF (T.GT.TS) GO TO 50
C MAIN ITERATION LOOP
   20 CALL MPMUL(X,R(I2),R(I3))
      CALL MPADDI(R(I3),-1,R(I3))
C TEMPORARILY REDUCE T
      TS3 = T
      T = (T+IT0)/2
      CALL MPMUL(R(I2),R(I3),R(I3))
C RESTORE T
      T = TS3
      CALL MPSUB(R(I2),R(I3),R(I2))
      IF (T.GE.TS) GO TO 40
C FOLLOWING LOOP ALMOST DOUBLES T (POSSIBLE
C BECAUSE NEWTONS METHOD HAS 2ND ORDER CONVERGENCE).
      T = TS
   30 TS2 = T
      T = (T+IT0)/2
      IF (T.GT.TS3) GO TO 30
      T = MIN0(TS,TS2)
      GO TO 20
C RETURN IF NEWTON ITERATION WAS CONVERGING
   40 IF ((R(I3).EQ.0) .OR. ((2* (R(I2+1)-R(I3+1))).GE.
     +    (T-IT0))) GO TO 50
      WRITE (LUN,FMT=9010)
C THE FOLLOWING MESSAGE MAY INDICATE THAT B**(T-1) IS TOO SMALL,
C OR THAT THE STARTING APPROXIMATION IS NOT ACCURATE ENOUGH.
      CALL MPERR
C MOVE RESULT TO Y AND RETURN AFTER RESTORING T
   50 T = TS
      CALL MPSTR(R(I2),Y)
C RESTORE M AND CHECK FOR OVERFLOW (UNDERFLOW IMPOSSIBLE)
      M = M - 2
      IF (Y(2).LE.M) RETURN
      WRITE (LUN,FMT=9020)
      CALL MPOVFL(Y)
      RETURN

 9000 FORMAT (' *** ATTEMPTED DIVISION BY ZERO IN CALL TO MREC ***')
 9010 FORMAT (' *** ERROR OCCURRED IN MPREC, NEWTON ITERATION',' NOT C',
     +       'ONVERGING PROPERLY ***')
 9020 FORMAT (' *** OVERFLOW OCCURRED IN MPREC ***')
      END
      SUBROUTINE MPROOT(X,N,Y)
C RETURNS Y = X**(1/N) FOR INTEGER N, ABS(N) .LE. MAX (B, 64).
C AND MP NUMBERS X AND Y,
C USING NEWTONS METHOD WITHOUT DIVISIONS.   SPACE = 4T+10
C (BUT Y(1) MAY BE R(3T+9))
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      REAL RX
      INTEGER EX,I2,I3,IT0,NP,TS,TS2,TS3,XES
C     ..
C     .. Local Arrays ..
      INTEGER IT(9)
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADDI,MPCHK,MPCMR,MPCRM,MPDIVI,MPERR,MPMUL,MPPWR,MPSTR,
     +         MPSUB
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ALOG,EXP,FLOAT,IABS,MAX0,MIN0,MOD
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
C     .. Data statements ..
      DATA IT(1),IT(2),IT(3),IT(4),IT(5)/0,8,6,5,4/
      DATA IT(6),IT(7),IT(8),IT(9)/4,4,4,4/
C     ..
C CHECK LEGALITY OF B, T, M, LUN AND MXR
      CALL MPCHK(4,10)
      IF (N.NE.1) GO TO 10
      CALL MPSTR(X,Y)
      RETURN

   10 I2 = 2*T + 7
      I3 = I2 + T + 2
      IF (N.NE.0) GO TO 20
      WRITE (LUN,FMT=9000)
      GO TO 30

   20 NP = IABS(N)
C LOSS OF ACCURACY IF NP LARGE, SO ONLY ALLOW NP .LE. MAX (B, 64)
      IF (NP.LE.MAX0(B,64)) GO TO 40
      WRITE (LUN,FMT=9010)
   30 CALL MPERR
      Y(1) = 0
      RETURN
C LOOK AT SIGN OF X
   40 IF (X(1)) 60,50,70
C X = 0 HERE, RETURN 0 IF N POSITIVE, ERROR IF NEGATIVE
   50 Y(1) = 0
      IF (N.GT.0) RETURN
      WRITE (LUN,FMT=9020)
      GO TO 30
C X NEGATIVE HERE, SO ERROR IF N IS EVEN
   60 IF (MOD(NP,2).NE.0) GO TO 70
      WRITE (LUN,FMT=9030)
      GO TO 30
C GET EXPONENT AND DIVIDE BY NP
   70 XES = X(2)
      EX = XES/NP
C REDUCE EXPONENT SO RX NOT TOO LARGE OR SMALL.
      X(2) = 0
      CALL MPCMR(X,RX)
C USE SINGLE-PRECISION ROOT FOR FIRST APPROXIMATION
      CALL MPCRM(EXP((FLOAT(NP*EX-XES)*ALOG(FLOAT(B))-
     +           ALOG(ABS(RX)))/FLOAT(NP)),R(I2))
C SIGN OF APPROXIMATION SAME AS THAT OF X
      R(I2) = X(1)
C RESTORE EXPONENT
      X(2) = XES
C CORRECT EXPONENT OF FIRST APPROXIMATION
      R(I2+1) = R(I2+1) - EX
C SAVE T (NUMBER OF DIGITS)
      TS = T
C START WITH SMALL T TO SAVE TIME
      T = 3
C ENSURE THAT B**(T-1) .GE. 100
      IF (B.LT.10) T = IT(B)
      IF (T.GT.TS) GO TO 110
C IT0 IS A NECESSARY SAFETY FACTOR
      IT0 = (T+4)/2
C MAIN ITERATION LOOP
   80 CALL MPPWR(R(I2),NP,R(I3))
      CALL MPMUL(X,R(I3),R(I3))
      CALL MPADDI(R(I3),-1,R(I3))
C TEMPORARILY REDUCE T
      TS3 = T
      T = (T+IT0)/2
      CALL MPMUL(R(I2),R(I3),R(I3))
      CALL MPDIVI(R(I3),NP,R(I3))
C RESTORE T
      T = TS3
      CALL MPSUB(R(I2),R(I3),R(I2))
C FOLLOWING LOOP ALMOST DOUBLES T (POSSIBLE BECAUSE
C NEWTONS METHOD HAS 2ND ORDER CONVERGENCE).
      IF (T.GE.TS) GO TO 100
      T = TS
   90 TS2 = T
      T = (T+IT0)/2
      IF (T.GT.TS3) GO TO 90
      T = MIN0(TS,TS2)
      GO TO 80
C NOW R(I2) IS X**(-1/NP)
C CHECK THAT NEWTON ITERATION WAS CONVERGING
  100 IF ((R(I3).EQ.0) .OR. ((2* (R(I2+1)-R(I3+1))).GE.
     +    (T-IT0))) GO TO 110
      WRITE (LUN,FMT=9040)
C THE FOLLOWING MESSAGE MAY INDICATE THAT B**(T-1) IS TOO SMALL,
C OR THAT THE INITIAL APPROXIMATION OBTAINED USING ALOG AND EXP
C IS NOT ACCURATE ENOUGH.
      CALL MPERR
C RESTORE T
  110 T = TS
      IF (N.LT.0) GO TO 120
      CALL MPPWR(R(I2),N-1,R(I2))
      CALL MPMUL(X,R(I2),Y)
      RETURN

  120 CALL MPSTR(R(I2),Y)
      RETURN

 9000 FORMAT (' *** N = 0 IN CALL TO MPROOT ***')
 9010 FORMAT (' *** ABS(N) TOO LARGE IN CALL TO MPROOT ***')
 9020 FORMAT (' *** X = 0 AND N NEGATIVE IN CALL TO MPROOT ***')
 9030 FORMAT (' *** X NEGATIVE AND N EVEN IN CALL TO MPROOT ***')
 9040 FORMAT (' *** ERROR OCCURRED IN MPROOT, NEWTON ITERATION',' NOT ',
     +       'CONVERGING PROPERLY ***')
      END
      SUBROUTINE MPSIGB(I,X)
C SETS SIGN OF MP NUMBER X TO I.
C I SHOULD BE 0, +1 OR -1.
C EXPONENT AND DIGITS OF X ARE UNCHANGED,
C BUT RESULT MUST BE A VALID MP NUMBER.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER I
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. External Subroutines ..
      EXTERNAL MPERR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IABS
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      X(1) = I
C CHECK FOR VALID SIGN
      IF (IABS(I).LE.1) GO TO 10
      WRITE (LUN,FMT=9000)
      GO TO 20
C RETURN IF X ZERO
   10 IF (I.EQ.0) RETURN
C CHECK FOR VALID EXPONENT AND LEADING DIGIT
      IF ((IABS(X(2)).LE.M) .AND. (X(3).GT.0) .AND. (X(3).LT.B)) RETURN
      WRITE (LUN,FMT=9010)
   20 CALL MPERR
      X(1) = 0
      RETURN

 9000 FORMAT (' *** INVALID SIGN IN CALL TO MPSIGB ***')
 9010 FORMAT (' *** X NOT VALID MP NUMBER IN CALL TO MPSIGB ***')
      END
      SUBROUTINE MPSIN(X,Y)
C RETURNS Y = SIN(X) FOR MP X AND Y,
C METHOD IS TO REDUCE X TO (-1, 1) AND USE MPSIN1, SO
C TIME IS O(M(T)T/LOG(T)).
C DIMENSION OF R IN CALLING PROGRAM MUST BE AT LEAST 5T+12
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      REAL RX,RY
      INTEGER I2,I3,IE,XS
C     ..
C     .. External Functions ..
      INTEGER MPCMPI
      EXTERNAL MPCMPI
C     ..
C     .. External Subroutines ..
      EXTERNAL MPABS,MPADDI,MPADDQ,MPART1,MPCHK,MPCMF,MPCMR,MPDIV,
     +         MPDIVI,MPERR,MPMUL,MPMULI,MPSIN1,MPSUB
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,IABS,SIN
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(5,12)
      I2 = 4*T + 11
      IF (X(1).NE.0) GO TO 20
   10 Y(1) = 0
      RETURN

   20 XS = X(1)
      IE = IABS(X(2))
      IF (IE.LE.2) CALL MPCMR(X,RX)
      CALL MPABS(X,R(I2))
C USE MPSIN1 IF ABS(X) .LE. 1
      IF (MPCMPI(R(I2),1).GT.0) GO TO 30
      CALL MPSIN1(R(I2),Y,1)
      GO TO 50
C FIND ABS(X) MODULO 2PI (IT WOULD SAVE TIME IF PI WERE
C PRECOMPUTED AND SAVED IN COMMON).
C FOR INCREASED ACCURACY COMPUTE PI/4 USING MPART1
   30 I3 = 2*T + 7
      CALL MPART1(5,R(I3))
      CALL MPMULI(R(I3),4,R(I3))
      CALL MPART1(239,Y)
      CALL MPSUB(R(I3),Y,Y)
      CALL MPDIV(R(I2),Y,R(I2))
      CALL MPDIVI(R(I2),8,R(I2))
      CALL MPCMF(R(I2),R(I2))
C SUBTRACT 1/2, SAVE SIGN AND TAKE ABS
      CALL MPADDQ(R(I2),-1,2,R(I2))
      XS = -XS*R(I2)
      IF (XS.EQ.0) GO TO 10
      R(I2) = 1
      CALL MPMULI(R(I2),4,R(I2))
C IF NOT LESS THAN 1, SUBTRACT FROM 2
      IF (R(I2+1).GT.0) CALL MPADDI(R(I2),-2,R(I2))
      IF (R(I2).EQ.0) GO TO 10
      R(I2) = 1
      CALL MPMULI(R(I2),2,R(I2))
C NOW REDUCED TO FIRST QUADRANT, IF LESS THAN PI/4 USE
C POWER SERIES, ELSE COMPUTE COS OF COMPLEMENT
      IF (R(I2+1).GT.0) GO TO 40
      CALL MPMUL(R(I2),Y,R(I2))
      CALL MPSIN1(R(I2),Y,1)
      GO TO 50

   40 CALL MPADDI(R(I2),-2,R(I2))
      CALL MPMUL(R(I2),Y,R(I2))
      CALL MPSIN1(R(I2),Y,0)
   50 Y(1) = XS
      IF (IE.GT.2) RETURN
C CHECK THAT ABSOLUTE ERROR LESS THAN 0.01 IF ABS(X) .LE. 100
C (IF ABS(X) IS LARGE THEN SINGLE-PRECISION SIN INACCURATE)
      IF (ABS(RX).GT.100.0) RETURN
      CALL MPCMR(Y,RY)
      IF (ABS(RY-SIN(RX)).LT.0.01) RETURN
      WRITE (LUN,FMT=9000)
C THE FOLLOWING MESSAGE MAY INDICATE THAT
C B**(T-1) IS TOO SMALL.
      CALL MPERR
      RETURN

 9000 FORMAT (' *** ERROR OCCURRED IN MPSIN, RESULT INCORRECT ***')
      END
      SUBROUTINE MPSINH(X,Y)
C RETURNS Y = SINH(X) FOR MP NUMBERS X AND Y, X NOT TOO LARGE.
C METHOD IS TO USE MPEXP OR MPEXP1, SPACE = 5T+12
C SAVE SIGN OF X AND CHECK FOR ZERO, SINH(0) = 0
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I2,I3,XS
C     ..
C     .. External Subroutines ..
      EXTERNAL MPABS,MPADDI,MPCHK,MPDIV,MPDIVI,MPEXP,MPEXP1,MPMUL,MPREC,
     +         MPSUB
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      XS = X(1)
      IF (XS.NE.0) GO TO 10
      Y(1) = 0
      RETURN
C CHECK LEGALITY OF B, T, M, LUN AND MXR
   10 CALL MPCHK(5,12)
      I3 = 4*T + 11
C WORK WITH ABS(X)
      CALL MPABS(X,R(I3))
      IF (R(I3+1).LE.0) GO TO 20
C HERE ABS(X) .GE. 1, IF TOO LARGE MPEXP GIVES ERROR MESSAGE
C INCREASE M TO AVOID OVERFLOW IF SINH(X) REPRESENTABLE
      M = M + 2
      CALL MPEXP(R(I3),R(I3))
      CALL MPREC(R(I3),Y)
      CALL MPSUB(R(I3),Y,Y)
C RESTORE M.  IF RESULT OVERFLOWS OR UNDERFLOWS, MPDIVI AT
C STATEMENT 30 WILL ACT ACCORDINGLY.
      M = M - 2
      GO TO 30
C HERE ABS(X) .LT. 1 SO USE MPEXP1 TO AVOID CANCELLATION
   20 I2 = I3 - (T+2)
      CALL MPEXP1(R(I3),R(I2))
      CALL MPADDI(R(I2),2,R(I3))
      CALL MPMUL(R(I3),R(I2),Y)
      CALL MPADDI(R(I2),1,R(I3))
      CALL MPDIV(Y,R(I3),Y)
C DIVIDE BY TWO AND RESTORE SIGN
   30 CALL MPDIVI(Y,2*XS,Y)
      RETURN

      END
      SUBROUTINE MPSIN1(X,Y,IS)
C COMPUTES Y = SIN(X) IF IS.NE.0, Y = COS(X) IF IS.EQ.0,
C USING TAYLOR SERIES.   ASSUMES ABS(X) .LE. 1.
C X AND Y ARE MP NUMBERS, IS AN INTEGER.
C TIME IS O(M(T)T/LOG(T)).   THIS COULD BE REDUCED TO
C O(SQRT(T)M(T)) AS IN MPEXP1, BUT NOT WORTHWHILE UNLESS
C T IS VERY LARGE.  ASYMPTOTICALLY FASTER METHODS ARE
C DESCRIBED IN THE REFERENCES GIVEN IN COMMENTS
C TO MPATAN AND MPPIGL.
C DIMENSION OF R IN CALLING PROGRAM MUST BE AT LEAST 3T+8
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER IS
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER B2,I,I2,I3,TS
C     ..
C     .. External Functions ..
      INTEGER MPCMPI
      EXTERNAL MPCMPI
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADD2,MPADDI,MPCHK,MPCIM,MPDIVI,MPERR,MPMUL,MPSTR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX0,MIN0
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(3,8)
      IF (X(1).NE.0) GO TO 20
C SIN(0) = 0, COS(0) = 1
   10 Y(1) = 0
      IF (IS.EQ.0) CALL MPCIM(1,Y)
      RETURN

   20 I2 = T + 5
      I3 = I2 + T + 2
      B2 = 2*MAX0(B,64)
      CALL MPMUL(X,X,R(I3))
      IF (MPCMPI(R(I3),1).LE.0) GO TO 30
      WRITE (LUN,FMT=9000)
      CALL MPERR
      GO TO 10

   30 IF (IS.EQ.0) CALL MPCIM(1,R(I2))
      IF (IS.NE.0) CALL MPSTR(X,R(I2))
      Y(1) = 0
      I = 1
      TS = T
      IF (IS.EQ.0) GO TO 40
      CALL MPSTR(R(I2),Y)
      I = 2
C POWER SERIES LOOP.  REDUCE T IF POSSIBLE
   40 T = R(I2+1) + TS + 2
      IF (T.LE.2) GO TO 70
      T = MIN0(T,TS)
C PUT R(I3) FIRST IN CASE ITS DIGITS ARE MAINLY ZERO
      CALL MPMUL(R(I3),R(I2),R(I2))
C IF I*(I+1) IS NOT REPRESENTABLE AS AN INTEGER, THE FOLLOWING
C DIVISION BY I*(I+1) HAS TO BE SPLIT UP.
      IF (I.GT.B2) GO TO 50
      CALL MPDIVI(R(I2),-I* (I+1),R(I2))
      GO TO 60

   50 CALL MPDIVI(R(I2),-I,R(I2))
      CALL MPDIVI(R(I2),I+1,R(I2))
   60 I = I + 2
      T = TS
      CALL MPADD2(R(I2),Y,Y,Y,0)
      IF (R(I2).NE.0) GO TO 40
   70 T = TS
      IF (IS.EQ.0) CALL MPADDI(Y,1,Y)
      RETURN

 9000 FORMAT (' *** ABS(X) .GT. 1 IN CALL TO MPSIN1 ***')
      END
      SUBROUTINE MPSQRT(X,Y)
C RETURNS Y = SQRT(X), USING SUBROUTINE MPROOT IF X .GT. 0.
C DIMENSION OF R IN CALLING PROGRAM MUST BE AT LEAST 4T+10
C (BUT Y(1) MAY BE R(3T+9)).  X AND Y ARE MP NUMBERS.
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I,I2,IY3
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK,MPERR,MPEXT,MPMUL,MPROOT
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(4,10)
C MPROOT NEEDS 4T+10 WORDS, BUT CAN OVERLAP SLIGHTLY.
      I2 = 3*T + 9
      IF (X(1)) 10,20,30
   10 WRITE (LUN,FMT=9000)
      CALL MPERR
   20 Y(1) = 0
      RETURN

   30 CALL MPROOT(X,-2,R(I2))
      I = R(I2+2)
      CALL MPMUL(X,R(I2),Y)
      IY3 = Y(3)
      CALL MPEXT(I,IY3,Y)
      RETURN

 9000 FORMAT (' *** X NEGATIVE IN CALL TO SUBROUTINE MPSQRT ***')
      END
      SUBROUTINE MPSTR(X,Y)
C SETS Y = X FOR MP X AND Y.
C SEE IF X AND Y HAVE THE SAME ADDRESS (THEY OFTEN DO)
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,T2
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      J = X(1)
      Y(1) = J + 1
      IF (J.EQ.X(1)) GO TO 10
C HERE X(1) AND Y(1) MUST HAVE THE SAME ADDRESS
      X(1) = J
      RETURN
C HERE X(1) AND Y(1) HAVE DIFFERENT ADDRESSES
   10 Y(1) = J
C NO NEED TO MOVE X(2), ... IF X(1) = 0
      IF (J.EQ.0) RETURN
      T2 = T + 2
      DO 20 I = 2,T2
          Y(I) = X(I)
   20 CONTINUE
      RETURN

      END
      SUBROUTINE MPSUB(X,Y,Z)
C SUBTRACTS Y FROM X, FORMING RESULT IN Z, FOR MP X, Y AND Z.
C FOUR GUARD DIGITS ARE USED, AND THEN R*-ROUNDING
C     .. Array Arguments ..
      INTEGER X(*),Y(*),Z(*)
C     ..
C     .. Local Arrays ..
      INTEGER Y1(1)
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADD2
C     ..
      Y1(1) = -Y(1)
      CALL MPADD2(X,Y,Z,Y1,0)
      RETURN

      END
      SUBROUTINE MPTAN(X,Y)
C SETS Y = TAN(X) FOR MP X AND Y
C USES SUBROUTINE MPSIN1 SO TIME IS O(M(T)T/LOG(T)).
C DIMENSION OF R IN CALLING PROGRAM MUST BE AT LEAST 6T+20
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER FL,I2,I3,TS,XS
C     ..
C     .. External Subroutines ..
      EXTERNAL MPABS,MPADDI,MPCHK,MPCLR,MPCMF,MPDIV,MPDIVI,MPMUL,MPMULI,
     +         MPOVFL,MPPI,MPREC,MPROOT,MPSIN1,MPSTR
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(6,20)
      TS = T
      I2 = 4*T + 15
      I3 = I2 + T + 3
      IF (X(1).NE.0) GO TO 20
   10 Y(1) = 0
      T = TS
      RETURN
C SAVE SIGN AND WORK WITH ABS(X)
   20 XS = X(1)
      CALL MPCLR(R(I2),T+1)
      CALL MPABS(X,R(I2))
C USE ONE GUARD DIGIT THROUGHOUT
      T = T + 1
      CALL MPPI(R(I3))
C COMPUTE ABS(X) MODULO PI
      CALL MPDIV(R(I2),R(I3),R(I2))
      CALL MPDIVI(R(I3),4,R(I3))
      CALL MPCMF(R(I2),R(I2))
      CALL MPMULI(R(I2),2,R(I2))
      IF (R(I2).EQ.0) GO TO 10
C NOW IN (0, 2), MAKE IT (-1, 1)
      IF (R(I2+1).GT.0) CALL MPADDI(R(I2),-2,R(I2))
      CALL MPMULI(R(I2),2,R(I2))
      IF (R(I2).EQ.0) GO TO 10
      XS = XS*R(I2)
      R(I2) = 1
C METHODS DEPEND ON WHETHER ABS(TAN(X)) .LT. 1 OR NOT
      FL = R(I2+1)
      IF (FL.LE.0) GO TO 30
      R(I2) = -1
      CALL MPADDI(R(I2),2,R(I2))
   30 CALL MPMUL(R(I2),R(I3),R(I2))
      CALL MPSIN1(R(I2),R(I2),1)
      CALL MPMUL(R(I2),R(I2),R(I3))
      R(I3) = -R(I3)
      CALL MPADDI(R(I3),1,R(I3))
      CALL MPROOT(R(I3),-2,R(I3))
      CALL MPMUL(R(I3),R(I2),R(I3))
      IF (FL.LE.0) GO TO 50
      IF (R(I3).NE.0) GO TO 40
      T = TS
C HERE X IS TOO CLOSE TO AN ODD MULTIPLE OF PI/2
C TREAT AS OVERFLOW THOUGH NOT QUITE THE SAME
      WRITE (LUN,FMT=9000)
      CALL MPOVFL(Y)
      RETURN

   40 CALL MPREC(R(I3),R(I3))
C RESTORE T AND MOVE RESULT NOW
   50 T = TS
      CALL MPSTR(R(I3),Y)
      Y(1) = XS*Y(1)
      RETURN

 9000 FORMAT (' *** TAN(X) TOO LARGE IN CALL TO MPTAN ***')
      END
      SUBROUTINE MPTANH(X,Y)
C RETURNS Y = TANH(X) FOR MP NUMBERS X AND Y,
C USING MPEXP OR MPEXP1, SPACE = 5T+12
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I2,XS
C     ..
C     .. External Functions ..
      INTEGER MPCOMP
      EXTERNAL MPCOMP
C     ..
C     .. External Subroutines ..
      EXTERNAL MPABS,MPADDI,MPCHK,MPCIM,MPCRM,MPDIV,MPEXP,MPEXP1,MPMULI
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ALOG,FLOAT
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      IF (X(1).NE.0) GO TO 10
C TANH(0) = 0
      Y(1) = 0
      RETURN
C CHECK LEGALITY OF B, T, M, LUN AND MXR
   10 CALL MPCHK(5,12)
      I2 = 4*T + 11
C SAVE SIGN AND WORK WITH ABS(X)
      XS = X(1)
      CALL MPABS(X,R(I2))
C SEE IF ABS(X) SO LARGE THAT RESULT IS +-1
      CALL MPCRM(0.5E0*FLOAT(T)*ALOG(FLOAT(B)),Y)
      IF (MPCOMP(R(I2),Y).LE.0) GO TO 20
C HERE ABS(X) IS VERY LARGE
      CALL MPCIM(XS,Y)
      RETURN
C HERE ABS(X) NOT SO LARGE
   20 CALL MPMULI(R(I2),2,R(I2))
      IF (R(I2+1).LE.0) GO TO 30
C HERE ABS(X) .GE. 1/2 SO USE MPEXP
      CALL MPEXP(R(I2),R(I2))
      CALL MPADDI(R(I2),-1,Y)
      CALL MPADDI(R(I2),1,R(I2))
      CALL MPDIV(Y,R(I2),Y)
      GO TO 40
C HERE ABS(X) .LT. 1/2, SO USE MPEXP1 TO AVOID CANCELLATION
   30 CALL MPEXP1(R(I2),R(I2))
      CALL MPADDI(R(I2),2,Y)
      CALL MPDIV(R(I2),Y,Y)
C RESTORE SIGN
   40 Y(1) = XS*Y(1)
      RETURN

      END
      SUBROUTINE MPUNFL(X)
C CALLED ON MULTIPLE-PRECISION UNDERFLOW, IE WHEN THE
C EXPONENT OF MP NUMBER X WOULD BE LESS THAN -M.
C SINCE M MAY HAVE BEEN OVERWRITTEN, CHECK B, T, M ETC.
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. External Subroutines ..
      EXTERNAL MPCHK
C     ..
      CALL MPCHK(1,4)
C THE UNDERFLOWING NUMBER IS SET TO ZERO
C AN ALTERNATIVE WOULD BE TO CALL MPMINR (X) AND RETURN,
C POSSIBLY UPDATING A COUNTER AND TERMINATING EXECUTION
C AFTER A PRESET NUMBER OF UNDERFLOWS.  ACTION COULD EASILY
C BE DETERMINED BY A FLAG IN LABELLED COMMON.
      X(1) = 0
      RETURN

      END
      SUBROUTINE MPUNPK(Y,X)
C RESTORES THE MP NUMBER X WHICH IS STORED IN COMPRESSED
C FORMAT IN THE INTEGER ARRAY Y.  FOR FURTHER DETAILS SEE
C SUBROUTINE MPPACK.
C CHECK FOR ZERO
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Array Arguments ..
      INTEGER X(*),Y(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I,IB,IS,J,K,K1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IABS
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      IF (Y(2).NE.0) GO TO 10
      X(1) = 0
      RETURN
C HERE Y IS NONZERO. GET SIGN THEN UNPACK DIGITS.
   10 IS = 1
      IF (Y(2).LT.0) IS = -1
      J = T/2
C DO LAST DIGIT SEPARATELY IF T ODD.
      IF ((2*J).LT.T) X(T+2) = Y(J+2)/B
C WORK BACKWARDS IN CASE X AND Y ARE THE SAME ARRAY.
      DO 20 IB = 1,J
          I = J - IB
          K = IABS(Y(I+2))
          K1 = K/B
          X(2*I+3) = K1
          X(2*I+4) = K - B*K1
   20 CONTINUE
C FINALLY MOVE EXPONENT AND SIGN TO X.
      X(2) = Y(1)
      X(1) = IS
      RETURN

      END
      SUBROUTINE MPZETA(N,X)
C RETURNS MP X = ZETA(N) FOR INTEGER N .GT. 1, WHERE
C ZETA(N) IS THE RIEMANN ZETA FUNCTION (THE SUM FROM
C I = 1 TO INFINITY OF I**(-N)).
C USES THE EULER-MACLAURIN SERIES UNLESS N = 2, 4, 6 OR 8.
C IN WORST CASE SPACE IS 8T+18+NL*((T+3)/2),
C WHERE NL IS THE NUMBER OF TERMS USED IN THE EULER-
C MACLAURIN SERIES, NL .LE. 1 + AL*T*LN(B), WHERE
C AL IS GIVEN IN A DATA STATEMENT BELOW.
C TIME IS O(T**3).
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      REAL AL,AL2,ALBT,FN,FNL
      INTEGER I,I2,I3,I4,IB,IP,IQ,J,N2,NL,NM,P,TS
C     ..
C     .. Local Arrays ..
      INTEGER ID(4)
C     ..
C     .. External Subroutines ..
      EXTERNAL MPADD,MPADD2,MPADDI,MPADDQ,MPBERN,MPCHK,MPCIM,MPCQM,
     +         MPDIVI,MPERR,MPMUL,MPMULQ,MPPI,MPPWR,MPUNPK
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ALOG,EXP,FLOAT,INT,MAX0,MOD
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
C     .. Data statements ..
C AL AND AL2 ARE EMPIRICALLY DETERMINED CONSTANTS, CHOSEN TO
C APPROXIMATELY MINIMIZE EXECUTION TIME.  IF SPACE IS CRITICAL,
C AL MAY BE REDUCED.
C ZETA(N) KNOWN IN TERMS OF BERNOULLI NUMBERS AND
C PI**N IF N IS EVEN.  ID DEFINES ZETA(2), ... , ZETA(8).
      DATA AL/0.1E0/,AL2/5.0E0/
      DATA ID(1),ID(2),ID(3),ID(4)/6,90,945,9450/
C     ..
C CHECK B, T ETC AND ENSURE ENOUGH SPACE FOR MPPI
      CALL MPCHK(3,8)
      X(1) = 0
      IF (N.GT.1) GO TO 10
      WRITE (LUN,FMT=9000) N
      GO TO 110
C HERE N .GT. 1.  SEE IF N = 2, 4, 6 OR 8.
   10 IF ((N.GT.8) .OR. (MOD(N,2).NE.0)) GO TO 20
C HERE ZETA(N) = (PI**N)/ID(N/2)
      CALL MPPI(X)
      CALL MPPWR(X,N,X)
      N2 = N/2
      CALL MPDIVI(X,ID(N2),X)
      GO TO 100
C SEE IF N IS VERY LARGE.  CAN RETURN WITH ZETA(N) = 1 TO
C REQUIRED PRECISION IF 2**N .GE. 2*B**(T-1)
   20 ALBT = ALOG(FLOAT(B))*FLOAT(T-1) + ALOG(2E0)
      FN = FLOAT(N)
      IF ((FN*ALOG(2E0)).GE.ALBT) GO TO 90
C HERE WE MAY USE EULER-MACLAURIN SERIES. FOR N = 3
C THE SERIES OF GOSPER (AS USED IN PROGRAM TEST) IS FASTER.
C SERIES FOR OTHER ODD N ARE GIVEN BY H. RIESEL IN
C BIT 13 (1973), 97-113.
C ESTIMATE NUMBER OF TERMS IN ASYMPTOTIC EXPANSION
      NL = 1 + INT(AL*ALBT)
C ESTIMATE NUMBER OF TERMS REQUIRED IN FINITE SUM.
C CONSTANTS ARE 2(1+LN(2*PI)) = 5.675 AND 1+LN(2*PI**2) = 3.982
      FNL = FLOAT(NL)
      NM = 2 + INT(EXP((ALBT+ (FN+2E0*FNL+0.5E0)*ALOG(FN+2E0*FNL)- ((FN-
     +     0.5E0)*ALOG(FN-1E0)+5.675E0*FNL+3.982E0))/ (FN+2E0*FNL+1E0)))
      P = (T+3)/2
      I2 = 6*T + 15
      I3 = I2 + T + 2
      I4 = I3 + T + 2
C SEE IF IT WOULD BE BETTER NOT TO USE ASYMPTOTIC EXPANSION
      IF (((FN-1E0)*ALOG(AL2*FNL+FLOAT(NM))).LT.ALBT) GO TO 30
C DONT USE ASYMPTOTIC EXPANSION, BUT RECOMPUTE NM
      NM = 2 + INT(EXP(ALBT/ (FN-1E0)))
      GO TO 60
C CHECK THAT SPACE IS SUFFICIENT
   30 CALL MPCHK(8,NL*P+18)
C COMPUTE REQUIRED BERNOULLI NUMBERS (IF ZETA(N) IS REQUIRED
C FOR SEVERAL N, IT WOULD SAVE TIME TO PRECOMPUTE THESE).
      CALL MPBERN(NL,P,R(I4))
      CALL MPCQM(N,2*NM,R(I2))
      CALL MPDIVI(R(I2),NM,R(I2))
C SUM EULER-MACLAURIN ASYMPTOTIC SERIES FIRST
      DO 40 I = 1,NL
          IP = I4 + (I-1)*P
          CALL MPUNPK(R(IP),R(I3))
          CALL MPMUL(R(I2),R(I3),R(I3))
          CALL MPADD(X,R(I3),X)
          CALL MPMULQ(R(I2),N+2*I-1,2*I+1,R(I2))
          CALL MPMULQ(R(I2),N+2*I,2*I+2,R(I2))
          CALL MPDIVI(R(I2),NM,R(I2))
          CALL MPDIVI(R(I2),NM,R(I2))
   40 CONTINUE
C ADD INTEGRAL APPROXIMATION AND MULTIPLY BY NM**(1-N)
      CALL MPADDQ(X,1,N-1,X)
      DO 50 I = 2,N
          CALL MPDIVI(X,NM,X)
   50 CONTINUE
C ADD FINITE SUM IN FORWARD DIRECTION SO CAN REDUCE T
C MORE EASILY THAN IF BACKWARD DIRECTION WERE USED.
   60 R(I2+1) = 0
      TS = T
      IQ = 1
      DO 80 I = 2,NM
C REDUCE T FOR I**(-N) COMPUTATION IF POSSIBLE
          T = MAX0(2,TS+R(I2+1))
          IB = MAX0(7*B*B,32767)/I
          CALL MPCIM(1,R(I2))
C DO SINGLE-PRECISION OPERATIONS WHERE POSSIBLE,
C MULTIPLE-PRECISION ONLY WHERE NECESSARY.
          DO 70 J = 1,N
              IQ = I*IQ
              IF ((IQ.LE.IB) .AND. (IQ.NE.B) .AND. (J.LT.N)) GO TO 70
              CALL MPDIVI(R(I2),IQ,R(I2))
              IQ = 1
   70     CONTINUE
C NOW R(I2) IS I**(-N).  HALVE LAST TERM IN FINITE SUM.
          IF (I.EQ.NM) CALL MPDIVI(R(I2),2,R(I2))
C RESTORE T FOR ADDITION
          T = TS
C LEAVE FINITE SUM LOOP IF MP UNDERFLOW OCCURRED
          IF (R(I2).EQ.0) GO TO 90
          CALL MPADD2(R(I2),X,X,X,0)
   80 CONTINUE
   90 CALL MPADDI(X,1,X)
C CHECK THAT 1 .LE. ZETA(N) .LT. 2
  100 IF ((X(1).EQ.1) .AND. (X(2).EQ.1) .AND. (X(3).EQ.1)) RETURN
      WRITE (LUN,FMT=9010)
  110 CALL MPERR
      RETURN

 9000 FORMAT (' *** N =',I10,' .LE. 1 IN CALL TO MPZETA ***')
 9010 FORMAT (' *** ERROR OCCURRED IN MPZETA, RESULT INCORRECT ***')
      END
      SUBROUTINE MP40D(N,X)
C PRINTING ROUTINE CALLED BY TEST PROGRAM, PRINTS MP X TO
C N DECIMAL PLACES, ASSUMING -10 .LT. X .LT. 100.
C SPACE = 3T + N + 14.
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I2
C     ..
C     .. External Subroutines ..
      EXTERNAL MP40E,MPCHK,MPOUT
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(3,N+14)
C DO NOTHING IF N NONPOSITIVE
      IF (N.LE.0) RETURN
      I2 = 3*T + 12
C CONVERT TO CHARACTER FORMAT AND CALL MP40E TO PRINT
      CALL MPOUT(X,R(I2),N+3,N)
      CALL MP40E(N+3,R(I2))
      RETURN

      END
      SUBROUTINE MP40E(N,X)
C WRITES X(1), ... , X(N) ON UNIT LUN, WHERE X IS AN
C INTEGER ARRAY OF DIMENSION AT LEAST N.  CALLED BY MP40D.
C FOLLOWING IS USUALLY FASTER THAN USING AN IMPLIED DO LOOP
C BUT DOES NOT WORK WITH RALPH COMPILER ON UNIVAC 1100
C COMPUTERS, WHEN IT SHOULD BE REPLACED BY IMPLIED DO LOOP.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      INTEGER X(N)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      WRITE (LUN,FMT=9000) X
      RETURN

 9000 FORMAT (8X,13A1,4 (1X,10A1),/ (11X,10A1,1X,10A1,1X,10A1,1X,10A1,
     +       1X,10A1))
      END
      SUBROUTINE MP40F(N,X)
C PRINTING ROUTINE CALLED BY TEST2 PROGRAM, PRINTS X TO
C N SIGNIFICANT FIGURES, N .GE. 2.
C DIM OF R IN COMMON AT LEAST 6T+N+17
C CHECK LEGALITY OF B, T, M, LUN AND MXR
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      INTEGER X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I2,J
C     ..
C     .. External Subroutines ..
      EXTERNAL MP40G,MPCHK,MPOUTE
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      CALL MPCHK(6,N+17)
C DO NOTHING IF N.LT.2
      IF (N.LT.2) RETURN
      I2 = 6*T + 15
C CONVERT TO PRINTABLE FORM AND CALL MP40G TO PRINT
      CALL MPOUTE(X,R(I2+1),J,N+2)
      R(I2) = J
      CALL MP40G(N+3,R(I2))
      RETURN

      END
      SUBROUTINE MP40G(N,X)
C WRITES X(1), ... , X(N) ON UNIT LUN, CALLED BY MP40F
C FOLLOWING IS USUALLY FASTER THAN USING AN IMPLIED DO LOOP
C BUT DOES NOT WORK WITH RALPH COMPILER ON UNIVAC 1100
C COMPUTERS, WHEN IT SHOULD BE REPLACED BY IMPLIED DO LOOP.
C     .. Parameters ..
      INTEGER RMAX
      PARAMETER (RMAX=15000)
C     ..
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      INTEGER X(N)
C     ..
C     .. Scalars in Common ..
      INTEGER B,LUN,M,MXR,T
C     ..
C     .. Arrays in Common ..
      INTEGER R(RMAX)
C     ..
C     .. Common blocks ..
      COMMON B,T,M,LUN,MXR,R
C     ..
      WRITE (LUN,FMT=9000) X
      RETURN
C
C
 9000 FORMAT (1X,I6,1X,13A1,4 (1X,10A1),
     +       / (11X,10A1,1X,10A1,1X,10A1,1X,10A1,1X,10A1))
      END
