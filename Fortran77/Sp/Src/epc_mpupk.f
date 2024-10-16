C--**--CH3973--524--C:U--13:7:2000
C--**--CH184--524--Fix--27:5:1999
      SUBROUTINE MPUPK(SOURCE,DEST,LDEST,LFIELD)

C     .. Parameters ..
      INTEGER NCHAR
      PARAMETER (NCHAR=4)
C     ..
C     .. Scalar Arguments ..
      INTEGER LDEST,LFIELD
C     ..
C     .. Array Arguments ..
      INTEGER DEST(LDEST),SOURCE(*)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,NWORDS
      CHARACTER*4 CH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C     ..
      IF (LDEST.LE.0) RETURN
c maximum of ldest characters contained in nwords words
      NWORDS = LDEST/NCHAR
      LFIELD = 0
      DO 20 I = 1,NWORDS
          WRITE (CH,FMT='(a4)') SOURCE(I)
          DO 10 J = 1,NCHAR
              IF (CH(J:J).EQ.'$') RETURN
              LFIELD = LFIELD + 1
              READ (CH(J:J),FMT='(a1)') DEST(LFIELD)
   10     CONTINUE
   20 CONTINUE

      I = MOD(LDEST,NCHAR)
      IF (I.EQ.0) RETURN
c and the remaining part word
      WRITE (CH,FMT='(a4)') SOURCE(NWORDS+1)
      DO 30 J = 1,I
          IF (CH(J:J).EQ.'$') RETURN
          LFIELD = LFIELD + 1
          READ (CH(J:J),FMT='(a1)') DEST(LFIELD)
   30 CONTINUE
      END
