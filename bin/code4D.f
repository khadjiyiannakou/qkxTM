      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=4,LX2=4,LX3=4,LX4=4)
      PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
      PARAMETER(NCOL=3,NCOL2=NCOL*NCOL)
      PARAMETER(LMAX=LX4/2+1,MAXDTLS=LMAX-1)
c     
      PARAMETER(NUMBIN=10)
      PARAMETER(PARBS=0.30,PARBDS=0.12)
C
      COMMON/ARRAYS/U11(NCOL,NCOL,LX1,LX2,LX3,LX4,4)
      COMPLEX*16 U11
      REAL*8 rndnum
C
      INQUIRE(IOLENGTH=l_rec) U11(1,1,1,1,1,1,1)
C
      OPEN(40,FILE="conf.0000.dat",access="direct",form='unformatted', 
     &CONVERT='BIG_ENDIAN', RECL=l_rec)
C
      CALL STARTF
      ISEED=137152
      CALL RINI(ISEED)
C
      CALL SETUP
C
      BETAG=4.0625d0
C
      IHEAT=1000
      ITOT=NITER
C
      WRITE(6,90)
90    FORMAT(' ****************************************************')
      WRITE(6,91)
91    FORMAT(' *')
      WRITE(6,92) LX1,LX2,LX3,LX4,BETAG,NCOL
92    FORMAT('  LX = ',4I6,'    BETA =',F8.4,'   NCOLOR = ',I4)
      WRITE(6,91)
      WRITE(6,93) IHEAT,ITOT
93    FORMAT('  NUM HEATS =',I6,'   NUM ITER =',I6)
      WRITE(6,91)
      WRITE(6,90)
C
      ACTLL=0.0d0
      IBIN=IHEAT/10
      DO 2 ITER=1,IHEAT
         IBIAS=2
         IF(ITER.LE.200) IBIAS=1
         IF(ITER.EQ.(ITER/5 )*5 ) IBIAS=1
 45      CALL UPDATG(IBIAS,ACTL,BETAG)
         CALL RENORM
         ACTLL=ACTLL+ACTL
         IF(ITER.EQ.(ITER/IBIN)*IBIN)THEN
         NBIN=ITER/IBIN
         ACTLL=ACTLL/IBIN
         WRITE(6,75) NBIN,1.0d0-ACTLL
75       FORMAT('  NBIN =',I4,' ACTLL HEATS =',F16.10)
         ACTLL=0.0d0
         ENDIF
2     CONTINUE
c        goto99
C
      CALL PLAQUETTE
99    ISEED = INT(rndnum()*float(259200))

C     
      ICALL=1
      DO IJ=1, NCOL
         DO IK=1, NCOL
            DO L1=1, LX1
               DO L2=1, LX2
                  DO L3=1, LX3
                     DO L4=1, LX4
                        DO MU=1, 4
C     
       WRITE(40,REC=ICALL) U11(IJ,IK,L1,L2,L3,L4,MU)
C     
                           ICALL=ICALL+1
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

C
      STOP
      END
C *******************************************************************
C *************         SUBROUTINE RUNGRB           *****************
C *******************************************************************
      SUBROUTINE PLAQUETTE
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=4,LX2=4,LX3=4,LX4=4)
      PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
      PARAMETER(NCOL=3,NCOL2=NCOL*NCOL)

      COMMON/NEXT/IUP(LSIZE,4),IDN(LSIZE,4)
      COMMON/ARRAYS/U11(NCOL2,LSIZE,4)

      COMPLEX*16 A11(NCOL2)
      COMPLEX*16 B11(NCOL2)
      COMPLEX*16 C11(NCOL2)
      COMPLEX*16 DUM11(NCOL2)
      COMPLEX*16 AVG, AKT
      COMPLEX*16 U11

      AVG=(0.0D0,0.0D0)

      DO NU=1, 4
         DO MU=NU+1,4
            DO I1=1,LSIZE
               I2=IUP(I1,NU)

               DO IJ=1, NCOL2
                  A11(IJ)=U11(IJ,I1,NU)
                  B11(IJ)=U11(IJ,I2,MU)
               ENDDO

               CALL VMX(1,A11,B11,C11,1)

               I2=IUP(I1,MU)

               DO IJ=1, NCOL2
                  A11(IJ)=U11(IJ,I2,NU)
               ENDDO

               CALL HERM(1,A11,DUM11,1)

               CALL VMX(1,C11,A11,B11,1)

               DO IJ=1, NCOL2
                  A11(IJ)=U11(IJ,I1,MU)
               ENDDO

               CALL HERM(1,A11,DUM11,1)
               CALL TRVMX(1,B11,A11,AKT,1)
     
               AVG=AVG+AKT/(LSIZE*6.0d0*NCOL)
            ENDDO
         ENDDO
      ENDDO
      WRITE(6,90) DREAL(AVG)
      WRITE(6,91) DIMAG(AVG)
 90   FORMAT('[Info][Plaquette]','              Real[Plaquette] = ', 
     &E15.8)
 91   FORMAT('[Info][Plaquette]','              Imag[Plaquette] = ', 
     &E15.8)
c      WRITE(6,*) "Average plaquette:", DREAL(AVG), DIMAG(AVG)
      RETURN
      END
C *******************************************************************
C *************         SUBROUTINE RUNGRB           *****************
C *******************************************************************
      SUBROUTINE RUNGRB(B11,C11)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NCOL=3,NCOL2=NCOL*NCOL)
C     
      COMPLEX*16 B11(NCOL2),C11(NCOL2)
      COMMON/SUBB/II1,JJ1,II2,JJ2,II3,JJ3,II4,JJ4
C     
      DO 10 LDU=1,NCOL-1
         DO 15 LDL=LDU+1,NCOL
            II1=LDU
            JJ1=LDU
            II2=LDL
            JJ2=LDL
            II3=LDL
            JJ3=LDU
            II4=LDU
            JJ4=LDL
            CALL SUBGRB(LDU,LDL,B11,C11)
 15      CONTINUE
 10   CONTINUE
C     
      RETURN
      END
C*********************************************************************
C*******************      SUBROUTINE SUBGRB      *********************
C*********************************************************************
      SUBROUTINE SUBGRB(LDU,LDL,B11,C11)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NCOL=3,NCOL2=NCOL*NCOL)
C     
      COMMON/SUBB/II1,JJ1,II2,JJ2,II3,JJ3,II4,JJ4
C     
      COMPLEX*16 B11(NCOL,NCOL),C11(NCOL,NCOL)
      COMPLEX*16 F11,F12,A11(NCOL,NCOL),S11(NCOL,NCOL)
      COMPLEX*16 T11(NCOL,NCOL)
C     
      F11=(B11(II1,JJ1)+CONJG(B11(II2,JJ2)))*0.5
      F12=(B11(II3,JJ3)-CONJG(B11(II4,JJ4)))*0.5
      UMAG=SQRT(F11*CONJG(F11)+F12*CONJG(F12))
      UMAG=1./UMAG
      F11=F11*UMAG
      F12=F12*UMAG
C     
      A11(II1,JJ1)=CONJG(F11)
      A11(II4,JJ4)=CONJG(F12)
      A11(II3,JJ3)=-F12
      A11(II2,JJ2)=F11
C     
      DO 3 IJ1=1,NCOL
         S11(IJ1,LDU)=C11(IJ1,LDU)*(A11(LDU,LDU)-1.0)
     &        +C11(IJ1,LDL)*A11(LDL,LDU)
         S11(IJ1,LDL)=C11(IJ1,LDU)*A11(LDU,LDL)
     &        +C11(IJ1,LDL)*(A11(LDL,LDL)-1.0)
         T11(IJ1,LDU)=B11(IJ1,LDU)*(A11(LDU,LDU)-1.0)
     &        +B11(IJ1,LDL)*A11(LDL,LDU)
         T11(IJ1,LDL)=B11(IJ1,LDU)*A11(LDU,LDL)
     &        +B11(IJ1,LDL)*(A11(LDL,LDL)-1.0)
 3    CONTINUE
      DO 4 IJ1=1,NCOL
         C11(IJ1,LDU)=C11(IJ1,LDU)+S11(IJ1,LDU)
         C11(IJ1,LDL)=C11(IJ1,LDL)+S11(IJ1,LDL)
         B11(IJ1,LDU)=B11(IJ1,LDU)+T11(IJ1,LDU)
         B11(IJ1,LDL)=B11(IJ1,LDL)+T11(IJ1,LDL)
 4    CONTINUE
C     
      RETURN
      END
C**************************************************************
C**************************************************************
C       HERE WE CONSTRUCT A FROZEN GAUGE CONFIGURATION
C**************************************************************
C**************************************************************
      SUBROUTINE STARTF
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=4,LX2=4,LX3=4,LX4=4)
      PARAMETER(LSIZE=LX1*LX2*LX3*LX4)
      PARAMETER(NCOL=3,NCOL2=NCOL*NCOL)
C
      COMMON/ARRAYS/U11(NCOL2,LSIZE,4)
      COMPLEX*16 U11,A11(NCOL2)
C
      DO 4 IJ=1,NCOL2
         A11(IJ)=(0.0,0.0)
4     CONTINUE
      DO 6 N1=1,NCOL
         NC=N1+NCOL*(N1-1)
         A11(NC)=(1.0,0.0)
6     CONTINUE
C
      DO 1 MU=1,4
         DO 1 NN=1,LSIZE
            DO 2 IJ=1,NCOL2
               U11(IJ,NN,MU)=A11(IJ)
2           CONTINUE
1     CONTINUE
C
      RETURN
      END
C**********************************************************
C  THIS ROUTINE REIMPOSES THE UNITARITY CONSTRAINTS ON
C  OUR SU(3) MATRICES
C**********************************************************
      SUBROUTINE RENORM
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=4,LX2=4,LX3=4,LX4=4)
      PARAMETER(LSIZE=LX1*LX2*LX3*LX4)
      PARAMETER(NCOL=3,NCOL2=NCOL*NCOL)
C
      COMMON/ARRAYS/U11(NCOL,NCOL,LSIZE,4)
      COMPLEX*16 ADUM(NCOL,NCOL),U11,CSUM
C
      DO 100 MU=1,4
         DO 100 NN=1,LSIZE
C
            DO 1 J=1,NCOL
               DO 1 I=1,NCOL
                  ADUM(I,J)=U11(I,J,NN,MU)
1           CONTINUE
C
            DO 20 N2=1,NCOL
C
               DO 10 N3=1,N2-1
C
                  CSUM=(0.0,0.0)
                  DO 5 N1=1,NCOL
                     CSUM=CSUM+ADUM(N1,N2)*CONJG(ADUM(N1,N3))
 5                CONTINUE
                  DO 6 N1=1,NCOL
                     ADUM(N1,N2)=ADUM(N1,N2)-CSUM*ADUM(N1,N3)
 6                CONTINUE
C     
 10            CONTINUE
C
               SUM=0.0
               DO 7 N1=1,NCOL
                  SUM=SUM+ADUM(N1,N2)*CONJG(ADUM(N1,N2))
 7             CONTINUE
                  ANORM=1.0/SQRT(SUM)
                  DO 8 N1=1,NCOL
                     ADUM(N1,N2)=ADUM(N1,N2)*ANORM
8                 CONTINUE
C
20          CONTINUE
C
            DO 2 J=1,NCOL
               DO 2 I=1,NCOL
                  U11(I,J,NN,MU)=ADUM(I,J)
2           CONTINUE
C
100   CONTINUE
      RETURN
      END
C****************************************************************
C****************************************************************
C     HERE WE SET UP  LINK POINTERS FOR FULL LATTICE            *
C****************************************************************
C****************************************************************
      SUBROUTINE SETUP
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=4,LX2=4,LX3=4,LX4=4)
      PARAMETER(LSIZE=LX1*LX2*LX3*LX4)
C
      COMMON/NEXT/IUP(LSIZE,4),IDN(LSIZE,4)
C
      LX12=LX1*LX2
      LX123=LX12*LX3
      LX1234=LX123*LX4
C
      NN=0
      DO 2 L4=1,LX4
         DO 2 L3=1,LX3
            DO 2 L2=1,LX2
               DO 2 L1=1,LX1
                  NN=NN+1
C     
               NU1=NN+1
               IF((L1+1).GT.LX1) NU1=NU1-LX1
               IUP(NN,1)=NU1
               ND1=NN-1
               IF((L1-1).LT.1) ND1=ND1+LX1
               IDN(NN,1)=ND1
C     
               NU2=NN+LX1
               IF((L2+1).GT.LX2) NU2=NU2-LX12
               IUP(NN,2)=NU2
               ND2=NN-LX1
               IF((L2-1).LT.1) ND2=ND2+LX12
               IDN(NN,2)=ND2
C     
               NU3=NN+LX12
               IF((L3+1).GT.LX3) NU3=NU3-LX123
               IUP(NN,3)=NU3
               ND3=NN-LX12
               IF((L3-1).LT.1) ND3=ND3+LX123
               IDN(NN,3)=ND3
C
               NU4=NN+LX123
               IF((L4+1).GT.LX4) NU4=NU4-LX1234
               IUP(NN,4)=NU4
               ND4=NN-LX123
               IF((L4-1).LT.1) ND4=ND4+LX1234
               IDN(NN,4)=ND4
C     
2     CONTINUE
C
      RETURN
      END
C****************************************************************
C****************************************************************
C*                  HERE WE UPDATE THE LINKS                    *
C****************************************************************
C****************************************************************
      SUBROUTINE UPDATG(IBIAS,ACTL,BETAG)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=4,LX2=4,LX3=4,LX4=4)
      PARAMETER(LSIZE=LX1*LX2*LX3*LX4)
      PARAMETER(NCOL=3,NCOL2=NCOL*NCOL)
C
      COMMON/ARRAYS/U11(NCOL2,LSIZE,4)
      COMMON/NEXT/IUP(LSIZE,4),IDN(LSIZE,4)
      COMMON/HB1/BBETAG,IIBIAS
      DIMENSION UINT11(NCOL2),DUM11(NCOL2)
     &,A11(NCOL2),B11(NCOL2),C11(NCOL2)
C
      COMPLEX*16 U11,A11,B11,C11,UINT11,DUM11
C
      BBETAG=BETAG
      IIBIAS=IBIAS
      SUM=0.
      DO 1 NN=1,LSIZE
         M1=NN
C
         DO 10 MU=1,4
C
            DO 11 IJ=1,NCOL2
                UINT11(IJ)=(0.,0.)
11          CONTINUE
            M5=IUP(M1,MU)
C
            DO 20 NU=1,4
               IF(MU.EQ.NU) GO TO 20
C
               DO 15 IJ=1,NCOL2
                  A11(IJ)=U11(IJ,M1,NU)
15             CONTINUE
               M3=IUP(M1,NU)
               DO 17 IJ=1,NCOL2
                  B11(IJ)=U11(IJ,M3,MU)
17             CONTINUE
               CALL VMX(1,A11,B11,C11,1)
               DO 19 IJ=1,NCOL2
                  A11(IJ)=U11(IJ,M5,NU)
19             CONTINUE
               CALL HERM(1,A11,DUM11,1)
               CALL VMX(1,C11,A11,B11,1)
               DO 23 IJ=1,NCOL2
                  UINT11(IJ)=UINT11(IJ)+B11(IJ)
23             CONTINUE
               M3=IDN(M1,NU)
               DO 25 IJ=1,NCOL2
                  B11(IJ)=U11(IJ,M3,MU)
25             CONTINUE
               DO 27 IJ=1,NCOL2
                  A11(IJ)=U11(IJ,M3,NU)
27             CONTINUE
               CALL HERM(1,A11,DUM11,1)
               CALL VMX(1,A11,B11,C11,1)
               M3=IDN(M5,NU)
               DO 31 IJ=1,NCOL2
                  A11(IJ)=U11(IJ,M3,NU)
31             CONTINUE
               CALL VMX(1,C11,A11,B11,1)
               DO 32 IJ=1,NCOL2
                  UINT11(IJ)=UINT11(IJ)+B11(IJ)
32             CONTINUE
20          CONTINUE
C
            CALL HERM(1,UINT11,DUM11,1)
C
            DO 36 IJ=1,NCOL2
               C11(IJ)=U11(IJ,M1,MU)
36          CONTINUE
            CALL VMX(1,UINT11,C11,B11,1)
C
            CALL RUNGRH(B11,C11)
C
            DO 40 IJ=1,NCOL2
               U11(IJ,M1,MU)=C11(IJ)
40          CONTINUE
            DO 42 I=1,NCOL
               NC=I+NCOL*(I-1)
               SUM=SUM+REALPART(B11(NC))
42          CONTINUE
C
10       CONTINUE
1     CONTINUE
C
      APQ=1.0-SUM/(24.d0*NCOL*LSIZE)
      ACTL=APQ
C
      RETURN
      END
C*********************************************************************
C     SUBROUTINE RUNGRH
C*********************************************************************
      SUBROUTINE RUNGRH(B11,C11)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NCOL=3,NCOL2=NCOL*NCOL)
C
      COMPLEX*16 B11(NCOL2),C11(NCOL2)
      COMMON/SUB/II1,JJ1,II2,JJ2,II3,JJ3,II4,JJ4
C
      DO 10 LDU=1,NCOL-1
         DO 15 LDL=LDU+1,NCOL
            II1=LDU
            JJ1=LDU
            II2=LDL
            JJ2=LDL
            II3=LDL
            JJ3=LDU
            II4=LDU
            JJ4=LDL
C
            CALL SUBGRH(LDU,LDL,B11,C11)
15       CONTINUE
10    CONTINUE
C
      RETURN
      END
C***********************************************************************
C     SUBROUTINE SUBGRH
C***********************************************************************
      SUBROUTINE SUBGRH(LDU,LDL,B11,C11)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NCOL=3,NCOL2=NCOL*NCOL)
C     
      COMMON/HB1/BETAG,IBIAS
      COMMON/HB2/UMAG,A0,A1,A2,A3
      COMMON/SUB/II1,JJ1,II2,JJ2,II3,JJ3,II4,JJ4
C     
      COMPLEX*16 B11(NCOL,NCOL),C11(NCOL,NCOL)
      COMPLEX*16 F11,F12,A11(NCOL,NCOL),S11(NCOL,NCOL)
      COMPLEX*16 E11,E12,G11,G12,T11(NCOL,NCOL)
C     
      F11=(B11(II1,JJ1)+CONJG(B11(II2,JJ2)))*0.5
      F12=(B11(II3,JJ3)-CONJG(B11(II4,JJ4)))*0.5
      UMAG=SQRT(F11*CONJG(F11)+F12*CONJG(F12))
      UMAG=1./UMAG
      F11=F11*UMAG
      F12=F12*UMAG
      IF(IBIAS.EQ.1)THEN
         CALL HEATB
         E11=CMPLX(A0,A1)
         E12=CMPLX(A2,A3)
         G11=CONJG(F11)*E11+CONJG(F12)*E12
         G12=-F12*E11+F11*E12
C     
         A11(II1,JJ1)=G11
         A11(II4,JJ4)=-CONJG(G12)
         A11(II3,JJ3)=G12
         A11(II2,JJ2)=CONJG(G11)
      ENDIF
      IF(IBIAS.EQ.2)THEN
         G11= F11*F11-CONJG(F12)*F12
         G12= F12*F11+CONJG(F11)*F12
C     
         A11(II1,JJ1)=CONJG(G11)
         A11(II4,JJ4)=CONJG(G12)
         A11(II3,JJ3)=-G12
         A11(II2,JJ2)=G11
      ENDIF
      IF(IBIAS.EQ.3)THEN
C     
         A11(II1,JJ1)=CONJG(F11)
         A11(II4,JJ4)=CONJG(F12)
         A11(II3,JJ3)=-F12
         A11(II2,JJ2)=F11
      ENDIF
C     
      DO 3 IJ1=1,NCOL
         S11(IJ1,LDU)=C11(IJ1,LDU)*(A11(LDU,LDU)-1.0)
     &        +C11(IJ1,LDL)*A11(LDL,LDU)     
         S11(IJ1,LDL)=C11(IJ1,LDU)*A11(LDU,LDL)
     &        +C11(IJ1,LDL)*(A11(LDL,LDL)-1.0)     
         T11(IJ1,LDU)=B11(IJ1,LDU)*(A11(LDU,LDU)-1.0)
     &        +B11(IJ1,LDL)*A11(LDL,LDU)     
         T11(IJ1,LDL)=B11(IJ1,LDU)*A11(LDU,LDL)
     &        +B11(IJ1,LDL)*(A11(LDL,LDL)-1.0)     
 3    CONTINUE
      DO 4 IJ1=1,NCOL
         C11(IJ1,LDU)=C11(IJ1,LDU)+S11(IJ1,LDU)
         C11(IJ1,LDL)=C11(IJ1,LDL)+S11(IJ1,LDL)
         B11(IJ1,LDU)=B11(IJ1,LDU)+T11(IJ1,LDU)
         B11(IJ1,LDL)=B11(IJ1,LDL)+T11(IJ1,LDL)
 4    CONTINUE
C     
      RETURN
      END
C*****************************************************
C SU2 SUBGROUP HEATBATH
C*****************************************************
      SUBROUTINE HEATB
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=4,LX2=4,LX3=4,LX4=4)
      PARAMETER(LSIZE=LX1*LX2*LX3*LX4)
      PARAMETER(NCOL=3,NCOL2=NCOL*NCOL)
C     
      COMMON/HB1/BETAG,IBIAS
      COMMON/HB2/UMAG,A0,A1,A2,A3
      REAL*8 RNDNUM
C     
      PI=4.0*DATAN(1.0d0)
      BTG=BETAG*2.0d0/NCOL
      TEMP=1.d0/BTG
      DD=RNDNUM()
C     
      UMAG=UMAG*TEMP
      A1=-DLOG(RNDNUM())*UMAG
      A2=-DLOG(RNDNUM())*UMAG
      A3=DCOS(2.0d0*PI*RNDNUM())
      A3=A3*A3
      A1=A1*A3
      A2=A2+A1
      A0=1.0d0-A2
      A3=RNDNUM()
      A3=A3*A3-1+A2*0.5d0
      IF(A3.GT.0.0d0)THEN
 45      X11=-DLOG(RNDNUM())*UMAG
         X22=-DLOG(RNDNUM())*UMAG
         CCC=DCOS(2.0d0*PI*RNDNUM())
         CCC=CCC*CCC
         AA=X11*CCC
         DBB=X22+AA
         XX=RNDNUM()
         XX=XX*XX-1.0d0+DBB*0.5d0
         IF(XX.GT.0.0d0)GOTO45
         A0=1.0d0-DBB
      ENDIF
C     
 54   CONTINUE
      X1=2.0*RNDNUM()-1.0d0
      X2=2.0*RNDNUM()-1.0d0
      X3=2.0*RNDNUM()-1.0d0
      CC=X1*X1+X2*X2+X3*X3-1.0d0
      IF(CC.LE.0.0)THEN
         A1=X1
         A2=X2
         A3=X3
      ELSE
         GOTO54
      ENDIF
      RAD=(1.d0-A0*A0)
      RAD=(A1*A1+A2*A2+A3*A3)/RAD
      RAD=1.d0/DSQRT(RAD)
      A1=A1*RAD
      A2=A2*RAD
      A3=A3*RAD
C     
      RETURN
      END
C***********************************************************************
C                 VECTOR MATRIX MULTIPLY ... 5*5 COMPLEX               *
C***********************************************************************
      SUBROUTINE VMX(NNN1,A,B,C,NNN2)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NCOL=3,NCOL2=NCOL*NCOL)
C
      COMPLEX*16 A(NCOL,NCOL),B(NCOL,NCOL),C(NCOL,NCOL),CSUM
C
      DO 1 J=1,NCOL
         DO 1 I=1,NCOL
            CSUM=(0.0,0.0)
            DO 2 K=1,NCOL
               CSUM=CSUM+A(I,K)*B(K,J)
2           CONTINUE
            C(I,J)=CSUM
1     CONTINUE
C
      RETURN
      END
C***********************************************************************
C***********************************************************************
C                   TRACE PRODUCT .. Ncol*Ncol COMPLEX                 *
C***********************************************************************
C***********************************************************************
      SUBROUTINE TRVMX(NNN1,A,B,CC,NNN2)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NCOL=3,NCOL2=NCOL*NCOL)
C
      COMPLEX*16 A(NCOL,NCOL),B(NCOL,NCOL),CC
C
      CC=(0.0,0.0)
      DO 1 I=1,NCOL
         DO 1 K=1,NCOL
            CC=CC+A(I,K)*B(K,I)
1     CONTINUE
C
      RETURN
      END
C***********************************************************************
C***********************************************************************
C                           HERMITIAN CONJUGATE                        *
C***********************************************************************
C***********************************************************************
      SUBROUTINE HERM(NNN1,A11,DUM11,NNN2)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NCOL=3,NCOL2=NCOL*NCOL)
C
      COMPLEX*16 A11(NCOL,NCOL),DUM11(NCOL,NCOL)
C
      DO 1  I=1,NCOL
         DO 1  J=1,NCOL
            DUM11(I,J)=CONJG(A11(J,I))
1     CONTINUE
      DO 4  I=1,NCOL
         DO 4  J=1,NCOL
            A11(I,J)=DUM11(I,J)
4     CONTINUE
C
      RETURN
      END
