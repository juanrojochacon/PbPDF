      SUBROUTINE DSSZr(X,Q,A,RUV,RDV,RUB,RDB,RS,RC,RB,RG)
**************************************************************************
*
*       DSSZ NLO GLOBAL ABALYSIS OF NUCLEAR PARTON DISTRIBUTIONS                                                                     
*       D. de Florian, R. Sassot, M. Stratmann, P. Zurita                                          
*                             
*  REFERENCE:   arXiv:....
*
*  USAGE:
*                                               
*       The distributions are obtained by calling the subroutine                
*                                                                            
*               DSSZ(X,Q,A,RUV,RDV,RUB,RDB,RS,RC,RB,RG)                       
*                                                                            
*  INPUT:                                                     
*            X= x_bjorken                                                    
*            Q= factorization scale                                          
*            A= atomic number  
*
*  OUTPUT: the subroutine returns the NUCLEAR RATIOS
*                R_f(x,Q,A) = f_A(x,Q)/f_p(x,Q)
*          where f_A corresponds the pdf of a parton of flavour f in a proton 
*          of a nucleus A, and f_p is the corresponding pdf in the free proton. 
*
*                    RUV :    U VALENCE DISTRIBUTION                     
*                    RDV :    D VALENCE DISTRIBUTION                     
*                    RUB :    UBAR DISTRIBUTION                          
*                    RDB :    DBAR DISTRIBUTION                          
*                    RS  :    STR(=STRBAR) DISTRIBUTION
*                    RC  :    CHARM DISTRIBUTION 
*                    RB  :    BOTTOM DISTRIBUTION
*                    RG  :    GLUON DISTRIBUTION           
*              
*          The "nuclear ratio" in the neutron can be obtained by isospin symmetry,
*          i.e., RUV_proton=RDV_neutron, etc
*                                                                            
*  INIT:  the routine DSSZINI(ISET) has to be called once to initialize a new grid
*                   ISET : 0 CENTRAL FIT
*                   ISET : 1 TO 25 (AND -1 TO -25) TO CALL HESSIAN UNCERTAINTY SETS
*                   
*  RANGE OF VALIDITY OF THE INTERPOLATION:                              
*                           10^(-4) < X < 1.0                     
*                           1 < Q**2 < 10^4                    
*                           A = 9, 56, 197, 208      
*
*  IN CASE OF PROBLEMS, DOUBTS, ETC, PLEASE REPORT TO :                 
*           deflo@df.uba.ar                                                     
*           sassot@df.uba.ar     
*           marco@bnl.gov
*           pia@df.uba.ar
*                                                                            
**************************************************************************
      IMPLICIT NONE
C...
      INTEGER NX, NQ, NARG, NNA
      PARAMETER (NX=49, NQ=23, NARG=3, NNA=4)
C...
      DOUBLE PRECISION XUVF(NX,NQ,NNA), XDVF(NX,NQ,NNA), 
     1     XUBF(NX,NQ,NNA), XDBF(NX,NQ,NNA), XSF(NX,NQ,NNA),
     2     XCF(NX,NQ,NNA), XBF(NX,NQ,NNA), XGF(NX,NQ,NNA)
      DOUBLE PRECISION RUV,RDV,RUB,RDB,RS,RC,RB,RG
      DOUBLE PRECISION XT(NARG)
      INTEGER NA(NARG),A
      DOUBLE PRECISION ARRF(NX+NQ+NNA) 
      DOUBLE PRECISION X,Q
      DOUBLE PRECISION DSSZFINT
C...
      COMMON/ DSSZGRID / XUVF, XDVF, XUBF, XDBF, XSF, XCF, XBF, XGF,
     1                   ARRF, NA
C---INTERPOLATAION AND OUTPUT
       xt(1)=x
       xt(2)=dlog(Q*Q)
       xt(3)=A
C....
       RUV = DSSZFINT(NARG,XT,NA,ARRF,XUVF) 
       RDV = DSSZFINT(NARG,XT,NA,ARRF,XDVF) 
       RUB = DSSZFINT(NARG,XT,NA,ARRF,XUBF) 
       RDB = DSSZFINT(NARG,XT,NA,ARRF,XDBF) 
       RS = DSSZFINT(NARG,XT,NA,ARRF,XSF)
       RC = DSSZFINT(NARG,XT,NA,ARRF,XCF)
       RB = DSSZFINT(NARG,XT,NA,ARRF,XBF)
       RG = DSSZFINT(NARG,XT,NA,ARRF,XGF)       
C...
       RETURN
       END
C
C-----------------------------
      SUBROUTINE DSSZINI(ISET)
C-----------------------------
      IMPLICIT NONE
C...
      INTEGER NX, NQ, NARG, NNA, ISET
      PARAMETER (NX=49, NQ=23, NARG=3, NNA=4)
C...
      DOUBLE PRECISION XUVF(NX,NQ,NNA), XDVF(NX,NQ,NNA), 
     1     XUBF(NX,NQ,NNA), XDBF(NX,NQ,NNA), XSF(NX,NQ,NNA),
     2     XCF(NX,NQ,NNA), XBF(NX,NQ,NNA), XGF(NX,NQ,NNA)
      DOUBLE PRECISION QS(NQ), XB(NX), XA(NNA), XT(NARG)
      DOUBLE PRECISION ARRF(NX+NQ+NNA) 
C...
      INTEGER NA(NARG)
      INTEGER I, J, K
C...
      CHARACTER *11 FNAME1(2)
      CHARACTER *2  FNAME2, FNAME2A
      CHARACTER *5  FNAME3
      CHARACTER *20 FNAME, FILE
      CHARACTER *1 DUMMY0
C...
      DATA FNAME1 / 'dssz-plus-','dssz-minus-'/
C...
      COMMON/ DSSZGRID / XUVF, XDVF, XUBF, XDBF, XSF, XCF, XBF, XGF,
     1                   ARRF, NA
C...
C...TABLE OF (X,Q2) SUPPORT POINTS 
C...
      DATA QS / 1.D0, 1.5D0, 2.D0, 3.D0, 4.5D0, 6.D0, 8.D0, 10.D0,
     1     15.D0, 20.D0, 30.D0, 45.d0, 60.D0, 80.D0, 100.D0,
     2     200.D0, 300.D0, 450.D0, 700.D0, 1000.D0,
     3     2500.D0, 6000.D0, 10000.D0/
      DATA XB / 1.0E-4, 1.4E-4, 2.0E-4, 3.0E-4, 4.5E-4, 6.7E-4,
     3     1.0E-3, 1.4E-3, 2.0E-3, 3.0E-3, 4.5E-3, 6.7E-3,
     4     0.01D0, 0.014D0, 0.02D0, 0.03D0, 0.045D0, 0.06D0, 0.08D0,
     5     0.1D0, 0.125D0, 0.15D0, 0.175D0, 0.2D0, 0.225D0, 0.25D0,
     6     0.3D0, 0.325D0, 0.35D0, 0.375D0, 0.4D0, 0.425D0, 0.45D0,
     7     0.475D0, 0.5D0, 0.525D0, 0.55D0, 0.0575D0, 0.6D0,
     8     0.625D0, 0.65D0, 0.7D0, 0.75D0, 0.8D0, 0.85D0,
     9     0.90D0, 0.925D0, 0.95D0, 1.D0/
C...
C...SUPPORTED NUCLEI: Be, Fe, Au, Pb
C...
      DATA XA/ 9.D0, 56.D0, 197.D0, 208.D0/
C...
      FNAME3='.grid'
      DUMMY0='0'
C...
      IF(ISET.EQ.0) THEN
         FILE='dssz-central.grid'
      ELSE IF((ISET.LT.0).AND.(ISET.GE.-9)) THEN
         WRITE (FNAME2,1000) ABS(ISET)
 1000    FORMAT(I1)
         call strcat(dummy0,fname2,fname2a)
         call strcat(fname1(2),fname2a,fname)
         call strcat(fname,fname3,file)
      ELSE IF(ISET.LT.-9) THEN
         write (fname2,1001) abs(iset)
 1001    format(I2)
         call strcat(fname1(2),fname2,fname)
         call strcat(fname,fname3,file)
      ELSE IF((iset.gt.0).and.(iset.le.9)) then
         write (fname2,1002) abs(iset)
 1002    format(I1)
         call strcat(dummy0,fname2,fname2a)
         call strcat(fname1(1),fname2a,fname)
         call strcat(fname,fname3,file)
      ELSE IF(iset.gt.9) then
         write (fname2,1003) abs(iset)
 1003    format(I2)
         call strcat(fname1(1),fname2,fname)
         call strcat(fname,fname3,file)
      ENDIF
C...
      OPEN(UNIT=71,FILE=file,STATUS='OLD')
C...
      do i=1,nx-1
         do j=1,nq
            do k=1,nna
               read(71,90) xuvf(i,j,k),xdvf(i,j,k),xubf(i,j,k),
     1              xdbf(i,j,k),xsf(i,j,k),
     2              xcf(i,j,k),xbf(i,j,k),xgf(i,j,k)
            enddo
         enddo
      enddo
      CLOSE(71)
 90   FORMAT (8(1PE15.7))
C... SET VALUE AT X=1 FOR ALL Q**2 AND ALL A VALUES
      DO J=1,NQ
         DO K=1,NNA
            XUVF(NX,J,K)  = 2.D0
            XDVF(NX,J,K)  = 2.D0
            XUBF(NX,J,K) = 2.D0
            XDBF(NX,J,K) = 2.D0
            XSF(NX,J,K) = 2.D0
            XCF(NX,J,K) = 2.D0
            XBF(NX,J,K) = 2.D0
            XGF(NX,J,K) = 2.D0
         ENDDO
      ENDDO
C...
      NA(1)=NX
      NA(2)=NQ
      NA(3)=NNA
C...      
      do i=1,NX
         ARRF(I) = XB(I)     
      enddo
      do j=1,NQ
         ARRF(J+NX) = DLOG(QS(J))     
      enddo
      do K=1,NNA
         ARRF(k+NX+NQ) = XA(k)     
      enddo
C...
      RETURN
      END
C...      


c----------------------------------------------------------------------------------------
      subroutine strcat(str1,str2,str)
c concatenates str1 and str2 into str. Ignores trailing blanks of str1,str2
      character *(*) str1,str2,str
      l1=istrl(str1)
      l2=istrl(str2)
      l =len(str)
      if(l.lt.l1+l2) then
          write(*,*) 'error: l1+l2>l in strcat'
          write(*,*) 'l1=',l1,' str1=',str1
          write(*,*) 'l2=',l2,' str2=',str2
          write(*,*) 'l=',l
          stop
      endif
      if(l1.ne.0) str(1:l1)=str1(1:l1)
      if(l2.ne.0) str(l1+1:l1+l2)=str2(1:l2)
      if(l1+l2+1.le.l) str(l1+l2+1:l)= ' '
      end
c
      function istrl(string)
c returns the position of the last non-blank character in string
      character * (*) string
      i = len(string)
      dowhile(i.gt.0.and.string(i:i).eq.' ')
         i=i-1
      enddo
      istrl = i
      end
c----------------------------------------------------------------------------------------
      double precision FUNCTION DSSZFINT(NARG,ARG,NENT,ENT,TABLE)
*********************************************************************
*                                                                   *
*   THE INTERPOLATION ROUTINE (CERN LIBRARY ROUTINE E104)           *
*                                                                   *
*********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION ARG(5),NENT(5),ENT(10),TABLE(10)
      DIMENSION D(5),NCOMB(5),IENT(5)
      KD=1
      M=1
      JA=1
         DO 5 I=1,NARG
      NCOMB(I)=1
      JB=JA-1+NENT(I)
         DO 2 J=JA,JB
      IF (ARG(I).LE.ENT(J)) GO TO 3
    2 CONTINUE
      J=JB
    3 IF (J.NE.JA) GO TO 4
      J=J+1
    4 JR=J-1
      D(I)=(ENT(J)-ARG(I))/(ENT(J)-ENT(JR))
      IENT(I)=J-JA
      KD=KD+IENT(I)*M
      M=M*NENT(I)
    5 JA=JB+1
      DSSZFINT=0.
   10 FAC=1.
      IADR=KD
      IFADR=1
         DO 15 I=1,NARG
      IF (NCOMB(I).EQ.0) GO TO 12
      FAC=FAC*(1.-D(I))
      GO TO 15
   12 FAC=FAC*D(I)
      IADR=IADR-IFADR
   15 IFADR=IFADR*NENT(I)
      DSSZFINT=DSSZFINT+FAC*TABLE(IADR)
      IL=NARG
   40 IF (NCOMB(IL).EQ.0) GO TO 80
      NCOMB(IL)=0
      IF (IL.EQ.NARG) GO TO 10
      IL=IL+1
         DO 50  K=IL,NARG
   50 NCOMB(K)=1
      GO TO 10
   80 IL=IL-1
      IF(IL.NE.0) GO TO 40
      RETURN
      END



