#include "copyright.h"
#include "dprec.fh"

MODULE EREPMOD

  IMPLICIT NONE

CONTAINS

  SUBROUTINE VDW_REP_IJ( &
       &     QI,SI,Z0I,ZQI,    &
       &     QJ,SJ,Z0J,ZQJ,    &
       &     RIJ,              &
       &     EREP,dEdQI,dEdQJ)


    IMPLICIT NONE
    
    _REAL_,INTENT(IN) :: QI,SI,Z0I,ZQI
    _REAL_,INTENT(IN) :: QJ,SJ,Z0J,ZQJ
    _REAL_,INTENT(IN) :: RIJ
    _REAL_,INTENT(OUT) :: EREP,dEdQI,dEdQJ

    !==========================================
    
    _REAL_ :: ZERO
    _REAL_ :: ZI,ZJ
    _REAL_ :: dZIdQI,dZJdQJ
    _REAL_ :: S
    _REAL_ :: dSdZI,dSdZJ
    
    !==========================================
    
    
    ZERO  = 0.D0
    !C INTENT OUT VARIABLES INITIALIZED
    EREP  = ZERO
    dEdQI = ZERO
    dEdQJ = ZERO

    ZI = Z0I * EXP(-ZQI*QI)
    ZJ = Z0J * EXP(-ZQJ*QJ)
    
    dZIdQI = ZI * (-ZQI)
    dZJdQJ = ZJ * (-ZQJ)
    
    CALL SS_Overlap(ZI,ZJ,ZERO,ZERO,ZERO,ZERO,ZERO,RIJ,S)
    EREP = SI*SJ*S
    
    CALL SS_Overlap_DZA(ZI,ZJ,ZERO,ZERO,ZERO,ZERO,ZERO,RIJ,dSdZI)
    dEdQI = SI*SJ * dSdZI * dZIdQI

    CALL SS_Overlap_DZB(ZI,ZJ,ZERO,ZERO,ZERO,ZERO,ZERO,RIJ,dSdZJ)
    dEdQJ = SI*SJ * dSdZJ * dZJdQJ

  END SUBROUTINE VDW_REP_IJ



  SUBROUTINE VDW_REP_IJ_DRI(         &
       & QI,SI,Z0I,ZQI,              &
       & QJ,SJ,Z0J,ZQJ,              &
       & CXI,CYI,CZI,  CXJ,CYJ,CZJ,  &
       & GXI,GYI,GZI)
    
    
    IMPLICIT NONE

    _REAL_,INTENT(IN) :: QI,SI,Z0I,ZQI
    _REAL_,INTENT(IN) :: QJ,SJ,Z0J,ZQJ
    _REAL_,INTENT(IN) :: CXI,CYI,CZI,CXJ,CYJ,CZJ
    _REAL_,INTENT(OUT) :: GXI,GYI,GZI

    !==============================================

    _REAL_ :: ZERO
    _REAL_ :: DSDA(3)
    _REAL_ :: ZI,ZJ

    !==============================================


    ZERO  = 0.D0
!C INTENT OUT VARIABLES INITIALIZED
    GXI = ZERO
    GYI = ZERO
    GZI = ZERO
    
    ZI = Z0I * EXP(-ZQI*QI)
    ZJ = Z0J * EXP(-ZQJ*QJ)
    
    CALL SS_Overlap_DRA(ZI,ZJ,CXI,CYI,CZI,CXJ,CYJ,CZJ,DSDA)
    GXI = SI*SJ*DSDA(1)
    GYI = SI*SJ*DSDA(2)
    GZI = SI*SJ*DSDA(3)

  END SUBROUTINE VDW_REP_IJ_DRI


!C====================================================================
!C====================================================================
!C====================================================================



  SUBROUTINE SS_Overlap(za,zb,XI,YI,ZI,XJ,YJ,ZJ,SSS)
    
    IMPLICIT NONE

    _REAL_,INTENT(IN) :: za,zb,XI,YI,ZI,XJ,YJ,ZJ
    _REAL_,INTENT(OUT) :: SSS

    !==============================================

    _REAL_ :: PI
    _REAL_ :: r2,r,d,a,ra,ra2,si,sj,nij

    !==============================================


    PI     = 3.141592653589793238462643383279502884197D0

    r2   = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
    r    = SQRT(r2)
    
    d    = (za*za-zb*zb)**3
    
    
    IF ( r < 1.D-3 ) THEN
       SSS = (za*zb)**3 / ( 8*PI * (za+zb)**3 )
    ELSE
       IF ( ABS(d) < 1.D-11 ) THEN
          a  = (za+zb)/2.d0
          ra = r*a
          ra2= ra*ra
          SSS = EXP(-ra)*a*a*a*(3.d0+3.d0*ra+ra2)/(192.D0*PI)
       ELSE
          si  = DELTAIJ(za,zb,r)
          sj  = DELTAIJ(zb,za,r)
          nij = (za*zb)**3 / ( 8.D0*PI * d )
          SSS = nij * (si-sj)
       END IF
    END IF

  END SUBROUTINE SS_Overlap



  SUBROUTINE SS_Overlap_DRA(za,zb,XI,YI,ZI,XJ,YJ,ZJ,SSP)
    
    IMPLICIT NONE

    _REAL_,INTENT(IN) :: za,zb,XI,YI,ZI,XJ,YJ,ZJ
    _REAL_,INTENT(OUT) :: SSP(3)

    CALL SS_Overlap_DRB(zb,za,XJ,YJ,ZJ,XI,YI,ZI,SSP)
    
  END SUBROUTINE SS_Overlap_DRA

      
  SUBROUTINE SS_Overlap_DRB(za,zb,XI,YI,ZI,XJ,YJ,ZJ,SSP)
    
    IMPLICIT NONE

    _REAL_,INTENT(IN) :: za,zb,XI,YI,ZI,XJ,YJ,ZJ
    _REAL_,INTENT(OUT) :: SSP(3)

    !==============================================

    _REAL_ :: PI
    _REAL_ :: CRD(3),r2,r
    _REAL_ :: d,a,ra,tmp,siP,sjP,nij

    !==============================================


    PI     = 3.141592653589793238462643383279502884197D0

    CRD(1) = xi-xj
    CRD(2) = yi-yj
    CRD(3) = zi-zj
    
    r2   = CRD(1)*CRD(1) + CRD(2)*CRD(2) + CRD(3)*CRD(3)
    r    = SQRT(r2)
    
    d    = (za*za-zb*zb)**3
    
    IF ( r < 1.D-3 ) THEN
       SSP(1)=0.0D0
       SSP(2)=0.0D0
       SSP(3)=0.0D0
    ELSE
       IF ( ABS(d) < 1.D-11 ) THEN
          a  = (za+zb)/2.D0
          ra = r*a
          tmp = EXP(-ra)*(a**5)*(1+ra) / (192.D0*PI)
          SSP(1) = CRD(1) * tmp
          SSP(2) = CRD(2) * tmp
          SSP(3) = CRD(3) * tmp
       ELSE
          siP = DELTAIJ_PRIME(za,zb,r)
          sjP = DELTAIJ_PRIME(zb,za,r)
          nij = (za*zb)**3 / ( 8.D0*PI * d )
          tmp = nij * (sjP-siP) / r
          SSP(1) = CRD(1) * tmp
          SSP(2) = CRD(2) * tmp
          SSP(3) = CRD(3) * tmp
       END IF
    END IF

  END SUBROUTINE SS_Overlap_DRB







  SUBROUTINE SS_Overlap_dza(za,zb,XI,YI,ZI,XJ,YJ,ZJ,SSS)
    
    IMPLICIT NONE

    _REAL_,INTENT(IN) :: za,zb,XI,YI,ZI,XJ,YJ,ZJ
    _REAL_,INTENT(OUT) :: SSS

    !==============================================

    _REAL_ :: PI
    _REAL_ :: r2,r,d,dt,ddza
    _REAL_ :: denom,zazb,zazb2,zazb3
    _REAL_ :: a,ra,ra2,dadza,dradza,dra2dza,fact,t1,dt1dza
    _REAL_ :: si,sj,nij,dsidza,dsjdza,ddenomdza,dnijdza

    !==============================================
    
    PI     = 3.141592653589793238462643383279502884197D0
    
    r2   = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
    r    = SQRT(r2)
    
    d    = (za*za-zb*zb)**3
    dt   = za*za-zb*zb
    ddza= (3 * dt*dt) * 2 * za
    
    IF ( r < 1.D-3 ) THEN
       denom = ( 8*PI * (za+zb)**3 )
       zazb  = za*zb
       zazb2 = zazb  * zazb
       zazb3 = zazb2 * zazb
       SSS = 3.D0 * (zazb2 / denom) * zb
       SSS = SSS - (zazb3/denom**2) * ( 3*8*PI * (za+zb)**2 )
    ELSE
       IF ( ABS(d) < 1.D-9 ) THEN
          a  = (za+zb)/2.D0
          ra = r*a
          ra2= ra*ra
          dadza   = 0.5D0
          dradza  = r*dadza
          dra2dza = 2.D0*ra * dradza
          fact = 1.D0 / (192.D0*PI)
          t1 = (3.D0+3.D0*ra+ra2)
          dt1dza = 3.D0*dradza + dra2dza
          
          SSS = fact * (-dradza * EXP(-ra)) * (a**3) * t1
          SSS = SSS + fact * EXP(-ra) * (3.D0*a*a*dadza) * t1
          SSS = SSS + fact * EXP(-ra) * (a**3) * (dt1dza)
       ELSE
          si  = DELTAIJ(za,zb,r)
          sj  = DELTAIJ(zb,za,r)

          denom = ( 8.D0*PI * d )
          nij = (za*zb)**3 / denom
          
          dsidza = DELTAIJ_dza(za,zb,r)
          dsjdza = DELTAIJ_dzb(zb,za,r)
          ddenomdza = 8.D0*PI * ddza
          
            dnijdza = (3.D0*(za*zb)**2 / denom) * zb
            dnijdza = dnijdza - ((za*zb)**3 / denom**2) * ddenomdza
            
            SSS = dnijdza*(si-sj) + nij*(dsidza-dsjdza)
         END IF
      END IF
      
    END SUBROUTINE SS_Overlap_dza


    SUBROUTINE SS_Overlap_dzb(za,zb,XI,YI,ZI,XJ,YJ,ZJ,SSS)
      
      IMPLICIT NONE

      _REAL_,INTENT(IN) :: za,zb,XI,YI,ZI,XJ,YJ,ZJ
      _REAL_,INTENT(OUT) :: SSS
      
      CALL SS_Overlap_dza(zb,za,XJ,YJ,ZJ,XI,YI,ZI,SSS)
      
    END SUBROUTINE SS_Overlap_dzb


!C==============================================
!C  AUXILIARY FUNCTIONS

    FUNCTION DELTAIJ(a,b,r) RESULT(x)
      
      IMPLICIT NONE
      _REAL_,INTENT(IN) :: a,b,r
      _REAL_ :: x
      
      _REAL_ :: tmp

      tmp = b * (4.D0*a + r * (a*a-b*b)) / r
      x = EXP(-a*r) * tmp

    END FUNCTION DELTAIJ



    FUNCTION DELTAIJ_DR(a,b,r) RESULT(x)
      
      IMPLICIT NONE
      _REAL_,INTENT(IN) :: a,b,r
      _REAL_ :: x
      
      _REAL_ :: tmp,dtmp

      tmp  = b * (4.D0*a + r * (a*a-b*b)) / r
      dtmp = -(b * (4.D0*a + r * (a*a-b*b)))/r**2 &
           & + (b/r) * (a*a-b*b)
      x = (-a)*EXP(-a*r) * tmp + EXP(-a*r) * dtmp

    END FUNCTION DELTAIJ_DR



    FUNCTION DELTAIJ_DZA(a,b,r) RESULT(x)
      
      IMPLICIT NONE
      _REAL_,INTENT(IN) :: a,b,r
      _REAL_ :: x
      
      _REAL_ :: tmp

      tmp = (b/r) * (4.D0*a + r * (a*a-b*b))
      x = (-r*EXP(-a*r)) * tmp
      x = x + EXP(-a*r) * (b/r)*(4.D0 + r*(2*a))

    END FUNCTION DELTAIJ_DZA



    FUNCTION DELTAIJ_DZB(a,b,r) RESULT(x)
      
      IMPLICIT NONE
      _REAL_,INTENT(IN) :: a,b,r
      _REAL_ :: x
      
      x = EXP(-a*r)*((4.D0*a + r * (a*a-b*b)) / r)
      x = x + EXP(-a*r)* (b/r)*( r*(-2.D0*b) )

    END FUNCTION DELTAIJ_DZB

!C-----------------------------------------------

    FUNCTION DELTAIJ_PRIME(a,b,r) RESULT(x)
      
      IMPLICIT NONE
      _REAL_,INTENT(IN) :: a,b,r
      _REAL_ :: x

      _REAL_ :: S,c,d

      S = DELTAIJ(a,b,r)
      c = b * (a*a-b*b) / r
      d = S * (1.0D0/r + a)
      x = EXP(-a*r) * c - d

    END FUNCTION DELTAIJ_PRIME



    FUNCTION DELTAIJ_PRIME_DR(a,b,r) RESULT(x)
      
      IMPLICIT NONE
      _REAL_,INTENT(IN) :: a,b,r
      _REAL_ :: x

      _REAL_ :: S,dS,c,dc,d,dd

      S  = DELTAIJ(a,b,r)
      dS = DELTAIJ_DR(a,b,r)
      c  =   b * (a*a-b*b) / r
      dc = - b * (a*a-b*b) / r**2
      d  =  S * (1.0D0/r + a)
      dd = dS * (1.0D0/r + a) + S * (-1.D0/r**2)
      x = (-a)*EXP(-a*r) * c + EXP(-a*r)*dc - dd

    END FUNCTION DELTAIJ_PRIME_DR



    FUNCTION DELTAIJ_PRIME_DZA(a,b,r) RESULT(x)
      
      IMPLICIT NONE
      _REAL_,INTENT(IN) :: a,b,r
      _REAL_ :: x

      _REAL_ :: S,dSda,c,dCda,d,dDda

      S = DELTAIJ(a,b,r)
      dSda = DELTAIJ_DZA(a,b,r)
      c = b * (a*a-b*b) / r
      dCda = (b/r) * (2*a)
      d = S * (1.0D0/r + a)
      dDda = dSda * (1.0D0/r + a) + S
      x = (-r*EXP(-a*r)) * c + EXP(-a*r)*dCda - dDda

    END FUNCTION DELTAIJ_PRIME_DZA



    FUNCTION DELTAIJ_PRIME_DZB(a,b,r) RESULT(x)
      
      IMPLICIT NONE
      _REAL_,INTENT(IN) :: a,b,r
      _REAL_ :: x

      _REAL_ :: S,dSdb,c,dCdb,d,dDdb
      
      S = DELTAIJ(a,b,r)
      dSdb = DELTAIJ_DZB(a,b,r)
      c = b * (a*a-b*b) / r
      dCdb = ( (a*a-b*b) / r ) + (b * (-2*b) / r)
      d = S * (1.0D0/r + a)
      dDdb = dSdb * (1.0D0/r + a)
      x = EXP(-a*r)*dCdb - dDdb

    END FUNCTION DELTAIJ_PRIME_DZB

!C-----------------------------------------------


    SUBROUTINE SS_RadDer(za,zb,r,SSP)
      
      IMPLICIT NONE

      _REAL_,INTENT(IN) :: za,zb,r
      _REAL_,INTENT(OUT) :: SSP

      _REAL_ :: PI,r2,d,a,ra,tmp,siP,sjP,nij

      PI     = 3.141592653589793238462643383279502884197D0
      r2   = r*r
      d    = (za*za-zb*zb)**3

      IF ( r < 1.D-3 ) THEN
         SSP = 0.0D0
      ELSE
         IF ( ABS(d) < 1.D-11 ) THEN
            a  = (za+zb)/2.D0
            ra = r*a
            tmp = EXP(-ra)*(a**5)*(1+ra) / (192.D0*PI)
            SSP = r*tmp
         ELSE
            siP = DELTAIJ_PRIME(za,zb,r)
            sjP = DELTAIJ_PRIME(zb,za,r)
            nij = (za*zb)**3 / ( 8.D0*PI * d )
            tmp = nij * (sjP-siP)
            SSP = tmp
         END IF
      END IF

    END SUBROUTINE SS_RadDer




    SUBROUTINE SS_RadDer_DR(za,zb,r,SSP)
      
      IMPLICIT NONE

      _REAL_,INTENT(IN) :: za,zb,r
      _REAL_,INTENT(OUT) :: SSP

      _REAL_ :: PI,r2,d,a,ra,tmp,dtmp,siP,sjP,nij,dra

      PI     = 3.141592653589793238462643383279502884197D0
      r2   = r*r
      d    = (za*za-zb*zb)**3

      IF ( r < 1.D-3 ) THEN
         SSP = 0.0D0
      ELSE
         IF ( ABS(d) < 1.D-11 ) THEN
            a  = (za+zb)/2.D0
            ra = r*a
            dra = a
            tmp = EXP(-ra)*(a**5)*(1+ra) / (192.D0*PI)
            dtmp = (-dra)*EXP(-ra)*(a**5)*(1+ra) / (192.D0*PI) &
                 & + EXP(-ra)*(a**5)*(dra) / (192.D0*PI)
            SSP = r*dtmp + tmp
         ELSE
            siP = DELTAIJ_PRIME_DR(za,zb,r)
            sjP = DELTAIJ_PRIME_DR(zb,za,r)
            nij = (za*zb)**3 / ( 8.D0*PI * d )
            tmp = nij * (sjP-siP)
            SSP = tmp
         END IF
      END IF

    END SUBROUTINE SS_RadDer_DR



    SUBROUTINE SS_RadDer_DZA(za,zb,r,SSP)
      
      IMPLICIT NONE

      _REAL_,INTENT(IN) :: za,zb,r
      _REAL_,INTENT(OUT) :: SSP

      _REAL_ :: PI,r2,d,dt,dddza,denom
      _REAL_ :: a,dadza,ra,dradza,fact,a5,da5dza
      _REAL_ :: siP,sjP,nij,dSiPdza,dSjPdza,dnijdza

      PI     = 3.141592653589793238462643383279502884197D0
      r2   = r*r
      d    = (za*za-zb*zb)**3
      dt   = za*za-zb*zb
      dddza= (3*dt*dt)*(2*za)


      IF ( r < 1.D-3 ) THEN
         SSP = 0.0D0
      ELSE
         IF ( ABS(d) < 1.D-9 ) THEN
            a  = (za+zb)/2
            dadza = 0.5D0
            ra = r*a
            dradza = r*dadza
            fact = 1.D0 / (192.D0*PI)
            a5   = a**5
            da5dza = (5*a**4) * dadza

            SSP =       fact * (-dradza*EXP(-ra))*a5*(1+ra)
            SSP = SSP + fact * EXP(-ra) * da5dza * (1+ra)
            SSP = SSP + fact * EXP(-ra) * a5 * dradza
            SSP = SSP * r
         ELSE
            siP = DELTAIJ_PRIME(za,zb,r)
            sjP = DELTAIJ_PRIME(zb,za,r)
            nij = (za*zb)**3 / ( 8.D0*PI * d )

            dSiPdza = DELTAIJ_PRIME_DZA(za,zb,r)
            dSjPdza = DELTAIJ_PRIME_DZB(zb,za,r)
            denom = ( 8.D0*PI * d )
            dnijdza = (3.D0*(za*zb)**2 / denom)*zb &
                 & - ((za*zb)**3 / denom**2) * ( 8.D0*PI * dddza )

            SSP = dnijdza * (sjP-siP) + nij * (dSjPdza-dSiPdza)

         END IF
      END IF

    END SUBROUTINE SS_RadDer_DZA



    SUBROUTINE SS_RadDer_DZB(za,zb,r,SSP)
      
      IMPLICIT NONE

      _REAL_,INTENT(IN) :: za,zb,r
      _REAL_,INTENT(OUT) :: SSP

      CALL SS_RadDer_DZA(zb,za,r,SSP)

    END SUBROUTINE SS_RadDer_DZB


  END MODULE EREPMOD
