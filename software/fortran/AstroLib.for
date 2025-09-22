*     -------------------------------------------------------------------------
*
*                                AstroLib.cs
*
* this library contains various astrodynamic routines.
*
*                            companion code for
*               fundamentals of astrodynamics and applications
*                                    2013
*                              by david vallado
*
*               email dvallado@comspoc.com, davallado@gmail.com
*
*    current :
*              17 jun 20  david vallado
*                           misc fixes
*    changes :
*              15 jan 19  david vallado
*                           combine with astiod etc
*              11 jan 18  david vallado
*                           misc cleanup
*              30 sep 15  david vallado
*                           fix jd, jdfrac
*              19 mar 14  david vallado
*                           original baseline
*       ----------------------------------------------------------------      



* -----------------------------------------------------------------------------
*
*                           FUNCTION GSTIME
*
*  this function finds the Greenwich sidereal time (iau-82). 
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    JD          - Julian Date                    days from 4713 BC
*
*  OutPuts       :
*    GSTIME      - Greenwich SIDEREAL Time        0 to 2Pi rad
*
*  Locals        :
*    Temp        - Temporary variable for reals   rad
*    TUT1        - Julian Centuries from the
*                  Jan 1, 2000 12 h epoch (UT1)
*
*  Coupling      :
*
*  References    :
*    vallado       2007, 193, Eq 3-43
* -----------------------------------------------------------------------------

      REAL*8 FUNCTION GSTIME ( JD )
        IMPLICIT NONE
        REAL*8 JD
* ----------------------------  Locals  -------------------------------
        REAL*8 Temp, TUT1

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------

        TUT1= ( JD - 2451545.0D0 ) / 36525.0D0
        Temp= - 6.2D-6*TUT1*TUT1*TUT1
     &        + 0.093104D0*TUT1*TUT1
     &        + (876600.0D0*3600.0D0 + 8640184.812866D0)*TUT1
     &        + 67310.54841D0
        Temp= DMOD( Temp*Deg2Rad/240.0D0,TwoPi ) ! 360/86400 = 1/240, to deg, to rad

        ! ------------------------ Check quadrants --------------------
        IF ( Temp .lt. 0.0D0 ) THEN
            Temp= Temp + TwoPi
          ENDIF

        GSTIME= Temp

      RETURN
      END

* -----------------------------------------------------------------------------
*
*                           FUNCTION GSTIM0
*
*  this function finds the Greenwich sidereal time at the beginning of a year.
*    This formula is derived from the Astronomical Almanac and is good only
*    0 hr UT1, Jan 1 of a year.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Year        - Year                           1998, 1999, etc.
*
*  OutPuts       :
*    GSTIM0      - Greenwich sidereal Time        0 to 2Pi rad
*
*  Locals        :
*    JD          - Julian Date                    days from 4713 BC
*    Temp        - Temporary variable for Reals   rad
*    TUT1        - Julian Centuries from the
*                  Jan 1, 2000 12 h epoch (UT1)
*
*  Coupling      :
*
*  References    :
*    Vallado       2007, 195, Eq 3-46
*
* -----------------------------------------------------------------------------

      REAL*8 FUNCTION GSTIM0 ( Year )
        IMPLICIT NONE
        INTEGER Year
* ----------------------------  Locals  -------------------------------
        REAL*8 JD, Temp, TUT1

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        JD  = 367.0D0 * Year - ( INT((7*(Year+INT(10/12)))*0.25D0) )+
     &                           ( INT(275/9) ) + 1721014.5D0
        TUT1= ( JD - 2451545.0D0 ) / 36525.0D0

        Temp= - 6.2D-6*TUT1*TUT1*TUT1
     &        + 0.093104D0*TUT1*TUT1
     &        + (876600.0D0*3600.0D0 + 8640184.812866D0)*TUT1
     &        + 67310.54841D0

        ! ------------------------ Check quadrants --------------------
        Temp= DMOD( Temp,TwoPi )
        IF ( Temp .lt. 0.0D0 ) THEN
            Temp= Temp + TwoPi
          ENDIF

        GSTIM0= Temp

      RETURN
      END

* -----------------------------------------------------------------------------
*
*                           SUBROUTINE LSTIME
*
*  this subroutine finds the Local sidereal time at a given location.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Lon         - SITE longitude (WEST -)        -2Pi to 2Pi rad
*    JD          - Julian Date                    days from 4713 BC
*
*  OutPuts       :
*    LST         - Local sidereal Time            0.0D0 to 2Pi rad
*    GST         - Greenwich sidereal Time        0.0D0 to 2Pi rad
*
*  Locals        :
*    None.
*
*  Coupling      :
*    GSTIME        Finds the Greenwich sidereal Time
*
*  References    :
*    Vallado       2007, 194, alg 15, ex 3-5
*
* -----------------------------------------------------------------------------

      SUBROUTINE LSTIME      ( Lon, JD, LST, GST )
        IMPLICIT NONE
        REAL*8 Lon, JD, LST, GST

* ----------------------------  Locals  -------------------------------
        EXTERNAL GSTIME
        REAL*8 GSTIME

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        GST = GSTIME( JD )
        LST = Lon + GST

        ! ----------------------- Check quadrants ---------------------
        LST = DMOD( LST,TwoPi )
        IF ( LST .lt. 0.0D0 ) THEN
            LST= LST + TwoPi
          ENDIF

      RETURN
      END
      
      
* ----------------------------------------------------------------------------
*
*                           SUBROUTINE INITREDUC
*
*  This procedure initializes the nutation matricies needed for reduction
*    calculations. The routine needs the filename of the files as input.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    FileN1      - Name for nutation file         nutation.dat
*    FileN2      - Name for planetary nutation    nut85.dat
*
*  Outputs       :
*    IAr80       - Integers for FK5 1980
*    RAr80       - Reals for FK5 1980
*    IAr00       - Integers for IAU 2000
*    RAr00       - Reals for IAU 2000
*    IAr96       - Integers for IAU 1996
*    RAr96       - Reals for IAU 1996
*    pIAr96      - Integers for planetary IAU 1996
*    pRAr96      - Reals for planetary IAU 1996
*
*  Locals        :
*    convrt      - conversion factor to degrees
*    i           - Index
*
*  Coupling      :
*    None        -
*
*  References    :
*
* -------------------------------------------------------------------------------

      SUBROUTINE InitReduc   ( FileN1 )
c      SUBROUTINE InitReduc   ( FileN1, FileN2, IAr00, RAr00, IAr80,
c     &                         RAr80, IAr96, RAr96, pIAr96, pRAr96 )
       IMPLICIT NONE
       CHARACTER*64 FileN1
c       INTEGER IAr00(5,106), IAr96(5,263), pIAr96(10,112)
c       REAL*8  RAr00(4,106), RAr96(4,263), pRAr96(4,112)


* ----------------------------  Locals  -------------------------------
       INTEGER i
       REAL*8 Convrt

       INCLUDE 'astrolib.cmn'

       ! --------------------  Implementation   ----------------------
c       Convrt= 0.0001D0/3600.0D0 ! 0.0001" to deg

       ! Open this for future use in the program
       OPEN( UNIT=44, FILE='timerec.rec',FORM='UNFORMATTED',
     &       ACCESS='DIRECT',RECL=80,STATUS='UNKNOWN')

       OPEN( 42, FILE = FileN1, STATUS = 'OLD' )
       READ( 42,*)
       READ( 42,*)
       ! ---------------- Read in IAU 2000 Theroy coefficients --------
       DO i = 1, 106
         Read( 42,*)
c         Read( 42,*) IAr00(1,i),IAr00(2,i),IAr00(3,i),IAr00(4,i),
c     &               IAr00(5,i),RAr00(1,i),RAr00(2,i),RAr00(3,i),
c     &               RAr00(4,i),RAr00(5,i),RAr00(6,i)
c         RAr00(1,i)= RAr00(1,i) * Convrt
c         RAr00(2,i)= RAr00(2,i) * Convrt
c         RAr00(3,i)= RAr00(3,i) * Convrt
c         RAr00(4,i)= RAr00(4,i) * Convrt
c         RAr00(5,i)= RAr00(5,i) * Convrt
c         RAr00(6,i)= RAr00(6,i) * Convrt
       ENDDO

       ! ---------------- Read in 1980 IAU Theory of Nutation ---------
       READ( 42,*)
       READ( 42,*)
       Convrt= 0.0001D0/3600.0D0 ! 0.0001" to deg
       DO i = 1, 106
         Read( 42,*) IAr80(1,i),IAr80(2,i),IAr80(3,i),IAr80(4,i),
     &               IAr80(5,i),RAr80(1,i),RAr80(2,i),RAr80(3,i),
     &               RAr80(4,i)
         RAr80(1,i)= RAr80(1,i) * Convrt
         RAr80(2,i)= RAr80(2,i) * Convrt
         RAr80(3,i)= RAr80(3,i) * Convrt
         RAr80(4,i)= RAr80(4,i) * Convrt
       ENDDO

       ! ---------------- Read in 1996 IAU Theory of Nutation ---------
c       READ( 42,*)
c       READ( 42,*)
c       Convrt= 0.000001D0/3600.0D0 ! 0.000001" to deg
c       DO i = 1, 263
c         Read( 42,*) IAr96(1,i),IAr96(2,i),IAr96(3,i),IAr96(4,i),
c     &               IAr96(5,i),RAr96(1,i),RAr96(2,i),RAr96(3,i),
c     &               RAr96(4,i),RAr96(5,i),RAr96(6,i)
c         RAr96(1,i)= RAr96(1,i) * Convrt
c         RAr96(2,i)= RAr96(2,i) * Convrt
c         RAr96(3,i)= RAr96(3,i) * Convrt
c         RAr96(4,i)= RAr96(4,i) * Convrt
c         RAr96(5,i)= RAr96(5,i) * Convrt
c         RAr96(6,i)= RAr96(6,i) * Convrt
c       ENDDO
c
c       CLOSE( 42 )

       ! ------------- Read in 1996 planetary correction terms --------
c       OPEN( 42, FILE = FileN2, STATUS = 'OLD' )
c       READ( 42,*)
c       READ( 42,*)
c       Convrt= 0.000001D0/3600.0D0 ! 0.000001" to deg
c       DO i= 1, 112
c         Read( 42,*) pIAr96(1,i),pIAr96(2,i),pIAr96(3,i),pIAr96(4,i),
c     &               pIAr96(5,i),pIAr96(6,i),pIAr96(7,i),pIAr96(8,i),
c     &               pIAr96(9,i),pIAr96(10,i),pRAr96(1,i),pRAr96(2,i),
c     &               pRAr96(3,i),pRAr96(4,i)
c         pRAr96(1,i)= pRAr96(1,i) * Convrt
c         pRAr96(2,i)= pRAr96(2,i) * Convrt
c         pRAr96(3,i)= pRAr96(3,i) * Convrt
c         pRAr96(4,i)= pRAr96(4,i) * Convrt
c        ENDDO
c
       CLOSE( 42 )

      RETURN
      END ! Subroutine InitReduc


* ----------------------------------------------------------------------------
*
*                           SUBROUTINE FK4
*
*  this subroutine converts vectors between the B1950 and J2000 epochs.  Be
*    aware that this process is not exact. There are different secular rates
*    for each system, and there are differences in the central location. The
*    matrices are multiplied directly for speed.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    rJ2000      - Initial J2000 GCRF Position vec ER, km, etc
*    vJ2000      - Initial J2000 GCRF Velocity vec
*    Direction   - Which set of vars to output    FROM  TOO
*
*  Outputs       :
*    rFK4        - Position vector FK4
*    vFK4        - Velocity vector FK4
*
*  Locals        :
*    r11,r12,r13 - Components of rotation matrix
*    r21,r22,r23 - Components of rotation matrix
*    r31,r32,r33 - Components of rotation matrix
*
*  Coupling      :
*
*
*  References    :
*    Vallado       2007, 240
*
* ----------------------------------------------------------------------------

      SUBROUTINE FK4         ( rJ2000, vJ2000, Direction, rFK4, vFK4 )
        IMPLICIT NONE
        CHARACTER*4 Direction
        REAL*8 rJ2000(3), vJ2000(3), rFK4(3), vFK4(3)

* ----------------------------  Locals  -------------------------------
        REAL*8 r11, r12, r13, r21, r22, r23, r31, r32, r33

        r11=  0.9999256794956877D0
        r12= -0.0111814832204662D0 
        r13= -0.0048590038153592D0 

        r21=  0.0111814832391717D0
        r22=  0.9999374848933135D0 
        r23= -0.0000271625947142D0 

        r31=  0.0048590037723143D0 
        r32= -0.0000271702937440D0 
        r33=  0.9999881946043742D0 

        IF ( Direction .eq. 'TOO ' ) THEN
            rFK4(1) = r11*rJ2000(1) + r21*rJ2000(2) + r31*rJ2000(3)
            rFK4(2) = r12*rJ2000(1) + r22*rJ2000(2) + r32*rJ2000(3)
            rFK4(3) = r13*rJ2000(1) + r23*rJ2000(2) + r33*rJ2000(3)
            vFK4(1) = r11*vJ2000(1) + r21*vJ2000(2) + r31*vJ2000(3)
            vFK4(2) = r12*vJ2000(1) + r22*vJ2000(2) + r32*vJ2000(3)
            vFK4(3) = r13*vJ2000(1) + r23*vJ2000(2) + r33*vJ2000(3)
          ELSE
            rJ2000(1) = r11*rFK4(1) + r12*rFK4(2) + r13*rFK4(3)
            rJ2000(2) = r21*rFK4(1) + r22*rFK4(2) + r23*rFK4(3)
            rJ2000(3) = r31*rFK4(1) + r32*rFK4(2) + r33*rFK4(3)
            vJ2000(1) = r11*vFK4(1) + r12*vFK4(2) + r13*vFK4(3)
            vJ2000(2) = r21*vFK4(1) + r22*vFK4(2) + r23*vFK4(3)
            vJ2000(3) = r31*vFK4(1) + r32*vFK4(2) + r33*vFK4(3)
          ENDIF

      RETURN
      END   ! SUBROUTINE FK4

* ----------------------------------------------------------------------------
*
*                           SUBROUTINE PRECESSION
*
*  this function calulates the transformation matrix for precession.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    TTT         - Julian Centuries of TT         centuries
*
*  Outputs       :
*    Prec        - Precession transformation (eci-mod)
*
*  Locals        :
*    TTT2        - TTT squared
*    TTT3        - TTT cubed
*    Zeta        - PRECESSION ANGLE               rad
*    z           - PRECESSION ANGLE               rad
*    Theta       - PRECESSION ANGLE               rad
*
*  Coupling      :
*    none
*
*  References    :
*    Vallado       2007, 228
*
* ----------------------------------------------------------------------------

      SUBROUTINE PRECESSION  ( TTT, Prec )
        IMPLICIT NONE
        REAL*8 TTT, Prec(3,3)

* ----------------------------  Locals  -------------------------------
        REAL*8 zeta, z, theta
        REAL*8 TTT2, TTT3
        REAL*8 coszeta, sinzeta, costheta, sintheta, cosz, sinz

        INCLUDE 'astmath.cmn'

        ! --------------------- PRECESSION angles ---------------------
        TTT2= TTT * TTT
        TTT3= TTT2 * TTT
        Zeta = 2306.2181D0*TTT + 0.30188D0*TTT2 + 0.017998D0*TTT3
        Theta= 2004.3109D0*TTT - 0.42665D0*TTT2 - 0.041833D0*TTT3
        z    = 2306.2181D0*TTT + 1.09468D0*TTT2 + 0.018203D0*TTT3

        Zeta = Zeta  * Deg2Rad / 3600.0D0
        Theta= Theta * Deg2Rad / 3600.0D0
        Z    = Z     * Deg2Rad / 3600.0D0

        coszeta   = DCOS(zeta)
        sinzeta   = DSIN(zeta)
        costheta = DCOS(theta)
        sintheta = DSIN(theta)
        cosz    = DCOS(z)
        sinz    = DSIN(z)

        ! ----------------- form matrix  J2000 to MOD -----------------
        prec(1,1) =  coszeta * costheta * cosz - sinzeta * sinz
        prec(1,2) = -sinzeta * costheta * cosz - coszeta * sinz
        prec(1,3) = -sintheta * cosz

        prec(2,1) =  coszeta * costheta * sinz + sinzeta * cosz
        prec(2,2) = -sinzeta * costheta * sinz + coszeta * cosz
        prec(2,3) = -sintheta * sinz

        prec(3,1) =  coszeta * sintheta
        prec(3,2) = -sinzeta * sintheta
        prec(3,3) =  costheta

      RETURN
      END   ! SUBROUTINE PRECESSION

* ----------------------------------------------------------------------------
*
*                           SUBROUTINE GCRF_MOD
*
*  this function transfroms a vector between the mean equator mean equinox of
*    epoch (eci) and the mean equator mean equinox of date (mod).
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    rJ2000      - Initial J2000 GCRF Position     ER, km, etc
*    vJ2000      - Initial J2000 GCRF Velocity
*    TTT         - Julian Centuries of TT         centuries
*    Direction   - Which set of vars to output    FROM  TOO
*
*  Outputs       :
*    rMOD        - Position vector of date        ER, km, etc
*                    Mean Equator, Mean Equinox
*    vMOD        - Velocity vector of date
*                    Mean Equator, Mean Equinox
*
*  Locals        :
*    Prec        - Precession transformation (eci-mod)
*    tmt         - transformation matrix transpose
*    Zeta        - PRECESSION ANGLE               rad
*    z           - PRECESSION ANGLE               rad
*    Theta       - PRECESSION ANGLE               rad
*
*  Coupling      :
*    PRECESSION  - Find precession matrix
*    ROT2        - Rotation about the second axis
*    ROT3        - Rotation about the third axis
*
*  References    :
*    Vallado       2007, 228
*
* ----------------------------------------------------------------------------

      SUBROUTINE GCRF_MOD  ( rJ2000, vJ2000, Direction, rMOD, vMOD, TTT)
        IMPLICIT NONE
        CHARACTER*4 Direction
        REAL*8 rJ2000(3), vJ2000(3), rMOD(3), vMOD(3), TTT

* ----------------------------  Locals  -------------------------------
        REAL*8 Prec(3,3),tmt(3,3)

        ! ---------------- get transformation matrices ----------------
        CALL Precession (TTT,Prec)

        ! ------------------- Perform matrix mmpy ---------------------
        IF ( Direction .eq. 'TOO ' ) THEN
            CALL MATMULT     ( Prec, rJ2000, 3,3,1,3,3,3, rmod )
            CALL MATMULT     ( Prec, vJ2000, 3,3,1,3,3,3, vmod )
          ELSE
            CALL MATTRANS( prec, 3,3, 3,3, tmt )
            CALL MATMULT     ( tmt, rmod  , 3,3,1,3,3,3, rj2000 )
            CALL MATMULT     ( tmt, vmod  , 3,3,1,3,3,3, vj2000 )
          ENDIF

c       ! --------------------- Perform rotations ---------------------
c        IF ( Direction .eq. 'TOO ' ) THEN
c            CALL ROT3( rJ2000, -Zeta, Temp1 )
c            CALL ROT2( Temp1 , Theta, Temp  )
c            CALL ROT3( Temp  ,  -z  , rMOD  )
c            CALL ROT3( vJ2000, -Zeta, Temp1 )
c            CALL ROT2( Temp1 , Theta, Temp  )
c            CALL ROT3( Temp  ,  -z  , vMOD  )
c          ELSE
c            CALL ROT3( rMOD  ,   z  , Temp   )
c            CALL ROT2( Temp  ,-Theta, Temp1  )
c            CALL ROT3( Temp1 ,  Zeta, rJ2000 )
c            CALL ROT3( vMOD  ,   z  , Temp   )
c            CALL ROT2( Temp  ,-Theta, Temp1  )
c            CALL ROT3( Temp1 ,  Zeta, vJ2000 )
c          ENDIF

      RETURN
      END   ! SUBROUTINE GCRF_MOD

* ----------------------------------------------------------------------------
*
*                           SUBROUTINE NUTATION
*
*  this function calulates the transformation matrix for nutation.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    TTT         - Julian Centuries of TT
*
*  Outputs       :
*    Nut         - Nutation transformation (mod-tod)
*    DeltaPsi    - NUTATION ANGLE                 rad
*    TrueEps     - True obliquity of the ecliptic rad
*    Omega       -                                rad
*    MeanEps     - Mean obliquity of the ecliptic rad
*
*  Locals        :
*    TTT2        - TTT squared
*    TTT3        - TTT cubed
*    MeanEps     - Mean obliquity of the ecliptic rad
*    l           -                                rad
*    ll          -                                rad
*    F           -                                rad
*    D           -                                rad
*    DeltaEps    - Change in obliquity            rad
*
*  Coupling      :
*    none
*
*  References    :
*    Vallado       2007, 228
*
* ----------------------------------------------------------------------------

      SUBROUTINE Nutation    ( TTT,
     &                         DeltaPsi, TrueEps, MeanEps, Omega, Nut )
        IMPLICIT NONE
        REAL*8  Omega,
     &          DeltaPsi, TrueEps, TTT, Nut(3,3), meaneps

* ----------------------------  Locals  -------------------------------
        INTEGER i
        REAL*8 cospsi, sinpsi, coseps, sineps, costrueeps, sintrueeps
        REAL*8 DeltaEps, Tempval, TTT4, TTT2, TTT3, l, l1, f, d

        INCLUDE 'astmath.cmn'
        INCLUDE 'astrolib.cmn'

        ! ---- Determine coefficients for IAU 1980 NUTATION Theory ----
        TTT2= TTT*TTT
        TTT3= TTT2*TTT
        TTT4= TTT2*TTT2

        MeanEps = -46.8150D0*TTT - 0.00059D0*TTT2 + 0.001813D0*TTT3 +
     &            84381.448D0
        MeanEps = DMOD( MeanEps/3600.0D0,360.0D0 )
        MeanEps = MeanEps * Deg2Rad

c     ---- Old values ----
c        rr   = 360.0
c        l    =  134.9629814 + (1325*rr  + 198.8673981)*TTT + 0.0086972 *TTT2 + 0.00001778*TTT3
c        l1   =  357.5277233 + (  99*rr  + 359.05034  )*TTT - 0.00016028*TTT2 - 0.00000333*TTT3
c        F    =   93.2719103 + (1342*rr  +  82.0175381)*TTT - 0.0036825 *TTT2 + 0.00000306*TTT3
c        D    =  297.8503631 + (1236*rr  + 307.111480 )*TTT - 0.00191417*TTT2 + 0.00000528*TTT3
c        Omega=  125.0445222 - (   5*rr  + 134.1362608)*TTT + 0.0020708 *TTT2 + 0.00000222*TTT3

        l    =  134.96340251D0 + ( 1717915923.2178D0*TTT +
     &          31.8792D0*TTT2 + 0.051635D0*TTT3 - 0.00024470D0*TTT4 )
     &          / 3600.0D0
        l1   =  357.52910918D0 + (  129596581.0481D0*TTT -
     &           0.5532D0*TTT2 - 0.000136D0*TTT3 - 0.00001149*TTT4 )
     &          / 3600.0D0
        F    =   93.27209062D0 + ( 1739527262.8478D0*TTT -
     &          12.7512D0*TTT2 + 0.001037D0*TTT3 + 0.00000417*TTT4 )
     &          / 3600.0D0
        D    =  297.85019547D0 + ( 1602961601.2090D0*TTT -
     &           6.3706D0*TTT2 + 0.006593D0*TTT3 - 0.00003169*TTT4 )
     &          / 3600.0D0
        Omega=  125.04455501D0 + (   -6962890.2665D0*TTT +
     &           7.4722D0*TTT2 + 0.007702D0*TTT3 - 0.00005939*TTT4 )
     &          / 3600.0D0

        l    = DMOD( l,360.0D0 )     * Deg2Rad
        l1   = DMOD( l1,360.0D0 )    * Deg2Rad
        F    = DMOD( F,360.0D0 )     * Deg2Rad
        D    = DMOD( D,360.0D0 )     * Deg2Rad
        Omega= DMOD( Omega,360.0D0 ) * Deg2Rad

        DeltaPsi= 0.0D0
        DeltaEps= 0.0D0

        DO i= 106, 1, -1
            Tempval= IAr80(1,i)*l + IAr80(2,i)*l1 + IAr80(3,i)*F +
     &               IAr80(4,i)*D + IAr80(5,i)*Omega
            DeltaPsi= DeltaPsi + (RAr80(1,i)+RAr80(2,i)*TTT) *
     &                 DSIN( TempVal )
            DeltaEps= DeltaEps + (RAr80(3,i)+RAr80(4,i)*TTT) *
     &                 DCOS( TempVal )
        ENDDO

        ! --------------- Find NUTATION Parameters --------------------
        DeltaPsi = DMOD( DeltaPsi,360.0D0 ) * Deg2Rad
        DeltaEps = DMOD( DeltaEps,360.0D0 ) * Deg2Rad
        TrueEps  = MeanEps + DeltaEps

        cospsi  = DCOS(deltapsi)
        sinpsi  = DSIN(deltapsi)
        coseps  = DCOS(meaneps)
        sineps  = DSIN(meaneps)
        costrueeps = DCOS(trueeps)
        sintrueeps = DSIN(trueeps)

        nut(1,1) =  cospsi
        nut(1,2) = -coseps * sinpsi
        nut(1,3) = -sineps * sinpsi

        nut(2,1) =  costrueeps * sinpsi
        nut(2,2) =  costrueeps * coseps * cospsi + sintrueeps * sineps
        nut(2,3) =  costrueeps * sineps * cospsi - sintrueeps * coseps

        nut(3,1) =  sintrueeps * sinpsi
        nut(3,2) =  sintrueeps * coseps * cospsi - sineps * costrueeps
        nut(3,3) =  sintrueeps * sineps * cospsi + costrueeps * coseps

      RETURN
      END   ! SUBROUTINE NUTATION
      
* ----------------------------------------------------------------------------
*
*                           SUBROUTINE GCRF_TOD
*
*  this function transforms vectors between the mean equator mean equinox of
*    date (j2000) and the true equator true equinox of date (tod).
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    rJ2000      - Initial J2000 GCRF Position     ER, km, etc
*    vJ2000      - Initial J2000 GCRF Velocity
*    TTT         - Julian Centuries of TT
*    Direction   - Which set of vars to output    FROM  TOO
*
*  Outputs       :
*    rTOD        - Position vector of date
*                    True Equator, True Equinox
*    vTOD        - Velocity vector of date
*                    True Equator, True Equinox
*
*  Locals        :
*    Prec        - Precession transformation (eci-mod)
*    Nut         - Nutation transformation (mod-tod)
*    tmt         - transformation matrix transpose
*    Omega       -                                rad
*    DeltaEps    - Change in obliquity            rad
*    DeltaPsi    - NUTATION ANGLE                 rad
*    TrueEps     - True obliquity of the ecliptic rad
*    MeanEps     - Mean obliquity of the ecliptic rad
*
*  Coupling      :
*    PRECESSION  - Find precession matrix
*    NUTATION    - Find nutation matrix
*
*  References    :
*    Vallado       2007, 228
*
* ----------------------------------------------------------------------------
      SUBROUTINE GCRF_TOD    ( rj2000, vj2000, Direction, rTOD, vTOD,
     &                        TTT )
        IMPLICIT NONE
        CHARACTER*4 Direction
        REAL*8 rj2000(3),vj2000(3),rtod(3),vtod(3), TTT

* ----------------------------  Locals  -------------------------------
        REAL*8  meaneps,omega,DeltaPsi, TrueEps,Prec(3,3),Nut(3,3),
     &          tm(3,3),tmt(3,3)

        INCLUDE 'astrolib.cmn'

        ! ---------------- get transformation matrices ----------------
        CALL Precession ( TTT,Prec)

        CALL nutation( ttt,  deltapsi,trueeps,meaneps,omega,nut)

        ! ------------------- Perform matrix mmpy ---------------------
        CALL MATMULT     ( Nut , Prec  , 3,3,3,3,3,3, tm )
        IF ( Direction .eq. 'TOO ' ) THEN
            CALL MATMULT     ( tm  , rJ2000, 3,3,1,3,3,3, rtod )
            CALL MATMULT     ( tm  , vJ2000, 3,3,1,3,3,3, vtod )
          ELSE
            CALL MATTRANS( tm, 3,3, 3,3, tmt )
            CALL MATMULT     ( tmt , rtod  , 3,3,1,3,3,3, rj2000 )
            CALL MATMULT     ( tmt , vtod  , 3,3,1,3,3,3, vj2000 )
          ENDIF

c        ! --------------------- Perform rotations ---------------------
c        IF ( Direction .eq. 'TOO ' ) THEN
c            CALL ROT1( rMOD ,   Eps       , Temp  )
c            CALL ROT3( Temp , -DeltaPsi, Temp1 )
c            CALL ROT1( Temp1, -TrueEps , rTOD  )
c            CALL ROT1( vMOD ,   Eps       , Temp  )
c            CALL ROT3( Temp , -DeltaPsi, Temp1 )
c            CALL ROT1( Temp1, -TrueEps , vTOD  )
c          ELSE
c            CALL ROT1( rTOD ,  TrueEps , Temp  )
c            CALL ROT3( Temp ,  DeltaPsi, Temp1 )
c            CALL ROT1( Temp1,  -Eps       , rMOD  )
c            CALL ROT1( vTOD ,  TrueEps , Temp  )
c            CALL ROT3( Temp ,  DeltaPsi, Temp1 )
c            CALL ROT1( Temp1,  -Eps       , vMOD  )
c          ENDIF

      RETURN
      END   ! SUBROUTINE GCRF_TOD

* ----------------------------------------------------------------------------
*
*                           SUBROUTINE SIDEREAL
*
*  this function calulates the transformation matrix that accounts for the
*    effects of nutation. Notice that deltaspi should not be moded to a
*    positive number because it is multiplied rather than used in a
*    trigonometric argument.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    JDUT1       - Julian Date of UT1             days from 4713 BC
*    DeltaPsi    - NUTATION ANGLE                 rad
*    MeanEps     - Mean obliquity of the ecliptic rad
*    Omega       -                                rad
*    LOD         - Excess length of day           sec
*    terms       - number of terms to include with ast 0, 2
*
*  Outputs       :
*    St          - Sidereal Time transformation (tod-pef)
*    StDot       - Sidereal Time rate transformation (tod-pef)
*    Omegaearth  - rotation of the earth          rad
*
*  Locals        :
*    GST         - Mean Greenwich SIDEREAL Time   0 to 2Pi rad
*    AST         - Apparent GST                   0 to 2Pi rad
*    Hr          - hour                           hr
*    minute         - minutes                        minute
*    SEC         - seconds                        SEC
*
*  Coupling      :
*    none
*
*  References    :
*    Vallado       2007, 228
*
* ----------------------------------------------------------------------------

      SUBROUTINE SIDEREAL    ( JDUT1,DeltaPsi,MeanEps,Omega,LOD,
     &                         ST,STDot,OmegaEarth,terms )
        IMPLICIT NONE
        REAL*8 JDUT1, DeltaPsi, GSTIME, Omega, LOD,
     &         ThetaSa,St(3,3),stdot(3,3),Omegaearth(3)
        EXTERNAL GSTIME
        INTEGER terms

* ----------------------------  Locals  -------------------------------
        REAL*8 GST, AST, Conv1, meaneps

        INCLUDE 'astmath.cmn'

        OmegaEarth(1) = 0.0D0
        OmegaEarth(2) = 0.0D0
        OmegaEarth(3) = 0.0D0

        Conv1 = pi / (180.0D0*3600.0D0)

        ! ------------------------ Find Mean GST ----------------------
        GST= GSTIME( JDUT1 )

        ! ------------------------ Find Mean AST ----------------------
        IF ((terms.gt.0) .and. (JDUT1.gt.2450449.5D0)) THEN
            AST= GST + DeltaPsi* DCOS(MeanEps)
     &           + 0.00264D0*Conv1*DSIN(Omega)
     &           + 0.000063D0*Conv1*DSIN(2.0D0*Omega)
          ELSE
            AST= GST + DeltaPsi* DCOS(MeanEps)
          ENDIF

        st(1,1) =  DCOS(ast)
        st(1,2) =  DSIN(ast)
        st(1,3) =  0.0

        st(2,1) = -DSIN(ast)
        st(2,2) =  DCOS(ast)
        st(2,3) =  0.0

        st(3,1) =  0.0
        st(3,2) =  0.0
        st(3,3) =  1.0

        ! ------------ compute sidereal time rate matrix --------------
        ThetaSA   =  7.29211514670698D-05 * (1.0D0 - LOD/86400.0D0)
        omegaearth(3) = ThetaSa

        stdot(1,1) = -omegaearth(3) * DSIN(ast)
        stdot(1,2) =  omegaearth(3) * DCOS(ast)
        stdot(1,3) =  0.0

        stdot(2,1) = -omegaearth(3) * DCOS(ast)
        stdot(2,2) = -omegaearth(3) * DSIN(ast)
        stdot(2,3) =  0.0

        stdot(3,1) =  0.0
        stdot(3,2) =  0.0
        stdot(3,3) =  0.0

      RETURN
      END   ! SUBROUTINE SIDEREAL Time

* ----------------------------------------------------------------------------
*
*                           SUBROUTINE GCRF_PEF
*
*  this function transforms a vector between the mean equator, mean equinox of epoch
*    frame (j2000), and the pseudo earth fixed frame (pef).
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    rJ2000      - Initial J2000 GCRF Position     ER, km, etc
*    vJ2000      - Initial J2000 GCRF Velocity
*    Direction   - Which set of vars to output    FROM  TOO
*    TTT         - Julian Centuries of TT
*    JDUT1       - Julian Date of UT1             days from 4713 BC
*    LOD         - Excess length of day           sec
*    terms       - number of terms to include with ast 0, 2
*
*  Outputs       :
*    rPEF        - Position vector of date        ER
*                    True Equator, True Equinox
*    vPEF        - Velocity vector of date        ER/TU
*                    True Equator, True Equinox
*
*  Locals        :
*    Prec        - Precession transformation (eci-mod)
*    Nut         - Nutation transformation (mod-tod)
*    St          - Sidereal Time transformation (tod-pef)
*    StDot       - Sidereal Time rate transformation (tod-pef)
*    tmt         - transformation matrix transpose
*    DeltaPsi    - NUTATION ANGLE                 rad
*    TrueEps     - True obliquity of the ecliptic rad
*    MeanEps     - Mean obliquity of the ecliptic rad
*    Omegaearth  - rotation of the earth          rad
*    Omega       -                                rad
*    GST         - Mean Greenwich SIDEREAL Time   0 to 2Pi rad
*    AST         - Apparent GST                   0 to 2Pi rad
*
*  Coupling      :
*    PRECESSION  - Find precession matrix
*    NUTATION    - Find nutation matrix
*    SIDEREAL    - Find sidereal time matrix
*
*
*  References    :
*    Vallado       2007, 228
*
* ----------------------------------------------------------------------------

      SUBROUTINE GCRF_PEF    ( rj2000, vj2000, Direction, rPEF, vPEF,
     &                        TTT, JDUT1, LOD, terms )
        IMPLICIT NONE
        CHARACTER*4 Direction
        REAL*8 rj2000(3), vj2000(3), rPEF(3), vPEF(3), JDUT1, LOD
        INTEGER Terms
        EXTERNAL GSTIME

* ----------------------------  Locals  -------------------------------
        REAL*8 ttt,meaneps,omegaearth(3),wcrossr(3),
     &         TrueEps, GSTIME, Omega, tm(3,3),tmt(3,3),
     &         Prec(3,3),Nut(3,3),st(3,3), DeltaPsi,
     &         stdot(3,3),tm1(3,3)

        INCLUDE 'astrolib.cmn'
     
        ! ---------------- get transformation matrices ----------------
        CALL Precession ( TTT,Prec)

        CALL nutation( ttt,  deltapsi,trueeps,meaneps,omega,nut)

        CALL sidereal( jdut1,deltapsi,meaneps,omega,lod,  st,stdot,
     &                 OmegaEarth,terms )

        ! ------------------- Perform matrix mmpy ---------------------
        CALL MATMULT     ( Nut , Prec  , 3,3,3,3,3,3, tm1 )
        CALL MATMULT     ( St  , tm1    , 3,3,3,3,3,3, tm )
        IF ( Direction .eq. 'TOO ' ) THEN
            CALL MATMULT     ( tm  , rJ2000, 3,3,1,3,3,3, rPEF )
            CALL MATMULT     ( tm  , vJ2000, 3,3,1,3,3,3, vPEF )
            CALL CROSS( OmegaEarth,rPEF,  wcrossr )
            CALL SUBVEC( vPEF, wcrossr,  vPEF )
          ELSE
            CALL MATTRANS( tm, 3,3, 3,3, tmt )
            CALL MATMULT     ( tmt , rPEF  , 3,3,1,3,3,3, rj2000 )
            CALL CROSS( OmegaEarth,rPEF,  wcrossr )
            CALL ADDVEC( vPEF, wcrossr,  vPEF )
            CALL MATMULT     ( tmt , vPEF  , 3,3,1,3,3,3, vj2000 )
          ENDIF

c       ! --------------------- Perform rotations ---------------------
c       IF ( Direction .eq. 'TOO ' ) THEN
c           CALL ROT3( rTOD,   AST  , rPEF )
c           CALL ROT3( vTOD,   AST  , Temp  )
c           vPEF(1)= Temp(1) + rPEF(2)*OmegaEarth(3)
c           vPEF(2)= Temp(2) - rPEF(1)*OmegaEarth(3)
c           vPEF(3)= Temp(3)
c         ELSE
c           CALL ROT3( rPEF, -AST  , rTOD )
c           DO i = 1,4
c               Temp(i)   = vPEF(i)
c             ENDDO
c           Temp(1)= vPEF(1) - rPEF(2)*OmegaEarth(3)
c           Temp(2)= vPEF(2) + rPEF(1)*OmegaEarth(3)
c           CALL ROT3( Temp,  -AST  , vTOD )
c         ENDIF

      RETURN
      END   ! SUBROUTINE GCRF_PEF

* ----------------------------------------------------------------------------
*
*                           SUBROUTINE POLARM
*
*  this function calulates the transformation matrix for polar motion.
*    the units for polar motion are input in rad because it's more
*    efficient to do this in the main routine.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    xp          - Polar motion coefficient       rad
*    yp          - Polar motion coefficient       rad
*
*  Outputs       :
*    PM          - Polar motion transformation (pef-ecef)
*
*  Locals        :
*    None.
*
*  Coupling      :
*    ROT2        - Rotation about the second axis
*    ROT3        - Rotation about the third axis
*
*  References    :
*    Vallado       2007, 228
*
* ----------------------------------------------------------------------------

      SUBROUTINE POLARM      ( xp, yp,  PM )
        IMPLICIT NONE
        REAL*8 xp, yp, PM(3,3)

* ----------------------------  Locals  -------------------------------
        REAL*8 cosxp,cosyp,sinxp,sinyp

        cosxp = DCOS(xp)
        sinxp = DSIN(xp)
        cosyp = DCOS(yp)
        sinyp = DSIN(yp)

        pm(1,1) =  cosxp
        pm(1,2) =  sinxp * sinyp
        pm(1,3) =  sinxp * cosyp

        pm(2,1) =  0.0
        pm(2,2) =  cosyp
        pm(2,3) = -sinyp

        pm(3,1) = -sinxp
        pm(3,2) =  cosxp * sinyp
        pm(3,3) =  cosxp * cosyp

        ! Approximate matrix using small angle approximations
c       pm(1,1) =  1.0
c       pm(1,2) =  0.0
c       pm(1,3) =  xp

c       pm(2,1) =  0.0
c       pm(2,2) =  1.0
c       pm(2,3) = -yp

c       pm(3,1) = -xp
c       pm(3,2) =  yp
c       pm(3,3) =  1.0

      RETURN
      END   ! SUBROUTINE PolarM

* ----------------------------------------------------------------------------
*
*                           SUBROUTINE GCRF_ITRF
*
*  this function calulates the reduction formula do a vector going from the
*    earth centered inertial frame (eci, eci) to the earth fixed (ITRF)
*    frame.  the results take into account the effects of precession, nutation,
*    sidereal time, and polar motion.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    rJ2000      - Initial J2000 GCRF Position     ER, km, etc
*    vJ2000      - Initial J2000 GCRF Velocity
*    Direction   - Which set of vars to output    FROM  TOO
*    TTT         - Julian Centuries of TT         centuries
*    JDUT1       - Julian Date of UT1             days from 4713 BC
*    LOD         - Excess length of day           sec
*    xp          - Polar motion coefficient       rad
*    yp          - Polar motion coefficient       rad
*    terms       - number of terms to include with ast 0, 2
*
*  Outputs       :
*    rITRF       - Position vector of date
*                    True Equator, True Equinox
*    vITRF       - Velocity vector of date
*                    True Equator, True Equinox
*
*  Locals        :
*    Prec        - Precession transformation (eci-mod)
*    Nut         - Nutation transformation (mod-tod)
*    St          - Sidereal Time transformation (tod-pef)
*    StDot       - Sidereal Time rate transformation (tod-pef)
*    PM          - Polar motion transformation (pef-ecef)
*    tmt         - transformation matrix transpose
*    DeltaPsi    - NUTATION ANGLE                 rad
*    TrueEps     - True obliquity of the ecliptic rad
*    MeanEps     - Mean obliquity of the ecliptic rad
*    Omegaearth  - rotation of the earth          rad
*
*  Coupling      :
*    PRECESSION  - Find precession matrix
*    NUTATION    - Find nutation matrix
*    SIDEREAL    - Find sidereal time matrix
*    POLARM      - Find polar motion matrix
*    ROT2        - Rotation about the second axis
*    ROT3        - Rotation about the third axis
*
*  References    :
*    Vallado       2007, 228
*
* ----------------------------------------------------------------------------

      SUBROUTINE GCRF_ITRF   ( rj2000,vj2000, Direction, rITRF,vITRF,
     &                        TTT, JDUT1, LOD, xp, yp, terms )
        IMPLICIT NONE
        CHARACTER*4 Direction
        REAL*8 rPEF(3), vPEF(3), rITRF(3), vITRF(3), xp, yp,
     &         TTT, Jdut1, LOD
        INTEGER terms

* ----------------------------  Locals  -------------------------------
        REAL*8 Deltapsi,meaneps,trueeps,omega,
     &         nut(3,3), st(3,3), stdot(3,3), pm(3,3), tm(3,3),
     &         tmt(3,3),  rj2000(3),vj2000(3),prec(3,3),
     &         omegaearth(3),wcrossr(3),tm1(3,3)

        INCLUDE 'astrolib.cmn'
     
        CALL precession(ttt,  prec)

        CALL nutation( ttt,  deltapsi,trueeps,meaneps,omega,nut)

        CALL sidereal(jdut1,deltapsi,meaneps,omega,lod,  st,stdot,
     &                 OmegaEarth,terms )

        CALL polarm(xp,yp, pm)

        ! ------------------- Perform matrix mmpy ---------------------
        CALL MATMULT     ( Nut , Prec  , 3,3,3,3,3,3, tm1  )
        CALL MATMULT     ( St  , tm1    , 3,3,3,3,3,3, tm  )
        IF ( Direction .eq. 'TOO ' ) THEN
            CALL MATMULT     ( tm  , rJ2000, 3,3,1,3,3,3, rPEF  )
            CALL MATMULT     ( pm  , rPEF  , 3,3,1,3,3,3, ritrf )

            CALL MATMULT     ( tm  , vJ2000, 3,3,1,3,3,3, vPEF  )
            CALL CROSS       ( omegaearth, rPEF, wcrossr )
            CALL SUBVEC      ( vPEF, wcrossr, vPEF )
            CALL MATMULT     ( pm  , vPEF  , 3,3,1,3,3,3, vitrf )
          ELSE
            CALL MATTRANS( tm, 3,3, 3,3, tmt )
            CALL MATMULT     ( tmt , ritrf  , 3,3,1,3,3,3, rj2000 )
            CALL MATMULT     ( tmt , vitrf  , 3,3,1,3,3,3, vj2000 )
          ENDIF

chk these accel
c        CALL cross(omegaearth,rPEF,  wcrossr)
c        aecef = pm*(st*nut*prec*aeci - stdot*rPEF - cross(omegaearth,temp) ...
c                - 2.0*cross(omegaearth,vPEF))

c        Write(20,*) 'xy,yp  ',xp*3600*57.29577,' ',yp*3600*57.29577
c        IF ( Direction .eq. 'TOO ' ) THEN
c            rITRF(1)= rPEF(1) + xp*rPEF(3)
c            rITRF(2)= rPEF(2) - yp*rPEF(3)
c            rITRF(3)= rPEF(3) - xp*rPEF(1) + yp*rPEF(2)
c            vITRF(1)= vPEF(1) + xp*vPEF(3)
c            vITRF(2)= vPEF(2) - yp*vPEF(3)
c            vITRF(3)= vPEF(3) - xp*vPEF(1) + yp*vPEF(2)
c          ELSE
c            rPEF(1)= rITRF(1) - xp*rITRF(3)
c            rPEF(2)= rITRF(2) + yp*rITRF(3)
c            rPEF(3)= rITRF(3) + xp*rITRF(1) - yp*rITRF(2)
c            vPEF(1)= vITRF(1) - xp*vITRF(3)
c            vPEF(2)= vITRF(2) + yp*vITRF(3)
c            vPEF(3)= vITRF(3) + xp*vITRF(1) - yp*vITRF(2)
c          ENDIF

      RETURN
      END   ! SUBROUTINE GCRF_ITRF

* ----------------------------------------------------------------------------
*
*                           PROCEDURE TRUEMEAN
*
*  this subroutine calculates the transformation matrix to teme.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    rJ2000      - Initial J2000 GCRF Position     ER, km, etc
*    vJ2000      - Initial J2000 GCRF Velocity
*    TTT         - Julian Centuries of TT
*    Order       - Number of coefficients used in nutation
*    terms       - number of terms to include with ast 0, 2
*
*  Outputs       :
*    rTm         - Position vector of date
*                    Mean Equator, True Equinox
*    vTm         - Velocity vector of date
*                    Mean Equator, True Equinox
*    DeltaPsiSp  - NUTATION ANGLE                 rad
*    TrueEpsSp   - True obliquity of the ecliptic rad
*
*  Locals        :
*    TTT2        - TTT squared
*    TTT3        - TTT cubed
*    Eps         - Mean obliquity of the ecliptic rad
*    l           -                                rad
*    ll          -                                rad
*    F           -                                rad
*    D           -                                rad
*    Omega       -                                rad
*    DeltaEps    - Change in obliquity            rad
*
*  Coupling      :
*    ROT2        - Rotation about the second axis
*    ROT3        - Rotation about the third axis
*
*  References    :
*    Vallado       2007, 236
*
* ----------------------------------------------------------------------------

      SUBROUTINE TrueMean    ( Order, TTT, Terms, tm )
        IMPLICIT NONE
        REAL*8  TTT, Nutteme(3,3)
        INTEGER Order, Terms

* ----------------------------  Locals  -------------------------------
        INTEGER ii, i
        REAL*8  DeltaEps,Tempval, TTT2, TTT3, rr, l, l1, f, d,
     &          Omega, DeltaPsi,TrueEps,Conv1,
     &          Prec(3,3),meaneps,st(3,3),eqe,tm(3,3),tm1(3,3)
        REAL*8 cospsi,sinpsi,coseps,sineps,costrueeps,sintrueeps

        INCLUDE 'astmath.cmn'
        INCLUDE 'astrolib.cmn'
        
        ! ----- Determine coefficients for IAU 1980 NUTATION Theory ----
        Conv1 = pi / (3600.0D0*180.0D0)

        ! ---------------- get transformation matrices ----------------
        CALL Precession ( TTT,Prec)

       TTT2= TTT*TTT
       TTT3= TTT2*TTT
       MeanEps  = 23.439291 - 0.0130042*TTT - 0.000000164*TTT2 +
     &             0.000000504*TTT3
       MeanEps = DMOD( MeanEps,360.0D0 )
       MeanEps = MeanEps * Deg2Rad

       rr   = 360.0  !deg
       l    =  134.9629814 + (1325*rr  + 198.8673981)*TTT +
     &           0.0086972 *TTT2 + 0.00001778*TTT3
       l1   =  357.5277233 + (  99*rr  + 359.05034  )*TTT -
     &           0.00016028*TTT2 - 0.00000333*TTT3
       F    =   93.2719103 + (1342*rr  +  82.0175381)*TTT -
     &           0.0036825 *TTT2 + 0.00000306*TTT3
       D    =  297.8503631 + (1236*rr  + 307.111480 )*TTT -
     &           0.00191417*TTT2 + 0.00000528*TTT3
       Omega=  125.0445222 - (   5*rr  + 134.1362608)*TTT +
     &           0.0020708 *TTT2 + 0.00000222*TTT3
       l    = DMOD( l,360.0D0 )     * Deg2Rad
       l1   = DMOD( l1,360.0D0 )    * Deg2Rad
       F    = DMOD( F,360.0D0 )     * Deg2Rad
       D    = DMOD( D,360.0D0 )     * Deg2Rad
       Omega= DMOD( Omega,360.0D0 ) * Deg2Rad

       DeltaPsi= 0.0
       DeltaEps= 0.0
       DO ii= 1, Order ! make sure the datafile is in the correct order
          i = ii

         Tempval= IAr80(1,i)*l + IAr80(2,i)*l1 + IAr80(3,i)*F +
     &            IAr80(4,i)*D + IAr80(5,i)*Omega
         DeltaPsi= DeltaPsi + (RAr80(1,i)+RAr80(2,i)*TTT) *
     &             DSIN( TempVal )
         DeltaEps= DeltaEps + (RAr80(3,i)+RAr80(4,i)*TTT) *
     &             DCOS( TempVal )
       ENDDO

       ! --------------- Find Approx Nutation Parameters --------------
       DeltaPsi = DMOD( DeltaPsi,360.0D0 ) * Deg2Rad
       DeltaEps = DMOD( DeltaEps,360.0D0 ) * Deg2Rad
       TrueEps  = MeanEps + DeltaEps

        cospsi  = DCOS(deltapsi)
        sinpsi  = DSIN(deltapsi)
        coseps  = DCOS(meaneps)
        sineps  = DSIN(meaneps)
        costrueeps = DCOS(trueeps)
        sintrueeps = DSIN(trueeps)

        nutteme(1,1) =  cospsi
        nutteme(1,2) = -coseps * sinpsi
        nutteme(1,3) = -sineps * sinpsi

        nutteme(2,1) =  costrueeps * sinpsi
        nutteme(2,2) =  costrueeps * coseps * cospsi +
     &                   sintrueeps * sineps
        nutteme(2,3) =  costrueeps * sineps * cospsi -
     &                   sintrueeps * coseps

        nutteme(3,1) =  sintrueeps * sinpsi
        nutteme(3,2) =  sintrueeps * coseps * cospsi -
     &                   sineps * costrueeps
        nutteme(3,3) =  sintrueeps * sineps * cospsi +
     &                   costrueeps * coseps

        IF (terms.gt.0) THEN
            eqe= DeltaPsi* DCOS(MeanEps)
     &           + 0.00264D0*Conv1*DSIN(Omega)
     &           + 0.000063D0*Conv1*DSIN(2.0D0*Omega)
          ELSE
            eqe= DeltaPsi* DCOS(MeanEps)
          ENDIF

        st(1,1) =  DCOS(eqe)
        st(1,2) =  DSIN(eqe)
        st(1,3) =  0.0

        st(2,1) = -DSIN(eqe)
        st(2,2) =  DCOS(eqe)
        st(2,3) =  0.0

        st(3,1) =  0.0
        st(3,2) =  0.0
        st(3,3) =  1.0

        CALL MATMULT     ( nutteme, prec , 3,3,3,3,3,3, tm1 )
        CALL MATMULT     ( st , tm1      , 3,3,3,3,3,3, tm )

      RETURN
      END   ! Subroutine TrueMean

* ----------------------------------------------------------------------------
*
*                           SUBROUTINE GCRF_TEME
*
*  this function transforms a vector from the mean equator, mean equinox of date
*    frame (j2000), to the true equator mean equinox frame (teme).
*
*  Author        : David Vallado                  719-573-2600   26 Sep 2002
*
*  Inputs          Description                    Range / Units
*    rJ2000      - Initial J2000 GCRF Position     ER, km, etc
*    vJ2000      - Initial J2000 GCRF Velocity
*    Direction   - Which set of vars to output    FROM  TOO
*    DeltaPsi    - NUTATION ANGLE                 rad
*    TrueEps     - True obliquity of the ecliptic rad
*    Omega       -                                rad
*    terms       - number of terms to include with ast 0, 2
*
*  Outputs       :
*    rTm         - Position vector of date
*                    True Equator, Mean Equinox
*    vTm         - Velocity vector of date
*                    True Equator, Mean Equinox
*
*  Locals        :
*    GST         - Mean Greenwich SIDEREAL Time   0 to 2Pi rad
*    AST         - Apparent GST                   0 to 2Pi rad
*    Hr          - hour                           hr
*    minute         - minutes                        minute
*    SEC         - seconds                        SEC
*    Temp        - Temporary vector
*    TempVal     - Temporary variable
*
*  Coupling      :
*
*
*  References    :
*    Vallado       2007, 236
*
* ----------------------------------------------------------------------------

      SUBROUTINE GCRF_TEME    ( rj2000, vj2000, Direction, rTM, vTM,
     &                         Order, TTT, Terms )
        IMPLICIT NONE
        CHARACTER*4 Direction
        REAL*8 rj2000(3), vj2000(3), rTM(3), vTM(3), TTT
        INTEGER Order, Terms

* ----------------------------  Locals  -------------------------------
        REAL*8  nutteme(3,3), tmt(3,3)

        INCLUDE 'astrolib.cmn'
        
        CALL TrueMean ( Order, TTT, Terms, Nutteme )

        ! ------------------- Perform matrix mmpy ---------------------
        IF ( Direction .eq. 'TOO ' ) THEN
            CALL MATMULT     ( nutteme , rJ2000, 3,3,1,3,3,3, rtm )
            CALL MATMULT     ( nutteme , vJ2000, 3,3,1,3,3,3, vtm )
          ELSE
            CALL MATTRANS( nutteme, 3,3, 3,3, tmt )
            CALL MATMULT     ( tmt , rtm  , 3,3,1,3,3,3, rj2000 )
            CALL MATMULT     ( tmt , vtm  , 3,3,1,3,3,3, vj2000 )
          ENDIF


c       IF (Direction .eq. 'TOO ') THEN
c           rtm(1)= rmod(1) - DeltaPsi*DSIN(eps)*rmod(3)
c           rtm(2)= rmod(2) - DeltaEps*rmod(3)
c           rtm(3)= rmod(3) +DeltaPsi*DSIN(eps)*rmod(1) +
c     &             Deltaeps*rmod(2)
c           vtm(1)= vmod(1) - DeltaPsi*DSIN(eps)*vmod(3)
c           vtm(2)= vmod(2) - DeltaEps*vmod(3)
c           vtm(3)= vmod(3) +DeltaPsi*DSIN(eps)*vmod(1) +
c     &             Deltaeps*vmod(2)
c         ELSE
c           rmod(1)= rtm(1) + DeltaPsi*DSIN(eps)*rtm(3)
c           rmod(2)= rtm(2) + Deltaeps*rtm(3)
c           rmod(3)= rtm(3) - DeltaPsi*DSIN(eps)*rtm(1) -
c     &              Deltaeps*rtm(2)
c           vmod(1)= vtm(1) + DeltaPsi*DSIN(eps)*vtm(3)
c           vmod(2)= vtm(2) + Deltaeps*vtm(3)
c           vmod(3)= vtm(3) - DeltaPsi*DSIN(eps)*vtm(1) -
c     &              Deltaeps*vtm(2)
c         ENDIF

      RETURN
      END   ! Subroutine GCRF_TEME

      
* -----------------------------------------------------------------------------------------
*                                       2body functions
* -----------------------------------------------------------------------------------------
      

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE rv2coe
*
*  this subroutine finds the classical orbital elements given the Geocentric
*    Equatorial Position and Velocity vectors.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R           - IJK Position vector            km
*    V           - IJK Velocity vector            km / s
*
*  Outputs       :
*    P           - SemiLatus rectum               km
*    A           - semimajor axis                 km
*    Ecc         - Eccentricity
*    Incl        - inclination                    0.0D0 to Pi rad
*    Omega       - Longitude of Ascending Node    0.0D0 to 2Pi rad
*    Argp        - Argument of Perigee            0.0D0 to 2Pi rad
*    Nu          - True anomaly                   0.0D0 to 2Pi rad
*    M           - Mean anomaly                   0.0D0 to 2Pi rad
*    ArgLat      - Argument of Latitude      (CI) 0.0D0 to 2Pi rad
*    LamTrue     - True Longitude            (CE) 0.0D0 to 2Pi rad
*    LonPer      - Longitude of Periapsis    (EE) 0.0D0 to 2Pi rad
*
*  Locals        :
*    HBar        - Angular Momentum H Vector      km2 / s
*    EBar        - Eccentricity     E Vector
*    NBar        - Line of Nodes    N Vector
*    c1          - V**2 - u/R
*    RDotV       - R DOT V
*    Hk          - Hk norm vector
*    SME         - Specfic Mechanical Energy      km2 / s2
*    i           - index
*    E           - Eccentric, Parabolic,
*                  Hyperbolic Anomaly             rad
*    Temp        - Temporary variable
*    TypeOrbit   - Type of orbit                  EE, EI, CE, CI
*
*  Coupling      :
*    MAG         - Magnitude of a vector
*    CROSS       - CROSS product of two vectors
*    DOT         - DOT product of two vectors
*    ANGLE       - Find the ANGLE between two vectors
*    NEWTONNU    - Find the mean anomaly
*
*  References    :
*    Vallado       2007, 121, Alg 9, Ex 2-5
*
* ------------------------------------------------------------------------------

      SUBROUTINE rv2coe      ( R, V, P, A, Ecc, Incl, Omega, Argp, Nu,
     &                         M, ArgLat, TrueLon, LonPer )
        IMPLICIT NONE
        REAL*8 R(3), V(3), P, A, Ecc, Incl, Omega, Argp, Nu, M, ArgLat,
     &         TrueLon, LonPer
        EXTERNAL DOT, MAG
* -----------------------------  Locals  ------------------------------
        REAL*8 c1, RDotV, hk, SME, Hbar(3), Ebar(3), Nbar(3),
     &         Dot, E, Temp, MAG, maghbar, magnbar, magr, magv
        INTEGER i
        CHARACTER*2 TypeOrbit

        INCLUDE 'astmath.cmn'
        INCLUDE 'astrolib.cmn'

        ! --------------------  Implementation   ----------------------
        magr = MAG( R )
        magv = MAG( V )
        ! ------------------  Find H N and E vectors   ----------------
        CALL CROSS( R, V, HBar )
        maghbar = MAG(Hbar)
        IF ( maghbar .ge. 0.0D0 ) THEN
            NBar(1)= -HBar(2)
            NBar(2)=  HBar(1)
            NBar(3)=   0.0D0
            magnbar = MAG( Nbar )
            c1 = magv**2 - mu/magr
            RDotV= DOT( R, V )
            DO i= 1 , 3
                EBar(i)= (c1*R(i) - RDotV*V(i))/mu
              ENDDO

            Ecc = MAG( EBar )

            ! ------------  Find a e and semi-Latus rectum   ----------
            SME= ( magv*magv*0.5D0 ) - ( mu/magr )
            IF ( DABS( SME ) .gt. Small ) THEN
                A= -mu / (2.0D0*SME)
              ELSE
                A= Infinite
              ENDIF
            P = maghbar*maghbar/mu

            ! -----------------  Find inclination   -------------------
            Hk= HBar(3)/maghbar
c            IF ( DABS( DABS(Hk) - 1.0D0 ) .lt. Small ) THEN
c                ! -------------  Equatorial Orbits   ------------------
c                IF ( DABS(HBar(3)) .gt. 0.0D0 ) THEN
c                    Hk= DSIGN(1.0D0, HBar(3))
c                  ENDIF
c              ENDIF
            Incl= DACOS( Hk ) 

            ! --------  Determine type of orbit for Later use  --------
            ! ------ Elliptical, Parabolic, Hyperbolic Inclined -------
            TypeOrbit= 'EI' 
            IF ( Ecc .lt. Small ) THEN
                ! ----------------  Circular Equatorial ---------------
                IF ( (Incl.lt.Small).or.(DABS(Incl-Pi).lt.Small) ) THEN
                    TypeOrbit= 'CE'
                  ELSE
                    ! --------------  Circular Inclined ---------------
                    TypeOrbit= 'CI'
                  ENDIF
              ELSE
                ! - Elliptical, Parabolic, Hyperbolic Equatorial --
                IF ( (Incl.lt.Small).or.(DABS(Incl-Pi).lt.Small) ) THEN
                    TypeOrbit= 'EE'
                  ENDIF
              ENDIF

            ! ----------  Find Longitude of Ascending Node ------------
            IF ( magnbar .gt. Small ) THEN
                Temp= NBar(1) / magnbar
                IF ( DABS(Temp) .gt. 1.0D0 ) THEN
                    Temp= DSIGN(1.0D0, Temp)
                  ENDIF
                Omega= DACOS( Temp ) 
                IF ( NBar(2) .lt. 0.0D0 ) THEN
                    Omega= TwoPi - Omega
                  ENDIF
              ELSE
                Omega= Undefined 
              ENDIF

            ! ---------------- Find Argument of perigee ---------------
            IF ( TypeOrbit .eq. 'EI' ) THEN
                CALL ANGLE( NBar, EBar, Argp )
                IF ( EBar(3) .lt. 0.0D0 ) THEN
                    Argp= TwoPi - Argp 
                  ENDIF
              ELSE
                Argp= Undefined 
              ENDIF

            ! ------------  Find True Anomaly at Epoch    -------------
            IF ( TypeOrbit(1:1) .eq. 'E' ) THEN
                CALL ANGLE( EBar, r, Nu )
                IF ( RDotV .lt. 0.0D0 ) THEN
                    Nu= TwoPi - Nu 
                  ENDIF
              ELSE
                Nu= Undefined 
              ENDIF

            ! ----  Find Argument of Latitude - Circular Inclined -----
            IF ( TypeOrbit .eq. 'CI' ) THEN
                CALL ANGLE( NBar, R, ArgLat )
                IF ( R(3) .lt. 0.0D0 ) THEN
                    ArgLat= TwoPi - ArgLat
                  ENDIF
              ELSE
                ArgLat= Undefined 
              ENDIF

            ! -- Find Longitude of Perigee - Elliptical Equatorial ----
            IF ( ( Ecc.gt.Small ) .and. (TypeOrbit.eq.'EE') ) THEN
                Temp= EBar(1)/Ecc
                IF ( DABS(Temp) .gt. 1.0D0 ) THEN
                    Temp= DSIGN(1.0D0, Temp)
                  ENDIF
                LonPer= DACOS( Temp ) 
                IF ( EBar(2) .lt. 0.0D0 ) THEN
                    LonPer= TwoPi - LonPer 
                  ENDIF
                IF ( Incl .gt. HalfPi ) THEN
                    LonPer= TwoPi - LonPer
                  ENDIF
              ELSE
                LonPer= Undefined
              ENDIF

            ! -------- Find True Longitude - Circular Equatorial ------
            IF ( ( magr.gt.Small ) .and. ( TypeOrbit.eq.'CE' ) ) THEN
                Temp= R(1)/magr
                IF ( DABS(Temp) .gt. 1.0D0 ) THEN
                    Temp= DSIGN(1.0D0, Temp)
                  ENDIF
                TrueLon= DACOS( Temp )
                IF ( R(2) .lt. 0.0D0 ) THEN
                    TrueLon= TwoPi - TrueLon
                  ENDIF
                IF ( Incl .gt. HalfPi ) THEN
                    TrueLon= TwoPi - TrueLon
                  ENDIF
              ELSE
                TrueLon= Undefined
              ENDIF

            ! ------------ Find Mean Anomaly for all orbits -----------
            CALL NEWTONNU(Ecc, Nu, E, M )

         ELSE
           P    = Undefined
           A    = Undefined
           Ecc  = Undefined
           Incl = Undefined
           Omega= Undefined 
           Argp = Undefined 
           Nu   = Undefined 
           M    = Undefined 
           ArgLat  = Undefined 
           TrueLon= Undefined 
           LonPer = Undefined 
         ENDIF 

      RETURN
      END  ! end rv2coe

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE coe2rv
*
*  this subroutine finds the position and velocity vectors in Geocentric
*    Equatorial (IJK) system given the classical orbit elements.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    P           - SemiLatus rectum               km
*    Ecc         - Eccentricity
*    Incl        - inclination                    0.0D0 to Pi rad
*    Omega       - Longitude of Ascending Node    0.0D0 to 2Pi rad
*    Argp        - Argument of Perigee            0.0D0 to 2Pi rad
*    Nu          - True anomaly                   0.0D0 to 2Pi rad
*    ArgLat      - Argument of Latitude      (CI) 0.0D0 to 2Pi rad
*    LamTrue     - True Longitude            (CE) 0.0D0 to 2Pi rad
*    LonPer      - Longitude of Periapsis    (EE) 0.0D0 to 2Pi rad
*
*  Outputs       :
*    R           - IJK Position vector            km
*    V           - IJK Velocity vector            km / s
*
*  Locals        :
*    Temp        - Temporary REAL*8 value
*    Rpqw        - PQW Position vector            km
*    Vpqw        - PQW Velocity vector            km / s
*    SinNu       - Sine of Nu
*    CosNu       - Cosine of Nu
*    TempVec     - PQW Velocity vector
*
*  Coupling      :
*    ROT3        - Rotation about the 3rd axis
*    ROT1        - Rotation about the 1st axis
*
*  References    :
*    Vallado       2007, 126, Alg 10, Ex 2-5
*
* ------------------------------------------------------------------------------

      SUBROUTINE coe2rv      ( P, Ecc, Incl, Omega, Argp, Nu, ArgLat,
     &                         TrueLon, LonPer, R, V )
        IMPLICIT NONE
        REAL*8 R(3), V(3), P, Ecc, Incl, Omega, Argp, Nu, ArgLat,
     &         TrueLon, LonPer
* -----------------------------  Locals  ------------------------------
        REAL*8 Rpqw(3), Vpqw(3), TempVec(3), Temp, SinNu, CosNu

        INCLUDE 'astmath.cmn'
        INCLUDE 'astrolib.cmn'

        ! --------------------  Implementation   ----------------------
*       Determine what type of orbit is involved and set up the
*       set up angles for the special cases.
        ! -------------------------------------------------------------
        IF ( Ecc .lt. Small ) THEN
            ! ----------------  Circular Equatorial  ------------------
            IF ( (Incl.lt.Small).or.( DABS(Incl-Pi).lt. Small ) ) THEN
                Argp = 0.0D0
                Omega= 0.0D0 
                Nu   = TrueLon 
              ELSE
                ! --------------  Circular Inclined  ------------------
                Argp= 0.0D0
                Nu  = ArgLat 
              ENDIF
          ELSE
            ! ---------------  Elliptical Equatorial  -----------------
            IF ( ( Incl.lt.Small) .or. (DABS(Incl-Pi).lt.Small) ) THEN
                Argp = LonPer
                Omega= 0.0D0 
              ENDIF 
          ENDIF

        ! ----------  Form PQW position and velocity vectors ----------
        CosNu= DCOS(Nu)
        SinNu= DSIN(Nu)
        Temp = P / (1.0D0 + Ecc*CosNu)
        Rpqw(1)= Temp*CosNu
        Rpqw(2)= Temp*SinNu
        Rpqw(3)=     0.0D0
        IF ( DABS(p) .lt. 0.00000001D0 ) THEN
            p= 0.00000001D0
          ENDIF
        Vpqw(1)=    -SinNu    * DSQRT(mu/P)
        Vpqw(2)=  (Ecc + CosNu) * DSQRT(mu/P)
        Vpqw(3)=      0.0D0

        ! ----------------  Perform transformation to IJK  ------------
        CALL ROT3( Rpqw   , -Argp , TempVec )
        CALL ROT1( TempVec, -Incl , TempVec )
        CALL ROT3( TempVec, -Omega,  R     )

        CALL ROT3( Vpqw   , -Argp , TempVec )
        CALL ROT1( TempVec, -Incl , TempVec )
        CALL ROT3( TempVec, -Omega, V     )

      RETURN
      END

* ----------------------------------------------------------------------------
*
*                           function flt2rv.m
*
*  this function transforms  the flight elements - latgc, lon, fpav, az,
*    position and velocity magnitude into an eci position and velocity vector.
*
*  author        : david vallado                  719-573-2600   17 jun 2002
*
*  revisions
*    vallado     - fix extra terms in rtasc calc                  8 oct 2002
*
*  inputs          description                    range / units
*    rmag        - eci position vector magnitude  km
*    vmag        - eci velocity vector magnitude  km/sec
*    latgc       - geocentric latitude            rad
*    lon         - longitude                      rad
*    fpa         - sat flight path angle          rad
*    az          - sat flight path az             rad
*    ttt         - julian centuries of tt         centuries
*    jdut1       - julian date of ut1             days from 4713 bc
*    lod         - excess length of day           sec
*    xp          - polar motion coefficient       arc sec
*    yp          - polar motion coefficient       arc sec
*    terms       - number of terms for ast calculation 0,2
*
*  outputs       :
*    r           - eci position vector            km
*    v           - eci velocity vector            km/s
*
*  locals        :
*    fpav        - sat flight path anglefrom vert rad
*
*  coupling      :
*    none        -
*
*  references    :
*    vallado       2007, 117
*    chobotov            67
*
* ----------------------------------------------------------------------------

      SUBROUTINE flt2rv      ( rmag, vmag, latgc, lon, fpa, az,
     &                         ttt, jdut1, lod, xp, yp, terms, r, v )
        IMPLICIT NONE
        REAL*8 R(3), V(3), rmag, vmag, latgc, lon, fpa, az, ttt, jdut1,
     &         lod, xp, yp
        INTEGER terms

* -----------------------------  Locals  ------------------------------
        REAL*8 recef(3), vecef(3), rtasc, decl, temp, fpav

        INCLUDE 'astmath.cmn'
        INCLUDE 'astrolib.cmn'

        ! --------------------  Implementation   ----------------------
        ! -------- form position vector
        recef(1) = rmag*dcos(latgc)*dcos(lon)
        recef(2) = rmag*dcos(latgc)*dsin(lon)
        recef(3) = rmag*dsin(latgc)

        ! -------- convert r to eci
        vecef(1) = 0.0D0
        vecef(2) = 0.0D0
        vecef(3) = 0.0D0
c        CALL ECI_ECEF( r,v, 'FROM', rECEF,vECEF,
c     &                 TTT, JDUT1, LOD, xp, yp, terms )

        ! ------------- calculate rtasc and decl ------------------
        temp= dsqrt( r(1)*r(1) + r(2)*r(2) )

* v needs to be defined herexxxxxxxxx
        if ( temp .lt. small ) THEN
            rtasc= datan2( v(2) , v(1) )
          else
            rtasc= datan2( r(2) , r(1) )
          ENDIF
        decl= asin( r(3)/rmag )

        ! -------- form velocity vector
        fpav = halfpi - fpa
        v(1)= vmag*( dcos(rtasc)*(-dcos(az)*dsin(fpav)*dsin(decl) +
     &                dcos(fpav)*dcos(decl)) - dsin(az)*dsin(fpav)*
     &                dsin(rtasc) )
        v(2)= vmag*( dsin(rtasc)*(-dcos(az)*dsin(fpav)*dsin(decl) +
     &                dcos(fpav)*dcos(decl)) + dsin(az)*dsin(fpav)*
     &                dcos(rtasc) )
        v(3)= vmag*( dcos(az)*dcos(decl)*dsin(fpav) + dcos(fpav)*
     &                dsin(decl) )

      RETURN
      END

* ----------------------------------------------------------------------------
*
*                           function rv2flt.m
*
*  this function transforms a position and velocity vector into the flight
*    elements - latgc, lon, fpa, az, position and velocity magnitude.
*
*  author        : david vallado                  719-573-2600   17 jun 2002
*
*  revisions
*    vallado     - add terms for ast calculation                 30 sep 2002
*
*  inputs          description                    range / units
*    r           - eci position vector            km
*    v           - eci velocity vector            km/s
*    ttt         - julian centuries of tt         centuries
*    jdut1       - julian date of ut1             days from 4713 bc
*    lod         - excess length of day           sec
*    xp          - polar motion coefficient       arc sec
*    yp          - polar motion coefficient       arc sec
*    terms       - number of terms for ast calculation 0,2
*
*  outputs       :
*    rmag        - eci position vector magnitude  km
*    vmag        - eci velocity vector magnitude  km/sec
*    latgc       - geocentric latitude            rad
*    lon         - longitude                      rad
*    fpa         - sat flight path angle          rad
*    az          - sat flight path az             rad
*
*  locals        :
*    fpav        - sat flight path anglefrom vert rad
*
*  coupling      :
*    none        -
*
*  references    :
*    vallado       2007, 117
*
* ----------------------------------------------------------------------------

      SUBROUTINE rv2flt      ( R, V, ttt, jdut1, lod, xp, yp, terms,
     &                         rmag, vmag, latgc, lon, fpa, az)
        IMPLICIT NONE
        REAL*8 R(3), V(3), rmag, vmag, latgc, lon, fpa, az, ttt, jdut1,
     &         lod, xp, yp, Mag, hmag, DOT
        INTEGER terms
        EXTERNAL DOT, MAG

* -----------------------------  Locals  ------------------------------
        REAL*8 recef(3), vecef(3), temp, fpav, rdotv, hcrossr(3), h(3)

        INCLUDE 'astmath.cmn'
        INCLUDE 'astrolib.cmn'

        ! --------------------  Implementation   ----------------------
        rmag = mag(r)
        vmag = mag(v)

        ! -------- convert r to ecef for lat/lon calculation
        CALL ECI_ECEF( r,v, 'TOO ', rECEF,vECEF,
     &                 TTT, JDUT1, LOD, xp, yp, terms )

        ! ----------------- find longitude value  ----------------- uses ecef
        temp = Dsqrt( recef(1)*recef(1) + recef(2)*recef(2) )
        if ( temp .lt. small ) THEN
            lon= Datan2( vecef(2), vecef(1) )
          else
            lon= Datan2( recef(2), recef(1) )
          ENDIF

*        latgc = Datan2( recef(3) , dsqrt(recef(1)**2 + recef(2)**2) )
        latgc = Dasin( recef(3) / rmag )

        CALL cross(r, v, h)
        hmag = mag(h)
        rdotv= dot(r, v)
        fpav= Datan2(hmag, rdotv)
        fpa = halfpi - fpav

        CALL cross(h, r, hcrossr)

        az = Datan2( r(1)*hcrossr(2) - r(2)*hcrossr(1),
     &               hcrossr(3)*rmag )

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           function eq2rv
*
*  this function finds the classical orbital elements given the equinoctial
*    elements.
*
*  author        : david vallado                  719-573-2600    9 jun 2002
*
*  revisions
*    vallado     - fix elliptical equatorial orbits case         19 oct 2002
*
*  inputs          description                    range / units
*    af          -
*    ag          -
*    n           - mean motion                    rad
*    meanlon     - mean longitude                 rad
*    chi         -
*    psi         -
*
*  outputs       :
*    r           - eci position vector            km
*    v           - eci velocity vector            km/s
*
*  locals        :
*    temp        - temporary variable
*    p           - semilatus rectum               km
*    a           - semimajor axis                 km
*    ecc         - eccentricity
*    incl        - inclination                    0.0  to pi rad
*    omega       - longitude of ascending node    0.0  to 2pi rad
*    argp        - argument of perigee            0.0  to 2pi rad
*    nu          - true anomaly                   0.0  to 2pi rad
*    m           - mean anomaly                   0.0  to 2pi rad
*    arglat      - argument of latitude      (ci) 0.0  to 2pi rad
*    truelon     - true longitude            (ce) 0.0  to 2pi rad
*    lonper      - longitude of periapsis    (ee) 0.0  to 2pi rad
*
*  coupling      :
*
*  references    :
*    vallado       2007, 116
*    chobotov            30
*
* ------------------------------------------------------------------------------

      SUBROUTINE eq2rv       ( af, ag, meanlon, n, chi, psi, r, v)
        IMPLICIT NONE
        REAL*8 R(3), V(3), af, ag, meanlon, n, chi, psi

* -----------------------------  Locals  ------------------------------
        REAL*8 p, a, ecc, incl, omega, argp, nu, m,
     &         arglat, truelon, lonper, e0
        INTEGER fr

        INCLUDE 'astmath.cmn'
        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
        arglat  = 999999.1D0
        lonper  = 999999.1D0
        truelon = 999999.1D0

        a = (mu/n**2)**(1.0D0/3.0D0)

        ecc = dsqrt (af**2 + ag**2)

        p = a * (1.0D0 - ecc*ecc)

        incl = 2.0D0 * datan( dsqrt(chi**2 + psi**2) )

        ! -------- setup retrograde factor ----------------------------
        fr = 1
        ! -------- set this so it only affects i = 180 deg orbits!! ---
        if (Dabs(incl-pi) .lt. small) THEN
            fr = -1
          ENDIF

        omega = Datan2( chi, psi)

        argp = Datan2( fr*ag, af ) - datan2( chi, psi )

        if ( ecc .lt. small ) THEN
            ! ----------------  circular equatorial  ------------------
            if ( (incl.lt.small) .or. ( dabs(incl-pi).lt.small) ) THEN
                argp = 0.0D0
                omega= 0.0D0
              else
                ! --------------  circular inclined  ------------------
                argp= 0.0D0
              ENDIF
          else
            ! ---------------  elliptical equatorial  -----------------
            if ( (incl.lt.small) .or. (dabs(incl-pi).lt.small) ) THEN
                omega= 0.0D0
              ENDIF
          ENDIF

        m = meanlon - omega - argp
        m = dmod(m+twopi, twopi)

        CALL newtonm ( ecc, m, e0, nu )

        ! ----------  fix for elliptical equatorial orbits ------------
        if ( ecc .lt. small ) THEN
            ! ----------------  circular equatorial  ------------------
            if ((incl.lt.small) .or. ( dabs(incl-pi).lt. small )) THEN
                argp    = undefined
                omega   = undefined
                truelon = nu
              else
                ! --------------  circular inclined  ------------------
                argp  = undefined
                arglat= nu
              ENDIF
            nu   = undefined
          else
            ! ---------------  elliptical equatorial  -----------------
            if ( ( incl.lt.small) .or. (dabs(incl-pi).lt.small) ) THEN
                lonper = argp
                argp    = undefined
                omega   = undefined
              ENDIF
          ENDIF

        ! -------- now convert back to position and velocity vectors
        CALL coe2rv(p, ecc, incl, omega, argp, nu, arglat, truelon,
     &       lonper, r, v)

      RETURN
      END

* ----------------------------------------------------------------------------
*
*                           function rv2eq.m
*
*  this function transforms a position and velocity vector into the flight
*    elements - latgc, lon, fpa, az, position and velocity magnitude.
*
*  author        : david vallado                  719-573-2600    7 jun 2002
*
*  revisions
*    vallado     - fix special orbit types (ee, hyper)           15 mar 2003
*
*  inputs          description                    range / units
*    r           - eci position vector            km
*    v           - eci velocity vector            km/s
*
*  outputs       :
*    af          -
*    ag          -
*    meanlon     - mean longitude                 rad
*    n           - mean motion                    rad/s
*    chi         -
*    psi         -
*
*  locals        :
*    none        -
*
*  coupling      :
*    none        -
*
*  references    :
*    vallado       2007, 116
*    chobotov            30
*
* ----------------------------------------------------------------------------

      SUBROUTINE rv2eq       ( R, V, af, ag, meanlon, n, chi, psi)
        IMPLICIT NONE
        REAL*8 R(3), V(3), af, ag, meanlon, n, chi, psi, DCOT
        EXTERNAL DCOT

* -----------------------------  Locals  ------------------------------
        REAL*8 p, a, ecc, incl, omega, argp, nu, m,
     &         arglat, truelon, lonper
        INTEGER fr

        INCLUDE 'astmath.cmn'
        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
        ! -------- convert to classical elements ----------------------
        CALL rv2coe( r, v,
     &       p, a, ecc, incl, omega, argp, nu, m, arglat, truelon,
     &       lonper )

        ! -------- setup retrograde factor ----------------------------
        fr = 1
        ! -------- set this so it only affects i = 180 deg orbits!! ---
        if (dabs(incl-pi) .lt. small) THEN
            fr = -1
          ENDIF

        if ( ecc .lt. small ) THEN
            ! ----------------  circular equatorial  ------------------
            if ((incl.lt.small) .or. ( dabs(incl-pi).lt. small )) THEN
                argp = 0.0D0
                omega= 0.0D0
              else
                ! --------------  circular inclined  ------------------
                argp= 0.0D0
              ENDIF
          else
            ! ---------------  elliptical equatorial  -----------------
            if ( ( incl.lt.small) .or. (dabs(incl-pi).lt.small) ) THEN
                argp = lonper
                omega= 0.0D0
              ENDIF
          ENDIF

        af = ecc * dcos(fr*omega + argp)
        ag = ecc * dsin(fr*omega + argp)

        if (fr .gt. 0  ) THEN
            chi = dtan(incl*0.5D0) * dsin(omega)
            psi = dtan(incl*0.5D0) * dcos(omega)
          else
            chi = dcot(incl*0.5D0) * dsin(omega)
            psi = dcot(incl*0.5D0) * dcos(omega)
          ENDIF

c        IF (DABS(ecc-1.0D0).le.small) THEN
c            n  = 2.0D0 * dsqrt(mu/(p*p*p))
c          ELSE
            IF (a.gt.0.0D0) THEN
                n  = dsqrt(mu/(a*a*a))
              ELSE
                n  = dsqrt(-mu/(a*a*a))
              ENDIF
c          ENDIF

        meanlon = fr*omega + argp + m
        meanlon = dmod(meanlon, twopi)

      RETURN
      END

* ----------------------------------------------------------------------------
*
*                           function adbar2rv.m
*
*  this function transforms the adbarv elements (rtasc, decl, fpa, azimuth,
*    position and velocity magnitude) into eci position and velocity vectors.
*
*  author        : david vallado                  719-573-2600    9 jun 2002
*
*  revisions
*                -
*
*  inputs          description                    range / units
*    rmag        - eci position vector magnitude  km
*    vmag        - eci velocity vector magnitude  km/sec
*    rtasc       - right ascension of sateillite  rad
*    decl        - declination of satellite       rad
*    fpav        - sat flight path angle from vertrad
*    az          - sat flight path azimuth        rad
*
*  outputs       :
*    r           - eci position vector            km
*    v           - eci velocity vector            km/s
*
*  locals        :
*    none        -
*
*  coupling      :
*    none        -
*
*  references    :
*    vallado       2007, 117
*    chobotov            70
*
* ----------------------------------------------------------------------------

      SUBROUTINE adbar2rv    ( rmag, vmag, rtasc, decl, fpav, az, r, v )
        IMPLICIT NONE
        REAL*8 R(3), V(3), rmag, vmag, rtasc, decl, fpav, az

        ! --------------------  Implementation   ----------------------
        ! -------- form position vector
        r(1)= rmag*dcos(decl)*dcos(rtasc)
        r(2)= rmag*dcos(decl)*dsin(rtasc)
        r(3)= rmag*dsin(decl)

        ! -------- form velocity vector
        v(1)= vmag*( dcos(rtasc)*(-dcos(az)*dsin(fpav)*dsin(decl) +
     &                dcos(fpav)*dcos(decl)) - dsin(az)*dsin(fpav)*
     &                dsin(rtasc) )
        v(2)= vmag*( dsin(rtasc)*(-dcos(az)*dsin(fpav)*dsin(decl) +
     &                dcos(fpav)*dcos(decl)) + dsin(az)*dsin(fpav)*
     &                dcos(rtasc) )
        v(3)= vmag*( dcos(az)*dcos(decl)*dsin(fpav) + dcos(fpav)*
     &               dsin(decl) )

      RETURN
      END

* ----------------------------------------------------------------------------
*
*                           function rv2adbar.m
*
*  this function transforms a position and velocity vector into the adbarv
*    elements - rtasc, decl, fpav, azimuth, position and velocity magnitude.
*
*  author        : david vallado                  719-573-2600    9 jun 2002
*
*  revisions
*                -
*
*  inputs          description                    range / units
*    r           - eci position vector            km
*    v           - eci velocity vector            km/s
*
*  outputs       :
*    rmag        - eci position vector magnitude  km
*    vmag        - eci velocity vector magnitude  km/sec
*    rtasc       - right ascension of sateillite  rad
*    decl        - declination of satellite       rad
*    fpav        - sat flight path angle from vertrad
*    az          - sat flight path azimuth        rad
*
*  locals        :
*    none        -
*
*  coupling      :
*    none        -
*
*  references    :
*    vallado       2007, 117
*    chobotov            70
*
* ----------------------------------------------------------------------------

      SUBROUTINE rv2adbar    ( R, V, rmag, vmag, rtasc, decl, fpav, az)
        IMPLICIT NONE
        REAL*8 R(3), V(3), rmag, vmag, rtasc, decl, fpav, az, Dot, MAG
        EXTERNAL DOT, MAG

* -----------------------------  Locals  ------------------------------
        REAL*8 temp, temp1, h(3), hcrossr(3), rdotv, hmag

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        rmag = mag(r)
        vmag = mag(v)

        ! ---------------- calculate rtasc and decl -------------------
        temp= dsqrt( r(1)*r(1) + r(2)*r(2) )
        if ( temp .lt. small ) THEN
            temp1= dsqrt( v(1)*v(1) + v(2)*v(2) )
            if ( dabs(temp1) .gt. small ) THEN
                rtasc= datan2( v(2) , v(1) )
              else
                rtasc= 0.0D0
              ENDIF
          else
            rtasc= datan2( r(2), r(1) )
          ENDIF
        decl= asin( r(3)/rmag )

        CALL cross(r, v, h)
        hmag = mag(h)
        rdotv= dot(r, v)
        fpav = datan2(hmag, rdotv)

        CALL cross(h, r, hcrossr)
        az = datan2( r(1)*hcrossr(2) - r(2)*hcrossr(1),
     &               hcrossr(3)*rmag )

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           function rv2rsw
*
*  this function converts position and velocity vectors into radial, along-
*    track, and cross-track coordinates. note that sometimes the middle vector
*    is called in-track.
*
*  author        : david vallado                  719-573-2600    9 jun 2002
*
*  revisions
*                -
*
*  inputs          description                    range / units
*    r           - position vector                km
*    v           - velocity vector                km/s
*
*  outputs       :
*    rrsw        - position vector                km
*    vrsw        - velocity vector                km/s
*    transmat    - transformation matrix 
*
*  locals        :
*    tempvec     - temporary vector
*    rvec,svec,wvec - direction cosines
*
*  coupling      :
*
*
*  references    :
*    vallado       2007, 163
*
* ------------------------------------------------------------------------------

      SUBROUTINE rv2rsw      ( r, v, rrsw, vrsw, transmat )
        IMPLICIT NONE
        REAL*8 r(3), v(3), rrsw(3), vrsw(3), transmat(3,3)

* -----------------------------  Locals  ------------------------------
        REAL*8 rvec(3), svec(3), wvec(3), tempvec(3)

        ! --------------------  Implementation   ----------------------
        ! in order to work correctly each of the components must be
        ! unit vectors
        ! radial component
        CALL norm( r, rvec)

        ! cross-track component
        CALL cross(r, v, tempvec)
        CALL norm( tempvec,wvec )

        ! along-track component
        CALL cross(wvec, rvec, tempvec)
        CALL norm( tempvec, svec )

        ! assemble transformation matrix from to rsw frame (individual
        !  components arranged in row vectors)
        transmat(1,1) = rvec(1)
        transmat(1,2) = rvec(2)
        transmat(1,3) = rvec(3)
        transmat(2,1) = svec(1)
        transmat(2,2) = svec(2)
        transmat(2,3) = svec(3)
        transmat(3,1) = wvec(1)
        transmat(3,2) = wvec(2)
        transmat(3,3) = wvec(3)

        CALL MatVecMult  ( transmat, r, 3,3,3,1, rrsw )
        CALL MatVecMult  ( transmat, v, 3,3,3,1, vrsw )

*   alt approach
*       rrsw(1) = mag(r)
*       rrsw(2) = 0.0D0
*       rrsw(3) = 0.0D0
*       vrsw(1) = dot(r, v)/rrsw(1)
*       vrsw(2) = dsqrt(v(1)**2 + v(2)**2 + v(3)**2 - vrsw(1)**2)
*       vrsw(3) = 0.0D0

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           function rv2ntw
*
*  this function converts position and velocity vectors into in-radial,
*    velocity, and cross-track coordinates. note that sometimes the first
*    vector is called along-radial.
*
*  author        : david vallado                  719-573-2600    5 jul 2002
*
*  revisions
*                -
*
*  inputs          description                    range / units
*    r           - position vector                km
*    v           - velocity vector                km/s
*
*  outputs       :
*    rntw        - position vector                km
*    vntw        - velocity vector                km/s
*    transmat    - transformation matrix
*
*  locals        :
*    tempvec     - temporary vector
*    tvec,nvec,wvec - direction cosines
*
*  coupling      :
*
*
*  references    :
*    vallado       2007, 164
*
* ------------------------------------------------------------------------------

      SUBROUTINE rv2ntw      ( R, V, rntw, vntw, transmat )
        IMPLICIT NONE
        REAL*8 R(3), V(3), rntw(3), vntw(3), transmat(3,3)

* -----------------------------  Locals  ------------------------------
        REAL*8 nvec(3), tvec(3), wvec(3), tempvec(3)

        ! --------------------  Implementation   ----------------------
        ! in order to work correctly each of the components must be
        ! unit vectors
        ! in-velocity component
        CALL norm( v, nvec)

        ! cross-track component
        CALL cross(r, v, tempvec)
        CALL norm( tempvec,wvec )

        ! along-radial component
        CALL cross(nvec, wvec, tempvec)
        CALL norm( tempvec, tvec )

        ! assemble transformation matrix from to ntw frame (individual
        !  components arranged in row vectors)
        transmat(1,1) = tvec(1)
        transmat(1,2) = tvec(2)
        transmat(1,3) = tvec(3)
        transmat(2,1) = nvec(1)
        transmat(2,2) = nvec(2)
        transmat(2,3) = nvec(3)
        transmat(3,1) = wvec(1)
        transmat(3,2) = wvec(2)
        transmat(3,3) = wvec(3)

        CALL MatVecMult  ( transmat, r, 3,3,3,1, rntw )
        CALL MatVecMult  ( transmat, v, 3,3,3,1, vntw )

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE FINDC2C3
*
*  this subroutine calculates the C2 and C3 functions for use in the Universal
*    Variable calcuLation of z.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    ZNew        - Z variable                     rad2
*
*  Outputs       :
*    C2New       - C2 FUNCTION value
*    C3New       - C3 FUNCTION value
*
*  Locals        :
*    SqrtZ       - Square root of ZNew
*
*  Coupling      :
*    SINH        - Hyperbolic Sine
*    COSH        - Hyperbolic Cosine
*
*  References    :
*    Vallado       2007, 71, Alg 1
*
* ------------------------------------------------------------------------------

      SUBROUTINE FINDC2C3    ( ZNew, C2New, C3New )
        IMPLICIT NONE
        REAL*8 ZNew, C2New, C3New
* -----------------------------  Locals  ------------------------------
        REAL*8 SqrtZ

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        IF ( ZNew .gt. Small ) THEN
            SqrtZ = DSQRT( ZNew )
            C2New = (1.0D0-DCOS( SqrtZ )) / ZNew
            C3New = (SqrtZ-DSIN( SqrtZ )) / ( SqrtZ**3 )
          ELSE
            IF ( ZNew .lt. -Small ) THEN
                SqrtZ = DSQRT( -ZNew )
                C2New = (1.0D0-COSH( SqrtZ )) / ZNew 
                C3New = (SINH( SqrtZ ) - SqrtZ) / ( SqrtZ**3 )
              ELSE
                C2New = 0.5D0
                C3New = 1.0D0/6.0D0 
              ENDIF 
          ENDIF 
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE NEWTONE
*
*  this subroutine solves Keplers equation when the Eccentric, paraboic, .or.
*    Hyperbolic anomalies are known. The Mean anomaly and true anomaly are
*    calculated.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Ecc         - Eccentricity                   0.0D0 to
*    E0          - Eccentric Anomaly              -2Pi to 2Pi rad
*
*  Outputs       :
*    M           - Mean Anomaly                   0.0D0 to 2Pi rad
*    Nu          - True Anomaly                   0.0D0 to 2Pi rad
*
*  Locals        :
*    Sinv        - Sine of Nu
*    Cosv        - Cosine of Nu
*
*  Coupling      :
*    SINH        - Hyperbolic Sine
*    COSH        - Hyperbolic Cosine
*
*  References    :
*    Vallado       2007, 85, Alg 6
*
* ------------------------------------------------------------------------------

      SUBROUTINE NEWTONE     ( Ecc, E0, M, Nu )
        IMPLICIT NONE
        REAL*8 Ecc, E0, M, Nu
* -----------------------------  Locals  ------------------------------
        Real*8 Sinv, Cosv

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        ! ------------------------- Circular --------------------------
        IF ( DABS( Ecc ) .lt. Small ) THEN
            M = E0
            Nu= E0 
          ELSE

            ! ----------------------- Elliptical ----------------------
            IF ( Ecc .lt. 0.999D0 ) THEN
                M= E0 - Ecc*DSIN(E0)
                Sinv= ( DSQRT( 1.0D0-Ecc*Ecc ) * DSIN(E0) ) /
     &                ( 1.0D0-Ecc*DCOS(E0) )
                Cosv= ( DCOS(E0)-Ecc ) / ( 1.0D0 - Ecc*DCOS(E0) ) 
                Nu  = DATAN2( Sinv, Cosv )
              ELSE

                ! ---------------------- Hyperbolic  ------------------
                IF ( Ecc .gt. 1.0001D0 ) THEN
                    M= Ecc*SINH(E0) - E0
                    Sinv= ( DSQRT( Ecc*Ecc-1.0D0 ) * SINH(E0) ) /
     &                    ( 1.0D0 - Ecc*COSH(E0) )
                    Cosv= ( COSH(E0)-Ecc ) / ( 1.0D0 - Ecc*COSH(E0) ) 
                    Nu  = DATAN2( Sinv, Cosv )
                  ELSE

                    ! -------------------- Parabolic ------------------
                    M= E0 + (1.0D0/3.0D0)*E0*E0*E0
                    Nu= 2.0D0*DATAN(E0) 
                  ENDIF 
              ENDIF
          ENDIF
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE NEWTONM
*
*  this subroutine performs the Newton Rhapson iteration to find the
*    Eccentric Anomaly given the Mean anomaly.  The True Anomaly is also
*    calculated.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Ecc         - Eccentricity                   0.0D0 to
*    M           - Mean Anomaly                   -2Pi to 2Pi rad
*
*  Outputs       :
*    E0          - Eccentric Anomaly              0.0D0 to 2Pi rad
*    Nu          - True Anomaly                   0.0D0 to 2Pi rad
*
*  Locals        :
*    E1          - Eccentric Anomaly, next value  rad
*    Sinv        - Sine of Nu
*    Cosv        - Cosine of Nu
*    Ktr         - Index
*    R1r         - CUBIC roots - 1 to 3
*    R1i         - imaginary component
*    R2r         -
*    R2i         -
*    R3r         -
*    R3i         -
*    S           - Variables for parabolic solution
*    W           - Variables for parabolic solution
*
*  Coupling      :
*    CUBIC       - Solves a CUBIC polynomial
*    SINH        - Hyperbolic Sine
*    COSH        - Hyperbolic Cosine
*
*  References    :
*    Vallado       2001, 73, Alg 2, Ex 2-1
*
* ------------------------------------------------------------------------------

      SUBROUTINE NEWTONM     ( Ecc, M, E0, Nu )
        IMPLICIT NONE
        REAL*8 Ecc, M, E0, Nu
        EXTERNAL DCot
* -----------------------------  Locals  ------------------------------
        INTEGER Ktr, NumIter
        REAL*8 E1, Sinv, Cosv, R1r, R1i, R2r, R2i, R3r, R3i, DCot
        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        NumIter =    50
        ! -------------------------- Hyperbolic  ----------------------
        IF ( (Ecc-1.0D0) .gt. Small ) THEN
           ! -------------------  Initial Guess -----------------------
            IF ( Ecc .lt. 1.6D0 ) THEN
                IF ( ((M.lt.0.0D0).and.(M.gt.-Pi)).or.(M.gt.Pi) ) THEN
                    E0= M - Ecc
                  ELSE
                    E0= M + Ecc
                  ENDIF
              ELSE
                IF ( (Ecc .lt. 3.6D0) .and. (DABS(M) .gt. Pi) ) THEN
                    E0= M - DSIGN(1.0D0, M)*Ecc
                  ELSE
                    E0= M/(Ecc-1.0D0)
                  ENDIF
              ENDIF
            Ktr= 1
            E1 = E0 + ( (M-Ecc*SINH(E0)+E0) / (Ecc*COSH(E0) - 1.0D0) )
            DO WHILE ((DABS(E1-E0).gt.Small ) .and. ( Ktr.le.NumIter ))
                E0= E1
                E1= E0 + ( ( M - Ecc*SINH(E0) + E0 ) /
     &                     ( Ecc*COSH(E0) - 1.0D0 ) )
                Ktr = Ktr + 1
              ENDDO
            ! ----------------  Find True Anomaly  --------------------
            Sinv= -( DSQRT( Ecc*Ecc-1.0D0 ) * SINH(E1) ) /
     &             ( 1.0D0 - Ecc*COSH(E1) )
            Cosv= ( COSH(E1) - Ecc ) / ( 1.0D0 - Ecc*COSH(E1) )
            Nu  = DATAN2( Sinv, Cosv )
          ELSE
            ! --------------------- Parabolic -------------------------
            IF ( DABS( Ecc-1.0D0 ) .lt. Small ) THEN
                CALL CUBIC( 1.0D0/3.0D0, 0.0D0, 1.0D0, -M, R1r, R1i,
     &                      R2r, R2i, R3r, R3i )
                E0= R1r
*                 S = 0.5D0 * (HalfPi - DATAN( 1.5D0*M ) )
*                 W = DATAN( DTAN( S )**(1.0D0/3.0D0) )
*                 E0= 2.0D0*DCOT(2.0D0*W)
                Ktr= 1
                Nu = 2.0D0 * DATAN(E0)
              ELSE
                ! -------------------- Elliptical ----------------------
                IF ( Ecc .gt. Small ) THEN
                    ! -----------  Initial Guess -------------
                    IF ( ((M .lt. 0.0D0) .and. (M .gt. -Pi)) .or.
     &                   (M .gt. Pi) ) THEN
                        E0= M - Ecc
                      ELSE
                        E0= M + Ecc
                      ENDIF
                    Ktr= 1
                    E1 = E0 + ( M - E0 + Ecc*DSIN(E0) ) /
     &                        ( 1.0D0 - Ecc*DCOS(E0) )
                    DO WHILE (( DABS(E1-E0) .gt. Small ) .and.
     &                       ( Ktr .le. NumIter ))
                        Ktr = Ktr + 1
                        E0= E1
                        E1= E0 + ( M - E0 + Ecc*DSIN(E0) ) /
     &                           ( 1.0D0 - Ecc*DCOS(E0) )
                      ENDDO
                    ! -------------  Find True Anomaly  ---------------
                    Sinv= ( DSQRT( 1.0D0-Ecc*Ecc ) * DSIN(E1) ) /
     &                    ( 1.0D0-Ecc*DCOS(E1) )
                    Cosv= ( DCOS(E1)-Ecc ) / ( 1.0D0 - Ecc*DCOS(E1) )
                    Nu  = DATAN2( Sinv, Cosv )
                  ELSE
                    ! -------------------- Circular -------------------
                    Ktr= 0
                    Nu= M
                    E0= M
                  ENDIF
              ENDIF
          ENDIF
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE NEWTONNU
*
*  this subroutine solves Keplers equation when the true anomaly is known.
*    The Mean and Eccentric, parabolic, or hyperbolic anomaly is also found.
*    The parabolic limit at 168ø is arbitrary. The hyperbolic anomaly is also
*    limited. The hyperbolic sine is used because it's not double valued.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Ecc         - Eccentricity                   0.0D0 to
*    Nu          - True Anomaly                   -2Pi to 2Pi rad
*
*  Outputs       :
*    E0          - Eccentric Anomaly              0.0D0 to 2Pi rad       153.02D0ø
*    M           - Mean Anomaly                   0.0D0 to 2Pi rad       151.7425D0ø
*
*  Locals        :
*    E1          - Eccentric Anomaly, next value  rad
*    SinE        - Sine of E
*    CosE        - Cosine of E
*    Ktr         - Index
*
*  Coupling      :
*    ASINH     - Arc hyperbolic sine
*    SINH        - Hyperbolic Sine
*
*  References    :
*    Vallado       2007, 85, Alg 5
*
* ------------------------------------------------------------------------------

      SUBROUTINE NEWTONNU    ( Ecc, Nu, E0, M )
        IMPLICIT NONE
        REAL*8 Ecc, Nu, E0, M
        EXTERNAL ASINH
* -----------------------------  Locals  ------------------------------
        REAL*8 SinE, CosE, ASINH

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        E0= 999999.9D0
        M = 999999.9D0
        ! --------------------------- Circular ------------------------
        IF ( DABS( Ecc ) .lt. 0.000001D0 ) THEN
            M = Nu
            E0= Nu 
          ELSE
            ! ---------------------- Elliptical -----------------------
            IF ( Ecc .lt. 0.999D0 ) THEN
                SinE= ( DSQRT( 1.0D0-Ecc*Ecc ) * DSIN(Nu) ) /
     &                ( 1.0D0+Ecc*DCOS(Nu) )
                CosE= ( Ecc + DCOS(Nu) ) / ( 1.0D0 + Ecc*DCOS(Nu) )
                E0  = DATAN2( SinE, CosE )
                M   = E0 - Ecc*DSIN(E0) 
              ELSE
                ! -------------------- Hyperbolic  --------------------
                IF ( Ecc .gt. 1.0001D0 ) THEN
                    IF ( ((Ecc .gt. 1.0D0) .and. (DABS(Nu)+0.00001D0
     &                     .lt. Pi-DACOS(1.0D0/Ecc)) ) ) THEN
                        SinE= ( DSQRT( Ecc*Ecc-1.0D0 ) * DSIN(Nu) ) /
     &                        ( 1.0D0 + Ecc*DCOS(Nu) )
                        E0  = ASINH( SinE )
                        M   = Ecc*DSINH(E0) - E0
                      ENDIF 
                  ELSE
                    ! ----------------- Parabolic ---------------------
                    IF ( DABS(Nu) .lt. 168.0D0/57.29578D0 ) THEN
                        E0= DTAN( Nu*0.5D0 )
                        M = E0 + (E0*E0*E0)/3.0D0 
                      ENDIF
                  ENDIF
              ENDIF
          ENDIF

        IF ( Ecc .lt. 1.0D0 ) THEN
            M = DMOD( M, 2.0D0*Pi )
            IF ( M .lt. 0.0D0 ) THEN
                M= M + 2.0D0*Pi 
              ENDIF
            E0 = DMOD( E0, 2.0D0*Pi )
          ENDIF 
      RETURN
      END  ! end newtonnu

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE KEPLER
*
*  this subroutine solves Keplers problem for orbit determination and returns a
*    future Geocentric Equatorial (IJK) position and velocity vector.  The
*    solution SUBROUTINE uses Universal variables.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*                                                               15 Jan 1993
*                     29 Nov 95 - needs repeat until - else it can't for 10.0D0TUs
*  Inputs          Description                    Range / Units
*    Ro          - IJK Position vector - initial  km
*    Vo          - IJK Velocity vector - initial  km / s
*    dtsec        - Length of time to propagate    s
*
*  OutPuts       :
*    R           - IJK Position vector            km
*    V           - IJK Velocity vector            km / s
*    Error       - Error flag                     'ok', ...
*
*  Locals        :
*    F           - f expression
*    G           - g expression
*    FDot        - f DOT expression
*    GDot        - g DOT expression
*    XOld        - Old Universal Variable X
*    XOldSqrd    - XOld squared
*    XNew        - New Universal Variable X
*    XNewSqrd    - XNew squared
*    ZNew        - New value of z
*    C2New       - C2(psi) FUNCTION
*    C3New       - C3(psi) FUNCTION
*    dtsec        - change in time                 s
*    TimeNew     - New time                       s
*    RDotV       - Result of Ro DOT Vo
*    A           - Semi .or. axis                 km
*    Alpha       - Reciprocol  1/a
*    SME         - Specific Mech Energy           km2 / s2
*    Period      - Time period for satellite      s
*    S           - Variable for parabolic case
*    W           - Variable for parabolic case
*    H           - Angular momentum vector
*    Temp        - Temporary REAL*8 value
*    i           - index
*
*  Coupling      :
*    MAG         - Magnitude of a vector
*    DOT         - DOT product of two vectors
*    COT         - Cotangent FUNCTION
*    FINDC2C3    - Find C2 and C3 functions
*    CROSS       - CROSS product of two vectors
*
*  References    :
*    Vallado       2007, 101, Alg 8, Ex 2-4
*
* ------------------------------------------------------------------------------

      SUBROUTINE KEPLER      ( Ro, Vo, dtsec, R, V, Error )
        IMPLICIT NONE
        REAL*8 Ro(3), Vo(3), dtsec, R(3), V(3), Dot
        CHARACTER*12 Error
        EXTERNAL Dot, DCot, Mag
* -----------------------------  Locals  ------------------------------
        INTEGER Ktr, i, NumIter
        REAL*8 H(3), F, G, FDot, GDot, Rval, XOld, XOldSqrd, XNew,
     &      XNewSqrd, ZNew, p, C2New, C3New, DtNew, RDotV, A,
     &      DCot, mag, Alpha, SME, Period, S, W, Temp,
     &      magro, magvo, magh, magr, tmp

        INCLUDE 'astmath.cmn'
        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
        NumIter    =    35
        ! --------------------  Initialize values   -------------------
        ktr = 0
        xOld= 0.0D0
        ZNew= 0.0D0
        Error= 'ok' 

        IF ( DABS( dtsec ) .gt. Small ) THEN
            magro = MAG( Ro )
            magvo = MAG( Vo )
            RDotV= DOT( Ro, Vo )

            ! -------------  Find SME, Alpha, and A  ------------------
            SME= ( magvo*magvo*0.5D0 ) - ( mu/magro )
            Alpha= -SME*2.0D0 / mu

            IF ( DABS( SME ) .gt. Small ) THEN
                A= -mu / ( 2.0D0*SME )
              ELSE
                A= Infinite
              ENDIF
            IF ( DABS( Alpha ) .lt. Small ) THEN ! Parabola
                Alpha= 0.0D0
              ENDIF

            ! ------------   Setup initial guess for x  ---------------
            ! -----------------  Circle and Ellipse -------------------
            IF ( Alpha .ge. Small ) THEN
                Period= TwoPi * DSQRT( DABS(A)**3/mu )
                ! ------- Next IF needed for 2body multi-rev ----------
                IF ( DABS( dtsec ) .gt. DABS( Period ) ) THEN
                    dtsec= DMOD( dtsec, Period )
                  ENDIF
                IF ( DABS(Alpha-1.0D0) .gt. Small ) THEN
                     XOld = DSQRT(mu) * dtsec * Alpha
                  ELSE
                     ! 1st guess can't be too close. ie a circle, r=a
                     XOld= DSQRT(mu) * dtsec*Alpha*0.97D0
                  ENDIF
              ELSE
                ! --------------------  Parabola  ---------------------
                IF ( DABS( Alpha ) .lt. Small ) THEN
                    CALL CROSS( ro, vo, h )
                    magh = MAG(h)
                    p= magh*magh/mu
                    S= 0.5D0 * (HalfPi - DATAN( 3.0D0*DSQRT( mu/
     &                         (p*p*p) )* dtsec ) )
                    W= DATAN( DTAN( S )**(1.0D0/3.0D0) )
                    XOld = DSQRT(p) * ( 2.0D0*DCOT(2.0D0*W) )
                    Alpha= 0.0D0 
                  ELSE
                    ! ------------------  Hyperbola  ------------------
                    Temp= -2.0D0*mu*dtsec /
     &                   ( A*( RDotV + DSIGN(1.0D0, dtsec)*DSQRT(-mu*A)*
     &                   (1.0D0-magro*Alpha) ) )
                    XOld= DSIGN(1.0D0, dtsec) * DSQRT(-A) *DLOG(Temp)
                  ENDIF

              ENDIF
     
            tmp = 1.0 / DSQRT(mu)  
            Ktr= 1
            DtNew = -10.0D0
            DO WHILE ( (DABS(DtNew-DSQRT(mu)*dtsec).ge.Small).and.
     &                 (Ktr.lt.NumIter) )
                XOldSqrd = XOld*XOld 
                ZNew     = XOldSqrd * Alpha

                ! ------------- Find C2 and C3 functions --------------
                CALL FINDC2C3( ZNew, C2New, C3New )

                ! ------- Use a Newton iteration for New values -------
                DtNew= XOldSqrd*XOld*C3New + (RDotV*tmp)*XOldSqrd
     &                   *C2New + magro*XOld*( 1.0D0 - ZNew*C3New )
                Rval = XOldSqrd*C2New + (RDotV*tmp)*XOld*(1.0D0-
     &                   ZNew*C3New) + magro*( 1.0D0 - ZNew*C2New )

                ! ------------- Calculate New value for x -------------
                XNew = XOld + ( DSQRT(mu)*dtsec - DtNew ) / Rval

				! ----- check if the univ param goes negative. if so, use bissection
				IF (XNew < 0.0D0 .and. dtsec > 0.0D0) THEN
					XNew = XOld*0.5D0
    			ENDIF

                  Ktr = Ktr + 1
                XOld = XNew 
              ENDDO

            IF ( Ktr .ge. NumIter ) THEN
                Error= 'KNotConv'
c               Write(*,*) ' Not converged in ', NumIter:2,' iterations '
                DO i= 1 , 3
                    V(i)= 0.0D0
                    R(i)= V(i)
                  ENDDO
              ELSE
                ! --- Find position and velocity vectors at New time --
                XNewSqrd = XNew*XNew
                F = 1.0D0 - ( XNewSqrd*C2New / magro )
                G = dtsec - XNewSqrd*XNew*C3New*tmp
                DO i= 1 , 3
                    R(i)= F*Ro(i) + G*Vo(i)
                  ENDDO
                magr = MAG( R )
                GDot = 1.0D0 - ( XNewSqrd*C2New / magr )
                FDot = ( DSQRT(mu)*XNew / ( magro*magr ) ) *
     &                 ( ZNew*C3New-1.0D0 )
                DO i= 1 , 3
                    V(i)= FDot*Ro(i) + GDot*Vo(i)
                  ENDDO
                Temp= F*GDot - FDot*G 
                IF ( DABS(Temp-1.0D0) .gt. 0.00001D0 ) THEN
                    Error= 'FandG'
                  ENDIF 
              ENDIF  ! IF (
          ELSE
            ! ----------- Set vectors to incoming since 0 time --------
            DO i=1, 3
                r(i)= ro(i)
                v(i)= vo(i)
              ENDDO
          ENDIF 

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE FINDTOF
*
*  this subroutine finds the time of flight given the initial position vectors,
*    Semi-parameter, and the sine and cosine values for the change in true
*    anomaly.  The result uses p-iteration theory to analytically find the result.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Ro          - Interceptor position vector    km
*    R           - TARGET position vector         km
*    p           - Semiparameter                  km
*
*  Outputs       :
*    Tof         - Time for transfer              s
*
*  Locals        :
*    SinDNu      - Sine of change in Nu           rad
*    CosDNu      - Cosine of change in Nu         rad
*    DeltaE      -
*    DeltaH      -
*    k           -
*    l           -
*    m           -
*    a           -
*    f           -
*    g           -
*    FDot        -
*    SinDeltaE   - Sine value
*    CosDeltaE   - Cosine value
*    RcrossR     - CROSS product of two positions
*
*  Coupling      :
*    CROSS       - CROSS product of two vectors
*    SINH        - Hyperbolic Sine
*    ACOSH     - Arc hyperbolic cosine
*
*  References    :
*    Vallado       2007, 134, Alg 11
*
* ------------------------------------------------------------------------------

      SUBROUTINE FINDTOF     ( Ro, R, p, Tof )
        IMPLICIT NONE
        REAL*8 Ro(3), R(3), p, Tof, MAG
        EXTERNAL Dot, ACOSH, MAG
* -----------------------------  Locals  ------------------------------
        REAL*8 RCrossR(3), CosDNu, SinDNu, Small, c , s, alpha, DeltaE,
     &    DeltaH, DNu, k, l, m, a, f, g, FDot, SinDeltaE, CosDeltaE,
     &    Dot, ACOSH, magro, magr, magrcrossr

        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
        Small = 0.00001D0 

        magro = MAG(ro)
        magr = MAG(r)
        CosDNu= DOT(Ro, R)/(magro*magr)
        CALL CROSS( Ro, R, RCrossR )
        magrcrossr = MAG(rcrossr)
        SinDNu= magRCrossR/(magro*magr)

        k= magro * magr*( 1.0D0-CosDNu )
        l= magro + magr
        m= magro * magr*( 1.0D0+CosDNu )
        a= (m*k*p) / ((2.0D0*m-l*l)*p*p + 2.0D0*k*l*p - k*k) 

        ! ------  Use F and G series to find Velocity Vectors  --------
        F = 1.0D0 - ( magr/p )*(1.0D0-CosDNu)
        G = magro*magr*SinDNu/DSQRT(mu*p)
        Alpha= 1.0D0/a 

        IF ( alpha .gt. Small ) THEN
            ! ------------------------ Elliptical ---------------------
            DNu  = DATAN2( SinDNu, CosDNu )
            FDot = DSQRT(mu/p) * DTAN(DNu*0.5D0)*
     &              ( ((1.0D0-CosDNu)/p)-(1.0D0/magro)-(1.0D0/magr) )
            COSDeltaE= 1.0D0-(magro/a)*(1.0D0-f)
            SinDeltaE= (-magro*magr*FDot)/DSQRT(mu*a)
            DeltaE   = DATAN2( SinDeltaE, CosDeltaE )
            Tof      = G + DSQRT(a*a*a/mu)*(DeltaE-SinDeltaE)
          ELSE
            ! ------------------------ Hyperbolic ---------------------
            IF ( alpha .lt. -Small ) THEN
                DeltaH = ACOSH( 1.0D0-(magro/a)*(1.0D0-F) )
                Tof    = G + DSQRT(-a*a*a/mu)*(SINH(DeltaH)-DeltaH)
              ELSE
                ! -------------------- Parabolic ----------------------
                DNu= DATAN2( SinDNu, CosDNu )
                c  = DSQRT( magr*magr+magro*magro -
     &                     2.0D0*magr*magro*DCOS(DNu) )
                s  = (magro+magr+c ) * 0.5D0
                Tof= ( 2.0D0/3.0D0 ) * DSQRT(s*s*s*0.5D0/mu) *
     &                     (1.0D0 -  ((s-c)/s)**1.5D0 )
              ENDIF 
          ENDIF

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE ijk2ll
*
*  These SUBROUTINEs convert a Geocentric Equatorial (IJK) position vector into
*    latitude and longitude.  Geodetic and Geocentric latitude are found.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Recef       - IJK position vector            km
*
*  OutPuts       :
*    Latgc       - Geocentric Latitude            -Pi to Pi rad
*    Latgd       - Geodetic Latitude              -Pi to Pi rad
*    Lon         - Longitude (WEST -)             -2Pi to 2Pi rad
*    Hellp       - Height above the ellipsoid     km
*
*  Locals        :
*  Escobal:
*    Rc          - Range of SITE wrt earth center km
*    Height      - Height above earth wrt SITE    km
*    Alpha       - ANGLE from Iaxis to point, LST rad
*    OldDelta    - Previous value of DeltaLat     rad
*    DeltaLat    - Diff between Delta and
*                  Geocentric lat                 rad
*    Delta       - Declination ANGLE of R in IJK  rad
*    RSqrd       - Magnitude of r squared         km2
*    SinTemp     - Sine of Temp                   rad
*    c           -
*
*  Almanac:
*    Temp        - Diff between Geocentric/
*                  Geodetic lat                   rad
*    SinTemp     - Sine of Temp                   rad
*    OldDelta    - Previous value of DeltaLat     rad
*    RtAsc       - Right ascension                rad
*    Decl        - Declination                    rad
*    i           - index
*
*  Borkowski:
*
*
*  Coupling      :
*    MAG         - Magnitude of a vector
*    gc2gd    - Converts between geocentric and geodetic latitude
*
*  References    :
*    Vallado       2007, 179, Alg 12 and Alg 13, Ex 3-3
*
* ------------------------------------------------------------------------------

      SUBROUTINE ijk2llA ( Recef, Latgc, Latgd, Lon, Hellp )
        IMPLICIT NONE
        REAL*8 Recef(3), JD, Latgc, Latgd, Lon, Hellp
        EXTERNAL MAG
* -----------------------------  Locals  ------------------------------
        INTEGER i
        REAL*8 RtAsc, OldDelta, c, Decl,
     &         Temp, SinTemp, MAG, magr

        INCLUDE 'astmath.cmn'
        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
        magr = MAG( Recef )

        ! ----------------- Find Longitude value  ---------------------
        Temp = DSQRT( Recef(1)*Recef(1) + Recef(2)*Recef(2) )
        IF ( DABS( Temp ) .lt. Small ) THEN
            RtAsc= DSIGN(1.0D0, Recef(3))*Pi*0.5D0
          ELSE
            RtAsc= DATAN2( Recef(2) / Temp , Recef(1) / Temp )
          ENDIF
        Lon  = RtAsc
        IF ( DABS(Lon) .ge. Pi ) THEN ! Mod it ?
            IF ( Lon .lt. 0.0D0 ) THEN
                Lon= TwoPi + Lon
              ELSE
                Lon= Lon - TwoPi 
              ENDIF
          ENDIF
        Decl = DASIN( Recef(3) / magr )
        Latgd= Decl 

        ! ------------- Iterate to find Geodetic Latitude -------------
        i= 1 
        OldDelta = Latgd + 10.0D0

        DO WHILE ((DABS(OldDelta-Latgd).ge.Small).and.(i.lt.10))
            OldDelta= Latgd 
            SinTemp = DSIN( Latgd ) 
            c       = rekm / (DSQRT( 1.0D0-EESqrd*SinTemp*SinTemp ))
            Latgd= DATAN( (Recef(3)+c*EESqrd*SinTemp)/Temp )
            i = i + 1
          ENDDO

        Hellp   = (Temp/DCOS(Latgd)) - c

        CALL gc2gd( Latgc, 'FROM', Latgd )

      RETURN
      END

      SUBROUTINE ijk2llE ( Recef, Latgc, Latgd, Lon, Hellp )
        IMPLICIT NONE
        REAL*8 Recef(3), Latgc, Latgd, Lon, Hellp
        EXTERNAL MAG

* -----------------------------  Locals  ------------------------------
        INTEGER i
        Real*8 rsite, DeltaLat, RSqrd
        REAL*8 RtAsc, OldDelta, Decl,
     &         Temp, SinTemp, OneMinusE2, MAG, magr
        CHARACTER Show

        INCLUDE 'astmath.cmn'
        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
        Show = 'N'

       ! -------------------  Initialize values   --------------------
        magr = MAG( Recef )
        OneMinuse2 = 1.0D0 - EeSqrd

       ! ---------------- Find Longitude value  ----------------------
        Temp = DSQRT( Recef(1)*Recef(1) + Recef(2)*Recef(2) )
        IF ( DABS( Temp ) .lt. Small ) THEN
            RtAsc= DSIGN(1.0D0, Recef(3))*Pi*0.5D0
          ELSE
            RtAsc= DATAN2( Recef(2) / Temp , Recef(1) / Temp )
          ENdif
        Lon  = RtAsc 

        IF ( DABS(Lon) .ge. Pi ) THEN
            IF ( Lon .lt. 0.0D0 ) THEN
                Lon= TwoPi + Lon
              ELSE
                Lon= Lon - TwoPi 
              ENDIF
          ENDIF
       ! -------------- Set up initial latitude value  ---------------  
        Decl    = DASIN( Recef(3) / magr )
        Latgc= Decl 
        DeltaLat= 100.0D0 
        RSqrd   = magr**2

       ! ---- Iterate to find Geocentric .and. Geodetic Latitude  -----  
        i= 1 
        DO WHILE ( ( DABS( OldDelta - DeltaLat ) .ge. Small ) .and.
     &             ( i .lt. 10 ))
            OldDelta = DeltaLat 
            rsite    = DSQRT( OneMinuse2 / (1.0D0 -
     &                 EeSqrd*(DCOS(Latgc))**2 ) )
            Latgd = DATAN( DTAN(Latgc) / OneMinuse2 ) 
            Temp     = Latgd-Latgc 
            SinTemp  = DSIN( Temp ) 
            Hellp    = DSQRT( RSqrd - rsite*rsite*SinTemp*SinTemp ) -
     &                 rsite*DCOS(Temp)
            DeltaLat = DASIN( Hellp*SinTemp / magr )
            Latgc = Decl - DeltaLat 
            i = i + 1
            IF ( Show .eq. 'Y' ) THEN
                write(*, *) 'E loops gc gd ', Latgc*57.29578D0,
     &                     Latgd*57.29578D0
              ENDIF
          ENDDO

        IF ( i .ge. 10 ) THEN
            Write(*, *) 'ijk2ll did NOT converge '
          ENDIF

      RETURN
      END

* ------------------------------- Borkowski method  --------------------------
      SUBROUTINE ijk2llB ( Recef , Latgc, Latgd, Lon, Hellp )
        IMPLICIT NONE
        REAL*8 Recef(3), Latgc, Latgd, Lon, Hellp
        EXTERNAL MAG
        REAL*8 MAG

* -----------------------------  Locals  ------------------------------
        Real*8 a, b, RtAsc, sqrtp, third, e, f, p, q, d, nu, g, t,
     &         aTemp, Temp

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------

        ! ---------------- Find Longitude value  ----------------------
        Temp = DSQRT( Recef(1)*Recef(1) + Recef(2)*Recef(2) )
        IF ( DABS( Temp ) .lt. Small ) THEN
            RtAsc= DSIGN(1.0D0, Recef(3))*Pi*0.5D0
          ELSE
            RtAsc= DATAN2( Recef(2) / Temp , Recef(1) / Temp )
          ENDIF
        Lon  = RtAsc 
        IF ( DABS(Lon) .ge. Pi ) THEN
            IF ( Lon .lt. 0.0D0 ) THEN
                Lon= TwoPi + Lon
              ELSE
                Lon= Lon - TwoPi
              ENDIF
          ENDIF

        a= 1.0D0 
        b= DSIGN(1.0D0, Recef(3))*6356.75160056D0/6378.137D0
       ! -------------- Set up initial latitude value  ---------------  
        aTemp= 1.0D0/(a*Temp) 
        e= (b*Recef(3)-a*a+b*b)*atemp
        f= (b*Recef(3)+a*a-b*b)*atemp
        third= 1.0D0/3.0D0 
        p= 4.0D0*Third*(e*f + 1.0D0 ) 
        q= 2.0D0*(e*e - f*f) 
        d= p*p*p + q*q 

        IF ( d .gt. 0.0D0 ) THEN
            nu= (DSQRT(d)-q)**third - (DSQRT(d)+q)**third
          ELSE
            SqrtP= DSQRT(-p)
            nu= 2.0D0*SqrtP*DCOS( third*DACOS(q/(p*SqrtP)) ) 
          ENDIF 
        g= 0.5D0*(DSQRT(e*e + nu) + e) 
        t= DSQRT(g*g + (f-nu*g)/(2.0D0*g-e)) - g 

        Latgd= DATAN(a*(1.0D0-t*t)/(2.0D0*b*t)) 
        hellp= (temp-a*t)*DCOS( Latgd) + (Recef(3)-b)*DSIN(Latgd)

        CALL gc2gd( Latgc, 'FROM', Latgd )
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE GC2GD
*
*  this subroutine converts from Geodetic to Geocentric Latitude for positions
*    on the surface of the Earth.  Notice that (1-f) squared = 1-eSqrd.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Latgd       - Geodetic Latitude              -Pi to Pi rad
*    Direction   - Which set of vars to output    FROM  TOO
*
*  Outputs       :
*    Latgc       - Geocentric Latitude            -Pi to Pi rad
*
*  Locals        :
*    None.
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 148, Eq 3-11
*
* ------------------------------------------------------------------------------

      SUBROUTINE gc2gd       ( Latgc, Direction, Latgd )
        IMPLICIT NONE
        REAL*8 Latgc, Latgd
        CHARACTER*4 Direction

        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
        IF ( Direction .eq. 'FROM' ) THEN
            Latgc= DATAN( (1.0D0 - EESqrd)*DTAN(Latgd) )
          ELSE
            Latgd= DATAN( DTAN(Latgc)/(1.0D0 - EESqrd) )
          ENDIF
      RETURN
      END


* ------------------------------------------------------------------------------
*
*                           SUBROUTINE SIGHT
*
*  this subroutine takes the position vectors of two satellites and determines
*    If there is line-of-SIGHT between the two satellites.  An oblate Earth
*    with radius of 1 ER is assumed.  The process forms the equation of
*    a line between the two vectors.  Differentiating and setting to zero finds
*    the minimum value, and when plugged back into the original line equation,
*    gives the minimum distance.  The parameter TMin is allowed to range from
*    0.0D0 to 1.0D0.  Scale the K-component to account for oblate Earth because it's
*    the only qunatity that changes.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R1          - Position vector of the 1st sat km
*    R2          - Position vector of the 2nd sat km
*    WhichKind   - Spherical .or. Ellipsoidal Earth 'S', 'E'*default
*
*  Outputs       :
*    LOS         - Line of SIGHT                  'YES', 'NO '
*
*  Locals        :
*    TR1         - Scaled R1 vector               km
*    TR2         - Scaled R2 vector               km
*    ADotB       - DOT product of a DOT b
*    TMin        - Minimum value of t from a to b
*    DistSqrd    - minute Distance squared to Earth  km
*    ASqrd       - Magnitude of A squared
*    BSqrd       - Magnitude of B squared
*
*  Coupling:
*    DOT         - DOT product of two vectors
*
*  References    :
*    Vallado       2007, 310, Alg 35, Ex 5-3
* ------------------------------------------------------------------------------

      SUBROUTINE SIGHT       ( R1, R2, WhichKind, LOS )
        IMPLICIT NONE
        REAL*8 R1(3), R2(3)
        INTEGER i
        CHARACTER WhichKind
        CHARACTER*3 LOS
        EXTERNAL Dot, MAG
* -----------------------------  Locals  ------------------------------
        REAL*8 TR1(3), TR2(3), ADotB, TMin, DistSqrd, ASqrd,
     &         BSqrd, Temp, Dot, Mag, magtr1, magtr2

        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
        DO i=1 , 3
            TR1(i)= R1(i)
            TR2(i)= R2(i)
          ENDDO
        ! --------------------- Scale z component ---------------------
        IF ( WhichKind .eq. 'E' ) THEN
            Temp= 1.0D0/DSQRT(1.0D0-EESqrd)
          ELSE
            Temp= 1.0D0
          ENDIF
        TR1(3)= TR1(3)*Temp
        TR2(3)= TR2(3)*Temp
        magtr1 = MAG(tr1)
        magtr2 = MAG(tr2)
        BSqrd= magTR2**2
        ASqrd= magTR1**2
        ADotB= DOT( TR1, TR2 )
        ! ---------------------- Find TMin ----------------------------
        DistSqrd= 0.0D0
        IF ( DABS(ASqrd + BSqrd - 2.0D0*ADotB) .lt. 0.0001D0 ) THEN
            TMin= 0.0D0
          ELSE
            TMin = ( ASqrd - ADotB ) / ( ASqrd + BSqrd - 2.0D0*ADotB )
          ENDIF
        ! ----------------------- Check LOS ---------------------------
        IF ( (TMin .lt. 0.0D0) .or. (TMin .gt. 1.0D0) ) THEN
            LOS= 'YES'
          ELSE
            DistSqrd= ( (1.0D0-TMin)*ASqrd + ADotB*TMin ) / rekm**2
            IF ( DistSqrd .gt. rekm ) THEN
                LOS= 'YES'
              ELSE
                LOS= 'NO '
              ENDIF
          ENDIF
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE SUN
*
*  this subroutine calculates the Geocentric Equatorial position vector
*    the SUN given the Julian Date.  This is the low precision formula .and.
*    is valid for years from 1950 to 2050.  Accuaracy of apparent coordinates
*    is 0.01D0 degrees.  Notice many of the calcuLations are performed in
*    degrees, and are not changed until Later.  This is due to the fact that
*    the Almanac uses degrees exclusively in their formuLations.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    JD          - Julian Date                    days from 4713 BC
*
*  Outputs       :
*    RSun        - IJK Position vector of the SUN AU
*    RtAsc       - Right Ascension                rad
*    Decl        - Declination                    rad
*
*  Locals        :
*    MeanLong    - Mean Longitude
*    MeanAnomaly - Mean anomaly
*    EclpLong    - Ecliptic Longitude
*    Obliquity   - Mean Obliquity of the Ecliptic
*    TUT1        - Julian Centuries of UT1 from
*                  Jan 1, 2000 12h
*    TTDB        - Julian Centuries of TDB from
*                  Jan 1, 2000 12h
*    Hr          - Hours                          0 .. 24              10
*    minute         - MiNutes                        0 .. 59              15
*    SEC         - Seconds                        0.0D0 .. 59.99D0         30.00D0
*    Temp        - Temporary variable
*    deg         - Degrees
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 281, Alg 29, Ex 5-1
*
* ------------------------------------------------------------------------------

      SUBROUTINE SUN         ( JD, RSun, RtAsc, Decl )
        IMPLICIT NONE
        REAL*8 JD, RSun(3), RtAsc, Decl
* -----------------------------  Locals  ------------------------------
        REAL*8 MeanLong, MeanAnomaly,
     &        EclpLong, Obliquity, TUT1, TTDB, magrsun

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        ! -------------------  Initialize values   --------------------
        TUT1= ( JD - 2451545.0D0 )/ 36525.0D0 

        MeanLong= 280.4606184D0 + 36000.77005361D0*TUT1 
        MeanLong= DMOD( MeanLong, 360.0D0 )  !deg

        TTDB= TUT1 
        MeanAnomaly= 357.5277233D0 + 35999.05034D0*TTDB 
        MeanAnomaly= DMOD( MeanAnomaly*Deg2Rad, TwoPi )  !rad
        IF ( MeanAnomaly .lt. 0.0D0 ) THEN
            MeanAnomaly= TwoPi + MeanAnomaly
          ENDIF

        EclpLong= MeanLong + 1.914666471D0*DSIN(MeanAnomaly)
     &              + 0.019994643D0*DSIN(2.0D0*MeanAnomaly) !deg

        Obliquity= 23.439291D0 - 0.0130042D0*TTDB  !deg

        MeanLong = MeanLong*Deg2Rad 
        IF ( MeanLong .lt. 0.0D0 ) THEN
            MeanLong= TwoPi + MeanLong
          ENDIF
        EclpLong = EclpLong *Deg2Rad 
        Obliquity= Obliquity *Deg2Rad 

        ! ------- Find magnitude of SUN vector, ) THEN components -----
        magRSun= 1.000140612D0 - 0.016708617D0*DCOS( MeanAnomaly )
     &                         - 0.000139589D0*DCOS( 2.0D0*MeanAnomaly )    ! in AU's

        RSun(1)= magRSun*DCOS( EclpLong )
        RSun(2)= magRSun*DCOS(Obliquity)*DSIN(EclpLong)
        RSun(3)= magRSun*DSIN(Obliquity)*DSIN(EclpLong)

        RtAsc= DATAN( DCOS(Obliquity)*DTAN(EclpLong) )
        ! --- Check that RtAsc is in the same quadrant as EclpLong ----
        IF ( EclpLong .lt. 0.0D0 ) THEN
            EclpLong= EclpLong + TwoPi    ! make sure it's in 0 to 2pi range
          ENDIF
        IF ( DABS( EclpLong-RtAsc ) .gt. Pi*0.5D0 ) THEN
            RtAsc= RtAsc + 0.5D0*Pi*DNINT( (EclpLong-RtAsc)/(0.5D0*Pi))
          ENDIF
        Decl = DASIN( DSIN(Obliquity)*DSIN(EclpLong) )

      RETURN
      END

*  References    :
*    Vallado       2007, 311, Eq 5-9

      SUBROUTINE SunIll   ( JD, Lat, Lon, SunIllum, SunAz, SunEl )
        IMPLICIT NONE
        REAL*8 SunIllum, Lat, Lon, JD, SunAz, SunEl

* -----------------------------  Locals  ------------------------------
        Real*8 RSun(3)
        Real*8 lst, gst, x, LHA, sinv, cosv
        Real*8 l0, l1, l2, l3, sRtAsc, sdecl

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        CALL SUN( JD, RSun, sRtAsc, sDecl ) ! AU's needed for Sun ill

        CALL LSTIME( Lon, JD, LST, GST )

        LHA = LST - sRtAsc

        SunEl  = DASIN( dsin(sDecl)*dsin(Lat) +
     &           dcos(sDecl)*dcos(Lat)*dcos(LHA) )

        Sinv= -dsin(LHA)*dcos(sDecl)*dcos(Lat)/(dcos(SunEl)*dcos(Lat))
        Cosv= ( dsin(sDecl)-dsin(SunEl)*dsin(Lat) )/
     &        ( dcos(SunEl)*dcos(Lat) )
        SunAz  = DATAN2( Sinv, Cosv )

        SunEl= SunEl/Deg2Rad

        IF (SunEl .gt. -18.01D0) THEN
            x= SunEl/90.0D0

         IF (SunEl .ge. 20) THEN
             l0=  3.74
             l1=  3.97
             l2= -4.07
             l3=  1.47
           ELSEIF ((SunEl .ge. 5.0).and.(SunEl .lt. 20.0)) THEN
                 l0=   3.05
                 l1=  13.28
                 l2= -45.98
                 l3=  64.33
               ELSEIF ((SunEl .ge. -0.8).and.(SunEl .lt. 5.0)) THEN
                     l0=    2.88
                     l1=   22.26
                     l2= -207.64
                     l3= 1034.30
                   ELSEIF ((SunEl .ge. -5.0).and.(SunEl .lt. -0.8))
     &                THEN
                         l0=    2.88
                         l1=   21.81
                         l2= -258.11
                         l3= -858.36
                       ELSEIF ((SunEl .ge. -12.0).and.(SunEl .lt.
     &                          -5.0)) THEN
                             l0=    2.70
                             l1=   12.17
                             l2= -431.69
                             l3=-1899.83
                           ELSEIF ((SunEl .ge. -18.0).and.(SunEl .lt.
     &                              -12.0)) THEN
                                 l0=   13.84
                                 l1=  262.72
                                 l2= 1447.42
                                 l3= 2797.93
                               ELSE
                                 l0= 0.0
                                 l1= 0.0
                                 l2= 0.0
                                 l3= 0.0
                               ENDIF

         l1= l0 + l1*x + l2*x*x + l3*x*x*x
         SunIllum= 10.0** l1
         IF ((SunIllum .lt. -1D+36).or.(SunIllum .gt. 999.999D0)) THEN
             SunIllum= 0.0
           ENDIF
       ELSE
         SunIllum= 0.0D0
       ENDIF

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE MOON
*
*  this subroutine calculates the Geocentric Equatorial (IJK) position vector
*    for the MOON given the Julian Date.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    JD          - Julian Date                    days from 4713 BC
*
*  Outputs       :
*    RMoon       - IJK Position vector of MOON    ER
*    RtAsc       - Right Ascension                rad
*    Decl        - Declination                    rad
*
*  Locals        :
*    EclpLong    - Ecliptic Longitude
*    EclpLat     - Eclpitic Latitude
*    HzParal     - Horizontal Parallax
*    l           - Geocentric Direction Cosines
*    m           -             "     "
*    n           -             "     "
*    TTDB        - Julian Centuries of TDB from
*                  Jan 1, 2000 12h
*    Hr          - Hours                          0 .. 24
*    minute         - MiNutes                        0 .. 59
*    SEC         - Seconds                        0.0D0 .. 59.99D0
*    deg         - Degrees
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 290, Alg 31, Ex 5-3
*
* ------------------------------------------------------------------------------

      SUBROUTINE MOON        ( JD, RMoon, RtAsc, Decl )
        IMPLICIT NONE
        REAL*8 JD, RMoon(3), RtAsc, Decl

* -----------------------------  Locals  ------------------------------
        REAL*8 TTDB, l, m, n, Obliquity, magrmoon, EclpLong, EclpLat,
     &         HzParal

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        TTDB = ( JD - 2451545.0D0 ) / 36525.0D0

        EclpLong= 218.32D0 + 481267.883D0*TTDB
     &              + 6.29D0*DSIN( (134.9D0+477198.85D0*TTDB)*Deg2Rad )
     &              - 1.27D0*DSIN( (259.2D0-413335.38D0*TTDB)*Deg2Rad )
     &              + 0.66D0*DSIN( (235.7D0+890534.23D0*TTDB)*Deg2Rad )
     &              + 0.21D0*DSIN( (269.9D0+954397.70D0*TTDB)*Deg2Rad )
     &              - 0.19D0*DSIN( (357.5D0+ 35999.05D0*TTDB)*Deg2Rad )
     &              - 0.11D0*DSIN( (186.6D0+966404.05D0*TTDB)*Deg2Rad )  ! Deg

        EclpLat =   5.13D0*DSIN( ( 93.3D0+483202.03D0*TTDB)*Deg2Rad )
     &              + 0.28D0*DSIN( (228.2D0+960400.87D0*TTDB)*Deg2Rad )
     &              - 0.28D0*DSIN( (318.3D0+  6003.18D0*TTDB)*Deg2Rad )
     &              - 0.17D0*DSIN( (217.6D0-407332.20D0*TTDB)*Deg2Rad )  ! Deg

        HzParal =  0.9508D0 + 0.0518D0*DCOS( (134.9D0+477198.85D0*TTDB)
     &              *Deg2Rad )
     &            + 0.0095D0*DCOS( (259.2D0-413335.38D0*TTDB)*Deg2Rad )
     &            + 0.0078D0*DCOS( (235.7D0+890534.23D0*TTDB)*Deg2Rad )
     &            + 0.0028D0*DCOS( (269.9D0+954397.70D0*TTDB)*Deg2Rad )  ! Deg

        EclpLong = DMOD( EclpLong*Deg2Rad, TwoPi )
        EclpLat  = DMOD( EclpLat*Deg2Rad, TwoPi )
        HzParal  = DMOD( HzParal*Deg2Rad, TwoPi )

        Obliquity= 23.439291D0 - 0.0130042D0*TTDB  !deg
        Obliquity= Obliquity *Deg2Rad

        ! ------------ Find the geocentric direction cosines ----------
        l= DCOS( EclpLat ) * DCOS( EclpLong )
        m= DCOS(Obliquity)*DCOS(EclpLat)*DSIN(EclpLong)
     &       - DSIN(Obliquity)*DSIN(EclpLat)
        n= DSIN(Obliquity)*DCOS(EclpLat)*DSIN(EclpLong)
     &       + DCOS(Obliquity)*DSIN(EclpLat)

        ! ------------- Calculate MOON position vector ----------------
        magRMoon= 1.0D0/DSIN( HzParal )
        RMoon(1)= magRMoon*l
        RMoon(2)= magRMoon*m
        RMoon(3)= magRMoon*n

        ! -------------- Find Rt Ascension and Declination ------------
        RtAsc= DATAN2( m, l )
        Decl = DASIN( n )

      RETURN
      END

*  References    :
*    Vallado       2007, 311, Eq 5-9

      SUBROUTINE MoonIll     ( MoonEl, f, MoonIllum )
        IMPLICIT NONE
        REAL*8 MoonEl, f, MoonIllum

* -----------------------------  Locals  ------------------------------
        REAL*8  x, l0, l1, l2, l3

        ! --------------------  Implementation   ----------------------
        x= MoonEl/90.0D0
c        g= 1.0

       IF (MoonEl .ge. 20) THEN
         l0= -1.95
         l1=  4.06
         l2= -4.24
         l3=  1.56
        ELSEIF ((MoonEl .ge. 5.0).and.(MoonEl .lt. 20.0)) THEN
             l0=  -2.58
             l1=  12.58
             l2= -42.58
             l3=  59.06
           ELSEIF ((MoonEl .gt. -0.8).and.(MoonEl .lt. 5.0)) THEN
                 l0=   -2.79
                 l1=   24.27
                 l2= -252.95
                 l3= 1321.29
               ELSE
                 l0= 0.0
                 l1= 0.0
                 l2= 0.0
                 l3= 0.0
                 f= 0.0
c                 g= 0.0
               ENDIF

       l1= l0 + l1*x + l2*x*x + l3*x*x*x
       l2= (-0.00868D0*f - 2.2D-9*f*f*f*f)

*       HzParal =   0.9508 + 0.0518*dcos( (134.9+477198.85*TTDB)*Deg2Rad )
*                + 0.0095*dcos( (259.2-413335.38*TTDB)*Deg2Rad )
*                + 0.0078*dcos( (235.7+890534.23*TTDB)*Deg2Rad )
*                + 0.0028*dcos( (269.9+954397.70*TTDB)*Deg2Rad )   { Deg }
*       HzParal  = REALMOD( HzParal*Deg2Rad, TwoPi )
*       l3= (2.0* POWER(10.0, (HzParal*rad / 0.951))*g ) { use g to eliminate neg el passes }

        MoonIllum= 10.0D0 ** ( l1 + l2 )
        IF ((MoonIllum .lt. -1d+36).or.(MoonIllum .gt. 0.999D0)) THEN
            MoonIllum= 0.0D0
          ENDIF
       RETURN
       END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE LIGHT
*
*  this subroutine determines If a spacecraft is sunlit .or. in the dark at a
*    particular time.  An oblate Earth and cylindrical shadow is assumed.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R           - Position vector of sat         km
*    JD          - Julian Date at desired time    Days from 4713 BC
*    WhichKind   - Spherical .or. Ellipsoidal Earth 'S', 'E'*default
*
*  OutPuts       :
*    Vis         - Visibility Flag                'YES', 'NO '
*
*  Locals        :
*    RtAsc       - Suns Right ascension           rad
*    Decl        - Suns Declination               rad
*    RSun        - SUN vector                     AU
*    AUER        - Conversion from AU to ER
*
*  Coupling      :
*    SUN         - Position vector of SUN
*    LNCOM1      - Multiple a vector by a constant
*    SIGHT       - Does Line-of-SIGHT exist beteen vectors
*
*  References    :
*    Vallado       2007, 310, Alg 35, Ex 5-6
*
* ------------------------------------------------------------------------------

      SUBROUTINE LIGHT       ( R, JD, WhichKind, LIT )
        IMPLICIT NONE
        REAL*8 R(3), JD
        CHARACTER WhichKind, Lit(3)
* -----------------------------  Locals  ------------------------------
        REAL*8 RSun(3), RtAsc, Decl

        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------

        CALL SUN( JD, RSun, RtAsc, Decl )
        CALL LNCOM1( AUER, RSun, RSun )

        ! ------------ Is the satellite in the shadow? ----------------
        CALL SIGHT( RSun, R, WhichKind, Lit )

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE CHECKHITEARTH
*
*  this subroutine checks to see If the trajectory hits the earth during the
*    transfer.  The first check determines If the satellite is initially
*    heading towards perigee, and finally heading away from perigee.  IF (
*    this is the case, the radius of perigee is calculated.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    RInt        - Initial Position vector of Int km
*    V1t         - Initial Velocity vector of trnskm/TU
*    RTgt        - Initial Position vector of Tgt km
*    V2t         - Final Velocity vector of trns  km/TU
*
*  Outputs       :
*    HitEarth    - Is Earth was impacted          'Y' 'N'
*
*  Locals        :
*    SME         - Specific mechanical energy
*    rp          - Radius of Perigee              km
*    TransA      - Semi-.or. axis of transfer     km
*    TransE      - Eccentricity of transfer
*    TransP      - Semi-paramater of transfer     km
*    HBar        - Angular momentum vector of
*                  transfer orbit
*
*  Coupling      :
*    DOT         - DOT product of vectors
*    MAG         - Magnitude of a vector
*    CROSS       - CALL CROSS product of vectors
*
*  References    :
*    Vallado       2007, 500, Alg 60
*
* ------------------------------------------------------------------------------

      SUBROUTINE CHECKHITEARTH ( Rint, V1t, Rtgt, V2t, HitEarth )
        IMPLICIT NONE
        REAL*8 RInt(3), V1t(3), RTgt(3), V2t(3)
        CHARACTER HitEarth
        EXTERNAL DOT, MAG
* -----------------------------  Locals  ------------------------------
        REAL*8 HBar(3), SME, rp, TransP, TransA, TransE, Dot, Mag,
     &         magrint, magv1t, maghbar

        ! --------------------  Implementation   ----------------------
        HitEarth= 'N'

        ! ---------- Find If trajectory intersects Earth --------------
        IF ((DOT(Rint,V1t).lt.0.0D0).and.(DOT(RTgt,V2t).gt.0.0D0)) THEN

            ! ---------------  Find H N and E vectors   ---------------
            CALL CROSS( RInt, V1t, HBar )
            magrint = MAG( rint )
            magv1t  = MAG( v1t )
            maghbar = MAG( HBar )

            IF ( maghbar .gt. 0.00001D0 ) THEN
                ! ---------  Find a e and semi-Latus rectum   ---------
                SME    = magV1t**2*0.5D0 - ( 1.0D0/magRInt )
                TransP = maghbar*maghbar
                TransE = 1.0D0
                IF ( DABS( SME ) .gt. 0.00001D0 ) THEN
                    TransA= -1.0D0 / (2.0D0*SME)
                    TransE= DSQRT( (TransA - TransP)/TransA )
                    rp= TransA*(1.0D0-TransE) 
                  ELSE
                    rp= TransP*0.5D0   ! Parabola
                  ENDIF

                IF ( DABS( rp ) .lt. 1.0D0 ) THEN
                    HitEarth= 'Y' 
                  ENDIF
              ELSE
                Write(*, *) 'The orbit does not exist '
              ENDIF
          ENDIF

      RETURN
      END





* -----------------------------------------------------------------------------
*
*                           SUBROUTINE SUNRISESET
*
*  this subroutine finds the Universal time for Sunrise and Sunset given the
*    day and SITE location.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    JD          - Julian Date                    days from 4713 BC
*    Latgd       - SITE latitude (SOUTH -)        -65ø to 65ø rad
*    Lon         - SITE longitude (WEST -)        -2Pi to 2Pi rad
*    WhichKind   - Character for which rise/set   'S' 'C' 'N' 'A'
*
*  OutPuts       :
*    UTSunRise   - Universal time of sunrise      hrs
*    UTSunSet    - Universal time of sunset       hrs
*    Error       - Error Parameter
*
*  Locals        :
*    SunAngle    - ANGLE between the SUN vector
*                  and a point on the Earth     rad
*    JDTemp      - Julian date for sunrise/set    days from 4713 BC
*    UTTemp      - Temporary UT time              days
*    TUT1        - Julian Centuries from the
*                  Jan 1, 2000 12 h epoch (UT1)
*    Ra          - Right ascension                rad
*    Decl        - Declination                    rad
*    MeanLonSun  -                                rad
*    MeanAnomalySun                               rad
*    LonEcliptic - Longitude of the ecliptic      rad
*    Obliquity   - Obliquity of the ecliptic      rad
*    GST         - for 0 h UTC of each day        rad
*    LHA         - Local hour ANGLE               rad
*    Year        - Year                           1900 .. 2100
*    Mon         - Month                          1 .. 12
*    Day         - Day                            1 .. 28,29,30,31
*    Hr          - Hour                           0 .. 23
*    minute         - Minute                         0 .. 59
*    Sec         - Second                         0.0D0 .. 59.999D0
*    Opt         - Idx to for rise and set calc    1,2
*
*  Coupling      :
*    INVJDay- Finds the Year day mon hr minute Sec from the Julian Date
*    JDay   - Finds the Julian date given Year, mon day, hr, minute, Sec
*
*  References    :
*    Vallado       2007, 283, Alg 30, Ex 5-2
*
* -----------------------------------------------------------------------------

      SUBROUTINE SUNRISESET  ( JD, Latgd,Lon, WhichKind, UTSunRise,
     &                         UTSunSet, Error )
        IMPLICIT NONE
        REAL*8 JD, Latgd, Lon, UTSunRise, UTSunSet
        CHARACTER WhichKind
        CHARACTER*12 Error
* ----------------------------  Locals  -------------------------------
        INTEGER Opt, Year, Month, Day, Hr, minute
        REAL*8 MeanAnomalySun, JDF, 
     &         JDTemp, UTTemp, SunAngle, TUT1, Ra, Sec, MeanLonSun,
     &         LonEcliptic, Decl, Obliquity, GST, LHA

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        Error= 'ok'

        ! -------------- Make sure lon is within +- 180 deg -----------
        IF ( Lon .gt. Pi ) THEN
            Lon= Lon - 2.0D0*Pi
          ENDIF
        IF ( Lon .lt. -Pi ) THEN
            Lon= Lon + 2.0D0*Pi
          ENDIF
        IF (WhichKind .eq. 'S') THEN
            SunAngle= (90.0D0+50.0D0/60.0D0 )*Deg2Rad
          ENDIF
        IF (WhichKind .eq. 'C') THEN
            SunAngle=  96.0D0 *Deg2Rad
          ENDIF
        IF (WhichKind .eq. 'N') THEN
            SunAngle= 102.0D0 *Deg2Rad
          ENDIF
        IF (WhichKind .eq. 'A') THEN
            SunAngle= 108.0D0 *Deg2Rad
          ENDIF
        CALL INVJDay( JD, 0.0D0, Year,Month,Day,Hr,minute,Sec )
        DO Opt= 1 , 2
            IF ( Opt .eq. 1 ) THEN
                CALL JDay( Year,Month,Day, 6,0,0.0D0, JDTemp, JDF )
              ELSE
                CALL JDay( Year,Month,Day,18,0,0.0D0, JDTemp, JDF )
              ENDIF
            JDTemp= JDTemp + JDF - Lon*Rad2Deg/15.0D0/24.0D0

            TUT1 = (JDTemp - 2451545.0D0)/36525.0D0 
            MeanLonSun = 280.4606184D0 + 36000.77005361D0*TUT1
            MeanAnomalySun= 357.5277233D0+35999.05034D0*TUT1
            MeanAnomalySun= DMOD( MeanAnomalySun*Deg2Rad,TwoPi )
            IF ( MeanAnomalySun .lt. 0.0D0 ) THEN
                MeanAnomalySun= MeanAnomalySun + TwoPi
              ENDIF
            LonEcliptic= MeanLonSun + 1.914666471D0*DSIN(MeanAnomalySun)
     &                     + 0.019994643D0*DSIN(2.0D0*MeanAnomalySun)
            LonEcliptic= DMOD( LonEcliptic*Deg2Rad,TwoPi )
            IF ( LonEcliptic .lt. 0.0D0 ) THEN
                LonEcliptic= LonEcliptic + TwoPi
              ENDIF
            Obliquity= 23.439291D0 - 0.0130042D0*TUT1
            Obliquity= Obliquity *Deg2Rad
            Ra  = DATAN( DCOS(Obliquity) * DTAN(LonEcliptic) )
            Decl= DASIN( DSIN(Obliquity) * DSIN(LonEcliptic) )
            IF ( Ra .lt. 0.0D0 ) THEN
                Ra= Ra + TwoPi
              ENDIF
            IF ( (LonEcliptic .gt. Pi) .and. (Ra .lt. Pi) ) THEN
                Ra= Ra + Pi
              ENDIF
            IF ( (LonEcliptic .lt. Pi) .and. (Ra .gt. Pi) ) THEN
                Ra= Ra - Pi
              ENDIF
            LHA= (DCOS(SunAngle) - DSIN(Decl)*DSIN(Latgd)) /
     &           (DCOS(Decl)*DCOS(Latgd) )
            IF ( DABS(LHA) .le. 1.0D0 ) THEN
                LHA= DACOS( LHA )
              ELSE
                Error= 'Not ok'
              ENDIF
            IF ( Error .eq. 'ok' ) THEN
                IF ( Opt .eq. 1 ) THEN
                    LHA= TwoPi - LHA
                  ENDIF
                GST= 1.75336855923327D0 + 628.331970688841D0*TUT1
     &                 + 6.77071394490334D-06*TUT1*TUT1
     &                 - 4.50876723431868D-10*TUT1*TUT1*TUT1
                GST= DMOD( GST,TwoPi )
                IF ( GST .lt. 0.0D0 ) THEN
                    GST= GST + TwoPi
                  ENDIF
                UTTemp= LHA + Ra  - GST
                UTTemp= UTTemp * Rad2Deg/15.0D0
                UTTemp= DMOD( UTTemp,24.0D0 )
                UTTemp= UTTemp - Lon*Rad2Deg/15.0D0
                IF ( UTTemp .lt. 0.0D0 ) THEN
                    UTTemp= UTTemp + 24.0D0
                    Error= 'Day before'
                  ENDIF
                IF ( UTTemp .gt. 24.0D0 ) THEN
                    UTTemp= UTTemp - 24.0D0
                    Error= 'Day After'
                  ENDIF
              ELSE
                UTTemp= 99.99D0
              ENDIF
            IF ( Opt .eq. 1 ) THEN
                UTSunRise= UTTemp
              ELSE
                UTSunSet = UTTemp
              ENDIF
          ENDDO
      RETURN
      END

* -----------------------------------------------------------------------------
*
*                           SUBROUTINE MOONRISESET
*
*  this subroutine finds the Universal time for Moonrise and Moonset given the
*    day and SITE location.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    JD          - Julian Date                    days from 4713 BC
*    Latgd       - SITE latitude (SOUTH -)        -65ø to 65ø rad
*    Lon         - SITE longitude (WEST -)        -2Pi to 2Pi rad
*
*  OutPuts       :
*    UTMoonRise  - Universal time of Moonrise     hrs
*    UTMoonSet   - Universal time of Moonset      hrs
*    MoonPhaseAng- Phase angle of the Moon        deg
*    Error       - Error Parameter
*
*  Locals        :
*    MoonAngle   - ANGLE between the Moon vector
*                  and a point on the Earth       rad
*    JDTemp      - Julian date for Moonrise/set   days from 4713 BC
*    UTTemp      - Temporary UT time              days
*    TUT1        - Julian Centuries from the
*                  Jan 1, 2000 12 h epoch (UT1)
*    RtAsc       - Right ascension                rad
*    Decl        - Declination                    rad
*    MeanLonMoon -                                rad
*    MeanAnomaly -                                rad
*    EclpLong    - Longitude of the ecliptic      rad
*    Obliquity   - Obliquity of the ecliptic      rad
*    RMoon
*    RMoonRS
*    RV
*    RhoSat
*    Try
*    l, m, n     - Direction cosines
*    EclpLat
*    MoonGHA, MoonGHAn
*    DGHA, DGHAn
*    LHAn
*    LST
*    DeltaUT, DeltaUTn
*    t, tn
*    HzParal
*    LonEclSun
*    LonEclMoon
*    TTDB
*    GST         - for 0 h UTC of each day        rad
*    LHA         - Local hour ANGLE               rad
*    Year        - Year                           1900 .. 2100
*    Mon         - Month                          1 .. 12
*    Day         - Day                            1 .. 28,29,30,31
*    Hr          - Hour                           0 .. 23
*    minute         - Minute                         0 .. 59
*    Sec         - Second                         0.0D0 .. 59.999D0
*    Opt         - Idx to for rise and set calc    1,2
*
*  Coupling      :
*    INVJDay- Finds the Year day mon hr minute Sec from the Julian Date
*    JDay   - Finds the Julian date given Year, mon day, hr, minute, Sec
*
*  References    :
*    Vallado       2007, 292, Alg 32, Ex 5-4
*
* -----------------------------------------------------------------------------

      SUBROUTINE MOONRISESET ( JD,Latgd,Lon, UTMoonRise, UTMoonSet,
     &                         MoonPhaseAng, Error )
        IMPLICIT NONE
        REAL*8 JD, Latgd, Lon, UTMoonRise, UTMoonSet, MoonPhaseAng
        CHARACTER*12 Error
* ----------------------------  Locals  -------------------------------
        INTEGER Opt, i, Year, Month, Day, Hr, minute, Try
        REAL*8 DeltaUT, tn, GST, t,
     &    l,m,n, EclpLong, EclpLat, Obliquity, MoonGHA, DGHA, LHAn,
     &    MoonGHAn, LHA, LST, JDF,
     &    LonEclSun, LonEclMoon, MeanAnomaly, MeanLong, ttdb,
     &    Sec, JDTemp, UTTemp, RtAsc, Decl

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        Error = 'ok'

        ! ---------- for once for MoonRise (1), ) THEN set (2) --------
        ! -------------- Make sure lon is within +- 180 deg -----------
        IF ( Lon .gt. Pi ) THEN
            Lon= Lon - 2.0D0*Pi
          ENDIF
        IF ( Lon .lt. -Pi ) THEN
            Lon= Lon + 2.0D0*Pi
          ENDIF

        Try= 1
        Opt= 1
        DO WHILE (Opt .le. 2)
            CALL INVJDay( JD, 0.0D0, Year,Month,Day,Hr,minute,Sec )
            CALL JDay( Year,Month,Day,0,0,0.0D0, JDTemp, JDF )
            UTTemp= 0.5D0

            IF ( Try .eq. 2 ) THEN
                IF ( Opt .eq. 1 ) THEN
                    UTTemp= 0.25D0
                  ELSE
                    UTTemp= 0.75D0
                  ENDIF
              ENDIF

            i = 0
            tn= UTTemp
            t = tn + 10.0D0
            JDTemp= JDTemp + JDF + UTTemp

            DO WHILE ( (DABS(tn-t).ge.0.008D0) .and. (i .le. 5) )
                TTDB = ( JDTemp + JDF - 2451545.0D0 ) / 36525.0D0
                EclpLong= 218.32D0 + 481267.883D0*TTDB
     &              + 6.29D0*DSIN( (134.9D0+477198.85D0*TTDB)*Deg2Rad )
     &              - 1.27D0*DSIN( (259.2D0-413335.38D0*TTDB)*Deg2Rad )
     &              + 0.66D0*DSIN( (235.7D0+890534.23D0*TTDB)*Deg2Rad )
     &              + 0.21D0*DSIN( (269.9D0+954397.70D0*TTDB)*Deg2Rad )
     &              - 0.19D0*DSIN( (357.5D0+ 35999.05D0*TTDB)*Deg2Rad )
     &              - 0.11D0*DSIN( (186.6D0+966404.05D0*TTDB)*Deg2Rad )
            EclpLat = 5.13D0*DSIN( ( 93.3D0+483202.03D0*TTDB)*Deg2Rad )
     &              + 0.28D0*DSIN( (228.2D0+960400.87D0*TTDB)*Deg2Rad )
     &              - 0.28D0*DSIN( (318.3D0+  6003.18D0*TTDB)*Deg2Rad )
     &              - 0.17D0*DSIN( (217.6D0-407332.20D0*TTDB)*Deg2Rad )
                EclpLong = DMOD( EclpLong*Deg2Rad, TwoPi )
                EclpLat  = DMOD( EclpLat*Deg2Rad, TwoPi )
                Obliquity= 23.439291D0 - 0.0130042D0*TTDB
                Obliquity= Obliquity *Deg2Rad
                ! ------- Find the geocentric direction cosines -------
                l= DCOS( EclpLat ) * DCOS( EclpLong )
                m= DCOS(Obliquity)*DCOS(EclpLat)*DSIN(EclpLong)
     &               - DSIN(Obliquity)*DSIN(EclpLat)
                n= DSIN(Obliquity)*DCOS(EclpLat)*DSIN(EclpLong)
     &               + DCOS(Obliquity)*DSIN(EclpLat)
                RtAsc= DATAN2( m,l )
                ! - Check that RtAsc is in the same quadrant as EclpLong
                IF ( EclpLong .lt. 0.0D0 ) THEN
                    EclpLong= EclpLong + TwoPi
                  ENDIF   
                IF ( DABS( EclpLong - RtAsc ) .gt. Pi*0.5D0 ) THEN
                    RtAsc= RtAsc + 0.5D0*Pi*
     &                     DINT( 0.5D0 + (EclpLong-RtAsc) / (0.5D0*Pi) )
                  ENDIF
                Decl = DASIN( n )
                CALL LSTIME( Lon,JDTemp + JDF,LST,GST )
                MoonGHAn= LST - Lon - RtAsc
                IF ( i .eq. 0 ) THEN
                    LHA = MoonGHAn + Lon
                    DGHA= 347.8D0 * Deg2Rad
                  ELSE
                    DGHA= (MoonGHAn - MoonGHA) / DeltaUT 
                  ENDIF
                IF ( DGHA .lt. 0.0D0 ) THEN
                    DGHA= DGHA + TwoPi/DABS(DeltaUT)
                  ENDIF
                LHAn= 0.00233D0 - (DSIN(Latgd)*DSIN(Decl)) /
     &                            (DCOS(Latgd)*DCOS(Decl))
                IF ( LHAn .gt. 1.0D0 ) THEN
                    LHAn= 0.0D0
                  ENDIF
                IF ( LHAn .lt. -1.0D0 ) THEN
                    LHAn= -1.0D0
                  ENDIF
                LHAn= DACOS( LHAn )
                IF ( Opt .eq. 1 ) THEN
                    LHAn= TwoPi - LHAn 
                  ENDIF 
                IF ( DABS( DGHA ) .gt. 0.0001D0 ) THEN
                    DeltaUT= (LHAn - LHA ) / DGHA
                  ELSE
                    DeltaUT= (LHAn - LHA )
                    DeltaUT= 1.0D0
                    Write( *,*)  'Fileout,x'
                  ENDIF
                t= tn 
                IF ( DABS( DeltaUT ) .gt. 0.5D0 ) THEN
                    IF ( DABS( DGHA ) .gt. 0.001D0 ) THEN
                        IF ( DeltaUT .lt. 0.0D0 ) THEN
                            DeltaUT= DeltaUT + TwoPi/DGHA
                            IF ( DABS( DeltaUT ) .gt. 0.51D0 ) THEN
                                i= 6 
                              ENDIF
                          ELSE
                            DeltaUT= DeltaUT - TwoPi/DGHA
                            IF ( DABS( DeltaUT ) .gt. 0.51D0 ) THEN
                                i= 6 
                              ENDIF
                          ENDIF
                      ELSE
                        DeltaUT= DeltaUT
                        Write(*,*) 'Fileout,y'
                      ENDIF 
                  ENDIF
                tn     = UTTemp + DeltaUT
                JDTemp = JDTemp + JDF - UTTemp + tn
                i = i + 1
                MoonGHA= MoonGHAn 

              ENDDO

            UTTemp= tn*24.0D0
            IF ( i .gt. 5 ) THEN
                UTTemp= 9999.99D0 
              ENDIF
            IF ( UTTemp .lt. 9999.0D0 ) THEN
                UTTemp= DMOD( UTTemp,24.0D0)
              ENDIF   
            IF ( UTTemp .lt. 0.0D0 ) THEN
                UTTemp= UTTEmp + 24.0D0
              ENDIF   
            IF ( UTTemp .gt. 900 ) THEN
                UTTemp= 24.0D0
              ENDIF

            IF (Opt .eq. 1 ) THEN
                UTMoonRise= UTTemp
              ENDIF
            IF (Opt .eq. 2 ) THEN
                UTMoonSet = UTTemp
              ENDIF

            Try= Try + 1 
            IF ( (i .gt. 5) .and. (Try .lt. 3) ) THEN
                Write(*,*) 'try #2 ',opt
              ELSE
                IF ( (i.gt.5) .and. (Try.gt.2) ) THEN
                    IF (Opt .eq. 1 ) THEN
                        Error = 'No Rise'
                      ENDIF
                    IF (Opt .eq. 2 ) THEN
                        Error = 'No Set'
                      ENDIF
                  ENDIF
                Opt= Opt + 1
                Try= 1
              ENDIF

          ENDDO

        ! ------------- determine phase ANGLE of the MOON --------------
        MeanLong= 280.4606184D0 + 36000.77005361D0*TTDB
        MeanLong= DMOD( MeanLong,360.0D0 )

        MeanAnomaly= 357.5277233D0 + 35999.05034D0*TTDB
        MeanAnomaly= DMOD( MeanAnomaly*Deg2Rad,TwoPi )
        IF ( MeanAnomaly .lt. 0.0D0 ) THEN
            MeanAnomaly= TwoPi + MeanAnomaly
          ENDIF

        LonEclSun= MeanLong + 1.914666471D0*DSIN(MeanAnomaly)
     &              + 0.019994643D0*DSIN(2.0D0*MeanAnomaly)

        LonEclMoon=   218.32D0 + 481267.883D0*TTDB
     &              + 6.29D0*DSIN( (134.9D0+477198.85D0*TTDB)*Deg2Rad )
     &              - 1.27D0*DSIN( (259.2D0-413335.38D0*TTDB)*Deg2Rad )
     &              + 0.66D0*DSIN( (235.7D0+890534.23D0*TTDB)*Deg2Rad )
     &              + 0.21D0*DSIN( (269.9D0+954397.70D0*TTDB)*Deg2Rad )
     &              - 0.19D0*DSIN( (357.5D0+ 35999.05D0*TTDB)*Deg2Rad )
     &              - 0.11D0*DSIN( (186.6D0+966404.05D0*TTDB)*Deg2Rad )
        LonEclMoon= DMOD( LonEclMoon, 360.0D0 )

        MoonPhaseAng= LonEclMoon - LonEclSun

        IF ( MoonPhaseAng .lt. 0.0D0 ) THEN
            MoonPhaseAng= 360.0D0 + MoonPhaseAng
          ENDIF

      RETURN
      END

      
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE SATFOV
*
*  this subroutine finds parameters reLating to a satellite's FOV.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Incl        - Inclination                    rad
*    Az          - Azimuth                        rad
*    SLatgd      - Geodetic Latitude of sat       rad
*    SLon        - Longitude of sat               rad
*    SAlt        - Altitudeof satellite           km
*    TFOV        - Total field of view            rad
*    EtaCtr      - Ctr where sensor looks         rad
*
*  Outputs       :
*    FovMax      - Maximum field of view          rad
*    TotalRng    -
*    RhoMax      -
*    RhoMin      -
*    TgtLat      -
*    TgtLon      -
*
*  Locals        :
*    r           -
*    etaHopriz   -
*    RhoHoriz    -
*    gamma       -
*    rho         -
*    FovMin      -
*    Lat         -
*    Lon         -
*    MaxLat      -
*    MinLKat     -
*    i           - Index
*
*  Coupling      :
*    PATH        - Finds tgt location given initial location, range, and az
*
*  References    :
*    Vallado       2007, 845, Eq 11-8 to Eq 11-13, Ex 11-1
*
* ------------------------------------------------------------------------------

      SUBROUTINE SATFOV      ( Incl, Az, SLatgd, SLon, SAlt, tFOV,
     &                         EtaCtr, FovMax, TotalRng, RhoMax, RhoMin,
     &                         TgtLat, TgtLon )
        IMPLICIT NONE
        REAL*8 Incl, Az, SLatgd, SLon, SAlt, tFOV, EtaCtr, FovMAx,
     &     TotalRng, RhoMax, RhoMin, TgtLat, TgtLon
* -----------------------------  Locals  ------------------------------
        INTEGER i
        REAL*8 r, EtaHoriz, rhoHoriz, gamma, rho, FovMin, Lat,
     &     Lon, maxLat, minLat

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        ! ------- Find satellite parameters and limiting cases --------
        r       = 1.0D0 + SAlt 
        EtaHoriz= DASIN(1.0D0/r) 
        RhoHoriz= r*DCOS(EtaHoriz) 

        ! ---------------- Find Ground range ANGLE --------------------
        FovMax= tFOV*0.5D0 + EtaCtr 
        Gamma = Pi - DASIN( r*DSIN(FovMax) )   ! must use larger ANGLE
        Rho   = DCOS( gamma ) + r*DCOS(FovMax) 
        RhoMax= DASIN( Rho*DSIN(FovMax) ) 

        ! -------- for minimum, If the sensor looks off axis ----------
        IF ( DABS(EtaCtr) .gt. 0.00001D0 ) THEN
            FovMin  = EtaCtr - tFOV*0.5D0
            Gamma   = Pi - DASIN( r*DSIN(FovMin) )  ! use larger
            Rho     = DCOS( gamma ) + r*DCOS(FovMin) 
            RhoMin  = DASIN( Rho*DSIN(FovMin) ) 
            TotalRng= RhoMax - RhoMin 
          ELSE
            ! --------------------- Nadir pointing --------------------
            FovMin  = 0.0D0
            RhoMin  = 0.0D0 
            TotalRng= 2.0D0*RhoMax  ! equal sided
          ENDIF

        ! -------------- Find location of center of FOV ---------------
        IF ( DABS(EtaCtr) .gt. 0.00001D0 ) THEN
            CALL PATH( SLatgd, SLon, RhoMin + TotalRng*0.5D0, Az,
     &                 Lat, Lon )
          ELSE
            Lat= SLatgd
            Lon= SLon 
          ENDIF 

        ! ----- Loop around the New circle with the sensor range ------
        DO i= 0 , 72
            Az= i*5.0D0/Rad2Deg
            CALL PATH( Lat, Lon, TotalRng*0.5D0, Az,  TgtLat, TgtLon )
            IF ( i .eq. 0 ) THEN
                MaxLat= TgtLat
              ENDIF 
            IF ( i .eq. 36 ) THEN
                MinLat= TgtLat
              ENDIF 
          ENDDO

      RETURN
      END

      
      
* ---------------------------------------------------------------------------
*
*                           SUBROUTINE SITE
*
*  this subroutine finds the position and velocity vectors for a SITE.  The
*    answer is returned in the Geocentric Equatorial (IJK) coordinate system.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Latgd       - Geodetic Latitude              -Pi/2 to Pi/2 rad
*    Alt         - Altitude                       km
*    LST         - Local SIDEREAL Time            -2Pi to 2Pi rad
*
*  OutPuts       :
*    RSecef      - ecef SITE position vector      km
*    VSecef      - ecef SITE velocity vector      km/s
*
*  Locals        :
*    EarthRate   - IJK Earth's rotation rate      rad/s
*    SinLat      - Variable containing  DSIN(Lat) rad
*    Temp        - Temporary Real value
*    Rdel        - Rdel component of SITE vector  km
*    Rk          - Rk component of SITE vector    km
*    CEarth      -
*
*  Coupling      :
*    CROSS       - CROSS product of two vectors
*
*  References    :
*    Vallado       2001, 404-407, Alg 47, Ex 7-1
*
* -----------------------------------------------------------------------------  

      SUBROUTINE SITE               ( Latgd,Alt,Lon, RSecef,VSecef )
        IMPLICIT NONE
        REAL*8 LatGd,Alt,Lon,RSecef(3),VSecef(3)
* -----------------------------  Locals  ------------------------------
        REAL*8 SinLat, CEarth, Rdel, Rk

        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
        SinLat      = DSIN( Latgd ) 

        ! ------  Find Rdel and Rk components of SITE vector  ---------
        CEarth= rekm / DSQRT( 1.0D0 - ( EESqrd*SinLat*SinLat ) )
        Rdel  = ( CEarth + Alt )*DCOS( Latgd ) 
        Rk    = ( (1.0D0-EESqrd)*CEarth + Alt )*SinLat

        ! ---------------  Find SITE position vector  -----------------
        RSecef(1) = Rdel * DCOS( Lon )
        RSecef(2) = Rdel * DSIN( Lon )
        RSecef(3) = Rk

        ! ---------------  Find SITE velocity vector  ------------------
        VSecef(1) = 0.0D0
        VSecef(2) = 0.0D0
        VSecef(3) = 0.0D0
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE ANGLESLAPLACE
*
*  this subroutine solves the problem of orbit determination using three
*    optical sightings and the method of Laplace.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Alpha1       - Right Ascension #1            rad
*    Alpha2       - Right Ascension #2            rad
*    Alpha3       - Right Ascension #3            rad
*    Delta1       - Declination #1                rad
*    Delta2       - Declination #2                rad
*    Delta3       - Declination #3                rad
*    JD1          - Julian Date of 1st sighting   Days from 4713 BC
*    JD2          - Julian Date of 2nd sighting   Days from 4713 BC
*    JD3          - Julian Date of 3rd sighting   Days from 4713 BC
*    RS1          - IJK SITE position vector #1   km
*    RS2          - IJK SITE position vector #2   km
*    RS3          - IJK SITE position vector #3   km
*
*  OutPuts        :
*    R            - IJK position vector           km
*    V            - IJK velocity vector           km / s
*
*  Locals         :
*    L1           - Line of SIGHT vector for 1st
*    L2           - Line of SIGHT vector for 2nd
*    L3           - Line of SIGHT vector for 3rd
*    LDot         - 1st derivative of L2
*    LDDot        - 2nd derivative of L2
*    RS2Dot       - 1st Derivative of RS2 - vel
*    RS2DDot      - 2nd Derivative of RS2
*    t12t13       - (t1-t2) * (t1-t3)
*    t21t23       - (t2-t1) * (t2-t3)
*    t31t32       - (t3-t1) * (t3-t2)
*    i            - index
*    D            -
*    D1           -
*    D2           -
*    D3           -
*    D4           -
*    OldR         - Previous iteration on r
*    Rho          - Range from SITE to satellite at t2
*    RhoDot       -
*    DMat         -
*    D1Mat        -
*    D2Mat        -
*    D3Mat        -
*    D4Mat        -
*    EarthRate    - Angular rotation of the earth
*    L2DotRS      - Vector L2 Dotted with RSecef
*    Temp         - Temporary vector
*    Temp1        - Temporary vector
*    Small        - Tolerance
*    Roots        -
*
*  Coupling       :
*    MAG          - Magnitude of a vector
*    DETERMINANT  - Evaluate the determinant of a matrix
*    CROSS        - CROSS product of two vectors
*    NORM         - Normlize a matrix
*    FACTOR       - Find the roots of a polynomial
*
*  References     :
*    Vallado       2001, 413-417
*
* ------------------------------------------------------------------------------  

      SUBROUTINE ANGLESLAPLACE ( Delta1,Delta2,Delta3,Alpha1,Alpha2,
     &                      Alpha3,JD1,JD2,JD3,RS1,RS2,RS3, r2,v2 )
        IMPLICIT NONE
        REAL*8 Delta1,Delta2,Delta3,Alpha1,Alpha2,Alpha3,JD1,JD2,JD3,
     &         RS1(3),RS2(3),RS3(3),r2(3),v2(3)
        EXTERNAL DETERMINANT, Dot, Mag
* -----------------------------  Locals  ------------------------------
        INTEGER i, j, k
        REAL*8 Small, Poly(16),Roots(15,2), MAG
        REAL*8 DMat(3,3), DMat1(3,3), DMat2(3,3), DMat3(3,3),DMat4(3,3)
        REAL*8 L1(3), L2(3), L3(3), LDot(3), LDDot(3), RS2Dot(3),
     &         RS2DDot(3), EarthRate(3), Temp(3), Temp1(3),magr2,
     &         D, D1, D2, D3, D4, Rho, RhoDot, t1t13, t1t3, t31t3,
     &         tau1, tau3, BigR2, L2DotRS, Determinant, Dot, magtemp,
     &         magtemp1, magrs2

        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
c        TUDay        =     0.00933809017716D0
*        TUDay        =    58.132440906D0
        Small        =     0.0000001D0
        EarthRate(1)= 0.0D0
        EarthRate(2)= 0.0D0
        EarthRate(3)= OmegaEarth

        JD1= JD1*86400.0D0    ! days to sec
        JD2= JD2*86400.0D0
        JD3= JD3*86400.0D0

        ! ---------- set middle to 0, find deltas to others -----------
        tau1= JD1-JD2 
        tau3= JD3-JD2 

        ! --------------- Find Line of SIGHT vectors ------------------
        L1(1)= DCOS(Delta1)*DCOS(Alpha1)
        L1(2)= DCOS(Delta1)*DSIN(Alpha1)
        L1(3)= DSIN(Delta1)
        L2(1)= DCOS(Delta2)*DCOS(Alpha2)
        L2(2)= DCOS(Delta2)*DSIN(Alpha2)
        L2(3)= DSIN(Delta2)
        L3(1)= DCOS(Delta3)*DCOS(Alpha3)
        L3(2)= DCOS(Delta3)*DSIN(Alpha3)
        L3(3)= DSIN(Delta3)

        ! -------------------------------------------------------------
*       Using Lagrange Interpolation formula to derive an expression
*       for L(t), substitute t=t2 and differentiate to obtain the
*       derivatives of L.
        ! -------------------------------------------------------------
        t1t13= 1.0D0 / (tau1*(tau1-tau3)) 
        t1t3 = 1.0D0 / (tau1*tau3) 
        t31t3= 1.0D0 / ((tau3-tau1)*tau3) 
        DO i= 1 , 3
            LDot(i)=      ( -tau3 * t1t13 )*L1(i) +
     &               ( (-tau1-tau3) * t1t3  )*L2(i) +
     &                      ( -tau1 * t31t3 )*L3(i)
            LDDot(i)= ( 2.0D0 * t1t13 )*L1(i) +
     &                  ( 2.0D0 * t1t3  )*L2(i) +
     &                  ( 2.0D0 * t31t3 )*L3(i)
          ENDDO
        CALL NORM( LDot,  LDot )
        CALL NORM( LDDot, LDDot )

        ! ------------------- Find 2nd derivative of RSecef ---------------
        CALL CROSS( RS1,RS2, Temp )
        magtemp = MAG(Temp)
        CALL CROSS( RS2,RS3, Temp1 )
        magtemp1 = MAG(Temp1)

*      needs a different test xxxx!!  
        IF ( ( DABS(magtemp) .gt. Small ) .and.
     &     ( DABS( magtemp1) .gt. Small )  ) THEN
           ! ------------ All sightings from one SITE -----------------
*          fix this testhere  
            DO i= 1 , 3
                RS2Dot(i)=      ( -tau3 * t1t13 )*RS1(i) +
     &                     ( (-tau1-tau3) * t1t3  )*RS2(i) +
     &                            ( -tau1 * t31t3 )*RS3(i)
                RS2DDot(i)= ( 2.0D0 * t1t13 )*RS1(i) +
     &                        ( 2.0D0 * t1t3  )*RS2(i) +
     &                        ( 2.0D0 * t31t3 )*RS3(i)
              ENDDO

            CALL CROSS( EarthRate,RS2,     RS2Dot )
            CALL CROSS( EarthRate,RS2Dot,  RS2DDot )
          ELSE
            ! ---------- Each sighting from a different SITE ----------
            DO i= 1 , 3
                RS2Dot(i)=      ( -tau3 * t1t13 )*RS1(i) +
     &                     ( (-tau1-tau3) * t1t3  )*RS2(i) +
     &                            ( -tau1 * t31t3 )*RS3(i)
                RS2DDot(i)= ( 2.0D0 * t1t13 )*RS1(i) +
     &                        ( 2.0D0 * t1t3  )*RS2(i) +
     &                        ( 2.0D0 * t31t3 )*RS3(i)
              ENDDO
          ENDIF 

        DO i= 1 , 3
            DMat(i,1) =2.0D0 * L2(i)
            DMat(i,2) =2.0D0 * LDot(i)
            DMat(i,3) =2.0D0 * LDDot(i)

            ! ----------------  Position determinants -----------------
            DMat1(i,1) =L2(i)
            DMat1(i,2) =LDot(i)
            DMat1(i,3) =RS2DDot(i)
            DMat2(i,1) =L2(i)
            DMat2(i,2) =LDot(i)
            DMat2(i,3) =RS2(i)

            ! ------------  Velocity determinants ---------------------
            DMat3(i,1) =L2(i)
            DMat3(i,2) =RS2DDot(i)
            DMat3(i,3) =LDDot(i)
            DMat4(i,1) =L2(i)
            DMat4(i,2) =RS2(i)
            DMat4(i,3) =LDDot(i)
          ENDDO

        D = DETERMINANT(DMat,3) 
        D1= DETERMINANT(DMat1,3) 
        D2= DETERMINANT(DMat2,3) 
        D3= DETERMINANT(DMat3,3) 
        D4= DETERMINANT(DMat4,3) 
* 
      ! ---------------  Iterate to find Rho magnitude ----------------
*     magr= 1.5D0   ! First Guess
*     Write( 'Input initial guess for magr ' )
*     Read( magr )
*     i= 1 
*     REPEAT
*         OldR= magr
*         Rho= -2.0D0*D1/D - 2.0D0*D2/(magr**3*D)
*         magr= DSQRT( Rho*Rho + 2.0D0*Rho*L2DotRS + magRS2**2 )
*         INC(i) 
*         magr= (OldR - magr ) / 2.0D0             ! Simple bissection
*         WriteLn( FileOut,'Rho guesses ',i:2,'Rho ',Rho:14:7,' magr ',magr:14:7,oldr:14:7 )
*! seems to converge, but wrong Numbers
*         INC(i) 
*     UNTIL ( DABS( OldR-magR ) .lt. Small ) .or. ( i .ge. 30 )
   

        IF ( DABS(D) .gt. 0.000001D0 ) THEN
            ! --------------- Solve eighth order poly -----------------
            L2DotRS= DOT( L2,RS2 ) 
            magrs2 = MAG(rs2)
            Poly( 1)=  1.0D0  ! r2^8th variable!!!!!!!!!!!!!!
            Poly( 2)=  0.0D0
            Poly( 3)=  (L2DotRS*4.0D0*D1/D - 4.0D0*D1*d1/(D*D)
     &                 - magRS2**2 )
            Poly( 4)=  0.0D0
            Poly( 5)=  0.0D0
            Poly( 6)=  Mu*(L2DotRS*4.0D0*D2/D - 8.0D0*D1*D2/(D*D) )
            Poly( 7)=  0.0D0
            Poly( 8)=  0.0D0
            Poly( 9)=  -4.0D0*Mu*D2*D2/(D*D)
            Poly(10)=  0.0D0
            Poly(11)=  0.0D0
            Poly(12)=  0.0D0
            Poly(13)=  0.0D0
            Poly(14)=  0.0D0
            Poly(15)=  0.0D0
            Poly(16)=  0.0D0
            CALL FACTOR( Poly,8,  Roots )

           ! ------------------ Find correct (xx) root ----------------
            BigR2= 0.0D0 
            DO j= 1 , 8
*                IF ( DABS( Roots(j,2) ) .lt. Small ) THEN
*                    WriteLn( 'Root ',j,Roots(j,1),' + ',Roots(j,2),'j' )
*        temproot= roots(j,1)*roots(j,1)
*        temproot= Temproot*TempRoot*TempRoot*TempRoot +
*                  Poly(3)*TempRoot*TempRoot*TempRoot + Poly(6)*roots(j,1)*Temproot + Poly(9)
*                    WriteLn( FileOut,'Root ',j,Roots(j,1),' + ',Roots(j,2),'j  value = ',temproot )
                    IF ( Roots(j,1) .gt. BigR2 ) THEN
                        BigR2= Roots(j,1)
                      ENDIF
*                  ENDIF
              ENDDO
        Write(*,*) 'BigR2 ',BigR2
        Write(*,*) 'Keep this root ? '
        READ(*,*) BigR2

            Rho= -2.0D0*D1/D - 2.0D0*Mu*D2 / (BigR2*BigR2*BigR2*D) 

            ! --------- Find the middle position vector ---------------
            DO k= 1 , 3
                r2(k)= Rho*L2(k) + RS2(k)
              ENDDO
            magr2 = MAG( r2 )
            ! ---------------- Find RhoDot magnitude ------------------
            RhoDot= -D3/D - Mu*D4/(magr2**3*D)
*        WriteLn( FileOut,'Rho ',Rho:14:7 )
*        WriteLn( FileOut,'RhoDot ',RhoDot:14:7 )

            ! -------------- Find middle velocity vector --------------
            DO i= 1 , 3
                V2(i)= RhoDot*L2(i) + Rho*LDot(i) + RS2Dot(i)
              ENDDO
         ELSE
           Write(*,*) 'Determinant value was zero ',D
         ENDIF

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE ANGLESGAUSS
*
*  this subroutine solves the problem of orbit determination using three
*    optical sightings.  The solution SUBROUTINE uses the Gaussian technique.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Alpha1       - Right Ascension #1            rad
*    Alpha2       - Right Ascension #2            rad
*    Alpha3       - Right Ascension #3            rad
*    Delta1       - Declination #1                rad
*    Delta2       - Declination #2                rad
*    Delta3       - Declination #3                rad
*    JD1          - Julian Date of 1st sighting   Days from 4713 BC
*    JD2          - Julian Date of 2nd sighting   Days from 4713 BC
*    JD3          - Julian Date of 3rd sighting   Days from 4713 BC
*    RSecef           - IJK SITE position vector      km
*
*  OutPuts        :
*    R            - IJK position vector at t2     km
*    V            - IJK velocity vector at t2     km / s
*
*  Locals         :
*    L1           - Line of SIGHT vector for 1st
*    L2           - Line of SIGHT vector for 2nd
*    L3           - Line of SIGHT vector for 3rd
*    Tau          - Taylor expansion series about
*                   Tau ( t - to )
*    TauSqr       - Tau squared
*    t21t23       - (t2-t1) * (t2-t3)
*    t31t32       - (t3-t1) * (t3-t2)
*    i            - index
*    D            -
*    Rho          - Range from SITE to sat at t2  km
*    RhoDot       -
*    DMat         -
*    RS1          - SITE vectors
*    RS2          -
*    RS3          -
*    EarthRate    - Velocity of Earth rotation
*    P            -
*    Q            -
*    OldR         -
*    OldV         -
*    F1           - F coefficient
*    G1           -
*    F3           -
*    G3           -
*    L2DotRS      -
*
*  Coupling       :
*    Detrminant   - Evaluate the determinant of a matrix
*    FACTOR       - Find roots of a polynomial
*    MATMULT      - Multiply two matrices together
*    GIBBS        - GIBBS method of orbit determination
*    HGIBBS       - Herrick GIBBS method of orbit determination
*    ANGLE        - ANGLE between two vectors
*
*  References     :
*    Vallado       2001, 417-421, Alg 49, Ex 7-2 (425-427)
*
* ------------------------------------------------------------------------------  

      SUBROUTINE ANGLESGAUSS ( Delta1,Delta2,Delta3,Alpha1,Alpha2,
     &                      Alpha3,JD1,JD2,JD3,RS1,RS2,RS3, r2,v2 )
        IMPLICIT NONE
        REAL*8 Delta1,Delta2,Delta3,Alpha1,Alpha2,Alpha3,JD1,JD2,JD3,
     &         RS1(3),RS2(3),RS3(3),r2(3),v2(3)
        EXTERNAL Determinant, Dot, MAG
* -----------------------------  Locals  ------------------------------
        INTEGER i, ll, j
        REAL*8 small, Roots(15,2), Poly(16),
     &         r1(3), r3(3), L1(3), L2(3), L3(3)
        CHARACTER*12 Error
        REAL*8 LMatIi(3,3), CMat(3,1),RhoMat(3,1), LMatI(3,3),
     &         RSMat(3,3),LIR(3,3), Determinant, Dot, magrs2
        REAL*8 rDot, tau1, tau3, u, uDot, p, MAG, magr2, magr1, magr3,
     &         f1, g1, f3, g3, a, ecc, incl, omega, argp,
     &         Nu, m, l, ArgPer, BigR2, a1, a1u, a3, a3u, d, d1,
     &         d2, c1, c3, L2DotRS, rhoold1, rhoold2, rhoold3,
     &         rad,  theta, theta1, copa, TauSqr

        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
*        TUDay        =    58.132440906D0
c        TUDay        =     0.00933809017716D0
        Small        =     0.0000001D0
        rad    = 57.29577951308D0 

        ! ---------- set middle to 0, find deltas to others -----------
        tau1= (JD1-JD2)*86400.0D0
        tau3= (JD3-JD2)*86400.0D0

        ! ----------------  Find Line of SIGHT vectors  ---------------
        L1(1)= DCOS(Delta1)*DCOS(Alpha1)
        L1(2)= DCOS(Delta1)*DSIN(Alpha1)
        L1(3)= DSIN(Delta1)

        L2(1)= DCOS(Delta2)*DCOS(Alpha2)
        L2(2)= DCOS(Delta2)*DSIN(Alpha2)
        L2(3)= DSIN(Delta2)

        L3(1)= DCOS(Delta3)*DCOS(Alpha3)
        L3(2)= DCOS(Delta3)*DSIN(Alpha3)
        L3(3)= DSIN(Delta3)

        ! ------------- Find L matrix and determinant -----------------
        ! --------- Called LMatI since it is only used for determ -----

        DO i= 1 , 3
              LMatIi(i,1) =L1(i)
              LMatIi(i,2) =L2(i)
              LMatIi(i,3) =L3(i)
              RSMat(i,1) =RS1(i)
              RSMat(i,2) =RS2(i)
              RSMat(i,3) =RS3(i)
          ENDDO

        D= DETERMINANT(LMatIi,3) 
        ! ------------------ Now assign the inverse -------------------
        LMatI(1,1) = ( L2(2)*L3(3)-L2(3)*L3(2)) / D
        LMatI(2,1) = (-L1(2)*L3(3)+L1(3)*L3(2)) / D
        LMatI(3,1) = ( L1(2)*L2(3)-L1(3)*L2(2)) / D
        LMatI(1,2) = (-L2(1)*L3(3)+L2(3)*L3(1)) / D
        LMatI(2,2) = ( L1(1)*L3(3)-L1(3)*L3(1)) / D
        LMatI(3,2) = (-L1(1)*L2(3)+L1(3)*L2(1)) / D
        LMatI(1,3) = ( L2(1)*L3(2)-L2(2)*L3(1)) / D
        LMatI(2,3) = (-L1(1)*L3(2)+L1(2)*L3(1)) / D
        LMatI(3,3) = ( L1(1)*L2(2)-L1(2)*L2(1)) / D

        CALL MATMULT( LMatI,RSMat,3,3,3, 3,3,3,   LIR )

        ! ------------ Find f and g series at 1st and 3rd obs ---------
*      speed by assuming circ sat vel for uDot here ??
*      some similartities in 1/6t3t1 ...  
        ! --- keep separated this time ----
        a1 =  Tau3 / (Tau3 - Tau1) 
        a1u=  (Tau3*((tau3-Tau1)*(Tau3-Tau1) - Tau3*Tau3 )) /
     &        (6.0D0*(Tau3 - Tau1))
        a3 = -Tau1 / (Tau3 - Tau1) 
        a3u= -(Tau1*((tau3-Tau1)*(Tau3-Tau1) - Tau1*Tau1 )) /
     &        (6.0D0*(Tau3 - Tau1))

        ! --- Form initial guess of r2 ----
        d1=  LIR(2,1)*a1 - LIR(2,2) + LIR(2,3)*a3
        d2=  LIR(2,1)*a1u + LIR(2,3)*a3u

        ! ------- Solve eighth order poly NOT same as LAPLACE ---------
        L2DotRS= DOT( L2,RS2 ) 
        magrs2 = MAG(rs2)
        Poly( 1)=  1.0D0  ! r2^8th variable!!!!!!!!!!!!!!
        Poly( 2)=  0.0D0
        Poly( 3)=  -(D1*D1 + 2.0D0*D1*L2DotRS + magRS2**2)
        Poly( 4)=  0.0D0
        Poly( 5)=  0.0D0
        Poly( 6)=  -2.0D0*Mu*(L2DotRS*D2 + D1*D2)
        Poly( 7)=  0.0D0
        Poly( 8)=  0.0D0
        Poly( 9)=  -Mu*Mu*D2*D2
        Poly(10)=  0.0D0
        Poly(11)=  0.0D0
        Poly(12)=  0.0D0
        Poly(13)=  0.0D0
        Poly(14)=  0.0D0
        Poly(15)=  0.0D0
        Poly(16)=  0.0D0
        CALL FACTOR( Poly,8,  Roots )

        ! ------------------ Select the correct root ------------------
        BigR2= 0.0D0 
        DO j= 1 , 8
*            IF ( DABS( Roots(j,2) ) .lt. Small ) THEN
*     temproot= roots(j,1)*roots(j,1)
*     temproot= Temproot*TempRoot*TempRoot*TempRoot +
*              Poly(3)*TempRoot*TempRoot*TempRoot + Poly(6)*roots(j,1)*Temproot + Poly(9)
*                WriteLn( FileOut,'Root ',j,Roots(j,1),' + ',Roots(j,2),'j  value = ',temproot )
                IF ( Roots(j,1) .gt. BigR2 ) THEN
                    BigR2= Roots(j,1)
                ENDIF  ! IF (
*             ENDIF  ! IF (
          ENDDO
        ! ------------ Solve matrix with u2 better known --------------
        u= Mu / ( BigR2*BigR2*BigR2 ) 

        c1= a1+a1u*u 
        c3= a3+a3u*u 
          CMat(1,1)= -c1
          CMat(2,1)= 1.0D0
          CMat(3,1)= -c3
        CALL MATMULT( LIR,CMat,3,3,1, 3,3,1,  RhoMat )

        Rhoold1=  RhoMat(1,1)/c1
        Rhoold2= -RhoMat(2,1)
        Rhoold3=  RhoMat(3,1)/c3


      ! -------- Loop through the refining process ------------  for WHILE () DO
      DO ll= 1 , 3
            Write( *,*) ' Iteration # ',ll
            ! ---------- Now form the three position vectors ----------
            DO i= 1 , 3
                R1(i)=  RhoMat(1,1)*L1(i)/c1 + RS1(i)
                R2(i)= -RhoMat(2,1)*L2(i)    + RS2(i)
                R3(i)=  RhoMat(3,1)*L3(i)/c3 + RS3(i)
              ENDDO

            CALL GIBBS(r1,r2,r3,  v2,theta,theta1,copa,error )

            IF ( (Error .ne. 'ok') .and. (copa .lt. 1.0D0/Rad) ) THEN
                ! --- HGibbs to get middle vector ----
                CALL HERRGIBBS(r1,r2,r3,JD1,JD2,JD3,
     &                     v2,theta,theta1,copa,error )
*                WriteLn( FileOut,'hgibbs ' )
              ENDIF 

            CALL rv2coe( r2,v2, p,a,ecc,incl,omega,argp,Nu,m,u,l,ArgPer)
            magr2 = MAG(r2)

        IF ( ll .le. 2 ) THEN
            ! --- Now get an improved estimate of the f and g series --
*       .or. can the analytic functions be found now??  
            u= Mu / ( magr2**3 )
            rDot= DOT(r2,v2)/magr2
            uDot= (-3.0D0*Mu*RDot) / (magr2**4)

            TauSqr= Tau1*Tau1 
            f1=  1.0D0 - 0.5D0*u*TauSqr -(1.0D0/6.0D0)*UDot*TauSqr*Tau1
     &                 + (1.0D0/24.0D0) * u*u*TauSqr*TauSqr
     &                 + (1.0D0/30.0D0)*U*UDot*TauSqr*TauSqr*Tau1
            g1= Tau1 - (1.0D0/6.0D0)*u*Tau1*TauSqr - (1.0D0/12.0D0) *
     &                 UDot*TauSqr*TauSqr
     &                 + (1.0D0/120.0D0)*u*u*TauSqr*TauSqr*Tau1
     &                 + (1.0D0/120.0D0)*u*UDot*TauSqr*TauSqr*TauSqr
            TauSqr= Tau3*Tau3 
            f3=  1.0D0 - 0.5D0*u*TauSqr -(1.0D0/6.0D0)*UDot*TauSqr*Tau3
     &                 + (1.0D0/24.0D0) * u*u*TauSqr*TauSqr
     &                 + (1.0D0/30.0D0)*U*UDot*TauSqr*TauSqr*Tau3
            g3= Tau3 - (1.0D0/6.0D0)*u*Tau3*TauSqr - (1.0D0/12.0D0) *
     &                 UDot*TauSqr*TauSqr
     &                 + (1.0D0/120.0D0)*u*u*TauSqr*TauSqr*Tau3
     &                 + (1.0D0/120.0D0)*u*UDot*TauSqr*TauSqr*TauSqr
          ELSE
            ! -------- Now use exact method to find f and g -----------
            CALL ANGLE( R1,R2, Theta )
            CALL ANGLE( R2,R3, Theta1 )
            magr1 = MAG(r1)
            magr3 = MAG(r3)

            f1= 1.0D0 - ( (magR1*(1.0D0 - DCOS(Theta)) / p ) )
            g1= ( magR1*magR2*DSIN(-theta) ) / DSQRT( p )  ! - ANGLE because backwards!!
            f3= 1.0D0 - ( (magR3*(1.0D0 - DCOS(Theta1)) / p ) )
            g3= ( magR3*magR2*DSIN(theta1) ) / DSQRT( p )

         ENDIF
            c1=  g3 / (f1*g3 - f3*g1)
            c3= -g1 / (f1*g3 - f3*g1) 
            ! ----- Solve for all three ranges via matrix equation ----
            CMat(1,1)= -c1
            CMat(2,1)= 1.0D0
            CMat(3,1)= -c3
            CALL MATMULT( LIR,CMat,3,3,1, 3,3,1,  RhoMat )

            ! ----------------- Check for convergence -----------------

          ENDDO   ! DO WHILE the ranges are still changing

        ! ---------------- Find all three vectors ri ------------------
        DO i= 1 , 3
            R1(i)=  RhoMat(1,1)*L1(i)/c1 + RS1(i)
            R2(i)= -RhoMat(2,1)*L2(i)    + RS2(i)
            R3(i)=  RhoMat(3,1)*L3(i)/c3 + RS3(i)
          ENDDO
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE RV_RADEC
*
*  this subroutine converts the right ascension and declination values with
*    position and velocity vectors of a satellite. Uses velocity vector to
*    find the solution of singular cases.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Rijk        - IJK position vector            km
*    Vijk        - IJK velocity vector            km/s
*    Direction   - Which set of vars to output    FROM  TOO
*
*  OutPuts       :
*    Rr          - Radius of the satellite        km
*    RtAsc       - Right Ascension                rad
*    Decl        - Declination                    rad
*    DRr         - Radius of the satellite rate   km/s
*    DRtAsc      - Right Ascension rate           rad/s
*    DDecl       - Declination rate               rad/s
*
*  Locals        :
*    Temp        - Temporary position vector
*    Temp1       - Temporary variable
*
*  Coupling      :
*    DOT         - DOT product of two vectors
*
*  References    :
*    Vallado       2001, 246-248, Alg 25
*
* ------------------------------------------------------------------------------  

      SUBROUTINE RV_RADEC    ( Rijk,Vijk, Direction, rr,RtAsc,Decl,
     &                         DRr,DRtAsc,DDecl )
        IMPLICIT NONE
        REAL*8 Rijk(3),Vijk(3),rr,RtAsc,Decl,DRr,DRtAsc,DDecl
        CHARACTER*4 Direction
        EXTERNAL Dot, MAG
* -----------------------------  Locals  ------------------------------
        INTEGER Small
        REAL*8 Temp, Temp1, Dot, MAG

        ! --------------------  Implementation   ----------------------
        Small        = 0.00000001D0
        IF ( Direction .eq. 'FROM' ) THEN
            Rijk(1)= rr*DCOS(Decl)*DCOS(RtAsc)
            Rijk(2)= rr*DCOS(Decl)*DSIN(RtAsc)
            Rijk(3)= rr*DSIN(Decl)
            Vijk(1)= DRr*DCOS(Decl)*DCOS(RtAsc) -
     &               rr*DSIN(Decl)*DCOS(RtAsc)*DDecl
     &               - rr*DCOS(Decl)*DSIN(RtAsc)*DRtAsc
            Vijk(2)= DRr*DCOS(Decl)*DSIN(RtAsc) -
     &               rr*DSIN(Decl)*DSIN(RtAsc)*DDecl
     &               + rr*DCOS(Decl)*DCOS(RtAsc)*DRtAsc
            Vijk(3)= DRr*DSIN(Decl) + rr*DCOS(Decl)*DDecl
          ELSE
            ! ------------- Calculate Angles and Rates ----------------
            rr = MAG(Rijk)
            Temp= DSQRT( Rijk(1)*Rijk(1) + Rijk(2)*Rijk(2) )
            IF ( Temp .lt. Small ) THEN
                RtAsc= DATAN2( Vijk(2), Vijk(1) )
              ELSE
                RtAsc= DATAN2( Rijk(2), Rijk(1) )
              ENDIF
            Decl= DASIN( Rijk(3)/rr )

            Temp1= -Rijk(2)*Rijk(2) - Rijk(1)*Rijk(1)  ! different now
            DRr= DOT(Rijk,Vijk)/rr 
            IF ( DABS(Temp1) .gt. Small ) THEN
                DRtAsc= ( Vijk(1)*Rijk(2) - Vijk(2)*Rijk(1) ) / Temp1
              ELSE
                DRtAsc= 0.0D0
              ENDIF
            IF ( DABS( Temp ) .gt. Small ) THEN
                DDecl= ( Vijk(3) - DRr*DSIN( Decl ) ) / Temp
              ELSE
                DDecl= 0.0D0
              ENDIF
          ENDIF
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE RV_TRADEC
*
*  this subroutine converts topocentric right-ascension declination with
*    position and velocity vectors. Uses velocity vector to find the
*    solution of singular cases.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Rijk        - IJK position vector            km
*    Vijk        - IJK velocity vector            km/s
*    RSecef          - IJK SITE position vector       km
*    Direction   - Which set of vars to output    FROM  TOO
*
*  OutPuts       :
*    Rho         - Top Radius of the sat          km
*    TRtAsc      - Top Right Ascension            rad
*    TDecl       - Top Declination                rad
*    DRho        - Top Radius of the sat rate     km/s
*    TDRtAsc     - Top Right Ascension rate       rad/s
*    TDDecl      - Top Declination rate           rad/s
*
*  Locals        :
*    RhoV        - IJK Range Vector from SITE     km
*    DRhoV       - IJK Velocity Vector from SITE  km / s
*    Temp        - Temporary REAL*8 value
*    Temp1       - Temporary REAL*8 value
*    i           - Index
*
*  Coupling      :
*    MAG         - Magnitude of a vector
*    LNCOM2      - Linear combination of 2 vectors
*    ADDVEC      - Add two vectors
*    DOT         - DOT product of two vectors
*
*  References    :
*    Vallado       2001, 248-250, Alg 26
*
* ------------------------------------------------------------------------------  

      SUBROUTINE RV_TRADEC   ( Rijk,Vijk,RSecef, Direction, Rho,TRtAsc,
     &                         TDecl,DRho,DTRtAsc,DTDecl )
        IMPLICIT NONE
        REAL*8 Rijk(3),VIjk(3),RSecef(3),Rho,TRtAsc,TDecl,DRho,DTRtAsc,
     &         DTDecl, MAG
        CHARACTER*4 Direction
        EXTERNAL DOT, MAG
* -----------------------------  Locals  ------------------------------
        INTEGER i
        REAL*8 Small, temp, temp1, RhoV(3),DRhoV(3), Dot, Magrhov

        ! --------------------  Implementation   ----------------------
        Small        = 0.00000001D0
        IF ( Direction .eq. 'FROM' ) THEN
            ! --------  Calculate Topocentric Vectors -----------------
            RhoV(1)= Rho*DCOS(TDecl)*DCOS(TRtAsc)
            RhoV(2)= Rho*DCOS(TDecl)*DSIN(TRtAsc)
            RhoV(3)= Rho*DSIN(TDecl)

            DRhoV(1)= DRho*DCOS(TDecl)*DCOS(TRtAsc)
     &                  - Rho*DSIN(TDecl)*DCOS(TRtAsc)*DTDecl
     &                  - Rho*DCOS(TDecl)*DSIN(TRtAsc)*DTRtAsc
            DRhoV(2)= DRho*DCOS(TDecl)*DSIN(TRtAsc)
     &                  - Rho*DSIN(TDecl)*DSIN(TRtAsc)*DTDecl
     &                  + Rho*DCOS(TDecl)*DCOS(TRtAsc)*DTRtAsc
            DRhoV(3)= DRho*DSIN(TDecl) + Rho*DCOS(TDecl)*DTDecl

            ! ------ Find IJK range vector from SITE to satellite -----
            CALL ADDVEC( RhoV,RSecef,  Rijk )
            DO i=1,3
                Vijk(i)= DRhoV(i)
              ENDDO

          ELSE
            ! ------ Find IJK range vector from SITE to satellite -----
            CALL LNCOM2( 1.0D0,-1.0D0, Rijk,RSecef,  RhoV )
            DO i=1,3
                DRhoV(i)= Vijk(i)  ! Same for topocentric
              ENDDO
            magrhov = MAG(rhoV)

            ! ------- Calculate Topocentric ANGLE and Rate Values -----
            Rho= MAG(Rhov)
            Temp= DSQRT( RhoV(1)*RhoV(1) + RhoV(2)*RhoV(2) )
            IF ( Temp .lt. Small ) THEN
                TRtAsc= DATAN2( DRhoV(2), DRhoV(1) )
              ELSE
                TRtAsc= DATAN2( RhoV(2), RhoV(1) )
              ENDIF

            TDecl= DASIN( RhoV(3)/magRhoV )

            Temp1= -RhoV(2)*RhoV(2) - RhoV(1)*RhoV(1)  ! different now
            DRho= DOT(RhoV,DRhoV)/Rho
            IF ( DABS(Temp1) .gt. Small ) THEN
                DTRtAsc= ( DRhoV(1)*RhoV(2) - DRhoV(2)*RhoV(1) ) / Temp1
              ELSE
                DTRtAsc= 0.0D0
              ENDIF
            IF ( DABS( Temp ) .gt. Small ) THEN
                DTDecl= ( DRhoV(3) - DRho*DSIN( TDecl ) ) / Temp
              ELSE
                DTDecl= 0.0D0
              ENDIF

          ENDIF
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE RV_RAZEL
*
*  this subroutine converts Range, Azimuth, and Elevation and their rates with
*    the Geocentric Equatorial (IJK) Position and Velocity vectors.  Notice the
*    value of small as it can affect rate term calculations. Uses velocity
*    vector to find the solution of singular cases.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Rijk        - IJK position vector            km
*    Vijk        - IJK velocity vector            km/s
*    RSecef          - IJK SITE Position Vector       km
*    Latgd       - Geodetic Latitude              -Pi/2 to Pi/2 rad
*    Lon         - Site longitude                 -Pi to Pi rad
*    Direction   - Which set of vars to output    FROM  TOO
*
*  Outputs       :
*    Rho         - Satellite Range from SITE      km
*    Az          - Azimuth                        0.0D0 to 2Pi rad
*    El          - Elevation                      -Pi/2 to Pi/2 rad
*    DRho        - Range Rate                     km / s
*    DAz         - Azimuth Rate                   rad / s
*    DEl         - Elevation rate                 rad / s
*
*  Locals        :
*    RhoVijk     - IJK Range Vector from SITE     km
*    DRhoVijk    - IJK Velocity Vector from SITE  km / s
*    Rhosez      - SEZ Range vector from SITE     km
*    DRhosez     - SEZ Velocity vector from SITE  km
*    WCrossR     - CALL CROSS product result      km / s
*    EarthRate   - IJK Earth's rotation rate vec  rad / s
*    TempVec     - Temporary vector
*    Temp        - Temporary REAL*8 value
*    Temp1       - Temporary REAL*8 value
*    i           - Index
*
*  Coupling      :
*    MAG         - Magnitude of a vector
*    ADDVEC      - Add two vectors
*    CROSS       - CROSS product of two vectors
*    ROT3        - Rotation about the 3rd axis
*    ROT2        - Rotation about the 2nd axis
*    DOT         - DOT product of two vectors
*    RVSEZ_RAZEL - Find R and V from SITE in Topocentric Horizon (SEZ) system
*    LNCOM2      - Combine two vectors and constants
*
*  References    :
*    Vallado       2001, 250-255, Alg 27
*
* ------------------------------------------------------------------------------

      SUBROUTINE RV_RAZEL    ( Reci,Veci,Latgd,Lon,alt,TTT,jdut1,lod,
     &                         xp,yp,terms, Direction,
     &                         Rho,Az,El,DRho,DAz,DEl )
        IMPLICIT NONE
        REAL*8 Reci(3),Veci(3),Latgd,Lon,Alt,Rho,Az,El,DRho,DAz,DEl,
     &         TTT,jdut1,lod,xp,yp
        INTEGER terms
        CHARACTER*4 Direction
        EXTERNAL Dot, MAG
* -----------------------------  Locals  ------------------------------
        REAL*8 Temp, Temp1, Rhoecef(3), magrhosez, RSecef(3),VSecef(3),
     &         DRhoecef(3), Rhosez(3), DRhosez(3), WCrossR(3),
     &         TempVec, Dot, MAG, Lat,
     &         recef(3), vecef(3)
        INTEGER i

        INCLUDE 'astmath.cmn'
        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
        CALL SITE ( Latgd,Alt,Lon, RSecef,VSecef )

        IF ( Direction .eq. 'FROM' ) THEN
            ! --------  Find SEZ range and velocity vectors -----------
            CALL RVSEZ_RAZEL( Rhosez,DRhosez, 'FROM',
     &                        Rho,Az,El,DRho,DAz,DEl )

            ! ---------  Perform SEZ to ECEF transformation -----------
            CALL ROT2( Rhosez ,Latgd-HalfPi, TempVec )
            CALL ROT3( TempVec,   -Lon   , Rhoecef )
            CALL ROT2( DRhosez,Latgd-HalfPi, TempVec )
            CALL ROT3( TempVec,   -Lon   , vecef )

            ! -----------  Find range and velocity vectors ------------
            CALL ADDVEC( Rhoecef,RSecef,Recef )

            CALL GCRF_ITRF ( reci,veci, 'FROM', rECEF,vECEF,
     &                      TTT, JDUT1, LOD, xp, yp, terms )
          ELSE

            ! ---------------- convert eci to ecef --------------------
            CALL GCRF_ITRF  ( reci,veci, 'TOO ', rECEF,vECEF,
     &                      TTT, JDUT1, LOD, xp, yp, terms )

            ! ----- find ecef range vector from site to satellite -----
            CALL SUBVEC( recef, rsecef,  rhoecef)
            rho = mag(rhoecef)

            ! ----------- Convert to SEZ for calculations -------------
            CALL ROT3( Rhoecef,    Lon   ,  TempVec )
            CALL ROT2( TempVec,HalfPi-Latgd,   Rhosez   )
            CALL ROT3( vecef  ,    Lon   ,  TempVec )
            CALL ROT2( TempVec,HalfPi-Latgd,  DRhosez   )

            ! ----------- Calculate Azimuth and Elevation -------------
            Temp= DSQRT( Rhosez(1)*Rhosez(1) + Rhosez(2)*Rhosez(2) )
            IF ( Temp .lt. Small ) THEN
                Az = DATAN2( DRhosez(2) , -DRhosez(1) )
              ELSE
                Az = DATAN2( Rhosez(2) , -Rhosez(1) )
              ENDIF

            IF ( ( Temp .lt. Small ) ) THEN   ! directly over the north pole
                El= DSIGN(1.0D0, Rhosez(3))*HalfPi ! +- 90
              ELSE
                magrhosez = MAG(rhosez)
                El= DASIN( Rhosez(3) / magRhosez )
              ENDIF

            ! ---- Calculate Range, Azimuth and Elevation rates -------
            DRho= DOT(Rhosez,DRhosez)/Rho
            IF ( DABS( Temp*Temp ) .gt. Small ) THEN
                DAz= ( DRhosez(1)*Rhosez(2) - DRhosez(2)*Rhosez(1) ) /
     &                 ( Temp*Temp )
              ELSE
                DAz= 0.0D0
              ENDIF

            IF ( DABS( Temp ) .gt. 0.00000001D0 ) THEN
                DEl= ( DRhosez(3) - DRho*DSIN( El ) ) / Temp
              ELSE
                DEl= 0.0D0
              ENDIF

          ENDIF  ! If
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE RV_ELATLON
*
*  this subroutine converts ecliptic latitude and longitude with position .and.
*    velocity vectors. Uses velocity vector to find the solution of singular
*    cases.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Rijk        - IJK position vector            km
*    Vijk        - IJK velocity vector            km/s
*    Direction   - Which set of vars to output    FROM  TOO
*
*  OutPuts       :
*    Rr          - Radius of the sat              km
*    EclLat      - Ecliptic Latitude              -Pi/2 to Pi/2 rad
*    EclLon      - Ecliptic Longitude             -Pi/2 to Pi/2 rad
*    DRr         - Radius of the sat rate         km/s
*    DEclLat     - Ecliptic Latitude rate         -Pi/2 to Pi/2 rad
*    EEclLon     - Ecliptic Longitude rate        -Pi/2 to Pi/2 rad
*
*  Locals        :
*    Obliquity   - Obliquity of the ecliptic      rad
*    Temp        -
*    Temp1       -
*    Re          - Position vec in eclitpic frame
*    Ve          - Velocity vec in ecliptic frame
*
*  Coupling      :
*    ROT1        - Rotation about 1st axis
*    DOT         - DOT product
*
*  References    :
*    Vallado       2001, 257-259, Eq 4-15
*
* ------------------------------------------------------------------------------

      SUBROUTINE RV_ELATLON  ( Rijk,Vijk, Direction, rr,EclLat,EclLon,
     &                         DRr,DEclLat,DEclLon )
        IMPLICIT NONE
        REAL*8 Rijk(3), Vijk(3), rr,EclLat,EclLon,DRr,DEclLat,DEclLon
        CHARACTER*4 Direction
        EXTERNAL DOT, MAG
* -----------------------------  Locals  ------------------------------
        REAL*8 Dot, Small, Re(3), Ve(3), Obliquity, Temp, Temp1, Mag

        ! --------------------  Implementation   ----------------------
        Small    = 0.00000001D0
        Obliquity= 0.40909280D0  !23.439291D0/rad
        IF ( Direction .eq. 'FROM' ) THEN
            Re(1)= rr*DCOS(EclLat)*DCOS(EclLon)
            Re(2)= rr*DCOS(EclLat)*DSIN(EclLon)
            Re(3)= rr*DSIN(EclLat)

            Ve(1)= DRr*DCOS(EclLat)*DCOS(EclLon)
     &               - rr*DSIN(EclLat)*DCOS(EclLon)*DEclLat
     &               - rr*DCOS(EclLat)*DSIN(EclLon)*DEclLon
            Ve(2)= DRr*DCOS(EclLat)*DSIN(EclLon)
     &               - rr*DSIN(EclLat)*DSIN(EclLon)*DEclLat
     &               + rr*DCOS(EclLat)*DCOS(EclLon)*DEclLon
            Ve(3)= DRr*DSIN(EclLat) + rr*DCOS(EclLat)*DEclLat

            CALL ROT1( Re, -Obliquity, Rijk )
            CALL ROT1( Ve, -Obliquity, Vijk )
          ELSE
            CALL ROT1( Rijk, Obliquity, Re )
            CALL ROT1( Vijk, Obliquity, Ve )

            ! ------------- Calculate Angles and Rates ----------------
            rr= MAG(Re)
            Temp= DSQRT( Re(1)*Re(1) + Re(2)*Re(2) )
            IF ( Temp .lt. Small ) THEN
                Temp1= DSQRT( Ve(1)*Ve(1) + Ve(2)*Ve(2) )
                IF ( DABS(Temp1) .gt. Small ) THEN
                    EclLon= DATAN2( Ve(2) , Ve(1) )
                  ELSE
                    EclLon= 0.0D0
                  ENDIF
              ELSE
                EclLon= DATAN2( Re(2) , Re(1) )
              ENDIF
            EclLat= DASIN( Re(3)/rr )

            Temp1= -Re(2)*Re(2) - Re(1)*Re(1)  ! different now
            DRr= DOT(re,Ve)/rr
            IF ( DABS( Temp1 ) .gt. Small ) THEN
                DEclLon= ( Ve(1)*Re(2) - Ve(2)*Re(1) ) / Temp1
              ELSE
                DEclLon= 0.0D0
              ENDIF
            IF ( DABS( Temp ) .gt. Small ) THEN
                DEclLat= ( Ve(3) - DRr*DSIN( EclLat ) ) / Temp
              ELSE
                DEclLat= 0.0D0
              ENDIF
          ENDIF  ! IF (

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE RVSEZ_RAZEL
*
*  this subroutine converts range, azimuth, and elevation values with slant
*    range and velocity vectors for a satellite from a radar SITE in the
*    Topocentric Horizon (SEZ) system.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    RhoVec      - SEZ Satellite range vector     km
*    DRhoVec     - SEZ Satellite velocity vector  km / s
*    Direction   - Which set of vars to output    FROM  TOO
*
*  OutPuts       :
*    Rho         - Satellite range from SITE      km
*    Az          - Azimuth                        0.0D0 to 2Pi rad
*    El          - Elevation                      -Pi/2 to Pi/2 rad
*    DRho        - Range Rate                     km / s
*    DAz         - Azimuth Rate                   rad / s
*    DEl         - Elevation rate                 rad / s
*
*  Locals        :
*    SinEl       - Variable for DSIN( El )
*    CosEl       - Variable for DCOS( El )
*    SinAz       - Variable for DSIN( Az )
*    CosAz       - Variable for DCOS( Az )
*    Temp        -
*    Temp1       -
*
*  Coupling      :
*    DOT         - DOT product
*
*  References    :
*    Vallado       2001, 250-251, Eq 4-4, Eq 4-5
*
* ------------------------------------------------------------------------------

      SUBROUTINE RVSEZ_RAZEL ( Rhosez,DRhosez,Direction, Rho,Az,El,
     &           DRho,DAz,DEl )
        IMPLICIT NONE
        CHARACTER*4 Direction
        REAL*8 RhoSez(3), DRhoSez(3),Rho,Az,El,DRho,DAz, DEl, MAG
        EXTERNAL Dot, MAG

        INCLUDE 'astmath.cmn'

* -----------------------------  Locals  ------------------------------
        REAL*8 Temp1, Temp, SinEl, CosEl, SinAz,CosAz, Dot, magrhosez

        ! --------------------  Implementation   ----------------------
        IF ( Direction .eq. 'FROM' ) THEN
            ! -------------------- Initialize values ------------------
            SinEl= DSIN(El)
            CosEl= DCOS(El)
            SinAz= DSIN(Az)
            CosAz= DCOS(Az)

            ! ----------------- Form SEZ range vector -----------------
            Rhosez(1) = -Rho*CosEl*CosAz
            Rhosez(2) =  Rho*CosEl*SinAz
            Rhosez(3) =  Rho*SinEl

            ! --------------- Form SEZ velocity vector ----------------
            DRhosez(1) = -DRho*CosEl*CosAz + Rhosez(3)*DEl*CosAz +
     &                     Rhosez(2)*DAz
            DRhosez(2) =  DRho*CosEl*SinAz - Rhosez(3)*DEl*SinAz -
     &                     Rhosez(1)*DAz
            DRhosez(3) =  DRho*SinEl       + Rho*DEl*CosEl
          ELSE
            ! ----------- Calculate Azimuth and Elevation -------------
            Temp= DSQRT( Rhosez(1)*Rhosez(1) + Rhosez(2)*Rhosez(2) )
            IF ( DABS( Rhosez(2) ) .lt. Small ) THEN
                IF ( Temp .lt. Small ) THEN
                    Temp1= DSQRT( DRhosez(1)*DRhosez(1) +
     &                     DRhosez(2)*DRhosez(2) )
                    Az   =  DATAN2( DRhosez(2)/Temp1 ,
     &                     -DRhosez(1)/Temp1 )
                  ELSE
                    IF ( Rhosez(1) .gt. 0.0D0 ) THEN
                        Az= Pi
                      ELSE
                        Az= 0.0D0
                      ENDIF
                  ENDIF
              ELSE
                Az= DATAN2( Rhosez(2)/Temp , -Rhosez(1)/Temp )
              ENDIF

            IF ( ( Temp .lt. Small ) ) THEN   ! directly over the north pole
                El= DSIGN(1.0D0,Rhosez(3))*HalfPi ! +- 90
              ELSE
                El= DASIN( Rhosez(3) / mag(Rhosez) )
              ENDIF

            ! -----  Calculate Range, Azimuth and Elevation rates -----
            DRho= DOT(Rhosez,DRhosez)/Rho 
            IF ( DABS( Temp*Temp ) .gt. Small ) THEN
                DAz= ( DRhosez(1)*Rhosez(2) - DRhosez(2)*Rhosez(1) ) /
     &                 ( Temp*Temp )
              ELSE
                DAz= 0.0D0
              ENDIF

            IF ( DABS( Temp ) .gt. Small ) THEN
                DEl= ( DRhosez(3) - DRho*DSIN( El ) ) / Temp
              ELSE
                DEl= 0.0D0 
              ENDIF
          ENDIF
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE RADEC_ELATLON
*
*  this subroutine converts right-ascension declination values with ecliptic
*    latitude and longitude values.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    RtAsc       - Right Ascension                rad
*    Decl        - Declination                    rad
*    Direction   - Which set of vars to output    FROM  TOO
*
*  OutPuts       :
*    EclLat      - Ecliptic Latitude              -Pi/2 to Pi/2 rad
*    EclLon      - Ecliptic Longitude             -Pi/2 to Pi/2 rad
*
*  Locals        :
*    Obliquity   - Obliquity of the ecliptic      rad
*    Sinv        -
*    Cosv        -
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2001, 259, Eq 4-19, Eq 4-20
*
* ------------------------------------------------------------------------------  

      SUBROUTINE RADEC_ELATLON ( RtAsc,Decl,Direction, EclLat, EclLon )
        IMPLICIT NONE
        REAL*8 RtAsc,Decl,EclLat,EclLon
        CHARACTER*4 Direction
* -----------------------------  Locals  ------------------------------
        REAL*8 Sinv, Cosv, Obliquity

        ! --------------------  Implementation   ----------------------
        Obliquity= 0.40909280D0  !23.439291D0/rad
        IF ( Direction .eq. 'FROM' ) THEN
            Decl = DASIN( DSIN(EclLat)*DCOS(Obliquity)
     &                    + DCOS(EclLat)*DSIN(Obliquity)*DSIN(EclLon) )
            Sinv = ( -DSIN(EclLat)*DSIN(Obliquity)
     &                 + DCOS(EclLat)*DCOS(Obliquity)*DSIN(EclLon) ) /
     &                 DCOS(Decl)
            Cosv = DCOS(EclLat)*DCOS(EclLon) / DCOS(Decl) 
            RtAsc= DATAN2( Sinv,Cosv ) 
          ELSE
            EclLat= DASIN( -DCOS(Decl)*DSIN(RtAsc)*DSIN(Obliquity)
     &                        + DSIN(Decl)*DCOS(Obliquity) )
            Sinv  = ( DCOS(Decl)*DSIN(RtAsc)*DCOS(Obliquity)
     &                  + DSIN(Decl)*DSIN(Obliquity) ) / DCOS(EclLat)
            Cosv  = DCOS(Decl)*DCOS(RtAsc) / DCOS(EclLat) 
            EclLon= DATAN2( Sinv,Cosv ) 
          ENDIF 
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE RADEC_AZEL
*
* this subroutine converts right ascension declination values with
*   azimuth, and elevation.  Notice the range is not defined because
*   Right ascension declination only allows a unit vector to be formed.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    RtAsc       - Right Ascension                0.0D0 to 2Pi rad
*    Decl        - Declination                    -Pi/2 to Pi/2 rad
*    LST         - Local SIDEREAL Time            -2Pi to 2Pi rad
*    Latgd       - Geodetic Latitude              -Pi/2 to Pi/2 rad
*    Direction   - Which set of vars to output    FROM  TOO
*
*  Outputs       :
*    Az          - Azimuth                        0.0D0 to 2Pi rad
*    El          - Elevation                      -Pi/2 to Pi/2 rad
*
*  Locals        :
*    LHA         - Local Hour ANGLE               -2Pi to 2Pi rad
*    Sinv        - Sine value
*    Cosv        - Cosine value
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2001, 255-257, Alg 28
*
* ------------------------------------------------------------------------------  

      SUBROUTINE RADEC_AZEL  ( RtAsc,Decl,LST,Latgd, Direction, Az,El )
        IMPLICIT NONE
        REAL*8 RtAsc,Decl,LST,Latgd,Az,El
        CHARACTER*4 Direction
* -----------------------------  Locals  ------------------------------
        REAL*8 Sinv, Cosv, LHA

        ! --------------------  Implementation   ----------------------
        IF ( Direction .eq. 'FROM' ) THEN
            Decl = DASIN( DSIN(El)*DSIN(Latgd) +
     &                 DCOS(el)*DCOS(Latgd)*DCOS(Az) )

            Sinv = -(DSIN(az)*DCOS(el)*DCOS(Latgd)) /
     &              (DCOS(Latgd)*DCOS(Decl))
            Cosv = (DSIN(el) - DSIN(Latgd)*DSIN(decl)) /
     &              (DCOS(Latgd)*DCOS(Decl))
            LHA  = DATAN2( Sinv,Cosv ) 
            RtAsc= LST - LHA 
          ELSE
            LHA = LST - RtAsc

            El  = DASIN( DSIN(Decl)*DSIN(Latgd) +
     &            DCOS(Decl)*DCOS(Latgd)*DCOS(LHA) )

            Sinv= -DSIN(LHA)*DCOS(Decl)*DCOS(Latgd)/
     &                (DCOS(el)*DCOS(Latgd))
            Cosv= ( DSIN(Decl)-DSIN(el)*DSIN(Latgd) )/
     &             (DCOS(el)*DCOS(Latgd))
            Az  = DATAN2( Sinv,Cosv ) 
          ENDIF 

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE GIBBS
*
*  this subroutine performs the GIBBS method of orbit determination.  This
*    method determines the velocity at the middle point of the 3 given position
*    vectors.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R1          - IJK Position vector #1         km
*    R2          - IJK Position vector #2         km
*    R3          - IJK Position vector #3         km
*
*  OutPuts       :
*    V2          - IJK Velocity Vector for R2     km / s
*    Theta       - ANGLE between vectors          rad
*    Error       - Flag indicating success        'ok',...
*
*  Locals        :
*    tover2      -
*    l           -
*    Small       - Tolerance for roundoff errors
*    r1mr2       - Magnitude of r1 - r2
*    r3mr1       - Magnitude of r3 - r1
*    r2mr3       - Magnitude of r2 - r3
*    p           - P Vector     r2 x r3
*    q           - Q Vector     r3 x r1
*    w           - W Vector     r1 x r2
*    d           - D Vector     p + q + w
*    n           - N Vector (r1)p + (r2)q + (r3)w
*    s           - S Vector
*                    (r2-r3)r1+(r3-r1)r2+(r1-r2)r3
*    b           - B Vector     d x r2
*    Theta1      - Temp ANGLE between the vectors rad
*    Pn          - P Unit Vector
*    R1N         - R1 Unit Vector
*    dn          - D Unit Vector
*    Nn          - N Unit Vector
*    i           - index
*
*  Coupling      :
*    MAG         - Magnitude of a vector
*    CROSS       - CALL CROSS product of two vectors
*    DOT         - DOT product of two vectors
*    ADD3VEC     - Add three vectors
*    LNCOM2      - Multiply two vectors by two constants
*    LNCOM3      - Add three vectors each multiplied by a constant
*    NORM        - Creates a Unit Vector
*    ANGLE       - ANGLE between two vectors
*
*  References    :
*    Vallado       2001, 432-445, Alg 52, Ex 7-5
*
* ------------------------------------------------------------------------------  

      SUBROUTINE GIBBS       ( R1,R2,R3, V2, Theta,Theta1,Copa, Error )
        IMPLICIT NONE
        REAL*8 R1(3), R2(3), R3(3), V2(3), Theta, Theta1, Copa
        CHARACTER*12 Error
        EXTERNAL DOT, MAG
* -----------------------------  Locals  ------------------------------
        INTEGER i
        REAL*8 tover2, l, Small, r1mr2, r3mr1, r2mr3, p(3), q(3), w(3),
     &         d(3), n(3), s(3), b(3), Pn(3), R1N(3), Dn(3), Nn(3),Dot,
     &         magr1, magr2, magr3, mag, magd, magn

        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
        Small= 0.000001D0
        Theta= 0.0D0
        Error = 'ok'
        Theta1= 0.0D0
        magr1 = MAG( R1 )
        magr2 = MAG( R2 )
        magr3 = MAG( R3 )
        DO i= 1 , 3
            V2(i)= 0.0D0
          ENDDO

        CALL CROSS( R2,R3,P )
        CALL CROSS( R3,R1,Q )
        CALL CROSS( R1,R2,W )
        CALL NORM( P,Pn )
        CALL NORM( R1,R1N )
        Copa=  DASIN( DOT( Pn,R1n ) ) 

        IF ( DABS( Copa ) .gt. 0.017452406D0 ) THEN
            Error= 'not coplanar'
          ENDIF

        ! --------------- .or. don't contiNue processing --------------
        CALL ADD3VEC( P,Q,W,D )
        CALL LNCOM3( magr1,magr2,magr3,P,Q,W,N )
        CALL NORM( N,Nn )
        CALL NORM( D,DN )
        magd = MAG(d)
        magn = MAG(n)

        ! -------------------------------------------------------------
*       Determine If  the orbit is possible.  Both D and N must be in
*         the same direction, and non-zero.
        ! -------------------------------------------------------------
        IF ( ( DABS(magd).lt.Small ) .or. ( DABS(magn).lt.Small ) .or.
     &      ( DOT(Nn,dn) .lt. Small ) ) THEN
            Error= 'impossible'
          ELSE
              CALL ANGLE( R1,R2, Theta )
              CALL ANGLE( R2,R3, Theta1 )

              ! ----------- Perform GIBBS method to find V2 -----------
              R1mr2= magr1-magr2
              R3mr1= magr3-magr1
              R2mr3= magr2-magr3
              CALL LNCOM3(R1mr2,R3mr1,R2mr3,R3,R2,R1,S)
              CALL CROSS( d,r2,b )
              L    = DSQRT( mu / (magd*magn) )
              Tover2= L / magr2
              CALL LNCOM2(Tover2,L,B,S,V2)
            ENDIF

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE HERRGIBBS
*
*  this subroutine implements the Herrick-GIBBS approximation for orbit
*    determination, and finds the middle velocity vector for the 3 given
*    position vectors.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R1          - IJK Position vector #1         km
*    R2          - IJK Position vector #2         km
*    R3          - IJK Position vector #3         km
*    JD1         - Julian Date of 1st sighting    days from 4713 BC
*    JD2         - Julian Date of 2nd sighting    days from 4713 BC
*    JD3         - Julian Date of 3rd sighting    days from 4713 BC
*
*  OutPuts       :
*    V2          - IJK Velocity Vector for R2     km / s
*    Theta       - ANGLE between vectors          rad
*    Error       - Flag indicating success        'ok',...
*
*  Locals        :
*    Dt21        - time delta between r1 and r2   s
*    Dt31        - time delta between r3 and r1   s
*    Dt32        - time delta between r3 and r2   s
*    p           - P vector    r2 x r3
*    Pn          - P Unit Vector
*    R1N         - R1 Unit Vector
*    Theta1      - temporary ANGLE between vec    rad
*    TolAngle    - Tolerance ANGLE  (1 deg)       rad
*    Term1       - 1st Term for HGibbs expansion
*    Term2       - 2nd Term for HGibbs expansion
*    Term3       - 3rd Term for HGibbs expansion
*    i           - Index
*
*  Coupling      :
*    MAG         - Magnitude of a vector
*    CROSS       - CALL CROSS product of two vectors
*    DOT         - DOT product of two vectors
*    NORM        - Creates a Unit Vector
*    LNCOM3      - Combination of three scalars and three vectors
*    ANGLE       - ANGLE between two vectors
*
*  References    :
*    Vallado       2001, 439-445, Alg 52, Ex 7-4
*
* ------------------------------------------------------------------------------  

      SUBROUTINE HERRGIBBS   ( R1,R2,R3,JD1,JD2,JD3, V2, Theta,Theta1,
     &                         Copa, Error )
        IMPLICIT NONE
        REAL*8 R1(3), R2(3), R3(3), JD1, JD2, JD3, V2(3), Theta,
     &         Theta1, Copa
        CHARACTER*12 Error
        EXTERNAL Dot
* -----------------------------  Locals  ------------------------------
        INTEGER i
        REAL*8 p(3), Pn(3), R1n(3),Dot, magr1, magr2, magr3,
     &         Dt21, Dt31, Dt32, Term1, Term2, Term3, TolAngle, mag

        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
        Error =  'ok'
        Theta = 0.0D0
        Theta1= 0.0D0
        magr1 = MAG( R1 )
        magr2 = MAG( R2 )
        magr3 = MAG( R3 )
        DO i= 1 , 3
            V2(i)= 0.0D0
          ENDDO
        TolAngle= 0.01745329251994D0
        Dt21= (JD2-JD1)*86400.0D0
        Dt31= (JD3-JD1)*86400.0D0   ! differences in times
        Dt32= (JD3-JD2)*86400.0D0

        CALL CROSS( R2,R3,P )
        CALL NORM( P,Pn )
        CALL NORM( R1,R1N )
        Copa=  DASIN( DOT( Pn,R1n ) )
        IF ( DABS( Copa ) .gt. 0.017452406D0 ) THEN
            Error= 'not coplanar'
          ENDIF

        ! --------------------------------------------------------------
*       Check the size of the angles between the three position vectors.
*       Herrick GIBBS only gives "reasonable" answers when the
*       position vectors are reasonably close.  10 deg is only an estimate.
        ! --------------------------------------------------------------
        CALL ANGLE( R1,R2, Theta )
        CALL ANGLE( R2,R3, Theta1 )
        IF ( (Theta .gt. TolAngle) .or. (Theta1 .gt. TolAngle) ) THEN
            Error= 'ANGLE > 1ø'
          ENDIF

        ! ----------- Perform Herrick-GIBBS method to find V2 ---------
        Term1= -Dt32*( 1.0D0/(Dt21*Dt31) +
     &        mu/(12.0D0*magr1*magr1*magr1) )
        Term2= (Dt32-Dt21)*( 1.0D0/(Dt21*Dt32) +
     &        mu/(12.0D0*magr2*magr2*magr2) )
        Term3=  Dt21*( 1.0D0/(Dt32*Dt31) +
     &         mu/(12.0D0*magr3*magr3*magr3) )
        CALL LNCOM3( Term1,Term2,Term3,R1,R2,R3, V2 )

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE LAMBERTUNIV
*
*  this subroutine solves the Lambert problem for orbit determination and returns
*    the velocity vectors at each of two given position vectors.  The solution
*    uses Universal Variables for calculation and a bissection technique
*    updating psi.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R1          - IJK Position vector 1          km
*    R2          - IJK Position vector 2          km
*    DM          - direction of motion            'L','S'
*    Dtsec        - Time between R1 and R2         s
*
*  OutPuts       :
*    V1          - IJK Velocity vector            km / s
*    V2          - IJK Velocity vector            km / s
*    Error       - Error flag                     'ok', ...
*
*  Locals        :
*    VarA        - Variable of the iteration,
*                  NOT the semi .or. axis!
*    Y           - Area between position vectors
*    Upper       - Upper bound for Z
*    Lower       - Lower bound for Z
*    CosDeltaNu  - Cosine of true anomaly change  rad
*    F           - f expression
*    G           - g expression
*    GDot        - g DOT expression
*    XOld        - Old Universal Variable X
*    XOldCubed   - XOld cubed
*    ZOld        - Old value of z
*    ZNew        - New value of z
*    C2New       - C2(z) FUNCTION
*    C3New       - C3(z) FUNCTION
*    TimeNew     - New time                       s
*    Small       - Tolerance for roundoff errors
*    i, j        - index
*
*  Coupling      :
*    MAG         - Magnitude of a vector
*    DOT         - DOT product of two vectors
*    FINDC2C3    - Find C2 and C3 functions
*
*  References    :
*    Vallado       2001, 459-464, Alg 55, Ex 7-5
*
* ------------------------------------------------------------------------------

      SUBROUTINE LAMBERTUNIV ( ro,r, dm,OverRev, Dtsec, vo,v, Error )
        IMPLICIT NONE
        REAL*8 Ro(3), r(3), Dtsec, Vo(3),v(3)
        CHARACTER Dm, OverRev
        CHARACTER*12 Error
        EXTERNAL DOT, MAG

        INCLUDE 'astmath.cmn'
        INCLUDE 'astconst.cmn'

* -----------------------------  Locals  ------------------------------
        INTEGER i, Loops, YNegKtr, NumIter
        REAL*8 VarA, Y, Upper, Lower,CosDeltaNu, magro, magr,
     &         F, G, GDot, XOld, XOldCubed, Dot, mag,
     &         PsiOld, PsiNew, C2New, C3New, dtNew

        ! --------------------  Implementation   ----------------------
        NumIter= 40
        Error  = 'ok' 
        PsiNew = 0.0D0 
        magro = MAG(ro)
        magr = MAG(r)
        DO i= 1 , 3
            vo(i)= 0.0D0
            v(i)= 0.0D0
          ENDDO

        CosDeltaNu= DOT(ro,r)/(magro*magr)
        IF ( Dm .eq. 'L' ) THEN
            VarA = -DSQRT( magro*magr*(1.0D0+CosDeltaNu) )
          ELSE
            VarA =  DSQRT( magro*magr*(1.0D0+CosDeltaNu) )
          ENDIF

        ! ---------------  Form Initial guesses   ---------------------
        PsiOld = 0.0D0 
        PsiNew = 0.0D0 
        xOld   = 0.0D0 
        C2New  = 0.5D0 
        C3New  = 1.0D0/6.0D0

        ! --------- Set up initial bounds for the bissection ----------
        IF ( OverRev .eq. 'N' ) THEN
            Upper= TwoPi*TwoPi
            Lower= -4.0D0*TwoPi 
          ELSE
            Upper= -0.001D0+4.0D0*TwoPi*TwoPi ! at 4, not alw work
            Lower=  0.001D0+TwoPi*TwoPi       ! 2.0D0, makes orbit bigger
          ENDIF                               ! how about 2 revs??xx

        ! -------  Determine If  the orbit is possible at all ---------
        IF ( DABS( VarA ) .gt. Small ) THEN

            Loops  = 0 
            YNegKtr= 1  ! y neg ktr
            DtNew = -10.0D0
            DO WHILE ((DABS(dtNew-Dtsec) .ge. 0.000001D0) .and.
     &              (Loops .lt. NumIter) .and. (YNegKtr .le. 10))
                IF ( DABS(C2New) .gt. Small ) THEN
                    Y= magro + magr -
     &                      ( VarA*(1.0D0-PsiOld*C3New)/DSQRT(C2New) )
                  ELSE
                    Y= magro + magr
                  ENDIF
                ! ----------- Check for negative values of y ----------
                IF ( ( VarA .gt. 0.0D0 ) .and. ( Y .lt. 0.0D0 ) ) THEN
                    YNegKtr= 1
                    DO WHILE (( Y.lt.0.0D0 ) .and. ( YNegKtr .lt. 10 ))
                        PsiNew= 0.8D0*(1.0D0/C3New)*( 1.0D0
     &                           - (magro+magr)*DSQRT(C2New)/VarA  )
                        ! -------- Find C2 and C3 functions -----------
                        CALL FINDC2C3( PsiNew, C2New,C3New )
                        PsiOld= PsiNew 
                        Lower= PsiOld 
                        IF ( DABS(C2New) .gt. Small ) THEN
                            Y= magro + magr - ( VarA*(1.0D0-
     &                               PsiOld*C3New)/DSQRT(C2New) )
                          ELSE
                            Y= magro + magr
                          ENDIF
                        YNegKtr = YNegKtr + 1
                      ENDDO ! while
                  ENDIF  ! If  y neg

                IF ( YNegKtr .lt. 10 ) THEN
                    IF ( DABS(C2New) .gt. Small ) THEN
                        XOld= DSQRT( Y/C2New )
                      ELSE
                        XOld= 0.0D0
                      ENDIF
                    XOldCubed= XOld*XOld*XOld 
                    dtNew    = (XOldCubed*C3New + VarA*DSQRT(Y)) /
     &                          DSQRT(mu)

                    ! --------  Readjust upper and lower bounds -------
                    IF ( dtNew .lt. Dtsec ) THEN
                        Lower= PsiOld 
                      ENDIF
                    IF ( dtNew .gt. Dtsec ) THEN
                        Upper= PsiOld 
                      ENDIF
                    PsiNew= (Upper+Lower) * 0.5D0 

                    ! ------------- Find c2 and c3 functions ----------
                    CALL FINDC2C3( PsiNew, C2New,C3New )
                    PsiOld = PsiNew
                    Loops = Loops + 1

                    ! --- Make sure the first guess isn't too close ---
                    IF ( (DABS(dtNew - Dtsec) .lt. Small)
     &                        .and. (Loops .eq. 1) ) THEN
                        dtNew= Dtsec-1.0D0
                      ENDIF
                  ENDIF  ! If  YNegKtr .lt. 10
c               write(20,'(4f14.6)') y,xold,dtnew,psinew
              ENDDO ! Do While Loop

            IF ( (Loops .ge. NumIter) .or. (YNegKtr .ge. 10) ) THEN
                Error= 'GNotConv'
                IF ( YNegKtr .ge. 10 ) THEN
                    Error= 'Y negative' 
                  ENDIF
              ELSE
                ! --- Use F and G series to find Velocity Vectors -----
                F   = 1.0D0 - Y/magro
                GDot= 1.0D0 - Y/magr
                G   = 1.0D0 / (VarA*DSQRT( Y/mu ))  ! 1 over G
                DO i= 1 , 3
                    vo(i)= ( r(i) - F*ro(i) )*G
                    v(i) = ( GDot*r(i) - ro(i) )*G
                  ENDDO
              ENDIF   ! If  the answer has converged
          ELSE
            Error= 'impos 180ø' 
          ENDIF  ! If  Var A .gt. 0.0D0

      RETURN
      END

* --------- Two recursion algorithms needed by the LambertBattin routine
*
      REAL*8 FUNCTION SEE    ( v )
        IMPLICIT NONE
        REAL*8 v
* -----------------------------  Locals  ------------------------------
        REAL*8 c(0:20),term, termold, del, delold, sum1, eta, SQRTopv
        INTEGER i

        ! --------------------  Implementation   ----------------------
          c(0) =    0.2D0
          c(1) =    9.0D0 /  35.0D0
          c(2) =   16.0D0 /  63.0D0
          c(3) =   25.0D0 /  99.0D0
          c(4) =   36.0D0 / 143.0D0
          c(5) =   49.0D0 / 195.0D0
          c(6) =   64.0D0 / 255.0D0
          c(7) =   81.0D0 / 323.0D0
          c(8) =  100.0D0 / 399.0D0
          c(9) =  121.0D0 / 483.0D0
          c(10)=  144.0D0 / 575.0D0
          c(11)=  169.0D0 / 675.0D0
          c(12)=  196.0D0 / 783.0D0
          c(13)=  225.0D0 / 899.0D0
          c(14)=  256.0D0 /1023.0D0
          c(15)=  289.0D0 /1155.0D0
          c(16)=  324.0D0 /1295.0D0
          c(17)=  361.0D0 /1443.0D0
          c(18)=  400.0D0 /1599.0D0
          c(19)=  441.0D0 /1763.0D0
          c(20)=  484.0D0 /1935.0D0
          SQRTOpv= DSQRT(1.0D0 + v) 
          eta    = v / ( ( 1.0D0+SQRTOpv )**2 )

          ! ------------------- Process Forwards ----------------------
          delold = 1.0D0
          termold= c(0)   ! * eta}
          sum1   = termold 
          i= 1 
          DO WHILE ((i.le.20) .and. (DABS(Termold) .gt. 0.000001D0 ))
              del  = 1.0D0 / ( 1.0D0 + c(i)*eta*delold )
              term = termold * (del - 1.0D0) 
              sum1 = sum1 + term 
              i    = i + 1
              delold = del
              termold= term 
            ENDDO

c          See= (1.0D0 / (8.0D0*(1.0D0+SQRTOpv))) *
c     &         ( 3.0D0 + Sum1 / ( 1.0D0+eta*sum1 ) )
          See= 1.0D0 / ((1.0D0/(8.0D0*(1.0D0+sqrtopv))) *
     &         ( 3.0D0 + sum1 / ( 1.0D0+eta*sum1 ) ) )
      RETURN
      END   ! Internal FUNCTION See


      REAL*8 FUNCTION k      ( v )
        IMPLICIT NONE
        REAL*8 v
* -----------------------------  Locals  ------------------------------
        REAL*8 d(0:20), del,delold,term,termold, sum1
        INTEGER i

        ! --------------------  Implementation   ----------------------
          d(0) =     1.0D0 /    3.0D0
          d(1) =     4.0D0 /   27.0D0
          d(2) =     8.0D0 /   27.0D0
          d(3) =     2.0D0 /    9.0D0
          d(4) =    22.0D0 /   81.0D0
          d(5) =   208.0D0 /  891.0D0
          d(6) =   340.0D0 / 1287.0D0
          d(7) =   418.0D0 / 1755.0D0
          d(8) =   598.0D0 / 2295.0D0
          d(9) =   700.0D0 / 2907.0D0
          d(10)=   928.0D0 / 3591.0D0
          d(11)=  1054.0D0 / 4347.0D0
          d(12)=  1330.0D0 / 5175.0D0
          d(13)=  1480.0D0 / 6075.0D0
          d(14)=  1804.0D0 / 7047.0D0
          d(15)=  1978.0D0 / 8091.0D0
          d(16)=  2350.0D0 / 9207.0D0
          d(17)=  2548.0D0 /10395.0D0
          d(18)=  2968.0D0 /11655.0D0
          d(19)=  3190.0D0 /12987.0D0
          d(20)=  3658.0D0 /14391.0D0

          ! ----------------- Process Forwards ------------------------
          sum1   = d(0)
          delold = 1.0D0 
          termold= d(0)
          i      = 1
          DO WHILE ((i.le.20) .and. (DABS(Termold) .gt. 0.000001D0 ))
              del  = 1.0D0 / ( 1.0D0 + d(i)*v*delold )
              term = termold * ( del - 1.0D0 )
              sum1 = sum1 + term 
              i    = i + 1
              delold = del 
              termold= term 
            ENDDO

          k= Sum1
      RETURN
      END  ! Internal SUBROUTINE K

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE LAMBERBATTIN
*
*  this subroutine solves Lambert's problem using Battins method. The method is
*    developed in Battin (1987). It uses contiNued fractions to speed the
*    solution and has several parameters that are defined differently than
*    the traditional Gaussian technique.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Ro          - IJK Position vector 1          km
*    R           - IJK Position vector 2          km
*    DM          - direction of motion            'L','S'
*    Dtsec        - Time between R1 and R2         s
*
*  OutPuts       :
*    Vo          - IJK Velocity vector            km / s
*    V           - IJK Velocity vector            km / s
*    Error       - Error flag                     'ok',...
*
*  Locals        :
*    i           - Index
*    Loops       -
*    u           -
*    b           -
*    Sinv        -
*    Cosv        -
*    rp          -
*    x           -
*    xn          -
*    y           -
*    l           -
*    m           -
*    CosDeltaNu  -
*    SinDeltaNu  -
*    DNu         -
*    a           -
*    Tan2w       -
*    RoR         -
*    h1          -
*    h2          -
*    Tempx       -
*    eps         -
*    denom       -
*    chord       -
*    k2          -
*    s           -
*    f           -
*    g           -
*    gDot        -
*    am          -
*    ae          -
*    be          -
*    tm          -
*    arg1        -
*    arg2        -
*    tc          -
*    AlpE        -
*    BetE        -
*    AlpH        -
*    BetH        -
*    DE          -
*    DH          -
*
*  Coupling      :
*    MAG         - Magnitude of a vector
*    ASINH     - Inverse hyperbolic sine
*    ARCCOSH     - Inverse hyperbolic cosine
*    SINH        - Hyperbolic sine
*
*  References    :
*    Vallado       2001, 464-467, Ex 7-5
*
* ------------------------------------------------------------------------------  

      SUBROUTINE LAMBERTBATTIN ( ro,r, dm,OverRev, Dtsec, vo,v, Error )
        IMPLICIT NONE
        REAL*8 Ro(3), r(3), Dtsec, Vo(3),v(3)
        CHARACTER Dm, OverRev
        CHARACTER*12 Error
        EXTERNAL Dot, See, K, ASINH, MAG
* -----------------------------  Locals  ------------------------------
        INTEGER i, Loops
        REAL*8 RCrossR(3),Dot, See, K, ASINH, MAG,
*         lambda,bigt,     testamt,
     &   u, b, Sinv,Cosv, rp, x, xn, y, l, m, CosDeltaNu, SinDeltaNu,
     &   DNu, a, tan2w, RoR, h1, h2, tempx, eps, denom, chord, k2, s,
     &   f, g, am, ae, be, tm, gDot, arg1, arg2, AlpE, BetE,
     &   AlpH, BetH, DE, DH, magr, magro, magrcrossr,y1, lim1

        INCLUDE 'astmath.cmn'
        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
        Error = 'ok'
        magr = MAG(r)
        magro = MAG(ro)
        CosDeltaNu= DOT(ro,r)/(magro*magr)
        CALL CROSS( ro,r, RCrossR )
        magrcrossr = MAG(rcrossr)
        SinDeltaNu= magRCrossr/(magro*magr)
        DNu   = DATAN2( SinDeltaNu,CosDeltaNu )

        RoR   = magr/magro
        eps   = RoR - 1.0D0
        tan2w = 0.25D0*eps*eps / ( DSQRT( RoR ) + RoR *
     &                             ( 2.0D0 + DSQRT( RoR ) ) )
        rp    = DSQRT( magro*magr )*( (DCOS(DNu*0.25D0))**2 + tan2w )

        IF ( DNu .lt. Pi ) THEN
            L = ( (DSIN(DNu*0.25D0))**2 + tan2w ) /
     &          ( (DSIN(DNu*0.25D0))**2 + tan2w + DCOS( DNu*0.5D0 ) )
          ELSE
            L = ( (DCOS(DNu*0.25D0))**2 + tan2w - DCOS( DNu*0.5D0 ) ) /
     &            ( (DCOS(DNu*0.25D0))**2 + tan2w )
          ENDIF

        m    = mu*Dtsec*Dtsec / ( 8.0D0*rp*rp*rp )
        x    = 10.0D0
        xn   = L  !0.0D0   !L    ! 0 for par and hyp
        chord= DSQRT( magro*magro + magr*magr - 2.0D0*magro*magr*
     &         DCOS( DNu ) )
        s    = ( magro + magr + chord )*0.5D0
        lim1 = dsqrt(m/l)

        Loops= 1
        DO WHILE ((DABS(xn-x) .ge. Small) .and. (Loops .le. 30))
            x    = xn 
            Tempx= See(x) 
            Denom= 1.0D0 / ( (1.0D0+2.0D0*x+L) * (4.0D0*x +
     &             tempx*(3.0D0+x) ) )
            h1   = ( L+x )**2 * ( 1.0D0+ 3.0D0*x + Tempx )*Denom
            h2   = m*( x - L + Tempx )*Denom

            ! ----------------------- Evaluate CUBIC ------------------
            b = 0.25D0*27.0D0*h2 / ((1.0D0+h1)**3 )
            if (b .lt. -1.0D0) THEN ! reset the initial condition
                xn = 1.0D0 - 2.0D0*l
             else
                if (y1 .gt. lim1) THEN
                    xn = xn * (lim1/y1)
                else
                  u = 0.5D0*b / ( 1.0D0 + DSQRT( 1.0D0 + b ) ) 
                  K2= K(u) 

                   y = ( ( 1.0D0+h1 ) / 3.0D0 )*
     &                 ( 2.0D0 + DSQRT( 1.0D0+b ) /
     &                 ( 1.0D0+2.0D0*u*k2*k2 ) )
                   xn= DSQRT( ( (1.0D0-L)*0.5D0 )**2 + m/(y*y) ) -
     &                 ( 1.0D0+L )*0.5D0
                endif
             endif

            Loops = Loops + 1
        y1=dsqrt(m/((L+x)*(1.0D0+x)) )
        write(*,'(i3,6f11.7)') loops,y,x,k2,b,u,y1
          ENDDO

        a=  mu*Dtsec*Dtsec / (16.0D0*Rp*rp*xn*y*y )

        ! ------------------ Find Eccentric anomalies -----------------
        ! ------------------------ Hyperbolic -------------------------
        IF ( a .lt. -Small ) THEN
            arg1 = DSQRT( s / ( -2.0D0*a ) )
            arg2 = DSQRT( ( s-chord ) / ( -2.0D0*a ) )
            ! ------- Evaluate f and g functions --------
            AlpH = 2.0D0 * ASINH( arg1 )
            BetH = 2.0D0 * ASINH( arg2 )
            DH   = AlpH - BetH
            F    = 1.0D0 - (a/magro)*(1.0D0 - COSH(DH) )
            GDot = 1.0D0 - (a/magr) *(1.0D0 - COSH(DH) )
            G    = Dtsec - DSQRT(-a*a*a/mu)*(SINH(DH)-DH)
          ELSE
            ! ------------------------ Elliptical ---------------------
            IF ( a .gt. small ) THEN
                arg1 = DSQRT( s / ( 2.0D0*a ) )
                arg2 = DSQRT( ( s-chord ) / ( 2.0D0*a ) )
                Sinv = arg2
                Cosv = DSQRT( 1.0D0 - (magro+magr-chord)/(4.0D0*a) )
                BetE = 2.0D0*DACOS(Cosv)
                BetE = 2.0D0*DASIN(Sinv)
                IF ( DNu .gt. Pi ) THEN
                    BetE= -BetE
                  ENDIF

                Cosv= DSQRT( 1.0D0 - s/(2.0D0*a) ) 
                Sinv= arg1 

                am  = s*0.5D0 
                ae  = Pi 
                be  = 2.0D0*DASIN( DSQRT( (s-chord)/s ) ) 
                tm  = DSQRT(am*am*am/mu)*(ae - (be-DSIN(be)))
                IF ( Dtsec .gt. tm ) THEN
                    AlpE= 2.0D0*pi-2.0D0*DASIN( Sinv )
                  ELSE
                    AlpE= 2.0D0*DASIN( Sinv )
                  ENDIF
                DE  = AlpE - BetE 
                F   = 1.0D0 - (a/magro)*(1.0D0 - DCOS(DE) )
                GDot= 1.0D0 - (a/magr)* (1.0D0 - DCOS(DE) )
                G   = Dtsec - DSQRT(a*a*a/mu)*(DE - DSIN(DE))
              ELSE
                ! --------------------- Parabolic ---------------------
                arg1 = 0.0D0
                arg2 = 0.0D0 
                Error= 'a = 0 '
                Write(10,*) ' a parabolic orbit '
              ENDIF
          ENDIF

        DO i= 1 , 3
            vo(i)= ( r(i) - F*ro(i) )/G
            v(i) = ( GDot*r(i) - ro(i) )/G
          ENDDO

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE TARGET
*
*  this subroutine accomplishes the targeting problem using KEPLER/PKEPLER .and.
*    Lambert.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    RInt        - Initial Position vector of Int km
*    VInt        - Initial Velocity vector of Int km/s
*    RTgt        - Initial Position vector of Tgt km
*    VTgt        - Initial Velocity vector of Tgt km/s
*    dm          - Direction of Motion for Gauss  'L','S'
*    Kind        - Type of propagator             'K','P'
*    Dtsec        - Time of flight to the int      s
*
*  Outputs       :
*    V1t         - Initial Transfer Velocity vec  km/s
*    V2t         - Final Transfer Velocity vec    km/s
*    DV1         - Initial Change Velocity vec    km/s
*    DV2         - Final Change Velocity vec      km/s
*    Error       - Error flag from Gauss          'ok', ...
*
*  Locals        :
*    TransNormal - CROSS product of trans orbit   km
*    IntNormal   - CROSS product of int orbit     km
*    R1Tgt       - Position vector after Dt, Tgt  km
*    V1Tgt       - Velocity vector after Dt, Tgt  km/s
*    RIRT        - RInt(4) * R1Tgt(4)
*    CosDeltaNu  - Cosine of DeltaNu              rad
*    SinDeltaNu  - Sine of DeltaNu                rad
*    DeltaNu     - ANGLE between position vectors rad
*    i           - Index
*
*  Coupling      :
*    KEPLER      - Find R and V at future time
*    LAMBERTUNIV - Find velocity vectors at each ENDIF of transfer
*    LNCOM2      - Linear combination of two vectors and constants
*
*  References    :
*    Vallado       2001, 468-474, Alg 58
*
* ------------------------------------------------------------------------------

      SUBROUTINE TARGET      ( RInt,VInt,RTgt,VTgt, Dm,Kind, Dtsec,
     &                         V1t,V2t,DV1,DV2, Error  )
        IMPLICIT NONE
        REAL*8 RInt(3),VInt(3),RTgt(3),VTgt(3),Dtsec,V1t(3),V2t(3),
     &         DV1(3),DV2(3)
        CHARACTER Kind, Dm
        CHARACTER*12 Error
* -----------------------------  Locals  ------------------------------
        INTEGER i
        REAL*8 R1Tgt(3), V1Tgt(3)

        ! --------------------  Implementation   ----------------------
        ! ----------- Propogate TARGET forward by time ----------------
        IF (Kind.eq.'K') THEN
            CALL KEPLER ( RTgt,VTgt,Dtsec,  R1Tgt,V1Tgt,Error )
          ENDIF
*        IF (Kind.eq.'P') THEN
*            CALL PKEPLER( RTgt,VTgt,Dtsec,  R1Tgt,V1Tgt )
*          ENDIF

        ! ----------- Calculate transfer orbit between r's ------------
        IF ( Error .eq. 'ok' ) THEN
            CALL LAMBERTUNIV( RInt,R1Tgt,dm,'N',Dtsec,  V1t,V2t,Error )

            IF ( Error .eq. 'ok' ) THEN
                CALL LNCOM2( -1.0D0, 1.0D0,VInt, V1t,  DV1 )
                CALL LNCOM2(  1.0D0,-1.0D0,V1Tgt,V2t,  DV2 )
              ELSE
                DO i= 1 , 3
                    V1t(i)= 0.0D0
                    V2t(i)= 0.0D0
                    DV1(i)= 0.0D0
                    dV2(i)= 0.0D0
                  ENDDO
              ENDIF 
          ENDIF 
      RETURN
      END

            
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE RNGAZ
*
*  this subroutine calculates the Range and Azimuth between two specified
*    ground points on a spherical Earth.  Notice the range will ALWAYS be
*    within the range of values listed since you do not know the direction of
*    firing, Long .or. short.  The SUBROUTINE will calculate rotating Earth ranges
*    If the Tof is passed in other than 0.0D0. Range is calulated in rad .and.
*    converted to ER by s = rO, but the radius of the Earth = 1 ER, so it's
*    s = O.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    LLat        - Start Geocentric Latitude      -Pi/2 to  Pi/2 rad
*    LLon        - Start Longitude (WEST -)       0.0D0 to 2Pi rad
*    TLat        - ENDIF Geocentric Latitude        -Pi/2 to  Pi/2 rad
*    TLon        - ENDIF Longitude (WEST -)         0.0D0 to 2Pi rad
*    Tof         - Time of Flight If ICBM, .or. 0.0D0 TU
*
*  OutPuts       :
*    Range       - Range between points           ER
*    Az          - Azimuth                        0.0D0 to 2Pi rad
*
*  Locals        :
*    None.
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 843, Eq 11-3, Eq 11-4, Eq 11-5
*
* ------------------------------------------------------------------------------

      SUBROUTINE RNGAZ       ( LLat, LLon, TLat, TLon, Tof, Range, Az )
        IMPLICIT NONE
        REAL*8 LLat, LLon, TLat, TLon, Tof, Range, Az

        INCLUDE 'astmath.cmn'
        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
        Range= DACOS( DSIN(LLat)*DSIN(TLat) +
     &         DCOS(LLat)*DCOS(TLat)*DCOS(TLon-LLon + OmegaEarth*Tof) )

        ! ------ Check If the Range is 0 .or. half the earth  ---------
        IF ( DABS( DSIN(Range)*DCOS(LLat) ) .lt. Small ) THEN
            IF ( DABS( Range - Pi ) .lt. Small ) THEN
                Az= Pi
              ELSE
                Az= 0.0D0
              ENDIF
          ELSE
            Az= DACOS( ( DSIN(TLat) - DCOS(Range) * DSIN(LLat)) /
     &                 ( DSIN(Range) * DCOS(LLat)) )
          ENDIF

        ! ------ Check If the Azimuth is grt than Pi ( 180deg ) -------
        IF ( DSIN( TLon - LLon + OmegaEarth*Tof ) .lt. 0.0D0 ) THEN
            Az= TwoPi - Az
          ENDIF
      RETURN
      END

      
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE HOHMANN
*
*  this subroutine calculates the delta v's DO a Hohmann transfer DO either
*    circle to circle, or ellipse to ellipse.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    RInit       - Initial position magnitude     ER
*    RFinal      - Final position magnitude       ER
*    eInit       - Eccentricity of first orbit
*    eFinal      - Eccentricity of final orbit
*    NuInit      - True Anomaly of first orbit    0 .or. Pi rad
*    NuFinal     - True Anomaly of final orbit    0 .or. Pi rad
*
*  OutPuts       :
*    DeltaVa     - Change in velocity at point a  ER / TU
*    DeltaVb     - Change in velocity at point b  ER / TU
*    DtTU        - Time of Flight DO the trans   TU
*
*  Locals        :
*    SME1        - Mech Energy of first orbit     ER2 / TU
*    SME2        - Mech Energy of transfer orbit  ER2 / TU
*    SME3        - Mech Energy of final orbit     ER2 / TU
*    VInit       - Velocity of first orbit at a   ER / TU
*    vTransa     - Velocity of trans orbit at a   ER / TU
*    vTransb     - Velocity of trans orbit at b   ER / TU
*    VFinal      - Velocity of final orbit at b   ER / TU
*    aInit       - Semimajor axis of first orbit  ER
*    aTrans      - Semimajor axis of Trans orbit  ER
*    aFinal      - Semimajor axis of final orbit  ER
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 327, Alg 36, Ex 6-1
*
* ------------------------------------------------------------------------------

      SUBROUTINE Hohmann     ( RInit,RFinal,eInit,eFinal,NuInit,
     &                          NuFinal,Deltava,Deltavb,DtTU )
        IMPLICIT NONE
        REAL*8 RInit,RFinal,eInit,eFinal,NuInit,NuFinal,
     &          Deltava,Deltavb,DtTU
* ----------------------------  Locals  -------------------------------
        REAL*8 VInit,VTrana,VTranb,VFinal, aInit,aTran,aFinal

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        aInit  = (rInit*(1.0D0+eInit*DCOS(NuInit)))
     &           / (1.0D0 - eInit*eInit )
        aTran  = ( RInit + RFinal ) / 2.0D0
        aFinal = (rFinal*(1.0D0+eFinal*DCOS(NuFinal)))
     &           / (1.0D0 - eFinal*eFinal )
        DeltaVa= 0.0D0
        DeltaVb= 0.0D0
        DtTU   = 0.0D0

        IF ( ( eInit .lt. 1.0D0 ) .or. ( eFinal .lt. 1.0D0 ) ) THEN
            ! ----------------  Find Delta v at point a  --------------
            VInit  = DSQRT( 2.0D0/rInit - (1.0D0/aInit) )
            VTrana = DSQRT( 2.0D0/rInit - (1.0D0/aTran) )
            DeltaVa= DABS( VTrana - VInit )

            ! ----------------  Find Delta v at point b  --------------
            VFinal = DSQRT( 2.0D0/rFinal - (1.0D0/aFinal) )
            VTranb = DSQRT( 2.0D0/rFinal - (1.0D0/aTran) )
            DeltaVb= DABS( VFinal - VTranb )

            ! ---------------  Find Transfer Time of Flight  ----------
            DtTU= Pi * DSQRT( aTran*aTran*aTran ) 
          ENDIF

      RETURN
      END  ! SUBROUTINE Hohmann

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE BIELLIPTIC
*
*  this subroutine calculates the delta v's DO a Bi-elliptic transfer DO either
*    circle to circle, .or. ellipse to ellipse.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    RInit       - Initial position magnitude     ER
*    R2          - Interim orbit magnitude        ER
*    RFinal      - Final position magnitude       ER
*    eInit       - Eccentricity of first orbit
*    eFinal      - Eccentricity of final orbit
*    NuInit      - True Anomaly of first orbit    0 .or. Pi rad
*    NuFinal     - True Anomaly of final orbit    0 .or. Pi rad, Opp of NuInit
*
*  OutPuts       :
*    DeltaVa     - Change in velocity at point a  ER / TU
*    DeltaVb     - Change in velocity at point b  ER / TU
*    DtTU        - Time of Flight DO the trans   TU
*
*  Locals        :
*    SME1        - Mech Energy of first orbit     ER2 / TU
*    SME2        - Mech Energy of transfer orbit  ER2 / TU
*    SME3        - Mech Energy of final orbit     ER2 / TU
*    VInit       - Velocity of first orbit at a   ER / TU
*    vTransa     - Velocity of trans orbit at a   ER / TU
*    vTransb     - Velocity of trans orbit at b   ER / TU
*    VFinal      - Velocity of final orbit at b   ER / TU
*    aInit       - Semimajor axis of first orbit  ER
*    aTrans      - Semimajor axis of Trans orbit  ER
*    aFinal      - Semimajor axis of final orbit  ER
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 327, Alg 37, Ex 6-2
*
* ------------------------------------------------------------------------------  

      SUBROUTINE BiElliptic  ( RInit,Rb,RFinal,eInit,eFinal,
     &                         NuInit,NuFinal,
     &                         Deltava,Deltavb,DeltaVc,DtTU )
        IMPLICIT NONE
        REAL*8 RInit,Rb,RFinal,eInit,eFinal,NuInit,NuFinal,
     &         Deltava,Deltavb,DeltaVc,DtTU
* ----------------------------  Locals  -------------------------------
        REAL*8 VInit,VTran1a,VTran1b,VTran2b,VTran2c,VFinal,
     &         aInit,aTran1,aTran2,aFinal

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        aInit  = (rInit*(1.0D0+eInit*DCOS(NuInit)))
     &           / (1.0D0 - eInit*eInit )
        aTran1 = (RInit + Rb) * 0.5D0
        aTran2 = (Rb + RFinal) * 0.5D0
        aFinal = (rFinal*(1.0D0+eFinal*DCOS(NuFinal)))
     &           / (1.0D0 - eFinal*eFinal )

        DeltaVa= 0.0D0
        DeltaVb= 0.0D0
        DeltaVc= 0.0D0
        DtTU   = 0.0D0

        IF ( ( eInit .lt. 1.0D0 ) .and. ( eFinal .lt. 1.0D0 ) ) THEN
            ! ----------------  Find Delta v at point a  --------------
            VInit  = DSQRT( 2.0D0/rInit - (1.0D0/aInit) )
            VTran1a= DSQRT( 2.0D0/rInit - (1.0D0/aTran1) )
            DeltaVa= DABS( VTran1a - VInit )

            ! ----------------  Find Delta v at point b  --------------
            VTran1b= DSQRT( 2.0D0/rb - (1.0D0/aTran1) )
            VTran2b= DSQRT( 2.0D0/rb - (1.0D0/aTran2) )
            DeltaVb= DABS( VTran1b - VTran2b )

            ! ----------------  Find Delta v at point c  --------------
            VTran2c= DSQRT( 2.0D0/rFinal - (1.0D0/aTran2) )
            VFinal = DSQRT( 2.0D0/rFinal - (1.0D0/aFinal) )
            DeltaVc= DABS( VFinal - VTran2c )

            ! ---------------  Find Transfer Time of Flight  ----------
            DtTU= Pi * DSQRT( aTran1*aTran1*aTran1 ) +
     &              Pi * DSQRT( aTran2*aTran2*aTran2 )
          ENDIF
      RETURN
      END   ! SUBROUTINE BiElliptic

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE ONETANGENT
*
*  this subroutine calculates the delta v's DO a One Tangent transfer DO either
*    circle to circle, .or. ellipse to ellipse.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    RInit       - Initial position magnitude     ER
*    RFinal      - Final position magnitude       ER
*    eInit       - Eccentricity of first orbit
*    eFinal      - Eccentricity of final orbit
*    NuInit      - True Anomaly of first orbit    0 .or. Pi rad
*    Nu2         - True Anomaly of second orbit   Same quad as NuInit, rad
*    NuFinal     - True Anomaly of final orbit    0 .or. Pi rad
*
*  OutPuts       :
*    DeltaVa     - Change in velocity at point a  ER / TU
*    DeltaVb     - Change in velocity at point b  ER / TU
*    DtTU        - Time of Flight DO the transf  TU
*
*  Locals        :
*    SME1        - Mech Energy of first orbit     ER2 / TU
*    SME2        - Mech Energy of transfer orbit  ER2 / TU
*    SME3        - Mech Energy of final orbit     ER2 / TU
*    VInit       - Velocity of first orbit at a   ER / TU
*    vTransa     - Velocity of trans orbit at a   ER / TU
*    vTransb     - Velocity of trans orbit at b   ER / TU
*    VFinal      - Velocity of final orbit at b   ER / TU
*    aInit       - Semimajor axis of first orbit  ER
*    aTrans      - Semimajor axis of Trans orbit  ER
*    aFinal      - Semimajor axis of final orbit  ER
*    E           - Ecc anomaly of trans at b      rad
*    Ratio       - Ratio of initial to final
*                    orbit radii
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 335, Alg 38, Ex 6-3
*
* ------------------------------------------------------------------------------

      SUBROUTINE OneTangent  ( RInit,RFinal,eInit,eFinal,NuInit,NuTran,
     &                         Deltava,Deltavb,DtTU )
        IMPLICIT NONE
        REAL*8 RInit,RFinal,eInit,eFinal,NuInit,NuTran,
     &         Deltava,Deltavb,DtTU
* ----------------------------  Locals  -------------------------------
        REAL*8 EAInit,VInit,VTrana,VTranb,VFinal, eTran,aInit,aTran,
     &         aFinal, fpaTranb,fpaFinal, E, Sinv,Cosv,Ratio

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        DeltaVa= 0.0D0 
        DeltaVb= 0.0D0
        DtTU   = 0.0D0
        Ratio  = rInit/rFinal 
        IF ( DABS(NuInit) .lt. 0.01D0 ) THEN ! check 0 .or. 180
            eTran  = ( Ratio-1.0D0 ) / ( DCOS(NuTran)-Ratio )   ! init at perigee
            EAInit= 0.0D0 
          ELSE
            eTran  = ( Ratio-1.0D0 ) / ( DCOS(NuTran)+Ratio )  ! init at apogee
            EAInit= Pi 
          ENDIF 
        IF ( eTran .ge. 0.0D0 ) THEN
            aInit = (rInit*(1.0D0+eInit*DCOS(NuInit)))
     &              / (1.0D0 - eInit*eInit )
            aFinal= (rFinal*(1.0D0+eFinal*DCOS(NuTran)))
     &              / (1.0D0 - eFinal*eFinal )
*                       nutran is used since it = nufinal!!  
*aInit = rinit
*afinal= rfinal
            IF ( DABS( eTran-1.0D0 ) .gt. 0.000001D0 ) THEN
                IF ( DABS(NuInit) .lt. 0.01D0 ) THEN ! check 0 .or. 180
                    aTran = (rInit*(1.0D0+eTran*DCOS(NuInit)))
     &                      / (1.0D0 - eTran*eTran ) ! per
                  ELSE
*                 aTran = (rInit*(1.0D0+eTran*DCOS(NuInit))) / (1.0D0 + eTran*eTran )  apo  
                    aTran= RInit/(1.0D0 + eTran)
                  ENDIF
              ELSE
                aTran = 999999.9D0   ! Infinite DO Parabolic orbit
              ENDIF 

            ! ----------------  Find Delta v at point a  --------------
            VInit  = DSQRT( 2.0D0/rInit - (1.0D0/aInit) )
            VTrana = DSQRT( 2.0D0/rInit - (1.0D0/aTran) )
            DeltaVa= DABS( VTrana - VInit )

            ! ----------------  Find Delta v at point b  --------------
            VFinal  = DSQRT( 2.0D0/rFinal - (1.0D0/aFinal) )
            VTranb  = DSQRT( 2.0D0/rFinal - (1.0D0/aTran) )
            fpaTranb= DATAN( ( eTran*DSIN(NuTran) )
     &                / ( 1.0D0 + eTran*DCOS(NuTran) ) )
            fpaFinal= DATAN( ( eFinal*DSIN(NuTran) )
     &                / ( 1.0D0 + eFinal*DCOS(NuTran) ) )
            DeltaVb = DSQRT( VTranb*VTranb + VFinal*VFinal
     &                       - 2.0D0*VTranb*VFinal
     &                       *DCOS( fpaTranb-fpaFinal ) )

            ! ---------------  Find Transfer Time of Flight  ----------
            IF ( eTran .lt. 0.99999D0 ) THEN
                Sinv= ( DSQRT( 1.0D0-eTran*eTran )*DSIN(NuTran) )
     &                 / ( 1.0D0 + eTran*DCOS(NuTran) )
                Cosv= (eTran+DCOS(NuTran))/(1.0D0+eTran*DCOS(NuTran) )
                E   = DATAN2( Sinv,Cosv ) 
                DtTU= DSQRT( aTran*aTran*aTran ) *
     &                    ( E - eTran*DSIN(E)
     &                    - (EAInit - ETran*DSIN(EAInit)) )
              ELSE
                IF ( DABS( eTran-1.0D0 ) .lt. 0.000001D0 ) THEN
*                  Parabolic DtTU  
                  ELSE
*                  Hyperbolic DtTU  
                  ENDIF
              ENDIF

          ELSE
            Write(*,*) 'one tangent burn is not possible DO this case '
          ENDIF
      RETURN
      END   ! SUBROUTINE OneTangent

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE IONLYCHG
*
*  this subroutine calculates the delta v's DO a change in inclination only.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    DeltaI      - Change in inclination          rad
*    VInit       - Initial velocity vector        ER/TU
*    fpa         - Flight path angle              rad
*
*  OutPuts       :
*    DeltaVionly - answer
*
*  Locals        :
*    None.
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 346, Alg 39, Ex 6-4
*
* ------------------------------------------------------------------------------  

      SUBROUTINE IOnlyChg    ( Deltai,VInit,fpa, DeltaViOnly )
        IMPLICIT NONE
        REAL*8 Deltai,VInit,fpa, DeltaViOnly

        ! --------------------  Implementation   ----------------------
        DeltaViOnly = 2.0D0*VInit*DCOS(fpa)*DSIN(0.5D0*Deltai) 
      RETURN
      END   ! SUBROUTINE IOnlyChg
    

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE NODEONLYCHG
*
*  this subroutine calculates the delta v's for a change in longitude of
*    ascending node only.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    DeltaOmega  - Change in Node                 Rad
*    ecc         - Ecc of first orbit
*    VInit       - Initial velocity vector        ER/TU
*    fpa         - Flight path angle              rad
*    Incl        - Inclination                    rad
*
*
*  OutPuts       :
*    iFinal      - Final inclination              rad
*    Deltav      - Change in velocity             ER/TU
*
*  Locals        :
*    VFinal      - Final velocity vector          ER/TU
*    ArgLat      - Argument of latitude           rad
*    ArgLat1     - Final Argument of latitude     rad
*    NuInit      - Initial true anomaly           rad
*    Theta       -
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 349, Alg 40, Ex 6-5
*
* ------------------------------------------------------------------------------  

      SUBROUTINE NodeOnlyChg ( iInit,ecc,DeltaOmega,VInit,fpa,incl,
     &                         iFinal,DeltaV )
        IMPLICIT NONE
        REAL*8 iInit,ecc,DeltaOmega,VInit,fpa,incl,iFinal,DeltaV
* ----------------------------  Locals  -------------------------------
        REAL*8 ArgLat,ArgLat1,Theta

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        IF ( DABS(ecc) .gt. 0.00000001D0 ) THEN
            ! ------------------------ Elliptical ---------------------
            Theta = DATAN( DSIN(iInit)*DTAN(DeltaOmega) ) 
            iFinal= DASIN( DSIN(Theta)/DSIN(DeltaOmega) ) 
            DeltaV= 2.0D0*VInit*DCOS(fpa)*DSIN(0.5D0*Theta)

            ArgLat = Pi*0.5D0  ! set at 90 deg
            ArgLat1= DACOS( DCOS(incl)*DSIN(incl)*
     &               (1.0D0-DCOS(DeltaOmega))/ DSIN(Theta) )
          ELSE
            ! ------------------------- Circular ----------------------
            theta = DACOS( DCOS(iinit)**2
     &                +DSIN(iinit)**2*DCOS(DeltaOmega) )
            DeltaV= 2.0D0*VInit*DSIN(0.5D0*Theta)

            ArgLat = DACOS( DTAN(iInit)*(DCOS(DeltaOmega)-DCOS(Theta))
     &                      / DSIN(Theta) )
            ArgLat1= DACOS( DCOS(incl)*DSIN(incl)*
     &               (1.0D0-DCOS(DeltaOmega))/ DSIN(Theta) )

          ENDIF

      RETURN
      END   ! SUBROUTINE NodeOnlyChg

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE IandNodeChg
*
*  this subroutine calculates the delta v's for a change in inclination .and.
*    longitude of ascending node.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    VInit       - Initial velocity vector        ER/TU
*    iInit       - Initial inclination            rad
*    fpa         - Flight path angle              rad
*    DeltaOmega  - Change in Node                 Rad
*    DeltaI      - Change in inclination          Rad
*    RFinal      - Final position magnitude       ER
*
*  OutPuts       :
*    iFinal      - Final inclination              rad
*    Deltav      - Change in velocity             ER/TU
*
*  Locals        :
*    ArgLat      - Argument of latitude           rad
*    ArgLat1     - Final Argument of latitude     rad
*    Theta       -
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 350, Alg 41, Ex 6-6
*
* ------------------------------------------------------------------------------

      SUBROUTINE IandNodeChg ( iInit,DeltaOmega,Deltai,VInit,fpa,
     &                         DeltaV,iFinal )
        IMPLICIT NONE
        REAL*8 iInit,DeltaOmega,Deltai,VInit,fpa, DeltaV,iFinal
* ----------------------------  Locals  -------------------------------
        REAL*8 ArgLat, ArgLat1, Theta

        ! --------------------  Implementation   ----------------------
        iFinal= iInit - Deltai
        theta = DACOS( DCOS(iinit)*DCOS(ifinal) +
     &                    DSIN(iinit)*DSIN(ifinal)*DCOS(DeltaOmega) )
        DeltaV= 2.0D0*VInit*DCOS(fpa)*DSIN(0.5D0*Theta)

        ArgLat = DACOS( (DSIN(ifinal)*DCOS(DeltaOmega) -
     &                  DCOS(Theta)*DSIN(iinit))
     &                  / (DSIN(Theta)*DCOS(iinit)) )
        ArgLat1= DACOS( (DCOS(iInit)*DSIN(iFinal) -
     &                  DSIN(IInit)*DCOS(iFinal)*DCOS(DeltaOmega))
     &                  / DSIN(Theta) )

      RETURN
      END   ! SUBROUTINE IandNodeChg

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE MINCOMBINEDPLANECHG
*
*  this subroutine calculates the delta v's .and. the change in inclination
*    necessary DO the minimum change in velocity when traveling between two
*    non-coplanar orbits.  The notation used is from the initial orbit (1) at
*    point a, transfer is made to the transfer orbit (2), .and. to the final
*    orbit (3) at point b.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    RInit       - Initial position magnitude     ER
*    RFinal      - Final position magnitude       ER
*    eInit       - Ecc of first orbit
*    e2          - Ecc of trans orbit
*    eFinal      - Ecc of final orbit
*    NuInit      - True Anomaly of first orbit    0 .or. Pi rad
*    NuFinal     - True Anomaly of final orbit    0 .or. Pi rad
*    iInit       - Incl of the first orbit        rad
*    iFinal      - Incl of the second orbit       rad
*
*  OutPuts       :
*    Deltai1     - Amount of incl chg req at a    rad
*    DeltaVa     - Change in velocity at point a  ER / TU
*    DeltaVb     - Change in velocity at point b  ER / TU
*    DtTU        - Time of Flight DO the trans   TU
*    NumIter     - Number of iterations
*
*  Locals        :
*    SME1        - Mech Energy of first orbit     ER2 / TU
*    SME2        - Mech Energy of transfer orbit  ER2 / TU
*    SME3        - Mech Energy of final orbit     ER2 / TU
*    VInit       - Velocity of first orbit at a   ER / TU
*    vTransa     - Velocity of trans orbit at a   ER / TU
*    vTransb     - Velocity of trans orbit at b   ER / TU
*    VFinal      - Velocity of final orbit at b   ER / TU
*    aInit       - Semimajor axis of first orbit  ER
*    aTrans      - Semimajor axis of Trans orbit  ER
*    aFinal      - Semimajor axis of final orbit  ER
*    e2          - Eccentricity of second orbit
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 355, Alg 42, Table 6-3
*
* ------------------------------------------------------------------------------

      SUBROUTINE MinCombinedPlaneChg( RInit,RFinal,eInit,eFinal,NuInit,
     &                         NuFinal,iInit,iFinal,
     &                         Deltai,Deltai1,DeltaVa,DeltaVb,DtTU )
        IMPLICIT NONE
        REAL*8 RInit,RFinal,eInit,eFinal,NuInit,NuFinal,iInit,iFinal,
     &         Deltai,Deltai1,DeltaVa,DeltaVb,DtTU
* ----------------------------  Locals  -------------------------------
        INTEGER numiter
        REAL*8 deltainew,temp,TDi, SME1,SME2,SME3, VInit,
     &         V1t,V3t,VFinal, a1,a2,a3

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        ! -------------------  Initialize values   --------------------
        a1  = (RInit*(1.0D0+eInit*DCOS(NuInit))) / (1.0D0 - eInit**2 )
        a2  = 0.5D0 * (RInit+RFinal) 
        a3  = (RFinal*(1.0D0+eFinal*DCOS(NuFinal))) / (1.0D0-eFinal**2 )
        SME1= -1.0D0 / (2.0D0*a1) 
        SME2= -1.0D0 / (2.0D0*a2) 
        SME3= -1.0D0 / (2.0D0*a3) 

        ! ----------------------- Find velocities ---------------------
        VInit = DSQRT( 2.0D0*( (1.0D0/RInit) + SME1 ) ) 
        V1t   = DSQRT( 2.0D0*( (1.0D0/RInit) + SME2 ) ) 

        VFinal= DSQRT( 2.0D0*( (1.0D0/RFinal) + SME3 ) ) 
        V3t   = DSQRT( 2.0D0*( (1.0D0/RFinal) + SME2 ) ) 

        ! ---------- Find the optimum change of inclination -----------
        TDi = iFinal-iInit 

        Temp= (1.0D0/TDi) * DATAN( (RFinal/RInit**1.5D0 - DCOS(TDi))
     &        / DSIN(TDi) )
        Temp= (1.0D0/TDi) * DATAN( DSIN(TDi)
     &        / (RFinal/RInit**1.5D0 + DCOS(TDi)) )

        DeltaVa= DSQRT( V1t*V1t + VInit*VInit
     &                 - 2.0D0*V1t*VInit*DCOS(Temp*Tdi) )
        DeltaVb= DSQRT( V3t*V3t + VFinal*VFinal
     &                 - 2.0D0*V3t*VFinal*DCOS(TDi*(1.0D0-Temp)) )

        Deltai = Temp*TDi 
        Deltai1= TDi*(1.0D0-Temp) 

        ! ---------------  Find Transfer Time of Flight  --------------
        DtTU= Pi * DSQRT( A2*A2*A2 ) 

        ! ---- Iterate to find the optimum change of inclination ------
        DeltaiNew  = Deltai    ! 1st guess, 0.01D0 to 0.025D0 seems good
        Deltai1    = 100.0D0   ! IF ( going to smaller orbit, should be
        NumIter    = 0         ! 1.0D0 - 0.025D0!

        DO WHILE (DABS(DeltaiNew-Deltai1) .gt. 0.000001D0)
            Deltai1= DeltaiNew
            DeltaVa= DSQRT( V1t*V1t + VInit*VInit
     &               - 2.0D0*V1t*VInit* DCOS(Deltai1) )

            DeltaVb= DSQRT( V3t*V3t + VFinal*VFinal
     &               - 2.0D0*V3t*VFinal* DCOS(TDi-Deltai1) )

            DeltaiNew= DASIN( (DeltaVa*VFinal*V3t*DSIN(TDi-Deltai1))
     &                 / (VInit*V1t*DeltaVb) )
            NumIter = NumIter + 1
          ENDDO   ! DO WHILE DABS()

      RETURN
      END   ! SUBROUTINE MinCombinedPlaneChg

*    Vallado       2007, 355, Ex 6-7

      SUBROUTINE CombinedPlaneChg ( RInit,RFinal,eInit,e2,eFinal,NuInit,
     &                              Nu2a,Nu2b,NuFinal,Deltai,
     &                              Deltava,Deltavb,DtTU )
        IMPLICIT NONE
        REAL*8 RInit,RFinal,eInit,e2,eFinal,NuInit,Nu2a,Nu2b,
     &         NuFinal,Deltai, Deltava,Deltavb,DtTU
* ----------------------------  Locals  -------------------------------
        REAL*8 SME1,SME2,SME3, VInit,vTransa,vTransb,VFinal, a1,a2,a3,
     &         fpa1,fpa2a,fpa2b,fpa3,E,Eo,Sinv,Cosv

        ! --------------------  Implementation   ----------------------
        a1 = (RInit*(1.0D0+eInit*DCOS(NuInit))) / (1.0D0-eInit*eInit )
        IF ( DABS( e2-1.0D0 ) .gt. 0.000001D0 ) THEN
            a2  = (RInit*(1.0D0+e2*DCOS(Nu2a))) / (1.0D0 - e2*e2 )
            SME2= -1.0D0 / (2.0D0*a2)
          ELSE
            a2  = 999999.9D0   ! Undefined DO Parabolic orbit
            SME2= 0.0D0 
          ENDIF 
        a3 = (RFinal*(1.0D0+eFinal*DCOS(NuFinal)))
     &        / (1.0D0 - eFinal*eFinal )
        SME1= -1.0D0 / (2.0D0*a1) 
        SME3= -1.0D0 / (2.0D0*a3)

        ! ----------------  Find Delta v at point a  ------------------
        VInit = DSQRT( 2.0D0*( (1.0D0/RInit) + SME1 ) ) 
        vTransa= DSQRT( 2.0D0*( (1.0D0/RInit) + SME2 ) )
        fpa2a= DATAN( ( e2*DSIN(Nu2a) ) / ( 1.0D0 + e2*DCOS(Nu2a) ) ) 
        fpa1 = DATAN( ( eInit*DSIN(NuInit) )
     &                 / ( 1.0D0 + eInit*DCOS(NuInit) ) )
        DeltaVa= DSQRT( vTransa*vTransa + VInit*VInit
     &           - 2.0D0*vTransa*VInit*( DSIN(fpa2a)*DSIN(fpa1)
     &           + DCOS(fpa2a)*DCOS(fpa1)*DCOS(Deltai)) )

        ! ----------------  Find Delta v at point b  ------------------
        VFinal = DSQRT( 2.0D0*( (1.0D0/RFinal) + SME3 ) ) 
        vTransb= DSQRT( 2.0D0*( (1.0D0/RFinal) + SME2 ) ) 
        fpa2b= DATAN( ( e2*DSIN(Nu2b) ) / ( 1.0D0 + e2*DCOS(Nu2b) ) ) 
        fpa3 = DATAN( ( eFinal*DSIN(NuFinal) )
     &                / ( 1.0D0 + eFinal*DCOS(NuFinal) ) )
        DeltaVb= DSQRT( vTransb*vTransb + VFinal*VFinal
     &           - 2.0D0*vTransb*VFinal*( DSIN(fpa2b)*DSIN(fpa3)
     &           + DCOS(fpa2b)*DCOS(fpa3)*DCOS(Deltai)) )

        ! ---------------  Find Transfer Time of Flight  --------------
        Sinv= ( DSQRT( 1.0D0-e2*e2 )*DSIN(Nu2b) )
     &        / ( 1.0D0 + e2*DCOS(Nu2b) )
        Cosv= ( e2+DCOS(Nu2b) ) / ( 1.0D0 + e2*DCOS(Nu2b) )
        E= DATAN2( Sinv,Cosv ) 
        Sinv= ( DSQRT( 1.0D0-e2*e2 )*DSIN(Nu2a) )
     &        / ( 1.0D0 + e2*DCOS(Nu2a) )
        Cosv= ( e2+DCOS(Nu2a) ) / ( 1.0D0 + e2*DCOS(Nu2a) ) 
        Eo= DATAN2( Sinv,Cosv ) 
        DtTU= DSQRT( A2**3 ) * ( (E - e2*DSIN(E)) - (Eo-e2*DSIN(Eo)) )

      RETURN
      END   ! SUBROUTINE CombinedPlaneChg

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE RENDEZVOUS
*
*  this subroutine calculates parameters for a Hohmann transfer rendezvous.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Rcs1        - Radius of circular orbit int   ER
*    Rcs2        - Radius of circular orbit tgt   ER
*    eInit       - Ecc of first orbit
*    eFinal      - Ecc of final orbit
*    NuInit      - True Anomaly of first orbit    0 .or. Pi rad
*    NuFinal     - True Anomaly of final orbit    0 .or. Pi rad
*    PhaseI      - Initial phase angle (Tgt-Int)  +(ahead) .or. -(behind) rad
*    NumRevs     - Number of revs to wait
*    kTgt        -
*    kInt        -
*
*  OutPuts       :
*    PhaseF      - Final Phase Angle              rad
*    WaitTime    - Wait before next intercept opp TU
*    Deltav      - Change in velocity             ER/TU
*
*  Locals        :
*    DtTUTrans   - Time of flight of trans orbit  TU
*    ATrans      - Semimajor axis of trans orbit  ER
*    AngVelTgt   - Angular velocity of target     rad / TU
*    AngVelInt   - Angular velocity of int        rad / TU
*    LeadAng     - Lead Angle                     rad
*
*  Coupling      :
*    None
*
*  References    :
*    Vallado       2007, 364, Alg 44, Alg 45, Ex 6-8, Ex 6-9
*
* ------------------------------------------------------------------------------  

      SUBROUTINE Rendezvous  ( Rcs1,Rcs3,PhaseI,eInit,eFinal,NuInit,
     &                         NuFinal, kTgt,kInt, PhaseF,WaitTime,
     &                         DeltaV )
        IMPLICIT NONE
        REAL*8 Rcs1,Rcs3,PhaseI,eInit,eFinal,NuInit,NuFinal,
     &         PhaseF,WaitTime,DeltaV
        INTEGER kTgt,kInt
* ----------------------------  Locals  -------------------------------
        REAL*8 PeriodTrans,Rp,DtTUTrans,LeadAng,aTrans,
     &         AngVelTgt,AngVelInt, a1,a2,a3,VInit,vTransa,VFinal,
     &         vTransb,SME1,SME2,SME3,DeltaVa,DeltaVb,VInt,VTrans

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        ATrans    = (Rcs1 + Rcs3) / 2.0D0 
        DtTUTrans = Pi*DSQRT( ATrans*ATrans*ATrans ) 
        AngVelInt = 1.0D0 / ( DSQRT(Rcs1*Rcs1*Rcs1) ) 
        AngVelTgt = 1.0D0 / ( DSQRT(Rcs3*Rcs3*Rcs3) ) 
        VInt      = DSQRT( 1.0D0/Rcs1 ) 

        ! --------- Check DO satellites in the same orbits ------------
        IF ( DABS( AngVelInt - AngVelTgt ) .lt. 0.000001D0 ) THEN
            PeriodTrans= ( kTgt*TwoPi + PhaseI ) / AngVelTgt
            aTrans     = (PeriodTrans/(TwoPi*kInt))**(2.0D0/3.0D0)
            Rp         = 2.0D0*aTrans - Rcs1 
            IF ( Rp .lt. 1.0D0 ) THEN
                Write(*,*) 'Error - transfer orbit intersects Earth'
              ENDIF
            VTrans  = DSQRT( (2.0D0/Rcs1) - (1.0D0/aTrans) ) 
            DeltaV  = 2.0D0*(VTrans-VInt) 
            WaitTime= 0.0D0 
      leadang= 0.0D0
          ELSE
            LeadAng = AngVelTgt * DtTUTRans
            PhaseF  = LeadAng - Pi 
            WaitTime= ( PhaseF - PhaseI + 2.0D0*Pi*kTgt )
     &            / ( AngVelInt - AngVelTgt )

            a1  = (rcs1*(1.0D0+eInit*DCOS(NuInit)))
     &            / (1.0D0 - eInit*eInit )
            a2  = ( Rcs1 + Rcs3 ) / 2.0D0 
            a3  = (rcs3*(1.0D0+eFinal*DCOS(NuFinal)))
     &            / (1.0D0 - eFinal*eFinal )
            SME1= -1.0D0 / (2.0D0*a1) 
            SME2= -1.0D0 / (2.0D0*a2) 
            SME3= -1.0D0 / (2.0D0*a3) 
            ! ----------------  Find Delta v at point a  --------------
            VInit  = DSQRT( 2.0D0*( (1.0D0/Rcs1) + SME1 ) )
            vTransa= DSQRT( 2.0D0*( (1.0D0/Rcs1) + SME2 ) ) 
            DeltaVa= DABS( vTransa - VInit ) 

            ! ----------------  Find Delta v at point b  --------------
            VFinal = DSQRT( 2.0D0*( (1.0D0/Rcs3) + SME3 ) ) 
            vTransb= DSQRT( 2.0D0*( (1.0D0/Rcs3) + SME2 ) ) 
            DeltaVb= DABS( VFinal - vTransb )
            DeltaV = DeltaVa + DeltaVb
          ENDIF 

      RETURN
      END   ! SUBROUTINE Rendezvous

* -----
*    Vallado       2007, 370, Alg 46, Ex 6-10
*
* ------
      SUBROUTINE NonCoplanarRendz ( PhaseNew,Deltai,Delta2Node,LonTrue,
     &                              RInt,RTgt,kTgt,kInt,TTrans,TPhase,
     &                              DVPhase,DVTrans1,DVTrans2 )
        IMPLICIT NONE
        REAL*8 PhaseNew,Deltai,Delta2Node,LonTrue,RInt,RTgt,
     &         TTrans,TPhase,DVPhase, DVTrans1, DVTrans2
        INTEGER kTgt, kInt
* ----------------------------  Locals  -------------------------------
        REAL*8 AngVelInt,AngVelTgt,atrans,aphase,lead,
     &         leadNew,TNode,LonTrueNew,VInt,VTgt,VPhase,VTrans1,
     &         VTrans2

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        AngVelInt= DSQRT( 1.0D0/(RInt*RInt*RInt) ) 
        AngVelTgt= DSQRT( 1.0D0/(RTgt*RTgt*RTgt) ) 
        ATrans   = (RInt + RTgt) * 0.5D0
        TTrans = Pi*DSQRT( ATrans*ATrans*ATrans ) 

        Lead = AngVelTgt * TTRans 

        TNode= Delta2Node/AngVelInt 

        LonTrueNew= LonTrue + AngVelTgt*TNode 
*fix     PhaseNew= 13.5D0/Rad  
        LeadNew= Pi + PhaseNew

        TPhase= (LeadNew - Lead + TwoPi*kTgt) / AngVelTgt 

        aPhase = (TPhase/(TwoPi*kInt))**(2.0D0/3.0D0)

        ! ----------------  Find Deltav's  -----------------
        VInt= DSQRT(1.0D0/RInt) 
        VPhase= DSQRT(2.0D0/RInt - 1.0D0/aPhase)
        DVPhase= VPhase - VInt 

        VTrans1= DSQRT(2.0D0/RInt - 1.0D0/aTrans) 
        DVTrans1= VTrans1 - VPhase 

        VTrans2= DSQRT(2.0D0/RTgt - 1.0D0/aTrans) 
        VTgt= DSQRT(1.0D0/RTgt) 
        DVTrans2= DSQRT(VTgt*VTgt + VTrans2*VTrans2
     &            - 2.0D0*VTgt*VTrans2*DCOS(Deltai))

      RETURN
      END   ! SUBROUTINE NonCoplanarRendz

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE HILLSR
*
*  this subroutine calculates various position information for Hills equations.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R           - Initial Position vector of INT ER
*    V           - Initial Velocity Vector of INT ER / TU
*    Alt         - Altitude of TGT satellite      ER
*    DtTU        - Desired Time                   TU
*
*  Outputs       :
*    RInit       - Final Position vector of INT   ER
*    VInit       - Final Velocity Vector of INT   ER / TU
*
*  Locals        :
*    nt          - Angular velocity times time    rad
*    Omega       -
*    Sinnt       - Sine of nt
*    Cosnt       - Cosine of nt
*    Radius      - Magnitude of range vector      ER
*
*  Coupling      :
*
*
*  References    :
*    Vallado       2007, 397, Alg 47, Ex 6-14
*
* ------------------------------------------------------------------------------

      SUBROUTINE HillsR      ( R,V, Alt,DtTU, RInit,VInit )
        IMPLICIT NONE
        REAL*8 R(3), V(3),  Alt,DtTU, RInit(3), VInit(3)
* ----------------------------  Locals  -------------------------------
        REAL*8 SinNt,CosNt,Omega,nt,Radius

        ! --------------------  Implementation   ----------------------
        Radius= 1.0D0 + Alt
        Omega = DSQRT( 1.0D0 / (Radius*Radius*Radius) ) 
        nt    = Omega*DtTU 
        CosNt = DCOS( nt ) 
        SinNt = DSIN( nt ) 

        ! --------------- Determine new positions  --------------------
        RInit(1)= ( V(1)/Omega ) * SinNt -
     &           ( (2.0D0*V(2)/Omega) + 3.0D0*R(1) ) * CosNt +
     &           ( (2.0D0*V(2)/Omega) + 4.0D0*R(1) )
        RInit(2)= ( 2.0D0*V(1)/Omega ) * CosNt +
     &           ( (4.0D0*V(2)/Omega) + 6.0D0*R(1) ) * SinNt +
     &           ( R(2) - (2.0D0*V(1)/Omega) ) -
     &           ( 3.0D0*V(2) + 6.0D0*Omega*R(1) )*DtTU
        RInit(3)= R(3)*CosNt + (V(3)/Omega)*SinNt

        ! --------------- Determine new velocities  -------------------
        VInit(1)= V(1)*CosNt + (2.0D0*V(2)+3.0D0*Omega*R(1))*SinNt
        VInit(2)= -2.0D0*V(1)*SinNt + (4.0D0*V(2)
     &         +6.0D0*Omega*R(1))*CosNt - (3.0D0*V(2)+6.0D0*Omega*R(1))
        VInit(3)= -R(3)*Omega*SinNt + V(3)*CosNt

      RETURN
      END   ! SUBROUTINE HillsR

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE HILLSV
*
*  this subroutine calculates initial velocity DO Hills equations.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R           - Initial Position vector of INT ER
*    Alt         - Altitude of TGT satellite      ER
*    DtTU        - Desired Time                   TU
*
*  Outputs       :
*    V           - Initial Velocity Vector of INT ER / TU
*
*  Locals        :
*    Numer       -
*    Denom       -
*    nt          - Angular velocity times time    rad
*    Omega       -
*    Sinnt       - Sine of nt
*    Cosnt       - Cosine of nt
*    Radius      - Magnitude of range vector      ER
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 410, Eq 6-60, Ex 6-15
*
* ------------------------------------------------------------------------------  

      SUBROUTINE HillsV      ( R, Alt,DtTU, V )
        IMPLICIT NONE
        REAL*8 R(3), Alt,DtTU, V(3)
* ----------------------------  Locals  -------------------------------
        REAL*8 Numer,Denom,SinNt,CosNt,Omega,nt,Radius, DCot
        EXTERNAL DCOT

        ! --------------------  Implementation   ----------------------
        Radius= 1.0D0 + Alt
        Omega = DSQRT( 1.0D0 / (Radius*Radius*Radius) ) 
        nt    = Omega*DtTU 
        CosNt = DCOS( nt ) 
        SinNt = DSIN( nt )

        ! --------------- Determine initial Velocity ------------------
        Numer= ( (6.0D0*r(1)*(nt-SinNt)-r(2))*Omega*SinNt
     &          - 2.0D0*Omega*r(1)*(4.0D0-3.0D0*CosNt)*(1.0D0-CosNt) )
        Denom= (4.0D0*SinNt-3.0D0*nt)*SinNt + 4.0D0*( 1.0D0-CosNt )
     &          *( 1.0D0-CosNt )

        IF ( DABS( Denom ) .gt. 0.000001D0 ) THEN
            V(2)= Numer / Denom
          ELSE
            V(2)= 0.0D0
          ENDIF
        IF ( DABS( SinNt ) .gt. 0.000001D0 ) THEN
            V(1)= -( Omega*r(1)*(4.0D0-3.0D0*CosNt)
     &            +2.0D0*(1.0D0-CosNt)*v(2) ) / ( SinNt )
          ELSE
            V(1)= 0.0D0
          ENDIF
        V(3)= -R(3)*Omega*DCot(nt)

      RETURN
      END   ! SUBROUTINE HillsV

      SUBROUTINE IJK_RSW     ( Rijk,Vijk,R,S,W, Direction,Rrsw,Vrsw )
        IMPLICIT NONE
        REAL*8 Rijk(3),Vijk(3),R(3),S(3),W(3),Rrsw(3),Vrsw(3)
        Character*4 Direction

        REAL*8 ROTMat(3,3)
        INTEGER i

        IF ( Direction .eq. 'FROM' ) THEN
            ! ---------------- Form Rotation Matrix -------------------
            DO i= 1, 3
                ROTMat(i,1)= R(i)
                ROTMat(i,2)= S(i)
                ROTMat(i,3)= W(i)
              ENDDO

            ! ----------------- Do multiplication ---------------------
            DO i= 1, 4
                rijk(i)= 0.0D0
                vijk(i)= 0.0D0
              ENDDO
            CALL MatVecMult( ROTMat,Rrsw,3,3,3,1, Rijk )
            CALL MatVecMult( ROTMat,Vrsw,3,3,3,1, Vijk )
          ELSE
            ! --------------- Form Rotation Matrix --------------------
            DO i= 1, 3
                ROTMat(1,i)= R(i)
                ROTMat(2,i)= S(i)
                ROTMat(3,i)= W(i)
              ENDDO

            ! ----------------- Do multiplication ---------------------
            DO i= 1, 4
                Rrsw(i)= 0.0D0
                Vrsw(i)= 0.0D0
              ENDDO
            CALL MatVecMult( ROTMat,Rijk,3,3, 3,1, Rrsw )
            CALL MatVecMult( ROTMat,Vijk,3,3, 3,1, Vrsw )
          ENDIF 

      RETURN
      END   ! SUBROUTINE IJKRSW

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE COW2Hill
*
*  this subroutine finds the equivalent relative motion vector given a geocentric
*    vector.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Rtgt        - Position vector of tgt         ER
*    Vtgt        - Velocity Vector of tgt         ER / TU
*    Rint        - Position vector of int         ER
*    Vint        - Velocity Vector of int         ER / TU
*
*  Outputs       :
*    RHill       - Position vector of int rel to
*                  target                         ER
*    VHill       - Velocity Vector of int rel to
*                  target                         ER / TU
*
*  Locals        :
*    None.
*
*  Coupling      :
*    CROSS       - Cross product of two vectors
*    NORM        - Unit vector
*    IJK_RSW     -
*    LNCOM2      - Linear combination of two scalars .and. two vectors
*    ROT3        - Rotation about the 3rd axis
*
*  References    :
*    Vallado       2007, 413
*
* ------------------------------------------------------------------------------

      SUBROUTINE Cow2Hill    ( rtgt,vtgt,rint,vint, RHill,VHill )
        IMPLICIT NONE
        REAL*8 rtgt(3), vtgt(3), rint(3), vint(3), RHill(3), VHill,
     &         MAG

        EXTERNAL MAG

* ----------------------------  Locals  -------------------------------
        REAL*8 r(3), s(3), w(3), rv(3),  rt(3), magrt, magri,
     &         vt(3), ri(3), vi(3), angly, anglz

        ! --------------------  Implementation   ----------------------
        ! ---------- Form RSW unit vectors DO transformation ----------
        CALL NORM ( rtgt,  R )
        CALL CROSS( rtgt,vtgt, RV )
        CALL NORM ( RV,  W )
        CALL CROSS( W,R, S )

        CALL IJK_RSW( RTgt,VTgt,R,S,W,'TOO', Rt,Vt )   ! IJK to RSW
        magrt = MAG(rt)
        CALL IJK_RSW( RInt,VInt,R,S,W,'TOO', Ri,Vi )   ! IJK to RSW

       ! --- Determine z offset to correct vector ----
        IF ( DABS(ri(3)) .gt. 0.0000001D0 ) THEN
            anglz= DATAN(ri(3)/magrt)  ! coord sys based at tgt
            CALL ROT2( Ri, -anglz, Ri )
            CALL ROT2( Vi, -anglz, Vi )  ! should be ROT2(a), but opp
          ELSE
            anglz= 0.0D0 
          ENDIF

        ! -------------- Determine y offset to correct vector ---------
        IF ( DABS(ri(2)) .gt. 0.0000001D0 ) THEN
            angly= DATAN(ri(2)/magrt)  ! should be ROT3(-a), but opp, but sign DO later
            CALL ROT3( Ri, angly, Ri )
            CALL ROT3( Vi, angly, Vi )
          ELSE
            angly= 0.0D0
          ENDIF

        ! ---------------------------- Do all 3 here ------------------
*     LNCOM2( 1.0D0,-1.0D0,ri,rt, rHill )
        CALL LNCOM2( 1.0D0,-1.0D0,vi,vt, VHill )

        ! ------------------- Now add in corrections ------------------
        magri = MAG(ri)
        RHill(1)= ri(1) - magrt  !4if not do rotri or 1
        RHill(2)= angly*magri
        RHill(3)= anglz*magri

      RETURN
      END   ! SUBROUTINE Cow2Hill

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE HILL2COW
*
*  this subroutine finds the equivalent geocentric vector given the target
*    and relative motion vectors.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Rtgt        - Position vector of tgt         ER
*    Vtgt        - Velocity Vector of tgt         ER / TU
*    RHill       - Position vector of int rel to
*                  target                         ER
*    VHill       - Velocity Vector of int rel to
*                  target                         ER / TU
*
*  Outputs       :
*    Rint        - Position vector of int         ER
*    Vint        - Velocity Vector of int         ER / TU
*
*  Locals        :
*    None.
*
*  Coupling      :
*    IJKRSW      - Form translation matrix given position and velocity vectors
*    ROT2        - Rotation about the 2nd axis
*    ROT3        - Rotation about the 3rd axis
*    MAG         - Magnitude of a vector
*    NORM        - Unit vector
*    CROSS       - Cross product of two vectors
*
*  References    :
*    Vallado       2007, 414
*
* ------------------------------------------------------------------------------

      SUBROUTINE Hill2Cow    ( rTgtijk,vTgtijk,RHill,VHill,
     &                         rIntijk,vIntijk )
        IMPLICIT NONE
        REAL*8 rTgtijk(3), vTgtijk(3), RHill(3), VHill(3),
     &         rIntijk(3), vIntijk(3), MAG
        EXTERNAL MAG

* ----------------------------  Locals  -------------------------------
        REAL*8 angly, anglz, rtem(3), vtem(3),  RTgtrsw(3), VTgtrsw(3),
     &         rv(3), R(3), S(3), W(3), magrtem

        ! --------------------  Implementation   ----------------------
        ! --- Form RSW unit vectors DO transformation ----
        CALL NORM ( rTgtijk,  R )
        CALL CROSS( rTgtijk,vTgtijk, RV )
        CALL NORM ( RV,  W )
        CALL CROSS( W,R, S )

        ! IJK to RSW
        CALL IJK_RSW( RTgtijk,VTgtijk,R,S,W,'TOO', RTgtrsw,VTgtrsw )

        RTem(1)= RTgtrsw(1)+RHill(1)   ! in RSW
        RTem(2)= RTgtrsw(2)
        RTem(3)= RTgtrsw(3)
        vTem(1)= vTgtrsw(1)+VHill(1)
        vTem(2)= vTgtrsw(2)+VHill(2)
        vTem(3)= vTgtrsw(3)+VHill(3)
        magrtem = MAG( rTem)

        ! --- Now perform rotation to fix y ----
        IF ( DABS(RHill(2)) .gt. 0.0000001D0 ) THEN
            angly= DATAN(RHill(2)/magrTem)  ! rtgt, but IF ( x non-zero, needs extra
            CALL ROT3( rTem,-angly,  rTem )  ! should be ROT3(a) but opp
            CALL ROT3( vTem,-angly,  vTem )
          ELSE
            angly= 0.0D0
          ENDIF

        ! --- Now perform rotation to fix z ----
        IF ( DABS(RHill(3)) .gt. 0.0000001D0 ) THEN
            anglz= DATAN(RHill(3)/magrTem)  ! should be ROT2(-a), but opp
            CALL ROT2( rTem,anglz,  rTem )
            CALL ROT2( vTem,anglz,  vTem )
          ELSE
            anglz= 0.0D0
          ENDIF

        ! --- Now RSW to IJK via MatMult!! ----
        CALL IJK_RSW( RIntijk,VIntijk,R,S,W, 'FROM', RTem,VTem )
      RETURN
      END   ! SUBROUTINE Hill2Cow
      
      
      
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE PATH
*
*  this subroutine determines the ENDIF position for a given range and azimuth
*    from a given point.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    LLat        - Start Geocentric Latitude      -Pi/2 to  Pi/2 rad
*    LLon        - Start Longitude (WEST -)       0.0D0 to 2Pi rad
*    Range       - Range between points           km
*    Az          - Azimuth                        0.0D0 to 2Pi rad
*
*  OutPuts       :
*    TLat        - ENDIF Geocentric Latitude        -Pi/2 to  Pi/2 rad
*    TLon        - ENDIF Longitude (WEST -)         0.0D0 to 2Pi rad
*
*  Locals        :
*    SinDeltaN   - Sine of Delta N                rad
*    CosDeltaN   - Cosine of Delta N              rad
*    DeltaN      - ANGLE between the two points   rad
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 843, Eq 11-6, Eq 11-7
*
* ------------------------------------------------------------------------------

      SUBROUTINE PATH        ( LLat, LLon, Range, Az, TLat, TLon )
        IMPLICIT NONE
        REAL*8 LLat, LLon, Range, Az, TLat, TLon
* -----------------------------  Locals  ------------------------------
        REAL*8 SinDN, CosDN, DeltaN

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        Az= DMOD( Az, TwoPi )
        IF ( LLon .lt. 0.0D0 ) THEN
            LLon= TwoPi + LLon 
          ENDIF   
        IF ( Range .gt. TwoPi ) THEN
            Range= DMOD( Range, TwoPi )
          ENDIF

        ! ----------------- Find Geocentric Latitude  -----------------
        TLat = DASIN( DSIN(LLat)*DCOS(Range) +
     &         DCOS(LLat)*DSIN(Range)*DCOS(Az) )

        ! ---- Find Delta N, the ANGLE between the points -------------
        IF ( (DABS(DCOS(TLat)) .gt. Small) .and.
     &        (DABS(DCOS(LLat)) .gt. Small) ) THEN
            SinDN = DSIN(Az)*DSIN(Range) / DCOS(TLat)
            CosDN = ( DCOS(Range)-DSIN(TLat)*DSIN(LLat) ) /
     &                 ( DCOS(TLat)*DCOS(LLat) )
            DeltaN= DATAN2(SinDN, CosDN)
          ELSE
            ! ------ Case where launch is within 3nm of a Pole --------
            IF ( DABS(DCOS(LLat)) .le. Small ) THEN
                IF ( (Range .gt. Pi) .and. (Range .lt. TwoPi) ) THEN
                    DeltaN= Az + Pi
                  ELSE
                    DeltaN= Az 
                  ENDIF
              ENDIF
            ! ----- Case where ENDIF point is within 3nm of a pole ----
            IF ( DABS( DCOS(TLat) ) .le. Small ) THEN
                DeltaN= 0.0D0 
              ENDIF
          ENDIF 

        TLon= LLon + DeltaN
        IF ( DABS(TLon) .gt. TwoPi ) THEN
            TLon= DMOD( TLon, TwoPi )
          ENDIF
        IF ( TLon .lt. 0.0D0 ) THEN
            TLon= TwoPi + TLon 
          ENDIF   
      RETURN
      END

* ----------------------------------------------------------------------------
*
*                           SUBROUTINE PLANETRV
*
*  this subroutine calculate the planetary ephemerides using the Epoch J2000.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    PlanetNum   - Number of planet               1..9D0
*    JD          - Julian Date                    days from 4713 BC
*
*  OutPuts       :
*    R           - XYZ position vector            AU
*    V           - XYZ velocity vector            AU / day
*
*  Locals        :
*    ArgLat      -
*    TrueLon     -
*    LonPer      -
*    TUUT1       - Julian cenuries from Jan 1, 2000
*    TU2         - Tu squared
*    TU3         - TU Cubed
*    N           -
*    obliquity   - angle between ecliptic .and.
*                    Earth equator                rad
*    a           - Semi .or. axis
*    ecc         - eccentricity
*    p           - semi-parameter
*    incl        - inclination
*    omega       - ecliptic long of ascending node
*    argp        - ecliptic arg of perihelion
*    nu          - ecliptic true anomaly
*    m           - ecliptic mean anomaly
*    LLong       - True longitude
*    LongP       - longitude of perihelion
*    e0          -
*
*  Coupling      :
*    LnCom1      -
*    NewtonM     -
*    RandV       -
*
*  References    :
*    Vallado       2007, 995
*
* ----------------------------------------------------------------------------

      SUBROUTINE PlanetRV ( PlanetNum,WhichEpoch,WhichCoord,JD, R,V )
        IMPLICIT NONE
        REAL*8 JD, R(4), V(4)
        INTEGER PlanetNum
        CHARACTER*3 WhichEpoch, WhichCoord

        INCLUDE 'astmath.cmn'

* ----------------------------  Locals  -------------------------------
        REAL*8 llong,longp,m,e0, TUDaySun,ArgLat,
     &         TrueLon,LonPer,Tut1,tut12,Tut13,eps, a,ecc,p,
     &         incl,omega,argp,nu

       ! ---------------- Approximate TTDB with TUT1 -----------------
        Tut1 = ( JD - 2451545.0D0 ) / 36525.0D0
        Tut12= Tut1*Tut1
        Tut13= Tut12*Tut1

        IF ( WhichEpoch .eq. 'J20' ) THEN
           ! --------- Mean equinox of J2000 in radians (XYZ) --------
            IF (PlanetNum.eq.1) THEN
                   a    =  0.387098310D0
                   ecc  =  0.20563175D0 +    0.000020406D0*Tut1
     &                  - 0.0000000284D0*Tut12 - 0.00000000017D0*Tut13
                   incl =  0.12226007D0 -    0.000103875D0*TUT1
     &                  + 0.0000000141D0*TUT12 + 0.00000000072D0*TUT13
                   longp=  1.35186430D0 +    0.002772705D0*TUT1
     &                  - 0.0000002344D0*TUT12 + 0.00000000068D0*TUT13
                   Omega=  0.84353321D0 -    0.002189043D0*TUT1
     &                  - 0.0000015416D0*TUT12 - 0.00000000342D0*TUT13
                   LLong=  4.40260885D0 + 2608.790314157D0*TUT1
     &                  - 0.0000000934D0*TUT12 + 0.00000000003D0*TUT13
              ELSEIF (PlanetNum.eq.2) THEN
                   a    =  0.723329820D0
                   ecc  =  0.00677188D0 -    0.000047766D0*Tut1
     &                  + 0.0000000975D0*Tut12 + 0.00000000044D0*TUT13
                   incl =  0.05924803D0 -    0.000014954D0*TUT1
     &                  - 0.0000005662D0*TUT12 + 0.00000000017D0*TUT13
                   longp=  2.29621986D0 +    0.000084903D0*TUT1
     &                  - 0.0000241260D0*TUT12 - 0.00000009306D0*TUT13
                   Omega=  1.33831707D0 -    0.004852155D0*TUT1
     &                  - 0.0000024881D0*TUT12 - 0.00000000346D0*TUT13
                   LLong=  3.17614670D0 + 1021.328554621D0*TUT1
     &                  + 0.0000000288D0*TUT12 - 0.00000000003D0*TUT13
              ELSEIF (PlanetNum.eq.3) THEN
                   a    =  1.000001018D0
                   ecc  =  0.01670862D0 -    0.000042037D0*Tut1
     &                  - 0.0000001236D0*TUT12 + 0.00000000004D0*TUT13
                   incl =  0.00000000D0 +    0.000227846D0*TUT1
     &                  - 0.0000001625D0*TUT12 - 0.00000000059D0*TUT13
                   LongP=  1.79659565D0 +    0.005629659D0*TUT1
     &                  + 0.0000026225D0*TUT12 + 0.00000000834D0*TUT13
                   Omega=  0.0D0 
                   LLong=  1.75347032D0 +  628.307584919D0*TUT1
     &                  - 0.0000000991D0*TUT12 + 0.00000000000D0*TUT13
              ELSEIF (PlanetNum.eq.4) THEN
                   a    =  1.523679342D0
                   ecc  =  0.09340062D0 +    0.000090483D0*Tut1
     &                  - 0.0000000806D0*Tut12 - 0.00000000035D0*TUT13
                   incl =  0.03228381D0 -    0.000142208D0*TUT1
     &                  - 0.0000003936D0*TUT12 - 0.00000000047D0*TUT13
                   longp=  5.86535757D0 +    0.007747339D0*TUT1
     &                  - 0.0000030231D0*TUT12 + 0.00000000524D0*TUT13
                   Omega=  0.86495189D0 -    0.005148453D0*TUT1
     &                  - 0.0000111689D0*TUT12 - 0.00000003740D0*TUT13
                   LLong=  6.20348092D0 +  334.061243148D0*TUT1
     &                  + 0.0000000456D0*TUT12 - 0.00000000005D0*TUT13
              ELSEIF (PlanetNum.eq.5) THEN
                   a    =  5.202603191D0 +   0.0000001913D0*TUT1
                   ecc  =  0.04849485D0 +    0.000163244D0*Tut1
     &                  - 0.0000004719D0*Tut12 - 0.00000000197D0*TUT13
                   incl =  0.02274635D0 -    0.000034683D0*TUT1
     &                  + 0.0000005791D0*TUT12 + 0.00000000161D0*TUT13
                   longp=  0.25012853D0 +    0.003762101D0*TUT1
     &                  + 0.0000126104D0*TUT12 - 0.00000008011D0*TUT13
                   Omega=  1.75343528D0 +    0.003083697D0*TUT1
     &                  + 0.0000157755D0*TUT12 - 0.00000012273D0*TUT13
                   LLong=  0.59954650D0 +   52.969096509D0*TUT1
     &                  - 0.0000014837D0*TUT12 + 0.00000000007D0*TUT13
              ELSEIF (PlanetNum.eq.6) THEN
                   a    =  9.554909596D0 -   0.0000021389D0*TUT1
                   ecc  =  0.05550862D0 -    0.000346818D0*Tut1
     &                  - 0.0000006456D0*Tut12 + 0.00000000338D0*TUT13
                   incl =  0.04343912D0 +    0.000044532D0*TUT1
     &                  - 0.0000008557D0*TUT12 + 0.00000000031D0*TUT13
                   longp=  1.62414732D0 +    0.009888156D0*TUT1
     &                  + 0.0000092169D0*TUT12 + 0.00000008521D0*TUT13
                   Omega=  1.98383764D0 -    0.004479648D0*TUT1
     &                  - 0.0000032018D0*TUT12 + 0.00000000623D0*TUT13
                   LLong=  0.87401675D0 +   21.329909545D0*TUT1
     &                  + 0.0000036659D0*TUT12 - 0.00000000033D0*TUT13
              ELSEIF (PlanetNum.eq.7) THEN
                   a    = 19.218446062D0 -   0.0000000372D0*TUT1
     &                  + 0.00000000098D0*TUT12
                   ecc  =  0.04629590D0 -    0.000027337D0*Tut1
     &                  + 0.0000000790D0*Tut12 + 0.00000000025D0*TUT13
                   incl =  0.01349482D0 -    0.000029442D0*TUT1
     &                  + 0.0000000609D0*TUT12 + 0.00000000028D0*TUT13
                   longp=  3.01950965D0 +    0.001558939D0*TUT1
     &                  - 0.0000016528D0*TUT12 + 0.00000000721D0*TUT13
                   Omega=  1.29164744D0 +    0.001294094D0*TUT1
     &                  + 0.0000070756D0*TUT12 + 0.00000000182D0*TUT13
                   LLong=  5.48129387D0 +    7.478159856D0*TUT1
     &                  - 0.0000000848D0*TUT12 + 0.00000000010D0*TUT13
              ELSEIF (PlanetNum.eq.8) THEN
                   a    = 30.110386869D0 -   0.0000001663D0*TUT1
     &                  + 0.00000000069D0*TUT12
                   ecc  =  0.00898809D0 +    0.000006408D0*Tut1
     &                  - 0.0000000008D0*TUT12
                   incl =  0.03089149D0 +    0.000003939D0*TUT1
     &                  + 0.0000000040D0*TUT12 - 0.00000000000D0*TUT13
                   longp=  0.83991686D0 +    0.000508915D0*TUT1
     &                  + 0.0000012306D0*TUT12 - 0.00000000000D0*TUT13
                   Omega=  2.30006570D0 -    0.000107601D0*TUT1
     &                  - 0.0000000382D0*TUT12 - 0.00000000136D0*TUT13
                   LLong=  5.31188628D0 +    3.813303564D0*TUT1
     &                  + 0.0000000103D0*TUT12 - 0.00000000003D0*TUT13
              ELSEIF (PlanetNum.eq.9) THEN
                   a    = 39.53758D0
                   ecc  =  0.250877D0 
                   incl =  0.2990156D0 
                   LongP=  3.920268D0 
                   Omega=  1.926957D0
                   LLong=  3.8203049D0
              ENDIF
          ENDIF

        IF ( WhichEpoch .eq. 'ODA' ) THEN
           ! ----------- Mean equinox of date in radians (XYZ) -----------
            IF (PlanetNum.eq.1) THEN
                   a    =  0.387098310D0
                   ecc  =  0.20563175D0 +    0.000020406D0*Tut1
     &                  - 0.0000000284D0*Tut12 - 0.00000000017D0*Tut13
                   incl =  0.12226007D0 +    0.000031791D0*TUT1
     &                  - 0.0000003157D0*TUT12 + 0.00000000093D0*TUT13
                   longp=  1.35186430D0 +    0.027165657D0*TUT1
     &                  + 0.0000051643D0*TUT12 + 0.00000000098D0*TUT13
                   Omega=  0.84353321D0 +    0.020702904D0*TUT1
     &                  + 0.0000030695D0*TUT12 + 0.00000000368D0*TUT13
                   LLong=  4.40260885D0 + 2608.814707111D0*TUT1
     &                  + 0.0000053053D0*TUT12 + 0.00000000031D0*TUT13
              ELSEIF (PlanetNum.eq.2) THEN
                   a    =  0.723329820D0
                   ecc  =  0.00677188D0 -    0.000047766D0*Tut1
     &                  + 0.0000000975D0*TUT12 + 0.00000000044D0*TUT13
                   incl =  0.05924803D0 +    0.000017518D0*TUT1
     &                  - 0.0000000154D0*TUT12 - 0.00000000012D0*TUT13
                   longp=  2.29621986D0 +    0.024473335D0*TUT1
     &                  - 0.0000187338D0*TUT12 - 0.00000009276D0*TUT13
                   Omega=  1.33831707D0 +    0.015727494D0*TUT1
     &                  + 0.0000070974D0*TUT12 - 0.00000000140D0*TUT13
                   LLong=  3.17614670D0 + 1021.352943053D0*TUT1
     &                  + 0.0000054210D0*TUT12 + 0.00000000026D0*TUT13
              ELSEIF (PlanetNum.eq.3) THEN
                   a    =  1.000001018D0
                   ecc  =  0.01670862D0 -    0.000042037D0*Tut1
     &                  - 0.0000001236D0*TUT12 + 0.00000000004D0*TUT13
                   incl =  0.000D0 
                   Omega=  0.0D0 
                   Longp=  1.79659565D0 +    0.030011406D0*TUT1
     &                  + 0.0000080219D0*TUT12 + 0.00000000871D0*TUT13
                   LLong=  1.75347032D0 +  628.331966666D0*TUT1
     &                  + 0.0000053002D0*TUT12 + 0.00000000037D0*TUT13
              ELSEIF (PlanetNum.eq.4) THEN
                   a    =  1.523679342D0
                   ecc  =  0.09340062D0 +    0.000090483D0*Tut1
     &                  - 0.0000000806D0*Tut12 - 0.00000000035D0*TUT13
                   incl =  0.03228381D0 -    0.000010489D0*TUT1
     &                  + 0.0000002227D0*TUT12 - 0.00000000010D0*TUT13
                   longp=  5.86535757D0 +    0.032132089D0*TUT1
     &                  + 0.0000023588D0*TUT12 + 0.00000000555D0*TUT13
                   Omega=  0.86495189D0 +    0.013475553D0*TUT1
     &                  + 0.0000002801D0*TUT12 + 0.00000004058D0*TUT13
                   LLong=  6.20348092D0 +  334.085627899D0*TUT1
     &                  + 0.0000054275D0*TUT12 + 0.00000000026D0*TUT13
              ELSEIF (PlanetNum.eq.5) THEN
                   a    =  5.202603191D0 +   0.0000001913D0*TUT1
                   ecc  =  0.04849485D0 +    0.000163244D0*Tut1
     &                  - 0.0000004719D0*Tut12 - 0.00000000197D0*TUT13
                   incl =  0.02274635D0 -    0.000095934D0*TUT1
     &                  + 0.0000000812D0*TUT12 - 0.00000000007D0*TUT13
                   longp=  0.25012853D0 +    0.028146345D0*TUT1
     &                  + 0.0000179991D0*TUT12 - 0.00000007974D0*TUT13
                   Omega=  1.75343528D0 +    0.017819026D0*TUT1
     &                  + 0.0000070017D0*TUT12 + 0.00000000993D0*TUT13
                   LLong=  0.59954650D0 +   52.993480754D0*TUT1
     &                  + 0.0000039050D0*TUT12 + 0.00000000044D0*TUT13
              ELSEIF (PlanetNum.eq.6) THEN
                   a    =  9.554909596D0 -   0.0000021389D0*Tut1
                   ecc  =  0.05550862D0 -    0.000346818D0*Tut1
     &                  - 0.0000006456D0*Tut12 + 0.00000000338D0*TUT13
                   incl =  0.04343912D0 -    0.000065211D0*TUT1
     &                  - 0.0000002646D0*TUT12 + 0.00000000155D0*TUT13
                   longp=  1.62414732D0 +    0.034274242D0*TUT1
     &                  + 0.0000146184D0*TUT12 + 0.00000008550D0*TUT13
                   Omega=  1.98383764D0 +    0.015308246D0*TUT1
     &                  - 0.0000021061D0*TUT12 - 0.00000004154D0*TUT13
                   LLong=  0.87401675D0 +   21.354295630D0*TUT1
     &                  + 0.0000090673D0*TUT12 - 0.00000000005D0*TUT13
              ELSEIF (PlanetNum.eq.7) THEN
                   a    = 19.218446062D0 -   0.0000000372D0*Tut1
     &                   + 0.00000000098D0*TUT12
                   ecc  =  0.04629590D0 -    0.000027337D0*Tut1
     &                  + 0.0000000790D0*Tut12 + 0.00000000025D0*TUT13
                   incl =  0.01349482D0 +    0.000013516D0*TUT1
     &                  + 0.0000006543D0*TUT12 - 0.00000000161D0*TUT13
                   longp=  3.01950965D0 +    0.025942197D0*TUT1
     &                  + 0.0000037437D0*TUT12 + 0.00000000756D0*TUT13
                   Omega=  1.29164744D0 +    0.009095361D0*TUT1
     &                  + 0.0000233843D0*TUT12 + 0.00000032317D0*TUT13
                   LLong=  5.48129387D0 +    7.502543115D0*TUT1
     &                  + 0.0000053117D0*TUT12 + 0.00000000045D0*TUT13
              ELSEIF (PlanetNum.eq.8) THEN
                   a    = 30.110386869D0 -   0.0000001663D0*Tut1
     &                  + 0.00000000069D0*TUT12
                   ecc  =  0.00898809D0 +    0.000006408D0*Tut1
     &                  - 0.0000000008D0*TUT12
                   incl =  0.03089149D0 -    0.000162459D0*TUT1
     &                  - 0.0000001236D0*TUT12 + 0.00000000049D0*TUT13
                   longp=  0.83991686D0 +    0.024893067D0*TUT1
     &                  + 0.0000066179D0*TUT12 - 0.00000000005D0*TUT13
                   Omega=  2.30006570D0 +    0.019237118D0*TUT1
     &                  + 0.0000045389D0*TUT12 - 0.00000001110D0*TUT13
                   LLong=  5.31188628D0 +    3.837687716D0*TUT1
     &                  + 0.0000053976D0*TUT12 + 0.00000000031D0*TUT13
              ELSEIF (PlanetNum.eq.9) THEN
                   a    = 39.53758D0
                   ecc  =  0.250877D0 
                   incl =  0.2990156D0 
                   LongP=  3.920268D0
                   Omega=  1.926957D0 
                   LLong=  3.8203049D0 
                 ENDIF 
          ENDIF 

        LLong= DMOD( LLong,TwoPI ) 
        LongP= DMOD( LongP,TwoPI ) 
        Omega= DMOD( Omega,TwoPI ) 
        Argp = LongP - Omega 
        M    = LLong - LongP 

        CALL NewTonM( ecc,M,  E0,Nu )

       ! ------------ Find Heliocentric ecliptic r .and. v -------------  
        p      = a*(1.0D0-ecc*ecc)
        ArgLat = Argp+Nu
        TrueLon= Omega+Argp+Nu
        LonPer = Omega+Argp+Pi
        CALL COE2RV(P,ecc,Incl,Omega,Argp,Nu,ArgLat,TrueLon,LonPer,R,V)

       ! --- Correct the velocity because we used TTdb - days!! ------
        TUDaySun= 1.0D0/58.1324409D0   ! 1.0D0 / days per sun TU
        CALL LNCOM1( TUDaySun,V,  V )

        IF ( WhichCoord .eq. 'GEO' ) THEN
           ! ----------- Find obliquity of the ecliptic angle --------
            Eps = 23.439291D0 - 0.0130042D0*TUT1 - 0.000000164D0*TUT12
     &              + 0.000000504D0*TUT13
            Eps = DMOD( Eps,360.0D0 ) 
            Eps = Eps * Deg2Rad

           ! ------------- Rotate to Geocentric coordinates ----------  
            CALL ROT1( R ,-eps, R )
            CALL ROT1( V ,-eps, V )
          ENDIF

       RETURN
       END

* ----------------------------------------------------------------------------
*
*                           SUBROUTINE INTERPLANETARY
*
*  this subroutine calculates the delta v's for an interplanetary mission.
*    The transfer assumes circular orbits for each of the planets.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R1          - Radius of planet 1 from sun    km
*    R2          - Radius of planet 2 from sun    km
*    Rbo         - Radius at burnout at planet 1  km
*    Rimpact     - Radius at impact on planet 2   km
*    Mu1         - Grav parameter of planet 1     km3/s2
*    Mut         - Grav parameter of planet Sun   km3/s2
*    Mu2         - Grav parameter of planet 2     km3/s2
*
*  OutPuts       :
*    DeltaV1     - Hyperb Exc vel at planet 1 SOI km/s
*    DeltaV2     - Hyperb Exc vel at planet 2 SOI km/s
*    Vbo         - Burnout vel at planet 1        km/s
*    Vretro      - Retro vel at surface planet 2  km/s
*
*  Locals        :
*    SME1        - Spec Mech Energy of 1st orbit  Km2/s
*    SMEt        - Spec Mech Energy of trans orbitKm2/s
*    SME2        - Spec Mech Energy of 2nd orbit  Km2/s
*    Vcs1        - Vel of 1st orbit at dv 1 point Km/s
*    Vcs2        - Vel of 2nd orbit at dv 2 point Km/s
*    Vt1         - Vel of Trans orbit at dv 1 pnt Km/s
*    Vt2         - Vel of Trans orbit at dv 2 pnt Km/s
*    A           - Semimajor Axis of Trans orbit  Km
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 47-48
*
* ----------------------------------------------------------------------------  

       SUBROUTINE Interplanetary ( R1,R2,Rbo,Rimpact,Mu1,Mut,Mu2,
     &                             Deltav1,Deltav2,Vbo,Vretro )
        IMPLICIT NONE
        REAL*8 R1,R2,Rbo,Rimpact,Mu1,Mut,Mu2,Deltav1,Deltav2,vbo,
     &         vretro

        INCLUDE 'astmath.cmn'

* ----------------------------  Locals  -------------------------------
        REAL*8  SME1,SME2,SMEt, Vcs1, Vcs2, Vt1, Vt2, A

       !  Find a, SME, apogee .and. perigee velocities of trans orbit --
        A   = (R1+R2) * 0.5D0
        SMEt= -Mut/ (2.0D0*A)
        Vt1 = DSQRT( 2.0D0*( (Mut/R1) + SMEt ) )
        Vt2 = DSQRT( 2.0D0*( (Mut/R2) + SMEt ) )

       ! --- Find circular velocities of launch .and. target planet ----
        Vcs1= DSQRT( Mut/R1 )
        Vcs2= DSQRT( Mut/R2 )

       ! --- Find delta velocities DO Hohmann transfer portion  -----
        DeltaV1= DABS( Vt1 - Vcs1 ) 
        DeltaV2= DABS( Vcs2 - Vt2 ) 

       !  Find SME .and. burnout/impact vel of launch / target planets -  
        SME1  = Deltav1*DeltaV1 * 0.5D0 
        SME2  = Deltav2*DeltaV2 * 0.5D0 
        Vbo   = DSQRT( 2.0D0*( (Mu1/Rbo) + SME1 ) ) 
        Vretro= DSQRT( 2.0D0*( (Mu2/Rimpact) + SME2 ) ) 

c        TP= Pi*DSQRT( a*a*a/Mut )   ! Transfer Period in secs

       RETURN
       END

       
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE PKEPLER
*
*  this subroutine propagates a satellite's position and velocity vector over
*    a given time period accounting for perturbations caused by J2.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Ro          - original position vector       km
*    Vo          - original velocity vector       km/s
*    NDot        - Time rate of change of n       rad/s
*    NDDot       - Time accel of change of n      rad/s2
*    DtSec        - Change in time                 s
*
*  Outputs       :
*    R           - updated position vector        km
*    V           - updated velocity vector        km/s
*
*  Locals        :
*    P           - Semi-paramter                  km
*    A           - semior axis                    km
*    Ecc         - eccentricity
*    incl        - inclination                    rad
*    Argp        - argument of periapsis          rad
*    ArgpDot     - change in argument of perigee  rad/s
*    Omega       - longitude of the asc node      rad
*    OmegaDot    - change in Omega                rad
*    E0          - eccentric anomaly              rad
*    E1          - eccentric anomaly              rad
*    M           - mean anomaly                   rad/s
*    MDot        - change in mean anomaly         rad/s
*    ArgLat      - argument of latitude           rad
*    ArgLatDot   - change in argument of latitude rad/s
*    TrueLon     - true longitude of vehicle      rad
*    TrueLonDot  - change in the true longitude   rad/s
*    LonPerg     - longitude of periapsis         rad
*    LonPeroDot  - longitude of periapsis change  rad/s
*    N           - mean angular motion            rad/s
*    NUo         - true anomaly                   rad
*    J2oP2       - J2 over p sqyared
*    Sinv,Cosv   - Sine .and. Cosine of Nu
*
*  Coupling:
*    rv2coe       - Orbit Elements from position .and. Velocity vectors
*    coe2rv       - Position .and. Velocity Vectors from orbit elements
*    NEWTONM     - Newton Rhapson to find Nu .and. Eccentric anomaly
*
*  References    :
*    Vallado       2007, 687, Alg 64
*
* ------------------------------------------------------------------------------

      SUBROUTINE PKEPLER       ( Ro,Vo,nDot,Nddot,DtSec,     R,V     )
        IMPLICIT NONE
        REAL*8 Ro(3), Vo(3), nDot,nDDot,DtSec, R(3), V(3)
* -----------------------------  Locals  ------------------------------
        REAL*8 Tndto3,p, a, Ecc, Incl, Omega, Argp, Nu, M, ArgLat,
     &         TrueLon, LonPerg, OmegaDot, E0, ArgpDot, MDot,ArgLatDot,
     &         TrueLonDot, LonPerDot, n, J2oP2, J2

        INCLUDE 'astmath.cmn'
        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
        J2    =  0.00108263D0

        CALL rv2coe( Ro,Vo,  p,a,Ecc,incl,Omega,Argp,Nu,M,ArgLat,
     &              TrueLon,LonPerg )
        n= DSQRT(mu/(A*A*A))

        ! ------------ Find the value of J2 perturbations -------------
        J2oP2   = (n*rekm*rekm*1.5D0*J2) / (p*p)
*     NBar    = n*( 1.0D0 + J2oP2*DSQRT(1.0D0-Ecc*Ecc)* (1.0D0 - 1.5D0*DSIN(Incl)*DSIN(Incl)) )
        OmegaDot= -J2oP2 * DCOS(Incl)
        ArgpDot =  J2oP2 * (2.0D0-2.5D0*DSIN(Incl)*DSIN(Incl))
        MDot    =  N

        Tndto3= 2.0D0*NDot*DtSec / (3.0D0*n)
        a     = a - Tndto3 * a
*     edot  = -Tndto3 * (1.0D0-Ecc)/DtSec
        Ecc   = Ecc - Tndto3 * (1.0D0-Ecc)
        p     = a*(1.0D0 - Ecc*Ecc)

        ! ---- Update the orbital elements DO each orbit type ---------
        IF ( Ecc .lt. Small ) THEN
           ! --------------  Circular Equatorial  ---------------------
           IF ( (incl .lt. Small) .or. (DABS(Incl-Pi).lt.Small) ) THEN
               TrueLonDot= OmegaDot + ArgpDot + MDot
               TrueLon   = TrueLon  + TrueLonDot * DtSec
               TrueLon   = DMOD(TrueLon, TwoPi)
             ELSE
               ! ---------------  Circular Inclined    ----------------
               Omega    = Omega + OmegaDot * DtSec
               Omega    = DMOD(Omega, TwoPi)
               ArgLatDot= ArgpDot + MDot
               ArgLat   = ArgLat + ArgLatDot * DtSec
               ArgLat   = DMOD(ArgLat, TwoPi)
             ENDIF
           ELSE
             ! ----- Elliptical, Parabolic, Hyperbolic Equatorial -----
             IF ( ( incl .lt. Small ) .or.
     &            ( DABS(Incl-Pi) .lt. Small ) ) THEN
                  LonPerDot= OmegaDot + ArgpDot
                  LonPerg  = LonPerg + LonPerDot * DtSec
                  LonPerg  = DMOD(LonPerg, TwoPi)
                  M        = M + MDOT*DtSec +NDot*DtSec**2 +
     &                       NDdot*DtSec**3
                  M        = DMOD(M, TwoPi)
                  CALL NEWTONM( Ecc,M,  e0,Nu )
                ELSE
                  ! ---- Elliptical, Parabolic, Hyperbolic Inclined ---
                  Omega= Omega + OmegaDot * DtSec
                  Omega= DMOD(Omega, TwoPi)
                  Argp = Argp  + ArgpDot  * DtSec
                  Argp = DMOD(Argp, TwoPi)
                  M    = M + MDOT*DtSec + NDot*DtSec**2 + NDdot*DtSec**3
                  M    = DMOD(M, TwoPi)
                  CALL NEWTONM( Ecc,M,  e0,Nu )
                ENDIF
           ENDIF

         ! ------------ Use coe2rv to find new vectors -----------------
         CALL coe2rv(P,Ecc,Incl,Omega,Argp,Nu,ArgLat,TrueLon,LonPerg,
     &               R,V)
      RETURN
      END  ! SUBROUTINE PKEPLER

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE J2DRAGPERT
*
*  this subroutine calculates the perturbations for the PREDICT problem
*    involving secular rates of change resulting from J2 .and. Drag only.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    incl        - Inclination                    rad
*    Ecc         - Eccentricity
*    N           - Mean Motion                    rad/s
*    NDot        - Mean Motion rate               rad / 2TU2
*
*  Outputs       :
*    OmegaDot    - Long of Asc Node rate          rad / s
*    ArgpDot     - Argument of perigee rate       rad / s
*    EDot        - Eccentricity rate              / s
*
*  Locals        :
*    P           - Semiparameter                  km
*    A           - Semimajor axis                 km
*    NBar        - Mean Mean motion               rad / s
*
*  Coupling      :
*    None
*
*  References    :
*    Vallado       2007, 645, Eq 9-37, 646, Eq 9-39, 648, Eq 9-50
*
* ------------------------------------------------------------------------------

      SUBROUTINE J2DragPert ( Incl,Ecc,N,NDot,  OmegaDOT,ArgpDOT,EDOT )
        IMPLICIT NONE
        REAL*8 Incl,Ecc,N,NDot, OmegaDOT,ArgpDOT,EDOT

* -----------------------------  Locals  ------------------------------
        REAL*8 P,A,J2,NBar

        ! --------------------  Implementation   ----------------------
        J2  =  0.00108263D0

        a   = (1.0D0/n) ** (2.0D0/3.0D0)
        p   = a*(1.0D0 - ecc**2)
        NBar= n*( 1.0D0+1.5D0*J2*(DSQRT(1.0D0-Ecc*Ecc)/(p*p))*
     &          ( 1.0D0-1.5D0*DSIN(Incl)**2 ))

* ------------------------- Find dot Terms  ---------------------------
        OmegaDot = -1.5D0*( J2/(p*p) ) * DCOS(Incl) * NBar
        ArgpDot  =  1.5D0*( J2/(p*p) ) * (2.0D0-2.5D0*DSIN(Incl)**2) *
     &                    Nbar
        EDot     = -(4.0D0/3.0D0) * (1.0D0-Ecc) * (NDot/Nbar)

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE PREDICT
*
*  this subroutine determines the azimuth .and. elevation DO the viewing
*    of a satellite from a known ground SITE.  Notice the Julian Date is left
*    in it's usual DAYS format because the DOT terms are input as radians per
*    day, thus no extra need DO conversion. Setup with vectors to simplify the
*    use with any propagator.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    JD          - Julian Date of desired obs     Day
*    Latgd       - Geodetic Latitude of SITE      -Pi to Pi rad
*    LST         - Local SIDEREAL Time            -2Pi to 2Pi rad
*    R           - updated position vector        km
*    V           - updated velocity vector        km/s
*    RS          - IJK SITE Vector                km
*    WhichKind   - Type of Sunrise                'S''C''N''A'
*
*  OutPuts       :
*    Rho         - Range from SITE to satellite   km
*    Az          - Azimuth                        rad
*    El          - Elevation                      rad
*    TRtAsc      - Topo Right ascension           rad
*    TDecl       - Topo Declination               rad
*    Vis         - Visibility
*                  'Radar SUN'   - both in SUN
*                  'Eye'  - SITE dark, sat in SUN
*                  'Radar Nite'  - both dark
*                  'Not Visible' - sat below horz
*  Locals        :
*    Temp        - Temporary Real value
*    SRtAsc      - Suns Right ascension           rad
*    SDecl       - Suns Declination               rad
*    SatAngle    - ANGLE between IJK SUN .and. Sat  rad
*    Dist        - Ppdculr dist of sat from RSun  km
*    rr          - Range rate
*    Drr         - Range acceleration
*    Dtrtasc     - Topocentric rtasc rate
*    DRho        - Slant range rate
*    DAz         - Azimuth rate
*    Del         - Elevation rate
*    SunAngle    - ANGLE between SUN .and. SITE     rad
*    AngleLimit  - ANGLE DO twilight conditions  rad
*    RhoVec      - SITE to sat vector in SEZ      km
*    TempVec     - Temporary vector
*    RHoV        - SITE to sat vector in IJK      km
*    RSun        - SUN vector                     AU
*    C           - Temporary Vector
*
*  Coupling      :
*    SUN         - Position vector of SUN
*    CROSS       - CROSS Product of two vectors
*    ROT2,ROT3   - Rotations about 2nd .and. 3rd axis
*    LNCOM1      - Combination of a vector .and. a scalar
*    LNCOM2      - Combination of two vectors .and. two scalars
*    RV_RAZEL    - Conversion with vectors .and. range azimuth elevation
*    RV_TRADEC   - Conversion with topocentric right ascension declination
*    ANGLE       - ANGLE between two vectors
*
*  References    :
*    Vallado       2007, 900, Alg 73, Ex 11-6
*
* ------------------------------------------------------------------------------

      SUBROUTINE PREDICT       ( JD,latgd,LST, r,v,rs, WhichKind,
     &                           Rho,Az,El,tRtasc,tDecl, Vis )
        IMPLICIT NONE
        REAL*8 JD, Latgd, LST, r(3),v(3),RS(3), Rho, Az, El, trtasc,
     &         tdecl
        CHARACTER WhichKind
        CHARACTER*11 Vis
        EXTERNAL MAG
* -----------------------------  Locals  ------------------------------
        REAL*8 RhoVec(3), TempVec(3), RhoV(3), RSun(3), C(3), MAG,
     &         rr, drr, dtrtasc, dtdecl, SRtAsc, SDecl, Dist,
     &         drho, daz, del, SunAngle, SatAngle, AngleLimit,
     &         magc, magrsun,magr

        INCLUDE 'astmath.cmn'
        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
        Az    =  0.0D0
        El    =  0.0D0
        Rho   =  0.0D0
        TRtAsc=  0.0D0
        TDecl =  0.0D0

        ! ------ Find IJK range vector from SITE to satellite ---------
        CALL LNCOM2( 1.0D0,-1.0D0,R,RS,  RhoV )
        Rho= MAG(RhoV)

        ! ------- Calculate Topocentric Rt Asc .and. Declination ------
        CALL RV_TRADEC(r,v,rs,'TOO',rr,trtasc,tdecl,Drr,Dtrtasc,Dtdecl)

        ! ---------------------- Rotate to SEZ ------------------------
        CALL ROT3( RhoV,       LST   ,  TempVec )
        CALL ROT2( TempVec,HalfPi-Latgd,   RhoVec )

        ! --------------- Check visibility constraints ----------------
        ! ------------------ Is it above the Horizon ------------------
        IF ( RhoVec(3) .gt. 0.0D0 ) THEN
            ! --------- Is the SITE in the LIGHT, or the dark? --------
            CALL SUN( JD,RSun,SRtAsc,SDecl )
            CALL LNCOM1( AUER,RSun, RSun )
            CALL ANGLE( RSun,RS, SunAngle )
            IF ( WHICHKind .eq.'S') THEN
               AngleLimit= (90.0D0 + 50.0D0/60.0D0)*Deg2Rad
              ELSE
                IF ( WHICHKind .eq.'C') THEN
                    AngleLimit=  96.0D0*Deg2Rad
                  ELSE
                    IF ( WHICHKind .eq.'N') THEN
                        AngleLimit= 102.0D0*Deg2Rad
                      ELSE
                        IF ( WHICHKind .eq.'A') THEN
                            AngleLimit= 108.0D0*Deg2Rad
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF

            IF ( SunAngle .lt. AngleLimit ) THEN
                Vis= 'Day        '
              ELSE
                ! ---------- This assumes a conical shadow ------------
                ! ----- Is the satellite in the shadow .or. not? ------
                CALL CROSS( RSun, R, C )
                Magc = MAG(c)
                Magr = MAG(r)
                Magrsun = MAG(rsun)
                SatAngle= DASIN( magc/ (magrsun*magr) )
                Dist= magr*DCOS( SatAngle - HalfPi )
                IF ( Dist .gt. 1.0D0 ) THEN
                    Vis= 'Terminator '
                  ELSE
                    Vis= 'Night      '
                  ENDIF
              ENDIF
          ELSE
            Vis= 'not visible'
          ENDIF      

        ! -----------  Calculate Azimuth .and. Elevation  -------------
c        CALL RV_RAZEL( Reci,Veci,Latgd,Lon,alt,TTT,jdut1,lod,
c     &                         xp,yp,terms, 'TOO',
c     &                         Rho,Az,El,DRho,DAz,DEl )
c need to define transformation coordinates.


c old way        CALL RV_RAZEL( r,v,rs,latgd,LST,'TOO', rho,az,el,drho,daz,del )
      RETURN
      END ! Predict

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE DERIV
*
*  this subroutine calculates the derivative of the two-body state vector
*    use with the Runge-Kutta algorithm.  Note time is not needed.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    X           - State Vector                   km, km/s
*
*  Outputs       :
*    XDot        - Derivative of State Vector     km/s,  km/TU2
*
*  Locals        :
*    RCubed      - Cube of R
*
*  Coupling      :
*    None.
*
*  References    :
*    None.
*
* ------------------------------------------------------------------------------

      SUBROUTINE Deriv ( X,  XDot )
        IMPLICIT NONE
        REAL*8 X(6), XDot(6)

* -----------------------------  Locals  ------------------------------
        Real*8 RCubed

        ! --------------------  Implementation   ----------------------
        RCubed= ( DSQRT( X(1)**2 + X(2)**2 + X(3)**2 ) )**3

        write(*,*) rcubed,'   rcubed',x(4)
        write(*,*) x(5),'  ',x(6)
        ! -----------------  Velocity Terms  --------------------------
        XDot(1)= X(4)
        XDot(2)= X(5)
        XDot(3)= X(6)

        ! --------------  Acceleration Terms   ------------------------
        XDot(4)= -X(1) / RCubed
        XDot(5)= -X(2) / RCubed
        XDot(6)= -X(3) / RCubed

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE InitGravityField
*
*  this subroutine reads .and. stores the gravity field DO use in the program.
*    coefficients. The routine can be configured DO either normalized .or.
*    unnormalized values.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Order       - Size of gravity field          1..70D0
*
*  Outputs       :
*    C           - Gravitational Coefficients
*    S           - Gravitational Coefficients
*
*  Locals        :
*    RCubed      - Cube of R
*
*  Coupling      :
*    None.
*
*  References    :
*    None.
*
* ------------------------------------------------------------------------------

      SUBROUTINE InitGravityField   ( Order, C,S )
        IMPLICIT NONE
        INTEGER Order
        REAL*8 C(70,70), S(70,70)

        INTEGER l, m, cexp, sexp
        REAL*8 cnor, snor, Cval, Sval
        CHARACTER*8 NOTEOF

        ! --------------------  Implementation   ----------------------
        OPEN( UNIT=25,File='c:/tplib/wgs84.dat',STATUS='OLD' )

* --------------- Set up Loop to READ through Input File --------------
        NOTEOF = 'TRUE'
        DO WHILE (NOTEOF.eq.'TRUE')
            Read(25,*,END=999) l,m,Cnor,Snor
            C(l,m)= Cnor  ! unnormalized values
            S(l,m)= Snor
          ENDDO  ! While Not EOF ( END=999 in Read )

 999    CLOSE( 25 )

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE LegPoly
*
*  this subroutine finds the Legendre polynomials DO the gravity field.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Latgc       - Geocentric Latitude of SITE    -Pi to Pi rad
*    Order       - Size of gravity field          1..70D0
*
*  Outputs       :
*    LArr        - Array of Legendre Polynomials
*
*  Locals        :
*    L,m         - Indices of gravitational potential
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 593, Eq 8-56
*
* ------------------------------------------------------------------------------

      SUBROUTINE LegPoly     ( Latgc, Order, LArr )
        IMPLICIT NONE
        REAL*8 Latgc, LArr(0:70,0:70)
        INTEGER Order

        INTEGER L,m

        ! --------------------  Implementation   ----------------------
        LArr(0,1)= 0.0D0
        LArr(0,0)= 1.0D0
        LArr(1,0)= DSIN(Latgc)
        LArr(1,1)= DCOS(Latgc)

        ! ------------------- Perform Recursions ----------------------
        DO L= 2,Order
            LArr(0,L-1)= 0.0D0
            DO m= 0,L
                IF ( m .eq. 0 ) THEN
                      LArr(L,0)= ( (2*L-1)* LArr(1,0) * LArr(L-1,0)
     &                           - (L-1)* LArr(L-2,0) )/L
                  ELSE
                    IF ( m .eq. L ) THEN
                        LArr(L,m)= (2*L-1) * LArr(1,1) * LArr(L-1,m-1)
                      ELSE
                        LArr(L,m)= LArr(L-2,m)
     &                           + (2*L-1) * LArr(1,1) * LArr(L-1,m-1)
                      ENDIF
                  ENDIF  
              ENDDO   ! DO m
          ENDDO   ! DO L

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE FullGeop
*
*  this subroutine finds the Legendre polynomial value DO the gravity field
*    DO a given order.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Order       - Size of gravity field          1..70D0
*
*  Outputs       :
*    C           - Gravitational Coefficients
*    S           - Gravitational Coefficients
*
*  Locals        :
*    RCubed      - Cube of R
*
*  Coupling      :
*    IJKTOLATLONA- Find sub satellite point
*
*  References    :
*    Vallado       2007,
*
* ------------------------------------------------------------------------------

      SUBROUTINE FullGeop    ( R,V, ITime,WhichOne,BC,Order,C,S,APert )
        IMPLICIT NONE
        REAL*8 R(3), V(3), ITime,BC,C(70,70),S(70,70),APert(6)
        INTEGER Order, WhichOne
        EXTERNAL MAG
        INTEGER L, m
        REAL*8 LArr(0:70,0:70), OORDelta, Temp, OOr,  SumM1, SumM2,MAG,
     &         SumM3, DistPartr, DistPartPhi, DistPartLon, RDelta,
     &         Latgc, latgd, hellp, Lon, LastOOr, SumL1,SumL2,magr
        ! --------------------  Implementation   ----------------------
        CALL IJK2llA( R, latgc,latgd,lon,hellp )

        ! -------------------- Find Legendre polynomials --------------
        CALL LegPoly( Latgc,Order, LArr )

        ! --------- Partial derivatives of disturbing potential -------
        magr = MAG(r)
        OOr= 1.0D0/magr
        LastOOr= 1.0D0/magr
        SumM1= 0.0D0
        SumM2= 0.0D0
        SumM3= 0.0D0
        DO L= 2,Order
            DO m= 0,L
                SumM1= SumM1 + LArr(L,m) * C(L,m)*DCOS(m*Lon)
     &                       + S(L,m)*DSIN(m*Lon)
                SumM2= SumM2 + C(L,m)*DCOS(m*Lon)
     &                       + S(L,m)*DSIN(m*Lon) *
     &                       ( LArr(L,m+1) - LArr(L,m)*m*DTAN(Latgc) )
                SumM3= SumM3 + m*LArr(L,m) * (S(L,m)*DCOS(m*Lon)
     &                       - C(L,m)*DSIN(m*Lon))
                SumL1 = 0.0D0 ! fixnjnjnjjjjnjnjnjn
                SumL2 = 0.0D0 ! fix
              ENDDO
          ENDDO

        DistPartR  = -OOr*OOr*SumL1 * SumM1
        DistPartPhi=  OOr*SumL2     * SumM2
        DistPartLon=  OOr*SumL2     * SumM3

        ! --------- Non-spherical pertubative acceleration ------------
        RDelta  = DSQRT( r(1)*r(1) + r(2)*r(2) )
        OORdelta= 1.0D0/RDelta
        Temp    = OOr*DistPartR - r(3)*OOr*OOr*OORDelta*DistPartPhi

        APert(1)= Temp*r(1) - OORDelta*DistPartLon*r(2)
        APert(2)= Temp*r(2) + OORDelta*DistPartLon*r(1)
        APert(3)= OOr*DistPartR*r(3) + OOR*OOr*RDelta*DistPartPhi

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE PERTACCEL
*
*  this subroutine calculates the actual value of the perturbing acceleration.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R           - Radius vector                  km
*    V           - Velocity vector                km/s
*    Time        - Initial time (Julian Date)     Days from 4713 BC
*    WhichOne    - Which perturbation to calc     1 2 3 4 5 ...
*    BC          - Ballistic Coefficient          kg/m2
*
*  Outputs       :
*    APert       - Perturbing acceleration        km/TU2
*
*  Locals        :
*    rs2         - SUN radius vector **2
*    rs3         - SUN radius vector **3
*    rm2         - MOON radius vector **2
*    rm3         - MOON radius vector **3
*    r32         - "z" component of Radius vec **2
*    r33         - "z" component of Radius vec **3
*    r34         - "z" component of Radius vec **4
*    r2          - Radius vector **2
*    r3          - Radius vector **3
*    r4          - Radius vector **4
*    r5          - Radius vector **5
*    r7          - Radius vector **7
*    Beta        -
*    Temp        - Temporary Real Value
*    rho         - Atmospheric Density
*    Va          - Relative Velocity Vector       km / s
*    RSun        - Radius Vector to SUN           AU
*    RMoon       - Radius Vector to MOON          km
*    RtAsc       - Right Ascension                rad
*    Decl        - Declination                    rad
*    AUER        - Conversion from AU to km
*    Temp1       -
*    Temp2       -
*    i           - Index
*
*  Coupling      :
*    MAG         - Magnitude of a vector
*    SUN         - SUN vector
*    MOON        - MOON vector
*    ATMOS       - Atmospheric density
*
*  References    :
*
* ------------------------------------------------------------------------------

      SUBROUTINE PERTACCEL   ( R,V, ITime, WhichOne, BC, APert )
        IMPLICIT NONE
        REAL*8 R(3), V(3), ITime,BC,APert(6)
        INTEGER WhichOne
        EXTERNAL MAG
* -----------------------------  Locals  ------------------------------
        INTEGER i
* fix the c and s vars
        REAL*8 C(70), S(70)

        REAL*8 Va(3), RSun(3), RMoon(3), rs2, rm2, rs3, rm3, r32, r33,
     &         r34, r2, r3, r4, r5, r7, Beta, Temp, rho, srtasc, magva,
     &         sdecl, mrtasc, mdecl, Temp1, Temp2, J2, J3, magapert,
     &         J4, GMS, GMM, DOT, Mag, magr, magv, magrsun, magrmoon
        EXTERNAL DOT

        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
        J2         =    0.00108263D0
        J3         =   -0.00000254D0
        J4         =   -0.00000161D0
        GMS        =    3.329529364D5
        GMM        =    0.01229997D0
        magr = MAG( R )
        magv = MAG( V )
        R2 = magr*magr
        R3 = R2*magr
        R4 = R2*R2
        R5 = R2*R3
        R7 = R5*R2
        R32= r(3)*r(3)
        R33= R32*r(3)
        R34= R32*R32

        ! -----------------   J2 Acceleration   -----------------------
        IF ( WhichOne .eq. 1 ) THEN
                Temp1=  (-1.5D0*J2) / R5
                Temp2=  1.0D0 - (5.0D0*R32) / R2
                APert(1)= Temp1*r(1) * Temp2  ! recheck with formulae
                APert(2)= Temp1*r(2) * Temp2
                APert(3)= Temp1*r(3) * ( 3.0D0-(5.0D0*R32) / R2 )
              ENDIF

        ! ------------------   J3 Acceleration   ----------------------
        IF ( WhichOne .eq. 2 ) THEN
                Temp1=  (-2.5D0*J3) / R7
                Temp2=  3.0D0*r(3)-(7.0D0*R33) / R2
                APert(1)= Temp1*r(1) * Temp2
                APert(2)= Temp1*r(2) * Temp2
                IF ( DABS( r(3) ) .gt. 0.0000001D0 ) THEN
                    APert(3)= Temp1*r(3) * ((6.0D0*r(3))-((7.0D0*R33)
     &                           / R2) - ( (3.0D0*r2) / (5.0D0*r(3)) ))
                  ELSE
                    APert(3)= 0.0D0
                  ENDIF
              ENDIF

        ! -----------------    J4 Acceleration   ----------------------
        IF ( WhichOne .eq. 3 ) THEN
                Temp1=  (-1.875D0*J4) / R7
                Temp2=  1.0D0-((14.0D0*R32)/R2)+((21.0D0*R34) / R4)
                APert(1)= Temp1*r(1) * Temp2
                APert(2)= Temp1*r(2) * Temp2
                APert(3)= Temp1*r(3) * (5.0D0-((70.0D0*R32)/(3.0D0*R2))
     &                    +((21.0D0*R34) / R4 ))
              ENDIF

        ! -----------------   SUN Acceleration   ----------------------
        IF ( WhichOne .eq. 4 ) THEN
                CALL SUN( ITime,RSun,SRtAsc,SDecl )
                DO i= 1,3
                    RSun(i)= RSun(i)*AuER    ! chg AU's to km's
                  ENDDO

                magrsun = MAG(rsun)

                RS2= magrsun*magrsun
                RS3= RS2*magrsun
                Temp= DOT( R,RSun )
                Temp1= -GMS/RS3
                Temp2= 3.0D0*Temp/RS2
                APert(1)= Temp1 * (r(1) - Temp2*RSun(1))
                APert(2)= Temp1 * (r(2) - Temp2*RSun(2))
                APert(3)= Temp1 * (r(3) - Temp2*RSun(3))
              ENDIF

        ! -----------------  MOON Acceleration   ----------------------
        IF ( WhichOne .eq. 5 ) THEN
                CALL MOON( ITime,RMoon,MRtAsc,MDecl )
                magrmoon = MAG(rmoon)
                RM2= magRMoon**2
                RM3= RM2*magRMoon
                Temp= DOT( R,RMoon )
                Temp1= -GMM/RM3
                Temp2= 3.0D0*Temp/RM2
                APert(1)= Temp1 * (r(1) - Temp2*RMoon(1))
                APert(2)= Temp1 * (r(2) - Temp2*RMoon(2))
                APert(3)= Temp1 * (r(3) - Temp2*RMoon(3))
              ENDIF

        ! -----------------  Drag Acceleration   ----------------------
        IF ( WhichOne .eq. 6 ) THEN
                Va(1)= V(1) + (OmegaEarth*r(2))   ! km/s
                Va(2)= V(2) - (OmegaEarth*r(1))
                Va(3)= V(3)
                magva = MAG( Va )

                CALL ATMOS( R, Rho )

                Temp= -1000.0D0 * magVa * 0.5D0*Rho* ( 1.0D0/BC )*
     &                 6378137.0D0
                APert(1)= Temp*Va(1)
                APert(2)= Temp*Va(2)
                APert(3)= Temp*Va(3)
              ENDIF

        ! ----------------- Solar Acceleration   ----------------------
        IF ( WhichOne .eq. 7 ) THEN
                CALL SUN( ITime,RSun,SRtAsc,SDecl )

                Beta= 0.4D0                          ! reflectivity
                magAPert= (Beta*2.0D0*4.51D-06)/BC   ! assume Csr = 2.0D0
                Temp= -magAPert/magrsun
                APert(1)= Temp*RSun(1)
                APert(2)= Temp*RSun(2)
                APert(3)= Temp*RSun(3)
              ENDIF

        ! ------------------- Square Gravity Field --------------------
        IF ( WhichOne .eq. 10 ) THEN

            CALL FullGeop( R,V,ITime,WhichOne,BC,50,C,S,APert )
          ENDIF

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE PDERIV
*
*  this subroutine calculates the derivative of the state vector DO use with
*    the Runge-Kutta algorithm.  The DerivType string is used to determine
*    which perturbation equations are used.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    ITime       - Initial Time (Julian Date)     Days from 4713 BC
*    X           - State Vector                   km  ,  km/s
*    DerivType   - String of which perts to incl  'Y' .and. 'N'
*                  Options are in order : J2, J3,
*                  J4, SUN, MOON, Drag, Solarrad
*    BC          - Ballistic Coefficient          kg/m2
*
*  Outputs       :
*    XDot        - Derivative of State Vector     km/s, km/TU2
*
*  Locals        :
*    RCubed      - Radius vector cubed            ER3
*    Ro          - Radius vector                  km
*    Vo          - Velocity vector                km/s
*    APert       - Perturbing acceleration        km/TU2
*    TempPert    - Temporary acceleration         km/TU2
*    i           - Index
*
*  Coupling      :
*    PERTACCEL   - Calculates the actual values of each perturbing acceleration
*    ADDVEC      - Adds two vectors together
*    MAG         - Magnitude of a vector
*
*  References    :
*
* ------------------------------------------------------------------------------

      SUBROUTINE PDERIV( ITime,X,DerivType,BC,  XDot )
        IMPLICIT NONE
        REAL*8 X(6),XDot(6),ITime,BC
        CHARACTER*10 DerivType

* -----------------------------  Locals  ------------------------------
        REAL*8 RCubed,Ro(3),Vo(3),APert(3),TempPert(3), magr, magv, MAG
        INTEGER i

        EXTERNAL MAG
        
        ! --------------------  Implementation   ----------------------
        DO i= 1, 3
            APert(i)= 0.0D0
            Ro(i)   = X(i)
            Vo(i)   = X(i+3)
        ENDDO
        magr = MAG( Ro )
        magv = MAG( Vo )
*        APert(4)= 0.0D0
        RCubed = magr**3

* -------------------------  Velocity Terms  --------------------------
        XDot(1)= X(4)
        XDot(2)= X(5)
        XDot(3)= X(6)

* ----------------------  Acceleration Terms  -------------------------
        IF ( DerivType(1:1).eq.'Y' ) THEN
            CALL PertAccel( Ro,Vo,ITime,1,BC, APert )
          ENDIF
        IF ( DerivType(2:2).eq.'Y' ) THEN
            CALL PertAccel( Ro,Vo,ITime,2,BC, TempPert )
            CALL AddVec( TempPert,APert,APert )
          ENDIF
        IF ( DerivType(3:3).eq.'Y' ) THEN
            CALL PertAccel( Ro,Vo,ITime,3,BC, TempPert )
            CALL AddVec( TempPert,APert,APert )
          ENDIF
        IF ( DerivType(4:4).eq.'Y' ) THEN
            CALL PertAccel( Ro,Vo,ITime,4,BC, TempPert )
            CALL AddVec( TempPert,APert,APert )
          ENDIF
        IF ( DerivType(5:5).eq.'Y' ) THEN
            CALL PertAccel( Ro,Vo,ITime,5,BC, TempPert )
            CALL AddVec( TempPert,APert,APert )
          ENDIF
        IF ( DerivType(6:6).eq.'Y' ) THEN
            CALL PertAccel( Ro,Vo,ITime,6,BC, TempPert )
            CALL AddVec( TempPert,APert,APert )
          ENDIF
        IF ( DerivType(7:7).eq.'Y' ) THEN
            CALL PertAccel( Ro,Vo,ITime,7,BC, TempPert )
            CALL AddVec( TempPert,APert,APert )
          ENDIF
        ! -------------------- new full gravity field -----------------
        IF ( DerivType(10:10) .eq. 'Y' ) THEN
            CALL PERTACCEL( Ro,Vo,ITime,10,BC, TempPert )
            CALL ADDVEC( TempPert,APert,APert )
          ENDIF

        XDot(4)= (-X(1) / RCubed) + APert(1)
        XDot(5)= (-X(2) / RCubed) + APert(2)
        XDot(6)= (-X(3) / RCubed) + APert(3)

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                                SUBROUTINE RK4
*
*  this subroutine is a fourth order Runge-Kutta integrator DO a 6 dimension
*    First Order differential equation.  The intended use is DO a satellite
*    equation of motion.  The user must provide an external SUBROUTINE containing
*    the system Equations of Motion.  Notice time is included since some
*    applications in PDERIV may need this.  The LAST position in DerivType is a
*    flag DO two-body motion.  Two-Body motion is used IF ( the 10th element is
*    set to '2', otherwise the Yes .and. No values determine which perturbations
*    to use. Be careful with the units. The ITime parameter comes as a JD because
*    it's used DO SUN/MOON calcs later on. The stepsize is changed to s
*    operations involving the state vector because it is canonical.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*                                 12 Nov 1993 - fix DO s, DT, etc.
*  Inputs          Description                    Range / Units
*    ITime       - Initial Time (Julian Date)     Days from 4713 BC
*    DtDay       - Step size                      Day
*    XDot        - Derivative of State Vector
*    DerivType   - String of which perts to incl  'Y' .and. 'N'
*                  Options are in order : J2, J3,
*                  J4, SUN, MOON, Drag, Solarrad
*    BC          - Ballistic Coefficient          kg/m2
*    X           - State vector at initial time   km, km/s
*
*  Outputs       :
*    X           - State vector at new time       km, km/s
*
*  Locals        :
*    K           - Storage DO values of state
*                   vector at different times
*    Temp        - Storage DO state vector
*    TempTime    - Temporary time storage half
*                   way between DtDay             Day
*    J           - Index
*    DtSec        - Step size                      s
*
*  Coupling      :
*    DERIV       - SUBROUTINE for Derivatives of EOM
*    PDeriv      - SUBROUTINE for Perturbed Derivatives of EOM
*
*  References    :
*    Vallado       2007, 526
*
* ------------------------------------------------------------------------------

      SUBROUTINE RK4           ( ITime,DtDay,XDot,DerivType,BC,  X )
        IMPLICIT NONE
        REAL*8  DtDay, ITime, X(6), XDot(6), BC
        CHARACTER*10 DerivType

* -----------------------------  Locals  ------------------------------
        REAL*8 K(6,3), TEMP(6), TempTime, DtSec, TUDay
        INTEGER J

        TUDay = 0.00933809017716D0
          DtSec= DtDay/TUDay
          ! --------- Evaluate 1st Taylor Series Term -----------------
          IF (DerivType(10:10).eq.'2') THEN
              CALL DERIV( X,XDot )
            ELSE
              CALL PDERIV( ITime,X,DerivType,BC,XDot )
            ENDIF

          TempTime = ITime + DtDay*0.5D0

          ! -------- Evaluate 2nd Taylor Series Term ------------------
          DO J = 1,6
             K(J,1)  = DtSec * XDot(J)
             TEMP(J) = X(J) + 0.5D0*K(J,1)
          ENDDO

          IF (DerivType(10:10).eq.'2') THEN
              CALL DERIV( Temp,XDot )
            ELSE
              CALL PDERIV( TempTime,Temp,DerivType,BC,XDot )
            ENDIF

          ! ---------- Evaluate 3rd Taylor Series Term ----------------
          DO J = 1,6
             K(J,2)  = DtSec * XDot(J)
             TEMP(J) = X(J) + 0.5D0*K(J,2)
          ENDDO

          IF (DerivType(10:10).eq.'2') THEN
              CALL DERIV( Temp,XDot )
            ELSE
              CALL PDERIV( TempTime,Temp,DerivType,BC,XDot )
            ENDIF

          ! ---------- Evaluate 4th Taylor Series Term ----------------
          DO J = 1,6
             K(J,3)  = DtSec * XDot(J)
             TEMP(J) = X(J) + K(J,3)
          ENDDO

          IF (DerivType(10:10).eq.'2') THEN
              CALL DERIV( Temp,XDot )
            ELSE
              CALL PDERIV( ITime+DtDay,Temp,DerivType,BC,XDot )
            ENDIF

          ! ---- Update the state vector, perform integration  --------
          DO J = 1,6
             X(J) = X(J) + ( K(J,1) + 2.0D0*(K(J,2) + K(J,3)) +
     &                       DtSec*XDot(J) ) / 6.0D0
          ENDDO

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                                SUBROUTINE RKF45
*
*  this subroutine is a fourth order Runge-Kutta-Fehlberg integrator DO a 6-D
*    First Order differential equation.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    ITime       - Initial Time (Julian Date)     Days from 4713 BC
*    DtDay       - Step size                      Day
*    XDot        - Derivative of State Vector
*    DerivType   - String of which perts to incl  'Y' .and. 'N'
*                  Options are in order : J2, J3,
*                  J4, SUN, MOON, Drag, Solarrad
*    BC          - Ballistic Coefficient          kg/m2
*    X           - State vector at initial time   km, km/s
*
*  Outputs       :
*    X           - State vector at new time       km, km/s
*
*  Locals        :
*    K           - Storage DO values of state
*                    vector at different times
*    Temp        - Storage DO state vector
*    TempTime    - Temporary time storage half
*                    way between DtDay            Day
*    J           - Index
*    DtSec        - Step size                      s
*
*  Coupling      :
*    DERIV       - SUBROUTINE DO Derivatives of EOM
*    PDeriv      - SUBROUTINE DO Perturbed Derivatives of EOM
*
*  References    :
*    Vallado       2007, 526
*
* ------------------------------------------------------------------------------  

      SUBROUTINE RKF45       ( ITime, DtDay, XDot,DerivType, BC, X )
        IMPLICIT NONE
        REAL*8 ITime,DtDay, XDot(6),X(6), BC
        CHARACTER*10 Derivtype

* -----------------------------  Locals  ------------------------------
        INTEGER Ktr, J
        REAL*8 K(6,6), Temp(6,1), DtSec, HMin, HMax, TStop, Time, Err,
     &         S, TempTime, Small, TUDay

        Small =     0.000001D0  ! this is pretty sensitive for RF45
        TUDay =     0.0093380913806D0
        HMin = DtDay/64.0D0
        HMax = DtDay*64.0D0
        Time = ITime
        TStop= ITime + DtDay
        DtSec = DtDay/TUDay

        Ktr= 1
        DO WHILE (Time .lt. TStop)
            IF ( Time + DtDay .gt. TStop ) THEN  ! Make sure you END exactly on the step
                DtDay= TStop - Time
              ENDIF

            ! ------------- Evaluate 1st Taylor Series Term -----------
            IF ( DerivType(10:10) .eq. '2' ) THEN
                CALL DERIV( X, XDot )
              ELSE
                CALL PDeriv( Time,X,DerivType,BC, XDot )
              ENDIF

            TempTime= Time + DtDay*0.25D0
            ! ------------- Evaluate 2nd Taylor Series Term -----------
            DO j= 1,6
                K(J,1)   = DtSec * XDot(J)   ! set # 1
                Temp(J,1)= X(J) + 0.25D0 * K(J,1)   !get ready DO 2
              ENDDO
            IF ( DerivType(10:10) .eq. '2' ) THEN
                CALL DERIV( Temp, XDot )
              ELSE
                CALL PDeriv( TempTime,Temp,DerivType,BC,  XDot )
              ENDIF

            TempTime= Time + DtDay*0.375D0
            ! ------------- Evaluate 3rd Taylor Series Term -----------
            DO j= 1,6
                K(J,2)   = DtSec * XDot(J)
                Temp(J,1)= X(J) + 0.09375D0 * K(J,1)
     &                     + 0.28125D0 * K(J,2)
              ENDDO
            IF ( DerivType(10:10) .eq. '2' ) THEN
                CALL DERIV( Temp, XDot )
              ELSE
                CALL PDeriv( TempTime,Temp,DerivType,BC,  XDot )
              ENDIF

            TempTime= Time + DtDay*12.0D0/13.0D0

            ! ------------- Evaluate 4th Taylor Series Term -----------
            DO j= 1,6
                  K(J,3)    = DtSec * XDot(J)
                  Temp(J,1) = X(J) + K(J,1) * 1932.0D0/2197.0D0
     &                        - K(J,2)*7200.0D0/2197.0D0
     &                        + K(J,3)*7296.0D0/2197.0D0
              ENDDO
            IF ( DerivType(10:10) .eq. '2' ) THEN
                CALL DERIV( Temp, XDot )
              ELSE
                CALL PDeriv( TempTime,Temp,DerivType,BC,  XDot )
              ENDIF

            ! ------------- Evaluate 5th Taylor Series Term -----------
            DO j= 1,6
                  K(J,4)    = DtSec * XDot(J)
                  Temp(J,1) = X(J) + K(J,1)* 439.0D0/ 216.0D0
     &                      - K(J,2) * 8.0D0 + K(J,3)*3680.0D0/ 513.0D0
     &                      - K(J,4) * 845.0D0/4104.0D0
              ENDDO
            IF ( DerivType(10:10) .eq. '2' ) THEN
                CALL DERIV( Temp, XDot )
              ELSE
                CALL PDeriv( Time+DtDay,Temp,DerivType,BC,  XDot )
              ENDIF 

            TempTime= Time + DtDay*0.5D0 

            ! ------------- Evaluate 6th Taylor Series Term -----------
            DO j= 1,6
                  K(J,5)    = DtSec * XDot(J)
                  Temp(J,1) = X(J) - K(J,1)*8.0D0/27.0D0
     &                      + K(J,2)* 2.0D0 - K(J,3)*3544.0D0/2565.0D0
     &                      + K(J,4)*1859.0D0/4104.0D0 - K(J,5)*0.275D0
              ENDDO
            IF ( DerivType(10:10) .eq. '2' ) THEN
                CALL DERIV( Temp, XDot )
              ELSE
                CALL PDeriv( TempTime,Temp,DerivType,BC,  XDot )
              ENDIF 

            DO j= 1,6
                  K(J,6)=  DtSec * XDot(J)
              ENDDO

            ! ------------------- Check DO convergence ---------------
            Err= 0.0D0 
            DO j= 1,6
                Err= DABS( K(J,1)*1.0D0/360.0D0
     &                 - K(J,3)*128.0D0/4275.0D0
     &                 - K(J,4)*2197.0D0/75240.0D0
     &                 + K(J,5)*0.02D0 + K(J,6)*2.0D0/55.0D0 )
              ENDDO

            ! ----- Update the State vector, perform integration ------
            IF ( ( Err .lt. Small ) .or.
     &           ( DtDay .le. 2.0D0*HMin+Small ) ) THEN
                DO j= 1,6
                    X(J)= X(J) + K(J,1)*25.0D0/216.0D0
     &                      + K(J,3)*1408.0D0/2565.0D0
     &                      + K(J,4)*2197.0D0/4104.0D0 - K(J,5)*0.2D0
                  ENDDO
                Time= Time + DtDay
                s   = 0.0D0
                Ktr = 1
              ELSE
                S= 0.84D0* (Small*DtDay/Err)**0.25D0
                IF ( ( S.lt.0.75D0 ).and.(DtDay.gt.2.0D0*HMin ) ) THEN  ! Reduce  Step  Size
                    DtDay= DtDay * 0.5D0
                  ENDIF
                IF ( ( S.gt.1.5D0 ).and.(2.0D0*DtDay.lt.HMax ) ) THEN    ! Increase Step Size
                    DtDay= DtDay * 2.0D0
                  ENDIF
                Ktr = Ktr + 1
              ENDIF

          ENDDO   ! WHILE

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                                SUBROUTINE COWELL
*
*  this subroutine uses a fourth order Runge-Kutta integrator on a 6 dimension
*    First Order differential equation.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R           - Initial position vector        km
*    V           - Initial velocity vector        km/s
*    ITime       - Initial Time (Julian Date)     Days from 4713 BC
*    FTime       - Final Time (Julian Date)       Days from 4713 BC
*    DtDay       - Step size                      Day
*    DerivType   - String of which perts to incl  'Y' .and. 'N'
*                  Options are in order : J2, J3,
*                  J4, SUN, MOON, Drag, Solarrad
*    BC          - Ballistic Coefficient          kg/m2
*
*  Outputs       :
*    R1          - Final position vector          km
*    V1          - Final velocity vector          km/s
*
*  Locals        :
*    Time        - Current time during the loop   Days from 4713 BC
*    X           - State vector at each time      km, km/s
*
*  Coupling      :
*    RK4         - Runge-Kutta algorithm
*    MAG         - Magnitude of a vector
*
*  References    :
*
* ------------------------------------------------------------------------------

      SUBROUTINE Cowell   ( R,V,ITime,FTime,DtDay,DerivType,BC, R1,V1 )
        IMPLICIT NONE
        REAL*8 R(3),V(3),ITime,FTime,DtDay,BC,R1(3),V1(3)
        CHARACTER*10 DerivType

* -----------------------------  Locals  ------------------------------
        REAL*8 Time, X(6), XDot(6),FT
        INTEGER i

        DO i= 1, 6
            IF ( i .le. 3 ) THEN
                X(i)= r(i)
              ELSE
                X(i)= v(i-3)
              ENDIF
          ENDDO

        ! --------- Loop through the time interval desired ------------
        Time= ITime
        Ft = FTime
        DO WHILE (Time .le. Ft)
            IF ( Time+DtDay .gt. Ft ) THEN
                DtDay = Ft - Time
                Ft = FTime - 1.0D0
              write(*,*) 'asifgy'
              ENDIF

            CALL RK4( Time,DtDay,XDot,DerivType,BC, X )
         write(*,'(6(f8.4))') (x(i), i=1,6)
*            CALL RKF45( Time,DtDay,XDot,DerivType,BC, X )

            Time = Time + DtDay
         write(*,*) Ft,' ',time,' ',dtday
          ENDDO

        ! ----------------- Update the state vector -------------------
        DO i= 1, 6
            IF ( i .le. 3 ) THEN
                r1(i)  =   X(i)
              ELSE
                v1(i-3)=   X(i)
              ENDIF
          ENDDO
*        CALL MAG( R1 )
*        CALL MAG( V1 )

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE ATMOS
*
*  this subroutine finds the atmospheric density at an altitude above an
*    oblate earth given the position vector in the Geocentric Equatorial
*    frame.  The position vector is in km's .and. the density is in gm/cm**3.
*    DO certain applications, it may not be necessary to find the Hellp
*    exact height difference.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R           - IJK Position vector            km
*
*  Outputs       :
*    Rho         - Density                        kg/m**3
*
*  Locals        :
*    Hellp       - Height above ellipsoid         km
*    OldDelta    - Previous value of DeltaLat     rad
*    Latgd       - Geodetic Latitude              -Pi/2 to Pi/2 rad
*    SinTemp     - Sine of Temp
*    RhoNom      - Nominal density at particular alt      gm/cm**3
*    NextBaseAlt - Next Base Altitude
*    LastBaseAlt - Last Base Altitude
*    H           - Scale Height                   km
*    i           - index
*    AtmosFile   - File of data DO the
*                    exponential atmosphere
*
*  Coupling      :
*    MAG         - Magnitude of a vector
*
*  References    :
*    Vallado       2007, 562, Ex 8-4
*
* ------------------------------------------------------------------------------

      SUBROUTINE ATMOS       ( R,Rho )
        IMPLICIT NONE
        Real*8 R(3),Rho
        EXTERNAL MAG
        INTEGER i
        REAL*8 Hellp, OldDelta, Latgd, SinTemp, c, Decl, Temp, H,
     &         RhoNom, BaseAlt, LastBaseAlt, Small,
     &         LastH, LastRhoNom, MAG, magr

        INCLUDE 'astconst.cmn'

        ! -------------------  Initialize values   --------------------
        Small      =     0.0000001D0
        OPEN( UNIT=25,File='atmosexp.dat',STATUS='OLD' )

        magr = MAG( R )
        Decl = DASIN( R(3) / magr )
        Latgd= Decl

        ! ---- Iterate to find Geocentric .and. Geodetic Latitude  ----
        Temp = DSQRT( R(1)*R(1) + R(2)*R(2) )
        i= 1
        OldDelta = Latgd*2.0D0
        DO WHILE ( ( DABS(OldDelta-Latgd).ge.Small ) .and. (i.lt.10) )
            OldDelta= Latgd
            SinTemp = DSIN( Latgd )
            c       = 1.0D0 / (DSQRT( 1.0D0-eeSqrd*SinTemp*SinTemp ))
            Latgd   = DATAN( (r(3)+c*eeSqrd*SinTemp)/Temp )
            i = i + 1
          ENDDO
        Hellp = ( (Temp/DCOS(Latgd)) - c ) * rekm

        IF ( i .ge. 10 ) THEN
            Write(*,*)  'IJKtoLatLon did NOT converge '
          ENDIF

        ! ---------- Determine density based on altitude --------------

        ! ---------- Set up Loop to READ through Input File -----------
        READ( 25,*,END=999) LastBaseAlt, LastRhoNom, LastH

        DO WHILE (LastBaseAlt .lt. Hellp )
            READ( 25,*,END=999) BaseAlt, RhoNom, H
            IF (BaseAlt .lt. Hellp ) THEN
                LastBaseAlt= BaseAlt
                LastRhoNom = RhoNom
                LastH      = H
              ENDIF
          ENDDO

        RHO   = LastRHONOM * EXP((LastBaseAlt-Hellp)/LastH)

 999    CLOSE( 25 )
      RETURN
      END
       
       
       
       
  
            