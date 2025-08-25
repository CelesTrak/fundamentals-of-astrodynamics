c     -------------------------------------------------------------------------
c
c                                MathTimeLib.for
c
c this library contains various math and basic time routines.
c
c                            companion code for
c               fundamentals of astrodynamics and applications
c                                    2013
c                              by david vallado
c
c               email dvallado@comspoc.com, davallado@gmail.com
c
c    current :
c               9 sep 2022 david vallado
c                           combine from previous versions
c    changes :
c               9 sep 2022 david vallado
c                           original baseline 
c                           
c       ----------------------------------------------------------------      





* ------------------------------------------------------------------------------
*
*                           FUNCTION DCOT
*
*  this function finds the Cotangent of an ANGLE in radians.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    XVal        - ANGLE to take Cotangent of                   rad
*
*  OutPuts       :
*    COT         - Result
*
*  Locals        :
*    Temp        - Temporary Real variable
*
* ------------------------------------------------------------------------------

      REAL*8 FUNCTION DCOT    ( XVal )
        IMPLICIT NONE
        REAL*8 XVal

        REAL*8 Temp

        INCLUDE 'astmath.cmn'

        Temp= DTAN( XVal )

        IF (DABS( Temp ) .lt. Small) THEN
            DCOT = Infinite
          ELSE
            DCOT = 1.0 / Temp
          ENDIF
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           FUNCTION ACOSH
*
*  this function evaluates the inverse hyperbolic cosine FUNCTION.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    XVal        - ANGLE Value                                  1.0D0 to Infinity
*
*  OutPuts       :
*    ACOSH     - Result                                       any real
*
*  Locals        :
*    Temp        - Temporary REAL*8 Value
*
*  Coupling      :
*    None.
*
* ------------------------------------------------------------------------------  

      REAL*8 FUNCTION ACOSH( XVal )
        IMPLICIT NONE
        REAL*8 XVal
* -----------------------------  Locals  ------------------------------
        REAL*8 Undefined, Temp

        ! --------------------  Implementation   ----------------------
        Undefined  = 999999.1D0
        IF ( XVal*XVal - 1.0D0 .lt. 0.0D0 ) THEN
            Temp= Undefined
            Write(*,*) 'Error in ACOSH FUNCTION '
          ELSE
            Temp= DLOG( XVal + DSQRT( XVal*XVal - 1.0D0 ) )
          ENDIF

        ACOSH= Temp 
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           FUNCTION SINH
*
*  this function evaluates the hyperbolic sine FUNCTION.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    XVal        - ANGLE Value                                  any real
*
*  OutPuts       :
*    SINH        - Result                                       any real
*
*  Locals        :
*    None.
*
*  Coupling      :
*    None.
*
* ------------------------------------------------------------------------------  

      REAL*8 FUNCTION SINH   ( XVal )
        IMPLICIT NONE
        REAL*8 XVal
        ! --------------------  Implementation   ----------------------
        SINH= 0.5D0*( EXP(XVal) - EXP(-XVal) ) 
      RETURN
      END


* ------------------------------------------------------------------------------
*
*                           FUNCTION ASINH
*
*  this function evaluates the inverse hyperbolic sine FUNCTION.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    XVal        - ANGLE Value                                  any real
*
*  OutPuts       :
*    ASINH     - Result                                       any real
*
*  Locals        :
*    None.
*
*  Coupling      :
*    None.
*
* ------------------------------------------------------------------------------  

      REAL*8 FUNCTION ASINH( XVal )
        IMPLICIT NONE
        REAL*8 XVal

        ! --------------------  Implementation   ----------------------
        ASINH= DLOG( XVal + DSQRT( XVal*XVal + 1.0D0 ) )

      RETURN
      END  ! end asinh

* ------------------------------------------------------------------------------
*
*                           FUNCTION ATANH
*
*  this function evaluates the inverse hyperbolic tangent FUNCTION.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    XVal        - ANGLE Value                                  -1.0D0 to 1.0D0
*
*  OutPuts       :
*    ATANH     - Result                                       any real
*
*  Locals        :
*    Temp        - Temporary REAL*8 Value
*
*  Coupling      :
*    None.
*
* ------------------------------------------------------------------------------  

      REAL*8 FUNCTION ATANH( XVal )
        IMPLICIT NONE
        REAL*8 XVal
* -----------------------------  Locals  ------------------------------
        REAL*8 Temp

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        IF ( 1.0D0 - DABS(XVal) .lt. Small ) THEN
            Temp= Undefined
            Write(*,*) 'Error in ATANH FUNCTION '
          ELSE
            Temp= 0.5D0 * DLOG( (1.0D0 + XVal) / (1.0D0 - XVal) )
          ENDIF

        ATANH= Temp 
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           FUNCTION DOT
*
*  this function finds the DOT product of two vectors.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Vec1        - Vector number 1
*    Vec2        - Vector number 2
*
*  OutPuts       :
*    DOT         - Result
*
*  Locals        :
*    None.
*
*  Coupling      :
*    None.
*
* ------------------------------------------------------------------------------  

      REAL*8 FUNCTION DOT    ( Vec1,Vec2 )
        IMPLICIT NONE
        REAL*8 Vec1(3), Vec2(3)

        ! --------------------  Implementation   ----------------------
        DOT= Vec1(1)*Vec2(1) + Vec1(2)*Vec2(2) + Vec1(3)*Vec2(3)
      RETURN
      END

      

      
* ------------------------------------------------------------------------------
*
*                           FUNCTION MAG
*
*  this subroutine finds the magnitude of a vector.  The tolerance is set to
*    0.00000001D0, thus the 1.0D0E-16 for the squared test of underflows.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Vec       - Vector
*
*  OutPuts       :
*    Mag       - Answer 
*
*  Locals        :
*    None.
*
*  Coupling      :
*    None.
*
* ------------------------------------------------------------------------------  

      REAL*8 FUNCTION MAG    ( Vec )
        IMPLICIT NONE
        REAL*8 Vec(3)
* -----------------------------  Locals  ------------------------------
        Real*8 Temp

        ! --------------------  Implementation   ----------------------
        Temp= Vec(1)*Vec(1) + Vec(2)*Vec(2) + Vec(3)*Vec(3)

        IF ( DABS( Temp ) .ge. 1.0D-16 ) THEN
            MAG = DSQRT( Temp )
          ELSE
            MAG = 0.0D0
          ENDIF
      RETURN
      END

            

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE CROSS
*
*  this subroutine crosses two vectors.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Vec1        - Vector number 1
*    Vec2        - Vector number 2
*
*  OutPuts       :
*    OutVec      - Vector result of A x B
*
*  Locals        :
*    None.
*
*  Coupling      :
*
* ------------------------------------------------------------------------------  

      SUBROUTINE CROSS       ( Vec1,Vec2, OutVec )
        IMPLICIT NONE
        REAL*8 Vec1(3), Vec2(3), OutVec(3)

        ! --------------------  Implementation   ----------------------
        OutVec(1)= Vec1(2)*Vec2(3) - Vec1(3)*Vec2(2)
        OutVec(2)= Vec1(3)*Vec2(1) - Vec1(1)*Vec2(3)
        OutVec(3)= Vec1(1)*Vec2(2) - Vec1(2)*Vec2(1)

      RETURN
      END


* ------------------------------------------------------------------------------
*
*                           SUBROUTINE NORM
*
*  this subroutine calculates a unit vector given the original vector.  If  a
*    zero vector is input, the vector is set to zero.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Vec         - Vector
*
*  OutPuts       :
*    OutVec      - Unit Vector
*
*  Locals        :
*    MagVec      - Magnitude of vector
*    i           - Index
*
*  Coupling      :
*    MAG           Magnitude of a vector
*
* ------------------------------------------------------------------------------  

      SUBROUTINE NORM        ( Vec, OutVec )
        IMPLICIT NONE
        REAL*8 Vec(3), OutVec(3)
        EXTERNAL MAG

* -----------------------------  Locals  ------------------------------
        REAL*8 MAG, MagVec
        INTEGER i

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        MagVec = MAG( Vec )
        IF ( MagVec .gt. Small ) THEN
            DO i= 1, 3
                OutVec(i)= Vec(i)/MagVec
              ENDDO
          ELSE
            DO i= 1, 3
                OutVec(i)= 0.0D0
              ENDDO
          ENDIF

      RETURN
      END
      

      
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE ANGLE
*
*  this subroutine calculates the ANGLE between two vectors.  The output is
*    set to 999999.1D0 to indicate an undefined value.  Be SURE to check  
*    this at the output phase.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Vec1        - Vector number 1
*    Vec2        - Vector number 2
*
*  OutPuts       :
*    Theta       - ANGLE between the two vectors  -Pi to Pi
*
*  Locals        :
*    Temp        - Temporary REAL variable
*
*  Coupling      :
*    DOT           DOT Product of two vectors
*    DACOS        Arc Cosine FUNCTION
*
* ------------------------------------------------------------------------------  

      SUBROUTINE ANGLE       ( Vec1,Vec2, Theta )
        IMPLICIT NONE
        REAL*8 Vec1(3), Vec2(3), Theta, magvec1, magvec2
        EXTERNAL Dot, Mag
* -----------------------------  Locals  ------------------------------
        REAL*8 Temp, Dot, Mag
        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        magvec1 = MAG(vec1)
        magvec2 = MAG(vec2)
        IF ( magVec1*magVec2 .gt. Small**2 ) THEN
            Temp= DOT(Vec1,Vec2) / (magVec1*magVec2)
            IF ( DABS( Temp ) .gt. 1.0D0 ) THEN
                Temp= DSIGN(1.0D0, Temp)
              ENDIF
            Theta= DACOS( Temp ) 
          ELSE
            Theta= Undefined
          ENDIF
      RETURN
      END  ! end angle

            

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE ROTi
*
*  this subroutine performs a rotation about the ith axis. i is specified
*    for each operation.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Vec         - Input vector
*    XVal        - ANGLE of rotation              rad
*
*  OutPuts       :
*    OutVec      - Vector result
*
*  Locals        :
*    c           - Cosine of the ANGLE XVal
*    s           - Sine of the ANGLE XVal
*    Temp        - Temporary REAL*8 value
*
*  Coupling      :
*    None.
*
* ------------------------------------------------------------------------------

      SUBROUTINE ROT1        ( Vec, XVal, OutVec )
        IMPLICIT NONE
        REAL*8 Vec(3), XVal, OutVec(3)
* -----------------------------  Locals  ------------------------------
        REAL*8 c, s, Temp

        Temp= Vec(3)
        c= DCOS( XVal )
        s= DSIN( XVal )

        OutVec(3)= c*Vec(3) - s*Vec(2)
        OutVec(2)= c*Vec(2) + s*Temp
        OutVec(1)= Vec(1)
      RETURN
      END


      SUBROUTINE ROT2        ( Vec, XVal, OutVec )
        IMPLICIT NONE
        REAL*8 Vec(3), XVal, OutVec(3)
* -----------------------------  Locals  ------------------------------
        REAL*8 c, s, Temp

        ! --------------------  Implementation   ----------------------
        Temp= Vec(3)
        c= DCOS( XVal )
        s= DSIN( XVal )

        OutVec(3)= c*Vec(3) + s*Vec(1)
        OutVec(1)= c*Vec(1) - s*Temp
        OutVec(2)= Vec(2)
      RETURN
      END


      SUBROUTINE ROT3        ( Vec, XVal, OutVec )
        IMPLICIT NONE
        REAL*8 Vec(3), XVal, OutVec(3)
* -----------------------------  Locals  ------------------------------
        REAL*8 c, s, Temp

        ! --------------------  Implementation   ----------------------
        Temp= Vec(2)
        c= DCOS( XVal )
        s= DSIN( XVal )

        OutVec(2)= c*Vec(2) - s*Vec(1)
        OutVec(1)= c*Vec(1) + s*Temp
        OutVec(3)= Vec(3)
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                                  rot1mat
*
*  this function sets up a rotation matrix for an input angle about the first
*    axis.
*
*  author        : david vallado                  719-573-2600   10 jan 2003
*
*  revisions
*                -
*
*  inputs          description                    range / units
*    xval        - angle of rotation              rad
*
*  outputs       :
*    outmat      - matrix result
*
*  locals        :
*    c           - cosine of the angle xval
*    s           - sine of the angle xval
*
*  coupling      :
*    none.
*
* ----------------------------------------------------------------------------- }

      SUBROUTINE ROT1MAT     ( XVal, OutMat )
        IMPLICIT NONE
        REAL*8 XVal, OutMat(3,3)

* -----------------------------  Locals  ------------------------------
        REAL*8 c, s

        c= Dcos( xval )
        s= Dsin( xval )

        outmat(1,1)= 1.0D0
        outmat(1,2)= 0.0D0
        outmat(1,3)= 0.0D0

        outmat(2,1)= 0.0D0
        outmat(2,2)= c
        outmat(2,3)= s

        outmat(3,1)= 0.0D0
        outmat(3,2)=  c
        outmat(3,3)= -s

      RETURN
      END

      SUBROUTINE ROT2MAT     ( XVal, OutMat )
        IMPLICIT NONE
        REAL*8 XVal, OutMat(3,3)

* -----------------------------  Locals  ------------------------------
        REAL*8 c, s

        c= Dcos( xval )
        s= Dsin( xval )

        outmat(1,1)= c
        outmat(1,2)= 0.0D0
        outmat(1,3)= -s

        outmat(2,1)= 0.0D0
        outmat(2,2)= 1.0D0
        outmat(2,3)= 0.0D0

        outmat(3,1)= c
        outmat(3,2)= 0.0D0
        outmat(3,3)= s

      RETURN
      END

      SUBROUTINE ROT3MAT     ( XVal, OutMat )
        IMPLICIT NONE
        REAL*8 XVal, OutMat(3,3)

* -----------------------------  Locals  ------------------------------
        REAL*8 c, s

        c= Dcos( xval )
        s= Dsin( xval )

        outmat(1,1)= c
        outmat(1,2)= s
        outmat(1,3)= 0.0D0

        outmat(2,1)= c
        outmat(2,2)= -s
        outmat(2,3)= 0.0D0

        outmat(3,1)= 0.0D0
        outmat(3,2)= 0.0D0
        outmat(3,3)= 1.0D0

      RETURN
      END
      
            

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE ADDVEC
*
*  this subroutine adds two vectors.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Vec1        - Vector number 1
*    Vec2        - Vector number 2
*
*  OutPuts       :
*    OutVec      - Vector result of A + B
*
*  Locals        :
*    i           - Index
*
*  Coupling      :
*
*
* ------------------------------------------------------------------------------  

      SUBROUTINE ADDVEC      ( Vec1,Vec2, OutVec )
        IMPLICIT NONE
        REAL*8 Vec1(3), Vec2(3), OutVec(3)
* -----------------------------  Locals  ------------------------------
        INTEGER i

        ! --------------------  Implementation   ----------------------
        DO i=1 , 3
            OutVec(i)= Vec1(i) + Vec2(i)
          ENDDO
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE SUBVEC
*
*  this subroutine subtracts two vectors.
*
*  Author        : David Vallado                  719-573-2600   13 Feb 2003
*
*  Inputs          Description                    Range / Units
*    Vec1        - Vector number 1
*    Vec2        - Vector number 2
*
*  OutPuts       :
*    OutVec      - Vector result of A - B
*
*  Locals        :
*    i           - Index
*
*  Coupling      :
*
*
* ------------------------------------------------------------------------------  

      SUBROUTINE SUBVEC      ( Vec1,Vec2, OutVec )
        IMPLICIT NONE
        REAL*8 Vec1(3), Vec2(3), OutVec(3)
* -----------------------------  Locals  ------------------------------
        INTEGER i

        ! --------------------  Implementation   ----------------------
        DO i=1 , 3
            OutVec(i)= Vec1(i) - Vec2(i)
          ENDDO
      RETURN
      END


* ------------------------------------------------------------------------------
*
*                           SUBROUTINE ADD3VEC
*
*  this subroutine adds three vectors.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Vec1        - Vector number 1
*    Vec2        - Vector number 2
*    Vec3        - Vector number 3
*
*  OutPuts       :
*    OutVec      - Vector result of Vec1 + Vec2 + Vec3
*
*  Locals        :
*    i           - Index
*
*  Coupling      :
*
*
* ------------------------------------------------------------------------------

      SUBROUTINE ADD3VEC      ( Vec1,Vec2,Vec3, OutVec )
        IMPLICIT NONE
        REAL*8 Vec1(3), Vec2(3), Vec3(3), OutVec(3)
* -----------------------------  Locals  ------------------------------
        INTEGER i

        ! --------------------  Implementation   ----------------------
        DO i = 1, 3
            OutVec(i)= Vec1(i) + Vec2(i) + Vec3(i)
          ENDDO
      RETURN
      END
      

* ------------------------------------------------------------------------------
      
      SUBROUTINE MatVecMult (Mat1, Vec2, Mat1r,Mat1c,Max1r,Max1c, Vec3 )
        IMPLICIT NONE
        INTEGER Mat1r,Mat1c,Max1r,Max1c
        REAL*8 Mat1(Max1r,Max1c), Vec2(Max1r), Vec3(Max1r)
        ! -- local vars ---
        INTEGER Row, ktr

        DO Row=1, mat1r
            Vec3(Row) = 0.0D0
            DO ktr= 1, mat1c
                Vec3(Row) = Vec3(Row) + Mat1(Row,ktr) * Vec2(ktr)
              ENDDO
          ENDDO
      RETURN
      END   ! SUBROUTINE MatVecMult
      
      
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE LNCOM1
*
*  this subroutine calculates the linear combination of a vector
*    multiplied by a constants.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    A1          - constant number
*    Vec         - Vector number
*
*  OutPuts       :
*    OutVec      - Vector result of A1*Vec1 + A2*Vec2
*
*  Locals        :
*    i           - Index
*
*  Coupling      :
*
*
* ------------------------------------------------------------------------------  

      SUBROUTINE LNCOM1      ( A, Vec,  OutVec )
        IMPLICIT NONE
        REAL*8 A, Vec(3), OutVec(3)
* -----------------------------  Locals  ------------------------------
        INTEGER i

        ! --------------------  Implementation   ----------------------
        DO i= 1, 3
            OutVec(i)= A*Vec(i)
          ENDDO
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE LNCOM2
*
*  this subroutine calculates the linear combination of two vectors
*    multiplied by two different constants.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    A1          - constant number 1
*    A2          - constant number 2
*    Vec1        - Vector number 1
*    Vec2        - Vector number 2
*
*  OutPuts       :
*    OutVec      - Vector result of A1*Vec1 + A2*Vec2
*
*  Locals        :
*    i           - Index
*
*  Coupling      :
*
*
* ------------------------------------------------------------------------------  

      SUBROUTINE LNCOM2      ( A1, A2, Vec1, Vec2,  OutVec )
        IMPLICIT NONE
        REAL*8 A1, A2, Vec1(3), Vec2(3), OutVec(3)
* -----------------------------  Locals  ------------------------------
        INTEGER i

        ! --------------------  Implementation   ----------------------
        DO i= 1, 3
            OutVec(i)= a1*Vec1(i) + a2*Vec2(i)
          ENDDO
      RETURN
      END

*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE LNCOM3
*
*  this subroutine calculates the linear combination of three vectors
*    multiplied by three different constants.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    A1          - constant number 1
*    A2          - constant number 2
*    A3          - constant number 3
*    Vec1        - Vector number 1
*    Vec2        - Vector number 2
*    Vec3        - Vector number 3
*
*  OutPuts       :
*    OutVec      - Vector result of A1*Vec1 + A2*Vec2 + A3*Vec3
*
*  Locals        :
*    i           - Index
*
*  Coupling      :
*
*
* ------------------------------------------------------------------------------  

      SUBROUTINE LNCOM3      ( A1, A2, A3, Vec1, Vec2, Vec3,  OutVec )
        IMPLICIT NONE
        REAL*8 A1, A2, A3, Vec1(3), Vec2(3), Vec3(3), OutVec(3)
* -----------------------------  Locals  ------------------------------
        INTEGER i

        ! --------------------  Implementation   ----------------------
        DO i= 1, 3
            OutVec(i)= a1*Vec1(i) + a2*Vec2(i) + a3*Vec3(i)
          ENDDO
      RETURN
      END
      
     

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE MATSCALE
*
*  this subroutine multiplies a matrix by a constant.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Mat1        - Matrix number 1
*    Scale       - Constant value to multiply by
*
*  OutPuts       :
*    Mat2        - Matrix result of Mat1 * Scale
*
*  Locals        :
*    Row         - Row Index
*    Col         - Column Index
*
*  Coupling      :
*
* ------------------------------------------------------------------------------

      SUBROUTINE MATSCALE   ( Mat1,Scale,Mat1r,Mat1c,Max1r,Max1c, Mat2 )
        IMPLICIT NONE
        INTEGER Mat1r,Mat1c, Max1r, Max1c
        REAL*8 Mat1(Max1r,Max1c),Scale,Mat2(Max1r,Max1c)
* -----------------------------  Locals  ------------------------------
        INTEGER Row,Col

        ! --------------------  Implementation   ----------------------
        DO Row=1, Mat1r
            DO Col= 1, Mat1c
                  Mat2(Row,Col)= Mat1(Row,Col)*Scale
              ENDDO
          ENDDO
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE MATADD
*
*  this subroutine adds two matricies together.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Mat1        - Matrix number 1
*    Mat2        - Matrix number 2
*    Mat1r       - Matrix number 1 rows
*    Mat1c       - Matrix number 1 columns
*
*  OutPuts       :
*    Mat3        - Matrix result of Mat1 + Mat2
*                    of size mat1r x mat1c
*
*  Locals        :
*    Row         - Row Index
*    Col         - Column Index
*
*  Coupling      :
*
* ------------------------------------------------------------------------------  

      SUBROUTINE MATADD    ( Mat1,Mat2, Mat1r,Mat1c,Max1r,Max1c, Mat3 )
        IMPLICIT NONE
        INTEGER Mat1r,Mat1c,Max1r,Max1c
        REAL*8 Mat1(Max1r,Max1c),Mat2(Max1r,Max1c),Mat3(Max1r,Max1c)
* -----------------------------  Locals  ------------------------------
        INTEGER Row,Col

        ! --------------------  Implementation   ----------------------
        DO Row = 1, Mat1r

            DO Col =1, Mat1c

                  Mat3(Row,Col)= Mat1(Row,Col) + Mat2(Row,Col)
              ENDDO

          ENDDO

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE MATSUB
*
*  this subroutine subtracts two matricies together.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Mat1        - Matrix number 1
*    Mat2        - Matrix number 2
*    Mat1r       - Matrix number 1 rows
*    Mat1c       - Matrix number 1 columns
*
*  OutPuts       :
*    Mat3        - Matrix result of Mat1 + Mat2
*                    of size mat1r x mat1c
*
*  Locals        :
*    Row         - Row Index
*    Col         - Column Index
*
*  Coupling      :
*
* ------------------------------------------------------------------------------

      SUBROUTINE MATSUB    ( Mat1,Mat2, Mat1r,Mat1c,Max1r,Max1c, Mat3 )
        IMPLICIT NONE
        INTEGER Mat1r,Mat1c,Max1r,Max1c
        REAL*8 Mat1(Max1r,Max1c),Mat2(Max1r,Max1c),Mat3(Max1r,Max1c)
* -----------------------------  Locals  ------------------------------
        INTEGER Row,Col

        ! --------------------  Implementation   ----------------------
        DO Row= 1, Mat1r
            DO Col= 1, Mat1c
                  Mat3(Row,Col)= Mat1(Row,Col) - Mat2(Row,Col)
              ENDDO

          ENDDO

      RETURN
      END
      
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE MATMULT
*
*  this subroutine multiplies two matricies together.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Mat1        - Matrix number 1
*    Mat2        - Matrix number 2
*    Mat1r       - Matrix number 1 rows
*    Mat1c       - Matrix number 1 columns
*    Mat2c       - Matrix number 2 columns
*
*  OutPuts       :
*    Mat3        - Matrix result of Mat1 * Mat2 of size mat1r x mat2c
*
*  Locals        :
*    Row         - Row Index
*    Col         - Column Index
*    ktr         - Index
*
*  Coupling      :
*
* ------------------------------------------------------------------------------  

      SUBROUTINE MATMULT     ( Mat1,Mat2, Mat1r,Mat1c,Mat2c,
     &                         Max1r,Max1c,Max2c, Mat3 )
        IMPLICIT NONE
        INTEGER Mat1r,Mat1c,Mat2c, Max1r,Max1c,Max2c
        REAL*8 Mat1(Max1r,Max1c),Mat2(Max1c,Max2c),Mat3(Max1r,Max2c)
* -----------------------------  Locals  ------------------------------
        INTEGER Row,Col,ktr

        ! --------------------  Implementation   ----------------------
        DO Row=1, Mat1r
            DO Col= 1, Mat2c
                Mat3(Row,Col) = 0.0D0
                DO ktr= 1, Mat1c
                    Mat3(Row,Col)= Mat3(Row,Col)+Mat1(Row,ktr)*
     &                                    Mat2(ktr,Col)
                  ENDDO
              ENDDO

          ENDDO
      RETURN
      END



* ------------------------------------------------------------------------------
*
*                           SUBROUTINE MATTRANS
*
*  this subroutine finds the transpose of a matrix.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Mat1        - Matrix number 1
*    Mat1r       - Matrix number 1 rows
*    Mat1c       - Matrix number 1 columns
*
*  OutPuts       :
*    Mat2        - Matrix result of transpose Mat2
*
*  Locals        :
*    Row         - Row Index
*    Col         - Column Index
*
*  Coupling      :
*
* ------------------------------------------------------------------------------  

      SUBROUTINE MATTRANS      ( Mat1, Mat1r,Mat1c,Max1r,Max1c, Mat2 )
        IMPLICIT NONE
        INTEGER Mat1r,Mat1c,Max1r,Max1c
        REAL*8 Mat1(Max1r,Max1c),Mat2(Max1c,Max1r)
* -----------------------------  Locals  ------------------------------
        INTEGER Row,Col

        ! --------------------  Implementation   ----------------------
        DO Row =1, Mat1r
            DO Col = 1, Mat1c
                  Mat2(Col,Row)= Mat1(Row,Col)
              ENDDO
          ENDDO
      RETURN
      END
*


* ------------------------------------------------------------------------------
*
*                           SUBROUTINE MAKEMAT
*
*  this subroutine forms a rotation matrix for a given axis of rotation.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Angl        - ANGLE through which to rotate
*    Numbr       - Number of axis for rotation
*
*  OutPuts       :
*    Matr        - Matrix of rotation
*
*  Locals        :
*    None.
*
*  Coupling      :
*
* ------------------------------------------------------------------------------  

      SUBROUTINE MAKEMAT     ( Angl, Numbr, Matr )
        IMPLICIT NONE
        REAL*8 Angl, Matr(3,3)
        INTEGER Numbr

        ! --------------------  Implementation   ----------------------
        IF (Numbr .eq.1 ) THEN
            Matr(1,1)=  1.0D0
            Matr(2,2)=  DCOS(Angl)
            Matr(2,3)=  DSIN(Angl)
            Matr(3,2)= -DSIN(Angl)
            Matr(3,3)=  DCOS(Angl)
          ENDIF
        IF (Numbr .eq.2 ) THEN
            Matr(2,2)=  1.0D0
            Matr(1,1)=  DCOS(Angl)
            Matr(1,3)= -DSIN(Angl)
            Matr(3,1)=  DSIN(Angl)
            Matr(3,3)=  DCOS(Angl)
          ENDIF
        IF (Numbr .eq.3 ) THEN
            Matr(1,1)=  DCOS(Angl)
            Matr(1,2)=  DSIN(Angl)
            Matr(2,1)= -DSIN(Angl)
            Matr(2,2)=  DCOS(Angl)
            Matr(3,3)=  1.0D0
          ENDIF
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE LUDECOMP
*
*  this subroutine decomposes a matrix into an LU form.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*                  Dave Vallado                     719-554-3638    1 Aug 1989
*
*  Inputs          Description                    Range / Units
*    Order       - Order of matrix
*
*  OutPuts       :
*    LU          - LU decomposition matrix
*    Index       - Index vector for pivoting
*
*  Locals        :
*    i           - Index
*    j           - Index
*    k           - Index
*    IMax        - Pivot row pointer
*    Scale       - Scale FACTOR vector
*    Sum         - Temporary Variables
*    AMax        - Temporary Variables
*    Dum         - Temporary Variables
*
*  Coupling      :
*
*  References    :
*    Numerical Recipes - Flannery
*
* ------------------------------------------------------------------------------
*
      SUBROUTINE LUDeComp    ( LU, Index, Order, MaxRow )
        IMPLICIT NONE
        INTEGER Order, Index(Order), MaxRow
        REAL*8 LU(MaxRow,MaxRow)
* -----------------------------  Locals  ------------------------------
        REAL*8 Small, Scale(30), Sum, AMax, Dum ! Hardwire Scale to 30
        INTEGER I, J, K, IMax

        ! --------------------  Implementation   ----------------------
        Small= 0.000001D0
        IMax = 0 
        DO I = 1, Order
            AMax = 0.0D0
            DO J = 1, Order
                IF ( (DABS(  LU(i,j)) .gt. AMax) ) THEN
                    AMax = DABS( LU(i,j))
                  ENDIF
              ENDDO
            IF ( DABS(AMax) .lt. Small ) THEN
                Write(*,*) ' Singular Matrix AMax '
              ENDIF 
            Scale(i)= 1.0D0 / AMax
          ENDDO

        DO j = 1, Order
            DO i = 1, j - 1
                Sum =   LU(i,j)
                DO k = 1, i - 1
                    Sum = Sum -   LU(i,k)*  LU(k,j)
                  ENDDO
                LU(i,j)= Sum
              ENDDO
            AMax = 0.0D0 
            DO i = j, Order
                Sum =   LU(i,j)
                DO k = 1, j - 1
                    Sum = Sum -   LU(i,k)*  LU(k,j)
                  ENDDO
                LU(i,j)=  Sum
                Dum =   Scale(i )*DABS(Sum)
                IF ( (Dum .ge. AMax) ) THEN
                    IMax = i
                    AMax = Dum 
                  ENDIF 
              ENDDO
            IF ( (j .ne. iMax) ) THEN
                DO k = 1, Order
                    Dum =   LU(iMax,K)
                    LU(iMax,k)=   LU(j,k )
                    LU(j,k)=   Dum
                  ENDDO
                Scale(iMax)=   Scale(j)
              ENDIF 
            Index(j) = IMax
            IF ( DABS(  LU(j,j)) .lt. Small ) THEN
                Write(*,*) ' Matrix is Singular LU '
              ENDIF 
            IF ( (j .ne. Order) ) THEN
                Dum = 1.0D0 /   LU(j,j)
                DO i = j + 1, Order
                    LU(i,j)=  Dum*  LU(i,j)
                  ENDDO
              ENDIF 
          ENDDO
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE LUBkSUB
*
*  this subroutine finds the inverse of a matrix using LU decomposition.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Order       - Order of matrix
*    LU          - LU decomposition matrix
*    Index       - Index vector for pivoting
*
*  OutPuts       :
*    B           - Solution Vector
*
*  Locals        :
*    i           - Index
*    j           - Index
*    I0          - Pointer to non-zero element
*    IPtr        - Pivot Rwo Pointer
*    Sum         - Temporary Variables
*
*  Coupling      :
*
*  References    :
*    Numerical Recipes - Flannery
*
* ------------------------------------------------------------------------------  

      SUBROUTINE LUBkSub     ( LU, Index, Order, B, MaxRow )
        IMPLICIT NONE
        INTEGER Order, MaxRow, Index(Order)
        REAL*8 LU(MaxRow,MaxRow), B(MaxRow)
* -----------------------------  Locals  ------------------------------
        INTEGER I, J, IPtr, I0
        REAL*8 Sum

        ! --------------------  Implementation   ----------------------
        I0 = 0 
           Write(*,'(A)') ' index'
            Write(*,'(3(I4,1X))') Index(1),Index(2),Index(3)
        DO i = 1, Order
            IPtr = Index(i)
            Sum =   B(IPtr)
            B(Iptr)=   B(i)
            IF ( (I0 .ne. 0) ) THEN
                DO j = I0, i - 1
                    Sum = Sum -   LU(i,j)*  B(j)
                  ENDDO
              ELSE
                IF ( (Sum .ne. 0.0D0) ) THEN
                   I0 = I
                 ENDIF
              ENDIF
            B(i)=  Sum
          ENDDO

        B(Order)=   B(Order)/  LU(Order,Order)
           Write(*,'(A)') '  in lubks b'
           Write(*,'(3(F12.7,1X))') B(1),B(2),B(3)
        DO i = (Order - 1), 1, -1
            Sum =   B(i)
            DO j = i + 1, Order
                Sum = Sum -   LU(i,j)*  B(j)
              ENDDO
            B(i)= Sum /   LU(i,i)
          ENDDO
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE MATINVERSE
*
*  this subroutine finds the inverse of a matrix using LU decomposition.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Mat         - Matrix to invert
*    Order       - Order of matrix
*
*  OutPuts       :
*    MatInv      - Inverted matrix
*
*  Locals        :
*    i           - Index
*    j           - Index
*    Index       - Index vector for pivoting
*    LU          - LU decomposition matrix
*    B           - Operational vector to form MatInv
*
*  Coupling      :
*    LUDeComp    -
*    LUBkSub     -
*
*  References    :
*    Numerical Recipes - Flannery
*
* ------------------------------------------------------------------------------  

      SUBROUTINE MATINVERSE  ( Mat, Order,MaxRow, MatInv )
        IMPLICIT NONE
        INTEGER Order,MaxRow
        REAL*8 Mat(MaxRow,MaxRow), MatInv(MaxRow,MaxRow)
* -----------------------------  Locals  ------------------------------
        INTEGER MaxR
        PARAMETER (MaxR = 10)
        INTEGER I, J, Index( MaxR )
        REAL*8 LU(MaxR,MaxR),B(MaxR)

        ! --------------------  Implementation   ----------------------
         DO i = 1, Order
             Index(i) = i
             DO j = 1, Order
                 LU(i,j)=   Mat(i,j)
               ENDDO
           ENDDO 

           Write(*,'(A)') ' lu'
             Write(*,'(3(F12.7,1X))') LU(1,1),LU(1,2),LU(1,3)
             Write(*,'(3(F12.7,1X))') LU(2,1),LU(2,2),LU(2,3)
             Write(*,'(3(F12.7,1X))') LU(3,1),LU(3,2),LU(3,3)

          CALL LUDeComp(LU, Index, Order, MaxR )
           Write(*,'(A)') ' lu'
             Write(*,'(3(F12.7,1X))') LU(1,1),LU(1,2),LU(1,3)
             Write(*,'(3(F12.7,1X))') LU(2,1),LU(2,2),LU(2,3)
             Write(*,'(3(F12.7,1X))') LU(3,1),LU(3,2),LU(3,3)
             Write(*,*)
         DO j = 1, Order
             DO i = 1, Order
                IF ( (i .eq. j) ) THEN
                    B(i)= 1.0D0
                  ELSE
                    B(i)= 0.0D0
                  ENDIF

               ENDDO

           Write(*,'(A)') ' b'
           Write(*,'(3(F12.7,1X))') B(1),B(2),B(3)
           CALL LUBkSub(LU, Index, Order, B, MaxR)
           Write(*,'(A)') ' b after lubsub in inv'
           Write(*,'(3(F12.7,1X))') B(1),B(2),B(3)
           Write(*,'(A)') ' lu'
             Write(*,'(3(F12.7,1X))') LU(1,1),LU(1,2),LU(1,3)
             Write(*,'(3(F12.7,1X))') LU(2,1),LU(2,2),LU(2,3)
             Write(*,'(3(F12.7,1X))') LU(3,1),LU(3,2),LU(3,3)


             DO i = 1, Order
                 MatInv(i,j)=    B(i)
               ENDDO
           ENDDO
           Write(*,'(A)') ' lu'
       Write(*,'(3(F12.7,1X))') MatInv(1,1),MatInv(1,2),MatInv(1,3)
       Write(*,'(3(F12.7,1X))') MatInv(2,1),MatInv(2,2),MatInv(2,3)
       Write(*,'(3(F12.7,1X))') MatInv(3,1),MatInv(3,2),MatInv(3,3)

      RETURN
      END

      
* ------------------------------------------------------------------------------
*
*                           FUNCTION DETERMINANT
*
*  this function calculates the determinant value using L-U decompisition.
*  The formula must have a NON-ZERO number in the 1,1 position.  If  the
*  FUNCTION senses a NON-ZERO number in row 1, it exchanges row1 for a row
*  WITH a NON-ZERO number.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*                                                                16 Jun 1993
*  Inputs          Description                    Range / Units
*    Order       - Order of determinaent (# of rows)
*    Mat1        - Matrix to find determinant of
*
*  OutPuts       :
*    Determinant - Result
*
*  Locals        :
*    i           - Index
*    j           - Index
*    k           - Index
*    n           - Index
*    Temp        -
*    D           -
*    Sum         -
*    L           -
*    U           -
*    Small       - Tolarance for comparing to 0.0D0
*
*  Coupling      :
*
*  References    :
*    Marion        pg. 168-172, 126-127
*
* ------------------------------------------------------------------------------  

      REAL*8 FUNCTION DETERMINANT ( Mat1, Order,MaxRow )
        IMPLICIT NONE
        INTEGER Order,MaxRow
        REAL*8 Mat1( MaxRow,MaxRow )
* -----------------------------  Locals  ------------------------------
        REAL*8 Small, Temp, D, Sum
        INTEGER i,j,k,MaxR
        PARAMETER (MaxR = 10)
        REAL*8 L(MaxR,MaxR), U(MaxR,MaxR)

        ! --------------------  Implementation   ----------------------
        Small = 0.00000001D0
        Sum   = 0.0D0
        ! ---------- Switch a non zero row to the first row -----------
        IF ( DABS(   Mat1(1,1 ) ) .lt. Small ) THEN
            j= 1
            DO WHILE (j .le. Order)
                IF ( DABS(   Mat1(j,1 ) ) .gt. Small ) THEN
                    DO k= 1, Order
                        Temp=    Mat1(1,k )
                          Mat1(1,k)=   Mat1(j,k )
                          Mat1(j,k)= Temp
                      ENDDO
                    j= Order + 1 
                  ENDIF 
                j= j+1 
              ENDDO
          ENDIF  ! If  DABS(Mat1(1,1)) .lt. Small

        DO i= 1, Order
              L(i,1)=  Mat1(i,1 )
          ENDDO
        DO j= 1, Order
              U(1,j)=  Mat1(1,j ) /   L(1,1 )
          ENDDO
        DO j= 2, Order
            DO i= j, Order
                Sum= 0.0D0
                DO k= 1, j-1
                    Sum= Sum+   L(i,k )*   U(k,j )
                  ENDDO
                L(i,j)=  Mat1(i,j )- Sum
              ENDDO   ! DO i
              U(j,j)=  1.0D0
            IF ( j .ne. Order ) THEN
                DO i= j+1, Order
                    Sum= 0.0D0
                    DO k= 1, j-1
                        Sum= Sum +    L(j,k)*   U(k,i )
                      ENDDO
                      U(j,i)=  (Mat1(j,i)-Sum)/   L(j,j )
                  ENDDO   ! DO i
              ENDIF   ! If  j
          ENDDO   ! DO j
        D= 1.0D0 
        DO i= 1, Order
            D= D*   L(i,i)
          ENDDO
        Determinant= D 

      RETURN
      END
      
      
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE WRITEMAT
*
*  this subroutine WRITEs a matrix.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Mat1        - Matrix to WRITE out
*    DecNum      - Number of decimals to write
*
*  OutPuts       :
*    None.
*
*  Locals        :
*    Row         - Row Index
*    Col         - Column Index
*
*  Coupling      :
*    GETVAL      - Gets a value from the matrix
*
* ------------------------------------------------------------------------------

      SUBROUTINE WRITEMAT    ( Mat1, Mat1r,Mat1c,Max1r,Max1c,Title )
        IMPLICIT NONE
        INTEGER Mat1r,Mat1c,Max1r,Max1c
        REAL*8 Mat1( Max1r,Max1c )
        CHARACTER*64 Title
* -----------------------------  Locals  ------------------------------
        INTEGER  Row,Col

        ! --------------------  Implementation   ----------------------
        Write(*,'(A)') Title
        DO Row= 1, Mat1r
            Write(*,20) (Mat1(Row,Col),Col= 1, Mat1c)
          ENDDO
   20   FORMAT( 3(F12.7,1X) )
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE READMAT
*
*  this subroutine reads a matrix.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Mat1        - Matrix to WRITE out
*    Mat1r       - Matrix number 1 rows
*    Mat1c       - Matrix number 1 columns
*
*  OutPuts       :
*    None.
*
*  Locals        :
*    Row         - Row Index
*    Col         - Column Index
*
*  Coupling      :
*    None.
*
* ------------------------------------------------------------------------------  

      SUBROUTINE READMAT    ( Mat1, Mat1r,Mat1c,Max1r,Max1c )
        IMPLICIT NONE
        INTEGER Mat1r,Mat1c,Max1r,Max1c
        REAL*8 Mat1( Max1r,Max1c )
* -----------------------------  Locals  ------------------------------
        INTEGER  Row,Col

        ! --------------------  Implementation   ----------------------
        DO Row= 1, Mat1r
            Read(10,*) (Mat1(Row,Col),Col= 1, Mat1c)
          ENDDO
      RETURN
      END
*
      SUBROUTINE FILEWRITEMAT ( Mat1, Mat1r,Mat1c,Max1r,Max1c,Title )
        IMPLICIT NONE
        INTEGER Mat1r,Mat1c,Max1r,Max1c
        REAL*8 Mat1( Max1r,Max1c )
        CHARACTER*64 Title
* -----------------------------  Locals  ------------------------------
        INTEGER row,col

        ! --------------------  Implementation   ----------------------
        Write(20,'(A)') Title
        DO Row= 1, Mat1r
            Write(20,*) (Mat1(Row,Col),Col= 1, Mat1c)
          ENDDO
      RETURN
      END


* ----- WRITE a matrix using exponential formats -------  

      SUBROUTINE FILEEXPWRITEMAT ( Mat1,Mat1r,Mat1c,Max1r,Max1c,Title )
        IMPLICIT NONE
        INTEGER Mat1r, Mat1c,Max1r,Max1c
        REAL*8 Mat1(Max1r,Max1c)
        CHARACTER*64 Title
* -----------------------------  Locals  ------------------------------
        INTEGER  row,col

        ! --------------------  Implementation   ----------------------
        Write(20,'(A)') Title
        DO Row= 1, Mat1r
            Write(20,*) (' ',Mat1(Row,Col),Col= 1,Mat1c)
          ENDDO
      RETURN
      END





*
* ---- this subroutine takes an arbitrary number of points .and.
*  fits a polynomial of degree <degree>. The datapoints are stored in a
*  matrix having rows <NumPts> by 2 <one for x, and y>.
*  See Matthews pg 230-232.
*
*    Assumes the form of
*       ax2 + bx + c = y
*    coeff indices are
*       1 = a, 2 = b, 3 = c
*    MIN is at
*       2ax + b = 0
*    y is found by plugging a back into the polynomial
*
*  dav  3 Apr 97
* ---------------------------------------------------------------------

      SUBROUTINE POLYFIT     ( Degree,NumPts,DataPoints,Coeff,
     &                         MaxRow,MaxDeg,MinX,MinY )
        IMPLICIT NONE
        INTEGER Degree, NumPts,MaxRow,MaxDeg
        REAL*8 DataPoints(MaxRow,2),Coeff(MaxDeg,1), MinX, MinY
* -----------------------------  Locals  ------------------------------
        INTEGER j,k,r,c,MaxR
        PARAMETER (MaxR = 20)
        Real*8 p, ai(MaxR,MaxR), a(MaxR,MaxR),b(MaxR,1), parr(MaxR,1)

        ! --------------------  Implementation   ----------------------
        ! ----------- Find the sum of the product terms (x*y) ---------
        DO k= 1, NumPts
            p= 1
            DO r= 1, Degree+1
                b(r,1)=  b(r,1) + DataPoints(k,2)*p
                p= p*DataPoints(k,1)
              ENDDO
          ENDDO

        ! ------------ Find the sum of powers for x (x**) -------------
        DO k= 1, NumPts
            p= DataPoints(k,1)
            DO j= 1, 2*Degree
                  parr(j,1)= parr(j,1) + p
                p= p*DataPoints(k,1)
              ENDDO
          ENDDO

        ! ---------- Find the matrix entries for the equations --------
        DO r= 1, Degree+1
            DO c= 1, Degree+1
                IF ( r+c .ne. 2 ) THEN
                      a(r,c)= parr(r+c-2,1)
                  ELSE
                      a(r,c)= NumPts
                    ENDIF
              ENDDO
          ENDDO

        ! --------- Solve linear equations for coeffeicients ----------
        CALL MATINVERSE( A      , 3,Maxr         , AI )
        CALL MATMULT   ( AI,b, 3,3,1,MaxR,MaxR,MaxR     , coeff )

        ! ----------------- Find minimum values -----------------------
        MinX= -coeff(2,1) / ( 2*coeff(3,1) )
        MinY= coeff(3,1)*MinX*MinX + coeff(2,1)*MinX + coeff(1,1)

      RETURN
      END
*

*
      SUBROUTINE DMulRSub    ( ALPR,ALPI,BETR,BETI )
        IMPLICIT NONE
        REAL*8 Alpr(4), Alpi(4), Betr(4), Beti(4)
* -----------------------------  Locals  ------------------------------
        REAL*8 Te1,Te2,Te3,Te4,Te5,Te6,Te7,Te8,Te9,Te10,Te11,Te12,
     &         Te13,Te14,Te15,Te16,TeM,DE15,DE16,TemTe7,TemTe8

        ! --------------------  Implementation   ----------------------
         Te1= ALPR(1)-ALPR(3)
         Te2= ALPI(1)-ALPI(3)
         Te5= ALPR(3)-ALPR(2)
         Te6= ALPI(3)-ALPI(2)
         TeM= Te5*Te5+Te6*Te6 

         ! --------------- Check for zero values of tem ---------------
         IF ( DABS( Tem ) .gt. 1.0D-20 ) THEN
             Te3= (Te1*Te5+Te2*Te6)/TeM
             Te4= (Te2*Te5-Te1*Te6)/TeM 
           ELSE
             Te3= 0.0D0
             Te4= 0.0D0 
           ENDIF 

         Te7 = Te3+1.0D0 
         Te9 = Te3*Te3-Te4*Te4 
         Te10= 2.0D0*Te3*Te4 
         DE15= Te7*BETR(3)-Te4*BETI(3)
         DE16= Te7*BETI(3)+Te4*BETR(3)
         Te11= Te3*BETR(2)-Te4*BETI(2)+BETR(1)-DE15
         Te12= Te3*BETI(2)+Te4*BETR(2)+BETI(1)-DE16

         Te7 = Te9-1.0D0 
         Te1 = Te9*BETR(2)-Te10*BETI(2)
         Te2 = Te9*BETI(2)+Te10*BETR(2)
         Te13= Te1-BETR(1)-Te7*BETR(3)+Te10*BETI(3)
         Te14= Te2-BETI(1)-Te7*BETI(3)-Te10*BETR(3)
         Te15= DE15*Te3-DE16*Te4 
         Te16= DE15*Te4+DE16*Te3 

         Te1= Te13*Te13-Te14*Te14-4.0D0*(Te11*Te15-Te12*Te16)
         Te2= 2.0D0*Te13*Te14-4.0D0*(Te12*Te15+Te11*Te16) 

* ---------------------------------------------------------------------
*   Sometimes, for stiff systems (the roots vary widely in order
*   of magnitude), Te1 and Te2 get large enough to have their
*   squares overflow the floating point range.  To prevent this,
*   when either one is large, they are scaled by 10**10 for the
*   purpose of finding TeM.  The scale is restored when the
*   magnitude computation is completed.  This should not affect
*   the accuracy of the computations, since the mantissa is not
*   affected, only the exponent.
* ---------------------------------------------------------------------

         IF ( ( Te1 .gt. 1.0D15 ) .or. ( Te2 .gt. 1.0D15 ) ) THEN
             Te1= Te1*1.0D-10
             Te2= Te2*1.0D-10
             TeM= 1.0D10*DSQRT(Te1*Te1+Te2*Te2)
           ELSE
             TeM= DSQRT(Te1*Te1+Te2*Te2)
           ENDIF

         IF ( Te1 .gt. 0.0D0 ) THEN
             Te3= DSQRT(0.5D0*(TeM+Te1))
             IF ( Te2 .lt. 0.0D0 ) THEN
                 Te3= -Te3
               ENDIF
             ! ---------------- check for zero values of te3 ----------
             IF ( DABS( Te3 ) .lt. 1.0D-15 ) THEN
                 Te4= 0.0D0
               ELSE
                 Te4= 0.5D0*Te2/Te3 
               ENDIF
           ELSE
             Te4= DSQRT(0.5D0*(TeM-Te1))
             ! -------------------- Check for underflows --------------
             IF ( DABS( Te4 ) .lt. 1.0D-15 ) THEN
                 Te3= 0.0D0
               ELSE
                 Te3= 0.5D0*Te2/Te4 
               ENDIF
           ENDIF 

         Te7 = Te13+Te3
         Te8 = Te14+Te4 
         Te9 = Te13-Te3
         Te10= Te14-Te4 
         Te1 = 2.0D0*Te15 
         Te2 = 2.0D0*Te16 

         IF ( (Te7*Te7+Te8*Te8-Te9*Te9-Te10*Te10) .le. 0.0D0 ) THEN
             Te7= Te9
             Te8= Te10 
           ENDIF 
         TeM= Te7*Te7+Te8*Te8 

         TemTe7= Tem*Te7 
         TemTe8= Tem*Te8 

         ! ------------ Check for values of almost zero ---------------
         IF ( DABS( TemTe7 ) .lt. 1.0D-20 ) THEN
             Te3= 0.0D0
             Te4= 0.0D0 
           ELSE
             Te3= Te1/TeM*Te7
             Te4= Te2/TeM*Te7 
           ENDIF 
         IF ( DABS( TemTe8 ) .gt. 1.0D-20 ) THEN
             Te3= Te3 + Te2/TeM*Te8
             Te4= Te4 - Te1/TeM*Te8 
           ENDIF

         ALPR(4)= ALPR(3)+Te3*Te5-Te4*Te6
         ALPI(4)= ALPI(3)+Te3*Te6+Te4*Te5
      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE FACTOR
*
*  this subroutine is a root finding algorithm.  It takes the polynomial .and.
*    returns the roots (real and imaginary) in the Roots array.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Poly        - Array of 16 coefficients
*                    representing the polynomial
*                    (1) is x^8th, (2) is x^7th, ...
*                    others are zero
*    NRoots      - Number of roots
*
*  OutPuts       :
*    Roots       - Array containing roots (real,imag)
*
*  Locals        :
*                -
*                -
*                -
*
*  Coupling      :
*    DMulRSub    -
*
*  References    :
*    This is Bairstows method?
*
* ------------------------------------------------------------------------------  

      SUBROUTINE FACTOR      ( Poly, NRootS, RootS )
        IMPLICIT NONE
        REAL*8 Poly(16), Roots(15,2)
        INTEGER NRootS
* -----------------------------  Locals  ------------------------------
        REAL*8 Alpr(4), Alpi(4), Betr(4), Beti(4), Small, DPoly(16)
        REAL*8 Temp1,    Temp2,    AXR, AXi, PMax, Tem1, Tem2,
     &   TempReal, TempImag, Temp7
        LOGICAL Skip
        INTEGER Mode,LoopCnt,kk,i,j,l,RootCnt

        ! --------------------  Implementation   ----------------------
         Small = 0.000001D0
         PMax= 0.0D0 
         DO KK= 1, NRootS+1
             IF ( DABS(Poly(KK)) .gt. PMax ) THEN
                 PMax= Poly(KK)
               ENDIF
           ENDDO

         IF ( DABS(PMax) .lt. Small ) THEN
             PMax= 1.0D0
           ENDIF
      
         DO KK= 1, NRootS+1
             DPoly(kk)= Poly(kk)/PMax
           ENDDO

         IF ( NRootS .gt. 0 ) THEN
             RootCnt= 0
             i= NRootS+1 

             DO WHILE ((DABS(DPoly(i)) .lt. Small) .and.
     &                 (RootCnt .ne. NRoots))
                 RootCnt= RootCnt+1
                 RootS(RootCnt,1)= 0.0D0
                 RootS(RootCnt,2)= 0.0D0
                 i= i-1
               ENDDO

             IF ( RootCnt .ne. NRootS ) THEN
                 AXR= 0.8D0
                 AXi= 0.0D0
                 L  = 1 
                 LOOPCNT= 1 
                 AlpR(1)= AXR
                 Alpi(1)= AXi
                 MODE= 1 

                 DO WHILE (RootCnt .lt. NRoots)
                     BetR(4)= DPoly(1)
                     Beti(4)= 0.0D0
                     DO i= 1, NRootS
                         J= i+1
                         Temp1= BetR(4)*AXR - Beti(4)*AXi
                         Beti(4)= Beti(4)*AXR + BetR(4)*AXi
                         BetR(4)= Temp1 + DPoly(J)
                        ENDDO

                     TempReal= BetR(4)
                     Tempimag= Beti(4)

                     IF ( RootCnt .ne. 0 ) THEN
                         DO i= 1, RootCnt
                             Tem1   = AXR - RootS(i,1)
                             Tem2   = AXi - RootS(i,2)
                             Temp1  = Tem1*Tem1 + Tem2*Tem2 
                             Temp2  = (BetR(4)*Tem1+Beti(4)*Tem2)/Temp1
                             Beti(4)= (Beti(4)*Tem1-BetR(4)*Tem2)/Temp1
                             BetR(4)= Temp2
                           ENDDO
                       ENDIF 

                     IF (Mode .eq. 1 ) THEN
                             BetR(1)= BetR(4)
                             Beti(1)= Beti(4)
                             AXR    = 0.85D0 
                             AlpR(2)= AXR
                             Alpi(2)= AXi
                             MODE   = 2 
                           ENDIF
                     IF (Mode .eq. 2 ) THEN
                             BetR(2)= BetR(4)
                             Beti(2)= Beti(4)
                             AXR    = 0.9D0 
                             AlpR(3)= AXR
                             Alpi(3)= AXi
                             MODE   = 3
                           ENDIF
                     IF (Mode .eq. 3 ) THEN
                             BetR(3)= BetR(4)
                             Beti(3)= Beti(4)
                             CALL DMULRSUB( AlpR,Alpi,BetR,Beti )
                             AXR    = AlpR(4)
                             AXi    = Alpi(4)
                             MODE   = 4 
                           ENDIF
                     IF (Mode .eq. 5 ) THEN
                             BetR(1)=  BetR(4)
                             Beti(1)=  Beti(4)
                             AXR    =  AlpR(2)
                             AXi    = -Alpi(2)
                             Alpi(2)= -Alpi(2)
                             MODE   =  6 
                           ENDIF
                     IF (Mode .eq. 6 ) THEN
                             BetR(2)=  BetR(4)
                             Beti(2)=  Beti(4)
                             AXR    =  AlpR(3)
                             AXi    = -Alpi(3)
                             Alpi(3)= -Alpi(3)
                             L      =  2 
                             MODE   =  3 
                           ENDIF 

                     ! --------------- the convergence mode -----------
                     IF (Mode .eq. 4 ) THEN
                         Skip= .FALSE.
                         IF ( DABS(TempReal)+DABS(Tempimag) .gt.
     &                          1.0D-20 ) THEN
                             Temp7= DABS(AlpR(3)-AXR)+DABS(Alpi(3)-AXi)
                         IF ( Temp7/(DABS(AXR)+DABS(AXi)) .gt.
     &                          1.0D-7 ) THEN
                             LoopCnt = LoopCnt + 1
                             DO i= 1, 3
                                 AlpR(i)= AlpR(i+1)
                                 Alpi(i)= Alpi(i+1)
                                 BetR(i)= BetR(i+1)
                                 Beti(i)= Beti(i+1)
                               ENDDO
                             IF ( LOOPCNT .lt. 100 ) THEN
                                 CALL DMULRSUB( AlpR,Alpi,BetR,Beti )
                                 AXR = AlpR(4)
                                 AXi = Alpi(4)
                                 MODE= 4
                                 Skip= .TRUE.
                               ENDIF
                           ENDIF

                         IF ( Skip .eqv. .FALSE. ) THEN
                             RootCnt= RootCnt+1
                             RootS(RootCnt,1)= AlpR(4)
                             RootS(RootCnt,2)= Alpi(4)
                             LOOPCNT= 0

                             IF ( RootCnt .lt. NRootS ) THEN
                                 IF ( DABS(RootS(RootCnt,2)) .gt.
     &                                 Small ) THEN
                                     IF ( L .eq. 1 ) THEN
                                         AXR    =  AlpR(1)
                                         AXi    = -Alpi(1)
                                         Alpi(1)= -Alpi(1)
                                         MODE   =  5
                                       ELSE
                                         AXR= 0.8D0
                                         AXi= 0.0D0
                                         L  = 1
                                         LOOPCNT= 1
                                         AlpR(1)= AXR
                                         Alpi(1)= AXi
                                         MODE   = 1
                                       ENDIF
                                   ELSE
                                     AXR= 0.8D0
                                     AXi= 0.0D0
                                     L  = 1
                                     LOOPCNT= 1
                                     AlpR(1)= AXR
                                     Alpi(1)= AXi
                                     MODE   = 1
                                   ENDIF
                               ENDIF
                           ENDIF   ! IF ( Skip
                           ENDIF 
                     ENDIF   ! Case

                   ENDDO   ! DO WHILE

             ENDIF   ! If  RootCnt .ne. NRoots

           ENDIF   ! If  nRoots .gt. 0

      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE QUADRATIC
*
*  this subroutine solves for the two Roots of a QUADRATIC equation.  There are
*    no restrictions on the coefficients, and imaginary results are passed
*    out as separate values.  The general form is y = ax2 + bx + c.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    a           - Coefficient of x squared term
*    b           - Coefficient of x term
*    c           - Constant
*
*  OutPuts       :
*    R1r         - Real portion of Root 1
*    R1i         - Imaginary portion of Root 1
*    R2r         - Real portion of Root 2
*    R2i         - Imaginary portion of Root 2
*
*  Locals        :
*    Discrim     - Discriminate b2 - 4ac
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, pg 974
*
* ------------------------------------------------------------------------------  

      SUBROUTINE QUADRATIC( a,b,c,opt, R1r,R1i,R2r,R2i )
        IMPLICIT NONE
        REAL*8 R1r,R1i,R2r,R2i,a,b,c
        CHARACTER opt
* -----------------------------  Locals  ------------------------------
    
        REAL*8 discrim, small
        DATA Small      /0.00000001D0/
        ! --------------------  implementation   ----------------------
        r1r = 0.0D0
        r1i = 0.0D0
        r2r = 0.0D0
        r2i = 0.0D0

        discrim = b * b - 4.0D0 * a * c

        ! ---------------------  real roots  --------------------------
        if (DAbs(discrim) .lt. small) THEN
                r1r = -b / (2.0D0 * a)
                r2r = r1r
                if (opt .eq.'U') THEN
                    r2r = 99999.9D0
                ENDIF
            else
                if (DAbs(a) .lt. small) THEN
                    r1r = -c / b
                else
                    if (discrim .gt. 0.0D0) THEN
                        r1r = (-b + DSqrt(discrim)) / (2.0D0 * a)
                        r2r = (-b - DSqrt(discrim)) / (2.0D0 * a)
                    else
                        ! ------------------ complex roots --------------------
                        if (opt .eq. 'I') THEN
                            r1r = -b / (2.0D0 * a)
                            r2r = r1r;
                            r1i = DSqrt(-discrim) / (2.0 * a)
                            r2i = -DSqrt(-discrim) / (2.0 * a)
                        else
                            r1r = 99999.9D0
                            r2r = 99999.9D0
                        ENDIF
                    ENDIF
                ENDIF
            ENDIF
        RETURN
        END  ! QUADRATIC
        
        
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE CUBIC
*
*  this subroutine solves for the three Roots of a CUBIC equation.  There are
*    no restrictions on the coefficients, and imaginary results are passed
*    out as separate values.  The general form is y = ax3 + bx2 + cx + d.  Note
*    that R1i will ALWAYS be ZERO since there is ALWAYS at least one REAL Root.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    a           - Coefficient of x cubed term
*    b           - Coefficient of x squared term
*    c           - Coefficient of x term
*    d           - Constant
*
*  OutPuts       :
*    R1r         - Real portion of Root 1
*    R1i         - Imaginary portion of Root 1
*    R2r         - Real portion of Root 2
*    R2i         - Imaginary portion of Root 2
*    R3r         - Real portion of Root 3
*    R3i         - Imaginary portion of Root 3
*
*  Locals        :
*    Temp1       - Temporary value
*    Temp2       - Temporary value
*    Root1       - Temporary value of the Root
*    Root2       - Temporary value of the Root
*    Root3       - Temporary value of the Root
*    P           - Coefficient of x squared term where x cubed term is 1.0D0
*    Q           - Coefficient of x term where x cubed term is 1.0D0
*    R           - Coefficient of constant term where x cubed term is 1.0D0
*    Delta       - Discriminator for use with Cardans formula
*    E0          - ANGLE holder for trigonometric solution
*    Phi         - ANGLE used in trigonometric solution
*    CosPhi      - Cosine of Phi
*    SinPhi      - Sine of Phi
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, pg 975
*
* ------------------------------------------------------------------------------
*
      SUBROUTINE CUBIC       ( a,b,c,d, R1r,R1i,R2r,R2i,R3r,R3i )
        IMPLICIT NONE
        REAL*8 R1r,R1i,R2r,R2i,R3r,R3i,a,b,c,d
* -----------------------------  Locals  ------------------------------
        REAL*8 OneThird, aa, bb, cc, dd,
     &   Temp1, Temp2, Root1, Root2, Root3, P, Q, R, Delta,
     &   E0, CosPhi, SinPhi, Phi

        REAL*8 Small,    Rad2Deg,  Deg2Rad, Pi
        DATA Small      /0.00000001D0/
c       use machine precision instead
        Pi      =  4.0D0 * DATAN(1.0D0)
        Rad2Deg = 180.0D0/pi
        Deg2Rad = pi/180.0D0

        ! --------------------  Implementation   ----------------------
        OneThird  = 1.0D0/3.0D0
        R1r  = 0.0D0
        R1i  = 0.0D0
        R2r  = 0.0D0
        R2i  = 0.0D0
        R3r  = 0.0D0
        R3i  = 0.0D0
        Root1= 0.0D0 
        Root2= 0.0D0 
        Root3= 0.0D0 

        ! ----------- Force coefficients into std form ----------------
        IF (DABS(a) .gt. small) THEN
         P= b/a 
         Q= c/a 
         R= d/a 

         aa= OneThird*( 3.0D0*Q - P*P ) 
         bb= (1.0D0/27.0D0)*( 2.0D0*P*P*P - 9.0D0*P*Q + 27.0D0*R )

         Delta= (aa*aa*aa/27.0D0) + (bb*bb*0.25D0) 
         
        ! ------------------ Use Cardans formula ----------------------
        IF ( Delta .gt. Small ) THEN
            Temp1= (-bb*0.5D0)+DSQRT(Delta)
            Temp2= (-bb*0.5D0)-DSQRT(Delta) 
            IF ( DABS(Temp1) .gt. Small ) THEN
                Temp1= Temp1**OneThird
              ENDIF
            IF ( DABS(Temp2) .gt. Small ) THEN
                Temp2= Temp2**OneThird
              ENDIF
            Root1= Temp1 + Temp2 
            Root2= -0.5D0*(Temp1 + Temp2) 
            Root3= -0.5D0*(Temp1 + Temp2) 
            R2i= -0.5D0*DSQRT( 3.0D0 )*(Temp1 - Temp2) 
            R3i= -R2i 
          ELSE
            ! --------------- Evaluate zero point ---------------------
            IF ( DABS( Delta ) .lt. Small  ) THEN
                IF ( DABS(bb) .gt. Small ) THEN
                    Root1= -2.0D0 * (bb*0.5D0)**OneThird
                    Root2=  (bb*0.5D0)**OneThird
                    Root3= Root2 
                  ENDIF 
              ELSE
                ! ------------ Use trigonometric identities -----------
                E0     = 2.0D0*DSQRT(-aa*OneThird)
                CosPhi = (-bb/(2.0D0*DSQRT(-aa*aa*aa/27.0D0)) ) 
                SinPhi = DSQRT( 1.0D0-CosPhi*CosPhi ) 
                Phi    = DATAN2( SinPhi,CosPhi ) 
                Root1= E0*DCOS( Phi*OneThird ) 
                Root2= E0*DCOS( Phi*OneThird + 120.0D0*deg2Rad )
                Root3= E0*DCOS( Phi*OneThird + 240.0D0*deg2Rad )
                R1r= Root1 - P*OneThird 
                R2r= Root2 - P*OneThird 
                R3r= Root3 - P*OneThird 
              ENDIF 
          ENDIF 
          ELSE
          CALL quadratic ( b, c, d, 'R', r1r,r1i,r2r,r2i )
          r3r  = 99999.9D0
          r3i  = 99999.9D0
        
          ENDIF

      RETURN
      END   !  CUBIC

*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE QUARTIC
*
*  this subroutine solves for the four Roots of a QUARTIC equation.  There are
*    no restrictions on the coefficients, and imaginary results are passed
*    out as separate values.  The general form is y = ax4 + bx3 + cx2 + dx + e.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    a           - Coeficient of x fourth term
*    b           - Coefficient of x cubed term
*    c           - Coefficient of x squared term
*    d           - Coefficient of x term
*    e           - Constant
*
*  OutPuts       :
*    R1r         - Real portion of Root 1
*    R1i         - Imaginary portion of Root 1
*    R2r         - Real portion of Root 2
*    R2i         - Imaginary portion of Root 2
*    R3r         - Real portion of Root 3
*    R3i         - Imaginary portion of Root 3
*    R4r         - Real portion of Root 4
*    R4i         - Imaginary portion of Root 4
*
*  Locals        :
*    Temp1       - Temporary value
*    Temp2       - Temporary value
*    Root1       - Temporary value of the Root
*    Root2       - Temporary value of the Root
*    Root3       - Temporary value of the Root
*    s           - alternate variable
*    h           - Temporary value
*    hSqr        - h squared
*    hCube       - h Cubed
*    P           - Term in auxillary equation
*    Q           - Term in auxillary equation
*    R           - Term in auxillary equation
*    Delta       - Discriminator for use with Cardans formula
*    E0          - ANGLE holder for trigonometric solution
*    Phi         - ANGLE used in trigonometric solution
*    CosPhi      - Cosine of Phi
*    SinPhi      - Sine of Phi
*    RPrime      - Values of Roots before final work
*    Temp        - Temporary variable in finding MAX RPrime
*    Eta         - Constant coefficient in QUADRATIC solutions
*    Beta        - Constant coefficient in QUADRATIC solutions
*
*  Coupling      :
*    QUADRATIC     Find roots of a QUADRATIC polynomial
*
*  References    :
*    Vallado       2007, pg 976
*
* ------------------------------------------------------------------------------  

      SUBROUTINE QUARTIC     ( a,b,c,d,e, R1r,R1i,R2r,R2i,
     &                         R3r,R3i,R4r,R4i )
        IMPLICIT NONE
        REAL*8 R1r,R1i,R2r,R2i,R3r,R3i,R4r,R4i,a,b,c,d,e
* -----------------------------  Locals  ------------------------------
        REAL*8 OneThird,
     &   Temp1, Temp2, Root1, Root2, Root3, s, h, P, Q, R, Delta, E0,
     &   CosPhi, SinPhi, Phi, RPrime, hSqr, HCube, Eta, Beta
        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        OneThird  = 1.0D0/3.0D0
        R1r  = 0.0D0 
        R1i  = 0.0D0 
        R2r  = 0.0D0 
        R2i  = 0.0D0 
        R3r  = 0.0D0
        R3i  = 0.0D0 
        R4r  = 0.0D0 
        R4i  = 0.0D0
        Root1= 0.0D0 
        Root2= 0.0D0 
        Root3= 0.0D0 
        ! ----------- Force coefficients into std form ----------------
        b= B/A 
        c= C/A 
        d= D/A 
        e= E/A 

        H= -b/4 
        HSqr= H**2
        HCube= HSqr * H 

        P=                      6.0D0*HSqr   + 3.0D0*b*h + c 
        Q=          4.0D0*HCube + 3.0D0*b*HSqr + 2.0D0*c*h + d 
        R= h*HCube +  b*HCube +   c*HSqr   +    d*h  + e 

        a=  OneThird*( -P*P-12.0D0*R )
        b=  (1.0D0/27.0D0)*( -2.0D0*P*P*P+72.0D0*P*R-27.0D0*Q*Q )
        s= -2.0D0*OneThird*P 

        Delta= (a*a*a/27.0D0) + (b*b*0.25D0) 

        IF ( DABS(Q) .gt. Small ) THEN
            ! ------------------ Use Cardans formula ------------------
            IF ( Delta .gt. Small ) THEN
                Temp1= (-b*0.5D0)+DSQRT(Delta)
                Temp2= (-b*0.5D0)-DSQRT(Delta) 
                IF ( DABS(Temp1) .gt. Small ) THEN
                    Temp1= Temp1**OneThird
                  ENDIF
                IF ( DABS(Temp2) .gt. Small ) THEN
                    Temp2= Temp2**OneThird
                  ENDIF   
                Root1= Temp1 + Temp2 
                Root2= -0.5D0*(Temp1 + Temp2) 
                Root3= -0.5D0*(Temp1 + Temp2) 
                R2i= -0.5D0*DSQRT( 3.0D0 )*(Temp1 - Temp2) 
                R3i= -R2i 
              ELSE
                ! --------------- Evaluate zero point -----------------
                IF ( DABS( Delta ) .lt. Small  ) THEN
                    IF ( DABS(b) .gt. Small ) THEN
                        Root1= -2.0D0* (b*0.25D0)**OneThird
                        Root2=  (b*0.25D0)**OneThird
                        Root3= Root2 
                      ENDIF 
                  ELSE
                    ! ------ Use trigonometric identities -------------
                    E0    = 2.0D0*DSQRT(-a*OneThird)
                    CosPhi= (-b/(2.0D0*DSQRT(-a*a*a/27.0D0)) ) 
                    SinPhi= DSQRT( 1.0D0-CosPhi*CosPhi ) 
                    Phi   = DATAN2( SinPhi,CosPhi ) 
                    Root1= E0*DCOS( Phi*OneThird ) 
                    Root2= E0*DCOS( Phi*OneThird + 120.0D0*deg2Rad )
                    Root3= E0*DCOS( Phi*OneThird + 240.0D0*deg2Rad )
                  ENDIF 
              ENDIF 

             ! --------------- Find largest value of Root -------------
             RPrime= Root1+s 
             IF ( (RPrime .lt. Root2+s) .and.
     &            (DABS(R2i).lt.0.0001D0) ) THEN
                 RPrime= Root2+s
              ENDIF
             IF ( (RPrime .lt. Root3+s) .and.
     &            (DABS(R3i).lt.0.0001D0) ) THEN
                 RPrime= Root3+s 
              ENDIF

             ! -- Evaluate coefficients of two resulting Quadratics ---
             IF ( RPrime .gt. Small ) THEN
                 Eta = 0.5D0*( P + RPrime - Q/DSQRT(RPrime) )
                 Beta= 0.5D0*( P + RPrime + Q/DSQRT(RPrime) ) 
               ELSE
                 Eta = 0.5D0*P
                 Beta= 0.5D0*P 
                 IF ( RPrime .lt. 0.0D0 ) THEN
                     RPrime= -RPrime 
                   ENDIF
               ENDIF 

             CALL QUADRATIC( 1.0D0, DSQRT(RPrime),Eta,
     &                       R1r,R1i,R2r,R2i )
             CALL QUADRATIC( 1.0D0,-DSQRT(RPrime),Beta,
     &                       R3r,R3i,R4r,R4i )

           ELSE
             ! ------ Case where solution reduces to a QUADRATIC -------
             CALL QUADRATIC( 1.0D0,P,R,   R1r,R1i,R3r,R3i )
             R  = DSQRT( R1r*R1r + R1i*R1i ) 
             IF ( DABS(R) .gt. Small ) THEN
                 Phi= DATAN2( R1i/R,R1r/R )
               ELSE
                 Phi= 0.0D0
               ENDIF
             R1r= DSQRT(R) * DCOS(Phi*0.5D0)
             R1i= DSQRT(R) * DSIN(Phi*0.5D0) 
             IF ( DABS(R1i) .gt. Small ) THEN
                 R2r= R1r
               ELSE
                 R2r= -R1r
               ENDIF
             R2i= -R1i 

             R  = DSQRT( R3r*R3r + R3i*R3i ) 
             IF ( DABS(R) .gt. Small ) THEN
                 Phi= DATAN2( R3i/R,R3r/R )
               ELSE
                 Phi= 0.0D0
               ENDIF
             R3r= DSQRT(R) * DCOS(Phi*0.5D0) 
             R3i= DSQRT(R) * DSIN(Phi*0.5D0) 
             IF ( DABS(R3i) .gt. Small ) THEN
                 R4r= R3r
               ELSE
                 R4r= -R3r
               ENDIF
             R4i= -R3i 
           ENDIF 

        R1r= R1r + h 
        R2r= R2r + h
        R3r= R3r + h 
        R4r= R4r + h 
      RETURN
      END
      
* -----------------------------------------------------------------------------
*
*                           function cubicspl
*
*  this function performs cubic splining of an input zero crossing
*  function in order to find function values.
*
*  author        : david vallado           davallado@gmail.com     7 aug 2005
*
*  revisions
*                -
*  inputs          description                    range / units
*    p0,p1,p2,p3 - function values used for splining
*    t0,t1,t2,t3 - time values used for splining
*
*  outputs       :
*    acu0..acu3  - splined polynomial coefficients. acu3 t^3, etc
*
*  locals        : none
*
*  coupling      :
*    none
*
*  references    :
*    vallado       2013, 1034
* ---------------------------------------------------------------------------

      SUBROUTINE cubicspl( p1, p2, p3, p4,  acu0, acu1, acu2, acu3 )
        IMPLICIT NONE
        REAL*8 p1, p2, p3, p4, acu0, acu1, acu2, acu3
        
        ! --------------------  Implementation   ----------------------
        acu0 = p2
        acu1 = -p1 / 3.0D0 - 0.5D0 * p2 + p3 - p4 / 6.0D0
        acu2 = 0.5D0 * p1 - p2 + 0.5D0 * p3
        acu3 = -p1 / 6.0D0 + 0.5D0 * p2 - 0.5D0 * p3 + p4 / 6.0D0
      RETURN
      END



* -----------------------------------------------------------------------------
*
*                           function cubicinterp
*
*  this function performs a cubic spline. four points are needed.
*
*  author        : david vallado           davallado@gmail.com   1 dec  2005
*
*  inputs          description                    range / units
*    valuein     - kp
*
*  outputs       :
*    out         - ap
*
*  locals        :
*                -
*
*  coupling      :
*    cubicspl
*
*  references    :
*    vallado       2013, 1027
* --------------------------------------------------------------------------- 

      REAL*8 FUNCTION cubicinterp( p1a, p1b, p1c, p1d, p2a, p2b,
     &       p2c, p2d, valuein )
        IMPLICIT NONE
        REAL*8 p1a, p1b, p1c, p1d, p2a, p2b, p2c, p2d, valuein
* -----------------------------  Locals  ------------------------------
        REAL*8 kc0, kc1, kc2, kc3, ac0, ac1, ac2, ac3,
     &         r1r, r1i, r2r, r2i, r3r, r3i, val

        ! --------------------  Implementation   ----------------------

        
        ! -------- assign function points ---------
        CALL cubicspl(p1a, p1b, p1c, p1d, ac0, ac1, ac2, ac3)
        CALL cubicspl(p2a, p2b, p2c, p2d, kc0, kc1, kc2, kc3)

        ! recover the original function values
        ! use the normalized time first, but at an arbitrary interval
        CALL cubic(kc3, kc2, kc1, kc0 - valuein,  
     &             r1r, r1i, r2r, r2i, r3r, r3i)

        if ((r1r.ge.-0.000001D0).and.(r1r.le.1.001D0)) THEN
            val = r1r
        else
            if ((r2r.ge.-0.000001D0).and.(r2r.le.1.001D0)) THEN
                val = r2r
            else
                if ((r3r.ge.-0.000001D0).and.(r3r.le.1.001D0)) THEN
                    val = r3r
                else
                    val = 0.0D0
                    !Console.Write("error in cubicinterp root {0} {1} {2} {3} \n", valuein, r1r, r2r, r3r)
                ENDIF
            ENDIF
        ENDIF

        cubicinterp = ac3*val**3 + ac2*val**2 + ac1*val + ac0
      RETURN
      END   !  cubicinterp
      
  

* ------------------------------------------------------------------------------
*
*                           FUNCTION FACTORIAL
*
*  this function finds the value of a FACTORIAL.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    x           - Input value
*
*  Outputs       :
*    FUNCTION    - answer
*
*  Locals        :
*    Temp        - Temporary variable
*    i           - Index
*
*  Coupling      :
*    None.
*
* ------------------------------------------------------------------------------

      REAL*8 FUNCTION  FACTORIAL ( x )
        IMPLICIT NONE
        INTEGER x
* -----------------------------  Locals  ------------------------------
        REAL*8 Temp
        INTEGER i

        ! --------------------  Implementation   ----------------------
        Temp= 1.0D0
        DO i= 2, x
            Temp= Temp * i
          ENDDO
        FACTORIAL= Temp
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           FUNCTION BINOMIAL
*
*  this function finds the value of a BINOMIAL coefficient
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    i           -
*    j           -
*
*  Outputs       :
*    FUNCTION    - answer
*
*  Locals        :
*    None.
*
*  Coupling      :
*    FACTORIAL   - Finds the FACTORIAL of a number
*
* ------------------------------------------------------------------------------

      REAL*8 FUNCTION  BINOMIAL ( i,j )
        IMPLICIT NONE
        INTEGER i,j
        EXTERNAL Factorial
* -----------------------------  Locals  ------------------------------
        REAL*8 FACTORIAL

        ! --------------------  Implementation   ----------------------
        BINOMIAL= FACTORIAL(j) / ( FACTORIAL(i)*FACTORIAL(j-i) )
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           FUNCTION MIN
*
*  this function determines the minimum of 2 values.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    X           - Value number 1
*    Y           - Value number 2
*
*  OutPuts       :
*    MIN         - Minimum of x or y
*
*  Locals        :
*    None.
*
*  Coupling      :
*    None.
*
* ------------------------------------------------------------------------------

      REAL*8 FUNCTION MIN ( X, Y )
        IMPLICIT NONE
        REAL*8 X,Y

        ! --------------------  Implementation   ----------------------
        IF ( X .lt. Y ) THEN
            MIN = X
          ELSE
            MIN = Y
          ENDIF
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           FUNCTION MAX
*
*  this function determines the maximum of 2 values.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    X           - Value number 1
*    Y           - Value number 2
*
*  OutPuts       :
*    MAX         - Minimum of x or y
*
*  Locals        :
*    None.
*
*  Coupling      :
*    None.
*
* ------------------------------------------------------------------------------

      REAL*8 FUNCTION MAX    ( X, Y )
        IMPLICIT NONE
        REAL*8 X,Y

        ! ----------------------  Implementation   --------------------
        IF ( X .gt. Y ) THEN
            MAX = X
          ELSE
            MAX = Y
          ENDIF
      RETURN
      END

* ------------------------------------------------------------------------------
*
*                           SUBROUTINE PLANE
*
* this subroutine calculates the equation of a PLANE given 3 points
*   pt1 - x1,y1,z1, pt2 - x2,y2,z2, pt3 - x3,y3,z3 , and outputs the
*   a b c d variables describing the PLANE. NOTE that the general equation
*   of a PLANE is defined here as:  ax + by + cz + d = 0  and the values
*   are obtained by solving the ordered determinant  x  y  z   1
*                                                    x1 y1 z1  1   =0
*                                                    x2 y2 z2  1
*                                                    x3 y3 z3  1
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    x1,y1,z1    - point # 1
*    x2,y2,z2    - point # 2
*    x3,y3,z3    - point # 3
*
*  OutPuts       :
*    a,b,c,d     - constants for the equation of the PLANE
*
*  Locals        :
*    z23
*
*  Coupling      :
*    None.
*
* ------------------------------------------------------------------------------  

      SUBROUTINE PLANE       ( x1,y1,z1,x2,y2,z2,x3,y3,z3, a,b,c,d )
        IMPLICIT NONE
        REAL*8 x1,y1,z1,x2,y2,z2,x3,y3,z3,a,b,c,d
* -----------------------------  Locals  ------------------------------
        REAL*8 z23,y23,x23,yz23,yz32,xz23,xz32,xy23,xy32

        z23= z2-z3 
        y23= y2-y3 
        x23= x2-x3 

        yz23= y2*z3 
        yz32= y3*z2 
        xz23= x2*z3
        xz32= x3*z2 
        xy23= x2*y3 
        xy32= x3*y2 

        a= y1*z23 - z1*y23 + yz23 - yz32 
        b= -(x1*z23 - z1*x23 + xz23 - xz32)
        c= x1*y23 - y1*x23 + xy23 - xy32 
        d= -(x1*(yz23 - yz32) - y1*(xz23-xz32) + z1*(xy23 -xy32))
      RETURN
      END

      SUBROUTINE INITTIME
        IMPLICIT NONE
* ----------------------------  Locals  -------------------------------
        INTEGER i, LMonth(12)
        CHARACTER*3 MonthTitle(12), DayTitle(12)

        ! --------------------  Implementation   ----------------------
        DO i = 1,12
            LMonth(i) = 31
          ENDDO
        LMonth( 2) = 28
        LMonth( 4) = 30
        LMonth( 6) = 30
        LMonth( 9) = 30
        LMonth(11) = 30

        MonthTitle( 1)= 'Jan'
        MonthTitle( 2)= 'Feb'
        MonthTitle( 3)= 'Mar'
        MonthTitle( 4)= 'Apr'
        MonthTitle( 5)= 'May'
        MonthTitle( 6)= 'Jun'
        MonthTitle( 7)= 'Jul'
        MonthTitle( 8)= 'Aug'
        MonthTitle( 9)= 'Sep'
        MonthTitle(10)= 'Oct'
        MonthTitle(11)= 'Nov'
        MonthTitle(12)= 'Dec'

        DayTitle(1)= 'Sun'
        DayTitle(2)= 'Mon'
        DayTitle(3)= 'Tue'
        DayTitle(4)= 'Wed'
        DayTitle(5)= 'Thr'
        DayTitle(6)= 'Fri'
        DayTitle(7)= 'Sat'

      RETURN
      END
*
* -----------------------------------------------------------------------------
*
*                           FUNCTION GETINTMON
*
*  this function finds the INTEGER equivalent of the 3 character string
*    representation of month.
*
*  Inputs          Description                    Range / Units
*    MonStr      - Month name                     'Jan','Feb' ...
*
*  OutPuts       :
*    GETINTMON   - INTEGER Month equivalent       1 .. 12
*
*  Locals        :
*    i           - Index
*
*  Coupling      :
*
* -----------------------------------------------------------------------------

      INTEGER FUNCTION GETINTMON       ( MonStr )
        CHARACTER*3 MonStr
* ----------------------------  Locals  -------------------------------
        INTEGER i
        CHARACTER*3 MonthTitle(12),MonthlTitle(12)

        ! --------------------  Implementation   ----------------------
        MonthTitle( 1)= 'JAN'
        MonthTitle( 2)= 'FEB'
        MonthTitle( 3)= 'MAR'
        MonthTitle( 4)= 'APR'
        MonthTitle( 5)= 'MAY'
        MonthTitle( 6)= 'JUN'
        MonthTitle( 7)= 'JUL'
        MonthTitle( 8)= 'AUG'
        MonthTitle( 9)= 'SEP'
        MonthTitle(10)= 'OCT'
        MonthTitle(11)= 'NOV'
        MonthTitle(12)= 'DEC'

        MonthlTitle( 1)= 'Jan'
        MonthlTitle( 2)= 'Feb'
        MonthlTitle( 3)= 'Mar'
        MonthlTitle( 4)= 'Apr'
        MonthlTitle( 5)= 'May'
        MonthlTitle( 6)= 'Jun'
        MonthlTitle( 7)= 'Jul'
        MonthlTitle( 8)= 'Aug'
        MonthlTitle( 9)= 'Sep'
        MonthlTitle(10)= 'Oct'
        MonthlTitle(11)= 'Nov'
        MonthlTitle(12)= 'Dec'

        i = 1
        DO WHILE (( MonthTitle(i).ne.MonStr ).and.( i .lt. 12 ))
            i= i+1
          ENDDO

        IF ( (i.eq.12) .and. (MonthTitle(i).ne.MonStr) ) THEN
            i = 1
            DO WHILE (( MonthlTitle(i).ne.MonStr ).and.( i .lt. 12 ))
                i= i+1
              ENDDO
          ENDIF

        GETINTMON= i
      RETURN
      END
*
* -----------------------------------------------------------------------------
*
*                           FUNCTION GETINTDAY
*
*  this function finds the INTEGER equivalent of the 3 character string
*    representation of the day of the week.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    DayStr      - Day name string                'Sun','Mon' ...
*
*  OutPuts       :
*    GETINTDAY   - INTEGER Day equivalent         1 .. 7
*
*  Locals        :
*    i           - Index
*
*  Coupling      :
*
* -----------------------------------------------------------------------------

      INTEGER FUNCTION GETINTDAY       ( DayStr )
        CHARACTER*3 DayStr
* ----------------------------  Locals  -------------------------------
        INTEGER i
        CHARACTER*3 DayTitle(12)

        ! --------------------  Implementation   ----------------------
        DayTitle(1)= 'Sun'
        DayTitle(2)= 'Mon'
        DayTitle(3)= 'Tue'
        DayTitle(4)= 'Wed'
        DayTitle(5)= 'Thr'
        DayTitle(6)= 'Fri'
        DayTitle(7)= 'Sat'

        i= 1
        DO WHILE (( (DayTitle(i)) .ne. (DayStr) ) .and.
     &           ( i .lt. 12 ))
            i= i+1
          ENDDO
        GETINTDAY= i
      RETURN
      END

* -----------------------------------------------------------------------------
*
*                           SUBROUTINE DAYOFWEEK
*
*  this function finds the day of the week. Integers are used for the days,
*    1 = 'SUN', 2 = 'Mon', ... 7 = 'Sat'.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    JD          - Julian date of interest        days from 4713 BC
*
*  OutPuts       :
*    DAYOFWEEK   - answer                         1 to 7
*
*  Locals        :
*    None.
*
*  Coupling      :
*    None.
*
*  References    :
*    vallado       2007, 188, Eq 3-39
*
* -----------------------------------------------------------------------------

      INTEGER FUNCTION DAYOFWEEK       ( JD )
        REAL*8 JD

        ! --------------------  Implementation   ----------------------
        ! ------- Be sure JD is at 0.0D0 h on the day of interest -----
        JD       = DINT(JD+0.5D0)
        DAYOFWEEK= IDINT( JD - 7 * IDINT( (JD+1)/7 ) + 2 )
      RETURN
      END
*
* -----------------------------------------------------------------------------
*
*                           SUBROUTINE DAYLIGHTST
*
*  this subroutine finds the dates for switiching to daylight savings time in
*    a given year. The date is set as the 2nd sunday in march and the first
*    sunday in november. The DST dates are adjusted -10hr to get 0200, -Zone to get
*    the local time zone, and -1.0hr to process the stop because the local time is
*    on DST before a stop.
*
*  Author        : David Vallado                  719-573-2600   17 Mar 2007
*
*  Inputs          Description                    Range / Units
*    Year        - Year                           1900 .. 2100
*    Lon         - SITE longitude (WEST -)        -2Pi to 2Pi rad
*
*  Outputs       :
*    StartDay    - Day in April when DST begins   1 .. 28,29,30,31
*    StopDay     - Day in October when DST ends   1 .. 28,29,30,31
*
*  Locals        :
*    DW          - Day of the week                1 .. 7
*    JDStartDST  - Julian date of start           Days from 4713 BC
*    JDStopDST   - Julian date of stop            Days from 4713 BC
*    Zone        - Time zone of site. Default
*                  of 0.0 gives Greenwich         hrs
*
*  Coupling      :
*    JDAY   - Find the Julian Date
*
*  References    :
*    Vallado       2007, 188
*
* -----------------------------------------------------------------------------

      SUBROUTINE DAYLIGHTST ( Year, Lon, StartDay, StopDay, JDStartDST,
     &                        JDStopDST )
        IMPLICIT NONE
        INTEGER Year, StartDay, StopDay
        REAL*8  Lon, JDStartDST, JDStopDST
* ----------------------------  Locals  -------------------------------
        REAL*8 Zone
        INTEGER DW
        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        ! ------- Find time zone information to adjust to a site ------ 
        Zone= DINT( Lon*Rad2Deg/15.0D0 )
        IF (Zone .gt. 0.0D0) THEN
            Zone= Zone - 24.0D0
          ENDIF

        IF (year .lt. 2007) THEN
             StartDay= 0
            Dw = 0
            DO WHILE ((DW .ne. 1) .and. (StartDay .ne. 7)) ! 1 is Sunday
                StartDay= StartDay + 1
                CALL JDay( Year,4,StartDay,12,0,0.0D0, JDStartDST )
                DW= IDINT( JDStartDST - 7* IDINT( (JDStartDST+1)/7 )+2)
              ENDDO
            JDStartDST= JDStartDST - 10.0D0/24.0D0   ! set to 0200 UTC

            StopDay= 32
            DO WHILE ((DW .ne. 1) .and. (StopDay .ne. 25)) ! 1 is Sunday
                StopDay= StopDay - 1
                CALL JDay( Year,10,StopDay,12,0,0.0D0, JDStopDST )
                DW= IDINT( JDStopDST - 7* IDINT( (JDStopDST+1)/7 ) + 2 )
              ENDDO
            JDStopDST= JDStopDST - 10.0D0/24.0D0   ! set to 0200 UTC
          ELSE
           StartDay = 7
            Dw = 0
            DO WHILE ((DW .ne. 1) .and. (StartDay .ne. 15)) ! 1 is Sunday
                StartDay= StartDay + 1
                CALL JDay( Year,3,StartDay,12,0,0.0D0, JDStartDST )
                DW= IDINT( JDStartDST - 7* IDINT( (JDStartDST+1)/7 )+2)
              ENDDO
            JDStartDST= JDStartDST - 10.0D0/24.0D0   ! set to 0200 UTC

            StopDay = 0
            DO WHILE ((DW .ne. 1) .and. (StopDay .ne. 8)) ! 1 is Sunday
                StopDay= StopDay + 1
                CALL JDay( Year,11,StopDay,12,0,0.0D0, JDStopDST )
                DW= IDINT( JDStopDST - 7* IDINT( (JDStopDST+1)/7 ) + 2 )
              ENDDO
            JDStopDST= JDStopDST - 10.0D0/24.0D0   ! set to 0200 UTC
          ENDIF

      RETURN
      END
      
* -----------------------------------------------------------------------------
*
*                           SUBROUTINE JDay
*
*  this subroutine finds the Julian date given the Year, Month, Day, and Time.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Year        - Year                           1900 .. 2100
*    Mon         - Month                          1 .. 12
*    Day         - Day                            1 .. 28,29,30,31
*    Hr          - Universal Time Hour            0 .. 23
*    minute         - Universal Time minute             0 .. 59
*    Sec         - Universal Time Sec             0.0D0 .. 59.999D0
*    WhichType   - Julian .or. Gregorian calender   'J' .or. 'G'
*
*  Outputs       :
*    jd          - julian date (days only)           days from 4713 bc
*    jdFrac      - julian date (fraction of a day)   days from 0 hr of the day
*
*  Locals        :
*    B           - Var to aid Gregorian dates
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 189, Alg 14, Ex 3-14
* -----------------------------------------------------------------------------

      SUBROUTINE JDay       ( Year,Mon,Day,Hr,minute, Sec, JD, JDFrac )
        IMPLICIT NONE
        INTEGER Year, Mon, Day, Hr, minute
        REAL*8  Sec, JD, dtt, JDFrac

        ! --------------------  Implementation   ----------------------
        JD= 367.0D0 * Year
     &        - INT( (7* (Year+INT ( (Mon+9)/12) ) ) * 0.25D0 )
     &        + INT( 275*Mon / 9 )
     &        + Day + 1721013.5D0  ! use - 678987.0 to go to mjd directly
        JDFrac = ( (Sec/60.0D0 + minute ) / 60.0D0 + Hr ) / 24.0D0
     
*       check that the day and fractional day are correct
		IF (DABS(JDFrac).gt.1.0D0) THEN
			dtt = INT(JDFrac)
			JD = JD + dtt
			JDFrac = JDFrac - dtt
		ENDIF	
*     &      - 0.5D0*DSIGN(1.0D0, 100.0D0*Year + Mon - 190002.5D0) + 0.5D0
      RETURN
      END


*  References    :
*    Vallado       2007, 190

      SUBROUTINE JDayALL   ( Year,Mon,Day,Hr,minute,Sec, WhichType, JD )
        IMPLICIT NONE
        INTEGER Year, Mon, Day, Hr, minute
        CHARACTER WhichType
        REAL*8  Sec, JD
* ----------------------------  Locals  -------------------------------
        REAL*8 B

        ! --------------------  Implementation   ----------------------
        IF ( Mon .le. 2 ) THEN
            Year= Year - 1
            Mon = Mon + 12
          ENDIF
        IF ( WhichType .eq. 'J' ) THEN
            ! -------- Use for Julian calender, every 4 years ---------
            B= 0.0D0
          ELSE
            ! --------------------- Use for Gregorian -----------------
            B= 2 - DINT(Year*0.01D0) + DINT(DINT(Year*0.01D0)*0.25D0)
          ENDIF
        JD= DINT( 365.25D0*(Year + 4716) )
     &       + DINT( 30.6001D0*(Mon+1) )
     &       + Day + B - 1524.5D0
     &       + ( (Sec/60.0D0 + minute ) / 60.0D0 + Hr ) / 24.0D0
      RETURN
      END

* -----------------------------------------------------------------------------
*
*                           SUBROUTINE DAYS2MDHMS
*
*  this subroutine converts the day of the year, days, to the equivalent month
*    day, hour, Minute and second.
*
*  Algorithm     : Set up array for the Number of days per month
*                  Find Leap Year - be sure to account for the 400 years
*                  Loop through a Temp value for WHILE the value is .lt. the days
*                  Perform INTEGER conversions to the correct day and month
*                  Convert remainder into H M S using type conversions
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Year        - Year                          +1900 .. 2100+
*    Days        - Julian Day of the year         0.0D0  .. 366.0D0
*
*  OutPuts       :
*    Mon         - Month                          1 .. 12
*    Day         - Day                            1 .. 28,29,30,31
*    Hr          - Hour                           0 .. 23
*    minute         - Minute                         0 .. 59
*    Sec         - Second                         0.0D0 .. 59.999D0
*
*  Locals        :
*    DayofYr     - Day of year
*    Temp        - Temporary REAL*8 values
*    IntTemp     - Temporary INTEGER value
*    i           - Index
*    LMonth[12]  - INTEGER Array containing the Number of days per month
*
*  Coupling      :
*    None.
* -----------------------------------------------------------------------------

      SUBROUTINE DAYS2MDHMS  ( Year,Days,  Mon,Day,Hr,minute,Sec )
        IMPLICIT NONE
        REAL*8 Days,Sec
        INTEGER Year, Mon, Day, Hr, minute
* ----------------------------  Locals  -------------------------------
        INTEGER IntTemp,i,DayofYr, LMonth(12)
        REAL*8 Temp

        ! --------------------  Implementation   ----------------------
        ! -------------- Set up array of days in month  ---------------
        DO i = 1,12
            LMonth(i) = 31
          ENDDO
        LMonth( 2) = 28
        LMonth( 4) = 30
        LMonth( 6) = 30
        LMonth( 9) = 30
        LMonth(11) = 30

        DayofYr= IDINT(Days )

        ! ---------------- Find month and Day of month ----------------
        IF (MOD(Year,4).eq.0) THEN
            LMonth(2)= 29
          ENDIF
        i= 1
        IntTemp= 0
        DO WHILE ( (DayofYr.gt.IntTemp + LMonth(i) ) .and. ( i.lt.12 ))
            IntTemp= IntTemp + LMonth(i)
            i= i+1
          ENDDO
        Mon= i
        Day= DayofYr - IntTemp

        ! ---------------- Find hours Minutes and seconds -------------
        Temp= (Days - DayofYr )*24.0D0
        Hr  = IDINT( Temp )
        Temp= (Temp-Hr) * 60.0D0
        minute = IDINT( Temp )
        Sec = (Temp-minute) * 60.0D0

      RETURN
      END
      
* -----------------------------------------------------------------------------
*
*                           SUBROUTINE INVJDay
*
*  this subroutine finds the Year, month, day, hour, Minute and second
*  given the Julian date. TU can be UT1, TDT, TDB, etc.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    jd          - julian date (days only)           days from 4713 bc
*    jdFrac      - julian date (fraction of a day)   days from 0 hr of the day
*
*  OutPuts       :
*    Year        - Year                           1900 .. 2100
*    Mon         - Month                          1 .. 12
*    Day         - Day                            1 .. 28,29,30,31
*    Hr          - Hour                           0 .. 23
*    minute         - Minute                         0 .. 59
*    Sec         - Second                         0.0D0 .. 59.999D0
*
*  Locals        :
*    Days        - Day of year plus fractional
*                  portion of a day               days
*    Tu          - Julian Centuries from 0 h
*                  Jan 0, 1900
*    Temp        - Temporary real values
*    LeapYrs     - Number of Leap years from 1900
*
*  Coupling      :
*    DAYS2MDHMS  - Finds MD HMS given Days and Year
*
*  References    :
*    Vallado       2007, 208, Alg 22, Ex 3-13
* -----------------------------------------------------------------------------

      SUBROUTINE INVJDay    ( jd, jdFrac, Year,Mon,Day,Hr,minute, Sec )
        IMPLICIT NONE
        INTEGER Year, Mon, Day, Hr, minute
        REAL*8  Sec, jd, jdFrac
* ----------------------------  Locals  -------------------------------
        INTEGER LeapYrs
        REAL*8  Days, Tu, Temp, dt, dtt

        ! --------------------  Implementation   ----------------------
		! check jdfrac for multiple days
		if (ABS(jdFrac) .ge. 1.0D0) THEN
			dtt = IDINT(jdFrac);
			jd = jd + dtt
			jdFrac = jdFrac - dtt
		ENDIF

		! check for fraction of a day included in the jd
		dt = jd - IDINT(jd) - 0.5
        dt = jd - real(IDINT(jd)) - 0.5d0  !! needed real??
		if (ABS(dt) .gt. 0.00000001D0) THEN
			jdFrac = jdFrac + dt
			jd = jd - dt;
        ENDIF 
        

        ! ---------------- Find Year and Days of the year -------------
        Temp   = jd-2415019.5D0
        Tu     = Temp / 365.25D0
        Year   = 1900 + IDINT( Tu )
        LeapYrs= IDINT( ( Year-1901 )*0.25D0 )
        Days   = IDINT(Temp - ((Year-1900)*365.0D0 + LeapYrs ))

        ! -------------- Check for case of beginning of a year --------
        IF ( Days + jdFrac .lt. 1.0D0 ) THEN
            Year   = Year - 1
            LeapYrs= IDINT( ( Year-1901 )*0.25D0 )
            Days   = IDINT(Temp - ((Year-1900)*365.0D0 + LeapYrs ))
          ENDIF

        ! ------------------ Find remaing data  -----------------------
        CALL DAYS2MDHMS( Year,Days + jdFrac, Mon,Day,Hr,minute,Sec )

      RETURN
      END

* -----------------------------------------------------------------------------
*
*                           SUBROUTINE FINDDAYS
*
*  this subroutine finds the fractional days through a year given the year,
*    month, day, hour, Minute and second.
*
*  Algorithm     : Set up array for the Number of days per month
*                  Find Leap Year - be sure to account for the 400 years 
*                  Check for a leap year
*                  Loop to find the elapsed days in the year
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Year        - Year                           1900 .. 2100
*    Mon         - Month                          1 .. 12
*    Day         - Day                            1 .. 28,29,30,31
*    Hr          - Hour                           0 .. 23
*    minute         - Minute                         0 .. 59
*    Sec         - Second                         0.0D0 .. 59.999D0
*
*  OutPuts       :
*    Days        - Day of year plus fraction of a
*                    day                          days
*
*  Locals        :
*    LMonth      - Length of months of year
*    i           - Index
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 207, Ex 3-12
*
* -----------------------------------------------------------------------------

      SUBROUTINE FINDDAYS    ( Year,Month,Day,Hr,minute, Sec,  Days )
        IMPLICIT NONE
        INTEGER Year, Month, Day, Hr, minute
        REAL*8 Sec, Days
* ----------------------------  Locals  -------------------------------
        INTEGER i, LMonth(12)

        ! --------------------  Implementation   ----------------------
        DO i = 1,12
            LMonth(i) = 31
          ENDDO
        LMonth( 2) = 28
        LMonth( 4) = 30
        LMonth( 6) = 30
        LMonth( 9) = 30
        LMonth(11) = 30
        IF (MOD(Year,4).eq.0) THEN
            LMonth(2)= 29
          ENDIF

        i   = 1
        Days= 0.0D0
        DO WHILE ((i .lt. Month) .and. ( i .lt. 12 ))
            Days= Days + LMonth(i)
            i= i + 1
          ENDDO

        Days= Days + Day + Hr/24.0D0 + minute/1440.0D0 + Sec/86400.0D0
      RETURN
      END


      