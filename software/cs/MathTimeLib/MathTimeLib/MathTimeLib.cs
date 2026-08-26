// ----------------------------------------------------------------------------
//
//                                   MathTimeLib.cs
//
// this library contains various mathematical routines.
//
//                               companion code for
//                  fundamentals of astrodynamics and applications
//                                      2022
//                                by david vallado
//
//                  email dvallado@comspoc.com, davallado@gmail.com
//
//   current :
//              20 jan 2025  david vallado
//                           original baseline
//                           
//   changes :
//   
//   defines
//     infinite     - infinite value                                 999999.9
//     undefined    - undefined value                                999999.1  
//     edirection   - direction for operation in a method            efrom, eto
//
// ----------------------------------------------------------------------------      

using System;
using System.Globalization;
using System.Linq;
using System.Numerics;  


namespace MathTimeMethods
{

    // -----------------------------------------------------------------------------
    public class MathTimeLib
    {
        public string MathTimeLibVersion = "MathTimeLib Version 2025-01-20";

        // if needed...
        public Complex ComplexLibr = new Complex();

        // ------------------------------- constants -----------------------------------
        public static class globals
        {
            public static double infinite = 999999.9;
            public static double undefined = 999999.1;
        }

        public enum Edirection { efrom, eto };


        // ==============================================================================
        //                              Trigonometric routines  
        // ==============================================================================

        // -------------------------------------------------------------------------------
        //
        //                               function cot
        //
        //  this function finds the cotangent of an input in radians.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    xval        - input to take cotangent of        rad
        //
        //  outputs       :
        //    cot         - result
        //
        //  locals        :
        //    temp        - temporary real variable
        // ----------------------------------------------------------------------------

        public double cot (double xval)
        {
            double temp;

            temp = Math.Tan(xval);
            if (Math.Abs(temp) < 0.00000001)
                return 999999.9;
            else
                return 1.0 / temp;
        }  //  cot


        // ------------------------------------------------------------------------------
        //
        //                               function acosh
        //
        //  this function evaluates the inverse hyperbolic cosine function.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    xval        - angle value                                  1.0 to infinity
        //
        //  outputs       :
        //    acosh       - result                                       any real
        //
        //  locals        :
        //    temp        - temporary value
        //
        //  coupling      :
        //    none.
        // -----------------------------------------------------------------------------

        public double acosh (double xval)
        {
            return Math.Log(xval + Math.Sqrt(xval * xval - 1.0));
        }  // end acosh

        // ------------------------------------------------------------------------------
        //
        //                               function asinh
        //
        //  this function evaluates the inverse hyperbolic sine function.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    xval        - angle value                                  any real
        //
        //  outputs       :
        //    asinh       - result                                       any real
        //
        //  locals        :
        //    none.
        //
        //  coupling      :
        //    none.
        // -----------------------------------------------------------------------------

        public double asinh (double xval)
        {
            return Math.Log(xval + Math.Sqrt(xval * xval + 1.0));
        }  //  asinh


        // ==============================================================================
        //                                 Vector routines 
        // ==============================================================================

        // ------------------------------------------------------------------------------
        //
        //                               function dot
        //
        //  this function finds the dot product of two vectors.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    vec1        - vector number 1
        //    vec2        - vector number 2
        //
        //  outputs       :
        //    dot         - result
        //
        //  locals        :
        //    none.
        // -----------------------------------------------------------------------------

        public double dot (double[] x, double[] y)
        {
            return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
        }  //  dot


        // ------------------------------------------------------------------------------
        //
        //                               function mag
        //
        //  this procedure finds the magnitude of a vector.  
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    vec         - vector
        //
        //  outputs       :
        //    vec         - answer stored in fourth component
        //
        //  locals        :
        //    none.
        //
        //  coupling      :
        //    none.
        // -----------------------------------------------------------------------------

        public double mag (double[] x)
        {
            return Math.Sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
        }  // mag



        // ------------------------------------------------------------------------------
        //
        //                               procedure cross
        //
        //  this procedure crosses two vectors.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    vec1        - vector number 1
        //    vec2        - vector number 2
        //
        //  outputs       :
        //    outvec      - vector result of a x b
        //
        //  locals        :
        //    none.
        //
        //  coupling      :
        //    none
        // ----------------------------------------------------------------------------

        public void cross (double[] vec1, double[] vec2, out double[] outvec)
        {
            double[] tempvec = new double[3];
            tempvec[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
            tempvec[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
            tempvec[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
            outvec = tempvec;
        }  //  cross


        // ------------------------------------------------------------------------------
        //
        //                               procedure norm
        //
        //  this procedure calculates a unit vector given the original vector.  if a
        //  zero vector is input, the vector is set to zero.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    vec        - vector
        //
        // outputs       :
        //    outvec     - unit vector
        //
        //  locals       :
        //    i          - index
        //    small      - value defining a small value
        //    magv       - magnitude of the vector
        //
        //  coupling     :
        //    mag        magnitude of a vector
        // ---------------------------------------------------------------------------*/

        public double[] norm (double[] vec1)
        {
            double[] norm = new double[3];
            double magr;

            magr = mag(vec1);
            norm[0] = vec1[0] / magr;
            norm[1] = vec1[1] / magr;
            norm[2] = vec1[2] / magr;
            return norm;
        }  // norm 


        // ------------------------------------------------------------------------------
        //
        //                               procedure angle
        //
        //  this procedure calculates the angle between two vectors.  the output is
        //    set to 999999.1 to indicate an undefined value.  be sure to check for
        //    this at the output phase.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    vec1        - vector number 1
        //    vec2        - vector number 2
        //
        //  outputs       :
        //    theta       - angle between the two vectors  -Math.PI to Math.PI
        //    magv1       - magnitude of vec1
        //    magv2       - magnitude of vec2
        //    small       - value defining a small value
        //    undefined   - large number to use in place of a not defined number 
        //
        //  locals        :
        //    temp        - temporary real variable
        //
        //  coupling      :
        //    dot           dot product of two vectors
        // -----------------------------------------------------------------------------

        public double angle (double[] vec1, double[] vec2)
        {
            double small, undefined, magv1, magv2, temp;
            small = 0.00000001;
            undefined = 999999.1;

            magv1 = mag(vec1);
            magv2 = mag(vec2);

            if (magv1 * magv2 > small * small)
            {
                temp = dot(vec1, vec2) / (magv1 * magv2);
                if (Math.Abs(temp) > 1.0)
                    temp = Math.Sign(temp) * 1.0;
                return Math.Acos(temp);
            }
            else
                return undefined;
        }  //  angle


        // ------------------------------------------------------------------------------
        //
        //                               procedure roti
        //
        //  this procedure performs a rotation about the ith axis. i is specified
        //    for each operation.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    vec         - input vector
        //    xval        - angle of rotation              rad
        //
        //  outputs       :
        //    outvec      - vector result
        //
        //  locals        :
        //    c           - cosine of the angle xval
        //    s           - sine of the angle xval
        //    temp        - temporary extended value
        //
        //  coupling      :
        //    none.
        //
        // -----------------------------------------------------------------------------

        public double[] rot1 (double[] vec, double xval)
        {
            double c, s, temp;
            double[] outvec = new double[3];

            temp = vec[2];
            c = Math.Cos(xval);
            s = Math.Sin(xval);

            outvec[2] = c * vec[2] - s * vec[1];
            outvec[1] = c * vec[1] + s * temp;
            outvec[0] = vec[0];

            return outvec;
        }  // rot1  


        public double[] rot2(double[] vec, double xval)
        {
            double c, s, temp;
            double[] outvec = new double[3];

            temp = vec[2];
            c = Math.Cos(xval);
            s = Math.Sin(xval);

            outvec[2] = c * vec[2] + s * vec[0];
            outvec[0] = c * vec[0] - s * temp;
            outvec[1] = vec[1];

            return outvec;
        }   // rot2  


        public double[] rot3(double[] vec, double xval)
        {
            double c, s, temp;
            double[] outvec = new double[3];

            temp = vec[1];
            c = Math.Cos(xval);
            s = Math.Sin(xval);

            outvec[1] = c * vec[1] - s * vec[0];
            outvec[0] = c * vec[0] + s * temp;
            outvec[2] = vec[2];

            return outvec;
        }  // rot3 


        // ------------------------------------------------------------------------------
        //
        //                                      rotimat
        //
        //  this function sets up a rotation matrix for an input angle about the first
        //    axis.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    xval        - angle of rotation              rad
        //
        //  outputs       :
        //    outmat      - matrix result
        //
        //  locals        :
        //    c           - cosine of the angle xval
        //    s           - sine of the angle xval
        //
        //  coupling      :
        //    none.
        // -----------------------------------------------------------------------------

        public double[,] rot1mat(double xval)
        {
            double c, s;
            double[,] outmat = new double[3, 3];

            c = Math.Cos(xval);
            s = Math.Sin(xval);

            outmat[0, 0] = 1.0;
            outmat[0, 1] = 0.0;
            outmat[0, 2] = 0.0;

            outmat[1, 0] = 0.0;
            outmat[1, 1] = c;
            outmat[1, 2] = s;

            outmat[2, 0] = 0.0;
            outmat[2, 1] = -s;
            outmat[2, 2] = c;
            return outmat;
        }  //  rot1mat


        public double[,] rot2mat(double xval)
        {
            double c, s;
            double[,] outmat = new double[3, 3];

            c = Math.Cos(xval);
            s = Math.Sin(xval);

            outmat[0, 0] = c;
            outmat[0, 1] = 0.0;
            outmat[0, 2] = -s;

            outmat[1, 0] = 0.0;
            outmat[1, 1] = 1.0;
            outmat[1, 2] = 0.0;

            outmat[2, 0] = s;
            outmat[2, 1] = 0.0;
            outmat[2, 2] = c;
            return outmat;
        }  // rot2mat


        public double[,] rot3mat(double xval)
        {
            double c, s;
            double[,] outmat = new double[3, 3];

            c = Math.Cos(xval);
            s = Math.Sin(xval);

            outmat[0, 0] = c;
            outmat[0, 1] = s;
            outmat[0, 2] = 0.0;

            outmat[1, 0] = -s;
            outmat[1, 1] = c;
            outmat[1, 2] = 0.0;

            outmat[2, 0] = 0.0;
            outmat[2, 1] = 0.0;
            outmat[2, 2] = 1.0;
            return outmat;
        }  // rot3mat


        // ------------------------------------------------------------------------------
        //
        //                               procedure addvec
        //
        //  this procedure adds two vectors possibly multiplied by a constant.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    a1          - constant multiplier
        //    a2          - constant multiplier
        //    vec1        - vector number 1
        //    vec2        - vector number 2
        //
        //  outputs       :
        //    outvec      - vector result of a + b
        //
        //  locals        :
        //    row         - index
        //
        //  coupling      :
        //     none
        // -----------------------------------------------------------------------------

        public void addvec(double a1, double[] vec1, double a2, double[] vec2, out double[] vec3)
        {
            vec3 = new double[] { 0.0, 0.0, 0.0 };
            int row;
            double[] tempvec = new double[3];

            for (row = 0; row <= 2; row++)
            {
                vec3[row] = 0.0;
                vec3[row] = a1 * vec1[row] + a2 * vec2[row];
            }
        } // addvec


        // ------------------------------------------------------------------------------
        //
        //                               procedure addvec3
        //
        //  this procedure adds three vectors possibly multiplied by a constant.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    a1          - constant multiplier
        //    a2          - constant multiplier
        //    a3          - constant multiplier
        //    vec1        - vector number 1
        //    vec2        - vector number 2
        //    vec3        - vector number 3
        //
        //  outputs       :
        //    outvec      - vector result of a + b + c
        //
        //  locals        :
        //    row         - index
        //
        //  coupling      :
        //     none
        // -----------------------------------------------------------------------------

        public void addvec3(double a1, double[] vec1, double a2, double[] vec2, double a3, double[] vec3, out double[] vec4)
        {
            vec4 = new double[] { 0.0, 0.0, 0.0 };
            int row;
            double[] tempvec = new double[3];

            for (row = 0; row <= 2; row++)
            {
                vec4[row] = 0.0;
                vec4[row] = a1 * vec1[row] + a2 * vec2[row] + a3 * vec3[row];
            }
        } // addvec3

        // ------------------------------------------------------------------------------
        //
        //                               procedure matvecmult
        //
        //  this procedure multiplies a 3x3 matrix and a 3x1 vector together.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    mat         - order x order matrix
        //    vec         - vector
        //    order       - square size of the mat and the size of the vector which must be the same
        //
        //  outputs       :
        //    vecout      - vector result of mat * vec
        //
        //  locals        :
        //    row         - row index
        //    col         - column index
        //    ktr         - index
        //
        //  coupling      :
        // -----------------------------------------------------------------------------

        public double[] matvecmult(double[,] mat, double[] vec, int order)
        {
            int row, ktr;
            double[] vecout = new double[order];
            for (row = 0; row < order; row++)
            {
                vecout[row] = 0.0;
                for (ktr = 0; ktr < order; ktr++)
                    vecout[row] = vecout[row] + mat[row, ktr] * vec[ktr];
            }
            return vecout;
        } // matvecmult  



        // -----------------------------------------------------------------------------
        // perform vector outer product
        // makes a order x order matrix from multiplying 2 [order] dimensioned vectors

        public double[,] vecouter(double[] vec1, double[] vec2, int order)
        {
            double[,] mat = new double[order, order];
            int row, col;
            for (row = 0; row <= order - 1; row++)
            {
                for (col = 0; col <= order - 1; col++)
                {
                    mat[row, col] = vec1[row] * vec2[col];
                }
            }
            return mat;
        } // vecouter  


        // ==============================================================================
        //                                   Matrix routines 
        // ==============================================================================

        // -----------------------------------------------------------------------------
        // form a copy of a square matrix

        public double[,] matequal(double[,] mat1, int matr)
        {
            int row, col;
            double[,] mat2 = new double[matr, matr];
            for (row = 0; row < matr; row++)
            {
                for (col = 0; col < matr; col++)
                {
                    mat2[row, col] = mat1[row, col];
                }
            }
            return mat2;
        }  // matequal 


        // -----------------------------------------------------------------------------
        // multiply a matrix * a scalar

        public double[,] matscale(double[,] mat1, int matr, int matc, double scale)
        {
            int row, col;
            double[,] mat2 = new double[matr, matr];
            for (row = 0; row < matr; row++)
            {
                for (col = 0; col < matc; col++)
                {
                    mat2[row, col] = mat1[row, col] * scale;
                }
            }
            return mat2;
        }  // matscale 


        // -----------------------------------------------------------------------------
        // matrixes must be the same size

        public double[,] matadd(double[,] mat1, double[,] mat2, int mat1r, int mat1c)
        {
            int row, col;
            double[,] mat3 = new double[mat1r, mat1c];
            for (row = 0; row < mat1r; row++)
            {
                for (col = 0; col < mat1c; col++)
                    mat3[row, col] = mat1[row, col] + mat2[row, col];
            }
            return mat3;
        } // matadd  


        // -----------------------------------------------------------------------------
        // matrixes must be the same size

        public double[,] matsub(double[,] mat1, double[,] mat2, int mat1r, int mat1c)
        {
            int row, col;
            double[,] mat3 = new double[mat1r, mat1c];
            for (row = 0; row < mat1r; row++)
            {
                for (col = 0; col < mat1c; col++)
                    mat3[row, col] = mat1[row, col] - mat2[row, col];
            }
            return mat3;
        } // matsub  


        // ------------------------------------------------------------------------------
        //
        //                               procedure matmult
        //
        //  this procedure multiplies two matricies up to 10x10 together.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    mat1        - matrix number 1
        //    mat2        - matrix number 2
        //    mat1r       - matrix number 1 rows
        //    mat1c       - matrix number 1 columns
        //    mat2c       - matrix number 2 columns
        //
        //  outputs       :
        //    mat3        - matrix result of mat1 * mat2 of size mat1r x mat2c
        //
        //  locals        :
        //    row         - row index
        //    col         - column index
        //    ktr         - index
        //
        //  coupling      :
        // -----------------------------------------------------------------------------

        public double[,] matmult(double[,] mat1, double[,] mat2, int mat1r, int mat1c, int mat2c)
        {
            int row, col, ktr;
            double[,] mat3 = new double[mat1r, mat2c];
            for (row = 0; row < mat1r; row++)
            {
                for (col = 0; col < mat2c; col++)
                {
                    mat3[row, col] = 0.0;
                    for (ktr = 0; ktr < mat1c; ktr++)
                        mat3[row, col] = mat3[row, col] + mat1[row, ktr] * mat2[ktr, col];
                }
            }
            return mat3;
        } // matmult  


        // -----------------------------------------------------------------------------
        // form the transpose of a square matrix

        public double[,] mattrans(double[,] mat1, int matr)
        {
            int row, col;
            double[,] mat2 = new double[matr, matr];
            for (row = 0; row < matr; row++)
            {
                for (col = 0; col < matr; col++)
                {
                    mat2[row, col] = mat1[col, row];
                }
            }
            return mat2;
        }  // mattrans 



        // -----------------------------------------------------------------------------
        // form the transponse of a non-square matrix

        public double[,] mattransx(double[,] mat1, int matr, int matc)
        {
            int row, col;
            double[,] mat2 = new double[matc, matr];
            for (row = 0; row < matr; row++)
            {
                for (col = 0; col < matc; col++)
                {
                    mat2[col, row] = mat1[row, col];
                }
            }
            return mat2;
        }  // mattransx


        // -------------------------------------------------------------------------------
        //
        //                               procedure ludecomp
        //
        //  this procedure decomposes a matrix into an lu form.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    order       - order of matrix
        //
        //  outputs       :
        //    lu          - lu decomposition matrix
        //    index       - index vector for pivoting
        //
        //  locals        :
        //    i           - index
        //    j           - index
        //    k           - index
        //    imax        - pivot row pointer
        //    scale       - scale factor vector
        //    sum         - temporary variables
        //    amax        - temporary variables
        //    dum         - temporary variables
        //
        //  coupling      :
        //    none
        //
        //  references    :
        //    numerical recipes - flannery
        // -----------------------------------------------------------------------------

        private void ludecomp(ref double[,] lu, out int[] indexx, int order)
        {
            const double small = 0.000001;
            int i, j, k, imax;
            double[] scale = new double[order];
            double sum, amax, dum;

            indexx = new int[] { 0, 0, 0, 0, 0, 0 };

            imax = 0;
            for (i = 0; i < order; i++)
            {
                amax = 0.0;
                for (j = 0; j < order; j++)
                    if (Math.Abs(lu[i, j]) > amax)
                        amax = Math.Abs(lu[i, j]);
                if (Math.Abs(amax) < small)
                {
                    // " singular matrix amax "
                }
                scale[i] = 1.0 / amax;
            }
            for (j = 0; j < order; j++)
            {
                for (i = 0; i <= j - 1; i++)
                {
                    sum = lu[i, j];
                    for (k = 0; k <= i - 1; k++)
                        sum = sum - lu[i, k] * lu[k, j];
                    lu[i, j] = sum;
                }
                amax = 0.0;
                for (i = j; i < order; i++)
                {
                    sum = lu[i, j];
                    for (k = 0; k <= j - 1; k++)
                        sum = sum - lu[i, k] * lu[k, j];
                    lu[i, j] = sum;
                    dum = scale[i] * Math.Abs(sum);
                    if (dum >= amax)
                    {
                        imax = i;
                        amax = dum;
                    }
                }
                if (j != imax)
                {
                    for (k = 0; k < order; k++)
                    {
                        dum = lu[imax, k];
                        lu[imax, k] = lu[j, k];
                        lu[j, k] = dum;
                    }
                    scale[imax] = scale[j];
                }
                indexx[j] = imax;
                if (Math.Abs(lu[j, j]) < small)
                {
                    // " matrix is singular lu "
                }
                if (j != order - 1)
                {
                    dum = 1.0 / lu[j, j];
                    for (i = j + 1; i < order; i++)
                        lu[i, j] = dum * lu[i, j];
                }
            }
        } // ludecmp



        // -------------------------------------------------------------------------------
        //
        //                               procedure lubksub
        //
        //  this procedure finds the inverse of a matrix using lu decomposition.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    order       - order of matrix
        //    lu          - lu decomposition matrix
        //    index       - index vector for pivoting
        //
        //  outputs       :
        //    b           - solution vector
        //
        //  locals        :
        //    i           - index
        //    j           - index
        //    i0          - pointer to non-zero element
        //    iptr        - pivot rwo pointer
        //    sum         - temporary variables
        //
        //  coupling      :
        //    none
        //
        //  references    :
        //    numerical recipes - flannery
        // -----------------------------------------------------------------------------

        private void lubksub(double[,] lu, int[] indexx, int order, ref double[] b)
        {
            int i, j, iptr, ii;
            double sum;

            ii = -1;
            for (i = 0; i < order; i++)
            {
                iptr = indexx[i];
                sum = b[iptr];
                b[iptr] = b[i];
                if (ii != -1)
                    for (j = ii; j <= i - 1; j++)
                        sum = sum - lu[i, j] * b[j];
                else
                {
                    if (sum != 0.0)
                        ii = i;
                }
                b[i] = sum;
            }

            //b[order - 1] = b[order - 1] / lu[order - 1, order - 1];

            for (i = order - 1; i >= 0; i--)
            {
                sum = b[i];
                for (j = i + 1; j < order; j++)
                    sum = sum - lu[i, j] * b[j];
                b[i] = sum / lu[i, i];
            }
        }  // lubksub


        // -------------------------------------------------------------------------------
        //
        //                               procedure matinverse
        //
        //  this procedure finds the inverse of a matrix using lu decomposition.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    mat         - matrix to invert, 0 array starts
        //    order       - order of matrix
        //
        //  outputs       :
        //    matinv      - inverted matrix, 0 array starts
        //
        //  locals        :
        //    i           - index
        //    j           - index
        //    index       - index vector for pivoting
        //    lu          - lu decomposition matrix
        //    b           - operational vector to form matinv
        //
        //  coupling      :
        //    ludecomp    -
        //    lubksub     -
        //
        //  references    :
        //    numerical recipes - flannery
        // -----------------------------------------------------------------------------

        public void matinverse(double[,] mat, int order, out double[,] matinv)
        {
            int i, j;
            int[] indexx = new int[order];
            double[,] lu = new double[order, order];
            double[] b = new double[order];
            matinv = new double[,] { { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 },
               { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 } };

            for (i = 0; i < order; i++)
            {
                indexx[i] = i;
                for (j = 0; j < order; j++)
                    lu[i, j] = mat[i, j];
            }

            ludecomp(ref lu, out indexx, order);

            for (j = 0; j < order; j++)
            {
                for (i = 0; i < order; i++)
                {
                    if (i == j)
                        b[i] = 1.0;
                    else
                        b[i] = 0.0;
                }

                lubksub(lu, indexx, order, ref b);

                for (i = 0; i < order; i++)
                    matinv[i, j] = b[i];
            }
        } // matinverse



        // ----------------------------------------------------------------------------- 
        //
        //                               procedure determinant
        //
        //  This function calculates the determinant value using L - U decompisition.
        //    The formula must have a non-zero number in the 0, 0 position. if the
        //    function senses a non-zero number in row 0, it exchanges row0 for a row
        //    with a non-zero number.
        //    has trouble with this??? 
        //    double[,] mat =  { { 3, 0, 2.0},
        //                           { 2, 0, -2},
        //                           { 0, 1.0, 1} };
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    mat1        - matrix to find determinant of
        //    order       - order of matrix
        //
        //  outputs       :
        //    determinant - determinant of mat1
        //
        //  locals        :
        //    i, j, k, n  - index
        //    index       - index vector for pivoting
        //    lu          - lu decomposition matrix
        //    b           - operational vector to form matinv
        //
        //  coupling      :
        //    ludecomp    -
        //    lubksub     -
        //
        //  references    :
        //    Marion        pg. 168 - 172, 126 - 127
        // -------------------------------------------------------------------------- - */

        public double determinant(double[,] mat1, int order)
        {
            double small = 0.00000001;
            int i, j, k;
            double temp, d, sum;
            double[,] l = new double[order, order];
            double[,] u = new double[order, order];

            sum = 0.0;
            // ----------- Switch a non zero row to the first row---------- 
            if (Math.Abs(mat1[0, 0]) < small)
            {
                j = 0;
                while (j < order)
                {
                    if (Math.Abs(mat1[j, 0]) > small)
                    {
                        for (k = 0; k < order; k++)
                        {
                            temp = mat1[0, k];
                            mat1[0, k] = mat1[j, k];
                            mat1[j, k] = temp;
                        }
                        j = order + 1;
                    }
                    j = j + 1;
                } // if abs
            }

            for (i = 0; i < order; i++)
                l[i, 0] = mat1[i, 0];
            for (j = 0; j < order; j++)
                u[0, j] = mat1[0, j] / l[0, 0];
            for (j = 1; j < order; j++)
            {
                for (i = j; i < order; i++)
                {
                    sum = 0.0;
                    for (k = 0; k <= j - 1; k++)
                        sum = sum + l[i, k] * u[k, j];
                    l[i, j] = mat1[i, j] - sum;
                } // for i 
                u[j, j] = 1.0;
                if (j != order)
                {
                    for (i = j + 1; i < order; i++)
                    {
                        sum = 0.0;
                        for (k = 0; k <= j - 1; k++)
                            sum = sum + l[j, k] * u[k, i];
                        //if (Math.Abs(l[j, j]) < small)
                        //    u[j, i] = 1.0;
                        //else
                            u[j, i] = (mat1[j, i] - sum) / l[j, j];
                    } // for i 
                } // if j 
            } //  for j 
            d = 1.0;
            for (i = 0; i < order; i++)
                d = d * l[i, i];
            return d;
        } // determinant


        // -----------------------------------------------------------------------------
        // https://www.geeksforgeeks.org/adjoint-inverse-matrix/
        // Function to get cofactor of A[p,q] in [,]temp. n is current
        // dimension of [,]mat

        public void getCofactor(double[,] mat, out double[,] temp, int p, int q, int n)
        {
            temp = new double[,] { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
            int i = 0, j = 0;

            // Looping for each element of the matrix
            for (int row = 0; row < n; row++)
            {
                for (int col = 0; col < n; col++)
                {
                    // Copying into temporary matrix only those element
                    // which are not in given row and column
                    if (row != p && col != q)
                    {
                        temp[i, j++] = mat[row, col];

                        // Row is filled, so increase row index and
                        // reset col index
                        if (j == n - 1)
                        {
                            j = 0;
                            i++;
                        }
                    }
                }
            }
        }  // getCofactor


        // -----------------------------------------------------------------------------
        // https://www.geeksforgeeks.org/adjoint-inverse-matrix/
        // function to get adjoint of 3x3 mat[order, order] in adj[order, order].

        public void adjoint(double[,] mat, out double[,] adj)
        {
            int order = 3;
            adj = new double[,] { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
            if (order == 1)
            {
                adj[0, 0] = 1;
                return;
            }

            // temp is used to store cofactors of [,]mat
            int sign = 1;
            double[,] temp = new double[order, order];

            for (int i = 0; i < order; i++)
            {
                for (int j = 0; j < order; j++)
                {
                    // Get cofactor of mat[i,j]
                    getCofactor(mat, out temp, i, j, order);

                    // sign of adj[j,i] positive if sum of row
                    // and column indexes is even.
                    sign = ((i + j) % 2 == 0) ? 1 : -1;

                    // Interchanging rows and columns to get the
                    // transpose of the cofactor matrix
                    adj[j, i] = sign * determ(temp, order - 1);
                }
            }
        }   // adjoint


        // -----------------------------------------------------------------------------
        // https://www.geeksforgeeks.org/adjoint-inverse-matrix/
        // Recursive function for finding determinant of matrix.
        //   n is current dimension of [,]mat. 

        public double determ(double[,] mat, int n)
        {
            int order = 3;
            double D = 0; // Initialize result

            // Base case : if matrix contains single element
            if (n == 1)
                return mat[0, 0];

            // store cofactors
            double[,] temp = new double[order, order];

            // store sign multiplier
            int sign = 1; 

            // Iterate for each element of first row
            for (int f = 0; f < n; f++)
            {
                // Getting Cofactor of mat[0,f]
                getCofactor(mat, out temp, 0, f, n);
                D = D + sign * mat[0, f] * determ(temp, n - 1);

                // terms are to be added with alternate sign
                sign = -sign;
            }
            return D;
        }  // determ


        // -------------------------------------------------------------------------------
        //
        //                               procedure mat33inverse
        //
        //  this procedure finds the inverse of a 3x3 matrix using determinants.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    mat         - matrix to invert, 0 array starts
        //
        //  outputs       :
        //    matinv      - inverted matrix, 0 array starts
        //
        //  locals        :
        //    i           - index
        //    j           - index
        //
        //  coupling      :
        //     determinant- determinant of a matrix
        //     adjoint    - adjoint of a matrix
        //     
        //  references    :
        //    https://www.geeksforgeeks.org/adjoint-inverse-matrix/
        // -----------------------------------------------------------------------------

        public bool mat33inverse(double[,] mat, out double[,] matinv)
        {
            double detval, temp;
            int order = 3;
            // find determinant of mat
            detval = determ(mat, order);
            matinv = new double[,] { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };

            if (detval == 0)
            {
                Console.Write("Singular matrix, can't find its inverse");
                return false;
            }

            // find adjoint
            double[,] adj = new double[order, order];
            adjoint(mat, out adj);

            // find Inverse using formula inverse(mat) = adj(mat)/det(mat)
            temp = 1.0 / detval;
            for (int i = 0; i < order; i++)
                for (int j = 0; j < order; j++)
                    matinv[i, j] = adj[i, j] * temp;
            
            return true;
        }    // mat33inverse

        // -----------------------------------------------------------------------------

        public void writemat(string matname, double[,] mat, int row, int col, out string outstr)
        {
            string fmt = "+#.#########;-#.#########";
            int r, c;

            outstr = "matrix " + matname + "\n";
            for (r = 0; r < row; r++)
            {
                for (c = 0; c < col; c++)
                    outstr = outstr + mat[r, c].ToString(fmt).PadLeft(16) + " ";
                outstr = outstr + "\n";
            }
        }  // writemat


        // -----------------------------------------------------------------------------

        public void writeexpmat(string matname, double[,] mat, int row, int col, out string outstr)
        {
            // format strings to show signs "and" to not round off if trailing 0!!
            string fmt = "+#.#########0E+00;-#.#########0E+00";
            int r, c;

            outstr = "matrix " + matname + "\n";
            for (r = 0; r < row; r++)
            {
                for (c = 0; c < col; c++)
                    outstr = outstr + mat[r, c].ToString(fmt).PadLeft(16) + " ";
                outstr = outstr + "\n";
            }
        }  // writeexpmat



        // -----------------------------------------------------------------------------
        //
        // https://rosettacode.org/wiki/Cholesky_decomposition#C.23
        // Content is available under GNU Free Documentation License 1.2 unless otherwise noted.
        // Returns the lower Cholesky Factor, L, of input matrix A. 
        // Satisfies the equation: L*L^T = A.
        // </summary>
        // <param name="a">Input matrix must be square, symmetric, 
        // and positive definite. This method does not check for these properties,
        // and may produce unexpected results of those properties are not met.</param>
        // <returns></returns>
        public double[,] cholesky (double[,] a)
        {
            int n = (int)Math.Sqrt(a.Length);

            double[,] ret = new double[n, n];
            for (int r = 0; r < n; r++)
                for (int c = 0; c <= r; c++)
                {
                    if (c == r)
                    {
                        double sum = 0;
                        for (int j = 0; j < c; j++)
                        {
                            sum = sum + ret[c, j] * ret[c, j];
                        }
                        ret[c, c] = Math.Sqrt(a[c, c] - sum);
                    }
                    else
                    {
                        double sum = 0;
                        for (int j = 0; j < c; j++)
                            sum = sum + ret[r, j] * ret[c, j];
                        ret[r, c] = 1.0 / ret[c, c] * (a[r, c] - sum);
                    }
                }

            return ret;
        }   // cholesky

        // ==============================================================================
        //                                 Polynomial routines 
        // ==============================================================================

        // ------------------------------------------------------------------------------
        //
        //                               function DMulRSub
        //
        // called by factor to find roots of a polynomial.
        //
        // the original code had several places with a double = 0.0
        // The values here seem particularly touchy, values of 1e-45 may not even be enough
        // so implement a try-catch in several places.
        // 
        // 
        // -------------------------------------------------------------------------------

        public void DMulRSub(ref double[] ALPR, ref double[] ALPI, double[] BETR, double[] BETI)
        {
            double Te1, Te2, Te3, Te4, Te5, Te6, Te7, Te8, Te9, Te10, Te11, Te12,
            Te13, Te14, Te15, Te16, tem, DE15, DE16;
            Te8 = 0.0;

            Te1 = ALPR[1] - ALPR[3];
            Te2 = ALPI[1] - ALPI[3];
            Te5 = ALPR[3] - ALPR[2];
            Te6 = ALPI[3] - ALPI[2];
            tem = Te5 * Te5 + Te6 * Te6;

            // ---------------- Check for zero values of tem -------------- }
            try
            {
                Te3 = (Te1 * Te5 + Te2 * Te6) / tem;
                Te4 = (Te2 * Te5 - Te1 * Te6) / tem;
            }
            catch
            {
                Te3 = 0.0; // 9999999.0;   // account for small denominator
                Te4 = 0.0; // 9999999.0;
                Console.WriteLine("te3te4 chk 0.0" + tem.ToString());
            }

            Te7 = Te3 + 1.0;
            Te9 = Te3 * Te3 - Te4 * Te4;
            Te10 = 2.0 * Te3 * Te4;
            DE15 = Te7 * BETR[3] - Te4 * BETI[3];
            DE16 = Te7 * BETI[3] + Te4 * BETR[3];
            Te11 = Te3 * BETR[2] - Te4 * BETI[2] + BETR[1] - DE15;
            Te12 = Te3 * BETI[2] + Te4 * BETR[2] + BETI[1] - DE16;

            Te7 = Te9 - 1.0;
            Te1 = Te9 * BETR[2] - Te10 * BETI[2];
            Te2 = Te9 * BETI[2] + Te10 * BETR[2];
            Te13 = Te1 - BETR[1] - Te7 * BETR[3] + Te10 * BETI[3];
            Te14 = Te2 - BETI[1] - Te7 * BETI[3] - Te10 * BETR[3];
            Te15 = DE15 * Te3 - DE16 * Te4;
            Te16 = DE15 * Te4 + DE16 * Te3;

            Te1 = Te13 * Te13 - Te14 * Te14 - 4.0 * (Te11 * Te15 - Te12 * Te16);
            Te2 = 2.0 * Te13 * Te14 - 4.0 * (Te12 * Te15 + Te11 * Te16);

            //Console.WriteLine(Te1 + " " + Te2 + " " + Te3 + " " + Te4 + " " + Te5 + " " + Te6);
            //Console.WriteLine(Te7 + " " + Te8 + " " + Te9 + " " + Te10 + " " + Te11);
            //Console.WriteLine(Te12 + " " + Te13 + " " + Te14 + " " + Te15 + " " + Te16);

            // ---------------------------------------------------------------------
            //   Sometimes, for stiff systems(the roots vary widely in order
            //   of magnitude), Te1 and Te2 get large enough to have their
            //   squares overflow the floating point range.To prevent this,
            //   when either one is large, they are scaled by 10**10 for the
            //   purpose of finding TeM.  The scale is restored when the
            //   magnitude computation is completed.This should not affect
            //   the accuracy of the computations, since the mantissa is not
            //   affected, only the exponent.
            //   this doesn't appear to be necessary with "modern" computers where e308
            //   is the limit of doubles and the check could move to e150 or so.
            // ----------------------------------------------------------------------

            if ((Te1 > 1.0e15) || (Te2 > 1.0e15))
            {
                Te1 = Te1 * 1.0E-10;
                Te2 = Te2 * 1.0E-10;
                tem = 1.0E10 * Math.Sqrt(Te1 * Te1 + Te2 * Te2);
            }
            else
                tem = Math.Sqrt(Te1 * Te1 + Te2 * Te2);

            if (Te1 > 0.0)
            {
                Te3 = Math.Sqrt(0.5 * (tem + Te1));
                if (Te2 < 0.0)
                    Te3 = -Te3;
                // ----------------- check for zero values of te3 -------------- 
                try
                {
                    Te4 = 0.5 * Te2 / Te3;
                }
                catch
                {
                    Te4 = 0.0; // 9999999.0;   // account for small denominator
                    Console.WriteLine("te3 0.0" + Te3.ToString());
                }
            }
            else
            {
                Te4 = Math.Sqrt(0.5 * (tem - Te1));
                // --------------------- Check for underflows ------------------ 
                try
                {
                    Te3 = 0.5 * Te2 / Te4;
                }
                catch
                {
                    Te3 = 0.0; // 9999999.0;    // account for small denominator
                    Console.WriteLine("te4 0.0" + Te4.ToString());
                }
            }

            Te7 = Te13 + Te3;
            Te8 = Te14 + Te4;
            Te9 = Te13 - Te3;
            Te10 = Te14 - Te4;
            Te1 = 2.0 * Te15;
            Te2 = 2.0 * Te16;

            if (Te7 * Te7 + Te8 * Te8 - Te9 * Te9 - Te10 * Te10 <= 0.0)
            {
                Te7 = Te9;
                Te8 = Te10;
            }
            tem = Te7 * Te7 + Te8 * Te8;

            // ------------- Check for values of almost zero -------------- 
            try
            {
                Te3 = Te1 / tem * Te7 + Te2 / tem * Te8;
                Te4 = Te2 / tem * Te7 - Te1 / tem * Te8;
            }
            catch
            {
                Te3 = 0.0; // 9999999.0;   // account for small denominator
                Te4 = 0.0; // 9999999.0;   
                Console.WriteLine("tem chk 0.0" + tem.ToString());
            }

            ALPR[4] = ALPR[3] + Te3 * Te5 - Te4 * Te6;
            ALPI[4] = ALPI[3] + Te3 * Te6 + Te4 * Te5;
            //Console.WriteLine(Te1 + " " + Te2 + " " + Te3 + " " + Te4 + " " + Te5 + " " + Te6);
            //Console.WriteLine(Te7 + " " + Te8 + " " + Te9 + " " + Te10 + " " + Te11);
            //Console.WriteLine(Te12 + " " + Te13 + " " + Te14 + " " + Te15 + " " + Te16);
            //Console.WriteLine(ALPR[4] + " " + ALPI[4]);
            //Console.WriteLine("end dmlrsub");

        }   // DMulRSub



        // -------------------------------------------------------------------------------
        //
        //                                  factor
        //
        //  This method is a root finding algorithm. It takes the polynomial and
        //    returns the roots (real and imaginary) in the RootS array. it works in "most" 
        //    cases however it misses the correct postive root in some of the angles-only 
        //    cases so don't use it there, use a Halley iteration.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    Poly        - Array of 16 coefficients
        //                        representing the polynomial
        //                        [1] is x^8th, [2] is x^7th, ...
        //                        others are zero
        //    nRoots      - Number of roots
        //
        //  OutPuts       :
        //    Roots       - Array[,] containing roots(real, imag)
        //
        //  Locals        :
        //                    -
        //                    -
        //                    -
        //
        //  Coupling      :
        //    DMulRSub    -
        //
        //  References    :
        //    Original FORTRAN code from USAFA/DFAS MiniTotal program, author unknown
        //    This is Bairstows method?
        //
        // -----------------------------------------------------------------------------

        public void factor(double[] Poly, int nRoots, out double[,] RootS)
        {
            double small = 1.0e-30;
            double[] DPoly = new double[17];
            double[] AlpR = new double[5];
            double[] Alpi = new double[5];
            double[] BetR = new double[5];
            double[] Beti = new double[5];
            bool skip = false;

            int moder, loopcnt, kk, i, j, L, rootcnt;
            double temp1, temp2, AXR, AXi, pmax, tem1, tem2,
            tempreal, tempimag, temp7;

            RootS = new double[,] {  { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 },
                 { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, {0.0, 0.0, 0.0 },
                 { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, {0.0, 0.0, 0.0 },
                 { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, {0.0, 0.0, 0.0 },
                 { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, {0.0, 0.0, 0.0 } };

            temp7 = 0.0;
            pmax = 0.0;
            for (kk = 1; kk <= nRoots + 1; kk++)
                if (Math.Abs(Poly[kk]) > pmax)
                    pmax = Poly[kk];

            if (Math.Abs(pmax) < small)
                pmax = 1.0;

            for (kk = 1; kk <= nRoots + 1; kk++)
                DPoly[kk] = Poly[kk] / pmax;
            //Console.WriteLine("Dpoly");
            //Console.WriteLine(DPoly[1] + " " + DPoly[2] + " " + DPoly[3] + " " + DPoly[4]);
            //Console.WriteLine(DPoly[5] + " " + DPoly[6] + " " + DPoly[7] + " " + DPoly[8]);
            //Console.WriteLine(DPoly[9]);


            if (nRoots > 0)
            {
                rootcnt = 0;
                i = nRoots + 1;

                while ((Math.Abs(DPoly[i]) < small) && (rootcnt != nRoots))
                {
                    rootcnt = rootcnt + 1;
                    RootS[rootcnt, 1] = 0.0;
                    RootS[rootcnt, 2] = 0.0;
                    i = i - 1;
                } //{ While }

                if (rootcnt != nRoots)
                {
                    AXR = 0.8;
                    AXi = 0.0;
                    L = 1;
                    loopcnt = 1;
                    AlpR[1] = AXR;
                    Alpi[1] = AXi;
                    moder = 1;

                    while (rootcnt < nRoots)
                    {
                        BetR[4] = DPoly[1];
                        Beti[4] = 0.0;
                        for (i = 1; i <= nRoots; i++)
                        {
                            j = i + 1;
                            temp1 = BetR[4] * AXR - Beti[4] * AXi;
                            Beti[4] = Beti[4] * AXR + BetR[4] * AXi;
                            BetR[4] = temp1 + DPoly[j];
                        }

                        tempreal = BetR[4];
                        tempimag = Beti[4];

                        if (rootcnt != 0)
                        {
                            for (i = 1; i <= rootcnt; i++)
                            {
                                tem1 = AXR - RootS[i, 1];
                                tem2 = AXi - RootS[i, 2];
                                temp1 = tem1 * tem1 + tem2 * tem2;
                                temp2 = (BetR[4] * tem1 + Beti[4] * tem2) / temp1;
                                Beti[4] = (Beti[4] * tem1 - BetR[4] * tem2) / temp1;
                                BetR[4] = temp2;
                            }
                        }
                        // yes, the moder value gets changed in the middle - as that was in the original FOR
                        // poor programming, but it seems to work. 
                        if (moder == 1)
                        {
                            BetR[1] = BetR[4];
                            Beti[1] = Beti[4];
                            AXR = 0.85;
                            AlpR[2] = AXR;
                            Alpi[2] = AXi;
                            moder = 2;
                            //Console.WriteLine("BetR case 1");
                            //Console.WriteLine(BetR[1] + " " + BetR[2] + " " + BetR[3] + " " + BetR[4]);
                        }
                        else if (moder == 2)
                        {
                            BetR[2] = BetR[4];
                            Beti[2] = Beti[4];
                            AXR = 0.9;
                            AlpR[3] = AXR;
                            Alpi[3] = AXi;
                            moder = 3;
                            //Console.WriteLine("BetR case 2");
                            //Console.WriteLine(BetR[1] + " " + BetR[2] + " " + BetR[3] + " " + BetR[4]);
                        }

                        else if (moder == 5)
                        {
                            BetR[1] = BetR[4];
                            Beti[1] = Beti[4];
                            AXR = AlpR[2];
                            AXi = -Alpi[2];
                            Alpi[2] = -Alpi[2];
                            moder = 6;
                            //Console.WriteLine("BetR case 5");
                            //Console.WriteLine(BetR[1] + " " + BetR[2] + " " + BetR[3] + " " + BetR[4]);
                        }
                        else if (moder == 6)
                        {
                            BetR[2] = BetR[4];
                            Beti[2] = Beti[4];
                            AXR = AlpR[3];
                            AXi = -Alpi[3];
                            Alpi[3] = -Alpi[3];
                            L = 2;
                            moder = 3;
                            //Console.WriteLine("BetR case 6");
                            //Console.WriteLine(BetR[1] + " " + BetR[2] + " " + BetR[3] + " " + BetR[4]);
                        }
                        else if (moder == 4)
                        {
                            // -------------------  the convergence moder ------------------- 
                            skip = false;
                            if (Math.Abs(tempreal) + Math.Abs(tempimag) > 1.0e-40)
                            {
                                temp7 = Math.Abs(AlpR[3] - AXR) + Math.Abs(Alpi[3] - AXi);
                                if (temp7 / (Math.Abs(AXR) + Math.Abs(AXi)) > 1.0e-7)
                                {
                                    loopcnt = loopcnt + 1;
                                    for (i = 1; i <= 3; i++)
                                    {
                                        AlpR[i] = AlpR[i + 1];
                                        Alpi[i] = Alpi[i + 1];
                                        BetR[i] = BetR[i + 1];
                                        Beti[i] = Beti[i + 1];
                                    }
                                    if (loopcnt < 100)
                                    {
                                        DMulRSub(ref AlpR, ref Alpi, BetR, Beti);
                                        AXR = AlpR[4];
                                        AXi = Alpi[4];
                                        moder = 4;
                                        // Console.WriteLine("goto 80");
                                        skip = true;  // goto 80
                                    }
                                }
                            }
                            if (!skip)
                            {
                                rootcnt = rootcnt + 1;
                                RootS[rootcnt, 1] = AlpR[4];
                                RootS[rootcnt, 2] = Alpi[4];
                                Console.WriteLine(" found a root " + RootS[rootcnt, 1].ToString() + " " + RootS[rootcnt, 2].ToString());
                                loopcnt = 0;

                                if (rootcnt < nRoots)
                                {
                                    if (Math.Abs(RootS[rootcnt, 2]) > 1.0e-5)
                                    {
                                        if (L == 1)
                                        {
                                            AXR = AlpR[1];
                                            AXi = -Alpi[1];
                                            Alpi[1] = -Alpi[1];
                                            moder = 5;
                                        }
                                        else
                                        {
                                            AXR = 0.8;
                                            AXi = 0.0;
                                            L = 1;
                                            loopcnt = 1;
                                            AlpR[1] = AXR;
                                            Alpi[1] = AXi;
                                            moder = 1;
                                        }
                                    }
                                    else
                                    {
                                        AXR = 0.8;
                                        AXi = 0.0;
                                        L = 1;
                                        loopcnt = 1;
                                        AlpR[1] = AXR;
                                        Alpi[1] = AXi;
                                        moder = 1;
                                    }
                                }

                                //Console.WriteLine("BetR case 4");
                                //Console.WriteLine(BetR[1] + " " + BetR[2] + " " + BetR[3] + " " + BetR[4]);
                                //Console.WriteLine(AlpR[1] + " " + AlpR[2] + " " + AlpR[3] + " " + AlpR[4]);
                            } // if Skip 
                        }
                        else if (!skip)  //moder == 3 &&
                        {
                            BetR[3] = BetR[4];
                            Beti[3] = Beti[4];
                            DMulRSub(ref AlpR, ref Alpi, BetR, Beti);
                            AXR = AlpR[4];
                            AXi = Alpi[4];
                            moder = 4;
                            //Console.WriteLine("BetR case 3");
                            //Console.WriteLine(BetR[1] + " " + BetR[2] + " " + BetR[3] + " " + BetR[4]);
                        }
                    } // while 

                } // if rootcnt<> nRoots

            } // if nRoots > 0 

        } // factor 


        // ------------------------------------------------------------------------------
        //
        //                               function quadratic
        //
        //  this function solves for the two roots of a quadratic equation.  there are
        //    no restrictions on the coefficients, and imaginary results are passed
        //    out as separate values.  the general form is y = ax2 + bx + c.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    a           - coefficient of x squared term
        //    b           - coefficient of x term
        //    c           - constant
        //    opt         - option for output              I all roots including imaginary
        //                                                     R only real roots
        //                                                     U only unique real roots (no repeated)
        //
        //  outputs       :
        //    r1r         - real portion of root 1
        //    r1i         - imaginary portion of root 1
        //    r2r         - real portion of root 2
        //    r2i         - imaginary portion of root 2
        //
        //  locals        :
        //    discrim     - discriminate b2 - 4ac
        //
        //  coupling      :
        //    none.
        //
        //  references    :
        //    vallado       2013, 1027
        // ----------------------------------------------------------------------------*/

        public void quadratic(double a, double b, double c, char opt,
                              out double r1r, out double r1i, out double r2r, out double r2i)
        {
            const double small = 0.0000001;
            double discrim;
            // --------------------  implementation   // --------------------
            r1r = 0.0;
            r1i = 0.0;
            r2r = 0.0;
            r2i = 0.0;

            discrim = b * b - 4.0 * a * c;

            // ---------------------  real roots  // ------------------------
            if (Math.Abs(discrim) < small)
            {
                r1r = -b / (2.0 * a);
                r2r = r1r;
                if (opt == 'U')
                    r2r = 99999.9;
            }
            else
            {
                if (Math.Abs(a) < small)
                    r1r = -c / b;
                else
                {
                    if (discrim > 0.0)
                    {
                        r1r = (-b + Math.Sqrt(discrim)) / (2.0 * a);
                        r2r = (-b - Math.Sqrt(discrim)) / (2.0 * a);
                    }
                    else
                    {
                        // ------------------ complex roots --------------------
                        if (opt == 'I')
                        {
                            r1r = -b / (2.0 * a);
                            r2r = r1r;
                            r1i = Math.Sqrt(-discrim) / (2.0 * a);
                            r2i = -Math.Sqrt(-discrim) / (2.0 * a);
                        }
                        else
                        {
                            r1r = 99999.9;
                            r2r = 99999.9;
                        }
                    }
                }
            }
        }  // quadratic


        // ------------------------------------------------------------------------------
        //
        //                               function parabolicspl
        //
        //  this function performs parabolic splining of an 3 input data points. 
        //  the points do not need to be equally spaced. 
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    p0,p1,p2    - function values used for splining
        //    t0,t1,t2    - time values used for splining
        //
        //  outputs       :
        //    ap0..ap2    - splined polynomial coefficients.     ap2 * t1 * t1 + ap1 * t1 + ap0
        //
        //  locals        : none
        //
        //  coupling      :
        //    none
        //
        //  references    :
        //    vallado       2013, 1034
        // -----------------------------------------------------------------------------

        public void parabolicspl(double p1, double p2, double p3, double t1, double t2, double t3,
                                 out double ap0, out double ap1, out double ap2)
        {
            // be sure time is normalized for the first value
            t3 = t3 - t1;
            t2 = t2 - t1;
            t1 = 0.0;

            // parabolic coefficients
            ap0 = p1;
            ap1 = (-1.0 / t2 - 1.0 / t3) * p1 + (-t3 / (t2 * t2 - t2 * t3)) * p2 + (-t2 / (t3 * t3 - t2 * t3)) * p3;
            ap2 = (-(t2 - t3) / (t2 * t3 * t3 - t2 * t2 * t3)) * p1 + (1.0 / (t2 * t2 - t2 * t3)) * p2 +
                (t2 / (t2 * t3 * t3 - t2 * t2 * t3)) * p3;

            //// check
            //double ans1 = ap2 * t1 * t1 + ap1 * t1 + ap0;
            //double ans2 = ap2 * t2 * t2 + ap1 * t2 + ap0;
            //double ans3 = ap2 * t3 * t3 + ap1 * t3 + ap0;
            //Console.WriteLine((ans1-p1).ToString() + " " + (ans2 - p2).ToString() + " "+ (ans3 - p3).ToString() );
        }  // parabolicspl



        // ------------------------------------------------------------------------------
        //
        //                               function cubicspl
        //
        //  this function performs cubic splining of an input zero crossing
        //  function in order to find function values.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    p0,p1,p2,p3 - function values used for splining
        //    t0,t1,t2,t3 - time values used for splining
        //
        //  outputs       :
        //    acu0..acu3  - splined polynomial coefficients. acu3 t^3, etc
        //
        //  locals        : none
        //
        //  coupling      :
        //    none
        //
        //  references    :
        //    vallado       2013, 1034
        // -----------------------------------------------------------------------------

        public void cubicspl(double p1, double p2, double p3, double p4,
                             out double acu0, out double acu1, out double acu2, out double acu3)
        {
            acu0 = p2;
            acu1 = -p1 / 3.0 - 0.5 * p2 + p3 - p4 / 6.0;
            acu2 = 0.5 * p1 - p2 + 0.5 * p3;
            acu3 = -p1 / 6.0 + 0.5 * p2 - 0.5 * p3 + p4 / 6.0;
        }  // cubicspl


        // ------------------------------------------------------------------------------
        //
        //                               function cubic
        //
        //  this function solves for the three roots of a cubic equation.  there are
        //    no restrictions on the coefficients, and imaginary results are passed
        //    out as separate values.  the general form is y = a3x3 + b2x2 + c1x + d0.  note
        //    that r1i will always be zero Math.Since there is always at least one real root.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    a3          - coefficient of x cubed term
        //    b2          - coefficient of x squared term
        //    c1          - coefficient of x term
        //    d0          - constant
        //    opt         - option for output              I all roots including imaginary
        //                                                     R only real roots
        //                                                     U only unique real roots (no repeated)
        //
        //  outputs       :
        //    r1r         - real portion of root 1
        //    r1i         - imaginary portion of root 1
        //    r2r         - real portion of root 2
        //    r2i         - imaginary portion of root 2
        //    r3r         - real portion of root 3
        //    r3i         - imaginary portion of root 3
        //
        //  locals        :
        //    temp1       - temporary value
        //    temp2       - temporary value
        //    p           - coefficient of x squared term where x cubed term is 1.0
        //    q           - coefficient of x term where x cubed term is 1.0
        //    r           - coefficient of constant term where x cubed term is 1.0
        //    delta       - discriminator for use with cardans formula
        //    e0          - angle holder for trigonometric solution
        //    phi         - angle used in trigonometric solution
        //    cosphi      - cosine of phi
        //    sinphi      - sine of phi
        //
        //  coupling      :
        //    quadric     - quadratic roots
        //
        //  references    :
        //    vallado       2013, 1027
        // -----------------------------------------------------------------------------

        public void cubic(double a3, double b2, double c1, double d0, char opt,
                          out double r1r, out double r1i, out double r2r, out double r2i, out double r3r, out double r3i)
        {
            const double rad = 57.29577951308230;
            const double onethird = 1.0 / 3.0;
            const double small = 0.00000001;
            double temp1, temp2, p, q, r, delta, e0, cosphi, sinphi, phi;
            // ------------------------  implementation   // ------------------------
            r1r = 0.0;
            r1i = 0.0;
            r2r = 0.0;
            r2i = 0.0;
            r3r = 0.0;
            r3i = 0.0;

            if (Math.Abs(a3) > small)
            {
                // ------------- force coefficients into std form -------------------
                p = b2 / a3;
                q = c1 / a3;
                r = d0 / a3;

                a3 = onethird * (3.0 * q - p * p);
                b2 = (1.0 / 27.0) * (2.0 * p * p * p - 9.0 * p * q + 27.0 * r);

                delta = (a3 * a3 * a3 / 27.0) + (b2 * b2 * 0.25);

                // -------------------- use cardans formula ------------------------
                if (delta > small)
                {
                    temp1 = (-b2 * 0.5) + Math.Sqrt(delta);
                    temp2 = (-b2 * 0.5) - Math.Sqrt(delta);
                    temp1 = Math.Sign(temp1) * Math.Pow(Math.Abs(temp1), onethird);
                    temp2 = Math.Sign(temp2) * Math.Pow(Math.Abs(temp2), onethird);
                    r1r = temp1 + temp2 - p * onethird;

                    if (opt == 'I')
                    {
                        r2r = -0.5 * (temp1 + temp2) - p * onethird;
                        r2i = -0.5 * Math.Sqrt(3.0) * (temp1 - temp2);
                        r3r = -0.5 * (temp1 + temp2) - p * onethird;
                        r3i = -r2i;
                    }
                    else
                    {
                        r2r = 99999.9;
                        r3r = 99999.9;
                    }
                }
                else
                {
                    // -------------------- evaluate zero point ------------------------
                    if (Math.Abs(delta) < small)
                    {
                        r1r = -2.0 * Math.Sign(b2) * Math.Pow(Math.Abs(b2 * 0.5), onethird) - p * onethird;
                        r2r = Math.Sign(b2) * Math.Pow(Math.Abs(b2 * 0.5), onethird) - p * onethird;
                        if (opt == 'U')
                            r3r = 99999.9;
                        else
                            r3r = r2r;
                    }
                    else
                    {
                        // --------------- use trigonometric identities ----------------
                        e0 = 2.0 * Math.Sqrt(-a3 * onethird);
                        cosphi = (-b2 / (2.0 * Math.Sqrt(-a3 * a3 * a3 / 27.0)));
                        sinphi = Math.Sqrt(1.0 - cosphi * cosphi);
                        phi = Math.Atan2(sinphi, cosphi);
                        if (phi < 0.0)
                            phi = phi + 2.0 * Math.PI;
                        r1r = e0 * Math.Cos(phi * onethird) - p * onethird;
                        r2r = e0 * Math.Cos(phi * onethird + 120.0 / rad) - p * onethird;
                        r3r = e0 * Math.Cos(phi * onethird + 240.0 / rad) - p * onethird;
                    } // if fabs(delta) 
                }  // if delta > small
            }  // if fabs > small
            else
            {
                quadratic(b2, c1, d0, opt, out r1r, out r1i, out r2r, out r2i);
                r3r = 99999.9;
                r3i = 99999.9;
            }
        }  // cubic


        // ------------------------------------------------------------------------------
        //
        //                               function cubicinterp
        //
        //  this function performs a cubic spline. four points are needed.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    valuein     - kp
        //
        //  outputs       :
        //    out         - ap
        //
        //  locals        :
        //                    -
        //
        //  coupling      :
        //    cubicspl
        //
        //  references    :
        //    vallado       2013, 1027
        // -----------------------------------------------------------------------------

        public double cubicinterp(double p1a, double p1b, double p1c, double p1d, double p2a, double p2b,
                                  double p2c, double p2d, double valuein)
        {
            double kc0, kc1, kc2, kc3, ac0, ac1, ac2, ac3,
                   r1r, r1i, r2r, r2i, r3r, r3i, value;

            // -------- assign function points ---------
            cubicspl(p1a, p1b, p1c, p1d, out ac0, out ac1, out ac2, out ac3);
            cubicspl(p2a, p2b, p2c, p2d, out kc0, out kc1, out kc2, out kc3);

            // recover the original function values
            // use the normalized time first, but at an arbitrary interval
            cubic(kc3, kc2, kc1, kc0 - valuein, 'R', out r1r, out r1i, out r2r, out r2i, out r3r, out r3i);

            if ((r1r >= -0.000001) && (r1r <= 1.001))
            {
                value = r1r;
            }
            else
            {
                if ((r2r >= -0.000001) && (r2r <= 1.001))
                {
                    value = r2r;
                }
                else
                {
                    if ((r3r >= -0.000001) && (r3r <= 1.001))
                    {
                        value = r3r;
                    }
                    else
                    {
                        value = 0.0;
                        Console.Write("error in cubicinterp root {0} {1} {2} {3} \n", valuein, r1r, r2r, r3r);
                    }
                }
            }

            return (ac3 * Math.Pow(value, 3) + ac2 * value * value + ac1 * value + ac0);
        } // cubicinterp   


        // ==============================================================================
        //                                 Misc routines  
        // ==============================================================================

        // -----------------------------------------------------------------------------
        //
        //                               function percentile
        //
        //  this function performs a percentile calculation. the interpolation method will
        //  vary between routines, so matlab, python, etc will not all behave the same. samira
        //  found that matlab is likely using a Hyndman & Fan’s R-6 interpolation while this
        //  approach and python are using an R7 interpolation. 
        //  30th percentile, first array(sorted: 0.03, 5.13, 5.63, 45.3, 45.3, 345.3, 3445.3) :
        //  c# Type 7 rank = (7-1)×0.30+1 = 2.8 → interpolate 80% between sorted[2] = 5.13 and
        //     sorted[3] = 5.63 → 5.53
        //  MATLAB Type 5: breakpoints at(k-0.5)/7; p=30% falls 60% between the k=2 point(21.4%,
        //     value 5.13) and k = 3 point(35.7%, value 5.63) → 5.43 
        //  42nd percentile, same thing c# -0.042352, MATLAB's Type 5 -0.102544.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    inArr       - array of values
        //    excelpercentile - percentile desired                   0.0 to 1.0
        //    arrSize     - size of array
        //
        //  outputs       :
        //    percentile  - result
        //
        // -----------------------------------------------------------------------------

        public double Percentile(double[] inArr, double excelPercentile, Int32 arrSize)
        {
            double b;
            b = 0.0;
            int k;
            k = 0;
            Array.Sort(inArr, 0, arrSize);  // needs sort for fractional % calcs
            if (arrSize > 0)
            {
                if (arrSize == 1)
                    b = inArr.Min(); 
                else
                {
                    double n = (arrSize - 1) * excelPercentile + 1;
                    for (int i = 0; i < arrSize; i++)
                    {
                        if (n == arrSize)
                            b = inArr.Max();
                        else
                        {
                            k = (int)n;
                            double d = n - k;
                            b = inArr[k - 1] + d * (inArr[k] - inArr[k - 1]);
                        }
                    }  // for
                }
            }
            else
                b = 0.0;
            
            return b;
        }  // percentile



        // ------------------------------------------------------------------------------
        //
        //                               function factorial
        //
        //  this function performs a factorial. note the use of double in the return as the 
        //  numbers get huge quickly. this is good to about n = 170. 
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    n           - order in
        //
        //  outputs       :
        //    factorial   - result
        //
        //  locals        :
        //                    -
        //
        //  coupling      :
        //    none
        // -----------------------------------------------------------------------------

        public double factorial (int n)
        {
            if (n == 0)
                return 1.0;
            if (n >= 1)
                return n * factorial(n - 1);
            else
                return 0.0;
        }  // factorial 


        // ------------------------------------------------------------------------------
        //
        //                               function linearfit
        //
        //  this function fits a line to 2 input data points. 
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    p0,p1       - function values used for fit
        //    t0,t1       - time values used for fit
        //
        //  outputs       :
        //    al0, al1    - lienar coefficients.                 al0 + al1 * t
        //
        //  locals        : none
        //
        //  coupling      :
        //    none
        //
        //  references    :
        //    vallado       2013, 1034
        // -----------------------------------------------------------------------------

        public void linearfit(double p1, double p2, double t1, double t2, out double al0, out double al1)
        {
            // linear coefficients
            al1 = (p2 - p1) / (t2 - t1);
            al0 = -al1 * t1 + p1;

            //// check
            //double ans1 = al1 * t1 + al0;
            //double ans2 = al1 * t2 + al0;
            //Console.WriteLine(al0 + " " + al1 + " " + (ans1 - p1).ToString() + " a " + (ans2 - p2).ToString() );
        }  // linearfit




        // ==============================================================================
        //                                   time functions
        // ==============================================================================

        // -----------------------------------------------------------------------------------------
        // convert a string month to integer

        public int getIntMonth(string monthstr)
        {
            int strNumber, strIndex;
            string[] monstr = new string[13] { "", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC" };

            strIndex = 0;
            for (strNumber = 0; strNumber < monstr.Length; strNumber++)
            {
                strIndex = monstr[strNumber].IndexOf(monthstr.ToUpper());
                if (strIndex >= 0)
                    break;
            }
            return strNumber;
        } // getIntMonth


        // ------------------------------------------------------------------------------
        //                                    JD to STK time 
        //
        // converts a julian date into STK time string. use 2 factor julian date to 
        // preserve milli-second accuracy. remove UTCG first if it's there
        // 
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    JD          - julian date (integer part)                    24516800. days
        //    JDF         - julian date (day fraction part)               0.37184856 days
        //
        //  outputs       :
        //    epoch       - epoch time                                    14 Jul 2012 18:00:00.000 UTCG    
        //    
         // ------------------------------------------------------------------------------

        public void JD2STKtime(double JD, double JDF, out string epoch)
        {
            string[] monstr = new string[13] { "", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC" };

            int year, mon, day, hr, minute;
            double second;
            invjday(JD, JDF, out year, out mon, out day, out hr, out minute, out second);
            epoch = day.ToString("00") + " " + monstr[mon] + " " + year.ToString() + " " +
                 hr.ToString("00") + ":" + minute.ToString("00") + ":" + second.ToString("00.000000000");
        }  // JD2STKtime


        // ------------------------------------------------------------------------------
         //                                    STK time to JD
         //
         // converts a STK string epoch time into a julian date. use 2 factor julian date to 
         // preserve milli-second accuracy. remove UTCG first if it's there
         // 
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
         //    epoch       - epoch time                                    14 Jul 2012 18:00:00.000 UTCG    
         //
         //  outputs       :
         //    JD          - julian date (integer part)                    24516800. days
         //    JDF         - julian date (day fraction part)               0.37184856 days
         //    
          // ------------------------------------------------------------------------------

        public void STKtime2JD(string epoch, out double JD, out double JDF)
        {
            int year, mon, day, hr, minute;
            double second;

            if (epoch.Contains("UTCG"))
                epoch = epoch.Replace("UTCG", "");
            epoch = epoch.Trim();
            string[] hms = epoch.Split(':');
            minute = Convert.ToInt16(hms.GetValue(1).ToString());
            second = Convert.ToDouble(hms.GetValue(2).ToString());
            string[] dmy = hms.GetValue(0).ToString().Split(' ');
            year = Convert.ToInt16(dmy.GetValue(2).ToString());
            string monthstr = dmy.GetValue(1).ToString();
            mon = getIntMonth(monthstr);
            day = Convert.ToInt16(dmy.GetValue(0).ToString());
            hr = Convert.ToInt16(dmy.GetValue(3).ToString());

            jday(year, mon, day, hr, minute, second, out JD, out JDF);
            DateTime startday = new DateTime(year, getIntMonth(monthstr), day, hr, minute, Convert.ToInt16(second));
        }  // STKtime2JD


        // -----------------------------------------------------------------------------
        // one approach using raw variables

        public void getTimeFromGPS(Int32 weeknumber, Int32 gpstime, Int32 numrollover,
                                   out Int32 year, out Int32 mon, out Int32 day, out Int32 hr, out Int32 minute, out double second)
        {
            double jd80, jd80f, jd, jdf;
            jday(1980, 1, 6, 0, 0, 0.0, out jd80, out jd80f);

            jd = jd80 + weeknumber * 7 + 1024 * 7 * numrollover;
            jdf = jd80f + gpstime / 86400.0;
            invjday(jd, jdf, out year, out mon, out day, out hr, out minute, out second);
        }


        // -----------------------------------------------------------------------------
        // another approach using datetime method

        public DateTime getTimeFromGPS1(Int32 weeknumber, Int32 gpstime, Int32 numrollover)
        {
            //double jd80, jd80f, jd, jdf;
            // gps epoch
            DateTime datum = new DateTime(1980, 1, 6, 0, 0, 0);
            DateTime week = datum.AddDays(Math.Floor(weeknumber / 10.0) * 7 + weeknumber - Math.Floor(weeknumber / 10.0) * 10.0);
            DateTime time = week.AddSeconds(gpstime);
            return time;
        }

        // -----------------------------------------------------------------------------
        // get gps weeks from year month day

        public void getGPSFromTime(Int32 year, Int32 month, Int32 day, out Int32 GPSWeek, out Int32 GPSWeekF)
        {
            // int is the same as Int32
            // gps epoch
            DateTime datum = new DateTime(1980, 1, 6, 0, 0, 0);
            DateTime datum1 = new DateTime(year, month, day, 0, 0, 0);
            TimeSpan days = datum1.Subtract(datum);
            double weeks = Math.Floor(days.Days / 7.0);

            GPSWeek = Convert.ToInt32(weeks);
            GPSWeekF = Convert.ToInt32(weeks) * 10 + Convert.ToInt32(days.Days % 7);
        }


        // ------------------------------------------------------------------------------
        //
        //                               procedure jday
        //
        //  this procedure finds the julian date given the year, month, day, and time.
        //    the julian date is defined by each elapsed day since noon, jan 1, 4713 bc.
        //    two values are passed back for improved accuracy
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    year        - year                           1900 .. 2100
        //    mon         - month                          1 .. 12
        //    day         - day                            1 .. 28,29,30,31
        //    hr          - universal time hour            0 .. 23
        //    minute         - universal time minute             0 .. 59
        //    second        - universal time second            0.0 .. 59.999
        //
        //  outputs       :
        //    jd          - julian date (days only)           days from 4713 bc
        //    jdFrac      - julian date (fraction of a day)   days from 0 hr of the day
        //
        //  locals        :
        //    none.
        //
        //  coupling      :
        //    none.
        //
        //  references    :
        //    vallado       2013, 183, alg 14, ex 3-14
        // -----------------------------------------------------------------------------

        public void jday(int year, int mon, int day, int hr, int minute, double second, out double jd, out double jdFrac)
        {
            jd = 367.0 * year -
                 Math.Floor((7 * (year + Math.Floor((mon + 9) / 12.0))) * 0.25) +
                 Math.Floor(275 * mon / 9.0) +
                 day + 1721013.5;  // use - 678987.0 to go to mjd directly
            jdFrac = (second+ minute * 60.0 + hr * 3600.0) / 86400.0;

            // check that the day and fractional day are correct
            if (Math.Abs(jdFrac) >= 1.0)
            {
                double dtt = Math.Floor(jdFrac);
                jd = jd + dtt;
                jdFrac = jdFrac - dtt;
            }

            // - 0.5*Math.Sign(100.0*year + mon - 190002.5) + 0.5;
        }  //  jday


        // ------------------------------------------------------------------------------
        //
        //                               procedure days2mdhms
        //
        //  this procedure converts the day of the year, days, to the equivalent month
        //    day, hour, minute and second.
        //
        //  algorithm     : set up array for the number of days per month
        //                      find leap year - use 1900 because 2000 is a leap year
        //                      loop through a temp value while the value is < the days
        //                      perform int conversions to the correct day and month
        //                      convert remainder into h m s using type conversions
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    year        - year                           1900 .. 2100
        //    days        - julian day of the year         0.0  .. 366.0
        //
        //  outputs       :
        //    mon         - month                          1 .. 12
        //    day         - day                            1 .. 28,29,30,31
        //    hr          - hour                           0 .. 23
        //    minute         - minute                         0 .. 59
        //    second      - second                         0.0 .. 59.999
        //
        //  locals        :
        //    dayofyr     - day of year
        //    temp        - temporary extended values
        //    inttemp     - temporary int value
        //    i           - index
        //    lmonth[12]  - int array containing the number of days per month
        //
        //  coupling      :
        //    none.
        // -----------------------------------------------------------------------------

        public void days2mdhms(int year, double days,
                               out int mon, out int day, out int hr, out int minute, out double second)
        {
            int i, inttemp, dayofyr;
            double temp;
            int[] lmonth = new int[13] { 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

            dayofyr = (int)Math.Floor(days);
            // ------------------ find month and day of month ------------------
            if ((year % 4) == 0)
                lmonth[2] = 29;

            i = 1;
            inttemp = 0;
            while ((dayofyr > inttemp + lmonth[i]) && (i < 12))
            {
                inttemp = inttemp + lmonth[i];
                i = i + 1;
            }
            mon = i;
            day = dayofyr - inttemp;

            // ------------------ find hours minutes and seconds ---------------
            temp = (days - dayofyr) * 24.0;
            hr = Convert.ToInt16(Math.Floor(temp));
            temp = (temp - hr) * 60.0;
            minute = Convert.ToInt16(Math.Floor(temp));
            second= (temp - minute) * 60.0;
        }  //  days2mdhms


        // -------------------------------------------------------------------------------
        //
        //                               function finddays
        //
        //  this function finds the fractional days through a year given the year,
        //    month, day, hour, minute and second.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    year        - year                                     1900 .. 2100
        //    month       - month                                    1 .. 12
        //    day         - day                                      1 .. 28,29,30,31
        //    hr          - hour                                     0 .. 23
        //    min         - minute                                   0 .. 59
        //    second      - second                                   0.0 .. 59.999
        //
        //  outputs       :
        //    days        - day of year plus fraction of a
        //                        day days
        //
        //  locals        :
        //    lmonth      - length of months of year
        //
        //  references    :
        //    vallado       2007, 207, ex 3-12
        //
        // -------------------------------------------------------------------------------

        public void findDays(int year, int month, int day, int hr, int min, double second, out int dayofyr)
        {
            int i;
            // shift index values to be 1-12
            int[] lmonth = new int[13] {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

            for (i = 1; i <= 12; i++)
            {
                lmonth[i] = 31;
                if (i == 2)
                    lmonth[i] = 28;
                if (i == 4 || i == 6 || i == 9 || i == 11)
                    lmonth[i] = 30;
            }

            if ((year - 1900) % 4 == 0)
            {
                lmonth[2] = 29;
                if (year % 100 == 0 && !(year % 400 == 0))
                    lmonth[2] = 28;
            }

            i = 1;
            dayofyr = 0;
            while (i < month && i < 12)
                    {
                dayofyr = dayofyr + lmonth[i];
                i = i + 1;
            }

            dayofyr = dayofyr + day;  // + hr / 24.0 + min / 1440.0 + second / 86400.0;
        }



        // ------------------------------------------------------------------------------
        //
        //                               procedure invjday
        //
        //  this procedure finds the year, month, day, hour, minute and second
        //  given the julian date. tu can be ut1, tdt, tdb, etc. jd is input
        //  with two arguments for additional accuracy
        //
        //  algorithm     : set up starting values
        //                      find leap year - use 1900 because 2000 is a leap year
        //                      find the elapsed days through the year in a loop
        //                      call routine to find each individual value
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    jd          - julian date (days only)           days from 4713 bc
        //    jdFrac      - julian date (fraction of a day)   days from 0 hr of the day
        //
        //  outputs       :
        //    year        - year                           1900 .. 2100
        //    mon         - month                          1 .. 12
        //    day         - day                            1 .. 28,29,30,31
        //    hr          - hour                           0 .. 23
        //    minute         - minute                         0 .. 59
        //    second        - second                         0.0 .. 59.999
        //
        //  locals        :
        //    days        - day of year plus fractional
        //                      portion of a day               days
        //    tu          - julian centuries from 0 h
        //                      jan 0, 1900
        //    temp        - temporary double values
        //    leapyrs     - number of leap years from 1900
        //
        //  coupling      :
        //    days2mdhms  - finds month, day, hour, minute and second given days and year
        //
        //  references    :
        //    vallado       2013, 202, alg 22, ex 3-13
        // -----------------------------------------------------------------------------

        public void invjday(double jd, double jdFrac,
                            out int year, out int mon, out int day, out int hr, out int minute, out double second)
        {
            int leapyrs;
            double dt, days, tu, temp;

            // check jdfrac for multiple days
            if (Math.Abs(jdFrac) >= 1.0)
            {
                jd = jd + Math.Floor(jdFrac);
                jdFrac = jdFrac - Math.Floor(jdFrac);
            }

            // check for fraction of a day included in the jd
            dt = jd - Math.Floor(jd) - 0.5;
            if (Math.Abs(dt) > 0.00000001)
            {
                jd = jd - dt;
                jdFrac = jdFrac + dt;
            }

            // ---------------- find year and days of the year -----------------
            temp = jd - 2415019.5;
            tu = temp / 365.25;
            year = 1900 + Convert.ToInt16(Math.Floor(tu));
            leapyrs = Convert.ToInt16(Math.Floor((year - 1901) * 0.25));

            days = Math.Floor(temp - ((year - 1900) * 365.0 + leapyrs));

            // ------------- check for case of beginning of a year -------------
            if (days + jdFrac < 1.0)
            {
                year = year - 1;
                leapyrs = Convert.ToInt16(Math.Floor((year - 1901) * 0.25));
                days = Math.Floor(temp - ((year - 1900) * 365.0 + leapyrs));
            }

            // ------------------ find remaining data  // -----------------------
            // now add the daily time in to preserve accuracy
            days2mdhms(year, days + jdFrac, out mon, out day, out hr, out minute, out second);
        }  // invjday



        // ------------------------------------------------------------------------------
        //
        //                               procedure hms_sec
        //
        //  this procedure converts hours, minutes and seconds to seconds from the
        //  beginning of the day.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    utsecond    - seconds                        0.0 .. 86400.0
        //
        //  outputs       :
        //    hr          - hours                          0 .. 24
        //    minute      - minutes                        0 .. 59
        //    second      - seconds                        0.0 .. 59.99
        //
        //  locals        :
        //    temp        - temporary variable
        //
        //  coupling      :
        //    none.
        // -----------------------------------------------------------------------------

        public void hms_sec(ref int hr, ref int minute, ref double second, Enum direct, ref double utsec)
        {
            double temp;

            // ------------------------  implementation   // ----------------
            if (direct.Equals(Edirection.eto))
                utsec= hr * 3600.0 + minute * 60.0 + second;
            else
            {
                temp = utsec/ 3600.0;
                hr = (int)Math.Truncate(temp);
                minute = (int)Math.Truncate((temp - hr) * 60.0);
                second= (temp - hr - minute / 60.0) * 3600.0;
            }
        }  // hms_sec


        // ------------------------------------------------------------------------------
        //
        //                               procedure hms_ut
        //
        //  this procedure converts hours, minutes and seconds into universal time.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    hr          - hours                          0 .. 24
        //    minute      - minutes                        0 .. 59
        //    second      - seconds                        0.0 .. 59.99
        //    direction   - which set of vars to output    from  too
        //
        //  outputs       :
        //    ut          - universal time                 hrmin.sec
        //
        //  locals        :
        //    none.
        //
        //  coupling      :
        //    none.
        //
        //  references    :
        //    vallado       2013, 199, alg 21, ex 3-10
        // -----------------------------------------------------------------------------

        public void hms_ut (ref int hr, ref int minute, ref double second, Enum direct, ref double ut)
        {
            // ------------------------  implementation   // ----------------
            if (direct.Equals(Edirection.eto))
                ut = hr * 100.0 + minute + second* 0.01;
            else
            {
                hr = (int)Math.Truncate(ut * 0.01);
                minute = (int)Math.Truncate(ut - hr * 100.0);
                second= (ut - hr * 100.0 - minute) * 100.0;
            }
        }  // hms_ut


        // ------------------------------------------------------------------------------
        //
        //                               procedure hms_rad
        //
        //  this procedure converts hours, minutes and seconds into radians.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    hr          - hours                          0 .. 24
        //    minute      - minutes                        0 .. 59
        //    second      - seconds                        0.0 .. 59.99
        //    direction   - which set of vars to output    from  too
        //
        //  outputs       :
        //    hms         - result                         rad
        //
        //  locals        :
        //    temp        - conversion from hours to rad   0.261799
        //
        //  coupling      :
        //    none.
        //
        //  references    :
        //    vallado       2013, 198, alg 19 alg 20, ex 3-9
        // -----------------------------------------------------------------------------

        public void hms_rad (ref int hr, ref int minute, ref double second, Enum direct, ref double hms)
        {
            const double rad2deg = 57.29577951308230;
            double temp;

            // ------------------------  implementation   // ----------------
            temp = 15.0 / rad2deg;
            if (direct.Equals(Edirection.eto))
                hms = (hr + minute / 60.0 + second/ 3600.0) * temp;
            else
            {
                temp = hms / temp;
                hr = Convert.ToInt32(Math.Truncate(temp));
                minute = Convert.ToInt32(Math.Truncate((temp - hr) * 60.0));
                second= (temp - hr - minute / 60.0) * 3600.0;
            }
        }  // hms_rad


        // ------------------------------------------------------------------------------
        //
        //                               procedure dms_rad
        //
        //  this procedure converts degrees, minutes and seconds into radians.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    deg         - degrees                        0 .. 360
        //    minute      - minutes                        0 .. 59
        //    second      - seconds                        0.0 .. 59.99
        //    direction   - which set of vars to output    from  too
        //
        //  outputs       :
        //    dms         - result                         rad
        //
        //  locals        :
        //    temp        - temporary variable
        //
        //  coupling      :
        //    none.
        //
        //  references    :
        //    vallado       2013, 197, alg 17 alg 18, ex 3-8
        // -----------------------------------------------------------------------------

        public void dms_rad (ref int deg, ref int minute, ref double second, Enum direct, ref double dms)
        {
            const double rad2deg = 57.29577951308230;
            double temp;

            // ------------------------  implementation   // ----------------
            if (direct.Equals(Edirection.eto))
                dms = (deg + minute / 60.0 + second/ 3600.0) / rad2deg;
            else
            {
                temp = dms * rad2deg;
                deg = (int)Math.Truncate(temp);
                minute = (int)Math.Truncate((temp - deg) * 60.0);
                second= (temp - deg - minute / 60.0) * 3600.0;
            }
        }  // dms_rad


        // -------------------------------------------------------------------------------
        //
        //                               function convtime
        //
        //  this function finds the time parameters and julian century values for inputs
        //    of utc or ut1.numerous outputs are found as shown in the local variables.
        //    because calucations are in utc, you must include timezone if (you enter a
        //    local time, otherwise it should be zero.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    year        - year                                         1900 .. 2100
        //    mon         - month                                        1 .. 12
        //    day         - day                                          1 .. 28,29,30,31
        //    hr          - universal time hour                          0 .. 23
        //    minute      - universal time minute                           0 .. 59
        //    second      - universal time second(utc)                     0.0  .. 59.999
        //    timezone    - offset to utc from local site                0 .. 23 hr
        //    dut1        - delta of ut1 - utc sec
        //    dat         - delta of tai - utc sec
        //
        //  outputs       :
        //    ut1         - universal time                               sec
        //    tut1        - julian centuries of ut1
        //    jdut1       - julian date (days only) days from 4713 bc
        //    jdut1Frac   - julian date (fraction of a day)              days from 0 hr of the day
        //    utc         - coordinated universal time                   sec
        //    tai         - atomic time                                  sec
        //    tdt         - terrestrial dynamical time                   sec
        //    ttdt        - julian centuries of tdt
        //    jdtt        - julian date(days only)                       days from 4713 bc
        //    jdttFrac    - julian date(fraction of a day)               days from 0 hr of the day
        //    tdb         - terrestrial barycentric time                 sec
        //    ttdb        - julian centuries of tdb
        //    jdtdb       - julian date of tdb                           days from 4713 bc
        //    tcb         - celestial barycentric time                   sec
        //    tcg         - celestial geocentric time                    sec
        //    jdtdb       - julian date(days only)                       days from 4713 bc
        //    jdtdbFrac   - julian date(fraction of a day)               days from 0 hr of the day
        //
        //  locals        :
        //    hrtemp      - temporary hours                              hr
        //    mintemp     - temporary minutes                            minute
        //    sectemp     - temporary seconds                            sec
        //    localhr     - difference to local time                     hr
        //    jd          - julian date of request                       days from 4713 bc
        //    me          - mean anomaly of the earth                    rad
        //
        //  coupling      :
        //    hms_2_second  - conversion between hr-minute-sec.and.seconds
        //    jday        - find the julian date
        //
        //  references    :
        //    vallado       2013, 201, alg 16, ex 3-7
        // ------------------------------------------------------------------------------*/

        public void convtime(int year, int mon, int day, int hr, int minute, double second, int timezone, double dut1, int dat,
                             out double ut1, out double tut1, out double jdut1, out double jdut1frac, out double utc, out double tai,
                             out double tt, out double ttt, out double jdtt, out double jdttfrac,
                             out double tdb, out double ttdb, out double jdtdb, out double jdtdbfrac)
        {
            const double rad2deg = 57.29577951308230;
            double mjd, mfme, jd, jdfrac, jdtai, jdtaifrac, sectemp, me, tdb2, ttdb2, jdtdb2, jdtdb2frac, dlje;
                       //tcbmtdb, tcb, ttcb, jdtcb, jdtcbfrac;
            int localhr, mintemp;

            // initialize
            utc = 0.0;
            mintemp = 0;
            sectemp = 0.0;

            // ------------------------  implementation   // ----------------
            jday(year, mon, day, hr + timezone, minute, second, out jd, out jdfrac);
            mjd = jd + jdfrac - 2400000.5;
            mfme = hr * 60.0 + minute + second/ 60.0;

            // ------------------ start if (ut1 is known ------------------
            localhr = timezone + hr;
            hms_sec(ref localhr, ref minute, ref second, Edirection.eto, ref utc);

            ut1 = utc + dut1;
            hms_sec(ref localhr, ref mintemp, ref sectemp, Edirection.efrom, ref ut1);
            jday(year, mon, day, localhr, mintemp, sectemp, out jdut1, out jdut1frac);
            tut1 = ((jdut1 - 2451545.0) + jdut1frac) / 36525.0;

            tai = utc + dat;
            hms_sec(ref localhr, ref mintemp, ref sectemp, Edirection.efrom, ref tai);
            jday(year, mon, day, localhr, mintemp, sectemp, out jdtai, out jdtaifrac);

            tt = tai + 32.184;   // sec
            hms_sec(ref localhr, ref mintemp, ref sectemp, Edirection.efrom, ref tt);
            jday(year, mon, day, localhr, mintemp, sectemp, out jdtt, out jdttfrac);
            ttt = (jdtt + jdttfrac - 2451545.0) / 36525.0;

            // tdb
            // std approach(digits)
            //         me= 357.53  + 0.9856003 * (jdtt - 2451545.0);   
            //         me= mod(me,360.0  );
            //         me= me* deg2rad;
            //         tdb1= tt + 0.001658  * sin(me) + 0.000014 *sin(2.0 *me);
            //         [hrtemp, mintemp, sectemp] = sec2hms(tdb1 );
            //         [jdtdb1, jdtdb1frac] = jday(year, mon, day, hrtemp, mintemp, sectemp );
            //         ttdb1= (jdtdb1 + jdtdb1frac - 2451545.0  )/ 36525.0;
            //         fprintf(1,'std  tdb %8.6f ttdb  %16.12f jdtdb  %18.11f %18.11f \n', tdb1, ttdb1, jdtdb1, jdtdb1frac );
            // ast alm approach(2012) bradley email
            me = 357.53 + 0.98560028 * (jdtt - 2451545.0);
            me = me % 360.0;
            me = me / rad2deg;
            dlje = 246.11 + 0.90251792 * (jdtt - 2451545.0);
            tdb2 = tt + 0.001657 * Math.Sin(me) + 0.000022 * Math.Sin(dlje);
            hms_sec(ref localhr, ref mintemp, ref sectemp, Edirection.efrom, ref tdb2);
            jday(year, mon, day, localhr, mintemp, sectemp, out jdtdb2, out jdtdb2frac);
            ttdb2 = (jdtdb2 + jdtdb2frac - 2451545.0) / 36525.0;
            //fprintf(1,'asta tdb %8.6f ttdb  %16.12f jdtdb  %18.11f %18.11f \n', tdb2, ttdb2, jdtdb2, jdtdb2frac );
            // usno circular approach
            tdb = tt + 0.001657 * Math.Sin(628.3076 * ttt + 6.2401)
                   + 0.000022 * Math.Sin(575.3385 * ttt + 4.2970)
                   + 0.000014 * Math.Sin(1256.6152 * ttt + 6.1969)
                   + 0.000005 * Math.Sin(606.9777 * ttt + 4.0212)
                   + 0.000005 * Math.Sin(52.9691 * ttt + 0.4444)
                   + 0.000002 * Math.Sin(21.3299 * ttt + 5.5431)
                   + 0.000010 * ttt * Math.Sin(628.3076 * ttt + 4.2490);  // USNO circ(14)
            hms_sec(ref localhr, ref mintemp, ref sectemp, Edirection.efrom, ref tdb);
            jday(year, mon, day, localhr, mintemp, sectemp, out jdtdb, out jdtdbfrac);
            ttdb = (jdtdb + jdtdbfrac - 2451545.0) / 36525.0;

            //fprintf(1,'usno tdb %8.6f ttdb  %16.12f jdtdb  %18.11f %18.11f \n', tdb, ttdb, jdtdb, jdtdbfrac );
            //hms_sec(ref localhr, ref minute, ref second, Edirection.efrom, ref tdb);
            //        fprintf(1,'hms %3i %3i %8.6f \n', h, m, s);


            // tcg
            // approx with tai
            //tcg = tt + 6.969290134e-10 * (jdtai - 2443144.5003725) * 86400.0;  // AAS 05-352 (10) and IERS TN(104)
            //hms_sec(ref localhr, ref minute, ref second, Edirection.efrom, ref tcg);
            //jday(year, mon, day, hrtemp, mintemp, sectemp, out jdtcg, out jdtcgfrac);
            //tt2 = tcg - 6.969290134e-10 * (jdtcg + jdtcgfrac - 2443144.5003725) * 86400.0;

            //fprintf(1,'tcg %8.6f jdtcg  %18.11f ', tcg, jdtcg );
            //hms_sec(ref localhr, ref minute, ref second, Edirection.efrom, ref tcg);
            //        fprintf(1,'hms %3i %3i %8.6f \n', h, m, s);        

            // binomial approach with days
            //        lg=6.969290134e-10*86400.0;
            //        tcg1 = tt + (jdtt - 2443144.5003725)*(lg + lg* lg + lg* lg*lg);
            // days from 77
            //        jdttx = jday(year, mon, day, 0, 0, 0.0); 
            //        ttx = tt/86400.0 + jdttx-2443144.5003725  % days from the 1977 epoch
            //        tcg2 = (jdttx - 6.969290134e-10*2443144.5003725) / (1.0 - 6.969290134e-10) % days
            //        tcg2 = (tcg2 - jdttx)*86400*86400;
            // secondfrom 77
            //        ttx = tt + (jdttx-2443144.5003725)*86400.0;  % s from the 1977 epoch
            //        tcg3 = ttx / (1.0 - 6.969290134e-10); % s
            //        tcg3 = tcg3 -(jdttx-2443144.5003725)*86400.0;
            // check with tcg
            //        tcg4 = tt + 6.969290134e-10*(jdtcg - 2443144.5003725)*86400.0;  % AAS 05-352 (10) and IERS TN(104)
            //        [hrtemp, mintemp, sectemp] = sec2hms(tcg4 );
            //        jdtcg4 = jday(year, mon, day, hrtemp, mintemp, sectemp );
            //        tt2 = tcg4-6.969290134e-10*(jdtcg4-2443144.5003725)*86400.0;
            //        difchk = tt2-tt


            //tcbmtdb = -1.55051976772e-8 * (jdtai + jdtaifrac - 2443144.5003725) * 86400.0 - 6.55e-5;  // second, value for de405 AAS 05-352 (10) and IERS TN(104)?
            //tcb = tdb + tcbmtdb;
            //hms_sec(ref localhr, ref minute, ref second, Edirection.efrom, ref tcb);
            //jday(year, mon, day, hrtemp, mintemp, sectemp, out jdtcb, out jdtcbfrac);
            //ttcb = (jdtcb + jdtcbfrac - 2451545.0) / 36525.0;
            //fprintf(1,'     tcb %8.6f ttcb  %16.12f jdtcb  %18.11f %18.11f \n', tcb, ttcb, jdtcb, jdtcbfrac );
        }  // convtime

    }  //  class MathTimeLib

}  //  MathTimeMethods

