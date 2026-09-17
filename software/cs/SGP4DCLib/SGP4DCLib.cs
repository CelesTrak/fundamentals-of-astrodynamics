/* -----------------------------------------------------------------------------
*
*                                   SGP4DCLib.cs
*
*    this file contains the differential corrections routines using the sgp4
*    analytical propagation code.a detailed discussion of the theory and history
*    may be found in the 2008 aiaa paper by vallado and crawford.
*
*                                companion code for
*                      fundamentals of astrodynamics and applications
*                                      2022
*                                 by david vallado
*
*                       (w) dvallado@comspoc.com, email davallado@gmail.com
*
*    current :
*              11 nov 15  david vallado
*                           convert to c#
*    changes :
*              30 sep 15  david vallado
*                           fix jd, jdfrac
*               3 nov 14  david vallado
*                           update to msvs2013 c++
*               8 aug 14 alek lidtke
*                          change the constructor of std::vector< std::vector<double> > nut (lines 206 and 377), also other vectors in lines 802, 1267, and 2121
*                          started using fopen_s in line 1345, fscanf_s in line 1356 with MSVS compiler
*               6 aug 08  david vallado
*                           add operationmode for afspc (a) or improved (i)
*              18 jun 08  david vallado
*                           original version
* -------------------------------------------------------------------------- - */

using System;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

using MathTimeMethods;     // Edirection, globals
using EOPSPWMethods;       // EOPDataClass, SPWDataClass, iau80Class, iau06Class
using AstroLibMethods;     // EOpt, gravityConst, astroConst, xysdataClass, jpldedataClass
using SGP4Methods;


namespace SGP4DCMethods
{
    public class SGP4DCLib
    {

        public string SGP4DCLibVersion = "SGP4DCLib Version 2023-11-15";

    // setup the classes so methods can be called
    public MathTimeLib MathTimeLibr = new MathTimeLib();

    public EOPSPWLib EOPSPWLibr = new EOPSPWLib();

    public AstroLib AstroLibr = new AstroLib();

    public SGP4Lib SGP4Libr = new SGP4Lib();

    public StringBuilder strbuild = new StringBuilder();

        public class obsrec
        {
            public int sennum;
            public long satnum;
            public int year, mon, day, hr, min;
            public double jd, jdf, sec, dtmin, lst;
            public int error;
            public char init, method;
            public double[] rsecef = new double[3];
            public double[] vsecef = new double[3];
            public int obstype;
            public double x, y, z, xdot, ydot, zdot, bstar,
                rng, az, el, drng, daz, del,
                rtasc, decl, trtasc, tdecl;
        } //  obsrec

        public class senrec
        {
            public int sennum;
            public string senname;
            public double senlat, senlon, senalt,
                rngmin, rngmax, azmin, azmax, elmin, elmax,
                biasrng, biasaz, biasel, biasdrng, biasdaz, biasdel,
                biastrtasc, biastdecl,
                noisex, noisey, noisez, noisexdot, noiseydot, noisezdot, noisebstar,
                noiserng, noiseaz, noiseel, noisedrng, noisedaz, noisedel,
                noisetrtasc, noisetdecl;
        }  // senrec

        public static double twopi = 2.0 * Math.PI;

        public static double deg2rad = Math.PI / 180.0;  //   0.0174532925199433

        public static double xpdotp = 1440.0 / twopi;  // 229.1831180523293

        static double DMAX(double a, double b)
        {
            if (a > b) return a;
            else return b;
        }

        static float FMAX(float a, float b)
        {
            if (a > b) return a;
            else return b;
        }

        static int IMAX(int a, int b)
        {
            if (a > b) return a;
            else return b;
        }

        static double DMIN(double a, double b)
        {
            if (a < b) return a;
            else return b;
        }

        static float FMIN(float a, float b)
        {
            if (a < b) return a;
            else return b;
        }

        static int IMIN(int a, int b)
        {
            if (a < b) return a;
            else return b;
        }

        static double DSQR(double a)
        { return (a * a); }

        static float FSQR(float a)
        { return (a * a); }

        static double sum_f(double a, double b)
        { return a + b; }

        //#define SIGN(a, b)  ((b) >= 0 ? Math.Abs(a) : -Math.Abs(a))
        static double SIGN(double a, double b)
        {
            if (b >= 0)
                return Math.Abs(a);
            else
                return -Math.Abs(a);
        }

        //#define NO_DIFF(a, b) sum_f(a,b)==(b)


        static double dpythag(double a, double b)
        {
            double absa, absb;

            absa = Math.Abs(a);
            absb = Math.Abs(b);

            if (absa == 0.0)
                return absb;
            else if (absb == 0.0)
                return absa;
            else if (absa > absb)
                return absa * Math.Sqrt(1.0 + DSQR(absb / absa));
            else
                return absb * Math.Sqrt(1.0 + DSQR(absa / absb));
        }


        /* ==================================================================
        This is zero-offset (C style) version of NR routine for doing the
        SVD back-substitution action.

        Solve A*x = b for x in least-squares manner by using the SVD
        T
        decomposition of matrix A = U * w * V

        u[0..m-1, 0..n-1] these are constants from A by SVD decomposition.
        w[0..n-1]
        v[0..n-1, 0..n-1]

        b[0..m-1] is input
        x[0..n-1] is output

        m & n are dimensions, normally m > n for true (over determined)
        solution.
        ================================================================== */

        public void dsvbksb
            (
            double[,] u,
            double[,] w,
            double[,] v,
            int m, int n,
            double[,] b,
            out double[,] dx
            )
        {
            int jj, j, i;
            double s;
            double[,] tmp = new double[n, 1];
            dx = new double[n, 1];

            for (j = 0; j < n; j++)
            {
                s = 0.0;
                // Here we understand that if w[j] == 0 then treat X/0 as 0, not infinity.
                if (w[j, 0] == 0.0)
                {
                    for (i = 0; i < m; i++)
                        s += u[i, j] * b[i, 0];
                    s /= w[j, 0];
                }
                tmp[j, 0] = s;
            }

            for (j = 0; j < n; j++)
            {
                s = 0.0;
                for (jj = 0; jj < n; jj++)
                    s += v[j, jj] * tmp[jj, 0];
                dx[j, 0] = s;
            }
        } // dsvbksb


        /* -----------------------------------------------------------------------------
        *   Perform Math.Singular Value Decomposition of a matrix A to get:
        *
        *				 T
        *	A = U * w * V
        *
        *   This can then be used to robustly (in a Math.Singular-matrix sense) solve
        *   the least-squares problem:
        *
        *   A * x = b
        *
        *   With known A and b and dim b > dim x for an over-determined system by
        *   calling the dsvbksb() function with U, w and V. For near-Math.Singular cases
        *   you can zero and small w[] terms to prevent cancling near-infinities.
        *
        *   Here the calling arguments are:
        *
        *	A[0..m-1, 0..n-1]	input as A
        *
        *	m, n				dimensions of matrix A (and for w & V)
        *
        *	U[0..m-1, 0..n-1]	output
        *	w[0..n-1		output
        *	V[0..n-1, 0..n-1]	output
        *
        *   Return value is 0 if OK, or -1 if failed (rare).
        *
        --------------------------------------------------------------------------- */

        // -----------------------------------------------------------------------------
        public int dsvdcmp
            (
            ref double[,] a, int m, int n,
            out double[,] w, out double[,] v
            )
        {
            bool flag = false;
            int i, its, j, jj, k, l = 0, nm = 0;
            double anorm, c, f, g, h, s, scale, x, y, z;
            const int SVD_ITMAX = 30;
            //const double SVD_EPS = std::numeric_limits<double>::epsilon();
            const double SVD_EPS = 0.000000001;
            // a is already sized [m,n] by the caller (it's the accumulated ATWA matrix) - this used
            // to be `out double[,] a` with `a = new double[m - 1, n - 1]` here, which both (a) threw
            // away the caller's real matrix contents (an `out` parameter's incoming value is never
            // visible to the method - it decomposed a fresh all-zero array instead of the real ATWA),
            // and (b) under-allocated by one in each dimension while every loop below indexes up to
            // m-1/n-1 inclusive (0-based, e.g. `for (k = i; k < m; k++) ... a[k, i]`), which is what
            // actually threw the IndexOutOfRangeException - unrelated to statesize, this would crash
            // at any size since it's off by one relative to the loop bounds regardless of m/n's value.
            v = new double[n, n];
            w = new double[n, 1];
            double[] rv1 = new double[n];

            g = scale = anorm = 0.0;
            for (i = 0; i < n; i++)
            {
                l = i + 2;
                rv1[i] = scale * g;
                g = s = scale = 0.0;
                if (i < m)
                {
                    for (k = i; k < m; k++)
                        scale += Math.Abs(a[k, i]);
                    if (scale != 0.0)
                    {
                        for (k = i; k < m; k++)
                        {
                            a[k, i] /= scale;
                            s += a[k, i] * a[k, i];
                        }
                        f = a[i, i];

                        g = -SIGN(Math.Sqrt(s), f);
                        h = f * g - s;
                        a[i, i] = f - g;
                        for (j = l - 1; j < n; j++)
                        {
                            for (s = 0.0, k = i; k < m; k++) s += a[k, i] * a[k, j];
                            f = s / h;
                            for (k = i; k < m; k++) a[k, j] += f * a[k, i];
                        }
                        for (k = i; k < m; k++) a[k, i] *= scale;
                    }
                }
                w[i, 0] = scale * g;
                g = s = scale = 0.0;
                if (i + 1 <= m && i + 1 != n)
                {
                    for (k = l - 1; k < n; k++) scale += Math.Abs(a[i, k]);
                    if (scale != 0.0)
                    {
                        for (k = l - 1; k < n; k++)
                        {
                            a[i, k] /= scale;
                            s += a[i, k] * a[i, k];
                        }
                        f = a[i, l - 1];
                        g = -SIGN(Math.Sqrt(s), f);
                        h = f * g - s;
                        a[i, l - 1] = f - g;
                        for (k = l - 1; k < n; k++) rv1[k] = a[i, k] / h;
                        for (j = l - 1; j < m; j++)
                        {
                            for (s = 0.0, k = l - 1; k < n; k++) s += a[j, k] * a[i, k];
                            for (k = l - 1; k < n; k++) a[j, k] += s * rv1[k];
                        }
                        for (k = l - 1; k < n; k++) a[i, k] *= scale;
                    }
                }
                anorm = DMAX(anorm, (Math.Abs(w[i, 0]) + Math.Abs(rv1[i])));
            }
            for (i = n - 1; i >= 0; i--)
            {
                if (i < n - 1)
                {
                    if (g != 0.0)
                    {
                        for (j = l; j < n; j++)
                            v[j, i] = (a[i, j] / a[i, l]) / g;
                        for (j = l; j < n; j++)
                        {
                            for (s = 0.0, k = l; k < n; k++) s += a[i, k] * v[k, j];
                            for (k = l; k < n; k++) v[k, j] += s * v[k, i];
                        }
                    }
                    for (j = l; j < n; j++) v[i, j] = v[j, i] = 0.0;
                }
                v[i, i] = 1.0;
                g = rv1[i];
                l = i;
            }
            for (i = IMIN(m, n) - 1; i >= 0; i--)
            {
                l = i + 1;
                g = w[i, 0];
                for (j = l; j < n; j++) a[i, j] = 0.0;
                if (g != 0.0)
                {
                    g = 1.0 / g;
                    for (j = l; j < n; j++)
                    {
                        for (s = 0.0, k = l; k < m; k++) s += a[k, i] * a[k, j];
                        f = (s / a[i, i]) * g;
                        for (k = i; k < m; k++) a[k, j] += f * a[k, i];
                    }
                    for (j = i; j < m; j++) a[j, i] *= g;
                }
                else for (j = i; j < m; j++) a[j, i] = 0.0;
                ++a[i, i];
            }
            for (k = n - 1; k >= 0; k--)
            {
                for (its = 1; its <= SVD_ITMAX; its++)
                {
                    flag = true;
                    for (l = k; l >= 0; l--)
                    {
                        nm = l - 1;
                        if (Math.Abs(rv1[l]) <= SVD_EPS * anorm)
                        {
                            flag = false;
                            break;
                        }
                        if (Math.Abs(w[nm, 0]) <= SVD_EPS * anorm) break;
                    }
                    if (flag)
                    {
                        c = 0.0;
                        s = 1.0;
                        for (i = l; i < k + 1; i++)
                        {
                            f = s * rv1[i];
                            rv1[i] = c * rv1[i];
                            if (Math.Abs(f) <= SVD_EPS * anorm) break;
                            g = w[i, 0];
                            h = dpythag(f, g);
                            w[i, 0] = h;
                            h = 1.0 / h;
                            c = g * h;
                            s = -f * h;
                            for (j = 0; j < m; j++)
                            {
                                y = a[j, nm];
                                z = a[j, i];
                                a[j, nm] = y * c + z * s;
                                a[j, i] = z * c - y * s;
                            }
                        }
                    }
                    z = w[k, 0];
                    if (l == k)
                    {
                        if (z < 0.0)
                        {
                            w[k, 0] = -z;
                            for (j = 0; j < n; j++) v[j, k] = -v[j, k];
                        }
                        break;
                    }

                    if (its == SVD_ITMAX)
                    {
                        //strBuild.AppendLine("dsvdcmp: No convergence in " + SVD_ITMAX + " iterations");
                        return -1;
                    }

                    x = w[l, 0];
                    nm = k - 1;
                    y = w[nm, 0];
                    g = rv1[nm];
                    h = rv1[k];
                    f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
                    g = dpythag(f, 1.0);
                    f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
                    c = s = 1.0;
                    for (j = l; j <= nm; j++)
                    {
                        i = j + 1;
                        g = rv1[i];
                        y = w[i, 0];
                        h = s * g;
                        g = c * g;
                        z = dpythag(f, h);
                        rv1[j] = z;
                        c = f / z;
                        s = h / z;
                        f = x * c + g * s;
                        g = g * c - x * s;
                        h = y * s;
                        y *= c;
                        for (jj = 0; jj < n; jj++)
                        {
                            x = v[jj, j];
                            z = v[jj, i];
                            v[jj, j] = x * c + z * s;
                            v[jj, i] = z * c - x * s;
                        }
                        z = dpythag(f, h);
                        w[j, 0] = z;
                        // BUG FIX: was "if (z < 0.0)". The original algorithm renormalizes c/s by 1/z
                        // whenever z is nonzero (dpythag's result is a magnitude, so it's essentially
                        // always positive) - "< 0.0" meant this renormalization almost never ran, so
                        // the rotation applied to the U matrix (the loop right below) kept using stale
                        // c/s values left over from the V-matrix rotation above instead of z's own.
                        // U and V each stayed individually orthogonal (verified independently), but
                        // U*W*V^T stopped reconstructing the original matrix at all - confirmed by
                        // testing this exact logic standalone against a plain SPD test matrix outside
                        // this codebase before touching this file.
                        if (z != 0.0)
                        {
                            z = 1.0 / z;
                            c = f * z;
                            s = h * z;
                        }
                        f = c * g + s * y;
                        x = c * y - s * g;
                        for (jj = 0; jj < m; jj++)
                        {
                            y = a[jj, j];
                            z = a[jj, i];
                            a[jj, j] = y * c + z * s;
                            a[jj, i] = z * c - y * s;
                        }
                    }
                    rv1[l] = 0.0;
                    rv1[k] = f;
                    w[k, 0] = x;
                }
            }

            return 0;
        }  // dsvdcmp


        /* -----------------------------------------------------------------------------
        *                           procedure printtle
        *
        *  this procedure prints out the tle in a tle format.
        *
        *  author        : david vallado                  719-573-2600    6 aug 2008
        *
        *  inputs          description                    range / units
        *    satrec      - record of satellite parameters for TLE
        *
        *  outputs       :
        *    none.
        *
        *  locals        :
        *    none.
        *
        *  coupling      :
        *    none.
        * --------------------------------------------------------------------------- */

        public void printtle
            (
                SGP4Lib.elsetrec satrec
            )
        {
            StringBuilder strBuild = new StringBuilder();

            string longstr1, longstr1a, longstr2;
            string intldesg;  // classification
            string bstarstr, bstarexp;
            double rad;
            const double xpdotp = 1440.0 / (2.0 * Math.PI);  // 229.1831180523293
            rad = 180.0 / Math.PI;
            intldesg = "testsat";

            // -------------------- write out tle format of data ----------------------
            // may need to do some of these...
            // added %1c and classification here was N before, 
            // changed intldesg to 7 from 10 characters, added 3 spaces after %7s to align columns, 
            // changed epochdays to %11.8lf and added a leading 0
            //strBuild.AppendLine(longstr1 + "1 %5sU%8s  %02i%12.8lf  .00000000  00000-0 ",
            //    satrec.satnum, intldesg, satrec.epochyr - 2000, satrec.epochdays);

            //strBuild.AppendLine(longstr1a + "%8.4e 0   1\n", satrec.bstar);
            //ptr = strchr(longstr1a, 'e');
            //bstarexp = longstr1a[ptr - longstr1a + 3];
            //strBuild.AppendLine(longstr1a + "%8.4e 0   1\n", satrec.bstar * 1000.0);
            //bstarstr = longstr1a;

            // now do the second line
            if (satrec.inclo < 0.0)
                satrec.inclo = 2.0 * Math.PI + satrec.inclo;
            if (satrec.nodeo < 0.0)
                satrec.nodeo = 2.0 * Math.PI + satrec.nodeo;
            if (satrec.argpo < 0.0)
                satrec.argpo = 2.0 * Math.PI + satrec.argpo;
            if (satrec.mo < 0.0)
                satrec.mo = 2.0 * Math.PI + satrec.mo;
            // use the '0' in the format to get left filling zeros
            //strBuild.AppendLine(longstr2 + "2 %5s %08.4lf%9.4lf %07.0lf %8.4lf %8.4lf %11.8lf 00000\n",
            //    satrec.satnum, satrec.inclo * rad, satrec.nodeo * rad,
            //    satrec.ecco * 10000000.0, satrec.argpo * rad, satrec.mo * rad, satrec.no_unkozai * xpdotp);

            //     printf("%s%s%s", longstr1,bstarstr, bstarexp);
            //strBuild.AppendLine(longstr1 + longstr1a);
            //strBuild.AppendLine(longstr2);
        }  //  printtle


        /* -----------------------------------------------------------------------------
        *                           procedure getsensorparams
        *
        *  this procedure gets the sensor parameters. note that the values in here are
        *  arbirtrary at this point, but may be filled in and added as appropriate.
        *
        *  author        : david vallado                  719-573-2600    1 dec 2007
        *
        *  inputs          description                    range / units
        *    sennum      - sensor number
        *
        *  outputs       :
        *    currsenrec  - structure containing sensor information
        *
        *  locals        :
        *    none.
        *
        *  coupling      :
        *    none.
        * --------------------------------------------------------------------------- */

        public void getsensorparams
            (
            int sennum,
            out senrec currsenrec
            )
        {
            double rad = 180.0 / Math.PI;

            currsenrec = new senrec();

            switch (sennum)
            {
                case 1:
                    {
                        ///*
                        currsenrec.noisex = 0.01;    // 10 m
                        currsenrec.noisey = 0.01;
                        currsenrec.noisez = 0.01;
                        currsenrec.noisexdot = 0.01;  // 10 cm/s
                        currsenrec.noiseydot = 0.01;
                        currsenrec.noisezdot = 0.01;
                        currsenrec.noisebstar = 0.0001;
                        //*/
                        /*
                        currsenrec.noisex = 1.0;    // 1 km
                        currsenrec.noisey = 1.0;
                        currsenrec.noisez = 1.0;
                        currsenrec.noisexdot = 0.001;  // 1 m/s
                        currsenrec.noiseydot = 0.001;
                        currsenrec.noisezdot = 0.001;
                        currsenrec.noisebstar = 1.0;
                        */
                    }
                    break;
                case 2:
                    {
                        currsenrec.noisex = 0.001;
                    }
                    break;
                case 3:
                    {
                        currsenrec.noisex = 0.001;
                    }
                    break;
                case 344:
                    {
                        currsenrec.sennum = 344;
                        //strcpy_s(currsenrec.senname, "fylingdales");
                        currsenrec.senlat = 54.37 / rad;
                        currsenrec.senlon = -0.67 / rad;
                        currsenrec.senalt = 0.3389;  // km
                        currsenrec.noiserng = 0.09;   // km
                        currsenrec.noiseaz = 0.02 / rad;  // rad
                        currsenrec.noiseel = 0.01 / rad;  // rad
                    }
                    break;
                case 932:
                    {
                        currsenrec.sennum = 932;
                        //strcpy_s(currsenrec.senname, "kaena point");
                        // updated to the exact AAS-paper values used in Ex10_4.m / geos6a.inp
                        // (was 21.572 / 201.733 / 0.3005 with 1.0 placeholder noise values)
                        currsenrec.senlat = 21.572056 / rad;
                        currsenrec.senlon = -158.266578 / rad;
                        currsenrec.senalt = 0.3002;  // km

                        currsenrec.noiserng = 0.0925;   // km
                        currsenrec.noiseaz = 0.0224 / rad;  // rad
                        currsenrec.noiseel = 0.0139 / rad;  // rad
                    }
                    break;
            }  // case

        } // getsensorparams


        /* -----------------------------------------------------------------------------
        *                           procedure state2satrec
        *
        *  this procedure converts the state to the satrec structure, and back. be careful 
        *  of the units as the variables are passed back and forth. 
        *
        *  author        : david vallado                  719-573-2600   15 jan 2008
        *
        *  inputs          description                    range / units
        *    xnom        - state vector                   varied
        *    direct      - direction of conversion        eTo, eFrom
        *    satrec      - structure of satellite parameters for TLE
        *
        *  outputs       :
        *    satrec      - structure of satellite parameters for TLE
        *    xnom        - state vector                   varied
        *
        *  locals        :
        *    rnom        - nom position vector at epoch   km
        *    vnom        - nom velocity vector at epoch   km/s
        *
        *  coupling      :
        *    none.
        *
        *  references    :
        * --------------------------------------------------------------------------- */

        public void state2satrec
            (
            ref double[] xnom,
            double jdepoch, double jdepochf,
            char statetype, int statesize,
            AstroLib.EOpt opt,
            Enum direct,
            ref SGP4Lib.elsetrec satrec
            )
        {
            double[] rnom = new double[3];
            double[] vnom = new double[3];
            double p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper;
            int month, day, hr, minute;
            double seconds;
            SGP4Lib.elsetrec satrecorig;


            satrecorig = satrec;
            if (direct.Equals(MathTimeLib.Edirection.eto))
            {
                switch (statetype)
                {
                    case 'v':
                        rnom[0] = xnom[0];
                        rnom[1] = xnom[1];
                        rnom[2] = xnom[2];
                        vnom[0] = xnom[3];
                        vnom[1] = xnom[4];
                        vnom[2] = xnom[5];

                        AstroLibr.rv2coe(rnom, vnom, out p, out a, out ecc, out incl, out raan, out argp, out nu, out m,
                                        out arglat, out truelon, out lonper);  // mu should be 3.986008e5,
                        satrec.no_unkozai = 60.0 * Math.Sqrt(3.986008e5 / (a * a * a));  // rad per min, sgp4 wgs-72 mu value
                        satrec.a = a / 6378.135;  // er
                        satrec.ecco = ecc;
                        satrec.inclo = incl;      // rad
                        satrec.nodeo = raan;
                        satrec.argpo = argp;
                        satrec.mo = m;
                        break;
                    case 't':  // tle (keplerian) elements
                        //  satrec.no    = xnom[0];  // rad/min
                        satrec.a = xnom[0];  // er
                        satrec.no_unkozai = 1.0 / 13.446839 * Math.Sqrt(1.0 / (satrec.a * satrec.a * satrec.a));  // rad / min
                        satrec.ecco = xnom[1];
                        satrec.inclo = xnom[2];
                        satrec.nodeo = xnom[3];
                        satrec.argpo = xnom[4];
                        satrec.mo = xnom[5];
                        break;
                    case 'e':  // equinoctial elements
                        // satrec.no    = xnom[2];    // rad / min
                        // satrec.a     = pow( 1.0 / (satrec.no * satrec.no * 13.446839 * 13.446839) , 1.0 / 3.0 );  // er
                        satrec.a = xnom[2];    // er
                        if (satrec.a < 1.0)  // can't be less than the radius of the earth
                        {
                            strbuild.AppendLine("changed a " + satrec.a);
                            satrec.a = satrecorig.a * 1.01;   // er   1.05, .9
                            xnom[2] = satrec.a;  // remember to change both!!!
                            strbuild.AppendLine(@" to satrec.a \n");
                        }
                        satrec.no_unkozai = (1.0 / 13.446839) * Math.Sqrt(1.0 / (satrec.a * satrec.a * satrec.a));  // rad / min
                        if (satrec.no_unkozai < 1.0e-5)
                        {
                            strbuild.AppendLine("changed no " + satrec.no_kozai);
                            satrec.no_unkozai = satrecorig.no_unkozai * 0.9;
                            satrec.a = Math.Pow(1.0 / (satrec.no_unkozai * satrec.no_unkozai * 13.446839 * 13.446839), 1.0 / 3.0);  // er
                            xnom[2] = satrec.a;  // remember to change both!!!
                            strbuild.AppendLine(" to " + satrec.no_unkozai + @" \n");
                        }
                        satrec.ecco = (Math.Sqrt(xnom[0] * xnom[0] + xnom[1] * xnom[1]));
                        if (satrec.ecco > 1.0)
                        {
                            strbuild.AppendLine("changed ecco " + satrec.ecco);
                            satrec.ecco = satrecorig.ecco * 0.9;
                            // be careful of order here...
                            // could iterate, but Math.Since we're just changing values, it's not needed
                            xnom[0] = (satrec.ecco * Math.Cos(satrec.argpo + satrec.nodeo));    // ke, af
                            xnom[1] = (satrec.ecco * Math.Sin(satrec.argpo + satrec.nodeo));    // he, ag
                            strbuild.AppendLine(" to " + satrec.ecco + @" \n");
                        }
                        satrec.inclo = 2.0 * Math.Atan(Math.Sqrt(xnom[4] * xnom[4] + xnom[5] * xnom[5]));
                        satrec.nodeo = Math.Atan2(xnom[4], xnom[5]);

                        // make sure nodeo stays within 0  to 2PI as it gets used without trigonmetric functions
                        if (satrec.nodeo < 0.0)
                        {
                            satrec.nodeo = twopi + satrec.nodeo;
                            xnom[4] = (Math.Tan(satrec.inclo * 0.5) * Math.Sin(satrec.nodeo));  // pe
                            xnom[5] = (Math.Tan(satrec.inclo * 0.5) * Math.Cos(satrec.nodeo));  // qe
                        }

                        satrec.argpo = Math.Atan2(xnom[1], xnom[0]) - satrec.nodeo;
                        satrec.mo = xnom[3] - satrec.nodeo - satrec.argpo;
                        // fmod doesn't seem to make a diff in the results
                        satrec.mo = satrec.mo % twopi;
                        break;
                } // case

                // set misc variables
                //whichconst = wgs72;
                //opsmode = 'a';
                satrec.satnumStr = "88888";
                satrec.nddot = 0.00000e0;
                satrec.bstar = 0.001;
                satrec.ndot = 0.000001;
                satrec.elnum = 1;
                satrec.revnum = 100;
                satrec.classification = 'U';
                satrec.intldesg = "          ";
                satrec.ephtype = 0;

                // no_kozai is the converted mean motion in revs/day in the TLE
                // no_unkozai is the converted semimajor axis in er
                //satrec.no_kozai = satrec.no_unkozai;
                double ak, adel, d1, del, del1, tumin, xke, xpdotp, j2x, cosio2;
                xke = 60.0 / Math.Sqrt(AstroLibr.astroConsts.re * AstroLibr.astroConsts.re * AstroLibr.astroConsts.re
                    / AstroLibr.astroConsts.mu);  // min / tu
                tumin = 1.0 / xke;   // tu / min
                j2x = 0.001082616;
                xpdotp = 1440.0 / (2.0 * Math.PI);   // 229.1831180523293   rad / min ?
                                                     // assume for 1st guess that no_kozai = no_unkozai
                ak = Math.Pow(xke / satrec.no_unkozai, 2.0 / 3.0);          // er
                cosio2 = Math.Cos(satrec.inclo) * Math.Cos(satrec.inclo);
                d1 = 0.75 * j2x * (3.0 * cosio2 - 1.0) / (Math.Pow(1.0 - satrec.ecco * satrec.ecco, 1.5));  // same both ways
                del1 = d1 / (ak * ak);
                adel = ak * (1.0 - del1 * del1 - del1 * (1.0 / 3.0 + 134.0 * del1 * del1 / 81.0));
                del = d1 / (adel * adel);
                satrec.no_kozai = satrec.no_unkozai * (1.0 + del);

                MathTimeLibr.invjday(satrec.jdsatepoch, satrec.jdsatepochF, out satrec.epochyr, out month, out day,
                    out hr, out minute, out seconds);
                Int32 days;
                MathTimeLibr.findDays(satrec.epochyr, month, day, hr, minute, seconds, out days);
                satrec.epochdays = days;

                // convert units 
                //satrec.no_kozai = satrec.no_kozai / xpdotp; //* rev / day to rad/min
                satrec.ndot = satrec.ndot / (xpdotp * 1440.0);  //* ? * minperday
                satrec.nddot = satrec.nddot / (xpdotp * 1440.0 * 1440);
                //satrec.inclo = satrec.inclo  * deg2rad;   // deg to rad if needed
                //satrec.nodeo = satrec.nodeo  * deg2rad;
                //satrec.argpo = satrec.argpo  * deg2rad;
                //satrec.mo = satrec.mo     * deg2rad;
                if (statesize > 6)
                    satrec.bstar = xnom[6];
                satrec.jdsatepoch = jdepoch;
                satrec.jdsatepochF = jdepochf;
            }
            else  // find xnom from satrec -----------------------------------------
            {
                switch (statetype)
                {
                    case 'v':  // vectors
                               // xnom contains the nominal pos and vel vector
                        break;
                    case 't':  // tle (keplerian) elements
                               //  xnom[0] = satrec.no;
                        satrec.a = Math.Pow(satrec.no_unkozai * satrec.tumin, (-2.0 / 3.0));
                        xnom[0] = satrec.a;  // er
                        xnom[1] = satrec.ecco;
                        xnom[2] = satrec.inclo;
                        xnom[3] = satrec.nodeo;
                        xnom[4] = satrec.argpo;
                        xnom[5] = satrec.mo;
                        break;
                    case 'e':  // equinoctial elements
                        xnom[0] = (satrec.ecco * Math.Cos(satrec.argpo + satrec.nodeo));    // ke, af
                        xnom[1] = (satrec.ecco * Math.Sin(satrec.argpo + satrec.nodeo));    // he, ag
                        xnom[2] = satrec.a;  // a is in er
                        xnom[3] = (satrec.mo + satrec.argpo + satrec.nodeo) % twopi;    // L
                        xnom[4] = (Math.Tan(satrec.inclo * 0.5) * Math.Sin(satrec.nodeo));  // pe
                        xnom[5] = (Math.Tan(satrec.inclo * 0.5) * Math.Cos(satrec.nodeo));  // qe
                        break;
                } // case
                if (statesize > 6)
                    xnom[6] = satrec.bstar;

            }
        } // state2satrec


        /* -----------------------------------------------------------------------------
        *
        *                           procedure findatwbatwa
        *
        * this procedure finds the a and b matrices for the differential correction
        *   problem.  remember that it isn't critical for the propagations to use
        *   the highest fidelity techniques because we're only trying to find the
        *   "slope". k is an index that allows us to do multiple rows at once. it's
        *   used for both the b and a matrix calculations.
        *
        *  algorithm     : find the a and b matrices by accumulation to reduce matrix
        *                  sizes calculate the matrix combinations
        *                  atw is found without matrix operations to avoid large matrices
        *
        *  author        : david vallado                  719-573-2600    6 aug 2008
        *
        *  inputs          description                    range / units
        *    firstob     - number of observations
        *    lastob      - number of observations
        *    statesize   - size of state                  6 , 7
        *    percentchg  - amount to modify the vectors
        *                  by in finite differencing
        *    deltaamtchg - tolerance for small value in
        *                  finite differencing            0.0000001
        *    whichconst  - parameter for sgp4 constants   wgs72, wgs721, wgs84
        *    satrec      - structure of satellite parameters for TLE
        *    obsrecfile    - array of records containing:
        *                  senum, jd, rsvec, obstype,
        *                  rng, az, el, drng, daz, del,
        *                  trtasc, tdecl data
        *    statetype   - type of elements (equinoctial, etc)  'e', 't'
        *    proptype    - type of propagation            's', 'n', 'k'
        *    xnom        - state vector                   varied
        *
        *  outputs       :
        *    atwa        - atwa matrix
        *    atwb        - atwb matrix
        *    atw         - atw matrix
        *    b           - b matrix, the residuals
        *    drng2       - range residual squared
        *    daz2        - azimuth residual squared
        *    del2        - elevation residual squared
        *    ddrng2      - range rate residual squared
        *    ddaz2       - azimuth rate residual squared
        *    ddel2       - elevation rate residual squared
        *    dtrtasc2    - topocentric right ascension residual squared
        *    dtdecl2     - topocentric declination residual squared
        *    dx2         - x position residual squared
        *    dy2         - y position residual squared
        *    dz2         - z position residual squared
        *    dxdot2      - xdot position residual squared
        *    dydot2      - ydot position residual squared
        *    dzdot2      - zdot position residual squared
        *
        *  locals        :
        *    rnom        - nom position vector at epoch   km
        *    vnom        - nom velocity vector at epoch   km/s
        *    a           - a matrix
        *    indobs      -
        *    at          -
        *    w1          -
        *    w2          -
        *    w3          -
        *    lst         -
        *    gst         -
        *    dtsec        -
        *    deltaamt    -
        *    rngpert     -
        *    azpert      - modified azimuth               -2Math.PI to 2Math.PI
        *    elpert      - modified azimuth               -Math.PI/2 to Math.PI/2
        *    drng        -
        *    daz         -
        *    del         -
        *    error       -
        *    i, j, k     -
        *
        *  coupling      :
        *    findsenptr  - find sensor data
        *    rv_razel    - find r and v given range, az, el, and rates
        *    rv_tradec   - find r and v given topocentric rtasc and decl
        *
        *  references    :
        *    vallado       2007, 753-765
        * --------------------------------------------------------------------------- */

        public void findatwaatwb
        (
            int firstob, int lastob, int statesize,
            double percentchg, double deltaamtchg, SGP4Lib.gravconsttype whichconst, char interp, double jdeopstart,
            double jd, double jdf, ref SGP4Lib.elsetrec satrec, SGP4DCLib.obsrec[] obsrecarr,
            char statetype, char proptype, double[] xnom, double rmag, int derivType, double cdam, double cram,
            AstroLib.EOpt opt, EOPSPWLib.iau80Class iau80arr, EOPSPWLib.iau06Class iau06arr,
            double[] eoparr, Int32 mjdeopstart, AstroLib.jpldedataClass[] jpldearr, double jdjpldestart,
            out double drng2, out double daz2, out double del2,
            out double ddrng2, out double ddaz2, out double ddel2,
            out double dtrtasc2, out double dtdecl2,
            out double dx2, out double dy2, out double dz2,
            out double dxdot2, out double dydot2, out double dzdot2,
            out double[,] atwa, out double[,] atwb, out double[,] atw, out double[,] b
        )
        {
            obsrec currobsrec;
            senrec currsenrec;
            drng2 = daz2 = del2 = ddrng2 = ddaz2 = ddel2 = dtrtasc2 = dtdecl2 = dx2 = dy2 = dz2 = dxdot2 = dydot2
                = dzdot2 = 0.0;

            // ------------ zero out and dimension matrices
            atwa = new double[statesize, statesize];
            b = new double[statesize, 1];
            atwb = new double[statesize, 1];
            atw = new double[statesize, statesize];

            int rowc, colc, r, c, i, j, ii; //k, sennum;
            double[] rs = new double[3];
            double[] r3 = new double[3];
            double[] v3 = new double[3];
            double[] rteme = new double[3];
            double[] vteme = new double[3];
            double[] ateme = new double[3];
            double[] ritrf = new double[3];
            double[] vitrf = new double[3];
            double[] aitrf = new double[3];
            double[] rnom = new double[3];
            double[] vnom = new double[3];
            double[] recinom = new double[3];
            double[] vecinom = new double[3];
            double[] recipert = new double[3];
            double[] vecipert = new double[3];
            double[] rpert = new double[3];
            double[] vpert = new double[3];
            double[] aeci = new double[3];
            double dut1, tut1, jdut1, jdtt, jdftt, ttt, lod, xp, yp, ddpsi, ddeps, ddx, ddy, jdxysstart,
                dtt, dayconv;  // jde
            int dat;

            int indobs = 7;  // size larger than most number in the state

            double[,] a = new double[indobs, statesize];
            double[,] at = new double[statesize, indobs];
            double[,] atwaacc = new double[statesize, statesize];
            double[,] atwbacc = new double[statesize, 1];
            double[,] X = new double[7, 1];
            double[,] XCurrn = new double[7, 1];
            double[,] XCurrp = new double[7, 1];
            double[,] Xf = new double[7, 1];
            double[,] xpert = new double[7, 1];
            double[,] trans = new double[3, 3];
            double deltaamt;
            double w1, w2, w3, w4, w5, w6, w7,
                rngnom, aznom, elnom, drngnom, daznom, delnom, rngpert, azpert, elpert, drngpert, dazpert, delpert,
                trtascnom, tdeclnom, dtrtascnom, dtdeclnom, trtascpert, tdeclpert, dtrtascpert, dtdeclpert;  // lst, gst, dtsec,
            double weight, bstarnom, bstarpert, rng2, dbstar2;

            // --------------------- initialize parameters ------------------
            jdtt = 0.0;
            jdftt = 0.0;
            dat = 0;
            dayconv = 1.0 / 86400.0;
            double jdyconv = 1.0 / 36525.0;

            for (ii = 0; ii < 6; ii++)
            {
                X[ii, 0] = xnom[ii];
            }

            currobsrec = obsrecarr[0];
            switch (currobsrec.obstype)
            {
                // range
                case 0:
                    indobs = 1;
                    break;
                // az-el
                case 1:
                    indobs = 2;
                    break;
                // rng t rtasc decl
                case 3:
                    indobs = 2;
                    break;
                // rng az el
                case 2:
                    indobs = 3;
                    break;
                // pos/vel vectors
                case 4:
                    indobs = statesize;  // add one for bstar possibility
                    break;
            }  // case

            // ------------- reset these since they will accumulate ---------
            ddpsi = ddeps = 0.0;
            rngnom = trtascnom = tdeclnom = drngnom = dtrtascnom = dtdeclnom = 0.0;
            rngpert = trtascpert = tdeclpert = drngpert = dtrtascpert = dtdeclpert = 0.0;
            aznom = elnom = rngnom = daznom = delnom = bstarnom = 0.0;
            azpert = elpert = rngpert = bstarpert = dazpert = delpert = 0.0;
            w1 = w2 = w3 = w4 = w5 = w6 = w7 = 0;

            SGP4Lib.elsetrec satrecp;

            dbstar2 = 0.0;
            rng2 = 0.0;

            //for (r = 0; r < statesize; r++)
            //{
            //    for (c = 0; c < statesize; c++)
            //        atwa[r, c] = 0.0;
            //    atwb[r, 0] = 0.0;
            //}

            // ------------------- loop through all the observations ------------------
            for (i = firstob; i <= lastob; i++)
            {
                currobsrec = obsrecarr[i];
                //printf( "ob %2i rsecef0 %11.5f jd %11.5f dtmin %8.3f hr %3i rng %11.7f \n",
                //         i, currobsrec.rsecef[0], currobsrec.jd, currobsrec.dtmin, currobsrec.hr, currobsrec.rng );  // ritrf
                //printf( "ob %2i jd %11.5f dtmin %9.5f hr %3i currobs %11.7f  %11.7f  %11.7f \n",
                //	i, currobsrec.jd, currobsrec.dtmin, currobsrec.hr, currobsrec.x, currobsrec.y, currobsrec.z);  // ritrf

                // --------- propagate the nominal vector to the epoch time -----------
                // -------------- and find the nominal observations -------------------
                // use IAU-76/FK5. so iau06 not needed
                //AstroLib.EOpt opt = AstroLib.EOpt.e80;
                jdxysstart = 0.0;
                satrecp = null;
                // minutes from the first observation to this observation's time - computed directly from
                // jd/jdf (like the Kepler branch's dtsec below) rather than trusting a caller-populated
                // currobsrec.dtmin, which nothing currently sets (defaults to 0.0) and was silently making
                // every SGP4 propagation evaluate at tsince=0 regardless of the observation's real time.
                double tsinceMin = (currobsrec.jd + currobsrec.jdf - obsrecarr[1].jd - obsrecarr[1].jdf) * 1440.0;
                switch (proptype)
                {
                    // dtmin is the minutes from epoch (1st obs in .e file) to currobs time
                    case 's':  // sgp4
                        SGP4Libr.sgp4(ref satrec, tsinceMin, rteme, vteme);
                        EOPSPWLibr.findeopparam(jd, jdf, 's', EOPSPWLibr.eopdata,
                            out dut1, out dat, out lod, out xp, out yp, out ddpsi, out ddeps, out ddx, out ddy);
                        jdftt = jdf + (dat + 32.184);   // sec 
                        ttt = (jd + jdf + (dat + 32.184 + jdftt) * dayconv - 2451545.0) * jdyconv;

                        AstroLibr.teme_eci(ref rteme, ref vteme, iau80arr, MathTimeLib.Edirection.eto, ttt, 0.0, 0.0,
                            ref recinom, ref vecinom);
                        break;
                    case 'a': // semianalytical
                              //for (ii = 0; ii < 3; ii++)
                              //{
                              //	rnom[ii] = X[ii, 0];
                              //	vnom[ii] = X[ii + 3, 0];
                              //}
                              //tsince = currobsrec.dtmin*60.0;  // should be 0.0 if using first point in .e file
                              //fullsun = 'y';
                              //astPert::srpperts(currobsrec.jd, currobsrec.jdf, rnom, vnom, interp, opt, fullsun, tsince, 
                              //	jpldearr, jdjpldestart, tut1, cram, recinom, vecinom, rmean1, vmean1);
                              //for (ii = 0; ii < 3; ii++)
                              //{
                              //	Xf[ii, 0] = recinom[ii];
                              //	Xf[ii + 3, 0] = vecinom[ii];
                              //}
                        break;
                    case 'n':  // numerical 
                               //// save first time through to each obs to save prop time for a matrix calcs
                               //// no, this needs to have the same perts from the start time, not form part way through
                               ////if (i == firstob)
                               ////{
                               ////	// reset to be sure it's the nominal value
                               //	for (ii = 0; ii < 6; ii++)
                               //		XCurrn[ii, 0] = X[ii, 0];
                               //	tsince = currobsrec.dtmin*60.0;
                               ////}
                               //// assumes all the obs are sequential
                               ////tsince = (currobsrec.dtmin - oldtn) * 60.0;  // sec

                        //astPert::rk4(currobsrec.jd, currobsrec.jdf, tsince, XCurrn, rmag, derivType, cdam, cram, interp, opt, 
                        //	jpldearr, jdjpldestart, iau80rec, jdeopstart, eoparr, Xf);
                        //   // update for next time step for speed 
                        //oldtn = currobsrec.dtmin; // min
                        //for (ii = 0; ii < 3; ii++)
                        //{
                        //	recinom[ii] = Xf[ii, 0];
                        //	vecinom[ii] = Xf[ii + 3, 0];
                        //}
                        break;
                    default:  // Kepler 2-body
                        double dtsec = (currobsrec.jd + currobsrec.jdf - obsrecarr[1].jd - obsrecarr[1].jdf) * 86400.0;  // s
                        for (ii = 0; ii < 3; ii++)
                        {
                            rnom[ii] = X[ii, 0];
                            vnom[ii] = X[ii + 3, 0];
                        }
                        AstroLibr.kepler(rnom, vnom, dtsec, out recinom, out vecinom, out string keplerOutTextNom);
                        for (ii = 0; ii < 3; ii++)
                        {
                            Xf[ii, 0] = recinom[ii];
                            Xf[ii + 3, 0] = vecinom[ii];
                        }
                        break;
                }  // case proptype

                // eci to itrf if observation type
                // change obs vector state to ecef so obs are created in the right coord system
                if (currobsrec.obstype != 4)
                {
                    dtt = (currobsrec.hr * 60 + currobsrec.min + currobsrec.sec / 60.0) * dayconv;
                    EOPSPWLibr.findeopparam(currobsrec.jd, dtt, 's', EOPSPWLibr.eopdata,
                        out dut1, out dat, out lod, out xp, out yp, out ddpsi, out ddeps, out ddx, out ddy);
                    // BUG FIX: this used to be "jdut1 = currobsrec.jdf + dut1" - missing currobsrec.jd
                    // (the ~2.45 million "big" part of the Julian date) entirely, and adding dut1 (in
                    // SECONDS) directly to a fractional-day value without converting to days first. Since
                    // jdut1 alone (not tut1) is what actually feeds eci_ecef's sidereal-time calculation,
                    // this was producing a wildly wrong Earth-rotation angle (~126 degrees off, confirmed
                    // by comparing against an independently-computed GMST for this exact epoch).
                    jdut1 = currobsrec.jd + currobsrec.jdf + dut1 / 86400.0;
                    tut1 = (jdut1 - 2451545.0) * jdyconv;
                    jdtt = currobsrec.jd;
                    jdftt = jdf + dat + 32.184;   // sec 
                    ttt = (currobsrec.jd + jdf + (dat + 32.184 + jdftt) * dayconv - 2451545.0) * jdyconv;
                    // eci_ecef takes (reci,veci,aeci, direct, recef,vecef,aecef, iau80arr, ttt,jdut1,lod,xp,yp,ddpsi,ddeps).
                    // converting the nominal ECI state (recinom/vecinom) to ECEF (ritrf/vitrf) to compare
                    // against the ground-station observations is ECI->ECEF, i.e. Edirection.eto (confirmed
                    // against AstroLib.eci_ecef's own direct.Equals(Edirection.eto) branch). aeci/aitrf are
                    // unused zero-filled scratch accel vectors (two-body/Kepler mode has no real accel here).
                    AstroLibr.eci_ecef(ref recinom, ref vecinom, ref aeci, MathTimeLib.Edirection.eto,
                        ref ritrf, ref vitrf, ref aitrf, iau80arr, ttt, jdut1, lod, xp, yp, ddpsi, ddeps);
                } // if obstype

                // ----------------- determine sensor characteristics -----------------
                rs[0] = currobsrec.rsecef[0];
                rs[1] = currobsrec.rsecef[1];
                rs[2] = currobsrec.rsecef[2];

                getsensorparams(currobsrec.sennum, out currsenrec);

                // ------------------------- find b matrix ----------------------------
                if (currobsrec.obstype == 3)
                    AstroLibr.rv_tradec(ref ritrf, ref vitrf, rs, MathTimeLib.Edirection.eto,
                        ref rngnom, ref trtascnom, ref tdeclnom, ref drngnom, ref dtrtascnom, ref dtdeclnom);
                else
                if (currobsrec.obstype == 4)
                    bstarnom = satrec.bstar;
                else
                    AstroLibr.rv_razel(ref ritrf, ref vitrf, currsenrec.senlat, currsenrec.senlon, currsenrec.senalt,
                        MathTimeLib.Edirection.eto,
                        ref rngnom, ref aznom, ref elnom, ref drngnom, ref daznom, ref delnom);

                switch (currobsrec.obstype)
                {
                    case 0:
                        b[0, 0] = currobsrec.rng - rngnom;
                        break;
                    case 1:
                        b[0, 0] = currobsrec.az - aznom;
                        //fix for 0-360...
                        if (Math.Abs(b[0, 0]) > Math.PI)
                            b[0, 0] = b[0, 0] - Math.Sign(b[0, 0]) * 2.0 * Math.PI;
                        b[1, 0] = currobsrec.el - elnom;
                        break;
                    case 2:
                        b[0, 0] = currobsrec.rng - rngnom;
                        b[1, 0] = currobsrec.az - aznom;
                        // fix for 0-360...
                        if (Math.Abs(b[1, 0]) > Math.PI)
                            b[1, 0] = b[1, 0] - Math.Sign(b[1, 0]) * 2.0 * Math.PI;
                        b[2, 0] = currobsrec.el - elnom;
                        break;
                    case 3:
                        b[0, 0] = currobsrec.trtasc - trtascnom;
                        // fix for 0-360...
                        if (Math.Abs(b[0, 0]) > Math.PI)
                            b[0, 0] = b[0, 0] - Math.Sign(b[0, 0]) * 2.0 * Math.PI;
                        b[1, 0] = currobsrec.tdecl - tdeclnom;
                        break;
                    case 4:
                        b[0, 0] = currobsrec.x - recinom[0];
                        b[1, 0] = currobsrec.y - recinom[1];
                        b[2, 0] = currobsrec.z - recinom[2];
                        b[3, 0] = currobsrec.xdot - vecinom[0];
                        b[4, 0] = currobsrec.ydot - vecinom[1];
                        b[5, 0] = currobsrec.zdot - vecinom[2];
                        // BUG FIX: b is sized [indobs,1], and indobs==statesize for obstype 4 (6 when
                        // bstar isn't being solved for, 7 when it is) - writing b[6,0] unconditionally
                        // threw IndexOutOfRangeException for any 6-state (position+velocity only) fit,
                        // which is exactly what processxyz/processtle use since state2satrec's 'v' case
                        // doesn't sync bstar from xnom anyway (see the statesize note where those are
                        // built). Only write the 7th residual when there's actually a 7th state to fit.
                        if (statesize == 7)
                            b[6, 0] = currobsrec.bstar - bstarnom;
                        break;
                }  // case

                //printf( "rnom %11.5f %11.5f %11.5f %8.3f %8.3f %8.3f %8.3f \n",
                //              rteme[0], rteme[1], rteme[2], rngnom, aznom*rad, elnom*rad, currobsrec.rng );  // ritrf
                //printf("recinom %11.5f %11.5f %11.5f curr %8.3f %8.3f %8.3f \n",
                //	recinom[0], recinom[1], recinom[2], currobsrec.x, currobsrec.y, currobsrec.z);  // ritrf


                // ------------------------ find a matrix -----------------------------
                // ------------- reset the perturbed vector to the nominal ------------
                if (proptype == 's')
                    satrecp = satrec;
                drng2 = daz2 = del2 = dtrtasc2 = dtdecl2 = dx2 = dy2 = dz2 = dxdot2 = dydot2 = dzdot2 = 0.0;
                // ----- perturb each element in the state individually (elements or vectors) ------
                for (j = 0; j < statesize; j++)
                {
                    // use the current nominal propagated results from previous section for numerical
                    // turn off for 2 body???????
                    //if (proptype == 'n')
                    //{
                    // initialize the entire vector
                    // finitediff changes only 1 component
                    for (ii = 0; ii < 6; ii++)
                        xpert[ii, 0] = xnom[ii];
                    //}

                    finitediff(whichconst, j, percentchg, deltaamtchg, statetype, statesize, proptype, ref satrecp,
                        opt, xnom, jd, jdf, out xpert, out deltaamt);

                    // --------- propagate the perturbed vector to the epoch time -----------
                    // ---------------- and find the perturbed observations -----------------
                    //AstroLib.EOpt opt = AstroLib.EOpt.e80;
                    switch (proptype)
                    {
                        case 's':  // sgp4
                            SGP4Libr.sgp4(ref satrecp, tsinceMin, r3, v3);
                            jdftt = jdf + dat + 32.184;   // sec
                            ttt = (jd + jdf + (dat + 32.184 + jdftt) * dayconv - 2451545.0) * jdyconv;
                            AstroLibr.teme_eci(ref r3, ref v3, iau80arr, MathTimeLib.Edirection.eto, ttt, ddpsi, ddeps,
                                ref recipert, ref vecipert);
                            break;
                        case 'a': // semianalytical
                            /*	for (ii = 0; ii < 3; ii++)
                                {
                                    rpert[ii] = xpert[ii, 0];
                                    vpert[ii] = xpert[ii + 3, 0];
                                }
                                tMath.Since = currobsrec.dtmin*60.0;
                                fullsun = 'y';
                                astPert::srpperts(currobsrec.jd, currobsrec.jdf, rpert, vpert, interp, opt, fullsun, tMath.Since, jpldearr, jdjpldestart, tut1, cram, recipert, vecipert, rmean1, vmean1);
                                for (ii = 0; ii < 3; ii++)
                                {
                                    Xf[ii, 0] = recipert[ii];
                                    Xf[ii + 3, 0] = vecipert[ii];
                                }*/
                            break;
                        case 'n':  // numerical use previous states to avoid long runtimes
                                   //// the tMath.Since is the same as the previous step, no I don't think so
                                   //tMath.Since = currobsrec.dtmin*60.0;
                                   //astPert::rk4(currobsrec.jd, currobsrec.jdf, tMath.Since, xpert, rmag, derivType, cdam, cram, interp, opt,
                                   //	jpldearr, jdjpldestart, iau80rec, jdeopstart, eoparr, Xf);
                                   //// update for next time step for speed 
                                   //for (ii = 0; ii < 3; ii++)
                                   //{
                                   //	recipert[ii] = Xf[ii, 0];
                                   //	vecipert[ii] = Xf[ii + 3, 0];
                                   //}
                            break;
                        default:  // Kepler 2-body
                            for (ii = 0; ii < 3; ii++)
                            {
                                rpert[ii] = xpert[ii, 0];
                                vpert[ii] = xpert[ii + 3, 0];
                            }
                            // BUG FIX: this used to propagate for currobsrec.dtmin*60.0 seconds - dtmin is
                            // never populated by any current caller (defaults to 0.0), and even when set
                            // it's meant for SGP4's satellite-epoch-relative time convention, not "seconds
                            // since the first observation" (which is what xnom/rnom/vnom represent here).
                            // The NOMINAL branch above computes the correct elapsed time directly from
                            // jd/jdf relative to obsrecarr[1] (dtsec, line ~1150) - the perturbed branch
                            // needs that exact same elapsed time, or it's comparing a correctly-propagated
                            // nominal against a perturbed state frozen at dt=0. That mismatch grows with
                            // how far an observation is from the first one, and hits velocity perturbations
                            // hardest (a velocity change barely affects position at dt=0 by construction),
                            // which is exactly the pattern of huge, badly-conditioned velocity columns in
                            // atwa that we were seeing.
                            double dtsecPert = (currobsrec.jd + currobsrec.jdf - obsrecarr[1].jd - obsrecarr[1].jdf) * 86400.0;  // s
                            AstroLibr.kepler(rpert, vpert, dtsecPert, out recipert, out vecipert, out string keplerOutTextPert);
                            for (ii = 0; ii < 3; ii++)
                            {
                                Xf[ii, 0] = recipert[ii];
                                Xf[ii + 3, 0] = vecipert[ii];
                            }
                            break;
                    }

                    // now get vectors to ITRF for obs procesMath.Sing
                    // eci to itrf if observation type
                    if (currobsrec.obstype != 4)
                    {
                        dtt = (currobsrec.hr * 60 + currobsrec.min + currobsrec.sec / 60.0) * dayconv;
                        EOPSPWLibr.findeopparam(currobsrec.jd, dtt, 's', EOPSPWLibr.eopdata,
                           out dut1, out dat, out lod, out xp, out yp, out ddpsi, out ddeps, out ddx, out ddy);
                        // same jdut1 fix as the nominal branch above - see comment there.
                        jdut1 = currobsrec.jd + currobsrec.jdf + dut1 / 86400.0;
                        tut1 = (jdut1 - 2451545.0) * jdyconv;
                        jdftt = currobsrec.jdf + dat + 32.184;   // sec
                        ttt = (currobsrec.jd + currobsrec.jdf + (dat + 32.184 + jdftt) * dayconv - 2451545.0) * jdyconv;
                        opt = AstroLib.EOpt.e80;
                        // Same fix as the nominal-state call above: convert the perturbed ECI state to
                        // ECEF with Edirection.eto. Also note the *source* here was wrong before - it read
                        // r3/v3 (only ever set on the 's' SGP4 branch above) instead of recipert/vecipert,
                        // which is what the Kepler/two-body branch (the `default:` case a few lines up)
                        // actually populates. For proptype 's' this would have used a stale/zero r3,v3.
                        AstroLibr.eci_ecef(ref recipert, ref vecipert, ref aeci, MathTimeLib.Edirection.eto,
                            ref ritrf, ref vitrf, ref aitrf, iau80arr, ttt, jdut1, lod, xp, yp, ddpsi, ddeps);
                    } // if obstype

                    if (currobsrec.obstype == 3)
                        AstroLibr.rv_tradec(ref ritrf, ref vitrf, rs, MathTimeLib.Edirection.eto, ref rngpert, ref trtascpert, ref tdeclpert,
                            ref drngpert, ref dtrtascpert, ref dtdeclpert);
                    else
                        if (currobsrec.obstype == 4)
                        bstarpert = satrec.bstar * (1.0 + percentchg);
                    else  // currobsrec.obstype = 0 or 1 or 2
                        AstroLibr.rv_razel(ref ritrf, ref vitrf, currsenrec.senlat, currsenrec.senlon, currsenrec.senalt, MathTimeLib.Edirection.eto,
                            ref rngpert, ref azpert, ref elpert, ref drngpert, ref dazpert, ref delpert);
                    switch (currobsrec.obstype)
                    {
                        case 0:
                            a[0, j] = (rngpert - rngnom) / deltaamt;
                            break;
                        case 1:
                            a[0, j] = (azpert - aznom) / deltaamt;
                            a[1, j] = (elpert - elnom) / deltaamt;
                            break;
                        case 2:
                            a[0, j] = (rngpert - rngnom) / deltaamt;
                            a[1, j] = (azpert - aznom) / deltaamt;
                            a[2, j] = (elpert - elnom) / deltaamt;
                            break;
                        case 3:
                            a[0, j] = (trtascpert - trtascnom) / deltaamt;
                            a[1, j] = (tdeclpert - tdeclnom) / deltaamt;
                            break;
                        case 4:
                            a[0, j] = (recipert[0] - recinom[0]) / deltaamt;
                            a[1, j] = (recipert[1] - recinom[1]) / deltaamt;
                            a[2, j] = (recipert[2] - recinom[2]) / deltaamt;
                            a[3, j] = (vecipert[0] - vecinom[0]) / deltaamt;
                            a[4, j] = (vecipert[1] - vecinom[1]) / deltaamt;
                            a[5, j] = (vecipert[2] - vecinom[2]) / deltaamt;
                            // same fix as the b-matrix/dbstar2 bounds issue - a is sized [indobs,statesize]
                            // and row 6 doesn't exist for a 6-state (no bstar) fit.
                            if (statesize == 7)
                                a[6, j] = (bstarpert - bstarnom) / deltaamt;
                            break;
                    } // case

                    //printf( "rpert %11.5f %11.5f %11.5f %8.3f %8.3f %8.3f \n",
                    //              r3[0], r3[1], r3[2], rngpert, azpert*rad, elpert*rad );  // ritrf
                    //printf("recipert %11.5f %11.5f %11.5f \n", recipert[0], recipert[1], recipert[2]);  // ritrf
                    //printf("xpert %11.5f %11.5f %11.5f \n", xpert[0, 0], xpert[1, 0], xpert[2, 0]);  // ritrf

                    // ----------------- reset the modified vector --------------------
                    satrecp = satrec;

                }  // for j = 0 to statesize

                // reset state for numerical cases to avoid starting from scratch each time
                if (proptype == 'n')
                {
                    for (ii = 0; ii < 3; ii++)
                    {
                        XCurrn[ii, 0] = recinom[ii];
                        XCurrn[ii + 3, 0] = vecinom[ii];
                    }
                }

                // ----------------- now form the matrix combinations -----------------
                at = MathTimeLibr.mattrans(a, statesize);  // indobs, 

                // ------------------------- assign weights ---------------------------
                switch (currobsrec.obstype)
                {
                    case 0:
                        {
                            w1 = 1.0 / (currsenrec.noiserng * currsenrec.noiserng);
                            rng2 = rng2 + b[0, 0] * b[0, 0] * w1;
                        }
                        break;
                    case 1:
                        {
                            w1 = 1.0 / (currsenrec.noiseaz * currsenrec.noiseaz);
                            w2 = 1.0 / (currsenrec.noiseel * currsenrec.noiseel);
                            daz2 = daz2 + b[0, 0] * b[0, 0] * w1;
                            del2 = del2 + b[1, 0] * b[1, 0] * w2;
                        }
                        break;
                    case 2:
                        {
                            w1 = 1.0 / (currsenrec.noiserng * currsenrec.noiserng);
                            w2 = 1.0 / (currsenrec.noiseaz * currsenrec.noiseaz);
                            w3 = 1.0 / (currsenrec.noiseel * currsenrec.noiseel);
                            drng2 = drng2 + b[0, 0] * b[0, 0] * w1;
                            daz2 = daz2 + b[1, 0] * b[1, 0] * w2;
                            del2 = del2 + b[2, 0] * b[2, 0] * w3;
                        }
                        break;
                    case 3:
                        {
                            w1 = 1.0 / (currsenrec.noisetrtasc * currsenrec.noisetrtasc);
                            w2 = 1.0 / (currsenrec.noisetdecl * currsenrec.noisetdecl);
                            dtrtasc2 = dtrtasc2 + b[0, 0] * b[0, 0] * w1;
                            dtdecl2 = dtdecl2 + b[1, 0] * b[1, 0] * w2;
                        }
                        break;
                    case 4:
                        {
                            w1 = 1.0 / (currsenrec.noisex * currsenrec.noisex);
                            w2 = 1.0 / (currsenrec.noisey * currsenrec.noisey);
                            w3 = 1.0 / (currsenrec.noisez * currsenrec.noisez);
                            w4 = 1.0 / (currsenrec.noisexdot * currsenrec.noisexdot);
                            w5 = 1.0 / (currsenrec.noiseydot * currsenrec.noiseydot);
                            w6 = 1.0 / (currsenrec.noisezdot * currsenrec.noisezdot);
                            w7 = 1.0 / (currsenrec.noisebstar * currsenrec.noisebstar);
                            dx2 = dx2 + b[0, 0] * b[0, 0] * w1;
                            dy2 = dy2 + b[1, 0] * b[1, 0] * w2;
                            dz2 = dz2 + b[2, 0] * b[2, 0] * w3;
                            dxdot2 = dxdot2 + b[3, 0] * b[3, 0] * w4;
                            dydot2 = dydot2 + b[4, 0] * b[4, 0] * w5;
                            dzdot2 = dzdot2 + b[5, 0] * b[5, 0] * w6;
                            // same fix as the b-matrix construction above - b[6,0] doesn't exist for a
                            // 6-state (no bstar) fit.
                            if (statesize == 7)
                                dbstar2 = dbstar2 + b[6, 0] * b[6, 0] * w7;
                        }
                        break;
                } // case

                for (rowc = 0; rowc < statesize; rowc++)
                {
                    for (colc = 0; colc < indobs; colc++)
                    {
                        switch (colc)
                        {
                            case 0:
                                weight = w1;
                                break;
                            case 1:
                                weight = w2;
                                break;
                            case 2:
                                weight = w3;
                                break;
                            case 3:
                                weight = w4;
                                break;
                            case 4:
                                weight = w5;
                                break;
                            case 5:
                                weight = w6;
                                break;
                            case 6:
                                weight = w7;
                                break;
                            default:
                                weight = 1.0;
                                break;
                        }  // case
                        atw[rowc, colc] = at[rowc, colc] * weight;
                    }  // for colc
                } // for rowc

                // ----------------- find the atwa / atwb matrices --------------------
                atwaacc = MathTimeLibr.matmult(atw, a, statesize, indobs, statesize);
                atwbacc = MathTimeLibr.matmult(atw, b, statesize, indobs, 1);

                // ------------------- accumulate the matricies -----------------------
                for (r = 0; r < statesize; r++)
                    for (c = 0; c < statesize; c++)
                        atwa[r, c] = atwaacc[r, c] + atwa[r, c];

                c = 0;
                for (r = 0; r < statesize; r++)
                    atwb[r, c] = atwbacc[r, c] + atwb[r, c];

            } // for i through the observations
            string outstr;
            MathTimeLibr.writeexpmat("atwa ", atwa, statesize, statesize, out outstr);
            MathTimeLibr.writeexpmat("atwb ", atwb, statesize, 1, out outstr);
            MathTimeLibr.writeexpmat("a ", a, indobs, statesize, out outstr);
            MathTimeLibr.writemat("b ", b, indobs, 1, out outstr);
        }  // findatwaatwb


        /* -----------------------------------------------------------------------------
        *
        *                           procedure finitediff
        *
        * this procedure perturbs the components of the state vector for processing
        * with the finite differencing for the a matrix.
        *
        *  author        : david vallado                  719-573-2600   15 jan 2008
        *
        *  inputs          description                    range / units
        *    whichconst  - parameter for sgp4 constants   wgs72, wgs721, wgs84
        *    pertelem    - which element to perturb
        *    percentchg  - amount to modify the vectors   0.001
        *                  by in finite differencing
        *    deltaamtchg - tolerance for small value in
        *                  finite differencing            0.0000001
        *    statetype   - type of elements (equinoctial, etc)  'e', 't'
        *    xnom        - state vector                   varied
        *
        *  outputs       :
        *    deltaamt    - amount each elemnt is perturbed
        *    satrec      - satellite record
        *
        *  locals        :
        *    jj          - index
        *
        *  coupling      :
        *    state2satrec- conversion between state and satellite structure
        *    sgp4init    - intiialize sgp4 constants
        *
        *  references    :
        *    vallado       2007, 753-765
        * --------------------------------------------------------------------------- */

        public void finitediff
            (
            SGP4Lib.gravconsttype whichconst,
            int pertelem, double percentchg, double deltaamtchg, char statetype, int statesize, char proptype,
            ref SGP4Lib.elsetrec satrec,
            AstroLib.EOpt opt,
            double[] xnom,
            double jd, double jdf,
            out double[,] xpert,
            out double deltaamt
            )
        {
            double[] xtemp = new double[statesize];
            xpert = new double[7, 1];
            int jj;
            double newpct;
            newpct = 1.0;

            for (jj = 0; jj < 6; jj++)
            {
                xtemp[jj] = xpert[jj, 0];
            }

            // chk if perturbing amt is too small. if so, up the percentchg and try again
            // this will execute 5 times, but leaves percentchg the same after each run
            jj = 1;
            do
            {
                deltaamt = xnom[pertelem] * percentchg;
                xtemp[pertelem] = xnom[pertelem] + deltaamt;

                // need entire state here...
                if (proptype == 's')
                    state2satrec(ref xtemp, jd, jdf, statetype, statesize, opt, MathTimeLib.Edirection.eto, ref satrec);

                if (Math.Abs(deltaamt) < deltaamtchg)       // 0.00001
                {
                    newpct = newpct * percentchg * 1.4;  // increase by 40% and try again
                                                         // printf(" %2i percentchg chgd %11.8f ",jj,newpct);
                }
                jj = jj + 1;
            } while ((Math.Abs(deltaamt) < deltaamtchg) && (jj < 5));

            // printf(" \n");
            // ---- obtain various parameters ----
            //      satrec.a    = pow( satrec.no * tumin , (-2.0/3.0) );
            //      satrec.alta = satrec.a * (1.0 + satrec.ecco) - 1.0;
            //      satrec.altp = satrec.a * (1.0 - satrec.ecco) - 1.0;

            // ---------------- initialize the orbit at sgp4epoch --------------------
            if (proptype == 's')
                SGP4Libr.sgp4init(whichconst, satrec.operationmode, satrec.satnumStr, satrec.jdsatepoch - 2433281.5, satrec.bstar,
                    satrec.ndot, satrec.nddot, satrec.ecco, satrec.argpo, satrec.inclo, satrec.mo, satrec.no_kozai,
                    satrec.nodeo, ref satrec);

            // assign perturbed vector
            for (jj = 0; jj < 6; jj++)
            {
                if (jj == pertelem)
                    xpert[jj, 0] = xtemp[jj];
                else
                    xpert[jj, 0] = xnom[jj];
            }

        }  // finitediff


        /* -----------------------------------------------------------------------------
        *
        *                           procedure leastsquares
        *
        * this procedure performs orbit determination using least squares differential
        *   correction method. a variety of observation types are possible.
        *
        *  algorithm     : find the atwa and atwb matrices by accumulation
        *                  calculate the matrix combinations
        *
        *  author        : david vallado                  719-573-2600    6 aug 2008
        *
        *  inputs          description                    range / units
        *    percentchg  - amount to modify the vectors   0.001
        *                  by in finite differencing
        *    deltaamtchg - tolerance for small value in
        *                  finite differencing            0.0000001
        *    epsilon     - tolerance for calculations     0.0002
        *    whichconst  - parameter for sgp4 for constants wgs72, wgs721, wgs84
        *    typeans     - type of dc (bksub or svd)      'b','s'
        *    firstob     - first obs record to use        0
        *    lastob      - last obs record to use         100
        *    statesize   - size of state                  6 , 7
        *    obsrecarr   - 10000 array structure containing:
        *                  senum, jd, rsvec, obstype,
        *                  rng, az, el, drng, daz, del,
        *                  trtasc, tdecl data
        *    loops       - number of iterations
        *    satrec      - record of satellite parameters for TLE
        *                  this exists for test purposes only
        *    statetype   - type of elements (equinoctial, etc)  'e', 't'
        *    outfile, outfile1 - output files
        *    derivType = 1 = 2-body
        *    derivType = 2 = J2
        *    derivType = 3 = J3
        *    derivType = 4 = J4
        *    derivType = 5 = Jx
        *    derivType = 6 = drag
        *    derivType = 7 = 3-body
        *    derivType = 8 = srp
        *    derivType = 10 = all
        *    derivtype = 11 - srp semianalytical
        *    derivtype = 12 - 3body semianalytical
        *    derivtype = 13 - zonal semianalytical
        *    interp      - interpolation for sun/moon vectors 'l', 's'
        *    opt         - method of calc sun/moon        'a', 'j'
        *
        *  outputs       :
        *    xnom        - nominal state vector
        *    x           - state vector
        *    dx          - state correction
        *    atwai       - covariance matrix
        *    atwa        - atwa matrix - needed for svd processing
        *    atwb        - atwab matrix
        *    satrec      - record of satellite parameters for TLE
        *
        *  locals        :
        *    numwork     - number of observations working with
        *    sigmaold    -
        *    sigmanew    -
        *    drng2       -
        *    daz2        -
        *    del2        -
        *    dtrtasc2    -
        *    dtdecl2     -
        *
        *  coupling      :
        *    findatwaatwb- find combination matrices for use in final equations
        *    matinverse  - matrix inverse routine
        *    matmult     - matrix multiply routine
        *
        *  references    :
        *    vallado       2007, 753-765
        * --------------------------------------------------------------------------- */

        public void leastsquares
            (
                double percentchg, double deltaamtchg, double epsilon, SGP4Lib.gravconsttype whichconst,
                char interp, double jdeopstart, ref SGP4Lib.elsetrec satrec,
                char typeans, char statetype, char proptype, int firstob, int lastob, int statesize,
                ref obsrec[] obsrecarr, out int loops,
                ref double[] xnom, double jd, double jdf, int derivType, double cdam, double cram, char opt,

                EOPSPWLib.iau80Class iau80arr, EOPSPWLib.iau06Class iau06arr, ref double[] eoparr,
                out double[,] x, out double[,] dx, out double[,] atwai, out double[,] atwa, out double[,] atwb
            )
        {
            int i, numwork, ok;  // j,
            double sigmanew, sigmaold, sigmaold2, drng2, daz2, del2, ddrng2, ddaz2, ddel2, dtrtasc2,
                dtdecl2, dx2, dy2, dz2, dxdot2, dydot2, dzdot2, wtest, wmax;
            SGP4Lib.elsetrec satrecorig;  // for restarting cases
            char restart;
            // NOTE: xnom is now `ref` - it carries the caller's initial nominal state guess in,
            // and the converged solution out. (Previously this was `out` and got zeroed out here,
            // which meant leastsquares always started the DC from a zero state vector - a bug in
            // the original conversion, since nothing currently calls this method.)
            x = new double[7, 1];  // rows
            dx = new double[7, 1];  // rows
            atwai = new double[7, 7];  // rows
            atwa = new double[7, 7];  // rows
            atwb = new double[7, 1];  // rows
            obsrec currobsrec;

            // should state size be 7 when just doing xyz processing, and no bstar?

            double rad, rmag;
            int indobs = 7;
            double[,] tmmp = new double[statesize, statesize];
            double[,] b = new double[indobs, 1];
            double[,] atw = new double[statesize, indobs];
            double[,] u = new double[statesize, indobs];
            double[,] w = new double[indobs, 1];
            double[,] wutatwb = new double[indobs, 1];
            double[,] ww = new double[indobs, indobs];
            double[,] v = new double[indobs, indobs];
            double[,] ut = new double[indobs, statesize];
            double[,] vw = new double[indobs, indobs];
            double[,] vwut = new double[indobs, 1];
            double[,] atwao = new double[statesize, statesize];

            double[] limitdx = new double[indobs];
            StringBuilder strBuild = new StringBuilder();

            // ---------------------- initialize parameters ---------------------------
            // ---- Initialize EOP coordinate params
            string nutLoc;
            nutLoc = @"D:\Codes\LIBRARY\DataLib\nut80.dat";
            EOPSPWLibr.iau80in(nutLoc, out iau80arr);
            //nutLoc = @"D:\Codes\LIBRARY\DataLib\";
            //EOPSPWLibr.iau06in(nutLoc, out iau06arr);

            // ---- Initialize EOP data
            string EOPupdate;
            string eopFileName = @"D:\Codes\LIBRARY\DataLib\EOP-All-v1.1_2020-02-12.txt";
            Int32 ktrActObs, mjdeopstart;
            // readeop only takes (ref eopdata, inFile, out ktrActualObs, out updDate) - it doesn't
            // return a start-day value at all (findeopparam works out its own start day internally
            // from eopdata[0].mjd). mjdeopstart is still threaded through to findatwaatwb below as a
            // formal parameter, but findatwaatwb no longer actually reads it - leaving it at 0 is fine.
            EOPSPWLibr.readeop(ref EOPSPWLibr.eopdata, eopFileName, out ktrActObs, out EOPupdate);
            mjdeopstart = 0;

            // ---- Initialize jplde data
            // readjplde only takes (ref jpldearr, jplLoc, infilename) - no start-date outputs; each
            // record already carries its own .mjd, and jdjpldestart below is only ever consumed by
            // the commented-out semianalytical/numerical branches in findatwaatwb, so 0.0 is fine.
            AstroLib.jpldedataClass[] jpldearr = AstroLibr.jpldearr;
            AstroLibr.readjplde(ref jpldearr, "D:/Codes/LIBRARY/DataLib/", "sunmooneph_430t12.txt");
            double jdjpldestart = 0.0;


            restart = 'n';
            satrecorig = satrec;
            rad = 180.0 / Math.PI;

            sigmaold = 20000.0;
            sigmaold2 = 20000.0;
            sigmanew = 10000.0;

            numwork = lastob - firstob;   // old +2

            limitdx[0] = 0.80 * xnom[0];
            limitdx[1] = 0.80 * xnom[1];
            limitdx[2] = 0.80 * xnom[2];
            limitdx[3] = 0.40 * xnom[3];
            limitdx[4] = 0.40 * xnom[4];
            limitdx[5] = 0.40 * xnom[5];
            if (statesize == 7)
                limitdx[6] = 0.10 * xnom[6];

            //fprintf(outfile, "nom coe %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f \n",
            //    xnom[0], xnom[1], xnom[2], xnom[3], xnom[4], xnom[5]);
            strbuild.AppendLine("nom coe " + xnom[0] + " " + xnom[1] + " " + xnom[2] + " " +
                xnom[3] + " " + xnom[4] + " " + xnom[5] + @"\n");
            rmag = Math.Sqrt(xnom[0] * xnom[0] + xnom[1] * xnom[1] + xnom[2] * xnom[2]);



            // -------------------------- loop through the iterations --------------------------- 
            // several conditions apply for a stopMath.PIng condition:
            // changes become small,
            // max iterations met,
            // overall change is negligible,
            // iterations are diverging

            // while through all the iterations
            // NOTE: bumped from 5 to 25. The vector-state (statetype=='v') correction step below caps
            // each component's per-iteration change at 1% of its current value, unlike Ex10_4.m's
            // reference implementation which takes a full, unclamped Newton step every iteration - so
            // this needs more iterations than MATLAB's fixed 5 to cover the same total correction from
            // a badly-off seed. Once things are actually converging (post the dtsec/dtmin timing fix),
            // this is a "needs more steps" issue, not a "stuck" one - so raising the cap is the safe
            // fix rather than touching the clamp itself, which may be relied on elsewhere for stability.
            loops = 0;
            while ((Math.Abs((sigmanew - sigmaold) / sigmaold) >= epsilon && sigmanew >= epsilon) &&
                !(sigmanew > sigmaold && sigmaold > sigmaold2 && sigmanew > 500000.0) && loops < 25)
            {
                sigmaold2 = sigmaold;
                sigmaold = sigmanew;
                drng2 = 0.0;
                daz2 = 0.0;
                del2 = 0.0;
                ddrng2 = 0.0;
                ddaz2 = 0.0;
                ddel2 = 0.0;
                dtrtasc2 = 0.0;
                dtdecl2 = 0.0;
                dx2 = 0.0;
                dy2 = 0.0;
                dz2 = 0.0;
                dxdot2 = 0.0;
                dydot2 = 0.0;
                dzdot2 = 0.0;

                // ---- place current nominal value into the state vector ----

                //    SGP4DCLib.obsrec[] obsrecarr,  //obsrec
                //double[] eoparr,   // <eopdata>
                //Int32 mjdeopstart,
                //double[] jpldearr,  // <jpldedata>
                //double jdjpldestart,

                findatwaatwb(firstob, lastob, statesize, percentchg, deltaamtchg,
                    whichconst, interp, jdeopstart, jd, jdf, ref satrec, obsrecarr, statetype, proptype, xnom, rmag,
                    derivType, cdam, cram, AstroLib.EOpt.e80, iau80arr, iau06arr, eoparr, mjdeopstart, jpldearr, jdjpldestart,
                    out drng2, out daz2, out del2, out ddrng2, out ddaz2, out ddel2, out dtrtasc2, out dtdecl2,
                    out dx2, out dy2, out dz2, out dxdot2, out dydot2, out dzdot2, out atwa, out atwb, out atw, out b);

                // ---- Ex10_4.m-style diagnostic dump, so a run's output can be diffed line-for-line
                // against the MATLAB reference. MATLAB only prints atwa/atwb on the first pass (j==1) -
                // capture them here, before either inversion path below has a chance to mutate atwa
                // (dsvdcmp does, in place, for the SVD path).
                if (loops == 0)
                {
                    strbuild.AppendLine("atwa = ");
                    for (int rr = 0; rr < statesize; rr++)
                    {
                        StringBuilder rowStr = new StringBuilder();
                        for (int cc = 0; cc < statesize; cc++)
                            rowStr.Append(string.Format(CultureInfo.InvariantCulture, "{0,10:F1} ", atwa[rr, cc]));
                        strbuild.AppendLine(rowStr.ToString().TrimEnd());
                    }
                    strbuild.AppendLine("atwb = ");
                    for (int rr = 0; rr < statesize; rr++)
                        strbuild.AppendLine(string.Format(CultureInfo.InvariantCulture, "{0,10:F1} ", atwb[rr, 0]));
                }

                // matrix inversion approach
                if (typeans == 'b')
                {
                    MathTimeLibr.matinverse(atwa, statesize, out atwai);
                    dx = MathTimeLibr.matmult(atwai, atwb, statesize, statesize, 1);
                    // writeexpmat("atwai1",atwai,statesize,statesize);
                    // matmult   ( atwai,atwa,tmmp, statesize,statesize,statesize );
                    // writeexpmat("tmmp ",tmmp,statesize,statesize);
                }

                // svd approach for matrix inversion
                if (typeans == 's')
                {
                    atwao = (double[,])atwa.Clone(); // save an actual independent copy for later -
                                                      // dsvdcmp now mutates atwa in place (ref, not out),
                                                      // so a plain reference assignment here would have
                                                      // let atwao silently track the mutated U matrix too
                    ok = dsvdcmp(ref atwa, statesize, statesize, out w, out v);
                    wmax = 0.0;
                    for (i = 0; i < statesize; i++)
                    {
                        if (w[i, 0] > wmax)
                            wmax = w[i, 0];
                    }
                    wtest = 1.0e-14 * wmax;  // -10 seems to be a limit
                    for (i = 0; i < statesize; i++)
                    {
                        if (w[i, 0] <= wtest)
                            w[i, 0] = 1.0e-10; // set to 0.0 with dsvbksb 1e-10 else, delete bad term
                        ww[i, i] = 1.0 / w[i, 0];
                    }

                    // dsvbksb(atwa, w, v, statesize, indobs, b, dx);
                    vw = MathTimeLibr.matmult(v, ww, statesize, statesize, statesize);
                    ut = MathTimeLibr.mattrans(atwa, statesize);
                    atwai = MathTimeLibr.matmult(vw, ut, statesize, statesize, statesize);
                    dx = MathTimeLibr.matmult(atwai, atwb, statesize, statesize, 1);
                    // writeexpmat("atwai2",atwai,statesize,statesize);
                    // matmult   ( atwai,atwao,tmmp, statesize,statesize,statesize );
                    // writeexpmat("tmmp ",tmmp,statesize,statesize);
                }

                // ---- guard against a singular/near-singular atwa producing NaN/Infinity in dx ----
                // IEEE-754 NaN comparisons are always false, so the "> 0.01" style clamps just below
                // silently let a NaN straight through into xnom instead of catching it - stop here
                // instead of corrupting the state and continuing to iterate on garbage.
                bool dxIsBad = false;
                for (i = 0; i < statesize; i++)
                    if (double.IsNaN(dx[i, 0]) || double.IsInfinity(dx[i, 0]))
                        dxIsBad = true;
                if (dxIsBad)
                {
                    strbuild.AppendLine("leastsquares: atwa was singular/near-singular (dx came back " +
                        "NaN/Infinity) at loop " + loops + " - stopping and keeping the last good xnom " +
                        "rather than corrupting it.");
                    break;
                }

                // ---- update state vector ----
                for (i = 0; i < statesize; i++)
                {
                    if (statetype != 'v')
                    {
                        // update elements and limit
                        if ((loops > -1) && (Math.Abs(dx[i, 0] / xnom[i]) > 1000.0))   // 100
                        {
                            dx[i, 0] = 0.10 * xnom[i] * Math.Sign(dx[i, 0]);   // 0.30 try leaving the same
                        }
                        else
                        if ((loops > 0) && (Math.Abs(dx[i, 0] / xnom[i]) > 200.0))   // 100
                        {
                            dx[i, 0] = 0.30 * xnom[i] * Math.Sign(dx[i, 0]);   // 0.30 try leaving the same
                        }
                        else
                        {
                            if ((loops > 0) && (Math.Abs(dx[i, 0] / xnom[i]) > 100.0))  // 20
                            {
                                dx[i, 0] = 0.70 * xnom[i] * Math.Sign(dx[i, 0]);   // 0.70 - 0.80 about same
                            }
                            else
                            if ((loops > 0) && (Math.Abs(dx[i, 0] / xnom[i]) > 10.0))      // 5
                            {
                                dx[i, 0] = 0.90 * xnom[i] * Math.Sign(dx[i, 0]);   // 0.90 try leaving the same
                            }
                        }
                    }
                    else
                    {
                        // Vector-state correction clamp intentionally left disabled (not just a temporary
                        // test): with the real bugs fixed (dsvdcmp, finitediff, the dtmin/dtsec timing
                        // mismatch, and the jdut1 Julian-date bug), the unclamped Newton step converges
                        // in ~7 loops and matches Ex10_4.m's MATLAB reference to sub-km position and
                        // ~0.001 km/s velocity - the clamp was never guarding against a real divergence
                        // here, just slowing/preventing convergence. Left commented out below rather than
                        // deleted - if a future test case (messier data, worse initial guess, etc.) shows
                        // a real need for damping, a much looser bound than the original 1% would be the
                        // starting point, not this exact line.
                        // if (Math.Abs(dx[i, 0] / xnom[i]) > 0.01)
                        //     dx[i, 0] = 0.01 * xnom[i] * Math.Sign(dx[i, 0]);
                    }  // limit corrections if vectors

                    xnom[i] = xnom[i] + dx[i, 0];
                }  // for i

                // -------------- re-initialize state with new values -----------------
                if (proptype == 's')
                {
                    state2satrec(ref xnom, jd, jdf, statetype, statesize, AstroLib.EOpt.e80, MathTimeLib.Edirection.eto, ref satrec);
                    SGP4Libr.sgp4init(whichconst, satrec.operationmode, satrec.satnumStr, satrec.jdsatepoch - 2433281.5, satrec.bstar,
                        satrec.ndot, satrec.nddot, satrec.ecco, satrec.argpo, satrec.inclo, satrec.mo, satrec.no_kozai,
                        satrec.nodeo, ref satrec);
                }

                currobsrec = obsrecarr[0];
                switch (currobsrec.obstype)
                {
                    case 0:
                        sigmanew = Math.Sqrt(drng2 / numwork);
                        break;
                    case 1:
                        sigmanew = Math.Sqrt((daz2 + del2) / numwork);
                        break;
                    case 2:
                        sigmanew = Math.Sqrt((drng2 + daz2 + del2) / numwork);
                        break;
                    case 3:
                        sigmanew = Math.Sqrt((dtrtasc2 + dtdecl2) / numwork);
                        break;
                    case 4:
                        sigmanew = Math.Sqrt((dx2 + dy2 + dz2 + dxdot2 + dydot2 + dzdot2) / numwork);
                        break;
                } // case

                // ---- Ex10_4.m-style per-iteration dump (dx / output r1,v1 / rms), every pass ----
                StringBuilder dxLine = new StringBuilder("dx ");
                for (int ii = 0; ii < statesize; ii++)
                    dxLine.Append(string.Format(CultureInfo.InvariantCulture, "{0,16:F8} ", dx[ii, 0]));
                strbuild.AppendLine(dxLine.ToString());
                strbuild.AppendLine("output: ");
                strbuild.AppendLine(string.Format(CultureInfo.InvariantCulture,
                    "r1 {0,16:F8} {1,16:F8} {2,16:F8} km v1 {3,16:F8} {4,16:F8} {5,16:F8} km ",
                    xnom[0], xnom[1], xnom[2], xnom[3], xnom[4], xnom[5]));
                strbuild.AppendLine(string.Format(CultureInfo.InvariantCulture, "rms  {0,16:F8}  ", sigmanew));

                // ----------------------- write out data -----------------------------
                //fprintf(outfile, "corrections  %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f ",
                //    dx[0, 0], dx[1, 0], dx[2, 0], dx[3, 0], dx[4, 0], dx[5, 0]);
                //fprintf(outfile, "loop %2i  nominal %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f ",
                //    loops, xnom[0], xnom[1], xnom[2], xnom[3], xnom[4], xnom[5]);
                //fprintf(outfile, "sigold %12.6f %1c  %12.6f %12.6f %12.6f \n", sigmanew,
                //    limited, sigmaold2, sigmaold, sigmanew);

                //printf("corrections  %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f \n",
                //    dx[0, 0], dx[1, 0], dx[2, 0], dx[3, 0], dx[4, 0], dx[5, 0]);
                //if (proptype == 's')
                //    SGP4DCLibr.printtle(satrec);
                // printf( "sigmaold %12.6f  sigmanew %12.6f \n",sigmaold,sigmanew );
                // fprintf( outfile," ---------------------------------------------------- \n" );
                loops = loops + 1;

                // try not solving for bstar if the iterations are diverging rapidly
                if ((sigmanew > sigmaold) && (sigmaold > sigmaold2) && (sigmanew > 500000.0) && (restart == 'n'))   // 1000
                {
                    restart = 'y';
                    sigmaold = 20000.0;
                    sigmanew = 10000.0;
                    loops = 0;
                    //printf("restarted\n");
                    //fprintf(outfile, "restarted\n");
                    //fprintf(outfile1, "restarted\n");
                    satrec = satrecorig;
                    state2satrec(ref xnom, jd, jdf, statetype, statesize, AstroLib.EOpt.e80, MathTimeLib.Edirection.efrom, ref satrec);
                }

            }  // do through iterations


            if (proptype == 's')
            {
                strBuild.AppendLine("loop " + loops + " nom coe " + satrec.a * 6378.135 + " "
                    + satrec.no_kozai * 1440.0 / (2.0 * Math.PI) + " "
                    + satrec.ecco + " " + (satrec.inclo * rad) + " " + satrec.nodeo * rad + " "
                    + (satrec.argpo * rad) + " " + (satrec.mo * rad) + " " + satrec.bstar + @"\n");
                strBuild.AppendLine("loop " + loops + " nom coe " + satrec.no_kozai * 1440.0 / (2.0 * Math.PI) + " "
                    + satrec.ecco + " " + (satrec.inclo * rad) + " " + satrec.nodeo * rad + " "
                    + (satrec.argpo * rad) + " " + (satrec.mo * rad) + " " + satrec.bstar + @"\n");
                strBuild.AppendLine("3 pos component uncertainty \n");
            }
            else
            {
                strBuild.AppendLine("loop " + loops + " nom coe " + xnom[0] + " " + xnom[1] + " " + xnom[2]
                    + " " + xnom[3] + " " + xnom[4] + " " + xnom[5] + @"\n");
            }

            // printtle(satrec);

        } // leastsquares

    }
}
