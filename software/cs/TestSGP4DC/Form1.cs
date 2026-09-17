/* ---------------------------------------------------------------------
*
*                              testSGP4DC.cs
*
*  this program tests the sgp4 differential correction. 
*
*                          companion code for
*             fundamentals of astrodynamics and applications
*                                  2022
*                            by david vallado
*
*     (w) 719-573-2600, email dvallado@agi.com, davallado@gmail.com
*     *****************************************************************
*  current :
*             3 nov 14  david vallado
*                        conversion to msvs c++
*  changes :
*             8 jul 08  david vallado
*                        original version for 2008 paper
*
*   typeans    = *_argv[1];
*   typeodrun  = *_argv[2];
*   statetype  = *_argv[3];
*   percentchg = atof(_argv[4]);    // 0.001
*   deltaamtchg= atof(_argv[5]);    // 0.0000001
*   rmsepsilon = atof(_argv[6]);    // 0.001
*   strcpy( fname, _argv[7]);
*   lastob     = atof(_argv[8]);    // 200
*   statesize  = atof(_argv[9]);    // 6 or 7
*
*       ----------------------------------------------------------------      */


using System;
using System.Collections.Generic;
using System.Data;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Windows.Forms;

using MathTimeMethods;  // Edirection, globals
using EOPSPWMethods;    // EOPDataClass, SPWDataClass, iau80Class, iau06Class
using AstroLibMethods;     // EOpt, gravityConst, astroConst, xysdataClass, jpldedataClass
using SGP4Methods;
using SGP4DCMethods;



namespace TestSGP4DC
{
    public partial class Form1 : Form
    {
        // setup the classes so methods can be called
        public MathTimeLib MathTimeLibr = new MathTimeLib();

        public EOPSPWLib EOPSPWLibr = new EOPSPWLib();

        public AstroLib AstroLibr = new AstroLib();

        public SGP4Lib SGP4Libr = new SGP4Lib();

        public SGP4DCLib SGP4DCLibr = new SGP4DCLib();

        public StringBuilder strbuild = new StringBuilder();


        public Form1()
        {
            InitializeComponent();
            
            if (this.cbTypeODRun.SelectedItem == null)
                this.cbTypeODRun.SelectedIndex = 1;
            if (this.cbmatinv.SelectedItem == null)
                this.cbmatinv.SelectedIndex = 0;
            if (this.cbproptype.SelectedItem == null)
                this.cbproptype.SelectedIndex = 2;
        }

        /*       ----------------------------------------------------------------     */
        // ---- Combo box option -> internal code mappings ----
        // cbTypeODRun: "Range - az - el" | "TLE" | "Ephemeris"   (kept as the literal string)
        // cbmatinv:    "SVD" -> 's'   (SVD-based inverse, dsvdcmp path in SGP4DCLib.leastsquares)
        //              "LuBkSub" -> 'b' (direct matinverse/back-substitution path)
        // cbproptype:  "Two-body" -> 't' (falls through to the default/Kepler case in findatwaatwb)
        //              "SGP4" -> 's', "Numerical" -> 'n', "Semianalytical" -> 'a'
        private char MapMatInvChoice(string choice)
        {
            switch (choice)
            {
                case "SVD": return 's';
                case "LuBkSub": return 'b';
                default: return 'b';
            }
        }

        private char MapPropTypeChoice(string choice)
        {
            switch (choice)
            {
                case "SGP4": return 's';
                case "Numerical": return 'n';
                case "Semianalytical": return 'a';
                case "Two-body":
                default: return 't';  // anything other than 's','a','n' takes the Kepler/2-body path
            }
        }

        /*       ----------------------------------------------------------------     */
        public void button4_Click(object sender, EventArgs e)
        {
            string typeodrun = "Range - az - el";
            if (this.cbTypeODRun.SelectedItem != null)
                typeodrun = this.cbTypeODRun.SelectedItem.ToString();

            string matinvChoice = null;
            if (this.cbmatinv.SelectedItem != null)
                matinvChoice = this.cbmatinv.SelectedItem.ToString();
            char typeans = MapMatInvChoice(matinvChoice);

            string proptypeChoice = null;
            if (this.cbproptype.SelectedItem != null)
                proptypeChoice = this.cbproptype.SelectedItem.ToString();
            char proptype = MapPropTypeChoice(proptypeChoice);

            int obsskip = 1;
            int obsread = 20;
            int.TryParse(this.obs2skip.Text, out obsskip);
            int.TryParse(this.obs2read.Text, out obsread);
            if (obsskip < 1) obsskip = 1;
            if (obsread < 1) obsread = 20;

            string outputFolder = this.OFLoc.Text;

            switch (typeodrun)
            {
                // ---------------------- obs range az and el (matches ex10_4.m / findatwaatwb.m) ---------------------
                case "Range - az - el":
                    // book example nominal (deliberately off, per Vallado ex10_4.m, to show DC convergence)
                    double[] xnomSeed = new double[6]
                    {
                        5975.29040000,  2568.64000000,  3120.58450000,   // reci  km
                        3.98384600,    -2.07115900,    -5.91709500       // veci  km/s
                    };

                    processrazel(typeans, proptype, obsskip, obsread, xnomSeed, outputFolder);
                    break;

                // ---------------------------- ephemeris -------------------------
                case "Ephemeris":
                    processxyz(typeans, proptype, obsskip, obsread, outputFolder);
                    break;

                // ---------------------------- tle, form ephemeris, then try to refind original ----------------------------
                case "TLE":
                    processtle(typeans, outputFolder, obsread);
                    break;
            }
        }




        // ---- UNCONVERTED LEGACY CODE (readobs) commented out below pending C# port ----
//         /*       ----------------------------------------------------------------     */
//         // read obs from an ephemeris file
// 
//         void readobs
//         (
//             string filename,
//             iau80data iau80arr,
//             iau06data iau06arr,
//      double[] eoparr,  //<eopdata>
//     double jdeopstart, double jdeopstartf,
//     char interp,
//     char proptype,
//     out double jdepoch, out double jdepochf,
//     out double[] obsarr,  //<obsrec>
//     int firstob, int obsskip, int obsread,
//     out double[] xnom = new double[7]
// )
//         {
//             //FILE* infile;
//             int j, ii, kk;
//             int year, mon, day, hr, min, timezone, dat;
//             string longstr1, longstr, coords;
// 
//             string strr, monstr;
//             char strr1, strr2, strr3;
//             string tmp, tmpc, tmpc1, units;
//             double dtsec, dtseczero, conv, sec, mfme;
//             double dut1, lod, xp, yp, ddpsi, ddeps, ddx, ddy, icrsx, y, s, jdxysstart,
//                 ut1, tut1, jdut1, jdut1f, utc, tai, tt, ttt, jdtt, jdftt, tcg, tdb, ttdb, jdtdb, jdtdbf, tcb,
//                 p, a, ecc, incl, raan, argp, nu, m, eccanom, arglat, truelon, lonper;
//             double[] rpos = new double[3];
//             double[] vpos = new double[3];
//             double[] apos = new double[3];
//             double[] reci = new double[3];
//             double[] veci = new double[3];
//             double[] aeci = new double[3];
//             double dtmin, periodmin;
//             SGP4Lib.elsetrec satrec;
//             SGP4Lib.elsetrec satrecx;
//             double[,] trans = new double[3, 3];
//             int sttrlen = 185;
//             StringBuilder strBuild = new StringBuilder();
// 
//             SGP4DCLib.obsrec currobsrec;
//             obsarr.resize(5000);
// 
// 
//             //conv = 0.001;
//             //timezone = 0;
// 
// 
//             // --------------------------------------------------------------------------------------
//             int err = fopen_s(&infile, filename, "r");
// 
//             if (infile == NULL)
//             {
//                 strBuild.AppendLine("Failed to open " + filename);
//             }
//             else
//             {
//                 conv = 1.0;  // set for km in case sscan doesn't work
//                 conv = 0.001;
//                 units= "Meters"; // default
// 
//                 // ---------------------- search the header for information -----------------------------
//                 do
//                 {
//                     fgets(longstr1, sttrlen, infile);
//                     //removeWS(longstr1);
//                     strr1 = char(longstr1[0]);
//                     strr3 = char(longstr1[4]);
//                     if (strr1 == 'D' || strr3 == 'D')   // search for units
//                     {
//                         std::string ExampleString(longstr1);
//                         std::istringstream yearDayMonthHourStringStream(ExampleString);
//                         yearDayMonthHourStringStream >> tmp >> units;
// 
//                         if (units.Contains("Kilometers"))
//                             conv = 1.0;
//                         else
//                             conv = 0.001;
//                     }
//                     if (strr1 == 'S' || strr3 == 'S')   // search for epoch time position
//                     {
//                         scanf(longstr1, "%14s %i %4s %i %3i %1c %2i %1c %lf ",
//                             &tmp, &currobsrec.day, &monstr, &currobsrec.year, &currobsrec.hr, &tmpc,
//                             &currobsrec.min, &tmpc1, &currobsrec.sec);
//                         currobsrec.mon = MathTimeLibr.getIntMonth(monstr);
//                         MathTimeLibr.jday(currobsrec.year, currobsrec.mon, currobsrec.day,
//                             currobsrec.hr, currobsrec.min, currobsrec.sec, satrec.jdsatepoch, satrec.jdsatepochF);
//                         currobsrec.jd = satrec.jdsatepoch;
//                         currobsrec.jdf = satrec.jdsatepochF;
//                         jdepoch = satrec.jdsatepoch;
//                         jdepochf = satrec.jdsatepochF;
//                     }
//                     if (strr1 == 'C' || strr3 == 'C')    // search for coordinate system
//                     {
//                         strr2 = longstr1[5];
//                         if ((longstr1[2] == 'o') || (strr2 == 'o'))
//                         {
//                             //std::string coord;
//                             scanf(longstr1, "%s %s ", &tmp, &coords);
//                             // dav move the eop initialization
//                             //              findatmosparam( jde, mfme, interp, fluxtype, f81type, inputtype, spwarr, jdspwstart,
//                             //                f107, f107bar, ap, avgap, aparr, kp, sumkp, kparr );
//                             //                              sethelp(iauhelp,'n');
//                             timezone = 0;
//                         }
//                     }
// 
//                 } while ((strr1 != 'E' && strr3 != 'E') && (feof(infile) == 0));
// 
//                 fgets(longstr1, sttrlen, infile);  // blank line
// 
//             }
// 
// 
//             // ---------------------- search the file for observations -----------------------------
//             // ---- read more than the points for the fit span to permit prediction analysis if possible
//             ii = 0;
//             while ((feof(infile) == 0) && (ii < obsread))   // read more in for comparison later
//             {
//                 fgets(longstr, sttrlen, infile);
//                 strncpy_s(strr, &longstr[0], 1);
//                 strr[1] = '\0';
//                 if (strcmp(strr, "E") != 0)
//                 {
//                     sscanf_s(longstr, "%lf %lf %lf %lf %lf %lf %lf ",
//                         &dtsec, &currobsrec.x, &currobsrec.y, &currobsrec.z, &currobsrec.xdot, &currobsrec.ydot, &currobsrec.zdot);
// 
//                     rpos[0] = currobsrec.x * conv;   // get inputs in km
//                     rpos[1] = currobsrec.y * conv;
//                     rpos[2] = currobsrec.z * conv;
//                     vpos[0] = currobsrec.xdot * conv;
//                     vpos[1] = currobsrec.ydot * conv;
//                     vpos[2] = currobsrec.zdot * conv;
// 
//                     if (ii == 0) // fix the jd epoch in case it starts later than 0.0 sec
//                     {
//                         dtseczero = dtsec;  // assumes epoch point is the first in the file
//                         MathTimeLibr.invjday(currobsrec.jd, currobsrec.jdf + dtsec / 86400.0, 
//                             out year, out mon, out day, out hr, out min, out sec);
//                         MathTimeLibr.jday(year, mon, day, hr, min, sec, out jdepoch, out jdepochf);
//                     }
// 
//                     MathTimeLibr.invjday(jdepoch, jdepochf + (dtsec - dtseczero) / 86400.0,
//                         out year, out mon, out day, out hr, out min, out sec);
//                     currobsrec.year = year;
//                     currobsrec.mon = mon;
//                     currobsrec.day = day;
//                     currobsrec.hr = hr;
//                     currobsrec.min = min;
//                     currobsrec.sec = sec;
//                     MathTimeLibr.jday(year, mon, day, 0, 0, 0.0, out currobsrec.jd, out currobsrec.jdf);
//                     mfme = hr * 60 + min + sec / 60.0;
// 
//                     // convert coordinate system as necessary
//                     // obs need to be in eci
//                     // do in iau-76/fk5 so not need iau06 parameters
//                     jdxysstart = 0.0;
//                     AstroLib.eOpt opt = e80;
//                     if (coords.Equals("J2000") || coords.Equals("MeanOfEpoch"))
//                     {
//                         EOPSPWLibr.findeopparam(currobsrec.jd, currobsrec.jdf, 's', eoparr, jdeopstart, 
//                             out dut1, out dat, out lod, out xp, out yp, out ddpsi, out ddeps,
//                         out ddx, out ddy, out x, out y, out s);
//                         MathTimeLibr.convtime(year, mon, day, hr, min, sec, timezone, dut1, dat,
//                             out ut1, out tut1, out jdut1, out jdut1f, out utc, out tai, out tt, out ttt, 
//                             out jdtt, out jdftt, out tcg, out tdb, out ttdb, out jdtdb, out jdtdbf, out tcb);
//                         //coordFK5::itrf_gcrf(ritrf, vitrf, aitrf, eFrom, rpos, vpos, apos, iau80rec, ttt, jdut1, lod, xp, yp, 2, ddpsi, ddeps, deltapsi, deltaeps, trans);
//                         //coordFK5::teme_ecef(rteme, vteme, ateme, eTo, ritrf, vitrf, aitrf, ttt, jdut1, lod, xp, yp, eqeterms);
//                         AstroLibr.eci_mod(ref rpos, ref vpos, MathTimeLib.eto, ref reci, ref veci, opt, ttt, jdut1 + jdut1f);
//                     }
//                     else
//                         if (coords.Equals("Fixed"))
//                     {
//                         EOPSPWLibr.findeopparam(currobsrec.jd, currobsrec.jdf, 's', EOPSPWLibr.eopdata, mjdeopstart + 2400000.5,
//                            out dut1, out dat, out lod, out xp, out yp, out ddpsi, out ddeps, out ddx, out ddy);
// 
//                         MathTimeLibr.convtime(year, mon, day, hr, min, sec, timezone, dut1, dat,
//                             out ut1, out tut1, out jdut1, out jdut1f, out utc, out tai, out tt, out ttt, 
//                             out jdtt, out jdftt, out tcg, out tdb, out ttdb, out jdtdb, out jdtdbf, out tcb);
//                         AstroLibr.eci_ecef(ref rpos, ref vpos, MathTimeLibr.eto, ref reci, ref veci, opt, iau80rec, iau06rec, jdtt, jdftt,
//                             jdut1 + jdut1f, jdxysstart, lod, xp, yp, ddpsi, ddeps, ddx, ddy, trans);
//                         //coordFK5::teme_ecef( rteme, vteme, ateme, eTo,rpos, vpos, apos, ttt, jdut1, lod, xp, yp, eqeterms);
//                     }
//                     else
//                     {
//                         for (j = 0; j < 3; j++)
//                         {
//                             reci[j] = rpos[j];
//                             veci[j] = vpos[j];
//                         }
//                     }
// 
//                     // now reset the obs vector to be eci, (teme?)
//                     currobsrec.x = reci[0];
//                     currobsrec.y = reci[1];
//                     currobsrec.z = reci[2];
//                     currobsrec.xdot = veci[0];
//                     currobsrec.ydot = veci[1];
//                     currobsrec.zdot = veci[2];
// 
//                     // setup the initial satrec TLE estimate if sgp4 used 
//                     // this will be the nominal 
//                     if (ii == 0 && proptype == 's')
//                     {
//                         AstroLibr.rv2coe(reci, veci, out p, out a, out ecc, out incl, out raan, out argp, out nu, 
//                             out m, out eccanom, out arglat, out truelon, out lonper);
//                         satrec.no_kozai = 60.0 * Math.Sqrt(3.986008e5 / (a * a * a));  // rad per min, sgp4 wgs-72 mu value
// 
//                         satrec.satnum = "55555";  // set some number here
//                         satrec.intldesg= "testsat";
//                         satrec.classification = 'U';
// 
//                         // set the time span based on the period size,
//                         // but have to calculate on second iteration when dtsec is known
//                         // so calc here before satrec.a changes on the 1st iteration
//                         periodmin = (2.0 * Math.PI * Math.Sqrt(Math.Pow(a, 3) / 398600.8)) / 60.0;  // minutes
//                         satrec.a = a / 6378.135;  // ER
// 
//                         satrec.ecco = ecc;
//                         satrec.inclo = incl;
//                         satrec.nodeo = raan;
//                         satrec.argpo = argp;
//                         satrec.mo = m;
//                         if (year >= 2000)
//                             satrec.epochyr = year - 2000;
//                         else
//                             satrec.epochyr = year - 1900;
//                         MathTimeLibr.finddays(year, mon, day, hr, min, sec, out satrec.epochdays);
//                         satrecx = satrec;    // store the answer for later
//                     }  // if ii = 0
//                     else
//                     {
//                         if (ii == 0)
//                         {
//                             for (j = 0; j < 3; j++)
//                             {
//                                 xnom[j] = reci[j];
//                                 xnom[j + 3] = veci[j];
//                             }
//                         }
//                     }
// 
//                     if (ii == 1)
//                     {
//                         // lastob = int(lastob / (dtsec / 60.0) * .5 * int( periodmin/ lastob )); // reset this based on the ephemeris file spacing
//                         // lastob = 200; // force this/. 901 is at 43 sec, 906 is at 15 min, ice is at 30 sec, gps is at 15 min
//                         //strBuild.AppendLine(" %3i with spacing of %lf min is %lf min \n", lastob, (dtsec - dtseczero) / 60.0, lastob*(dtsec - dtseczero) / 60.0);
//                     }
// 
//                     dtmin = (dtsec - dtseczero) / 60.0;  // in min
//                     currobsrec.dtmin = dtmin;
//                     currobsrec.sennum = 1;
//                     currobsrec.satnum = 55555;
//                     currobsrec.obstype = 4;  // xyz obs
//                     obsarr[ii] = currobsrec;
// 
//                     ii = ii + 1;   // increase the obs by one even if obsstep is different
// 
//                     if (ii % 500 == 0)
//                         strBuild.AppendLine(@"over 500 obs \n");
//                 }   // if strr
// 
//                 // reread any additional steps if needed
//                 if (obsskip > 1)
//                 {
//                     for (kk = 1; kk < obsskip; kk++)
//                     {
//                         if (feof(infile) == 0)
//                             fgets(longstr, sttrlen, infile);
//                     }
//                 }
// 
//             } // while ii through producing the observation vectors
// 
//             fclose(infile);
// 
//         }  // readobs






        /* ----------------------------------------------------------------------------------------
        *  processrazel
        *
        *  New C# port of the MATLAB Ex10_4.m / findatwaatwb.m range-az-el differential
        *  correction example (Vallado, Fundamentals of Astrodynamics 2022, sec 10.6).
        *  Reads a GEOS-format range/az/el obs file, seeds a (deliberately off) nominal
        *  state vector, and calls SGP4DCLib.leastsquares (proptype -> Kepler two-body
        *  path) to iterate to a converged solution, matching the MATLAB test case so it
        *  can become a TestAll regression case.
        *
        *  input file format assumed (matches geos6a.inp / ex10_4.m filedat columns):
        *      col1-3  (unused id fields)
        *      col4    year            col5   month      col6   day
        *      col7    hour            col8   minute     col9   sec
        *      col10   range (km)      col11  az (deg)   col12  el (deg)
        *  whitespace-delimited, one observation per line.
        *
        *  NOTE: verify this column layout against the actual geos6a.inp on disk before
        *  trusting results - I don't have that file to check against directly.
        * ---------------------------------------------------------------------------------------- */
        // vestigial - SGP4DCLib.leastsquares takes this by ref but never actually uses it for the
        // Kepler/two-body path (EOP data is (re)loaded internally from a hardcoded path instead).
        private double[] eoparrUnused = new double[1];

        // The exe runs from <ProjectRoot>\bin\Debug\ (or bin\Release\); BaseDirectory always reflects
        // the exe's actual location regardless of how the process was launched (unlike the ambiguous
        // "current directory", which differs between running under the VS debugger and double-clicking
        // the exe directly) - so this reliably resolves to <ProjectRoot>\tests either way.
        private string GetTestsDir()
        {
            string exeDir = AppDomain.CurrentDomain.BaseDirectory;
            string testsDir = Path.GetFullPath(Path.Combine(exeDir, "..", "..", "tests"));
            return testsDir;
        }

        // Classical Jacobi eigenvalue algorithm for a real symmetric matrix - used to reproduce
        // MATLAB's eig(p) output for the covariance matrix (eigenaxes/eigenvalues) so a run's full
        // report can be diffed line-for-line against the Ex10_4.m reference. Returns eigenvalues
        // ascending and eigenvectors as columns, matching MATLAB's default eig() ordering/convention.
        private void JacobiEigenSymmetric(double[,] aIn, int n, out double[] eigenvalues, out double[,] eigenvectors)
        {
            double[,] a = (double[,])aIn.Clone();
            double[,] v = new double[n, n];
            for (int i = 0; i < n; i++) v[i, i] = 1.0;

            for (int sweep = 0; sweep < 100; sweep++)
            {
                double off = 0.0;
                for (int p = 0; p < n; p++)
                    for (int q = p + 1; q < n; q++)
                        off += a[p, q] * a[p, q];
                if (off < 1e-24) break;

                for (int p = 0; p < n - 1; p++)
                {
                    for (int q = p + 1; q < n; q++)
                    {
                        if (Math.Abs(a[p, q]) < 1e-300) continue;
                        double theta = (a[q, q] - a[p, p]) / (2.0 * a[p, q]);
                        double t = (theta == 0.0) ? 1.0 :
                            Math.Sign(theta) / (Math.Abs(theta) + Math.Sqrt(theta * theta + 1.0));
                        double c = 1.0 / Math.Sqrt(t * t + 1.0);
                        double s = t * c;
                        double app = a[p, p], aqq = a[q, q], apq = a[p, q];
                        a[p, p] = c * c * app - 2 * s * c * apq + s * s * aqq;
                        a[q, q] = s * s * app + 2 * s * c * apq + c * c * aqq;
                        a[p, q] = 0.0; a[q, p] = 0.0;
                        for (int i = 0; i < n; i++)
                        {
                            if (i != p && i != q)
                            {
                                double aip = a[i, p], aiq = a[i, q];
                                a[i, p] = c * aip - s * aiq;
                                a[p, i] = a[i, p];
                                a[i, q] = s * aip + c * aiq;
                                a[q, i] = a[i, q];
                            }
                            double vip = v[i, p], viq = v[i, q];
                            v[i, p] = c * vip - s * viq;
                            v[i, q] = s * vip + c * viq;
                        }
                    }
                }
            }

            double[] evals = new double[n];
            for (int i = 0; i < n; i++) evals[i] = a[i, i];
            int[] order = Enumerable.Range(0, n).OrderBy(i => evals[i]).ToArray();

            eigenvalues = new double[n];
            eigenvectors = new double[n, n];
            for (int newIdx = 0; newIdx < n; newIdx++)
            {
                int oldIdx = order[newIdx];
                eigenvalues[newIdx] = evals[oldIdx];
                for (int r = 0; r < n; r++)
                    eigenvectors[r, newIdx] = v[r, oldIdx];
            }
        }

        private void ReportDiff(StringBuilder strBuild, string label, double[] refr, double[] rFinal)
        {
            double dx0 = refr[0] - rFinal[0];
            double dy0 = refr[1] - rFinal[1];
            double dz0 = refr[2] - rFinal[2];
            double mag = Math.Sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0);
            // Matches Ex10_4.m's exact wording/format: "diff <label>  %16.8f  %16.8f  %16.8f  %16.8f  km"
            strBuild.AppendLine(string.Format(CultureInfo.InvariantCulture,
                "diff {0}  {1,16:F8}  {2,16:F8}  {3,16:F8}  {4,16:F8}  km", label, dx0, dy0, dz0, mag));
        }

        private void processrazel(char typeans, char proptype, int obsskip, int obsread,
            double[] xnomSeed, string outputFolder)
        {
            StringBuilder strBuild = new StringBuilder();
            strBuild.AppendLine("=============== processrazel (Ex10_4 / GEOS range-az-el DC) ===============");

            // ---- locate the obs file - check the tests folder first, then a couple of legacy spots ----
            string outputFolderSafe = outputFolder == null ? "" : outputFolder;
            string testsDir = GetTestsDir();
            string[] candidates = new string[]
            {
                Path.Combine(testsDir, "geos6a.inp"),
                Path.Combine(testsDir, "GEOS6a.INP"),
                "geos6a.inp",
                Path.Combine("..", "geos6a.inp"),
                Path.Combine(outputFolderSafe, "geos6a.inp"),
            };
            string fname = null;
            for (int c = 0; c < candidates.Length; c++)
            {
                if (File.Exists(candidates[c]))
                {
                    fname = candidates[c];
                    break;
                }
            }
            if (fname == null)
            {
                MessageBox.Show("Could not find geos6a.inp in any of:\n" + string.Join("\n", candidates) +
                    "\nPlace the obs file in the tests folder (" + testsDir + ") and try again.");
                return;
            }

            // ---- Kaena Point sensor, AAS-paper values used directly in Ex10_4.m ----
            //      (SGP4DCLib.getsensorparams sennum 932 has been updated to match these;
            //       GEOS6a.INP's obstype/sennum columns (fields 0 and 2) confirm 2/932)
            double latgd = 21.572056 * Math.PI / 180.0;
            double lon = -158.266578 * Math.PI / 180.0;
            double alt = 0.3002;  // km
            double[] rsecef;
            double[] vsecef;
            AstroLibr.site(latgd, lon, alt, out rsecef, out vsecef);

            // ---- read the obs file, honoring obsskip / obsread from the GUI ----
            string[] allLines = File.ReadAllLines(fname);
            List<string> rawLines = new List<string>();
            for (int li = 0; li < allLines.Length; li++)
            {
                string t = allLines[li].Trim();
                if (t.Length > 0)
                    rawLines.Add(t);
            }

            List<string> picked = new List<string>();
            for (int idx = 0; idx < rawLines.Count && picked.Count < obsread; idx += obsskip)
                picked.Add(rawLines[idx]);

            int numobs = picked.Count;
            if (numobs < 3)
            {
                MessageBox.Show("Only found " + numobs + " usable observations in " + fname +
                    " - need at least 3 (and generally 10+) for a useful fit.");
                return;
            }

            // obsrecarr is sized numobs+1 and left 1-based (index 0 unused by the fit loop) to
            // match the 1-based firstob..lastob loop convention already used inside
            // SGP4DCLib.findatwaatwb / leastsquares.
            SGP4DCLib.obsrec[] obsrecarr = new SGP4DCLib.obsrec[numobs + 1];
            for (int j = 1; j <= numobs; j++)
            {
                string[] f = picked[j - 1].Split((char[])null, StringSplitOptions.RemoveEmptyEntries);
                int obstype = int.Parse(f[0], CultureInfo.InvariantCulture);
                int sennum = int.Parse(f[2], CultureInfo.InvariantCulture);
                int year = int.Parse(f[3], CultureInfo.InvariantCulture);
                int mon = int.Parse(f[4], CultureInfo.InvariantCulture);
                int day = int.Parse(f[5], CultureInfo.InvariantCulture);
                int hr = int.Parse(f[6], CultureInfo.InvariantCulture);
                int min = int.Parse(f[7], CultureInfo.InvariantCulture);
                double sec = double.Parse(f[8], CultureInfo.InvariantCulture);
                // Ex10_4.m subtracts known calibration biases from the raw file values before use -
                // matching that exactly (confirmed against a direct obsrecarr(1).az/obsrecarr(2).az
                // comparison: both were off from MATLAB by the exact same 0.0081 deg, i.e. this bias,
                // which is why it showed up as a constant offset rather than something time-varying).
                double rng = double.Parse(f[9], CultureInfo.InvariantCulture) - 0.08;
                double az = double.Parse(f[10], CultureInfo.InvariantCulture) * Math.PI / 180.0 - 0.0081 * Math.PI / 180.0;
                double el = double.Parse(f[11], CultureInfo.InvariantCulture) * Math.PI / 180.0 - 0.0045 * Math.PI / 180.0;

                double jd, jdf;
                MathTimeLibr.jday(year, mon, day, hr, min, sec, out jd, out jdf);

                SGP4DCLib.obsrec rec = new SGP4DCLib.obsrec();
                rec.sennum = sennum;
                rec.obstype = obstype;   // range - az - el (2, per GEOS6a.INP field 0)
                rec.year = year;
                rec.mon = mon;
                rec.day = day;
                rec.hr = hr;
                rec.min = min;
                rec.sec = sec;
                rec.jd = jd;
                rec.jdf = jdf;
                rec.rng = rng;
                rec.az = az;
                rec.el = el;
                rec.rsecef = (double[])rsecef.Clone();
                rec.vsecef = (double[])vsecef.Clone();
                obsrecarr[j] = rec;
            }

            // index 0 is peeked at (for obstype) by both findatwaatwb and leastsquares even though
            // the main loop is 1-based (firstob..lastob) - alias it to the first real observation
            // so that peek doesn't hit a null obsrec.
            obsrecarr[0] = obsrecarr[1];

            // ---- seed the nominal state and run the least-squares differential correction ----
            double[] xnom = new double[7];
            Array.Copy(xnomSeed, xnom, 6);

            strBuild.AppendLine("input: ");
            strBuild.AppendLine(string.Format(CultureInfo.InvariantCulture,
                "rnom {0,16:F8} {1,16:F8} {2,16:F8} km vnom {3,16:F8} {4,16:F8} {5,16:F8} km ",
                xnom[0], xnom[1], xnom[2], xnom[3], xnom[4], xnom[5]));

            SGP4Lib.elsetrec satrec = new SGP4Lib.elsetrec();
            double jd0 = obsrecarr[1].jd;
            double jdf0 = obsrecarr[1].jdf;
            int loops;
            double[,] x, dx, atwai, atwa, atwb;

            SGP4DCLibr.leastsquares(
                0.01, 0.01, 0.0001,
                SGP4Lib.gravconsttype.wgs84,
                's', 0.0,
                ref satrec,
                typeans, 'v', proptype,
                1, numobs, 6,
                ref obsrecarr, out loops,
                ref xnom, jd0, jdf0,
                1, 0.044, 0.02, 'a',
                new EOPSPWLib.iau80Class(), new EOPSPWLib.iau06Class(),
                ref eoparrUnused,
                out x, out dx, out atwai, out atwa, out atwb);

            if (SGP4DCLibr.strbuild.Length > 0)
            {
                strBuild.AppendLine("--- leastsquares diagnostics ---");
                strBuild.Append(SGP4DCLibr.strbuild.ToString());
                SGP4DCLibr.strbuild.Clear();
            }

            double[] rFinal = new double[] { xnom[0], xnom[1], xnom[2] };
            double[] vFinal = new double[] { xnom[3], xnom[4], xnom[5] };

            double p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper;
            AstroLibr.rv2coe(rFinal, vFinal, out p, out a, out ecc, out incl,
                out raan, out argp, out nu, out m, out arglat, out truelon, out lonper);
            double rad = 180.0 / Math.PI;

            // ---- Ex10_4.m-format final report, so this can be diffed line-for-line against a saved
            // MATLAB run. Labels, ordering, and the (mislabeled - it's actually meters, not km, per
            // the *1000 below) "cov r, 1 sigma" units are all kept exactly as the reference script
            // prints them, rather than "corrected", specifically so future runs stay diffable.
            strBuild.AppendLine("output: ");
            strBuild.AppendLine(string.Format(CultureInfo.InvariantCulture,
                "r1 {0,16:F8} {1,16:F8} {2,16:F8} km v1 {3,16:F8} {4,16:F8} {5,16:F8} km ",
                rFinal[0], rFinal[1], rFinal[2], vFinal[0], vFinal[1], vFinal[2]));
            strBuild.AppendLine("          p km       a km      ecc      incl deg     raan deg     argp deg" +
                "      nu deg      m deg      arglat   truelon    lonper");
            strBuild.AppendLine(string.Format(CultureInfo.InvariantCulture,
                "coes {0,11:F4} {1,11:F4} {2,13:F9} {3,13:F7} {4,11:F5} {5,11:F5} {6,11:F5} {7,11:F5} {8,11:F5} {9,11:F5} {10,11:F5}",
                p, a, ecc, incl * rad, raan * rad, argp * rad, nu * rad, m * rad, arglat * rad, truelon * rad, lonper * rad));

            // ---- reference answers from Ex10_4.m, for a quick pass/fail sanity check ----
            double[] rGtds = new double[] { 5748.6070, 2679.9324, 3442.7902 };
            double[] rOdtkShort = new double[] { 5746.469300, 2679.944362, 3442.416160 };
            double[] rOdtkLong = new double[] { 5748.398438, 2680.135601, 3442.296315 };
            ReportDiff(strBuild, "odtk short bad init ", rOdtkShort, rFinal);
            ReportDiff(strBuild, "gtds w bad init     ", rGtds, rFinal);
            ReportDiff(strBuild, "odtk long w bad init", rOdtkLong, rFinal);

            if (atwai != null)
            {
                strBuild.AppendLine("covariance p = ");
                for (int rr = 0; rr < 6; rr++)
                {
                    StringBuilder rowStr = new StringBuilder();
                    for (int cc = 0; cc < 6; cc++)
                        rowStr.Append(string.Format(CultureInfo.InvariantCulture, "{0,14:E6} ", atwai[rr, cc]));
                    strBuild.AppendLine(rowStr.ToString().TrimEnd());
                }
                strBuild.AppendLine(string.Format(CultureInfo.InvariantCulture,
                    "cov r, 1 sigma  {0,16:F8}  {1,16:F8}  {2,16:F8}  km",
                    Math.Sqrt(Math.Abs(atwai[0, 0])) * 1000.0,
                    Math.Sqrt(Math.Abs(atwai[1, 1])) * 1000.0,
                    Math.Sqrt(Math.Abs(atwai[2, 2])) * 1000.0));

                double[] eigenvalues;
                double[,] eigenvectors;
                JacobiEigenSymmetric(atwai, 6, out eigenvalues, out eigenvectors);
                for (int col = 0; col < 6; col++)
                {
                    StringBuilder axRow = new StringBuilder("eigenaxes  ");
                    for (int rr = 0; rr < 6; rr++)
                        axRow.Append(string.Format(CultureInfo.InvariantCulture, "{0,16:F8}  ", eigenvectors[rr, col]));
                    axRow.Append("km^2");
                    strBuild.AppendLine(axRow.ToString());
                }
                double[] eigenvaluesM = new double[6];
                for (int i = 0; i < 6; i++)
                    eigenvaluesM[i] = Math.Sqrt(Math.Abs(eigenvalues[i])) * 1000.0;  // m, per Ex10_4.m
                StringBuilder evRow = new StringBuilder("eigenvalues   ");
                for (int i = 0; i < 6; i++)
                    evRow.Append(string.Format(CultureInfo.InvariantCulture, "{0,16:F8}  ", eigenvaluesM[i]));
                evRow.Append("m");
                strBuild.AppendLine(evRow.ToString());
                double magPos = Math.Sqrt(eigenvaluesM[0] * eigenvaluesM[0] + eigenvaluesM[1] * eigenvaluesM[1] + eigenvaluesM[2] * eigenvaluesM[2]);
                double magVel = Math.Sqrt(eigenvaluesM[3] * eigenvaluesM[3] + eigenvaluesM[4] * eigenvaluesM[4] + eigenvaluesM[5] * eigenvaluesM[5]);
                strBuild.AppendLine(string.Format(CultureInfo.InvariantCulture, "mag eigenvalues  {0,16:F8} m", magPos));
                strBuild.AppendLine(string.Format(CultureInfo.InvariantCulture, "mag eigenvalues  {0,16:F8} m/s", magVel));
            }

            strBuild.AppendLine("");
            strBuild.AppendLine("obs file    " + fname);
            strBuild.AppendLine("numobs used " + numobs + "   loops " + loops);

            string outText = strBuild.ToString();
            try
            {
                Directory.CreateDirectory(testsDir);
                string outPath = Path.Combine(testsDir, "testdc_razel.out");
                File.WriteAllText(outPath, outText);
                MessageBox.Show(outText + "\n\nFull results written to:\n" + outPath);
            }
            catch (Exception ex)
            {
                MessageBox.Show(outText + "\n\n(could not write output file: " + ex.Message + ")");
            }
        }

        /* ----------------------------------------------------------------------------------------
        *  processtle
        *
        *  Self-consistency test for the SGP4 (proptype 's') differential-correction path.
        *  For each TLE in a file: propagate it with SGP4 to synthesize a set of "truth"
        *  position/velocity observations (obstype 4) spanning one orbit, perturb the state
        *  slightly, and run leastsquares to see whether it converges back to the truth
        *  ephemeris. This is the first thing in this codebase to actually exercise the
        *  proptype=='s' branch of findatwaatwb/leastsquares end-to-end.
        *
        *  KNOWN LIMITATION: statesize is 6 (position+velocity only), not 7. state2satrec's
        *  'v' (vector) case - the one leastsquares actually uses here - never copies xnom[6]
        *  into satrec.bstar, so a bstar unknown would just silently ride along at whatever
        *  satrec.bstar happened to already be instead of being fit. Recovering bstar would
        *  need that gap closed in SGP4DCLib.cs first; flagging rather than papering over it.
        * ---------------------------------------------------------------------------------------- */
        private void processtle(char typeans, string outputFolder, int lastob)
        {
            StringBuilder strBuild = new StringBuilder();
            strBuild.AppendLine("=============== processtle (SGP4/TLE self-consistency DC test) ===============");

            if (lastob < 4) lastob = 36;  // points across one orbit; the razel obsread default (20) is a bit coarse but usable

            string outputFolderSafe = outputFolder == null ? "" : outputFolder;
            string testsDir = GetTestsDir();
            string[] candidateDirs = new string[] { testsDir, "", outputFolderSafe };
            string fname = null;
            for (int c = 0; c < candidateDirs.Length && fname == null; c++)
            {
                if (!Directory.Exists(string.IsNullOrEmpty(candidateDirs[c]) ? "." : candidateDirs[c]))
                    continue;
                string[] hits = Directory.GetFiles(string.IsNullOrEmpty(candidateDirs[c]) ? "." : candidateDirs[c], "*.tle");
                if (hits.Length > 0)
                    fname = hits[0];
            }
            if (fname == null)
            {
                MessageBox.Show("Could not find a *.tle file in the tests folder (" + testsDir +
                    "), next to the exe, or in the Output file location folder.\n" +
                    "Drop a file with standard 2-line (or 3-line with a name line) TLEs, extension .tle, there.");
                return;
            }

            // ---- parse every TLE pair in the file (name line, if present, is just skipped) ----
            string[] allLines = File.ReadAllLines(fname);
            List<string[]> tlePairs = new List<string[]>();
            for (int li = 0; li < allLines.Length - 1; li++)
            {
                string a = allLines[li].TrimEnd();
                string b = allLines[li + 1].TrimEnd();
                if (a.StartsWith("1 ") && b.StartsWith("2 "))
                {
                    tlePairs.Add(new string[] { a, b });
                    li++;  // consume both lines of the pair
                }
            }
            if (tlePairs.Count == 0)
            {
                MessageBox.Show("Found " + fname + " but no valid TLE pairs in it (looked for a '1 ...' line " +
                    "immediately followed by a '2 ...' line).");
                return;
            }
            int maxSats = 50;
            bool truncated = tlePairs.Count > maxSats;
            if (truncated)
                tlePairs = tlePairs.GetRange(0, maxSats);

            // nutation coefficients only - real EOP time-series data isn't needed here since
            // teme_eci only takes ttt/ddpsi/ddeps (no UT1/polar motion), and using a plain
            // JD-based ttt (skipping the ~tens-of-seconds UTC-TT leap-second offset) is well
            // within precession/nutation precision needs for this test - the same simplification
            // findatwaatwb itself already makes by hardcoding ddpsi=ddeps=0.0 for this call.
            string nutLoc = @"D:\Codes\LIBRARY\DataLib\nut80.dat";
            EOPSPWLib.iau80Class iau80arr;
            EOPSPWLibr.iau80in(nutLoc, out iau80arr);

            int nTested = 0;
            int nConverged = 0;

            foreach (string[] pair in tlePairs)
            {
                nTested++;
                double startmfe, stopmfe, deltamin;
                SGP4Lib.elsetrec satrecTruth;
                SGP4Libr.twoline2rv(pair[0], pair[1], 'c', 'd', 'a', SGP4Lib.gravconsttype.wgs72,
                    out startmfe, out stopmfe, out deltamin, out satrecTruth);

                double periodmin = 2.0 * Math.PI / satrecTruth.no_kozai;
                double dtmin = periodmin / lastob;

                SGP4DCLib.obsrec[] obsrecarr = new SGP4DCLib.obsrec[lastob + 2];
                for (int ii = 1; ii <= lastob + 1; ii++)
                {
                    double tsince = dtmin * (ii - 1);
                    double[] rteme = new double[3];
                    double[] vteme = new double[3];
                    SGP4Libr.sgp4(ref satrecTruth, tsince, rteme, vteme);

                    double jdSample = satrecTruth.jdsatepoch;
                    double jdfSample = satrecTruth.jdsatepochF + tsince / 1440.0;
                    double ttt = (jdSample + jdfSample - 2451545.0) / 36525.0;
                    double[] reci = new double[3];
                    double[] veci = new double[3];
                    AstroLibr.teme_eci(ref rteme, ref vteme, iau80arr, MathTimeLib.Edirection.eto, ttt, 0.0, 0.0,
                        ref reci, ref veci);

                    int yy, mo, dd, hh, mi;
                    double ss;
                    MathTimeLibr.invjday(jdSample, jdfSample, out yy, out mo, out dd, out hh, out mi, out ss);
                    double jdRec, jdfRec;
                    MathTimeLibr.jday(yy, mo, dd, hh, mi, ss, out jdRec, out jdfRec);

                    SGP4DCLib.obsrec rec = new SGP4DCLib.obsrec();
                    rec.sennum = 1;   // getsensorparams case 1 carries the x/y/z/xdot/ydot/zdot/bstar noise values
                    rec.obstype = 4;  // position/velocity vector obs
                    rec.year = yy; rec.mon = mo; rec.day = dd; rec.hr = hh; rec.min = mi; rec.sec = ss;
                    rec.jd = jdRec; rec.jdf = jdfRec;
                    rec.x = reci[0]; rec.y = reci[1]; rec.z = reci[2];
                    rec.xdot = veci[0]; rec.ydot = veci[1]; rec.zdot = veci[2];
                    rec.bstar = satrecTruth.bstar;
                    obsrecarr[ii] = rec;
                }
                obsrecarr[0] = obsrecarr[1];
                int numobs = lastob + 1;

                // seed the nominal a little off truth so the DC has something to converge on
                double[] xnom = new double[7];
                xnom[0] = obsrecarr[1].x * 1.0002;
                xnom[1] = obsrecarr[1].y * 1.0002;
                xnom[2] = obsrecarr[1].z * 0.9998;
                xnom[3] = obsrecarr[1].xdot * 0.999;
                xnom[4] = obsrecarr[1].ydot * 1.001;
                xnom[5] = obsrecarr[1].zdot * 0.999;

                SGP4Lib.elsetrec satrecFit = new SGP4Lib.elsetrec();
                satrecFit.satnumStr = satrecTruth.satnumStr;
                satrecFit.operationmode = 'a';

                int loops;
                double[,] x, dx, atwai, atwa, atwb;
                SGP4DCLibr.leastsquares(
                    0.001, 0.0000001, 0.0002,
                    SGP4Lib.gravconsttype.wgs72,
                    's', 0.0,
                    ref satrecFit,
                    typeans, 'v', 's',
                    1, numobs, 6,
                    ref obsrecarr, out loops,
                    ref xnom, obsrecarr[1].jd, obsrecarr[1].jdf,
                    1, 0.044, 0.02, 'a',
                    new EOPSPWLib.iau80Class(), new EOPSPWLib.iau06Class(),
                    ref eoparrUnused,
                    out x, out dx, out atwai, out atwa, out atwb);

                if (SGP4DCLibr.strbuild.Length > 0)
                {
                    strBuild.AppendLine("--- leastsquares diagnostics ---");
                    strBuild.Append(SGP4DCLibr.strbuild.ToString());
                    SGP4DCLibr.strbuild.Clear();
                }

                double dxp = xnom[0] - obsrecarr[1].x;
                double dyp = xnom[1] - obsrecarr[1].y;
                double dzp = xnom[2] - obsrecarr[1].z;
                double magm = Math.Sqrt(dxp * dxp + dyp * dyp + dzp * dzp) * 1000.0;  // m

                strBuild.AppendLine(string.Format(CultureInfo.InvariantCulture,
                    "sat {0,-10} loops {1,2}  final pos residual {2,10:F3} m", satrecTruth.satnumStr, loops, magm));
                if (magm < 1.0)
                    nConverged++;
            }

            strBuild.AppendLine("");
            strBuild.AppendLine(nConverged + " / " + nTested + " satellites converged to < 1 m position residual");
            if (truncated)
                strBuild.AppendLine("(file had more than " + maxSats + " TLEs - only the first " + maxSats + " were tested)");

            string outText = strBuild.ToString();
            try
            {
                Directory.CreateDirectory(testsDir);
                string outPath = Path.Combine(testsDir, "testdc_tle.out");
                File.WriteAllText(outPath, outText);
                MessageBox.Show(outText + "\n\nFull results written to:\n" + outPath);
            }
            catch (Exception ex)
            {
                MessageBox.Show(outText + "\n\n(could not write output file: " + ex.Message + ")");
            }
        }

        // ---- UNCONVERTED LEGACY CODE (processrazel_OLD) commented out below pending C# port ----
//         /*       ----------------------------------------------------------------     */
//         void processrazel
//         (
//             char batchmode, char proptype, int derivType, double cdam, double cram,
//             EOPSPWLib.iau80Class iau80arr,
//             EOPSPWLib.iau06Class iau06arr,
//             double[] eoparr,
//             double jdeopstart, double jdeopstartf,
//             double[] spwarr,
//             double jdspwstart, double jdspwstartf,
//             double[] jpldearr,
//         char jplopt,
//             double jdjpldestart,
//             double jdepoch, double jdepochf, char interp,
//             int firstob, int obsskip, int obsread,
//             double[] obsarr,
//             out double[] xnom
//         )
//         {
//            // FILE* infile, *outfile, *outfile1;
//             double[] obsrecarr;
//             SGP4Lib.elsetrec satrec;
//             long search_satno = 0;
//             string longstr1, longstr2;
//             string strr;
//             string coords;
//             char typerun;
//             char statetype, typeans, opsmode;
//             int statesize;
//             string fname;
//             double percentchg, deltaamtchg, rmsepsilon;
// 
//             int day, mon, year, hr, min, loops, i;
//             double conv, sec;
//             obsrecarr.resize(5000);
// 
//             double dut1, ut1, tut1, jdut1, jdut1f, utc, tai, tt, ttt, jdtt, jdftt, tcg, tdb, jdxysstart,
//                 ttdb, jdtdb, jdtdbf, tcb, mfme, lod, xp, yp, ddpsi, ddeps, ddx, ddy, icrsx, y, s;
//             char keepit;
//             int timezone, dat, ktr;
//             double[,] trans = new double[3, 3];
// 
//             double rngest, zrng, zaz, zel, zrtasc, zdecl, theta, theta1, copa;
//             string error;
//             int ii, j, skipped, eqeterms;
//             double[] rpos = new double[3];
//             double[] vpos = new double[3];
//             double[] rnomTEME = new double[3];
//             double[] vnomTEME = new double[3];
//             double[] dr = new double[3];
//             double[] dv = new double[3];
//             double[] r1 = new double[3];
//             double[] r2 = new double[3];
//             double[] r3 = new double[3];
//             double[] v1 = new double[3];
//             double[] v2 = new double[3];
//             double[] v3 = new double[3];
//             double[] rteme = new double[3];
//             double[] vteme = new double[3];
//             double dtmin, jd, jdf;
//             double[] r1teme = new double[3];
//             double[] v1teme = new double[3];
//             double[] r2teme = new double[3];
//             double[] v2teme = new double[3];
//             double[] r3teme = new double[3];
//             double[] v3teme = new double[3];
//             double[] rsum = new double[3];
//             double[] vsum = new double[3];
// 
//             SGP4DCLib.obsrec currobsrec;
//             SGP4DCLib.senrec currsenrec;
// 
//             double rad = 180.0 / Math.PI;
//             double dayconv = 1.0 / 86400.0;
// 
//             StringBuilder strBuild = new StringBuilder();
// 
//             eqeterms = 2;
// 
//             statesize = 6;
//             if (batchmode == 'n') // needs a default value here in batch mode
//                 statesize = 7;
//             // set these variables up with intial dimensions
//             // these must be set here
//             // function matrices can be resized
//             double[,] x = new double[statesize, 1];
//             double[,] dx = new double[statesize, 1];
//             double[,] atwai = new double[statesize, statesize];
//             double[,] atwa = new double[statesize, statesize];
//             double[,] atwb = new double[statesize, 1];
// 
//             SGP4Lib.gravconsttype whichconst = SGP4Lib.gravconsttype.wgs72; // type of run for sgp4
//             typerun = 'c';   // c catalog (just reads tle), v verification, m manual
//             opsmode = 'a';
//             typeans = 's';
//             interp = 'l'; // set eop interpolation to linear
//                           // strBuild.AppendLine("input interp n l s \n\n");
// 
//             //     sethelp(iauhelp,'n');
//             timezone = 0;
// 
//             jdxysstart = 0.0;
// 
//             if (batchmode == 'n')
//             {
//                 statetype = 't';   // v vectors, t tle vars, e equinoctial
//                 strBuild.AppendLine("input type of elements tle or equin \n\n");
//                 statetype = getchar();
//                 fflush(stdin);
//                 strcpy_s(fname, "../geos6a.inp");
//                 strBuild.AppendLine("input obs filename: \n");
//                 //scanf_s("%s", &fname, sizeof(fname));
// 
//                 percentchg = 0.001;  // percent to change components for finite differencing
//                 strBuild.AppendLine("input percentchg for finite differencing (.001) \n");
//                 //scanf_s("%lf", &percentchg);
// 
//                 deltaamtchg = 0.0000001;  // delta amt to change components for finite differencing
//                 strBuild.AppendLine("input deltaamtchg to chk in finite differencing (.0000001) \n");
//                 //scanf_s("%lf", &deltaamtchg);
// 
//                 rmsepsilon = 0.0002;  // percent to change components for finite differencing
//                 strBuild.AppendLine("input rmsepsilon to quit (.02) \n");
//                 //scanf_s("%lf", &rmsepsilon);
// 
//                 statesize = 7;        // solve for bstar or not
//                 strBuild.AppendLine("input statesize (6 or 7) \n");
//                 //scanf_s("%i", &statesize);
//             }
// 
//             strBuild.AppendLine("typerun     %c c catalog (just reads tle), v verification, m manual \n", interp);
//             strBuild.AppendLine("interp type %c eop interpolation, l linear, s spline \n", interp);
//             strBuild.AppendLine("state type  %c v vectors, t tle vars, e equinoctial \n", statetype);
//             strBuild.AppendLine("fname       %s \n", fname);
//             strBuild.AppendLine("percentchg  %lf \n", percentchg);
//             strBuild.AppendLine("deltaamtchg %11.8lf \n", deltaamtchg);
//             strBuild.AppendLine("rmsepsilon  %lf \n", rmsepsilon);
//             strBuild.AppendLine("state size  %i \n", statesize);
// 
//             int err = fopen_s(&infile, fname, "r");
//             err = fopen_s(&outfile, "testdc.out", "w");
//             err = fopen_s(&outfile1, "testdc1.out", "w");
// 
//             fstrBuild.AppendLine(outfile, "typeans %c \n", typeans);
// 
//             if (infile == NULL)
//             {
//                 strBuild.AppendLine("Failed to open %s ", fname);  // , strerror(errno)
//             }
//             else
//             {
//                 conv = 0.001;
//                 // read header lines
//                 do
//                 {
//                     fgets(longstr1, 91, infile);
//                     strr = longstr1[0];
//                     strr[1] = '\0';
//                 } while ((strcmp(strr, "#") == 0) && (feof(infile) == 0));
// 
//                 // ----------------------- read obs from the input file -------------------------
//                 coords = "Fixed";  // assume earth fixed observations unless specified otherwise
//                 ii = 0;
//                 // limit to obsread?/?????????
//                 while (feof(infile) == 0)
//                 {
//                     if (ii != 0)
//                         fgets(longstr2, 120, infile);
//                     else
//                         longstr2 = longstr1;  // data for the 1st point
//                     strr = longstr2[0];
//                     strr[1] = '\0';
// 
//                     if (!strr.Equals('#') == 0)
//                     {
//                         sscanf_s(longstr2, "%d", &currobsrec.obstype);
//                         switch (currobsrec.obstype)
//                         {
//                             case 0:
//                                 {
//                                     sscanf_s(longstr2, "%d %d %d %d %d %d %d %d %lf %lf ",
//                                         &currobsrec.obstype, &currobsrec.satnum, &currobsrec.sennum,
//                                         &currobsrec.year, &currobsrec.mon, &currobsrec.day,
//                                         &currobsrec.hr, &currobsrec.min, &currobsrec.sec,
//                                         &currobsrec.rng);
//                                 }
//                                 break;
//                             case 1:
//                                 {
//                                     sscanf_s(longstr2, "%d %d %d %d %d %d %d %d %lf %lf %lf ",
//                                         &currobsrec.obstype, &currobsrec.satnum, &currobsrec.sennum,
//                                         &currobsrec.year, &currobsrec.mon, &currobsrec.day,
//                                         &currobsrec.hr, &currobsrec.min, &currobsrec.sec,
//                                         &currobsrec.az, &currobsrec.el);
//                                     currobsrec.az = currobsrec.az / rad;
//                                     currobsrec.el = currobsrec.el / rad;
//                                 }
//                                 break;
//                             case 2:
//                                 {
//                                     sscanf_s(longstr2, "%d %d %d %d %d %d %d %d %lf %lf %lf %lf ",
//                                         &currobsrec.obstype, &currobsrec.satnum, &currobsrec.sennum,
//                                         &currobsrec.year, &currobsrec.mon, &currobsrec.day,
//                                         &currobsrec.hr, &currobsrec.min, &currobsrec.sec,
//                                         &currobsrec.rng, &currobsrec.az, &currobsrec.el);
//                                     currobsrec.az = currobsrec.az / rad;
//                                     currobsrec.el = currobsrec.el / rad;
//                                 }
//                                 break;
//                             case 3:
//                                 {
//                                     sscanf_s(longstr2, "%d %d %d %d %d %d %d %d %lf %lf %lf ",
//                                         &currobsrec.obstype, &currobsrec.satnum, &currobsrec.sennum,
//                                         &currobsrec.year, &currobsrec.mon, &currobsrec.day,
//                                         &currobsrec.hr, &currobsrec.min, &currobsrec.sec,
//                                         &currobsrec.trtasc, &currobsrec.tdecl);
//                                     currobsrec.trtasc = currobsrec.trtasc / rad;
//                                     currobsrec.tdecl = currobsrec.tdecl / rad;
//                                 }
//                                 break;
//                         } // case obstype
// 
//                         year = currobsrec.year;
//                         mon = currobsrec.mon;
//                         day = currobsrec.day;
//                         hr = currobsrec.hr;
//                         min = currobsrec.min;
//                         sec = currobsrec.sec;
//                         mfme = hr * 60 + min + sec / 60.0;
//                         // not 0.0 for utc
//                         MathTimeLibr.jday(year, mon, day, hr, min, sec, currobsrec.jd, currobsrec.jdf);
//                         if (ii == 1) // need the second one since it uses gibbs and thus the middle vector
//                         {
//                             MathTimeLibr.jday(year, mon, day, hr, min, sec, satrec.jdsatepoch, satrec.jdsatepochF);
//                             obsrecarr[0].dtmin = ((obsrecarr[0].jd +
//                                 (obsrecarr[0].hr * 60 + obsrecarr[0].min + obsrecarr[0].sec / 60.0) / 1440.0)
//                                 - satrec.jdsatepoch) * 1440.0;
//                         }
//                         // calculate the sensor each time because the obs might be from a different site
//                         SGP4DCLibr.getsensorparams(currobsrec.sennum, out currsenrec);
//                         AstroLibr.site(currsenrec.senlat, currsenrec.senlon, currsenrec.senalt,
//                             currobsrec.rsecef, currobsrec.vsecef);
// 
//                         dtmin = ((currobsrec.jd + currobsrec.jdf) - (jdepoch + jdepochf)) * 1440.0;  // satrec.jdepoch
//                         if (ii == 1)
//                             currobsrec.dtmin = 0.0;  // second state is the epoch
//                         else
//                             currobsrec.dtmin = dtmin;
//                         obsrecarr[ii] = currobsrec;
//                         ii = ii + 1;
//                     }   // if strr
//                 } // while ii through producing the observation vectors
// 
// 
//                 strBuild.AppendLine("over 30 obs %3i \n", ii - 1);
// 
//                 // setup nominal vector -----------------------------------------
//                 satrec.bstar = 0.0001;  // set a deafult
// 
//                 // need to get some kind of range estimate for these data types...
//                 if ((currobsrec.obstype == 1) || (currobsrec.obstype == 3))
//                 {
//                     strBuild.AppendLine("input range estimate in er (km?) \n");
//                     scanf_s("%lf", &rngest);
//                 }
// 
//                 zrng = 0.0;
//                 zaz = 0.0;
//                 zel = 0.0;
//                 zrtasc = 0.0;
//                 zdecl = 0.0;
//                 rsum[0] = 0.0;
//                 rsum[1] = 0.0;
//                 rsum[2] = 0.0;
//                 vsum[0] = 0.0;
//                 vsum[1] = 0.0;
//                 vsum[2] = 0.0;
//                 skipped = 0;
//                 keepit = 'y';
// 
//                 // try averaging the first 10 observations to get a value
//                 for (ktr = 1; ktr < 10; ktr++)     // 10
//                 {
//                     switch (currobsrec.obstype)
//                     {
//                         case 1:
//                             {
//                                 AstroLibr.rv_radec(ref r1, ref v1, MathTimeLib.efrom,
//                                     ref rngest, ref obsrecarr[ktr].az, ref obsrecarr[ktr].el, zrng, zaz, zel);
//                             }
//                             break;
//                         case 2:  // very simplistic iod here just gets 3 position vectors in ecef,
//                                  // obs are in topocentric (sez) coordinates (ecef) range az el obs
//                             {
//                                 AstroLibr.rv_razel(r1, v1, currsenrec.senlat, currsenrec.senlon, currsenrec.senalt, MathTimeLibr.eFrom,
//                                     obsrecarr[ktr].rng, obsrecarr[ktr].az, obsrecarr[ktr].el, zrng, zaz, zel);
//                                 AstroLibr.rv_razel(r2, v2, currsenrec.senlat, currsenrec.senlon, currsenrec.senalt, MathTimeLibr.eFrom,
//                                     obsrecarr[ktr + 1].rng, obsrecarr[ktr + 1].az, obsrecarr[ktr + 1].el, zrng, zaz, zel);
//                                 AstroLibr.rv_razel(r3, v3, currsenrec.senlat, currsenrec.senlon, currsenrec.senalt, MathTimeLibr.eFrom,
//                                     obsrecarr[ktr + 2].rng, obsrecarr[ktr + 2].az, obsrecarr[ktr + 2].el, zrng, zaz, zel);
//                             }
//                             break;
//                         case 3:
//                             {
//                                 // but the rs vector should be in eci here - see code...
//                                 AstroLibr.rv_tradec(r1, v1, obsrecarr[ktr].rsecef, MathTimeLibr.eFrom,
//                                     rngest, obsrecarr[ktr].trtasc, obsrecarr[ktr].tdecl, zrng, zrtasc, zdecl);
//                             }
//                             break;
//                     } // case obstype
// 
//                     // convert to eci  perform iod (in earth inertial, use teme as approx eci)
//                     jd = obsrecarr[ktr].jd;
//                     jdf = (obsrecarr[ktr].hr * 60.0 + obsrecarr[ktr].min + obsrecarr[ktr].sec / 60.0) * dayconv;
//                     EOPSPWLibr.findeopparam(jd, jdf, 's', eoparr, jdeopstart, dut1, dat, lod, xp, yp, ddpsi, ddeps,
//                         ddx, ddy, icrsx, y, s);
//                     MathTimeLibr.convtime(obsrecarr[ktr].year, obsrecarr[ktr].mon, obsrecarr[ktr].day, obsrecarr[ktr].hr, obsrecarr[ktr].min, obsrecarr[ktr].sec,
//                         timezone, dut1, dat, ut1, tut1, jdut1, jdut1f, utc, tai, tt, ttt, jdtt, jdftt, tcg, tdb, ttdb, jdtdb, jdtdbf, tcb);
//                     // itrf to teme
//                     //AstroLibr.teme_ecef(r1teme, v1teme, ateme, MathTimeLibr.eFrom, r1, v1, aitrf, ttt, xp, yp, jdut1, lod, eqeterms);
//                     // these will actually be regular eci
//                     AstroLibr.eci_ecef(r1, v1, aitrf, MathTimeLibr.eTo, r1teme, v1teme, ateme, e80,
//                         iau80rec, iau06rec, jdtt, jdftt, jdut1 + jdut1f, jdxysstart, lod, xp, yp, ddpsi, ddeps, ddx, ddy, trans);
// 
//                     jd = obsrecarr[ktr + 1].jd;
//                     jdf = (obsrecarr[ktr + 1].hr * 60.0 + obsrecarr[ktr + 1].min + obsrecarr[ktr + 1].sec / 60.0) * dayconv;
//                     EOPSPWLibr.findeopparam(jd, jdf, 's', eoparr, jdeopstart, dut1, dat, lod, xp, yp, ddpsi, ddeps,
//                         ddx, ddy, icrsx, y, s);
//                     MathTimeLibr.convtime(obsrecarr[ktr + 1].year, obsrecarr[ktr + 1].mon, obsrecarr[ktr + 1].day, obsrecarr[ktr + 1].hr, obsrecarr[ktr + 1].min, obsrecarr[ktr + 1].sec,
//                         timezone, dut1, dat, ut1, tut1, jdut1, jdut1f, utc, tai, tt, ttt, jdtt, jdftt, tcg, tdb, ttdb, jdtdb, jdtdbf, tcb);
//                     // itrf to teme
//                     //AstroLibr.teme_ecef(r2teme, v2teme, ateme, MathTimeLibr.eFrom, r2, v2, aitrf, ttt, xp, yp, jdut1, lod, eqeterms);
//                     AstroLibr.eci_ecef(r2, v2, aitrf, MathTimeLibr.eTo, r2teme, v2teme, ateme, e80,
//                         iau80rec, iau06rec, jdtt, jdftt, jdut1 + jdut1f, jdxysstart, lod, xp, yp, ddpsi, ddeps, ddx, ddy, trans);
// 
//                     jd = obsrecarr[ktr + 2].jd;
//                     jdf = (obsrecarr[ktr + 2].hr * 60.0 + obsrecarr[ktr + 2].min + obsrecarr[ktr + 2].sec / 60.0) * dayconv;
//                     EOPSPWLibr.findeopparam(jd, jdf, 's', eoparr, jdeopstart, dut1, dat, lod, xp, yp, ddpsi, ddeps,
//                         ddx, ddy, icrsx, y, s);
//                     MathTimeLibr.convtime(obsrecarr[ktr + 2].year, obsrecarr[ktr + 2].mon, obsrecarr[ktr + 2].day, obsrecarr[ktr + 2].hr, obsrecarr[ktr + 2].min, obsrecarr[ktr + 2].sec,
//                         timezone, dut1, dat, ut1, tut1, jdut1, jdut1f, utc, tai, tt, ttt, jdtt, jdftt, tcg, tdb, ttdb, jdtdb, jdtdbf, tcb);
//                     // itrf to teme
//                     //AstroLibr.teme_ecef(r3teme, v3teme, ateme, MathTimeLibr.eFrom, r3, v3, aitrf, ttt, xp, yp, jdut1, lod, eqeterms);
//                     AstroLibr.eci_ecef(r3, v3, aitrf, MathTimeLibr.eTo, r3teme, v3teme, ateme, e80,
//                         iau80rec, iau06rec, jdtt, jdftt, jdut1 + jdut1f, jdxysstart, lod, xp, yp, ddpsi, ddeps, ddx, ddy, trans);
// 
//                     // these are really eci here for now
//                     AstroLibr.gibbs(r1teme, r2teme, r3teme, v2teme, theta, theta1, copa, error);
//                     fstrBuild.AppendLine(outfile1, " gibbs  %11.5lf  %11.5lf  %11.5lf %11.7lf %11.7lf %11.7lf \n", r2teme[0], r2teme[1], r2teme[2], v2teme[0], v2teme[1], v2teme[2]);
//                     // most of these angles are just 0.72 deg apart
//                     // the times in herrgibbs need to be in increasing order...
//                     AstroLibr.herrgibbs(r1teme, r2teme, r3teme, obsrecarr[ktr - 1].jd + obsrecarr[ktr - 1].jdf,
//                         obsrecarr[ktr].jd + obsrecarr[ktr].jdf, obsrecarr[ktr + 1].jd + obsrecarr[ktr + 1].jdf,
//                         v2teme, theta, theta1, copa, error);   // this is earth fixed
//                     fstrBuild.AppendLine(outfile1, "hgibbs  %11.5lf  %11.5lf  %11.5lf %11.7lf %11.7lf %11.7lf ", r2[0], r2[1], r2[2], v2[0], v2[1], v2[2]);
//                     fstrBuild.AppendLine(outfile1, "hgibbs  %11.5lf  %11.5lf  %11.5lf %11.7lf %11.7lf %11.7lf \n", obsrecarr[ktr].dtmin - obsrecarr[ktr].dtmin,
//                         obsrecarr[ktr + 1].dtmin - obsrecarr[ktr].dtmin,
//                         obsrecarr[ktr + 2].dtmin - obsrecarr[ktr].dtmin,
//                         theta * rad, theta1 * rad, copa * rad);
// 
//                     // now sum up the vector in teme to get an initial guess
//                     // use two-body for simplicity and since short time span for iod
//                     // move vectors back to middle vector time which is one ahead
//                     AstroLibr.kepler(r2teme, v2teme, obsrecarr[ktr + 1].dtmin * 60.0, r1, v1);  // these will be all in eci/teme
// 
//                     //                  if (ktr > 0) // don't ask for the first one
//                     //                    {
//                     //                      strBuild.AppendLine("%11.7lf %11.7lf %11.7lf  keepit? y, n \n",r1[0], r1[1], r1[2] );
//                     //                      keepit = getchar();
//                     //                      fflush(stdin);
//                     //                    }
//                     //                  if (keepit == 'y')
//                     {
//                         rsum[0] = rsum[0] + r1[0];
//                         rsum[1] = rsum[1] + r1[1];
//                         rsum[2] = rsum[2] + r1[2];
//                         vsum[0] = vsum[0] + v1[0];
//                         vsum[1] = vsum[1] + v1[1];
//                         vsum[2] = vsum[2] + v1[2];
//                     }
//                     //                    else
//                     //                      skipped = skipped + 1;
// 
//                 }  // for ktr to sum up vectors
//                 fstrBuild.AppendLine(outfile1, "res eci  %11.5lf  %11.5lf  %11.5lf %11.7lf %11.7lf %11.7lf %3i %3i \n", rsum[0], rsum[1], rsum[2], vsum[0], vsum[1], vsum[2], ktr, skipped);
// 
//                 //              if (ktr > 1)
//                 //                  ktr = ktr - skipped;
// 
//                 rteme[0] = rsum[0] / ktr;
//                 rteme[1] = rsum[1] / ktr;
//                 rteme[2] = rsum[2] / ktr;
//                 vteme[0] = vsum[0] / ktr;
//                 vteme[1] = vsum[1] / ktr;
//                 vteme[2] = vsum[2] / ktr;
// 
//                 // set it outright for testing =================================================
//                 rteme[0] = 5975.2904;
//                 rteme[1] = 2568.6400;
//                 rteme[2] = 3120.5845;
//                 vteme[0] = 3.983846;
//                 vteme[1] = -2.071159;
//                 vteme[2] = -5.917095;
// 
//                 for (i = 0; i < 6; i++)
//                 {
//                     if (i < 3)
//                         xnom[i] = rteme[i];
//                     else
//                         xnom[i] = vteme[i - 3];
//                 }
// 
//                 if (statesize > 6)
//                     xnom[6] = 0.0001; // set initial guess as it changes in state2satrec
//                 else
//                     satrec.bstar = 0.0001;
//                 SGP4DC::state2satrec(xnom, jdepoch, jdepochf, 'v', statesize, MathTimeLibr.eTo, satrec); // force this method to assign values
//                 if (obsrecarr[1].year >= 2000)
//                     satrec.epochyr = obsrecarr[1].year - 2000;
//                 else
//                     satrec.epochyr = obsrecarr[1].year - 1900;
//                 MathTimeLibr.finddays(obsrecarr[1].year, obsrecarr[1].mon, obsrecarr[1].day, obsrecarr[1].hr,
//                     obsrecarr[1].min, obsrecarr[1].sec, satrec.epochdays);
//                 //               satrec.jdsatepoch = satrec.jdsatepoch;  // set earlier
// 
//                 SGP4Funcs::sgp4init(whichconst, opsmode, satrec.satnum, satrec.jdsatepoch - 2433281.5, satrec.bstar,
//                     satrec.ndot, satrec.nddot, satrec.ecco, satrec.argpo, satrec.inclo, satrec.mo, satrec.no_kozai,
//                     satrec.nodeo, satrec);
// 
// 
//                 strBuild.AppendLine("rteme  %11.5lf  %11.5lf  %11.5lf %11.7lf %11.7lf %11.7lf \n", rteme[0], rteme[1], rteme[2], vteme[0], vteme[1], vteme[2]);
//                 strBuild.AppendLine("coe  %11.5lf  %11.5lf  %11.5lf %11.7lf %11.7lf %11.7lf %11.7lf \n", satrec.a * 6378.135, satrec.ecco, satrec.inclo * rad,
//                     satrec.nodeo * rad, satrec.argpo * rad, satrec.mo * rad, satrec.no_kozai);
//                 fstrBuild.AppendLine(outfile1, "rteme  %11.5lf  %11.5lf  %11.5lf %11.7lf %11.7lf %11.7lf \n", rteme[0], rteme[1], rteme[2], vteme[0], vteme[1], vteme[2]);
//                 fstrBuild.AppendLine(outfile1, "coe  %11.5lf  %11.5lf  %11.5lf %11.7lf %11.7lf %11.7lf %11.7lf \n", satrec.a * 6378.135, satrec.ecco, satrec.inclo * rad,
//                     satrec.nodeo * rad, satrec.argpo * rad, satrec.mo * rad, satrec.no_kozai);
// 
//                 // nominal vector in teme
//                 //              satrec.a = 7203.20025180/6378.135;  //er
//                 //              satrec.ecco = 0.0015278327;
//                 //              satrec.inclo = 114.954226 / rad;
//                 //              satrec.nodeo = 190.308229 / rad;
//                 //              satrec.argpo = 345.4910919 / rad;
//                 //              satrec.mo = 165.996981 / rad;
//                 // answer??
//                 satrec.a = 7214.3395180 / 6378.135;  //er
//                 satrec.ecco = 0.0006009327;
//                 satrec.inclo = 114.962396 / rad;
//                 satrec.nodeo = 190.297730 / rad;
//                 satrec.argpo = 123.8855619 / rad;
//                 satrec.mo = 123.828381 / rad;
//                 satrec.no_kozai = (1.0 / 13.446839) * sqrt(1.0 / (satrec.a * satrec.a * satrec.a));  // rad / min
// 
//                 SGP4Funcs::sgp4init(whichconst, opsmode, satrec.satnum, satrec.jdsatepoch - 2433281.5, satrec.bstar,
//                     satrec.ndot, satrec.nddot, satrec.ecco, satrec.argpo, satrec.inclo, satrec.mo, satrec.no_kozai,
//                     satrec.nodeo, satrec);
//                 // call the propagator to get the initial state vector value
//                 SGP4Funcs::sgp4(satrec, 0.0, rnomTEME, vnomTEME);
// 
//                 SGP4DC::state2satrec(xnom, jdepoch, jdepochf, statetype, statesize, MathTimeLibr.eFrom, satrec);
//                 //rnomTEME[0] = 6157.4121;
//                 //rnomTEME[1] = 2455.16193;
//                 //rnomTEME[2] = 2823.89107;
//                 //vnomTEME[0] = 3.667402;
//                 //vnomTEME[1] = -2.204994;
//                 //vnomTEME[2] = -6.0647111;
//                 strBuild.AppendLine("running new satellite \n");
// 
//                 // --------- now determine the orbit for the satellite ----------
//                 SGP4DC::leastsquares(percentchg, deltaamtchg, rmsepsilon, whichconst, interp, jdeopstart, satrec, typeans,
//                     statetype, proptype, firstob, obsread, statesize, obsrecarr, loops,
//                     xnom, jdepoch, jdepochf, derivType, cdam, cram, jplopt, jpldearr, jdjpldestart, iau80rec, iau06rec, eoparr, x, dx, atwai, atwa, atwb, outfile, outfile1);
// 
//                 SGP4Funcs::sgp4(satrec, 0.0, rnomTEME, vnomTEME);
//                 rpos[0] = obsrecarr[0].x;
//                 rpos[1] = obsrecarr[0].y;
//                 rpos[2] = obsrecarr[0].z;
//                 vpos[0] = obsrecarr[0].xdot;
//                 vpos[1] = obsrecarr[0].ydot;
//                 vpos[2] = obsrecarr[0].zdot;
//                 for (j = 0; j < 3; j++)
//                 {
//                     dr[j] = rnomTEME[j] - rpos[j];
//                     dv[j] = vnomTEME[j] - vpos[j];
//                 }
//                 fstrBuild.AppendLine(outfile, " difference %11lf  %11lf  %11lf %11lf \n", dr[0], dr[1], dr[2], MathTimeLibr.mag(dr));
//                 fstrBuild.AppendLine(outfile1, " difference %11lf  %11lf  %11lf %11lf \n", dr[0], dr[1], dr[2], MathTimeLibr.mag(dr));
//                 strBuild.AppendLine(" difference %11lf  %11lf  %11lf %11lf \n", dr[0], dr[1], dr[2], MathTimeLibr.mag(dr));
// 
//                 fclose(infile);
//                 fclose(outfile);
//                 fclose(outfile1);
//             }
//         }  // processrazel
// 
// 
//         //// trim from start (in place)
//         //static inline void ltrim(std::string &s)
//         //{
//         //    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch)
// 
//         //    {
//         //        return !std::isspace(ch);
//         //    }));
//         //}
// 
//         //void removeWS(char* c)
//         //{
//         //    char* s = c;
//         //    while (s++, *(c++))
//         //	{
//         //    if (isspace(*s)) ++s;
//         //    *s = *c;
//         //}
//         //}
// 
//         //void whitespace(char str[])
//         //{
//         //    int i;
//         //    int j;
//         //    for (i = 0; str[i] != 0; i++)
//         //    {
//         //        if (isspace(str[i]))
//         //        {
//         //            for (j = i; str[j] != 0; j++)
//         //            {
//         //                str[j] = str[j + 1];
//         //            }
//         //        }
//         //    }
//         //}
// 
// 

        /* ----------------------------------------------------------------------------------------
        *  processxyz
        *
        *  New C# port of the "read a .e file and run LS OD on it" path. Reads a standard AGI/STK
        *  ephemeris (.e) file's EphemerisTimePosVel block (time-from-epoch-seconds, x,y,z,vx,vy,vz
        *  in km/km-s), builds obstype-4 (position/velocity) observations from it, seeds the nominal
        *  state from the first observation (as the original design intended), and runs leastsquares.
        *
        *  ASSUMPTION I COULDN'T VERIFY: I don't have one of your actual .e files to check the exact
        *  header keywords/spacing against, so this uses the standard documented STK .e layout
        *  (ScenarioEpoch / CoordinateSystem / EphemerisTimePosVel / END Ephemeris keywords, one
        *  header keyword per line). It also assumes CoordinateSystem is an inertial frame close
        *  enough to this codebase's ECI (J2000/ICRF, not Fixed/ECF or TEME) and does NOT do any
        *  frame conversion on the read values - if your files use a different coordinate system
        *  or header layout (e.g. what ConvertEphem writes), send me a sample and I'll fix this
        *  against the real thing rather than my best guess.
        * ---------------------------------------------------------------------------------------- */
        private void processxyz(char typeans, char proptype, int obsskip, int obsread, string outputFolder)
        {
            StringBuilder strBuild = new StringBuilder();
            strBuild.AppendLine("=============== processxyz (.e file position/velocity DC) ===============");

            string outputFolderSafe = outputFolder == null ? "" : outputFolder;
            string testsDir = GetTestsDir();
            string[] candidateDirs = new string[] { testsDir, "", outputFolderSafe };
            string fname = null;
            for (int c = 0; c < candidateDirs.Length && fname == null; c++)
            {
                string dir = string.IsNullOrEmpty(candidateDirs[c]) ? "." : candidateDirs[c];
                if (!Directory.Exists(dir))
                    continue;
                string[] hits = Directory.GetFiles(dir, "*.e");
                if (hits.Length > 0)
                    fname = hits[0];
            }
            if (fname == null)
            {
                MessageBox.Show("Could not find a *.e ephemeris file in the tests folder (" + testsDir +
                    "), next to the exe, or in the Output file location folder.");
                return;
            }

            string[] allLines = File.ReadAllLines(fname);

            // ---- pull the scenario epoch and coordinate system out of the header ----
            DateTime? scenarioEpoch = null;
            string coordSystem = null;
            int dataStart = -1;
            for (int li = 0; li < allLines.Length; li++)
            {
                string line = allLines[li].Trim();
                if (line.StartsWith("ScenarioEpoch", StringComparison.OrdinalIgnoreCase))
                {
                    string dateStr = line.Substring("ScenarioEpoch".Length).Trim();
                    DateTime parsed;
                    // STK's epoch text form has more fractional-second digits (commonly 9) than .NET's
                    // DateTime format strings support (max 7, "fffffff") - so exact-format parsing isn't
                    // reliable here. General TryParse handles arbitrary fractional-second precision fine.
                    if (DateTime.TryParse(dateStr, CultureInfo.InvariantCulture,
                            DateTimeStyles.AllowWhiteSpaces, out parsed))
                        scenarioEpoch = parsed;
                }
                else if (line.StartsWith("CoordinateSystem", StringComparison.OrdinalIgnoreCase))
                    coordSystem = line.Substring("CoordinateSystem".Length).Trim();
                else if (line.StartsWith("EphemerisTimePosVel", StringComparison.OrdinalIgnoreCase))
                {
                    dataStart = li + 1;
                    break;
                }
            }
            if (scenarioEpoch == null || dataStart < 0)
            {
                MessageBox.Show("Couldn't find the expected \"ScenarioEpoch\" and \"EphemerisTimePosVel\" " +
                    "header lines in " + fname + " - this parser assumes the standard STK .e layout. " +
                    "Send me a sample of the actual file and I'll fix the parsing against it.");
                return;
            }
            if (coordSystem != null && coordSystem.IndexOf("J2000", StringComparison.OrdinalIgnoreCase) < 0 &&
                coordSystem.IndexOf("ICRF", StringComparison.OrdinalIgnoreCase) < 0)
            {
                strBuild.AppendLine("WARNING: CoordinateSystem is \"" + coordSystem + "\", not J2000/ICRF - " +
                    "this code reads the vectors as-is with no frame conversion, so results may be off if " +
                    "that's actually a rotating (Fixed) or TEME frame.");
            }

            double epochJd, epochJdf;
            MathTimeLibr.jday(scenarioEpoch.Value.Year, scenarioEpoch.Value.Month, scenarioEpoch.Value.Day,
                scenarioEpoch.Value.Hour, scenarioEpoch.Value.Minute,
                scenarioEpoch.Value.Second + scenarioEpoch.Value.Millisecond / 1000.0, out epochJd, out epochJdf);

            List<double[]> rows = new List<double[]>();  // [tsec, x,y,z,vx,vy,vz] - converted to km, km/s below
            for (int li = dataStart; li < allLines.Length; li++)
            {
                string line = allLines[li].Trim();
                if (line.Length == 0) continue;
                if (line.StartsWith("END", StringComparison.OrdinalIgnoreCase)) break;
                string[] f = line.Split((char[])null, StringSplitOptions.RemoveEmptyEntries);
                if (f.Length < 7) break;  // not a data row (e.g. a following keyword) - stop
                double[] row = new double[7];
                bool ok = true;
                for (int k = 0; k < 7; k++)
                    ok &= double.TryParse(f[k], NumberStyles.Float, CultureInfo.InvariantCulture, out row[k]);
                if (!ok) break;
                // STK .e files are in meters/meters-per-second; this codebase works in km/km-s throughout
                row[1] /= 1000.0; row[2] /= 1000.0; row[3] /= 1000.0;
                row[4] /= 1000.0; row[5] /= 1000.0; row[6] /= 1000.0;
                rows.Add(row);
            }
            if (rows.Count < 3)
            {
                MessageBox.Show("Only found " + rows.Count + " data rows under EphemerisTimePosVel in " + fname +
                    " - need at least 3.");
                return;
            }

            List<double[]> picked = new List<double[]>();
            for (int idx = 0; idx < rows.Count && picked.Count < obsread; idx += obsskip)
                picked.Add(rows[idx]);
            int numobs = picked.Count;

            SGP4DCLib.obsrec[] obsrecarr = new SGP4DCLib.obsrec[numobs + 1];
            for (int j = 1; j <= numobs; j++)
            {
                double[] row = picked[j - 1];
                // BUG FIX: this used to recompute jd/jdf per row via a fresh jday() call that dropped
                // scenarioEpoch's fractional seconds (DateTime.Millisecond was never included, unlike
                // the epochJd/epochJdf computed above - which were then never actually used). Reusing
                // that one correct epoch value and just adding each row's offset avoids both the
                // dropped-milliseconds bug and the redundant/wasteful re-parsing every iteration.
                double jd = epochJd;
                double jdf = epochJdf + row[0] / 86400.0;  // add the row's time-from-epoch offset (sec -> days)

                SGP4DCLib.obsrec rec = new SGP4DCLib.obsrec();
                rec.sennum = 1;
                rec.obstype = 4;
                rec.jd = jd; rec.jdf = jdf;
                rec.x = row[1]; rec.y = row[2]; rec.z = row[3];
                rec.xdot = row[4]; rec.ydot = row[5]; rec.zdot = row[6];
                rec.bstar = 0.0;
                obsrecarr[j] = rec;
            }
            obsrecarr[0] = obsrecarr[1];

            // seed the nominal from the first observation (matches the original design's intent)
            double[] xnom = new double[7];
            xnom[0] = obsrecarr[1].x; xnom[1] = obsrecarr[1].y; xnom[2] = obsrecarr[1].z;
            xnom[3] = obsrecarr[1].xdot; xnom[4] = obsrecarr[1].ydot; xnom[5] = obsrecarr[1].zdot;

            SGP4Lib.elsetrec satrec = new SGP4Lib.elsetrec();
            int loops;
            double[,] x, dx, atwai, atwa, atwb;
            SGP4DCLibr.leastsquares(
                0.001, 0.0000001, 0.0002,
                SGP4Lib.gravconsttype.wgs84,
                's', 0.0,
                ref satrec,
                typeans, 'v', proptype,
                1, numobs, 6,
                ref obsrecarr, out loops,
                ref xnom, obsrecarr[1].jd, obsrecarr[1].jdf,
                1, 0.044, 0.02, 'a',
                new EOPSPWLib.iau80Class(), new EOPSPWLib.iau06Class(),
                ref eoparrUnused,
                out x, out dx, out atwai, out atwa, out atwb);

            if (SGP4DCLibr.strbuild.Length > 0)
            {
                strBuild.AppendLine("--- leastsquares diagnostics ---");
                strBuild.Append(SGP4DCLibr.strbuild.ToString());
                SGP4DCLibr.strbuild.Clear();
            }

            double[] rFinal = new double[] { xnom[0], xnom[1], xnom[2] };
            double[] vFinal = new double[] { xnom[3], xnom[4], xnom[5] };

            // BUG FIX: this used to subtract rFinal/vFinal (the fitted state AT THE FIRST OBSERVATION'S
            // EPOCH) directly from the last observation's raw position - but the last observation is
            // ~(numobs-1) sample intervals LATER in time, so that "diff" was really just measuring the
            // satellite's actual orbital displacement over that whole span (tens of thousands of km for
            // an hour-long LEO arc), not fit quality. Propagate the fitted epoch state forward to the
            // last observation's actual time via the same two-body Kepler propagator the DC itself uses
            // internally, then compare - that's an apples-to-apples residual check.
            double dtsecLast = (obsrecarr[numobs].jd + obsrecarr[numobs].jdf - obsrecarr[1].jd - obsrecarr[1].jdf) * 86400.0;
            double[] rAtLast, vAtLast;
            AstroLibr.kepler(rFinal, vFinal, dtsecLast, out rAtLast, out vAtLast, out string keplerOutTextLast);
            double dxp = rAtLast[0] - obsrecarr[numobs].x;
            double dyp = rAtLast[1] - obsrecarr[numobs].y;
            double dzp = rAtLast[2] - obsrecarr[numobs].z;

            strBuild.AppendLine("obs file    " + fname);
            strBuild.AppendLine("coord sys   " + (coordSystem ?? "(not found in header)"));
            strBuild.AppendLine("numobs used " + numobs + "   loops " + loops);
            strBuild.AppendLine(string.Format(CultureInfo.InvariantCulture,
                "nominal r  {0,16:F8} {1,16:F8} {2,16:F8} km", rFinal[0], rFinal[1], rFinal[2]));
            strBuild.AppendLine(string.Format(CultureInfo.InvariantCulture,
                "nominal v  {0,16:F8} {1,16:F8} {2,16:F8} km/s", vFinal[0], vFinal[1], vFinal[2]));
            strBuild.AppendLine(string.Format(CultureInfo.InvariantCulture,
                "diff vs last obs read (propagated to its time)  {0,12:F8} {1,12:F8} {2,12:F8} km", dxp, dyp, dzp));

            string outText = strBuild.ToString();
            try
            {
                Directory.CreateDirectory(testsDir);
                string outPath = Path.Combine(testsDir, "testdc_xyz.out");
                File.WriteAllText(outPath, outText);
                MessageBox.Show(outText + "\n\nFull results written to:\n" + outPath);
            }
            catch (Exception ex)
            {
                MessageBox.Show(outText + "\n\n(could not write output file: " + ex.Message + ")");
            }
        }

        // ---- UNCONVERTED LEGACY CODE (processxyz) commented out below pending C# port ----
//         /*  --------------------------------------------------------------------------    
//         // this procedure reads a .e file and performs LS OD on it
//         // there are options to read a span of data, at given timesteps
//         // the process use numerical partial derivatives for the A matrix
//         // the nominal vector is formed at the timeof the first measurement
//         //  -------------------------------------------------------------------------- */
// 
//         public void processxyz
//         (
//                     char batchmode, char proptype, int derivType, double cdam, double cram,
//                     EOPSPWLib.iau80Class iau80arr,
//                     EOPSPWLib.iau06Class iau06arr,
//                     double[] eoparr,
//                     double jdeopstart, double jdeopstartf,
//                     double[] spwarr,
//                     double jdspwstart, double jdspwstartf,
//                     double[] jpldearr,
//                 char jplopt,
//                     double jdjpldestart,
//                     double jdepoch, double jdepochf, char interp,
//                     int firstob, int obsskip, int obsread,
//                     double[] obsarr,
//                     out double[] xnom
//         )
//         {
//             FILE* outfile, *outfile1;
//             SGP4Lib.elsetrec satrec;
//             long search_satno = 0;
//             char typerun;
//             int loops, eqeterms, j;
//             double ttt;
//             double[]<double[]< double >> trans = double[] < double[] < double >> (3, double[] < double > (3));
//             char statetype, typeans, opsmode;
//             int statesize;
//             double percentchg, deltaamtchg, rmsepsilon;
//             double magcov;
//             double rpos[3], vpos[3], rnomTEME[3], vnomTEME[3], dr[3], dv[3], magdr1, magdr, reci[3], veci[3], aeci[3],
//                 ateme[3], recin[3], vecin[3];
// 
//             SGP4Lib.gravconsttype whichconst = SGP4Lib.gravconsttype.wgs72; // type of run for sgp4
//             typerun = 'c';  // c catalog, v verification, m manual
//             opsmode = 'a';
//             typeans = 'b';  // mat inv or svd
//             eqeterms = 2;
//             interp = 's';
// 
//             if (batchmode == 'n')
//             {
//                 statetype = 'v';   // v vectors, t tle vars, e equinoctial
//                                    //strBuild.AppendLine("input type of elements tle or equin \n\n");
//                                    //statetype = getchar();
//                                    //fflush(stdin);
//                 percentchg = 0.001;  // percent to change components for finite differencing
//                                      //strBuild.AppendLine("input percentchg for finite differencing (.001) \n");
//                                      //scanf_s("%lf", &percentchg, 5 * sizeof(char));
//                 deltaamtchg = 0.0000001;  // delta amt to change components for finite differencing
//                                           //strBuild.AppendLine("input deltaamtchg to chk in finite differencing (.0000001) \n");
//                                           //scanf_s("%lf", &deltaamtchg, 10 * sizeof(char));
//                 rmsepsilon = 0.0002;  // percent to change components for finite differencing
//                                       //strBuild.AppendLine("input rmsepsilon to quit (.02) \n");
//                                       //scanf_s("%lf", &rmsepsilon, 5 * sizeof(char));
//                 statesize = 6;        // solve for bstar or not
//                                       //strBuild.AppendLine("input statesize (6 or 7) \n");
//                                       //scanf_s("%i", &statesize, sizeof(char));
//             }
// 
//             // set these variables up with initial dimensions
//             // these must be set here
//             // function matrices can be resized
//             double[,] x = new double[statesize, 1];
//             double[,] dx = new double[statesize, 1];
//             double[,] atwai = new double[statesize, statesize];
//             double[,] atwa = new double[statesize, statesize];
//             double[,] atwb = new double[statesize, 1];
// 
//             strBuild.AppendLine("state type  %c \n", statetype);
//             strBuild.AppendLine("percentchg  %lf pct to change components for finite diff \n", percentchg);
//             strBuild.AppendLine("deltaamtchg %11.8lf static change amount for each case \n", deltaamtchg);
//             strBuild.AppendLine("rmsepsilon  %lf \n", rmsepsilon);
//             strBuild.AppendLine("state size  %i \n", statesize);
// 
// 
//             // open files for reading and writing results
//             int err = fopen_s(&outfile, "testdc.out", "w");
//             err = fopen_s(&outfile1, "testdc1.out", "w");
// 
//             //	fstrBuild.AppendLine(outfile, "typeans %c \n", typeans);
//             fstrBuild.AppendLine(outfile, "state type  %c \n", statetype);
//             fstrBuild.AppendLine(outfile, "percentchg  %lf \n", percentchg);
//             fstrBuild.AppendLine(outfile, "deltaamtchg %11.8lf \n", deltaamtchg);
//             fstrBuild.AppendLine(outfile, "rmsepsilon  %lf \n", rmsepsilon);
//             fstrBuild.AppendLine(outfile, "state size  %i \n", statesize);
// 
// 
//             // setup nominal vector -----------------------------------------
//             eOpt opt = e80;
//             if (proptype == 's')
//             {
//                 if (proptype == 's')
//                 {
//                     recin[0] = obsarr[0].x;
//                     recin[1] = obsarr[0].y;
//                     recin[2] = obsarr[0].z;
//                     vecin[0] = obsarr[0].xdot;
//                     vecin[1] = obsarr[0].ydot;
//                     vecin[2] = obsarr[0].zdot;
// 
//                     // convert eci to teme for sgp4 tle creation!
//                     ttt = (jdepoch + jdepochf - 2451545.0) / 36525.0;
//                     AstroLibr.teme_eci(rnomTEME, vnomTEME, ateme, MathTimeLibr.eFrom, recin, vecin, aeci, iau80rec, opt, ttt, 0.0, 0.0);
// 
//                     xnom[0] = rnomTEME[0];
//                     xnom[1] = rnomTEME[1];
//                     xnom[2] = rnomTEME[2];
//                     xnom[3] = vnomTEME[0];
//                     xnom[4] = vnomTEME[1];
//                     xnom[5] = vnomTEME[2];
//                     SGP4DC::state2satrec(xnom, jdepoch, jdepochf, statetype, statesize, MathTimeLibr.eTo, satrec);
//                 }
//                 satrec.bstar = 0.0001;  // set a deafult for opsmode
//                 SGP4Funcs::sgp4init(whichconst, 'a', satrec.satnum, satrec.jdsatepoch - 2433281.5, satrec.bstar,
//                     satrec.ndot, satrec.nddot, satrec.ecco, satrec.argpo, satrec.inclo, satrec.mo, satrec.no_kozai,
//                     satrec.nodeo, satrec);
//                 // call the propagator to get the initial state vector value
//                 SGP4Funcs::sgp4(satrec, 0.0, rnomTEME, vnomTEME);
//                 ttt = (satrec.jdsatepoch + satrec.jdsatepochF - 2451545.0) / 36525.0;
//                 AstroLibr.teme_eci(rnomTEME, vnomTEME, ateme, MathTimeLibr.eTo, recin, vecin, aeci, iau80rec, opt, ttt, 0.0, 0.0);
//                 for (j = 0; j < 3; j++)
//                 {
//                     xnom[j] = recin[j];
//                     xnom[j + 3] = vecin[j];
//                 }
//             }
//             else
//             {
//                 for (j = 0; j < 3; j++)
//                 {
//                     recin[j] = xnom[j];
//                     vecin[j] = xnom[j + 3];
//                 }
//             }
// 
//             strBuild.AppendLine(" init nom vector %11lf  %11lf  %11lf  %11lf  %11lf  %11lf \n", recin[0], recin[1], recin[2], vecin[0], vecin[1], vecin[2]);
// 
//             // print initial difference
//             reci[0] = obsarr[0].x;
//             reci[1] = obsarr[0].y;
//             reci[2] = obsarr[0].z;
//             veci[0] = obsarr[0].xdot;
//             veci[1] = obsarr[0].ydot;
//             veci[2] = obsarr[0].zdot;
//             for (j = 0; j < 3; j++)
//             {
//                 dr[j] = recin[j] - reci[j];
//                 dv[j] = vecin[j] - veci[j];
//             }
//             fstrBuild.AppendLine(outfile, " init diff %11lf  %11lf  %11lf %11lf \nNote the this assumes the .e file starts at satellite epoch\n", dr[0], dr[1], dr[2], MathTimeLibr.mag(dr));
//             magdr1 = 1000.0 * MathTimeLibr.mag(dr); // initial difference in m
//             strBuild.AppendLine(" init diff %11lf  %11lf  %11lf %11lf \nNote the this assumes the .e file starts at satellite epoch\n", dr[0], dr[1], dr[2], MathTimeLibr.mag(dr));
//             magdr1 = 1000.0 * MathTimeLibr.mag(dr); // initial difference in m
// 
//             strBuild.AppendLine("running new satellite \n");
// 
//             // --------- now determine the orbit for the satellite ----------
//             // note xnom is changed in here - it should come out in x though?....
//             SGP4DC::leastsquares(percentchg, deltaamtchg, rmsepsilon, whichconst, interp, jdeopstart, satrec, typeans,
//                 statetype, proptype, firstob, obsread, statesize, obsarr, loops,
//                 xnom, jdepoch, jdepochf, derivType, cdam, cram, jplopt, jpldearr, jdjpldestart, iau80rec, iau06rec, eoparr, x, dx, atwai, atwa, atwb, outfile, outfile1);
// 
//             if (proptype == 's')
//             {
//                 SGP4Funcs::sgp4(satrec, 0.0, rnomTEME, vnomTEME);
//                 AstroLibr.teme_eci(rnomTEME, vnomTEME, ateme, MathTimeLibr.eTo, recin, vecin, aeci, iau80rec, opt, ttt, 0.0, 0.0);
//             }
//             else
//             {
//                 for (j = 0; j < 3; j++)
//                 {
//                     recin[j] = xnom[j];
//                     vecin[j] = xnom[j + 3];
//                 }
//             }
// 
//             rpos[0] = obsarr[0].x;
//             rpos[1] = obsarr[0].y;
//             rpos[2] = obsarr[0].z;
//             vpos[0] = obsarr[0].xdot;
//             vpos[1] = obsarr[0].ydot;
//             vpos[2] = obsarr[0].zdot;
//             for (j = 0; j < 3; j++)
//             {
//                 dr[j] = recin[j] - rpos[j];
//                 dv[j] = vecin[j] - vpos[j];
//             }
//             magdr = MathTimeLibr.mag(dr) * 1000.0;
//             //               magcov =  sqrt( xnom[0] / satrec.ecco * atwai[1][0] +
//             //                              -xnom[1] / (satrec.ecco * satrec.ecco) * atwai[4][0] +
//             //                               xnom[1] / (satrec.ecco * satrec.ecco) * atwai[5][0]);
// 
//             magcov = sqrt(atwai[2][2]) * 6378135.0;  // in m
// 
//             fstrBuild.AppendLine(outfile, " difference %11lf  %11lf  %11lf %11lf %10.2lf \n", dr[0], dr[1], dr[2], MathTimeLibr.mag(dr), magcov);
//             fstrBuild.AppendLine(outfile1, " difference %11lf  %11lf  %11lf km %12lf m %10.2lf \n\n", dr[0], dr[1], dr[2], 1000.0 * MathTimeLibr.mag(dr), magcov);
//             if (proptype == 's')
//                 strBuild.AppendLine("%4i %8.1lf %8.4lf m a %10.5lf e %10.7lf i %8.3lf %11.6lf m\n", loops, magdr1, magdr, satrec.a, satrec.ecco, satrec.inclo * 57.2957795,
//                     magcov);
// 
//             // ------------- now take the converged TLE and propagate (and ----------------
//             // ---------- convert coordinate system) over ephemeris and compare -----------
//             //for (j = 0; j < obsread; j++)
//             //{
//             //	if (proptype == 's')
//             //		SGP4Funcs::sgp4(satrec, obsrecarr[j].dtmin, rteme, vteme);
//             //	dr[0] = 1000.0 * (obsrecarr[j].x - reci[0]);
//             //	dr[1] = 1000.0 * (obsrecarr[j].y - reci[1]);
//             //	dr[2] = 1000.0 * (obsrecarr[j].z - reci[2]);
//             //	dv[0] = 1000.0 * (obsrecarr[j].xdot - veci[0]);
//             //	dv[1] = 1000.0 * (obsrecarr[j].ydot - veci[1]);
//             //	dv[2] = 1000.0 * (obsrecarr[j].zdot - veci[2]);
//             //	fstrBuild.AppendLine(outfile, "diff %11.3lf  %11.4lf %11.4lf %11.4lf  %11.3lf m \n", obsrecarr[j].dtmin, dr[0], dr[1], dr[2], MathTimeLibr.mag(dr));
//             //}
//             fstrBuild.AppendLine(outfile, " init diff %11lf  %11lf  %11lf %11lf \n", dr[0], dr[1], dr[2], MathTimeLibr.mag(dr));
//             magdr1 = 1000.0 * MathTimeLibr.mag(dr); // initial difference in m
// 
//             if (proptype == 's')
//                 SGP4DC::printtle(satrec);
// 
//             fclose(outfile);
//             fclose(outfile1);
// 
//         }  // processxyz



        // ---- UNCONVERTED LEGACY CODE (processtle) commented out below pending C# port ----
// /*       ---------------------------------------------------------------      */
// void processtle
// (
//             char batchmode, char proptype, int derivType, double cdam, double cram,
//             EOPSPWLib.iau80Class iau80arr,
//             EOPSPWLib.iau06Class iau06arr,
//             double[] eoparr,
//             double jdeopstart, double jdeopstartf,
//             double[] spwarr,
//             double jdspwstart, double jdspwstartf,
//             double[] jpldearr,
//         char jplopt,
//             double jdjpldestart,
//             double jdepoch, double jdepochf, char interp,
//             int firstob, int obsskip, int obsread,
//             double[] obsarr,
//             out double[] xnom
// )
//         {
//     FILE* infile, *outfile, *outfile1, *outfile2;
//     double[]<obsrec> obsrecarr;
//     SGP4Lib.elsetrec satrec;
//     long search_satno = 0;
//     char longstr1[130], longstr2[130];
//     char str[2];
//     char typerun;
//             double startmfe, stopmfe, deltamin, magdr;
//     int ll1m, l1m, l10m, l50m, l1km, l10km;
//     double totalnum, periodmin, avgiter, magcov;
//     int i, ii, j, loops;
//     double rnomTEME[3], vnomTEME[3], rpos[3], vpos[3], dr[3], dv[3], dtmin, jd, jdf, magdr1;
//     SGP4DCLib.obsrec currobsrec;
//     double rad, apalt, pralt;
//     char statetype, typeans, opsmode;
//     int statesize, lastob;
//     char fname[160];
//     double percentchg, deltaamtchg, rmsepsilon;
//     obsrecarr.resize(5000);
// 
//     interp = 'l'; // set eop interpolation to linear
// 
//     if (batchmode == 'n') // needs a default value here in batch mode
//         statesize = 7;
//             // set these variables up with intial dimensions
//             // these must be set here
//             // function matrices can be resized
//             double[,] x = new double[statesize, 1];
//             double[,] dx = new double[statesize, 1];
//             double[,] atwai = new double[statesize, statesize];
//             double[,] atwa = new double[statesize, statesize];
//             double[,] atwb = new double[statesize, 1];
// 
//             rad = 180.0 / pi;
//             SGP4Lib.gravconsttype whichconst = SGP4Lib.gravconsttype.wgs72; // type of run for sgp4
//             typerun = 'c';  // c catalog, v verification, m manual
//     if (batchmode == 'n')
//     {
//         statetype = 't';   // v vectors, t tle vars, e equinoctial
//         strBuild.AppendLine("input type of elements tle or equin \n\n");
// 
//         statetype = getchar();
//         fflush(stdin);
// 
//         strBuild.AppendLine("input TLE filename: \n");
//         scanf_s("%s", fname, 85 * sizeof(char));
// 
//         percentchg = 0.001;  // percent to change components for finite differencing
//                              //strBuild.AppendLine("input percentchg for finite differencing (.001) \n");
//                              //scanf_s("%lf", &percentchg, 5 * sizeof(char));
// 
//         deltaamtchg = 0.0000001;  // delta amt to change components for finite differencing
//                                   //strBuild.AppendLine("input deltaamtchg to chk in finite differencing (.0000001) \n");
//                                   //scanf_s("%lf", &deltaamtchg, 10 * sizeof(char));
// 
//         rmsepsilon = 0.0002;  // percent to change components for finite differencing
//                               //strBuild.AppendLine("input rmsepsilon to quit (.02) \n");
//                               //scanf_s("%lf", &rmsepsilon, 5 * sizeof(char));
// 
//         lastob = 72;        // how many obs per rev to consider
//         strBuild.AppendLine("input lastob to quit (72) \n");  // this will be 72 points per rev, but calculate later for 2 periods
//         scanf_s("%i", &lastob);
// 
//         statesize = 6;        // solve for bstar or not
//                               //strBuild.AppendLine("input statesize (6 for posVel or 7 for posVel and B*) \n");
//                               //scanf_s("%i", &statesize, 2 * sizeof(char));
//     }
//     strBuild.AppendLine("fname %s \n", fname);
// 
//     firstob = 0;
//     opsmode = 'a';
//     typeans = 's';
// 
//     int err = fopen_s(&outfile, "testdc.out", "w");
//     err = fopen_s(&outfile1, "testdc1.out", "w");
//     err = fopen_s(&outfile2, "testdc2.out", "w");
// 
//     fstrBuild.AppendLine(outfile, "typeans %c \n", typeans);
//     fstrBuild.AppendLine(outfile, "state type %c \n", statetype);
//     fstrBuild.AppendLine(outfile, "percentchg %lf \n", percentchg);
//     fstrBuild.AppendLine(outfile, "deltaamtchg %11.8lf \n", deltaamtchg);
//     fstrBuild.AppendLine(outfile, "rmsepsilon %lf \n", rmsepsilon);
//     fstrBuild.AppendLine(outfile, "state size %i \n", statesize);
//     fstrBuild.AppendLine(outfile, "pts per period %i \n", lastob);
// 
//     fstrBuild.AppendLine(outfile2, "typeans %c \n", typeans);
//     fstrBuild.AppendLine(outfile2, "state type %c \n", statetype);
//     fstrBuild.AppendLine(outfile2, "percentchg %lf \n", percentchg);
//     fstrBuild.AppendLine(outfile2, "deltaamtchg %11.8lf \n", deltaamtchg);
//     fstrBuild.AppendLine(outfile2, "rmsepsilon %lf \n", rmsepsilon);
//     fstrBuild.AppendLine(outfile2, "state size %i \n", statesize);
//     fstrBuild.AppendLine(outfile2, "pts per period %i \n", lastob);
//     fstrBuild.AppendLine(outfile2, " iter   norad      inti diff (m)  final diff    a           e          i     ap alt   pr alt\n");
// 
//     err = fopen_s(&infile, fname, "r");
//             if (infile == NULL)
//             {
//                 strBuild.AppendLine("Failed to open %s ", fname);
//             }
//             else
//             {
//                 // initialize EOP file
//                 double[] eoparr = new eoprec[];
//                 EOPSPWLibr.readeop(eoparr, "../EOP-All-v1.1_2018-01-04.txt", out jdeopstart, out jdeopstartf);
// 
//                 // --- setup counters for bins of entire catalog
//                 ll1m = 0;
//                 l1m = 0;
//                 l10m = 0;
//                 l50m = 0;
//                 l1km = 0;
//                 l10km = 0;
//                 avgiter = 0;
// 
//                 while (feof(infile) == 0)
//                 {
//                     do
//                     {
//                         fgets(longstr1, 130, infile);
//                         strncpy_s(str, &longstr1[0], 1);
//                         str[1] = '\0';
//                     } while ((strcmp(str, "#") == 0) && (feof(infile) == 0));
// 
//                     // ---- setup the sgp4 vectors from the input file ----
//                     if (feof(infile) == 0)
//                     {
//                         fgets(longstr2, 130, infile);
//                         // convert the char string to sgp4 elements
//                         // includes initialization of sgp4
//                         SGP4Funcs::twoline2rv(longstr1, longstr2, typerun, 'd', opsmode, whichconst,
//                             startmfe, stopmfe, deltamin, satrec);
//                         fstrBuild.AppendLine(outfile, "%5s xx \n", satrec.satnum);
//                         strBuild.AppendLine(" %5s ", satrec.satnum);
//                     } // if going through input file
// 
//                     // set the time span based on the period size
//                     //satrec.a = pow(satrec.no * satrec.tumin, (-2.0 / 3.0));
//                     periodmin = 2.0 * pi * sqrt(pow(satrec.radiusearthkm * satrec.a, 3) / 398600.8) / 60.0;  // in  min
//                     dtmin = int((periodmin) / lastob);
// 
//                     // ---- loop to create ephemeris points (observations) ----------
//                     // ----------------- and setup obsrecarr ------------------------
//                     for (ii = 0; ii <= lastob; ii++)
//                     {
//                         currobsrec.dtmin = dtmin;
//                         SGP4Funcs::sgp4(satrec, dtmin * ii, rpos, vpos);
//                         jd = satrec.jdsatepoch;
//                         jdf = (dtmin * ii) / 1440.0;
//                         MathTimeLibr.invjday(jd, jdf, currobsrec.year, currobsrec.mon, currobsrec.day,
//                             currobsrec.hr, currobsrec.min, currobsrec.sec);
//                         MathTimeLibr.jday(currobsrec.year, currobsrec.mon, currobsrec.day, 0, 0, 0.0, currobsrec.jd, currobsrec.jdf);
//                         currobsrec.sennum = 1;
//                         currobsrec.init = satrec.init;
//                         currobsrec.method = satrec.method;
//                         currobsrec.obstype = 4;  // xyz obs
//                         currobsrec.x = rpos[0];
//                         currobsrec.y = rpos[1];
//                         currobsrec.z = rpos[2];
//                         currobsrec.xdot = vpos[0];
//                         currobsrec.ydot = vpos[1];
//                         currobsrec.zdot = vpos[2];
//                         currobsrec.bstar = satrec.bstar;
//                         obsrecarr[ii] = currobsrec;
//                     } // for ii through producing the observation vectors
// 
//                     // setup nominal vector -----------------------------------------
//                     /*
//                     // create nominal to be a little off
//                     satrec.no = satrec.no*0.999;
//                     satrec.ecco = satrec.ecco*0.999;        //
//                     satrec.inclo = satrec.inclo + 0.01/rad;
//                     satrec.nodeo = satrec.nodeo*1.0002;     //
//                     satrec.argpo = satrec.argpo*1.0001;     //
//                     satrec.mo = satrec.mo;                  //
//                     satrec.bstar = 0.0001; // initial guess only 0.001
//                     */
//                     // -------- or use vector from first call from TLE...
//                     SGP4Funcs::sgp4(satrec, 0.0, rnomTEME, vnomTEME);
//                     for (i = 0; i < 6; i++)
//                     {
//                         if (i < 3)
//                             xnom[i] = rnomTEME[i];
//                         else
//                             xnom[i] = vnomTEME[i - 3];
//                     }
//                     if (statesize > 6)
//                         xnom[6] = 0.0001; // set initial guess as it changes in state2satrec
//                     else
//                         satrec.bstar = 0.0001;
//                     SGP4DC::state2satrec(xnom, jdepoch, jdepochf, 'v', statesize, MathTimeLibr.eTo, satrec); // force this method to assign values
// 
//                     SGP4Funcs::sgp4init(whichconst, opsmode, satrec.satnum, satrec.jdsatepoch - 2433281.5, satrec.bstar,
//                         satrec.ndot, satrec.nddot, satrec.ecco, satrec.argpo, satrec.inclo, satrec.mo, satrec.no_kozai,
//                         satrec.nodeo, satrec);
//                     // call the propagator to get the initial state vector value
//                     SGP4Funcs::sgp4(satrec, 0.0, rnomTEME, vnomTEME);
// 
//                     // print initial difference -----------------------
//                     rpos[0] = obsrecarr[0].x;
//                     rpos[1] = obsrecarr[0].y;
//                     rpos[2] = obsrecarr[0].z;
//                     vpos[0] = obsrecarr[0].xdot;
//                     vpos[1] = obsrecarr[0].ydot;
//                     vpos[2] = obsrecarr[0].zdot;
//                     for (j = 0; j < 3; j++)
//                     {
//                         dr[j] = rnomTEME[j] - rpos[j];
//                         dv[j] = vnomTEME[j] - vpos[j];
//                     }
//                     fstrBuild.AppendLine(outfile, " init diff %11lf  %11lf  %11lf %11lf \n", dr[0], dr[1], dr[2], MathTimeLibr.mag(dr));
//                     magdr1 = 1000.0 * MathTimeLibr.mag(dr); // initial difference in m
//                                                             //               strBuild.AppendLine( " init diff %11lf  %11lf  %11lf %11lf \n",dr[0], dr[1], dr[2], mag(dr) );
// 
//                     //               if (statetype == 'v')
//                     //                 {
//                     //                   // set some rnomTEME vnomTEME value...
// 
//                     //                   rv2coe ( rnomTEME, vnomTEME, p, a, ecc, incl, raan, argp, nu, m, eccanom, arglat, truelon, lonper);
//                     //                   satrec.no = 60.0 * sqrt(3.986008e5/(a * a * a));  // rad per min, sgp4 wgs-72 mu value
//                     //                   satrec.a = a/6378.135;  // er
//                     //                   satrec.ecco = ecc;
//                     //                   satrec.inclo = incl;
//                     //                   satrec.nodeo = raan;
//                     //                   satrec.argpo = argp;
//                     //                   satrec.mo = m;
//                     //                 }
// 
//                     SGP4DC::state2satrec(xnom, jdepoch, jdepochf, statetype, statesize, MathTimeLibr.eFrom, satrec);
// 
//                     //               strBuild.AppendLine( "running new satellite \n" );
// 
//                     // --------- now determine the orbit for the satellite ----------
//                     SGP4DC::leastsquares(percentchg, deltaamtchg, rmsepsilon, whichconst, interp, jdeopstart, satrec, typeans,
//                         statetype, proptype, firstob, lastob, statesize, obsrecarr, loops,
//                         xnom, jdepoch, jdepochf, derivType, cdam, cram, jplopt, jpldearr, jdjpldestart, iau80rec, iau06rec, eoparr, x, dx, atwai, atwa, atwb, outfile, outfile1);
//                     SGP4Funcs::sgp4(satrec, 0.0, rnomTEME, vnomTEME);  // the final answer from LS
//                     rpos[0] = obsrecarr[0].x;
//                     rpos[1] = obsrecarr[0].y;
//                     rpos[2] = obsrecarr[0].z;
//                     vpos[0] = obsrecarr[0].xdot;
//                     vpos[1] = obsrecarr[0].ydot;
//                     vpos[2] = obsrecarr[0].zdot;
//                     for (j = 0; j < 3; j++)
//                     {
//                         dr[j] = rnomTEME[j] - rpos[j];
//                         dv[j] = vnomTEME[j] - vpos[j];
//                     }
//                     magdr = MathTimeLibr.mag(dr) * 1000.0;
//                     fstrBuild.AppendLine(outfile, " difference %11lf  %11lf  %11lf %11lf \n", dr[0], dr[1], dr[2], MathTimeLibr.mag(dr));
//                     fstrBuild.AppendLine(outfile1, " %8s  difference %11lf  %11lf  %11lf km %12lf m \n\n", satrec.satnum, dr[0], dr[1], dr[2], 1000.0 * MathTimeLibr.mag(dr));
//                     apalt = 6378.135 * satrec.a * (1.0 + satrec.ecco) - 6378.135;
//                     pralt = 6378.135 * satrec.a * (1.0 - satrec.ecco) - 6378.135;
//                     //               magcov =  1000.0 * sqrt(atwai[0][0] + atwai[1][1] + atwai[2][2]);
//                     magcov = sqrt(atwai[2][2]) * 6378135.0;  // in m
// 
//                     if (1000.0 * MathTimeLibr.mag(dr) < 0.0001)
//                         fstrBuild.AppendLine(outfile2, " %4i %8s %12.0lf %14.6lf %10.2lf %12lf %8.4lf %7.0lf %7.0lf %10.2lf\n",
//                         loops, satrec.satnum, magdr1, 0.000001, satrec.a * 6378, satrec.ecco, satrec.inclo * rad, apalt, pralt, magcov);
//                     else
//                         fstrBuild.AppendLine(outfile2, " %4i %8s %12.0lf %14.6lf %10.2lf %12lf %8.4lf %7.0lf %7.0lf %10.2lf\n",
//                         loops, satrec.satnum, magdr1, 1000.0 * MathTimeLibr.mag(dr), satrec.a * 6378, satrec.ecco, satrec.inclo * rad, apalt, pralt, magcov);
// 
//                     strBuild.AppendLine("%4i %8.1lf %8.4lf m a %10.5lf e %9.6lf i %7.3lf %9.4lf m\n", loops, magdr1, 1000.0 * MathTimeLibr.mag(dr), satrec.a, satrec.ecco, satrec.inclo * rad,
//                         magcov);
// 
// 
//                     if (magdr < 1.0)
//                         ll1m = ll1m + 1;
//                     else
//                     {
//                         if (magdr < 10.0)
//                             l1m = l1m + 1;
//                         else
//                         {
//                             if (magdr < 100.0)
//                                 l10m = l10m + 1;
//                             else
//                             {
//                                 if (magdr < 1000.0)
//                                     l50m = l50m + 1;
//                                 else
//                                 {
//                                     if (magdr < 10000.0)
//                                         l1km = l1km + 1;
//                                     else
//                                         l10km = l10km + 1;
//                                 }
//                             }
//                         }
//                     }
//                     avgiter = avgiter + loops;
// 
//                 }  // while not eof
// 
//                 totalnum = ll1m + l1m + l10m + l50m + l1km + l10km;
//                 strBuild.AppendLine("%5i %8.2lf      < 1m  avgiter = %8.2lf \n", ll1m, 100.0 * (ll1m / totalnum), avgiter / totalnum);
//                 strBuild.AppendLine("%5i %8.2lf   1m < 10m \n", l1m, 100.0 * ((ll1m + l1m) / totalnum));
//                 strBuild.AppendLine("%5i %8.2lf  10m < 100m \n", l10m, 100.0 * ((ll1m + l1m + l10m) / totalnum));
//                 strBuild.AppendLine("%5i %8.2lf 100m < 1km \n", l50m, 100.0 * ((ll1m + l1m + l10m + l50m) / totalnum));
//                 strBuild.AppendLine("%5i %8.2lf  1km < 10km \n", l1km, 100.0 * ((ll1m + l1m + l10m + l50m + l1km) / totalnum));
//                 strBuild.AppendLine("%5i %8.2lf 10km > \n", l10km, 100.0 * ((ll1m + l1m + l10m + l50m + l1km + l10km) / totalnum));
// 
//                 fstrBuild.AppendLine(outfile, "%5i %8.2lf      < 1m  avgiter = %8.2lf \n", ll1m, 100.0 * (ll1m / totalnum), avgiter / totalnum);
//                 fstrBuild.AppendLine(outfile, "%5i %8.2lf      < 1m \n", ll1m, 100.0 * (ll1m / totalnum));
//                 fstrBuild.AppendLine(outfile, "%5i %8.2lf   1m < 10m \n", l1m, 100.0 * ((ll1m + l1m) / totalnum));
//                 fstrBuild.AppendLine(outfile, "%5i %8.2lf  10m < 100m \n", l10m, 100.0 * ((ll1m + l1m + l10m) / totalnum));
//                 fstrBuild.AppendLine(outfile, "%5i %8.2lf 100m < 1km \n", l50m, 100.0 * ((ll1m + l1m + l10m + l50m) / totalnum));
//                 fstrBuild.AppendLine(outfile, "%5i %8.2lf  1km < 10km \n", l1km, 100.0 * ((ll1m + l1m + l10m + l50m + l1km) / totalnum));
//                 fstrBuild.AppendLine(outfile, "%5i %8.2lf 10km > \n", l10km, 100.0 * ((ll1m + l1m + l10m + l50m + l1km + l10km) / totalnum));
// 
//                 fclose(infile);
//                 fclose(outfile);
//                 fclose(outfile1);
//             }
// 
// 
//     // Read an STK epehemris
//     // by Alek Lidtke
//     //void LoadEphemerisTable(std::string fullFileName)  // ExternalSpaceObject::
//     //{
//     //	/* Read a STK v 10 ephemeris file that contains the ephemeris of this object. The epochs in Julian Days of all the ephemeris points will
//     //	be stored in EphemerisEpochs attribute, the Cartesian positions in km in EphemerisPositions, and the Cartesian velocities
//     //	in km/sec in EphemerisVelocities. The state vectors are defined in the True Equator Mean Equinox reference frame.
//     //	@param - name of a file that contains Julian day epoch, Cartesian positions and velocities in meteres and metres per second in TEME reference frame.
//     //	*/
//     //	EphemerisEpochs = double[]<double>(); // Initialise the vectors that will hold the ephemeris points.
//     //	EphemerisPositions = double[]< double[]<double> >();
//     //	EphemerisVelocities = double[]< double[]<double> >();
// 
//     //	/* Parse the ephemeris file. */
//     //	std::ifstream TLEfileStream(fullFileName, std::ifstream::in);
// 
//     //	std::string EphemerisLine; // Currently read ephemeris point.
// 
//     //	int noEphemPoints; // Counter of the total number of ephemeris points read.
//     //	double JDAYepoch;
//     //	int EphemerisPointsRead = 0; // Counter.
//     //	double temp; // Temporary number read from the file.
// 
//     //	while (std::getline(TLEfileStream, EphemerisLine))
//     //	{
//     //		std::istringstream iss(EphemerisLine);
// 
//     //		/* Get information about the ephemeris file being parsed from the header and parse the ephemeris points. */
//     //		if (EphemerisLine.find("NumberOfEphemerisPoints") == 0){
//     //			std::string noEphemPointsString = EphemerisLine.substr(EphemerisLine.find(" ") + 1); // Don't want the white space before the number of points.
//     //			std::basic_istringstream<char> noEphemPointsSS(noEphemPointsString); // Convert fromstring to an int.
//     //			noEphemPointsSS >> noEphemPoints;
//     //			EphemerisEpochs.resize(noEphemPoints); // Allocate memory for the ephemeris points.
//     //			EphemerisPositions.resize(noEphemPoints);
//     //			EphemerisVelocities.resize(noEphemPoints);
//     //			for (double[]<int>::size_type i = 0; i<EphemerisEpochs.size(); i++){
//     //				EphemerisPositions.at(i).resize(3); // Resize every entry to be able to hold epoch, 3 Cartesian poisiton and 3 Cartesian velocity compontents (Julian Days, metres, and metres per second, respectively).
//     //				EphemerisVelocities.at(i).resize(3);
//     //			}
//     //		}
//     //		else if (EphemerisLine.find("# Time of first point:") == 0){
//     //			std::string JDAYepochString = EphemerisLine.substr(EphemerisLine.find("UTCG =") + 6, 14); // Only get this part of the string that contains the epoch of the first point.
//     //			std::basic_istringstream<char> JDAYepochSS(JDAYepochString);
//     //			JDAYepochSS >> JDAYepoch; // Find the epoch of the first point - will get offset times from this epoch in seconds for every ephemeris point.
//     //		}
//     //		else if (EphemerisLine.size() != 0){ // Don't try to call .at(0) for empty strings.
//     //			if (isdigit(EphemerisLine.at(0))){ // Only actual ephemeris points will start with digits.
//     //				std::basic_istringstream<char> ephemPointSS(EphemerisLine); // Parse this ephemeris point - always 7 numbers in scientific notation but don't preassume anyting.
//     //				for (int counter = 0; counter<7; counter++){ // Read the epoch, position, and velocity for this point.
//     //					if (counter == 0){
//     //						ephemPointSS >> temp; // Change the epoch to Julian Days from the time offset in seconds from the first point.
//     //						EphemerisEpochs.at(EphemerisPointsRead) = JDAYepoch + temp*SECOND_IN_JDAYS;
//     //					}
//     //					else if ((counter>0) && (counter<4)){ // Now parsing position components.
//     //						ephemPointSS >> temp;
//     //						EphemerisPositions.at(EphemerisPointsRead).at(counter - 1) = temp / 1000.0; // Change to km.
//     //					}
//     //					else if ((counter>3) && (counter<7)){ // Now parsing velocity compontents.
//     //						ephemPointSS >> temp;
//     //						EphemerisVelocities.at(EphemerisPointsRead).at(counter - 4) = temp / 1000.0;
//     //					}
//     //				}
//     //				EphemerisPointsRead += 1; // Make sure this ephemeris point won't be overwritten.
//     //			}
//     //		}
//     //	};
// 
//     //	int LastLineInfo = EphemerisLine.find("END Ephemeris");
//     //	if ((LastLineInfo == 0) || (EphemerisLine.empty()))
//     //	{ // Reached the last line of the file.
//     //		std::cout << "Parsed entre ephemeris file " << fullFileName << " and read " << EphemerisEpochs.size() << " ephemeris points." << std::endl;
//     //	}
//     //	else if ((LastLineInfo != 0) && (!EphemerisLine.empty()))
//     //	{ // Haven't reached the last line of the file. Last line may be empty.
//     //		std::cerr << "Haven't reached the end of the ephemeris file while reading the ephemeris from " << fullFileName << "." << std::endl;
//     //	}
//     //	else if (noEphemPoints != (int)EphemerisEpochs.size())
//     //	{
//     //		std::cerr << "Haven't read as many ephemeris points as quoted in the file header for " << fullFileName << "." << std::endl;
//     //	};
//     //}  // LoadEphemerisTable
// 
// 
// }  // processtle


    }
}
