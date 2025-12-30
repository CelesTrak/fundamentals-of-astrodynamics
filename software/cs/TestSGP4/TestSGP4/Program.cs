using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Text.RegularExpressions;
using System.IO;

using SGP4Methods;

namespace SGP4Program
{
    class SGP4Programs
    {

        /* ---------------------------------------------------------------------
        *
        *                              testSGP4.cs
        *
        *  this program tests the sgp4 propagator. an stk ephemeris file is generated
        *  along with the test output. the code for this is left justified for easy
        *  location.
        *
        *                          companion code for
        *             fundamentals of astrodynamics and applications
        *                                  2022
        *                            by david vallado
        *
        *     email dvallado@comspoc.com, davallado@gmail.com
        *     *****************************************************************
        *    current :
        *              30 dec 25  david vallado
        *                           add alpha5 support
        *    changes :
        *              29 aug 24  david vallado
        *                           add check for sgp4-xp tle
        *               7 dec 15  david vallado
        *                           fix jd, jdfrac
        *               3 nov 14  david vallado
        *                           update to msvs2013 c++
        *              11 nov 13  david vallado
        *                           conversion to msvs c++
        *                           misc fixes to constants and options
        *                           add singly averaged state elements to be exported
        *               3 sep 08  david vallado
        *                           add switch for afspc compatibility and improved operation
        *              14 may 08  david vallado
        *                           fixes for linux suggested by brian micek
        *                           misc fixes noted by the community - manual operation,
        *                           formats, char lengths
        *              14 aug 06  david vallado
        *                           update mfe for verification time steps, constants
        *              20 jul 05  david vallado
        *                           fixes for paper, corrections from paul crawford
        *               7 jul 04  david vallado
        *                           fix record file and get working
        *              14 may 01  david vallado
        *                           2nd edition baseline
        *                     80  norad
        *                           original baseline
        *       ----------------------------------------------------------------      */

        public static void Main(string[] args)
        {
            // setup the class so methods can be called
            SGP4Lib SGP4Libr = new SGP4Lib();
            
            string infilename, outfilename, outfilenameE;
            double[] ro = new double[3];
            double[] vo = new double[3];
            char typerun, typeinput, opsmode;
            SGP4Lib.gravconsttype whichconst;
            int whichcon;
            outfilenameE = "";

            // ----------------------------  locals  -------------------------------
            double p, a, ecc, incl, node, argp, nu, m, arglat, truelon, lonper;
            double sec, jd, jdfrac, rad, tsince, startmfe, stopmfe, deltamin;
            int year; int mon; int day; int hr; int min;
            string longstr1 = "";
            string longstr2 = "";
            SGP4Lib.elsetrec satrec = new SGP4Lib.elsetrec();

            //SGP4Libr.InitElset(ref satrec);

            whichconst = SGP4Lib.gravconsttype.wgs72;

            StringBuilder strbuild = new StringBuilder();
            StringBuilder strbuildE = new StringBuilder();   // for STK e file

            rad = 180.0 / Math.PI;
            // ------------------------  implementation   --------------------------
            string[] monstr = new string[13] { "", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC" };

            Console.Write("%s\n", SGP4Libr.SGP4Version);

            //opsmode = "a" best understanding of how afspc code works
            //opsmode = "i" improved sgp4 resulting in smoother behavior
            Console.Write("input operation mode (a), i \n");
            string linex = Console.ReadLine(); // Read string from console
            opsmode = linex[0];

            Console.Write("typerun = c compare 1 year of full satcat data \n",
                          "typerun = v verification run, requires modified elm file with \n",
                          "              start, stop, and delta times \n",
                          "typerun = m manual operation- either mfe, epoch, or day of yr \n");
            Console.Write("input type of run c, v, m \n");
            linex = Console.ReadLine(); // Read string from console
            typerun = linex[0];

            //typeinput = 'm' input start stop mfe
            //typeinput = 'e' input start stop ymd hms
            //typeinput = 'd' input start stop yr dayofyr
            if ((typerun != 'v') && (typerun != 'c'))
            {
                Console.Write("input mfe, epoch (YMDHMS), or dayofyr approach, m,e,d \n");
                linex = Console.ReadLine(); // Read string from console
                typeinput = linex[0];
            }
            else
                typeinput = 'e';

            //Console.Write("input which constants 721 (72) 84 \n");
            //linex = Console.ReadLine(); // Read string from console
            //if (int.TryParse(linex, out whichcon)) // Try to parse the string as an integer
            //{
            //}
            //else
            //{
            //    Console.WriteLine("Not an integer!");
            //}
            // standard operation 
            whichcon = 72;
            //if (whichcon == 721) whichconst = SGP4Lib.gravconsttype.wgs72old;
            if (whichcon == 72) whichconst = SGP4Lib.gravconsttype.wgs72;
            //if (whichcon == 84) whichconst = SGP4Lib.gravconsttype.wgs84;

            // ---------------- setup files for operation ------------------
            // input 2-line element set file
            Console.Write("input elset filename: \n");
            infilename = Console.ReadLine();

            string[] TLEFileData = File.ReadAllLines(infilename);

            if (typerun == 'c')
                outfilename = "tcsall.out";
            else
            {
                if (typerun == 'v')
                    outfilename = "tcsver.out";
                else
                    outfilename = "tcs.out";
            }

            //        dbgfile = fopen("sgp4test.dbg", "w");
            //        fprintf(dbgfile,"this is the debug output\n\n" );

            // ----------------- test simple propagation -------------------
            // loop through the file
            int ktr = 0;
            int ktr1 = 0;
            while (ktr < TLEFileData.Count())
            {
                // sgp4fix addiional parameters to store from the TLE
                satrec.classification = 'U';
                satrec.intldesg = "        ";
                satrec.ephtype = 0;
                satrec.elnum = 0;
                satrec.revnum = 0;

                longstr1 = TLEFileData[ktr];
                while ((TLEFileData[ktr].Contains('#')) && (ktr < TLEFileData.Count()))
                {
                    ktr = ktr + 1;
                    longstr1 = TLEFileData[ktr];
                } // while skipping comment lines

                ktr = ktr + 1;
                longstr2 = TLEFileData[ktr];
                while ((TLEFileData[ktr].Contains('#')) && (ktr < TLEFileData.Count()))
                {
                    ktr = ktr + 1;
                    longstr2 = TLEFileData[ktr];
                } // while skipping comment lines

                if (ktr < TLEFileData.Count())
                {
                    ktr1 = ktr1 + 1; // increment for the .e file
                    // convert the char string to sgp4 elements
                    // includes initialization of sgp4
                    SGP4Libr.twoline2rv(longstr1, longstr2, typerun, typeinput, opsmode, whichconst,
                                out startmfe, out stopmfe, out deltamin, out satrec);

                    // sgp4fix note that the ephtype must be 0 for SGP4. SGP4-XP uses 4.
                    if (satrec.ephtype == 0)
                    {
                        // call the propagator to get the initial state vector value
                        // no longer need gravconst since it is assigned in sgp4init
                        SGP4Libr.sgp4(ref satrec, 0.0, ro, vo);

                        // generate .e files for stk
                        jd = satrec.jdsatepoch;
                        outfilenameE = longstr2.Substring(2, 5).ToString() + ".e";
                        SGP4Libr.invjday(satrec.jdsatepoch, satrec.jdsatepochF, out year, out mon, out day, out hr, out min, out sec);
                        //string datestr = String.Format(" {0,2}{1,4}{2,5} {3,2}:{4,2}:{5} ", day, monstr[mon], year, hr, min, sec);
                        // add datetime support 29 Jul 2022 dav
                        DateTime dt = new DateTime(year, mon, day, hr, min, 0); // TS
                        dt = dt.AddSeconds(sec); // TS: Supports double seconds
                        string datestr = $"{dt:dd MMM yyyy HH:mm:ss.fff}"; // TS

                        strbuildE.AppendLine("stk.v.12.0\n");
                        strbuildE.AppendLine("BEGIN Ephemeris\n");
                        strbuildE.AppendLine("NumberOfEphemerisPoints xxxxx");
                        strbuildE.AppendLine("ScenarioEpoch  " + datestr);
                        strbuildE.AppendLine("InterpolationMethod		Lagrange");
                        strbuildE.AppendLine("InterpolationOrder		5");
                        strbuildE.AppendLine("CentralBody Earth");
                        strbuildE.AppendLine("CoordinateSystem TEME");
                        strbuildE.AppendLine("DistanceUnit Kilometers");   // Kilometers
                        strbuildE.AppendLine("EphemerisTimePosVel\n");

                        // .e file
                        tsince = 0.0;
                        strbuildE.AppendLine((tsince.ToString("0.00000000")).PadLeft(17) +
                                            (ro[0].ToString("0.00000000")).PadLeft(17) + " " + ro[1].ToString("0.00000000").PadLeft(17) + " " + ro[2].ToString("0.00000000").PadLeft(17) + " " +
                                            (vo[0].ToString("0.00000000")).PadLeft(17) + " " + vo[1].ToString("0.000000000").PadLeft(17) + " " + vo[2].ToString("0.000000000").PadLeft(17));

                        // output file
                        strbuild.AppendLine(satrec.satnumStr.ToString() + " xx");
                        strbuild.AppendLine((tsince.ToString("0.00000000")).PadLeft(17) +
                                            ro[0].ToString("0.00000000").PadLeft(17) + " " + ro[1].ToString("0.00000000").PadLeft(17) + " " + ro[2].ToString("0.00000000").PadLeft(17) + " " +
                                            vo[0].ToString("0.000000000").PadLeft(17) + " " + vo[1].ToString("0.000000000").PadLeft(17) + " " + vo[2].ToString("0.000000000").PadLeft(17));
                        tsince = startmfe;

                        // check so the first value isn't written twice
                        if (Math.Abs(tsince) > 1.0e-8)
                            tsince = tsince - deltamin;

                        // ----------------- loop to perform the propagation ----------------
                        while ((tsince < stopmfe) && (satrec.error == 0))
                        {
                            tsince = tsince + deltamin;

                            if (tsince > stopmfe)
                                tsince = stopmfe;

                            SGP4Libr.sgp4(ref satrec, tsince, ro, vo);

                            if (satrec.error > 0)
                                Console.Write(@"# *** error: t:= %f *** code = %3d\n",
                                        satrec.t, satrec.error);

                            if (satrec.error == 0)
                            {
                                if ((typerun != 'v') && (typerun != 'c'))
                                {
                                    jd = satrec.jdsatepoch;
                                    jdfrac = satrec.jdsatepochF + tsince / 1440.0;
                                    if (jdfrac < 0.0)
                                    {
                                        jd = jd - 1.0;
                                        jdfrac = jdfrac + 1.0;
                                    }
                                    SGP4Libr.invjday(jd, jdfrac, out year, out mon, out day, out hr, out min, out sec);
                                    strbuild.AppendLine(tsince.ToString("0.00000000").PadLeft(17) +
                                                        ro[0].ToString("0.00000000").PadLeft(17) + " " + ro[1].ToString("0.00000000").PadLeft(17) + " " + ro[2].ToString("0.00000000").PadLeft(17) + " " +
                                                        vo[0].ToString("0.000000000").PadLeft(17) + " " + vo[1].ToString("0.000000000").PadLeft(17) + " " + vo[2].ToString("0.000000000").PadLeft(17) + " " +
                                                        year.ToString() + " " + mon.ToString() + " " + day.ToString() + " " + hr.ToString() + " " + min.ToString() + " " + sec.ToString("0.000000")
                                                       );
                                }
                                else
                                {
                                    jd = satrec.jdsatepoch;
                                    jdfrac = satrec.jdsatepochF + tsince / 1440.0;
                                    if (jdfrac < 0.0)
                                    {
                                        jd = jd - 1.0;
                                        jdfrac = jdfrac + 1.0;
                                    }
                                    SGP4Libr.invjday(jd, jdfrac, out year, out mon, out day, out hr, out min, out sec);

                                    strbuildE.AppendLine((tsince * 60).ToString("0.00000000").PadLeft(17) +
                                                       ro[0].ToString("0.00000000").PadLeft(17) + " " + ro[1].ToString("0.00000000").PadLeft(17) + " " + ro[2].ToString("0.00000000").PadLeft(17) + " " +
                                                       vo[0].ToString("0.000000000").PadLeft(17) + " " + vo[1].ToString("0.000000000").PadLeft(17) + " " + vo[2].ToString("0.000000000").PadLeft(17));

                                    SGP4Libr.rv2coe(ro, vo, satrec.mu, out p, out a, out ecc, out incl, out node, out argp, out nu, out m, out arglat, out truelon, out lonper);
                                    strbuild.AppendLine(tsince.ToString("0.00000000").PadLeft(17) +
                                                        ro[0].ToString("0.00000000").PadLeft(17) + " " + ro[1].ToString("0.00000000").PadLeft(17) + " " + ro[2].ToString("0.00000000").PadLeft(17) + " " +
                                                        vo[0].ToString("0.000000000").PadLeft(17) + " " + vo[1].ToString("0.000000000").PadLeft(17) + " " + vo[2].ToString("0.000000000").PadLeft(17) + " " +
                                                        a.ToString("0.000000").PadLeft(17) + " " + ecc.ToString("0.000000") + " " + (incl * rad).ToString("0.00000") + " " +
                                                        (node * rad).ToString("0.00000") + " " + (argp * rad).ToString("0.00000") + " " + (nu * rad).ToString("0.00000") + " " + (m * rad).ToString("0.00000") + " " +
                                                        year.ToString() + " " + mon.ToString() + " " + day.ToString() + " " + hr.ToString() + ":" + min.ToString() + ":" + sec.ToString("0.000000")
                                                       );
                                }
                            } // if satrec.error == 0

                        } // while propagating the orbit

                        strbuildE.Replace("xxxxx", ktr1.ToString());
                        strbuildE.AppendLine("END Ephemeris\n");

                        File.WriteAllText(outfilenameE, strbuildE.ToString());
                        strbuildE.Clear();
                    }  // if ephtype == 0

                } // if not eof

                ktr = ktr + 1;
            } // while through the input file

            File.WriteAllText(outfilename, strbuild.ToString());
            strbuild.Clear();


            // sgp4fix demonstrate method of running SGP4 directly from orbital element values
            //1 08195U 75081A   06176.33215444  .00000099  00000-0  11873-3 0   813
            //2 08195  64.1586 279.0717 6877146 264.7651  20.2257  2.00491383225656
            const double deg2rad = Math.PI / 180.0;         //   0.0174532925199433
            const double xpdotp = 1440.0 / (2.0 * Math.PI);  // 229.1831180523293

            whichconst = SGP4Lib.gravconsttype.wgs72;
            opsmode = 'a';
            satrec.satnumStr = "08195";
            satrec.jdsatepoch = 2453911.8321544402;
            satrec.no_kozai = 2.00491383;
            satrec.ecco = 0.6877146;
            satrec.inclo = 64.1586;
            satrec.nodeo = 279.0717;
            satrec.argpo = 264.7651;
            satrec.mo = 20.2257;
            satrec.nddot = 0.00000e0;
            satrec.bstar = 0.11873e-3;
            satrec.ndot = 0.00000099;
            satrec.elnum = 813;
            satrec.revnum = 22565;
            satrec.classification = 'U';
            satrec.intldesg = "        ";
            satrec.ephtype = 0;

            // convert units and initialize
            satrec.no_kozai = satrec.no_kozai / xpdotp; //* rad/min
            satrec.ndot = satrec.ndot / (xpdotp * 1440.0);  //* ? * minperday
            satrec.nddot = satrec.nddot / (xpdotp * 1440.0 * 1440);
            satrec.inclo = satrec.inclo * deg2rad;
            satrec.nodeo = satrec.nodeo * deg2rad;
            satrec.argpo = satrec.argpo * deg2rad;
            satrec.mo = satrec.mo * deg2rad;

            // set start/stop times for propagation
            startmfe = 0.0;
            stopmfe = 2880.0;
            deltamin = 120.0;

            SGP4Libr.sgp4init(whichconst, opsmode, satrec.satnumStr, satrec.jdsatepoch - 2433281.5, satrec.bstar,
                   satrec.ndot, satrec.nddot, satrec.ecco, satrec.argpo, satrec.inclo, satrec.mo, satrec.no_kozai,
                   satrec.nodeo, ref satrec);

            tsince = startmfe;
            while ((tsince < stopmfe) && (satrec.error == 0))
            {
                tsince = tsince + deltamin;

                if (tsince > stopmfe)
                    tsince = stopmfe;

                SGP4Libr.sgp4(ref satrec, tsince, ro, vo);

                jd = satrec.jdsatepoch;
                jdfrac = satrec.jdsatepochF + tsince / 1440.0;
                if (jdfrac < 0.0)
                {
                    jd = jd - 1.0;
                    jdfrac = jdfrac + 1.0;
                }
                SGP4Libr.invjday(jd, jdfrac, out year, out mon, out day, out hr, out min, out sec);

                strbuild.AppendLine(tsince.ToString("0.00000000").PadLeft(5) +
                                    ro[0].ToString("0.0000000000").PadLeft(5) + " " + ro[1].ToString("0.0000000000").PadLeft(5) + " " + ro[2].ToString("0.0000000000").PadLeft(5) + " " +
                                    vo[0].ToString("0.00000000000").PadLeft(5) + " " + vo[1].ToString("0.00000000000").PadLeft(5) + " " + vo[2].ToString("0.00000000000").PadLeft(5));
            } // while propagating the orbit

        }  // testSGP4    

    } // class SGP4Programs

}  // namespace

