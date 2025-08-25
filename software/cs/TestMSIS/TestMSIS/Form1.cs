// TestMSIS.cpp : main project file.
//
// test driver for msis00 model
// dav 12 nov 2004 fix ap array set (to 8) and include sw tests for each
// dav 4 oct 2004
//
// The NRLMSISE-00 model was developed by Mike Picone, Alan Hedin, and
// Doug Drob. They also wrote a NRLMSISE-00 distribution package in
// FORTRAN which is available at
// http://uap-www.nrl.navy.mil/models_web/msis/msis_home.htm
//
// 
// remove .h from include do not do it for other libraries,
// remove clrscr(); everywhere and replace by system("cls");
// finally after writing all #include stuff add this-
// using namespace std; it allows you to do things like cin>> cout<< etc.
//

//using static MSIS_Methods.MSISLib;


using EOPSPWMethods;       // EOPDataClass, SPWDataClass, iau80Class, iau06Class
using MathTimeMethods;     // Edirection, globals
//using AstroLibMethods;     // EOpt, gravityConst, astroConst, xysdataClass, jpldedataClass
using MSIS_Methods;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Windows.Forms.VisualStyles;


namespace TestMSIS
{
    public partial class Form1 : Form
    {
        public string errstr;

        public MathTimeLib MathTimeLibr = new MathTimeLib();
        public EOPSPWLib EOPSPWLibr = new EOPSPWLib();
        public MSISLib MSISLibr = new MSISLib();

        public class msis_inputCl
        {
            public int year;      // year, currently ignored 
            public int doy;       // day of year 
            public double sec;    // seconds in day (UT) 
            public double alt;    // altitude in kilometers 
            public double glat;   // geodetic latitude 
            public double glong;  // geodetic longitude 
            public double stl;    // local apparent solar time (hours), see note below 
            public double f107a;  // 81 day average of F10.7 flux (centered on doy) 
            public double f107;   // daily F10.7 flux for previous day 
        }
        //public msis_inputCl msisinput = new msis_inputCl();


        public Form1()
        {
            InitializeComponent();
        }


        private void btnSearch_Click(object sender, EventArgs e)
        {
            msis_inputCl[] input = new msis_inputCl[17];
            int i, mass;
            MSISLib.msistype msis86r, msis90r, msis00r;
            MSISLib.lpolytype lpoly = new MSISLib.lpolytype();
            MSISLib.fittype fit = new MSISLib.fittype();
            MSISLib.lsqvtype lsqv = new MSISLib.lsqvtype();
            double[] d= new double[10];
            double[] t = new double[3];
            double[] aph = new double[8];
            double[] apar = new double[8];
            double[] apar9 = new double[8];

            /* input values */
            for (i = 0; i < 8; i++)
            {
                aph[i] = 100;
                apar[i] = 4;
                apar9[i] = 40;
            }

            mass = 48;
            Case 3 is off
            for (i = 0; i < 17; i++)
            {
                input[i] = new msis_inputCl();

                input[i].doy = 172;
                input[i].year = 0; /* without effect */
                input[i].sec = 29000;
                input[i].alt = 400;
                input[i].glat = 60;
                input[i].glong = -70;
                input[i].stl = 16;
                input[i].f107a = 150;
                input[i].f107 = 150;
            }

            // setup differences for each test
            input[1].doy = 81;
            input[2].sec = 75000;
            input[2].alt = 1000;
            input[3].alt = 100;
            input[4].glat = 0;
            input[5].glong = 0;
            input[6].stl = 4;
            input[7].f107a = 70;
            input[8].f107 = 180;
            input[10].alt = 0;
            input[11].alt = 10;
            input[12].alt = 30;
            input[13].alt = 50;
            input[14].alt = 70;
            input[16].alt = 100;

            errstr = errstr + ("now do the msis-86 test cases --------------------------------------\n");
            MSISLibr.msis86init(out msis86r);

            /* evaluate std test cases */
            for (i = 0; i < 17; i++)
            {
                errstr = errstr + "------case number " + i + "\n";
                if (i == 2)  // it's different for the msis00 test case
                    input[i].alt = 400.0;
                if (i == 9)
                    MSISLibr.gts5(ref msis86r, ref lpoly, ref lsqv,
                         input[i].doy, input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                         input[i].stl, input[i].f107a, input[i].f107, apar9, mass, d, t);
                else
                {
                    if (i >= 15)
                    {
                        // test switch settings
                        MSISLibr.tselec(ref msis86r.csw, svt);
                        MSISLibr.gts5(ref msis86r, ref lpoly, ref lsqv,
                             input[i].doy, input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                             input[i].stl, input[i].f107a, input[i].f107, aph, mass, d, t);
                    }
                    else
                        MSISLibr.gts5(ref msis86r, ref lpoly, ref lsqv,
                           input[i].doy, input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                           input[i].stl, input[i].f107a, input[i].f107, apar, mass, d, t);
                }
                errstr = errstr + input[i].doy + " " + input[i].sec + ".00 " +
                    input[i].alt + ".00 " +
                    input[i].glat + ".00 " + input[i].glong + ".00 " +
                    input[i].stl + ".00 " + input[i].f107a + ".00 " +
                    input[i].f107 + ".00\n";
                errstr = errstr + t[1] + " " + t[2] + " " 
                    + d[1].ToString("E6") + " " + d[2].ToString("E6") + " " + d[3].ToString("E6") + " " 
                    + d[4].ToString("E6") + " " 
                    + d[5].ToString("E6") + " " + d[6].ToString("E6") + " " + d[7].ToString("E6") + " " 
                    + d[8].ToString("E6") + "\n";
            }

            string directory = @"D:\Codes\LIBRARY\cs\TestMSIS\";
            //System.IO.File.WriteAllText(directory + "tmsiscs.out", errstr);

            errstr = errstr + ("now do the msis-90 test cases --------------------------------------\n");
            MSISLibr.msis90init(out msis90r);

            /* evaluate std test cases */
            for (i = 0; i < 17; i++)
            {
                errstr = errstr + "------case number " + i + "\n";
                if (i == 2)  // it's different for the msis00 test case
                    input[i].alt = 400.0;
                if (i == 9)
                    MSISLibr.gtd6(ref msis90r, ref lpoly, ref fit, ref lsqv,
                        input[i].doy, input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                        input[i].stl, input[i].f107a, input[i].f107, apar9, mass, d, t);
                else
                {
                    if (i >= 15)
                    {
                        // test switch settings
                        MSISLibr.tselec(ref msis90r.csw, svt);
                        MSISLibr.gtd6(ref msis90r, ref lpoly, ref fit, ref lsqv,
                            input[i].doy, input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                            input[i].stl, input[i].f107a, input[i].f107, aph, mass, d, t);
                    }
                    else
                        MSISLibr.gtd6(ref msis90r, ref lpoly, ref fit, ref lsqv,
                          input[i].doy, input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                          input[i].stl, input[i].f107a, input[i].f107, apar, mass, d, t);
                }
                errstr = errstr + input[i].doy + " " + input[i].sec + ".00 " + 
                    input[i].alt + ".00 " +
                    input[i].glat + ".00 " + input[i].glong + ".00 " +
                    input[i].stl + ".00 " + input[i].f107a + ".00 " +
                    input[i].f107 + ".00\n";
                errstr = errstr + t[1] + " " + t[2] + " "
                    + d[1].ToString("E6") + " " + d[2].ToString("E6") + " " + d[3].ToString("E6") + " "
                    + d[4].ToString("E6") + " "
                    + d[5].ToString("E6") + " " + d[6].ToString("E6") + " " + d[7].ToString("E6") + " "
                    + d[8].ToString("E6") + "\n";
            }

            System.IO.File.WriteAllText(directory + "tmsiscs.out", errstr);

            errstr = errstr + ("now do the msis-00 test cases --------------------------------------\n");
            MSISLibr.msis00init(out msis00r);

            /* evaluate std test cases */
            for (i = 0; i < 17; i++)
            {
                errstr = errstr + "------case number " + i + "\n";
                if (i == 2)   // alt above 500 km
                {
                    input[i].alt = 1000.0;
                    errstr = errstr + ("-- high alt, gtd7 results \n");
                    MSISLibr.gtd7(ref msis00r, ref lpoly, ref fit, ref lsqv,
                      input[i].doy, input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                      input[i].stl, input[i].f107a, input[i].f107, apar, mass, d, t);
                    errstr = errstr + input[i].doy + " " + input[i].sec + ".00 " +
                        input[i].alt + ".00 " +
                        input[i].glat + ".00 " + input[i].glong + ".00 " +
                        input[i].stl + ".00 " + input[i].f107a + ".00 " +
                        input[i].f107 + ".00\n";
                    errstr = errstr + t[1] + " " + t[2] + " "
                        + d[1].ToString("E6") + " " + d[2].ToString("E6") + " " + d[3].ToString("E6") + " "
                        + d[4].ToString("E6") + " "
                        + d[5].ToString("E6") + " " + d[6].ToString("E6") + " " + d[7].ToString("E6") + " "
                        + d[8].ToString("E6") + "\n";
                    errstr = errstr + ("-- high alt, gtd7d results \n");
                    MSISLibr.gtd7d(ref msis00r, ref lpoly, ref fit, ref lsqv,
                      input[i].doy, input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                      input[i].stl, input[i].f107a, input[i].f107, apar, mass, d, t);
                }
                else
                    if (i == 9)
                    MSISLibr.gtd7(ref msis00r, ref lpoly, ref fit, ref lsqv,
                      input[i].doy, input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                      input[i].stl, input[i].f107a, input[i].f107, apar9, mass, d, t);
                else
                {
                    if (i >= 15)
                    {
                        // test switch settings
                        MSISLibr.tselec(ref msis00r.csw, svt);
                        MSISLibr.gtd7(ref msis00r, ref lpoly, ref fit, ref lsqv,
                          input[i].doy, input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                          input[i].stl, input[i].f107a, input[i].f107, aph, mass, d, t);
                    }
                    else
                        MSISLibr.gtd7(ref msis00r, ref lpoly, ref fit, ref lsqv,
                        input[i].doy, input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                        input[i].stl, input[i].f107a, input[i].f107, apar, mass, d, t);
                }
                errstr = errstr + input[i].doy + " " + input[i].sec + ".00 " +
                    input[i].alt + ".00 " +
                    input[i].glat + ".00 " + input[i].glong + ".00 " +
                    input[i].stl + ".00 " + input[i].f107a + ".00 " +
                    input[i].f107 + ".00\n";
                errstr = errstr + t[1] + " " + t[2] + " "
                    + d[1].ToString("E6") + " " + d[2].ToString("E6") + " " + d[3].ToString("E6") + " "
                    + d[4].ToString("E6") + " "
                    + d[5].ToString("E6") + " " + d[6].ToString("E6") + " " + d[7].ToString("E6") + " "
                    + d[8].ToString("E6") + "\n";
            }

            errstr = errstr + ("\n");

            System.IO.File.WriteAllText(directory + "tmsiscs.out", errstr);

            this.opsStatusL.Text = "done";
            Refresh();
        }


        int[] sv = new int[26] { 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
        int[] svt = new int[26] { 0, 1, 1, 1, 1, 1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };


        // local functions
        // test the routine against the nrl test cases

        // interface for plugin to hpop
        //void mainplugin(double[], double[], double );



        //using namespace System;


        public void mainplugin
               (
                 double[] r, double[] v, double jd
               )
        {
            int i, mass;
            MSISLib.msistype msis00r;
            double[] apar = new double[8];
            double f107, f107a, mfme;
            char interp, inputtype, fluxtype, f81type;

            double f107ctr81;
            double ap;
            double apavg;
            double kp;
            double sumkp;
            double[] kparr = new double[8];
            double[] aparr = new double[8];

            mass = 48;

            //	doy   =   172;
            //   	year  =  2005;
            //  	sec   = 29000;
            //	alt   =   400;
            //	glat  =    60;
            //	glong =   -70;
            //	stl   =     6;

            Int32 ktrActObs;
            string EOPupdate;
            string fileLoc;
            int year, mon, day, hr, minute, dat;
            double dut1, lod, xp, yp, ddpsi, ddeps, ddx, ddy;
            double f107bar, avgap;
            double jdFrac, second;
            jdFrac = 0.0;
            fileLoc = @"D:\Codes\LIBRARY\DataLib\nut80.dat";
            EOPSPWLibr.iau80in(fileLoc, out EOPSPWLibr.iau80arr);

            string eopFileName = @"D:\Codes\LIBRARY\DataLib\EOP-All-v1.1_2025-01-10.txt";
            EOPSPWLibr.readeop(ref EOPSPWLibr.eopdata, eopFileName, out ktrActObs, out EOPupdate);

            MathTimeLibr.invjday(jd, jdFrac, out year, out mon, out day, out hr, out minute, out second);
            //MathTimeLibr.jday(year, mon, day, hr, minute, second, out jd, out jdFrac);
            EOPSPWLibr.findeopparam(jd, jdFrac, 's', EOPSPWLibr.eopdata, out dut1, out dat,
               out lod, out xp, out yp, out ddpsi, out ddeps, out ddx, out ddy);


            string spwFileName = @"D:\Codes\LIBRARY\DataLib\SpaceWeather-All-v1.2_2025-01-10.txt";
            string errstr;
            EOPSPWLibr.readspw(ref EOPSPWLibr.spwdata, spwFileName, out ktrActObs, out errstr);

            MathTimeLibr.jday(year, mon, day, hr, minute, second, out jd, out jdFrac);
            //                              interp, adj obs, last ctr, actu const
            EOPSPWLibr.findspwparam(jd, jdFrac, 's', 'a', 'l', 'a', EOPSPWLibr.spwdata, out f107, out f107bar,
               out ap, out avgap, aparr, out kp, out sumkp, kparr);



            /* get space weather data values */
            for (i = 0; i < 8; i++)
                apar[i] = 4;
            f107a = 150;
            f107 = 150;

            //    interp      - interpolation       n-none, a-ap only, f-f10.7 only, b-both
            interp = 'n';
            //    inputtype   - input type          a-actual, u - user   c - constant
            inputtype = 'a'; // actual

            fluxtype = 'a';  // adjusted
            f81type = 'c';   // centered

            //        initeop( eoparr, jdeopsatrt);
            //		EopSpw::initspw(spwarr, "D:/Codes/LIBRARY/CPP/TestMSIS/Data/SpaceWeather-All-v1.2_11-03-2014.txt",  jdspwstart);
            EOPSPWLibr.readspw(ref EOPSPWLibr.spwdata, spwFileName, out ktrActObs, out errstr);

            jd = 2452363.5;
            r[0] = 1;
            r[1] = 2;
            r[2] = 3;

            // initialize the msis routine - probably do in the initialize part of the code for plug-in
            MSISLibr.msis00init(out msis00r);

            mfme = 0.0;
            //        EopSpw::findatmosparam( jd, mfme, interp, fluxtype, f81type, inputtype, spwarr, jdspwstart, f107, f107ctr81,
            //                        ap, apavg, aparray, kp, sumkp, kparray );

            //       interfaceatmos( jd, mfme, r, interp, fluxtype, f81type, inputtype, msis00r, spwarr, jdspwstart );
        }


        //     test case for running jb2006
        //     input to jb2006:
        //           amjd   : date and time, in modified julian days
        //                    and fraction (mjd =3d jd-2400000.5)
        //           sun(1) : right ascension of sun (radians)
        //           sun(2) : declination of sun (radians)
        //           sat(1) : right ascension of position (radians)
        //           sat(2) : geocentric latitude of position (radians)
        //           sat(3) : height of position (km)
        //           geo(1) : 10.7-cm solar flux (1.0e-22*watt/(m**2*hertz))
        //                    (tabular time 1.0 day earlier)
        //           geo(2) : 10.7-cm solar flux, ave.
        //                    81-day centered on the input time
        //           geo(3) : geomagnetic planetary 3-hour index
        //                    a-sub-p for a tabular time 0.279 days earlier
        //                    (6.7 hours earlier)
        //           s10    : euv index (26-34 nm) scaled to f10
        //                    (tabular time 1.0 day earlier)
        //           s10b   : euv 81-day ave. centered index
        //           xm10   : mg2 index scaled to f10
        //           xm10b  : mg2 81-day ave. centered index
        //                    (tabular time 5.0 days earlier)
        //
        //     output from jb2006:
        //           temp(1): exospheric temperature above input position (deg k)
        //           temp(2): temperature at input position (deg k)
        //           rho    : total mass-desnity at input position (kg/m**3)
        //

        public void testjb06()
        {
            //const double pi       = 3.1415926535879323846;
            double[] sun = new double[3];
            double[] sat = new double[4];
            double[] geo = new double[4];
            double[] temp1 = new double[4];
            //      double  pi = 3.1415927;
            double s10, s10b, f10, f10b, xm10, xm10b, ap, d1950, amjd, rho;
            int iyr, iyy, iday, id1950;

            //    set solar indices
            //     use 1 day lag for euv and f10 influence
            s10 = 140;
            s10b = 100;
            f10 = 135;
            f10b = 95;

            //     use 5 day lag for mg fuv influence
            xm10 = 130;
            xm10b = 95;

            //     use 6.7 hr lag for ap influence
            ap = 30;
            geo[1] = f10;
            geo[2] = f10b;
            geo[3] = ap;

            //     set point of interest location (radians and km)
            sat[1] = 90.0 * Math.PI / 180.0;
            sat[2] = 45.0 * Math.PI / 180.0;
            sat[3] = 400.0;

            //     set sun location (radians)
            sun[1] = 90.0 * Math.PI / 180.0;
            sun[2] = 20.0 * Math.PI / 180.0;

            //    set time of interest
            iyr = 1;
            iday = 200;
            if (iyr < 50) iyr = iyr + 100;
            iyy = (iyr - 50) * 365 + ((iyr - 1) / 4 - 12);
            id1950 = iyy + iday;
            d1950 = id1950;
            amjd = d1950 + 33281.0;

            //    compute density kg/m3 rho
            //      jb2006 (amjd,sun,sat,geo,s10,s10b,xm10,xm10b,temp1, rho);

            rho = 0.0;
            errstr = errstr + ("%8.0f  %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %10.1f %10.1f %14.8g  \n",
                     d1950, s10, s10b, f10, f10b, xm10, xm10b, ap, temp1[1], temp1[2], rho);
            //     output results:
            // 18828.  140. 100. 135.  95. 130.  95.  30.  1145.8  1137.7 0.4066d-11

        }  // testjb06

    } // class

}  // namespace

