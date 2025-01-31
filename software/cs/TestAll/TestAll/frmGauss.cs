using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text.RegularExpressions;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Windows.Forms;

using MathTimeMethods;     // Edirection, globals
using EOPSPWMethods;       // EOPDataClass, SPWDataClass, iau80Class, iau06Class
using AstroLibMethods;     // EOpt, gravityConst, astroConst, xysdataClass, jpldedataClass

namespace TestAllTool
{
    public partial class frmGauss : Form
    {
        public MathTimeLib MathTimeLibr = new MathTimeLib();

        public EOPSPWLib EOPSPWLibr = new EOPSPWLib();

        public AstroLib AstroLibr = new AstroLib();

        public StringBuilder strbuild = new StringBuilder();
        public StringBuilder strbuildall = new StringBuilder();
        public StringBuilder strbuildObs = new StringBuilder();
        public StringBuilder strbuildallsum = new StringBuilder();
        public string methodType;

        public List<string> fileNames = new List<string>();


        // ------------------------- sensor data
        // location, weights, biases, etc
        public class senClass
        {
            public int sennum;
            public string senname;
            public double senlat, senlon, senalt;
            public double rngmin, rngmax, azmin, azmax, elmin, elmax;
            public double biasrng, biasaz, biasel, biasdrng, biasdaz, biasdel,
                biastrtasc, biastdecl;
            public double noisex, noisey, noisez, noisexdot, noiseydot, noisezdot, noisebstar;
            public double noiserng, noiseaz, noiseel, noisedrng, noisedaz, noisedel,
                noisetrtasc, noisetdecl;
            public double[] rsecef = new double[3];
            public double[] vsecef = new double[3];
        }
        public senClass[] senClassArr = new senClass[100];


        // ------------------------- obs data
        // overall idea is to store eveyrthing as tracks
        // initial assignment from sensor is generally pretty good for short tracks, so keep it. 
        // trackcleaner does a good job of stitching together short tracks into slightly longer ones (minutes)
        // this program then takes those tracks and does processing and OD on them to get longer (multi-rev) tracks processed
        public class obsClass
        {
            public double jd;
            public double jdf;
            public double mjd;
            public int year;
            public int month;
            public int day;
            public int hr;
            public int minute;
            public double second;
            public string datestr;
            public Int32 sennum;
            public string initsscnum;
            public string sscnum;
            public int passnum;   // may not have this usually, but perhaps with gobs
            public double rtasc;
            public double decl;
            public double rtascdot;
            public double decldot;
            public double az;
            public double el;
            public double azdot;
            public double eldot;
            public double lat;
            public double lon;
            public double alt;
            public double[] rsecef = new double[3];
            public double[] rseci = new double[3];
            public double dut1;
            public int dat;
            public double lod;
            public double xp;
            public double yp;
            public double ddpsi;
            public double ddeps;
            public double ddx;
            public double ddy;
        };
        obsClass[] obsClassArr = new obsClass[50000];


        // ------------------------- track data
        // store all the tracks
        // uses the obs class for individual data
        public class trackClass
        {
            public Int32 trackktr;
            public Int32 obsintrk;
            public double passlen;
            public Int32 sscnum;
            public obsClass[] tracks = new obsClass[50000];
        };
        trackClass[] trackClassArr = new trackClass[5000];



        public frmGauss()
        {
            InitializeComponent();

            if (this.cbGaussList.SelectedItem == null)
                this.cbGaussList.SelectedIndex = 4;
            if (this.comboBoxSen1.SelectedItem == null)
                this.comboBoxSen1.SelectedIndex = 2;
            if (this.comboBoxSen2.SelectedItem == null)
                this.comboBoxSen2.SelectedIndex = 1;
            if (this.comboBoxSen3.SelectedItem == null)
                this.comboBoxSen3.SelectedIndex = 0;

            // save values from time to time
            // http://www.dotnetperls.com/settings
            this.trtasc1.Text = TestAll.Properties.Settings.Default.trtasc1;
            this.tdecl1.Text = TestAll.Properties.Settings.Default.tdecl1;
            this.trtasc2.Text = TestAll.Properties.Settings.Default.trtasc2;
            this.tdecl2.Text = TestAll.Properties.Settings.Default.tdecl2;
            this.trtasc3.Text = TestAll.Properties.Settings.Default.trtasc3;
            this.tdecl3.Text = TestAll.Properties.Settings.Default.tdecl3;
            this.time1.Text = TestAll.Properties.Settings.Default.time1;
            this.time2.Text = TestAll.Properties.Settings.Default.time2;
            this.time3.Text = TestAll.Properties.Settings.Default.time3;
            this.lat1.Text = TestAll.Properties.Settings.Default.lat1;
            this.lon1.Text = TestAll.Properties.Settings.Default.lon1;
            this.alt1.Text = TestAll.Properties.Settings.Default.alt1;
            this.lat2.Text = TestAll.Properties.Settings.Default.lat2;
            this.lon2.Text = TestAll.Properties.Settings.Default.lon2;
            this.alt2.Text = TestAll.Properties.Settings.Default.alt2;
            this.lat3.Text = TestAll.Properties.Settings.Default.lat3;
            this.lon3.Text = TestAll.Properties.Settings.Default.lon3;
            this.alt3.Text = TestAll.Properties.Settings.Default.alt3;
        }


        // ------------------------------------------------------------------------------------------------------
        // collect sensor data for JTMA sites (and other if added)
        // use sensor name
        public void getSensorECI
            (
                string ecistr, string timestr, out double latgd, out double lon, out double alt
            )
        {
            double[] rsecef = new double[3];
            double[] vsecef = new double[3];
            double[] asecef = new double[3];
            double[] rseci = new double[3];
            double[] vseci = new double[3];
            double[] aseci = new double[3];
            double magr, latgc;
            double jd, jdf, second;
            double dut1, lod, xp, yp, ddpsi, ddeps, ddx, ddy, jdtt, jdftt, ttt, jdut1;
            int dat;
            int year, month, day, hr, minute;

            string tmpstr1 = Regex.Replace(ecistr, @"\s+", " ");
            string tmpstr = Regex.Replace(tmpstr1, ",", "");
            string[] linesplt = tmpstr.Split(' ');
            rseci[0] = Convert.ToDouble(linesplt[0]);
            rseci[1] = Convert.ToDouble(linesplt[1]);
            rseci[2] = Convert.ToDouble(linesplt[2]);
            magr = MathTimeLibr.mag(rseci);

            vseci[0] = 0.0;
            vseci[1] = 0.0;
            vseci[2] = 0.0;

            MathTimeLibr.STKtime2JD(timestr, out jd, out jdf);
            MathTimeLibr.invjday(jd, jdf, out year, out month, out day, out hr, out minute, out second);

            string nutLoc;
            Int32 ktrActObs;
            string EOPupdate;
            EOPSPWLib.iau80Class iau80arr;
            nutLoc = @"D:\Codes\LIBRARY\DataLib\nut80.dat";
            EOPSPWLibr.iau80in(nutLoc, out iau80arr);

            string eopFileName = @"D:\Codes\LIBRARY\DataLib\EOP-All-v1.1_2023-01-01.txt";
            EOPSPWLibr.readeop(ref EOPSPWLibr.eopdata, eopFileName, out ktrActObs, out EOPupdate);

            EOPSPWLibr.findeopparam(jd, jdf, 's', EOPSPWLibr.eopdata, 
                       out dut1, out dat, out lod, out xp, out yp, out ddpsi, out ddeps, out ddx, out ddy);
            jdtt = jd;
            jdftt = jdf + (dat + 32.184) / 86400.0;
            ttt = (jdtt + jdftt - 2451545.0) / 36525.0;
            // note you have to use tdb for time of interst AND j2000 (when dat = 32)
            //  ttt = (jd + jdFrac + (dat + 32.184) / 86400.0 - 2451545.0 - (32 + 32.184) / 86400.0) / 36525.0;
            jdut1 = jd + jdf + dut1 / 86400.0;

            // -------- convert r to ecef for lat/lon calculation
            AstroLibr.eci_ecef(ref rseci, ref vseci, ref aseci, MathTimeLib.Edirection.eto, ref rsecef, ref vsecef, ref asecef,
                EOPSPWLibr.iau80arr, ttt, jdut1, lod, xp, yp, ddpsi, ddeps);

            AstroLibr.ecef2ll(rsecef, out latgc, out latgd, out lon, out alt);
        }  // getSensorECI

        // ------------------------------------------------------------------------------------------------------
        // collect sensor data for JTMA sites (and other if added)
        // use sensor name
        public void getSensorLLA
            (
                int obsktr, string senname, out int sennum, ref double lat, ref double lon, ref double alt
            )
        {
            double rad = 180.0 / Math.PI;

            sennum = 0;
            switch (senname)
            {
                case "JTMA-001":
                    lat = 35.305774 / rad;
                    lon = -106.457277 / rad;
                    alt = 1.7882;
                    sennum = 6201;
                    break;
                case "JTMA-003":
                    lat = 35.305782 / rad;
                    lon = -106.457083 / rad;
                    alt = 1.7882;
                    sennum = 6203;
                    break;
                case "JTMA-004":
                    lat = 35.305715 / rad;
                    lon = -106.457170 / rad;
                    alt = 1.7971;
                    sennum = 6204;
                    break;
                case "JTMA-005":
                    lat = 40.037885 / rad;
                    lon = -75.596769 / rad;
                    alt = 0.078;
                    sennum = 6205;
                    break;
                case "JTMA-006":
                    lat = -32.007550 / rad;
                    lon = 116.136740 / rad;
                    alt = 0.374;
                    sennum = 6206;
                    break;
                case "JTMA-007":
                    lat = -32.294010 / rad;
                    lon = 148.585270 / rad;
                    alt = 0.280;
                    sennum = 6207;
                    break;
                case "JTMA-008":
                    lat = 38.165700 / rad;
                    lon = -2.326862 / rad;
                    alt = 1.654;
                    sennum = 6208;
                    break;
                case "JTMA-009":
                    lat = -25.743746 / rad;
                    lon = 29.298861 / rad;
                    alt = 1.480;
                    sennum = 6209;
                    break;
                case "JTMA-010":
                    lat = 30.597332 / rad;
                    lon = 34.763020 / rad;
                    alt = 0.872;
                    sennum = 6210;
                    break;
                case "JTMA-011":
                    lat = 20.707100 / rad;
                    lon = -156.257467 / rad;
                    alt = 3.0703;
                    sennum = 6211;
                    break;
                case "JTMA-012":
                    lat = -22.953317 / rad;
                    lon = -68.179997 / rad;
                    alt = 2.660;
                    sennum = 6212;
                    break;
                case "JTMA-013":
                    lat = -31.820550 / rad;
                    lon = 117.281500 / rad;
                    alt = 0.331933;
                    sennum = 6213;
                    break;
                case "JTMA-014":
                    lat = 34.291194 / rad;
                    lon = -115.382611 / rad;
                    alt = 0.884987;
                    sennum = 6214;
                    break;
                case "JTMA-015":
                    lat = -23.771000 / rad;
                    lon = 133.919000 / rad;
                    alt = 0.884987;
                    sennum = 6215;
                    break;
                case "JTMA-016":
                    lat = 20.900587 / rad;
                    lon = -156.504102 / rad;
                    alt = 0.227;
                    sennum = 6216;
                    break;
                case "JTMA-017":
                    lat = 30.236 / rad;
                    lon = -5.609 / rad;
                    alt = 0.674;
                    sennum = 6217;
                    break;
                case "Other":
                    // don't re-assign the lat lon values, leave as input
                    break;
            }  // switch


            lat = Convert.ToDouble(lat) / rad;      //rad
            lon = Convert.ToDouble(lon) / rad;      //rad
            alt = Convert.ToDouble(alt);            // km

            senClassArr[0] = new senClass();

            // -----------  Find ecef site pos/vel vector 
            AstroLibr.site(lat, lon, alt, out senClassArr[0].rsecef, out senClassArr[0].vsecef);

            senClassArr[0].senlat = lat;
            senClassArr[0].senlon = lon;
            senClassArr[0].senalt = alt;
            senClassArr[0].biastrtasc = 0.0000048482;
            senClassArr[0].biastdecl = 0.0000048482;
            senClassArr[0].noisetrtasc = 0.0000048482;  // 1 arcsec = 0.0002777 deg = 0.0000048482 rad
            senClassArr[0].noisetdecl = 0.0000048482;
            senClassArr[0].senname = senname;
        }  // getSensorLLA


        // sennum comes in with 2 left chars?
        public void writeJTMTrack(trackClass[] trackClassArr, Int32 trknumktr, out string tempstr, out string tempsumstr)
        {
            string senstr;
            //StringBuilder strbuild = new StringBuilder();
            double rad = 180.0 / Math.PI;

            tempstr = "";
            senstr = "JTMA-0" + (trackClassArr[trknumktr].tracks[1].sennum).ToString();
            //senstr = "JTMA-0" + (trackClassArr[trknumktr].tracks[1].sennum).ToString().Remove(0, 2);
            //matlab:
            //senstr = (trackClassArr[trknumktr].tracks[1].sennum).ToString().Remove(0, 2);

            for (Int32 ktrtk = 1; ktrtk <= trackClassArr[trknumktr].obsintrk; ktrtk++)
            {
                // jtm format 3,2021-11-15T00:10:50.951000,JTMA-017,31306,0,1,79.5186996,-4.8842573,10.862204
                tempstr = tempstr + "3," + trackClassArr[trknumktr].tracks[ktrtk].year + "-" + trackClassArr[trknumktr].tracks[ktrtk].month
                    + "-" + trackClassArr[trknumktr].tracks[ktrtk].day + "T" + trackClassArr[trknumktr].tracks[ktrtk].hr + ":"
                    + trackClassArr[trknumktr].tracks[ktrtk].minute + ":" + trackClassArr[trknumktr].tracks[ktrtk].second.ToString("00.######") + ","
                    + senstr + "," + trackClassArr[trknumktr].tracks[ktrtk].sscnum + ", 0,1,"
                    + (trackClassArr[trknumktr].tracks[ktrtk].rtasc * rad) + "," + (trackClassArr[trknumktr].tracks[ktrtk].decl * rad)
                    + ", 10.00 \n";

                // use for matlab testing - track number instead of 3
                //tempstr = tempstr + trknumktr + "," + trackClassArr[trknumktr].tracks[ktrtk].year + "," + trackClassArr[trknumktr].tracks[ktrtk].month
                //    + "," + trackClassArr[trknumktr].tracks[ktrtk].day + "," + trackClassArr[trknumktr].tracks[ktrtk].hr + ","
                //    + trackClassArr[trknumktr].tracks[ktrtk].minute + "," + trackClassArr[trknumktr].tracks[ktrtk].second.ToString("00.######") + ","
                //    + senstr + "," + trackClassArr[trknumktr].tracks[ktrtk].sscnum + ", 0,1,"
                //    + (trackClassArr[trknumktr].tracks[ktrtk].rtasc * rad) + "," + (trackClassArr[trknumktr].tracks[ktrtk].decl * rad)
                //    + ", 10.00 \n";
            }

            // write out summary for each track
            Int32 lastktr = trackClassArr[trknumktr].obsintrk;
            tempsumstr = trackClassArr[trknumktr].tracks[1].year + "-" + trackClassArr[trknumktr].tracks[1].month
                    + "-" + trackClassArr[trknumktr].tracks[1].day + "T" + trackClassArr[trknumktr].tracks[1].hr + ":"
                    + trackClassArr[trknumktr].tracks[1].minute + ":" + trackClassArr[trknumktr].tracks[1].second.ToString("00.######")
                    + "," + senstr + "," + trackClassArr[trknumktr].tracks[1].sscnum + "\n"
                    + trackClassArr[trknumktr].tracks[lastktr].year + "-" + trackClassArr[trknumktr].tracks[lastktr].month
                    + "-" + trackClassArr[trknumktr].tracks[lastktr].day + "T" + trackClassArr[trknumktr].tracks[lastktr].hr + ":"
                    + trackClassArr[trknumktr].tracks[lastktr].minute + ":" + trackClassArr[trknumktr].tracks[lastktr].second.ToString("00.######")
                    + "," + senstr + "," + trackClassArr[trknumktr].tracks[lastktr].sscnum + "," + lastktr;

        }  // writeJTMTrack



        // doall the tests
        // may be a single case, or all of them
        // input data is stored in the incoming list
        public void doangletests
            (
              double[] jd, double[] jdf,
              double[] trtasc, double[] tdecl,
              double[] latgd, double[] lon, double[] alt,
              double[] rseci1, double[] vseci1,
              double[] rseci2, double[] vseci2,
              double[] rseci3, double[] vseci3,
              int idx1, int idx2, int idx3,
              string ans
            )
        {
            double[] r2 = new double[3];
            double[] v2 = new double[3];
            double n, af, ag, chi, psi, meanlonM, meanlonNu;
            double[] h = new double[3];
            Int16 fr;
            double rad;
            double[] initguess = new double[30];
            string errstr;
            string coestr, eqstr, posstr, velstr, rsstr, timestr;
            double rng1, rng2, rng3;
            Int32 iyear1, imon1, iday1, ihr1, iminute1;
            Int32 iyear2, imon2, iday2, ihr2, iminute2;
            Int32 iyear3, imon3, iday3, ihr3, iminute3;
            double isecond1, isecond2, isecond3, bigr2x;
            Int32 numhalfrev;
            //conv = Math.PI / (180.0 * 3600.0);
            rad = 180.0 / Math.PI;
            errstr = "";
            char diffsites = 'n';
            string[] monstr = new string[13] { "", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC" };

            this.trtasc1.Text = (trtasc[idx1] * rad).ToString();
            this.tdecl1.Text = (tdecl[idx1] * rad).ToString();
            this.trtasc2.Text = (trtasc[idx2] * rad).ToString();
            this.tdecl2.Text = (tdecl[idx2] * rad).ToString();
            this.trtasc3.Text = (trtasc[idx3] * rad).ToString();
            this.tdecl3.Text = (tdecl[idx3] * rad).ToString();
            this.lat1.Text = (latgd[idx1] * rad).ToString();
            this.lon1.Text = (lon[idx1] * rad).ToString();
            this.alt1.Text = alt[idx1].ToString();
            this.lat2.Text = (latgd[idx2] * rad).ToString();
            this.lon2.Text = (lon[idx2] * rad).ToString();
            this.alt2.Text = alt[idx2].ToString();
            this.lat3.Text = (latgd[idx3] * rad).ToString();
            this.lon3.Text = (lon[idx3] * rad).ToString();
            this.alt3.Text = alt[idx3].ToString();

            // gooding tests cases from Gooding paper (1997 CMDA)
            double[] los1;
            double[] los2;
            double[] los3;

            los1 = new double[3];
            los2 = new double[3];
            los3 = new double[3];
            los1[0] = Math.Cos(tdecl[idx1]) * Math.Cos(trtasc[idx1]);
            los1[1] = Math.Cos(tdecl[idx1]) * Math.Sin(trtasc[idx1]);
            los1[2] = Math.Sin(tdecl[idx1]);

            los2[0] = Math.Cos(tdecl[idx2]) * Math.Cos(trtasc[idx2]);
            los2[1] = Math.Cos(tdecl[idx2]) * Math.Sin(trtasc[idx2]);
            los2[2] = Math.Sin(tdecl[idx2]);

            los3[0] = Math.Cos(tdecl[idx3]) * Math.Cos(trtasc[idx3]);
            los3[1] = Math.Cos(tdecl[idx3]) * Math.Sin(trtasc[idx3]);
            los3[2] = Math.Sin(tdecl[idx3]);

            strbuildall.AppendLine("los1 " + los1[0].ToString("0.00000000") + " " +
                los1[1].ToString("0.00000000") + " " + los1[2].ToString("0.00000000") +
                " " + MathTimeLibr.mag(los1).ToString("0.00000000"));
            strbuildall.AppendLine("los2 " + los2[0].ToString("0.00000000") + " " +
                los2[1].ToString("0.00000000") + " " + los2[2].ToString("0.00000000") +
                " " + MathTimeLibr.mag(los2).ToString("0.00000000"));
            strbuildall.AppendLine("los3 " + los3[0].ToString("0.00000000") + " " +
                los3[1].ToString("0.00000000") + " " + los3[2].ToString("0.000000") +
                " " + MathTimeLibr.mag(los3).ToString("0.00000000"));

            double p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper;

            // to get initial guess, take measurments (1/2 and 2/3), assume circular orbit
            // find velocity and compare - just distinguish between LEO, GPS and GEO for now
            double dtrtasc1, dtdecl1, dtrtasc2, dtdecl2, dt1, dt2;
            dt1 = (jd[idx2] - jd[idx1]) * 86400.0 + (jdf[idx2] - jdf[idx1]) * 86400.0;
            dt2 = (jd[idx3] - jd[idx2]) * 86400.0 + (jdf[idx3] - jdf[idx2]) * 86400.0;
            dtrtasc1 = (trtasc[idx2] - trtasc[idx1]) / dt1;
            dtrtasc2 = (trtasc[idx3] - trtasc[idx2]) / dt2;
            dtdecl1 = (tdecl[idx2] - tdecl[idx1]) / dt1;
            dtdecl2 = (tdecl[idx3] - tdecl[idx2]) / dt2;

            strbuildall.AppendLine("rtasc " + (trtasc[idx1] * rad).ToString() + " " + (trtasc[idx2] * rad).ToString()
                + " " + (trtasc[idx3] * rad).ToString());
            strbuildall.AppendLine("decl " + (tdecl[idx1] * rad).ToString() + " " + (tdecl[idx2] * rad).ToString()
                + " " + (tdecl[idx3] * rad).ToString());


            strbuildall.AppendLine("");
            strbuildallsum.AppendLine("Laplace -----------------------------------");
            strbuildall.AppendLine("Laplace -----------------------------------");

            rsstr = rseci1[0].ToString() + " " + rseci1[1].ToString() + " " + rseci1[2].ToString();
            this.rseci1t.Text = rsstr;
            strbuildall.AppendLine("\n\ninputs: \n");
            strbuildall.AppendLine("Site obs1 "
                + this.rseci1t.Text + " km  lat " + (latgd[idx1] * rad).ToString() + " lon " + (lon[idx1] * rad).ToString() + " "
                + alt[idx1].ToString());
            rsstr = rseci2[0].ToString() + " " + rseci2[1].ToString() + " " + rseci2[2].ToString();
            this.rseci2t.Text = rsstr;
            strbuildall.AppendLine("Site obs2 "
                + this.rseci2t.Text + " km  lat " + (latgd[idx2] * rad).ToString() + " lon " + (lon[idx2] * rad).ToString() + " "
                + alt[idx2].ToString());
            rsstr = rseci3[0].ToString() + " " + rseci3[1].ToString() + " " + rseci3[2].ToString();
            this.rseci3t.Text = rsstr;
            strbuildall.AppendLine("Site obs3 "
                + this.rseci3t.Text + " km  lat " + (latgd[idx3] * rad).ToString() + " lon " + (lon[idx3] * rad).ToString() + " "
                + alt[idx3].ToString());

            MathTimeLibr.invjday(jd[idx1], jdf[idx1], out iyear1, out imon1, out iday1, out ihr1, out iminute1, out isecond1);
            MathTimeLibr.JD2STKtime(jd[idx1], jdf[idx1], out timestr);
            this.time1.Text = timestr;
            strbuildall.AppendLine("obs#1 " + iyear1.ToString() + " " + imon1.ToString() + " " + iday1.ToString()
                + " " + ihr1.ToString("00") + " " + iminute1.ToString("00") + " " + isecond1.ToString("0.000")
                + " " + (trtasc[idx1] * rad).ToString() + " " + (tdecl[idx1] * rad).ToString().ToString());
            MathTimeLibr.invjday(jd[idx2], jdf[idx2], out iyear2, out imon2, out iday2, out ihr2, out iminute2, out isecond2);
            MathTimeLibr.JD2STKtime(jd[idx2], jdf[idx2], out timestr);
            this.time2.Text = timestr;
            strbuildall.AppendLine("obs#2 " + iyear2.ToString() + " " + imon2.ToString() + " " + iday2.ToString()
                + " " + ihr2.ToString("00") + " " + iminute2.ToString("00") + " " + isecond2.ToString("0.000")
                + " " + (trtasc[idx2] * rad).ToString() + " " + (tdecl[idx2] * rad).ToString().ToString());
            MathTimeLibr.invjday(jd[idx3], jdf[idx3], out iyear3, out imon3, out iday3, out ihr3, out iminute3, out isecond3);
            MathTimeLibr.JD2STKtime(jd[idx3], jdf[idx3], out timestr);
            this.time3.Text = timestr;
            strbuildall.AppendLine("obs#3 " + iyear3.ToString() + " " + imon3.ToString() + " " + iday3.ToString()
                + " " + ihr3.ToString("00") + " " + iminute3.ToString("00") + " " + isecond3.ToString("0.000")
                + " " + (trtasc[idx3] * rad).ToString() + " " + (tdecl[idx3] * rad).ToString().ToString());
            //if (caseopt == 2)
            //    diffsites = 'y';
            //else 
            //diffsites = 'n';

            MathTimeLibr.invjday(jd[idx1], jdf[idx1], out iyear1, out imon1, out iday1, out ihr1, out iminute1, out isecond1);
            MathTimeLibr.invjday(jd[idx2], jdf[idx2], out iyear2, out imon2, out iday2, out ihr2, out iminute2, out isecond2);
            MathTimeLibr.invjday(jd[idx3], jdf[idx3], out iyear3, out imon3, out iday3, out ihr3, out iminute3, out isecond3);

            // ------------- build tests for Der approach -------------
            strbuildObs.AppendLine("&A3dav_input\n   indata = 0, "
                + latgd[idx1] * rad + "," + lon[idx1] * rad + "," + alt[idx1] + ","
                + iyear1 + "," + imon1 + "," + iday1 + "," + ihr1 + "," + iminute1 + "," + isecond1 + ","
                + trtasc[idx1] * rad + "," + tdecl[idx1] * rad + ",\n"
                + latgd[idx2] * rad + "," + lon[idx2] * rad + "," + alt[idx2] + ","
                + iyear2 + "," + imon2 + "," + iday2 + "," + ihr2 + "," + iminute2 + "," + isecond2 + ","
                + trtasc[idx2] * rad + "," + tdecl[idx2] * rad + ",\n"
                + latgd[idx3] * rad + "," + lon[idx3] * rad + "," + alt[idx3] + ","
                + iyear3 + "," + imon3 + "," + iday3 + "," + ihr3 + "," + iminute3 + "," + isecond3 + ","
                + trtasc[idx3] * rad + "," + tdecl[idx3] * rad + ",\n/"
                );
            // alternate der format
            int dayofyr1, dayofyr2, dayofyr3, hm, mm;
            double s;
            string strrtasc1, strrtasc2, strrtasc3;
            string strrs1, strrs2, strrs3;
            hm = mm = 0;
            s = 0.0;
            MathTimeLibr.findDays(iyear1, imon1, iday1,
                ihr1, iminute1, isecond1, out dayofyr1);
            MathTimeLibr.hms_rad(ref hm, ref mm, ref s, MathTimeLib.Edirection.efrom, ref trtasc[idx1]);
            strrtasc1 = hm + " " + mm + " " + s;
            double inval = latgd[idx1];
            MathTimeLibr.dms_rad(ref hm, ref mm, ref s, MathTimeLib.Edirection.efrom, ref inval);
            strrs1 = hm + " " + mm + " " + s;
            inval = lon[idx1];
            MathTimeLibr.dms_rad(ref hm, ref mm, ref s, MathTimeLib.Edirection.efrom, ref inval);
            strrs1 = strrs1 + " " + hm + " " + mm + " " + s + " " + alt[idx1];

            MathTimeLibr.findDays(iyear2, imon2, iday2,
                ihr2, iminute2, isecond2, out dayofyr2);
            MathTimeLibr.hms_rad(ref hm, ref mm, ref s, MathTimeLib.Edirection.efrom, ref trtasc[idx2]);
            strrtasc2 = hm + " " + mm + " " + s;
            inval = latgd[idx2];
            MathTimeLibr.dms_rad(ref hm, ref mm, ref s, MathTimeLib.Edirection.efrom, ref inval);
            strrs2 = hm + " " + mm + " " + s;
            inval = lon[idx2];
            MathTimeLibr.dms_rad(ref hm, ref mm, ref s, MathTimeLib.Edirection.efrom, ref inval);
            strrs2 = strrs2 + " " + hm + " " + mm + " " + s + " " + alt[idx2];

            MathTimeLibr.findDays(iyear3, imon3, iday3,
                ihr3, iminute3, isecond3, out dayofyr3);
            MathTimeLibr.hms_rad(ref hm, ref mm, ref s, MathTimeLib.Edirection.efrom, ref trtasc[idx3]);
            strrtasc3 = hm + " " + mm + " " + s;
            inval = latgd[idx3];
            MathTimeLibr.dms_rad(ref hm, ref mm, ref s, MathTimeLib.Edirection.efrom, ref inval);
            strrs3 = hm + " " + mm + " " + s;
            inval = lon[idx3];
            MathTimeLibr.dms_rad(ref hm, ref mm, ref s, MathTimeLib.Edirection.efrom, ref inval);
            strrs3 = strrs3 + " " + hm + " " + mm + " " + s + " " + alt[idx3];

            strbuildObs.AppendLine("1 lat lon info " + latgd[idx1] * rad + "," + lon[idx1] * rad + ","
                + alt[idx1] + " " + strrs1 + "\n"
                + "37775 " + " 10  "
                + (iyear1 - 2000) + " " + dayofyr1 + " "
                + ihr1 + " " + iminute1 + " " + isecond1 + " "
                + tdecl[idx1] * rad + " " + strrtasc1 + "\n"
                + "2 lat lon info " + latgd[idx2] * rad + "," + lon[idx2] * rad + ","
                + alt[idx2] + " " + strrs2 + "\n"
                + "37775 " + " 10  "
                + (iyear2 - 2000) + " " + dayofyr2 + " "
                + ihr2 + " " + iminute2 + " " + isecond2 + " "
                + tdecl[idx2] * rad + " " + strrtasc2 + "\n"
                + "3 lat lon info " + latgd[idx3] * rad + "," + lon[idx3] * rad + ","
                + alt[idx3] + " " + strrs3 + "\n"
                + "37775 " + " 08  "
                + (iyear3 - 2000) + " " + dayofyr3 + " "
                + ihr3 + " " + iminute3 + " " + isecond3 + " "
                + tdecl[idx3] * rad + " " + strrtasc3 + "\n");
            // ------------- end build tests for Der approach -------------


            // ------------- write out the 3 obs tracks in jtm format to avoid geosc weirdness -------------
            strbuildObs.AppendLine("3," + iyear1 + "-" + imon1
                + "-" + iday1 + "T" + ihr1.ToString("00") + ":" + iminute1.ToString("00") + ":" + isecond1.ToString("00.######") + ","
                + "JTMA-010" + "," + "37775" + ", 0,1,"
                + (trtasc[idx1] * rad) + "," + (tdecl[idx1] * rad)
                + ",10.00 \n");
            strbuildObs.AppendLine("3," + iyear2 + "-" + imon2
                + "-" + iday2 + "T" + ihr2.ToString("00") + ":" + iminute2.ToString("00") + ":" + isecond2.ToString("00.######") + ","
                + "JTMA-010" + "," + "37775" + ", 0,2,"
                + (trtasc[idx2] * rad) + "," + (tdecl[idx2] * rad)
                + ",20.00 \n");
            strbuildObs.AppendLine("3," + iyear3 + "-" + imon3
                + "-" + iday3 + "T" + ihr3.ToString("00") + ":" + iminute3.ToString("00") + ":" + isecond3.ToString("00.######") + ","
                + "JTMA-010" + "," + "37775" + ", 0,3,"
                + (trtasc[idx3] * rad) + "," + (tdecl[idx3] * rad)
                + ",30.00 \n");
            //string tempstr1, tempstr2, tempstr3;
            //writeJTMTrack(trackClassArr, currTrack1, out tempstr1, out tempsumstr);
            //writeJTMTrack(trackClassArr, currTrack2, out tempstr2, out tempsumstr);
            //writeJTMTrack(trackClassArr, currTrack3, out tempstr3, out tempsumstr);
            //tempstr = tempstr1 + tempstr2 + tempstr3;
            //File.WriteAllText(@"D:\datafile\CSSI\COMSPOC\UCT\Nov21test\temp\tr"
            //    + currTrack1 + "-" + currTrack2 + "-" + currTrack3 + "_" + strnum + ".jtm", tempstr);
            // ------------- end of write -------------


            if (methodType == "Laplace" || methodType == "All")
            {
                AstroLibr.angleslaplace(tdecl[idx1], tdecl[idx2], tdecl[idx3], trtasc[idx1], trtasc[idx2], trtasc[idx3],
                            jd[idx1], jdf[idx1], jd[idx2], jdf[idx2], jd[idx3], jdf[idx3],
                            diffsites, rseci1, rseci2, rseci3, out r2, out v2, out bigr2x, out errstr);
                //this.bigrL.Text = errstr.Substring(errstr[errstr.IndexOf("root pick")+10], 10);
                this.bigrL.Text = bigr2x.ToString();
                strbuildall.AppendLine(errstr);
                posstr = "r2 " + r2[0].ToString("0.000000") + " " +
                    r2[1].ToString("0.000000") + " " + r2[2].ToString("0.000000");
                velstr = "v2 " + v2[0].ToString("0.000000") + " " +
                    v2[1].ToString("0.000000") + " " + v2[2].ToString("0.000000");
                this.rL.Text = posstr;
                this.vL.Text = velstr;
                strbuildall.AppendLine(posstr + velstr);
                AstroLibr.rv2coe(r2, v2, out p, out a, out ecc, out incl, out raan, out argp, out nu, out m,
                    out arglat, out truelon, out lonper);
                coestr = "a= " + a.ToString("0.0000") + " e= " + ecc.ToString("0.000000000") + " i= " + (incl * rad).ToString("0.0000")
                    + "  " + (raan * rad).ToString("0.0000") + "  " + (argp * rad).ToString("0.0000") + "  "
                    + (nu * rad).ToString("0.0000") + "  " + (m * rad).ToString("0.0000") + "  " + (arglat * rad).ToString("0.0000");
                // + " " + (truelon * rad).ToString("0.0000") + " " + (lonper * rad).ToString("0.0000"));
                this.coeL.Text = coestr;
                strbuildall.AppendLine("\nLaplace coes " + coestr);
                strbuildall.AppendLine(ans);
                strbuildallsum.AppendLine("               " + coestr);
                strbuildallsum.AppendLine(ans);
                AstroLibr.rv2eq(r2, v2, out a, out n, out af, out ag, out chi, out psi, out meanlonM, out meanlonNu, out fr);
                eqstr = "n= " + n.ToString("0.0000") + " af= " + af.ToString("0.000000000") + " ag= " + ag.ToString("0.0000")
                    + " chi= " + chi.ToString("0.0000") + " psi= " + psi.ToString("0.0000") + "  " + (meanlonNu * rad).ToString("0.0000")
                    + "  " + (meanlonM * rad).ToString("0.0000");
                this.eqL.Text = eqstr;
                MathTimeLibr.cross(r2, v2, out h);
                this.hL.Text = h[0].ToString("0.0000") + "  " + h[1].ToString("0.0000") + "  " + h[2].ToString("0.0000");
            }

            if (methodType == "Gauss" || methodType == "All")
            {
                strbuildallsum.AppendLine("Gauss  -----------------------------------");
                strbuildall.AppendLine("Gauss  -----------------------------------");

                AstroLibr.anglesgauss(tdecl[idx1], tdecl[idx2], tdecl[idx3], trtasc[idx1], trtasc[idx2], trtasc[idx3],
                     jd[idx1], jdf[idx1], jd[idx2], jdf[idx2], jd[idx3], jdf[idx3],
                     rseci1, rseci2, rseci3, out r2, out v2, out errstr);
                strbuildall.AppendLine(errstr);
                posstr = "r2 " + r2[0].ToString("0.000000") + " " +
                    r2[1].ToString("0.000000") + " " + r2[2].ToString("0.000000");
                velstr = "v2 " + v2[0].ToString("0.000000") + " " +
                    v2[1].ToString("0.000000") + " " + v2[2].ToString("0.000000");
                this.rG.Text = posstr;
                this.vG.Text = velstr;
                strbuildall.AppendLine(posstr + velstr);
                AstroLibr.rv2coe(r2, v2, out p, out a, out ecc, out incl, out raan, out argp, out nu, out m,
                    out arglat, out truelon, out lonper);
                coestr = "a= " + a.ToString("0.0000") + " e= " + ecc.ToString("0.000000000") + " i= " + (incl * rad).ToString("0.0000")
                    + "  " + (raan * rad).ToString("0.0000") + "  " + (argp * rad).ToString("0.0000") + "  "
                    + (nu * rad).ToString("0.0000") + "  " + (m * rad).ToString("0.0000") + "  " + (arglat * rad).ToString("0.0000");
                // + " " + (truelon * rad).ToString("0.0000") + " " + (lonper * rad).ToString("0.0000"));
                this.coeG.Text = coestr;
                strbuildall.AppendLine("Gauss coes " + coestr);
                strbuildall.AppendLine(ans);
                strbuildallsum.AppendLine("               " + coestr);
                strbuildallsum.AppendLine(ans);
                AstroLibr.rv2eq(r2, v2, out a, out n, out af, out ag, out chi, out psi, out meanlonM, out meanlonNu, out fr);
                eqstr = "n= " + n.ToString("0.0000") + " af= " + af.ToString("0.000000000") + " ag= " + ag.ToString("0.0000")
                    + " chi= " + chi.ToString("0.0000") + " psi= " + psi.ToString("0.0000") + "  " + (meanlonNu * rad).ToString("0.0000")
                    + "  " + (meanlonM * rad).ToString("0.0000");
                this.eqG.Text = eqstr;
                MathTimeLibr.cross(r2, v2, out h);
                this.hG.Text = h[0].ToString("0.0000") + "  " + h[1].ToString("0.0000") + "  " + h[2].ToString("0.0000");
            }

            int n12, n13, n23;
            n12 = n13 = n23 = 0;

            double pctchg = 0.05;  // 5% seems reasonable
            //changing the pctchg doens't seem to vary the results very much
            // BUT what about changing as the processing goes through - coarser to finer?????????????????????????????
            //for (int ii = 0; ii < 10; ii++)
            //{
            // didn't seem to work either - converges pretty quickly once it goes
            //pctchg = 0.005 * ii;
            strbuildallsum.AppendLine("Double-r -------------------------- pctchg -------------------------" + pctchg);
            strbuildall.AppendLine("Double-r --------------------------- pctchg -----------------------" + pctchg);

            if (methodType == "Double-r" || methodType == "All")
            {
                strbuildall.AppendLine("Double-r -----------------------------------");
                // initial guesses needed for double-r and Gooding
                // use result from Gauss as it's usually pretty good
                // this seems to really help Gooding!!
                AstroLibr.getGaussRoot(tdecl[idx1], tdecl[idx2], tdecl[idx3], trtasc[idx1], trtasc[idx2], trtasc[idx3],
                     jd[idx1], jdf[idx1], jd[idx2], jdf[idx2], jd[idx3], jdf[idx3],
                     rseci1, rseci2, rseci3, out bigr2x);
                this.bigrG.Text = bigr2x.ToString();
                // initguess[caseopt] = bigr2x;
                // bigr2x = AstroLibr.gravConst.re * bigr2x;



                #region MoultonApproach Not work yet
                // Moulton approach - after Battin (140)
                //double coef, xn, M, ax, sp3, sin_psi, cos_psi, magr2, magr23, d, d1c, d2c, s1, s2, s3, s4, s5, s6,
                //    tau12, tau13, tau32;
                //double[,] dmat = new double[3, 3];
                //double[,] dmat1c = new double[3, 3];
                //double[,] dmat2c = new double[3, 3];
                //double[] earthrate = new double[3];
                //double[] ldot = new double[3];
                //double[] lddot = new double[3];
                //double[] rs2cdot = new double[3];
                //double[] rs2cddot = new double[3];

                //// ---- set middle to 0, deltas to other times ----
                //// need to separate into jd and jdf portions for accuracy since it's often only seconds or minutes
                //tau12 = (jd[idx1] - jd[idx2]) * 86400.0 + (jdf[idx1] - jdf[idx2]) * 86400.0; // days to sec
                //tau13 = (jd[idx1] - jd[idx3]) * 86400.0 + (jdf[idx1] - jdf[idx3]) * 86400.0;
                //tau32 = (jd[idx3] - jd[idx2]) * 86400.0 + (jdf[idx3] - jdf[idx2]) * 86400.0;

                //// --------------------------------------------------------------
                //// using lagrange interpolation formula to derive an expression
                //// for l(t), substitute t=t2=0 and differentiate to obtain the
                //// derivatives of los.
                //// ---------------------------------------------------------------
                //s1 = -tau32 / (tau12 * tau13);
                //s2 = (tau12 + tau32) / (tau12 * tau32); // be careful here! it's -t1-t3 which can be 0
                //s3 = -tau12 / (-tau13 * tau32);
                //s4 = 2.0 / (tau12 * tau13);
                //s5 = 2.0 / (tau12 * tau32);
                //s6 = 2.0 / (-tau13 * tau32);

                //// Escobal says that for Earth orbiting satellites the ldot and lddot need additional terms
                //for (int i = 0; i < 3; i++)
                //{
                //    ldot[i] = s1 * los1[i] + s2 * los2[i] + s3 * los3[i];  // rad / s
                //    lddot[i] = s4 * los1[i] + s5 * los2[i] + s6 * los3[i];  // rad / s^2
                //}

                //lod = 0.0;
                ////double omegaearthc = gravConst.earthrot * (1.0 - lod / 86400.0) * tu;  // rad/tu
                //double omegaearth = AstroLibr.gravConst.earthrot * (1.0 - lod / 86400.0);  // rad/s
                //earthrate[0] = 0.0;
                //earthrate[1] = 0.0;
                //earthrate[2] = omegaearth;

                //MathTimeLibr.cross(earthrate, rseci2, out rs2cdot);
                //MathTimeLibr.cross(earthrate, rs2cdot, out rs2cddot);

                //for (int i = 0; i < 3; i++)
                //{
                //    dmat[i, 0] = los2[i];  // rad
                //    dmat[i, 1] = ldot[i];  // rad / TU
                //    dmat[i, 2] = lddot[i];  // rad / TU^2

                //    // ----  position determinants ----
                //    dmat1c[i, 0] = los2[i];
                //    dmat1c[i, 1] = ldot[i];
                //    dmat1c[i, 2] = rs2cddot[i];

                //    dmat2c[i, 0] = los2[i];
                //    dmat2c[i, 1] = ldot[i];
                //    dmat2c[i, 2] = rseci2[i];
                //}
                //d = 2.0 * MathTimeLibr.determinant(dmat, 3);
                //d1c = MathTimeLibr.determinant(dmat1c, 3);
                //d2c = MathTimeLibr.determinant(dmat2c, 3);

                //double axstar = d1c / d;
                //double bxstar = d2c / d;


                //magr2 = MathTimeLibr.mag(rseci2);
                //// determine elongation angle psi Escobal 271 eq 7.221
                //// note that above horizon, psi is between pi / 2 < psi < pi
                //// put position vector negative
                //// cos(psi) < 0, sin(psi) > 0
                //cos_psi = -(r2[0] * los2[0] + r2[1] * los2[1] + r2[2] * los2[2]) / magr2;
                //psi = Math.Acos(cos_psi);
                //sin_psi = Math.Sqrt(1.0 - cos_psi * cos_psi);
                //sp3 = sin_psi * sin_psi * sin_psi;

                //// find parameters M and n from the simaltaneous equations
                //M = Math.Atan2(magr2 * Math.Sin(psi), magr2 * Math.Cos(psi) - bxstar);
                //n = magr2 * sin_psi / Math.Sin(M);

                //magr23 = magr2 * magr2 * magr2;
                //coef = n * magr23 * sp3 / (AstroLibr.gravConst.mu * axstar);

                //// now either solve with a 4th observation (Escobal 271)


                //// or iterate through the possible psi angles (90 - 180 deg)
                //// search for roots at one degree per increment
                //double x, xi1, xitmp, sin_xi, cos_xi, sin_xim, cos_xim, y1, y2, fx, fx1,
                //    dfx, dfx2, del, dx0, dx1, dx2;
                //x = 0.0;
                //xi1 = 0.0;
                //fx1 = 0.0;
                //int iroot = 0;
                //int np = 90;
                //for (int psii = 0; psii <= np; psii++)
                //{
                //    xitmp = 0.50 * Math.PI * psii / np;
                //    //xi_deg = psitmp * rad;
                //    sin_xi = Math.Sin(xitmp);
                //    cos_xi = Math.Cos(xitmp);
                //    sin_xim = Math.Sin(xitmp + M);
                //    cos_xim = Math.Cos(xitmp + M);
                //    // solve simaltaneous equations
                //    // escobal 271, eq 7.228
                //    y1 = Math.Pow(sin_xi, 4);
                //    y2 = coef * sin_xim;
                //    fx = y1 + y2;

                //    if (fx == 0.0)
                //    {
                //        // writeline 'problem in root_cal  fx = 0'
                //        return;
                //    }
                //    if (psii == 0)
                //    {
                //        xi1 = xitmp;
                //        fx1 = fx;
                //    }
                //    else
                //    {
                //        if (fx1 * fx < 0.0)
                //        {
                //            x = 0.5 * (xi1 + xitmp);
                //            xi1 = xitmp;
                //            fx1 = fx;
                //            int ktr = 1;
                //            double dx = 100.0;
                //            while (Math.Abs(dx) > 1.0e-10 && ktr < 10)
                //            {
                //                sin_xi = Math.Sin(x);
                //                cos_xi = Math.Cos(x);
                //                sin_xim = Math.Sin(x + M);
                //                cos_xim = Math.Cos(x + M);
                //                y1 = Math.Pow(sin_xi, 4);
                //                y2 = coef * sin_xim;
                //                fx = y1 + y2;
                //                dfx = 4.0 * Math.Pow(sin_xi, 3) * cos_xi + coef * cos_xim;
                //                dfx2 = 12.0 * sin_xi * sin_xi * cos_xi * cos_xi - 4.0 * sin_xi * sin_xi * sin_xi * sin_xi - coef * sin_xim;
                //                del = dfx / (Math.Abs(dfx));
                //                //
                //                // Laguerre find order n where dx2 > 0
                //                int ii = 2;
                //                dx2 = 100.0;
                //                while (ii <= 9 && dx2 <= 0.0)
                //                {
                //                    dx0 = (ii - 1) * (ii - 1) * dfx * dfx;
                //                    dx1 = ii * (ii - 1) * fx * dfx2;
                //                    dx2 = dx0 - dx1;
                //                    ii = ii + 1;
                //                }
                //                if (dx2 > 0.0)
                //                    dx = ii * fx / (dfx + del * Math.Sqrt(dx2));
                //                else
                //                    dx = fx / (dfx);

                //                x = x - dx;
                //                ktr = ktr + 1;
                //            }  // while dx && ktr
                //        }  // if fx1 * fx < 0.0

                //        // converged x
                //        //xi_deg = xn * rad;
                //        double sin_xi2, cos_xi20, rho2, x3, x30, xi2, xi20, magr2n;
                //        rho2 = magr2 * sin_psi / Math.Sin(x);
                //        if (rho2 < 300000.0)  // km
                //        {
                //            //iroot = iroot + 1;
                //            //zr2(iroot) = rho2;
                //            //rho2 = magr2 * (sin_psi * Math.Cos(x) + cos_psi * Math.Sin(x)) / Math.Sin(x);
                //            //r23 = rho2 * rho2 * rho2;
                //            //rho = -ax / r23 - bx;

                //            //// determine slant range Escobal 283, eqn(7.280)
                //            //double c2losdr = 2.0 * (los2[0] * r2[0] + los2[1] * r2[1] + los2[2] * r2[2]);
                //            //double Rimag2 = r2[0] * r2[0] + r2[1] * r2[1] + r2[2] * r2[2];
                //            //double Rr2 = 4.0 * (Rimag2 - posi_range * posi_range);
                //            //csi2 = Math.Sqrt(c2losdr * c2losdr - Rr2);
                //            //rho2 = 0.50 * (-c2losdr + csi2);

                //            //for (int j = 0; j < 3; j++)
                //            //    rn[j] = magr2 * R2n[j] + rho2 * L2[j];
                //            //magr2n = Math.Sqrt(R2n(1) * R2n(1) + R2n(2) * R2n(2) + R2n(3) * R2n(3));
                //            //cos_xi20 = (R2n(1) * rn(1) + R2n(2) * rn(2) + R2n(3) * rn(3)) / xr2 / magr2;
                //            //xi20 = Math.Acos(cos_xi20);
                //            ////xi20_deg = xi20 * rad;
                //            ////
                //            //sin_xi2 = rho * Math.Sin(psi) / xr2;
                //            //xi2 = Math.Asin(sin_xi2);
                //            ////xi2_deg = xi2 * rad;  // angle opposite rho
                //            ////x30 = psi_deg + xi_deg + xi20_deg;
                //            ////x3 = psi_deg + xi_deg + xi2_deg;
                //            ////
                //            //rho2 = xr2 * Math.Sin(psi + x) / Math.Sin(psi);
                //            ////
                //            ////              write(63, 603) i,xi_deg,y1,y2,fx,xr2,rho2;
                //        }
                //        else // fx1 * fx
                //        {
                //            xi1 = xitmp;
                //            fx1 = fx;
                //        } // fx1* fx
                //    }
                //}  // for
                #endregion

                #region GTDS Approach not work yet
                // preliminary serach tries to find acceptable intiial guesses
                // unknown if this is better or not from simply using Gauss equation solution
                // it "may" be better in tougher cases...
                // gtds approach for finding initial range?
                double rho2, beta, d, delta, c1, c2, cmin,
                    abart, bbart, rho1min, rho1max, magr1min, magr1max;
                double magrs2, rng2u, rng2l, l2dotrs;
                double[] abarv = new double[3];
                double[] bbarv = new double[3];
                double[] r1min = new double[3];
                double[] r1max = new double[3];
                double[] r1x = new double[3];
                double[] r2x = new double[3];
                double[] r2xp = new double[3];
                double[] r2u = new double[3];
                double[] r2l = new double[3];
                double[,] r1u = new double[4, 2];
                double[,] r1l = new double[4, 2];

                // input initial slant range guess in km
                rho2 = bigr2x - 6378.1363;  // or 29000.0 etc ............

                // start section to find upper and lower bounds
                double K = 1.25;

                // level number 0, 1, 2, ...
                for (int L = 0; L < 3; L++)
                {
                    magrs2 = MathTimeLibr.mag(rseci2);
                    rng2u = magrs2 + rho2 * Math.Pow(K, L);
                    rng2l = magrs2 + rho2 * Math.Pow(K, -L);

                    // ---- now form two r2 position vectors ----- 
                    l2dotrs = MathTimeLibr.dot(los2, rseci2);
                    rng2u = 0.5 * (-2.0 * l2dotrs + Math.Sqrt(l2dotrs * l2dotrs - 4.0 * (magrs2 * magrs2 - rng2u * rng2u)));
                    for (int i = 0; i < 3; i++)
                        r2u[i] = rng2u * los2[i] + rseci2[i];
                    rng2l = 0.5 * (-2.0 * l2dotrs + Math.Sqrt(l2dotrs * l2dotrs - 4.0 * (magrs2 * magrs2 - rng2l * rng2l)));
                    for (int i = 0; i < 3; i++)
                        r2l[i] = rng2l * los2[i] + rseci2[i];

                    // ---- now form r1 guesses
                    magr1min = 6578.1363;  // assume 200 km alt 
                    magr1max = 100000.0;   // km

                    // in following, do for both 2u, 2l ================================================add in
                    for (int i = 0; i < 3; i++)
                    {
                        r2x[i] = r2u[i];
                        r2xp[i] = r2u[i] + rseci2[i];  // needed below for bbart
                    }

                    beta = Math.Acos(magrs2 / MathTimeLibr.mag(r2x));
                    abart = MathTimeLibr.dot(r2x, r2x) / MathTimeLibr.dot(los1, r2x);
                    bbart = MathTimeLibr.dot(r2x, r2xp) / MathTimeLibr.dot(los1, r2x);

                    for (int i = 0; i < 3; i++)
                    {
                        abarv[i] = abart * los1[i] - r2x[i];
                        bbarv[i] = rseci1[i] + bbart * los1[i] - r2x[i];
                    }

                    d = MathTimeLibr.dot(abarv, abarv) - Math.Pow((MathTimeLibr.mag(r2x) / Math.Tan(beta)), 2);
                    delta = MathTimeLibr.dot(abarv, bbarv) - d * MathTimeLibr.dot(bbarv, bbarv);

                    if (delta < 0.0)
                    {
                        //                      r1x = rs2eci + h2 * K ^ e2; // doesn't look like right formula gtds ref off?
                    }
                    else
                    {
                        // test the 2 assumed values for r2
                        c1 = -1.0 / d * (MathTimeLibr.dot(abarv, bbarv) - Math.Sqrt(delta));
                        c2 = -1.0 / d * (MathTimeLibr.dot(abarv, bbarv) + Math.Sqrt(delta));
                        cmin = Math.Cos(beta) * Math.Cos(beta) - 1.0;

                        // determine improved values of r1min and r1max
                        rho1min = 0.0;
                        rho1max = 0.0;
                        if (c1 < 0.0)
                            rho1min = abart * c1 + bbart;
                        else
                        if (c2 < 0.0)
                            rho1min = abart * c2 + bbart;

                        if (rho1min > 0.0)
                        {
                            for (int i = 0; i < 3; i++)
                                r1min[i] = rho1min * los1[i] + rseci1[i];
                            magr1min = MathTimeLibr.mag(r1min);
                        }

                        if (c1 >= 0.0)
                            rho1max = abart * c1 + bbart;
                        else
                            if (c2 >= 0.0)
                            rho1max = abart * c2 + bbart;

                        if (rho1max >= 0.0)
                        {
                            for (int i = 0; i < 3; i++)
                                r1max[i] = rho1max * los1[i] + rseci1[i];
                            magr1max = MathTimeLibr.mag(r1max);
                        }

                    }

                    // determine if any successes
                    // set the upper and lower for each value of L and r2
                //    if (magr1min < magr1mino && magr1max > magr1maxo)
                //    {
                //        r1p[L, r2] = magr1mino;
                //        r1u[L, r2] = magr1mino - magr1min;
                //        r1l[L, r2] = magr1max;
                //    }
                //    else
                //    {
                //        // but the parens = 1.0??????
                //        r1mid = magr1min + (magr2 - magr1min) * ((magr1mino - 6378.0) / (magr1mino - 6378.0));
                //        r1p[L, r2] = r1mid;
                //        r1u[L, r2] = magr1max - r1mid;
                //        r1l[L, r2] = r1mid - r1min;
                //    }

                //    if (r1mid > r1max)
                //        r1u[L, r2] = r1mid - magr1min;

                //    r1min[L, r2] = r1mid;
                //    rng1U[L, r2] = r1max - r1mid;
                //    rng1L[L, r2] = r1mid - r1min;

                //    if (r1max > rng1U)
                //        rng1U[L, r2] = r1mid - r1min;
                //    // end section to find upper and lower bounds

                //    // preliminary search for a given level
                //    r2 = rs2eci + h2 * K ^ e2;
                //    if (e1 > 0.0)
                //        r1 = r2(L, r2) + r1u(L, R2) * K ^ -e1;
                //    if (e1 < 0.0 && no upper)
                //          r1 = r1(L, r2) - r1e(L, R2) * (1.0 - J ^ -e1);
                //if (e1 < 0.0 && upper)
                //    r1 = r1(L, r2) - r1e(L, R2) * (1.0 - J ^ -e1);
                //// find Q quality valkues
                //q0 = q;
                //tau10 = r1;
                //tau20 = r2;

            }  // for L 0-3
            
                #endregion //===========================================================

                // try simply checking values
                // test case 4 escobal requires a larger initial guess - 63000? (ans)
                if (bigr2x < 0.0 || bigr2x > 50000.0)
                    bigr2x = 35000.0;  // simply set this to about GEO, allowing for less than that too. 
                                       // use Laplace guess
                                       // bigr2x = Convert.ToDouble(this.bigrL.Text);

                this.bigrD.Text = bigr2x.ToString();

                strbuildall.AppendLine("double-r initial - from Gauss " + bigr2x);

                double tau1 = (jd[idx1] - jd[idx2]) * 86400.0 + (jdf[idx1] - jdf[idx2]) * 86400.0; // days to sec
                double tau2 = (jd[idx3] - jd[idx1]) * 86400.0 + (jdf[idx3] - jdf[idx1]) * 86400.0;
                double tau3 = (jd[idx3] - jd[idx2]) * 86400.0 + (jdf[idx3] - jdf[idx2]) * 86400.0;

                // adding these in "seems" to help GEO cases where the obs are more than a day apart. 2014 JS An Improved...
                // haven't tested this for LEO orbits yet
                // use abs() ????? yes
                // for now, these are set in anglesdoubler, using abs
                n12 = (int)(tau1 / 86400.0);
                n13 = (int)(tau2 / 86400.0);
                n23 = (int)(tau3 / 86400.0);

                strbuildall.AppendLine("niis " + n12 + " " + n13 + " " + n23);
                // from jim notes on gooding for near/circular orbits
                //if (0.0 < tau3 - tau1 && tau3 - tau1 < Peroid *0.5)
                //    k = 0;
                //if (Peroid *0.5 < tau3 - tau1 && tau3-tau1 < Peroid)
                //    k = 1;
                //if (Peroid < tau3 - tau1 && tau3 - tau1 < 1.5 * Peroid)
                //    k = 2;
                //if (1.5 * Peroid < tau3 - tau1 && tau3 - tau1 < 2.0 * Peroid)
                //    k = 3;

                rng1 = bigr2x;  // old 12500 needs to be in km!! seems to do better when all the same? if too far off (*2) NAN
                rng2 = bigr2x * 1.02;  // 1.02 might be better? make the initial guess a bit different
                rng3 = bigr2x * 1.08;
                AstroLibr.anglesdoubler(tdecl[idx1], tdecl[idx2], tdecl[idx3], trtasc[idx1], trtasc[idx2], trtasc[idx3],
                     jd[idx1], jdf[idx1], jd[idx2], jdf[idx2], jd[idx3], jdf[idx3],
                     rseci1, rseci2, rseci3, rng1, rng2, out r2, out v2, out errstr, pctchg);
                string ltrstr;
                ltrstr = errstr.Substring(errstr.Length - 2, 2);
                if (ltrstr[0].Equals('\n'))
                    ltrstr = errstr.Substring(errstr.Length - 1, 1);
                else
                    if (ltrstr[1].Equals('\n'))
                        ltrstr = errstr.Substring(errstr.Length - 2, 1);

                strbuildall.AppendLine(errstr);
                posstr = "r2 " + r2[0].ToString("0.000000") + " " +
                    r2[1].ToString("0.000000") + " " + r2[2].ToString("0.000000");
                velstr = "v2 " + v2[0].ToString("0.000000") + " " +
                    v2[1].ToString("0.000000") + " " + v2[2].ToString("0.000000");
                this.rDR.Text = posstr;
                this.vDR.Text = velstr;
                strbuildall.AppendLine(posstr + velstr);
                AstroLibr.rv2coe(r2, v2, out p, out a, out ecc, out incl, out raan, out argp, out nu, out m,
                    out arglat, out truelon, out lonper);
                coestr = "a= " + a.ToString("0.0000") + " e= " + ecc.ToString("0.000000000") + " i= " + (incl * rad).ToString("0.0000")
                    + "  " + (raan * rad).ToString("0.0000") + "  " + (argp * rad).ToString("0.0000") + "  "
                    + (nu * rad).ToString("0.0000") + "  " + (m * rad).ToString("0.0000") + "  " + (arglat * rad).ToString("0.0000");
                // + " " + (truelon * rad).ToString("0.0000") + " " + (lonper * rad).ToString("0.0000"));
                this.coeDR.Text = coestr;
                strbuildall.AppendLine("Double-r " + ltrstr + " coes " + coestr);
                strbuildall.AppendLine(ans);
                strbuildallsum.AppendLine("Double-r " + ltrstr + " coes " + coestr);
                strbuildallsum.AppendLine(ans);
                AstroLibr.rv2eq(r2, v2, out a, out n, out af, out ag, out chi, out psi, out meanlonM, out meanlonNu, out fr);
                eqstr = "n= " + n.ToString("0.0000") + " af= " + af.ToString("0.000000000") + " ag= " + ag.ToString("0.0000")
                    + " chi= " + chi.ToString("0.0000") + " psi= " + psi.ToString("0.0000") + "  " + (meanlonNu * rad).ToString("0.0000")
                    + "  " + (meanlonM * rad).ToString("0.0000");
                this.eqDR.Text = eqstr;
                MathTimeLibr.cross(r2, v2, out h);
                this.hDR.Text = h[0].ToString("0.0000") + "  " + h[1].ToString("0.0000") + "  " + h[2].ToString("0.0000");
            }
            //}


            if (methodType == "Gooding" || methodType == "All")
            {
                strbuildallsum.AppendLine("Gooding -----------------------------------");
                strbuildall.AppendLine("Gooding -----------------------------------");
                numhalfrev = 0;

                AstroLibr.getGaussRoot(tdecl[idx1], tdecl[idx2], tdecl[idx3], trtasc[idx1], trtasc[idx2], trtasc[idx3],
                     jd[idx1], jdf[idx1], jd[idx2], jdf[idx2], jd[idx3], jdf[idx3],
                     rseci1, rseci2, rseci3, out bigr2x);

                // try simply checking values
                if (bigr2x < 0.0 || bigr2x > 75000.0)
                    bigr2x = 40000.0;  // simply set this to about GEO, allowing for less than that too. 

                rng1 = bigr2x;  // old 12500 needs to be in km!! seems to do better when all the same? if too far off (*2) NAN
                rng2 = rng1 * 1.02;  // 1.02 might be better? make the initial guess a bit different
                rng3 = rng1 * 1.08;
                this.bigrGoo.Text = bigr2x.ToString();

                if (n12 > 0 || n13 > 0 || n23 > 0)
                    numhalfrev = 2;


                AstroLibr.anglesgooding(tdecl[idx1], tdecl[idx2], tdecl[idx3], trtasc[idx1], trtasc[idx2], trtasc[idx3],
                    jd[idx1], jdf[idx1], jd[idx2], jdf[idx2], jd[idx3], jdf[idx3],
                    rseci1, rseci2, rseci3, numhalfrev, rng1, rng2, rng3, out r2, out v2, out errstr);
                strbuildall.AppendLine(errstr);
                posstr = "r2 " + r2[0].ToString("0.000000") + " " +
                    r2[1].ToString("0.000000") + " " + r2[2].ToString("0.000000");
                velstr = "v2 " + v2[0].ToString("0.000000") + " " +
                    v2[1].ToString("0.000000") + " " + v2[2].ToString("0.000000");
                this.rGoo.Text = posstr;
                this.vGoo.Text = velstr;
                strbuildall.AppendLine(posstr + velstr);
                AstroLibr.rv2coe(r2, v2, out p, out a, out ecc, out incl, out raan, out argp, out nu, out m,
                    out arglat, out truelon, out lonper);
                coestr = "a= " + a.ToString("0.0000") + " e= " + ecc.ToString("0.000000000") + " i= " + (incl * rad).ToString("0.0000")
                    + "  " + (raan * rad).ToString("0.0000") + "  " + (argp * rad).ToString("0.0000") + "  "
                    + (nu * rad).ToString("0.0000") + "  " + (m * rad).ToString("0.0000") + "  " + (arglat * rad).ToString("0.0000");
                // + " " + (truelon * rad).ToString("0.0000") + " " + (lonper * rad).ToString("0.0000"));
                this.coeGoo.Text = coestr;
                strbuildall.AppendLine("Gooding coes " + coestr);
                strbuildall.AppendLine(ans);
                strbuildallsum.AppendLine("               " + coestr);
                strbuildallsum.AppendLine(ans);
                AstroLibr.rv2eq(r2, v2, out a, out n, out af, out ag, out chi, out psi, out meanlonM, out meanlonNu, out fr);
                eqstr = "n= " + n.ToString("0.0000") + " af= " + af.ToString("0.000000000") + " ag= " + ag.ToString("0.0000")
                    + " chi= " + chi.ToString("0.0000") + " psi= " + psi.ToString("0.0000") + "  " + (meanlonNu * rad).ToString("0.0000")
                    + "  " + (meanlonM * rad).ToString("0.0000");
                this.eqGoo.Text = eqstr;
                MathTimeLibr.cross(r2, v2, out h);
                this.hGoo.Text = h[0].ToString("0.0000") + "  " + h[1].ToString("0.0000") + "  " + h[2].ToString("0.0000");
            }

        }  // doangletests


        private void btnSearch_Click(object sender, EventArgs e)
        {
            int dat;
            double jdut1, dut1, jdtt, jdftt, ttt, xp, yp, ddx, ddy, ddpsi, ddeps, lod;
            Int32 year, mon, day, hr, minute, sennum;
            double second;
            double[] rseci1 = new double[3];
            double[] rseci2 = new double[3];
            double[] rseci3 = new double[3];
            double[] vseci1 = new double[3];
            double[] vseci2 = new double[3];
            double[] vseci3 = new double[3];
            double[] aseci1 = new double[3];
            double[] aseci2 = new double[3];
            double[] aseci3 = new double[3];
            double[] rsecef1 = new double[3];
            double[] rsecef2 = new double[3];
            double[] rsecef3 = new double[3];
            double[] vsecef1 = new double[3];
            double[] vsecef2 = new double[3];
            double[] vsecef3 = new double[3];
            double[] asecef1 = new double[3];
            double[] asecef2 = new double[3];
            double[] asecef3 = new double[3];
            double[] jd = new double[50];
            double[] jdf = new double[50];
            double[] trtasc = new double[50];
            double[] tdecl = new double[50];
            double[] latgd = new double[50];
            double[] lon = new double[50];
            double[] alt = new double[50];
            int idx1, idx2, idx3;
            string ans;

            // read input data
            // note the input data has a # line between each case
            string infilename = this.testfilename.Text;
            string[] fileData = File.ReadAllLines(infilename);

            double rad = 180.0 / Math.PI;
            double itol = 5.0 / rad;   // rad

            this.opsStatusG.Text = "Test Angles - working";
            Refresh();

            // get EOP data in
            string nutLoc;
            Int32 ktrActObs;
            string EOPupdate;
            EOPSPWLib.iau80Class iau80arr;
            EOPSPWLib.iau06Class iau06arr;
            nutLoc = @"D:\Codes\LIBRARY\DataLib\nut80.dat";
            EOPSPWLibr.iau80in(nutLoc, out iau80arr);
            nutLoc = @"D:\Codes\LIBRARY\DataLib\";
            EOPSPWLibr.iau06in(nutLoc, out iau06arr);

            string eopFileName = @"D:\Codes\LIBRARY\DataLib\EOP-All-v1.1_2020-02-12.txt";
            EOPSPWLibr.readeop(ref EOPSPWLibr.eopdata, eopFileName, out ktrActObs, out EOPupdate);

            // now read it in
            AstroLib.xysdataClass[] xysarr = AstroLibr.xysarr;
            AstroLibr.readXYS(ref xysarr, nutLoc, "xysdata.dat");

            if (this.cbGaussList.SelectedItem == null)
                this.cbGaussList.SelectedIndex = 1;
            methodType = this.cbGaussList.SelectedItem.ToString();


            // -2 is gui, -1 is all, number runs a specific one
            int strt, stp;
            int tmpint = Convert.ToInt32(this.storedCase.Text);
            strt = 0;
            stp = 30;
            if (tmpint >= 0)
            {
                strt = tmpint;
                stp = tmpint;
            }
            if (tmpint == -1)
            {
                strt = 0;
                stp = 30;
            }

            dat = 35; // defaults
            dut1 = 0.0;
            xp = 0.0;
            yp = 0.0;
            ddpsi = 0.0;
            ddeps = 0.0;
            lod = 0.0;
            ddx = 0.0;
            ddy = 0.0;


            strbuildallsum.AppendLine("this test seems to work pretty good. Note\n"
               + "1. Gauss is used as initial guess, unless it's over 50000, then set to 35000 km\n"
               + "2. nij use abs() only case that's a little different are the jim tests - but still nearly hyprebolic.\n"
               + "3. dm is all 'S'. with 'L', only case 17 is better, and it started with ~43300 km\n"
               + "4. tol 0.1 km, chkamt 0.15 to limit changes, pct chg each of 3 iterations in double-r 0.05 \n"
               + "5. pct chg is pctchg*0.5 each iteration \n"
              // + "6. do knij/n corection  \n"
               + "7. the if w should also be just the first way (same as dm line above)" );

            if (tmpint >= -1)
            {
                for (int caseopt = strt; caseopt <= stp; caseopt++)
                {
                    strbuildall.AppendLine("caseopt " + caseopt.ToString());
                    int ktr = 1;     // skip header, go to next # comment line

                    string line = fileData[ktr];
                    line.Replace(@"\s+", " ");
                    string[] linesplt = line.Split(' ');
                    int tmpcase = Convert.ToInt32(linesplt[1]);
                    while (tmpcase != caseopt)
                    {
                        line = fileData[ktr];
                        line.Replace(@"\s+", " ");
                        linesplt = line.Split(' ');
                        if (line[0].Equals('#'))
                            tmpcase = Convert.ToInt32(linesplt[1]);

                        ktr = ktr + 1;
                    }

                    // get all the data for caseopt
                    int obsktr = 0;
                    // set the first case only
                    if (caseopt == 0)
                    {
                        ans = fileData[ktr];
                        ktr = 2;
                    }
                    else
                        ans = fileData[ktr - 1];

                    // read in the whole input file
                    while (ktr < fileData.Count() && !fileData[ktr][0].Equals('#'))
                    {
                        line = fileData[ktr];
                        linesplt = line.Split(',');
                        mon = Convert.ToInt32(linesplt[1]);
                        day = Convert.ToInt32(linesplt[0]);
                        year = Convert.ToInt32(linesplt[2]);
                        hr = Convert.ToInt32(linesplt[3]);
                        minute = Convert.ToInt32(linesplt[4]);
                        second = Convert.ToDouble(linesplt[5]);
                        MathTimeLibr.jday(year, mon, day, hr, minute, second, out jd[obsktr], out jdf[obsktr]);

                        latgd[obsktr] = Convert.ToDouble(linesplt[6]) / rad;
                        lon[obsktr] = Convert.ToDouble(linesplt[7]) / rad;
                        alt[obsktr] = Convert.ToDouble(linesplt[8]);

                        trtasc[obsktr] = Convert.ToDouble(linesplt[9]) / rad;
                        tdecl[obsktr] = Convert.ToDouble(linesplt[10]) / rad;

                        obsktr = obsktr + 1;
                        ktr = ktr + 1;
                    }

                    idx1 = 0;
                    idx2 = 1;
                    idx3 = 2;
                    strbuildallsum.AppendLine("\n\n ================================ case number " + caseopt.ToString() + " ================================");
                    strbuildall.AppendLine("\n\n ================================ case number " + caseopt.ToString() + " ================================");

                    double convrt = Math.PI / (180 * 3600.0);

                    // setup specific index where mutliple choices exist
                    switch (caseopt)
                    {
                        case 0:
                            idx1 = 2;
                            idx2 = 4;
                            idx3 = 5;

                            dut1 = -0.609641;
                            dat = 35;
                            xp = 0.137495 * convrt;
                            yp = 0.342416 * convrt;
                            lod = 0.0;
                            ddpsi = 0.0;
                            ddeps = 0.0;

                            break;
                        case 1:
                            // book example
                            dut1 = -0.609641;      // sec
                            dat = 35;              // sec
                            lod = 0.0;
                            xp = 0.137495 * convrt;  // " to rad
                            yp = 0.342416 * convrt;
                            ddpsi = 0.0;  // " to rad
                            ddeps = 0.0;
                            ddx = 0.0;    // " to rad
                            ddy = 0.0;
                            //latgd = 40.0 / rad;
                            //lon = -110.0 / rad;
                            //alt = 2.0;  // km
                            // ---- select points to use
                            idx1 = 2;
                            idx2 = 4;
                            idx3 = 5;

                            idx1 = 5;
                            idx2 = 9;
                            idx3 = 13;
                            break;
                    }  // end switch

                    // get new site info only if gui, which isn't this option
                    // if a value, set the lat lon alt
                    // if "Other", leave the values as they were input
                    //getSensorLLA(this.comboBoxSen1.SelectedItem.ToString(), ref latgd[idx1], ref lon[idx1], ref alt[idx1]);
                    //this.lat1.Text = (latgd[idx1] * rad).ToString();
                    //this.lon1.Text = (lon[idx1] * rad).ToString();
                    //this.alt1.Text = alt[idx1].ToString();

                    //jd1 = jd[idx1] + jdf[idx1];
                    //jd2 = jd[idx2] + jdf[idx2];
                    //jd3 = jd[idx3] + jdf[idx3];
                    AstroLibr.site(latgd[idx1], lon[idx1], alt[idx1], out rsecef1, out vsecef1);
                    if (caseopt != 0 && caseopt != 1)
                        EOPSPWLibr.findeopparam(jd[idx1], jdf[idx1], 's', EOPSPWLibr.eopdata, 
                            out dut1, out dat, out lod, out xp, out yp, out ddpsi, out ddeps, out ddx, out ddy);
                    //MathTimeLibr.convtime(year[idx1], mon[idx1], day[idx1], hr[idx1], minute[idx1], second[idx1], 0, dut1, dat,
                    //    out ut1, out tut1, out jdut1, out jdut1frac, out utc, out tai,
                    //    out tt, out ttt, out jdtt, out jdttfrac, out tdb, out ttdb, out jdtdb, out jdtdbfrac);
                    jdtt = jd[idx1];
                    jdftt = jdf[idx1] + (dat + 32.184) / 86400.0;
                    ttt = (jdtt + jdftt - 2451545.0) / 36525.0;
                    // note you have to use tdb for time of interst AND j2000 (when dat = 32)
                    //  ttt = (jd + jdFrac + (dat + 32.184) / 86400.0 - 2451545.0 - (32 + 32.184) / 86400.0) / 36525.0;
                    jdut1 = jd[idx1] + jdf[idx1] + dut1 / 86400.0;
                    AstroLibr.eci_ecef(ref rseci1, ref vseci1, ref aseci1, MathTimeLib.Edirection.efrom, ref rsecef1, ref vsecef1, ref asecef1,
                         iau80arr, ttt, jdut1, lod, xp, yp, ddpsi, ddeps);

                    // get new site info
                    // if a value, set the lat lon alt
                    // if "Other", leave the values as they were input
                    //getSensorLLA(this.comboBoxSen2.SelectedItem.ToString(), ref latgd[idx2], ref lon[idx2], ref alt[idx2]);
                    //this.lat2.Text = (latgd[idx2] * rad).ToString();
                    //this.lon2.Text = (lon[idx2] * rad).ToString();
                    //this.alt2.Text = alt[idx2].ToString();

                    AstroLibr.site(latgd[idx2], lon[idx2], alt[idx2], out rsecef2, out vsecef2);
                    if (caseopt != 0 && caseopt != 1)
                        EOPSPWLibr.findeopparam(jd[idx2], jdf[idx2], 's', EOPSPWLibr.eopdata, 
                        out dut1, out dat, out lod, out xp, out yp, out ddpsi, out ddeps, out ddx, out ddy);
                    //MathTimeLibr.convtime(year[idx2], mon[idx2], day[idx2], hr[idx2], minute[idx2], second[idx2], 0, dut1, dat,
                    //    out ut1, out tut1, out jdut1, out jdut1frac, out utc, out tai,
                    //    out tt, out ttt, out jdtt, out jdttfrac, out tdb, out ttdb, out jdtdb, out jdtdbfrac);
                    jdtt = jd[idx2];
                    jdftt = jdf[idx2] + (dat + 32.184) / 86400.0;
                    ttt = (jdtt + jdftt - 2451545.0) / 36525.0;
                    jdut1 = jd[idx2] + jdf[idx2] + dut1 / 86400.0;
                    AstroLibr.eci_ecef(ref rseci2, ref vseci2, ref aseci2, MathTimeLib.Edirection.efrom, ref rsecef2, ref vsecef2, ref asecef2,
                         iau80arr, ttt, jdut1, lod, xp, yp, ddpsi, ddeps);
                    double gst, lst;
                    AstroLibr.lstime(lon[idx2], jdut1, out lst, out gst);
                    strbuildall.AppendLine("\nlst " + lst.ToString() + " " + (lst * rad).ToString());


                    // get new site info
                    // if a value, set the lat lon alt
                    // if "Other", leave the values as they were input
                    //getSensorLLA(this.comboBoxSen3.SelectedItem.ToString(), ref latgd[idx3], ref lon[idx3], ref alt[idx3]);
                    //this.lat3.Text = (latgd[idx3] * rad).ToString();
                    //this.lon3.Text = (lon[idx3] * rad).ToString();
                    //this.alt3.Text = alt[idx3].ToString();

                    AstroLibr.site(latgd[idx3], lon[idx3], alt[idx3], out rsecef3, out vsecef3);
                    if (caseopt != 0 && caseopt != 1)
                        EOPSPWLibr.findeopparam(jd[idx3], jdf[idx3], 's', EOPSPWLibr.eopdata, 
                        out dut1, out dat, out lod, out xp, out yp, out ddpsi, out ddeps, out ddx, out ddy);
                    jdtt = jd[idx3];
                    jdftt = jdf[idx3] + (dat + 32.184) / 86400.0;
                    ttt = (jdtt + jdftt - 2451545.0) / 36525.0;
                    jdut1 = jd[idx3] + jdf[idx3] + dut1 / 86400.0;
                    AstroLibr.eci_ecef(ref rseci3, ref vseci3, ref aseci3, MathTimeLib.Edirection.efrom, ref rsecef3, ref vsecef3, ref asecef3,
                         iau80arr, ttt, jdut1, lod, xp, yp, ddpsi, ddeps);

                    //if (Math.Abs(latgd[idx1] - latgd[idx2]) < 0.001 && Math.Abs(latgd[idx1] - latgd[idx3]) < 0.001
                    //    && Math.Abs(lon[idx1] - lon[idx2]) < 0.001 && Math.Abs(lon[idx1] - lon[idx3]) < 0.001)
                    //    diffsites = 'n';
                    //else
                    //    diffsites = 'y';


                    // write output
                    strbuildall.AppendLine("rseci1 " + rseci1[0].ToString("0.000000") + " " +
                    rseci1[1].ToString("0.000000") + " " + rseci1[2].ToString("0.000000"));
                    strbuildall.AppendLine("rseci2 " + rseci2[0].ToString("0.000000") + " " +
                        rseci2[1].ToString("0.000000") + " " + rseci2[2].ToString("0.000000"));
                    strbuildall.AppendLine("rseci3 " + rseci3[0].ToString("0.000000") + " " +
                        rseci3[1].ToString("0.000000") + " " + rseci3[2].ToString("0.000000"));

                    if (caseopt == 19)
                    {  // curtis example -many mistakes!
                        rseci1 = new double[] { 3489.8, 3430.2, 4078.5 };
                        rseci2 = new double[] { 3460.1, 3460.1, 4078.5 };
                        rseci3 = new double[] { 3429.9, 3490.1, 4078.5 };
                    }

                    doangletests(jd, jdf, trtasc, tdecl, latgd, lon, alt, rseci1, vseci1, rseci2, vseci2, rseci3, vseci3,
                        idx1, idx2, idx3, ans);

                }  // caseopt cases

            } // if tmpint >= -1
            else
            {
                // do case from screen GUI
                // assign all the screen values
                idx1 = 1;
                idx2 = 2;
                idx3 = 3;
                latgd[idx1] = Convert.ToDouble(this.lat1.Text) / rad;
                lon[idx1] = Convert.ToDouble(this.lon1.Text) / rad;
                alt[idx1] = Convert.ToDouble(this.alt1.Text);
                latgd[idx2] = Convert.ToDouble(this.lat2.Text) / rad;
                lon[idx2] = Convert.ToDouble(this.lon2.Text) / rad;
                alt[idx2] = Convert.ToDouble(this.alt2.Text);
                latgd[idx3] = Convert.ToDouble(this.lat3.Text) / rad;
                lon[idx3] = Convert.ToDouble(this.lon3.Text) / rad;
                alt[idx3] = Convert.ToDouble(this.alt3.Text);

                MathTimeLibr.STKtime2JD(this.time1.Text, out jd[idx1], out jdf[idx1]);
                MathTimeLibr.STKtime2JD(this.time2.Text, out jd[idx2], out jdf[idx2]);
                MathTimeLibr.STKtime2JD(this.time3.Text, out jd[idx3], out jdf[idx3]);

                trtasc[idx1] = Convert.ToDouble(this.trtasc1.Text) / rad;
                tdecl[idx1] = Convert.ToDouble(this.tdecl1.Text) / rad;
                trtasc[idx2] = Convert.ToDouble(this.trtasc2.Text) / rad;
                tdecl[idx2] = Convert.ToDouble(this.tdecl2.Text) / rad;
                trtasc[idx3] = Convert.ToDouble(this.trtasc3.Text) / rad;
                tdecl[idx3] = Convert.ToDouble(this.tdecl3.Text) / rad;

                // save gui values for next time
                TestAll.Properties.Settings.Default.trtasc1 = this.trtasc1.Text;
                TestAll.Properties.Settings.Default.tdecl1 = this.tdecl1.Text;
                TestAll.Properties.Settings.Default.trtasc2 = this.trtasc2.Text;
                TestAll.Properties.Settings.Default.tdecl2 = this.tdecl2.Text;
                TestAll.Properties.Settings.Default.trtasc3 = this.trtasc3.Text;
                TestAll.Properties.Settings.Default.tdecl3 = this.tdecl3.Text;
                TestAll.Properties.Settings.Default.time1 = this.time1.Text;
                TestAll.Properties.Settings.Default.time2 = this.time2.Text;
                TestAll.Properties.Settings.Default.time3 = this.time3.Text;
                TestAll.Properties.Settings.Default.lat1 = this.lat1.Text;
                TestAll.Properties.Settings.Default.lon1 = this.lon1.Text;
                TestAll.Properties.Settings.Default.alt1 = this.alt1.Text;
                TestAll.Properties.Settings.Default.lat2 = this.lat2.Text;
                TestAll.Properties.Settings.Default.lon2 = this.lon2.Text;
                TestAll.Properties.Settings.Default.alt2 = this.alt2.Text;
                TestAll.Properties.Settings.Default.lat3 = this.lat3.Text;
                TestAll.Properties.Settings.Default.lon3 = this.lon3.Text;
                TestAll.Properties.Settings.Default.alt3 = this.alt3.Text;
                TestAll.Properties.Settings.Default.Save();

                // get new site info from gui
                // if a value, set the lat lon alt
                // if "Other", leave the values as they were input
                // -1 - don't assign the obsclassarr values
                getSensorLLA(-1, this.comboBoxSen1.SelectedItem.ToString(), out sennum, ref latgd[idx1], ref lon[idx1], ref alt[idx1]);
                this.lat1.Text = (latgd[idx1] * rad).ToString();
                this.lon1.Text = (lon[idx1] * rad).ToString();
                this.alt1.Text = alt[idx1].ToString();

                //jd1 = jd[idx1] + jdf[idx1];
                //jd2 = jd[idx2] + jdf[idx2];
                //jd3 = jd[idx3] + jdf[idx3];
                AstroLibr.site(latgd[idx1], lon[idx1], alt[idx1], out rsecef1, out vsecef1);
                EOPSPWLibr.findeopparam(jd[idx1], jdf[idx1], 's', EOPSPWLibr.eopdata, 
                    out dut1, out dat, out lod, out xp, out yp, out ddpsi, out ddeps, out ddx, out ddy);
                //MathTimeLibr.convtime(year[idx1], mon[idx1], day[idx1], hr[idx1], minute[idx1], second[idx1], 0, dut1, dat,
                //    out ut1, out tut1, out jdut1, out jdut1frac, out utc, out tai,
                //    out tt, out ttt, out jdtt, out jdttfrac, out tdb, out ttdb, out jdtdb, out jdtdbfrac);
                jdtt = jd[idx1];
                jdftt = jdf[idx1] + (dat + 32.184) / 86400.0;
                ttt = (jdtt + jdftt - 2451545.0) / 36525.0;
                // note you have to use tdb for time of interst AND j2000 (when dat = 32)
                //  ttt = (jd + jdFrac + (dat + 32.184) / 86400.0 - 2451545.0 - (32 + 32.184) / 86400.0) / 36525.0;
                jdut1 = jd[idx1] + jdf[idx1] + dut1 / 86400.0;
                AstroLibr.eci_ecef(ref rseci1, ref vseci1, ref aseci1, MathTimeLib.Edirection.efrom, ref rsecef1, ref vsecef1, ref asecef1,
                         iau80arr, ttt, jdut1, lod, xp, yp, ddpsi, ddeps);

                // get new site info
                // if a value, set the lat lon alt
                // if "Other", leave the values as they were input
                // -1 - don't assign the obsclassarr values
                getSensorLLA(-1, this.comboBoxSen2.SelectedItem.ToString(), out sennum, ref latgd[idx2], ref lon[idx2], ref alt[idx2]);
                this.lat2.Text = (latgd[idx2] * rad).ToString();
                this.lon2.Text = (lon[idx2] * rad).ToString();
                this.alt2.Text = alt[idx2].ToString();

                AstroLibr.site(latgd[idx2], lon[idx2], alt[idx2], out rsecef2, out vsecef2);
                EOPSPWLibr.findeopparam(jd[idx2], jdf[idx2], 's', EOPSPWLibr.eopdata, 
                    out dut1, out dat, out lod, out xp, out yp, out ddpsi, out ddeps, out ddx, out ddy);
                //MathTimeLibr.convtime(year[idx2], mon[idx2], day[idx2], hr[idx2], minute[idx2], second[idx2], 0, dut1, dat,
                //    out ut1, out tut1, out jdut1, out jdut1frac, out utc, out tai,
                //    out tt, out ttt, out jdtt, out jdttfrac, out tdb, out ttdb, out jdtdb, out jdtdbfrac);
                jdtt = jd[idx2];
                jdftt = jdf[idx2] + (dat + 32.184) / 86400.0;
                ttt = (jdtt + jdftt - 2451545.0) / 36525.0;
                jdut1 = jd[idx2] + jdf[idx2] + dut1 / 86400.0;
                AstroLibr.eci_ecef(ref rseci2, ref vseci2, ref aseci2, MathTimeLib.Edirection.efrom, ref rsecef2, ref vsecef2, ref asecef2,
                         iau80arr, ttt, jdut1, lod, xp, yp, ddpsi, ddeps);
                double gst, lst;
                AstroLibr.lstime(lon[idx2], jdut1, out lst, out gst);
                strbuildall.AppendLine("\nlst " + lst.ToString() + " " + (lst * rad).ToString());


                // get new site info
                // if a value, set the lat lon alt
                // if "Other", leave the values as they were input
                // -1 - don't assign the obsclassarr values
                getSensorLLA(-1, this.comboBoxSen3.SelectedItem.ToString(), out sennum, ref latgd[idx3], ref lon[idx3], ref alt[idx3]);
                this.lat3.Text = (latgd[idx3] * rad).ToString();
                this.lon3.Text = (lon[idx3] * rad).ToString();
                this.alt3.Text = alt[idx3].ToString();

                AstroLibr.site(latgd[idx3], lon[idx3], alt[idx3], out rsecef3, out vsecef3);
                EOPSPWLibr.findeopparam(jd[idx3], jdf[idx3], 's', EOPSPWLibr.eopdata,
                    out dut1, out dat, out lod, out xp, out yp, out ddpsi, out ddeps, out ddx, out ddy);
                jdtt = jd[idx3];
                jdftt = jdf[idx3] + (dat + 32.184) / 86400.0;
                ttt = (jdtt + jdftt - 2451545.0) / 36525.0;
                jdut1 = jd[idx3] + jdf[idx3] + dut1 / 86400.0;
                AstroLibr.eci_ecef(ref rseci3, ref vseci3, ref aseci3, MathTimeLib.Edirection.efrom, ref rsecef3, ref vsecef3, ref asecef3,
                         iau80arr, ttt, jdut1, lod, xp, yp, ddpsi, ddeps);

                //if (Math.Abs(latgd[idx1] - latgd[idx2]) < 0.001 && Math.Abs(latgd[idx1] - latgd[idx3]) < 0.001
                //    && Math.Abs(lon[idx1] - lon[idx2]) < 0.001 && Math.Abs(lon[idx1] - lon[idx3]) < 0.001)
                //    diffsites = 'n';
                //else
                //    diffsites = 'y';


                // write output
                strbuildall.AppendLine("rseci1 " + rseci1[0].ToString("0.000000") + " " +
                rseci1[1].ToString("0.000000") + " " + rseci1[2].ToString("0.000000"));
                strbuildall.AppendLine("rseci2 " + rseci2[0].ToString("0.000000") + " " +
                    rseci2[1].ToString("0.000000") + " " + rseci2[2].ToString("0.000000"));
                strbuildall.AppendLine("rseci3 " + rseci3[0].ToString("0.000000") + " " +
                    rseci3[1].ToString("0.000000") + " " + rseci3[2].ToString("0.000000"));
                ans = "";

                doangletests(jd, jdf, trtasc, tdecl, latgd, lon, alt, rseci1, vseci1, rseci2, vseci2, rseci3, vseci3,
                    idx1, idx2, idx3, ans);
            }  // GUI test

            string directory = @"D:\Codes\LIBRARY\cs\TestAll\";
            strbuild.AppendLine("angles only tests case results written to " + directory + "testall-Angles.out ");
            strbuild.AppendLine(@"geo data for chap 9 plot written to D:\faabook\current\excel\testgeo.out for ch9 plot ");

            File.WriteAllText(directory + "testall-Angles.out", strbuildall.ToString());
            File.WriteAllText(directory + "testall-Anglessum.out", strbuildallsum.ToString());
            File.WriteAllText(directory + "testall-Anglesobs.out", strbuildObs.ToString());

            this.opsStatusG.Text = "Done";
            Refresh();

        }  // test gauss

        private void comboBoxSen1_SelectedIndexChanged(object sender, EventArgs e)
        {
            double latgd, lon, alt, rad;
            int sennum;
            rad = 180.0 / Math.PI;

            latgd = Convert.ToDouble(this.lat1.Text) / rad;
            lon = Convert.ToDouble(this.lon1.Text) / rad;
            alt = Convert.ToDouble(this.alt1.Text);

            // get new site info
            if (this.comboBoxSen1.SelectedItem.ToString().Equals("OtherLLA"))
            {
                // leave the lat lon alt values as they were input
                // set eci vectors
                //xxx
            }
            else if (this.comboBoxSen1.SelectedItem.ToString().Equals("OtherECI"))
            {
                // get the lat lon alt from eci
                getSensorECI(this.rseci1t.Text, this.time1.Text, out latgd, out lon, out alt);
            }
            else
            {
                // if a sensor value given, get jtm location
                // also has the Other option in there which keeps it the same
                // -1 - don't assign the obsclassarr values
                getSensorLLA(-1, this.comboBoxSen1.SelectedItem.ToString(), out sennum, ref latgd, ref lon, ref alt);
            }

            this.lat1.Text = (latgd * rad).ToString();
            this.lon1.Text = (lon * rad).ToString();
            this.alt1.Text = alt.ToString();
        }

        private void comboBoxSen2_SelectedIndexChanged(object sender, EventArgs e)
        {
            double latgd, lon, alt, rad;
            int sennum;
            rad = 180.0 / Math.PI;

            latgd = Convert.ToDouble(this.lat2.Text) / rad;
            lon = Convert.ToDouble(this.lon2.Text) / rad;
            alt = Convert.ToDouble(this.alt2.Text);

            // get new site info
            if (this.comboBoxSen2.SelectedItem.ToString().Equals("OtherLLA"))
            {
                // leave the lat lon alt values as they were input
                // set eci vectors
                //xxx
            }
            else if (this.comboBoxSen2.SelectedItem.ToString().Equals("OtherECI"))
            {
                // get the lat lon alt from eci
                getSensorECI(this.rseci2t.Text, this.time2.Text, out latgd, out lon, out alt);
            }
            else
            {
                // if a sensor value given, get jtm location
                // also has the Other option in there which keeps it the same
                // -1 - don't assign the obsclassarr values
                getSensorLLA(-1, this.comboBoxSen2.SelectedItem.ToString(), out sennum, ref latgd, ref lon, ref alt);
            }

            this.lat2.Text = (latgd * rad).ToString();
            this.lon2.Text = (lon * rad).ToString();
            this.alt2.Text = alt.ToString();
        }

        private void comboBoxSen3_SelectedIndexChanged(object sender, EventArgs e)
        {
            double latgd, lon, alt, rad;
            int sennum;
            rad = 180.0 / Math.PI;

            latgd = Convert.ToDouble(this.lat3.Text) / rad;
            lon = Convert.ToDouble(this.lon3.Text) / rad;
            alt = Convert.ToDouble(this.alt3.Text);

            // get new site info
            if (this.comboBoxSen3.SelectedItem.ToString().Equals("OtherLLA"))
            {
                // leave the lat lon alt values as they were input
                // set eci vectors
                //xxx
            }
            else if (this.comboBoxSen3.SelectedItem.ToString().Equals("OtherECI"))
            {
                // get the lat lon alt from eci
                getSensorECI(this.rseci3t.Text, this.time3.Text, out latgd, out lon, out alt);
            }
            else
            {
                // if a sensor value given, get jtm location
                // also has the Other option in there which keeps it the same
                // -1 - don't assign the obsclassarr values
                getSensorLLA(-1, this.comboBoxSen3.SelectedItem.ToString(), out sennum, ref latgd, ref lon, ref alt);
            }

            this.lat3.Text = (latgd * rad).ToString();
            this.lon3.Text = (lon * rad).ToString();
            this.alt3.Text = alt.ToString();
        }


        // writeout input to file so it can be added later 
        private void button1_Click(object sender, EventArgs e)
        {
            double jd, jdf, second;
            int year, month, day, hr, minute;

            strbuildallsum.AppendLine("# 30 new test case ");
            MathTimeLibr.STKtime2JD(this.time1.Text, out jd, out jdf);
            MathTimeLibr.invjday(jd, jdf, out year, out month, out day, out hr, out minute, out second);
            strbuildallsum.AppendLine(day.ToString("00") + "," + month + "," + year + ","
                 + hr.ToString("00") + "," + minute.ToString("00") + "," + second.ToString("00.000000000") + ","
                 + this.lat1.Text + "," + this.lon1.Text + "," + this.alt1.Text + ","
                 + this.trtasc1.Text + "," + this.tdecl1.Text + ",");

            MathTimeLibr.STKtime2JD(this.time2.Text, out jd, out jdf);
            MathTimeLibr.invjday(jd, jdf, out year, out month, out day, out hr, out minute, out second);
            strbuildallsum.AppendLine(day.ToString("00") + "," + month + "," + year + ","
                 + hr.ToString("00") + "," + minute.ToString("00") + "," + second.ToString("00.000000000") + ","
                 + this.lat2.Text + "," + this.lon2.Text + "," + this.alt2.Text + ","
                 + this.trtasc2.Text + "," + this.tdecl2.Text + ",");

            MathTimeLibr.STKtime2JD(this.time3.Text, out jd, out jdf);
            MathTimeLibr.invjday(jd, jdf, out year, out month, out day, out hr, out minute, out second);
            strbuildallsum.AppendLine(day.ToString("00") + "," + month + "," + year + ","
                 + hr.ToString("00") + "," + minute.ToString("00") + "," + second.ToString("00.000000000") + ","
                 + this.lat3.Text + "," + this.lon3.Text + "," + this.alt3.Text + ","
                 + this.trtasc3.Text + "," + this.tdecl3.Text + ",");
        }



        // take jims results and build input test case
        private void button2_Click(object sender, EventArgs e)
        {
            string t1, t2, t3, t1h, t2h, t3h, rtasc1, rtasc2, rtasc3, decl1, decl2, decl3, rs1, rs2, rs3;
            string s241, s242, s243;

            s241 = "-7.411645,72.451909,-0.0630844,";
            s242 = "-7.411646,72.452308,-0.0633503,";
            s243 = "-7.411874,72.452251,-0.0633579,";

            string infilename = @"D:\datafile\CSSI\COMSPOC\UCT\DoubleR_Results.txt";
            string[] fileData = File.ReadAllLines(infilename);

            // Obs: 2009-02-09T15:39:00.841000, 66.4741666666667, 1.2877, 2009-02-09T19:16:19.334000, 120.923333333333, 1.3061, 2009-02-09T22:43:52.944000, 172.9275, 1.3149
            // Tracker numbers: 241, 242, 241

            this.opsStatusG.Text = "Status:  starting to do jf results";
            Refresh();

            int ktr = 0;
            int i = 0;
            ktr = 0;
            while (ktr < fileData.Count())
            {
                // skip comment lines
                if (fileData[ktr][0] != '/')
                {
                    if (fileData[ktr].Contains("Obs:"))
                    {
                        i = i + 1;
                        string line = fileData[ktr];
                        line.Replace(@"\s+", " ");
                        string[] linesplt = line.Split(',');

                        t1 = linesplt[0].Substring(5, linesplt[0].Length - 5);
                        t1h = t1.Split('T')[1];
                        t1 = t1.Split('T')[0];
                        rtasc1 = linesplt[1];
                        decl1 = linesplt[2];

                        t2 = linesplt[3].Trim();
                        t2h = t2.Split('T')[1];
                        t2 = t2.Split('T')[0];
                        rtasc2 = linesplt[4];
                        decl2 = linesplt[5];

                        t3 = linesplt[6].Trim();
                        t3h = t3.Split('T')[1];
                        t3 = t3.Split('T')[0];
                        rtasc3 = linesplt[7];
                        decl3 = linesplt[8];

                        ktr = ktr + 1;
                        line = fileData[ktr];
                        linesplt = line.Split(',');
                        rs1 = linesplt[0].Substring(linesplt[0].Length - 3, 3);
                        rs2 = linesplt[1].Trim();
                        rs3 = linesplt[2].Trim();

                        // get est from Jim's runs of answer
                        ktr = ktr + 1;
                        line = fileData[ktr];
                        linesplt = line.Split(':');

                        strbuildallsum.AppendLine("# " + (33 + i).ToString() + " test case " + linesplt[1]);
                        string line1 = t1.Substring(8, 2) + "," + t1.Substring(5, 2) + "," + t1.Substring(0, 4) + "," + t1h.Replace(":", ",");
                        if (rs1.Equals("241"))
                            line1 = line1 + "," + s241;
                        if (rs1.Equals("242"))
                            line1 = line1 + "," + s242;
                        if (rs1.Equals("243"))
                            line1 = line1 + "," + s243;
                        strbuildallsum.AppendLine(line1 + rtasc1 + "," + decl1);

                        line1 = t2.Substring(8, 2) + "," + t2.Substring(5, 2) + "," + t2.Substring(0, 4) + "," + t2h.Replace(":", ",");
                        if (rs2.Equals("241"))
                            line1 = line1 + "," + s241;
                        if (rs2.Equals("242"))
                            line1 = line1 + "," + s242;
                        if (rs2.Equals("243"))
                            line1 = line1 + "," + s243;
                        strbuildallsum.AppendLine(line1 + rtasc2 + "," + decl2);

                        line1 = t3.Substring(8, 2) + "," + t3.Substring(5, 2) + "," + t3.Substring(0, 4) + "," + t3h.Replace(":", ",");
                        if (rs3.Equals("241"))
                            line1 = line1 + "," + s241;
                        if (rs3.Equals("242"))
                            line1 = line1 + "," + s242;
                        if (rs3.Equals("243"))
                            line1 = line1 + "," + s243;
                        strbuildallsum.AppendLine(line1 + rtasc3 + "," + decl3);
                    }
                }

                ktr = ktr + 1;
            }  // while

            string directory = @"D:\Codes\LIBRARY\cs\TestAll\";
            File.WriteAllText(directory + "testall-Anglessum.out", strbuildallsum.ToString());

            this.opsStatusG.Text = "Status:  done";
            Refresh();
        }


        //
        // read JTM obs file (jtm csv) --------------------------------------------------------
        // produces obsClassArr
        // each file comes through trackcleaner so it is already somewhat assembled into tracks
        // run linear fit to get the fit coefficients on each track as it comes in
        //
        public void readJTMobs
            (
              string infilename, ref obsClass[] obsClassArr, ref Int32 numobs
            )
        {
           // obsClassArr = new obsClass[50000];
            double obsvalue, lon, alt, second;
            int year; int month; int day; int hr; int minute;
            string[] monthStr = new string[13] { "", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };
            Int32 i; //, obsktr;
            //Int32 obsintrk;
            string[] linedata;
            double[] rseci = new double[3];
            double[] vseci = new double[3];
            double[] aseci = new double[3];
            double[] rsecef = new double[3];
            double[] vsecef = new double[3];
            double[] asecef = new double[3];

            double rtasc0, decl0, t0;

            double rad = 180.0 / Math.PI;
            rtasc0 = 0.0;
            decl0 = 0.0;
            t0 = 0.0;

            // -------------------------------- read EOP data in ----------------------------------
            // xxxxxxxxdo this once outside....
            int dat;
            double dut1, lod, xp, yp, ddpsi, ddeps, ddx, ddy;
            double ttt, jdut1;

            string nutLoc;
            Int32 ktrActObs;
            string EOPupdate;
            EOPSPWLib.iau80Class iau80arr;
            EOPSPWLib.iau06Class iau06arr;
            nutLoc = @"D:\Codes\LIBRARY\DataLib\nut80.dat";
            EOPSPWLibr.iau80in(nutLoc, out iau80arr);
            nutLoc = @"D:\Codes\LIBRARY\DataLib\";
            EOPSPWLibr.iau06in(nutLoc, out iau06arr);

            string eopFileName = @"D:\Codes\LIBRARY\DataLib\EOP-All-v1.1_2022-01-15.txt";
            EOPSPWLibr.readeop(ref EOPSPWLibr.eopdata, eopFileName, out ktrActObs, out EOPupdate);

            AstroLib.xysdataClass[] xysarr = AstroLibr.xysarr;
            AstroLibr.readXYS(ref xysarr, nutLoc, "xysdata.dat");

            string[] fileData = File.ReadAllLines(infilename);

            //numobs = 0;
            // find the original assignment
            string tmpstr = Regex.Replace(fileData[3], @"\s+", " ");
            linedata = tmpstr.Split(' ');
            //string sscnum = linedata[0].Remove(0, 2);  // Convert.ToInt32(
            string sscnum = linedata[3];  // Convert.ToInt32(
            int sennum;

            // check for being past header info (#)
            for (i = 1; i < fileData.Count(); i++)
            {
                if (fileData[i].Length > 1 && fileData[i][0] != '#')
                {
                    // set new record as they are needed. obs will start at 1
                    numobs = numobs + 1;
                    obsClassArr[numobs] = new obsClass();

                    linedata = fileData[i].Split(',');

                    string[] linedata1 = new string[3];
                    string[] linedata2 = new string[3];
                    linedata1 = linedata[1].Split('T');
                    linedata2 = linedata1[0].Split('-');
                    day = Convert.ToInt32(linedata2[2]);
                    month = Convert.ToInt32(linedata2[1]);
                    year = Convert.ToInt32(linedata2[0]);

                    linedata2 = linedata1[1].Split(':');
                    hr = Convert.ToInt32(linedata2[0]);
                    minute = Convert.ToInt32(linedata2[1]);
                    second = Convert.ToDouble(linedata2[2]);
                    obsClassArr[numobs].datestr = day.ToString() + ' ' + month.ToString() + ' ' + year.ToString()
                        + ' ' + linedata2[0] + ' ' + linedata2[1] + ' ' + linedata2[2];
                    obsClassArr[numobs].year = year;
                    obsClassArr[numobs].month = month;
                    obsClassArr[numobs].day = day;
                    obsClassArr[numobs].hr = hr;
                    obsClassArr[numobs].minute = minute;
                    obsClassArr[numobs].second = second;

                    MathTimeLibr.jday(year, month, day, hr, minute, second,
                        out obsClassArr[numobs].jd, out obsClassArr[numobs].jdf);
                    obsClassArr[numobs].mjd = obsClassArr[numobs].jd + obsClassArr[numobs].jdf - 2400000.5;
                    string senname = linedata[2];
                    // obsClassArr[numobs].sennum = Convert.ToInt32(linedata[7]);

                    double latgd;
                    latgd = 0.0;
                    lon = 0.0;
                    alt = 0.0;
                    getSensorLLA(numobs, senname, out sennum, ref latgd, ref lon, ref alt);
                    obsClassArr[numobs].sennum = sennum;

                    //sscnum
                    //obsClassArr[numobs].sscnum = 99999;
                    //obsClassArr[numobs].sscnum = Convert.ToInt32(linedata[3]);
                    obsClassArr[numobs].sscnum = sscnum;

                    obsClassArr[numobs].lat = senClassArr[0].senlat;
                    obsClassArr[numobs].lon = senClassArr[0].senlon;
                    obsClassArr[numobs].alt = senClassArr[0].senalt;

                    AstroLibr.site(obsClassArr[numobs].lat, obsClassArr[numobs].lon, obsClassArr[numobs].alt, out rsecef, out vsecef);
                    EOPSPWLibr.findeopparam(obsClassArr[numobs].jd, obsClassArr[numobs].jdf, 's', EOPSPWLibr.eopdata,
                        out dut1, out dat, out lod, out xp, out yp, out ddpsi, out ddeps, out ddx, out ddy);
                    double jdtt = obsClassArr[numobs].jd;
                    double jdftt = obsClassArr[numobs].jdf + (dat + 32.184) / 86400.0;
                    ttt = (jdtt + jdftt - 2451545.0) / 36525.0;
                    jdut1 = obsClassArr[numobs].jd + obsClassArr[numobs].jdf + dut1 / 86400.0;
                    AstroLibr.eci_ecef(ref rseci, ref vseci, ref aseci, MathTimeLib.Edirection.efrom, ref rsecef, ref vsecef, ref asecef, 
                         iau80arr, ttt, jdut1, lod, xp, yp, ddpsi, ddeps);
                    obsClassArr[numobs].rsecef[0] = rsecef[0];
                    obsClassArr[numobs].rsecef[1] = rsecef[1];
                    obsClassArr[numobs].rsecef[2] = rsecef[2];
                    obsClassArr[numobs].rseci[0] = rseci[0];
                    obsClassArr[numobs].rseci[1] = rseci[1];
                    obsClassArr[numobs].rseci[2] = rseci[2];
                    obsClassArr[numobs].dut1 = dut1;
                    obsClassArr[numobs].lod = lod;
                    obsClassArr[numobs].dat = dat;
                    obsClassArr[numobs].xp = xp;
                    obsClassArr[numobs].yp = yp;
                    obsClassArr[numobs].ddpsi = ddpsi;
                    obsClassArr[numobs].ddeps = ddeps;
                    obsClassArr[numobs].ddx = ddx;
                    obsClassArr[numobs].ddy = ddy;

                    // for debuging only - leave in degrees
                    //obsvalue = (Convert.ToDouble(linedata[6])); 
                    obsvalue = (Convert.ToDouble(linedata[6])) / rad;
                    obsClassArr[numobs].rtasc = obsvalue;
                    //obsvalue = (Convert.ToDouble(linedata[7])); 
                    obsvalue = (Convert.ToDouble(linedata[7])) / rad;
                    obsClassArr[numobs].decl = obsvalue;

                    // set initial linear values
                    rtasc0 = obsClassArr[numobs].rtasc;
                    decl0 = obsClassArr[numobs].decl;
                    t0 = obsClassArr[numobs].mjd;
                }


                Array.Clear(trackClassArr, 0, trackClassArr.Length);
                int totalTracks = 1;
                trackClassArr[totalTracks] = new trackClass();

                int oldpassnum = 0;
                int obsintrk = 0;
                for (int jj = 1; jj < numobs; jj++)
                {
                    if (oldpassnum > 0 && obsClassArr[jj].passnum != oldpassnum)
                    {
                        trackClassArr[totalTracks].passlen =
                            (trackClassArr[totalTracks].tracks[trackClassArr[totalTracks].obsintrk].jd
                            + trackClassArr[totalTracks].tracks[trackClassArr[totalTracks].obsintrk].jdf
                            - trackClassArr[totalTracks].tracks[1].jd - trackClassArr[totalTracks].tracks[1].jdf) * 86400.0;

                        totalTracks = totalTracks + 1;
                        trackClassArr[totalTracks] = new trackClass();
                        obsintrk = 0;
                    }

                    // assign new track obs info
                    obsintrk = obsintrk + 1;
                    trackClassArr[totalTracks].tracks[obsintrk] = obsClassArr[jj];
                    trackClassArr[totalTracks].obsintrk = obsintrk;
                    trackClassArr[totalTracks].sscnum = Convert.ToInt32(obsClassArr[jj].sscnum);
                    trackClassArr[totalTracks].trackktr = obsClassArr[jj].passnum;
                    oldpassnum = obsClassArr[jj].passnum;
                }

            }  // for i reading all obs in

        }  // readJTMObs




        // process jTM tracks
        // read them in and run various iods on the tracks
        private void button3_Click(object sender, EventArgs e)
        {
            StringBuilder strbuild = new StringBuilder();
            string directory, strname;

            this.opsStatusG.Text = "Status: Working on jtm files";
            Refresh();

            strname = @"D:\SDAODEvaluation\DOC_Test\JTM\temptest\99940_99940_JTMA-004_2021-10-22T033631.522000000.jtm";

            OpenFileDialog openFileDialog1 = new OpenFileDialog();
            if (Directory.Exists(Path.GetDirectoryName(strname)))
            {
                openFileDialog1.InitialDirectory = Path.GetDirectoryName(strname);
            }
            else
            {
                openFileDialog1.InitialDirectory = "c:\\";
            }
            openFileDialog1.Filter = "JTM files (*.jtm)|*.jtm|(*.*)|*.*";
            openFileDialog1.FilterIndex = 1;
            // Allow the user to select multiple files 
            openFileDialog1.Multiselect = true;
            openFileDialog1.RestoreDirectory = true;
            openFileDialog1.ShowDialog();

            // ------- setup list of filenames
            fileNames.Clear();  // make sure it is empty
            foreach (var files in openFileDialog1.FileNames)
                fileNames.Add(files);

            // read in jtmdata
            // setup the numbers for each case
            int num1 = 16;
            int num2 = 16;
            int num3 = 7;
            int num4 = 48;

            int numobs = 0;
            for (int i = 0; i < 4; i++)
            {
                readJTMobs(fileNames[i], ref obsClassArr, ref numobs);
                if (i == 0)
                    num1 = numobs;
                else if (i == 1)
                    num2 = numobs - num1;
                else if (i == 2)
                    num3 = numobs - num2 - num1;
                else if (i == 3)
                    num4 = numobs - num3 - num2 - num1;
            }

            // process them
            // setup so each obs in a track goes to the middle point of each other track

            int ctr1 = Convert.ToInt32(Math.Floor(num1 * 0.5));
            int ctr2 = Convert.ToInt32(Math.Floor(num1 + num2 * 0.5));
            int ctr3 = Convert.ToInt32(Math.Floor(num1 + num2 + num3 * 0.5));
            int ctr4 = Convert.ToInt32(Math.Floor(num1 + num2 + num3 + num4 * 0.5));

            string errstr, coestr;
            double trtasc1, tdecl1, trtasc2, tdecl2, trtasc3, tdecl3, jd1, jdf1, jd2, jdf2, jd3, jdf3, rad;
            double trtasct, tdeclt, jdt, jdft;
            double p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper;
            int ktr, ktrstart, ktr1, ktr2, ktr3;
            double[] r2 = new double[3]; 
            double[] v2 = new double[3]; 
            double[] rsecit = new double[3]; 
            double[] rseci1 = new double[3]; 
            double[] rseci2 = new double[3];
            double[] rseci3 = new double[3];
            rad = 180.0 / (2.0 * Math.PI);
            ktrstart = 1;
            ktr1 = 1;
            ktr2 = 2;
            ktr3 = 3;
            // go through the 4 cases
            for (int zz = 0; zz < 4; zz++)
            {
                switch (zz)
                {
                    case 0:
                        ktrstart = num1;
                        ktr1 = ctr1;
                        ktr2 = ctr2;
                        ktr3 = ctr3;
                        break;
                    case 1:
                        ktrstart = num1;
                        ktr1 = ctr1;
                        ktr2 = ctr3;
                        ktr3 = ctr4;
                        break;
                    case 2:
                        ktrstart = num1 + num2;
                        ktr1 = ctr2;
                        ktr2 = ctr3;
                        ktr3 = ctr4;
                        break;
                    case 3:
                        ktrstart = num1;
                        ktr1 = ctr1;
                        ktr2 = ctr3;
                        ktr3 = ctr4;
                        break;
                }  // switch
                strbuild.AppendLine("Case " + zz);
                for (ktr = ktrstart; ktr < ktrstart + num1; ktr++)
                {
                    // this one will check each obs in the track
                    tdecl1 = obsClassArr[ktr].decl;
                    trtasc1 = obsClassArr[ktr].rtasc;
                    jd1 = obsClassArr[ktr].jd;
                    jdf1 = obsClassArr[ktr].jdf;
                    rseci1 = obsClassArr[ktr].rseci;

                    // these will just be the center ones
                    tdecl2 = obsClassArr[ktr2].decl;
                    trtasc2 = obsClassArr[ktr2].rtasc;
                    jd2 = obsClassArr[ktr2].jd;
                    jdf2 = obsClassArr[ktr2].jdf;
                    rseci2 = obsClassArr[ktr2].rseci;

                    tdecl3 = obsClassArr[ktr3].decl;
                    trtasc3 = obsClassArr[ktr3].rtasc;
                    jd3 = obsClassArr[ktr3].jd;
                    jdf3 = obsClassArr[ktr3].jdf;
                    rseci3 = obsClassArr[ktr3].rseci;

                    // but these need to be time ordered
                    // get earliest one into 1 position
                    if (jd2 + jdf2 < jd1 + jdf1)
                    {
                        tdeclt = tdecl1;
                        trtasct = trtasc1;
                        jdt = jd1;
                        jdft = jdf1;
                        rsecit = rseci1;

                        tdecl1 = tdecl2;
                        trtasc1 = trtasc2;
                        jd1 = jd2;
                        jdf1 = jdf2;
                        rseci1 = rseci2;

                        tdecl2 = tdeclt;
                        trtasc2 = trtasct;
                        jd2 = jdt;
                        jdf2 = jdft;
                        rseci2 = rsecit;
                    }
                    else if (jd3 + jdf3 < jd1 + jdf1)
                    {
                        tdeclt = tdecl1;
                        trtasct = trtasc1;
                        jdt = jd1;
                        jdft = jdf1;
                        rsecit = rseci1;

                        tdecl1 = tdecl3;
                        trtasc1 = trtasc3;
                        jd1 = jd3;
                        jdf1 = jdf3;
                        rseci1 = rseci3;

                        tdecl3 = tdeclt;
                        trtasc3 = trtasct;
                        jd3 = jdt;
                        jdf3 = jdft;
                        rseci3 = rsecit;
                    }
                    // get latest one into 3 position
                    if (jd2 + jdf2 > jd3 + jdf3)
                    {
                        tdeclt = tdecl3;
                        trtasct = trtasc3;
                        jdt = jd3;
                        jdft = jdf3;
                        rsecit = rseci3;

                        tdecl3 = tdecl2;
                        trtasc3 = trtasc2;
                        jd3 = jd2;
                        jdf3 = jdf2;
                        rseci3 = rseci2;

                        tdecl2 = tdeclt;
                        trtasc2 = trtasct;
                        jd2 = jdt;
                        jdf2 = jdft;
                        rseci2 = rsecit;
                    }


                    double rng1, rng2, pctchg;
                    string ltrstr;
                    rng1 = 25000.0;  // old 12500 needs to be in km!! seems to do better when all the same? if too far off (*2) NAN
                    rng2 = 25000.0 * 1.02;  // 1.02 might be better? make the initial guess a bit different
                    //rng3 = 25000.0 * 1.08;
                    pctchg = 0.10;

                    strbuild.AppendLine((jd1 + jdf1).ToString() + " " + (jd2 + jdf2).ToString() + " " + (jd3 + jdf3).ToString());
                    if (false)
                    {
                        AstroLibr.anglesgauss(tdecl1, tdecl2, tdecl3, trtasc1, trtasc2, trtasc3,
                             jd1, jdf1, jd2, jdf2, jd3, jdf3,
                             rseci1, rseci2, rseci3, out r2, out v2, out errstr);

                        ltrstr = errstr.Substring(errstr.Length - 2, 2);
                        AstroLibr.rv2coe(r2, v2, out p, out a, out ecc, out incl, out raan, out argp, out nu, out m,
                            out arglat, out truelon, out lonper);
                        coestr = "a= " + a.ToString("0.0000") + " e= " + ecc.ToString("0.000000000") + " i= " + (incl * rad).ToString("0.0000")
                            + "  " + (raan * rad).ToString("0.0000") + "  " + (argp * rad).ToString("0.0000") + "  "
                            + (nu * rad).ToString("0.0000") + "  " + (m * rad).ToString("0.0000") + "  " + (arglat * rad).ToString("0.0000");
                        this.coeDR.Text = coestr;
                        strbuild.AppendLine("Gauss case:" + ktr + " " + ktr2 + " " + ktr3 + " iters: " + ltrstr + " coes " + coestr);
                    }

                    if (true)
                    {
                        AstroLibr.anglesdoubler(tdecl1, tdecl2, tdecl3, trtasc1, trtasc2, trtasc3,
                             jd1, jdf1, jd2, jdf2, jd3, jdf3,
                             rseci1, rseci2, rseci3, rng1, rng2, out r2, out v2, out errstr, pctchg);

                        ltrstr = errstr.Substring(errstr.Length - 2, 2);
                        AstroLibr.rv2coe(r2, v2, out p, out a, out ecc, out incl, out raan, out argp, out nu, out m,
                            out arglat, out truelon, out lonper);
                        coestr = "a= " + a.ToString("0.0000") + " e= " + ecc.ToString("0.000000000") + " i= " + (incl * rad).ToString("0.0000")
                            + "  " + (raan * rad).ToString("0.0000") + "  " + (argp * rad).ToString("0.0000") + "  "
                            + (nu * rad).ToString("0.0000") + "  " + (m * rad).ToString("0.0000") + "  " + (arglat * rad).ToString("0.0000");
                        this.coeDR.Text = coestr;
                        strbuild.AppendLine("Double-r case:" + ktr + " " + ktr2 + " " + ktr3 + " iters: " + ltrstr + " coes " + coestr);
                    }
                }
            }


            directory = @"D:\Codes\LIBRARY\cs\TestAll\";
            File.WriteAllText(directory + "testall-Anglessum.out", strbuild.ToString());

            this.opsStatusG.Text = "Status:  done";
            Refresh();
        }



    }  // class

}  // namespace

