// ----------------------------------------------------------------------------
//
//                                   EOPSPWLib
//
// this library contains various eop and space weather routines. These
// methods read the combined files and find specific values given a particular 
// date. 
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
//     EOPdataClass - EOP class parameters
//     eopdata      - array for EOP data 
//     eopdataU     - array for EOP data from USNO
//     SPWdataClass - SPW class parameters
//     spwdata      - array for SPW data 
//     iau80Class   - array for iau80 parameters
//     iau06Class   - array for iau06 parameters
//
// ----------------------------------------------------------------------------

using System;
using System.Linq;
using System.Text.RegularExpressions;
using System.IO;

using MathTimeMethods;  // Edirection, globals

namespace EOPSPWMethods
{
    public class EOPSPWLib
    {
        public string EOPSPWLibVersion = "EOPSPWLib Version 2025-01-20";

        // assume max of 24000 days from 1957 and eop from 1962  
        // include extra for 10 years of daily predictions
        public const long numb = 60000;
        // actual number of records read in
        public long numbeop, numbspw;

        // this is just the EOP data setup so it can grow
        public class EOPdataClass
        {
            public double xp, xerr, yp, yerr, dut1, dut1err, lod, ddpsi, ddeps, dx, dy;
            public int year, mon, day, mjd, dat;
        }
        // use a class for all the eop data so it can be processed
        public EOPdataClass[] eopdata = new EOPdataClass[numb];
        // temp for USNO processing
        public EOPdataClass[] eopdataU = new EOPdataClass[numb];


        // set it up so it can grow
        public class SPWdataClass
        {
            public double adjf10, adjctrf81, adjlstf81, Dadjf10, obsf10, obsctrf81, obslstf81, cp;
            public double adjf30, obsf30, adjlstf30f81, obslstf30f81;
            public Int32 mjd;
            public Int32 yr, mon, day, avgap, sumkp, c9;
            public Int32 brsn, nd, isn, q;
            public Int32[] aparr = new int[8];
            public Int32[] kparr = new int[8];
        }
        // use a struct for all the space weather data so it can be processed
        public SPWdataClass[] spwdata = new SPWdataClass[numb];

        // this will be a fixed size
        public class iau80Class
        {
            public double[,] rar80 = new double[106, 4];
            public Int32[,] iar80 = new int[106, 5];
        }
        public iau80Class iau80arr = new iau80Class();

        // this will be a fixed size since the constants stay the same (mostly)
        public class iau06Class
        {
            public double[,] ax0 = new double[1600, 2];  // reals
            public Int32[,] ax0i = new int[1600, 14];  // integers
            public double[,] ay0 = new double[1275, 2];  // reals
            public Int32[,] ay0i = new int[1275, 14];  // integers
            public double[,] as0 = new double[66, 2];  // reals
            public Int32[,] as0i = new int[66, 14];  // integers

            public double[,] ag0 = new double[35, 2];  // reals
            public Int32[,] ag0i = new int[35, 14];  // integers

            public double[,] apn0 = new double[1358, 2];  // reals
            public Int32[,] apn0i = new int[1358, 14];  // integers
            public double[,] apl0 = new double[1056, 2];  // reals
            public Int32[,] apl0i = new int[1056, 14];  // integers

            public double[,] aapn0 = new double[678, 8];  // reals
            public Int32[,] aapn0i = new int[678, 5];  // integers
            public double[,] aapl0 = new double[687, 5];  // reals
            public Int32[,] aapl0i = new int[687, 14];  // integers
        }
        public iau06Class iau06arr = new iau06Class();


        // setup the class so methods can be called
        public MathTimeLib MathTimeLibr = new MathTimeLib();



        // ---------------------------------------------------------------------------
        //
        //                           function iau80in
        //
        //  this function initializes the nutation matricies needed for reduction
        //    calculations. the routine needs the filename of the files as input.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    eoploc      - file location of eop
        //    
        //  outputs       :
        //    iau80arr    - record containing the iau80 constants rad
        //
        //  locals        :
        //    convrt      - conversion factor to degrees
        //    i,j         - index
        //
        //  coupling      :
        //    none        -
        // --------------------------------------------------------------------------- 

        public void iau80in(string nutLoc, out iau80Class iau80arr)
        {
            double convrt;
            int i, j;
            //string pattern;
            string pattern;
            // string line;

            // ------------------------  implementation   -------------------
            convrt = 0.0001 * Math.PI / (180.0 * 3600.0); // 0.0001" to rad

            iau80arr = new iau80Class();

            string[] fileData = File.ReadAllLines(nutLoc);

            pattern = @"^\s?(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)";
            for (i = 0; i < 106; i++)
            {
                //line = Regex.Replace(fileData[i], @"\s+", "|");
                //string[] linedata = line.Split('|');

                string[] linedata = Regex.Split(fileData[i], pattern);

                iau80arr.iar80[i, 0] = Convert.ToInt32(linedata[1]);
                iau80arr.iar80[i, 1] = Convert.ToInt32(linedata[2]);
                iau80arr.iar80[i, 2] = Convert.ToInt32(linedata[3]);
                iau80arr.iar80[i, 3] = Convert.ToInt32(linedata[4]);
                iau80arr.iar80[i, 4] = Convert.ToInt32(linedata[5]);
                iau80arr.rar80[i, 0] = Convert.ToDouble(linedata[6]);
                iau80arr.rar80[i, 1] = Convert.ToDouble(linedata[7]);
                iau80arr.rar80[i, 2] = Convert.ToDouble(linedata[8]);
                iau80arr.rar80[i, 3] = Convert.ToDouble(linedata[9]);

                for (j = 0; j < 4; j++)
                    iau80arr.rar80[i, j] = iau80arr.rar80[i, j] * convrt;
            }
        }  // iau80in


        // ---------------------------------------------------------------------------
        //
        //                           function iau06in
        //
        //  this function initializes the matricies needed for iau 2006 reduction
        //    calculations. the routine uses the files listed as inputs, but they are
        //    are not input to the routine as they are static files. 
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    iau06xtab5.2.a.dat  - file for x coefficient
        //    iau06ytab5.2.b.dat  - file for y coefficient
        //    iau06stab5.2.d.dat  - file for s coefficient
        //    iau06ansofa.dat     - file for nutation coefficients
        //    iau06anpl.dat       - file for planetary nutation coefficients
        //    iau06nlontab5.3.a.dat - file for longitude coefficients
        //    iau06nobltab5.3.b.dat - file for obliquity coefficients
        //    iau06gsttab5.3.e.dat  - file for gmst coefficients
        //
        //  outputs       :
        //    a0x         - real coefficients for x                     rad
        //    a0xi        - integer coefficients for x
        //    a0y         - real coefficients for y                     rad
        //    a0yi        - integer coefficients for y
        //    a0s         - real coefficients for s                     rad
        //    a0si        - integer coefficients for s
        //    apn         - real coefficients for nutation              rad
        //    apni        - integer coefficients for nutation
        //    ape         - real coefficients for obliquity             rad
        //    apei        - integer coefficients for obliquity
        //    agst        - real coefficients for gst                   rad
        //    agsti       - integer coefficients for gst
        //
        //  locals        :
        //    convrt      - conversion factor to radians
        //    i           - index
        //
        //  coupling      :
        //    none        -
        //
        //  references    :
        //    vallado     2013, pg 205-219, 910-912
        // ----------------------------------------------------------------------------- */

        public void iau06in(string nutLoc, out iau06Class iau06arr)
        {
            string pattern;
            // string line;
            Int32 i, j;
            // Int32 k, numsegs, ktr;

            // ------------------------  implementation   -------------------
            // " to rad
            double convrtu = (0.000001 * Math.PI) / (180.0 * 3600.0);  // if micro arcsecond
            double convrtmu = (0.0000001 * Math.PI) / (180.0 * 3600.0);  // if 0.1 micro arcsecond
            //double convrtm = (0.001 * Math.PI) / (180.0 * 3600.0);     // if milli arcsecond

            // ------------------------------
            //  note that since all these coefficients have only a single decimal place, one could store them as integers, and then simply
            //  divide by one additional power of ten. it would make memory storage much smaller and potentially faster.
            iau06arr = new iau06Class();

            // xys values
            // tab5.2a.txt in IERS
            // find first data point, and number for the segment
            //ktr = 0;  // ktr going through the file
            //k = 0;    // ktr through array of values
            //while (ktr < fileData.Count())
            //{
            //    while (!(fileData[ktr][0] == 'j'))
            //        ktr = ktr + 1;
            //    numsegs = Convert.ToInt32(fileData[ktr].Split('=')[2]);
            //    ktr = ktr + 2;
            //    for (i = 0; i < numsegs; i++)
            //    {
            //        //line = Regex.Replace(fileData[i], @"\s+", "|");
            //        //string[] linedata = line.Split('|');
            //        // reals
            //        string[] linedata = Regex.Split(fileData[ktr], pattern);
            //        iau06arr.axs0[k, 0] = Convert.ToDouble(linedata[2]) * convrtu;  // rad
            //        iau06arr.axs0[k, 1] = Convert.ToDouble(linedata[3]) * convrtu;  // rad
            //                                                                        // integers
            //        for (j = 0; j < 13; j++)
            //            iau06arr.a0xi[k, j] = Convert.ToInt32(linedata[j + 4]);
            //        k = k + 1;
            //        ktr = ktr + 1;
            //    }
            //}

            pattern = @"^\s+?(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)";
            // iau06x.txt in IERS
            string[] fileData = File.ReadAllLines(nutLoc + "iau06xtab5.2.a.dat");
            // start on line 2 to be below header info
            for (i = 2; i < fileData.Count(); i++)
            {
                //line = Regex.Replace(fileData[i], @"\s+", "|");
                //string[] linedata = line.Split('|');
                // reals
                string[] linedata = Regex.Split(fileData[i], pattern);
                iau06arr.ax0[i - 2, 0] = Convert.ToDouble(linedata[2]) * convrtu;  // rad
                iau06arr.ax0[i - 2, 1] = Convert.ToDouble(linedata[3]) * convrtu;  // rad
                // integers
                for (j = 0; j <= 13; j++)
                    iau06arr.ax0i[i - 2, j] = Convert.ToInt32(linedata[j + 4]);
            }

            // tab5.2b.txt in IERS
            fileData = File.ReadAllLines(nutLoc + "iau06ytab5.2.b.dat");
            for (i = 2; i < fileData.Count(); i++)
            {
                //line = Regex.Replace(fileData[i], @"\s+", "|");
                //string[] linedata = line.Split('|');
                // reals
                string[] linedata = Regex.Split(fileData[i], pattern);
                iau06arr.ay0[i - 2, 0] = Convert.ToDouble(linedata[2]) * convrtu;  // rad
                iau06arr.ay0[i - 2, 1] = Convert.ToDouble(linedata[3]) * convrtu;  // rad
                // integers
                for (j = 0; j <= 13; j++)
                    iau06arr.ay0i[i - 2, j] = Convert.ToInt32(linedata[j + 4]);
            }

            // tab5.2d.txt in IERS
            fileData = File.ReadAllLines(nutLoc + "iau06stab5.2.d.dat");
            for (i = 2; i < fileData.Count(); i++)
            {
                //line = Regex.Replace(fileData[i], @"\s+", "|");
                //string[] linedata = line.Split('|');
                // reals
                string[] linedata = Regex.Split(fileData[i], pattern);
                // reals
                iau06arr.as0[i - 2, 0] = Convert.ToDouble(linedata[2]) * convrtu;  // rad
                iau06arr.as0[i - 2, 1] = Convert.ToDouble(linedata[3]) * convrtu;  // rad
                // integers
                for (j = 0; j <= 13; j++)
                    iau06arr.as0i[i - 2, j] = Convert.ToInt32(linedata[j + 4]);
            }


            string pattern1 = @"^\s+?(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)";
            //     2000a nutation values old approach iau2000a
            fileData = File.ReadAllLines(nutLoc + "iau06ansofa.dat");  // 678
            for (i = 2; i < fileData.Count(); i++)
            {
                //line = Regex.Replace(fileData[i], @"\s+", "|");
                //string[] linedata = line.Split('|');
                // reals
                string[] linedata = Regex.Split(fileData[i], pattern1);
                // reals
                for (j = 0; j < 6; j++)
                    iau06arr.aapn0[i - 2, j] = Convert.ToDouble(linedata[j + 6]) * convrtmu;  // rad
                // integers
                for (j = 0; j < 5; j++)
                    iau06arr.aapn0i[i - 2, j] = Convert.ToInt32(linedata[j + 1]);
            }

            //
            pattern1 = @"^\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)";
            //     2000a planetary nutation values
            fileData = File.ReadAllLines(nutLoc + "iau06anplsofa.dat");  // 687
            for (i = 2; i < fileData.Count(); i++)
            {
                //line = Regex.Replace(fileData[i], @"\s+", "|");
                //string[] linedata = line.Split('|');
                string[] linedata = Regex.Split(fileData[i], pattern1);
                // reals
                for (j = 0; j < 4; j++)
                    iau06arr.aapl0[i - 2, j] = Convert.ToDouble(linedata[j + 15]) * convrtmu;  // rad
                // integers
                for (j = 0; j < 14; j++)   // 13 ??
                    iau06arr.aapl0i[i - 2, j] = Convert.ToInt32(linedata[j + 1]);
            }

            // tab5.3a.txt in IERS
            // nutation values planetary now included new iau2006
            // nutation in longitude
            fileData = File.ReadAllLines(nutLoc + "iau06nlontab5.3.a.dat");  //1358
            for (i = 2; i < fileData.Count(); i++)
            {
                //line = Regex.Replace(fileData[i], @"\s+", "|");
                //string[] linedata = line.Split('|');
                // reals
                string[] linedata = Regex.Split(fileData[i], pattern);
                iau06arr.apn0[i - 2, 0] = Convert.ToDouble(linedata[2]) * convrtu;  // rad
                iau06arr.apn0[i - 2, 1] = Convert.ToDouble(linedata[3]) * convrtu;  // rad
                // integers
                for (j = 0; j <= 13; j++)
                    iau06arr.apn0i[i - 2, j] = Convert.ToInt32(linedata[j + 4]);
            }

            // tab5.3b.txt in IERS
            // nutation in obliquity
            fileData = File.ReadAllLines(nutLoc + "iau06nobltab5.3.b.dat");  // 1056
            for (i = 2; i < fileData.Count(); i++)
            {
                //line = Regex.Replace(fileData[i], @"\s+", "|");
                //string[] linedata = line.Split('|');
                // reals
                string[] linedata = Regex.Split(fileData[i], pattern);
                iau06arr.apl0[i - 2, 0] = Convert.ToDouble(linedata[2]) * convrtu;  // rad
                iau06arr.apl0[i - 2, 1] = Convert.ToDouble(linedata[3]) * convrtu;  // rad
                // integers
                for (j = 0; j <= 13; j++)
                    iau06arr.apl0i[i - 2, j] = Convert.ToInt32(linedata[j + 4]);
            }

            // tab5.2e.txt in IERS
            // gmst values
            // note - these are very similar to the first 34 elements of iau06s.dat,
            // but they are not the same.
            fileData = File.ReadAllLines(nutLoc + "iau06gsttab5.2.e.dat");
            for (i = 2; i < fileData.Count(); i++)
            {
                //line = Regex.Replace(fileData[i], @"\s+", "|");
                //string[] linedata = line.Split('|');
                // reals
                string[] linedata = Regex.Split(fileData[i], pattern);
                iau06arr.ag0[i - 2, 0] = Convert.ToDouble(linedata[2]) * convrtu;  // rad
                iau06arr.ag0[i - 2, 1] = Convert.ToDouble(linedata[3]) * convrtu;  // rad
                // integers
                for (j = 0; j <= 13; j++)
                    iau06arr.ag0i[i - 2, j] = Convert.ToInt32(linedata[j + 4]);
            }

        }    // iau06in



        // ---------------------------------------------------------------------------
        //
        //                           function readeop
        //
        //  this function initializes the earth orientation parameter data with an EOP file. the input
        //  data files are from CelesTrak and the eoppn.txt file contains the nutation
        //  daily values used for optimizing the speed of operation. note these nutation
        //  values do not have corrections applied (dx/dy/ddpsi/ddeps) so a single
        //  routine can process past and future data with these corrections. the
        //  corrections are not known in the future. the files could be combined, but it
        //  may make more sense to keep them separate because the values can be calculated
        //  long into the future.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    eopdata[0]  - array of eop data records
        //    infile      - name of eop file to read
        //     
        //  outputs       :
        //    eopdata[0]  - array of eop data records
        //    mjdeopstart - modified julian date of the start of the eopdata[0] data
        //    ktrActualObs- ktr of number of actual obs (no predicted)
        //    updDate     - date string of update in case it's not the current time
        //
        //  locals        :
        //                -
        //
        //  coupling      :
        // ---------------------------------------------------------------------------

        public void readeop(ref EOPdataClass[] eopdata,
                             string inFile, out Int32 ktrActualObs, out string updDate)
        {
            Int32 numrecsobs, i, ktr;
            string[] linedata;

            updDate = "";

            // eiter do the whole array at once, or as each line is read
            //            initEOPArray(ref eopdata);

            // read the whole file at once into lines of an array
            string[] EOParray = File.ReadAllLines(inFile);

            // BEGIN OBSERVED
            // #   Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
            // # (0h UTC)           "         "          s          s          "        "          "         "     s 
            // 1962 01 01 37665 -0.012700  0.213000  0.0326338  0.0017230  0.064041  0.006305  0.000000  0.000000   2 
            // 1962 01 02 37666 -0.015900  0.214100  0.0320547  0.0016690  0.063758  0.006529  0.000000  0.000000   2 
            // 1962 01 03 37667 -0.019000  0.215200  0.0315526  0.0015820  0.063649  0.006754  0.000000  0.000000   2 

            // find beginning of data
            i = 0;
            while (!EOParray[i].Contains("NUM_OBSERVED_POINTS"))
            {
                i = i + 1;
                if (EOParray[i].Contains("UPDATED"))
                {
                    updDate = EOParray[i].Substring(8, 17);
                }
            }

            linedata = EOParray[i].Split(' ');
            numrecsobs = Convert.ToInt32(linedata[1]);  // starts at 0

            ktr = 0;
            // ---- process observed records only
            for (ktr = 0; ktr < numrecsobs; ktr++)
            {
                // set new record as they are needed
                eopdata[ktr] = new EOPdataClass();

                // replace multiple spaces with just one
                string line3 = Regex.Replace(EOParray[ktr + i + 2].ToString(), @"\s+", " ");
                linedata = line3.Split(' ');
                // do all at once?
                //eopdata = new EOPdataClass();

                //eopdata[ktr].year = Convert.ToInt32(linedata[0]);  // starts at 0
                //eopdata[ktr].mon = Convert.ToInt32(linedata[1]);
                //eopdata[ktr].day = Convert.ToInt32(linedata[2]);
                eopdata[ktr].mjd = Convert.ToInt32(linedata[3]);
                eopdata[ktr].xp = Convert.ToDouble(linedata[4]);
                //eopdata[ktr].xerr = 0.0;
                eopdata[ktr].yp = Convert.ToDouble(linedata[5]);
                //eopdata[ktr].yerr = 0.0;
                eopdata[ktr].dut1 = Convert.ToDouble(linedata[6]);
                //eopdata[ktr].dut1err = 0.0;
                eopdata[ktr].lod = Convert.ToDouble(linedata[7]);
                eopdata[ktr].ddpsi = Convert.ToDouble(linedata[8]);
                eopdata[ktr].ddeps = Convert.ToDouble(linedata[9]);
                eopdata[ktr].dx = Convert.ToDouble(linedata[10]);
                eopdata[ktr].dy = Convert.ToDouble(linedata[11]);
                eopdata[ktr].dat = Convert.ToInt32(linedata[12]);

            }

            ktrActualObs = Convert.ToInt32(ktr);  // skip end processed as they can have errors (do later)

            // ---- process predicted records
            i = i + ktr;  // start from end of last read
            while (!EOParray[i].Contains("NUM_PREDICTED_POINTS"))
            {
                i = i + 1;
            }

            linedata = EOParray[i].Split(' ');
            numrecsobs = Convert.ToInt32(linedata[1]);  // starts at 0

            // ---- process observed records only
            for (ktr = 0; ktr < numrecsobs; ktr++)
            {
                // replace multiple spaces with just one
                string line3 = Regex.Replace(EOParray[ktr + i + 2].ToString(), @"\s+", " ");
                linedata = line3.Split(' ');
                Int32 ktr1 = ktrActualObs + ktr;
                // set new record as they are needed
                eopdata[ktr1] = new EOPdataClass();

                //eopdata[ktr1].year = Convert.ToInt32(linedata[0]);  // starts at 0
                //eopdata[ktr1].mon = Convert.ToInt32(linedata[1]);
                //eopdata[ktr1].day = Convert.ToInt32(linedata[2]);
                eopdata[ktr1].mjd = Convert.ToInt32(linedata[3]);
                eopdata[ktr1].xp = Convert.ToDouble(linedata[4]);
                //eopdata[ktr1].xerr = 0.0;
                eopdata[ktr1].yp = Convert.ToDouble(linedata[5]);
                //eopdata[ktr1].yerr = 0.0;
                eopdata[ktr1].dut1 = Convert.ToDouble(linedata[6]);
                //eopdata[ktr1].dut1err = 0.0;
                eopdata[ktr1].lod = Convert.ToDouble(linedata[7]);
                eopdata[ktr1].ddpsi = Convert.ToDouble(linedata[8]);
                eopdata[ktr1].ddeps = Convert.ToDouble(linedata[9]);
                eopdata[ktr1].dx = Convert.ToDouble(linedata[10]);
                eopdata[ktr1].dy = Convert.ToDouble(linedata[11]);
                eopdata[ktr1].dat = Convert.ToInt32(linedata[12]);
            }

            ktrActualObs = ktrActualObs + Convert.ToInt32(ktr);  // skip end proccessed as they can have errors (do later)
            numbeop = ktrActualObs;

        }   // readeop


        // ---------------------------------------------------------------------------
        //
        //                           function findeopparam
        //
        //  this routine finds the eop parameters for a given time. several types of
        //  interpolation are available. 
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    jd          - julian date of epoch (0 hrs utc)
        //    jdFrac      - fractional date of epoch
        //    interp      - interpolation                                n-none, l-linear, s-spline
        //    eopdata[0]  - array of eop data records
        //    jdeopstart  - julian date of the start of the eopdata[0] data (set in initeop)
        //
        //  outputs       :
        //    dut1        - delta ut1 (ut1-utc)                               sec
        //    dat         - number of leap seconds                            sec
        //    lod         - excess length of day                              sec
        //    xp          - x component of polar motion                       rad
        //    yp          - y component of polar motion                       rad
        //    ddpsi       - correction to delta psi (iau-76 theory)           rad
        //    ddeps       - correction to delta eps (iau-76 theory)           rad
        //    dx          - correction to x (cio theory)                      rad
        //    dy          - correction to y (cio theory)                      rad
        //
        //  locals        :
        //                -
        //
        //  coupling      :
        //    none        -
        // ---------------------------------------------------------------------------

        public void findeopparam(double jd, double jdFrac, char interp, EOPdataClass[] eopdata,
               out double dut1, out int dat, out double lod, out double xp, out double yp,
               out double ddpsi, out double ddeps, out double dx, out double dy)
        {
            Int32 recnum;
            int off1, off2;
            //eopdata eopdata[0], lasteoprec, tempeoprec;
            // EOPdataClass[] nexteoprec;
            double fixf, mfme, jd1, jdeopstarto;
            double convrt = Math.PI / (3600.0 * 180.0);  // " to rad

            // check if any whole days in jdF
            jd1 = Math.Floor(jd + jdFrac) + 0.5;  // want jd at 0 hr
            mfme = (jd + jdFrac - jd1) * 1440.0;
            if (mfme < 0.0)
                mfme = 1440.0 + mfme;

            // ---- read data for day of interest
            jdeopstarto = Math.Floor(jd + jdFrac - eopdata[0].mjd - 2400000.5); // needed to get correct start day
            recnum = Convert.ToInt32(jdeopstarto);

            // check for out of bound values
            if ((recnum >= 1) && (recnum < numbeop - 1))
            {
                // ---- set non-interpolated values
                dut1 = eopdata[recnum].dut1;
                dat = eopdata[recnum].dat;
                lod = eopdata[recnum].lod;
                xp = eopdata[recnum].xp;
                yp = eopdata[recnum].yp;
                ddpsi = eopdata[recnum].ddpsi;
                ddeps = eopdata[recnum].ddeps;
                dx = eopdata[recnum].dx;
                dy = eopdata[recnum].dy;

                // ---- do linear interpolation
                if (interp == 'l')
                {
                    fixf = mfme / 1440.0;

                    dut1 = eopdata[recnum].dut1 + (eopdata[recnum + 1].dut1 - eopdata[recnum].dut1) * fixf;
                    lod = eopdata[recnum].lod + (eopdata[recnum + 1].lod - eopdata[recnum].lod) * fixf;
                    xp = eopdata[recnum].xp + (eopdata[recnum + 1].xp - eopdata[recnum].xp) * fixf;
                    yp = eopdata[recnum].yp + (eopdata[recnum + 1].yp - eopdata[recnum].yp) * fixf;
                    ddpsi = eopdata[recnum].ddpsi + (eopdata[recnum + 1].ddpsi - eopdata[recnum].ddpsi) * fixf;
                    ddeps = eopdata[recnum].ddeps + (eopdata[recnum + 1].ddeps - eopdata[recnum].ddeps) * fixf;
                    dx = eopdata[recnum].dx + (eopdata[recnum + 1].dx - eopdata[recnum].dx) * fixf;
                    dy = eopdata[recnum].dy + (eopdata[recnum + 1].dy - eopdata[recnum].dy) * fixf;
                }

                // ---- do spline interpolations
                if (interp == 's')
                {
                    fixf = mfme / 1440.0;  // since using mjd as time argument
                    off1 = 1;
                    off2 = 2;
                    // setup so the interval is in between points 2 and 3
                    dut1 = MathTimeLibr.cubicinterp(eopdata[recnum - off1].dut1, eopdata[recnum].dut1, eopdata[recnum + off1].dut1, eopdata[recnum + off2].dut1,
                                         eopdata[recnum - off1].mjd, eopdata[recnum].mjd, eopdata[recnum + off1].mjd, eopdata[recnum + off2].mjd,
                                         eopdata[recnum].mjd + fixf);
                    lod = MathTimeLibr.cubicinterp(eopdata[recnum - off1].lod, eopdata[recnum].lod, eopdata[recnum + off1].lod, eopdata[recnum + off2].lod,
                                         eopdata[recnum - off1].mjd, eopdata[recnum].mjd, eopdata[recnum + off1].mjd, eopdata[recnum + off2].mjd,
                                         eopdata[recnum].mjd + fixf);
                    xp = MathTimeLibr.cubicinterp(eopdata[recnum - off1].xp, eopdata[recnum].xp, eopdata[recnum + off1].xp, eopdata[recnum + off2].xp,
                                         eopdata[recnum - off1].mjd, eopdata[recnum].mjd, eopdata[recnum + off1].mjd, eopdata[recnum + off2].mjd,
                                         eopdata[recnum].mjd + fixf);
                    yp = MathTimeLibr.cubicinterp(eopdata[recnum - off1].yp, eopdata[recnum].yp, eopdata[recnum + off1].yp, eopdata[recnum + off2].yp,
                                         eopdata[recnum - off1].mjd, eopdata[recnum].mjd, eopdata[recnum + off1].mjd, eopdata[recnum + off2].mjd,
                                         eopdata[recnum].mjd + fixf);
                    ddpsi = MathTimeLibr.cubicinterp(eopdata[recnum - off1].ddpsi, eopdata[recnum].ddpsi, eopdata[recnum + off1].ddpsi, eopdata[recnum + off2].ddpsi,
                                         eopdata[recnum - off1].mjd, eopdata[recnum].mjd, eopdata[recnum + off1].mjd, eopdata[recnum + off2].mjd,
                                         eopdata[recnum].mjd + fixf);
                    ddeps = MathTimeLibr.cubicinterp(eopdata[recnum - off1].ddeps, eopdata[recnum].ddeps, eopdata[recnum + off1].ddeps, eopdata[recnum + off2].ddeps,
                                         eopdata[recnum - off1].mjd, eopdata[recnum].mjd, eopdata[recnum + off1].mjd, eopdata[recnum + off2].mjd,
                                         eopdata[recnum].mjd + fixf);
                    dx = MathTimeLibr.cubicinterp(eopdata[recnum - off1].dx, eopdata[recnum].dx, eopdata[recnum + off1].dx, eopdata[recnum + off2].dx,
                                         eopdata[recnum - off1].mjd, eopdata[recnum].mjd, eopdata[recnum + off1].mjd, eopdata[recnum + off2].mjd,
                                         eopdata[recnum].mjd + fixf);
                    dy = MathTimeLibr.cubicinterp(eopdata[recnum - off1].dy, eopdata[recnum].dy, eopdata[recnum + off1].dy, eopdata[recnum + off2].dy,
                                         eopdata[recnum - off1].mjd, eopdata[recnum].mjd, eopdata[recnum + off1].mjd, eopdata[recnum + off2].mjd,
                                         eopdata[recnum].mjd + fixf);
                }
            }
            // set default values
            else
            {
                dut1 = 0.0;
                dat = 37;
                lod = 0.0;
                xp = 0.0;
                yp = 0.0;
                ddpsi = 0.0;
                ddeps = 0.0;
                dx = 0.0;
                dy = 0.0;
            }

            // now convert units for use in operations
            xp = xp * convrt;  // " to rad
            yp = yp * convrt;
            ddpsi = ddpsi * convrt;  // " to rad
            ddeps = ddeps * convrt;
            dx = dx * convrt;  // " to rad
            dy = dy * convrt;
        }  //  findeopparam




        // ---------------------------------------------------------------------------
        //
        //                           function readspw
        //
        //  this function initializes the space weather data from the cssi/SPW files. 
        //  It reads the actual data, and the predicted.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    inFile      - path and name of input file
        //    
        //  outputs       :
        //    spwdata     - array of spw data records
        //    jdspwstart  - julian date of the start of the spwarr data
        //
        //  locals        :
        //                -
        //
        //  coupling      :
        //    jday        - julian date
        // ---------------------------------------------------------------------------

        public void readspw(ref SPWdataClass[] spwdata, string inFile, out Int32 ktrActualObs, out string errstr)
        {
            double jd, jdFrac;
            Int32 oldjd;
            Int32 numrecsobs, i, ktr;

            //        initSPWArray(ref spwdata);
            oldjd = 0;
            errstr = "ok\n";
            string line3 = "";
            DateTime now = DateTime.UtcNow;
            MathTimeLibr.jday(now.Year, now.Month, now.Day, 0, 0, 0.0, out jd, out jdFrac);
            Int32 nowmjd = Convert.ToInt32(jd + jdFrac - 2400000.5);

            // read the whole file at once into lines of an array
            string[] SPWarray = File.ReadAllLines(inFile);

            // NUM_OBSERVED_POINTS 20591
            // BEGIN OBSERVED
            // 1957 10 01 1700 19 43 40 30 20 37 23 43 37 273  32  27  15   7  22   9  32  22  21 1.1 5 236 268.0 0 265.2 230.6 269.3 266.6 230.9 1.9 1.8 1.7 1.6
            // 1957 10 02 1700 20 37 37 17 17 27 23 17 30 203  22  22   6   6  12   9   6  15  12 0.7 3 234 252.0 0 266.0 231.4 253.3 267.4 231.7 1.9 1.8 1.7 1.6 
            //       Regex data = new Regex(@"(\d+):\s*(\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)");

            // find beginning of data
            i = 0;
            while (!SPWarray[i].Contains("NUM_OBSERVED_POINTS"))
            {
                i = i + 1;
            }

            string[] linedata = SPWarray[i].Split(' ');
            numrecsobs = Convert.ToInt32(linedata[1]);  // starts at 0

            ktr = 0;
            // ---- process observed records
            // use 15 to avoid last minute gaps in some data - better is to use regex and match success !!!!!!!!!!!!!!!!!
            // update to read from fixed format with all data included
            for (ktr = 0; ktr < numrecsobs; ktr++)  // use  - 15 if old format of spw files
            {
                // set new record as they are needed
                spwdata[ktr] = new SPWdataClass();

                // replace multiple spaces with just one
                try
                {
                    line3 = Regex.Replace(SPWarray[ktr + i + 2].ToString(), @"\s+", " ");
                    linedata = line3.Split(' ');
                    if (line3.Contains("-1"))
                        errstr = errstr + "Error " + line3 + "\n";

                    //spwdata[ktr].yr = Convert.ToInt32(linedata[0]);  // starts at 0
                    //spwdata[ktr].mon = Convert.ToInt32(linedata[1]);
                    //spwdata[ktr].day = Convert.ToInt32(linedata[2]);
                    //spwdata[ktr].brsn = Convert.ToInt32(linedata[3]);
                    //spwdata[ktr].nd = Convert.ToInt32(linedata[4]);
                    spwdata[ktr].kparr[0] = Convert.ToInt32(linedata[5]);
                    spwdata[ktr].kparr[1] = Convert.ToInt32(linedata[6]);
                    spwdata[ktr].kparr[2] = Convert.ToInt32(linedata[7]);
                    spwdata[ktr].kparr[3] = Convert.ToInt32(linedata[8]);
                    spwdata[ktr].kparr[4] = Convert.ToInt32(linedata[9]);
                    spwdata[ktr].kparr[5] = Convert.ToInt32(linedata[10]);
                    spwdata[ktr].kparr[6] = Convert.ToInt32(linedata[11]);
                    spwdata[ktr].kparr[7] = Convert.ToInt32(linedata[12]);
                    spwdata[ktr].sumkp = Convert.ToInt32(linedata[13]);
                    spwdata[ktr].aparr[0] = Convert.ToInt32(linedata[14]);
                    spwdata[ktr].aparr[1] = Convert.ToInt32(linedata[15]);
                    spwdata[ktr].aparr[2] = Convert.ToInt32(linedata[16]);
                    spwdata[ktr].aparr[3] = Convert.ToInt32(linedata[17]);
                    spwdata[ktr].aparr[4] = Convert.ToInt32(linedata[18]);
                    spwdata[ktr].aparr[5] = Convert.ToInt32(linedata[19]);
                    spwdata[ktr].aparr[6] = Convert.ToInt32(linedata[20]);
                    spwdata[ktr].aparr[7] = Convert.ToInt32(linedata[21]);
                    spwdata[ktr].avgap = Convert.ToInt32(linedata[22]);
                    //spwdata[ktr].cp = Convert.ToDouble(linedata[23]);
                    //spwdata[ktr].c9 = Convert.ToInt32(linedata[24]);
                    //spwdata[ktr].isn = Convert.ToInt32(linedata[25]);
                    spwdata[ktr].adjf10 = Convert.ToDouble(linedata[26]);
                    //spwdata[ktr].q = Convert.ToInt32(linedata[27]);
                    spwdata[ktr].adjctrf81 = Convert.ToDouble(linedata[28]);
                    spwdata[ktr].adjlstf81 = Convert.ToDouble(linedata[29]);
                    spwdata[ktr].obsf10 = Convert.ToDouble(linedata[30]);
                    spwdata[ktr].obsctrf81 = Convert.ToDouble(linedata[31]);
                    spwdata[ktr].obslstf81 = Convert.ToDouble(linedata[32]);
                    //spwdata[ktr].adjf30 = Convert.ToDouble(linedata[33]);
                    //spwdata[ktr].obsf30 = Convert.ToDouble(linedata[34]);
                    //spwdata[ktr].adjlstf30f81 = Convert.ToDouble(linedata[35]);
                    //spwdata[ktr].obslstf30f81 = Convert.ToDouble(linedata[36]);
                    MathTimeLibr.jday(Convert.ToInt32(linedata[0]), Convert.ToInt32(linedata[1]), 
                        Convert.ToInt32(linedata[2]), 0, 0, 0.0, out jd, out jdFrac);
                    if (ktr > 0)
                        oldjd = spwdata[ktr - 1].mjd;
                    spwdata[ktr].mjd = Convert.ToInt32(jd + jdFrac - 2400000.5);

                    // find last 3 days data
                    if (spwdata[ktr].mjd >= nowmjd - 3)
                    {
                        errstr = errstr + spwdata[ktr].yr + " " + spwdata[ktr].mon + " " + spwdata[ktr].day
                            + spwdata[ktr].aparr[0] + " " + spwdata[ktr].aparr[1] + " " + spwdata[ktr].aparr[2] + " " + spwdata[ktr].aparr[3]
                            + spwdata[ktr].aparr[4] + " " + spwdata[ktr].aparr[5] + " " + spwdata[ktr].aparr[6] + " " + spwdata[ktr].aparr[7]
                            + spwdata[ktr].adjf10 + " " + spwdata[ktr].obsf10;
                    }

                    if (ktr > 0 && spwdata[ktr].mjd != oldjd + 1)
                        errstr = errstr + "Error day before missing " + line3 + "\n";
                }
                catch
                {
                    errstr = errstr + "Error # act obs wrong " + numrecsobs + "\n";
                } // catch if num of recs doesn't match

            }  // for through observed obs

            ktrActualObs = Convert.ToInt32(ktr);  // skip end proccessed as they can have errors (do later)
            numbspw = ktrActualObs;

            try
            {
                Int32 numpredObs = Convert.ToInt32(SPWarray[ktrActualObs + i + 4].ToString().Split(' ')[1]);

                // ---- check file for 45 day predicted values
                ktr = ktrActualObs - 1;
                for (int ktr1 = ktrActualObs + 6 + i; ktr1 < ktrActualObs + 4 + numpredObs + i; ktr1++)  // use  - 15 if old format of spw files
                {
                    ktr = ktr + 1;
                    // set new record as they are needed
                    spwdata[ktr] = new SPWdataClass();

                    // replace multiple spaces with just one
                    try
                    {
                        line3 = Regex.Replace(SPWarray[ktr1].ToString(), @"\s+", " ");
                        linedata = line3.Split(' ');
                        if (line3.Contains("-1"))
                            errstr = errstr + "Error " + line3 + "\n";

                        //spwdata[ktr].yr = Convert.ToInt32(linedata[0]);  // starts at 0
                        //spwdata[ktr].mon = Convert.ToInt32(linedata[1]);
                        //spwdata[ktr].day = Convert.ToInt32(linedata[2]);
                        //spwdata[ktr].brsn = Convert.ToInt32(linedata[3]);
                        //spwdata[ktr].nd = Convert.ToInt32(linedata[4]);
                        spwdata[ktr].kparr[0] = Convert.ToInt32(linedata[5]);
                        spwdata[ktr].kparr[1] = Convert.ToInt32(linedata[6]);
                        spwdata[ktr].kparr[2] = Convert.ToInt32(linedata[7]);
                        spwdata[ktr].kparr[3] = Convert.ToInt32(linedata[8]);
                        spwdata[ktr].kparr[4] = Convert.ToInt32(linedata[9]);
                        spwdata[ktr].kparr[5] = Convert.ToInt32(linedata[10]);
                        spwdata[ktr].kparr[6] = Convert.ToInt32(linedata[11]);
                        spwdata[ktr].kparr[7] = Convert.ToInt32(linedata[12]);
                        spwdata[ktr].sumkp = Convert.ToInt32(linedata[13]);
                        spwdata[ktr].aparr[0] = Convert.ToInt32(linedata[14]);
                        spwdata[ktr].aparr[1] = Convert.ToInt32(linedata[15]);
                        spwdata[ktr].aparr[2] = Convert.ToInt32(linedata[16]);
                        spwdata[ktr].aparr[3] = Convert.ToInt32(linedata[17]);
                        spwdata[ktr].aparr[4] = Convert.ToInt32(linedata[18]);
                        spwdata[ktr].aparr[5] = Convert.ToInt32(linedata[19]);
                        spwdata[ktr].aparr[6] = Convert.ToInt32(linedata[20]);
                        spwdata[ktr].aparr[7] = Convert.ToInt32(linedata[21]);
                        spwdata[ktr].avgap = Convert.ToInt32(linedata[22]);
                        //spwdata[ktr].cp = Convert.ToDouble(linedata[23]);
                        //spwdata[ktr].c9 = Convert.ToInt32(linedata[24]);
                        //spwdata[ktr].isn = Convert.ToInt32(linedata[25]);
                        spwdata[ktr].adjf10 = Convert.ToDouble(linedata[26]);
                        //spwdata[ktr].q = Convert.ToInt32(linedata[27]);
                        spwdata[ktr].adjctrf81 = Convert.ToDouble(linedata[28]);
                        spwdata[ktr].adjlstf81 = Convert.ToDouble(linedata[29]);
                        spwdata[ktr].obsf10 = Convert.ToDouble(linedata[30]);
                        spwdata[ktr].obsctrf81 = Convert.ToDouble(linedata[31]);
                        spwdata[ktr].obslstf81 = Convert.ToDouble(linedata[32]);
                        //spwdata[ktr].adjf30 = Convert.ToDouble(linedata[33]);
                        //spwdata[ktr].obsf30 = Convert.ToDouble(linedata[34]);
                        //spwdata[ktr].adjlstf30f81 = Convert.ToDouble(linedata[35]);
                        //spwdata[ktr].obslstf30f81 = Convert.ToDouble(linedata[36]);
                        MathTimeLibr.jday(Convert.ToInt32(linedata[0]), Convert.ToInt32(linedata[1]),
                            Convert.ToInt32(linedata[2]), 0, 0, 0.0, out jd, out jdFrac);
                        oldjd = spwdata[ktr - 1].mjd;
                        spwdata[ktr].mjd = Convert.ToInt32(jd + jdFrac - 2400000.5);

                        if (spwdata[ktr].mjd != oldjd + 1 && ktr < ktrActualObs + 3 + numpredObs)
                            errstr = errstr + "Error day before missing " + line3 + "\n";

                        // find next 3 days data
                        if (spwdata[ktr].mjd <= nowmjd + 3)
                        {
                            errstr = errstr + spwdata[ktr].yr + " " + spwdata[ktr].mon + " " + spwdata[ktr].day
                                + spwdata[ktr].aparr[0] + " " + spwdata[ktr].aparr[1] + " " + spwdata[ktr].aparr[2] + " " + spwdata[ktr].aparr[3]
                                + spwdata[ktr].aparr[4] + " " + spwdata[ktr].aparr[5] + " " + spwdata[ktr].aparr[6] + " " + spwdata[ktr].aparr[7]
                                + spwdata[ktr].adjf10 + " " + spwdata[ktr].obsf10;
                        }
                    }
                    catch
                    {
                        errstr = errstr + "Error # act obs wrong " + numpredObs + "\n";
                    }

                }  // for

            }
            catch
            {
                errstr = errstr + "Error # act obs wrong " + numrecsobs + "\n";
            }

        }   // readspw



        // ---------------------------------------------------------------------------
        //
        //                           function findspwparam
        //
        //  this routine finds the atmospheric parameters for a given time.
        //    ap/kp 3 hourly data is valid at 0000, 0300 hrs, etc
        //    apavg and kpsum are valid at 1200 hrs
        //    f107 and f107bar values are valid at 1700/2000 hrs depending on the date
        //    ap arrays go 0-7, but msisarr goes 1-8 to match msis code and fortran source
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    jd          - julian date of epoch (0 hrs utc)              days from 4713 bc
        //    mfme        - minutes from midnight epoch                       mins
        //    interp      - interpolation                                 l - linear, s - spline
        //    fluxtype    - flux type                                     a-adjusted, o-observed
        //    f81type     - flux 81-day avg type                          l-last, c-centered
        //    inputtype   - input type                                    a-actual, c - constant
        //    spwarr      - array of space weather data
        //
        //  outputs       :
        //    f107        - f10.7 value (current day)
        //    f107bar     - f10.7 81-day avg value
        //    ap          - planetary effect array
        //    avgap       - daily average ap planetary value
        //    aparr       - last 8 values of 3-hourly ap
        //    kp          - planetary effect array
        //    sumkp       - daily kp sum of planetary values
        //    kparr       - last 8 values of 3-hourly kp
        //
        //  locals        :
        //    fluxtime    - minutes from midnight where f107 is valid (1020 or 1200)
        //
        //  coupling      :
        //    none        -
        // ---------------------------------------------------------------------------

        public void findspwparam(double jd, double jdFrac, char interp, char fluxtype, char f81type, char inputtype,
                                  SPWdataClass[] spwdata,
                                  out double f107, out double f107bar,
                                  out double ap, out double avgap, double[] aparr, out double kp, out double sumkp,
                                  double[] kparr)
        {
            double jd1, mfme, jdspwstarto, fluxtime, fixf, tf107, tf107bar;
            Int32 recnum, recnumt, idx, i, j, off1, off2;

            f107 = 0.0;
            f107bar = 0.0;
            ap = 0.0;
            avgap = 0.0;
            kp = 0.0;
            sumkp = 0.0;

            // check if any whole days in jdF
            jd1 = Math.Floor(jd + jdFrac) + 0.5;  // want jd at 0 hr
            mfme = (jd + jdFrac - jd1) * 1440.0;
            if (mfme < 0.0)
                mfme = 1440.0 + mfme;

            // ---- set flux time based on when measurments were taken
            // ---- before may 31, 1991, use 1700 hrs (1020 minutes)
            if (jd > 2448407.5)
                fluxtime = 1200.0;
            else
                fluxtime = 1020.0;

            // ---- read data for day of interest
            jdspwstarto = Math.Floor(jd + jdFrac - spwdata[0].mjd - 2400000.5);
            recnum = Convert.ToInt32(jdspwstarto);

            // --------------------  implementation   ----------------------
            // check for out of bound values
            if ((recnum >= 1) && (recnum < numbspw - 1))
            {
                // ---- set non-interpolated values
                if (fluxtype == 'a')  // adjusted
                {
                    f107 = spwdata[recnum].adjf10;
                    if (f81type == 'l')  // last
                        f107bar = spwdata[recnum].adjlstf81;
                    else
                        f107bar = spwdata[recnum].adjctrf81;
                }
                else // observed
                {
                    f107 = spwdata[recnum].obsf10;
                    if (f81type == 'l')  // last
                        f107bar = spwdata[recnum].obslstf81;
                    else
                        f107bar = spwdata[recnum].obsctrf81;
                }
                avgap = spwdata[recnum].avgap;
                sumkp = spwdata[recnum].sumkp;

                // ---- get last ap/kp array value from the current time value
                idx = (int)(Math.Floor(mfme / 180.0)); // values change at 0, 3, 6, ... hrs
                if (idx < 0) idx = 0;
                if (idx > 7) idx = 7;

                ap = spwdata[recnum].aparr[idx];
                kp = spwdata[recnum].kparr[idx] * 0.1;

                // some tricks to get the last 8 values
                j = idx;
                for (i = 1; i <= 8; i++)
                {
                    if (j >= 0)
                    {
                        aparr[8 - i] = spwdata[recnum].aparr[j];
                        kparr[8 - i] = spwdata[recnum].kparr[j];
                    }
                    else
                    {
                        aparr[8 - i] = spwdata[recnum - 1].aparr[8 + j];
                        kparr[8 - i] = spwdata[recnum - 1].kparr[8 + j];
                    }
                    j = j - 1;
                }

                // ---- do linear interpolation
                if (interp == 'l')
                {
                    if (mfme > fluxtime - 720.0) // go 12 hrs before...
                    {
                        if (mfme > fluxtime)
                            recnumt = recnum + 1;
                        else
                            recnumt = recnum - 1;
                        fixf = (fluxtime - mfme) / 1440.0;
                    }
                    else
                    {
                        recnumt = recnum - 1;
                        fixf = (mfme + (1440 - fluxtime)) / 1440.0;
                    }
                    if (fluxtype == 'a') // adjusted or observed values
                    {
                        tf107 = spwdata[recnumt].adjf10;
                        if (f81type == 'l')
                            tf107bar = spwdata[recnumt].adjlstf81;
                        else
                            tf107bar = spwdata[recnumt].adjctrf81;
                    }
                    else
                    {
                        tf107 = spwdata[recnumt].obsf10;
                        if (f81type == 'l')
                            tf107bar = spwdata[recnumt].obslstf81;
                        else
                            tf107bar = spwdata[recnumt].obsctrf81;
                    }

                    if (mfme <= fluxtime)
                    {
                        if (mfme > fluxtime - 720.0)
                        {
                            f107 = f107 + (tf107 - f107) * fixf;
                            f107bar = f107bar + (tf107bar - f107bar) * fixf;
                        }
                        else
                        {
                            f107 = tf107 + (f107 - tf107) * fixf;
                            f107bar = tf107bar + (f107bar - tf107bar) * fixf;
                        }
                    }
                    else
                    {
                        f107 = f107 + (tf107 - f107) * fixf;
                        f107bar = f107bar + (tf107bar - f107bar) * fixf;
                    }
                }

                // ---- do spline interpolations
                if (interp == 's')
                {
                    off1 = 1;
                    off2 = 2;
                    fixf = mfme / 1440.0;  // days for mjd
                    // setup so the interval is in between points 2 and 3
                    if (fluxtype == 'a') // adjusted 
                    {
                        f107 = MathTimeLibr.cubicinterp(
                               spwdata[recnum - off1].adjf10, spwdata[recnum].adjf10, spwdata[recnum + off1].adjf10, spwdata[recnum + off2].adjf10,
                               spwdata[recnum - off1].mjd, spwdata[recnum].mjd, spwdata[recnum + off1].mjd, spwdata[recnum + off2].mjd,
                               spwdata[recnum].mjd + fixf);
                        if (f81type == 'l')
                            f107bar = MathTimeLibr.cubicinterp(
                               spwdata[recnum - off1].adjlstf81, spwdata[recnum].adjlstf81, spwdata[recnum + off1].adjlstf81, spwdata[recnum + off2].adjlstf81,
                               spwdata[recnum - off1].mjd, spwdata[recnum].mjd, spwdata[recnum + off1].mjd, spwdata[recnum + off2].mjd,
                               spwdata[recnum].mjd + fixf);
                        else
                            f107bar = MathTimeLibr.cubicinterp(
                               spwdata[recnum - off1].adjctrf81, spwdata[recnum].adjctrf81, spwdata[recnum + off1].adjctrf81, spwdata[recnum + off2].adjctrf81,
                               spwdata[recnum - off1].mjd, spwdata[recnum].mjd, spwdata[recnum + off1].mjd, spwdata[recnum + off2].mjd,
                               spwdata[recnum].mjd + fixf);
                    }
                    else  // observed values
                    {
                        f107 = MathTimeLibr.cubicinterp(
                            spwdata[recnum - off1].obsf10, spwdata[recnum].obsf10, spwdata[recnum + off1].obsf10, spwdata[recnum + off2].obsf10,
                            spwdata[recnum - off1].mjd, spwdata[recnum].mjd, spwdata[recnum + off1].mjd, spwdata[recnum + off2].mjd,
                            spwdata[recnum].mjd + fixf);
                        if (f81type == 'l')
                            f107bar = MathTimeLibr.cubicinterp(
                                spwdata[recnum - off1].obslstf81, spwdata[recnum].obslstf81, spwdata[recnum + off1].obslstf81, spwdata[recnum + off2].obslstf81,
                                spwdata[recnum - off1].mjd, spwdata[recnum].mjd, spwdata[recnum + off1].mjd, spwdata[recnum + off2].mjd,
                                spwdata[recnum].mjd + fixf);
                        else
                            f107bar = MathTimeLibr.cubicinterp(
                                spwdata[recnum - off1].obsctrf81, spwdata[recnum].obsctrf81, spwdata[recnum + off1].obslstf81, spwdata[recnum + off2].obslstf81,
                                spwdata[recnum - off1].mjd, spwdata[recnum].mjd, spwdata[recnum + off1].mjd, spwdata[recnum + off2].mjd,
                                spwdata[recnum].mjd + fixf);
                    }

                }  // spline interpolation

                else
                     if (inputtype == 'u')  // ---- user input data
                {
                    // this is for data that may be simulated, or otherwise different from
                    // the current noaa data

                    // there could also be the interpolation stuff from above

                }
                else
                     if (inputtype == 'c')  // ---- constant data
                {
                    // this data is the same all the time
                    // leave the same as when it enters
                }
            }

        }  // findspwparam


        // ---------------------------------------------------------------------------
        //
        //                           function kp2ap
        //
        //  this function converts kp to ap using cubic splines. notice the arrays go
        //  beyond the range of values to permit endpoint evaluations without additional
        //  logic. the arrays have an extra element so they will start at 1. also, the normal
        //  cubic splining is usually between pts 2 and 3, but here it is between 3 and 4 
        //  because i've added 2 additional points at the start. 
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    kpin        - kp
        //
        //  outputs       :
        //    kp2ap       - ap
        //
        //  locals        :
        //    idx         - index of function value above the input value so the input
        //                  value is between the 2nd and 3rd point
        //
        //  coupling      :
        //    cubicspl    - perform the splining operation given 4 points
        //
        //  references    :
        //    vallado       2013, 558
        // ---------------------------------------------------------------------------

        public double kp2ap(double kpin)
        {
            double[] bap = new double[]
                 { 0, -0.00001, -0.001,
                   0, 2, 3, 4, 5, 6, 7, 9, 12, 15, 18, 22, 27, 32,
                   39, 48, 56, 67, 80, 94, 111, 132, 154, 179, 207, 236, 300, 400, 900
                 };

            double[] bkp = new double[]
           { 0, -6.6666667,  -3.3333  ,
             0, 3.3333, 6.6667, 10.0, 13.3333, 16.6667, 20.0, 23.3333,
             26.6667, 30.0, 33.3333, 36.6667, 40.0, 43.3333,
             46.6667, 50.0, 53.3333, 56.6667, 60.0, 63.3333, 66.6667, 70.0,
             73.3333, 76.6667, 80.0, 83.3333, 86.66667, 90.0, 93.3333
           };

            int idx;

            idx = 1;
            while ((idx < 33) && (kpin > bkp[idx]))
            {
                idx = idx + 1;
            }

            if (idx > 2)
            {
                return MathTimeLibr.cubicinterp(bap[idx - 2], bap[idx - 1], bap[idx], bap[idx + 1],
                                     bkp[idx - 2], bkp[idx - 1], bkp[idx], bkp[idx + 1],
                                     kpin);
            } // if idx > 3
            else
                return 0.0;
        }   // kp2ap


        // ---------------------------------------------------------------------------
        //
        //                           function ap2kp
        //
        //  this function converts ap to kp using cubic splines. notice the values go
        //  beyond the range of values to permit endpoint evaluations without additional
        //  logic. the arrays have an extra element so they will start at 1.
        //
        //  author        : david vallado             davallado@gmail.com      20 jan 2025
        //
        //  inputs          description                              range / units
        //    kpin        - kp
        //
        //  outputs       :
        //    ap2kp       - ap
        //
        //  locals        :
        //    idx         - index of function value above the input value so the input
        //                  value is between the 2nd and 3rd point
        //
        //  coupling      :
        //
        //  references    :
        //    vallado       2013, 558
        // ---------------------------------------------------------------------------

        public double ap2kp(double apin)
        {
            double[] bap = new double[]
           { 0, -0.00001, -0.001 ,
             0, 2, 3, 4, 5, 6, 7, 9, 12, 15, 18, 22, 27, 32,
             39, 48, 56, 67, 80, 94, 111, 132, 154, 179, 207, 236, 300, 400, 900
           };
            double[] bkp = new double[]
           { 0, -6.6666667,  -3.3333  ,
             0, 3.3333, 6.6667, 10.0, 13.3333, 16.6667, 20.0, 23.3333,
             26.6667, 30.0, 33.3333, 36.6667, 40.0, 43.3333,
             46.6667, 50.0, 53.3333, 56.6667, 60.0, 63.3333, 66.6667, 70.0,
             73.3333, 76.6667, 80.0, 83.3333, 86.66667, 90.0, 93.3333
           };

            int idx;

            idx = 1;
            while ((idx < 33) && (apin > bap[idx]))
            {
                idx = idx + 1;
            }

            if (idx > 2)
            {
                return MathTimeLibr.cubicinterp(bkp[idx - 2], bkp[idx - 1], bkp[idx], bkp[idx + 1],
                                     bap[idx - 2], bap[idx - 1], bap[idx], bap[idx + 1],
                                     apin);
            } // if idxs > 3
            else
                return 0.0;
        }  // ap2kp

    }  // EOPSPWLib



    // -----------------------------------------------------------------------------
    // alternate way to import the data
    //public class iau80Data
    //{
    //    public double[,] rar80;
    //    public int[,] iar80;

    //    /// <summary>
    //    /// constructor for the iau80 record
    //    /// </summary>
    //    public iau80Data
    //        ()
    //    {
    //        double convrt;

    //        // ------------------------  implementation   -------------------
    //        convrt = 0.0001 * Math.PI / (180.0 * 3600.0); // 0.0001" to rad

    //        // ---- Load terms (sorted) for IAU80 nutation
    //        iar80 = new int[,]
    //                { { 0,  0,  0,  0,  1 },
    //                  { 0,  0,  2, -2,  2 },
    //                  { 0,  0,  2,  0,  2 },
    //                  { 0,  0,  0,  0,  2 },
    //                  { 0,  1,  0,  0,  0 },
    //                  { 1,  0,  0,  0,  0 },
    //                  { 0,  1,  2, -2,  2 },
    //                  { 0,  0,  2,  0,  1 },
    //                  { 1,  0,  2,  0,  2 },
    //                  { 0, -1,  2, -2,  2 },
    //                  { 1,  0,  0, -2,  0 },
    //                  { 0,  0,  2, -2,  1 },
    //                  {-1,  0,  2,  0,  2 },
    //                  { 1,  0,  0,  0,  1 },
    //                  { 0,  0,  0,  2,  0 },
    //                  {-1,  0,  2,  2,  2 },
    //                  {-1,  0,  0,  0,  1 },
    //                  { 1,  0,  2,  0,  1 },
    //                  { 2,  0,  0, -2,  0 },
    //                  {-2,  0,  2,  0,  1 },
    //                  { 0,  0,  2,  2,  2 },
    //                  { 2,  0,  2,  0,  2 },
    //                  { 2,  0,  0,  0,  0 },
    //                  { 1,  0,  2, -2,  2 },
    //                  { 0,  0,  2,  0,  0 },
    //                  { 0,  0,  2, -2,  0 },
    //                  {-1,  0,  2,  0,  1 },
    //                  { 0,  2,  0,  0,  0 },
    //                  { 0,  2,  2, -2,  2 },
    //                  {-1,  0,  0,  2,  1 },
    //                  { 0,  1,  0,  0,  1 },
    //                  { 1,  0,  0, -2,  1 },
    //                  { 0, -1,  0,  0,  1 },
    //                  { 2,  0, -2,  0,  0 },
    //                  {-1,  0,  2,  2,  1 },
    //                  { 1,  0,  2,  2,  2 },
    //                  { 0, -1,  2,  0,  2 },
    //                  { 0,  0,  2,  2,  1 },
    //                  { 1,  1,  0, -2,  0 },
    //                  { 0,  1,  2,  0,  2 },
    //                  {-2,  0,  0,  2,  1 },
    //                  { 0,  0,  0,  2,  1 },
    //                  { 2,  0,  2, -2,  2 },
    //                  { 1,  0,  0,  2,  0 },
    //                  { 1,  0,  2, -2,  1 },
    //                  { 0,  0,  0, -2,  1 },
    //                  { 0, -1,  2, -2,  1 },
    //                  { 2,  0,  2,  0,  1 },
    //                  { 1, -1,  0,  0,  0 },
    //                  { 1,  0,  0, -1,  0 },
    //                  { 0,  0,  0,  1,  0 },
    //                  { 0,  1,  0, -2,  0 },
    //                  { 1,  0, -2,  0,  0 },
    //                  { 2,  0,  0, -2,  1 },
    //                  { 0,  1,  2, -2,  1 },
    //                  { 1,  1,  0,  0,  0 },
    //                  { 1, -1,  0, -1,  0 },
    //                  {-1, -1,  2,  2,  2 },
    //                  { 0, -1,  2,  2,  2 },
    //                  { 1, -1,  2,  0,  2 },
    //                  { 3,  0,  2,  0,  2 },
    //                  {-2,  0,  2,  0,  2 },
    //                  { 1,  0,  2,  0,  0 },
    //                  {-1,  0,  2,  4,  2 },
    //                  { 1,  0,  0,  0,  2 },
    //                  {-1,  0,  2, -2,  1 },
    //                  { 0, -2,  2, -2,  1 },
    //                  {-2,  0,  0,  0,  1 },
    //                  { 2,  0,  0,  0,  1 },
    //                  { 3,  0,  0,  0,  0 },
    //                  { 1,  1,  2,  0,  2 },
    //                  { 0,  0,  2,  1,  2 },
    //                  { 1,  0,  0,  2,  1 },
    //                  { 1,  0,  2,  2,  1 },
    //                  { 1,  1,  0, -2,  1 },
    //                  { 0,  1,  0,  2,  0 },
    //                  { 0,  1,  2, -2,  0 },
    //                  { 0,  1, -2,  2,  0 },
    //                  { 1,  0, -2,  2,  0 },
    //                  { 1,  0, -2, -2,  0 },
    //                  { 1,  0,  2, -2,  0 },
    //                  { 1,  0,  0, -4,  0 },
    //                  { 2,  0,  0, -4,  0 },
    //                  { 0,  0,  2,  4,  2 },
    //                  { 0,  0,  2, -1,  2 },
    //                  {-2,  0,  2,  4,  2 },
    //                  { 2,  0,  2,  2,  2 },
    //                  { 0, -1,  2,  0,  1 },
    //                  { 0,  0, -2,  0,  1 },
    //                  { 0,  0,  4, -2,  2 },
    //                  { 0,  1,  0,  0,  2 },
    //                  { 1,  1,  2, -2,  2 },
    //                  { 3,  0,  2, -2,  2 },
    //                  {-2,  0,  2,  2,  2 },
    //                  {-1,  0,  0,  0,  2 },
    //                  { 0,  0, -2,  2,  1 },
    //                  { 0,  1,  2,  0,  1 },
    //                  {-1,  0,  4,  0,  2 },
    //                  { 2,  1,  0, -2,  0 },
    //                  { 2,  0,  0,  2,  0 },
    //                  { 2,  0,  2, -2,  1 },
    //                  { 2,  0, -2,  0,  1 },
    //                  { 1, -1,  0, -2,  0 },
    //                  {-1,  0,  0,  1,  1 },
    //                  {-1, -1,  0,  2,  1 },
    //                  { 0,  1,  0,  1,  0 } };
    //        rar80 = new double[,]
    //                { { -171996.0, -174.2, 92025.0,  8.9 },
    //                  {  -13187.0,   -1.6,  5736.0, -3.1 },
    //                  {   -2274.0,   -0.2,   977.0, -0.5 },
    //                  {    2062.0,    0.2,  -895.0,  0.5 },
    //                  {    1426.0,   -3.4,    54.0, -0.1 },
    //                  {     712.0,    0.1,    -7.0,  0.0 },
    //                  {    -517.0,    1.2,   224.0, -0.6 },
    //                  {    -386.0,   -0.4,   200.0,  0.0 },
    //                  {    -301.0,    0.0,   129.0, -0.1 },
    //                  {     217.0,   -0.5,   -95.0,  0.3 },
    //                  {    -158.0,    0.0,    -1.0,  0.0 },
    //                  {     129.0,    0.1,   -70.0,  0.0 },
    //                  {     123.0,    0.0,   -53.0,  0.0 },
    //                  {      63.0,    0.1,   -33.0,  0.0 },
    //                  {      63.0,    0.0,    -2.0,  0.0 },
    //                  {     -59.0,    0.0,    26.0,  0.0 },
    //                  {     -58.0,   -0.1,    32.0,  0.0 },
    //                  {     -51.0,    0.0,    27.0,  0.0 },
    //                  {      48.0,    0.0,     1.0,  0.0 },
    //                  {      46.0,    0.0,   -24.0,  0.0 },
    //                  {     -38.0,    0.0,    16.0,  0.0 },
    //                  {     -31.0,    0.0,    13.0,  0.0 },
    //                  {      29.0,    0.0,    -1.0,  0.0 },
    //                  {      29.0,    0.0,   -12.0,  0.0 },
    //                  {      26.0,    0.0,    -1.0,  0.0 },
    //                  {     -22.0,    0.0,     0.0,  0.0 },
    //                  {      21.0,    0.0,   -10.0,  0.0 },
    //                  {      17.0,   -0.1,     0.0,  0.0 },
    //                  {     -16.0,    0.1,     7.0,  0.0 },
    //                  {      16.0,    0.0,    -8.0,  0.0 },
    //                  {     -15.0,    0.0,     9.0,  0.0 },
    //                  {     -13.0,    0.0,     7.0,  0.0 },
    //                  {     -12.0,    0.0,     6.0,  0.0 },
    //                  {      11.0,    0.0,     0.0,  0.0 },
    //                  {     -10.0,    0.0,     5.0,  0.0 },
    //                  {      -8.0,    0.0,     3.0,  0.0 },
    //                  {      -7.0,    0.0,     3.0,  0.0 },
    //                  {      -7.0,    0.0,     3.0,  0.0 },
    //                  {      -7.0,    0.0,     0.0,  0.0 },
    //                  {       7.0,    0.0,    -3.0,  0.0 },
    //                  {      -6.0,    0.0,     3.0,  0.0 },
    //                  {      -6.0,    0.0,     3.0,  0.0 },
    //                  {       6.0,    0.0,    -3.0,  0.0 },
    //                  {       6.0,    0.0,     0.0,  0.0 },
    //                  {       6.0,    0.0,    -3.0,  0.0 },
    //                  {      -5.0,    0.0,     3.0,  0.0 },
    //                  {      -5.0,    0.0,     3.0,  0.0 },
    //                  {      -5.0,    0.0,     3.0,  0.0 },
    //                  {       5.0,    0.0,     0.0,  0.0 },
    //                  {      -4.0,    0.0,     0.0,  0.0 },
    //                  {      -4.0,    0.0,     0.0,  0.0 },
    //                  {      -4.0,    0.0,     0.0,  0.0 },
    //                  {       4.0,    0.0,     0.0,  0.0 },
    //                  {       4.0,    0.0,    -2.0,  0.0 },
    //                  {       4.0,    0.0,    -2.0,  0.0 },
    //                  {      -3.0,    0.0,     0.0,  0.0 },
    //                  {      -3.0,    0.0,     0.0,  0.0 },
    //                  {      -3.0,    0.0,     1.0,  0.0 },
    //                  {      -3.0,    0.0,     1.0,  0.0 },
    //                  {      -3.0,    0.0,     1.0,  0.0 },
    //                  {      -3.0,    0.0,     1.0,  0.0 },
    //                  {      -3.0,    0.0,     1.0,  0.0 },
    //                  {       3.0,    0.0,     0.0,  0.0 },
    //                  {      -2.0,    0.0,     1.0,  0.0 },
    //                  {      -2.0,    0.0,     1.0,  0.0 },
    //                  {      -2.0,    0.0,     1.0,  0.0 },
    //                  {      -2.0,    0.0,     1.0,  0.0 },
    //                  {      -2.0,    0.0,     1.0,  0.0 },
    //                  {       2.0,    0.0,    -1.0,  0.0 },
    //                  {       2.0,    0.0,     0.0,  0.0 },
    //                  {       2.0,    0.0,    -1.0,  0.0 },
    //                  {       2.0,    0.0,    -1.0,  0.0 },
    //                  {      -1.0,    0.0,     0.0,  0.0 },
    //                  {      -1.0,    0.0,     1.0,  0.0 },
    //                  {      -1.0,    0.0,     0.0,  0.0 },
    //                  {      -1.0,    0.0,     0.0,  0.0 },
    //                  {      -1.0,    0.0,     0.0,  0.0 },
    //                  {      -1.0,    0.0,     0.0,  0.0 },
    //                  {      -1.0,    0.0,     0.0,  0.0 },
    //                  {      -1.0,    0.0,     0.0,  0.0 },
    //                  {      -1.0,    0.0,     0.0,  0.0 },
    //                  {      -1.0,    0.0,     0.0,  0.0 },
    //                  {      -1.0,    0.0,     0.0,  0.0 },
    //                  {      -1.0,    0.0,     0.0,  0.0 },
    //                  {      -1.0,    0.0,     0.0,  0.0 },
    //                  {      -1.0,    0.0,     1.0,  0.0 },
    //                  {      -1.0,    0.0,     0.0,  0.0 },
    //                  {      -1.0,    0.0,     0.0,  0.0 },
    //                  {      -1.0,    0.0,     0.0,  0.0 },
    //                  {       1.0,    0.0,     0.0,  0.0 },
    //                  {       1.0,    0.0,     0.0,  0.0 },
    //                  {       1.0,    0.0,    -1.0,  0.0 },
    //                  {       1.0,    0.0,     0.0,  0.0 },
    //                  {       1.0,    0.0,    -1.0,  0.0 },
    //                  {       1.0,    0.0,    -1.0,  0.0 },
    //                  {       1.0,    0.0,     0.0,  0.0 },
    //                  {       1.0,    0.0,     0.0,  0.0 },
    //                  {       1.0,    0.0,     0.0,  0.0 },
    //                  {       1.0,    0.0,     0.0,  0.0 },
    //                  {       1.0,    0.0,     0.0,  0.0 },
    //                  {       1.0,    0.0,    -1.0,  0.0 },
    //                  {       1.0,    0.0,     0.0,  0.0 },
    //                  {       1.0,    0.0,     0.0,  0.0 },
    //                  {       1.0,    0.0,     0.0,  0.0 },
    //                  {       1.0,    0.0,     0.0,  0.0 },
    //                  {       1.0,    0.0,     0.0,  0.0 } };
    //        for (int i = 0; i <= 105; i++)
    //        {
    //            for (int k = 0; k <= 3; k++)
    //            {
    //                rar80[i, k] = rar80[i, k] * convrt;
    //            }
    //        }
    //    }
    //}  // class iau80rec

}  //  EOPSPWMethods


