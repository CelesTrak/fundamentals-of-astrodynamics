// <copyright file="JulianDate.cs" company="Fundamentals Of Astrodynamics Contributors">
// Copyright (c) 2025 Fundamentals Of Astrodynamics Contributors.
// This code is released under the [GNU Affero General Public License v3.0](./LICENSE).
// You are free to use, modify, and distribute it under the terms of this license,
// provided that any modifications or derivative works are also made available
// under the same license.
// </copyright>

using System;
using System.Globalization;

namespace FundamentalsOfAstrodynamics.MathTime
{
    /// <summary>
    /// This class stores and converts a Julian Date (JD).
    /// The Julian Date is defined by each elapsed day since noon on January 1, 4713 BCE.
    /// </summary>
    public class JulianDate
    {
        /// <summary>
        /// Julian date (JD).
        /// </summary>
        /// <remarks>
        /// Per <see href="https://en.wikipedia.org/wiki/Julian_day">Wikipedia</see>:
        /// The Julian date (JD) of any instant is the Julian day number plus the fraction
        /// of a day since the preceding noon in Universal Time. Julian dates are expressed
        /// as a Julian day number with a decimal fraction added.
        /// 
        /// For example, the Julian Date for 00:30:00.0 UT January 1, 2013, is 2456293.520833.
        /// </remarks>
        public decimal JD { get; set; }

        /// <summary>
        /// Julian day number (JDN).
        /// </summary>
        /// <remarks>
        /// Per <see href="https://en.wikipedia.org/wiki/Julian_day">Wikipedia</see>:
        /// The Julian day number (JDN) has the same epoch as the Julian period, but counts
        /// the number of days since the epoch rather than the number of years since then.
        /// Specifically, Julian day number 0 is assigned to the day starting at noon
        /// Universal Time on Monday, January 1, 4713 BC, proleptic Julian calendar
        /// (November 24, 4714 BC, in the proleptic Gregorian calendar).
        /// 
        /// For example, the Julian day number for the day starting at 12:00 UT (noon)
        /// on January 1, 2000, was 2451545.
        /// </remarks>
        public long JDN
        {
            get
            {
                return (long)Math.Truncate(JD);
            }
            set
            {
                JD = value;
            }
        }

        /// <summary>
        /// Epoch of 12:00 November 16, 1858.
        /// </summary>
        public decimal Reduced
        {
            get
            {
                return JD - 2400000.0M;
            }
            set
            {
                JD = value + 2400000.0M;
            }
        }

        /// <summary>
        /// Epoch of 0:00 November 17, 1858.
        /// </summary>
        /// <remarks>
        /// Introduced by the Smithsonian Astrophysical Observatory (SAO) in 1957 to record
        /// the orbit of Sputnik via an IBM 704 (36-bit machine) and using only 18 bits until
        /// August 7, 2576. MJD is the epoch of VAX/VMS and its successor OpenVMS, using
        /// 63-bit date/time, which allows times to be stored up to July 31, 31086, 02:48:05.47.
        /// </remarks>
        public decimal Modified
        {
            get
            {
                return JD - 2400000.5M;
            }
            set
            {
                JD = value + 2400000.5M;
            }
        }

        /// <summary>
        /// Epoch of 0:00 May 24, 1968.
        /// </summary>
        /// <remarks>
        /// Introduced by the National Aeronautics and Space Administration (NASA) in 1979 as
        /// part of a parallel grouped binary time code (PB-5) "designed specifically, although
        /// not exclusively, for spacecraft applications". TJD was a 4-digit day count from
        /// MJD 40000, which was May 24, 1968, represented as a 14-bit binary number. Since this
        /// code was limited to four digits, TJD recycled to zero on MJD 50000, or October 10, 1995,
        /// "which gives a long ambiguity period of 27.4 years". (NASA codes PB-1–PB-4 used a
        /// 3-digit day-of-year count.) Only whole days are represented. Time of day is expressed
        /// by a count of seconds of a day, plus optional milliseconds, microseconds and nanoseconds
        /// in separate fields. Later PB-5J was introduced which increased the TJD field to 16 bits,
        /// allowing values up to 65535, which will occur in the year 2147. There are five digits
        /// recorded after TJD 9999.
        /// </remarks>
        public decimal Truncated
        {
            get
            {
                return Math.Floor(JD - 2440000.5M);
            }
            set
            {
                JD = Math.Floor(value + 2440000.5M);
            }
        }

        /// <summary>
        /// Epoch of 12:00 December 31, 1899
        /// </summary>
        /// <remarks>
        /// Introduced by the International Astronomical Union (IAU) in 1955.
        /// </remarks>
        public decimal Dublin
        {
            get
            {
                return JD - 2415020.0M;
            }
            set
            {
                JD = value + 2415020.0M;
            }
        }

        /// <summary>
        /// Epoch of 0:00 January 1, 1950.
        /// </summary>
        /// <remarks>
        /// Introduced by the Centre national d'études spatiales (CNES).
        /// </remarks>
        public decimal CNES
        {
            get
            {
                return JD - 2433282.5M;
            }
            set
            {
                JD = value + 2433282.5M;
            }
        }

        /// <summary>
        /// Epoch of 0:00 January 1, 1958.
        /// </summary>
        /// <remarks>
        /// Introduced by the Consultative Committee for Space Data Systems (CCSDS).
        /// </remarks>
        public decimal CCSDS
        {
            get
            {
                return JD - 2436204.5M;
            }
            set
            {
                JD = value + 2436204.5M;
            }
        }

        /// <summary>
        /// Epoch of 0:00 January 1, 2000.
        /// </summary>
        /// <remarks>
        /// Introduced by the European Space Agency (ESA).
        /// </remarks>
        public decimal ModifiedJD2000
        {
            get
            {
                return JD - 2451544.5M;
            }
            set
            {
                JD = value + 2451544.5M;
            }
        }

        /// <summary>
        /// Epoch of 0:00 October 15, 1582.
        /// </summary>
        /// <remarks>
        /// Count of days of the Gregorian calendar.
        /// It was invented by Bruce G. Ohms of IBM in 1986 and is named for Aloysius Lilius, who devised the Gregorian Calendar.
        /// </remarks>
        public decimal Lilian
        {
            get
            {
                return Math.Floor(JD - 2299159.5M);
            }
            set
            {
                JD = value + 2299159.5M;
            }
        }

        /// <summary>
        /// Epoch of 0:00 January 1, 1 (proleptic Gregorian calendar).
        /// </summary>
        /// <remarks>
        /// Count of days of the Common Era (CE).
        /// It was named (after the Latin ablative feminine singular for "from a fixed date") by Howard Jacobson.
        /// </remarks>
        public decimal RataDie
        {
            get
            {
                return Math.Floor(JD - 1721424.5M);
            }
            set
            {
                JD = value + 1721424.5M;
            }
        }

        /// <summary>
        /// Epoch of 12:00 December 29, 1873.
        /// </summary>
        /// <remarks>
        /// Count of Martian days.
        /// </remarks>
        public decimal MarsSol
        {
            get
            {
                return (JD - 2405522.0M) / 1.02749M;
            }
            set
            {
                JD = (value * 1.02749M) + 2405522.0M;
            }
        }

        /// <summary>
        /// Epoch of 12:00 November 7, 2035.
        /// </summary>
        /// <remarks>
        /// Count of Martian solar days, with the epoch adjusted
        /// to begin at the landing of Ares III on the surface of Mars.
        /// </remarks>
        public decimal AresIII
        {
            get
            {
                return (JD - 2464639.0M) / 1.02749M;
            }
            set
            {
                JD = (value * 1.02749M) + 2464639.0M;
            }
        }

        /// <summary>
        /// Epoch of 0:00 January 1, 1970.
        /// </summary>
        /// <remarks>
        /// Count of seconds, excluding leap seconds.
        /// Unix time is typically encoded as a signed integer.
        /// </remarks>
        public int UnixTime
        {
            get
            {
                // Perform all math first, then cast
                // from decimal to int.
                return (int)((JD - 2440587.5M) * 86400.0M);
            }
            set
            {
                // Cast from int to decimal first.
                decimal temp = value;
                // Then do the math.
                JD = (temp / 86400.0M) + 2440587.5M;
            }
        }

        /// <summary>
        /// Epoch of 0:00 January 1, 1970.
        /// </summary>
        /// <remarks>
        /// Count of milliseconds, excluding leap seconds.
        /// JavaScript Date is typically encoded as a double.
        /// </remarks>
        public double JavaScriptDate
        {
            get
            {
                // Cast _jd to double first, since double has
                // a larger number space than decimal.
                // Then perform all math using doubles.
                return (((double)JD - 2440587.5D) * 86400000.0D);
            }
            set
            {
                // Cast from double to decimal first.
                decimal temp = (decimal)value;
                // Then do the math.
                JD = (temp / 86400000.0M) + 2440587.5M;
            }
        }

        public DateTime Gregorian
        {
            get
            {
                // ---------------- find year and days of the year -----------------
                decimal T1900 = (JD - 2415019.5M) / 365.25M;
                int Year = 1900 + (int)Math.Truncate(T1900);
                int LeapYears = (int)Math.Truncate((Year - 1901) * 0.25M);
                decimal Days = (JD - 2415019.5M) - ((Year - 1900) * 365.0M + LeapYears);

                // ------------- check for case of beginning of a year -------------
                if (Days < 1.0M)
                {
                    Year--;
                    LeapYears = (int)Math.Truncate((Year - 1901) * 0.25);
                    Days = (JD - 2415019.5M) - ((Year - 1900) * 365.0M + LeapYears);
                }

                // ------------------ find remaining data  // -----------------------
                // now add the daily time in to preserve accuracy
                return Days2DateTime(Year,
                                      Days);
            }
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="JulianDate"/> class.
        /// </summary>
        public JulianDate()
        {
            JD = decimal.Zero;
        }

        public JulianDate(decimal jd)
        {
            JD = jd;
        }

        public JulianDate(double jd)
        {
            JD = (decimal)jd;
        }

        public JulianDate(DateTime jd)
        {
            if ((jd.Year >= 1801) && (jd.Year <= 2099))
            {
                JD = (367.0M * jd.Year) -
                    Math.Truncate(7.0M * (jd.Year + Math.Truncate((jd.Month + 9.0M) / 12.0M)) * 0.25M) +
                    Math.Truncate(275.0M * jd.Month / 9.0M) +
                    jd.Day +
                    1721013.5M -
                    0.5M * Math.Sign((100.0M * jd.Year) + jd.Month - 190002.5M) +
                    0.5M;
                JD += (jd.Second + jd.Minute * 60.0M + jd.Hour * 3600.0M) / 86400.0M;
            }
            else
            {
                // The formula above works for years between 1801 and 2099
                // because it doesn't account for leap days.
                throw new ArgumentOutOfRangeException(
                    nameof(jd),
                    "Year must be between 1801 and 2099."
                    );
            }
        }

        public void jday
                (
                  int year, int mon, int day, int hr, int minute, double sec,
                  out double jd, out double jdFrac
                )
        {
            jd = 367.0 * year -
                 Math.Floor(7 * (year + Math.Floor((mon + 9) / 12.0)) * 0.25) +
                 Math.Floor(275 * mon / 9.0) +
                 day + 1721013.5;  // use - 678987.0 to go to mjd directly
            jdFrac = (sec + minute * 60.0 + hr * 3600.0) / 86400.0;

            // check that the day and fractional day are correct
            if (Math.Abs(jdFrac) > 1.0)
            {
                double dtt = Math.Floor(jdFrac);
                jd = jd + dtt;
                jdFrac = jdFrac - dtt;
            }

            // - 0.5*sgn(100.0*year + mon - 190002.5) + 0.5;
        }  //  jday

        /// <summary>
        /// Finds the Julian Date given the year, month, day, and time.
        /// The julian date is defined by each elapsed day since noon, jan 1, 4713 bc.
        /// </summary>
        /// <param name="year"></param>
        /// <param name="mon"></param>
        /// <param name="day"></param>
        /// <param name="hr"></param>
        /// <param name="minute"></param>
        /// <param name="sec"></param>
        /// <param name="jd"></param>
        /// <param name="jdFrac"></param>
        /// <remarks>
        /// The formula used mirrors the one used at <see href="https://aa.usno.navy.mil/faq/JD_formula"></see>
        /// </remarks>
        public void jdayFixed
                (int year,
                 int mon,
                 int day,
                 int hr,
                 int minute,
                 double sec,
                 out double jd,
                 out double jdFrac)
        {
            jd = 367.0 * year -
                 Math.Truncate(7 * (year + Math.Truncate((mon + 9) / 12.0)) * 0.25) +
                 Math.Truncate(275 * mon / 9.0) +
                 day + 1721013.5 // use - 678987.0 to go to mjd directly
                 - 0.5 * Math.Sign(100 * year + mon - 190002.5) + 0.5;
            jdFrac = (sec + minute * 60.0 + hr * 3600.0) / 86400.0;

            // check that the day and fractional day are correct
            if (Math.Abs(jdFrac) > 1.0)
            {
                double dtt = Math.Truncate(jdFrac);
                jd += dtt;
                jdFrac -= dtt;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="year"></param>
        /// <param name="days"></param>
        /// <returns></returns>
        public DateTime Days2DateTime(int year,
                                      decimal days)
        {
            int i = 1, DayOfYear, SumDays = 0;
            int Month, Day, Hour, Minute, Second;
            decimal temp;

            // Leading 0 so that C/C++/C# array index matches month number.
            int[] lMonth = { 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

            // Adjust array for leap year.
            if ((year % 4) == 0)
            {
                lMonth[2] = 29;
            }

            DayOfYear = (int)Math.Truncate(days);

            while ((DayOfYear > SumDays + lMonth[i]) && (i < 12))
            {
                SumDays += lMonth[i++];
            }
            Month = i;
            Day = DayOfYear - SumDays;

            // ------------------ find hours minutes and seconds ---------------
            temp = (days - DayOfYear) * 24.0M;
            Hour = (int)(Math.Truncate(temp));
            temp = (temp - Hour) * 60.0M;
            Minute = (int)(Math.Truncate(temp));
            temp = (temp - Minute) * 60.0M;
            Second = (int)(Math.Truncate(temp)); ;

            Console.WriteLine($"{year}, {Month}, {Day}, {Hour}, {Minute}, {Second}");

            return new DateTime(year, Month, Day, Hour, Minute, Second);
        }  //  days2mdhms

        /// <summary>
        /// Technically "cheating", but it can be done this way too.
        /// </summary>
        /// <param name="year"></param>
        /// <param name="days"></param>
        /// <returns></returns>
        public static DateTime Days2DateTime2(int year,
                                             decimal days)
        {
            // Default constructor is 01/01/0001
            DateTime result = new DateTime();
            result = result.AddYears(year - 1);
            result = result.AddDays((double)(days - 1.0M));
            return result;
        }

        /// <inheritdoc/>
        public override string ToString() => JD.ToString(CultureInfo.CurrentCulture) + " JD";
    }
}
