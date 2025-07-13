// <copyright file="JulianDay.cs" company="Fundamentals Of Astrodynamics Contributors">
// Copyright (c) 2025 Fundamentals Of Astrodynamics Contributors.
// This code is released under the [GNU Affero General Public License v3.0](./LICENSE).
// You are free to use, modify, and distribute it under the terms of this license,
// provided that any modifications or derivative works are also made available
// under the same license.
// </copyright>

using System.ComponentModel;

using FundamentalsOfAstrodynamics.MathTime;

using Xunit.Abstractions;

namespace MathTime.xUnit
{
    public class JulianDateTests
    {
        private readonly ITestOutputHelper _o;

        public JulianDateTests(ITestOutputHelper output) => _o = output;

        [Fact]
        [Trait("Original", "DateTime Input")]
        public void OriginalJD1()
        {
            JulianDate jdDate = new JulianDate();
            jdDate.jday(1858, 11, 16, 12, 0, 0, out double jd, out double jdFrac);

            _o.WriteLine("    jd: " + jd.ToString());
            _o.WriteLine("jdFrac: " + jdFrac.ToString());

            Assert.Equal(2400000, jd + jdFrac);
        }

        [Fact]
        [Trait("Fixed", "DateTime Input")]
        public void FixedJD1()
        {
            JulianDate jdDate = new JulianDate();
            jdDate.jdayFixed(1858, 11, 16, 12, 0, 0, out double jd, out double jdFrac);

            _o.WriteLine("    jd: " + jd.ToString());
            _o.WriteLine("jdFrac: " + jdFrac.ToString());

            Assert.Equal(2400000, jd + jdFrac);
        }

        [Fact]
        [Trait("Original", "DateTime Input")]
        public void OriginalJD2()
        {
            JulianDate jdDate = new JulianDate();
            jdDate.jday(1877, 8, 11, 7, 30, 0, out double jd, out double jdFrac);

            _o.WriteLine("    jd: " + jd.ToString());
            _o.WriteLine("jdFrac: " + jdFrac.ToString());

            Assert.Equal(2406842.8125, jd + jdFrac);
        }

        [Fact]
        [Trait("Fixed", "DateTime Input")]
        public void FixedJD2()
        {
            JulianDate jdDate = new JulianDate();
            jdDate.jdayFixed(1877, 8, 11, 7, 30, 0, out double jd, out double jdFrac);

            _o.WriteLine("    jd: " + jd.ToString());
            _o.WriteLine("jdFrac: " + jdFrac.ToString());

            Assert.Equal(2406842.8125, jd + jdFrac);
        }

        [Fact]
        [Trait("Original", "DateTime Input")]
        public void OriginalJD3()
        {
            JulianDate jdDate = new JulianDate();
            jdDate.jday(1801, 1, 1, 12, 0, 0, out double jd, out double jdFrac);

            _o.WriteLine("    jd: " + jd.ToString());
            _o.WriteLine("jdFrac: " + jdFrac.ToString());

            Assert.Equal(2378862, jd + jdFrac);
        }

        [Fact]
        [Trait("Fixed", "DateTime Input")]
        public void FixedJD3()
        {
            JulianDate jdDate = new JulianDate();
            jdDate.jdayFixed(1801, 1, 1, 12, 0, 0, out double jd, out double jdFrac);

            _o.WriteLine("    jd: " + jd.ToString());
            _o.WriteLine("jdFrac: " + jdFrac.ToString());

            Assert.Equal(2378862, jd + jdFrac);
        }

        [Fact]
        [Trait("Original", "DateTime Input")]
        public void OriginalJD4()
        {
            JulianDate jdDate = new JulianDate();
            jdDate.jday(2099, 12, 31, 12, 0, 0, out double jd, out double jdFrac);

            _o.WriteLine("    jd: " + jd.ToString());
            _o.WriteLine("jdFrac: " + jdFrac.ToString());

            Assert.Equal(2488069, jd + jdFrac);
        }

        [Fact]
        [Trait("Fixed", "DateTime Input")]
        public void FixedJD4()
        {
            JulianDate jdDate = new JulianDate();
            jdDate.jdayFixed(2099, 12, 31, 12, 0, 0, out double jd, out double jdFrac);

            _o.WriteLine("    jd: " + jd.ToString());
            _o.WriteLine("jdFrac: " + jdFrac.ToString());

            Assert.Equal(2488069, jd + jdFrac);
        }

        [Fact]
        [Trait("Original", "DateTime Input")]
        public void OriginalJD5()
        {
            JulianDate jdDate = new JulianDate();
            jdDate.jday(1900, 2, 28, 12, 0, 0, out double jd, out double jdFrac);

            _o.WriteLine("    jd: " + jd.ToString());
            _o.WriteLine("jdFrac: " + jdFrac.ToString());

            Assert.Equal(2415079, jd + jdFrac);
        }

        [Fact]
        [Trait("Fixed", "DateTime Input")]
        public void FixedJD5()
        {
            JulianDate jdDate = new JulianDate();
            jdDate.jdayFixed(1900, 2, 28, 12, 0, 0, out double jd, out double jdFrac);

            _o.WriteLine("    jd: " + jd.ToString());
            _o.WriteLine("jdFrac: " + jdFrac.ToString());

            Assert.Equal(2415079, jd + jdFrac);
        }

        [Fact]
        [Trait("Original", "DateTime Input")]
        public void OriginalJD6()
        {
            JulianDate jdDate = new JulianDate();
            jdDate.jday(1900, 3, 1, 12, 0, 0, out double jd, out double jdFrac);

            _o.WriteLine("    jd: " + jd.ToString());
            _o.WriteLine("jdFrac: " + jdFrac.ToString());

            Assert.Equal(2415080, jd + jdFrac);
        }

        [Fact]
        [Trait("Fixed", "DateTime Input")]
        public void FixedJD6()
        {
            JulianDate jdDate = new JulianDate();
            jdDate.jdayFixed(1900, 3, 1, 12, 0, 0, out double jd, out double jdFrac);

            _o.WriteLine("    jd: " + jd.ToString());
            _o.WriteLine("jdFrac: " + jdFrac.ToString());

            Assert.Equal(2415080, jd + jdFrac);
        }

        [Fact]
        [Trait("JD", "Decimal Input")]
        public void JdDec()
        {
            JulianDate j = new JulianDate(2456293.520833M);

            _o.WriteLine("JdDec(): " + j.JD);

            Assert.Equal(2456293.520833M, j.JD);
        }

        [Fact]
        [Trait("JD", "DateTime Input")]
        public void JdDt()
        {
            JulianDate j = new JulianDate(new DateTime(2013, 1, 1, 0, 30, 0));

            _o.WriteLine("JdDt(): " + j.JD);

            Assert.Equal(2456293.520833M, j.JD, 6);
        }

        [Fact]
        [Trait("JDN", "Decimal Input")]
        public void JdnDec()
        {
            JulianDate j = new JulianDate(2456293.520833M);

            _o.WriteLine("JdnDec(): " + j.JDN);

            Assert.Equal(2456293, j.JDN);
        }

        [Fact]
        [Trait("JDN", "DateTime Input")]
        public void JdnDt()
        {
            JulianDate j = new JulianDate(new DateTime(2000, 1, 1, 12, 0, 0));

            _o.WriteLine("JdnDt(): " + j.JDN);

            Assert.Equal(2451545, j.JDN);
        }

        [Fact]
        [Trait("Reduced", "Decimal Input")]
        public void ReducedDec()
        {
            JulianDate j = new JulianDate(2400000.0M);

            _o.WriteLine("ReducedDec(): " + j.ToString());

            Assert.Equal(0, j.Reduced);

            j = new JulianDate();
            j.Reduced = 0;

            Assert.Equal(2400000.0M, j.JD);
        }

        [Fact]
        [Trait("Reduced", "DateTime Input")]
        public void ReducedDt()
        {
            JulianDate j = new JulianDate(new DateTime(1858, 11, 16, 12, 0, 0));

            _o.WriteLine("ReducedDt(): " + j.ToString());

            Assert.Equal(0, j.Reduced);
        }

        [Fact]
        [Trait("Modified", "Decimal Input")]
        public void ModifiedDec()
        {
            JulianDate j = new JulianDate(2400000.5M);

            Assert.Equal(0, j.Modified);
        }

        [Fact]
        [Trait("Modified", "DateTime Input")]
        public void ModifiedDt()
        {
            JulianDate j = new JulianDate(new DateTime(1858, 11, 17, 0, 0, 0));

            Assert.Equal(0, j.Modified);
        }

        [Fact]
        [Trait("Truncated", "Decimal Input")]
        public void TruncatedDec()
        {
            JulianDate j = new JulianDate(2440000.5M);

            Assert.Equal(0, j.Truncated);
        }

        [Fact]
        [Trait("Truncated", "DateTime Input")]
        public void TruncatedDt()
        {
            JulianDate j = new JulianDate(new DateTime(1968, 5, 24, 0, 0, 0));

            Assert.Equal(0, j.Truncated);
        }

        [Fact]
        [Trait("Dublin", "Decimal Input")]
        public void DublinDec()
        {
            JulianDate j = new JulianDate(2415020.0M);

            Assert.Equal(0, j.Dublin);
        }

        [Fact]
        [Trait("Dublin", "DateTime Input")]
        public void DublinDt()
        {
            JulianDate j = new JulianDate(new DateTime(1899, 12, 31, 12, 0, 0));

            Assert.Equal(0, j.Dublin);
        }

        [Fact]
        [Trait("CNES", "Decimal Input")]
        public void CNESDec()
        {
            JulianDate j = new JulianDate(2433282.5M);

            Assert.Equal(0, j.CNES);
        }

        [Fact]
        [Trait("CNES", "DateTime Input")]
        public void CNESDt()
        {
            JulianDate j = new JulianDate(new DateTime(1950, 1, 1, 0, 0, 0));

            Assert.Equal(0, j.CNES);
        }

        [Fact]
        [Trait("CCSDS", "Decimal Input")]
        public void CCSDSDec()
        {
            JulianDate j = new JulianDate(2436204.5M);

            Assert.Equal(0, j.CCSDS);
        }

        [Fact]
        [Trait("CCSDS", "DateTime Input")]
        public void CCSDSDt()
        {
            JulianDate j = new JulianDate(new DateTime(1958, 1, 1, 0, 0, 0));

            Assert.Equal(0, j.CCSDS);
        }

        [Fact]
        [Trait("ModifiedJD2000", "Decimal Input")]
        public void ModifiedJD2000Dec()
        {
            JulianDate j = new JulianDate(2451544.5M);

            Assert.Equal(0, j.ModifiedJD2000);
        }

        [Fact]
        [Trait("ModifiedJD2000", "DateTime Input")]
        public void ModifiedJD2000Dt()
        {
            JulianDate j = new JulianDate(new DateTime(2000, 1, 1, 0, 0, 0));

            Assert.Equal(0, j.ModifiedJD2000);
        }

        [Fact]
        [Trait("Lilian", "Decimal Input")]
        public void LilianDec()
        {
            JulianDate j = new JulianDate(2299159.5M);

            Assert.Equal(0, j.Lilian);
        }

        [Fact]
        [Trait("Lilian", "DateTime Input")]
        public void LilianDt()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => _ = new FundamentalsOfAstrodynamics.MathTime.JulianDate(new DateTime(1582, 10, 15, 0, 0, 0)));

            JulianDate j = new JulianDate(new DateTime(2025, 7, 11, 0, 0, 0));

            Assert.Equal(161708, j.Lilian);
        }

        [Fact]
        [Trait("RataDie", "Decimal Input")]
        public void RataDieDec()
        {
            JulianDate j = new JulianDate(1721424.5M);

            Assert.Equal(0, j.RataDie);
        }

        [Fact]
        [Trait("RataDie", "DateTime Input")]
        public void RataDieDt()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => _ = new FundamentalsOfAstrodynamics.MathTime.JulianDate(new DateTime(1, 1, 1, 0, 0, 0)));

            JulianDate j = new JulianDate(new DateTime(2025, 7, 11, 0, 0, 0));

            Assert.Equal(739443, j.RataDie);
        }

        [Fact]
        [Trait("MarsSol", "Decimal Input")]
        public void MarsSolDec()
        {
            JulianDate j = new JulianDate(2405522.0M);

            Assert.Equal(0, j.MarsSol);
        }

        [Fact]
        [Trait("MarsSol", "DateTime Input")]
        public void MarsSolDt()
        {
            JulianDate j = new JulianDate(new DateTime(1873, 12, 29, 12, 0, 0));

            Assert.Equal(0, j.MarsSol);

            // The Ares III crew step foot on Mars.
            j = new JulianDate(new DateTime(2035, 11, 7, 12, 0, 0));

            _o.WriteLine("    JD: " + j.JD);

            Assert.Equal(57535.353142123037693797506545M, j.MarsSol);
        }

        [Fact]
        [Trait("AresIII", "Decimal Input")]
        public void AresIIIDec()
        {
            // The Ares III crew step foot on Mars.
            JulianDate j = new JulianDate(2464639.0M);

            Assert.Equal(0, j.AresIII);
        }

        [Fact]
        [Trait("AresIII", "DateTime Input")]
        public void AresIIIDt()
        {
            // The Ares III crew step foot on Mars.
            JulianDate j = new JulianDate(new DateTime(2035, 11, 7, 12, 0, 0));

            Assert.Equal(0, j.AresIII);

            // The Ares III aborts their mission due to a sand storm (according to the book).
            j.AresIII = 6;

            Assert.Equal(2464645.16494M, j.JD);
        }

        [Fact]
        [Trait("UnixTime", "Decimal Input")]
        public void UnixTimeDec()
        {
            JulianDate j = new JulianDate(2440587.5M);

            Assert.Equal(0, j.UnixTime);
        }

        [Fact]
        [Trait("UnixTime", "DateTime Input")]
        public void UnixTimeDt()
        {
            JulianDate j = new JulianDate(new DateTime(1970, 1, 1, 0, 0, 0));

            Assert.Equal(0, j.UnixTime);
        }

        [Fact]
        [Trait("JavaScriptDate", "Decimal Input")]
        public void JavaScriptDateDec()
        {
            JulianDate j = new JulianDate(2440587.5M);

            Assert.Equal(0, j.JavaScriptDate);
        }

        [Fact]
        [Trait("JavaScriptDate", "DateTime Input")]
        public void JavaScriptDateDt()
        {
            JulianDate j = new JulianDate(new DateTime(1970, 1, 1, 0, 0, 0));

            Assert.Equal(0, j.JavaScriptDate);
        }

        [Fact]
        [Trait("Gregorian", "Decimal Input")]
        public void GregorianDec()
        {
            JulianDate j = new JulianDate(2460868.0M);

            Assert.Equal(new DateTime(2025, 7, 11, 12, 0, 0), j.Gregorian);
        }

        [Fact]
        [Trait("Gregorian", "DateTime Input")]
        public void GregorianDt()
        {
            // The Ares III crew step foot on Mars.
            JulianDate j = new JulianDate(new DateTime(2025, 7, 11, 12, 0, 0));

            Assert.Equal(new DateTime(2025, 7, 11, 12, 0, 0), j.Gregorian);
        }
    }
}
