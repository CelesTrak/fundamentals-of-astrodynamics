using System;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.Text.RegularExpressions;
using System.IO;
using System.Collections.Generic;


using MathTimeMethods;     // Edirection, globals
using EOPSPWMethods;       // EOPDataClass, SPWDataClass, iau80Class, iau06Class
using AstroLibMethods;     // EOpt, gravityConst, astroConst, xysdataClass, jpldedataClass
using SGP4Methods;

using AstroLambertkMethods;
using System.Data;
using System.Diagnostics.Eventing.Reader;
using static AstroLibMethods.AstroLib;
using System.Threading;

namespace TestAllTool
{
    public partial class Form1 : Form
    {
        public MathTimeLib MathTimeLibr = new MathTimeLib();

        public EOPSPWLib EOPSPWLibr = new EOPSPWLib();

        public AstroLib AstroLibr = new AstroLib();

        public AstroLambertkLib AstroLambertkLibr = new AstroLambertkLib();

        public SGP4Lib Sgp4Libr = new SGP4Lib();

        StringBuilder strbuild = new StringBuilder();

        // public enum MathTimeLib.Edirection { efrom, eto };

        //   public iau80Data iau80rec;
        // 000. gives leading 0's
        // ;+00.;-00. gives signs in front
        string fmt = "0.00000000000";
        string fmtE = "0.0000000000E0";
        string fmt1 = "0.000000000000";

        public Form1()
        {
            InitializeComponent();
        }


        public void testvecouter()
        {
            double[] vec1 = new double[3];
            double[] vec2 = new double[3];
            double[,] mat1 = new double[3, 3];
            vec1 = new double[3] { 2.3, 4.7, -1.6 };
            vec2 = new double[3] { 0.3, -0.7, 6.0 };

            mat1 = MathTimeLibr.vecouter(vec1, vec2, 3);

            strbuild.AppendLine("vecout = " + mat1[0, 0].ToString(fmt).PadLeft(4) + " " + mat1[0, 1].ToString(fmt).PadLeft(4) + " " + mat1[0, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("vecout = " + mat1[1, 0].ToString(fmt).PadLeft(4) + " " + mat1[1, 1].ToString(fmt).PadLeft(4) + " " + mat1[1, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("vecout = " + mat1[2, 0].ToString(fmt).PadLeft(4) + " " + mat1[2, 1].ToString(fmt).PadLeft(4) + " " + mat1[2, 2].ToString(fmt).PadLeft(4));
        }
        public void testmatadd()
        {
            double[,] mat1 = new double[3, 3];
            double[,] mat2 = new double[3, 3];
            double[,] mat3 = new double[3, 3];
            Int32 mat1r, mat1c;
            mat1 = new double[3, 3] {{ 1.0,   2.0,   3.0 },
                                     { -1.1,   0.5,   2.0 },
                                     { -2.00,  4.00,  7.0 }};
            mat2 = new double[3, 3]{ { 1.0,  1.4, 1.8 },
                                     { 0.0,  2.6, -0.6 },
                                     { 1.9,  0.1, 7.1 } };

            mat1r = 3;
            mat1c = 3;

            mat3 = MathTimeLibr.matadd(mat1, mat2, mat1r, mat1c);

            strbuild.AppendLine("matadd = " + mat3[0, 0].ToString(fmt).PadLeft(4) + " " + mat3[0, 1].ToString(fmt).PadLeft(4) + " " + mat3[0, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("matadd = " + mat3[1, 0].ToString(fmt).PadLeft(4) + " " + mat3[1, 1].ToString(fmt).PadLeft(4) + " " + mat3[1, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("matadd = " + mat3[2, 0].ToString(fmt).PadLeft(4) + " " + mat3[2, 1].ToString(fmt).PadLeft(4) + " " + mat3[2, 2].ToString(fmt).PadLeft(4));
        }
        public void testmatsub()
        {
            double[,] mat1 = new double[3, 3];
            double[,] mat2 = new double[3, 3];
            double[,] mat3 = new double[3, 3];
            Int32 mat1r, mat1c;
            mat1 = new double[3, 3] {{ 1.0,   2.0,   3.0 },
                                     { -1.1,   0.5,   2.0 },
                                     { -2.00,  4.00,  7.0 }};
            mat2 = new double[3, 3]{ { 1.0,  1.4, 1.8 },
                                     { 0.0,  2.6, -0.6 },
                                     { 1.9,  0.1, 7.1 } };

            mat1r = 3;
            mat1c = 3;

            mat3 = MathTimeLibr.matsub(mat1, mat2, mat1r, mat1c);

            strbuild.AppendLine("matsub = " + mat3[0, 0].ToString(fmt).PadLeft(4) + " " + mat3[0, 1].ToString(fmt).PadLeft(4) + " " + mat3[0, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("matsub = " + mat3[1, 0].ToString(fmt).PadLeft(4) + " " + mat3[1, 1].ToString(fmt).PadLeft(4) + " " + mat3[1, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("matsub = " + mat3[2, 0].ToString(fmt).PadLeft(4) + " " + mat3[2, 1].ToString(fmt).PadLeft(4) + " " + mat3[2, 2].ToString(fmt).PadLeft(4));
        }
        public void testmatmult()
        {
            double[,] mat1 = new double[3, 3];
            double[,] mat2 = new double[3, 3];
            double[,] mat3 = new double[3, 3];
            Int32 mat1r, mat1c, mat2c;
            mat1r = 3;
            mat1c = 3;
            mat2c = 2;

            mat1 = new double[3, 3] {{ 1.0,   2.0,   3.0 },
                                     { -1.1,   0.5,   2.0 },
                                     { -2.00,  4.00,  7.0 }};
            mat2 = new double[3, 2]{         { 1.0,  1.4 },
                                             { 0.0,  2.6 },
                                             { 1.9,  0.1 } };
            mat3 = MathTimeLibr.matmult(mat1, mat2, mat1r, mat1c, mat2c);

            strbuild.AppendLine("matmult = " + mat3[0, 0].ToString(fmt).PadLeft(4) + " " + mat3[0, 1].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("matmult = " + mat3[1, 0].ToString(fmt).PadLeft(4) + " " + mat3[1, 1].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("matmult = " + mat3[2, 0].ToString(fmt).PadLeft(4) + " " + mat3[2, 1].ToString(fmt).PadLeft(4));
        }

        public void testmattrans()
        {
            double[,] mat1 = new double[3, 3];
            double[,] mat3 = new double[3, 3];
            int matr;
            matr = 3;
            mat1 = new double[3, 3] {{ 1.0,   2.0,   3.0 },
                                     { -1.1,   0.5,   2.0 },
                                     { -2.00,  4.00,  7.0 }};
            mat3 = MathTimeLibr.mattrans(mat1, matr);

            strbuild.AppendLine("mattrans = " + mat3[0, 0].ToString(fmt).PadLeft(4) + " " + mat3[0, 1].ToString(fmt).PadLeft(4) + " " + mat3[0, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("mattrans = " + mat3[1, 0].ToString(fmt).PadLeft(4) + " " + mat3[1, 1].ToString(fmt).PadLeft(4) + " " + mat3[1, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("mattrans = " + mat3[2, 0].ToString(fmt).PadLeft(4) + " " + mat3[2, 1].ToString(fmt).PadLeft(4) + " " + mat3[2, 2].ToString(fmt).PadLeft(4));
        }

        public void testmattransx()
        {
            double[,] mat1 = new double[3, 3];
            double[,] mat3 = new double[3, 3];
            int matr, matc;
            matr = 3;
            matc = 3;
            mat1 = new double[3, 3] {{ 1.0,   2.0,   3.0 },
                                     { -1.1,   0.5,   2.0 },
                                     { -2.00,  4.00,  7.0 }};

            mat3 = MathTimeLibr.mattransx(mat1, matr, matc);
        }

        public void testmatinverse()
        {
            double[,] mat1 = new double[3, 3];
            double[,] matinv = new double[3, 3];

            // enter by COL!!!!!!!!!!!!!!
            mat1 = new double[,] { { 3, 5, 6 }, { 2, 0, 3 }, { 1, 2, 8 } };
            MathTimeLibr.matinverse(mat1, 3, out matinv);

            strbuild.AppendLine("matinv = " + matinv[0, 0].ToString(fmt).PadLeft(4) + " " + matinv[0, 1].ToString(fmt).PadLeft(4) + " " + matinv[0, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("matinv = " + matinv[1, 0].ToString(fmt).PadLeft(4) + " " + matinv[1, 1].ToString(fmt).PadLeft(4) + " " + matinv[1, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("matinv = " + matinv[2, 0].ToString(fmt).PadLeft(4) + " " + matinv[2, 1].ToString(fmt).PadLeft(4) + " " + matinv[2, 2].ToString(fmt).PadLeft(4));

            //Results: test before
            // 0.1016949    0.4745763 - 0.2542373
            // 0.2203390 - 0.3050847 - 0.0508475
            //- 0.0677966    0.0169492    0.1694915

            mat1 = new double[,] { { 1, 3, 3 }, { 1, 4, 3 }, { 1, 3, 4 } };
            MathTimeLibr.matinverse(mat1, 3, out matinv);

            strbuild.AppendLine("matinv = " + matinv[0, 0].ToString(fmt).PadLeft(4) + " " + matinv[0, 1].ToString(fmt).PadLeft(4) + " " + matinv[0, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("matinv = " + matinv[1, 0].ToString(fmt).PadLeft(4) + " " + matinv[1, 1].ToString(fmt).PadLeft(4) + " " + matinv[1, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("matinv = " + matinv[2, 0].ToString(fmt).PadLeft(4) + " " + matinv[2, 1].ToString(fmt).PadLeft(4) + " " + matinv[2, 2].ToString(fmt).PadLeft(4));


            double[,] ata = new double[,] {{264603537.493561, 206266447.729262, 274546062925.826, -282848493891885, 362835957483807, -4.3758299682612E+17 },
            { 206266447.729262, 160790924.64848, 214016946538.904, -220488942186083, 282841585443473, -3.41109159752805E+17 },
            { 274546062925.826, 214016946538.904, 284862180536440, -2.93476576836794E+17, 3.76469583735348E+17, -4.54025256502168E+20},
            { -282848493891885, -220488942186083, -2.93476576836794E+17, 3.02351477439543E+20, -3.87854240635812E+20, 4.67755241586584E+23},
            { 362835957483807, 282841585443473, 3.76469583735348E+17, -3.87854240635812E+20, 4.97536553328938E+20, -6.00032966815125E+23},
            { -4.3758299682612E+17, -3.41109159752805E+17, -4.54025256502168E+20, 4.67755241586584E+23, -6.00032966815125E+23, 7.23644441510866E+26} };

            MathTimeLibr.matinverse(ata, 6, out matinv);

            strbuild.AppendLine("matinv = " + matinv[0, 0].ToString(fmt1).PadLeft(4) + " " + matinv[0, 1].ToString(fmt1).PadLeft(4) + " " + matinv[0, 2].ToString(fmt1).PadLeft(4) + " " + matinv[0, 3].ToString(fmt1).PadLeft(4) + " " + matinv[0, 4].ToString(fmt1).PadLeft(4) + " " + matinv[0, 5].ToString(fmt1).PadLeft(4));
            strbuild.AppendLine("matinv = " + matinv[1, 0].ToString(fmt1).PadLeft(4) + " " + matinv[1, 1].ToString(fmt1).PadLeft(4) + " " + matinv[1, 2].ToString(fmt1).PadLeft(4) + " " + matinv[1, 3].ToString(fmt1).PadLeft(4) + " " + matinv[1, 4].ToString(fmt1).PadLeft(4) + " " + matinv[1, 5].ToString(fmt1).PadLeft(4));
            strbuild.AppendLine("matinv = " + matinv[2, 0].ToString(fmt1).PadLeft(4) + " " + matinv[2, 1].ToString(fmt1).PadLeft(4) + " " + matinv[2, 2].ToString(fmt1).PadLeft(4) + " " + matinv[2, 3].ToString(fmt1).PadLeft(4) + " " + matinv[2, 4].ToString(fmt1).PadLeft(4) + " " + matinv[2, 5].ToString(fmt1).PadLeft(4));
            strbuild.AppendLine("matinv = " + matinv[3, 0].ToString(fmt1).PadLeft(4) + " " + matinv[3, 1].ToString(fmt1).PadLeft(4) + " " + matinv[3, 2].ToString(fmt1).PadLeft(4) + " " + matinv[3, 3].ToString(fmt1).PadLeft(4) + " " + matinv[3, 4].ToString(fmt1).PadLeft(4) + " " + matinv[3, 5].ToString(fmt1).PadLeft(4));
            strbuild.AppendLine("matinv = " + matinv[4, 0].ToString(fmt1).PadLeft(4) + " " + matinv[4, 1].ToString(fmt1).PadLeft(4) + " " + matinv[4, 2].ToString(fmt1).PadLeft(4) + " " + matinv[4, 3].ToString(fmt1).PadLeft(4) + " " + matinv[4, 4].ToString(fmt1).PadLeft(4) + " " + matinv[4, 5].ToString(fmt1).PadLeft(4));
            strbuild.AppendLine("matinv = " + matinv[5, 0].ToString(fmt1).PadLeft(4) + " " + matinv[5, 1].ToString(fmt1).PadLeft(4) + " " + matinv[5, 2].ToString(fmt1).PadLeft(4) + " " + matinv[5, 3].ToString(fmt1).PadLeft(4) + " " + matinv[5, 4].ToString(fmt1).PadLeft(4) + " " + matinv[5, 5].ToString(fmt1).PadLeft(4));

            double[,] mat3 = MathTimeLibr.matmult(ata, matinv, 6, 6, 6);

            strbuild.AppendLine("mat3 = " + mat3[0, 0].ToString(fmt1).PadLeft(4) + " " + mat3[0, 1].ToString(fmt1).PadLeft(4) + " " + mat3[0, 2].ToString(fmt1).PadLeft(4) + " " + mat3[0, 3].ToString(fmt1).PadLeft(4) + " " + mat3[0, 4].ToString(fmt1).PadLeft(4) + " " + mat3[0, 5].ToString(fmt1).PadLeft(4));
            strbuild.AppendLine("mat3 = " + mat3[1, 0].ToString(fmt1).PadLeft(4) + " " + mat3[1, 1].ToString(fmt1).PadLeft(4) + " " + mat3[1, 2].ToString(fmt1).PadLeft(4) + " " + mat3[1, 3].ToString(fmt1).PadLeft(4) + " " + mat3[1, 4].ToString(fmt1).PadLeft(4) + " " + mat3[1, 5].ToString(fmt1).PadLeft(4));
            strbuild.AppendLine("mat3 = " + mat3[2, 0].ToString(fmt1).PadLeft(4) + " " + mat3[2, 1].ToString(fmt1).PadLeft(4) + " " + mat3[2, 2].ToString(fmt1).PadLeft(4) + " " + mat3[2, 3].ToString(fmt1).PadLeft(4) + " " + mat3[2, 4].ToString(fmt1).PadLeft(4) + " " + mat3[2, 5].ToString(fmt1).PadLeft(4));
            strbuild.AppendLine("mat3 = " + mat3[3, 0].ToString(fmt1).PadLeft(4) + " " + mat3[3, 1].ToString(fmt1).PadLeft(4) + " " + mat3[3, 2].ToString(fmt1).PadLeft(4) + " " + mat3[3, 3].ToString(fmt1).PadLeft(4) + " " + mat3[3, 4].ToString(fmt1).PadLeft(4) + " " + mat3[3, 5].ToString(fmt1).PadLeft(4));
            strbuild.AppendLine("mat3 = " + mat3[4, 0].ToString(fmt1).PadLeft(4) + " " + mat3[4, 1].ToString(fmt1).PadLeft(4) + " " + mat3[4, 2].ToString(fmt1).PadLeft(4) + " " + mat3[4, 3].ToString(fmt1).PadLeft(4) + " " + mat3[4, 4].ToString(fmt1).PadLeft(4) + " " + mat3[4, 5].ToString(fmt1).PadLeft(4));
            strbuild.AppendLine("mat3 = " + mat3[5, 0].ToString(fmt1).PadLeft(4) + " " + mat3[5, 1].ToString(fmt1).PadLeft(4) + " " + mat3[5, 2].ToString(fmt1).PadLeft(4) + " " + mat3[5, 3].ToString(fmt1).PadLeft(4) + " " + mat3[5, 4].ToString(fmt1).PadLeft(4) + " " + mat3[5, 5].ToString(fmt1).PadLeft(4));
        }

        public void testdeterminant()
        {
            double[,] mat1 = new double[3, 3];
            double det;
            int order;

            order = 3;

            mat1[0, 0] = 6.0;
            mat1[0, 1] = 1.0;
            mat1[0, 2] = 1.0;
            mat1[1, 0] = 4.0;
            mat1[1, 1] = -2.0;
            mat1[1, 2] = 5.0;
            mat1[2, 0] = 2.0;
            mat1[2, 1] = 8.0;
            mat1[2, 2] = 7.0;

            det = MathTimeLibr.determinant(mat1, order);

            strbuild.AppendLine("det = " + det.ToString(fmt).PadLeft(4) + " ans -306");
        }

        public void testcholesky()
        {
            double[,] mat1 = new double[3, 3];
            double[,] a = new double[3, 3];

            a[0, 0] = 1.0;
            a[0, 1] = 0.0;
            a[0, 2] = 1.0;
            a[1, 0] = 0.0;
            a[1, 1] = 2.0;
            a[1, 2] = 0.0;
            a[2, 0] = 1.0;
            a[2, 1] = 0.0;
            a[2, 2] = 3.0;

            mat1 = MathTimeLibr.cholesky(a);

            strbuild.AppendLine("matcho = " + mat1[0, 0].ToString(fmt1).PadLeft(4) + " " + mat1[0, 1].ToString(fmt1).PadLeft(4) + " " + mat1[0, 2].ToString(fmt1).PadLeft(4));
            strbuild.AppendLine("matcho = " + mat1[1, 0].ToString(fmt1).PadLeft(4) + " " + mat1[1, 1].ToString(fmt1).PadLeft(4) + " " + mat1[1, 2].ToString(fmt1).PadLeft(4));
            strbuild.AppendLine("matcho = " + mat1[2, 0].ToString(fmt1).PadLeft(4) + " " + mat1[2, 1].ToString(fmt1).PadLeft(4) + " " + mat1[2, 2].ToString(fmt1).PadLeft(4));
        }

        public void testposvelcov2pts()
        {
            double[] reci = new double[3];
            double[] veci = new double[3];
            double[,] cov = new double[6, 6];
            double[,] sigmapts = new double[6, 12];

            reci = new double[] { 5102.50895791863, 6123.01140072233, 6378.13692818659 };
            veci = new double[] { -4.74322014817, 0.79053648924, 5.53375572723 };

            cov[0, 0] = 12559.93762571587;
            cov[0, 1] = cov[1, 0] = 12101.56371305036;
            cov[0, 2] = cov[2, 0] = -440.3145384949657;
            cov[0, 3] = cov[3, 0] = -0.8507401236198346;
            cov[0, 4] = cov[4, 0] = 0.9383675791981778;
            cov[0, 5] = cov[5, 0] = -0.0318596430999798;
            cov[1, 1] = 12017.77368889201;
            cov[1, 2] = cov[2, 1] = 270.3798093532698;
            cov[1, 3] = cov[3, 1] = -0.8239662300032132;
            cov[1, 4] = cov[4, 1] = 0.9321640899868708;
            cov[1, 5] = cov[5, 1] = -0.001327326827629336;
            cov[2, 2] = 4818.009967057008;
            cov[2, 3] = cov[3, 2] = 0.02033418761460195;
            cov[2, 4] = cov[4, 2] = 0.03077663516695039;
            cov[2, 5] = cov[5, 2] = 0.1977541628188323;
            cov[3, 3] = 5.774758755889862e-005;
            cov[3, 4] = cov[4, 3] = -6.396031584925255e-005;
            cov[3, 5] = cov[5, 3] = 1.079960679599204e-006;
            cov[4, 4] = 7.24599391355188e-005;
            cov[4, 5] = cov[5, 4] = 1.03146660433274e-006;
            cov[5, 5] = 1.870413627417302e-005;

            AstroLibr.posvelcov2pts(reci, veci, cov, out sigmapts);

            strbuild.AppendLine("sigmapts = " + sigmapts[0, 0].ToString(fmt1).PadLeft(4) + " " + sigmapts[0, 1].ToString(fmt1).PadLeft(4) + " " + sigmapts[0, 2].ToString(fmt1).PadLeft(4) + " " + sigmapts[0, 3].ToString(fmt1).PadLeft(4) + " " + sigmapts[0, 4].ToString(fmt1).PadLeft(4) + " " + sigmapts[0, 5].ToString(fmt1).PadLeft(4));
            strbuild.AppendLine("sigmapts = " + sigmapts[1, 0].ToString(fmt1).PadLeft(4) + " " + sigmapts[1, 1].ToString(fmt1).PadLeft(4) + " " + sigmapts[1, 2].ToString(fmt1).PadLeft(4) + " " + sigmapts[1, 3].ToString(fmt1).PadLeft(4) + " " + sigmapts[1, 4].ToString(fmt1).PadLeft(4) + " " + sigmapts[1, 5].ToString(fmt1).PadLeft(4));
            strbuild.AppendLine("sigmapts = " + sigmapts[2, 0].ToString(fmt1).PadLeft(4) + " " + sigmapts[2, 1].ToString(fmt1).PadLeft(4) + " " + sigmapts[2, 2].ToString(fmt1).PadLeft(4) + " " + sigmapts[2, 3].ToString(fmt1).PadLeft(4) + " " + sigmapts[2, 4].ToString(fmt1).PadLeft(4) + " " + sigmapts[2, 5].ToString(fmt1).PadLeft(4));
            strbuild.AppendLine("sigmapts = " + sigmapts[3, 0].ToString(fmt1).PadLeft(4) + " " + sigmapts[3, 1].ToString(fmt1).PadLeft(4) + " " + sigmapts[3, 2].ToString(fmt1).PadLeft(4) + " " + sigmapts[3, 3].ToString(fmt1).PadLeft(4) + " " + sigmapts[3, 4].ToString(fmt1).PadLeft(4) + " " + sigmapts[3, 5].ToString(fmt1).PadLeft(4));
            strbuild.AppendLine("sigmapts = " + sigmapts[4, 0].ToString(fmt1).PadLeft(4) + " " + sigmapts[4, 1].ToString(fmt1).PadLeft(4) + " " + sigmapts[4, 2].ToString(fmt1).PadLeft(4) + " " + sigmapts[4, 3].ToString(fmt1).PadLeft(4) + " " + sigmapts[4, 4].ToString(fmt1).PadLeft(4) + " " + sigmapts[4, 5].ToString(fmt1).PadLeft(4));
            strbuild.AppendLine("sigmapts = " + sigmapts[5, 0].ToString(fmt1).PadLeft(4) + " " + sigmapts[5, 1].ToString(fmt1).PadLeft(4) + " " + sigmapts[5, 2].ToString(fmt1).PadLeft(4) + " " + sigmapts[5, 3].ToString(fmt1).PadLeft(4) + " " + sigmapts[5, 4].ToString(fmt1).PadLeft(4) + " " + sigmapts[5, 5].ToString(fmt1).PadLeft(4));
        }

        public void testposcov2pts()
        {
            double[,] cov2 = new double[6, 6];
            double[,] cov = new double[6, 6];
            double[,] sigmapts = new double[6, 12];
            double[] yu = new double[6];
            double[,] covout = new double[6, 6];
            double[] r1 = new double[3];
            double[] v1 = new double[3];
            double[] reci = new double[3];

            r1 = new double[] { 5102.50895791863, 6123.01140072233, 6378.13692818659 };
            v1 = new double[] { -4.74322014817, 0.79053648924, 5.53375572723 };

            cov2[0, 0] = 12559.93762571587;
            cov2[0, 1] = cov2[1, 0] = 12101.56371305036;
            cov2[0, 2] = cov2[2, 0] = -440.3145384949657;
            cov2[0, 3] = cov2[3, 0] = -0.8507401236198346;
            cov2[1, 1] = 12017.77368889201;
            cov2[1, 2] = cov2[2, 1] = 270.3798093532698;
            cov2[1, 3] = cov2[3, 1] = -0.8239662300032132;
            cov2[2, 2] = 4818.009967057008;
            cov2[2, 3] = cov2[3, 2] = 0.02033418761460195;
            cov2[3, 3] = 5.774758755889862e-005;

            // form sigmapts pos/vel
            AstroLibr.poscov2pts(r1, cov2, out sigmapts);
            strbuild.AppendLine("sigmapts = " + sigmapts[0, 0].ToString(fmt1).PadLeft(4) + " " + sigmapts[0, 1].ToString(fmt1).PadLeft(4) + " " + sigmapts[0, 2].ToString(fmt1).PadLeft(4) + " " + sigmapts[0, 3].ToString(fmt1).PadLeft(4) + " " + sigmapts[0, 4].ToString(fmt1).PadLeft(4) + " " + sigmapts[0, 5].ToString(fmt1).PadLeft(4));
            strbuild.AppendLine("sigmapts = " + sigmapts[1, 0].ToString(fmt1).PadLeft(4) + " " + sigmapts[1, 1].ToString(fmt1).PadLeft(4) + " " + sigmapts[1, 2].ToString(fmt1).PadLeft(4) + " " + sigmapts[1, 3].ToString(fmt1).PadLeft(4) + " " + sigmapts[1, 4].ToString(fmt1).PadLeft(4) + " " + sigmapts[1, 5].ToString(fmt1).PadLeft(4));
            strbuild.AppendLine("sigmapts = " + sigmapts[2, 0].ToString(fmt1).PadLeft(4) + " " + sigmapts[2, 1].ToString(fmt1).PadLeft(4) + " " + sigmapts[2, 2].ToString(fmt1).PadLeft(4) + " " + sigmapts[2, 3].ToString(fmt1).PadLeft(4) + " " + sigmapts[2, 4].ToString(fmt1).PadLeft(4) + " " + sigmapts[2, 5].ToString(fmt1).PadLeft(4));

            // reassemble covariance at each step and write out
            //AstroLibr.remakecovpv(sigmapts, out yu, out covout);
            //AstroLibr.poscov2pts(reci, cov, out sigmapts);
        }

        public void testremakecovpv()
        {
            double[,] sigmapts = new double[6, 12];
            double[] yu = new double[6];
            double[,] cov = new double[6, 6];
            double[,] cov2 = new double[6, 6];

            sigmapts = new double[3, 6] { { 5377.0260353, 4827.9918806, 5102.5089579, 5102.5089579, 5102.5089579, 5102.5089579 },
                { 6123.0114007, 6123.0114007, 6391.5382004, 5854.4846011, 6123.0114007, 6123.0114007},
                { 6378.1369282, 6378.1369282, 6378.1369282, 6378.1369282, 6548.1606318, 6208.1132245} };

     //       AstroLibr.remakecovpv(sigmapts, out yu, out cov);

            strbuild.AppendLine("cov = " + cov[0, 0].ToString(fmt1).PadLeft(4) + " " + cov[0, 1].ToString(fmt1).PadLeft(4) + " " + cov[0, 2].ToString(fmt1).PadLeft(4));
            strbuild.AppendLine("cov = " + cov[1, 0].ToString(fmt1).PadLeft(4) + " " + cov[1, 1].ToString(fmt1).PadLeft(4) + " " + cov[1, 2].ToString(fmt1).PadLeft(4));
            strbuild.AppendLine("cov = " + cov[2, 0].ToString(fmt1).PadLeft(4) + " " + cov[2, 1].ToString(fmt1).PadLeft(4) + " " + cov[2, 2].ToString(fmt1).PadLeft(4));
        }
        public void testremakecovp()
        {
            double[,] sigmapts = new double[6, 12];
            double[] yu = new double[6];
            double[,] cov = new double[6, 6];

            sigmapts = new double[3, 6] { { 5377.0260353, 4827.9918806, 5102.5089579, 5102.5089579, 5102.5089579, 5102.5089579 },
                { 6123.0114007, 6123.0114007, 6391.5382004, 5854.4846011, 6123.0114007, 6123.0114007},
                { 6378.1369282, 6378.1369282, 6378.1369282, 6378.1369282, 6548.1606318, 6208.1132245} };

//            AstroLibr.remakecovp(sigmapts, out yu, out cov);

            strbuild.AppendLine("cov = " + cov[0, 0].ToString(fmt1).PadLeft(4) + " " + cov[0, 1].ToString(fmt1).PadLeft(4) + " " + cov[0, 2].ToString(fmt1).PadLeft(4));
            strbuild.AppendLine("cov = " + cov[1, 0].ToString(fmt1).PadLeft(4) + " " + cov[1, 1].ToString(fmt1).PadLeft(4) + " " + cov[1, 2].ToString(fmt1).PadLeft(4));
            strbuild.AppendLine("cov = " + cov[2, 0].ToString(fmt1).PadLeft(4) + " " + cov[2, 1].ToString(fmt1).PadLeft(4) + " " + cov[2, 2].ToString(fmt1).PadLeft(4));
        }
        public void testmatequal()
        {
            double[,] mat1 = new double[3, 3];
            double[,] mat3 = new double[3, 3];
            int matr;
            matr = 3;
            mat1 = new double[,] { { 1, 3, 3 }, { 1, 4, 3 }, { 1, 3, 4 } };

            mat3 = MathTimeLibr.matequal(mat1, matr);

            strbuild.AppendLine("matequal = " + mat3[0, 0].ToString(fmt).PadLeft(4) + " " + mat3[0, 1].ToString(fmt).PadLeft(4) + " " + mat3[0, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("matequal = " + mat3[1, 0].ToString(fmt).PadLeft(4) + " " + mat3[1, 1].ToString(fmt).PadLeft(4) + " " + mat3[1, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("matequal = " + mat3[2, 0].ToString(fmt).PadLeft(4) + " " + mat3[2, 1].ToString(fmt).PadLeft(4) + " " + mat3[2, 2].ToString(fmt).PadLeft(4));
        }
        public void testmatscale()
        {
            double[,] mat1 = new double[3, 3];
            double[,] mat3 = new double[3, 3];
            int matr, matc;
            double scale;
            matr = 3;
            matc = 3;
            scale = 1.364;
            mat1 = new double[,] { { 1, 3, 3 }, { 1, 4, 3 }, { 1, 3, 4 } };

            mat3 = MathTimeLibr.matscale(mat1, matr, matc, scale);

            strbuild.AppendLine("matscale = " + mat3[0, 0].ToString(fmt).PadLeft(4) + " " + mat3[0, 1].ToString(fmt).PadLeft(4) + " " + mat3[0, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("matscale = " + mat3[1, 0].ToString(fmt).PadLeft(4) + " " + mat3[1, 1].ToString(fmt).PadLeft(4) + " " + mat3[1, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("matscale = " + mat3[2, 0].ToString(fmt).PadLeft(4) + " " + mat3[2, 1].ToString(fmt).PadLeft(4) + " " + mat3[2, 2].ToString(fmt).PadLeft(4));
        }
        public void testnorm()
        {
            double[] vec1 = new double[3];
            double[] vec2 = new double[3];
            vec1 = new double[3] { 2.3, 4.7, -1.6 };

            vec2 = MathTimeLibr.norm(vec1);

            strbuild.AppendLine("norm = " + vec2[0].ToString(fmt).PadLeft(4) + " " + vec2[1].ToString(fmt).PadLeft(4) + " " + vec2[2].ToString(fmt).PadLeft(4));
        }
        public void testmag()
        {
            double[] x = new double[3];
            double magx;
            x = new double[3] { 1.0, 2.0, 5.0 };

            magx = MathTimeLibr.mag(x);

            strbuild.AppendLine("mag = " + magx.ToString(fmt).PadLeft(4));
        }
        public void testcross()
        {
            double[] vec1 = new double[3];
            double[] vec2 = new double[3];
            double[] outvec = new double[3];
            vec1 = new double[3] { 1.0, 2.0, 5.0 };
            vec2 = new double[3] { 2.3, 4.7, -1.6 };

            MathTimeLibr.cross(vec1, vec2, out outvec);

            strbuild.AppendLine("cross = " + outvec[0].ToString(fmt).PadLeft(4) + " " + outvec[1].ToString(fmt).PadLeft(4) + " " + outvec[2].ToString(fmt).PadLeft(4));
        }
        public void testdot()
        {
            double[] x = new double[3];
            double[] y = new double[3];
            double dotp;
            x = new double[3] { 1, 2, 5 };
            y = new double[3] { 2.3, 4.7, -1.6 };

            dotp = MathTimeLibr.dot(x, y);

            strbuild.AppendLine("x " + x[0].ToString(fmt).PadLeft(4) + " " + x[1].ToString(fmt).PadLeft(4) + " " + x[2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("y " + y[0].ToString(fmt).PadLeft(4) + " " + y[1].ToString(fmt).PadLeft(4) + " " + y[2].ToString(fmt).PadLeft(4));

            strbuild.AppendLine("dot = " + dotp.ToString(fmt).PadLeft(4));
        }

        public void testangle()
        {
            double[] vec1 = new double[3];
            double[] vec2 = new double[3];
            double ang;
            vec1 = new double[3] { 1, 2, 5 };
            vec2 = new double[3] { 2.3, 4.7, -1.6 };

            ang = MathTimeLibr.angle(vec1, vec2);

            strbuild.AppendLine("angle = " + ang.ToString(fmt).PadLeft(4));
        }
        public void testasinh()
        {
            double xval, ans;
            xval = 1.45;

            ans = MathTimeLibr.asinh(xval);

            strbuild.AppendLine("asinh = " + ans.ToString(fmt).PadLeft(4));
        }
        public void testcot()
        {
            double xval, ans;
            xval = 0.47238734;

            ans = MathTimeLibr.cot(xval);

            strbuild.AppendLine("cot = " + ans.ToString(fmt).PadLeft(4));
        }
        public void testacosh()
        {
            double xval, ans;
            xval = 1.43;

            ans = MathTimeLibr.acosh(xval);

            strbuild.AppendLine("acosh = " + ans.ToString(fmt).PadLeft(4));
        }
        public void testaddvec()
        {
            double a1, a2;
            double[] vec1 = new double[3];
            double[] vec2 = new double[3];
            double[] vec3 = new double[3];
            vec1 = new double[3] { 1, 2, 5 };
            vec2 = new double[3] { 2.3, 4.7, -5.6 };
            a1 = 1.0;
            a2 = 2.0;

            MathTimeLibr.addvec(a1, vec1, a2, vec2, out vec3);

            strbuild.AppendLine("vec1 " + vec1[0].ToString(fmt).PadLeft(4) + " " + vec1[1].ToString(fmt).PadLeft(4) + " " + vec1[2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("vec2 " + vec2[0].ToString(fmt).PadLeft(4) + " " + vec2[1].ToString(fmt).PadLeft(4) + " " + vec2[2].ToString(fmt).PadLeft(4));

            strbuild.AppendLine("addvec = " + vec3[0].ToString(fmt).PadLeft(4) + " " + vec3[1].ToString(fmt).PadLeft(4) + " " + vec3[2].ToString(fmt).PadLeft(4));
        }


        public void testPercentile()
        {
            double ans;
            double[] sequence = new double[15];
            double excelPercentile;
            Int32 arrSize;
            excelPercentile = 0.3;
            arrSize = 7;
            sequence[0] = 45.3;
            sequence[1] = 5.63;
            sequence[2] = 5.13;
            sequence[3] = 345.3;
            sequence[4] = 45.3;
            sequence[5] = 3445.3;
            sequence[6] = 0.03;

            ans = MathTimeLibr.Percentile(sequence, excelPercentile, arrSize);

            strbuild.AppendLine("percentile = " + ans.ToString(fmt).PadLeft(4));
        }
        public void testrot1()
        {
            double[] vec = new double[3];
            double[] vec3 = new double[3];
            double xval, rad;
            rad = 180.0 / Math.PI;
            vec[0] = 23.4;
            vec[1] = 6723.4;
            vec[2] = -2.4;
            xval = 225.0 / rad;

            vec3 = MathTimeLibr.rot1(vec, xval);

            strbuild.AppendLine("testrot1 = " + vec3[0].ToString(fmt).PadLeft(4) + " " + vec3[1].ToString(fmt).PadLeft(4) + " " + vec3[2].ToString(fmt).PadLeft(4));
        }
        public void testrot2()
        {
            double[] vec3 = new double[3];
            double[] vec = new double[3];
            double xval, rad;
            rad = 180.0 / Math.PI;
            vec[0] = 23.4;
            vec[1] = 6723.4;
            vec[2] = -2.4;
            xval = 23.4 / rad;
            vec3 = MathTimeLibr.rot2(vec, xval);

            strbuild.AppendLine("testrot2 = " + vec3[0].ToString(fmt).PadLeft(4) + " " + vec3[1].ToString(fmt).PadLeft(4) + " " + vec3[2].ToString(fmt).PadLeft(4));
        }
        public void testrot3()
        {
            double[] vec3 = new double[3];
            double[] vec = new double[3];
            double xval, rad;
            rad = 180.0 / Math.PI;
            vec[0] = 23.4;
            vec[1] = 6723.4;
            vec[2] = -2.4;
            xval = 323.4 / rad;
            vec3 = MathTimeLibr.rot3(vec, xval);

            strbuild.AppendLine("testrot3 = " + vec3[0].ToString(fmt).PadLeft(4) + " " + vec3[1].ToString(fmt).PadLeft(4) + " " + vec3[2].ToString(fmt).PadLeft(4));
        }
        public void testfactorial()
        {
            Int32 n;
            double ans;
            n = 4;

            ans = MathTimeLibr.factorial(n);

            strbuild.AppendLine("factorial = " + ans.ToString().PadLeft(4));
        }
        public void testcubicspl()
        {
            double p1, p2, p3, p4, acu0, acu1, acu2, acu3;
            p1 = 1.0;
            p2 = 3.5;
            p3 = 5.6;
            p4 = 32.0;

            MathTimeLibr.cubicspl(p1, p2, p3, p4, out acu0, out acu1, out acu2, out acu3);

            strbuild.AppendLine("cubicspl = " + acu0.ToString(fmt).PadLeft(11) + acu1.ToString(fmt).PadLeft(11)
                + acu2.ToString(fmt).PadLeft(11) + acu3.ToString(fmt).PadLeft(11));
        }
        public void testcubic()
        {
            double a3, b2, c1, d0;
            char opt;
            double r1r, r1i, r2r, r2i, r3r, r3i;
            a3 = 1.7;
            b2 = 3.5;
            c1 = 5.6;
            d0 = 32.0;
            opt = 'I';  // all roots, unique, real

            MathTimeLibr.cubic(a3, b2, c1, d0, opt, out r1r, out r1i, out r2r, out r2i, out r3r, out r3i);

            strbuild.AppendLine("cubic = " + r1r.ToString(fmt).PadLeft(11) + r1i.ToString(fmt).PadLeft(11)
                + r2r.ToString(fmt).PadLeft(11) + r2i.ToString(fmt).PadLeft(11)
                + r3r.ToString(fmt).PadLeft(11) + r3i.ToString(fmt).PadLeft(11));
        }
        public void testcubicinterp()
        {
            double p1a, p1b, p1c, p1d, p2a, p2b, p2c, p2d, valuein;
            double ans;
            p1a = 1.7;
            p1b = 3.5;
            p1c = 5.6;
            p1d = 11.7;
            p2a = 21.7;
            p2b = 35.5;
            p2c = 57.6;
            p2d = 181.7;
            valuein = 4.0;

            ans = MathTimeLibr.cubicinterp(p1a, p1b, p1c, p1d, p2a, p2b, p2c, p2d, valuein);

            strbuild.AppendLine("cubicint = " + ans.ToString(fmt).PadLeft(11));
        }

        public void testquadratic()
        {
            double a, b, c;
            char opt;
            double r1r, r1i, r2r, r2i;
            a = 1.7;
            b = 3.5;
            c = 5.6;
            opt = 'I';  // all roots, unique, real

            MathTimeLibr.quadratic(a, b, c, opt, out r1r, out r1i, out r2r, out r2i);

            strbuild.AppendLine("quad = " + r1r.ToString(fmt).PadLeft(11) + r1i.ToString(fmt).PadLeft(11)
                + r2r.ToString(fmt).PadLeft(11) + r2i.ToString(fmt).PadLeft(11));
        }

        public void testconvertMonth()
        {
            string monstr;
            monstr = "Jan";
            int mon;

            mon = MathTimeLibr.getIntMonth(monstr);
        }

        public void testjday()
        {
            int year;
            double jd, jdFrac;
            int mon, day, hr, minute;
            double second;
            year = 2020;
            mon = 12;
            day = 15;
            hr = 16;
            minute = 58;
            second = 50.208;
            MathTimeLibr.jday(year, mon, day, hr, minute, second, out jd, out jdFrac);

            strbuild.AppendLine("jd " + jd.ToString(fmt).PadLeft(4) + " " + jdFrac.ToString(fmt).PadLeft(4));


            // alt tests
            MathTimeLibr.invjday(2450382.5, jdFrac, out year, out mon, out day, out hr, out minute, out second);
            strbuild.AppendLine("year" + year + " mon " + mon + " day " + day + " " + hr + ":" + minute + ":" + second);

            MathTimeLibr.invjday(2450382.5, jdFrac -0.2, out year, out mon, out day, out hr, out minute, out second);
            strbuild.AppendLine("year" + year + " mon " + mon + " day " + day + " " + hr + ":" + minute + ":" + second);

            MathTimeLibr.invjday(2450382.5 + 1.0, jdFrac + 1.5, out year, out mon, out day, out hr, out minute, out second);
            strbuild.AppendLine("year" + year + " mon " + mon + " day " + day + " " + hr + ":" + minute + ":" + second);

            MathTimeLibr.invjday(2450382.5, -0.5, out year, out mon, out day, out hr, out minute, out second);
            strbuild.AppendLine("year" + year + " mon " + mon + " day " + day+ " " + hr + ":" + minute + ":" + second);

            MathTimeLibr.invjday(2450382.5, 0.5, out year, out mon, out day, out hr, out minute, out second);
            strbuild.AppendLine("year" + year + " mon " + mon + " day " + day + " " + hr + ":" + minute + ":" + second);

        }
        public void testdays2mdhms()
        {
            int year;
            double days;
            int mon, day, hr, minute;
            double second;
            year = 2020;
            days = 237.456982367;

            MathTimeLibr.days2mdhms(year, days, out mon, out day, out hr, out minute, out second);

            strbuild.AppendLine("year" + year + " days " + days);
            strbuild.AppendLine("year" + year + " mon " + mon + " day " + day + hr + ":" + minute + ":" + second);
        }
        public void testinvjday()
        {
            int year;
            int mon, day, hr, minute;
            double jd, jdF, jdb, mfme, dt, second;
            jd = 2449877.0;
            jdF = 0.3458762;
            mfme = 0.0;
            dt = 0.0;

            MathTimeLibr.invjday(jd, jdF, out year, out mon, out day, out hr, out minute, out second);

            strbuild.AppendLine("year" + year + " mon " + mon + " day " + day + hr + ":" + minute + ":" + second);


            // some stressing cases
            int i;
            for (i = 1; i < 11; i++)
            {
                if (i == 1)
                {
                    jd = 2457884.5;
                    jdF = 0.160911640046296;
                }
                if (i == 2)
                {
                    jd = 2457884.5;
                    jdF = -0.160911640046296;
                }
                if (i == 3)
                {
                    jd = 2457884.5;
                    jdF = 0.660911640046296;
                }
                if (i == 4)
                {
                    jd = 2457884.5;
                    jdF = -0.660911640046296;
                }
                if (i == 5)
                {
                    jd = 2457884.5;
                    jdF = 2.160911640046296;
                }
                if (i == 6)
                {
                    jd = 2457884.660911640046296;
                    jdF = 0.0;
                }
                if (i == 7)
                {
                    jd = 2457884.0;
                    jdF = 2.160911640046296;
                }
                if (i == 8)
                {
                    jd = 2457884.5;
                    jdF = 0.0;
                }
                if (i == 9)
                {
                    jd = 2457884.5;
                    jdF = 0.5;
                }
                if (i == 10)
                {
                    jd = 2457884.5;
                    jdF = 1.0;
                }
                if (i == 11)
                {
                    jd = 2457884.3;
                    jdF = 1.0;
                }

                jdb = Math.Floor(jd + jdF) + 0.5;
                mfme = (jd + jdF - jdb) * 1440.0;
                if (mfme < 0.0)
                {
                    mfme = 1440.0 + mfme;
                }
                MathTimeLibr.invjday(jd, jdF, out year, out mon, out day, out hr, out minute, out second);
                if (Math.Abs(jdF) >= 1.0)
                {
                    jd = jd + Math.Floor(jdF);
                    jdF = jdF - Math.Floor(jdF);
                }
                dt = jd - Math.Floor(jd) - 0.5;
                if (Math.Abs(dt) > 0.00000001)
                {
                    jd = jd - dt;
                    jdF = jdF + dt;
                }
                // this gets it even to the day
                if (jdF < 0.0)
                {
                    jd = jd - 1.0;
                    jdF = 1.0 + jdF;
                }
                strbuild.AppendLine(i + " " + year + " " + mon + " " + day + " " + hr + ":" + minute + ":" + second + " " +
                     mfme +" " + hr * 60.0 + " "+ minute + " " + second / 60.0 + " " + jd + " " + jdF + " " + dt);
            }
        }   // through stressing cases


        // tests eop, spw, and fk5 iau80
        public void testiau80in()
        {
            Int32 year, mon, day, hr, minute, dat;
            double jd, jdFrac, second, dut1, lod, xp, yp, ddx, ddy, ddpsi, ddeps;
            int i;
            double f107, f107bar, ap, avgap, kp, sumkp;
            double[] aparr = new double[8];
            double[] kparr = new double[8];
            Int32 ktrActObs;
            string EOPupdate;

            year = 2017;
            mon = 4;
            day = 6;
            hr = 0;
            minute = 0;
            second = 0.0;
            // 2017 04 06 57849  0.009587  0.384603  0.4636260  0.0014635 -0.098066 -0.012276 -0.000063  0.000300  37 

            string fileLoc;
            fileLoc = @"D:\Codes\LIBRARY\DataLib\nut80.dat";
            EOPSPWLibr.iau80in(fileLoc, out EOPSPWLibr.iau80arr);

            string eopFileName = @"D:\Codes\LIBRARY\DataLib\EOP-All-v1.1_2018-01-04.txt";
            EOPSPWLibr.readeop(ref EOPSPWLibr.eopdata, eopFileName,  out ktrActObs, out EOPupdate);
            int y, m, d, h, mm;
            double ss;
            strbuild.AppendLine("EOP tests  mfme    dut1  dat    lod           xp                      yp               ddpsi                   ddeps               ddx                 ddy");
            for (i = 0; i < 90; i++)
            {
                MathTimeLibr.jday(year, mon, day, hr + i, minute, second, out jd, out jdFrac);
                EOPSPWLibr.findeopparam(jd, jdFrac, 's', EOPSPWLibr.eopdata, out dut1, out dat,
                   out lod, out xp, out yp, out ddpsi, out ddeps, out ddx, out ddy);
                MathTimeLibr.invjday(jd, jdFrac, out y, out m, out d, out h, out mm, out ss);
                strbuild.AppendLine(y.ToString("0000") + " " + m.ToString("00") + " " + d.ToString("00") + " " + (h * 60 + mm).ToString("0000") + " " +
                    dut1.ToString(fmt).PadLeft(4) + " " + dat.ToString("00").PadLeft(4) + " " + lod.ToString(fmt).PadLeft(4) + " " + xp.ToString(fmtE).PadLeft(4) + " " + yp.ToString(fmtE).PadLeft(4) + " " +
                    ddpsi.ToString(fmtE).PadLeft(4) + " " + ddeps.ToString(fmtE).PadLeft(4) + " " + ddx.ToString(fmtE).PadLeft(4) + " " + ddy.ToString(fmtE).PadLeft(4));
            }

            string spwFileName = @"D:\Codes\LIBRARY\DataLib\SpaceWeather-All-v1.2_2018-01-04.txt";
            string errstr;
            EOPSPWLibr.readspw(ref EOPSPWLibr.spwdata, spwFileName, out ktrActObs, out errstr);
            strbuild.AppendLine("SPW tests  mfme f107 f107bar ap apavg  kp sumkp aparr[]  ");
            for (i = 0; i < 45; i++)
            {
                MathTimeLibr.jday(year, mon, day, hr + i, minute, second, out jd, out jdFrac);
                //                              interp, adj obs, last ctr, actu const
                EOPSPWLibr.findspwparam(jd, jdFrac, 's', 'a', 'l', 'a', EOPSPWLibr.spwdata, out f107, out f107bar,
                   out ap, out avgap, aparr, out kp, out sumkp, kparr);
                MathTimeLibr.invjday(jd, jdFrac, out y, out m, out d, out h, out mm, out ss);
                strbuild.AppendLine("sala " + y.ToString("0000") + " " + m.ToString("00") + " " + d.ToString("00") + " " + (h * 60 + mm).ToString("0000") + " " +
                   f107.ToString(fmt).PadLeft(4) + " " + f107bar.ToString(fmt).PadLeft(4) + " " + ap.ToString(fmt).PadLeft(4) + " " + avgap.ToString(fmt).PadLeft(4) + " " + kp.ToString(fmt).PadLeft(4) + " " +
                   sumkp.ToString(fmt).PadLeft(4) + " " + aparr[0].ToString(fmt).PadLeft(4) + " " + aparr[1].ToString(fmt).PadLeft(4) + " " + aparr[2].ToString(fmt).PadLeft(4));

                EOPSPWLibr.findspwparam(jd, jdFrac, 's', 'o', 'c', 'a', EOPSPWLibr.spwdata, out f107, out f107bar,
                   out ap, out avgap, aparr, out kp, out sumkp, kparr);
                MathTimeLibr.invjday(jd, jdFrac, out y, out m, out d, out h, out mm, out ss);
                strbuild.AppendLine("soca " + y.ToString("0000") + " " + m.ToString("00") + " " + d.ToString("00") + " " + (h * 60 + mm).ToString("0000") + " " +
                   f107.ToString(fmt).PadLeft(4) + " " + f107bar.ToString(fmt).PadLeft(4) + " " + ap.ToString(fmt).PadLeft(4) + " " + avgap.ToString(fmt).PadLeft(4) + " " + kp.ToString(fmt).PadLeft(4) + " " +
                   sumkp.ToString(fmt).PadLeft(4) + " " + aparr[0].ToString(fmt).PadLeft(4) + " " + aparr[1].ToString(fmt).PadLeft(4) + " " + aparr[2].ToString(fmt).PadLeft(4));
            }
        }

        public void testfundarg()
        {
            double[] fArgs = new double[14];
            double ttt;
            AstroLib.EOpt opt = AstroLib.EOpt.e80;
            ttt = 0.042623631889;

            AstroLibr.fundarg(ttt, opt, out fArgs);
            strbuild.AppendLine("fundarg = " + fArgs[0].ToString(fmt).PadLeft(4) + " " + fArgs[1].ToString(fmt).PadLeft(4) + " "
                + fArgs[2].ToString(fmt).PadLeft(4) + " " + fArgs[3].ToString(fmt).PadLeft(4) + " " + fArgs[4].ToString(fmt).PadLeft(4) + " "
                + fArgs[5].ToString(fmt).PadLeft(4) + " " + fArgs[6].ToString(fmt).PadLeft(4) + " " + fArgs[7].ToString(fmt).PadLeft(4) + " "
                + fArgs[8].ToString(fmt).PadLeft(4) + " " + fArgs[9].ToString(fmt).PadLeft(4) + " " + fArgs[10].ToString(fmt).PadLeft(4) + " "
                + fArgs[11].ToString(fmt).PadLeft(4) + " " + fArgs[12].ToString(fmt).PadLeft(4) + " " + fArgs[13].ToString(fmt).PadLeft(4));

            // do in deg
            for (int i = 0; i < 14; i++)
                fArgs[i] = fArgs[i] * 180.0 / Math.PI;

            strbuild.AppendLine("fundarg = " + fArgs[0].ToString(fmt).PadLeft(4) + " " + fArgs[1].ToString(fmt).PadLeft(4) + " "
                + fArgs[2].ToString(fmt).PadLeft(4) + " " + fArgs[3].ToString(fmt).PadLeft(4) + " " + fArgs[4].ToString(fmt).PadLeft(4) + " "
                + fArgs[5].ToString(fmt).PadLeft(4) + " " + fArgs[6].ToString(fmt).PadLeft(4) + " " + fArgs[7].ToString(fmt).PadLeft(4) + " "
                + fArgs[8].ToString(fmt).PadLeft(4) + " " + fArgs[9].ToString(fmt).PadLeft(4) + " " + fArgs[10].ToString(fmt).PadLeft(4) + " "
                + fArgs[11].ToString(fmt).PadLeft(4) + " " + fArgs[12].ToString(fmt).PadLeft(4) + " " + fArgs[13].ToString(fmt).PadLeft(4));


            AstroLibr.fundarg(ttt, AstroLib.EOpt.e06cio, out fArgs);
            strbuild.AppendLine("fundarg = " + fArgs[0].ToString(fmt).PadLeft(4) + " " + fArgs[1].ToString(fmt).PadLeft(4) + " "
                + fArgs[2].ToString(fmt).PadLeft(4) + " " + fArgs[3].ToString(fmt).PadLeft(4) + " " + fArgs[4].ToString(fmt).PadLeft(4) + " "
                + fArgs[5].ToString(fmt).PadLeft(4) + " " + fArgs[6].ToString(fmt).PadLeft(4) + " " + fArgs[7].ToString(fmt).PadLeft(4) + " "
                + fArgs[8].ToString(fmt).PadLeft(4) + " " + fArgs[9].ToString(fmt).PadLeft(4) + " " + fArgs[10].ToString(fmt).PadLeft(4) + " "
                + fArgs[11].ToString(fmt).PadLeft(4) + " " + fArgs[12].ToString(fmt).PadLeft(4) + " " + fArgs[13].ToString(fmt).PadLeft(4));

            // do in deg
            for (int i = 0; i < 14; i++)
                fArgs[i] = fArgs[i] * 180.0 / Math.PI;

            strbuild.AppendLine("fundarg = " + fArgs[0].ToString(fmt).PadLeft(4) + " " + fArgs[1].ToString(fmt).PadLeft(4) + " "
                + fArgs[2].ToString(fmt).PadLeft(4) + " " + fArgs[3].ToString(fmt).PadLeft(4) + " " + fArgs[4].ToString(fmt).PadLeft(4) + " "
                + fArgs[5].ToString(fmt).PadLeft(4) + " " + fArgs[6].ToString(fmt).PadLeft(4) + " " + fArgs[7].ToString(fmt).PadLeft(4) + " "
                + fArgs[8].ToString(fmt).PadLeft(4) + " " + fArgs[9].ToString(fmt).PadLeft(4) + " " + fArgs[10].ToString(fmt).PadLeft(4) + " "
                + fArgs[11].ToString(fmt).PadLeft(4) + " " + fArgs[12].ToString(fmt).PadLeft(4) + " " + fArgs[13].ToString(fmt).PadLeft(4));


        }
        public void testprecess()
        {
            double ttt, psia, wa, epsa, chia;
            double[,] prec = new double[3, 3];
            AstroLib.EOpt opt = AstroLib.EOpt.e80;

            ttt = 0.042623631889;
            // ttt = 0.04262362174880504;
            prec = AstroLibr.precess(ttt, opt, out psia, out wa, out epsa, out chia);

            strbuild.AppendLine("prec = " + prec[0, 0].ToString(fmt).PadLeft(4) + " " + prec[0, 1].ToString(fmt).PadLeft(4) + " " + prec[0, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("prec = " + prec[1, 0].ToString(fmt).PadLeft(4) + " " + prec[1, 1].ToString(fmt).PadLeft(4) + " " + prec[1, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("prec = " + prec[2, 0].ToString(fmt).PadLeft(4) + " " + prec[2, 1].ToString(fmt).PadLeft(4) + " " + prec[2, 2].ToString(fmt).PadLeft(4));

            prec = AstroLibr.precess(ttt, AstroLib.EOpt.e06eq, out psia, out wa, out epsa, out chia);

            strbuild.AppendLine("prec00 = " + prec[0, 0].ToString(fmt).PadLeft(4) + " " + prec[0, 1].ToString(fmt).PadLeft(4) + " " + prec[0, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("prec00 = " + prec[1, 0].ToString(fmt).PadLeft(4) + " " + prec[1, 1].ToString(fmt).PadLeft(4) + " " + prec[1, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("prec00 = " + prec[2, 0].ToString(fmt).PadLeft(4) + " " + prec[2, 1].ToString(fmt).PadLeft(4) + " " + prec[2, 2].ToString(fmt).PadLeft(4));
        }

        public void testnutation()
        {
            double[] fArgs = new double[14];
            double ttt, ddpsi, ddeps;
            double[,] nut = new double[3, 3];
            double[,] nut00 = new double[3, 3];
            AstroLib.EOpt opt = AstroLib.EOpt.e80;
            string fileLoc;
            fileLoc = @"D:\Codes\LIBRARY\DataLib\nut80.dat";
            EOPSPWLibr.iau80in(fileLoc, out EOPSPWLibr.iau80arr);
            fileLoc = @"D:\Codes\LIBRARY\DataLib\";
            EOPSPWLibr.iau06in(fileLoc, out EOPSPWLibr.iau06arr);

            double deltapsi, deltaeps, trueeps, meaneps;
            ttt = 0.042623631889;
            ddpsi = -0.052195;
            ddeps = -0.003875;

            AstroLibr.fundarg(ttt, opt, out fArgs);

            nut = AstroLibr.nutation(ttt, ddpsi, ddeps, EOPSPWLibr.iau80arr, fArgs, out deltapsi, out deltaeps, out trueeps, out meaneps);

            strbuild.AppendLine("nut = " + nut[0, 0].ToString(fmt).PadLeft(4) + " " + nut[0, 1].ToString(fmt).PadLeft(4) + " " + nut[0, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("nut = " + nut[1, 0].ToString(fmt).PadLeft(4) + " " + nut[1, 1].ToString(fmt).PadLeft(4) + " " + nut[1, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("nut = " + nut[2, 0].ToString(fmt).PadLeft(4) + " " + nut[2, 1].ToString(fmt).PadLeft(4) + " " + nut[2, 2].ToString(fmt).PadLeft(4));

            AstroLibr.fundarg(ttt, AstroLib.EOpt.e06cio, out fArgs);
            nut00 = AstroLibr.precnutbias00a(ttt, ddpsi, ddeps, EOPSPWLibr.iau06arr, AstroLib.EOpt.e06cio, fArgs);
            strbuild.AppendLine("nut06 c= " + nut[0, 0].ToString(fmt).PadLeft(4) + " " + nut[0, 1].ToString(fmt).PadLeft(4) + " " + nut[0, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("nut06 c= " + nut[1, 0].ToString(fmt).PadLeft(4) + " " + nut[1, 1].ToString(fmt).PadLeft(4) + " " + nut[1, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("nut06 c= " + nut[2, 0].ToString(fmt).PadLeft(4) + " " + nut[2, 1].ToString(fmt).PadLeft(4) + " " + nut[2, 2].ToString(fmt).PadLeft(4));

            //AstroLibr.fundarg(ttt, AstroLib.EOpt.e00a, out fArgs);
            //nut00 = AstroLibr.nutation00a(ttt, ddpsi, ddeps, EOPSPWLibr.iau06arr, AstroLib.EOpt.e00a);
            //strbuild.AppendLine("nut06 a= " + nut[0, 0].ToString(fmt).PadLeft(4) + " " + nut[0, 1].ToString(fmt).PadLeft(4) + " " + nut[0, 2].ToString(fmt).PadLeft(4));
            //strbuild.AppendLine("nut06 a= " + nut[1, 0].ToString(fmt).PadLeft(4) + " " + nut[1, 1].ToString(fmt).PadLeft(4) + " " + nut[1, 2].ToString(fmt).PadLeft(4));
            //strbuild.AppendLine("nut06 a= " + nut[2, 0].ToString(fmt).PadLeft(4) + " " + nut[2, 1].ToString(fmt).PadLeft(4) + " " + nut[2, 2].ToString(fmt).PadLeft(4));
        }


        public void testnutationqmod()
        {
            double[] fArgs = new double[14];
            double[,] nutq = new double[3, 3];
            double ttt, deltapsi, deltaeps, meaneps;
            AstroLib.EOpt opt = AstroLib.EOpt.e80;
            string fileLoc;
            fileLoc = @"D:\Codes\LIBRARY\DataLib\nut80.dat";
            EOPSPWLibr.iau80in(fileLoc, out EOPSPWLibr.iau80arr);
            ttt = 0.042623631889;
            // ttt = 0.04262362174880504;

            AstroLibr.fundarg(ttt, opt, out fArgs);

            nutq = AstroLibr.nutationqmod(ttt, EOPSPWLibr.iau80arr, fArgs, out deltapsi, out deltaeps, out meaneps);
        }
        public void testsidereal()
        {
            double[] fArgs = new double[14];
            double[,] nut = new double[3, 3];
            double[,] st = new double[3, 3];
            double ttt, jdut1, deltapsi, deltaeps, meaneps, trueeps, ddpsi, ddeps, lod;
            int eqeterms = 2;
            AstroLib.EOpt opt = AstroLib.EOpt.e80;
            EOPSPWLib.iau80Class iau80arr;
            EOPSPWLib.iau06Class iau06arr;
            jdut1 = 2453101.82740678310;
            ttt = 0.042623631889;
            //ttt = 0.04262362174880504;
            lod = 0.001556;
            string fileLoc;
            fileLoc = @"D:\Codes\LIBRARY\DataLib\nut80.dat";
            EOPSPWLibr.iau80in(fileLoc, out iau80arr);
            fileLoc = @"D:\Codes\LIBRARY\DataLib\";
            EOPSPWLibr.iau06in(fileLoc, out iau06arr);

            ddpsi = -0.052195;
            ddeps = -0.003875;

            AstroLibr.fundarg(ttt, opt, out fArgs);
            nut = AstroLibr.nutation(ttt, ddpsi, ddeps, EOPSPWLibr.iau80arr, fArgs, out deltapsi, out deltaeps, out trueeps, out meaneps);
            st = AstroLibr.sidereal(jdut1, deltapsi, meaneps, fArgs, lod, eqeterms, opt);
            strbuild.AppendLine("st = " + st[0, 0].ToString(fmt).PadLeft(4) + " " + st[0, 1].ToString(fmt).PadLeft(4) + " " + st[0, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("st = " + st[1, 0].ToString(fmt).PadLeft(4) + " " + st[1, 1].ToString(fmt).PadLeft(4) + " " + st[1, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("st = " + st[2, 0].ToString(fmt).PadLeft(4) + " " + st[2, 1].ToString(fmt).PadLeft(4) + " " + st[2, 2].ToString(fmt).PadLeft(4));


            st = AstroLibr.sidereal(jdut1, deltapsi, meaneps, fArgs, lod, eqeterms, AstroLib.EOpt.e06eq);
            strbuild.AppendLine("st00 = " + st[0, 0].ToString(fmt).PadLeft(4) + " " + st[0, 1].ToString(fmt).PadLeft(4) + " " + st[0, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("st00 = " + st[1, 0].ToString(fmt).PadLeft(4) + " " + st[1, 1].ToString(fmt).PadLeft(4) + " " + st[1, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("st00 = " + st[2, 0].ToString(fmt).PadLeft(4) + " " + st[2, 1].ToString(fmt).PadLeft(4) + " " + st[2, 2].ToString(fmt).PadLeft(4));
        }
        public void testpolarm()
        {
            double[,] pm = new double[3, 3];
            double xp, yp, ttt;
            AstroLib.EOpt opt = AstroLib.EOpt.e80;

            ttt = 0.042623631889;
            //ttt = 0.04262363188899416;
            xp = 0.0;
            yp = 0.0;
            opt = AstroLib.EOpt.e80;

            pm = AstroLibr.polarm(xp, yp, ttt, opt);

            strbuild.AppendLine("pm = " + pm[0, 0].ToString(fmt).PadLeft(4) + " " + pm[0, 1].ToString(fmt).PadLeft(4) + " " + pm[0, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("pm = " + pm[1, 0].ToString(fmt).PadLeft(4) + " " + pm[1, 1].ToString(fmt).PadLeft(4) + " " + pm[1, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("pm = " + pm[2, 0].ToString(fmt).PadLeft(4) + " " + pm[2, 1].ToString(fmt).PadLeft(4) + " " + pm[2, 2].ToString(fmt).PadLeft(4));

            pm = AstroLibr.polarm(xp, yp, ttt, AstroLib.EOpt.e06eq);

            strbuild.AppendLine("pm06 = " + pm[0, 0].ToString(fmt).PadLeft(4) + " " + pm[0, 1].ToString(fmt).PadLeft(4) + " " + pm[0, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("pm06 = " + pm[1, 0].ToString(fmt).PadLeft(4) + " " + pm[1, 1].ToString(fmt).PadLeft(4) + " " + pm[1, 2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("pm06 = " + pm[2, 0].ToString(fmt).PadLeft(4) + " " + pm[2, 1].ToString(fmt).PadLeft(4) + " " + pm[2, 2].ToString(fmt).PadLeft(4));
        }
        public void testgstime()
        {
            double gst, jdut1;
            jdut1 = 2453101.82740678310;

            gst = AstroLibr.gstime(jdut1);

            strbuild.AppendLine("gst = " + gst.ToString(fmt).PadLeft(4) + " " + (gst * 180.0 / Math.PI).ToString(fmt).PadLeft(4));
        }

        public void testlstime()
        {
            double rad = 180.0 / Math.PI;
            double lon, jdut1, lst, gst;
            lon = -104.0 / rad;
            jdut1 = 2453101.82740678310;

            AstroLibr.lstime(lon, jdut1, out lst, out gst);

            strbuild.AppendLine("lst = " + lst.ToString(fmt).PadLeft(4) + " " + (lst * 180.0 / Math.PI).ToString(fmt).PadLeft(4));
        }

        public void testhms_sec()
        {
            int hr, min;
            double second;
            double utsec;
            hr = 12;
            min = 34;
            second = 56.233;
            utsec = 0.0;
            strbuild.AppendLine("hr minute second " + hr + " " + min + " " + second);

            MathTimeLibr.hms_sec(ref hr, ref min, ref second, MathTimeLib.Edirection.eto, ref utsec);
            strbuild.AppendLine("utsec = " + utsec.ToString(fmt).PadLeft(4));

            MathTimeLibr.hms_sec(ref hr, ref min, ref second, MathTimeLib.Edirection.efrom, ref utsec);
            strbuild.AppendLine("hr minute second " + hr + " " + min + " " + second);
        }

        public void testhms_ut()
        {
            int hr, min;
            double second;
            double ut;
            hr = 13;
            min = 22;
            second = 45.98;
            ut = 0.0;
            strbuild.AppendLine("hr minute second " + hr + " " + min + " " + second);

            MathTimeLibr.hms_ut(ref hr, ref min, ref second, MathTimeLib.Edirection.eto, ref ut);
            strbuild.AppendLine("ut = " + ut.ToString(fmt).PadLeft(4));

            MathTimeLibr.hms_ut(ref hr, ref min, ref second, MathTimeLib.Edirection.efrom, ref ut);
            strbuild.AppendLine("hr minute second " + hr + " " + min + " " + second);
        }

        public void testhms_rad()
        {
            int hr, min;
            double second;
            double hms;
            hr = 15;
            min = 15;
            second = 53.63;
            hms = 0.0;
            strbuild.AppendLine("hr minute second " + hr + " " + min + " " + second);

            MathTimeLibr.hms_rad(ref hr, ref min, ref second, MathTimeLib.Edirection.eto, ref hms);
            strbuild.AppendLine("hms = " + hms.ToString(fmt).PadLeft(4));

            MathTimeLibr.hms_rad(ref hr, ref min, ref second, MathTimeLib.Edirection.efrom, ref hms);
            strbuild.AppendLine("hr minute second " + hr + " " + min + " " + second);
        }

        public void testdms_rad()
        {
            int deg, min;
            double second;
            double dms;
            deg = -35;
            min = -15;
            second = -53.63;
            dms = 0.0;
            strbuild.AppendLine("deg minute second " + deg + " " + min + " " + second);

            MathTimeLibr.dms_rad(ref deg, ref min, ref second, MathTimeLib.Edirection.eto, ref dms);
            strbuild.AppendLine("dms = " + dms.ToString(fmt).PadLeft(4));

            MathTimeLibr.dms_rad(ref deg, ref min, ref second, MathTimeLib.Edirection.efrom, ref dms);
            strbuild.AppendLine("deg minute second " + deg + " " + min + " " + second);
        }


        public void testeci_ecef()
        {
            double[] fArgs = new double[14];
            double[] fArgs06 = new double[14];

            double[] reci = new double[3];
            double[] veci = new double[3];
            double[] recii = new double[3];
            double[] vecii = new double[3];
            double[] recef = new double[3];
            double[] vecef = new double[3];
            double conv;
            Int32 year, mon, day, hr, minute, dat;
            double jd, jdFrac, jdut1, second, dut1, ttt, lod, xp, yp, ddx, ddy, ddpsi, ddeps,
                jdtt, jdftt;
            double x, y, s;

            conv = Math.PI / (180.0 * 3600.0);  // arcsec to rad

            recef = new double[] { -1033.4793830, 7901.2952754, 6380.3565958 };
            vecef = new double[] { -3.225636520, -2.872451450, 5.531924446 };
            year = 2004;
            mon = 4;
            day = 6;
            hr = 7;
            minute = 51;
            second = 28.386009;

            // sofa example -------------
            //year = 2007;
            //mon = 4;
            //day = 5;
            //hr = 12;
            //minute = 0;
            //second = 0.0;

            MathTimeLibr.jday(year, mon, day, hr, minute, second, out jd, out jdFrac);

            dut1 = -0.4399619;      // second
            dat = 32;               // second
            xp = -0.140682 * conv;  // " to rad
            yp = 0.333309 * conv;
            lod = 0.0015563;
            ddpsi = -0.052195 * conv;  // " to rad
            ddeps = -0.003875 * conv;
            ddx = -0.000205 * conv;    // " to rad
            ddy = -0.000136 * conv;

            // sofa example -------------
            //dut1 = -0.072073685;      // second
            //dat = 33;                // second
            //xp = 0.0349282 * conv;  // " to rad
            //yp = 0.4833163 * conv;
            //lod = 0.0;
            //ddpsi = -0.0550655 * conv;  // " to rad
            //ddeps = -0.006358 * conv;
            //ddx = 0.000175 * conv;    // " to rad
            //ddy = -0.0002259 * conv;

            jdtt = jd;
            jdftt = jdFrac + (dat + 32.184) / 86400.0;
            ttt = (jdtt + jdftt - 2451545.0) / 36525.0;
            strbuild.AppendLine("ttt wo base (use this) " + ttt.ToString());
            jdut1 = jd + jdFrac + dut1 / 86400.0;

            strbuild.AppendLine("ITRF          IAU-76/FK5   " + recef[0].ToString(fmt).PadLeft(4) + " " + recef[1].ToString(fmt).PadLeft(4) + " " + recef[2].ToString(fmt).PadLeft(4) + " "
                + vecef[0].ToString(fmt).PadLeft(4) + " " + vecef[1].ToString(fmt).PadLeft(4) + " " + vecef[2].ToString(fmt).PadLeft(4));

            string fileLoc;
            fileLoc = @"D:\Codes\LIBRARY\DataLib\nut80.dat";
            EOPSPWLibr.iau80in(fileLoc, out EOPSPWLibr.iau80arr);

            fileLoc = @"D:\Codes\LIBRARY\DataLib\";
            EOPSPWLibr.iau06in(fileLoc, out EOPSPWLibr.iau06arr);

            // test creating xys file
            // done, works. :-)
            //fileLoc = @"D:\Codes\LIBRARY\DataLib\";
            //AstroLibr.createXYS(fileLoc, iau06arr, fArgs);

            // now read it in
            //AstroLib.xysdataClass[] xysarr = AstroLibr.xysarr;
            fileLoc = @"D:\Codes\LIBRARY\DataLib\";
            AstroLibr.readXYS(ref AstroLibr.xysarr, fileLoc, "xysdata.dat");

            // now test it for interpolation
            //jdtt = jd + jdFrac + (dat + 32.184) / 86400.0;
            AstroLibr.fundarg(ttt, AstroLib.EOpt.e06cio, out fArgs06);
            AstroLibr.iau06xysS(ttt, EOPSPWLibr.iau06arr, fArgs06, out x, out y, out s);
            strbuild.AppendLine("iau06xysS wo   x   " + x.ToString() + " y " + y.ToString() + " s " + s.ToString());
            strbuild.AppendLine("iau06xysS wo   x   " + (x / conv).ToString() + " y " + (y / conv).ToString() + " s " + (s / conv).ToString() + " arcsec");
            x = x + ddx;
            y = y + ddy;
            strbuild.AppendLine("iau06xysS      x   " + x.ToString() + " y " + y.ToString() + " s " + s.ToString());
            strbuild.AppendLine("iau06xysS      x   " + (x / conv).ToString() + " y " + (y / conv).ToString() + " s " + (s / conv).ToString() + " arcsec");
            AstroLibr.findxysparam(jdtt + jdftt, 0.0, 'n', AstroLibr.xysarr, out x, out y, out s);
            strbuild.AppendLine("findxysparam n x " + x.ToString() + "   y " + y.ToString() + "   s " + s.ToString());
            AstroLibr.findxysparam(jdtt + jdftt, 0.0, 'l', AstroLibr.xysarr, out x, out y, out s);
            strbuild.AppendLine("findxysparam l x " + x.ToString() + " y " + y.ToString() + " s " + s.ToString());
            AstroLibr.findxysparam(jdtt + jdftt, 0.0, 's', AstroLibr.xysarr, out x, out y, out s);
            strbuild.AppendLine("findxysparam s x " + x.ToString() + " y " + y.ToString() + " s " + s.ToString());

            // get more accurate vaules for now
            AstroLibr.iau06xysS(ttt, EOPSPWLibr.iau06arr, fArgs, out x, out y, out s);

            AstroLibr.eci_ecef(ref recii, ref vecii, MathTimeLib.Edirection.efrom, ref recef, ref vecef,
                 EOPSPWLibr.iau80arr, jdtt, jdftt, jdut1, lod, xp, yp, ddpsi, ddeps);
            strbuild.AppendLine("GCRF          IAU-76/FK5   " + recii[0].ToString(fmt).PadLeft(4) + " " + recii[1].ToString(fmt).PadLeft(4) + " " + recii[2].ToString(fmt).PadLeft(4) + " "
                + vecii[0].ToString(fmt).PadLeft(4) + " " + vecii[1].ToString(fmt).PadLeft(4) + " " + vecii[2].ToString(fmt).PadLeft(4));

            AstroLibr.eci_ecef06(ref reci, ref veci, MathTimeLib.Edirection.efrom, ref recef, ref vecef,
                 AstroLib.EOpt.e06cio, EOPSPWLibr.iau06arr, AstroLibr.xysarr, jdtt, jdftt, jdut1, lod, xp, yp, ddx, ddy);
            strbuild.AppendLine("GCRF          IAU-2006 CIO " + reci[0].ToString(fmt).PadLeft(4) + " " + reci[1].ToString(fmt).PadLeft(4) + " " + reci[2].ToString(fmt).PadLeft(4) + " "
                + veci[0].ToString(fmt).PadLeft(4) + " " + veci[1].ToString(fmt).PadLeft(4) + " " + veci[2].ToString(fmt).PadLeft(4));


            // try backwards
            AstroLibr.eci_ecef(ref recii, ref vecii, MathTimeLib.Edirection.eto, ref recef, ref vecef,
                EOPSPWLibr.iau80arr, jdtt, jdftt, jdut1, lod, xp, yp, ddpsi, ddeps);
            strbuild.AppendLine("ITRF rev       IAU-76/FK5   " + recef[0].ToString(fmt).PadLeft(4) + " " + recef[1].ToString(fmt).PadLeft(4) + " " + recef[2].ToString(fmt).PadLeft(4) + " "
                + vecef[0].ToString(fmt).PadLeft(4) + " " + vecef[1].ToString(fmt).PadLeft(4) + " " + vecef[2].ToString(fmt).PadLeft(4));

            AstroLibr.eci_ecef06(ref reci, ref veci, MathTimeLib.Edirection.eto, ref recef, ref vecef,
                 AstroLib.EOpt.e06cio, EOPSPWLibr.iau06arr, AstroLibr.xysarr, jdtt, jdftt, jdut1, lod, xp, yp, ddx, ddy);
            strbuild.AppendLine("ITRF rev       IAU-2006 CIO " + recef[0].ToString(fmt).PadLeft(4) + " " + recef[1].ToString(fmt).PadLeft(4) + " " + recef[2].ToString(fmt).PadLeft(4) + " "
               + vecef[0].ToString(fmt).PadLeft(4) + " " + vecef[1].ToString(fmt).PadLeft(4) + " " + vecef[2].ToString(fmt).PadLeft(4));


            // writeout data for table interpolation
            Int32 i;
            year = 1980;
            mon = 1;
            day = 1;
            hr = 0;
            minute = 0;
            second = 0.0;
            Int32 ktrActObs;
            string EOPupdate;
            char interp = 'x';  // full series


            fileLoc = @"D:\Codes\LIBRARY\DataLib\";
            EOPSPWLibr.iau06in(fileLoc, out EOPSPWLibr.iau06arr);

            // read interpolated one
            //EOPSPWLibr.initEOPArrayP(ref EOPSPWLibr.eopdataP);

            // read existing data - this does not find x, y, s!
            //getCurrEOPFileName(this.EOPSPWLoc.Text, out eopFileName);
            string eopFileName = @"D:\Codes\LIBRARY\DataLib\EOP-All-v1.1_2018-01-04.txt";
            EOPSPWLibr.readeop(ref EOPSPWLibr.eopdata, eopFileName, out ktrActObs, out EOPupdate);

            // now find table of CIO values
            double deltapsi, deltaeps, tempval;

            // rad to "
            double convrt = (180.0 * 3600.0) / Math.PI;
            strbuild.AppendLine("CIO tests      x           y           s          ddpsi           ddeps      ddx       ddy   mjd ");
            for (i = 0; i < 14; i++)   // 14500
            {
                MathTimeLibr.jday(year, mon, day + i, hr, minute, second, out jd, out jdFrac);
                EOPSPWLibr.findeopparam(jd, jdFrac, 's', EOPSPWLibr.eopdata, out dut1, out dat,
                   out lod, out xp, out yp, out ddpsi, out ddeps, out ddx, out ddy);
                jdtt = jd;
                jdftt = jdFrac + (dat + 32.184) / 86400.0;
                ttt = (jdtt + jdftt - 2451545.0) / 36525.0;

                AstroLibr.fundarg(ttt, AstroLib.EOpt.e80, out fArgs);

                //ddpsi = 0.0;
                //ddeps = 0.0;
                deltapsi = 0.0;
                deltaeps = 0.0;
                int ii;
                for (ii = 105; ii >= 0; ii--)
                {
                    tempval = EOPSPWLibr.iau80arr.iar80[ii, 0] * fArgs[0] + EOPSPWLibr.iau80arr.iar80[ii, 1] * fArgs[1] + EOPSPWLibr.iau80arr.iar80[ii, 2] * fArgs[2] +
                             EOPSPWLibr.iau80arr.iar80[ii, 3] * fArgs[3] + EOPSPWLibr.iau80arr.iar80[ii, 4] * fArgs[4];
                    deltapsi = deltapsi + (EOPSPWLibr.iau80arr.rar80[ii, 0] + EOPSPWLibr.iau80arr.rar80[ii, 1] * ttt) * Math.Sin(tempval);
                    deltaeps = deltaeps + (EOPSPWLibr.iau80arr.rar80[ii, 2] + EOPSPWLibr.iau80arr.rar80[ii, 3] * ttt) * Math.Cos(tempval);
                }

                // --------------- find nutation parameters --------------------
                deltapsi = (deltapsi + ddpsi) % (2.0 * Math.PI);
                deltaeps = (deltaeps + ddeps) % (2.0 * Math.PI);
                deltapsi = deltapsi * convrt;
                deltaeps = deltaeps * convrt;

                // CIO parameters
                AstroLibr.fundarg(ttt, AstroLib.EOpt.e06cio, out fArgs06);
                //ddx = 0.0;
                //ddy = 0.0;
                AstroLibr.iau06xys(jdtt, jdftt, ddx, ddy, interp, EOPSPWLibr.iau06arr, fArgs06, AstroLibr.xysarr, out x, out y, out s);
                x = x * convrt;
                y = y * convrt;
                s = s * convrt;
        
                strbuild.AppendLine(" " + x.ToString(fmt).PadLeft(4) + " " + y.ToString(fmt).PadLeft(4) + " " + s.ToString(fmt).PadLeft(4) + " " 
                    + deltapsi.ToString(fmt).PadLeft(4) + " " + deltaeps.ToString(fmt).PadLeft(4) + " " 
                    + ddx.ToString(fmt).PadLeft(4) + " " + ddy.ToString(fmt).PadLeft(4) + " "
                    + (jd + jdFrac - 2400000.5).ToString().PadLeft(4));
            }  // for loop
        }

        public void testtod2ecef()
        {
            double[] fArgs = new double[14];
            double[] reci = new double[3];
            double[] veci = new double[3];
            double[] recii = new double[3];
            double[] vecii = new double[3];
            double[] rmod = new double[3];
            double[] vmod = new double[3];
            double[] rtod = new double[3];
            double[] vtod = new double[3];
            double[] rtirs = new double[3];
            double[] vtirs = new double[3];
            double[] rcirs = new double[3];
            double[] vcirs = new double[3];
            double[] rpef = new double[3];
            double[] vpef = new double[3];
            double[] recef = new double[3];
            double[] vecef = new double[3];
            double[] recefi = new double[3];
            double[] vecefi = new double[3];
            double[] rtemp = new double[3];
            double[] vtemp = new double[3];
            double conv;
            Int32 year, mon, day, hr, minute, dat;
            double jd, jdFrac, jdut1, second, dut1, ttt, lod, xp, yp, ddx, ddy, ddpsi, ddeps,
                jdtt, jdftt;

            conv = Math.PI / (180.0 * 3600.0);

            recef = new double[] { -1033.4793830, 7901.2952754, 6380.3565958 };
            vecef = new double[] { -3.225636520, -2.872451450, 5.531924446 };
            year = 2004;
            mon = 4;
            day = 6;
            hr = 7;
            minute = 51;
            second = 28.386009;
            MathTimeLibr.jday(year, mon, day, hr, minute, second, out jd, out jdFrac);

            dut1 = -0.4399619;      // second
            dat = 32;               // second
            xp = -0.140682 * conv;  // " to rad
            yp = 0.333309 * conv;
            lod = 0.0015563;
            ddpsi = -0.052195 * conv;  // " to rad
            ddeps = -0.003875 * conv;
            ddx = -0.000205 * conv;    // " to rad
            ddy = -0.000136 * conv;

            reci = new double[] { 0.0, 0.0, 0.0 };
            veci = new double[] { 0.0, 0.0, 0.0 };

            string fileLoc;
            // can do it either way... with or without  EOPSPWLibr.
            EOPSPWLib.iau80Class iau80arr;
            EOPSPWLib.iau06Class iau06arr;
            fileLoc = @"D:\Codes\LIBRARY\DataLib\nut80.dat";
            EOPSPWLibr.iau80in(fileLoc, out iau80arr);
            fileLoc = @"D:\Codes\LIBRARY\DataLib\";
            EOPSPWLibr.iau06in(fileLoc, out iau06arr);
            jdtt = jd;
            jdftt = jdFrac + (dat + 32.184) / 86400.0;
            ttt = (jdtt + jdftt - 2451545.0) / 36525.0;

            // now read it in
            AstroLib.xysdataClass[] xysarr = AstroLibr.xysarr;
            fileLoc = @"D:\Codes\LIBRARY\DataLib\";
            AstroLibr.readXYS(ref xysarr, fileLoc, "xysdata.dat");

            jdut1 = jd + jdFrac + dut1 / 86400.0;

            strbuild.AppendLine("ITRF start    IAU-76/FK5   " + recef[0].ToString(fmt).PadLeft(4) + " " + recef[1].ToString(fmt).PadLeft(4) + " " + recef[2].ToString(fmt).PadLeft(4) + " "
                + vecef[0].ToString(fmt).PadLeft(4) + " " + vecef[1].ToString(fmt).PadLeft(4) + " " + vecef[2].ToString(fmt).PadLeft(4));

            // PEF
            AstroLibr.ecef_pef(ref recef, ref vecef, MathTimeLib.Edirection.eto, ref rpef, ref vpef, ttt, xp, yp);
            strbuild.AppendLine("PEF           IAU-76/FK5   " + rpef[0].ToString(fmt).PadLeft(4) + " " +
                rpef[1].ToString(fmt).PadLeft(4) + " " + rpef[2].ToString(fmt).PadLeft(4) + " "
                + vpef[0].ToString(fmt).PadLeft(4) + " " + vpef[1].ToString(fmt).PadLeft(4) + " " + vpef[2].ToString(fmt).PadLeft(4));
            AstroLibr.ecef_pef(ref recii, ref vecii, MathTimeLib.Edirection.efrom, ref rpef, ref vpef, ttt, xp, yp);
            strbuild.AppendLine("ITRF  rev     IAU-76/FK5   " + recii[0].ToString(fmt).PadLeft(4) + " " + recii[1].ToString(fmt).PadLeft(4) + " " + recii[2].ToString(fmt).PadLeft(4) + " "
                + vecii[0].ToString(fmt).PadLeft(4) + " " + vecii[1].ToString(fmt).PadLeft(4) + " " + vecii[2].ToString(fmt).PadLeft(4));

            AstroLibr.ecef_tirs(ref recef, ref vecef, MathTimeLib.Edirection.eto, ref rtemp, ref vtemp, ttt, xp, yp);
            strbuild.AppendLine("TIRS          IAU-2006 CIO " + rtemp[0].ToString(fmt).PadLeft(4) + " " +
                rtemp[1].ToString(fmt).PadLeft(4) + " " + rtemp[2].ToString(fmt).PadLeft(4) + " "
                + vtemp[0].ToString(fmt).PadLeft(4) + " " + vtemp[1].ToString(fmt).PadLeft(4) + " " + vtemp[2].ToString(fmt).PadLeft(4));
            AstroLibr.ecef_tirs(ref recefi, ref vecefi, MathTimeLib.Edirection.efrom, ref rtemp, ref vtemp, ttt, xp, yp);
            strbuild.AppendLine("ITRF rev      IAU-2006 CIO " + recefi[0].ToString(fmt).PadLeft(4) + " " + recefi[1].ToString(fmt).PadLeft(4) + " " + recefi[2].ToString(fmt).PadLeft(4) + " "
                + vecefi[0].ToString(fmt).PadLeft(4) + " " + vecefi[1].ToString(fmt).PadLeft(4) + " " + vecefi[2].ToString(fmt).PadLeft(4));

            // TOD
            AstroLibr.ecef_tod(ref recef, ref vecef, MathTimeLib.Edirection.eto, ref rtod, ref vtod, iau80arr, ttt, jdut1, lod, xp, yp, 0.0, 0.0);
            strbuild.AppendLine("TOD wo corr   IAU-76/FK5   " + rtod[0].ToString(fmt).PadLeft(4) + " " + rtod[1].ToString(fmt).PadLeft(4) + " " + rtod[2].ToString(fmt).PadLeft(4) + " "
                + vtod[0].ToString(fmt).PadLeft(4) + " " + vtod[1].ToString(fmt).PadLeft(4) + " " + vtod[2].ToString(fmt).PadLeft(4));
            AstroLibr.ecef_tod(ref recef, ref vecef, MathTimeLib.Edirection.eto, ref rtod, ref vtod, iau80arr, ttt, jdut1, lod, xp, yp, ddpsi, ddeps);
            strbuild.AppendLine("TOD w corr    IAU-76/FK5   " + rtod[0].ToString(fmt).PadLeft(4) + " " + rtod[1].ToString(fmt).PadLeft(4) + " " + rtod[2].ToString(fmt).PadLeft(4) + " "
                + vtod[0].ToString(fmt).PadLeft(4) + " " + vtod[1].ToString(fmt).PadLeft(4) + " " + vtod[2].ToString(fmt).PadLeft(4));
            AstroLibr.ecef_tod(ref recefi, ref vecefi, MathTimeLib.Edirection.efrom, ref rtod, ref vtod, iau80arr, ttt, jdut1, lod, xp, yp, ddpsi, ddeps);
            strbuild.AppendLine("ITRFi         IAU-76/FK5   " + recefi[0].ToString(fmt).PadLeft(4) + " " + recefi[1].ToString(fmt).PadLeft(4) + " " + recefi[2].ToString(fmt).PadLeft(4) + " "
                + vecefi[0].ToString(fmt).PadLeft(4) + " " + vecefi[1].ToString(fmt).PadLeft(4) + " " + vecefi[2].ToString(fmt).PadLeft(4));

            AstroLibr.ecef_cirs(ref recef, ref vecef, MathTimeLib.Edirection.eto, ref rtemp, ref vtemp,
                AstroLib.EOpt.e06cio, iau06arr, ttt, jdut1, lod, xp, yp, ddx, ddy);
            strbuild.AppendLine("CIRS          IAU-2006 CIO " + rtemp[0].ToString(fmt).PadLeft(4) + " " +
                rtemp[1].ToString(fmt).PadLeft(4) + " " + rtemp[2].ToString(fmt).PadLeft(4) + " "
                + vtemp[0].ToString(fmt).PadLeft(4) + " " + vtemp[1].ToString(fmt).PadLeft(4) + " " + vtemp[2].ToString(fmt).PadLeft(4));
            AstroLibr.ecef_cirs(ref recefi, ref vecefi, MathTimeLib.Edirection.efrom, ref rtemp, ref vtemp,
               AstroLib.EOpt.e06cio, iau06arr, ttt, jdut1, lod, xp, yp, ddx, ddy);
            strbuild.AppendLine("ITRF rev      IAU-2006 CIO " + recefi[0].ToString(fmt).PadLeft(4) + " " + recefi[1].ToString(fmt).PadLeft(4) + " " + recefi[2].ToString(fmt).PadLeft(4) + " "
                + vecefi[0].ToString(fmt).PadLeft(4) + " " + vecefi[1].ToString(fmt).PadLeft(4) + " " + vecefi[2].ToString(fmt).PadLeft(4));

            // MOD
            AstroLibr.ecef_mod(ref recef, ref vecef, MathTimeLib.Edirection.eto, ref rmod, ref vmod,
                 iau80arr, jdtt, jdftt, jdut1, lod, xp, yp, 0.0, 0.0);
            strbuild.AppendLine("MOD wo corr   IAU-76/FK5   " + rmod[0].ToString(fmt).PadLeft(4) + " " +
                rmod[1].ToString(fmt).PadLeft(4) + " " + rmod[2].ToString(fmt).PadLeft(4) + " "
                + vmod[0].ToString(fmt).PadLeft(4) + " " + vmod[1].ToString(fmt).PadLeft(4) + " " + vmod[2].ToString(fmt).PadLeft(4));
            AstroLibr.ecef_mod(ref recef, ref vecef, MathTimeLib.Edirection.eto, ref rmod, ref vmod,
               iau80arr, jdtt, jdftt, jdut1, lod, xp, yp, ddpsi, ddeps);
            strbuild.AppendLine("MOD  w corr   IAU-76/FK5   " + rmod[0].ToString(fmt).PadLeft(4) + " " +
                rmod[1].ToString(fmt).PadLeft(4) + " " + rmod[2].ToString(fmt).PadLeft(4) + " "
                + vmod[0].ToString(fmt).PadLeft(4) + " " + vmod[1].ToString(fmt).PadLeft(4) + " " + vmod[2].ToString(fmt).PadLeft(4));
            AstroLibr.ecef_mod(ref recefi, ref vecefi, MathTimeLib.Edirection.efrom, ref rmod, ref vmod,
               iau80arr, jdtt, jdftt, jdut1, lod, xp, yp, ddpsi, ddeps);
            strbuild.AppendLine("ITRF  rev     IAU-76/FK5   " + recefi[0].ToString(fmt).PadLeft(4) + " " + recefi[1].ToString(fmt).PadLeft(4) + " " + recefi[2].ToString(fmt).PadLeft(4) + " "
                + vecefi[0].ToString(fmt).PadLeft(4) + " " + vecefi[1].ToString(fmt).PadLeft(4) + " " + vecefi[2].ToString(fmt).PadLeft(4));


            // J2000
            AstroLibr.eci_ecef(ref recii, ref vecii, MathTimeLib.Edirection.efrom, ref recef, ref vecef,
                 iau80arr, jdtt, jdftt, jdut1, lod, xp, yp, 0.0, 0.0);
            strbuild.AppendLine("J2000 wo corr IAU-76/FK5   " + recii[0].ToString(fmt).PadLeft(4) + " " + recii[1].ToString(fmt).PadLeft(4) + " " + recii[2].ToString(fmt).PadLeft(4) + " "
                + vecii[0].ToString(fmt).PadLeft(4) + " " + vecii[1].ToString(fmt).PadLeft(4) + " " + vecii[2].ToString(fmt).PadLeft(4));

            // GCRF
            AstroLibr.eci_ecef(ref reci, ref veci, MathTimeLib.Edirection.efrom, ref recef, ref vecef,
                  iau80arr, jdtt, jdftt, jdut1, lod, xp, yp, ddpsi, ddeps);
            strbuild.AppendLine("GCRF w corr   IAU-76/FK5   " + reci[0].ToString(fmt).PadLeft(4) + " " + reci[1].ToString(fmt).PadLeft(4) + " " + reci[2].ToString(fmt).PadLeft(4) + " "
                + veci[0].ToString(fmt).PadLeft(4) + " " + veci[1].ToString(fmt).PadLeft(4) + " " + veci[2].ToString(fmt).PadLeft(4));
            AstroLibr.eci_ecef(ref reci, ref veci, MathTimeLib.Edirection.eto, ref recefi, ref vecefi,
                  iau80arr, jdtt, jdftt, jdut1, lod, xp, yp, ddpsi, ddeps);
            strbuild.AppendLine("ITRF rev      IAU-76/FK5   " + recefi[0].ToString(fmt).PadLeft(4) + " " + recefi[1].ToString(fmt).PadLeft(4) + " " + recefi[2].ToString(fmt).PadLeft(4) + " "
                + vecef[0].ToString(fmt).PadLeft(4) + " " + vecefi[1].ToString(fmt).PadLeft(4) + " " + vecefi[2].ToString(fmt).PadLeft(4));

            AstroLibr.eci_ecef06(ref recii, ref vecii, MathTimeLib.Edirection.efrom, ref recef, ref vecef,
                 AstroLib.EOpt.e06cio, iau06arr, AstroLibr.xysarr, jdtt, jdftt, jdut1, lod, xp, yp, ddx, ddy);
            strbuild.AppendLine("GCRF          IAU-2006 CIO " + recii[0].ToString(fmt).PadLeft(4) + " " + recii[1].ToString(fmt).PadLeft(4) + " " + recii[2].ToString(fmt).PadLeft(4) + " "
                + vecii[0].ToString(fmt).PadLeft(4) + " " + vecii[1].ToString(fmt).PadLeft(4) + " " + vecii[2].ToString(fmt).PadLeft(4));
            AstroLibr.eci_ecef06(ref recii, ref vecii, MathTimeLib.Edirection.eto, ref recefi, ref vecefi,
                 AstroLib.EOpt.e06cio, iau06arr, AstroLibr.xysarr, jdtt, jdftt, jdut1, lod, xp, yp, ddx, ddy);
            strbuild.AppendLine("ITRF rev      IAU-2006 CIO " + recefi[0].ToString(fmt).PadLeft(4) + " " + recefi[1].ToString(fmt).PadLeft(4) + " " + recefi[2].ToString(fmt).PadLeft(4) + " "
                + vecefi[0].ToString(fmt).PadLeft(4) + " " + vecefi[1].ToString(fmt).PadLeft(4) + " " + vecefi[2].ToString(fmt).PadLeft(4));

            // sofa
            strbuild.AppendLine("SOFA ECI CIO  5102.508959486507   6123.011392959787   6378.136934384333");
            strbuild.AppendLine("SOFA ECI 06a  5102.508965811828   6123.011397147659   6378.136925303720");
            strbuild.AppendLine("SOFA ECI 00a  5102.508965732931   6123.011397847143   6378.136924695331");


            // now reverses from eci
            strbuild.AppendLine("GCRF wco STARTIAU-76/FK5   " + reci[0].ToString(fmt).PadLeft(4) + " " + reci[1].ToString(fmt).PadLeft(4) + " " + reci[2].ToString(fmt).PadLeft(4) + " "
                + veci[0].ToString(fmt).PadLeft(4) + " " + veci[1].ToString(fmt).PadLeft(4) + " " + veci[2].ToString(fmt).PadLeft(4));

            // PEF
            AstroLibr.eci_pef(ref reci, ref veci, MathTimeLib.Edirection.eto, ref rpef, ref vpef,
                 iau80arr, jdtt, jdftt, jdut1, lod, ddpsi, ddeps);
            strbuild.AppendLine("PEF           IAU-76/FK5   " + rpef[0].ToString(fmt).PadLeft(4) + " " +
                rpef[1].ToString(fmt).PadLeft(4) + " " + rpef[2].ToString(fmt).PadLeft(4) + " "
                + vpef[0].ToString(fmt).PadLeft(4) + " " + vpef[1].ToString(fmt).PadLeft(4) + " " + vpef[2].ToString(fmt).PadLeft(4));
            AstroLibr.eci_pef(ref reci, ref veci, MathTimeLib.Edirection.efrom, ref rpef, ref vpef,
                 iau80arr, jdtt, jdftt, jdut1, lod, ddpsi, ddeps);
            strbuild.AppendLine("ECI rev       IAU-76/FK5   " + reci[0].ToString(fmt).PadLeft(4) + " " + reci[1].ToString(fmt).PadLeft(4) + " " + reci[2].ToString(fmt).PadLeft(4) + " "
                + veci[0].ToString(fmt).PadLeft(4) + " " + veci[1].ToString(fmt).PadLeft(4) + " " + veci[2].ToString(fmt).PadLeft(4));

            AstroLibr.eci_tirs(ref recii, ref vecii, MathTimeLib.Edirection.eto, ref rtirs, ref vtirs,
                AstroLib.EOpt.e06cio, iau06arr, jdtt, jdftt, jdut1, lod, ddx, ddy);
            strbuild.AppendLine("TIRS          IAU-2006 CIO  " + rpef[0].ToString(fmt).PadLeft(4) + " " +
                rtirs[1].ToString(fmt).PadLeft(4) + " " + rtirs[2].ToString(fmt).PadLeft(4) + " "
                + vtirs[0].ToString(fmt).PadLeft(4) + " " + vtirs[1].ToString(fmt).PadLeft(4) + " " + vtirs[2].ToString(fmt).PadLeft(4));
            AstroLibr.eci_tirs(ref recii, ref vecii, MathTimeLib.Edirection.efrom, ref rtirs, ref vtirs,
                AstroLib.EOpt.e06cio, iau06arr, jdtt, jdftt, jdut1, lod, ddx, ddy);
            strbuild.AppendLine("ECI rev       IAU-2006 CIO " + recii[0].ToString(fmt).PadLeft(4) + " " + recii[1].ToString(fmt).PadLeft(4) + " " + recii[2].ToString(fmt).PadLeft(4) + " "
                + vecii[0].ToString(fmt).PadLeft(4) + " " + vecii[1].ToString(fmt).PadLeft(4) + " " + vecii[2].ToString(fmt).PadLeft(4));

            // TOD
            AstroLibr.eci_tod(ref reci, ref veci, MathTimeLib.Edirection.eto, ref rtod, ref vtod,
                 iau80arr, jdtt, jdftt, jdut1, lod, ddpsi, ddeps);
            strbuild.AppendLine("TOD           IAU-76/FK5   " + rtod[0].ToString(fmt).PadLeft(4) + " " + rtod[1].ToString(fmt).PadLeft(4) + " " + rtod[2].ToString(fmt).PadLeft(4) + " "
                + vtod[0].ToString(fmt).PadLeft(4) + " " + vtod[1].ToString(fmt).PadLeft(4) + " " + vtod[2].ToString(fmt).PadLeft(4));
            AstroLibr.eci_tod(ref reci, ref veci, MathTimeLib.Edirection.efrom, ref rtod, ref vtod,
                 iau80arr, jdtt, jdftt, jdut1, lod, ddpsi, ddeps);
            strbuild.AppendLine("ECI rev       IAU-76/FK5   " + reci[0].ToString(fmt).PadLeft(4) + " " + reci[1].ToString(fmt).PadLeft(4) + " " + reci[2].ToString(fmt).PadLeft(4) + " "
                + veci[0].ToString(fmt).PadLeft(4) + " " + veci[1].ToString(fmt).PadLeft(4) + " " + veci[2].ToString(fmt).PadLeft(4));

            AstroLibr.eci_cirs(ref recii, ref vecii, MathTimeLib.Edirection.eto, ref rcirs, ref vcirs,
                 EOpt.e06cio, iau06arr, jdtt, jdftt, ddx, ddy);
            strbuild.AppendLine("CIRS          IAU-2006 CIO " + rcirs[0].ToString(fmt).PadLeft(4) + " " + rcirs[1].ToString(fmt).PadLeft(4) + " " + rcirs[2].ToString(fmt).PadLeft(4) + " "
                + vcirs[0].ToString(fmt).PadLeft(4) + " " + vcirs[1].ToString(fmt).PadLeft(4) + " " + vcirs[2].ToString(fmt).PadLeft(4));
            AstroLibr.eci_cirs(ref recii, ref vecii, MathTimeLib.Edirection.efrom, ref rcirs, ref vcirs,
                EOpt.e06cio, iau06arr, jdtt, jdftt, ddx, ddy);
            strbuild.AppendLine("ECI rev       IAU-2006 CIO " + recii[0].ToString(fmt).PadLeft(4) + " " + recii[1].ToString(fmt).PadLeft(4) + " " + recii[2].ToString(fmt).PadLeft(4) + " "
                + vecii[0].ToString(fmt).PadLeft(4) + " " + vecii[1].ToString(fmt).PadLeft(4) + " " + vecii[2].ToString(fmt).PadLeft(4));

            // MOD
            AstroLibr.eci_mod(ref reci, ref veci, MathTimeLib.Edirection.eto, ref rmod, ref vmod, iau80arr, ttt);
            strbuild.AppendLine("MOD           IAU-76/FK5   " + rmod[0].ToString(fmt).PadLeft(4) + " " +
                rmod[1].ToString(fmt).PadLeft(4) + " " + rmod[2].ToString(fmt).PadLeft(4) + " "
                + vmod[0].ToString(fmt).PadLeft(4) + " " + vmod[1].ToString(fmt).PadLeft(4) + " " + vmod[2].ToString(fmt).PadLeft(4));
            AstroLibr.eci_mod(ref recii, ref vecii, MathTimeLib.Edirection.efrom, ref rmod, ref vmod, iau80arr, ttt);
            strbuild.AppendLine("ECI rev       IAU-76/FK5   " + recii[0].ToString(fmt).PadLeft(4) + " " + recii[1].ToString(fmt).PadLeft(4) + " " + recii[2].ToString(fmt).PadLeft(4) + " "
                + vecii[0].ToString(fmt).PadLeft(4) + " " + vecii[1].ToString(fmt).PadLeft(4) + " " + vecii[2].ToString(fmt).PadLeft(4));
        }
        public void testteme_ecef()
        {
            double[] fArgs = new double[14];
            double[] rteme = new double[3];
            double[] vteme = new double[3];
            int eqeterms = 2;
            double[] recef = new double[3];
            double[] vecef = new double[3];
            double conv;
            Int32 year, mon, day, hr, minute, dat;
            double jd, jdFrac, jdut1, second, dut1, ttt, lod, xp, yp, ddx, ddy, ddpsi, ddeps;

            conv = Math.PI / (180.0 * 3600.0);

            recef = new double[] { -1033.4793830, 7901.2952754, 6380.3565958 };
            vecef = new double[] { -3.225636520, -2.872451450, 5.531924446 };
            year = 2004;
            mon = 4;
            day = 6;
            hr = 7;
            minute = 51;
            second = 28.386009;
            MathTimeLibr.jday(year, mon, day, hr, minute, second, out jd, out jdFrac);

            dut1 = -0.4399619;      // second
            dat = 32;               // second
            xp = -0.140682 * conv;  // " to rad
            yp = 0.333309 * conv;
            lod = 0.0015563;
            ddpsi = -0.052195 * conv;  // " to rad
            ddeps = -0.003875 * conv;
            ddx = -0.000205 * conv;    // " to rad
            ddy = -0.000136 * conv;
            eqeterms = 2;

            // note you have to use tdb for time of ineterst AND j2000 (when dat = 32)
            ttt = (jd + jdFrac + (dat + 32.184) / 86400.0 - 2451545.0 - (32 + 32.184) / 86400.0) / 36525.0;
            jdut1 = jd + jdFrac + dut1 / 86400.0;

            strbuild.AppendLine("ITRF          IAU-76/FK5   " + recef[0].ToString(fmt).PadLeft(4) + " " + recef[1].ToString(fmt).PadLeft(4) + " " + recef[2].ToString(fmt).PadLeft(4) + " "
                + vecef[0].ToString(fmt).PadLeft(4) + " " + vecef[1].ToString(fmt).PadLeft(4) + " " + vecef[2].ToString(fmt).PadLeft(4));

            AstroLibr.teme_ecef(ref rteme, ref vteme, MathTimeLib.Edirection.efrom, ttt, jdut1, lod, xp, yp,
                eqeterms, ref recef, ref vecef);
            strbuild.AppendLine("TEME          IAU-76/FK5   " + rteme[0].ToString(fmt).PadLeft(4) + " " + rteme[1].ToString(fmt).PadLeft(4) + " " + rteme[2].ToString(fmt).PadLeft(4) + " "
                + vteme[0].ToString(fmt).PadLeft(4) + " " + vteme[1].ToString(fmt).PadLeft(4) + " " + vteme[2].ToString(fmt).PadLeft(4));

            AstroLibr.teme_ecef(ref rteme, ref vteme, MathTimeLib.Edirection.eto, ttt, jdut1, lod, xp, yp,
                eqeterms, ref recef, ref vecef);
            strbuild.AppendLine("ITRF          IAU-76/FK5   " + recef[0].ToString(fmt).PadLeft(4) + " " + recef[1].ToString(fmt).PadLeft(4) + " " + recef[2].ToString(fmt).PadLeft(4) + " "
                + vecef[0].ToString(fmt).PadLeft(4) + " " + vecef[1].ToString(fmt).PadLeft(4) + " " + vecef[2].ToString(fmt).PadLeft(4));

        }
        public void testteme_eci()
        {
            double[] fArgs = new double[14];
            double[] rteme = new double[3];
            double[] vteme = new double[3];
            double[] reci = new double[3];
            double[] veci = new double[3];
            double conv;
            Int32 year, mon, day, hr, minute, dat;
            double jd, jdFrac, jdut1, second, dut1, ttt, xp, yp, ddx, ddy, ddpsi, ddeps;
            EOPSPWLib.iau80Class iau80arr;

            conv = Math.PI / (180.0 * 3600.0);

            reci = new double[] { 5102.5089579, 6123.0114007, 6378.1369282 };
            veci = new double[] { -4.743220157, 0.790536497, 5.533755727 };
            year = 2004;
            mon = 4;
            day = 6;
            hr = 7;
            minute = 51;
            second = 28.386009;
            MathTimeLibr.jday(year, mon, day, hr, minute, second, out jd, out jdFrac);

            dut1 = -0.4399619;      // second
            dat = 32;               // second
            xp = -0.140682 * conv;  // " to rad
            yp = 0.333309 * conv;
            ddpsi = -0.052195 * conv;  // " to rad
            ddeps = -0.003875 * conv;
            ddx = -0.000205 * conv;    // " to rad
            ddy = -0.000136 * conv;

            string fileLoc = @"D:\Codes\LIBRARY\DataLib\nut80.dat";
            EOPSPWLibr.iau80in(fileLoc, out iau80arr);

            // note you have to use tdb for time of ineterst AND j2000 (when dat = 32)
            ttt = (jd + jdFrac + (dat + 32.184) / 86400.0 - 2451545.0 - (32 + 32.184) / 86400.0) / 36525.0;
            jdut1 = jd + jdFrac + dut1 / 86400.0;

            strbuild.AppendLine("GCRF          IAU-76/FK5   " + reci[0].ToString(fmt).PadLeft(4) + " " + reci[1].ToString(fmt).PadLeft(4) + " " + reci[2].ToString(fmt).PadLeft(4) + " "
                + veci[0].ToString(fmt).PadLeft(4) + " " + veci[1].ToString(fmt).PadLeft(4) + " " + veci[2].ToString(fmt).PadLeft(4));

            AstroLibr.teme_eci(ref rteme, ref vteme, iau80arr, MathTimeLib.Edirection.efrom, ttt, ddpsi, ddeps, AstroLib.EOpt.e80, ref reci, ref veci);
            strbuild.AppendLine("TEME          IAU-76/FK5   " + rteme[0].ToString(fmt).PadLeft(4) + " " + rteme[1].ToString(fmt).PadLeft(4) + " " + rteme[2].ToString(fmt).PadLeft(4) + " "
                + vteme[0].ToString(fmt).PadLeft(4) + " " + vteme[1].ToString(fmt).PadLeft(4) + " " + vteme[2].ToString(fmt).PadLeft(4));

            AstroLibr.teme_eci(ref rteme, ref vteme, iau80arr, MathTimeLib.Edirection.eto, ttt, ddpsi, ddeps, AstroLib.EOpt.e80, ref reci, ref veci);
            strbuild.AppendLine("GCRF          IAU-76/FK5   " + reci[0].ToString(fmt).PadLeft(4) + " " + reci[1].ToString(fmt).PadLeft(4) + " " + reci[2].ToString(fmt).PadLeft(4) + " "
                + veci[0].ToString(fmt).PadLeft(4) + " " + veci[1].ToString(fmt).PadLeft(4) + " " + veci[2].ToString(fmt).PadLeft(4));
        }
        public void testqmod2ecef()
        {
            double[] fArgs = new double[14];
            double[] rqmod = new double[3];
            double[] vqmod = new double[3];
            double ttt, jdutc;
            AstroLib.EOpt opt = AstroLib.EOpt.e80;
            double[] recef = new double[3];
            double[] vecef = new double[3];
            EOPSPWLib.iau80Class iau80arr;

            ttt = 0.042623631889;
            jdutc = 2453101.82740678310;

            AstroLibr.fundarg(ttt, opt, out fArgs);

            string fileLoc = @"D:\Codes\LIBRARY\DataLib\nut80.dat";
            EOPSPWLibr.iau80in(fileLoc, out iau80arr);

            AstroLibr.qmod2ecef(rqmod, vqmod, ttt, jdutc, iau80arr, opt, out recef, out vecef);
        }
        public void testcsm2efg()
        {
            double[] fArgs = new double[14];
            double[] r1pef = new double[3];
            double[] v1pef = new double[3];
            double[] r1ecef = new double[3];
            double[] v1ecef = new double[3];
            double[] r2ecef = new double[3];
            double[] v2ecef = new double[3];
            double[] r2ric = new double[3];
            double[] v2ric = new double[3];
            double ttt, lod, xp, yp, jdut1, ddpsi, ddeps;
            int eqeterms;
            xp = 0.0;
            yp = 0.0;
            lod = 0.0;
            jdut1 = 2453101.82740678310;
            ttt = 0.042623631889;
            ddpsi = -0.052195;
            ddeps = -0.003875;
            eqeterms = 2;

            //AstroLibr.csm2efg(r1pef, v1pef, r2ric, v2ric, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps, AstroLib.EOpt.e80,
            //    out r1ecef, out v1ecef, out r2ecef, out v2ecef);
        }

        public void testrv_elatlon()
        {
            double rr, ecllat, ecllon, drr, decllat, decllon, rad;
            double[] rijk = new double[3];
            double[] vijk = new double[3];
            rad = 180.0 / Math.PI;
            rr = 12756.00;
            ecllat = 60.04570;
            ecllon = 256.004567345;
            drr = 0.045670;
            decllat = 6.798614;
            decllon = 0.00768;

            AstroLibr.rv_elatlon(ref rijk, ref vijk, MathTimeLib.Edirection.efrom, ref rr, ref ecllat, ref ecllon, ref drr, ref decllat, ref decllon);
            strbuild.AppendLine("rv ecllat " + rijk[0].ToString(fmt).PadLeft(4) + " " + rijk[1].ToString(fmt).PadLeft(4) + " " + rijk[2].ToString(fmt).PadLeft(4) + " " +
                                " " + vijk[0].ToString(fmt).PadLeft(4) + " " + vijk[1].ToString(fmt).PadLeft(4) + " " + vijk[2].ToString(fmt).PadLeft(4));

            AstroLibr.rv_elatlon(ref rijk, ref vijk, MathTimeLib.Edirection.eto, ref rr, ref ecllat, ref ecllon, ref drr, ref decllat, ref decllon);

            strbuild.AppendLine("ecllat  " + rr.ToString(fmt).PadLeft(4) + " " + (ecllat * rad).ToString(fmt).PadLeft(4) + " " + (ecllon * rad).ToString(fmt).PadLeft(4) + " " +
                              " " + drr.ToString(fmt).PadLeft(4) + " " + decllat.ToString(fmt).PadLeft(4) + " " + decllon.ToString(fmt).PadLeft(4));
        }

        public void testrv2radec()
        {
            double[] r = new double[3];
            double[] v = new double[3];
            double rr, rtasc, decl, drr, drtasc, ddecl;
            r = new double[] { -605.79221660, -5870.22951108, 3493.05319896 };
            v = new double[] { -1.56825429, -3.70234891, -6.47948395 };
            rr = rtasc = decl = drr = drtasc = ddecl = 0.0;

            AstroLibr.rv_radec(ref r, ref v, MathTimeLib.Edirection.eto, ref rr, ref rtasc, ref decl, ref drr, ref drtasc, ref ddecl);
            strbuild.AppendLine("rv radec " + r[0].ToString(fmt).PadLeft(4) + " " + r[1].ToString(fmt).PadLeft(4) + " " + r[2].ToString(fmt).PadLeft(4) + " " +
                                " " + v[0].ToString(fmt).PadLeft(4) + " " + v[1].ToString(fmt).PadLeft(4) + " " + v[2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("radec " + rr.ToString(fmt).PadLeft(4) + " " + rtasc.ToString(fmt).PadLeft(4) + " " + decl.ToString(fmt).PadLeft(4) + " " +
                                "  " + drr.ToString(fmt).PadLeft(4) + " " + drtasc.ToString(fmt).PadLeft(4) + " " + ddecl.ToString(fmt).PadLeft(4));

            AstroLibr.rv_radec(ref r, ref v, MathTimeLib.Edirection.efrom, ref rr, ref rtasc, ref decl, ref drr, ref drtasc, ref ddecl);

            strbuild.AppendLine("rv radec  " + r[0].ToString(fmt).PadLeft(4) + " " + r[1].ToString(fmt).PadLeft(4) + " " + r[2].ToString(fmt).PadLeft(4) + " " +
                                "  " + v[0].ToString(fmt).PadLeft(4) + " " + v[1].ToString(fmt).PadLeft(4) + " " + v[2].ToString(fmt).PadLeft(4));
        }
        public void testrv_razel()
        {
            double rad = 180.0 / Math.PI;
            double[] recef = new double[3];
            double[] vecef = new double[3];
            double[] rsecef = new double[3];
            double rho, az, el, drho, daz, del, latgd, lon, alt;
            recef = new double[] { -605.79221660, -5870.22951108, 3493.05319896 };
            vecef = new double[] { -1.56825429, -3.70234891, -6.47948395 };
            //rsecef = new double[] { -1605.79221660, -570.22951108, 193.05319896 };
            lon = -104.883 / rad;
            latgd = 39.883 / rad;
            alt = 2.102;
            rho = 0.0186569;
            az = -0.3501725;
            el = -0.5839385;
            drho = 0.6811410;
            daz = -0.4806057;
            del = 0.6284403;

            AstroLibr.rv_razel(ref recef, ref vecef, latgd, lon, alt, MathTimeLib.Edirection.eto, ref rho, ref az, ref el, ref drho, ref daz, ref del);
            strbuild.AppendLine("rv_razel  " + recef[0].ToString(fmt).PadLeft(4) + " " + recef[1].ToString(fmt).PadLeft(4) + " " + recef[2].ToString(fmt).PadLeft(4) + " " +
                                "  " + vecef[0].ToString(fmt).PadLeft(4) + " " + vecef[1].ToString(fmt).PadLeft(4) + " " + vecef[2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("razel " + rho.ToString(fmt).PadLeft(4) + " " + az.ToString(fmt).PadLeft(4) + " " + el.ToString(fmt).PadLeft(4) + " " +
                                "  " + drho.ToString(fmt).PadLeft(4) + " " + daz.ToString(fmt).PadLeft(4) + " " + del.ToString(fmt).PadLeft(4));

            AstroLibr.rv_razel(ref recef, ref vecef, latgd, lon, alt, MathTimeLib.Edirection.efrom, ref rho, ref az, ref el, ref drho, ref daz, ref del);
            strbuild.AppendLine("rv_razel  " + recef[0].ToString(fmt).PadLeft(4) + " " + recef[1].ToString(fmt).PadLeft(4) + " " + recef[2].ToString(fmt).PadLeft(4) + " " +
                                "  " + vecef[0].ToString(fmt).PadLeft(4) + " " + vecef[1].ToString(fmt).PadLeft(4) + " " + vecef[2].ToString(fmt).PadLeft(4));
        }

        public void testrv_tradec()
        {
            double[] reci = new double[3];
            double[] veci = new double[3];
            double[] rseci = new double[3];
            double rho, trtasc, tdecl, drho, dtrtasc, dtdecl;
            reci = new double[] { 4066.716, -2847.545, 3994.302 };
            veci = new double[] { -1.56825429, -3.70234891, -6.47948395 };
            rseci = new double[] { -1605.79221660, -570.22951108, 193.05319896 };
            rho = 0.2634728;
            trtasc = -0.1492353;
            tdecl = 0.0519525;
            drho = 0.3072265;
            dtrtasc = 0.2045751;
            dtdecl = -0.7510033;

            AstroLibr.rv_tradec(ref reci, ref veci, rseci, MathTimeLib.Edirection.eto, ref rho, ref trtasc, ref tdecl, ref drho, ref dtrtasc,
                ref dtdecl);
            strbuild.AppendLine("rv tradec  " + reci[0].ToString(fmt).PadLeft(4) + " " + reci[1].ToString(fmt).PadLeft(4) + " "
                + reci[2].ToString(fmt).PadLeft(4) + " " +
                "  " + veci[0].ToString(fmt).PadLeft(4) + " " + veci[1].ToString(fmt).PadLeft(4) + " "
                + veci[2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("tradec " + rho.ToString(fmt).PadLeft(4) + " " + trtasc.ToString(fmt).PadLeft(4) + " "
                + tdecl.ToString(fmt).PadLeft(4) + " " +
                "  " + drho.ToString(fmt).PadLeft(4) + " " + dtrtasc.ToString(fmt).PadLeft(4) + " "
                + dtdecl.ToString(fmt).PadLeft(4));

            AstroLibr.rv_tradec(ref reci, ref veci, rseci, MathTimeLib.Edirection.efrom, ref rho, ref trtasc, ref tdecl, ref drho,
                ref dtrtasc, ref dtdecl);
            strbuild.AppendLine("rv tradec  " + reci[0].ToString(fmt).PadLeft(4) + " " + reci[1].ToString(fmt).PadLeft(4) + " "
                + reci[2].ToString(fmt).PadLeft(4) + " " +
                "  " + veci[0].ToString(fmt).PadLeft(4) + " " + veci[1].ToString(fmt).PadLeft(4) + " "
                + veci[2].ToString(fmt).PadLeft(4));
        }
        public void testrvsez_razel()
        {
            double[] rhosez = new double[3];
            double[] drhosez = new double[3];
            double rho, az, el, drho, daz, del;
            rhosez = new double[] { -605.79221660, -5870.22951108, 3493.05319896 };
            drhosez = new double[] { -1.56825429, -3.70234891, -6.47948395 };
            rho = 0.0186569;
            az = -0.3501725;
            el = -0.5839385;
            drho = 0.6811410;
            daz = -0.4806057;
            del = 0.6284403;

            AstroLibr.rvsez_razel(ref rhosez, ref drhosez, MathTimeLib.Edirection.eto, ref rho, ref az, ref el, ref drho, ref daz, ref del);
            strbuild.AppendLine("rv rhosez  " + rhosez[0].ToString(fmt).PadLeft(4) + " " + rhosez[1].ToString(fmt).PadLeft(4) + " " + rhosez[2].ToString(fmt).PadLeft(4) + " " +
                             "  " + drhosez[0].ToString(fmt).PadLeft(4) + " " + drhosez[1].ToString(fmt).PadLeft(4) + " " + drhosez[2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("rhosez " + rho.ToString(fmt).PadLeft(4) + " " + az.ToString(fmt).PadLeft(4) + " " + el.ToString(fmt).PadLeft(4) + " " +
                                "  " + drho.ToString(fmt).PadLeft(4) + " " + daz.ToString(fmt).PadLeft(4) + " " + del.ToString(fmt).PadLeft(4));

            AstroLibr.rvsez_razel(ref rhosez, ref drhosez, MathTimeLib.Edirection.efrom, ref rho, ref az, ref el, ref drho, ref daz, ref del);

            strbuild.AppendLine("rv rhosez  " + rhosez[0].ToString(fmt).PadLeft(4) + " " + rhosez[1].ToString(fmt).PadLeft(4) + " " + rhosez[2].ToString(fmt).PadLeft(4) + " " +
                             "  " + drhosez[0].ToString(fmt).PadLeft(4) + " " + drhosez[1].ToString(fmt).PadLeft(4) + " " + drhosez[2].ToString(fmt).PadLeft(4));
        }
        public void testrv2rsw()
        {
            double[] r = new double[3];
            double[] v = new double[3];
            double[] rrsw = new double[3];
            double[] vrsw = new double[3];
            double[,] tm = new double[3, 3];
            r = new double[] { -605.79221660, -5870.22951108, 3493.05319896 };
            v = new double[] { -1.56825429, -3.70234891, -6.47948395 };

            tm = AstroLibr.rv2rsw(r, v, out rrsw, out vrsw);

            strbuild.AppendLine("rv2rsw " + rrsw[0].ToString(fmt) + " " + rrsw[1].ToString(fmt) + " " + rrsw[2].ToString(fmt) + " " +
                      vrsw[0].ToString(fmt) + " " + vrsw[1].ToString(fmt) + " " + vrsw[2].ToString(fmt));
        }

        public void testrv2pqw()
        {
            double[] r = new double[3];
            double[] v = new double[3];
            double[] rpqw = new double[3];
            double[] vpqw = new double[3];
            double[,] tm = new double[3, 3];
            r = new double[] { -605.79221660, -5870.22951108, 3493.05319896 };
            v = new double[] { -1.56825429, -3.70234891, -6.47948395 };

            AstroLibr.rv2pqw(r, v, out rpqw, out vpqw);

            strbuild.AppendLine("rv2pqw " + rpqw[0].ToString(fmt) + " " + rpqw[1].ToString(fmt) + " " + rpqw[2].ToString(fmt) + " " +
                  vpqw[0].ToString(fmt) + " " + vpqw[1].ToString(fmt) + " " + vpqw[2].ToString(fmt));
        }

        public void testrv2coe()
        {
            Int32 i;
            double[] r = new double[3];
            double[] v = new double[3];
            double[] r1 = new double[3];
            double[] v1 = new double[3];
            double rad = 180.0 / Math.PI;
            double p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper;

            for (i = 1; i <= 21; i++)
            {
                if (i == 1)
                {
                    r = new double[] { -605.79221660, -5870.22951108, 3493.05319896 };
                    v = new double[] { -1.56825429, -3.70234891, -6.47948395 };
                }
                if (i == 2)
                {
                    strbuild.AppendLine(" coe test ----------------------------");
                    r = new double[] { 6524.834, 6862.875, 6448.296 };
                    v = new double[] { 4.901327, 5.533756, -1.976341 };
                }

                // ------- elliptical orbit tests -------------------
                if (i == 3)
                {
                    strbuild.AppendLine(" coe test elliptical ----------------------------");
                    r = new double[] { 1.1372844 * AstroLibr.gravConst.re, -1.0534274 * AstroLibr.gravConst.re, -0.8550194 * AstroLibr.gravConst.re };
                    v = new double[] { 0.6510489 * AstroLibr.gravConst.velkmps, 0.4521008 * AstroLibr.gravConst.velkmps, 0.0381088 * AstroLibr.gravConst.velkmps };
                }
                if (i == 4)
                {
                    strbuild.AppendLine(" coe test elliptical ----------------------------");
                    r = new double[] { 1.056194 * AstroLibr.gravConst.re, -0.8950922 * AstroLibr.gravConst.re, -0.0823703 * AstroLibr.gravConst.re };
                    v = new double[] { -0.5981066 * AstroLibr.gravConst.velkmps, -0.6293575 * AstroLibr.gravConst.velkmps, 0.1468194 * AstroLibr.gravConst.velkmps };
                }

                // ------- circular inclined orbit tests -------------------
                if (i == 5)
                {
                    strbuild.AppendLine(" coe test near circular inclined ----------------------------");
                    r = new double[] { -0.422277 * AstroLibr.gravConst.re, 1.0078857 * AstroLibr.gravConst.re, 0.7041832 * AstroLibr.gravConst.re };
                    v = new double[] { -0.5002738 * AstroLibr.gravConst.velkmps, -0.5415267 * AstroLibr.gravConst.velkmps, 0.4750788 * AstroLibr.gravConst.velkmps };
                }
                if (i == 6)
                {
                    strbuild.AppendLine(" coe test near circular inclined ----------------------------");
                    r = new double[] { -0.7309361 * AstroLibr.gravConst.re, -0.6794646 * AstroLibr.gravConst.re, -0.8331183 * AstroLibr.gravConst.re };
                    v = new double[] { -0.6724131 * AstroLibr.gravConst.velkmps, 0.0341802 * AstroLibr.gravConst.velkmps, 0.5620652 * AstroLibr.gravConst.velkmps };
                }

                if (i == 7) // -- CI u = 45 deg
                {
                    strbuild.AppendLine(" coe test circular inclined ----------------------------");
                    r = new double[] { -2693.34555010128, 6428.43425355863, 4491.37782050409 };
                    v = new double[] { -3.95484712246016, -4.28096585381370, 3.75567104538731 };
                }
                if (i == 8) // -- CI u = 315 deg
                {
                    strbuild.AppendLine(" coe test circular inclined ----------------------------");
                    r = new double[] { -7079.68834483379, 3167.87718823353, -2931.53867301568 };
                    v = new double[] { 1.77608080328182, 6.23770933190509, 2.45134017949138 };
                }

                // ------- elliptical equatorial orbit tests -------------------
                if (i == 9)
                {
                    strbuild.AppendLine(" coe test elliptical near equatorial ----------------------------");
                    r = new double[] { 21648.6109280739, -14058.7723188698, -0.0003598029 };
                    v = new double[] { 2.16378060719980, 3.32694348486311, 0.00000004164788 };
                }
                if (i == 10)
                {
                    strbuild.AppendLine(" coe test elliptical near equatorial ----------------------------");
                    r = new double[] { 7546.9914487222, 24685.1032834356, -0.0003598029 };
                    v = new double[] { 3.79607016047138, -1.15773520476223, 0.00000004164788 };
                }

                if (i == 11) // -- EE w = 20 deg
                {
                    strbuild.AppendLine(" coe test elliptical equatorial ----------------------------");
                    r = new double[] { -22739.1086596208, -22739.1086596208, 0.0 };
                    v = new double[] { 2.48514004188565, -2.02004112073465, 0.0 };
                }
                if (i == 12) // -- EE w = 240 deg
                {
                    strbuild.AppendLine(" coe test elliptical equatorial ----------------------------");
                    r = new double[] { 28242.3662822040, 2470.8868808397, 0.0 };
                    v = new double[] { 0.66575215057746, -3.62533022188304, 0.0 };
                }

                // ------- circular equatorial orbit tests -------------------
                if (i == 13)
                {
                    strbuild.AppendLine(" coe test circular near equatorial ----------------------------");
                    r = new double[] { -2547.3697454933, 14446.8517254604, 0.000 };
                    v = new double[] { -5.13345156333487, -0.90516601477599, 0.00000090977789 };
                }
                if (i == 14)
                {
                    strbuild.AppendLine(" coe test circular near equatorial ----------------------------");
                    r = new double[] { 7334.858850000, -12704.3481945462, 0.000 };
                    v = new double[] { -4.51428154312046, -2.60632166411836, 0.00000090977789 };
                }

                if (i == 15) // -- CE l = 65 deg
                {
                    strbuild.AppendLine(" coe test circular equatorial ----------------------------");
                    r = new double[] { 6199.6905946008, 13295.2793851394, 0.0 };
                    v = new double[] { -4.72425923942564, 2.20295826245369, 0.0 };
                }
                if (i == 16) // -- CE l = 65 deg i = 180 deg
                {
                    strbuild.AppendLine(" coe test circular equatorial ----------------------------");
                    r = new double[] { 6199.6905946008, -13295.2793851394, 0.0 };
                    v = new double[] { -4.72425923942564, -2.20295826245369, 0.0 };
                }

                // ------- parabolic orbit tests -------------------
                if (i == 17)
                {
                    strbuild.AppendLine(" coe test parabolic ----------------------------");
                    r = new double[] { 0.5916109 * AstroLibr.gravConst.re, -1.2889359 * AstroLibr.gravConst.re, -0.3738343 * AstroLibr.gravConst.re };
                    v = new double[] { 1.1486347 * AstroLibr.gravConst.velkmps, -0.0808249 * AstroLibr.gravConst.velkmps, -0.1942733 * AstroLibr.gravConst.velkmps };
                }

                if (i == 18)
                {
                    strbuild.AppendLine(" coe test parabolic ----------------------------");
                    r = new double[] { -1.0343646 * AstroLibr.gravConst.re, -0.4814891 * AstroLibr.gravConst.re, 0.1735524 * AstroLibr.gravConst.re };
                    v = new double[] { 0.1322278 * AstroLibr.gravConst.velkmps, 0.7785322 * AstroLibr.gravConst.velkmps, 1.0532856 * AstroLibr.gravConst.velkmps };
                }

                if (i == 19)
                {
                    strbuild.AppendLine(" coe test hyperbolic ---------------------------");
                    r = new double[] { 0.9163903 * AstroLibr.gravConst.re, 0.7005747 * AstroLibr.gravConst.re, -1.3909623 * AstroLibr.gravConst.re };
                    v = new double[] { 0.1712704 * AstroLibr.gravConst.velkmps, 1.1036199 * AstroLibr.gravConst.velkmps, -0.3810377 * AstroLibr.gravConst.velkmps };
                }

                if (i == 20)
                {
                    strbuild.AppendLine(" coe test hyperbolic ---------------------------");
                    r = new double[] { 12.3160223 * AstroLibr.gravConst.re, -7.0604653 * AstroLibr.gravConst.re, -3.7883759 * AstroLibr.gravConst.re };
                    v = new double[] { -0.5902725 * AstroLibr.gravConst.velkmps, 0.2165037 * AstroLibr.gravConst.velkmps, 0.1628339 * AstroLibr.gravConst.velkmps };
                }

                if (i == 21)
                {
                    strbuild.AppendLine(" coe test rectilinear --------------------------");
                    r = new double[] { -1984.03023322569, 1525.27235370582, 6364.76955283447 };
                    v = new double[] { -1.60595491095, 1.23461759098, 5.15190381139 };
                    v = new double[] { -1.60936089585, 1.23723602618, 5.16283021192 };  // 196
                }

                strbuild.AppendLine(" start r " + r[0].ToString(fmt) + " " + r[1].ToString(fmt) + " " + r[2].ToString(fmt) + " " +
                "v " + v[0].ToString(fmt) + " " + v[1].ToString(fmt) + " " + v[2].ToString(fmt));
                // --------  coe2rv       - classical elements to posisiotn and velocity
                // --------  rv2coe       - position and velocity vectors to classical elements
                AstroLibr.rv2coe(r, v, out p, out a, out ecc, out incl, out raan, out argp, out nu, out m, out arglat, out truelon, out lonper);
                strbuild.AppendLine("           p km       a km      ecc      incl deg     raan deg     argp deg      nu deg      m deg      arglat   truelon    lonper");
                strbuild.AppendLine("ans coes " + p.ToString().PadLeft(17) + " " + a.ToString(fmt).PadLeft(17) + " " + ecc.ToString("0.000000000") + " " + (incl * rad).ToString(fmt) + " " +
                          (raan * rad).ToString(fmt) + " " + (argp * rad).ToString(fmt) + " " + (nu * rad).ToString(fmt) + " " + (m * rad).ToString(fmt) + " " +
                          (arglat * rad).ToString(fmt) + " " + (truelon * rad).ToString(fmt) + " " + (lonper * rad).ToString(fmt));

                // rectilinear orbits have sign(a) determines orbit type, arglat
                // is nu, but the magnitude is off...?
                if (Math.Abs(ecc - 1.0) < 0.0000001)
                    p = MathTimeLibr.mag(r) * 1.301;
                AstroLibr.coe2rv(p, ecc, incl, raan, argp, nu, arglat, truelon, lonper, out r1, out v1);
                strbuild.AppendLine(" end  r " + r1[0].ToString(fmt) + " " + r1[1].ToString(fmt) + " " + r1[2].ToString(fmt) + " " +
                "v " + v1[0].ToString(fmt) + " " + v1[1].ToString(fmt) + " " + v1[2].ToString(fmt));
            }  // through for
        }

        public void testfindc2c3()
        {
            double znew, c2new, c3new;

            // --------  findc2c3     - find c2 c3 parameters for f and g battins method
            znew = -39.47842;
            AstroLibr.findc2c3(znew, out c2new, out c3new);
            strbuild.AppendLine("findc2c3 z " + znew.ToString(fmt) + " " + c2new.ToString(fmt) + " " + c3new.ToString(fmt));

            znew = 0.0;
            AstroLibr.findc2c3(znew, out c2new, out c3new);
            strbuild.AppendLine("findc2c3 z " + znew.ToString(fmt) + " " + c2new.ToString(fmt) + " " + c3new.ToString(fmt));

            znew = 0.57483;
            AstroLibr.findc2c3(znew, out c2new, out c3new);
            strbuild.AppendLine("findc2c3 z " + znew.ToString(fmt) + " " + c2new.ToString(fmt) + " " + c3new.ToString(fmt));

            znew = 39.47842;
            AstroLibr.findc2c3(znew, out c2new, out c3new);
            strbuild.AppendLine("findc2c3 z " + znew.ToString(fmt) + " " + c2new.ToString(fmt) + " " + c3new.ToString(fmt));
        }


        public void testcoe2rv()
        {
            double[] r = new double[3];
            double[] v = new double[3];
            double p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper;

            r = new double[] { -605.79221660, -5870.22951108, 3493.05319896 };
            v = new double[] { -1.56825429, -3.70234891, -6.47948395 };
            AstroLibr.rv2coe(r, v, out p, out a, out ecc, out incl, out raan, out argp, out nu, out m, out arglat, out truelon, out lonper);

            AstroLibr.coe2rv(p, ecc, incl, raan, argp, nu, arglat, truelon, lonper, out r, out v);

            strbuild.AppendLine("coe2rv r " + r[0].ToString(fmt) + " " + r[1].ToString(fmt) + " " + r[2].ToString(fmt) + " " +
                "v " + v[0].ToString(fmt) + " " + v[1].ToString(fmt) + " " + v[2].ToString(fmt));
        }

        public void testlon2nu()
        {
            double rad = 180.0 / Math.PI;
            double jdut1, lon, incl, raan, argp;
            string strtext;
            jdut1 = 2449470.5;
            incl = 35.324598 / rad;
            lon = -121.3487 / rad;
            raan = 45.0 / rad;
            argp = 34.456798 / rad;

            lon = AstroLibr.lon2nu(jdut1, lon, incl, raan, argp, out strtext);

        }

        // faster version?
        public void testnewtonmx()
        {
            double rad = 180.0 / Math.PI;
            double ecc, eccanom, m, nu;
            ecc = 0.4;
            m = 334.566986 / rad;

            AstroLibr.newtonmx(ecc, m, out eccanom, out nu);

            strbuild.AppendLine(" newtonmx " + ecc.ToString() + " m " + (m * rad).ToString() +
                " eccanom " + (eccanom * rad).ToString() + " nu " + (nu * rad).ToString());
        }

        // --------  newtonm      - find eccentric and true anomaly given ecc and mean anomaly
        public void testnewtonm()
        {
            double rad = 180.0 / Math.PI;
            double ecc, eccanom, m, nu;
            ecc = 0.4;
            eccanom = 334.566986 / rad;
            AstroLibr.newtone(ecc, eccanom, out m, out nu);

            strbuild.AppendLine(" newtone ecc " + ecc.ToString(fmt) + " eccanom " + (eccanom * rad).ToString(fmt) +
                " m " + (m * rad).ToString(fmt) + " nu " + (nu * rad).ToString(fmt));

            ecc = 0.34;
            m = 235.4 / rad;
            AstroLibr.newtonm(ecc, m, out eccanom, out nu);
            strbuild.AppendLine(" newtonm ecc " + ecc.ToString(fmt) + " m " + (m * rad).ToString(fmt) +
                " eccanom " + (eccanom * rad).ToString(fmt) + " nu " + (nu * rad).ToString(fmt));
        }


        // --------  newtone      - find true and mean anomaly given ecc and eccentric anomaly
        public void testnewtone()
        {
            double rad = 180.0 / Math.PI;
            double ecc, eccanom, m, nu;
            ecc = 0.34;
            eccanom = 334.566986 / rad;
            AstroLibr.newtone(ecc, eccanom, out m, out nu);

            strbuild.AppendLine(" newtone ecc " + ecc.ToString(fmt) + " eccanom " + (eccanom * rad).ToString(fmt) +
                " m " + (m * rad).ToString(fmt) + " nu " + (nu * rad).ToString(fmt));
        }

        // --------  newtonnu     - find eccentric and mean anomaly given ecc and true anomaly
        public void testnewtonnu()
        {
            double rad = 180.0 / Math.PI;
            double ecc, eccanom, m, nu;
            ecc = 0.34;
            nu = 134.567001 / rad;

            AstroLibr.newtonnu(ecc, nu, out eccanom, out m);

            strbuild.AppendLine(" newtonnu ecc " + ecc.ToString(fmt) + " nu " + (nu * rad).ToString(fmt) +
                " eccanom " + (eccanom * rad).ToString(fmt) + " m " + (m * rad).ToString(fmt));
        }


        public void keplerc2c3
        (
            double[] r1, double[] v1, double dtseco, out double[] r2, out double[] v2,
            out double c2new, out double c3new, out double xnew, out double znew
        )
        {
            // -------------------------  implementation   -----------------
            int ktr, i, numiter, mulrev;
            double[] h = new double[3];
            double[] rx = new double[3];
            double[] vx = new double[3];
            double f, g, fdot, gdot, rval, xold, xoldsqrd,
                  xnewsqrd, p, dtnew, rdotv, a, dtsec,
                  alpha, sme, period, s, w, temp, magro, magvo, magh, magr, magv;
            char show;
            //string errork;

            show = 'n';
            double small, twopi, halfpi;

            for (int ii = 0; ii < 3; ii++)
            {
                rx[ii] = 0.0;
                vx[ii] = 0.0;
            }
            r2 = rx;  // seems to be the only way to get these variables out
            v2 = vx;
            xnew = 0.0;
            c2new = 0.0;
            c3new = 0.0;

            small = 0.000000001;
            twopi = 2.0 * Math.PI;
            halfpi = Math.PI * 0.5;

            // -------------------------  implementation   -----------------
            // set constants and intermediate printouts
            numiter = 50;

            if (show == 'y')
            {
                //            printf(" r1 %16.8f %16.8f %16.8f ER \n",r1[0]/AstroLibr.gravConst.re,r1[1]/AstroLibr.gravConst.re,r1[2]/AstroLibr.gravConst. );
                //            printf(" vo %16.8f %16.8f %16.8f ER/TU \n",vo[0]/velkmps, vo[1]/velkmps, vo[2]/velkmps );
            }

            // --------------------  initialize values   -------------------
            ktr = 0;
            xold = 0.0;
            znew = 0.0;
            //errork = "      ok";
            dtsec = dtseco;
            mulrev = 0;

            if (Math.Abs(dtseco) > small)
            {
                magro = MathTimeLibr.mag(r1);
                magvo = MathTimeLibr.mag(v1);
                rdotv = MathTimeLibr.dot(r1, v1);

                // -------------  find sme, alpha, and a  ------------------
                sme = ((magvo * magvo) * 0.5) - (AstroLibr.gravConst.mu / magro);
                alpha = -sme * 2.0 / AstroLibr.gravConst.mu;

                if (Math.Abs(sme) > small)
                    a = -AstroLibr.gravConst.mu / (2.0 * sme);
                else
                    a = 999999.9;
                if (Math.Abs(alpha) < small)   // parabola
                    alpha = 0.0;

                if (show == 'y')
                {
                    //           printf(" sme %16.8f  a %16.8f alp  %16.8f ER \n",sme/(AstroLibr.gravConst.mu/AstroLibr.gravConst.), a/AstroLibr.gravConst.re, alpha * AstroLibr.gravConst. );
                    //           printf(" sme %16.8f  a %16.8f alp  %16.8f km \n",sme, a, alpha );
                    //           printf(" ktr      xn        psi           r2          xn+1        dtn \n" );
                }

                // ------------   setup initial guess for x  ---------------
                // -----------------  circle and ellipse -------------------
                if (alpha >= small)
                {
                    period = twopi * Math.Sqrt(Math.Abs(a * a * a) / AstroLibr.gravConst.mu);
                    // ------- next if needed for 2body multi-rev ----------
                    if (Math.Abs(dtseco) > Math.Abs(period))
                        // including the truncation will produce vertical lines that are parallel
                        // (plotting chi vs time)
                        //                    dtsec = rem( dtseco,period );
                        mulrev = Convert.ToInt16(dtseco / period);
                    if (Math.Abs(alpha - 1.0) > small)
                        xold = Math.Sqrt(AstroLibr.gravConst.mu) * dtsec * alpha;
                    else
                        // - first guess can't be too close. ie a circle, r2=a
                        xold = Math.Sqrt(AstroLibr.gravConst.mu) * dtsec * alpha * 0.97;
                }
                else
                {
                    // --------------------  parabola  ---------------------
                    if (Math.Abs(alpha) < small)
                    {
                        MathTimeLibr.cross(r1, v1, out h);
                        magh = MathTimeLibr.mag(h);
                        p = magh * magh / AstroLibr.gravConst.mu;
                        s = 0.5 * (halfpi - Math.Atan(3.0 * Math.Sqrt(AstroLibr.gravConst.mu / (p * p * p)) * dtsec));
                        w = Math.Atan(Math.Pow(Math.Tan(s), (1.0 / 3.0)));
                        xold = Math.Sqrt(p) * (2.0 * MathTimeLibr.cot(2.0 * w));
                        alpha = 0.0;
                    }
                    else
                    {
                        // ------------------  hyperbola  ------------------
                        temp = -2.0 * AstroLibr.gravConst.mu * dtsec /
                              (a * (rdotv + Math.Sign(dtsec) * Math.Sqrt(-AstroLibr.gravConst.mu * a) *
                              (1.0 - magro * alpha)));
                        xold = Math.Sign(dtsec) * Math.Sqrt(-a) * Math.Log(temp);
                    }
                } // if alpha

                ktr = 1;
                dtnew = -10.0;
                double tmp = 1.0 / Math.Sqrt(AstroLibr.gravConst.mu);
                while ((Math.Abs(dtnew * tmp - dtsec) >= small) && (ktr < numiter))
                {
                    xoldsqrd = xold * xold;
                    znew = xoldsqrd * alpha;

                    // ------------- find c2 and c3 functions --------------
                    AstroLibr.findc2c3(znew, out c2new, out c3new);

                    // ------- use a newton iteration for new values -------
                    rval = xoldsqrd * c2new + rdotv * tmp * xold * (1.0 - znew * c3new) +
                             magro * (1.0 - znew * c2new);
                    dtnew = xoldsqrd * xold * c3new + rdotv * tmp * xoldsqrd * c2new +
                             magro * xold * (1.0 - znew * c3new);

                    // ------------- calculate new value for x -------------
                    xnew = xold + (dtsec * Math.Sqrt(AstroLibr.gravConst.mu) - dtnew) / rval;

                    // ----- check if the univ param goes negative. if so, use bissection
                    if (xnew < 0.0)
                        xnew = xold * 0.5;

                    if (show == 'y')
                    {
                        //  printf("%3i %11.7f %11.7f %11.7f %11.7f %11.7f \n", ktr,xold,znew,rval,xnew,dtnew);
                        //  printf("%3i %11.7f %11.7f %11.7f %11.7f %11.7f \n", ktr,xold/sqrt(AstroLibr.gravConst.),znew,rval/AstroLibr.gravConst.re,xnew/sqrt(AstroLibr.gravConst.),dtnew/sqrt(mu));
                    }

                    ktr = ktr + 1;
                    xold = xnew;
                }  // while

                if (ktr >= numiter)
                {
                    //errork = "knotconv";
                    //           printf("not converged in %2i iterations \n",numiter );
                    for (i = 0; i < 3; i++)
                    {
                        v2[i] = 0.0;
                        r2[i] = v2[i];
                    }
                }
                else
                {
                    // --- find position and velocity vectors at new time --
                    xnewsqrd = xnew * xnew;
                    f = 1.0 - (xnewsqrd * c2new / magro);
                    g = dtsec - xnewsqrd * xnew * c3new / Math.Sqrt(AstroLibr.gravConst.mu);

                    for (i = 0; i < 3; i++)
                        r2[i] = f * r1[i] + g * v1[i];
                    magr = MathTimeLibr.mag(r2);
                    gdot = 1.0 - (xnewsqrd * c2new / magr);
                    fdot = (Math.Sqrt(AstroLibr.gravConst.mu) * xnew / (magro * magr)) * (znew * c3new - 1.0);
                    for (i = 0; i < 3; i++)
                        v2[i] = fdot * r1[i] + gdot * v1[i];
                    magv = MathTimeLibr.mag(v2);
                    temp = f * gdot - fdot * g;
                    //if (Math.Abs(temp - 1.0) > 0.00001)
                    //    errork = "fandg";

                    if (show == 'y')
                    {
                        //           printf("f %16.8f g %16.8f fdot %16.8f gdot %16.8f \n",f, g, fdot, gdot );
                        //           printf("f %16.8f g %16.8f fdot %16.8f gdot %16.8f \n",f, g, fdot, gdot );
                        //           printf("r1 %16.8f %16.8f %16.8f ER \n",r2[0]/AstroLibr.gravConst.re,r2[1]/AstroLibr.gravConst.re,r2[2]/AstroLibr.gravConst. );
                        //           printf("v1 %16.8f %16.8f %16.8f ER/TU \n",v[0]/velkmps, v[1]/velkmps, v[2]/velkmps );
                    }
                }
            } // if fabs
            else
                // ----------- set vectors to incoming since 0 time --------
                for (i = 0; i < 3; i++)
                {
                    r2[i] = r1[i];
                    v2[i] = v1[i];
                }

            //       fprintf( fid,"%11.5f  %11.5f %11.5f  %5i %3i ",znew, dtseco/60.0, xold/(rad), ktr, mulrev );
        }  // keplerc2c3


        public void testfindfandg()
        {
            double[] r1;
            double[] v1;
            double[] r2;
            double[] v2;
            double dtsec, x, c2, c3, z, f, g, fdot, gdot;
            string opt = "pqw";  //  pqw, series, c2c3

            r1 = new double[] { 4938.49830042171, -1922.24810472241, 4384.68293292613 };
            v1 = new double[] { 0.738204644165659, 7.20989453238397, 2.32877392066299 };
            r2 = new double[] { -1105.78023519582, 2373.16130661458, 6713.89444816503 };
            v2 = new double[] { 5.4720951867079, -4.39299050886976, 2.45681739563752 };
            dtsec = 6372.69272563561; // 1ld
            dtsec = 60.0; // must be small step sizes!!

            //// dan hyperbolic test
            //r1 = new double[] { 4070.5942270000000, 3786.8271570000002, 4697.0576309999997 };
            ////v1 = new double[] { -32553.559100851671, -37563.543526937596, -37563.543526937596 };
            //// exact opposite from r velocity
            //v1 = new double[] { -34845.69531184976, -32416.55098811211, -40208.43885307875 };
            //r2 = new double[] { -1105.78023519582, 2373.16130661458, 6713.89444816503 };
            //v2 = new double[] { 5.4720951867079, -4.39299050886976, 2.45681739563752 };
            //dtsec = 0.25; // must be small step sizes!!

            strbuild.AppendLine(" r1 " + r1[0].ToString(fmt) + " " + r1[1].ToString(fmt) + " " + r1[2].ToString(fmt) + " " +
               "v1 " + v1[0].ToString(fmt) + " " + v1[1].ToString(fmt) + " " + v1[2].ToString(fmt));

            for (int i = 0; i <= 5; i++)
            {
                if (i == 0)
                    dtsec = 60.0;
                if (i == 1)
                    dtsec = 0.1;
                if (i == 2)
                    dtsec = 1.0;
                if (i == 3)
                    dtsec = 10.0;
                if (i == 4)
                    dtsec = 100.0;
                if (i == 5)
                    dtsec = 500.0;

                keplerc2c3(r1, v1, dtsec, out r2, out v2, out c2, out c3, out x, out z);
                strbuild.AppendLine(" r2 " + r2[0].ToString(fmt) + " " + r2[1].ToString(fmt) + " " + r2[2].ToString(fmt) + " " +
                    "v2 " + v2[0].ToString(fmt) + " " + v2[1].ToString(fmt) + " " + v2[2].ToString(fmt));

                strbuild.AppendLine("c2 " + c2.ToString(fmt) + " c3 " + c3.ToString(fmt) + " x " +
                    x.ToString(fmt) + " z " + z.ToString(fmt) + " dtsec " + dtsec.ToString(fmt));

                opt = "pqw";
                AstroLibr.findfandg(r1, v1, r2, v2, dtsec, x, c2, c3, z, opt, out f, out g, out fdot, out gdot);
                double ans = f * gdot - g * fdot;
                strbuild.AppendLine("f and g pqw    " + f.ToString(fmt) + " " + g.ToString(fmt) + " " +
                    fdot.ToString(fmt) + " " + gdot.ToString(fmt) + " " + ans.ToString(fmt));

                opt = "series";  //  pqw, series, c2c3
                AstroLibr.findfandg(r1, v1, r2, v2, dtsec, x, c2, c3, z, opt, out f, out g, out fdot, out gdot);
                ans = f * gdot - g * fdot;
                strbuild.AppendLine("f and g series " + f.ToString(fmt) + " " + g.ToString(fmt) + " " +
                    fdot.ToString(fmt) + " " + gdot.ToString(fmt) + " " + ans.ToString(fmt));

                opt = "c2c3";  //  pqw, series, c2c3
                AstroLibr.findfandg(r1, v1, r2, v2, dtsec, x, c2, c3, z, opt, out f, out g, out fdot, out gdot);
                ans = f * gdot - g * fdot;

                strbuild.AppendLine("f and g c2c3   " + f.ToString(fmt) + " " + g.ToString(fmt) + " " +
                    fdot.ToString(fmt) + " " + gdot.ToString(fmt) + " " + ans.ToString(fmt) + "\n");
            }

        }

        public void testcheckhitearth()
        {
            string hitearthstr = "";
            double[] r1 = new double[3];
            double[] v1t = new double[3];
            double[] r2 = new double[3];
            double[] v2t = new double[3];
            double ang, magr1, magr2, cosdeltanu, altpad, rp, a;
            Int32 nrev;
            char hitearth;

            nrev = 0;
            r1 = new double[] { 2.500000 * AstroLibr.gravConst.re, 0.000000, 0.000000 };
            r2 = new double[] { 1.9151111 * AstroLibr.gravConst.re, 1.6069690 * AstroLibr.gravConst.re, 0.000000 };
            // assume circular initial orbit for vel calcs
            v1t = new double[] { 0.0, Math.Sqrt(AstroLibr.gravConst.mu / r1[0]), 0.0 };
            ang = Math.Atan(r2[1] / r2[0]);
            v2t = new double[] { -Math.Sqrt(AstroLibr.gravConst.mu / r2[1]) * Math.Cos(ang), Math.Sqrt(AstroLibr.gravConst.mu / r2[0]) * Math.Sin(ang), 0.0 };
            altpad = 100.0; // km

            magr1 = MathTimeLibr.mag(r1);
            magr2 = MathTimeLibr.mag(r2);
            cosdeltanu = MathTimeLibr.dot(r1, r2) / (magr1 * magr2);

            AstroLibr.checkhitearth(altpad, r1, v1t, r2, v2t, nrev, out hitearth, out hitearthstr, out rp, out a);

            strbuild.AppendLine("hitearth? " + hitearthstr + " " + (Math.Acos(cosdeltanu) * 180.0 / Math.PI).ToString(fmt));
        }

        public void testcheckhitearthc()
        {
            string hitearthstr = "";
            double[] r1c = new double[3];
            double[] v1tc = new double[3];
            double[] r2c = new double[3];
            double[] v2tc = new double[3];
            double ang, magr1c, magr2c, cosdeltanu, altpadc, rp, a;
            Int32 nrev;
            char hitearth;

            nrev = 0;
            r1c = new double[] { 2.500000, 0.000000, 0.000000 };
            r2c = new double[] { 1.9151111, 1.6069690, 0.000000 };
            // assume circular initial orbit for vel calcs
            v1tc = new double[] { 0.0, Math.Sqrt(1.0 / r1c[0]), 0.0 };
            ang = Math.Atan(r2c[1] / r2c[0]);
            v2tc = new double[] { -Math.Sqrt(1.0 / r2c[1]) * Math.Cos(ang), Math.Sqrt(1.0 / r2c[0]) * Math.Sin(ang), 0.0 };
            altpadc = 100.0 / AstroLibr.gravConst.re; // er

            magr1c = MathTimeLibr.mag(r1c);
            magr2c = MathTimeLibr.mag(r2c);
            cosdeltanu = MathTimeLibr.dot(r1c, r2c) / (magr1c * magr2c);
            AstroLibr.checkhitearthc(altpadc, r1c, v1tc, r2c, v2tc, nrev, out hitearth, out hitearthstr, out rp, out a);

            strbuild.AppendLine("hitearth? " + hitearthstr + " " + (Math.Acos(cosdeltanu) * 180.0 / Math.PI).ToString(fmt));
        }


        public void testgibbs()
        {
            double[] r1 = new double[3];
            double[] r2 = new double[3];
            double[] r3 = new double[3];
            double[] v2 = new double[3];
            double copa, theta, theta1, rad;
            string errorstr;

            rad = 180.0 / Math.PI;

            r1 = new double[] { 0.0000000, 0.000000, AstroLibr.gravConst.re };
            r2 = new double[] { 0.0000000, -4464.696, -5102.509 };
            r3 = new double[] { 0.0000000, 5740.323, 3189.068 };

            AstroLibr.gibbs(r1, r2, r3, out v2, out theta, out theta1, out copa, out errorstr);

            strbuild.AppendLine("testgibbs " + v2[0].ToString(fmt) + " " + v2[1].ToString(fmt) + " " + v2[2].ToString(fmt));
            strbuild.AppendLine("testgibbs " + (theta * rad).ToString(fmt) + " " + (theta1 * rad).ToString(fmt) + " " + (copa * rad).ToString(fmt));
        }


        public void testhgibbs()
        {
            double[] r1 = new double[3];
            double[] r2 = new double[3];
            double[] r3 = new double[3];
            double[] v2 = new double[3];
            double copa, theta, theta1, rad, jd1, jd2, jd3;
            string errorstr;

            rad = 180.0 / Math.PI;

            r1 = new double[] { 0.0000000, 0.000000, AstroLibr.gravConst.re };
            r2 = new double[] { 0.0000000, -4464.696, -5102.509 };
            r3 = new double[] { 0.0000000, 5740.323, 3189.068 };
            jd1 = 2451849.5;
            jd2 = jd1 + 1.0 / 1440.0 + 16.48 / 86400.0;
            jd3 = jd1 + 2.0 / 1440.0 + 33.04 / 86400.0;
            AstroLibr.herrgibbs(r1, r2, r3, jd1, jd2, jd3, out v2, out theta, out theta1, out copa, out errorstr);

            strbuild.AppendLine("testherrgibbs " + v2[0].ToString(fmt) + " " + v2[1].ToString(fmt) + " " + v2[2].ToString(fmt));
            strbuild.AppendLine("testherrgibbs " + (theta * rad).ToString(fmt) + " " + (theta1 * rad).ToString(fmt) + " " + (copa * rad).ToString(fmt));
        }



        public void testgeo()
        {
            StringBuilder strbuildObs = new StringBuilder();
            double rad = 180.0 / Math.PI;

            // misc test
            double lona, londot, lons, lonp, z, j22, c22, s22, omegaearth, dt;
            dt = 86400.0;  // 1 day in second
            c22 = 1.57461532572292E-06;
            s22 = -9.03872789196567E-07;
            j22 = Math.Sqrt(c22 * c22 + s22 * s22);
            // stable longitude point
            lons = 0.5 * Math.Atan2(s22, c22);
            omegaearth = 0.000072921158553;   // rad /s
            z = 6.6017;  // rad
            // initial longitude with 0 initial drift
            lonp = 12.0 / rad;
            lonp = lons;

            lona = lonp - 1.0 / rad;
            lona = lons;
            londot = 0.0;
            for (int jj = 0; jj < 400; jj++)
            {
                lona = lona + londot * dt;
                londot = 3.0 * omegaearth / z * Math.Sqrt(2.0 * j22) *
                    Math.Sqrt(Math.Cos(2.0 * (lona - lons)) - Math.Cos(2.0 * (lonp - lons)));
                strbuildObs.AppendLine(jj + " " + lona * rad + " " + (lonp - lons) * rad + " " + londot * rad / 86400.0);
            } // for through all the tracks testing rtasc/decl rates
            //File.WriteAllText(@"D:\faabook\current\excel\testgeo.out", strbuildObs.ToString());
        }








        // test angles-only routines
        // output these results separately to the testall directory
        public void testangles()
        {
            double[] rseci1 = new double[3];
            double[] vseci1 = new double[3];
            double[] rseci2 = new double[3];
            double[] vseci2 = new double[3];
            double[] rseci3 = new double[3];
            double[] vseci3 = new double[3];
            double[] rsecef1 = new double[3];
            double[] rsecef2 = new double[3];
            double[] rsecef3 = new double[3];
            double[] vsecef1 = new double[3];
            double[] vsecef2 = new double[3];
            double[] vsecef3 = new double[3];
            double[] r2 = new double[3];
            double[] v2 = new double[3];
            double rad, lod;
            double[] latgd = new double[15];
            double[] lon = new double[15];
            double[] alt = new double[15];
            double[] trtasc = new double[15];
            double[] tdecl = new double[15];
            double[] initguess = new double[30];
            string errstr;
            Int32[] year = new int[15];
            Int32[] mon = new int[15];
            Int32[] day = new int[15];
            Int32[] hr = new int[15];
            Int32[] minute = new int[15];
            double[] second = new double[15];
            double p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper;

            int dat;
            double[] jd = new double[100];
            double[] jdf = new double[100];
            double jdut1, dut1, jdtt, jdftt, ttt, xp, yp, ddx, ddy, ddpsi, ddeps;
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
            StringBuilder strbuildall = new StringBuilder();
            
            this.opsStatus.Text = "Test Angles ";
            Refresh();

            string fileLoc;
            string ans;
            Int32 ktrActObs;
            string EOPupdate;
            EOPSPWLib.iau80Class iau80arr;
            fileLoc = @"D:\Codes\LIBRARY\DataLib\nut80.dat";
            EOPSPWLibr.iau80in(fileLoc, out iau80arr);

            string eopFileName = @"D:\Codes\LIBRARY\DataLib\EOP-All-v1.1_2020-02-12.txt";
            EOPSPWLibr.readeop(ref EOPSPWLibr.eopdata, eopFileName,  out ktrActObs, out EOPupdate);

            // gooding tests cases from Gooding paper (1997 CMDA)
            double[] los1;
            double[] los2;
            double[] los3;

            // read input data
            // note the input data has a # line between each case
            string infilename = @"D:\Codes\LIBRARY\DataLib\anglestest.dat";
            string[] fileData = File.ReadAllLines(infilename);

            // --- read obs data in
            //int caseopt = 0;  // set this for whichever case to run

            // find mins
            // orbits only need to be close
            //double rtol = 500.0; // km
            //double vtol = 0.1;   // km/s
            //double atol = 500.0; // km
            //double ptol = 500.0; // km
            //double etol = 0.1;  // 
            double itol = 5.0 / rad;   // rad

            for (int caseopt = 0; caseopt <= 24; caseopt++)  // 0  23
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
                while (ktr < fileData.Count() && !fileData[ktr][0].Equals('#'))
                {
                    line = fileData[ktr];
                    linesplt = line.Split(',');
                    mon[obsktr] = Convert.ToInt32(linesplt[1]);
                    day[obsktr] = Convert.ToInt32(linesplt[0]);
                    year[obsktr] = Convert.ToInt32(linesplt[2]);
                    hr[obsktr] = Convert.ToInt32(linesplt[3]);
                    minute[obsktr] = Convert.ToInt32(linesplt[4]);
                    second[obsktr] = Convert.ToDouble(linesplt[5]);
                    MathTimeLibr.jday(year[obsktr], mon[obsktr], day[obsktr], hr[obsktr], minute[obsktr], second[obsktr],
                        out jd[obsktr], out jdf[obsktr]);

                    latgd[obsktr] = Convert.ToDouble(linesplt[6]) / rad;
                    lon[obsktr] = Convert.ToDouble(linesplt[7]) / rad;
                    alt[obsktr] = Convert.ToDouble(linesplt[8]) / rad;

                    trtasc[obsktr] = Convert.ToDouble(linesplt[9]) / rad;
                    tdecl[obsktr] = Convert.ToDouble(linesplt[10]) / rad;
                    if (obsktr == 0)
                        initguess[tmpcase] = Convert.ToDouble(linesplt[11]);  // initial guess in km

                    obsktr = obsktr + 1;
                    ktr = ktr + 1;
                }

                int idx1, idx2, idx3;
                idx1 = 0;
                idx2 = 1;
                idx3 = 2;
                strbuild.AppendLine("/n/n ================================ case number " + caseopt.ToString() + " ================================");
                strbuildall.AppendLine("/n/n ================================ case number " + caseopt.ToString() + " ================================");
                strbuild.AppendLine(ans);
                switch (caseopt)
                {
                    case 0:
                        idx1 = 2;
                        idx2 = 4;
                        idx3 = 5;

                        break;
                    case 1:
                        // book example
                        //dut1 = -0.609641;      // second
                        //dat = 35;              // second
                        //lod = 0.0;
                        //xp = 0.137495 * conv;  // " to rad
                        //yp = 0.342416 * conv;
                        //ddpsi = 0.0;  // " to rad
                        //ddeps = 0.0;
                        //ddx = 0.0;    // " to rad
                        //ddy = 0.0;
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
                    case int n when (n >= 2 && n <= 12):
                        idx1 = 0;
                        idx2 = 1;
                        idx3 = 2;
                        break;
                }  // end switch

                //    // herrick interplantetary Herrick pp 384-5 & 418 (gibbs)
                //    // units of days in 1910 (?) and au
                //    // says mu = 1.0
                //    // k = 0, rho = 5.9 and 5.9 au?
                //    // some initi rho of 1000000 still works
                //    // days 7.8205, 26.7480, 48.6262 (in 1910 Nov)
                //    //rtasc = 3-50-24.3, decl = 25-11-10.5
                //    //rtasc = 3-13-3.0  decl = 22-29-31.3
                //    //rtasc = 4-54-19.5  decl = 20-14-51.9
                //    double tau12 = 0.325593;
                //    double tau13 = 0.701944;
                //    rseci1 = new double[] { 0.7000687 * AstroLibr.astroConsts.au,
                //        0.6429399 * AstroLibr.astroConsts.au, 0.2789211 * AstroLibr.astroConsts.au };
                //    rseci2 = new double[] { 0.4306907, 0.8143496, 0.3532745 };
                //    rseci3 = new double[] { 0.0628371, 0.9007098, 0.3907417 };
                //    los1 = new double[] { 0.9028975, 0.0606048, 0.4255621 };
                //    los2 = new double[] { 0.9224764, 0.0518570, 0.3825549 };
                //    los3 = new double[] { 0.9347684, 0.0802269, 0.3460802 };
                //    //(light-corrected times not used but aze 0.325578 and 0.701903)
                //    break;
                //case 2:
                //    // extreme case Lane ex 2
                //    tau12 = 1.570796327;
                //    tau13 = 3.141592654;
                //    rseci1 = new double[] { 0.0, -1.0, 0.0 };
                //    rseci2 = new double[] { 1.0, 0.0, 0.0 };
                //    rseci3 = new double[] { 0.0, 1.0, 0.0 };
                //    los1 = new double[] { 0.0, 1.0, 0.0 };
                //    los2 = new double[] { 0.0, 0.0, 1.0 };
                //    los3 = new double[] { 0.0, -1.0, 0.0 };

                //    break;
                //case 3:
                //    // rectilinear case
                //    tau12 = 1.047197552;
                //    tau13 = 2.960420507;
                //    rseci1 = new double[] { 1.0, 0.0, 0.0 };
                //    rseci2 = new double[] { 0.0, 1.0, 0.0 };
                //    rseci3 = new double[] { 0.0, 0.0, 0.0 };
                //    los1 = new double[] { -1.0, 0.0, 0.5 };
                //    los2 = new double[] { 0.0, -1.0, 1.5 };
                //    los3 = new double[] { 0.0, 0.0, 2.0 };
                //    break;
                //case 5:
                //    // revised escobal example 
                //    tau12 = 0.0381533;
                //    tau13 = 0.0399364;
                //    rseci1 = new double[] { 0.16606957, 0.84119785, -0.51291356 };
                //    rseci2 = new double[] { -0.73815134, -0.41528280, 0.53035336 };
                //    rseci3 = new double[] { -0.73343987, -0.42352540, 0.53037164 };
                //    los1 = new double[] { -0.92475472, -0.37382824, -0.07128226 };
                //    los2 = new double[] { 0.80904274, -0.55953385, 0.17992142 };
                //    los3 = new double[] { 0.85044131, -0.49106628, 0.18868888 };
                //    break;
                //case 7:
                //    // example from Lane
                //    tau12 = 2.2;
                //    tau13 = 0.35225232;
                //    rseci2 = new double[] { 0.89263524, 0.28086002, 0.35277012 };
                //    rseci3 = new double[] { -0.02703285, 0.93585748, 0.35152067 };
                //    los1 = new double[] { 0.76526944, 0.12314580, 0.63182102 };
                //    los2 = new double[] { -0.21266402, -0.54295751, 0.81238609 };
                //    los3 = new double[] { 0.52029946, -0.39083440, 0.75930030 };
                //    break;


                // state already exists use that period
                // TLE exists use that period
                // otherwise options, 95 min, 108 min, 150 min, 250 min, 7.2 hr, 12 hr, and 24 hr 

                //                for (int z = 0; z <= -10; z++)
                //{
                //    switch (z)
                //    {
                //        case 0:
                //            idx1 = 2;
                //            idx2 = 4;
                //            idx3 = 5;
                //            break;
                //        case 1:
                //            idx1 = 0;
                //            idx2 = 2;
                //            idx3 = 3;
                //            break;
                //        case 2:
                //            idx1 = 0;
                //            idx2 = 2;
                //            idx3 = 9;
                //            break;
                //        case 3:
                //            idx1 = 0;
                //            idx2 = 2;
                //            idx3 = 6;
                //            break;
                //        case 4:
                //            idx1 = 0;
                //            idx2 = 2;
                //            idx3 = 7;
                //            break;
                //        case 5:
                //            idx1 = 0;
                //            idx2 = 2;
                //            idx3 = 12;
                //            break;
                //        case 6:
                //            idx1 = 2;
                //            idx2 = 4;
                //            idx3 = 14;
                //            break;
                //        case 7:
                //            idx1 = 2;
                //            idx2 = 4;
                //            idx3 = 7;
                //            break;
                //        case 8:
                //            idx1 = 0;
                //            idx2 = 13;
                //            idx3 = 14;
                //            break;
                //        case 9:
                //            idx1 = 2;
                //            idx2 = 4;
                //            idx3 = 10;
                //            break;
                //        case 10:
                //            idx1 = 10;
                //            idx2 = 11;
                //            idx3 = 12;
                //            break;
                //    }
                //    strbuild.AppendLine("\nz " + z.ToString());




                //could put in separate doangles function




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
                // note you have to use tdb for time of ineterst AND j2000 (when dat = 32)
                //  ttt = (jd + jdFrac + (dat + 32.184) / 86400.0 - 2451545.0 - (32 + 32.184) / 86400.0) / 36525.0;
                jdut1 = jd[idx1] + jdf[idx1] + dut1 / 86400.0;
                AstroLibr.eci_ecef(ref rseci1, ref vseci1, MathTimeLib.Edirection.efrom, ref rsecef1, ref vsecef1,
                     iau80arr, jdtt, jdftt, jdut1, lod, xp, yp, ddpsi, ddeps);

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
                AstroLibr.eci_ecef(ref rseci2, ref vseci2, MathTimeLib.Edirection.efrom, ref rsecef2, ref vsecef2,
                    iau80arr, jdtt, jdftt, jdut1, lod, xp, yp, ddpsi, ddeps);
                double gst, lst;
                AstroLibr.lstime(lon[idx2], jdut1, out lst, out gst);
                strbuildall.AppendLine("\nlst " + lst.ToString() + " " + (lst * rad).ToString());


                AstroLibr.site(latgd[idx3], lon[idx3], alt[idx3], out rsecef3, out vsecef3);
                EOPSPWLibr.findeopparam(jd[idx3], jdf[idx3], 's', EOPSPWLibr.eopdata, 
                    out dut1, out dat, out lod, out xp, out yp, out ddpsi, out ddeps, out ddx, out ddy);
                jdtt = jd[idx3];
                jdftt = jdf[idx3] + (dat + 32.184) / 86400.0;
                ttt = (jdtt + jdftt - 2451545.0) / 36525.0;
                jdut1 = jd[idx3] + jdf[idx3] + dut1 / 86400.0;
                AstroLibr.eci_ecef(ref rseci3, ref vseci3, MathTimeLib.Edirection.efrom, ref rsecef3, ref vsecef3,
                   iau80arr, jdtt, jdftt, jdut1, lod, xp, yp, ddpsi, ddeps);

                if (Math.Abs(latgd[idx1] - latgd[idx2]) < 0.001 && Math.Abs(latgd[idx1] - latgd[idx3]) < 0.001
                    && Math.Abs(lon[idx1] - lon[idx2]) < 0.001 && Math.Abs(lon[idx1] - lon[idx3]) < 0.001)
                    diffsites = 'n';
                else
                    diffsites = 'y';


                // write output
                strbuildall.AppendLine("rseci1 " + rseci1[0].ToString("0.000000") + " " +
                rseci1[1].ToString("0.000000") + " " + rseci1[2].ToString("0.000000"));
                strbuildall.AppendLine("rseci2 " + rseci2[0].ToString("0.000000") + " " +
                    rseci2[1].ToString("0.000000") + " " + rseci2[2].ToString("0.000000"));
                strbuildall.AppendLine("rseci3 " + rseci3[0].ToString("0.000000") + " " +
                    rseci3[1].ToString("0.000000") + " " + rseci3[2].ToString("0.000000"));

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
                strbuild.AppendLine("Laplace -----------------------------------");
                strbuildall.AppendLine("Laplace -----------------------------------");

                strbuildall.AppendLine("\n\ninputs: \n");
                strbuildall.AppendLine("Site obs1 "
                    + rseci1[0].ToString() + " " + rseci1[1].ToString() + " " + rseci1[2].ToString()
                    + " km  lat " + (latgd[idx1] * rad).ToString() + " lon " + (lon[idx1] * rad).ToString() + " "
                    + alt[idx1].ToString());
                strbuildall.AppendLine("Site obs2 "
                    + rseci2[0].ToString() + " " + rseci2[1].ToString() + " " + rseci2[2].ToString()
                    + " km  lat " + (latgd[idx2] * rad).ToString() + " lon " + (lon[idx2] * rad).ToString() + " "
                    + alt[idx2].ToString());
                strbuildall.AppendLine("Site obs3 "
                    + rseci3[0].ToString() + " " + rseci3[1].ToString() + " " + rseci3[2].ToString()
                    + " km  lat " + (latgd[idx3] * rad).ToString() + " lon " + (lon[idx3] * rad).ToString() + " "
                    + alt[idx3].ToString());
                MathTimeLibr.invjday(jd[idx1], jdf[idx1], out iyear1, out imon1, out iday1, out ihr1, out iminute1, out isecond1);
                strbuildall.AppendLine("obs#1 " + iyear1.ToString() + " " + imon1.ToString() + " " + iday1.ToString()
                    + " " + ihr1.ToString("00") + " " + iminute1.ToString("00") + " " + isecond1.ToString("0.000")
                    + " " + (trtasc[idx1] * rad).ToString() + " " + (tdecl[idx1] * rad).ToString().ToString());
                MathTimeLibr.invjday(jd[idx2], jdf[idx2], out iyear2, out imon2, out iday2, out ihr2, out iminute2, out isecond2);
                strbuildall.AppendLine("obs#2 " + iyear2.ToString() + " " + imon2.ToString() + " " + iday2.ToString()
                    + " " + ihr2.ToString("00") + " " + iminute2.ToString("00") + " " + isecond2.ToString("0.000")
                    + " " + (trtasc[idx2] * rad).ToString() + " " + (tdecl[idx2] * rad).ToString().ToString());
                MathTimeLibr.invjday(jd[idx3], jdf[idx3], out iyear3, out imon3, out iday3, out ihr3, out iminute3, out isecond3);
                strbuildall.AppendLine("obs#3 " + iyear3.ToString() + " " + imon3.ToString() + " " + iday3.ToString()
                    + " " + ihr3.ToString("00") + " " + iminute3.ToString("00") + " " + isecond3.ToString("0.000")
                    + " " + (trtasc[idx3] * rad).ToString() + " " + (tdecl[idx3] * rad).ToString().ToString());
                //if (caseopt == 2)
                //    diffsites = 'y';
                //else 
                //diffsites = 'n';

                AstroLibr.angleslaplace(tdecl[idx1], tdecl[idx2], tdecl[idx3], trtasc[idx1], trtasc[idx2], trtasc[idx3],
                            jd[idx1], jdf[idx1], jd[idx2], jdf[idx2], jd[idx3], jdf[idx3],
                            diffsites, rseci1, rseci2, rseci3, out r2, out v2, out bigr2x, out errstr);
                strbuildall.AppendLine(errstr);
                strbuildall.AppendLine("r2 " + r2[0].ToString("0.000000") + " " +
                    r2[1].ToString("0.000000") + " " + r2[2].ToString("0.000000")
                    + "v2 " + v2[0].ToString("0.000000") + " " +
                    v2[1].ToString("0.000000") + " " + v2[2].ToString("0.000000"));
                AstroLibr.rv2coe(r2, v2, out p, out a, out ecc, out incl, out raan, out argp, out nu, out m, out arglat, out truelon, out lonper);
                strbuildall.AppendLine("\nlaplace coes a= " + a.ToString("0.0000") + " e= " + ecc.ToString("0.000000000") + " i= " + (incl * rad).ToString("0.0000") + " " +
                    (raan * rad).ToString("0.0000") + " " + (argp * rad).ToString("0.0000") + " " + (nu * rad).ToString("0.0000") + " " + (m * rad).ToString("0.0000") + " " +
                    (arglat * rad).ToString("0.0000")); // + " " + (truelon * rad).ToString("0.0000") + " " + (lonper * rad).ToString("0.0000"));
                strbuild.AppendLine("laplace coes a= " + a.ToString("0.0000") + " e= " + ecc.ToString("0.000000000") + " i= " + (incl * rad).ToString("0.0000") + " " +
                    (raan * rad).ToString("0.0000") + " " + (argp * rad).ToString("0.0000") + " " + (nu * rad).ToString("0.0000") + " " + (m * rad).ToString("0.0000") + " " +
                    (arglat * rad).ToString("0.0000")); // + " " + (truelon * rad).ToString("0.0000") + " " + (lonper * rad).ToString("0.0000"));

                strbuild.AppendLine("Gauss  -----------------------------------");
                strbuildall.AppendLine("Gauss  -----------------------------------");
                if (caseopt == 23)
                {  // curtis example -many mistakes!
                    rseci1 = new double[] { 3489.8, 3430.2, 4078.5 };
                    rseci2 = new double[] { 3460.1, 3460.1, 4078.5 };
                    rseci3 = new double[] { 3429.9, 3490.1, 4078.5 };
                }
                AstroLibr.anglesgauss(tdecl[idx1], tdecl[idx2], tdecl[idx3], trtasc[idx1], trtasc[idx2], trtasc[idx3],
                     jd[idx1], jdf[idx1], jd[idx2], jdf[idx2], jd[idx3], jdf[idx3],
                     rseci1, rseci2, rseci3, out r2, out v2, out errstr);
                strbuildall.AppendLine(errstr);
                strbuildall.AppendLine("r2 " + r2[0].ToString("0.000000") + " " +
                     r2[1].ToString("0.000000") + " " + r2[2].ToString("0.000000")
                     + "v2 " + v2[0].ToString("0.000000") + " " +
                     v2[1].ToString("0.000000") + " " + v2[2].ToString("0.000000"));
                AstroLibr.rv2coe(r2, v2, out p, out a, out ecc, out incl, out raan, out argp, out nu, out m, out arglat, out truelon, out lonper);
                strbuildall.AppendLine("gauss coes a= " + a.ToString("0.0000") + " e= " + ecc.ToString("0.000000000") + " i= " + (incl * rad).ToString("0.0000") + " " +
                     (raan * rad).ToString("0.0000") + " " + (argp * rad).ToString("0.0000") + " " + (nu * rad).ToString("0.0000") + " " + (m * rad).ToString("0.0000") + " " +
                     (arglat * rad).ToString("0.0000")); // + " " + (truelon * rad).ToString("0.0000") + " " + (lonper * rad).ToString("0.0000"));
                strbuild.AppendLine("gauss coes a= " + a.ToString("0.0000") + " e= " + ecc.ToString("0.000000000") + " i= " + (incl * rad).ToString("0.0000") + " " +
                     (raan * rad).ToString("0.0000") + " " + (argp * rad).ToString("0.0000") + " " + (nu * rad).ToString("0.0000") + " " + (m * rad).ToString("0.0000") + " " +
                     (arglat * rad).ToString("0.0000")); // + " " + (truelon * rad).ToString("0.0000") + " " + (lonper * rad).ToString("0.0000"));

                double pctchg = 0.05;
                strbuild.AppendLine("Double-r -----------------------------------");
                strbuildall.AppendLine("Double-r -----------------------------------");
                // initial guesses needed for double-r and Gooding
                // use result from Gauss as it's usually pretty good
                // this seems to really help Gooding!!
                AstroLibr.getGaussRoot(tdecl[idx1], tdecl[idx2], tdecl[idx3], trtasc[idx1], trtasc[idx2], trtasc[idx3],
                     jd[idx1], jdf[idx1], jd[idx2], jdf[idx2], jd[idx3], jdf[idx3],
                     rseci1, rseci2, rseci3, out bigr2x);
                initguess[caseopt] = bigr2x;

                rng1 = initguess[caseopt];  // old 12500 needs to be in km!! seems to do better when all the same? if too far off (*2) NAN
                rng2 = initguess[caseopt] * 1.02;  // 1.02 might be better? make the initial guess a bit different
                rng3 = initguess[caseopt] * 1.08;
                AstroLibr.anglesdoubler(tdecl[idx1], tdecl[idx2], tdecl[idx3], trtasc[idx1], trtasc[idx2], trtasc[idx3],
                     jd[idx1], jdf[idx1], jd[idx2], jdf[idx2], jd[idx3], jdf[idx3],
                     rseci1, rseci2, rseci3, rng1, rng2, out r2, out v2, out errstr, pctchg);
                strbuildall.AppendLine(errstr);
                strbuildall.AppendLine("r2 " + r2[0].ToString("0.000000") + " " +
     r2[1].ToString("0.000000") + " " + r2[2].ToString("0.000000")
     + "v2 " + v2[0].ToString("0.000000") + " " +
     v2[1].ToString("0.000000") + " " + v2[2].ToString("0.000000"));
                AstroLibr.rv2coe(r2, v2, out p, out a, out ecc, out incl, out raan, out argp, out nu, out m, out arglat, out truelon, out lonper);
                strbuildall.AppendLine("doubler coes a= " + a.ToString("0.0000") + " e= " + ecc.ToString("0.000000000") + " i= " + (incl * rad).ToString("0.0000") + " " +
                    (raan * rad).ToString("0.0000") + " " + (argp * rad).ToString("0.0000") + " " + (nu * rad).ToString("0.0000") + " " + (m * rad).ToString("0.0000") + " " +
                    (arglat * rad).ToString("0.0000")); // + " " + (truelon * rad).ToString("0.0000") + " " + (lonper * rad).ToString("0.0000"));
                strbuild.AppendLine("doubler coes a= " + a.ToString("0.0000") + " e= " + ecc.ToString("0.000000000") + " i= " + (incl * rad).ToString("0.0000") + " " +
                    (raan * rad).ToString("0.0000") + " " + (argp * rad).ToString("0.0000") + " " + (nu * rad).ToString("0.0000") + " " + (m * rad).ToString("0.0000") + " " +
                    (arglat * rad).ToString("0.0000")); // + " " + (truelon * rad).ToString("0.0000") + " " + (lonper * rad).ToString("0.0000"));


                strbuild.AppendLine("Gooding -----------------------------------");
                strbuildall.AppendLine("Gooding -----------------------------------");
                numhalfrev = 0;

                AstroLibr.getGaussRoot(tdecl[idx1], tdecl[idx2], tdecl[idx3], trtasc[idx1], trtasc[idx2], trtasc[idx3],
                     jd[idx1], jdf[idx1], jd[idx2], jdf[idx2], jd[idx3], jdf[idx3],
                     rseci1, rseci2, rseci3, out bigr2x);
                initguess[caseopt] = bigr2x;

                rng1 = initguess[caseopt];  // old 12500 needs to be in km!! seems to do better when all the same? if too far off (*2) NAN
                rng2 = initguess[caseopt] * 1.02;  // 1.02 might be better? make the initial guess a bit different
                rng3 = initguess[caseopt] * 1.08;

                AstroLibr.anglesgooding(tdecl[idx1], tdecl[idx2], tdecl[idx3], trtasc[idx1], trtasc[idx2], trtasc[idx3],
                    jd[idx1], jdf[idx1], jd[idx2], jdf[idx2], jd[idx3], jdf[idx3],
                    rseci1, rseci2, rseci3, numhalfrev, rng1, rng2, rng3, out r2, out v2, out errstr);
                strbuildall.AppendLine(errstr);
                strbuildall.AppendLine("r2 " + r2[0].ToString("0.000000") + " " +
     r2[1].ToString("0.000000") + " " + r2[2].ToString("0.000000")
     + "v2 " + v2[0].ToString("0.000000") + " " +
     v2[1].ToString("0.000000") + " " + v2[2].ToString("0.000000"));
                AstroLibr.rv2coe(r2, v2, out p, out a, out ecc, out incl, out raan, out argp, out nu, out m, out arglat, out truelon, out lonper);
                strbuildall.AppendLine("gooding coes  a= " + a.ToString("0.0000") + " e= " + ecc.ToString("0.000000000") + " i= " + (incl * rad).ToString("0.0000") + " " +
                    (raan * rad).ToString("0.0000") + " " + (argp * rad).ToString("0.0000") + " " + (nu * rad).ToString("0.0000") + " " + (m * rad).ToString("0.0000") + " " +
                    (arglat * rad).ToString("0.0000")); // + " " + (truelon * rad).ToString("0.0000") + " " + (lonper * rad).ToString("0.0000"));
                strbuild.AppendLine("gooding coes  a= " + a.ToString("0.0000") + " e= " + ecc.ToString("0.000000000") + " i= " + (incl * rad).ToString("0.0000") + " " +
                    (raan * rad).ToString("0.0000") + " " + (argp * rad).ToString("0.0000") + " " + (nu * rad).ToString("0.0000") + " " + (m * rad).ToString("0.0000") + " " +
                    (arglat * rad).ToString("0.0000")); // + " " + (truelon * rad).ToString("0.0000") + " " + (lonper * rad).ToString("0.0000"));

                //                }  // loop through cases of caseopt = 0

            } // caseopt

            string directory = @"D:\Codes\LIBRARY\cs\TestAll\";
            this.strbuild.AppendLine("angles only tests case results written to " + directory + "testall-Angles.out ");
            this.strbuild.AppendLine(@"geo data for chap 9 plot written to D:\faabook\current\excel\testgeo.out for ch9 plot ");

            File.WriteAllText(directory + "testall-Angles.out", strbuildall.ToString());
            File.WriteAllText(directory + "testall-Anglessum.out", strbuild.ToString());

            this.opsStatus.Text = "Test Angles - Done";
            Refresh();

        }   // testangles





        public void testlambertumins()
        {
            double[,] tbiSu = new double[10, 3];
            double[,] tbiLu = new double[10, 3];
            double[] r1 = new double[3];
            double[] r2 = new double[3];
            double tmin, tminp, tminenergy;
            int i;

            r1 = new double[] { 2.500000 * AstroLibr.gravConst.re, 0.000000, 0.000000 };
            r2 = new double[] { 1.9151111 * AstroLibr.gravConst.re, 1.6069690 * AstroLibr.gravConst.re, 0.000000 };
            char dm = 'S';
            //char de = 'L';
            Int32 nrev = 0;

            // timing of routines
            var watch = System.Diagnostics.Stopwatch.StartNew();
            for (i = 0; i < 1000; i++)
            {
                for (int j = 1; j < 5; j++)
                    AstroLibr.lambertminT(r1, r2, dm, 'L', nrev, out tmin, out tminp, out tminenergy);

                for (int j = 1; j < 5; j++)
                    AstroLibr.lambertminT(r1, r2, dm, 'H', nrev, out tmin, out tminp, out tminenergy);
            }
            watch.Stop();
            var elapsedMs = watch.ElapsedMilliseconds;
            strbuild.AppendLine("time for Lambert minT " + watch.ElapsedMilliseconds);

            double kbi, tof;

            // timing of routines
            watch = System.Diagnostics.Stopwatch.StartNew();
            for (i = 0; i < 1000; i++)
            {
                AstroLibr.lambertumins(r1, r2, 1, 'S', out kbi, out tof);
                tbiSu[1, 1] = kbi;
                tbiSu[1, 2] = tof;
                AstroLibr.lambertumins(r1, r2, 2, 'S', out kbi, out tof);
                tbiSu[2, 1] = kbi;
                tbiSu[2, 2] = tof;
                AstroLibr.lambertumins(r1, r2, 3, 'S', out kbi, out tof);
                tbiSu[3, 1] = kbi;
                tbiSu[3, 2] = tof;
                AstroLibr.lambertumins(r1, r2, 4, 'S', out kbi, out tof);
                tbiSu[4, 1] = kbi;
                tbiSu[4, 2] = tof;
                AstroLibr.lambertumins(r1, r2, 5, 'S', out kbi, out tof);
                tbiSu[5, 1] = kbi;
                tbiSu[5, 2] = tof;

                AstroLibr.lambertumins(r1, r2, 1, 'L', out kbi, out tof);
                tbiLu[1, 1] = kbi;
                tbiLu[1, 2] = tof;
                AstroLibr.lambertumins(r1, r2, 2, 'L', out kbi, out tof);
                tbiLu[2, 1] = kbi;
                tbiLu[2, 2] = tof;
                AstroLibr.lambertumins(r1, r2, 3, 'L', out kbi, out tof);
                tbiLu[3, 1] = kbi;
                tbiLu[3, 2] = tof;
                AstroLibr.lambertumins(r1, r2, 4, 'L', out kbi, out tof);
                tbiLu[4, 1] = kbi;
                tbiLu[4, 2] = tof;
                AstroLibr.lambertumins(r1, r2, 5, 'L', out kbi, out tof);
                tbiLu[5, 1] = kbi;
                tbiLu[5, 2] = tof;
            }
            watch.Stop();
            elapsedMs = watch.ElapsedMilliseconds;
            strbuild.AppendLine("time for Lambert umin " + watch.ElapsedMilliseconds);

            double[,] tbidk = new double[10, 3];
            double[,] tbirk = new double[10, 3];
            double tusec = 806.8111238242922;
            double ootusec = 1.0 / tusec;

            // timing of routines
            watch = System.Diagnostics.Stopwatch.StartNew();

            double s, tau;
            AstroLambertkLibr.lambertkmins1st(r1, r2, out s, out tau);

            // for general cases, use 'x' for dm to get the tof/kbi values
            for (i = 0; i < 1000; i++)
            {
                AstroLambertkLibr.lambertkmins(s, tau, 1, 'x', 'L', out kbi, out tof);
                tbidk[1, 1] = kbi;
                tbidk[1, 2] = tof * ootusec;
                AstroLambertkLibr.lambertkmins(s, tau, 2, 'x', 'L', out kbi, out tof);
                tbidk[2, 1] = kbi;
                tbidk[2, 2] = tof * ootusec;
                AstroLambertkLibr.lambertkmins(s, tau, 3, 'x', 'L', out kbi, out tof);
                tbidk[3, 1] = kbi;
                tbidk[3, 2] = tof * ootusec;
                AstroLambertkLibr.lambertkmins(s, tau, 4, 'x', 'L', out kbi, out tof);
                tbidk[4, 1] = kbi;
                tbidk[4, 2] = tof * ootusec;
                AstroLambertkLibr.lambertkmins(s, tau, 5, 'x', 'L', out kbi, out tof);
                tbidk[5, 1] = kbi;
                tbidk[5, 2] = tof * ootusec;

                AstroLambertkLibr.lambertkmins(s, tau, 1, 'x', 'H', out kbi, out tof);
                tbirk[1, 1] = kbi;
                tbirk[1, 2] = tof * ootusec;
                AstroLambertkLibr.lambertkmins(s, tau, 2, 'x', 'H', out kbi, out tof);
                tbirk[2, 1] = kbi;
                tbirk[2, 2] = tof * ootusec;
                AstroLambertkLibr.lambertkmins(s, tau, 3, 'x', 'H', out kbi, out tof);
                tbirk[3, 1] = kbi;
                tbirk[3, 2] = tof * ootusec;
                AstroLambertkLibr.lambertkmins(s, tau, 4, 'x', 'H', out kbi, out tof);
                tbirk[4, 1] = kbi;
                tbirk[4, 2] = tof * ootusec;
                AstroLambertkLibr.lambertkmins(s, tau, 5, 'x', 'H', out kbi, out tof);
                tbirk[5, 1] = kbi;
                tbirk[5, 2] = tof * ootusec;
            }

            watch.Stop();
            elapsedMs = watch.ElapsedMilliseconds;
            strbuild.AppendLine("time for Lambert kmin " + watch.ElapsedMilliseconds);
        }

        public void testlambertminT()
        {
            double[] r1 = new double[3];
            double[] r2 = new double[3];
            char dm, de;
            Int32 nrev;
            double tmin, tminp, tminenergy;
            r1 = new double[] { 2.500000 * AstroLibr.gravConst.re, 0.000000, 0.000000 };
            r2 = new double[] { 1.9151111 * AstroLibr.gravConst.re, 1.6069690 * AstroLibr.gravConst.re, 0.000000 };
            dm = 'S';
            de = 'L';
            nrev = 0;

            AstroLibr.lambertminT(r1, r2, dm, de, nrev, out tmin, out tminp, out tminenergy);
            strbuild.AppendLine("lambertminT tmin  s " + tmin.ToString(fmt) + " minp " + tminp.ToString(fmt) +
                " minener " + tminenergy.ToString(fmt));

            AstroLibr.lambertminT(r1, r2, dm, 'H', nrev, out tmin, out tminp, out tminenergy);
            strbuild.AppendLine("lambertminT tmin  s " + tmin.ToString(fmt) + " minp " + tminp.ToString(fmt) +
                " minener " + tminenergy.ToString(fmt));
        }

        public void testlambhodograph()
        {
            double rad = 180.0 / Math.PI;
            double[] r1 = new double[3];
            double[] v1 = new double[3];
            double[] r2 = new double[3];
            double p, ecc, dnu, dtsec;
            double[] v1t;
            double[] v2t;
            r1 = new double[] { 2.500000 * AstroLibr.gravConst.re, 0.000000, 0.000000 };
            r2 = new double[] { 1.9151111 * AstroLibr.gravConst.re, 1.6069690 * AstroLibr.gravConst.re, 0.000000 };
            // assume circular initial orbit for vel calcs
            v1 = new double[] { 0.0, Math.Sqrt(AstroLibr.gravConst.mu / r1[0]), 0.0 };
            p = 12345.235;  // km
            ecc = 0.023487;
            dnu = 34.349128 / rad;
            dtsec = 92854.234;


            AstroLibr.lambhodograph(r1, r2, v1, p, ecc, dnu, dtsec, out v1t, out v2t);

            strbuild.AppendLine("lamb hod " + v1t[0].ToString(fmt) + " " + v1t[1].ToString(fmt) + " " + v1t[2].ToString(fmt) + " \nlamb hod" +
                   v2t[0].ToString(fmt) + " " + v2t[1].ToString(fmt) + " " + v2t[2].ToString(fmt));
        }

        public void testlambertbattin()
        {
            double[] r1 = new double[3];
            double[] v1 = new double[3];
            double[] r2 = new double[3];
            double[] v1t = new double[3];
            double[] v2t = new double[3];
            double dtsec, altpadc;
            string errorsum = "";
            string errorout = "";
            char hitearth, dm, de;
            Int32 nrev;

            r1 = new double[] { 2.500000 * AstroLibr.gravConst.re, 0.000000, 0.000000 };
            r2 = new double[] { 1.9151111 * AstroLibr.gravConst.re, 1.6069690 * AstroLibr.gravConst.re, 0.000000 };
            v1 = new double[] { 0.0, Math.Sqrt(AstroLibr.gravConst.mu / r1[0]), 0.0 };
            dm = 'S';
            de = 'L';
            nrev = 0;
            dtsec = 76.0 * 60.0;
            altpadc = 100.0 / AstroLibr.gravConst.re;  //er

            AstroLibr.lambertbattin(r1, r2, v1, dm, de, nrev, dtsec, altpadc, 'y', out v1t, out v2t, out hitearth, out errorsum, out errorout);

            strbuild.AppendLine("lambertbattin " + v1t[0].ToString(fmt) + " " + v1t[1].ToString(fmt) + " " + v1t[2].ToString(fmt) + " \nlambertbattin " +
                v2t[0].ToString(fmt) + " " + v2t[1].ToString(fmt) + " " + v2t[2].ToString(fmt));
        }

        public void testeq2rv()
        {
            double[] r = new double[3];
            double[] v = new double[3];
            double a, af, ag, chi, psi, meanlon;
            Int16 fr;
            a = 7236.346;
            af = 0.23457;
            ag = 0.47285;
            chi = 0.23475;
            psi = 0.28374;
            meanlon = 2.230482378;
            fr = 1;

            AstroLibr.eq2rv(a, af, ag, chi, psi, meanlon, fr, out r, out v);

            strbuild.AppendLine("eq2rv " + r[0].ToString(fmt) + " " + r[1].ToString(fmt) + " " + r[2].ToString(fmt) + " " +
                 v[0].ToString(fmt) + " " + v[1].ToString(fmt) + " " + v[2].ToString(fmt));
        }

        public void testrv2eq()
        {
            double[] r, v;
            double a, n, af, ag, chi, psi, meanlonM, meanlonNu;
            Int16 fr;

            r = new double[] { 2.500000 * AstroLibr.gravConst.re, 0.000000, 0.000000 };
            // assume circular initial orbit for vel calcs
            v = new double[] { 0.0, Math.Sqrt(AstroLibr.gravConst.mu / r[0]), 0.0 };

            AstroLibr.rv2eq(r, v, out a, out n, out af, out ag, out chi, out psi, out meanlonM, out meanlonNu, out fr);

            strbuild.AppendLine("rv2eq   a " + a.ToString(fmt) + " n " + n.ToString(fmt) + " af " + af.ToString(fmt) + " ag "
                + ag.ToString(fmt) + " chi " + chi.ToString(fmt) + " psi " + psi.ToString(fmt) + " mm " +
                meanlonM.ToString(fmt) + " mnu " + meanlonNu);
        }



        // test building the lambert envelope
        private void testAll()
        {
            // done in lambert form

        }


        /* ------------------------------------------------------------------------------
        *                                 testAllMoving 
        *                                    
        * calc the values needed to graph the whole envelope response. this is the most 
        * generic construction for the lambert solutions.
        *  you need to entere which lambert technique to use
        *      output these results separately to the matlab directory
        * 
        *  author        : david vallado             davallado@gmail.com  10 oct 2019
        *
        *  inputs        description                                   range / units
        *
        *  outputs       :
        *    
        *  locals        :
        *    r1          - ijk position vector 1                km
        *    r2          - ijk position vector 2                km
        *    v1          - ijk velocity vector 1 if avail       km/s
        *    dm          - direction of motion                  'L', 'S'
        *    de          - orbital energy                       'L', 'H'
        *    dtsec       - time between r1 and r2               second
        *    dtwait      - time to wait before starting         second
        *    nrev        - number of revs to complete           0, 1, 2, 3,  
        *    tbi         - array of times for nrev              [,] 
        *    altpad      - altitude pad for hitearth calc       km
        *    v1t         - ijk transfer velocity vector         km/s
        *    v2t         - ijk transfer velocity vector         km/s
        *    hitearth    - flag if hit or not                   'y', 'n'
        *    error       - error flag                           1, 2, 3,   use numbers since c++ is so horrible at strings
        *  
          ------------------------------------------------------------------------------*/

        public void testAllMoving()
        {
            double dtwait, dtsec, dtseco;
            string detailSum, detailAll, errorout;
            char dm, de, hitearth, whichcase;
            Int32 ktr, ktr1, ktr2, ktr3, ktr4, i, iktr, nrev, numiter;
            double f, g, gdot;
            double[] v1t = new double[3];
            double[] v2t = new double[3];
            double tofsh, kbish, toflg, kbilg, toflo, kbilo, tofhi, kbihi;
            string outstr;
            double[] r1 = new double[3];
            double[] r2 = new double[3];
            double[] v1 = new double[3];
            double[] v2 = new double[3];
            double[] r1o = new double[3];
            double[] r2o = new double[3];
            double[] v1o = new double[3];
            double[] v2o = new double[3];
            double[] dv1 = new double[3];
            double[] dv2 = new double[3];
            StringBuilder strbuild = new StringBuilder();
            StringBuilder strbuildDV = new StringBuilder();

            // initialize variables
            detailSum = "";
            detailAll = "";
            outstr = "";
            char show = 'y';     // for test180, show = n, show180 = y
            //char show180 = 'n';  // for testlamb known show = y, show180 = n, n/n for envelope
            tofsh = 0.0;
            kbish = 0.0;
            toflg = 0.0;
            kbilg = 0.0;
            toflo = 0.0;
            kbilo = 0.0;
            tofhi = 0.0;
            kbihi = 0.0;

            whichcase = 'k';
            //    whichcase = 'u';
            //    whichcase = 'b';

            dm = 'S';
            de = 'H';
            double altpadc = 200.0 / 6378.137;  // set 200 km for altitude you set as the over limit. 
            ktr1 = 0;
            ktr2 = 0;
            ktr3 = 0;
            ktr4 = 0;

            //mu = 3.986004415e5;
            //tusec = 806.8111238242922;
            numiter = 17;
            dtwait = 0.0;

            this.opsStatus.Text = "working on all cases ";
            Refresh();

            // book fig
            //   r1o = new double[] { 2.500000 * 6378.137, 0.000000, 0.000000 };
            //   r2o = new double[] { 1.9151111 * 6378.137, 1.6069690 * 6378.137, 0.000000 };

            // do reverse case to check retro direct!
            //r2 = new double[] { 2.500000 * 6378.137, 0.000000, 0.000000 };
            //r1 = new double[] { 1.9151111 * 6378.137, 1.6069690 * 6378.137, 0.000000 };

            // assume circular initial orbit for vel calcs
            //    v1o = new double[] { 0.0, Math.Sqrt(mu / r1o[0]), 0.0 };
            //    double ang = Math.Atan(r2o[1] / r2o[0]);
            //    v2o = new double[] { -Math.Sqrt(mu / r2o[1]) * Math.Cos(ang), Math.Sqrt(mu / r2o[0]) * Math.Sin(ang), 0.0 };

            // book fig 7-17 case tgt/int fixed
            r1o = new double[] { -6518.1083, -2403.8479, -22.1722 };
            v1o = new double[] { 2.604057, -7.105717, -0.263218 };
            r2o = new double[] { 6697.4756, 1794.5832, 0.000 };
            v2o = new double[] { -1.962373, 7.323674, 0.000 };

            // book fig 7-18 case tgt/int moving
            r1o = new double[] { -6175.1034, 2757.0706, 1626.6556 };
            v1o = new double[] { 2.376641, 1.139677, 7.078097 };
            r2o = new double[] { -1078.007289, 8796.641859, 1890.7135 };
            v2o = new double[] { 2.654700, 1.018600, 7.015400 };


            // nathaniel test case (1 rev)
            r1o = new double[] { 7117.5156243154161, -4257.1942597683246, -583.88210887807986 };
            v1o = new double[] { -7.0417929290503141, -2.8681303717265290, 1.2224374606557487 };
            r2o = new double[] { 7172.3074180808417, -4123.2685183007470, -560.47505356742181 };
            v2o = new double[] { 5.9109725501309471, 2.0990367313315055, -1.1164901518784618 };  // xx

            //double tmin, tminp, tminenergy;
            //AstroLibr.lambertminT(r1o, r2o, 'S', 'L', 0, out tmin, out tminp, out tminenergy);
            //AstroLibr.lambertminT(r1o, r2o, 'S', 'L', 1, out tmin, out tminp, out tminenergy);

            //double tmaxrp;
            //AstroLibr.lambertTmaxrp(r1o, r2o, 'S', 0, out tmaxrp, out v1t);
            //AstroLibr.lambertTmaxrp(r1o, r2o, 'S', 2, out tmaxrp, out v1t);

            //char opt = 'f';  // fixed target through dtsec
            char opt = 'm';  // moving target through dtsec

            // be sure to set the correct S/L, L/H below as needed

            // dtwait
            int lktr1 = 250;  // 250
            double step1 = 120.0;  // 120
            // dtsec
            int lktr2 = 100;  // 100
            double step2 = 120.0;  // 120

            // do a loop for dtwait
            for (int loopktr = 0; loopktr < lktr1; loopktr++)  // 0
            {
                dtwait = loopktr * step1;   // secs

                this.opsStatus.Text = "working dtwait = " + dtwait.ToString();
                Refresh();

                // propagate through dtwait
                AstroLibr.kepler(r1o, v1o, dtwait, out r1, out v1);

                hitearth = '-';

                // -------------- do 0 rev cases first
                nrev = 0;
                ktr = 0;  // overall ktr
                for (iktr = 1; iktr <= 1; iktr++)   // 1    2
                {
                    if (iktr == 1)
                    {
                        dm = 'S';
                        de = 'L';   // or 'H'
                    }
                    if (iktr == 2)
                    {
                        dm = 'L';
                        de = 'H';    // or 'L'
                    }

                    dtseco = 0.0;
                    // calc the actual lambert values
                    for (i = 1; i <= lktr2; i++)
                    {
                        dtsec = dtseco + i * step2;

                        // propagate through dtsec and dtwait
                        if (opt.Equals('m'))
                            AstroLibr.kepler(r2o, v2o, dtwait + dtsec, out r2, out v2);
                        else
                            AstroLibr.kepler(r2o, v2o, dtwait, out r2, out v2);

                        if (nrev > 0)
                            getmins(1, 'u', nrev, r1, r2, 0.0, 0.0, dm, de, out tofsh, out kbish, out toflg, out kbilg, out outstr);
                        // strbuild.AppendLine(outstr);
                        hitearth = '-';

                        if (de == 'L')
                        {
                            if (whichcase == 'k')
                            {
                                AstroLambertkLibr.lambertK(r1, v1, r2, dm, de, nrev, dtwait, dtsec, toflo, kbilo, numiter, altpadc, 'n', show,
                                     out v1t, out v2t, out f, out g, out gdot, out hitearth, out errorout, out detailSum, out detailAll);
                                // writeout details if there's a problem...probably not needed
                                //if (detailAll.Contains("gnot"))
                                //{
                                //    AstroLambertkLibr.lambertK(r1, v1, r2, dm, de, nrev, dtwait, dtsec, tbilk, numiter, altpadc,
                                //         out v1t, out v2t, out f, out g, out gdot, out hitearth, out errorout, out detailSum, out detailAll, 'y', show180);
                                //}
                            }
                            if (whichcase == 'b')
                                AstroLibr.lambertbattin(r1, r2, v1, dm, de, nrev, dtsec, altpadc, 'y', out v1t, out v2t, out hitearth, out detailSum, out detailAll);
                            if (whichcase == 'u')
                            {
                                if (dm == 'S')
                                    AstroLibr.lambertuniv(r1, r2, v1, dm, de, nrev, dtsec, kbish, altpadc, 'y', out v1t, out v2t, out hitearth, out detailSum, out detailAll);
                                else
                                    AstroLibr.lambertuniv(r1, r2, v1, dm, de, nrev, dtsec, kbilg, altpadc, 'y', out v1t, out v2t, out hitearth, out detailSum, out detailAll);
                            }
                        }
                        else
                        {
                            if (whichcase == 'k')
                            {
                                AstroLambertkLibr.lambertK(r1, v1, r2, dm, de, nrev, dtwait, dtsec, tofhi, kbihi, numiter, altpadc, 'n', show,
                                     out v1t, out v2t, out f, out g, out gdot, out hitearth, out errorout, out detailSum, out detailAll);
                                //if (detailAll.Contains("gnot"))
                                //{
                                //    AstroLambertkLibr.lambertK(r1, v1, r2, dm, de, nrev, dtwait, dtsec, tbihk, numiter, altpadc,
                                //         out v1t, out v2t, out f, out g, out gdot, out hitearth, out errorout, out detailSum, out detailAll, 'y', show180);
                                //}
                            }
                            if (whichcase == 'b')
                                AstroLibr.lambertbattin(r1, r2, v1, dm, de, nrev, dtsec, altpadc, 'y', out v1t, out v2t, out hitearth, out detailSum, out detailAll);
                            if (whichcase == 'u')
                            {
                                if (dm == 'S')
                                    AstroLibr.lambertuniv(r1, r2, v1, dm, de, nrev, dtsec, kbish, altpadc, 'y', out v1t, out v2t, out hitearth, out detailSum, out detailAll);
                                else
                                    AstroLibr.lambertuniv(r1, r2, v1, dm, de, nrev, dtsec, kbilg, altpadc, 'y', out v1t, out v2t, out hitearth, out detailSum, out detailAll);
                            }
                        }
                        ktr = ktr + 1;
                        //strbuild.AppendLine(detailSum);
                        dv1 = new double[] { 0.0, 0.0, 0.0 };
                        dv2 = new double[] { 0.0, 0.0, 0.0 };
                        if (MathTimeLibr.mag(v1t) > 0.0000001)
                        {
                            for (int ii = 0; ii < 3; ii++)
                            {
                                dv1[ii] = v1t[ii] - v1[ii];
                                dv2[ii] = v2t[ii] - v2[ii];
                            }
                        }
                        strbuild.AppendLine(detailSum + " " + MathTimeLibr.mag(dv1).ToString() + " " + MathTimeLibr.mag(dv2).ToString());
                        double magdv1 = MathTimeLibr.mag(dv1);
                        double magdv2 = MathTimeLibr.mag(dv2);
                        strbuildDV.AppendLine(dtwait.ToString("0.0000000").PadLeft(12) + " " +
                              dtsec.ToString("0.0000000").PadLeft(15) + " " +
                              magdv1.ToString() + " " + magdv2.ToString());  // + " " + dm + "  " + de);
                                                                             //   strbuild.AppendLine(detailAll);
                    }  // for i through all the times

                    strbuild.AppendLine(" ");
                    strbuildDV.AppendLine(" ");
                    if (iktr == 1)
                        ktr1 = ktr;
                    if (iktr == 2)
                        ktr2 = ktr;
                } // for iktr through cases

                //strbuild.AppendLine(" ");
                this.opsStatus.Text = "working dtwait = " + dtwait.ToString() + " now the 4 cases";
                Refresh();

                // set this for doing just the nrev cases
                getmins(1, 'u', nrev, r1, r2, 0.0, 0.0, dm, de, out tofsh, out kbish, out toflg, out kbilg, out outstr);

                for (iktr = 5; iktr <= 1; iktr++)  //4
                {
                    if (iktr == 1)
                    {
                        dm = 'S';
                        de = 'L';
                    }
                    if (iktr == 2)
                    {
                        dm = 'S';
                        de = 'H';
                    }
                    if (iktr == 3)
                    {
                        dm = 'L';
                        de = 'L';
                    }
                    if (iktr == 4)
                    {
                        dm = 'L';
                        de = 'H';
                    }

                    for (nrev = 1; nrev <= 1; nrev++)   // 1  4
                    {
                        dtseco = 0.0;
                        if (nrev > 0)
                        {
                            // probably need to use kmins...
                            if (whichcase == 'k')
                            {
                                if (de == 'H')
                                    dtseco = toflo;// * tusec; // use univ for test
                                else
                                    dtseco = tofhi;// * tusec;  // use univ for test
                            }
                            else
                            {
                                // use univ for all these cases (univ and battin)
                                if (dm == 'L')
                                    dtseco = toflg;
                                else
                                    dtseco = tofsh;
                            }
                        }

                        // calc the actual lambert values
                        for (i = 1; i <= lktr2; i++)
                        {
                            dtsec = dtseco + i * step2;

                            // propagate through dtsec and dtwait
                            if (opt.Equals('m'))
                                AstroLibr.kepler(r2o, v2o, dtwait + dtsec, out r2, out v2);
                            else
                                AstroLibr.kepler(r2o, v2o, dtwait, out r2, out v2);

                            getmins(1, 'u', nrev, r1, r2, 0.0, 0.0, dm, de, out tofsh, out kbish, out toflg, out kbilg, out outstr);
                            // strbuild.AppendLine(outstr);
                            hitearth = '-';

                            if (de == 'L')
                            {
                                if (whichcase == 'k')
                                {
                                    AstroLambertkLibr.lambertK(r1, v1, r2, dm, de, nrev, dtwait, dtsec, toflo, kbilo, numiter, altpadc, 'n', show,
                                         out v1t, out v2t, out f, out g, out gdot, out hitearth, out errorout, out detailSum, out detailAll);
                                    //if (detailAll.Contains("gnot"))
                                    //{
                                    //    AstroLambertkLibr.lambertK(r1, v1, r2, dm, de, nrev, dtwait, dtsec, tbilk, numiter, altpadc,
                                    //         out v1t, out v2t, out f, out g, out gdot, out hitearth, out errorout, out detailSum, out detailAll, 'y', show180);
                                    //}
                                }
                                if (whichcase == 'b')
                                    AstroLibr.lambertbattin(r1, r2, v1, dm, de, nrev, dtsec, altpadc, 'y', out v1t, out v2t, out hitearth, out detailSum, out detailAll);
                                if (whichcase == 'u')
                                {
                                    if (dm == 'S')
                                        AstroLibr.lambertuniv(r1, r2, v1, dm, de, nrev, dtsec, kbish, altpadc, 'y', out v1t, out v2t, out hitearth, out detailSum, out detailAll);
                                    else
                                        AstroLibr.lambertuniv(r1, r2, v1, dm, de, nrev, dtsec, kbilg, altpadc, 'y', out v1t, out v2t, out hitearth, out detailSum, out detailAll);
                                }
                            }
                            else
                            {
                                if (whichcase == 'k')
                                {
                                    AstroLambertkLibr.lambertK(r1, v1, r2, dm, de, nrev, dtwait, dtsec, tofhi, kbihi, numiter, altpadc, 'n', show,
                                         out v1t, out v2t, out f, out g, out gdot, out hitearth, out errorout, out detailSum, out detailAll);
                                    //if (detailAll.Contains("gnot"))
                                    //{
                                    //    AstroLambertkLibr.lambertK(r1, v1, r2, dm, de, nrev, dtwait, dtsec, tbihk, numiter, altpadc,
                                    //         out v1t, out v2t, out f, out g, out gdot, out hitearth, out errorout, out detailSum, out detailAll, 'y', show180);
                                    //}
                                }
                                if (whichcase == 'b')
                                    AstroLibr.lambertbattin(r1, r2, v1, dm, de, nrev, dtsec, altpadc, 'y', out v1t, out v2t, out hitearth, out detailSum, out detailAll);
                                if (whichcase == 'u')
                                {
                                    if (dm == 'S')
                                        AstroLibr.lambertuniv(r1, r2, v1, dm, de, nrev, dtsec, kbish, altpadc, 'y', out v1t, out v2t, out hitearth, out detailSum, out detailAll);
                                    else
                                        AstroLibr.lambertuniv(r1, r2, v1, dm, de, nrev, dtsec, kbilg, altpadc, 'y', out v1t, out v2t, out hitearth, out detailSum, out detailAll);
                                }
                            }
                            ktr = ktr + 1;
                            //strbuild.AppendLine(detailSum);
                            dv1 = new double[] { 0.0, 0.0, 0.0 };
                            dv2 = new double[] { 0.0, 0.0, 0.0 };
                            if (MathTimeLibr.mag(v1t) > 0.0000001)
                            {
                                for (int ii = 0; ii < 3; ii++)
                                {
                                    dv1[ii] = v1t[ii] - v1[ii];
                                    dv2[ii] = v2t[ii] - v2[ii];
                                }
                            }
                            strbuild.AppendLine(detailSum + " " + MathTimeLibr.mag(dv1).ToString() + " " + MathTimeLibr.mag(dv2).ToString());
                            double magdv1 = MathTimeLibr.mag(dv1);
                            double magdv2 = MathTimeLibr.mag(dv2);
                            strbuildDV.AppendLine(dtwait.ToString("0.0000000").PadLeft(12) + " " +
                                  dtsec.ToString("0.0000000").PadLeft(15) + " " +
                                  magdv1.ToString() + " " + magdv2.ToString());  // + " " + dm + "  " + de);
                        }  // for i through all the times

                        strbuild.AppendLine(" ");
                        strbuildDV.AppendLine(" ");
                    }  // if nrev > 0

                    //                    strbuild.AppendLine(" ");
                    if (iktr == 2)
                        ktr3 = ktr;
                    if (iktr == 4)
                        ktr4 = ktr;
                } // for iktr through cases


            } // loop through dtwait

            strbuild.AppendLine("ktrs " + ktr1.ToString() + " " + ktr2.ToString() + " " + ktr3.ToString() + " " + ktr4.ToString());

            string directory = @"d:\codes\library\matlab\";
            File.WriteAllText(directory + "tlambertAllx.out", strbuild.ToString());
            File.WriteAllText(directory + "tlamb3dx.out", strbuildDV.ToString());

            this.opsStatus.Text = "Done ";
            Refresh();
        }  // testAllMoving


        /* ------------------------------------------------------------------------------
        *                                    getmins 
        *                                    
        * find tbi mins for k, lambert, etc. does either universal and k approach.
        *   universal seems to find them better though... could store in an array, but 
        *   faster to do 1 at a time. also note that the k lambert needs SH - LL, and LL - SH 
        *   reversal and it is done here. 
        *                                    
        *  author        : david vallado             davallado@gmail.com  10 oct 2019
        *
        *  inputs        description                                   range / units
        *    loopktr     - ktr for whether or not to write output           0 (writes)   
        *    app         - which approach to use                           'u', 'k'
        *    nrev        - number of revolultions 
        *    r1          - first position vector                            km
        *    r2          - second position vector                           km
        *    s           - parameter for k only, not needed for univ
        *    tau         - parameter for k only, not needed for univ
        *    dm          - parameter for k only, not needed for univ        'S', 'L'
        *    de          - parameter for k only, not needed for univ        'L', 'H'
        *    
        *  outputs       :
        *    tof         - min tof for the specified nrev                    s
        *    kbi         - min psi/k/etc for the given nrev 
        *    outstr      - output string if case 0
        *
        *  locals :
        *
        *  coupling      :
        *
         ------------------------------------------------------------------------------*/

        private void getmins
            (
                int loopktr, char app, int nrev, double[] r1, double[] r2, double s, double tau, char dm, char de,
                out double tof1, out double kbi1, out double tof2, out double kbi2, out string outstr
            )
        {
            StringBuilder strbuild = new StringBuilder();
            string detailSum;
            double tofsh, toflg, kbish, kbilg;
            double[] v1t = new double[3];
            tof1 = 0.0;
            kbi1 = 0.0;
            tof2 = 0.0;
            kbi2 = 0.0;

            //tusec = 806.8111238242922;
            if (nrev > 0)
            {

                if (app == 'u')
                {
                    // universal variable approach 
                    AstroLibr.lambertumins(r1, r2, nrev, 'S', out kbi1, out tof1);
                    AstroLibr.lambertumins(r1, r2, nrev, 'L', out kbi2, out tof2);
                }
                else
                {
                    // -----------do these calcs one time to save time
                    // call this outside getmins
                    // AstroLambertkLibr.lambertkmins1st(r1, r2, out s, out tau);

                    // k value approaches
                    AstroLambertkLibr.lambertkmins(s, tau, nrev, 'x', 'L', out kbi1, out tof1);
                    //tof1 = tof1 / tusec;  // in tu, tof is in secs
                    AstroLambertkLibr.lambertkmins(s, tau, nrev, 'x', 'H', out kbi2, out tof2);
                    //tof2 = tof2 / tusec;
                    // switch these here so it's not needed elsewhere
                    // no, switch in kmins
                    //if ((dm == 'S' && de == 'H') || ((dm == 'L' && de == 'L')))
                    //{
                    //    temp = kbi1;
                    //    kbi1 = kbi2;
                    //    kbi2 = temp;
                    //    temp = tof1;
                    //    tof1 = tof2;
                    //    tof2 = temp;
                    //}
                }

                // -----------------------------put min values etc points on plot ---------------------------------
                // writeout just one time
                if (loopktr == 0)
                {
                    strbuild.AppendLine("Lambertumin");
                    detailSum = "S   L   1  0.000 " + tof1.ToString("0.#######").PadLeft(15) + " psimin " + kbi1.ToString("0.#######").PadLeft(15) + " -0\n";
                    strbuild.AppendLine(detailSum);

                    strbuild.AppendLine("Lambertkmin");
                    detailSum = "S   L   1  0.000 " + (tof1).ToString("0.#######").PadLeft(15) + " kmin " + kbi1.ToString("0.#######").PadLeft(15) + " -0\n";
                    strbuild.AppendLine(detailSum);

                    // -----------------------------put min values etc points on plot ---------------------------------
                    double tmin, tminp, tminenergy;
                    strbuild.AppendLine("Lamberttmin");
                    AstroLibr.lambertminT(r1, r2, 'S', 'L', 1, out tmin, out tminp, out tminenergy);
                    detailSum = "S   L   tminp  0.000 " + tminp.ToString("0.#######").PadLeft(15) + "5.0000000".PadLeft(15) + " -0";
                    strbuild.AppendLine(detailSum);
                    detailSum = "S   L   tminp  0.000 " + tminp.ToString("0.#######").PadLeft(15) + "-5.0000000".PadLeft(15) + " -0\n";
                    strbuild.AppendLine(detailSum);

                    detailSum = "S   L   tminenergy  0.000 " + tminenergy.ToString("0.#######").PadLeft(15) + "5.0000000".PadLeft(15) + " -0";
                    strbuild.AppendLine(detailSum);
                    detailSum = "S   L   tminenergy  0.000 " + tminenergy.ToString("0.#######").PadLeft(15) + "-5.0000000".PadLeft(15) + " -0\n";
                    strbuild.AppendLine(detailSum);

                    AstroLibr.lambertumins(r1, r2, 1, 'S', out kbish, out tofsh);
                    detailSum = "S   L   1  0.000 " + tmin.ToString("0.#######").PadLeft(15) + (tofsh - 5).ToString("0.#######").PadLeft(15) + " -0";
                    strbuild.AppendLine(detailSum);
                    detailSum = "S   L   1  0.000 " + tmin.ToString("0.#######").PadLeft(15) + (tofsh + 5).ToString("0.#######").PadLeft(15) + " -0\n";
                    strbuild.AppendLine(detailSum);
                    AstroLibr.lambertumins(r1, r2, 2, 'S', out kbish, out tofsh);
                    AstroLibr.lambertminT(r1, r2, 'S', 'L', 2, out tmin, out tminp, out tminenergy);
                    detailSum = "S   L   2  0.000 " + tmin.ToString("0.#######").PadLeft(15) + (tofsh - 5).ToString("0.#######").PadLeft(15) + " -0";
                    strbuild.AppendLine(detailSum);
                    detailSum = "S   L   2  0.000 " + tmin.ToString("0.#######").PadLeft(15) + (tofsh + 5).ToString("0.#######").PadLeft(15) + " -0\n";
                    strbuild.AppendLine(detailSum);
                    AstroLibr.lambertumins(r1, r2, 3, 'S', out kbish, out tofsh);
                    AstroLibr.lambertminT(r1, r2, 'S', 'L', 3, out tmin, out tminp, out tminenergy);
                    detailSum = "S   L   3  0.000 " + tmin.ToString("0.#######").PadLeft(15) + (tofsh - 5).ToString("0.#######").PadLeft(15) + " -0";
                    strbuild.AppendLine(detailSum);
                    detailSum = "S   L   3  0.000 " + tmin.ToString("0.#######").PadLeft(15) + (tofsh + 5).ToString("0.#######").PadLeft(15) + " -0\n";
                    strbuild.AppendLine(detailSum);

                    AstroLibr.lambertminT(r1, r2, 'L', 'H', 1, out tmin, out tminp, out tminenergy);
                    detailSum = "L   H   tminp  0.000 " + tminp.ToString("0.#######").PadLeft(15) + "5.0000000".PadLeft(15) + " -0";
                    strbuild.AppendLine(detailSum);
                    detailSum = "L   H   tminp  0.000 " + tminp.ToString("0.#######").PadLeft(15) + "-5.0000000".PadLeft(15) + " -0\n";
                    strbuild.AppendLine(detailSum);

                    detailSum = "L   H   tminenergy  0.000 " + tminenergy.ToString("0.#######").PadLeft(15) + "5.0000000".PadLeft(15) + " -0";
                    strbuild.AppendLine(detailSum);
                    detailSum = "L   H   tminenergy  0.000 " + tminenergy.ToString("0.#######").PadLeft(15) + "-5.0000000".PadLeft(15) + " -0\n";
                    strbuild.AppendLine(detailSum);

                    AstroLibr.lambertumins(r1, r2, 1, 'L', out kbilg, out toflg);
                    detailSum = "L   H   1  0.000 " + tmin.ToString("0.#######").PadLeft(15) + (toflg - 5).ToString("0.#######").PadLeft(15) + " -0";
                    strbuild.AppendLine(detailSum);
                    detailSum = "L   H   1  0.000 " + tmin.ToString("0.#######").PadLeft(15) + (toflg + 5).ToString("0.#######").PadLeft(15) + " -0\n";
                    strbuild.AppendLine(detailSum);
                    AstroLibr.lambertumins(r1, r2, 2, 'L', out kbilg, out toflg);
                    AstroLibr.lambertminT(r1, r2, 'S', 'H', 2, out tmin, out tminp, out tminenergy);
                    detailSum = "L   H   2  0.000 " + tmin.ToString("0.#######").PadLeft(15) + (toflg - 5).ToString("0.#######").PadLeft(15) + " -0";
                    strbuild.AppendLine(detailSum);
                    detailSum = "L   H   2  0.000 " + tmin.ToString("0.#######").PadLeft(15) + (toflg + 5).ToString("0.#######").PadLeft(15) + " -0\n";
                    strbuild.AppendLine(detailSum);
                    AstroLibr.lambertumins(r1, r2, 3, 'L', out kbilg, out toflg);
                    AstroLibr.lambertminT(r1, r2, 'S', 'H', 3, out tmin, out tminp, out tminenergy);
                    detailSum = "L   H   3  0.000 " + tmin.ToString("0.#######").PadLeft(15) + (toflg - 5).ToString("0.#######").PadLeft(15) + " -0";
                    strbuild.AppendLine(detailSum);
                    detailSum = "L   H   3  0.000 " + tmin.ToString("0.#######").PadLeft(15) + (toflg + 5).ToString("0.#######").PadLeft(15) + " -0\n";
                    strbuild.AppendLine(detailSum);

                    // find max rp values for each nrev
                    double tmaxrp;
                    AstroLibr.lambertTmaxrp(r1, r2, 'S', 0, out tmaxrp, out v1t);
                    strbuild.AppendLine("x   x   0  0.000 " + tmaxrp.ToString() + " " +
                        v1t[0].ToString("0.0000000").PadLeft(15) + v1t[1].ToString("0.0000000").PadLeft(15) + v1t[2].ToString("0.0000000").PadLeft(15));
                    strbuild.AppendLine("x   x   0  0.000 " + tmaxrp.ToString() + " " +
                        v1t[0].ToString("0.0000000").PadLeft(15) + v1t[1].ToString("0.0000000").PadLeft(15) + v1t[2].ToString("0.0000000").PadLeft(15) + "\n");
                    AstroLibr.lambertTmaxrp(r1, r2, 'S', 1, out tmaxrp, out v1t);
                    strbuild.AppendLine("x   x   1  0.000 " + tmaxrp.ToString() + " " +
                        v1t[0].ToString("0.0000000").PadLeft(15) + v1t[1].ToString("0.0000000").PadLeft(15) + v1t[2].ToString("0.0000000").PadLeft(15));
                    strbuild.AppendLine("x   x   1  0.000 " + tmaxrp.ToString() + " " +
                        v1t[0].ToString("0.0000000").PadLeft(15) + v1t[1].ToString("0.0000000").PadLeft(15) + v1t[2].ToString("0.0000000").PadLeft(15) + "\n");
                    AstroLibr.lambertTmaxrp(r1, r2, 'S', 2, out tmaxrp, out v1t);
                    strbuild.AppendLine("x   x   2  0.000 " + tmaxrp.ToString() + " " +
                        v1t[0].ToString("0.0000000").PadLeft(15) + v1t[1].ToString("0.0000000").PadLeft(15) + v1t[2].ToString("0.0000000").PadLeft(15));
                    strbuild.AppendLine("x   x   2  0.000 " + tmaxrp.ToString() + " " +
                        v1t[0].ToString("0.0000000").PadLeft(15) + v1t[1].ToString("0.0000000").PadLeft(15) + v1t[2].ToString("0.0000000").PadLeft(15) + "\n");
                }
            }  // if nrev > 0

            outstr = strbuild.ToString();
        }   // getmins


        /* ------------------------------------------------------------------------------
        *                                      makesurf 
        *                                    
        *  make a surface from a fixdat result with numbers takes a text file of number of points, 
        *  then all the points, and cross-hatches it to get a surface. you run fixdat first 
        *  in most cases.
        *                                    
        *  author        : david vallado             davallado@gmail.com  10 oct 2019
        *
        *  inputs        description                                   range / units
        *    infilename  - in filename  
        *    outfilename - out filename  
        *    
        *  outputs       :
        *
        *  locals :
        *
        *  coupling      :
        *
          ------------------------------------------------------------------------------*/

        private void makesurf
            (
        string infilename,
        string outfilename
            )
        {
            Int32 ktr, numPerLine, NumLines, i, j;
            string line, line1, Restoflgine;
            string[] linesplt;
            StringBuilder strbuild = new StringBuilder();
            Restoflgine = "";

            string[] fileData = File.ReadAllLines(infilename);

            // process all the x lines 
            numPerLine = 0;
            ktr = 0;     // reset the file
            NumLines = 0;
            while (ktr < fileData.Count() - 1)  // not eof
            {
                line = fileData[ktr];
                linesplt = line.Split(' ');
                numPerLine = Convert.ToInt32(linesplt[0]);
                // matlab uses Inf or Nan to start a new line
                // needs to be in each col as well
                strbuild.AppendLine(" Nan Nan NaN NaN NaN NaN");  // numPerLine.ToString() +
                ktr = ktr + 1;
                NumLines = NumLines + 1;

                for (i = 1; i <= numPerLine; i++)
                {
                    line = fileData[ktr];
                    line1 = Regex.Replace(line, @"\s+", " ");
                    linesplt = line1.Split(' ');
                    //int posrest = line1.IndexOf(linesplt[2], 15); // start at position 3
                    //Restoflgine = line1.Substring(posrest -1, line1.Length -posrest);
                    Restoflgine = linesplt[2].ToString() + " " + linesplt[3].ToString() + " " +
                                 linesplt[4].ToString() + " " + linesplt[5].ToString();
                    strbuild.AppendLine(linesplt[0].ToString() + " " + linesplt[1].ToString() + " " + Restoflgine);
                    ktr = ktr + 1;
                }
            }

            // ------process the y lines-------
            // 'process y lines ');
            // the number of lines needs to be constant!! 
            //numPerLine = 0;
            int numinrow = 0;  // position of each y line
            // go through each point in initial line
            while (numinrow < numPerLine)
            {
                this.opsStatus.Text = "Done with line " + numinrow.ToString();
                Refresh();

                // ---get the nth point from the first row---
                ktr = 0;  // reset the file
                //line = fileData[ktr];
                //linesplt = line.Split(' ');
                //k = Convert.ToInt32(linesplt[0]);
                strbuild.AppendLine(" Nan Nan NaN NaN NaN NaN");  // k.ToString() +
                ktr = ktr + 1 + numinrow;  // get to first point of data
                line = fileData[ktr];
                line1 = Regex.Replace(line, @"\s+", " ");
                linesplt = line1.Split(' ');
                Restoflgine = linesplt[2].ToString() + " " + linesplt[3].ToString() + " " +
                             linesplt[4].ToString() + " " + linesplt[5].ToString();
                strbuild.AppendLine(linesplt[0].ToString() + " " + linesplt[1].ToString() + " " + Restoflgine);

                // ---get nth number from each other segment---
                // since they are all evenly spaced, simply add the delta until the end of file
                Int32 ktr0 = ktr + 1;
                for (j = 1; j < NumLines; j++)
                {
                    ktr = ktr0 + j * numPerLine;
                    line = fileData[ktr];
                    line1 = Regex.Replace(line, @"\s+", " ");
                    linesplt = line1.Split(' ');
                    Restoflgine = linesplt[2].ToString() + " " + linesplt[3].ToString() + " " +
                                 linesplt[4].ToString() + " " + linesplt[5].ToString();
                    strbuild.AppendLine(linesplt[0].ToString() + " " + linesplt[1].ToString() + " " + Restoflgine);
                    ktr0 = ktr0 + 1;
                }

                numinrow = numinrow + 1;
            }  // while

            string directory = @"d:\codes\library\matlab\";
            File.WriteAllText(directory + "surf.out", strbuild.ToString());
        }  // makesurf



        /* ------------------------------------------------------------------------------
        *                                      fixdat 
        *                                    
        *  fix the blank lines in a datafile. let 4 values be taken depending on the 
        *  intindex values inserts the number of points for each segment, then it can be 
        *  used in makesurf.
        *                                    
        *  author        : david vallado             davallado@gmail.com  10 oct 2019
        *
        *  inputs        description                                   range / units
        *    infilename  - in filename  
        *    outfilename - out filename  
        *    intindxx    - which indices to use, 1.2.3.4.5, 2.5.7.8 etc
        *    
        *  outputs       :
        *
        *  locals :
        *
          ------------------------------------------------------------------------------*/

        private void fixdat
            (
            string infilename,
            string outfilename,
            int intindx1, int intindx2, int intindx3, int intindx4, int intindx5, int intindx6
            )
        {
            Int32 ktr, i, j;
            string LongString;
            string[] DatArray = new string[2001];
            StringBuilder strbuild = new StringBuilder();
            LongString = "";

            string[] fileData = File.ReadAllLines(infilename);

            i = 0;
            ktr = 0;
            while (ktr < fileData.Count() - 1)  // not eof
            {
                LongString = fileData[ktr];

                if ((LongString.Contains("xx")) || LongString.Length < 10 || (i == 2000))
                {
                    this.opsStatus.Text = "Break " + ktr.ToString();
                    Refresh();

                    // ----Put a mandatory break at 2000----
                    if ((i == 2000) && (!LongString.Contains("xx")))
                        strbuild.AppendLine((i + 1).ToString() + " xx broken");
                    else
                    {
                        if (i > 0)
                            strbuild.AppendLine(i.ToString() + " xx ");
                    }

                    // ----Write out all the data to this point----
                    for (j = 1; j <= i; j++)
                        strbuild.AppendLine(DatArray[j].ToString());

                    // ----Write out crossover point for the mandatory break ---
                    if ((i == 2000) && ((!LongString.Contains('x')) || LongString.Length > 10))
                    {
                        strbuild.AppendLine(LongString);
                        // write out last one in loop ! 
                        i = 1;
                        DatArray[i] = LongString;
                    }
                    else
                    {
                        i = 0;
                        ktr = ktr + 1;   // file counter
                    }
                }
                else
                // ----Add new data ----
                {
                    i = i + 1;
                    string line1 = Regex.Replace(LongString, @"\s+", " ");
                    string[] linesplt = line1.Split(' ');

                    DatArray[i] = linesplt[intindx1] + " " + linesplt[intindx2] + " " + linesplt[intindx3] + " " + linesplt[intindx4]
                        + " " + linesplt[intindx5] + " " + linesplt[intindx6];
                    ktr = ktr + 1;   // file counter
                }

            }  // while through file

            // ----Write out the last section----
            if (i > 1)
            {
                if ((i == 2000) && (!LongString.Contains('x')) && LongString.Length > 10)
                    strbuild.AppendLine((i + 1).ToString() + " xx broken");
                else
                {
                    if (i > 0)
                        strbuild.AppendLine(i.ToString() + " xx ");
                }
                for (j = 1; j <= i; j++)
                    strbuild.AppendLine(DatArray[j]);
            }

            string directory = @"d:\codes\library\matlab\";
            File.WriteAllText(directory + "t.out", strbuild.ToString());
        }  // fixdat


        public void testlambertuniv()
        {
            double[] v1t = new double[3];
            double[] v2t = new double[3];
            double[] v1t1 = new double[3];
            double[] v2t1 = new double[3];
            double[] v1t2 = new double[3];
            double[] v2t2 = new double[3];
            double[] v1t3 = new double[3];
            double[] v2t3 = new double[3];
            double[] v1t4 = new double[3];
            double[] v2t4 = new double[3];
            double[,] tbiSu = new double[10, 3];
            double[,] tbiLu = new double[10, 3];
            double[] r1 = new double[3];
            double[] r2 = new double[3];
            double[] r3 = new double[3];
            double[] r4 = new double[3];
            double[] v1 = new double[3];
            double[] dv1 = new double[3];
            double[] dv11 = new double[3];
            double[] dv12 = new double[3];
            double[] dv13 = new double[3];
            double[] dv14 = new double[3];
            double[] v2 = new double[3];
            double[] v3 = new double[3];
            double[] v4 = new double[3];
            double[] dv2 = new double[3];
            double[] dv21 = new double[3];
            double[] dv22 = new double[3];
            double[] dv23 = new double[3];
            double[] dv24 = new double[3];
            double kbi, tof, dtwait, altpad, ang, dtsec;
            Int32 nrev, i, j;
            string errorsum = "";
            string errorout = "";
            char show = 'n';     // for test180, show = n, show180 = y
            //char show180 = 'n';  // for testlamb known show = y, show180 = n, n/n for envelope
            char hitearth, dm, de;
            nrev = 0;
            r1 = new double[] { 2.500000 * AstroLibr.gravConst.re, 0.000000, 0.000000 };
            r2 = new double[] { 1.9151111 * AstroLibr.gravConst.re, 1.6069690 * AstroLibr.gravConst.re, 0.000000 };
            // assume circular initial orbit for vel calcs
            v1 = new double[] { 0.0, Math.Sqrt(AstroLibr.gravConst.mu / r1[0]), 0.0 };
            ang = Math.Atan(r2[1] / r2[0]);
            v2 = new double[] { -Math.Sqrt(AstroLibr.gravConst.mu / r2[1]) * Math.Cos(ang), Math.Sqrt(AstroLibr.gravConst.mu / r2[0]) * Math.Sin(ang), 0.0 };
            dtsec = 76.0 * 60.0;
            altpad = 100.0;  // 100 km


            AstroLibr.lambertumins(r1, r2, 1, 'S', out kbi, out tof);
            tbiSu[1, 1] = kbi;
            tbiSu[1, 2] = tof;
            AstroLibr.lambertumins(r1, r2, 2, 'S', out kbi, out tof);
            tbiSu[2, 1] = kbi;
            tbiSu[2, 2] = tof;
            AstroLibr.lambertumins(r1, r2, 3, 'S', out kbi, out tof);
            tbiSu[3, 1] = kbi;
            tbiSu[3, 2] = tof;
            AstroLibr.lambertumins(r1, r2, 4, 'S', out kbi, out tof);
            tbiSu[4, 1] = kbi;
            tbiSu[4, 2] = tof;
            AstroLibr.lambertumins(r1, r2, 5, 'S', out kbi, out tof);
            tbiSu[5, 1] = kbi;
            tbiSu[5, 2] = tof;

            AstroLibr.lambertumins(r1, r2, 1, 'L', out kbi, out tof);
            tbiLu[1, 1] = kbi;
            tbiLu[1, 2] = tof;
            AstroLibr.lambertumins(r1, r2, 2, 'L', out kbi, out tof);
            tbiLu[2, 1] = kbi;
            tbiLu[2, 2] = tof;
            AstroLibr.lambertumins(r1, r2, 3, 'L', out kbi, out tof);
            tbiLu[3, 1] = kbi;
            tbiLu[3, 2] = tof;
            AstroLibr.lambertumins(r1, r2, 4, 'L', out kbi, out tof);
            tbiLu[4, 1] = kbi;
            tbiLu[4, 2] = tof;
            AstroLibr.lambertumins(r1, r2, 5, 'L', out kbi, out tof);
            tbiLu[5, 1] = kbi;
            tbiLu[5, 2] = tof;


            if (show == 'y')
                strbuild.AppendLine(" TEST ------------------ s/l  d  0 rev ------------------");
            hitearth = ' ';
            dm = 'S';
            de = 'L';
            nrev = 0;
            dtwait = 0.0;


            AstroLibr.lambertuniv(r1, r2, v1, dm, de, nrev, dtsec, kbi,
                          altpad, 'y', out v1t, out v2t, out hitearth, out errorsum, out errorout);

            AstroLibr.lambertuniv(r1, r2, v1, dm, 'H', nrev, dtsec, kbi,
              altpad, 'y', out v1t, out v2t, out hitearth, out errorsum, out errorout);

            dtsec = 21000.0;
            AstroLibr.lambertuniv(r1, r2, v1, 'S', 'H', 1, dtsec, kbi,
              altpad, 'y', out v1t, out v2t, out hitearth, out errorsum, out errorout);
            AstroLibr.lambertuniv(r1, r2, v1, 'S', 'L', 1, dtsec, kbi,
              altpad, 'y', out v1t, out v2t, out hitearth, out errorsum, out errorout);
            AstroLibr.lambertuniv(r1, r2, v1, 'L', 'H', 1, dtsec, kbi,
              altpad, 'y', out v1t, out v2t, out hitearth, out errorsum, out errorout);
            AstroLibr.lambertuniv(r1, r2, v1, 'L', 'L', 1, dtsec, kbi,
              altpad, 'y', out v1t, out v2t, out hitearth, out errorsum, out errorout);


            // test chap 7 fig 7-18 runs
            if (show == 'y')
                strbuild.AppendLine(" TEST ------------------ s/l  d  0 rev ------------------");
            hitearth = ' ';
            dm = 'S';
            de = 'L';
            nrev = 0;
            dtwait = 0.0;

            // fig 7-18 fixed tgt and int
            r1 = new double[] { -6518.1083, -2403.8479, -22.1722 };
            v1 = new double[] { 2.604057, -7.105717, -0.263218 };
            r2 = new double[] { 6697.4756, 1794.5832, 0.0 };
            v2 = new double[] { -1.962373, 7.323674, 0.0 };
            strbuild.AppendLine("dtwait  dtsec       dv1       dv2 ");
            this.opsStatus.Text = "Status:  on case 80a";
            Refresh();
            //           for (i = 0; i < 250; i++)
            {
                i = 0;
                dtsec = i * 60.0;
                AstroLibr.lambertuniv(r1, r2, v1, 'S', 'L', nrev, dtsec, kbi,
                              altpad, 'y', out v1t1, out v2t1, out hitearth, out errorsum, out errorout);
                AstroLibr.lambertuniv(r1, r2, v1, 'S', 'H', nrev, dtsec, kbi,
                              altpad, 'y', out v1t2, out v2t2, out hitearth, out errorsum, out errorout);
                AstroLibr.lambertuniv(r1, r2, v1, 'L', 'L', nrev, dtsec, kbi,
                              altpad, 'y', out v1t3, out v2t3, out hitearth, out errorsum, out errorout);
                AstroLibr.lambertuniv(r1, r2, v1, 'L', 'H', nrev, dtsec, kbi,
                              altpad, 'y', out v1t4, out v2t4, out hitearth, out errorsum, out errorout);

                if (errorout.Contains("ok"))
                {
                    MathTimeLibr.addvec(1.0, v1t1, -1.0, v1, out dv11);
                    MathTimeLibr.addvec(1.0, v2t1, -1.0, v2, out dv21);
                    MathTimeLibr.addvec(1.0, v1t2, -1.0, v1, out dv12);
                    MathTimeLibr.addvec(1.0, v2t2, -1.0, v2, out dv22);
                    MathTimeLibr.addvec(1.0, v1t3, -1.0, v1, out dv13);
                    MathTimeLibr.addvec(1.0, v2t3, -1.0, v2, out dv23);
                    MathTimeLibr.addvec(1.0, v1t4, -1.0, v1, out dv14);
                    MathTimeLibr.addvec(1.0, v2t4, -1.0, v2, out dv24);
                    strbuild.AppendLine(dtwait.ToString() + " " + dtsec.ToString() +
                        "  " + MathTimeLibr.mag(dv11).ToString() + "  " + MathTimeLibr.mag(dv21).ToString() +
                        "  " + MathTimeLibr.mag(dv12).ToString() + "  " + MathTimeLibr.mag(dv22).ToString() +
                        "  " + MathTimeLibr.mag(dv13).ToString() + "  " + MathTimeLibr.mag(dv23).ToString() +
                        "  " + MathTimeLibr.mag(dv14).ToString() + "  " + MathTimeLibr.mag(dv24).ToString());
                }
                else
                    strbuild.AppendLine(errorsum + " " + errorout);
            }


            // fig 7-19 moving tgt
            r1 = new double[] { 5328.7862, 4436.1273, 101.4720 };
            v1 = new double[] { -4.864779, 5.816486, 0.240163 };
            r2 = new double[] { 6697.4756, 1794.5831, 0.0 };
            v2 = new double[] { -1.962372, 7.323674, 0.0 };
            strbuild.AppendLine("dtwait  dtsec       dv1       dv2 ");
            // diff vectors, needs new umins

            this.opsStatus.Text = "Status:  on case 80b";
            Refresh();
            //           for (i = 0; i < 250; i++)
            {
                i = 0;
                dtsec = i * 60.0;
                AstroLibr.kepler(r2, v2, dtsec, out r3, out v3);
                //                for (j = 0; j < 25; j++)
                {
                    j = 0;
                    dtwait = j * 10.0;
                    dtwait = 0.0;  // set to 0 for now

                    AstroLibr.lambertuniv(r1, r3, v1, 'S', 'L', nrev, dtsec, kbi,
                                  altpad, 'y', out v1t1, out v2t1, out hitearth, out errorsum, out errorout);
                    AstroLibr.lambertuniv(r1, r3, v1, 'S', 'H', nrev, dtsec, kbi,
                                  altpad, 'y', out v1t2, out v2t2, out hitearth, out errorsum, out errorout);
                    AstroLibr.lambertuniv(r1, r3, v1, 'L', 'L', nrev, dtsec, kbi,
                                  altpad, 'y', out v1t3, out v2t3, out hitearth, out errorsum, out errorout);
                    AstroLibr.lambertuniv(r1, r3, v1, 'L', 'H', nrev, dtsec, kbi,
                                  altpad, 'y', out v1t4, out v2t4, out hitearth, out errorsum, out errorout);
                    if (errorout.Contains("ok"))
                    {
                        MathTimeLibr.addvec(1.0, v1t1, -1.0, v1, out dv11);
                        MathTimeLibr.addvec(1.0, v2t1, -1.0, v3, out dv21);
                        MathTimeLibr.addvec(1.0, v1t2, -1.0, v1, out dv12);
                        MathTimeLibr.addvec(1.0, v2t2, -1.0, v3, out dv22);
                        MathTimeLibr.addvec(1.0, v1t3, -1.0, v1, out dv13);
                        MathTimeLibr.addvec(1.0, v2t3, -1.0, v3, out dv23);
                        MathTimeLibr.addvec(1.0, v1t4, -1.0, v1, out dv14);
                        MathTimeLibr.addvec(1.0, v2t4, -1.0, v3, out dv24);
                        strbuild.AppendLine(dtwait.ToString() + " " + dtsec.ToString() +
                            "  " + MathTimeLibr.mag(dv11).ToString() + "  " + MathTimeLibr.mag(dv21).ToString() +
                            "  " + MathTimeLibr.mag(dv12).ToString() + "  " + MathTimeLibr.mag(dv22).ToString() +
                            "  " + MathTimeLibr.mag(dv13).ToString() + "  " + MathTimeLibr.mag(dv23).ToString() +
                            "  " + MathTimeLibr.mag(dv14).ToString() + "  " + MathTimeLibr.mag(dv24).ToString());
                    }
                    else
                        strbuild.AppendLine(errorsum + " " + errorout);
                }
            }


            // fig 7-21
            StringBuilder strbuildFig = new StringBuilder();
            r1 = new double[] { -6175.1034, 2757.0706, 1626.6556 };
            v1 = new double[] { 2.376641, 1.139677, 7.078097 };
            r2 = new double[] { -6078.007289, 2796.641859, 1890.7135 };
            v2 = new double[] { 2.654700, 1.018600, 7.015400 };

            strbuildFig.AppendLine("dtwait  dtsec       dv1       dv2 ");
            this.opsStatus.Text = "Status:  on case 80c";
            Refresh();
            int totaldts = 15000;
            int totaldtw = 30000;
            int step1 = 60;   // 60 orig
            int step2 = 600;  // 600 orig
            int stop1 = (int)(totaldts / step1);
            int stop2 = (int)(totaldtw / step2);
            for (i = 0; i < stop1; i++)  // orig 250, 15000 s total 
            {
                dtsec = i * step1;    // orig 60
                AstroLibr.kepler(r1, v1, dtsec, out r4, out v4);
                for (j = 0; j < stop2; j++)  // orig 25 600*25 = 15000 s total
                {
                    dtwait = j * step2;   // orig 600
                    AstroLibr.kepler(r2, v2, dtsec + dtwait, out r3, out v3);

                    AstroLibr.lambertuniv(r4, r3, v4, 's', 'd', nrev, dtsec, kbi,
                                  altpad, 'y', out v1t1, out v2t1, out hitearth, out errorsum, out errorout);
                    AstroLibr.lambertuniv(r4, r3, v4, 's', 'r', nrev, dtsec, kbi,
                                  altpad, 'y', out v1t2, out v2t2, out hitearth, out errorsum, out errorout);
                    AstroLibr.lambertuniv(r4, r3, v4, 'l', 'd', nrev, dtsec, kbi,
                                  altpad, 'y', out v1t3, out v2t3, out hitearth, out errorsum, out errorout);
                    AstroLibr.lambertuniv(r4, r3, v4, 'l', 'r', nrev, dtsec, kbi,
                                  altpad, 'y', out v1t4, out v2t4, out hitearth, out errorsum, out errorout);
                    if (errorout.Contains("ok"))
                    {
                        MathTimeLibr.addvec(1.0, v1t1, -1.0, v4, out dv11);
                        MathTimeLibr.addvec(1.0, v2t1, -1.0, v3, out dv21);
                        MathTimeLibr.addvec(1.0, v1t2, -1.0, v4, out dv12);
                        MathTimeLibr.addvec(1.0, v2t2, -1.0, v3, out dv22);
                        MathTimeLibr.addvec(1.0, v1t3, -1.0, v4, out dv13);
                        MathTimeLibr.addvec(1.0, v2t3, -1.0, v3, out dv23);
                        MathTimeLibr.addvec(1.0, v1t4, -1.0, v4, out dv14);
                        MathTimeLibr.addvec(1.0, v2t4, -1.0, v3, out dv24);
                        strbuildFig.AppendLine(dtwait.ToString() + " " + dtsec.ToString() +
                            "  " + MathTimeLibr.mag(dv11).ToString() + "  " + MathTimeLibr.mag(dv21).ToString() +
                            "  " + MathTimeLibr.mag(dv12).ToString() + "  " + MathTimeLibr.mag(dv22).ToString() +
                            "  " + MathTimeLibr.mag(dv13).ToString() + "  " + MathTimeLibr.mag(dv23).ToString() +
                            "  " + MathTimeLibr.mag(dv14).ToString() + "  " + MathTimeLibr.mag(dv24).ToString());
                    }
                    else
                        strbuildFig.AppendLine(errorsum + " " + errorout);
                }
            }

            // write data out
            string directory = @"D:\Codes\LIBRARY\Matlab\";
            File.WriteAllText(directory + "surfMovingSalltest.out", strbuildFig.ToString());
        }

        // test all the known problem cases for lambert
        // output these results separately to the testall directory
        private void testknowncases()
        {
            //  this file runs all the known problem cases. 
            double tusec = 806.8111238242922;
            Int16 numiter = 16;
            Int32 caseopt, nrev;
            double dtwait, dtsec;
            double[] r1 = new double[3];
            double[] r2 = new double[3];
            double[] v1 = new double[3];
            double[] v2 = new double[3];
            double[] v1tk = new double[3];
            double[] v2tk = new double[3];
            double[] v1tu = new double[3];
            double[] v2tu = new double[3];
            double[] v1tb = new double[3];
            double[] v2tb = new double[3];
            double[] v1tt = new double[3];
            double[] v2tt = new double[3];
            string detailSum, detailAll, errorout;
            double[] dv1 = new double[3];
            double[] dv2 = new double[3];
            double[] dv1t = new double[3];
            double[] dv2t = new double[3];
            double[] r3h = new double[3];
            double[] v3h = new double[3];
            double[] dr = new double[3];
            double ang, f, g, gdot, s, tau;
            double tmin, tminp, tminenergy;
            StringBuilder strbuildAll = new StringBuilder();
            detailSum = "";
            detailAll = "";
            //int i;
            char dm, de, hitearth;
            // for test180, show = n, show180 = y
            // for testlamb, show = y, show180 = n known cases
            // for envelope, show = n, show180 = n 

            double altpadc = 200.0 / AstroLibr.gravConst.re;  // set 200 km for altitude you set as the over limit. 

            dtsec = 0.0;
            nrev = 0;
            for (caseopt = 0; caseopt <= 85; caseopt++) // 74
            {
                this.opsStatus.Text = "working on lambert case " + caseopt.ToString();
                Refresh();

                dtwait = 0.0;
                // Problem cases during evaluation
                switch (caseopt)
                {
                    case 0:
                        strbuildAll.AppendLine("\n-------- lambert test book pg 497 short way \n");
                        nrev = 0;
                        r1 = new double[] { 2.500000 * AstroLibr.gravConst.re, 0.000000, 0.000000 };
                        r2 = new double[] { 1.9151111 * AstroLibr.gravConst.re, 1.6069690 * AstroLibr.gravConst.re, 0.000000 };
                        // assume circular initial orbit for vel calcs
                        v1 = new double[] { 0.0, Math.Sqrt(AstroLibr.gravConst.mu / r1[0]), 0.0 };
                        ang = Math.Atan(r2[1] / r2[0]);
                        v2 = new double[] { -Math.Sqrt(AstroLibr.gravConst.mu / r2[1]) * Math.Cos(ang), Math.Sqrt(AstroLibr.gravConst.mu / r2[0]) * Math.Sin(ang), 0.0 };
                        dtsec = 99900.3;
                        dtsec = 76.0 * 60.0;
                        dtsec = 21000.0;

                        strbuildAll.AppendLine("r1 " + " " + r1[0].ToString("0.00000000000") + " " + r1[1].ToString("0.00000000000") + " " + r1[2].ToString("0.00000000000"));
                        strbuildAll.AppendLine("r2 " + " " + r2[0].ToString("0.00000000000") + " " + r2[1].ToString("0.00000000000") + " " + r2[2].ToString("0.00000000000"));
                        strbuildAll.AppendLine("v1 " + " " + v1[0].ToString("0.00000000000") + " " + v1[1].ToString("0.00000000000") + " " + v1[2].ToString("0.00000000000"));
                        strbuildAll.AppendLine("v2 " + " " + v2[0].ToString("0.00000000000") + " " + v2[1].ToString("0.00000000000") + " " + v2[2].ToString("0.00000000000"));
                        break;
                    case 1:
                        nrev = 1;
                        r2 = new double[] { -1105.78023519582, 2373.16130661458, 6713.89444816503 };
                        v2 = new double[] { 5.4720951867079, -4.39299050886976, 2.45681739563752 };
                        r1 = new double[] { 4938.49830042171, -1922.24810472241, 4384.68293292613 };
                        v1 = new double[] { 0.738204644165659, 7.20989453238397, 2.32877392066299 };
                        dtsec = 6372.69272563561; // 1ld
                        break;
                }  // switch

                strbuildAll.AppendLine("===== Lambert Case " + caseopt.ToString() + " === ");

                ang = Math.Atan(r2[1] / r2[0]);
                double magr1 = MathTimeLibr.mag(r1);
                double magr2 = MathTimeLibr.mag(r2);

                // this value stays constant in all calcs, vara changes with df
                double cosdeltanu = MathTimeLibr.dot(r1, r2) / (magr1 * magr2);

                //strbuild.AppendLine("now do findtbi calcs");
                //strbuild.AppendLine("iter       y         dtnew          psiold      psinew   psinew-psiold   dtdpsi      dtdpsi2    lower    upper     ");

                AstroLambertkLibr.lambertkmins1st(r1, r2, out s, out tau);
                strbuildAll.AppendLine(" s " + s.ToString(fmt) + " tau " + tau.ToString(fmt));

                double kbi, tof;
                double[,] tbidk = new double[10, 3];
                AstroLambertkLibr.lambertkmins(s, tau, 1, 'x', 'L', out kbi, out tof);
                tbidk[1, 1] = kbi;
                tbidk[1, 2] = tof / tusec;
                AstroLambertkLibr.lambertkmins(s, tau, 2, 'x', 'L', out kbi, out tof);
                tbidk[2, 1] = kbi;
                tbidk[2, 2] = tof / tusec;
                AstroLambertkLibr.lambertkmins(s, tau, 3, 'x', 'L', out kbi, out tof);
                tbidk[3, 1] = kbi;
                tbidk[3, 2] = tof / tusec;
                AstroLambertkLibr.lambertkmins(s, tau, 4, 'x', 'L', out kbi, out tof);
                tbidk[4, 1] = kbi;
                tbidk[4, 2] = tof / tusec;
                AstroLambertkLibr.lambertkmins(s, tau, 5, 'x', 'L', out kbi, out tof);
                tbidk[5, 1] = kbi;
                tbidk[5, 2] = tof / tusec;
                AstroLambertkLibr.lambertkmins(s, tau, 6, 'x', 'L', out kbi, out tof);
                tbidk[6, 1] = kbi;
                tbidk[6, 2] = tof / tusec;

                double[,] tbirk = new double[10, 3];
                AstroLambertkLibr.lambertkmins(s, tau, 1, 'x', 'H', out kbi, out tof);
                tbirk[1, 1] = kbi;
                tbirk[1, 2] = tof / tusec;
                AstroLambertkLibr.lambertkmins(s, tau, 2, 'x', 'H', out kbi, out tof);
                tbirk[2, 1] = kbi;
                tbirk[2, 2] = tof / tusec;
                AstroLambertkLibr.lambertkmins(s, tau, 3, 'x', 'H', out kbi, out tof);
                tbirk[3, 1] = kbi;
                tbirk[3, 2] = tof / tusec;
                AstroLambertkLibr.lambertkmins(s, tau, 4, 'x', 'H', out kbi, out tof);
                tbirk[4, 1] = kbi;
                tbirk[4, 2] = tof / tusec;
                AstroLambertkLibr.lambertkmins(s, tau, 5, 'x', 'H', out kbi, out tof);
                tbirk[5, 1] = kbi;
                tbirk[5, 2] = tof / tusec;
                AstroLambertkLibr.lambertkmins(s, tau, 6, 'x', 'H', out kbi, out tof);
                tbirk[6, 1] = kbi;
                tbirk[6, 2] = tof / tusec;

                strbuildAll.AppendLine("From k variables ");
                strbuildAll.AppendLine(" " + tbidk[1, 1].ToString("#0.00000000000") + "  " + (tbidk[1, 2] * tusec).ToString("0.00000000000") + " s " + (tbidk[1, 2]).ToString("0.00000000000") + " tu ");
                strbuildAll.AppendLine(" " + tbidk[2, 1].ToString("#0.00000000000") + "  " + (tbidk[2, 2] * tusec).ToString("0.00000000000") + " s " + (tbidk[2, 2]).ToString("0.00000000000") + " tu ");
                strbuildAll.AppendLine(" " + tbidk[3, 1].ToString("#0.00000000000") + "  " + (tbidk[3, 2] * tusec).ToString("0.00000000000") + " s " + (tbidk[3, 2]).ToString("0.00000000000") + " tu ");
                strbuildAll.AppendLine(" " + tbidk[4, 1].ToString("#0.00000000000") + "  " + (tbidk[4, 2] * tusec).ToString("0.00000000000") + " s " + (tbidk[4, 2]).ToString("0.00000000000") + " tu ");
                strbuildAll.AppendLine(" " + tbidk[5, 1].ToString("#0.00000000000") + "  " + (tbidk[5, 2] * tusec).ToString("0.00000000000") + " s " + (tbidk[5, 2]).ToString("0.00000000000") + " tu ");
                strbuildAll.AppendLine("");
                strbuildAll.AppendLine(" " + tbirk[1, 1].ToString("#0.00000000000") + "  " + (tbirk[1, 2] * tusec).ToString("0.00000000000") + " s " + (tbirk[1, 2]).ToString("0.00000000000") + " tu ");
                strbuildAll.AppendLine(" " + tbirk[2, 1].ToString("#0.00000000000") + "  " + (tbirk[2, 2] * tusec).ToString("0.00000000000") + " s " + (tbirk[2, 2]).ToString("0.00000000000") + " tu ");
                strbuildAll.AppendLine(" " + tbirk[3, 1].ToString("#0.00000000000") + "  " + (tbirk[3, 2] * tusec).ToString("0.00000000000") + " s " + (tbirk[3, 2]).ToString("0.00000000000") + " tu ");
                strbuildAll.AppendLine(" " + tbirk[4, 1].ToString("#0.00000000000") + "  " + (tbirk[4, 2] * tusec).ToString("0.00000000000") + " s " + (tbirk[4, 2]).ToString("0.00000000000") + " tu ");
                strbuildAll.AppendLine(" " + tbirk[5, 1].ToString("#0.00000000000") + "  " + (tbirk[5, 2] * tusec).ToString("0.00000000000") + " s " + (tbirk[5, 2]).ToString("0.00000000000") + " tu ");


                //strbuild.AppendLine("lambertTest" + caseopt.ToString() + " " + r1[0].ToString("0.00000000000") + " " + r1[1].ToString("0.00000000000") + " " + r1[2].ToString("0.00000000000") +
                //    " " + v1[0].ToString("0.00000000000") + " " + v1[1].ToString("0.00000000000") + " " + v1[2].ToString("0.00000000000") +
                //    " " + r2[0].ToString("0.00000000000") + " " + r2[1].ToString("0.00000000000") + " " + r2[2].ToString("0.00000000000") +
                //    " " + v2[0].ToString("0.00000000000") + " " + v2[1].ToString("0.00000000000") + " " + v2[2].ToString("0.00000000000") + " " + dtsec.ToString());

                AstroLibr.lambertminT(r1, r2, 'S', 'L', 1, out tmin, out tminp, out tminenergy);
                strbuildAll.AppendLine("mint S " + tmin.ToString("0.0000") + " minp " + tminp.ToString("0.0000") + " minener " + tminenergy.ToString("0.0000"));
                AstroLibr.lambertminT(r1, r2, 'S', 'L', 2, out tmin, out tminp, out tminenergy);
                strbuildAll.AppendLine("mint S " + tmin.ToString("0.0000") + " minp " + tminp.ToString("0.0000") + " minener " + tminenergy.ToString("0.0000"));
                AstroLibr.lambertminT(r1, r2, 'S', 'L', 3, out tmin, out tminp, out tminenergy);
                strbuildAll.AppendLine("mint S " + tmin.ToString("0.0000") + " minp " + tminp.ToString("0.0000") + " minener " + tminenergy.ToString("0.0000"));

                AstroLibr.lambertminT(r1, r2, 'L', 'H', 1, out tmin, out tminp, out tminenergy);
                strbuildAll.AppendLine("mint L " + tmin.ToString("0.0000") + " minp " + tminp.ToString("0.0000") + " minener " + tminenergy.ToString("0.0000"));
                AstroLibr.lambertminT(r1, r2, 'L', 'H', 2, out tmin, out tminp, out tminenergy);
                strbuildAll.AppendLine("mint L " + tmin.ToString("0.0000") + " minp " + tminp.ToString("0.0000") + " minener " + tminenergy.ToString("0.0000"));
                AstroLibr.lambertminT(r1, r2, 'L', 'H', 3, out tmin, out tminp, out tminenergy);
                strbuildAll.AppendLine("mint L " + tmin.ToString("0.0000") + " minp " + tminp.ToString("0.0000") + " minener " + tminenergy.ToString("0.0000"));

                char modecon = 'n';  // 'c' to shortcut bad cases (hitearth) at iter 3 or 'n'

                strbuild.AppendLine(" TEST ------------------ s/l  H  0 rev ------------------");
                hitearth = ' ';
                dm = 'S';
                de = 'H';
                AstroLambertkLibr.lambertK(r1, v1, r2, dm, de, nrev, dtwait, dtsec, 0.0, 0.0, numiter, altpadc, modecon, 'n',
                    out v1tk, out v2tk, out f, out g, out gdot, out hitearth, out errorout, out detailSum, out detailAll);
                strbuildAll.AppendLine(detailAll);
                //strbuild.AppendLine("k#" + caseopt + " " + detailSum + " diffs " + MathTimeLib::mag(dr).ToString("0.00000000000"));
                strbuild.AppendLine("lamk v1t " + " " + v1tk[0].ToString("0.00000000000") + " " + v1tk[1].ToString("0.00000000000") + " " + v1tk[2].ToString("0.00000000000"));
                strbuild.AppendLine("lamk v2t " + " " + v2tk[0].ToString("0.00000000000") + " " + v2tk[1].ToString("0.00000000000") + " " + v2tk[2].ToString("0.00000000000"));
                //strbuild.AppendLine(magv1t.ToString("0.0000000").PadLeft(12) + " " + magv2t.ToString("0.0000000").PadLeft(12));

                AstroLibr.kepler(r1, v1tk, dtsec, out r3h, out v3h);
                strbuild.AppendLine("r3h " + " " + r3h[0].ToString("0.00000000000") + " " + r3h[1].ToString("0.00000000000") + " " + r3h[2].ToString("0.00000000000"));
                for (int j = 0; j < 3; j++)
                    dr[j] = r2[j] - r3h[j];
                if (MathTimeLibr.mag(dr) > 0.01)
                    strbuild.AppendLine("velk does not get to r2 position (km) " + MathTimeLibr.mag(dr).ToString() + "\n");

                AstroLibr.lambertuniv(r1, r2, v1, dm, de, nrev, dtsec, 0.0, altpadc * AstroLibr.gravConst.re, 'n', out v1tu, out v2tu, out hitearth, out detailSum, out detailAll);
                //strbuild.AppendLine(detailSum);
                strbuild.AppendLine("univ v1t " + " " + v1tu[0].ToString("0.00000000000") + " " + v1tu[1].ToString("0.00000000000") + " " + v1tu[2].ToString("0.00000000000"));
                strbuild.AppendLine("univ v2t " + " " + v2tu[0].ToString("0.00000000000") + " " + v2tu[1].ToString("0.00000000000") + " " + v2tu[2].ToString("0.00000000000"));
                AstroLibr.kepler(r1, v1tu, dtsec, out r3h, out v3h);
                for (int j = 0; j < 3; j++)
                    dr[j] = r2[j] - r3h[j];
                if (MathTimeLibr.mag(dr) > 0.01)
                    strbuild.AppendLine("velu does not get to r2 position (km) " + MathTimeLibr.mag(dr).ToString() + "\n");

                for (int j = 0; j < 3; j++)
                {
                    dv1[j] = v1tk[j] - v1tu[j];
                    dv2[j] = v2tk[j] - v2tu[j];
                }
                if (MathTimeLibr.mag(dv1) > 0.01 || MathTimeLibr.mag(dv2) > 0.01)
                    strbuild.AppendLine("velk does not match velu \n");

                AstroLibr.lambertbattin(r1, r2, v1, dm, de, nrev, dtsec, altpadc * AstroLibr.gravConst.re, 'n', out v1tb, out v2tb, out hitearth, out detailSum, out detailAll);
                //strbuild.AppendLine(detailSum);
                strbuild.AppendLine("batt v1t " + " " + v1tb[0].ToString("0.00000000000") + " " + v1tb[1].ToString("0.00000000000") + " " + v1tb[2].ToString("0.00000000000"));
                strbuild.AppendLine("batt v2t " + " " + v2tb[0].ToString("0.00000000000") + " " + v2tb[1].ToString("0.00000000000") + " " + v2tb[2].ToString("0.00000000000"));
                AstroLibr.kepler(r1, v1tb, dtsec, out r3h, out v3h);
                for (int j = 0; j < 3; j++)
                    dr[j] = r2[j] - r3h[j];
                if (MathTimeLibr.mag(dr) > 0.01)
                    strbuild.AppendLine("velb does not get to r2 position (km) " + MathTimeLibr.mag(dr).ToString() + "\n");
                //strbuild.AppendLine("diffs " + MathTimeLibr.mag(dr).ToString("0.00000000000"));

                for (int j = 0; j < 3; j++)
                {
                    dv1[j] = v1tk[j] - v1tb[j];
                    dv2[j] = v2tk[j] - v2tb[j];
                }
                if (MathTimeLibr.mag(dv1) > 0.01 || MathTimeLibr.mag(dv2) > 0.01)
                    strbuild.AppendLine("velk does not match velb \n");
                //strbuild.AppendLine("diffs " + MathTimeLib::mag(dr).ToString("0.00000000000"));

                //// teds approach
                //double[] r2ted = new double[] { r2[0], r2[1], r2[2] };
                //double[] v2ted = new double[] { v2[0], v2[1], v2[2] };
                //double[] r1ted = new double[] { r1[0], r1[1], r1[2] };
                //double[] v1ted = new double[] { v1[0], v1[1], v1[2] };
                //Cartesian r1com = new Cartesian(r1ted[0], r1ted[1], r1ted[2]);
                //Cartesian v1com = new Cartesian(v1ted[0], v1ted[1], v1ted[2]);
                //Cartesian r2com = new Cartesian(r2ted[0], r2ted[1], r2ted[2]);
                //Cartesian v2com = new Cartesian(v2ted[0], v2ted[1], v2ted[2]);
                //var result = LambertDeltaV.FindMinimumDeltaV(r2com, v2com, r1com, v1com, dtsec, Lambert.EngagementType.Prox, 0);  // .Intercept
                //double[] v1tr = { result.Velocities.Item1.X, result.Velocities.Item1.Y, result.Velocities.Item1.Z };  // LambertKMin/s
                //double[] v2tr = { result.Velocities.Item2.X, result.Velocities.Item2.Y, result.Velocities.Item2.Z };  // LambertKMin/s
                //for (int i = 0; i < 3; i++)
                //{
                //    dv1t[i] = v1[i] - v1tr[i];
                //    dv2t[i] = v2[i] - v2tr[i];
                //}
                //magv1tt = MathTimeLibr.mag(dv1);
                //magv2tt = MathTimeLibr.mag(dv2);
                ////strbuild.AppendLine(detailAll);  // dont do again
                //double knew = 1.1;
                //detailAll = ("T" + detailSum.PadLeft(2) + result.LambertCalculations.ToString().PadLeft(4) + 0.ToString().PadLeft(3) + "   " + dm + "  " + df + dtwait.ToString("0.00000000000").PadLeft(15) +
                //           dtsec.ToString("0.00000000000").PadLeft(15) + knew.ToString("0.00000000000").PadLeft(15) +
                //           v1tr[0].ToString("0.00000000000").PadLeft(15) + v1tr[1].ToString("0.00000000000").PadLeft(15) + v1tr[2].ToString("0.00000000000").PadLeft(15) +
                //           v2tr[0].ToString("0.00000000000").PadLeft(15) + v2tr[1].ToString("0.00000000000").PadLeft(15) + v2tr[2].ToString("0.00000000000").PadLeft(15) +
                //           (Math.Acos(cosdeltanu) * 180 / Math.PI).ToString("0.00000000000").PadLeft(15) + caseopt + hitearth);
                ////                strbuild.AppendLine(detailAll);

                ////                strbuild.AppendLine(magv1tt.ToString("0.00000000000") + "  " + magv2tt.ToString("0.00000000000") + " \n");
                //if ((Math.Abs(magv1t - magv1tt) > 0.01) || (Math.Abs(magv2t - magv2tt) > 0.01))
                //    strbuild.AppendLine("Error between the approaches");

                // timing of routines
                //var watch = System.Diagnostics.Stopwatch.StartNew();
                //int l = 0;
                //for (l = 1; l < 500; l++)
                {
                    strbuild.AppendLine(" TEST ------------------ s/l L 0 rev ------------------");
                    dm = 'L';
                    de = 'L';
                    // k near 180 is about 53017 while battin is 30324!
                    AstroLambertkLibr.lambertK(r1, v1, r2, dm, de, nrev, dtwait, dtsec, 0.0, 0.0, numiter, altpadc, modecon, 'n',
                        out v1tk, out v2tk, out f, out g, out gdot, out hitearth, out errorout, out detailSum, out detailAll);
                    strbuildAll.AppendLine(detailAll);
                    //strbuild.AppendLine("k#" + caseopt + " " + detailSum + " diffs " + MathTimeLibr.mag(dr).ToString("0.00000000000"));
                    strbuild.AppendLine("lamk v1t " + " " + v1tk[0].ToString("0.00000000000") + " " + v1tk[1].ToString("0.00000000000") + " " + v1tk[2].ToString("0.00000000000"));
                    strbuild.AppendLine("lamk v2t " + " " + v2tk[0].ToString("0.00000000000") + " " + v2tk[1].ToString("0.00000000000") + " " + v2tk[2].ToString("0.00000000000"));
                    //strbuild.AppendLine(magv1t.ToString("0.0000000").PadLeft(12) + " " + magv2t.ToString("0.0000000").PadLeft(12));

                    AstroLibr.kepler(r1, v1tk, dtsec, out r3h, out v3h);
                    for (int j = 0; j < 3; j++)
                        dr[j] = r2[j] - r3h[j];
                    if (MathTimeLibr.mag(dr) > 0.01)
                        strbuild.AppendLine("velk does not get to r2 (km) position " + MathTimeLibr.mag(dr).ToString() + "\n");

                    AstroLibr.lambertuniv(r1, r2, v1, dm, de, nrev, dtsec, 0.0, altpadc * AstroLibr.gravConst.re, 'n', out v1tu, out v2tu, out hitearth, out detailSum, out detailAll);
                    //strbuild.AppendLine(detailSum);
                    strbuild.AppendLine("univ v1t " + " " + v1tu[0].ToString("0.00000000000") + " " + v1tu[1].ToString("0.00000000000") + " " + v1tu[2].ToString("0.00000000000"));
                    strbuild.AppendLine("univ v2t " + " " + v2tu[0].ToString("0.00000000000") + " " + v2tu[1].ToString("0.00000000000") + " " + v2tu[2].ToString("0.00000000000"));
                    AstroLibr.kepler(r1, v1tu, dtsec, out r3h, out v3h);
                    for (int j = 0; j < 3; j++)
                        dr[j] = r2[j] - r3h[j];
                    if (MathTimeLibr.mag(dr) > 0.01)
                        strbuild.AppendLine("velu does not get to r2 (km) position " + MathTimeLibr.mag(dr).ToString() + "\n");

                    for (int j = 0; j < 3; j++)
                    {
                        dv1[j] = v1tk[j] - v1tu[j];
                        dv2[j] = v2tk[j] - v2tu[j];
                    }
                    if (MathTimeLibr.mag(dv1) > 0.01 || MathTimeLibr.mag(dv2) > 0.01)
                        strbuild.AppendLine("velk does not match velu \n");

                    AstroLibr.lambertbattin(r1, r2, v1, dm, de, nrev, dtsec, altpadc * AstroLibr.gravConst.re, 'n', out v1tb, out v2tb, out hitearth, out detailSum, out detailAll);
                    //strbuild.AppendLine(detailSum);
                    strbuild.AppendLine("batt v1t " + " " + v1tb[0].ToString("0.00000000000") + " " + v1tb[1].ToString("0.00000000000") + " " + v1tb[2].ToString("0.00000000000"));
                    strbuild.AppendLine("batt v2t " + " " + v2tb[0].ToString("0.00000000000") + " " + v2tb[1].ToString("0.00000000000") + " " + v2tb[2].ToString("0.00000000000"));
                    AstroLibr.kepler(r1, v1tb, dtsec, out r3h, out v3h);
                    for (int j = 0; j < 3; j++)
                        dr[j] = r2[j] - r3h[j];
                    if (MathTimeLibr.mag(dr) > 0.01)
                        strbuild.AppendLine("velb does not get to r2 (km) position " + MathTimeLibr.mag(dr).ToString() + "\n");
                    //strbuild.AppendLine("diffs " + MathTimeLibr.mag(dr).ToString("0.00000000000"));

                    for (int j = 0; j < 3; j++)
                    {
                        dv1[j] = v1tk[j] - v1tb[j];
                        dv2[j] = v2tk[j] - v2tb[j];
                    }
                    if (MathTimeLibr.mag(dv1) > 0.01 || MathTimeLibr.mag(dv2) > 0.01)
                        strbuild.AppendLine("velk does not match velb \n");
                    //strbuild.AppendLine("diffs " + MathTimeLib::mag(dr).ToString("0.00000000000"));

                    for (int j = 0; j < 3; j++)
                    {
                        dv1[j] = v1tk[j] - v1tb[j];
                        dv2[j] = v2tk[j] - v2tb[j];
                    }
                    if (MathTimeLibr.mag(dv1) > 0.01 || MathTimeLibr.mag(dv2) > 0.01)
                        strbuild.AppendLine("velk does not match velb \n");
                    //strbuild.AppendLine("diffs " + MathTimeLib::mag(dr).ToString("0.00000000000"));
                }
                //watch.Stop();
                //var elapsedMs = watch.ElapsedMilliseconds;
                //Console.WriteLine(watch.ElapsedMilliseconds); 

                // use random nrevs, but check if nrev = 0 and set to 1
                // but then you have to check that there is enough time for 1 rev
                int nnrev = nrev;
                if (nnrev == 0)
                    nnrev = 1;

                AstroLibr.lambertminT(r1, r2, 'S', 'L', nnrev, out tmin, out tminp, out tminenergy);
                strbuildAll.AppendLine("mint S " + tmin.ToString("0.0000") + " minp " + tminp.ToString("0.0000") + " minener " + tminenergy.ToString("0.0000"));
                AstroLibr.lambertminT(r1, r2, 'L', 'L', nnrev, out tmin, out tminp, out tminenergy);
                strbuildAll.AppendLine("mint L " + tmin.ToString("0.0000") + " minp " + tminp.ToString("0.0000") + " minener " + tminenergy.ToString("0.0000"));

                strbuild.AppendLine(" TEST ------------------ S  L " + nnrev.ToString() + " rev ------------------");
                //if (dtsec / tusec >= tbidk[nnrev, 2])
                // do inside lambertk now
                {
                    dm = 'S';
                    de = 'L';
                    AstroLambertkLibr.lambertkmins(s, tau, nnrev, dm, de, out kbi, out tof);
                    AstroLambertkLibr.lambertK(r1, v1, r2, dm, de, nnrev, dtwait, dtsec, tof, kbi, numiter, altpadc, modecon, 'n',
                        out v1tk, out v2tk, out f, out g, out gdot, out hitearth, out errorout, out detailSum, out detailAll);
                    strbuildAll.AppendLine(detailAll);
                    //strbuild.AppendLine("k#" + caseopt + " " + detailSum + " diffs " + MathTimeLibr.mag(dr).ToString("0.00000000000"));
                    strbuild.AppendLine("lamk v1t " + " " + v1tk[0].ToString("0.00000000000") + " " + v1tk[1].ToString("0.00000000000") + " " + v1tk[2].ToString("0.00000000000"));
                    strbuild.AppendLine("lamk v2t " + " " + v2tk[0].ToString("0.00000000000") + " " + v2tk[1].ToString("0.00000000000") + " " + v2tk[2].ToString("0.00000000000"));
                    //strbuild.AppendLine(magv1t.ToString("0.0000000").PadLeft(12) + " " + magv2t.ToString("0.0000000").PadLeft(12));

                    AstroLibr.kepler(r1, v1tk, dtsec, out r3h, out v3h);
                    for (int j = 0; j < 3; j++)
                        dr[j] = r2[j] - r3h[j];
                    if (MathTimeLibr.mag(dr) > 0.01)
                        strbuild.AppendLine("velk does not get to r2 (km) position " + MathTimeLibr.mag(dr).ToString() + "\n");

                    AstroLibr.lambertumins(r1, r2, nnrev, dm, out kbi, out tof);
                    AstroLibr.lambertuniv(r1, r2, v1, dm, de, nnrev, dtsec, kbi, altpadc * AstroLibr.gravConst.re, 'n', out v1tu, out v2tu, out hitearth, out detailSum, out detailAll);
                    //strbuild.AppendLine(detailSum);
                    strbuild.AppendLine("univ v1t " + " " + v1tu[0].ToString("0.00000000000") + " " + v1tu[1].ToString("0.00000000000") + " " + v1tu[2].ToString("0.00000000000"));
                    strbuild.AppendLine("univ v2t " + " " + v2tu[0].ToString("0.00000000000") + " " + v2tu[1].ToString("0.00000000000") + " " + v2tu[2].ToString("0.00000000000"));
                    AstroLibr.kepler(r1, v1tu, dtsec, out r3h, out v3h);
                    for (int j = 0; j < 3; j++)
                        dr[j] = r2[j] - r3h[j];
                    if (MathTimeLibr.mag(dr) > 0.01)
                        strbuild.AppendLine("velu does not get to r2 (km) position " + MathTimeLibr.mag(dr).ToString() + "\n");

                    for (int j = 0; j < 3; j++)
                    {
                        dv1[j] = v1tk[j] - v1tu[j];
                        dv2[j] = v2tk[j] - v2tu[j];
                    }
                    if (MathTimeLibr.mag(dv1) > 0.01 || MathTimeLibr.mag(dv2) > 0.01)
                        strbuild.AppendLine("velk does not match velu \n");

                    AstroLibr.lambertbattin(r1, r2, v1, dm, de, nnrev, dtsec, altpadc * AstroLibr.gravConst.re, 'n', out v1tb, out v2tb, out hitearth, out detailSum, out detailAll);
                    //strbuild.AppendLine(detailSum);
                    strbuild.AppendLine("batt v1t " + " " + v1tb[0].ToString("0.00000000000") + " " + v1tb[1].ToString("0.00000000000") + " " + v1tb[2].ToString("0.00000000000"));
                    strbuild.AppendLine("batt v2t " + " " + v2tb[0].ToString("0.00000000000") + " " + v2tb[1].ToString("0.00000000000") + " " + v2tb[2].ToString("0.00000000000"));
                    AstroLibr.kepler(r1, v1tb, dtsec, out r3h, out v3h);
                    for (int j = 0; j < 3; j++)
                        dr[j] = r2[j] - r3h[j];
                    if (MathTimeLibr.mag(dr) > 0.01)
                        strbuild.AppendLine("velb does not get to r2 (km) position " + MathTimeLibr.mag(dr).ToString() + "\n");
                    //strbuild.AppendLine("diffs " + MathTimeLibr.mag(dr).ToString("0.00000000000"));

                    for (int j = 0; j < 3; j++)
                    {
                        dv1[j] = v1tk[j] - v1tb[j];
                        dv2[j] = v2tk[j] - v2tb[j];
                    }
                    if (MathTimeLibr.mag(dv1) > 0.01 || MathTimeLibr.mag(dv2) > 0.01)
                        strbuild.AppendLine("velk does not match velb \n");
                    //strbuild.AppendLine("diffs " + MathTimeLib::mag(dr).ToString("0.00000000000"));
                }
                //else
                //    strbuild.AppendLine(" ------------------------- not enough time for 1 revs ");

                strbuild.AppendLine(" TEST ------------------ L  L " + nnrev.ToString() + " rev ------------------");
                //if (dtsec / tusec >= tbidk[nnrev, 2])
                // do inside lambertk now
                {
                    dm = 'L';
                    de = 'L';
                    // switch tdi!!  tdidk to tdirk 'H'
                    AstroLambertkLibr.lambertkmins(s, tau, nnrev, dm, de, out kbi, out tof);  // 'H'

                    //double tofk1, kbik2, tofk2, kbik1;
                    //string outstr;
                    //getmins(1, 'k', nrev, r1, r2, s, tau, dm, de, out tofk1, out kbik1, out tofk2, out kbik2, out outstr);

                    AstroLambertkLibr.lambertK(r1, v1, r2, dm, de, nnrev, dtwait, dtsec, tof, kbi, numiter, altpadc, modecon, 'n',
                        out v1tk, out v2tk, out f, out g, out gdot, out hitearth, out errorout, out detailSum, out detailAll);
                    strbuildAll.AppendLine(detailAll);
                    //strbuild.AppendLine("k#" + caseopt + " " + detailSum + " diffs " + MathTimeLibr.mag(dr).ToString("0.00000000000"));
                    strbuild.AppendLine("lamk v1t " + " " + v1tk[0].ToString("0.00000000000") + " " + v1tk[1].ToString("0.00000000000") + " " + v1tk[2].ToString("0.00000000000"));
                    strbuild.AppendLine("lamk v2t " + " " + v2tk[0].ToString("0.00000000000") + " " + v2tk[1].ToString("0.00000000000") + " " + v2tk[2].ToString("0.00000000000"));
                    //strbuild.AppendLine(magv1t.ToString("0.0000000").PadLeft(12) + " " + magv2t.ToString("0.0000000").PadLeft(12));

                    AstroLibr.kepler(r1, v1tk, dtsec, out r3h, out v3h);
                    for (int j = 0; j < 3; j++)
                        dr[j] = r2[j] - r3h[j];
                    if (MathTimeLibr.mag(dr) > 0.01)
                        strbuild.AppendLine("velk does not get to r2 (km) position " + MathTimeLibr.mag(dr).ToString() + "\n");

                    AstroLibr.lambertumins(r1, r2, nnrev, dm, out kbi, out tof);
                    AstroLibr.lambertuniv(r1, r2, v1, dm, de, nnrev, dtsec, kbi, altpadc * AstroLibr.gravConst.re, 'n', out v1tu, out v2tu, out hitearth, out detailSum, out detailAll);
                    //strbuild.AppendLine(detailSum);
                    strbuild.AppendLine("univ v1t " + " " + v1tu[0].ToString("0.00000000000") + " " + v1tu[1].ToString("0.00000000000") + " " + v1tu[2].ToString("0.00000000000"));
                    strbuild.AppendLine("univ v2t " + " " + v2tu[0].ToString("0.00000000000") + " " + v2tu[1].ToString("0.00000000000") + " " + v2tu[2].ToString("0.00000000000"));
                    AstroLibr.kepler(r1, v1tu, dtsec, out r3h, out v3h);
                    for (int j = 0; j < 3; j++)
                        dr[j] = r2[j] - r3h[j];
                    if (MathTimeLibr.mag(dr) > 0.01)
                        strbuild.AppendLine("velu does not get to r2 (km) position " + MathTimeLibr.mag(dr).ToString() + "\n");

                    for (int j = 0; j < 3; j++)
                    {
                        dv1[j] = v1tk[j] - v1tu[j];
                        dv2[j] = v2tk[j] - v2tu[j];
                    }
                    if (MathTimeLibr.mag(dv1) > 0.01 || MathTimeLibr.mag(dv2) > 0.01)
                        strbuild.AppendLine("velk does not match velu \n");

                    AstroLibr.lambertbattin(r1, r2, v1, dm, de, nnrev, dtsec, altpadc * AstroLibr.gravConst.re, 'n', out v1tb, out v2tb, out hitearth, out detailSum, out detailAll);
                    //strbuild.AppendLine(detailSum);
                    strbuild.AppendLine("batt v1t " + " " + v1tb[0].ToString("0.00000000000") + " " + v1tb[1].ToString("0.00000000000") + " " + v1tb[2].ToString("0.00000000000"));
                    strbuild.AppendLine("batt v2t " + " " + v2tb[0].ToString("0.00000000000") + " " + v2tb[1].ToString("0.00000000000") + " " + v2tb[2].ToString("0.00000000000"));
                    AstroLibr.kepler(r1, v1tb, dtsec, out r3h, out v3h);
                    for (int j = 0; j < 3; j++)
                        dr[j] = r2[j] - r3h[j];
                    if (MathTimeLibr.mag(dr) > 0.01)
                        strbuild.AppendLine("velb does not get to r2 (km) position " + MathTimeLibr.mag(dr).ToString() + "\n");
                    //strbuild.AppendLine("diffs " + MathTimeLibr.mag(dr).ToString("0.00000000000"));

                    for (int j = 0; j < 3; j++)
                    {
                        dv1[j] = v1tk[j] - v1tb[j];
                        dv2[j] = v2tk[j] - v2tb[j];
                    }
                    if (MathTimeLibr.mag(dv1) > 0.01 || MathTimeLibr.mag(dv2) > 0.01)
                        strbuild.AppendLine("velk does not match velb \n");
                    //strbuild.AppendLine("diffs " + MathTimeLib::mag(dr).ToString("0.00000000000"));
                }
                //else
                //    strbuild.AppendLine(" ------------------------- not enough time for 1 revs ");

                strbuild.AppendLine(" TEST ------------------ S  H " + nnrev.ToString() + " rev ------------------");
                //if (dtsec / tusec >= tbirk[nnrev, 2])
                // do inside lambertk now
                {
                    dm = 'S';
                    de = 'H';
                    // switch tdi!!  tdirk to tdidk  'L'
                    AstroLambertkLibr.lambertkmins(s, tau, nnrev, dm, de, out kbi, out tof);  // 'L'
                    AstroLambertkLibr.lambertK(r1, v1, r2, dm, de, nnrev, dtwait, dtsec, tof, kbi, numiter, altpadc, modecon, 'n',
                        out v1tk, out v2tk, out f, out g, out gdot, out hitearth, out errorout, out detailSum, out detailAll);
                    strbuildAll.AppendLine(detailAll);
                    //strbuild.AppendLine("k#" + caseopt + " " + detailSum + " diffs " + MathTimeLibr.mag(dr).ToString("0.00000000000"));
                    strbuild.AppendLine("lamk v1t " + " " + v1tk[0].ToString("0.00000000000") + " " + v1tk[1].ToString("0.00000000000") + " " + v1tk[2].ToString("0.00000000000"));
                    strbuild.AppendLine("lamk v2t " + " " + v2tk[0].ToString("0.00000000000") + " " + v2tk[1].ToString("0.00000000000") + " " + v2tk[2].ToString("0.00000000000"));
                    //strbuild.AppendLine(magv1t.ToString("0.0000000").PadLeft(12) + " " + magv2t.ToString("0.0000000").PadLeft(12));

                    AstroLibr.kepler(r1, v1tk, dtsec, out r3h, out v3h);
                    for (int j = 0; j < 3; j++)
                        dr[j] = r2[j] - r3h[j];
                    if (MathTimeLibr.mag(dr) > 0.01)
                        strbuild.AppendLine("velk does not get to r2 (km) position " + MathTimeLibr.mag(dr).ToString() + "\n");

                    AstroLibr.lambertumins(r1, r2, nnrev, dm, out kbi, out tof);
                    AstroLibr.lambertuniv(r1, r2, v1, dm, de, nnrev, dtsec, kbi, altpadc * AstroLibr.gravConst.re, 'n', out v1tu, out v2tu, out hitearth, out detailSum, out detailAll);
                    //strbuild.AppendLine(detailSum);
                    strbuild.AppendLine("univ v1t " + " " + v1tu[0].ToString("0.00000000000") + " " + v1tu[1].ToString("0.00000000000") + " " + v1tu[2].ToString("0.00000000000"));
                    strbuild.AppendLine("univ v2t " + " " + v2tu[0].ToString("0.00000000000") + " " + v2tu[1].ToString("0.00000000000") + " " + v2tu[2].ToString("0.00000000000"));
                    AstroLibr.kepler(r1, v1tu, dtsec, out r3h, out v3h);
                    for (int j = 0; j < 3; j++)
                        dr[j] = r2[j] - r3h[j];
                    if (MathTimeLibr.mag(dr) > 0.01)
                        strbuild.AppendLine("velu does not get to r2 (km) position " + MathTimeLibr.mag(dr).ToString() + "\n");

                    for (int j = 0; j < 3; j++)
                    {
                        dv1[j] = v1tk[j] - v1tu[j];
                        dv2[j] = v2tk[j] - v2tu[j];
                    }
                    if (MathTimeLibr.mag(dv1) > 0.01 || MathTimeLibr.mag(dv2) > 0.01)
                        strbuild.AppendLine("velk does not match velu \n");

                    AstroLibr.lambertbattin(r1, r2, v1, dm, de, nnrev, dtsec, altpadc * AstroLibr.gravConst.re, 'n', out v1tb, out v2tb, out hitearth, out detailSum, out detailAll);
                    //strbuild.AppendLine(detailSum);
                    strbuild.AppendLine("batt v1t " + " " + v1tb[0].ToString("0.00000000000") + " " + v1tb[1].ToString("0.00000000000") + " " + v1tb[2].ToString("0.00000000000"));
                    strbuild.AppendLine("batt v2t " + " " + v2tb[0].ToString("0.00000000000") + " " + v2tb[1].ToString("0.00000000000") + " " + v2tb[2].ToString("0.00000000000"));
                    AstroLibr.kepler(r1, v1tb, dtsec, out r3h, out v3h);
                    for (int j = 0; j < 3; j++)
                        dr[j] = r2[j] - r3h[j];
                    if (MathTimeLibr.mag(dr) > 0.01)
                        strbuild.AppendLine("velb does not get to r2 (km) position " + MathTimeLibr.mag(dr).ToString() + "\n");
                    //strbuild.AppendLine("diffs " + MathTimeLibr.mag(dr).ToString("0.00000000000"));

                    for (int j = 0; j < 3; j++)
                    {
                        dv1[j] = v1tk[j] - v1tb[j];
                        dv2[j] = v2tk[j] - v2tb[j];
                    }
                    if (MathTimeLibr.mag(dv1) > 0.01 || MathTimeLibr.mag(dv2) > 0.01)
                        strbuild.AppendLine("velk does not match velb \n");
                    //strbuild.AppendLine("diffs " + MathTimeLib::mag(dr).ToString("0.00000000000"));
                }
                //else
                //    strbuild.AppendLine(" ------------------------- not enough time for 1 revs ");

                strbuild.AppendLine(" TEST ------------------ L  H " + nnrev.ToString() + " rev ------------------");
                //if (dtsec / tusec >= tbirk[nnrev, 2])
                // do inside lambertk now
                {
                    dm = 'L';
                    de = 'H';
                    AstroLambertkLibr.lambertkmins(s, tau, nnrev, dm, de, out kbi, out tof);
                    AstroLambertkLibr.lambertK(r1, v1, r2, dm, de, nnrev, dtwait, dtsec, tof, kbi, numiter, altpadc, modecon, 'n',
                        out v1tk, out v2tk, out f, out g, out gdot, out hitearth, out errorout, out detailSum, out detailAll);
                    strbuildAll.AppendLine(detailAll);
                    //strbuild.AppendLine("k#" + caseopt + " " + detailSum + " diffs " + MathTimeLibr.mag(dr).ToString("0.00000000000"));
                    strbuild.AppendLine("lamk v1t " + " " + v1tk[0].ToString("0.00000000000") + " " + v1tk[1].ToString("0.00000000000") + " " + v1tk[2].ToString("0.00000000000"));
                    strbuild.AppendLine("lamk v2t " + " " + v2tk[0].ToString("0.00000000000") + " " + v2tk[1].ToString("0.00000000000") + " " + v2tk[2].ToString("0.00000000000"));
                    //strbuild.AppendLine(magv1t.ToString("0.0000000").PadLeft(12) + " " + magv2t.ToString("0.0000000").PadLeft(12));

                    AstroLibr.kepler(r1, v1tk, dtsec, out r3h, out v3h);
                    for (int j = 0; j < 3; j++)
                        dr[j] = r2[j] - r3h[j];
                    if (MathTimeLibr.mag(dr) > 0.01)
                        strbuild.AppendLine("velk does not get to r2 (km) position " + MathTimeLibr.mag(dr).ToString() + "\n");

                    AstroLibr.lambertumins(r1, r2, nnrev, dm, out kbi, out tof);
                    AstroLibr.lambertuniv(r1, r2, v1, dm, de, nnrev, dtsec, kbi, altpadc * AstroLibr.gravConst.re, 'n', out v1tu, out v2tu, out hitearth, out detailSum, out detailAll);
                    //strbuild.AppendLine(detailSum);
                    strbuild.AppendLine("univ v1t " + " " + v1tu[0].ToString("0.00000000000") + " " + v1tu[1].ToString("0.00000000000") + " " + v1tu[2].ToString("0.00000000000"));
                    strbuild.AppendLine("univ v2t " + " " + v2tu[0].ToString("0.00000000000") + " " + v2tu[1].ToString("0.00000000000") + " " + v2tu[2].ToString("0.00000000000"));
                    AstroLibr.kepler(r1, v1tu, dtsec, out r3h, out v3h);
                    for (int j = 0; j < 3; j++)
                        dr[j] = r2[j] - r3h[j];
                    if (MathTimeLibr.mag(dr) > 0.01)
                        strbuild.AppendLine("velu does not get to r2 (km) position " + MathTimeLibr.mag(dr).ToString() + "\n");

                    for (int j = 0; j < 3; j++)
                    {
                        dv1[j] = v1tk[j] - v1tu[j];
                        dv2[j] = v2tk[j] - v2tu[j];
                    }
                    if (MathTimeLibr.mag(dv1) > 0.01 || MathTimeLibr.mag(dv2) > 0.01)
                        strbuild.AppendLine("velk does not match velu \n");

                    AstroLibr.lambertbattin(r1, r2, v1, dm, de, nnrev, dtsec, altpadc * AstroLibr.gravConst.re, 'n', out v1tb, out v2tb, out hitearth, out detailSum, out detailAll);
                    //strbuild.AppendLine(detailSum);
                    strbuild.AppendLine("batt v1t " + " " + v1tb[0].ToString("0.00000000000") + " " + v1tb[1].ToString("0.00000000000") + " " + v1tb[2].ToString("0.00000000000"));
                    strbuild.AppendLine("batt v2t " + " " + v2tb[0].ToString("0.00000000000") + " " + v2tb[1].ToString("0.00000000000") + " " + v2tb[2].ToString("0.00000000000"));
                    AstroLibr.kepler(r1, v1tb, dtsec, out r3h, out v3h);
                    for (int j = 0; j < 3; j++)
                        dr[j] = r2[j] - r3h[j];
                    if (MathTimeLibr.mag(dr) > 0.01)
                        strbuild.AppendLine("velb does not get to r2 (km) position " + MathTimeLibr.mag(dr).ToString() + "\n");
                    //strbuild.AppendLine("diffs " + MathTimeLibr.mag(dr).ToString("0.00000000000"));

                    for (int j = 0; j < 3; j++)
                    {
                        dv1[j] = v1tk[j] - v1tb[j];
                        dv2[j] = v2tk[j] - v2tb[j];
                    }
                    if (MathTimeLibr.mag(dv1) > 0.01 || MathTimeLibr.mag(dv2) > 0.01)
                        strbuild.AppendLine("velk does not match velb \n");
                    //strbuild.AppendLine("diffs " + MathTimeLib::mag(dr).ToString("0.00000000000"));
                }
                //else
                //    strbuild.AppendLine(" ------------------------- not enough time for 1 revs ");

                strbuild.AppendLine(" --------------------------------end case " + caseopt + "------------------------------------------------ ");
                string resultStr = strbuildAll.ToString();
            }

            string directory = @"D:\Codes\LIBRARY\cs\TestAll\";
            File.WriteAllText(directory + "testall-lambertknown.out", strbuildAll.ToString());

            this.opsStatus.Text = "Done ";
            Refresh();
        }  // testknowncases


        public void testradecgeo2azel()
        {
            double rad = 180.0 / Math.PI;
            double rtasc, decl, rr, latgd, lon, alt, az, el;
            double ttt, jdut1, lod, xp, yp, ddpsi, ddeps;

            rtasc = 294.98914583 / rad;
            decl = -20.8234944 / rad;
            xp = 0.0;
            yp = 0.0;
            lod = 0.0;
            jdut1 = 2453101.82740678310;
            ttt = 0.042623631889;
            ddpsi = -0.052195;
            ddeps = -0.003875;
            rr = 12373.3546098;  // km
            latgd = 39.007 / rad;
            lon = -104.883 / rad;
            alt = 0.3253;

            AstroLibr.radecgeo2azel(rtasc, decl, rr, latgd, lon, alt, ttt, jdut1, lod, xp, yp, ddpsi, ddeps, out az, out el);
        }

        public void testijk2ll()
        {
            double[] r = new double[3];
            double latgc, latgd, lon, hellp, rad;
            rad = 180.0 / Math.PI;

            r = new double[] { 1.023 * AstroLibr.gravConst.re, 1.076 * AstroLibr.gravConst.re, 1.011 * AstroLibr.gravConst.re };

            AstroLibr.ecef2ll(r, out latgc, out latgd, out lon, out hellp);

            strbuild.AppendLine("ecef2ll " + (latgd * rad).ToString(fmt) + " " +
                (lon * rad).ToString(fmt) + " " + hellp.ToString(fmt) + "\n");

            AstroLibr.ecef2llb(r, out latgc, out latgd, out lon, out hellp);

            strbuild.AppendLine("ecef2llb " + (latgd * rad).ToString(fmt) + " " +
                (lon * rad).ToString(fmt) + " " + hellp.ToString(fmt) + "\n");
        }

        public void testgd2gc()
        {
            double rad = 180.0 / Math.PI;
            double latgd, ans;
            latgd = 34.173429 / rad;

            ans = AstroLibr.gd2gc(latgd);

            strbuild.AppendLine("gd2gc " + ans.ToString(fmt) + "\n");
        }

        public void testsite()
        {
            double rad = 180.0 / Math.PI;
            double latgd, lon, alt;
            double[] rsecef;
            double[] vsecef;
            latgd = 39.007 / rad;
            lon = -104.883 / rad;
            alt = 0.3253;

            AstroLibr.site(latgd, lon, alt, out rsecef, out vsecef);

            strbuild.AppendLine("site " + rsecef[0].ToString(fmt) + " " + rsecef[1].ToString(fmt) + " " + rsecef[2].ToString(fmt) + " " +
                        vsecef[0].ToString(fmt) + " " + vsecef[1].ToString(fmt) + " " + vsecef[2].ToString(fmt));
        }


        // --------  sun          - analytical sun ephemeris
        public void testsun()
        {
            double jd, rtasc, decl;
            double[] rsun;
            jd = 2449444.5;
            AstroLibr.sun(jd, out rsun, out rtasc, out decl);

            strbuild.AppendLine("sun " + rsun[0].ToString(fmt) + " " + rsun[1].ToString(fmt) + " " + rsun[2].ToString(fmt));
        }

        // --------  moon         - analytical moon ephemeris
        public void testmoon()
        {
            double jd, rtasc, decl;
            double[] rmoon;

            jd = 2449470.5;
            AstroLibr.moon(jd, out rmoon, out rtasc, out decl);

            strbuild.AppendLine("moon " + rmoon[0].ToString(fmt) + " " + rmoon[1].ToString(fmt) + " " + rmoon[2].ToString(fmt));
        }

        public void testkepler()
        {
            double[] r1;
            double[] v1;
            double[] r2;
            double[] v2;
            double dtsec;
            dtsec = 42397.344;  // s

            r1 = new double[] { 2.500000 * AstroLibr.gravConst.re, 0.000000, 0.000000 };
            // assume circular initial orbit for vel calcs
            v1 = new double[] { 0.0, Math.Sqrt(AstroLibr.gravConst.mu / r1[0]), 0.0 };
            strbuild.AppendLine("kepler " + r1[0].ToString(fmt) + " " + r1[1].ToString(fmt) + " " + r1[2].ToString(fmt) + " " +
                                v1[0].ToString(fmt) + " " + v1[1].ToString(fmt) + " " + v1[2].ToString(fmt) + " " + dtsec);

            AstroLibr.kepler(r1, v1, dtsec, out r2, out v2);

            strbuild.AppendLine("kepler " + r2[0].ToString(fmt) + " " + r2[1].ToString(fmt) + " " + r2[2].ToString(fmt) + " " +
                                v2[0].ToString(fmt) + " " + v2[1].ToString(fmt) + " " + v2[2].ToString(fmt));

            // test multi-rev case
            double p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper, period;
            AstroLibr.rv2coe(r1, v1, out p, out a, out ecc, out incl, out raan, out argp, out nu, out m, out arglat, out truelon, out lonper);
            period = 2.0 * Math.PI * Math.Sqrt(Math.Pow(MathTimeLibr.mag(r1), 3) / AstroLibr.gravConst.mu);

            AstroLibr.kepler(r1, v1, dtsec + 7.0 * period, out r2, out v2);

            strbuild.AppendLine("kepler " + r2[0].ToString(fmt) + " " + r2[1].ToString(fmt) + " " + r2[2].ToString(fmt) + " " +
                                v2[0].ToString(fmt) + " " + v2[1].ToString(fmt) + " " + v2[2].ToString(fmt) + " " + (dtsec + 7.0 * period));

        }



        // test in geoloc.sln
        //public void testcovct2cl()
        //{
        //    double[,] cartcov = new double[6, 6];
        //    double[] cartstate = new double[6];
        //    string anomclass;
        //    double[,] classcov = new double[6, 6];
        //    double[,] tm = new double[6, 6];

        //    MathTimeLibr.covct2cl(cartcov, cartstate, anomclass, out classcov, out tm);

        //}
        //public void testcovcl2ct()
        //{
        //    MathTimeLibr.covcl2ct
        //    (double[,] classcov, double[] classstate, string anomclass, out double[,] cartcov, out double[,] tm
        //            );
        //}
        //public void testcovct2eq()
        //{
        //    double[] classState = new double[6];
        //    double[] cartState = new double[6];
        //    double[] eqState = new double[6];
        //    double[] flState = new double[6];
        //    double[,] cartCov = new double[6, 6];
        //    double[,] classCov = new double[6, 6];
        //    double[,] eqCov = new double[6, 6];
        //    double[,] flCov = new double[6, 6];
        //    double[,] rswCov = new double[6, 6];
        //    double[,] ntwCov = new double[6, 6];
        //    double[,] tm = new double[6, 6];

        //    cartCov = new double[,] { { 1, 0, 0, 0, 0, 0 }, { 0, 1, 0, 0, 0, 0 }, { 0, 0, 1, 0, 0, 0 },
        //                         { 0, 0, 0, 1, 0, 0 }, { 0, 0, 0, 0, 1, 0 }, { 0, 0, 0, 0, 0, 1 } };


        //    MathTimeLibr.covct2eq
        //    (     double[,] cartcov, double[] cartstate, string anomeq, Int16 fr, out double[,] eqcov, out  tm                );
        //}
        //public void testcoveq2ct()
        //{
        //    MathTimeLibr.coveq2ct
        //     (                double[,] eqcov, double[] eqstate, string anomeq, Int16 fr, out double[,] cartcov, out  tm                );
        //}
        //public void testcovcl2eq()
        //{
        //    MathTimeLibr.covcl2eq
        //    (
        //            double[,] classcov, double[] classstate, string anomclass, string anomeq, Int16 fr, out double[,] eqcov, out  tm                );
        //}
        //public void testcoveq2cl()
        //{
        //    MathTimeLibr.coveq2cl(double[,] eqcov, double[] eqstate, string anomeq, string anomclass, Int16 fr, out double[,] classcov, out  tm);
        //}
        //public void testcovct2fl()
        //{
        //    MathTimeLibr.covct2fl
        //      (
        //            double[,] cartcov, double[] cartstate, string anomflt, double ttt, double jdut1, double lod,
        //            double xp, double yp, Int16 terms, double ddpsi, double ddeps, out double[,] flcov, out  tm
        //            );
        //}
        //public void testcovfl2ct()
        //{
        //    MathTimeLibr.covfl2ct(double[,] flcov, double[] flstate, string anomflt, double ttt, double jdut1, double lod,
        //            double xp, double yp, Int16 terms, double ddpsi, double ddeps, out double[,] cartcov, out  tm);

        //}
        //public void testcovct_rsw()
        //{
        //    MathTimeLibr.covct_rsw(ref double[,] cartcov, double[] cartstate, MathTimeLib.Edirection direct, ref double[,] rswcov, out  tm);
        //        direct = MathTimeLib.Edirection.eto;
        //        MathTimelibr.covct_ntw(ref cartCovo, cartState, direct, ref ntwCov, out tm);

        //    }
        //    public void testcovct_ntw()
        //{
        //        direct = MathTimeLib.Edirection.eto;
        //        MathTimelibr.covct_ntw(ref cartCovo, cartState, direct, ref ntwCov, out tm);

        //        MathTimeLibr.covct_ntw(ref double[,] cartcov, double[] cartstate, MathTimeLib.Edirection direct, ref double[,] ntwcov, out  tm);
        //}

        public void testsunmoonjpl()
        {
            AstroLib.jpldedataClass[] jpldearr = AstroLibr.jpldearr;
            double[] rsun = new double[3];
            double[] rmoon = new double[3];
            double rtascs, decls, rtascm, declm, rsmag, rmmag,  jd, jdF;

            MathTimeLibr.jday(2017, 5, 11, 3, 51, 42.7657, out jd, out jdF);

            strbuild.AppendLine(" =============================   test sun and moon ephemerides =============================\n");

            // read in jpl sun moon files
            // answers
            strbuild.AppendLine("2017  5 11  0  96576094.2145 106598001.2476 46210616.7776     151081093.9419  0.9804616 -252296.5509 -302841.7334 -93212.7720");
            strbuild.AppendLine("2017  5 11 12  95604355.9737 107353047.2919 46537942.1006     151098145.9151  0.9802403 -218443.5158 -325897.7785 -102799.8515");
            strbuild.AppendLine("2017  5 12  0  94625783.6875 108100430.4112 46861940.2387     151115133.0492  0.9800199 -182165.5046 -345316.4032 -111246.7742");

            // for 1 day centers, need to adjust the initjpl function
            //AstroLibr.initjplde(ref jpldearr, "D:/Codes/LIBRARY/DataLib/", "sunmooneph_430t.txt", out jdjpldestart, out jdjpldestartFrac);
            AstroLibr.readjplde(ref jpldearr, "D:/Codes/LIBRARY/DataLib/", "sunmooneph_430t12.txt");

            AstroLibr.findjpldeparam(jd, 0.0, 'l', jpldearr, out rsun, out rsmag, out rmoon, out rmmag);
            strbuild.AppendLine("findjpldeephem 0000 hrs l\n " + jd.ToString() + " 0.00000 " + rsun[0].ToString() + " " +
                rsun[1].ToString() + " " + rsun[2].ToString() + " " +
                rmoon[0].ToString() + " " + rmoon[1].ToString() + " " + rmoon[2].ToString());

            AstroLibr.findjpldeparam(jd, 0.0, 's', jpldearr, out rsun, out rsmag, out rmoon, out rmmag);
            strbuild.AppendLine("findjpldeephem 0000 hrs s\n " + jd.ToString() + " " + jdF.ToString() + " " + rsun[0].ToString() + " " +
                rsun[1].ToString() + " " + rsun[2].ToString() + " " +
                rmoon[0].ToString() + " " + rmoon[1].ToString() + " " + rmoon[2].ToString());

            AstroLibr.sunmoonjpl(jd, 0.0, 's', ref jpldearr, out rsun, out rtascs, out decls, out rmoon, out rtascm, out declm);
            strbuild.AppendLine("sunmoon 0000 hrs s\n " + jd.ToString() + " " + jdF.ToString() + " " + rsun[0].ToString() + " " +
                rsun[1].ToString() + " " + rsun[2].ToString() + " " +
                rmoon[0].ToString() + " " + rmoon[1].ToString() + " " + rmoon[2].ToString());


            AstroLibr.findjpldeparam(jd, jdF, 'l', jpldearr, out rsun, out rsmag, out rmoon, out rmmag);
            strbuild.AppendLine("findjpldeephem hrs l\n " + jd.ToString() + " " + jdF.ToString() + " " + rsun[0].ToString() + " " +
                rsun[1].ToString() + " " + rsun[2].ToString() + " " +
                rmoon[0].ToString() + " " + rmoon[1].ToString() + " " + rmoon[2].ToString());

            AstroLibr.sunmoonjpl(jd, jdF, 'l', ref jpldearr, out rsun, out rtascs, out decls, out rmoon, out rtascm, out declm);
            strbuild.AppendLine("sunmoon hrs l\n " + jd.ToString() + " " + jdF.ToString() + " " + rsun[0].ToString() + " " +
                rsun[1].ToString() + " " + rsun[2].ToString() + " " +
                rmoon[0].ToString() + " " + rmoon[1].ToString() + " " + rmoon[2].ToString());

            AstroLibr.findjpldeparam(jd, 1.0, 'l', jpldearr, out rsun, out rsmag, out rmoon, out rmmag);
            strbuild.AppendLine("findjpldeephem 2400 hrs s\n " + jd.ToString() + " " + jdF.ToString() + " " + rsun[0].ToString() + " " +
                rsun[1].ToString() + " " + rsun[2].ToString() + " " +
                rmoon[0].ToString() + " " + rmoon[1].ToString() + " " + rmoon[2].ToString());


            // ex 8.5 test
            MathTimeLibr.jday(2020, 2, 18, 15, 8, 47.23847, out jd, out jdF);
            AstroLibr.findjpldeparam(jd, 0.0, 's', jpldearr, out rsun, out rsmag, out rmoon, out rmmag);
            strbuild.AppendLine("ex findjpldeephem 0000 hrs s\n " + jd.ToString() + " " + jdF.ToString() + " " + rsun[0].ToString() + " " +
                rsun[1].ToString() + " " + rsun[2].ToString() + " " +
                rmoon[0].ToString() + " " + rmoon[1].ToString() + " " + rmoon[2].ToString());

            // test interpolation of vectors
            // shows spline is MUCH better - 3 km sun variation in mid day linear, 60m diff with spline. 
            MathTimeLibr.jday(2017, 5, 11, 3, 51, 42.7657, out jd, out jdF);
            MathTimeLibr.jday(2000, 1, 1, 0, 0, 0.0, out jd, out jdF);
            strbuild.AppendLine("findjplde  mfme     rsun x             y                 z             rmoon x             y                z      (km)");

            var watch = System.Diagnostics.Stopwatch.StartNew();
            // the code that you want to measure comes here
            // read in jpl sun moon files - seems to be the slowest part (800 msec)
            AstroLibr.readjplde(ref jpldearr, "D:/Codes/LIBRARY/DataLib/", "sunmooneph_430t12.txt");

            watch.Stop();
            var elapsedMs = watch.ElapsedMilliseconds;

            watch = System.Diagnostics.Stopwatch.StartNew();

            int ii;
            for (ii = 0; ii < 36500; ii++)
            {
                // seems pretty fast (45 msec)
                for (int jj = 0; jj < 24; jj++)
                {
                    AstroLibr.findjpldeparam(jd + ii, (jj * 1.0) / 24.0, 's', jpldearr, out rsun, out rsmag, out rmoon, out rmmag);
                    // the write takes some time (160 msec)
                    //strbuild.AppendLine(" " + jd.ToString() + " " + (ii * 60.0).ToString("0000") + " " +
                    //    rsun[0].ToString() + " " + rsun[1].ToString() + " " + rsun[2].ToString() + " " +
                    //    rmoon[0].ToString() + " " + rmoon[1].ToString() + " " + rmoon[2].ToString());
                }
            }

            watch.Stop();
            elapsedMs = watch.ElapsedMilliseconds;
        }

        public void testkp2ap()
        {
            int i;
            double kp, ap;

            for (i = 1; i <= 27; i++)
            {
                kp = 10.0 * i / 3.0;
                ap = EOPSPWLibr.kp2ap(kp);
                // get spacing correct, leading 0, front spaces
                strbuild.AppendLine(i.ToString("##") + (0.1 * kp).ToString("  0.######") + ap.ToString("  0.######"));
            }
        }


        public void testazel2radec()
        {
            double rad = 180.0 / Math.PI;
            double[] reci = new double[3];
            double[] veci = new double[3];
            double[] recef = new double[3];
            double[] vecef = new double[3];
            double[] rsecef = new double[3];
            double[] vsecef = new double[3];
            double[] rseci = new double[3];
            double[] vseci = new double[3];
            double rho, az, el, drho, daz, del, alt, latgd, lon;
            double rr, rtasc, decl, drr, drtasc, ddecl;
            double trtasc, tdecl, dtrtasc, dtdecl;
            double ttt, xp, yp, lod, jdut1, ddpsi, ddeps, dut1, lst, gst, jdtt, jdftt;
            int year, mon, day, hr, minute, dat;
            double second, jd, jdFrac;

            rr = trtasc = tdecl = rtasc = decl = drr = dtrtasc = dtdecl = drtasc = ddecl = 0.0;
            rho = az = el = drho = daz = del = 0.0;

            EOPSPWLib.iau80Class iau80arr;
            string fileLoc;
            fileLoc = @"D:\Codes\LIBRARY\DataLib\nut80.dat";
            EOPSPWLibr.iau80in(fileLoc, out iau80arr);

            xp = 0.0;
            yp = 0.0;
            lod = 0.0;
            ddpsi = -0.052195;
            ddeps = -0.003875;
            dut1 = -0.37816;

            year = 2015;
            mon = 12;
            day = 15;
            hr = 16;
            dat = 36;
            minute = 58;
            second = 50.208;
            MathTimeLibr.jday(year, mon, day, hr, minute, second, out jd, out jdFrac);

            // note you have to use tdb for time of ineterst AND j2000 (when dat = 32)
            jdtt = jd;
            jdftt = jdFrac + (dat + 32.184) / 86400.0;
            ttt = (jdtt + jdftt - 2451545.0) / 36525.0;
            jdut1 = jd + jdFrac + dut1 / 86400.0;

            recef = new double[] { -605.79221660, -5870.22951108, 3493.05319896 };
            //recef = new double[] { -100605.79221660, -1005870.22951108, 1003493.05319896 };
            vecef = new double[] { -1.56825429, -3.70234891, -6.47948395 };
            AstroLibr.eci_ecef(ref reci, ref veci, MathTimeLib.Edirection.efrom, ref recef, ref vecef,
                iau80arr, jdtt, jdftt, jdut1, lod, xp, yp, ddpsi, ddeps);

            lon = -104.883 / rad;
            latgd = 39.007 / rad;
            alt = 2.102;
            AstroLibr.site(latgd, lon, alt, out rsecef, out vsecef);

            AstroLibr.eci_ecef(ref rseci, ref vseci, MathTimeLib.Edirection.efrom, ref rsecef, ref vsecef,
                 iau80arr, jdtt, jdftt, jdut1, lod, xp, yp, ddpsi, ddeps);

            AstroLibr.lstime(lon, jdut1, out lst, out gst);

            // print out initial conditions
            strbuild.AppendLine("recef  " + recef[0].ToString(fmt).PadLeft(4) + " " + recef[1].ToString(fmt).PadLeft(4) + " " + recef[2].ToString(fmt).PadLeft(4) + " " +
                                "v  " + vecef[0].ToString(fmt).PadLeft(4) + " " + vecef[1].ToString(fmt).PadLeft(4) + " " + vecef[2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("rs ecef  " + rsecef[0].ToString(fmt).PadLeft(4) + " " + rsecef[1].ToString(fmt).PadLeft(4) + " " + rsecef[2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("reci  " + reci[0].ToString(fmt).PadLeft(4) + " " + reci[1].ToString(fmt).PadLeft(4) + " " + reci[2].ToString(fmt).PadLeft(4) + " " +
                                "v  " + veci[0].ToString(fmt).PadLeft(4) + " " + veci[1].ToString(fmt).PadLeft(4) + " " + veci[2].ToString(fmt).PadLeft(4));
            strbuild.AppendLine("rs eci  " + rseci[0].ToString(fmt).PadLeft(4) + " " + rseci[1].ToString(fmt).PadLeft(4) + " " + rseci[2].ToString(fmt).PadLeft(4));


            AstroLibr.rv_razel(ref recef, ref vecef, latgd, lon, alt, MathTimeLib.Edirection.eto, ref rho, ref az, ref el, ref drho, ref daz, ref del);

            AstroLibr.rv_radec(ref reci, ref veci, MathTimeLib.Edirection.eto, ref rr, ref rtasc, ref decl, ref drr, ref drtasc, ref ddecl);

            AstroLibr.rv_tradec(ref reci, ref veci, rseci, MathTimeLib.Edirection.eto, ref rho, ref trtasc, ref tdecl, ref drho, ref dtrtasc, ref dtdecl);

            // print out results
            strbuild.AppendLine("razel " + rho.ToString(fmt).PadLeft(4) + " " + az.ToString(fmt).PadLeft(4) + " " + el.ToString(fmt).PadLeft(4) + " " +
                                "  " + drho.ToString(fmt).PadLeft(4) + " " + daz.ToString(fmt).PadLeft(4) + " " + del.ToString(fmt).PadLeft(4));
            strbuild.AppendLine("radec " + rr.ToString(fmt).PadLeft(4) + " " + rtasc.ToString(fmt).PadLeft(4) + " " + decl.ToString(fmt).PadLeft(4) + " " +
                                "  " + drr.ToString(fmt).PadLeft(4) + " " + drtasc.ToString(fmt).PadLeft(4) + " " + ddecl.ToString(fmt).PadLeft(4));
            strbuild.AppendLine("tradec " + rho.ToString(fmt).PadLeft(4) + " " + trtasc.ToString(fmt).PadLeft(4) + " " + tdecl.ToString(fmt).PadLeft(4) + " " +
                                drho.ToString(fmt).PadLeft(4) + " " + dtrtasc.ToString(fmt).PadLeft(4) + " " + dtdecl.ToString(fmt).PadLeft(4));

            double rtasc1;
            AstroLibr.azel_radec(az, el, lst, latgd, out rtasc, out decl, out rtasc1);
            strbuild.AppendLine("radec " + rtasc.ToString(fmt).PadLeft(4) + " rtasc1 " + rtasc1.ToString(fmt).PadLeft(4) + " " + decl.ToString(fmt).PadLeft(4));
            strbuild.AppendLine("radec " + (Math.PI * 2 - rtasc).ToString(fmt).PadLeft(4) + " rtasc1 " + (Math.PI * 2 - rtasc1).ToString(fmt).PadLeft(4) + " " + decl.ToString(fmt).PadLeft(4));
        }


        /* ------------------------------------------------------------------------------
         *
         *                           function LegPolyEx
         *
         *   this function finds the exact (from equations) Legendre polynomials for the gravity field. 
         *   note that the arrays are indexed from 0 to coincide with the usual nomenclature (eq 8-21 
         *   in my text). fortran implementations will have indicies of 1 greater as they often 
         *   start at 1. these are exact expressions derived from mathematica. 
         *      
         *  author        : david vallado             davallado@gmail.com  16 dec 2019
         *
         *  inputs        description                                   range / units
         *    latgc       - Geocentric lat of satellite                   pi to pi rad          
         *    order       - size of gravity field                         1..2160..
         *
         *  outputs       :
         *    LegArr      - [,] array of Legendre polynomials
         *
         *  locals :
         *    L,m         - degree and order indices
         *    conv        - conversion to un-normalize
         *
         *  coupling      :
         *   none
         *
         *  references :
         *    vallado       2013, 597, Eq 8-57
          ------------------------------------------------------------------------------*/

        public void LegPolyEx
           (
           double latgc,
           Int32 order,
           out double[,] LegArrEx
           )
        {
            LegArrEx = new double[order + 2, order + 2];

            double s = LegArrEx[1, 0] = Math.Sin(latgc);
            double c = LegArrEx[1, 1] = Math.Cos(latgc);

            // -------------------- exact epxressions ---------------------- }
            LegArrEx[2, 0] = 0.5 * (3 * s * s - 1.0);
            LegArrEx[2, 1] = 3.0 * s * c;
            LegArrEx[2, 2] = 3.0 * c * c;

            // include (-1)^m for all the terms
            LegArrEx[3, 0] = -0.5 * s * (3 - 5 * s * s);
            LegArrEx[3, 1] = (3.0 / 2) * c * (-1 + 5 * s * s); // 15s^2 - 3
            LegArrEx[3, 2] = 15 * s * c * c;
            LegArrEx[3, 3] = 15 * (c * c * c);

            LegArrEx[4, 0] = 1.0 / 8.0 * (35.0 * Math.Pow(s, 4) - 30.0 * Math.Pow(s, 2) + 3.0);
            LegArrEx[4, 1] = 2.5 * c * (-3 * s + 7 * Math.Pow(s, 3));
            LegArrEx[4, 2] = 7.5 * Math.Pow(c, 2) * (-1 + 7 * Math.Pow(s, 2));
            LegArrEx[4, 3] = 105.0 * Math.Pow(c, 3) * s;
            LegArrEx[4, 4] = 105.0 * Math.Pow(c, 4);

            LegArrEx[5, 0] = (1.0 / 8) * s * (15 - 70 * Math.Pow(s, 2) + 63 * Math.Pow(s, 4));
            LegArrEx[5, 1] = (15.0 / 8) * c * (1 - 14 * Math.Pow(s, 2) + 21 * Math.Pow(s, 4));
            LegArrEx[5, 2] = (105.0 / 2) * c * c * (-s + 3 * Math.Pow(s, 3));
            LegArrEx[5, 3] = (105.0 / 2) * Math.Pow(c, 3) * (-1 + 9 * Math.Pow(s, 2));
            LegArrEx[5, 4] = 945.0 * s * Math.Pow(c, 4);
            LegArrEx[5, 5] = 945.0 * Math.Pow(c, 5);

            LegArrEx[6, 0] = 1.0 / 16 * (-5 + 105 * Math.Pow(s, 2) - 315 * Math.Pow(s, 4) + 231 * Math.Pow(s, 6));
            LegArrEx[6, 1] = (21.0 / 8) * c * (5 * s - 30 * Math.Pow(s, 3) + 33 * Math.Pow(s, 5));
            LegArrEx[6, 2] = (105.0 / 8) * c * c * (1 - 18 * Math.Pow(s, 2) + 33 * Math.Pow(s, 4));
            LegArrEx[6, 3] = (315.0 / 2) * Math.Pow(c, 3) * (-3 * s + 11 * Math.Pow(s, 3));
            LegArrEx[6, 4] = 945.0 / 2 * Math.Pow(c * c, 2) * (-1 + 11 * Math.Pow(s, 2));
            LegArrEx[6, 5] = 10395.0 * s * Math.Pow(c, 5);
            LegArrEx[6, 6] = 10395.0 * Math.Pow(c * c, 3);

            LegArrEx[7, 0] = 1.0 / 16 * (-35 * s + 315 * Math.Pow(s, 3) - 693 * Math.Pow(s, 5) + 429 * Math.Pow(s, 7));
            LegArrEx[7, 1] = (7.0 / 16) * c * (-5 + 135 * Math.Pow(s, 2) - 495 * Math.Pow(s, 4) + 429 * Math.Pow(s, 6));
            LegArrEx[7, 2] = (63.0 / 8) * c * c * (15 * s - 110 * Math.Pow(s, 3) + 143 * Math.Pow(s, 5));
            LegArrEx[7, 3] = (315.0 / 8) * Math.Pow(c, 3) * (3 - 66 * Math.Pow(s, 2) + 143 * Math.Pow(s, 4));
            LegArrEx[7, 4] = 3465.0 / 2 * Math.Pow(c * c, 2) * (-3 * s + 13 * Math.Pow(s, 3));
            LegArrEx[7, 5] = (10395.0 / 2) * Math.Pow(c, 5) * (-1 + 13 * Math.Pow(s, 2));
            LegArrEx[7, 6] = 135135.0 * s * Math.Pow(c * c, 3);
            LegArrEx[7, 7] = 135135.0 * Math.Pow(c, 7);

            LegArrEx[8, 0] = 1.0 / 128 * (35 - 1260 * Math.Pow(s, 2) + 6930 * Math.Pow(s, 4) - 12012 * Math.Pow(s, 6) + 6435 * Math.Pow(s, 8));
            LegArrEx[8, 1] = (9.0 / 16) * c * (-35 * s + 385 * Math.Pow(s, 3) - 1001 * Math.Pow(s, 5) + 715 * Math.Pow(s, 7));
            LegArrEx[8, 2] = (315.0 / 16) * c * c * (-1 + 33 * Math.Pow(s, 2) - 143 * Math.Pow(s, 4) + 143 * Math.Pow(s, 6));
            LegArrEx[8, 3] = (3465.0 / 8) * Math.Pow(c, 3) * (3 * s - 26 * Math.Pow(s, 3) + 39 * Math.Pow(s, 5));
            LegArrEx[8, 4] = 10395.0 / 8 * Math.Pow(c * c, 2) * (1 - 26 * Math.Pow(s, 2) + 65 * Math.Pow(s, 4));
            LegArrEx[8, 5] = (135135.0 / 2) * Math.Pow(c, 5) * (-s + 5 * Math.Pow(s, 3));
            LegArrEx[8, 6] = (135135.0 / 2) * Math.Pow(c * c, 3) * (-1 + 15 * Math.Pow(s, 2));
            LegArrEx[8, 7] = 2027025.0 * s * Math.Pow(c, 7);
            LegArrEx[8, 8] = 2027025.0 * Math.Pow(c * c, 4);

            LegArrEx[9, 0] = 1.0 / 128 * (315 * s - 4620 * Math.Pow(s, 3) + 18018 * Math.Pow(s, 5) - 25740 * Math.Pow(s, 7) + 12155 * Math.Pow(s, 9));
            LegArrEx[9, 1] = (45.0 / 128) * c * (7 - 308 * Math.Pow(s, 2) + 2002 * Math.Pow(s, 4) - 4004 * Math.Pow(s, 6) + 2431 * Math.Pow(s, 8));
            LegArrEx[9, 2] = (495.0 / 16) * c * c * (-7 * s + 91 * Math.Pow(s, 3) - 273 * Math.Pow(s, 5) + 221 * Math.Pow(s, 7));
            LegArrEx[9, 3] = (3465.0 / 16) * Math.Pow(c, 3) * (-1 + 39 * Math.Pow(s, 2) - 195 * Math.Pow(s, 4) + 221 * Math.Pow(s, 6));
            LegArrEx[9, 4] = 135135.0 / 8 * Math.Pow(c * c, 2) * (s - 10 * Math.Pow(s, 3) + 17 * Math.Pow(s, 5));
            LegArrEx[9, 5] = (135135.0 / 8) * Math.Pow(c, 5) * (1 - 30 * Math.Pow(s, 2) + 85 * Math.Pow(s, 4));
            LegArrEx[9, 6] = (675675.0 / 2) * Math.Pow(c * c, 3) * (-3 * s + 17 * Math.Pow(s, 3));
            LegArrEx[9, 7] = (2027025.0 / 2) * Math.Pow(c, 7) * (-1 + 17 * Math.Pow(s, 2));
            LegArrEx[9, 8] = 34459425.0 * s * Math.Pow(c * c, 4);
            LegArrEx[9, 9] = 34459425.0 * Math.Pow(c, 9);

            LegArrEx[10, 0] = 1.0 / 256 * (-63 + 3465 * Math.Pow(s, 2) - 30030 * Math.Pow(s, 4) + 90090 * Math.Pow(s, 6) - 109395 * Math.Pow(s, 8) + 46189 * Math.Pow(s, 10));
            LegArrEx[10, 1] = (55.0 / 128) * c * (63 * s - 1092 * Math.Pow(s, 3) + 4914 * Math.Pow(s, 5) - 7956 * Math.Pow(s, 7) + 4199 * Math.Pow(s, 9));
            LegArrEx[10, 2] = (495.0 / 128) * c * c * (7 - 364 * Math.Pow(s, 2) + 2730 * Math.Pow(s, 4) - 6188 * Math.Pow(s, 6) + 4199 * Math.Pow(s, 8));
            LegArrEx[10, 3] = (6435.0 / 16) * Math.Pow(c, 3) * (-7 * s + 105 * Math.Pow(s, 3) - 357 * Math.Pow(s, 5) + 323 * Math.Pow(s, 7));
            LegArrEx[10, 4] = 45045.0 / 16 * Math.Pow(c * c, 2) * (-1 + 45 * Math.Pow(s, 2) - 255 * Math.Pow(s, 4) + 323 * Math.Pow(s, 6));
            LegArrEx[10, 5] = (135135.0 / 8) * Math.Pow(c, 5) * (15 * s - 170 * Math.Pow(s, 3) + 323 * Math.Pow(s, 5));
            LegArrEx[10, 6] = (675675.0 / 8) * Math.Pow(c * c, 3) * (3 - 102 * Math.Pow(s, 2) + 323 * Math.Pow(s, 4));
            LegArrEx[10, 7] = (11486475.0 / 2) * Math.Pow(c, 7) * (-3 * s + 19 * Math.Pow(s, 3));
            LegArrEx[10, 8] = 34459425.0 / 2 * Math.Pow(c * c, 4) * (-1 + 19 * Math.Pow(s, 2));
            LegArrEx[10, 9] = 654729075.0 * s * Math.Pow(c, 9);
            LegArrEx[10, 10] = 654729075.0 * Math.Pow(c * c, 5);

            LegArrEx[11, 0] = 1.0 / 256 * (-693 * s + 15015 * Math.Pow(s, 3) - 90090 * Math.Pow(s, 5) + 218790 * Math.Pow(s, 7) - 230945 * Math.Pow(s, 9) + 88179 * Math.Pow(s, 11));
            LegArrEx[11, 1] = (33.0 / 256) * c * (-21 + 1365 * Math.Pow(s, 2) - 13650 * Math.Pow(s, 4) + 46410 * Math.Pow(s, 6) - 62985 * Math.Pow(s, 8) + 29393 * Math.Pow(s, 10));
            LegArrEx[11, 2] = (2145.0 / 128) * c * c * (21 * s - 420 * Math.Pow(s, 3) + 2142 * Math.Pow(s, 5) - 3876 * Math.Pow(s, 7) + 2261 * Math.Pow(s, 9));
            LegArrEx[11, 3] = (45045.0 / 128) * Math.Pow(c, 3) * (1 - 60 * Math.Pow(s, 2) + 510 * Math.Pow(s, 4) - 1292 * Math.Pow(s, 6) + 969 * Math.Pow(s, 8));
            LegArrEx[11, 4] = 135135.0 / 16 * Math.Pow(c * c, 2) * (-5 * s + 85 * Math.Pow(s, 3) - 323 * Math.Pow(s, 5) + 323 * Math.Pow(s, 7));
            LegArrEx[11, 5] = (135135.0 / 16) * Math.Pow(c, 5) * (-5 + 255 * Math.Pow(s, 2) - 1615 * Math.Pow(s, 4) + 2261 * Math.Pow(s, 6));
            LegArrEx[11, 6] = (2297295.0 / 8) * Math.Pow(c * c, 3) * (15 * s - 190 * Math.Pow(s, 3) + 399 * Math.Pow(s, 5));
            LegArrEx[11, 7] = (34459425.0 / 8) * Math.Pow(c, 7) * (1 - 38 * Math.Pow(s, 2) + 133 * Math.Pow(s, 4));
            LegArrEx[11, 8] = 654729075.0 / 2 * Math.Pow(c * c, 4) * (-s + 7 * Math.Pow(s, 3));
            LegArrEx[11, 9] = (654729075.0 / 2) * Math.Pow(c, 9) * (-1 + 21 * Math.Pow(s, 2));
            LegArrEx[11, 10] = 13749310575.0 * s * Math.Pow(c * c, 5);
            LegArrEx[11, 11] = 13749310575.0 * Math.Pow(c, 11);

            LegArrEx[12, 0] = (231.0 - 18018 * Math.Pow(s, 2) + 225225 * Math.Pow(s, 4) - 1021020 * Math.Pow(s, 6) + 2078505 * Math.Pow(s, 8) - 1939938 * Math.Pow(s, 10) + 676039 * Math.Pow(s, 12)) / 1024;
            LegArrEx[12, 1] = (39.0 / 256) * c * (-231 * s + 5775 * Math.Pow(s, 3) - 39270 * Math.Pow(s, 5) + 106590 * Math.Pow(s, 7) - 124355 * Math.Pow(s, 9) + 52003 * Math.Pow(s, 11));
            LegArrEx[12, 2] = (3003.0 / 256) * c * c * (-3 + 225 * Math.Pow(s, 2) - 2550 * Math.Pow(s, 4) + 9690 * Math.Pow(s, 6) - 14535 * Math.Pow(s, 8) + 7429 * Math.Pow(s, 10));
            LegArrEx[12, 3] = (15015.0 / 128) * Math.Pow(c, 3) * (45 * s - 1020 * Math.Pow(s, 3) + 5814 * Math.Pow(s, 5) - 11628 * Math.Pow(s, 7) + 7429 * Math.Pow(s, 9));
            LegArrEx[12, 4] = 135135.0 / 128 * Math.Pow(c * c, 2) * (5 - 340 * Math.Pow(s, 2) + 3230 * Math.Pow(s, 4) - 9044 * Math.Pow(s, 6) + 7429 * Math.Pow(s, 8));
            LegArrEx[12, 5] = (2297295.0 / 16) * Math.Pow(c, 5) * (-5 * s + 95 * Math.Pow(s, 3) - 399 * Math.Pow(s, 5) + 437 * Math.Pow(s, 7));
            LegArrEx[12, 6] = (2297295.0 / 16) * Math.Pow(c * c, 3) * (-5 + 285 * Math.Pow(s, 2) - 1995 * Math.Pow(s, 4) + 3059 * Math.Pow(s, 6));
            LegArrEx[12, 7] = (130945815.0 / 8) * Math.Pow(c, 7) * (5 * s - 70 * Math.Pow(s, 3) + 161 * Math.Pow(s, 5));
            LegArrEx[12, 8] = 654729075.0 / 8 * Math.Pow(c * c, 4) * (1 - 42 * Math.Pow(s, 2) + 161 * Math.Pow(s, 4));
            LegArrEx[12, 9] = (4583103525.0 / 2) * Math.Pow(c, 9) * (-3 * s + 23 * Math.Pow(s, 3));
            LegArrEx[12, 10] = (13749310575.0 / 2) * Math.Pow(c * c, 5) * (-1 + 23 * Math.Pow(s, 2));
            LegArrEx[12, 11] = 316234143225.0 * s * Math.Pow(c, 11);
            LegArrEx[12, 12] = 316234143225.0 * Math.Pow(c * c, 6);

            LegArrEx[13, 0] = (3003.0 * s - 90090 * Math.Pow(s, 3) + 765765 * Math.Pow(s, 5) - 2771340 * Math.Pow(s, 7) + 4849845 * Math.Pow(s, 9) - 4056234 * Math.Pow(s, 11) + 1300075 * Math.Pow(s, 13)) / 1024;
            LegArrEx[13, 1] = ((91.0 * c * (33 - 2970 * Math.Pow(s, 2) + 42075 * Math.Pow(s, 4) - 213180 * Math.Pow(s, 6) + 479655 * Math.Pow(s, 8) - 490314 * Math.Pow(s, 10) + 185725 * Math.Pow(s, 12)) / 1024));
            LegArrEx[13, 2] = (1365.0 / 256) * c * c * (-99 * s + 2805 * Math.Pow(s, 3) - 21318 * Math.Pow(s, 5) + 63954 * Math.Pow(s, 7) - 81719 * Math.Pow(s, 9) + 37145 * Math.Pow(s, 11));
            LegArrEx[13, 3] = (15015.0 / 256) * Math.Pow(c, 3) * (-9 + 765 * Math.Pow(s, 2) - 9690 * Math.Pow(s, 4) + 40698 * Math.Pow(s, 6) - 66861 * Math.Pow(s, 8) + 37145 * Math.Pow(s, 10));
            LegArrEx[13, 4] = 255255.0 / 128 * Math.Pow(c * c, 2) * (45 * s - 1140 * Math.Pow(s, 3) + 7182 * Math.Pow(s, 5) - 15732 * Math.Pow(s, 7) + 10925 * Math.Pow(s, 9));
            LegArrEx[13, 5] = (2297295.0 / 128) * Math.Pow(c, 5) * (5 - 380 * Math.Pow(s, 2) + 3990 * Math.Pow(s, 4) - 12236 * Math.Pow(s, 6) + 10925 * Math.Pow(s, 8));
            LegArrEx[13, 6] = (43648605.0 / 16) * Math.Pow(c * c, 3) * (-5 * s + 105 * Math.Pow(s, 3) - 483 * Math.Pow(s, 5) + 575 * Math.Pow(s, 7));
            LegArrEx[13, 7] = (218243025.0 / 16) * Math.Pow(c, 7) * (-1 + 63 * Math.Pow(s, 2) - 483 * Math.Pow(s, 4) + 805 * Math.Pow(s, 6));
            LegArrEx[13, 8] = 4583103525.0 / 8 * Math.Pow(c * c, 4) * (3 * s - 46 * Math.Pow(s, 3) + 115 * Math.Pow(s, 5));
            LegArrEx[13, 9] = (4583103525.0 / 8) * Math.Pow(c, 9) * (3 - 138 * Math.Pow(s, 2) + 575 * Math.Pow(s, 4));
            LegArrEx[13, 10] = (105411381075.0 / 2) * Math.Pow(c * c, 5) * (-3 * s + 25 * Math.Pow(s, 3));
            LegArrEx[13, 11] = (316234143225.0 / 2) * Math.Pow(c, 11) * (-1 + 25 * Math.Pow(s, 2));
            LegArrEx[13, 12] = 7905853580625.0 * s * Math.Pow(c * c, 6);
            LegArrEx[13, 13] = 7905853580625.0 * Math.Pow(c, 13);

            LegArrEx[14, 0] = (-429.0 + 45045 * Math.Pow(s, 2) - 765765 * Math.Pow(s, 4) + 4849845 * Math.Pow(s, 6) - 14549535 * Math.Pow(s, 8) + 22309287 * Math.Pow(s, 10) - 16900975 * Math.Pow(s, 12) + 5014575 * Math.Pow(s, 14)) / 2048;
            LegArrEx[14, 1] = ((105.0 * c * (429 * s - 14586 * Math.Pow(s, 3) + 138567 * Math.Pow(s, 5) - 554268 * Math.Pow(s, 7) + 1062347 * Math.Pow(s, 9) - 965770 * Math.Pow(s, 11) + 334305 * Math.Pow(s, 13)) / 1024));
            LegArrEx[14, 2] = ((1365.0 * c * c * (33 - 3366 * Math.Pow(s, 2) + 53295 * Math.Pow(s, 4) - 298452 * Math.Pow(s, 6) + 735471 * Math.Pow(s, 8) - 817190 * Math.Pow(s, 10) + 334305 * Math.Pow(s, 12)) / 1024));
            LegArrEx[14, 3] = (23205.0 / 256) * Math.Pow(c, 3) * (-99 * s + 3135 * Math.Pow(s, 3) - 26334 * Math.Pow(s, 5) + 86526 * Math.Pow(s, 7) - 120175 * Math.Pow(s, 9) + 58995 * Math.Pow(s, 11));
            LegArrEx[14, 4] = 2297295.0 / 256 * Math.Pow(c * c, 2) * (-1 + 95 * Math.Pow(s, 2) - 1330 * Math.Pow(s, 4) + 6118 * Math.Pow(s, 6) - 10925 * Math.Pow(s, 8) + 6555 * Math.Pow(s, 10));
            LegArrEx[14, 5] = (43648605.0 / 128) * Math.Pow(c, 5) * (5 * s - 140 * Math.Pow(s, 3) + 966 * Math.Pow(s, 5) - 2300 * Math.Pow(s, 7) + 1725 * Math.Pow(s, 9));
            LegArrEx[14, 6] = (218243025.0 / 128) * Math.Pow(c * c, 3) * (1 - 84 * Math.Pow(s, 2) + 966 * Math.Pow(s, 4) - 3220 * Math.Pow(s, 6) + 3105 * Math.Pow(s, 8));
            LegArrEx[14, 7] = (654729075.0 / 16) * Math.Pow(c, 7) * (-7 * s + 161 * Math.Pow(s, 3) - 805 * Math.Pow(s, 5) + 1035 * Math.Pow(s, 7));
            LegArrEx[14, 8] = 4583103525.0 / 16 * Math.Pow(c * c, 4) * (-1 + 69 * Math.Pow(s, 2) - 575 * Math.Pow(s, 4) + 1035 * Math.Pow(s, 6));
            LegArrEx[14, 9] = (105411381075.0 / 8) * Math.Pow(c, 9) * (3 * s - 50 * Math.Pow(s, 3) + 135 * Math.Pow(s, 5));
            LegArrEx[14, 10] = (316234143225.0 / 8) * Math.Pow(c * c, 5) * (1 - 50 * Math.Pow(s, 2) + 225 * Math.Pow(s, 4));
            LegArrEx[14, 11] = (7905853580625.0 / 2) * Math.Pow(c, 11) * (-s + 9 * Math.Pow(s, 3));
            LegArrEx[14, 12] = 7905853580625.0 / 2 * Math.Pow(c * c, 6) * (-1 + 27 * Math.Pow(s, 2));
            LegArrEx[14, 13] = 213458046676875.0 * s * Math.Pow(c, 13);
            LegArrEx[14, 14] = 213458046676875.0 * Math.Pow(c * c, 7);

            LegArrEx[15, 0] = (-6435.0 * s + 255255 * Math.Pow(s, 3) - 2909907 * Math.Pow(s, 5) + 14549535 * Math.Pow(s, 7) - 37182145 * Math.Pow(s, 9) + 50702925 * Math.Pow(s, 11) - 35102025 * Math.Pow(s, 13) + 9694845 * Math.Pow(s, 15)) / 2048;
            LegArrEx[15, 1] = ((15.0 * c * (-429 + 51051 * Math.Pow(s, 2) - 969969 * Math.Pow(s, 4) + 6789783 * Math.Pow(s, 6) - 22309287 * Math.Pow(s, 8) + 37182145 * Math.Pow(s, 10) - 30421755 * Math.Pow(s, 12) + 9694845 * Math.Pow(s, 14)) / 2048));
            LegArrEx[15, 2] = ((1785.0 * c * c * (429 * s - 16302 * Math.Pow(s, 3) + 171171 * Math.Pow(s, 5) - 749892 * Math.Pow(s, 7) + 1562275 * Math.Pow(s, 9) - 1533870 * Math.Pow(s, 11) + 570285 * Math.Pow(s, 13)) / 1024));
            LegArrEx[15, 3] = ((69615.0 * Math.Pow(c, 3) * (11 - 1254 * Math.Pow(s, 2) + 21945 * Math.Pow(s, 4) - 134596 * Math.Pow(s, 6) + 360525 * Math.Pow(s, 8) - 432630 * Math.Pow(s, 10) + 190095 * Math.Pow(s, 12)) / 1024));
            LegArrEx[15, 4] = 3968055.0 / 256 * Math.Pow(c * c, 2) * (-11 * s + 385 * Math.Pow(s, 3) - 3542 * Math.Pow(s, 5) + 12650 * Math.Pow(s, 7) - 18975 * Math.Pow(s, 9) + 10005 * Math.Pow(s, 11));
            LegArrEx[15, 5] = (43648605.0 / 256) * Math.Pow(c, 5) * (-1 + 105 * Math.Pow(s, 2) - 1610 * Math.Pow(s, 4) + 8050 * Math.Pow(s, 6) - 15525 * Math.Pow(s, 8) + 10005 * Math.Pow(s, 10));
            LegArrEx[15, 6] = (218243025.0 / 128) * Math.Pow(c * c, 3) * (21 * s - 644 * Math.Pow(s, 3) + 4830 * Math.Pow(s, 5) - 12420 * Math.Pow(s, 7) + 10005 * Math.Pow(s, 9));
            LegArrEx[15, 7] = (654729075.0 / 128) * Math.Pow(c, 7) * (7 - 644 * Math.Pow(s, 2) + 8050 * Math.Pow(s, 4) - 28980 * Math.Pow(s, 6) + 30015 * Math.Pow(s, 8));
            LegArrEx[15, 8] = 15058768725.0 / 16 * Math.Pow(c * c, 4) * (-7 * s + 175 * Math.Pow(s, 3) - 945 * Math.Pow(s, 5) + 1305 * Math.Pow(s, 7));
            LegArrEx[15, 9] = (105411381075.0 / 16) * Math.Pow(c, 9) * (-1 + 75 * Math.Pow(s, 2) - 675 * Math.Pow(s, 4) + 1305 * Math.Pow(s, 6));
            LegArrEx[15, 10] = (1581170716125.0 / 8) * Math.Pow(c * c, 5) * (5 * s - 90 * Math.Pow(s, 3) + 261 * Math.Pow(s, 5));
            LegArrEx[15, 11] = (7905853580625.0 / 8) * Math.Pow(c, 11) * (1 - 54 * Math.Pow(s, 2) + 261 * Math.Pow(s, 4));
            LegArrEx[15, 12] = 71152682225625.0 / 2 * Math.Pow(c * c, 6) * (-3 * s + 29 * Math.Pow(s, 3));
            LegArrEx[15, 13] = (213458046676875.0 / 2) * Math.Pow(c, 13) * (-1 + 29 * Math.Pow(s, 2));
            LegArrEx[15, 14] = 6190283353629375.0 * s * Math.Pow(c * c, 7);
            LegArrEx[15, 15] = 6190283353629375.0 * Math.Pow(c, 15);

            LegArrEx[16, 0] = (6435.0 - 875160 * Math.Pow(s, 2) + 19399380 * Math.Pow(s, 4) - 162954792 * Math.Pow(s, 6) + 669278610 * Math.Pow(s, 8) - 1487285800 * Math.Pow(s, 10) + 1825305300 * Math.Pow(s, 12) - 1163381400 * Math.Pow(s, 14) + 300540195 * Math.Pow(s, 16)) / 32768;
            LegArrEx[16, 1] = ((17.0 * c * (-6435 * s + 285285 * Math.Pow(s, 3) - 3594591 * Math.Pow(s, 5) + 19684665 * Math.Pow(s, 7) - 54679625 * Math.Pow(s, 9) + 80528175 * Math.Pow(s, 11) - 59879925 * Math.Pow(s, 13) + 17678835 * Math.Pow(s, 15)) / 2048));
            LegArrEx[16, 2] = ((765.0 * c * c * (-143 + 19019 * Math.Pow(s, 2) - 399399 * Math.Pow(s, 4) + 3062059 * Math.Pow(s, 6) - 10935925 * Math.Pow(s, 8) + 19684665 * Math.Pow(s, 10) - 17298645 * Math.Pow(s, 12) + 5892945 * Math.Pow(s, 14)) / 2048));
            LegArrEx[16, 3] = ((101745.0 * Math.Pow(c, 3) * (143 * s - 6006 * Math.Pow(s, 3) + 69069 * Math.Pow(s, 5) - 328900 * Math.Pow(s, 7) + 740025 * Math.Pow(s, 9) - 780390 * Math.Pow(s, 11) + 310155 * Math.Pow(s, 13)) / 1024));
            LegArrEx[16, 4] = (1322685.0 * Math.Pow(c * c, 2) * (11 - 1386 * Math.Pow(s, 2) + 26565 * Math.Pow(s, 4) - 177100 * Math.Pow(s, 6) + 512325 * Math.Pow(s, 8) - 660330 * Math.Pow(s, 10) + 310155 * Math.Pow(s, 12)) / 1024);
            LegArrEx[16, 5] = (3968055.0 / 256) * Math.Pow(c, 5) * (-231 * s + 8855 * Math.Pow(s, 3) - 88550 * Math.Pow(s, 5) + 341550 * Math.Pow(s, 7) - 550275 * Math.Pow(s, 9) + 310155 * Math.Pow(s, 11));
            LegArrEx[16, 6] = (43648605.0 / 256) * Math.Pow(c * c, 3) * (-21 + 2415 * Math.Pow(s, 2) - 40250 * Math.Pow(s, 4) + 217350 * Math.Pow(s, 6) - 450225 * Math.Pow(s, 8) + 310155 * Math.Pow(s, 10));
            LegArrEx[16, 7] = (5019589575.0 / 128) * Math.Pow(c, 7) * (21 * s - 700 * Math.Pow(s, 3) + 5670 * Math.Pow(s, 5) - 15660 * Math.Pow(s, 7) + 13485 * Math.Pow(s, 9));
            LegArrEx[16, 8] = 15058768725.0 / 128 * Math.Pow(c * c, 4) * (7 - 700 * Math.Pow(s, 2) + 9450 * Math.Pow(s, 4) - 36540 * Math.Pow(s, 6) + 40455 * Math.Pow(s, 8));
            LegArrEx[16, 9] = (75293843625.0 / 16) * Math.Pow(c, 9) * (-35 * s + 945 * Math.Pow(s, 3) - 5481 * Math.Pow(s, 5) + 8091 * Math.Pow(s, 7));
            LegArrEx[16, 10] = (527056905375.0 / 16) * Math.Pow(c * c, 5) * (-5 + 405 * Math.Pow(s, 2) - 3915 * Math.Pow(s, 4) + 8091 * Math.Pow(s, 6));
            LegArrEx[16, 11] = (14230536445125.0 / 8) * Math.Pow(c, 11) * (15 * s - 290 * Math.Pow(s, 3) + 899 * Math.Pow(s, 5));
            LegArrEx[16, 12] = 71152682225625.0 / 8 * Math.Pow(c * c, 6) * (3 - 174 * Math.Pow(s, 2) + 899 * Math.Pow(s, 4));
            LegArrEx[16, 13] = (2063427784543125.0 / 2) * Math.Pow(c, 13) * (-3 * s + 31 * Math.Pow(s, 3));
            LegArrEx[16, 14] = (6190283353629375.0 / 2) * Math.Pow(c * c, 7) * (-1 + 31 * Math.Pow(s, 2));
            LegArrEx[16, 15] = 191898783962510625.0 * s * Math.Pow(c, 15);
            LegArrEx[16, 16] = 191898783962510625.0 * Math.Pow(c * c, 8);

            LegArrEx[17, 0] = (109395.0 * s - 5542680 * Math.Pow(s, 3) + 81477396 * Math.Pow(s, 5) - 535422888 * Math.Pow(s, 7) + 1859107250 * Math.Pow(s, 9) - 3650610600 * Math.Pow(s, 11) + 4071834900 * Math.Pow(s, 13) - 2404321560 * Math.Pow(s, 15) + 583401555 * Math.Pow(s, 17)) / 32768;
            LegArrEx[17, 1] = ((153.0 * c * (715 - 108680 * Math.Pow(s, 2) + 2662660 * Math.Pow(s, 4) - 24496472 * Math.Pow(s, 6) + 109359250 * Math.Pow(s, 8) - 262462200 * Math.Pow(s, 10) + 345972900 * Math.Pow(s, 12) - 235717800 * Math.Pow(s, 14) + 64822395 * Math.Pow(s, 16)) / 32768));
            LegArrEx[17, 2] = ((2907.0 * c * c * (-715 * s + 35035 * Math.Pow(s, 3) - 483483 * Math.Pow(s, 5) + 2877875 * Math.Pow(s, 7) - 8633625 * Math.Pow(s, 9) + 13656825 * Math.Pow(s, 11) - 10855425 * Math.Pow(s, 13) + 3411705 * Math.Pow(s, 15)) / 2048));
            LegArrEx[17, 3] = ((14535.0 * Math.Pow(c, 3) * (-143 + 21021 * Math.Pow(s, 2) - 483483 * Math.Pow(s, 4) + 4029025 * Math.Pow(s, 6) - 15540525 * Math.Pow(s, 8) + 30045015 * Math.Pow(s, 10) - 28224105 * Math.Pow(s, 12) + 10235115 * Math.Pow(s, 14)) / 2048));
            LegArrEx[17, 4] = (305235.0 * Math.Pow(c * c, 2) * (1001 * s - 46046 * Math.Pow(s, 3) + 575575 * Math.Pow(s, 5) - 2960100 * Math.Pow(s, 7) + 7153575 * Math.Pow(s, 9) - 8064030 * Math.Pow(s, 11) + 3411705 * Math.Pow(s, 13)) / 1024);
            LegArrEx[17, 5] = ((43648605.0 * Math.Pow(c, 5) * (7 - 966 * Math.Pow(s, 2) + 20125 * Math.Pow(s, 4) - 144900 * Math.Pow(s, 6) + 450225 * Math.Pow(s, 8) - 620310 * Math.Pow(s, 10) + 310155 * Math.Pow(s, 12)) / 1024));
            LegArrEx[17, 6] = (1003917915.0 / 256) * Math.Pow(c * c, 3) * (-21 * s + 875 * Math.Pow(s, 3) - 9450 * Math.Pow(s, 5) + 39150 * Math.Pow(s, 7) - 67425 * Math.Pow(s, 9) + 40455 * Math.Pow(s, 11));
            LegArrEx[17, 7] = (3011753745.0 / 256) * Math.Pow(c, 7) * (-7 + 875 * Math.Pow(s, 2) - 15750 * Math.Pow(s, 4) + 91350 * Math.Pow(s, 6) - 202275 * Math.Pow(s, 8) + 148335 * Math.Pow(s, 10));
            LegArrEx[17, 8] = 75293843625.0 / 128 * Math.Pow(c * c, 4) * (35 * s - 1260 * Math.Pow(s, 3) + 10962 * Math.Pow(s, 5) - 32364 * Math.Pow(s, 7) + 29667 * Math.Pow(s, 9));
            LegArrEx[17, 9] = (75293843625.0 / 128) * Math.Pow(c, 9) * (35 - 3780 * Math.Pow(s, 2) + 54810 * Math.Pow(s, 4) - 226548 * Math.Pow(s, 6) + 267003 * Math.Pow(s, 8));
            LegArrEx[17, 10] = (2032933777875.0 / 16) * Math.Pow(c * c, 5) * (-35 * s + 1015 * Math.Pow(s, 3) - 6293 * Math.Pow(s, 5) + 9889 * Math.Pow(s, 7));
            LegArrEx[17, 11] = (14230536445125.0 / 16) * Math.Pow(c, 11) * (-5 + 435 * Math.Pow(s, 2) - 4495 * Math.Pow(s, 4) + 9889 * Math.Pow(s, 6));
            LegArrEx[17, 12] = 412685556908625.0 / 8 * Math.Pow(c * c, 6) * (15 * s - 310 * Math.Pow(s, 3) + 1023 * Math.Pow(s, 5));
            LegArrEx[17, 13] = (6190283353629375.0 / 8) * Math.Pow(c, 13) * (1 - 62 * Math.Pow(s, 2) + 341 * Math.Pow(s, 4));
            LegArrEx[17, 14] = (191898783962510625.0 / 2) * Math.Pow(c * c, 7) * (-s + 11 * Math.Pow(s, 3));
            LegArrEx[17, 15] = (191898783962510625.0 / 2) * Math.Pow(c, 15) * (-1 + 33 * Math.Pow(s, 2));
            LegArrEx[17, 16] = 6332659870762850625.0 * s * Math.Pow(c * c, 8);
            LegArrEx[17, 17] = 6332659870762850625.0 * Math.Pow(c, 17);

            LegArrEx[18, 0] = (-12155.0 + 2078505 * Math.Pow(s, 2) - 58198140 * Math.Pow(s, 4) + 624660036 * Math.Pow(s, 6) - 3346393050 * Math.Pow(s, 8) + 10039179150 * Math.Pow(s, 10) - 17644617900 * Math.Pow(s, 12) + 18032411700 * Math.Pow(s, 14) - 9917826435 * Math.Pow(s, 16) + 2268783825 * Math.Pow(s, 18)) / 65536;
            LegArrEx[18, 1] = ((171.0 * c * (12155 * s - 680680 * Math.Pow(s, 3) + 10958948 * Math.Pow(s, 5) - 78278200 * Math.Pow(s, 7) + 293543250 * Math.Pow(s, 9) - 619109400 * Math.Pow(s, 11) + 738168900 * Math.Pow(s, 13) - 463991880 * Math.Pow(s, 15) + 119409675 * Math.Pow(s, 17)) / 32768));
            LegArrEx[18, 2] = ((14535.0 * c * c * (143 - 24024 * Math.Pow(s, 2) + 644644 * Math.Pow(s, 4) - 6446440 * Math.Pow(s, 6) + 31081050 * Math.Pow(s, 8) - 80120040 * Math.Pow(s, 10) + 112896420 * Math.Pow(s, 12) - 81880920 * Math.Pow(s, 14) + 23881935 * Math.Pow(s, 16)) / 32768));
            LegArrEx[18, 3] = ((101745.0 * Math.Pow(c, 3) * (-429 * s + 23023 * Math.Pow(s, 3) - 345345 * Math.Pow(s, 5) + 2220075 * Math.Pow(s, 7) - 7153575 * Math.Pow(s, 9) + 12096045 * Math.Pow(s, 11) - 10235115 * Math.Pow(s, 13) + 3411705 * Math.Pow(s, 15)) / 2048));
            LegArrEx[18, 4] = (3357585.0 * Math.Pow(c * c, 2) * (-13 + 2093 * Math.Pow(s, 2) - 52325 * Math.Pow(s, 4) + 470925 * Math.Pow(s, 6) - 1950975 * Math.Pow(s, 8) + 4032015 * Math.Pow(s, 10) - 4032015 * Math.Pow(s, 12) + 1550775 * Math.Pow(s, 14)) / 2048);
            LegArrEx[18, 5] = ((77224455.0 * Math.Pow(c, 5) * (91 * s - 4550 * Math.Pow(s, 3) + 61425 * Math.Pow(s, 5) - 339300 * Math.Pow(s, 7) + 876525 * Math.Pow(s, 9) - 1051830 * Math.Pow(s, 11) + 471975 * Math.Pow(s, 13)) / 1024));
            LegArrEx[18, 6] = ((1003917915.0 * Math.Pow(c * c, 3) * (7 - 1050 * Math.Pow(s, 2) + 23625 * Math.Pow(s, 4) - 182700 * Math.Pow(s, 6) + 606825 * Math.Pow(s, 8) - 890010 * Math.Pow(s, 10) + 471975 * Math.Pow(s, 12)) / 1024));
            LegArrEx[18, 7] = (75293843625.0 / 256) * Math.Pow(c, 7) * (-7 * s + 315 * Math.Pow(s, 3) - 3654 * Math.Pow(s, 5) + 16182 * Math.Pow(s, 7) - 29667 * Math.Pow(s, 9) + 18879 * Math.Pow(s, 11));
            LegArrEx[18, 8] = 75293843625.0 / 256 * Math.Pow(c * c, 4) * (-7 + 945 * Math.Pow(s, 2) - 18270 * Math.Pow(s, 4) + 113274 * Math.Pow(s, 6) - 267003 * Math.Pow(s, 8) + 207669 * Math.Pow(s, 10));
            LegArrEx[18, 9] = (225881530875.0 / 128) * Math.Pow(c, 9) * (315 * s - 12180 * Math.Pow(s, 3) + 113274 * Math.Pow(s, 5) - 356004 * Math.Pow(s, 7) + 346115 * Math.Pow(s, 9));
            LegArrEx[18, 10] = (14230536445125.0 / 128) * Math.Pow(c * c, 5) * (5 - 580 * Math.Pow(s, 2) + 8990 * Math.Pow(s, 4) - 39556 * Math.Pow(s, 6) + 49445 * Math.Pow(s, 8));
            LegArrEx[18, 11] = (412685556908625.0 / 16) * Math.Pow(c, 11) * (-5 * s + 155 * Math.Pow(s, 3) - 1023 * Math.Pow(s, 5) + 1705 * Math.Pow(s, 7));
            LegArrEx[18, 12] = 2063427784543125.0 / 16 * Math.Pow(c * c, 6) * (-1 + 93 * Math.Pow(s, 2) - 1023 * Math.Pow(s, 4) + 2387 * Math.Pow(s, 6));
            LegArrEx[18, 13] = (191898783962510625.0 / 8) * Math.Pow(c, 13) * (s - 22 * Math.Pow(s, 3) + 77 * Math.Pow(s, 5));
            LegArrEx[18, 14] = (191898783962510625.0 / 8) * Math.Pow(c * c, 7) * (1 - 66 * Math.Pow(s, 2) + 385 * Math.Pow(s, 4));
            LegArrEx[18, 15] = (2110886623587616875.0 / 2) * Math.Pow(c, 15) * (-3 * s + 35 * Math.Pow(s, 3));
            LegArrEx[18, 16] = 6332659870762850625.0 / 2 * Math.Pow(c * c, 8) * (-1 + 35 * Math.Pow(s, 2));
            LegArrEx[18, 17] = 221643095476699771875.0 * s * Math.Pow(c, 17);
            LegArrEx[18, 18] = 221643095476699771875.0 * Math.Pow(c * c, 9);

            LegArrEx[19, 0] = (-230945.0 * s + 14549535 * Math.Pow(s, 3) - 267711444 * Math.Pow(s, 5) + 2230928700 * Math.Pow(s, 7) - 10039179150 * Math.Pow(s, 9) + 26466926850 * Math.Pow(s, 11) - 42075627300 * Math.Pow(s, 13) + 39671305740 * Math.Pow(s, 15) - 20419054425 * Math.Pow(s, 17) + 4418157975 * Math.Pow(s, 19)) / 65536;
            LegArrEx[19, 1] = ((95.0 * c * (-2431 + 459459 * Math.Pow(s, 2) - 14090076 * Math.Pow(s, 4) + 164384220 * Math.Pow(s, 6) - 951080130 * Math.Pow(s, 8) + 3064591530 * Math.Pow(s, 10) - 5757717420 * Math.Pow(s, 12) + 6263890380 * Math.Pow(s, 14) - 3653936055 * Math.Pow(s, 16) + 883631595 * Math.Pow(s, 18)) / 65536));
            LegArrEx[19, 2] = ((5985.0 * c * c * (7293 * s - 447304 * Math.Pow(s, 3) + 7827820 * Math.Pow(s, 5) - 60386040 * Math.Pow(s, 7) + 243221550 * Math.Pow(s, 9) - 548354040 * Math.Pow(s, 11) + 695987820 * Math.Pow(s, 13) - 463991880 * Math.Pow(s, 15) + 126233085 * Math.Pow(s, 17)) / 32768));
            LegArrEx[19, 3] = ((1119195.0 * Math.Pow(c, 3) * (39 - 7176 * Math.Pow(s, 2) + 209300 * Math.Pow(s, 4) - 2260440 * Math.Pow(s, 6) + 11705850 * Math.Pow(s, 8) - 32256120 * Math.Pow(s, 10) + 48384180 * Math.Pow(s, 12) - 37218600 * Math.Pow(s, 14) + 11475735 * Math.Pow(s, 16)) / 32768));
            LegArrEx[19, 4] = (25741485.0 * Math.Pow(c * c, 2) * (-39 * s + 2275 * Math.Pow(s, 3) - 36855 * Math.Pow(s, 5) + 254475 * Math.Pow(s, 7) - 876525 * Math.Pow(s, 9) + 1577745 * Math.Pow(s, 11) - 1415925 * Math.Pow(s, 13) + 498945 * Math.Pow(s, 15)) / 2048);
            LegArrEx[19, 5] = ((77224455.0 * Math.Pow(c, 5) * (-13 + 2275 * Math.Pow(s, 2) - 61425 * Math.Pow(s, 4) + 593775 * Math.Pow(s, 6) - 2629575 * Math.Pow(s, 8) + 5785065 * Math.Pow(s, 10) - 6135675 * Math.Pow(s, 12) + 2494725 * Math.Pow(s, 14)) / 2048));
            LegArrEx[19, 6] = ((1930611375.0 * Math.Pow(c * c, 3) * (91 * s - 4914 * Math.Pow(s, 3) + 71253 * Math.Pow(s, 5) - 420732 * Math.Pow(s, 7) + 1157013 * Math.Pow(s, 9) - 1472562 * Math.Pow(s, 11) + 698523 * Math.Pow(s, 13)) / 1024));
            LegArrEx[19, 7] = ((25097947875.0 * Math.Pow(c, 7) * (7 - 1134 * Math.Pow(s, 2) + 27405 * Math.Pow(s, 4) - 226548 * Math.Pow(s, 6) + 801009 * Math.Pow(s, 8) - 1246014 * Math.Pow(s, 10) + 698523 * Math.Pow(s, 12)) / 1024));
            LegArrEx[19, 8] = 225881530875.0 / 256 * Math.Pow(c * c, 4) * (-63 * s + 3045 * Math.Pow(s, 3) - 37758 * Math.Pow(s, 5) + 178002 * Math.Pow(s, 7) - 346115 * Math.Pow(s, 9) + 232841 * Math.Pow(s, 11));
            LegArrEx[19, 9] = (1581170716125.0 / 256) * Math.Pow(c, 9) * (-9 + 1305 * Math.Pow(s, 2) - 26970 * Math.Pow(s, 4) + 178002 * Math.Pow(s, 6) - 445005 * Math.Pow(s, 8) + 365893 * Math.Pow(s, 10));
            LegArrEx[19, 10] = (45853950767625.0 / 128) * Math.Pow(c * c, 5) * (45 * s - 1860 * Math.Pow(s, 3) + 18414 * Math.Pow(s, 5) - 61380 * Math.Pow(s, 7) + 63085 * Math.Pow(s, 9));
            LegArrEx[19, 11] = (2063427784543125.0 / 128) * Math.Pow(c, 11) * (1 - 124 * Math.Pow(s, 2) + 2046 * Math.Pow(s, 4) - 9548 * Math.Pow(s, 6) + 12617 * Math.Pow(s, 8));
            LegArrEx[19, 12] = 63966261320836875.0 / 16 * Math.Pow(c * c, 6) * (-s + 33 * Math.Pow(s, 3) - 231 * Math.Pow(s, 5) + 407 * Math.Pow(s, 7));
            LegArrEx[19, 13] = (63966261320836875.0 / 16) * Math.Pow(c, 13) * (-1 + 99 * Math.Pow(s, 2) - 1155 * Math.Pow(s, 4) + 2849 * Math.Pow(s, 6));
            LegArrEx[19, 14] = (2110886623587616875.0 / 8) * Math.Pow(c * c, 7) * (3 * s - 70 * Math.Pow(s, 3) + 259 * Math.Pow(s, 5));
            LegArrEx[19, 15] = (2110886623587616875.0 / 8) * Math.Pow(c, 15) * (3 - 210 * Math.Pow(s, 2) + 1295 * Math.Pow(s, 4));
            LegArrEx[19, 16] = 73881031825566590625.0 / 2 * Math.Pow(c * c, 8) * (-3 * s + 37 * Math.Pow(s, 3));
            LegArrEx[19, 17] = (221643095476699771875.0 / 2) * Math.Pow(c, 17) * (-1 + 37 * Math.Pow(s, 2));
            LegArrEx[19, 18] = 8200794532637891559375.0 * s * Math.Pow(c * c, 9);
            LegArrEx[19, 19] = 8200794532637891559375.0 * Math.Pow(c, 19);

            LegArrEx[20, 0] = (1.0 / 262144) * (46189 - 9699690 * Math.Pow(s, 2) + 334639305 * Math.Pow(s, 4) - 4461857400 * Math.Pow(s, 6) + 30117537450 * Math.Pow(s, 8) - 116454478140 * Math.Pow(s, 10) + 273491577450 * Math.Pow(s, 12) - 396713057400 * Math.Pow(s, 14) + 347123925225 * Math.Pow(s, 16) - 167890003050 * Math.Pow(s, 18) + 34461632205 * Math.Pow(s, 20));
            LegArrEx[20, 1] = ((105.0 * c * (-46189 * s + 3187041 * Math.Pow(s, 3) - 63740820 * Math.Pow(s, 5) + 573667380 * Math.Pow(s, 7) - 2772725670 * Math.Pow(s, 9) + 7814045070 * Math.Pow(s, 11) - 13223768580 * Math.Pow(s, 13) + 13223768580 * Math.Pow(s, 15) - 7195285845 * Math.Pow(s, 17) + 1641030105 * Math.Pow(s, 19)) / 65536));
            LegArrEx[20, 2] = ((21945.0 * c * c * (-221 + 45747 * Math.Pow(s, 2) - 1524900 * Math.Pow(s, 4) + 19213740 * Math.Pow(s, 6) - 119399670 * Math.Pow(s, 8) + 411265530 * Math.Pow(s, 10) - 822531060 * Math.Pow(s, 12) + 949074300 * Math.Pow(s, 14) - 585262485 * Math.Pow(s, 16) + 149184555 * Math.Pow(s, 18)) / 65536));
            LegArrEx[20, 3] = ((1514205.0 * Math.Pow(c, 3) * (663 * s - 44200 * Math.Pow(s, 3) + 835380 * Math.Pow(s, 5) - 6921720 * Math.Pow(s, 7) + 29801850 * Math.Pow(s, 9) - 71524440 * Math.Pow(s, 11) + 96282900 * Math.Pow(s, 13) - 67856520 * Math.Pow(s, 15) + 19458855 * Math.Pow(s, 17)) / 32768));
            LegArrEx[20, 4] = (77224455.0 * Math.Pow(c * c, 2) * (13 - 2600 * Math.Pow(s, 2) + 81900 * Math.Pow(s, 4) - 950040 * Math.Pow(s, 6) + 5259150 * Math.Pow(s, 8) - 15426840 * Math.Pow(s, 10) + 24542700 * Math.Pow(s, 12) - 19957800 * Math.Pow(s, 14) + 6486285 * Math.Pow(s, 16)) / 32768);
            LegArrEx[20, 5] = ((386122275.0 * Math.Pow(c, 5) * (-65 * s + 4095 * Math.Pow(s, 3) - 71253 * Math.Pow(s, 5) + 525915 * Math.Pow(s, 7) - 1928355 * Math.Pow(s, 9) + 3681405 * Math.Pow(s, 11) - 3492615 * Math.Pow(s, 13) + 1297257 * Math.Pow(s, 15)) / 2048));
            LegArrEx[20, 6] = ((25097947875.0 * Math.Pow(c * c, 3) * (-1 + 189 * Math.Pow(s, 2) - 5481 * Math.Pow(s, 4) + 56637 * Math.Pow(s, 6) - 267003 * Math.Pow(s, 8) + 623007 * Math.Pow(s, 10) - 698523 * Math.Pow(s, 12) + 299367 * Math.Pow(s, 14)) / 2048));
            LegArrEx[20, 7] = ((225881530875.0 * Math.Pow(c, 7) * (21 * s - 1218 * Math.Pow(s, 3) + 18879 * Math.Pow(s, 5) - 118668 * Math.Pow(s, 7) + 346115 * Math.Pow(s, 9) - 465682 * Math.Pow(s, 11) + 232841 * Math.Pow(s, 13)) / 1024));
            LegArrEx[20, 8] = (1581170716125.0 * Math.Pow(c * c, 4) * (3 - 522 * Math.Pow(s, 2) + 13485 * Math.Pow(s, 4) - 118668 * Math.Pow(s, 6) + 445005 * Math.Pow(s, 8) - 731786 * Math.Pow(s, 10) + 432419 * Math.Pow(s, 12)) / 1024);
            LegArrEx[20, 9] = (45853950767625.0 / 256) * Math.Pow(c, 9) * (-9 * s + 465 * Math.Pow(s, 3) - 6138 * Math.Pow(s, 5) + 30690 * Math.Pow(s, 7) - 63085 * Math.Pow(s, 9) + 44733 * Math.Pow(s, 11));
            LegArrEx[20, 10] = (137561852302875.0 / 256) * Math.Pow(c * c, 5) * (-3 + 465 * Math.Pow(s, 2) - 10230 * Math.Pow(s, 4) + 71610 * Math.Pow(s, 6) - 189255 * Math.Pow(s, 8) + 164021 * Math.Pow(s, 10));
            LegArrEx[20, 11] = (21322087106945625.0 / 128) * Math.Pow(c, 11) * (3 * s - 132 * Math.Pow(s, 3) + 1386 * Math.Pow(s, 5) - 4884 * Math.Pow(s, 7) + 5291 * Math.Pow(s, 9));
            LegArrEx[20, 12] = 63966261320836875.0 / 128 * Math.Pow(c * c, 6) * (1 - 132 * Math.Pow(s, 2) + 2310 * Math.Pow(s, 4) - 11396 * Math.Pow(s, 6) + 15873 * Math.Pow(s, 8));
            LegArrEx[20, 13] = (2110886623587616875.0 / 16) * Math.Pow(c, 13) * (-s + 35 * Math.Pow(s, 3) - 259 * Math.Pow(s, 5) + 481 * Math.Pow(s, 7));
            LegArrEx[20, 14] = (2110886623587616875.0 / 16) * Math.Pow(c * c, 7) * (-1 + 105 * Math.Pow(s, 2) - 1295 * Math.Pow(s, 4) + 3367 * Math.Pow(s, 6));
            LegArrEx[20, 15] = (14776206365113318125.0 / 8) * Math.Pow(c, 15) * (15 * s - 370 * Math.Pow(s, 3) + 1443 * Math.Pow(s, 5));
            LegArrEx[20, 16] = 221643095476699771875.0 / 8 * Math.Pow(c * c, 8) * (1 - 74 * Math.Pow(s, 2) + 481 * Math.Pow(s, 4));
            LegArrEx[20, 17] = (8200794532637891559375.0 / 2) * Math.Pow(c, 17) * (-s + 13 * Math.Pow(s, 3));
            LegArrEx[20, 18] = (8200794532637891559375.0 / 2) * Math.Pow(c * c, 9) * (-1 + 39 * Math.Pow(s, 2));
            LegArrEx[20, 19] = 319830986772877770815625.0 * s * Math.Pow(c, 19);
            LegArrEx[20, 20] = 319830986772877770815625.0 * Math.Pow(c * c, 10);

            LegArrEx[21, 0] = (1.0 / 262144) * (969969 * s - 74364290 * Math.Pow(s, 3) + 1673196525 * Math.Pow(s, 5) - 17210021400 * Math.Pow(s, 7) + 97045398450 * Math.Pow(s, 9) - 328189892940 * Math.Pow(s, 11) + 694247850450 * Math.Pow(s, 13) - 925663800600 * Math.Pow(s, 15) + 755505013725 * Math.Pow(s, 17) - 344616322050 * Math.Pow(s, 19) + 67282234305 * Math.Pow(s, 21));
            LegArrEx[21, 1] = (1.0 / 262144) * 231 * c * (4199 - 965770 * Math.Pow(s, 2) + 36216375 * Math.Pow(s, 4) - 521515800 * Math.Pow(s, 6) + 3780989550 * Math.Pow(s, 8) - 15628090140 * Math.Pow(s, 10) + 39070225350 * Math.Pow(s, 12) - 60108039000 * Math.Pow(s, 14) + 55599936075 * Math.Pow(s, 16) - 28345065450 * Math.Pow(s, 18) + 6116566755 * Math.Pow(s, 20));
            LegArrEx[21, 2] = ((26565.0 * c * c * (-4199 * s + 314925 * Math.Pow(s, 3) - 6802380 * Math.Pow(s, 5) + 65756340 * Math.Pow(s, 7) - 339741090 * Math.Pow(s, 9) + 1019223270 * Math.Pow(s, 11) - 1829375100 * Math.Pow(s, 13) + 1933910820 * Math.Pow(s, 15) - 1109154735 * Math.Pow(s, 17) + 265937685 * Math.Pow(s, 19)) / 65536));
            LegArrEx[21, 3] = ((504735.0 * Math.Pow(c, 3) * (-221 + 49725 * Math.Pow(s, 2) - 1790100 * Math.Pow(s, 4) + 24226020 * Math.Pow(s, 6) - 160929990 * Math.Pow(s, 8) + 590076630 * Math.Pow(s, 10) - 1251677700 * Math.Pow(s, 12) + 1526771700 * Math.Pow(s, 14) - 992401605 * Math.Pow(s, 16) + 265937685 * Math.Pow(s, 18)) / 65536));
            LegArrEx[21, 4] = (22713075.0 * Math.Pow(c * c, 2) * (1105 * s - 79560 * Math.Pow(s, 3) + 1615068 * Math.Pow(s, 5) - 14304888 * Math.Pow(s, 7) + 65564070 * Math.Pow(s, 9) - 166890360 * Math.Pow(s, 11) + 237497820 * Math.Pow(s, 13) - 176426952 * Math.Pow(s, 15) + 53187537 * Math.Pow(s, 17)) / 32768);
            LegArrEx[21, 5] = ((5019589575.0 * Math.Pow(c, 5) * (5 - 1080 * Math.Pow(s, 2) + 36540 * Math.Pow(s, 4) - 453096 * Math.Pow(s, 6) + 2670030 * Math.Pow(s, 8) - 8306760 * Math.Pow(s, 10) + 13970460 * Math.Pow(s, 12) - 11974680 * Math.Pow(s, 14) + 4091349 * Math.Pow(s, 16)) / 32768));
            LegArrEx[21, 6] = ((15058768725.0 * Math.Pow(c * c, 3) * (-45 * s + 3045 * Math.Pow(s, 3) - 56637 * Math.Pow(s, 5) + 445005 * Math.Pow(s, 7) - 1730575 * Math.Pow(s, 9) + 3492615 * Math.Pow(s, 11) - 3492615 * Math.Pow(s, 13) + 1363783 * Math.Pow(s, 15)) / 2048));
            LegArrEx[21, 7] = ((225881530875.0 * Math.Pow(c, 7) * (-3 + 609 * Math.Pow(s, 2) - 18879 * Math.Pow(s, 4) + 207669 * Math.Pow(s, 6) - 1038345 * Math.Pow(s, 8) + 2561251 * Math.Pow(s, 10) - 3026933 * Math.Pow(s, 12) + 1363783 * Math.Pow(s, 14)) / 2048));
            LegArrEx[21, 8] = (45853950767625.0 * Math.Pow(c * c, 4) * (3 * s - 186 * Math.Pow(s, 3) + 3069 * Math.Pow(s, 5) - 20460 * Math.Pow(s, 7) + 63085 * Math.Pow(s, 9) - 89466 * Math.Pow(s, 11) + 47027 * Math.Pow(s, 13)) / 1024);
            LegArrEx[21, 9] = ((45853950767625.0 * Math.Pow(c, 9) * (3 - 558 * Math.Pow(s, 2) + 15345 * Math.Pow(s, 4) - 143220 * Math.Pow(s, 6) + 567765 * Math.Pow(s, 8) - 984126 * Math.Pow(s, 10) + 611351 * Math.Pow(s, 12)) / 1024));
            LegArrEx[21, 10] = (4264417421389125.0 / 256) * Math.Pow(c * c, 5) * (-3 * s + 165 * Math.Pow(s, 3) - 2310 * Math.Pow(s, 5) + 12210 * Math.Pow(s, 7) - 26455 * Math.Pow(s, 9) + 19721 * Math.Pow(s, 11));
            LegArrEx[21, 11] = (4264417421389125.0 / 256) * Math.Pow(c, 11) * (-3 + 495 * Math.Pow(s, 2) - 11550 * Math.Pow(s, 4) + 85470 * Math.Pow(s, 6) - 238095 * Math.Pow(s, 8) + 216931 * Math.Pow(s, 10));
            LegArrEx[21, 12] = 234542958176401875.0 / 128 * Math.Pow(c * c, 6) * (9 * s - 420 * Math.Pow(s, 3) + 4662 * Math.Pow(s, 5) - 17316 * Math.Pow(s, 7) + 19721 * Math.Pow(s, 9));
            LegArrEx[21, 13] = (2110886623587616875.0 / 128) * Math.Pow(c, 13) * (1 - 140 * Math.Pow(s, 2) + 2590 * Math.Pow(s, 4) - 13468 * Math.Pow(s, 6) + 19721 * Math.Pow(s, 8));
            LegArrEx[21, 14] = (2110886623587616875.0 / 16) * Math.Pow(c * c, 7) * (-35 * s + 1295 * Math.Pow(s, 3) - 10101 * Math.Pow(s, 5) + 19721 * Math.Pow(s, 7));
            LegArrEx[21, 15] = (14776206365113318125.0 / 16) * Math.Pow(c, 15) * (-5 + 555 * Math.Pow(s, 2) - 7215 * Math.Pow(s, 4) + 19721 * Math.Pow(s, 6));
            LegArrEx[21, 16] = 1640158906527578311875.0 / 8 * Math.Pow(c * c, 8) * (5 * s - 130 * Math.Pow(s, 3) + 533 * Math.Pow(s, 5));
            LegArrEx[21, 17] = (8200794532637891559375.0 / 8) * Math.Pow(c, 17) * (1 - 78 * Math.Pow(s, 2) + 533 * Math.Pow(s, 4));
            LegArrEx[21, 18] = (106610328924292590271875.0 / 2) * Math.Pow(c * c, 9) * (-3 * s + 41 * Math.Pow(s, 3));
            LegArrEx[21, 19] = (319830986772877770815625.0 / 2) * Math.Pow(c, 19) * (-1 + 41 * Math.Pow(s, 2));
            LegArrEx[21, 20] = 13113070457687988603440625.0 * s * Math.Pow(c * c, 10);
            LegArrEx[21, 21] = 13113070457687988603440625.0 * Math.Pow(c, 21);

        } // LegPolyEx





        // ----------------------------------------------------------------------------
        // fukushima method (JG 2018)
        //   for very large spherical harmonic expansions and calcs of normalized associated 
        //   Legendre polynomials
        //   
        //   Plm are converted to X-numbers
        //   Clm, Slm treated as F-numbers
        // -----------------------------------------------------------------------------

        // initialize legendre function values
        public void pinit
            (
            Int32 n,
            Int32 m,
            ref double[] p
            )
        {
            p = new double[360];

            if (n == 0)
                p[0] = 1.0;
            else if (n == 1)
                p[0] = 1.7320508075688773;
            else if (n == 2)
            {
                if (m == 0)
                {
                    p[0] = 0.5590169943749474;
                    p[1] = 1.6770509831248423;
                }
                else if (m == 1)
                {
                    p[0] = 0.0;
                    p[1] = 1.9364916731037084;
                }
                else if (m == 2)
                {
                    p[0] = 0.9682458365518542;
                    p[1] = -0.9682458365518542;
                }
                else if (n == 3)
                {
                    if (m == 0)
                    {
                        p[0] = 0.9921567416492215;
                        p[1] = 1.6535945694153691;
                    }
                    else if (m == 1)
                    {
                        p[0] = 0.4050462936504913;
                        p[1] = 2.0252314682524563;
                    }
                    else if (m == 2)
                    {
                        p[0] = 1.2808688457449498;
                        p[1] = -1.2808688457449498;
                    }
                    else if (m == 3)
                    {
                        p[0] = 1.5687375497513917;
                        p[1] = -0.5229125165837972;
                    }
                }
                else if (n == 4)
                {
                    if (m == 0)
                    {
                        p[0] = 0.421875;
                        p[1] = 0.9375;
                        p[2] = 1.640625;
                    }
                    else if (m == 1)
                    {
                        p[0] = 0.0;
                        p[1] = 0.5929270612815711;
                        p[2] = 2.0752447144854989;
                    }
                    else if (m == 2)
                    {
                        p[0] = 0.6288941186718159;
                        p[1] = 0.8385254915624211;
                        p[2] = -1.4674196102342370;
                    }
                    else if (m == 3)
                    {
                        p[0] = 0.0;
                        p[1] = 1.5687375497513917;
                        p[2] = -0.7843687748756958;
                    }
                    else if (m == 4)
                    {
                        p[0] = 0.8319487194983835;
                        p[1] = -1.1092649593311780;
                        p[2] = 0.2773162398327945;
                    }
                }
            }
        }  // pinit


        /* ----------------------------------------------------------------------------
        *      x2f
        *      
        *      convert x to f number
        *      
           ------------------------------------------------------------------------------*/
        public double x2f
            (
              double x,
              Int32 ix
            )
        {
            double x2fv;
            Int32 IND;
            double BIG, BIGI;

            IND = 960;
            BIG = Math.Pow(2.0, IND);
            BIGI = Math.Pow(2.0, -IND);

            if (ix == 0)
                x2fv = x;
            else if (ix == -1)
                x2fv = x * BIGI;
            else if (ix == 1)
                x2fv = x * BIG;
            else if (ix < 0)
                x2fv = 0.0;
            else if (x < 0)
                x2fv = -BIG;
            else
                x2fv = BIG;

            return x2fv;
        }


        /* ----------------------------------------------------------------------------
         *                                  xnorm
         *                                  
         * uses the "x" factor approach - value and exponent
         * 
          ------------------------------------------------------------------------------*/

        public void xnorm
                (
                ref double x,
                ref Int32 ix
                )
        {
            Int32 IND = 960;
            double w;

            double BIG = Math.Pow(2, IND);
            double BIGI = Math.Pow(2, -IND);
            double BIGS = Math.Pow(2.0, 480);  // IND / 2
            double BIGSI = Math.Pow(2.0, -480);  // IND / 2

            w = Math.Abs(x);
            if (w >= BIGS)
            {
                x = x * BIGI;
                ix = ix + 1;
            }
            else if (w < BIGSI)
            {
                x = x * BIG;
                ix = ix - 1;
            }
        }  // xnorm


        /* ----------------------------------------------------------------------------
         *                                       xl2sum
         *                                       
         * routine to compute the two-term linear sum of X-numbers
         * with F-number coefficients
         * 
          ---------------------------------------------------------------------------- */
        public void xlsum2
                (
                   double f,
                   double g,
                   double x,
                   double y,
                   out double z,
                   Int32 ix,
                   Int32 iy,
                   out Int32 iz)
        {
            Int32 id;
            Int32 IND = 960;
            double BIGI = Math.Pow(2, -IND);

            id = ix - iy;
            if (id == 0)
            {
                z = f * x + g * y;
                iz = ix;
            }
            else if (id == 1)
            {
                z = f * x + g * (y * BIGI);
                iz = ix;
            }
            else if (id == -1)
            {
                z = g * y + f * (x * BIGI);
                iz = iy;
            }
            else if (id > 1)
            {
                z = f * x;
                iz = ix;
            }
            else
            {
                z = g * y;
                iz = iy;
            }

            xnorm(ref z, ref iz);
        }  // xlsum2


        /* ----------------------------------------------------------------------------
         *                                  dpeven
         *                                  
         * find Pnnj and Pn,n-1,j when degree n is even and n >= 6. the returned values are 
         * (xp[j], ip[j]) and (xp1[j], ip1[j]), double X-number vectors representing
         * Pnnj and Pn,n-1,j, respectively. required initial values for Pn-2,n-2,j as 
         * (xpold[j], ipold[j]) are needed. 
         * 
         *  inputs          description                                  range / units
         *    n           - degree 
         *
         *
         *
         *  outputs       :
         *    xp1         - value
         *    ip1         - exponent
         *
         *  locals        :
         *                -
         *
         *  coupling      :
         *    xnorm 
         *    xlsum2  
         *
         *  references    :
         * Fukushima (2012a)
           ---------------------------------------------------------------------------- */

        public void dpeven
            (
            Int32 n,
            double[] xpold,
            out double[] xp,
            out double[] xp1,
            Int32[] ipold,
            out Int32[] ip,
            out Int32[] ip1
            )
        {
            xpold = new double[360];
            xp = new double[360];
            xp1 = new double[360];
            ipold = new Int32[360];
            ip = new Int32[360];
            ip1 = new Int32[360];

            Int32 jx, jxm2, jxm1, n2, j, itemp, jm1, jp1;
            double gamma, gamma2, xtemp, alpha2;

            jx = n / 2;
            jxm2 = jx - 2;
            jxm1 = jx - 1;
            n2 = n * 2;
            gamma = Math.Sqrt(Convert.ToDouble(n2 + 1) * Convert.ToDouble(n2 - 1) / (Convert.ToDouble(n) * Convert.ToDouble(n - 1))) * 0.125;
            gamma2 = gamma * 2.0;
            xlsum2(gamma2, xpold[0], -gamma, xpold[1], out xp[0], ipold[0], ipold[1], out ip[0]);
            xlsum2(-gamma2, xpold[0], gamma2, xpold[1], out xtemp, ipold[0], ipold[1], out itemp);
            xlsum2(1.0, xtemp, -gamma, xpold[2], out xp[1], itemp, ipold[2], out ip[1]);
            j = 2;
            while (j <= jxm2)
            {
                jm1 = j - 1;
                jp1 = j + 1;
                xlsum2(-gamma, xpold[jm1], gamma2, xpold[j], out xtemp, ipold[jm1], ipold[j], out itemp);
                xlsum2(1.0, xtemp, -gamma, xpold[jp1], out xp[j], itemp, ipold[jp1], out ip[j]);
                j = j + 1;
            }
            xlsum2(-gamma, xpold[jxm2], gamma2, xpold[jxm1], out xp[jxm1], ipold[jxm2], ipold[jxm1], out ip[jxm1]);
            xp[jx] = -gamma * xpold[jxm1];
            ip[jx] = ipold[jxm1];
            xnorm(ref xp[jx], ref ip[jx]);
            alpha2 = Math.Sqrt(2.0 / Convert.ToDouble(n)) * 2.0;
            xp1[0] = 0.0;
            ip1[0] = 0;
            j = 1;
            while (j <= jx)
            {
                xp1[j] = -Convert.ToDouble(j) * alpha2 * xp[j];
                ip1[j] = ip[j];
                xnorm(ref xp1[j], ref ip1[j]);
                j = j + 1;
            }
        }   // dpeven


        /* ----------------------------------------------------------------------------
         *                                 dpodd
         *                                 
         *  routine to return Pnnj and Pn,n-1,j when n is odd and n >= 5. Same as Table 4 
         *  but when n is odd and n >= 5.
         ---------------------------------------------------------------------------- */
        public void dpodd
            (
            Int32 n,
            double[] xpold,
            out double[] xp,
            out double[] xp1,
            Int32[] ipold,
            out Int32[] ip,
            out Int32[] ip1
     )
        {
            xpold = new double[360];
            xp = new double[360];
            xp1 = new double[360];
            ipold = new Int32[360];
            ip = new Int32[360];
            ip1 = new Int32[360];

            Int32 jx, jxm2, jxm1, n2, j, itemp, jm1, jp1;
            double gamma, gamma2, xtemp, alpha;

            jx = (n - 1) / 2;
            jxm2 = jx - 2;
            jxm1 = jx - 1;
            n2 = n * 2;
            gamma = Math.Sqrt(Convert.ToDouble(n2 + 1) * Convert.ToDouble(n2 - 1) / (Convert.ToDouble(n) * Convert.ToDouble(n - 1))) * 0.125;
            gamma2 = gamma * 2.0;
            xlsum2(gamma * 3.0, xpold[0], -gamma, xpold[1], out xp[0], ipold[0], ipold[1], out ip[0]);
            j = 1;
            while (j <= jxm2)
            {
                jm1 = j - 1;
                jp1 = j + 1;
                xlsum2(-gamma, xpold[jm1], gamma2, xpold[j], out xtemp, ipold[jm1], ipold[j], out itemp);
                xlsum2(1.0, xtemp, -gamma, xpold[jp1], out xp[j], itemp, ipold[jp1], out ip[j]);
                j = j + 1;
            }
            xlsum2(-gamma, xpold[jxm2], gamma2, xpold[jxm1], out xp[jxm1], ipold[jxm2], ipold[jxm1], out ip[jxm1]);
            xp[jx] = -gamma * xpold[jxm1];
            ip[jx] = ipold[jxm1];
            xnorm(ref xp[jx], ref ip[jx]);
            alpha = Math.Sqrt(2.0 / Convert.ToDouble(n));
            j = 0;
            while (j <= jx)
            {
                xp1[j] = Convert.ToDouble(2 * j + 1) * alpha * xp[j];
                ip1[j] = ip[j];
                xnorm(ref xp1[j], ref ip1[j]);
                j = j + 1;
            }
        }   // dpodd


        /* ----------------------------------------------------------------------------
         *                             gpeven
         *                             
         * routine to return Pnmj when n is even. The returned values are (xp0[j], ip0[j]), 
         * a double X-number vector representing Pnmj. We assume that Pn,m+1,j and Pn,m+2,j are
         * externally provided as (xp1[j], ip1[j]) and (xp2[j], ip2[j]), respectively.
         * The routine internally calls xnorm and xlsum2 provided in Tables 7 and 8 of 
         * Fukushima (2012a).
          ---------------------------------------------------------------------------- */

        public void gpeven
            (
            Int32 jmax,
            Int32 n,
            Int32 m,
            double[] xp2,
            double[] xp1,
            out double[] xp0,
            Int32[] ip2,
            Int32[] ip1,
            out Int32[] ip0
            )
        {
            xp2 = new double[360];
            xp1 = new double[360];
            xp0 = new double[360];
            ip2 = new Int32[360];
            ip1 = new Int32[360];
            ip0 = new Int32[360];

            Int32 j, m1, m2, modd;
            double u, alpha2, beta;

            m1 = m + 1;
            m2 = m + 2;
            modd = m - Convert.ToInt32(m * 0.5) * 2;
            if (m == 0)
                u = Math.Sqrt(0.5 / (Convert.ToDouble(n) * Convert.ToDouble(n + 1)));
            else
                u = Math.Sqrt(1.0 / (Convert.ToDouble(n - m) * Convert.ToDouble(n + m1)));

            alpha2 = 4.0 * u;
            beta = Math.Sqrt(Convert.ToDouble(n - m1) * Convert.ToDouble(n + m2)) * u;
            xp0[0] = beta * xp2[0];
            ip0[0] = ip2[0];
            xnorm(ref xp0[0], ref ip0[0]);
            if (modd == 0)
            {
                j = 1;
                while (j <= jmax)
                {
                    xlsum2(Convert.ToDouble(j) * alpha2, xp1[j], beta, xp2[j], out xp0[j], ip1[j], ip2[j], out ip0[j]);
                    j = j + 1;
                }
            }
            else
            {
                j = 1;
                while (j <= jmax)
                {
                    xlsum2(-Convert.ToDouble(j) * alpha2, xp1[j], beta, xp2[j], out xp0[j], ip1[j], ip2[j], out ip0[j]);
                    j = j + 1;
                }
            }
        }  // gpeven


        /* ----------------------------------------------------------------------------
        *                                 gpodd
        *                                 
        * Table 7: Fortran routine to return Pnmj when n is odd. Same as
        * Table 6 but when n is odd.
         ---------------------------------------------------------------------------- */
        public void gpodd
            (
            Int32 jmax,
            Int32 n,
            Int32 m,
            double[] xp2,
            double[] xp1,
            out double[] xp0,
            Int32[] ip2,
            Int32[] ip1,
            out Int32[] ip0
            )
        {
            xp2 = new double[360];
            xp1 = new double[360];
            xp0 = new double[360];
            ip2 = new Int32[360];
            ip1 = new Int32[360];
            ip0 = new Int32[360];

            Int32 j, m1, m2, modd;
            double u, alpha, beta;

            m1 = m + 1;
            m2 = m + 2;
            modd = m - Convert.ToInt32(m * 0.5) * 2;
            if (m == 0)
                u = Math.Sqrt(0.50 / (Convert.ToDouble(n) * Convert.ToDouble(n + 1)));
            else
                u = Math.Sqrt(1.0 / (Convert.ToDouble(n - m) * Convert.ToDouble(n + m1)));

            alpha = 2.0 * u;
            beta = Math.Sqrt(Convert.ToDouble(n - m1) * Convert.ToDouble(n + m2)) * u;
            if (modd == 0)
            {
                j = 0;
                while (j <= jmax)
                {
                    xlsum2(Convert.ToDouble(2 * j + 1) * alpha, xp1[j], beta, xp2[j], out xp0[j], ip1[j], ip2[j], out ip0[j]);
                    j = j + 1;
                }
            }
            else
            {
                j = 0;
                while (j <= jmax)
                {
                    xlsum2(-Convert.ToDouble(2 * j + 1) * alpha, xp1[j], beta, xp2[j], out xp0[j], ip1[j], ip2[j], out ip0[j]);
                    j = j + 1;
                }
            }
        }  // gpodd



        // Fukushima combined approach to find matrix of normalized Legendre polynomials
        //
        //

        public void LegPolyFF
            (
                double[] recef,
                double latgc,
                Int32 order,
                char normalized,
                double[,,] normArr,
                AstroLib.gravityConst gravData,
                out double[,] ALFArr
            )
        {
            Int32 L, m;
            double x, y, z, f, g;

            ALFArr = new double[order + 2, order + 2];
            double magr = MathTimeLibr.mag(recef);

            // initial values
            ALFArr[0, 0] = 1.0;
            ALFArr[0, 1] = 0.0;
            ALFArr[1, 0] = Math.Sin(latgc);
            ALFArr[1, 1] = Math.Cos(latgc);
            m = 2;
            L = m + 1;
            x = ALFArr[1, 0];
            y = ALFArr[1, 1];

            // find zonal terms
            for (L = 2; L <= order + 1; L++)
            {
                // find tesseral and sectoral terms
                //               for (m = 0; m <= order + 1; m++)
                {
                    f = normArr[L, m, 0] * Math.Sin(latgc);
                    g = -normArr[L, m, 1];
                    z = f * x + g * y;
                    ALFArr[L, m] = z;
                    y = x;
                    x = z;
                }
            }

        }  // LegPolyFF



        /* ----------------------------------------------------------------------------
        *                                   xfsh2f
        *
        *  transfrom the Cnm, Snm, the 4 fully normalized spherical harmonic coefficients 
        *  of a given function depend on the spherical surface, to (Akm, Bkm), the 
        *  corresponding Fourier series coefficients of the function. in the program, 
        *  (i) x2f and xnorm are the Fortran function/routine to handle X-numbers
        *      (Fukushima, 2012a, tables 6 and 7), and 
        *  (ii) pinit, dpeven, dpodd, gpeven, and gpodd are the Fortran routines listed 
        *  in Tables 3–7, respectively.
        *  
        *    
        *  
        *    Fukushima 2018 table xx
          ---------------------------------------------------------------------------- */

        public void xfsh2f
            (
            Int32 nmax,
            AstroLib.gravityConst gravData,
            out double[,] a,
            out double[,] b
            )
        {
            Int32 NX = 100;  // 2200;
            Int32 j, m, k, L, jmax, n1;
            a = new double[NX, NX];
            b = new double[NX, NX];
            //     double[,] pja = new double[NX, NX];  // test to see the values
            jmax = 0;

            Int32[] ipold, ip, ip0;
            Int32[] ip1, ip2;
            double[] xpold, xp, xp0;
            double[] xp1, xp2;
            double pj;
            xpold = new double[NX];
            xp = new double[NX];
            xp0 = new double[NX];
            xp1 = new double[NX];
            xp2 = new double[NX];
            ipold = new Int32[NX];
            ip = new Int32[NX];
            ip0 = new Int32[NX];
            ip1 = new Int32[NX];
            ip2 = new Int32[NX];

            // initialize all to 0
            for (m = 0; m <= nmax; m++)
            {
                for (k = 0; k <= nmax; k++)
                {
                    a[k, m] = 0.0;
                    b[k, m] = 0.0;
                }
            }

            // initialize the first 4x4 values
            for (L = 0; L <= 4; L += 2)
            {
                jmax = Convert.ToInt32(L * 0.5);

                for (m = 0; m <= L; m++)
                {
                    pinit(L, m, ref xp);

                    for (j = 0; j <= jmax; j++)
                    {
                        k = 2 * j;
                        a[k, m] = a[k, m] + xp[j] * gravData.c[L, m];
                        b[k, m] = b[k, m] + xp[j] * gravData.s[L, m];
                    }
                }
            }

            // even terms
            for (j = 0; j <= jmax; j++)
            {
                ip[j] = 0;
            }

            for (L = 6; L <= nmax; L += 2)
            {

                for (j = 0; j <= jmax; j++)
                {
                    xpold[j] = xp[j];
                    ipold[j] = ip[j];
                }
                jmax = Convert.ToInt32(L * 0.5);
                n1 = L - 1;
                dpeven(L, xpold, out xp, out xp1, ipold, out ip, out ip1);

                for (j = 0; j <= jmax; j++)
                {
                    k = 2 * j;
                    pj = x2f(xp[j], ip[j]);
                    //       pja[k, n1] = pj;
                    a[k, L] = a[k, L] + pj * gravData.c[L, L];
                    b[k, L] = b[k, L] + pj * gravData.s[L, L];
                    pj = x2f(xp1[j], ip1[j]);
                    a[k, n1] = a[k, n1] + pj * gravData.c[L, n1];
                    b[k, n1] = b[k, n1] + pj * gravData.s[L, n1];
                    xp2[j] = xp[j];
                    ip2[j] = ip[j];
                }

                for (m = L - 2; m >= 0; m -= 1)
                {
                    gpeven(jmax, L, m, xp2, xp1, out xp0, ip2, ip1, out ip0);

                    for (j = 0; j <= jmax; j++)
                    {
                        k = 2 * j;
                        pj = x2f(xp0[j], ip0[j]);
                        //         pja[k, m] = pj;
                        a[k, m] = a[k, m] + pj * gravData.c[L, m];
                        b[k, m] = b[k, m] + pj * gravData.s[L, m];
                        xp2[j] = xp1[j];
                        ip2[j] = ip1[j];
                        xp1[j] = xp0[j];
                        ip1[j] = ip0[j];
                    }
                }
            }

            for (L = 1; L <= 3; L += 2)
            {
                jmax = Convert.ToInt32((L - 1) * 0.5);

                for (m = 0; m <= L; m++)
                {
                    pinit(L, m, ref xp);

                    for (j = 0; j <= jmax; j++)
                    {
                        k = 2 * j + 1;
                        a[k, m] = a[k, m] + xp[j] * gravData.c[L, m];
                        b[k, m] = b[k, m] + xp[j] * gravData.s[L, m];
                    }
                }
            }

            // odd terms
            for (j = 0; j <= jmax; j++)
            {
                ip[j] = 0;
            }

            for (L = 5; L <= nmax; L += 2)
            {

                for (j = 0; j <= jmax; j++)
                {
                    xpold[j] = xp[j];
                    ipold[j] = ip[j];
                }
                jmax = Convert.ToInt32((L - 1) * 0.5);
                n1 = L - 1;
                dpodd(L, xpold, out xp, out xp1, ipold, out ip, out ip1);

                for (j = 0; j <= jmax; j++)
                {
                    k = 2 * j + 1;
                    pj = x2f(xp[j], ip[j]);
                    a[k, L] = a[k, L] + pj * gravData.c[L, L];
                    b[k, L] = b[k, L] + pj * gravData.s[L, L];
                    //   pja[k, n] = pj;
                    pj = x2f(xp1[j], ip1[j]);
                    a[k, n1] = a[k, n1] + pj * gravData.c[L, n1];
                    b[k, n1] = b[k, n1] + pj * gravData.s[L, n1];
                    xp2[j] = xp[j];
                    ip2[j] = ip[j];
                }

                for (m = L - 2; m >= 0; m -= 1)
                {
                    gpodd(jmax, L, m, xp2, xp1, out xp0, ip2, ip1, out ip0);

                    for (j = 0; j <= jmax; j++)
                    {
                        k = 2 * j + 1;
                        pj = x2f(xp0[j], ip0[j]);
                        //     pja[k, m] = pj;
                        a[k, m] = a[k, m] + pj * gravData.c[L, m];
                        b[k, m] = b[k, m] + pj * gravData.s[L, m];
                        xp2[j] = xp1[j];
                        ip2[j] = ip1[j];
                        xp1[j] = xp0[j];
                        ip1[j] = ip0[j];
                    }
                }
            }
        }  // xfsh2f

        // from 2017 fukishima paper, ALFs
        //        j n = 2j Pn(1=2) d Pn(1=2)
        //1 2 􀀀2.7950849718747371205114670859140954E􀀀01 +4.46E􀀀16
        //2 4 􀀀8.6718750000000000000000000000000019E􀀀01 +1.28E􀀀16
        //3 8 􀀀3.0362102888840987874508856147660683E􀀀01 +4.89E􀀀16
        //4 16 􀀀8.6085221363787000197086857609730105E􀀀01 􀀀2.47E􀀀16
        //5 32 􀀀3.1119962497760147366174972709413071E􀀀01 􀀀1.14E􀀀15
        //6 64 􀀀8.5832418550243444685033693028444350E􀀀01 􀀀1.97E􀀀16
        //7 128 􀀀3.1316449107909472026965626279232704E􀀀01 􀀀2.13E􀀀15
        //8 256 􀀀8.5762286823362136598509949562358587E􀀀01 +4.73E􀀀16
        //9 512 􀀀3.1365884210353617421696699741828314E􀀀01 +1.94E􀀀15
        //10 1024 􀀀8.5744308441362078316451963835453876E􀀀01 􀀀1.12E􀀀15
        //11 2048 􀀀3.1378260213136415436393598399168296E􀀀01 +1.01E􀀀14


        // GMAT Pines approach
        //------------------------------------------------------------------------------
        public void FullGeopGMATPines
            (
            double jday,
            double[] pos,
            double latgc,
            Int32 nn, Int32 mm,
            AstroLib.gravityConst gravData,
            out double[] acc
            //out double[,] gradient
            )
        {
            acc = new double[3];

            // Int32 XS = fillgradient ? 2 : 1;
            // calculate vector components ----------------------------------
            double magr = Math.Sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);    // Naming scheme from ref [3]
            double s = pos[0] / magr;
            double t = pos[1] / magr;
            double u = pos[2] / magr; // sin(phi), phi = geocentric latitude

            // Calculate values for A -----------------------------------------
            int ord = 750;
            double[,] A = new double[ord, ord];
            double[,] N1 = new double[ord, ord];
            double[,] N2 = new double[ord, ord];
            double[,] V = new double[ord, ord];
            double[,] VR01 = new double[ord, ord];
            double[,] VR11 = new double[ord, ord];
            double[,] VR02 = new double[ord, ord];
            double[,] VR12 = new double[ord, ord];
            double[,] VR22 = new double[ord, ord];
            double[] Re = new double[ord];
            double[] Im = new double[ord];
            double sqrt2 = Math.Sqrt(2.0);
            Int32 XS = 2;
            u = Math.Sin(latgc);

            // get leg poly normalization numbers (do once)
            // all 0
            for (Int32 m = 0; m <= mm + 2; ++m)
            {
                for (Int32 L = m + 2; L <= nn + 2; ++L)
                {
                    N1[L, m] = Math.Sqrt(((2.0 * L + 1) * (2 * L - 1)) / ((L - m) * (L + m)));  // double in denom??
                    N2[L, m] = Math.Sqrt(((2.0 * L + 1) * (L - m - 1) * (L + m - 1)) / ((2.0 * L - 3) * (L + m) * (L - m)));
                }
            }

            // NANs
            for (Int32 L = 0; L <= nn + 2; ++L)
            {
                V[L, 0] = Math.Sqrt((2.0 * (2 * L + 1)));   // Temporary, to make following loop work
                for (Int32 m = 1; m <= L + 2 && m <= mm + 2; ++m)
                {
                    V[L, m] = V[L, m - 1] / Math.Sqrt(((L + m) * (L - m + 1)));  // need real on L-m?
                }
                V[L, 0] = Math.Sqrt((2.0 * L + 1));       // Now set true value
            }

            for (Int32 L = 0; L <= nn; ++L)
                for (Int32 m = 0; m <= L && m <= mm; ++m)
                {
                    //double nn = L;
                    VR01[L, m] = Math.Sqrt(((nn - m) * (nn + m + 1)));  // need real on L-m?
                    VR11[L, m] = Math.Sqrt(((2.0 * nn + 1) * (nn + m + 2) * (nn + m + 1)) / ((2.0 * nn + 3)));
                    VR02[L, m] = Math.Sqrt(((nn - m) * (nn - m - 1) * (nn + m + 1) * (nn + m + 2)));  // need real on L-m?
                    VR12[L, m] = Math.Sqrt((2.0 * nn + 1) / (2.0 * nn + 3) * ((nn - m) * (nn + m + 1) * (nn + m + 2) * (nn + m + 3)));  // need real on L-m?
                    VR22[L, m] = Math.Sqrt((2.0 * nn + 1) / (2.0 * nn + 5) * ((nn + m + 1.0) * (nn + m + 2) * (nn + m + 3) * (nn + m + 4)));
                    if (m == 0)
                    {
                        VR01[L, m] /= sqrt2;
                        VR11[L, m] /= sqrt2;
                        VR02[L, m] /= sqrt2;
                        VR12[L, m] /= sqrt2;
                        VR22[L, m] /= sqrt2;
                    }
                }

            // generate legendre polynomials - the off-diagonal elements
            A[1, 0] = u * Math.Sqrt(3.0);
            for (Int32 L = 1; L <= nn + XS; ++L)
                A[L + 1, L] = u * Math.Sqrt(2.0 * L + 3) * A[L, L];

            // apply column-fill recursion formula (Table 2, Row I, Ref.[1])
            for (Int32 m = 0; m <= mm + XS; ++m)
            {
                for (Int32 L = m + 2; L <= nn + XS; ++L)
                    A[L, m] = u * N1[L, m] * A[L - 1, m] - N2[L, m] * A[L - 2, m];  // uses anm bnm from fukushima eq 6, 7
                // Ref.[3], Eq.(24)
                if (m == 0)
                    Re[m] = 1;
                else
                    Re[m] = s * Re[m - 1] - t * Im[m - 1]; // real part of (s + i*t)^m
                if (m == 0)
                    Im[m] = 0;
                else
                    Im[m] = s * Im[m - 1] + t * Re[m - 1]; // imaginary part of (s + i*t)^m
            }

            // Now do summation ------------------------------------------------
            // initialize recursion
            double FieldRadius = AstroLibr.gravConst.re;
            double rho = FieldRadius / magr;
            double Factor = AstroLibr.gravConst.mu;
            double rho_np1 = -Factor / magr * rho;   // rho(0) ,Ref[3], Eq 26 , factor = mu for gravity
            double rho_np2 = rho_np1 * rho;
            double a1 = 0;
            double a2 = 0;
            double a3 = 0;
            double a4 = 0;
            for (Int32 L = 1; L <= nn; ++L)
            {
                rho_np1 *= rho;
                rho_np2 *= rho;
                double sum1 = 0;
                double sum2 = 0;
                double sum3 = 0;
                double sum4 = 0;

                for (Int32 m = 0; m <= L && m <= mm; ++m) // wcs - removed "m<=L"
                {
                    double Cval = gravData.c[L, m]; // Cnm(jday, L, m);
                    double Sval = gravData.s[L, m]; // Snm(jday, L, m);
                    // Pines Equation 27 (Part of)
                    double D = (Cval * Re[m] + Sval * Im[m]) * sqrt2;
                    double E, F;
                    if (m == 0)
                        E = 0;
                    else
                        E = (Cval * Re[m - 1] + Sval * Im[m - 1]) * sqrt2;
                    if (m == 0)
                        F = 0;
                    else
                        F = (Sval * Re[m - 1] - Cval * Im[m - 1]) * sqrt2;

                    // Correct for normalization
                    double Avv00 = A[L, m];
                    double Avv01 = VR01[L, m] * A[L, m + 1];
                    double Avv11 = VR11[L, m] * A[L + 1, m + 1];
                    // Pines Equation 30 and 30b (Part of)
                    sum1 += m * Avv00 * E;
                    sum2 += m * Avv00 * F;
                    sum3 += Avv01 * D;
                    sum4 += Avv11 * D;

                    // Truncate the gradient at GRADIENT_MAX x GRADIENT_MAX
                    //if (fillgradient)
                    //{
                    //    if ((m <= gradientlimit) && (L <= gradientlimit))
                    //    {
                    //        // Pines Equation 27 (Part of)
                    //        // 2015.09.18 GMT-5295 m<=2  -> m<=1
                    //        double G = m <= 1 ? 0 : (Cval * AstroLibr.gravConst.[m - 2] + Sval * Im[m - 2]) * sqrt2;
                    //        double H = m <= 1 ? 0 : (Sval * AstroLibr.gravConst.[m - 2] - Cval * Im[m - 2]) * sqrt2;
                    //        // Correct for normalization
                    //        double Avv02 = VR02[L][m] * A[L][m + 2];
                    //        double Avv12 = VR12[L][m] * A[L + 1][m + 2];
                    //        double Avv22 = VR22[L][m] * A[L + 2][m + 2];
                    //        if (GmatMathUtil::IsNaN(Avv02) || GmatMathUtil::IsInf(Avv02))
                    //            Avv02 = 0.0;  // ************** wcs added ****

                    //        // Pines Equation 36 (Part of)
                    //        sum11 += m * (m - 1) * Avv00 * G;
                    //        sum12 += m * (m - 1) * Avv00 * H;
                    //        sum13 += m * Avv01 * E;
                    //        sum14 += m * Avv11 * E;
                    //        sum23 += m * Avv01 * F;
                    //        sum24 += m * Avv11 * F;
                    //        sum33 += Avv02 * D;
                    //        sum34 += Avv12 * D;
                    //        sum44 += Avv22 * D;
                    //    }
                    //    else
                    //    {
                    //        if (matrixTruncationWasPosted == false)
                    //        {
                    //            MessageInterface::ShowMessage("*** WARNING *** Gradient data "

                    //                  "for the state transition matrix and A-matrix "

                    //                  "computations are truncated at degree and order "

                    //                  "<= %d.\L", gradientlimit);
                    //            matrixTruncationWasPosted = true;
                    //        }
                    //    }
                    //}
                }
                // Pines Equation 30 and 30b (Part of)
                double rr = rho_np1 / FieldRadius;
                a1 += rr * sum1;
                a2 += rr * sum2;
                a3 += rr * sum3;
                a4 -= rr * sum4;
                //if (fillgradient)
                //{
                //    // Pines Equation 36 (Part of)
                //    a11 += rho_np2 / FieldRadius / FieldRadius * sum11;
                //    a12 += rho_np2 / FieldRadius / FieldRadius * sum12;
                //    a13 += rho_np2 / FieldRadius / FieldRadius * sum13;
                //    a14 -= rho_np2 / FieldRadius / FieldRadius * sum14;
                //    a23 += rho_np2 / FieldRadius / FieldRadius * sum23;
                //    a24 -= rho_np2 / FieldRadius / FieldRadius * sum24;
                //    a33 += rho_np2 / FieldRadius / FieldRadius * sum33;
                //    a34 -= rho_np2 / FieldRadius / FieldRadius * sum34;
                //    a44 += rho_np2 / FieldRadius / FieldRadius * sum44;
                //}
            }

            // Pines Equation 31 
            acc[0] = a1 + a4 * s;
            acc[1] = a2 + a4 * t;
            acc[2] = a3 + a4 * u;
            //if (fillgradient)
            //{
            //    // Pines Equation 37
            //    gradient[0, 0] = a11 + s * s * a44 + a4 / magr + 2 * s * a14;
            //    gradient[1, 1] = -a11 + t * t * a44 + a4 / magr + 2 * t * a24;
            //    gradient[2, 2] = a33 + u * u * a44 + a4 / magr + 2 * u * a34;
            //    gradient[0, 1] = gradient[1, 0] = a12 + s * t * a44 + s * a24 + t * a14;
            //    gradient[0, 2] = gradient[2, 0] = a13 + s * u * a44 + s * a34 + u * a14;
            //    gradient[1, 2] = gradient[2, 1] = a23 + t * u * a44 + u * a24 + t * a34;
            //}
        }  // FullGeopGMATPines


        // -----------------------------------------------------------------------------------------------\
        // Gottlieb approach for acceleration
        // gotpot in his nomenclature
        //
        // -----------------------------------------------------------------------------------------------\

    //    public void FullGeopGot
    //(
    //    AstroLib.gravityConst gravData,
    //    double[] recef,
    //    double[,,] normArr,
    //    int order,
    //    out double[,] legarrGot,
    //    out double[] G,
    //    out string straccum
    //)
    //    {
    //        straccum = "";

    //        legarrGot = new double[order + 2, order + 2];
    //        G = new double[3];

    //        //normArr = new double[order + 2, order + 2, 7];
    //        double[] zeta, eta, xi;
    //        zeta = new double[order + 1];
    //        eta = new double[order + 1];
    //        xi = new double[order + 1];
    //        double[] ctrigArr = new double[order + 1];
    //        double[] strigArr = new double[order + 1];

    //        double Ri, Xovr, Yovr, Zovr, sinlat, magr;
    //        double muor, muor2, Reor, Reorn;
    //        double Sumh, Sumgam, Sumj, Sumk, Lambda;
    //        double Sumh_N, Sumgam_N, Sumj_N, Sumk_N, Mxpnm;
    //        double BnmVal, pnm, snm, cnm;
    //        Int32 mm1, mp1, nm1, nm2, npmp1, Lim, Sum_Init;
    //        double cn;

    //        magr = MathTimeLibr.mag(recef);
    //        Ri = 1.0 / magr;
    //        Xovr = recef[0] * Ri;
    //        Yovr = recef[1] * Ri;
    //        Zovr = recef[2] * Ri;
    //        sinlat = Zovr;
    //        double coslat = Math.Cos(Math.Asin(sinlat));
    //        Reor = AstroLibr.gravConst.re * Ri;
    //        Reorn = Reor;
    //        muor = AstroLibr.gravConst.mu * Ri;
    //        muor2 = muor * Ri;

    //        // include two-body or not
    //        //if (Want_Central_force == true)
    //        //    Sum_Init = 1;
    //        //else
    //        // note, 1 makes the two body pretty close, except for the 1st component
    //        Sum_Init = 0;

    //        // initial values
    //        // ctrigArr[0] = 1.0;    
    //        ctrigArr[1] = Xovr;
    //        //  strigArr[0] = 0.0;
    //        strigArr[1] = Yovr;
    //        Sumh = 0.0;
    //        Sumj = 0.0;
    //        Sumk = 0.0;
    //        Sumgam = Sum_Init;

    //        // normArr(L, m, 0) xi Gottlieb eta
    //        // normArr(L, m, 1) eta Gottlieb zeta
    //        // normArr(L, m, 2) alpha Gottlieb alpha
    //        // normArr(L, m, 3) beta Gottlieb beta
    //        // normArr(L, m, 5) delta Gottlieb zn
    //        legarrGot[0, 0] = 1.0;
    //        legarrGot[0, 1] = 0.0;
    //        legarrGot[1, 0] = Math.Sqrt(3) * sinlat;
    //        legarrGot[1, 1] = Math.Sqrt(3); // * coslat;

    //        for (int n = 2; n <= order; n++)
    //        {
    //            // get the power for each n
    //            Reorn = Reorn * Reor;
    //            //pn = legPoly[n, 0];
    //            cn = gravData.cNor[n, 0];
    //            //sn = gravData.sNor[n, 0];

    //            nm1 = n - 1;
    //            nm2 = n - 2;

    //            // eq 3-17, eq 7-14  alpha(n) beta(n)
    //            legarrGot[n, 0] = sinlat * normArr[n, 0, 2] * legarrGot[nm1, 0] - normArr[n, 0, 3] * legarrGot[nm2, 0];
    //            // inner diagonal eq 7-16
    //            // n-1,n-2, 6, not 5, no nm1
    //            legarrGot[n, nm1] = normArr[n - 1, nm2, 6] * sinlat * legarrGot[n, n];
    //            //      legPoly[n, nm1] = normArr[n, nm1, 7] * sinlat * legPoly[n, n];
    //            // diagonal eq 7-8
    //            legarrGot[n, n] = normArr[n, n, 4] * coslat * legarrGot[nm1, nm1];

    //            Sumh_N = normArr[1, 0, 6] * legarrGot[n, 0] * cn;  // 0 by 2016 paper
    //            Sumgam_N = legarrGot[n, 0] * cn * (n + 1);  // double

    //            if (order > 0)
    //            {
    //                for (int m = 1; m <= nm2; m++)
    //                {
    //                    // eq 3-18, eq 7-12   xin(m) eta(m)
    //                    legarrGot[n, m] = normArr[n, m, 0] * sinlat * legarrGot[nm1, m] - normArr[n, m, 1] * legarrGot[nm2, m];
    //                }
    //                // got all the Legendre functions now

    //                Sumj_N = 0.0;
    //                Sumk_N = 0.0;
    //                ctrigArr[n] = ctrigArr[1] * ctrigArr[nm1] - strigArr[1] * strigArr[nm1]; // mm1????
    //                strigArr[n] = strigArr[1] * ctrigArr[nm1] + ctrigArr[1] * strigArr[nm1];

    //                if (n < order)
    //                    Lim = n;
    //                else
    //                    Lim = order;

    //                for (int m = 1; m <= Lim; m++)
    //                {
    //                    mm1 = m - 1;
    //                    mp1 = m + 1;
    //                    npmp1 = (n + mp1);  // double
    //                    pnm = legarrGot[n, m];
    //                    cnm = gravData.cNor[n, m];
    //                    snm = gravData.sNor[n, m];
    //                    //ctmm1 = gravData.cNor[n, mm1];
    //                    //stmm1 = gravData.sNor[n, mm1];

    //                    Mxpnm = m * pnm;  // double
    //                    BnmVal = cnm * ctrigArr[m] + snm * strigArr[m];
    //                    Sumh_N = Sumh_N + legarrGot[n, mp1] * BnmVal * normArr[n, m, 6];  // zn(m)
    //                    Sumgam_N = Sumgam_N + npmp1 * pnm * BnmVal;
    //                    Sumj_N = Sumj_N + Mxpnm * (cnm * ctrigArr[m] + snm * strigArr[m]);
    //                    Sumk_N = Sumk_N - Mxpnm * (cnm * strigArr[m] - snm * ctrigArr[m]);
    //                }   // for through m

    //                Sumj = Sumj + Reorn * Sumj_N;
    //                Sumk = Sumk + Reorn * Sumk_N;
    //            }  // if order > 0

    //            // ---- Sums bleow here have values when N m = 0
    //            Sumh = Sumh + Reorn * Sumh_N;
    //            Sumgam = Sumgam + Reorn * Sumgam_N;
    //        }  // loop

    //        Lambda = Sumgam + sinlat * Sumh;
    //        G[0] = -muor2 * (Lambda * Xovr - Sumj);
    //        G[1] = -muor2 * (Lambda * Yovr - Sumk);
    //        G[2] = -muor2 * (Lambda * Zovr - Sumh);

    //        // if (show == 'y')
    //        {
    //            straccum = straccum + "Gottlieb case nonspherical, no two-body ---------- " + "\n";
    //            straccum = straccum + "legarrGot 4 0   " + legarrGot[4, 0].ToString() + "  4 1   "
    //               + legarrGot[4, 1].ToString() + "  4 4   " + legarrGot[4, 4].ToString() + "\n";
    //            straccum = straccum + "legarrGot 5 0   " + legarrGot[5, 0].ToString() + "  5 1   "
    //                + ctrigArr[2].ToString() + "  Tan   " + strigArr[2].ToString() + "\n";
    //            straccum = straccum + "legarrGot" + order + " 0   " + legarrGot[order, 0].ToString() + "  " + order + " 1   "
    //                + legarrGot[order, 1].ToString() + " + order + " + legarrGot[order, order].ToString() + "\n";
    //            //straccum = straccum + "trigarr " + order + " Sin  " + trigArr[order, 0].ToString() + "  Cos   "
    //            //    + trigArr[order, 1].ToString() + "  Tan   " + trigArr[order, 3].ToString() + "\n";
    //            straccum = straccum + "apertGot ecef " + order + " " + order + " " + G[0].ToString() + "     "
    //                    + G[1].ToString() + "     " + G[2].ToString() + "\n";
    //        }

    //    }  // FullGeopGot;



        public void testproporbit()
        {
            double[] fArgs = new double[14];
            double[] reci = new double[3];
            double[] veci = new double[3];
            double[] aeci = new double[3];
            double[] aeci2 = new double[3];
            double[] recef = new double[3];
            double[] vecef = new double[3];
            double[] aecef = new double[3];
            double[] rsecef = new double[3];
            double[] vsecef = new double[3];
            double[] rseci = new double[3];
            double[] vseci = new double[3];
            double psia, wa, epsa, chia;
            double meaneps, deltapsi, deltaeps, trueeps;
            double[] omegaearth = new double[3];
            double[] rpef = new double[3];
            double[] vpef = new double[3];
            double[] apef = new double[3];
            double[] crossr = new double[3];
            double[] tempvec1 = new double[3];
            double[,] tm = new double[3, 3];
            double[,] prec = new double[3, 3];
            double[,] nut = new double[3, 3];
            double[,] st = new double[3, 3];
            double[,] pm = new double[3, 3];
            double[,] precp = new double[3, 3];
            double[,] nutp = new double[3, 3];
            double[,] stp = new double[3, 3];
            double[,] pmp = new double[3, 3];
            double[,] temp = new double[3, 3];
            double[,] temp1 = new double[3, 3];
            double[,] transeci2ecef = new double[3, 3];
            double[,] transecef2eci = new double[3, 3];
            double[,] convArr = new double[152, 152];

            double[] adrag = new double[3];
            double[] vrel = new double[3];

            double[] rsun = new double[3];
            double[] rsatsun = new double[3];
            double[] rmoon = new double[3];
            double[] rsat3 = new double[3];
            double[] rearth3 = new double[3];
            double[] a3body = new double[3];
            double[] athirdbody = new double[3];
            double[] athirdbody1 = new double[3];
            double[] athirdbody2 = new double[3];
            double[] aPertG = new double[3];
            double[] aPertGot = new double[3];
            double[] aPertM = new double[3];
            double[] aPertM1 = new double[3];

            double[] asrp = new double[3];

            double ttt, ttdb, xp, yp, lod, jdut1, ddpsi, ddeps, ddx, ddy, dut1, jdtt, jdftt;
            int year, mon, day, hr, minute, dat, eqeterms;
            double second, jdF, jdutc, jdFutc, jdtdb, jdFtdb;

            double hellp, latgd, lon;
            double cd, cr, area, mass, q, tmptdb;

            double latgc;
            Int32 order;
            double[,] LegArrMU;  // montenbruck
            double[,] LegArrMN;
            double[,] LegArrGU;  // gtds
            double[,] LegArrGN;
            double[,] LegGottN;
            double[,] LegArrEx;  // exact 
            // 152 is arbitrary
            Int32 orderSize = 500;
            double[,,] normArr = new double[orderSize, orderSize, 7];

            AstroLib.gravityConst gravData;

            double rad = 180.0 / Math.PI;              // deg to rad
            double conv = Math.PI / (180.0 * 3600.0);  // " to rad

            StringBuilder strbuildall = new StringBuilder();
            //  strbuild.Clear();
            StringBuilder strbuildplot = new StringBuilder();
            strbuildplot.Clear();

            // ------------------ BOOK EXAMPLE ----------------------------------------------
            //STK answers
            //18 Feb 2020 15:08:47.238    -605.79079600    -5870.23042200    3493.05191600    -1.56825100    -3.70234800    -6.47948500    0.000748395693    0.007252325791    -0.004327586380
            //j2x0        0.000748395693    0.007252325791    -0.004327586380
            //j21x21      0.000748452078    0.007252229967    -0.004327507399
            //            0.000000056385   -0.000000095824    -0.000000078981
            //drag
            //            0.000748730855    0.007255346798    -0.004317258342
            //3body
            //            0.000748731063    0.007255346392    -0.004317259450
            //srp
            //            0.000748730594    0.007255346775    -0.004317258523
            // ------------------------------- initial state -------------------------------
            reci = new double[] { -605.79079600, -5870.23042200, 3493.05191600 };
            veci = new double[] { -1.568251000, -3.702348000, -6.479485000 };
            // print out initial conditions
            strbuildall.AppendLine("reci  " + reci[0].ToString(fmt).PadLeft(4) + " " + reci[1].ToString(fmt).PadLeft(4) + " " + reci[2].ToString(fmt).PadLeft(4) + " " +
                                "v  " + veci[0].ToString(fmt).PadLeft(4) + " " + veci[1].ToString(fmt).PadLeft(4) + " " + veci[2].ToString(fmt).PadLeft(4));
            cd = 2.2;
            cr = 1.2;
            area = 40.0;     // m^2 
            mass = 1000.0;   // kg

            // ------------------------------- establish time parameters -------------------------------
            // read in FK5 parameters
            EOPSPWLib.iau80Class iau80arr;
            string fileLoc;
            fileLoc = @"D:\Codes\LIBRARY\DataLib\nut80.dat";
            EOPSPWLibr.iau80in(fileLoc, out iau80arr);
            // read in ICRF parameters
            fileLoc = @"D:\Codes\LIBRARY\DataLib\";
            EOPSPWLib.iau06Class iau06arr;
            EOPSPWLibr.iau06in(fileLoc, out iau06arr);
            AstroLib.xysdataClass[] xysarr = AstroLibr.xysarr;
            fileLoc = @"D:\Codes\LIBRARY\DataLib\";
            AstroLibr.readXYS(ref xysarr, fileLoc, "xysdata.dat");

            // example date/time
            // stkeducationalfiles ex 8-5
            year = 2020;
            mon = 2;
            day = 18;
            hr = 15;
            minute = 8;
            second = 47.23847;
            MathTimeLibr.jday(year, mon, day, hr, minute, second, out jdutc, out jdFutc);  // utc

            // these values are not consistent wih the EOP files - are they interpolated already????
            // no. change to use actual at 0000 values, rerun problem. 
            //2020 02 18 58897  0.030655  0.336009 -0.1990725  0.0000639 -0.108041 -0.007459  0.000315 -0.000055  37
            //2020 02 19 58898  0.030313  0.337617 - 0.1991016  0.0000235 - 0.107939 - 0.007476  0.000324 - 0.000036  37

            xp = 0.030655 * conv;
            yp = 0.336009 * conv;
            lod = 0.0000639;
            ddpsi = -0.108041 * conv;  // " to rad
            ddeps = -0.007459 * conv;
            ddx = 0.000315 * conv;     // " to rad
            ddy = 0.000055 * conv;
            dut1 = -0.1990725;   // second
            dat = 37;            // second
            AstroLib.EOpt opt = AstroLib.EOpt.e80;
            eqeterms = 2;
            jdtt = jdutc;
            jdftt = jdFutc + (dat + 32.184) / 86400.0;

            // method to do calculations in
            char normalized = 'y';
            strbuildall.AppendLine("normalized = " + normalized.ToString());
            strbuildall.AppendLine(year.ToString("0000") + " " + mon.ToString("00") + " " + day.ToString("00") + " " + hr.ToString("00") + ":" +
                minute.ToString("00") + ":" + second.ToString());
            strbuildall.AppendLine("dat " + dat.ToString() + " lod " + lod.ToString());
            strbuildall.AppendLine("jdutc " + (jdutc + jdFutc).ToString());
            strbuildall.AppendLine("xp yp " + (xp / conv).ToString() + " " + (yp / conv).ToString() + " arcsec");
            strbuildall.AppendLine("dpsi deps " + (ddpsi / conv).ToString() + " " + (ddeps / conv).ToString() + " arcsec");
            strbuildall.AppendLine("dx dy " + (ddx / conv).ToString() + " " + (ddy / conv).ToString() + " arcsec \n");

            jdut1 = jdutc + jdFutc + dut1 / 86400.0;
            strbuildall.AppendLine("jdut1 " + jdut1.ToString());

            // watch if getting tdb that j2000 is also tdb
            ttt = (jdutc + jdFutc + (dat + 32.184) / 86400.0 - 2451545.0) / 36525.0;
            strbuildall.AppendLine("jdttt " + (jdutc + jdFutc + (dat + 32.184) / 86400.0).ToString() + " ttt " + ttt.ToString() + "\n");

            AstroLibr.fundarg(ttt, opt, out fArgs);

            tmptdb = (dat + 32.184 + 0.001657 * Math.Sin(628.3076 * ttt + 6.2401)
                 + 0.000022 * Math.Sin(575.3385 * ttt + 4.2970)
                 + 0.000014 * Math.Sin(1256.6152 * ttt + 6.1969)
                 + 0.000005 * Math.Sin(606.9777 * ttt + 4.0212)
                 + 0.000005 * Math.Sin(52.9691 * ttt + 0.4444)
                 + 0.000002 * Math.Sin(21.3299 * ttt + 5.5431)
                 + 0.000010 * ttt * Math.Sin(628.3076 * ttt + 4.2490)) / 86400.0;  // USNO circ(14)
            MathTimeLibr.jday(year, mon, day, hr, minute, second + tmptdb, out jdtdb, out jdFtdb);
            ttdb = (jdtdb + jdFtdb - 2451545.0) / 36525.0;
            strbuildall.AppendLine("jdttb " + (jdtdb + jdFtdb).ToString() + " ttdb " + ttdb.ToString());

            // get reduction matrices
            deltapsi = 0.0;
            meaneps = 0.0;

            omegaearth[0] = 0.0;
            omegaearth[1] = 0.0;
            omegaearth[2] = AstroLibr.gravConst.earthrot * (1.0 - lod / 86400.0);

            prec = AstroLibr.precess(ttt, opt, out psia, out wa, out epsa, out chia);
            nut = AstroLibr.nutation(ttt, ddpsi, ddeps, iau80arr, fArgs, out deltapsi, out deltaeps, out trueeps, out meaneps);
            st = AstroLibr.sidereal(jdut1, deltapsi, meaneps, fArgs, lod, eqeterms, opt);
            pm = AstroLibr.polarm(xp, yp, ttt, opt);

            //// ---- perform transformations eci to ecef
            pmp = MathTimeLibr.mattrans(pm, 3);
            stp = MathTimeLibr.mattrans(st, 3);
            nutp = MathTimeLibr.mattrans(nut, 3);
            precp = MathTimeLibr.mattrans(prec, 3);
            temp = MathTimeLibr.matmult(pmp, stp, 3, 3, 3);
            temp1 = MathTimeLibr.matmult(temp, nutp, 3, 3, 3);
            transeci2ecef = MathTimeLibr.matmult(temp1, precp, 3, 3, 3);
            recef = MathTimeLibr.matvecmult(transeci2ecef, reci, 3);
            strbuildall.AppendLine("recef  " + recef[0].ToString(fmt).PadLeft(4) + " " + recef[1].ToString(fmt).PadLeft(4) + " " + recef[2].ToString(fmt).PadLeft(4) + " " +
                                "v  " + vecef[0].ToString(fmt).PadLeft(4) + " " + vecef[1].ToString(fmt).PadLeft(4) + " " + vecef[2].ToString(fmt).PadLeft(4));

            //----perform transformations ecef to eci
            // note the rotations occur only for velocity so the full transformation is fine here
            transecef2eci = MathTimeLibr.mattrans(transeci2ecef, 3);
            reci = MathTimeLibr.matvecmult(transecef2eci, recef, 3);
            strbuildall.AppendLine("reci  " + reci[0].ToString(fmt).PadLeft(4) + " " + reci[1].ToString(fmt).PadLeft(4) + " " + reci[2].ToString(fmt).PadLeft(4) + " " +
                                "v  " + veci[0].ToString(fmt).PadLeft(4) + " " + veci[1].ToString(fmt).PadLeft(4) + " " + veci[2].ToString(fmt).PadLeft(4));

            AstroLibr.eci_ecef(ref reci, ref veci, MathTimeLib.Edirection.eto, ref recef, ref vecef,
                  iau80arr, jdtt, jdftt, jdut1, lod, xp, yp, ddpsi, ddeps);

            // print out initial conditions
            strbuildall.AppendLine("reci  " + reci[0].ToString(fmt).PadLeft(4) + " " + reci[1].ToString(fmt).PadLeft(4) + " " + reci[2].ToString(fmt).PadLeft(4) + " " +
                                "v  " + veci[0].ToString(fmt).PadLeft(4) + " " + veci[1].ToString(fmt).PadLeft(4) + " " + veci[2].ToString(fmt).PadLeft(4));
            strbuildall.AppendLine("recef  " + recef[0].ToString(fmt).PadLeft(4) + " " + recef[1].ToString(fmt).PadLeft(4) + " " + recef[2].ToString(fmt).PadLeft(4) + " " +
                                "v  " + vecef[0].ToString(fmt).PadLeft(4) + " " + vecef[1].ToString(fmt).PadLeft(4) + " " + vecef[2].ToString(fmt).PadLeft(4));

            AstroLibr.ecef2ll(recef, out latgc, out latgd, out lon, out hellp);
            //// or
            //latgc = 52.0 / rad;  // 52 and 34
            //lon = 5.0 / rad;
            //alt = 6880.0;
            //AstroLibr.site(latgd, lon, alt, out recef, out vecef);
            //strbuildall.AppendLine("new site loc");
            strbuildall.AppendLine("satellite lat " + latgc * rad + " lon " + lon * rad + "alt " + hellp + " km");
            //strbuildall.AppendLine("recef  " + recef[0].ToString(fmt).PadLeft(4) + " " + recef[1].ToString(fmt).PadLeft(4) + " " + recef[2].ToString(fmt).PadLeft(4) + " " +
            //                    "v  " + vecef[0].ToString(fmt).PadLeft(4) + " " + vecef[1].ToString(fmt).PadLeft(4) + " " + vecef[2].ToString(fmt).PadLeft(4));

            // ----------------------------------------------------------------------------------------------
            // --------------------------------------- GRAVITY FIELD ----------------------------------------
            this.opsStatus.Text = "Status: Reading gravity field EGM-08 test";
            Refresh();
            int startKtr;
            // get past text in each file
            //if (fname.Contains("GEM"))    // GEM10bunnorm36.grv, GEMT3norm50.grv
            //startKtr = 17;
            //if (fname.Contains("EGM-96")) // EGM-96norm70.grv
            //startKtr = 73;
            //if (fname.Contains("EGM-08")) // EGM-08norm100.grv
            // fully normalized, 4415, .1363, order 100
            //startKtr = 83;  // or 21 for the larger file... which has gfc in the first col too
            //string fname = "D:/Dataorig/Gravity/EGM-08norm100.grv";  // skip line 83
            // fully normalized, 4415, .1363, order 2190
            startKtr = 0;  // or 21 for the larger file... which has gfc in the first col too
            string fname = "D:/Dataorig/Gravity/EGM2008_to2190_TideFree.txt";
            // fully normalized, 4415, .1363, order 360
            //string fname = "D:/Dataorig/Gravity/GGM03C-Data.txt";  

            normalized = 'y';  // if file has normalized coefficients

            AstroLibr.readGravityField(fname, normalized, startKtr, out order, out gravData);
            strbuildall.AppendLine("\nread in gravity field " + fname + " " + order.ToString() + " --------------- ");
            strbuildall.AppendLine("\ncoefficients --------------- ");
            strbuildall.AppendLine("c  2  0  " + gravData.c[2, 0].ToString() + " s " + gravData.s[2, 0].ToString());
            strbuildall.AppendLine("c  4  0  " + gravData.c[4, 0].ToString() + " s " + gravData.s[4, 0].ToString());
            strbuildall.AppendLine("c  4  4  " + gravData.c[4, 4].ToString() + " s " + gravData.s[4, 4].ToString());
            strbuildall.AppendLine("c 21  1 " + gravData.c[21, 1].ToString() + " s " + gravData.s[21, 1].ToString());
            strbuildall.AppendLine("\nnormalized coefficients --------------- ");
            strbuildall.AppendLine("c  2  0  " + gravData.cNor[2, 0].ToString() + " s " + gravData.sNor[2, 0].ToString());
            strbuildall.AppendLine("c  4  0  " + gravData.cNor[4, 0].ToString() + " s " + gravData.sNor[4, 0].ToString());
            strbuildall.AppendLine("c  4  4  " + gravData.cNor[4, 4].ToString() + " s " + gravData.sNor[4, 4].ToString());
            strbuildall.AppendLine("c 21  1 " + gravData.cNor[21, 1].ToString() + " s " + gravData.sNor[21, 1].ToString());
            strbuildall.AppendLine("c 500  1 " + gravData.cNor[500, 1].ToString() + " s " + gravData.sNor[500, 1].ToString());

            // --------------------------------------------------------------------------------------------------
            // calculate legendre polynomials
            //this.opsStatus.Text = "Status: Calculate Legendre polynomial recursions Unn and Nor ";
            //Refresh();

            strbuildall.AppendLine("\nCalculate Legendre polynomial recursions Unn and Nor  --------------- ");
            order = 21;

            // GTDS version
            // does with  unnormalized elements, then normalized from there. But unnormalized only go to about 170
            AstroLibr.LegPolyGTDS(latgc, order, normalized, out LegArrGU, out LegArrGN);

            // Gottlieb version does large versions  ------------------------------------------------
            double[] norm1 = new double[order + 2];
            double[] norm2 = new double[order + 2];
            double[] norm10 = new double[order + 2];
            double[] norm11 = new double[order + 2];
            double[] normn10 = new double[order + 2];
            double[,] normn1 = new double[order + 2, order + 2];
            double[,] norm1m = new double[order + 2, order + 2];
            double[,] norm2m = new double[order + 2, order + 2];
            AstroLibr.LegPolyGottN(order, out norm1, out norm2, out norm11, out normn10, out norm1m, out norm2m, out normn1);
            double[] aPertGt = new double[3];
            AstroLibr.FullGeopGott(recef, order, gravData, norm1, norm2, norm11, normn10,
                  norm1m, norm2m, normn1, out LegGottN, out aPertGot);

            // Montenbruck version
            // does with  unnormalized elements, then normalized from there. But unnormalized only go to about 170
            AstroLibr.LegPolyMont(latgc, order, normalized, out LegArrMU, out LegArrMN);

            // Geodyn version
            //AstroLibr.geodynlegp(latgc, degree, order, out LegArrOU, out LegArrON);

            // Exact values
            LegPolyEx(latgc, order, out LegArrEx);

            // Fukushima approach do as 1-d arrays for now
            //AstroLibr.LegPolyFN(latgc, order, normArr, out LegArrFN);
            //double[] pmm = new double[8];
            //double[] psm = new double[8];
            //Int32[] ipsm = new Int32[8];
            // get the values in X-numbers
            //alfsx(Math.Cos(latgc), 6, normArr, out psm, out ipsm);
            //alfmx(Math.Sin(latgc), 3, 6, normArr, psm[3], ipsm[3], out pmm);

            //strbuildall.AppendLine("latgc = " + latgc * rad);
            //double[,] P = AstroLibr.legendre(degree, order, latgc);
            //double P = AstroLibr.legendreS(degree, order, latgc);
            //AstroLibr.legendreMa(degree, order, latgc);

            string errstr = " ";
            //double dr1, dr2, dr3, sumdr1, sumdr2, sumdr3;
            //sumdr1 = 0.0;
            //sumdr2 = 0.0;
            //sumdr3 = 0.0;
            strbuildall.AppendLine("\nwrite out normalized Legendre polynomials --------------- ");

            // order xxxxxxxxxxxxxxxxxx
            for (int L = 0; L < 21; L++)
            {
                string tempstr1 = "MN  ";  // montenbruck
                string tempstr2 = "GN  ";  // gtds
                //string tempstr3 = "MU  ";
                //string tempstr3a = "LU  ";  // exact
                //string tempstr4 = "OU  ";  // geodyn\
                string tempstr5 = "GtN ";  // gottlieb\
                //string tempstr6 = "GtU ";  // gottlieb
                //string tempstr7 = "FN  ";  // Fukushima, test ones
                int stopL = L;
                for (int m = 0; m <= stopL; m++)
                {
                    tempstr1 = tempstr1 + " " + L.ToString() + "  " + m.ToString() + "   " + LegArrMN[L, m].ToString();
                    tempstr2 = tempstr2 + " " + L.ToString() + "  " + m.ToString() + "   " + LegArrGN[L, m].ToString();
                    tempstr5 = tempstr5 + " " + L.ToString() + "  " + m.ToString() + "   " + LegGottN[L, m].ToString();
                    //tempstr7 = tempstr7 + " " + L.ToString() + "  " + m.ToString() + "   " + LegArrF[L, m].ToString();
                    //tempstr3 = tempstr3 + " " + L.ToString() + "  " + m.ToString() + "   " + LegArrMU[L, m].ToString();
                    //tempstr3a = tempstr3a + " " + L.ToString() + "  " + m.ToString() + "   " + LegArrEx[L, m].ToString();
                    //   tempstr6 = tempstr6 + " " + L.ToString() + "  " + m.ToString() + "   " + LegArrGotU[L, m].ToString();
                    //tempstr4 = tempstr4 + " " + L.ToString() + "  " + m.ToString() + "   "; // + LegArrOU[L + 1, m + 1].ToString();
                    // check error values
                    //dr1 = 100.0 * (LegArrF[L, m] - LegGottN[L, m]) / LegArrF[L, m];
                    //dr2 = 100.0 * (LegArrF[L, m] - LegArrGN[L, m]) / LegArrF[L, m];
                    //dr3 = 100.0 * (LegArrF[L, m] - LegArrMN[L, m]) / LegArrF[L, m];
                    //sumdr1 = sumdr1 + dr1;
                    //sumdr2 = sumdr2 + dr2;
                    //sumdr3 = sumdr3 + dr3;
                    //errstr = errstr + "\n" + L.ToString() + "  " + m.ToString() //+ "   " + dr1.ToString()
                    //    + " " + dr2.ToString() + " " + dr3.ToString();
                }
                // normalized ones
                strbuildall.AppendLine(tempstr2);
                strbuildall.AppendLine(tempstr1);
                strbuildall.AppendLine(tempstr5 + "\n");
                // strbuildall.AppendLine(tempstr7 + "\n");
                // unnormalized ones
                //strbuildall.AppendLine(tempstr3);
                //strbuildall.AppendLine(tempstr3a);
                //strbuildall.AppendLine(tempstr6);
                //strbuildall.AppendLine(tempstr4 + "\n");
            }
            strbuildplot.AppendLine(errstr);

            // ---------------------------------------- now accelerations ---------------------------------
            order = 21;
            strbuildall.AppendLine("\naccelerations --------------- ");
            string straccum = "";
            this.opsStatus.Text = "Status: Calculate Accelerations --------------- ";
            Refresh();

            if (order < 100)
            {
                // GTDS acceleration for non-spherical portion
                AstroLibr.FullGeopG(recef, order, normalized, convArr, normArr, gravData, out aPertG, 'y', out straccum);
                strbuildall.AppendLine(straccum);
                aeci = MathTimeLibr.matvecmult(transecef2eci, aPertG, 3);
                straccum = straccum + "apertG eci  " + order + " " + order + " " + aeci[0].ToString() + "     "
                        + aeci[1].ToString() + "     " + aeci[2].ToString() + "\n";
                strbuildall.AppendLine(straccum);

                AstroLibr.FullGeopG(recef, order, normalized, convArr, normArr, gravData, out aPertG, 'n', out straccum);
                strbuildall.AppendLine(straccum);

                // Montenbruck acceleration
                AstroLibr.FullGeopM(recef, order, normalized, convArr, gravData, out aPertM, 'y', out straccum);
                strbuildall.AppendLine(straccum);
                aeci = MathTimeLibr.matvecmult(transecef2eci, aPertM, 3);
                straccum = straccum + "apertM eci  " + order + " " + order + " " + aeci[0].ToString() + "     "
                        + aeci[1].ToString() + "     " + aeci[2].ToString() + "\n";
                strbuildall.AppendLine(straccum);

                AstroLibr.FullGeopM(recef, order, normalized, convArr, gravData, out aPertM, 'n', out straccum);
                strbuildall.AppendLine(straccum);

                // Montenbruck code acceleration
                AstroLibr.FullGeopMC(recef, order, normalized, convArr, gravData, out aPertM1, 'y', out straccum);
                strbuildall.AppendLine(straccum);
                aeci = MathTimeLibr.matvecmult(transecef2eci, aPertM1, 3);
                straccum = straccum + "apertM1 eci " + order + " " + order + " " + aeci[0].ToString() + "     "
                        + aeci[1].ToString() + "     " + aeci[2].ToString() + "\n";
                strbuildall.AppendLine(straccum);

                AstroLibr.FullGeopMC(recef, order, normalized, convArr, gravData, out aPertM1, 'n', out straccum);
                strbuildall.AppendLine(straccum);
            }

            // Pines acceleration
            //strbuildall.AppendLine("Pines acceleration ");
            //double[,] LegPineN = new double[order+2, order+2];
            //double[] aPertPi = new double[3];
            //AstroLibr.FullGeopPines(latgc, order, gravData, recef, out LegPineN, out aPertPi);
            //strbuildall.AppendLine(straccum);
            //strbuildall.AppendLine("0 " + LegPineN[0, 0].ToString());
            //strbuildall.AppendLine("1 " + LegPineN[1, 0].ToString() + " " + LegPineN[1, 1].ToString());
            //strbuildall.AppendLine("2 " + LegPineN[2, 0].ToString() + " " + LegPineN[2, 1].ToString()
            //    + " " + LegPineN[2, 2].ToString());
            //strbuildall.AppendLine("3 " + LegPineN[3, 0].ToString() + " " + LegPineN[3, 1].ToString()
            //    + " " + LegPineN[3, 2].ToString() + " " + LegPineN[3, 3].ToString());
            //strbuildall.AppendLine("acceleration " + order + " " + aPertPi[0] + " " + aPertPi[1] + " "
            //    + aPertPi[2]);

            // Lear acceleration
            //strbuildall.AppendLine("Lear acceleration ");
            ////double[] G = new double[3];
            //double[] aPertLr = new double[3];
            //AstroLibr.FullGeopLear(latgc, order, gravData, recef, out LegLearN, out aPertLr);
            //strbuildall.AppendLine(straccum);
            //strbuildall.AppendLine("0 " + LegLearN[0, 0].ToString());
            //strbuildall.AppendLine("1 " + LegLearN[1, 0].ToString() + " " + LegLearN[1, 1].ToString());
            //strbuildall.AppendLine("2 " + LegLearN[2, 0].ToString() + " " + LegLearN[2, 1].ToString()
            //    + " " + LegLearN[2, 2].ToString());
            //strbuildall.AppendLine("3 " + LegLearN[3, 0].ToString() + " " + LegLearN[3, 1].ToString()
            //    + " " + LegLearN[3, 2].ToString() + " " + LegLearN[3, 3].ToString());
            //strbuildall.AppendLine("acceleration " + order + " " + aPertLr[0] + " " + aPertLr[1] + " "
            //    + aPertLr[2]);

            // Gottlieb acceleration
            strbuildall.AppendLine("Gottlieb acceleration ");
            // this part can be done one time
            AstroLibr.LegPolyGottN(order, out norm1, out norm2, out norm11, out normn10, out norm1m,
                out norm2m, out normn1);
            AstroLibr.FullGeopGott(recef, order, gravData, norm1, norm2, norm11, normn10,
                  norm1m, norm2m, normn1, out LegGottN, out aPertGot);
            strbuildall.AppendLine(straccum);
            strbuildall.AppendLine("0 " + LegGottN[0, 0].ToString());
            strbuildall.AppendLine("1 " + LegGottN[1, 0].ToString() + " " + LegGottN[1, 1].ToString());
            strbuildall.AppendLine("2 " + LegGottN[2, 0].ToString() + " " + LegGottN[2, 1].ToString()
                + " " + LegGottN[2, 2].ToString());
            strbuildall.AppendLine("3 " + LegGottN[3, 0].ToString() + " " + LegGottN[3, 1].ToString()
                + " " + LegGottN[3, 2].ToString() + " " + LegGottN[3, 3].ToString());
            strbuildall.AppendLine("acceleration " + order + " " + aPertGot[0] + " " + aPertGot[1]
                + " " + aPertGot[2]);

            // Fukushima acceleration
            //strbuildall.AppendLine("Fukushima acceleration ");
            //LegPolyFF(recef, latgc, order, 'y', normArr, gravData, out LegArrF);
            //double[,] a = new double[360, 360];
            //double[,] b = new double[360, 360];
            //xfsh2f(80, gravData, out a, out b);
            //strbuildall.AppendLine(a[2, 0].ToString());
            //strbuildall.AppendLine("a  2  0  " + a[2, 0].ToString() + " b " + b[2, 0].ToString());
            //strbuildall.AppendLine("a  2  1  " + a[2, 1].ToString() + " b " + b[2, 1].ToString());
            //strbuildall.AppendLine("a  4  0  " + a[4, 0].ToString() + " b " + b[4, 0].ToString());
            //strbuildall.AppendLine("a  4  1  " + a[4, 0].ToString() + " b " + b[4, 1].ToString());
            //strbuildall.AppendLine("a  4  4  " + a[4, 4].ToString() + " b " + b[4, 4].ToString());
            //strbuildall.AppendLine("a 10 10  " + a[10, 0].ToString() + " b " + b[10, 0].ToString());
            //strbuildall.AppendLine("a 21  1 " + a[21, 1].ToString() + " b " + b[21, 1].ToString());

            // Pines approach
            //strbuildall.AppendLine("Pines acceleration ");
            //FullGeopPines(jdutc, recef, latgc, order, order, gravData, out aeci);
            //strbuildall.AppendLine("apertP    4 4   " + aeci[0].ToString() + "     " + aeci[1].ToString() + "     "
            //+ aeci[2].ToString());

            strbuildall.AppendLine(straccum);
            strbuildall.AppendLine("\ngravity field " + fname + " " + order.ToString() + " --------------- ");
            strbuildall.AppendLine(" summary accelerations ----------------------------------------------- ");
            strbuildall.AppendLine("apertG bf  " + order + " " + order + " " + aPertG[0].ToString() + "     "
                + aPertG[1].ToString() + "     " + aPertG[2].ToString());
            strbuildall.AppendLine("apertM bf  " + order + " " + order + " " + aPertM[0].ToString() + "     "
                + aPertM[1].ToString() + "     " + aPertM[2].ToString());
            strbuildall.AppendLine("apertMC bf " + order + " " + order + " " + aPertM1[0].ToString() + "     "
                + aPertM1[1].ToString() + "     " + aPertM1[2].ToString());
            strbuildall.AppendLine("apertGt bf " + order + " " + order + " " + aPertGot[0].ToString() + "     "
                + aPertGot[1].ToString() + "     " + aPertGot[2].ToString());

            aPertG = MathTimeLibr.matvecmult(transecef2eci, aPertG, 3);
            aPertM = MathTimeLibr.matvecmult(transecef2eci, aPertM, 3);
            aPertM1 = MathTimeLibr.matvecmult(transecef2eci, aPertM1, 3);
            aPertGt = MathTimeLibr.matvecmult(transecef2eci, aPertGot, 3);
            strbuildall.AppendLine("apertG  eci " + order + " " + order + " " + aPertG[0].ToString() + "     "
                + aPertG[1].ToString() + "     " + aPertG[2].ToString());
            strbuildall.AppendLine("apertM  eci " + order + " " + order + " " + aPertM[0].ToString() + "     "
                + aPertM[1].ToString() + "     " + aPertM[2].ToString());
            strbuildall.AppendLine("apertMC eci " + order + " " + order + " " + aPertM1[0].ToString() + "     "
                + aPertM1[1].ToString() + "     " + aPertM1[2].ToString());
            strbuildall.AppendLine("apertGt eci " + order + " " + order + " " + aPertGt[0].ToString() + "     "
                + aPertGt[1].ToString() + "     " + aPertGt[2].ToString());

            strbuildall.AppendLine("STK ans 4x4         -0.0000003723020	-0.0000031362090   	-0.0000102647170\n");  // no 2-body


            // -------------------------- add in two body term since full geop is only disturbing part
            AstroLibr.eci_ecef(ref reci, ref veci, MathTimeLib.Edirection.efrom, ref recef, ref vecef,
                  iau80arr, jdtt, jdftt, jdut1, lod, xp, yp, ddpsi, ddeps);
            aeci2[0] = -AstroLibr.gravConst.mu * reci[0] / (Math.Pow(MathTimeLibr.mag(reci), 3));
            aeci2[1] = -AstroLibr.gravConst.mu * reci[1] / (Math.Pow(MathTimeLibr.mag(reci), 3));
            aeci2[2] = -AstroLibr.gravConst.mu * reci[2] / (Math.Pow(MathTimeLibr.mag(reci), 3));
            strbuildall.AppendLine("a2body      " + aeci2[0].ToString() + "     " + aeci2[1].ToString() + "     "
                + aeci2[2].ToString());

            aPertG[0] = aPertG[0] + aeci2[0];
            aPertG[1] = aPertG[1] + aeci2[1];
            aPertG[2] = aPertG[2] + aeci2[2];

            double[] temm = new double[3];
            temm[0] = aPertM1[0];
            temm[1] = aPertM1[1];
            temm[2] = aPertM1[2];

            aPertM[0] = aPertM[0] + aeci2[0];
            aPertM[1] = aPertM[1] + aeci2[1];
            aPertM[2] = aPertM[2] + aeci2[2];

            aPertM1[0] = aPertM1[0] + aeci2[0];
            aPertM1[1] = aPertM1[1] + aeci2[1];
            aPertM1[2] = aPertM1[2] + aeci2[2];

            strbuildall.AppendLine(" now with two body included");
            strbuildall.AppendLine("apertG " + order + " " + order + " " + aPertG[0].ToString() + "     " + aPertG[1].ToString() + "     " + aPertG[2].ToString());
            strbuildall.AppendLine("apertM " + order + " " + order + " " + aPertM[0].ToString() + "     " + aPertM[1].ToString() + "     " + aPertM[2].ToString());
            strbuildall.AppendLine("apertMC " + order + " " + order + " " + aPertM1[0].ToString() + "     " + aPertM1[1].ToString() + "     " + aPertM1[2].ToString());
            strbuildall.AppendLine("STK ans 4x4 w2   0.0007483593980          0.0072522125910         -0.0043275195170\n");  // no 2-body
            //                    4x4 j2000   0.00074835849281         0.00725221243453        -0.00432751993509
            //                    4x4 icrf    0.00074835939828         0.00725221259059        -0.00432751951698
            //                    all j2000   0.00074845403274         0.00725223127396        -0.00432750265312
            //                    all icrf    0.00074845493821         0.00725223143002        -0.00432750223499

            this.opsStatus.Text = "Status: Other perts";
            Refresh();

            strbuildall.AppendLine("------------------ find drag acceleration");
            double density = 1.5e-12;  // kg / m3
            double magv = MathTimeLibr.mag(vecef);
            vrel[0] = vecef[0]; // vecef normal is veci to tod, then - wxr
            vrel[1] = vecef[1];
            vrel[2] = vecef[2];
            strbuildall.AppendLine(" vrel " + vrel[0].ToString() + " " + vrel[1].ToString() + " " + vrel[2].ToString());
            //                 kg / m3        m2  /  kg     km / s  km / s
            adrag[0] = -0.5 * density * cd * area / mass * magv * magv * vrel[0] / magv * 1000.0;  // simplify vel, get units to km/s2
            adrag[1] = -0.5 * density * cd * area / mass * magv * magv * vrel[1] / magv * 1000.0;  // simplify vel, get units to km/s2
            adrag[2] = -0.5 * density * cd * area / mass * magv * magv * vrel[2] / magv * 1000.0;  // simplify vel, get units to km/s2

            strbuildall.AppendLine(" adrag ecef" + adrag[0].ToString() + " " + adrag[1].ToString() + " " + adrag[2].ToString());

            strbuildall.AppendLine(" agrav + drag ecef" + (temm[0] + adrag[0]).ToString() + " " +
                (temm[1] + adrag[1]).ToString() + " " + (temm[2] + adrag[2]).ToString());


            transecef2eci = MathTimeLibr.matmult(temp1, pm, 3, 3, 3);
            aeci = MathTimeLibr.matvecmult(transecef2eci, adrag, 3);
            strbuildall.AppendLine(" adrag eci " + aeci[0].ToString() + " " + aeci[1].ToString() + " " + aeci[2].ToString());
            strbuildall.AppendLine("ans drag JR spline      0.0000000001040	0.0000000002090	0.0000000003550\n");
            strbuildall.AppendLine("ans drag JR daily       0.0000000000840	0.0000000001720	0.0000000002900\n");
            strbuildall.AppendLine("ans drag MSIS daily     0.0000000000730	0.0000000001510	0.0000000002530\n");

            double[] temmm = new double[3];
            temmm[0] = temm[0] + adrag[0];
            temmm[1] = temm[1] + adrag[1];
            temmm[2] = temm[2] + adrag[2];
            aeci = MathTimeLibr.matvecmult(transecef2eci, temmm, 3);
            strbuildall.AppendLine(" agrav+drag eci " + aeci[0].ToString() + " " + aeci[1].ToString() + " " + aeci[2].ToString());

            strbuildall.AppendLine(" ------------------ find third body acceleration");
            AstroLib.jpldedataClass[] jpldearr = AstroLibr.jpldearr;
            double musun, mumoon, rsmag, rmmag;
            musun = 1.32712428e11;    // km3 / s2
            mumoon = 4902.799;        // km3 / s2

            // make sure to fix mfme inside findjpl if using 12 hr data!!!!!!!!!!!!!!!!!!!!!!
            AstroLibr.readjplde(ref jpldearr, @"D:\Codes\LIBRARY\DataLib\", "sunmooneph_430t12.txt");
            //AstroLibr.initjplde(ref jpldearr, @"D:\Codes\LIBRARY\DataLib\", "sunmooneph_430t.txt");

            // sun
            AstroLibr.findjpldeparam(jdtdb, jdFtdb, 's', jpldearr, out rsun, out rsmag, out rmoon, out rmmag);
            // stk value (chk that tdb is argument)
            double[] rsuns = new double[] { 126916355.384390, -69567131.339884, -30163629.424510 };
            // 126054577.0753 18th       46117
            // 126743428.8631 18th 1200  46118 if 12 hr
            // 127422565.3381 19th       46118 if daily
            // JPL ans  18 Feb 2020 15:08:47.23847 M       0.6306    
            double[] rmoonj = new double[] { 14462.2967, -357096.9762, -151599.3021 };
            // JPL ans  2020  2 18 15:08:47.23847 S       0.6306  
            double[] rsunj = new double[] { 126921698.4134, -69564121.8695, -30156263.9220 };

            MathTimeLibr.addvec(1.0, rsuns, -1.0, rsun, out tempvec1); // km
            strbuildall.AppendLine(" diff rsun stk-mine " + tempvec1[0].ToString() + " " + tempvec1[1].ToString() + " " +
                tempvec1[2].ToString() + " " + MathTimeLibr.mag(tempvec1).ToString());
            MathTimeLibr.addvec(1.0, rsunj, -1.0, rsun, out tempvec1);
            strbuildall.AppendLine(" diff rsun jpl-mine " + tempvec1[0].ToString() + " " + tempvec1[1].ToString() + " " +
                tempvec1[2].ToString() + " " + MathTimeLibr.mag(tempvec1).ToString());
            MathTimeLibr.addvec(1.0, rsuns, -1.0, rsunj, out tempvec1);
            strbuildall.AppendLine(" diff rsun stk-jpl  " + tempvec1[0].ToString() + " " + tempvec1[1].ToString() + " " +
                tempvec1[2].ToString() + " " + MathTimeLibr.mag(tempvec1).ToString());
            MathTimeLibr.addvec(1.0, rmoonj, -1.0, rmoon, out tempvec1);
            strbuildall.AppendLine(" diff rmoon jpl-mine " + tempvec1[0].ToString() + " " + tempvec1[1].ToString() + " " +
                tempvec1[2].ToString() + " " + MathTimeLibr.mag(tempvec1).ToString());
            strbuildall.AppendLine(" rsun  " + rsun[0].ToString() + " " + rsun[1].ToString() + " " + rsun[2].ToString());
            strbuildall.AppendLine(" rmoon " + rmoon[0].ToString() + " " + rmoon[1].ToString() + " " + rmoon[2].ToString());

            double mu3 = musun;
            rsat3[0] = rsun[0] - reci[0];
            rsat3[1] = rsun[1] - reci[1];
            rsat3[2] = rsun[2] - reci[2];
            double magrsat3 = MathTimeLibr.mag(rsat3);
            rearth3[0] = rsun[0];
            rearth3[1] = rsun[1];
            rearth3[2] = rsun[2];
            double magrearth3 = MathTimeLibr.mag(rearth3);
            athirdbody[0] = mu3 * (rsat3[0] / Math.Pow(magrsat3, 3) - rearth3[0] / Math.Pow(magrearth3, 3));
            athirdbody[1] = mu3 * (rsat3[1] / Math.Pow(magrsat3, 3) - rearth3[1] / Math.Pow(magrearth3, 3));
            athirdbody[2] = mu3 * (rsat3[2] / Math.Pow(magrsat3, 3) - rearth3[2] / Math.Pow(magrearth3, 3));
            strbuildall.AppendLine(" a3bodyS  eci " + athirdbody[0].ToString() + " " + athirdbody[1].ToString() + " " + athirdbody[2].ToString());
            athirdbody2[0] = -mu3 / Math.Pow(magrearth3, 3) * (rearth3[0] - 3.0 * rearth3[0] * (MathTimeLibr.dot(reci, rearth3) / Math.Pow(magrearth3, 2))
                - 7.5 * rearth3[0] * Math.Pow(((MathTimeLibr.dot(reci, rearth3) / Math.Pow(magrearth3, 2))), 2));
            athirdbody2[1] = -mu3 / Math.Pow(magrearth3, 3) * (rearth3[1] - 3.0 * rearth3[1] * (MathTimeLibr.dot(reci, rearth3) / Math.Pow(magrearth3, 2))
                - 7.5 * rearth3[1] * Math.Pow(((MathTimeLibr.dot(reci, rearth3) / Math.Pow(magrearth3, 2))), 2));
            athirdbody2[2] = -mu3 / Math.Pow(magrearth3, 3) * (rearth3[2] - 3.0 * rearth3[2] * (MathTimeLibr.dot(reci, rearth3) / Math.Pow(magrearth3, 2))
                - 7.5 * rearth3[2] * Math.Pow(((MathTimeLibr.dot(reci, rearth3) / Math.Pow(magrearth3, 2))), 2));
            strbuildall.AppendLine(" a3bodyS2 eci" + athirdbody2[0].ToString() + " " + athirdbody2[1].ToString() + " " + athirdbody2[2].ToString());
            q = (Math.Pow(MathTimeLibr.mag(reci), 2) + 2.0 * MathTimeLibr.dot(reci, rsat3)) *
                (Math.Pow(magrearth3, 2) + magrearth3 * magrsat3 + Math.Pow(magrsat3, 2)) /
                (Math.Pow(magrearth3, 3) * Math.Pow(magrsat3, 3) * (magrearth3 + magrsat3));
            athirdbody1[0] = mu3 * (rsat3[0] * q - reci[0] / Math.Pow(magrearth3, 3));
            athirdbody1[1] = mu3 * (rsat3[1] * q - reci[1] / Math.Pow(magrearth3, 3));
            athirdbody1[2] = mu3 * (rsat3[2] * q - reci[2] / Math.Pow(magrearth3, 3));
            strbuildall.AppendLine(" a3bodyS1 eci" + athirdbody1[0].ToString() + " " + athirdbody1[1].ToString() + " " + athirdbody1[2].ToString());
            strbuildall.AppendLine("ans sun        0.0000000001820	0.0000000001620	-0.0000000001800\n");
            a3body[0] = athirdbody1[0];
            a3body[1] = athirdbody1[1];
            a3body[2] = athirdbody1[2];

            // moon
            mu3 = mumoon;
            rsat3[0] = rmoon[0] - reci[0];
            rsat3[1] = rmoon[1] - reci[1];
            rsat3[2] = rmoon[2] - reci[2];
            magrsat3 = MathTimeLibr.mag(rsat3);
            rearth3[0] = rmoon[0];
            rearth3[1] = rmoon[1];
            rearth3[2] = rmoon[2];
            magrearth3 = MathTimeLibr.mag(rearth3);
            athirdbody[0] = mu3 * (rsat3[0] / Math.Pow(magrsat3, 3) - rearth3[0] / Math.Pow(magrearth3, 3));
            athirdbody[1] = mu3 * (rsat3[1] / Math.Pow(magrsat3, 3) - rearth3[1] / Math.Pow(magrearth3, 3));
            athirdbody[2] = mu3 * (rsat3[2] / Math.Pow(magrsat3, 3) - rearth3[2] / Math.Pow(magrearth3, 3));
            strbuildall.AppendLine(" a3bodyM  eci " + athirdbody[0].ToString() + " " + athirdbody[1].ToString() + " " + athirdbody[2].ToString());
            athirdbody2[0] = -mu3 / Math.Pow(magrearth3, 3) * (rearth3[0] - 3.0 * rearth3[0] * (MathTimeLibr.dot(reci, rearth3) / Math.Pow(magrearth3, 2))
                - 7.5 * rearth3[0] * Math.Pow(((MathTimeLibr.dot(reci, rearth3) / Math.Pow(magrearth3, 2))), 2));
            athirdbody2[1] = -mu3 / Math.Pow(magrearth3, 3) * (rearth3[1] - 3.0 * rearth3[1] * (MathTimeLibr.dot(reci, rearth3) / Math.Pow(magrearth3, 2))
                - 7.5 * rearth3[1] * Math.Pow(((MathTimeLibr.dot(reci, rearth3) / Math.Pow(magrearth3, 2))), 2));
            athirdbody2[2] = -mu3 / Math.Pow(magrearth3, 3) * (rearth3[2] - 3.0 * rearth3[2] * (MathTimeLibr.dot(reci, rearth3) / Math.Pow(magrearth3, 2))
                - 7.5 * rearth3[2] * Math.Pow(((MathTimeLibr.dot(reci, rearth3) / Math.Pow(magrearth3, 2))), 2));
            strbuildall.AppendLine(" a3bodyM2 eci" + athirdbody2[0].ToString() + " " + athirdbody2[1].ToString() + " " + athirdbody2[2].ToString());
            q = (Math.Pow(MathTimeLibr.mag(reci), 2) + 2.0 * MathTimeLibr.dot(reci, rsat3)) *
                (Math.Pow(magrearth3, 2) + magrearth3 * magrsat3 + Math.Pow(magrsat3, 2)) /
                (Math.Pow(magrearth3, 3) * Math.Pow(magrsat3, 3) * (magrearth3 + magrsat3));
            athirdbody1[0] = mu3 * (rsat3[0] * q - reci[0] / Math.Pow(magrearth3, 3));
            athirdbody1[1] = mu3 * (rsat3[1] * q - reci[1] / Math.Pow(magrearth3, 3));
            athirdbody1[2] = mu3 * (rsat3[2] * q - reci[2] / Math.Pow(magrearth3, 3));
            strbuildall.AppendLine(" a3bodyM1 eci" + athirdbody1[0].ToString() + " " + athirdbody1[1].ToString() + " " + athirdbody1[2].ToString());
            strbuildall.AppendLine("ans moon        0.0000000000860	-0.0000000004210	-0.0000000006980\n");
            // total acceleration
            a3body[0] = a3body[0] + athirdbody1[0];
            a3body[1] = a3body[1] + athirdbody1[1];
            a3body[2] = a3body[2] + athirdbody1[2];
            strbuildall.AppendLine(" a3body S/M eci" + a3body[0].ToString() + " " + a3body[1].ToString() + " " + a3body[2].ToString());
            strbuildall.AppendLine("ans sun/moon    0.0000000002730	-0.0000000002680	-0.0000000008800\n");


            strbuildall.AppendLine(" ------------------ find srp acceleration\n");
            double psrp = 4.56e-6;  // N/m2 = kgm/s2 / m2 = kg/ms2
            rsatsun[0] = rsun[0] - reci[0];
            rsatsun[1] = rsun[1] - reci[1];
            rsatsun[2] = rsun[2] - reci[2];
            double magrsatsun = MathTimeLibr.mag(rsatsun);
            //           kg/ms2      m2      kg      km            km
            asrp[0] = -(psrp * cr * area / mass * rsatsun[0] / magrsatsun) / 1000.0;  // result in km/s
            asrp[1] = -(psrp * cr * area / mass * rsatsun[1] / magrsatsun) / 1000.0;
            asrp[2] = -(psrp * cr * area / mass * rsatsun[2] / magrsatsun) / 1000.0;
            strbuildall.AppendLine(" asrp eci " + asrp[0].ToString() + " " + asrp[1].ToString() + " " + asrp[2].ToString());
            strbuildall.AppendLine("ans srp        -0.0000000001970	0.0000000001150	0.0000000000480\n");

            strbuildall.AppendLine(" ------------------ add perturbing accelerations\n");
            aecef[0] = adrag[0];  // plus gravity xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            aecef[1] = adrag[1];
            aecef[2] = adrag[2];
            strbuildall.AppendLine(" aecef " + aecef[0].ToString() + " " + aecef[1].ToString() + " " + aecef[2].ToString());

            // ---- move acceleration from earth fixed coordinates to eci
            // there are no cross products here as normal
            aeci = MathTimeLibr.matvecmult(transecef2eci, aecef, 3);
            strbuildall.AppendLine(" aeci " + aeci[0].ToString() + " " + aeci[1].ToString() + " " + aeci[2].ToString());

            // find two body component of eci acceleration
            aeci2[0] = -AstroLibr.gravConst.mu * reci[0] / (Math.Pow(MathTimeLibr.mag(reci), 3));
            aeci2[1] = -AstroLibr.gravConst.mu * reci[1] / (Math.Pow(MathTimeLibr.mag(reci), 3));
            aeci2[2] = -AstroLibr.gravConst.mu * reci[2] / (Math.Pow(MathTimeLibr.mag(reci), 3));
            strbuildall.AppendLine(" aeci2body " + aeci2[0].ToString() + " " + aeci2[1].ToString() + " " + aeci2[2].ToString());

            // total acceleration
            aeci[0] = aeci2[0] + a3body[0] + asrp[0] + aeci[0];
            aeci[1] = aeci2[1] + a3body[1] + asrp[1] + aeci[1];
            aeci[2] = aeci2[2] + a3body[2] + asrp[2] + aeci[2];
            strbuildall.AppendLine("total aeci " + aeci[0].ToString() + " " + aeci[1].ToString() + " " + aeci[2].ToString());




            // ------------------------------------------- timing comparisons
            //strbuildall.AppendLine("\n ===================================== Timing Comparisons =====================================");
            //// timing of routines
            //var watch = System.Diagnostics.Stopwatch.StartNew();

            //for (int i = 0; i < 500; i++)
            //{
            //    straccum = "";
            //    order = 50;
            //    // normalized calcs, show
            //    AstroLibr.FullGeopM(recef, order, 'y', convArr, gravData, out aPertM, 'n', out straccum);
            //}
            ////  stop timer
            //watch.Stop();
            //var elapsedMs = watch.ElapsedMilliseconds;
            //strbuildall.AppendLine("Done with Montenbruck calcs " + (watch.ElapsedMilliseconds * 0.001).ToString() + " second  ");

            //watch = System.Diagnostics.Stopwatch.StartNew();
            //for (int i = 0; i < 500; i++)
            //{
            //    straccum = "";
            //    order = 50;
            //    // normalized calcs, show
            //    AstroLibr.FullGeopG(recef, order, 'y', convArr, normArr, gravData, out aPertG, 'n', out straccum);
            //}
            ////  stop timer
            //watch.Stop();
            //elapsedMs = watch.ElapsedMilliseconds;
            //strbuildall.AppendLine("Done with GTDS calcs " + (watch.ElapsedMilliseconds * 0.001).ToString() + " second  ");


            //watch = System.Diagnostics.Stopwatch.StartNew();
            //for (int i = 0; i < 500; i++)
            //{
            //    straccum = "";
            //    order = 100;
            //    // normalized calcs, show
            //    // GTDS version
            //    AstroLibr.LegPolyGTDS(latgc, order, 'y', out LegArrGU, out LegArrGN);
            //}
            ////  stop timer
            //watch.Stop();
            //elapsedMs = watch.ElapsedMilliseconds;
            //strbuildall.AppendLine("Done with GTDS ALF calcs " + (watch.ElapsedMilliseconds * 0.001).ToString() + " second  ");


            //watch = System.Diagnostics.Stopwatch.StartNew();
            //for (int i = 0; i < 500; i++)
            //{
            //    straccum = "";
            //    order = 100;
            //    // normalized calcs, show
            //    // Gottlieb version
            //    AstroLibr.LegPolyGottN(order, out norm1, out norm2, out norm11, out normn10, out norm1m, out norm2m, out normn1);
            //}
            ////  stop timer
            //watch.Stop();
            //elapsedMs = watch.ElapsedMilliseconds;
            //strbuildall.AppendLine("Done with Gott ALF calcs " + (watch.ElapsedMilliseconds * 0.001).ToString() + " second  ");

            //watch = System.Diagnostics.Stopwatch.StartNew();
            //for (int i = 0; i < 500; i++)
            //{
            //    straccum = "";
            //    order = 100;
            //    // normalized calcs, show
            //    // Montenbruck version
            //    AstroLibr.LegPolyMont(latgc, order, 'y', out LegArrMU, out LegArrMN);
            //}
            ////  stop timer
            //watch.Stop();
            //elapsedMs = watch.ElapsedMilliseconds;
            //strbuildall.AppendLine("Done with Mont ALF calcs " + (watch.ElapsedMilliseconds * 0.001).ToString() + " second  ");


            //watch = System.Diagnostics.Stopwatch.StartNew();
            //for (int i = 0; i < 500; i++)
            //{
            //    straccum = "";
            //    order = 100;
            //    // normalized calcs, show
            //    // Fukushima version
            //    //AstroLibr.LegPolyFN(latgc, order, 'y', out LegArrF);
            //}
            ////  stop timer
            //watch.Stop();
            //elapsedMs = watch.ElapsedMilliseconds;
            //strbuildall.AppendLine("Done with Fukushima ALF calcs " + (watch.ElapsedMilliseconds * 0.001).ToString() + " second  ");



            //// ------------------------------------------- pole test case comparisons
            //strbuildall.AppendLine("\n ===================================== Pole Test Comparisons =====================================");

            //rad = 180.0 / Math.PI;
            //for (int i = 0; i < 500; i++)
            //{
            //    lon = 154.0 / rad;
            //    latgc = (89.9 + (i / 1000.0)) / rad;
            //    double magr = 7378.382745;

            //    recef[0] = (magr * Math.Cos(latgc) * Math.Cos(lon));
            //    recef[1] = (magr * Math.Cos(latgc) * Math.Sin(lon));
            //    recef[2] = (magr * Math.Sin(latgc));

            //    straccum = "";
            //    order = 50;
            //    // normalized calcs, show
            //    AstroLibr.FullGeopG(recef, order, 'y', convArr, normArr, gravData, out aPertG, 'n', out straccum);

            //    strbuildall.AppendLine("test pole " + (latgc * rad).ToString() + " " + (lon * rad).ToString() + " " + aPertM[0].ToString() + "     " + aPertM[1].ToString() + "     " + aPertM[2].ToString());
            //}


            // available files:
            // GEM10Bunnorm36.grv
            // GEMT1norm36.grv
            // GEMT2norm36.grv
            // GEMT3norm36.grv
            // GEMT3norm50.grv
            // EGM-96norm70.grv
            // EGM-96norm254.grv
            // EGM-08norm100.grv
            // EGM-96norm70.grv
            // EGM-96norm70.grv
            // EGM-96norm70.grv
            // GGM01Cnorm90.grv
            // GGM02Cnorm90.grv
            // GGM03Cnorm100.grv
            // JGM2norm70.grv
            // JGM3norm70.grv
            // WGS-84_EGM96norm70.grv
            // WGS-84Enorm180.grv
            // WGS-84norm70.grv
            // 
            //string fname = "D:\Dataorig\Gravity\EGM96A.TXT";       // norm
            //string fname = "D:\Dataorig\Gravity\egm2008_gfc.txt";  // norm

            // --------------------gottlieb 1993 test
        //    strbuildall.AppendLine("===================================== Gottlieb 1993 test case ===================================== ");
        //    strbuildall.AppendLine("GEM-10B unnormalized 36x36 ");
        //    // get past text in each file
        //    //if (fname.Contains("GEM"))    // GEM10bunnorm36.grv, GEMT3norm50.grv
        //    //    startKtr = 17;
        //    //if (fname.Contains("EGM-96")) // EGM-96norm70.grv
        //    //    startKtr = 73;
        //    //if (fname.Contains("EGM-08")) // EGM-08norm100.grv
        //    //    startKtr = 83;  // or 21 for the larger file... which has gfc in the first col too
        //    fname = "D:/Dataorig/Gravity/GEM10Bunnorm36.grv";
        //    normalized = 'n';
        //    //double latgc;
        //    //Int32 degree, order;
        //    //double[,] LegArr;  // montenbruck
        //    //double[,] LegArrN;
        //    //double[,] LegArrG;  // gtds
        //    //double[,] LegArrGN;
        //    ////  double[,] LegArrEx;
        //    //double[,] LegArr1;  // geodyn

        //    //AstroLib.gravityModelData gravData;

        //    recef = new double[] { 5489.1500, 802.2220, 3140.9160 };  // km
        //    strbuildall.AppendLine("recef = " + recef[0].ToString() + " " + recef[1].ToString() + " " + recef[2].ToString());
        //    // these are from the vector
        //    latgc = Math.Asin(recef[2] / MathTimeLibr.mag(recef));
        //    double templ = Math.Sqrt(recef[0] * recef[0] + recef[1] * recef[1]);
        //    double rtasc;
        //    if (Math.Abs(templ) < 0.0000001)
        //        rtasc = Math.Sign(recef[2]) * Math.PI * 0.5;
        //    else
        //        rtasc = Math.Atan2(recef[1], recef[0]);
        //    lon = rtasc;
        //    strbuildall.AppendLine("latgc lon " + (latgc * rad).ToString() + " " + (lon * rad).ToString());

        //    this.opsStatus.Text = "Status: Reading gravity field Gottlieb test";
        //    Refresh();

        //    AstroLibr.readGravityField(fname, normalized, 17, out order, out gravData);
        //    strbuildall.AppendLine("\ncoefficents --------------- ");
        //    strbuildall.AppendLine("c  2  0  " + gravData.c[2, 0].ToString() + " s " + gravData.s[2, 0].ToString());
        //    strbuildall.AppendLine("c  4  0  " + gravData.c[4, 0].ToString() + " s " + gravData.s[4, 0].ToString());
        //    strbuildall.AppendLine("c  4  4  " + gravData.c[4, 4].ToString() + " s " + gravData.s[4, 4].ToString());
        //    strbuildall.AppendLine("c 21  1 " + gravData.c[21, 1].ToString() + " s " + gravData.s[21, 1].ToString());
        //    strbuildall.AppendLine("\nnormalized coefficents --------------- ");
        //    strbuildall.AppendLine("c  2  0  " + gravData.cNor[2, 0].ToString() + " s " + gravData.sNor[2, 0].ToString());
        //    strbuildall.AppendLine("c  4  0  " + gravData.cNor[4, 0].ToString() + " s " + gravData.sNor[4, 0].ToString());
        //    strbuildall.AppendLine("c  4  4  " + gravData.cNor[4, 4].ToString() + " s " + gravData.sNor[4, 4].ToString());
        //    strbuildall.AppendLine("c 21  1 " + gravData.cNor[21, 1].ToString() + " s " + gravData.sNor[21, 1].ToString());

        //    this.opsStatus.Text = "Status: Gottlieb test legpoly calcs";
        //    Refresh();

        //    order = 36;
        //    AstroLibr.LegPolyGTDS(latgc, order, normalized, out LegArrGU, out LegArrGN);
        //    // get geodyn version
        //    //AstroLibr.geodynlegp(latgc, degree, order, out LegArrOU, out LegArrON);
        //    // get exact values
        //    // LegPolyEx(latgc, order, out LegArrEx);

        //    errstr = " ";
        //    //sumdr1 = 0.0;
        //    //sumdr2 = 0.0;
        //    strbuildall.AppendLine("\nLegendre polynomials --------------- ");
        //    for (int L = 1; L <= 6; L++)  // order xxxxxxxxxxxxxxxxxx
        //    {
        //        string tempstr1 = "MN ";  // montenbruck
        //        string tempstr2 = "GN ";  // gtds
        //        //string tempstr3 = "MU ";
        //        string tempstr4 = "OU ";  // geodyn
        //        for (int m = 0; m <= L; m++)
        //        {
        //            tempstr1 = tempstr1 + " " + L.ToString() + "  " + m.ToString() + "   " + LegArrMN[L, m].ToString();
        //            tempstr2 = tempstr2 + " " + L.ToString() + "  " + m.ToString() + "   " + LegArrGN[L, m].ToString();
        //            //tempstr3 = tempstr3 + " " + L.ToString() + "  " + m.ToString() + "   " + LegArrMU[L, m].ToString();
        //            //tempstr4 = tempstr4 + " " + L.ToString() + "  " + m.ToString() + "   " + LegArrOU[L + 1, m + 1].ToString();
        //            //dr1 = 100.0 * (LegArr[L, m] - LegArrEx[L, m]) / LegArrEx[L, m];
        //            //dr2 = 100.0 * (LegArr1[L, m] - LegArrEx[L, m]) / LegArrEx[L, m];
        //            //sumdr1 = sumdr1 + dr1;
        //            //sumdr2 = sumdr2 + dr2;
        //            //errstr = errstr + "\n" + L.ToString() + "  " + m.ToString() + "   " + dr1.ToString()
        //            //    + " " + dr2.ToString();
        //        }
        //        strbuildall.AppendLine(tempstr1);
        //        strbuildall.AppendLine(tempstr2);
        //        //  strbuild.AppendLine(tempstr3);
        //        strbuildall.AppendLine(tempstr4 + "\n");
        //    }
        ////    strbuildall.AppendLine("totals gtds " + sumdr1.ToString() + " montenbruck " + sumdr2.ToString());
        //    strbuildplot.AppendLine(errstr);

        //    strbuildall.AppendLine("\naccelerations --------------- ");
        //    jdutc = 2451573.0;
        //    jdF = 0.1;
        //    straccum = "";
        //    order = 4;
        //    // normalized calcs, show
        //    AstroLibr.FullGeopM(recef, order, 'y', convArr, gravData, out aPertM, 'y', out straccum);
        //    // add in two body term since full geop is only disturbing part
        //    jdut1 = jdutc + jdF;
        //    //AstroLibr.eci_ecef(ref reci, ref veci, iau80arr, MathTimeLib.Edirection.efrom, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps, AstroLib.EOpt.e80, ref recef, ref vecef);
        //    // time is not given, so let ecef and eci be =
        //    reci[0] = recef[0];
        //    reci[1] = recef[1];
        //    reci[2] = recef[2];

        //    aeci2[0] = -398600.47 * reci[0] / (Math.Pow(MathTimeLibr.mag(reci), 3));
        //    aeci2[1] = -398600.47 * reci[1] / (Math.Pow(MathTimeLibr.mag(reci), 3));
        //    aeci2[2] = -398600.47 * reci[2] / (Math.Pow(MathTimeLibr.mag(reci), 3));
        //    //aPertG[0] = aPertG[0] + aeci2[0];
        //    //aPertG[1] = aPertG[1] + aeci2[1];
        //    //aPertG[2] = aPertG[2] + aeci2[2];
        //    aPertM[0] = aPertM[0] + aeci2[0];
        //    aPertM[1] = aPertM[1] + aeci2[1];
        //    aPertM[2] = aPertM[2] + aeci2[2];


        //    strbuildall.AppendLine(straccum);
        //    //   strbuild.AppendLine("apertG 4 4   " + aPertG[0].ToString() + "     " + aPertG[1].ToString() + "     " + aPertG[2].ToString());
        //    strbuildall.AppendLine("apertM 4 4   " + aPertM[0].ToString() + "     " + aPertM[1].ToString() + "     " + aPertM[2].ToString());
        //    strbuildall.AppendLine("ans          -0.00844269212018857E+00 -0.00123393633785485E+00 -0.00484659352346614E+00  km/s2  \n");

        //    straccum = "";
        //    order = 5;
        //    // normalized calcs, show
        //    AstroLibr.FullGeopG(recef, order, 'y', convArr, normArr, gravData, out aPertG, 'y', out straccum);
        //    strbuildall.AppendLine(straccum);
        //    strbuildall.AppendLine("apertG 5 5   " + aPertG[0].ToString() + "     " + aPertG[1].ToString() + "     " + aPertG[2].ToString());
        //    // strbuild.AppendLine("apertM 5 5   " + aPertM[0].ToString() + "     " + aPertM[1].ToString() + "     " + aPertM[2].ToString());
        //    strbuildall.AppendLine("ans          -0.00844260633555472E+00 -0.00123393243051834E+00 -0.00484652486332608E+00  km/s2  \n");



            // --------------------fonte 1993 test
            // strbuild.AppendLine("\n ===================================== Fonte 1993 test case =====================================");
            // strbuild.AppendLine("GEM-10B unnormalized 36x36 ");
            // fname = "D:/Dataorig/Gravity/GEM10Bunnorm36.grv";
            // normal = 'n';
            // recef = new double[] { 180.295260378399, -1145.13224944286, -6990.09446227757 }; // km
            // strbuild.AppendLine("recef = " + recef[0].ToString() + " " + recef[1].ToString() + " " + recef[2].ToString());
            // latgc = -1.40645188850273;
            // lon = -4.09449590512370;
            // strbuild.AppendLine("latgc lon " + (latgc * rad).ToString() + " " + (lon * rad).ToString());

            // this.opsStatus.Text = "Status: Reading gravity field Fonte test";
            // Refresh();

            // // Un-normalized Polynomial Validation GEM10B
            // // GTDS vs Lundberg Truth GTDS (21x21 GEM10B)
            // strbuild.AppendLine("\ncoefficients --------------- ");
            // AstroLibr.initGravityField(fname, normal, out gravData);
            // strbuild.AppendLine("c  4  0    " + gravData.c[4, 0].ToString() + " s " + gravData.s[4, 0].ToString());
            // strbuild.AppendLine("c 21  0   " + gravData.c[21, 0].ToString() + " s " + gravData.s[21, 0].ToString());
            // strbuild.AppendLine("c 21  5    " + gravData.c[21, 5].ToString() + " s " + gravData.s[21, 5].ToString());
            // strbuild.AppendLine("c 21 20   " + gravData.c[21, 20].ToString() + " s " + gravData.s[21, 20].ToString());
            // strbuild.AppendLine("c 21 21    " + gravData.c[21, 21].ToString() + " s " + gravData.s[21, 21].ToString());

            // // GTDS Emulation vs Lundberg Truth (21x21 GEM10B)
            // degree = 21;
            // order = 21;
            // AstroLibr.LegPoly(latgc, order, out LegArr, out LegArrG, out LegArrN, out LegArrGN);
            // // get geodyn version
            // AstroLibr.geodynlegp(latgc, degree, order, out LegArr1);
            // // get exact values
            // //   LegPolyEx(latgc, order, out LegArrEx);

            // dr1 = 0.0;
            // dr2 = 0.0;
            // sumdr1 = 0.0;
            // sumdr2 = 0.0;
            // strbuild.AppendLine("\nLegendre polynomials --------------- ");
            // for (int L = 1; L <= 6; L++)  // order
            // {
            //     string tempstr1 = "M ";
            //     string tempstr2 = "G ";
            //     string tempstr3 = "E ";
            //     string tempstr4 = "O ";
            //     for (int m = 0; m <= L; m++)
            //     {
            //         tempstr1 = tempstr1 + " " + L.ToString() + "  " + m.ToString() + "   " + LegArrN[L, m].ToString();
            //         tempstr2 = tempstr2 + " " + L.ToString() + "  " + m.ToString() + "   " + LegArrGN[L, m].ToString();
            //         // tempstr3 = tempstr3 + " " + L.ToString() + "  " + m.ToString() + "   " + LegArrEx[L, m].ToString();
            //         tempstr4 = tempstr4 + " " + L.ToString() + "  " + m.ToString() + "   " + LegArr1[L + 1, m + 1].ToString();
            //         //    dr1 = 100.0 * (LegArr[L, m] - LegArrEx[L, m]) / LegArrEx[L, m];
            //         //    dr2 = 100.0 * (LegArr1[L, m] - LegArrEx[L, m]) / LegArrEx[L, m];
            //         //sumdr1 = sumdr1 + dr1;
            //         //sumdr2 = sumdr2 + dr2;
            //         //errstr = errstr + "\n" + L.ToString() + "  " + m.ToString() + "   " + dr1.ToString()
            //         //    + " " + dr2.ToString();
            //     }
            //     strbuild.AppendLine(tempstr1);
            //     strbuild.AppendLine(tempstr2);
            //    // strbuild.AppendLine(tempstr3);
            //     strbuild.AppendLine(tempstr4 + "\n");
            // }
            //// strbuild.AppendLine("totals gtds " + sumdr1.ToString() + " montenbruck " + sumdr2.ToString());
            // strbuild.AppendLine("ans 21  0 0.385389365005720                                                                      21  5   354542.107743601  354542.1077435970657340");
            // strbuild.AppendLine("ans 21 20         -2442182686.11423  -2442182686.11409981594");
            // strbuild.AppendLine("ans 21 21          405012060.632803  405012060.6327805324689" + "\n");

            // strbuild.AppendLine("\naccelerations --------------- ");
            // AstroLibr.FullGeop(recef, jd, jdF, order, gravData, out aPert, out aPert1);

            // strbuild.AppendLine("apertG 21 21   " + aPert[0].ToString() + " " + aPert[1].ToString() + " " + aPert[2].ToString());
            // strbuild.AppendLine("apertM 21 21   " + aPert1[0].ToString() + "     " + aPert1[1].ToString() + "     " + aPert1[2].ToString());
            // strbuild.AppendLine("ans             8.653210294968294E-7  -6.515584998975128E-6  -1.931032474628621E-5 ");
            // strbuild.AppendLine("ans             8.653210294968E-7     -6.5155849989750E-6    -1.931032474628616E-5");

            // // --------------------fonte 1993 test
            // strbuild.AppendLine("\n===================================== Fonte 1993 test case =====================================");
            // strbuild.AppendLine("GEM-T3 normalized 50x50 ");
            // fname = "D:/Dataorig/Gravity/GEMT3norm50.grv";          // norm only released as 36x36 though...
            // normal = 'y';

            // this.opsStatus.Text = "Status: Reading gravity field fonte 93 test";
            // Refresh();

            // AstroLibr.initGravityField(fname, normal, out gravData);
            // strbuild.AppendLine("\ncoefficients --------------- ");
            // strbuild.AppendLine("c  4  0   " + gravData.c[4, 0].ToString() + " " + gravData.s[4, 0].ToString());
            // strbuild.AppendLine("c 21 20   " + gravData.c[21, 20].ToString() + " " + gravData.s[21, 20].ToString());
            // strbuild.AppendLine("c 50  0   " + gravData.c[50, 0].ToString() + " " + gravData.s[50, 0].ToString());
            // strbuild.AppendLine("c 50 50   " + gravData.c[50, 50].ToString() + " " + gravData.s[50, 50].ToString());
            // strbuild.AppendLine("c 50  5   " + gravData.c[50, 5].ToString() + " " + gravData.s[50, 5].ToString());

            // strbuild.AppendLine("\nLegendre polynomials --------------- ");
            // // GTDS Emulation vs Lundberg Truth (21x21 GEM10B)
            // degree = 50;
            // order = 50;

            // AstroLibr.LegPoly(latgc, order, out LegArr, out LegArrG, out LegArrN, out LegArrGN);
            // // get geodyn version
            // AstroLibr.geodynlegp(latgc, degree, order, out LegArr1);
            // // get exact
            // //  LegPolyEx(latgc, order, out LegArrEx);

            // strbuild.AppendLine("legarr4    0   " + LegArrN[4, 0].ToString() + " " + LegArrN[4, 1].ToString());

            // strbuild.AppendLine("50  0          " + LegArrN[50, 0].ToString());
            // strbuild.AppendLine("50  0 alt      " + LegArrGN[50, 0].ToString());
            // strbuild.AppendLine("ans 50  0      0.09634780379822722     9.634780379823085162E-02");
            // strbuild.AppendLine("50  0 geody    " + LegArr1[50, 0].ToString() + "\n");
            // //   strbuild.AppendLine("50  0 exact    " + LegArrEx[50, 0].ToString() + "\n");
            // //    strbuild.AppendLine("50  0 exact    " + LegArrEx[50, 0].ToString() + "\n");

            // strbuild.AppendLine("50 21       " + LegArrN[50, 21].ToString());
            // strbuild.AppendLine("50 21 alt   " + LegArrGN[50, 21].ToString());
            // strbuild.AppendLine("ans 50  21  -1.443200082785759E+28  -14432000827857661203015450149.6553");
            // strbuild.AppendLine("50 21 geody " + LegArr1[50, 21].ToString() + "\n");
            // //   strbuild.AppendLine("50 21 exact  " + LegArrEx[50, 21].ToString() + "\n");

            // strbuild.AppendLine("50 49       " + LegArrN[50, 49].ToString());
            // strbuild.AppendLine("50 49 alt   " + LegArrGN[50, 49].ToString());
            // strbuild.AppendLine("ans 50  49  -8.047341511222794E+39  -8.047341511222872818E+39");
            // strbuild.AppendLine("50 49 geody " + LegArr1[50, 49].ToString() + "\n");
            // // strbuild.AppendLine("50 49 exact " + ex5049.ToString() + "\n");

            // strbuild.AppendLine("50 50       " + LegArrN[50, 50].ToString());
            // strbuild.AppendLine("50 50 alt   " + LegArrGN[50, 50].ToString());
            // strbuild.AppendLine("ans 50 50      1.334572710963763E+39   1.334572710963775698E+39" + "\n");
            // strbuild.AppendLine("50 50 geody " + LegArr1[50, 50].ToString() + "\n");
            // //strbuild.AppendLine("50 50 exact " + ex550.ToString() + "\n");

            // strbuild.AppendLine("\naccelerations --------------- ");
            // normalized calcs, show
            // AstroLibr.FullGeop(recef, order, normalized, gravData, out aPert, out aPert1);
            // strbuild.AppendLine("apert 50 50   " + aPert[0].ToString() + " " + aPert[1].ToString() + " " + aPert[2].ToString());
            // strbuild.AppendLine("ans           8.683465146150188E-007    -6.519678538340073E-006   -1.931876804829165E-005");
            // strbuild.AppendLine("ans           8.68346514615019361E-07   -6.51967853834008023E-06  -1.93187680482916393E-05");

            // recef = new double[] { 487.0696937, -5330.5022406, 4505.7372146 };  // m
            // vecef = new double[] { -2.101083975, 4.624581986, 5.688300377 };

            // write out results
            string directory = @"d:\codes\library\matlab\";
            File.WriteAllText(directory + "legpoly.txt", strbuildall.ToString());

            File.WriteAllText(directory + "legendreAcc.txt", strbuildplot.ToString());

        }

        public void testhill()
        {
            double[] r, v, rh, vh, rint, vint;
            double alt, dts;
            r = new double[3];
            v = new double[3];
            rh = new double[3];
            vh = new double[3];
            // StringBuilder strbuild = new StringBuilder();

            dts = 1400.0; // second

            // circular orbit
            alt = 590.0;
            r[0] = AstroLibr.gravConst.re + alt;
            r[1] = 0.0;
            r[2] = 0.0;
            v[0] = 0.0;
            v[1] = Math.Sqrt(AstroLibr.gravConst.mu / MathTimeLibr.mag(r));
            v[2] = 0.0;

            rh[0] = 0.0;
            rh[1] = 0.0;
            rh[2] = 0.0;
            vh[0] = -0.1;
            vh[1] = -0.04;
            vh[2] = -0.02;

            for (int i = 1; i <= 50; i++)
            {
                dts = i * 60.0;  // second
                AstroLibr.hillsr(rh, vh, alt, dts, out rint, out vint);
                strbuild.AppendLine(dts.ToString() + " " + rint[0].ToString() + " " + rint[1].ToString() + " " + rint[2].ToString() +
                    " " + vint[0].ToString() + " " + vint[1].ToString() + " " + vint[2].ToString());
            }


            AstroLibr.hillsv(r, alt, dts, out vint);



        }  // test hill


        /* ----------------------------------------------------------------------------
*
*                                 function printcov
*
* this function prints a covariance matrix
*
* author        : david vallado                  719 - 573 - 2600   23 may 2003
*
* revisions
*
* inputs description range / units
*   covin     - 6x6 input covariance matrix
*   covtype   - type of covariance             'cl','ct','fl','sp','eq',
*   cu        - covariance units(deg or rad)  't' or 'm'
*   anom      - anomaly                        'mean' or 'true' or 'tau'
*
* outputs       :
*
* locals        :
*
* references    :
* none
*
* printcov(covin, covtype, cu, anom)
* ----------------------------------------------------------------------------*/

        public void printcov(double[,] covin, string covtype, char cu, string anom, out string strout)
        {
            int i;
            string semi = "";
            strout = "";

            if (anom.Equals("truea") || anom.Equals("meana"))
                semi = "a m  ";
            else
            {
                if (anom.Equals("truen") || anom.Equals("meann"))
                    semi = "n rad";
            }

            if (covtype.Equals("ct"))
            {
                strout = "cartesian covariance \n";
                strout = strout + "        x  m            y m             z  m           xdot  m/s       ydot  m/s       zdot  m/s  \n";
            }

            if (covtype.Equals("cl"))
            {
                strout = strout + "classical covariance \n";
                if (cu == 'm')
                {
                    strout = strout + "          " + semi + "          ecc           incl rad      raan rad         argp rad        ";
                    if (anom.Contains("mean")) // 'meana' 'meann'
                        strout = strout + "m rad \n";
                    else     // 'truea' 'truen'
                        strout = strout + " nu rad \n";
                }
                else
                {
                    strout = strout + "          " + semi + "           ecc           incl deg      raan deg         argp deg        ";
                    if (anom.Contains("mean")) // 'meana' 'meann'
                        strout = strout + " m deg \n";
                    else     // 'truea' 'truen'
                        strout = strout + " nu deg \n";
                }
            }

            if (covtype.Equals("eq"))
            {
                strout = strout + "equinoctial covariance \n";
                //            if (cu == 'm')
                if (anom.Contains("mean")) // 'meana' 'meann'
                    strout = strout + "         " + semi + "           af              ag           chi             psi         meanlonM rad\n";
                else     // 'truea' 'truen'
                    strout = strout + "         " + semi + "           af              ag           chi             psi         meanlonNu rad\n";
            }

            if (covtype.Equals("fl"))
            {
                strout = strout + "flight covariance \n";
                strout = strout + "       lon  rad      latgc rad        fpa rad         az rad           r  m           v  m/s  \n";
            }

            if (covtype.Equals("sp"))
            {
                strout = strout + "spherical covariance \n";
                strout = strout + "      rtasc deg       decl deg        fpa deg         az deg           r  m           v  m/s  \n";
            }

            // format strings to show signs "and" to not round off if trailing 0!!
            string fmt = "+#.#########0E+00;-#.#########0E+00";
            for (i = 0; i < 6; i++)
                strout = strout + covin[i, 0].ToString(fmt) + " " + covin[i, 1].ToString(fmt) + " " + covin[i, 2].ToString(fmt) + " " +
                 covin[i, 3].ToString(fmt) + " " + covin[i, 4].ToString(fmt) + " " + covin[i, 5].ToString(fmt) + "\n";
        }  // printcov


        /* ----------------------------------------------------------------------------
        *
        *                                  function printdiff
        *
        * this function prints a covariance matrix difference
        *
        * author        : david vallado                  719 - 573 - 2600   23 may 2003
        *
        * revisions
        *
        * inputs description range / units
        *   strin    - title
        *   mat1     - 6x6 input matrix
        *   mat2     - 6x6 input matrix
        *
        * outputs       :
        *
        * locals        :
        *
        * ----------------------------------------------------------------------------*/

        public void printdiff(string strin, double[,] mat1, double[,] mat2, out string strout)
        {
            double small = 1e-18;
            double[,] dr = new double[6, 6];
            double[,] diffmm = new double[6, 6];
            int i, j;

            // format strings to show signs "and" to not round off if trailing 0!!
            string fmt = "+#.#########0E+00;-#.#########0E+00";

            strout = "diff " + strin + "\n";
            for (i = 0; i < 6; i++)
            {
                for (j = 0; j < 6; j++)
                    dr[i, j] = mat1[i, j] - mat2[i, j];
                strout = strout + dr[i, 0].ToString(fmt) + " " + dr[i, 1].ToString(fmt) + " " + dr[i, 2].ToString(fmt) + " " +
                    dr[i, 3].ToString(fmt) + " " + dr[i, 4].ToString(fmt) + " " + dr[i, 5].ToString(fmt) + "\n";
            }

            strout = strout + "pctdiff % " + strin + " pct over 1e-18  \n";
            // fprintf(1, '%14.4f%14.4f%14.4f%14.4f%14.4f%14.4f \n', 100.0 * ((mat1' - mat2') / mat1'));
            // fprintf(1, 'Check consistency of both approaches tmct2cl-inv(tmcl2ct) diff pct over 1e-18 \n');
            // fprintf(1, '-------- accuracy of tm comparing ct2cl and cl2ct --------- \n');
            //tm1 = mat1';
            //tm2 = mat2';
            fmt = "+0.###0;-0.###0";
            for (i = 0; i < 6; i++)
            {
                for (j = 0; j < 6; j++)
                {
                    if (Math.Abs(dr[i, j]) < small || Math.Abs(mat1[i, j]) < small)
                        diffmm[i, j] = 0.0;
                    else
                        diffmm[i, j] = 100.0 * (dr[i, j] / mat1[i, j]);
                }
                strout = strout + diffmm[i, 0].ToString(fmt) + " " + diffmm[i, 1].ToString(fmt) + " " + diffmm[i, 2].ToString(fmt) + " " +
                     diffmm[i, 3].ToString(fmt) + " " + diffmm[i, 4].ToString(fmt) + " " + diffmm[i, 5].ToString(fmt) + "\n";
            }

        }  // printdiff



        public void testcovct2rsw()
        {
            int year, mon, day, hr, minute, timezone, dat;
            double second, dut1, ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac;
            double lod, xp, yp, ddpsi, ddeps;
            double[,] tm = new double[,] { { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 } };
            string anom = "meana";  // truea/n, meana/n
            double[] cartstate = new double[6];
            double[] recef = new double[3];
            double[] vecef = new double[3];
            double[] avec = new double[3];

            double[,] cartcovrsw = new double[6, 6];
            double[,] cartcovntw = new double[6, 6];
            double[,] tmct2cl = new double[6, 6];
            double[,] tmcl2ct = new double[6, 6];
            string strout;

            double[] reci = new double[3] { -605.79221660, -5870.22951108, 3493.05319896 };
            double[] veci = new double[3] { -1.56825429, -3.70234891, -6.47948395 };
            double[] aeci = new double[3] { 0.001, 0.002, 0.003 };

            // StringBuilder strbuild = new StringBuilder();
            // strbuild.Clear();

            EOPSPWLib.iau80Class iau80arr;
            string fileLoc;
            fileLoc = @"D:\Codes\LIBRARY\DataLib\nut80.dat";
            EOPSPWLibr.iau80in(fileLoc, out iau80arr);

            year = 2000;
            mon = 12;
            day = 15;
            hr = 16;
            minute = 58;
            second = 50.208;
            dut1 = 0.10597;
            dat = 32;
            timezone = 0;
            xp = 0.0;
            yp = 0.0;
            lod = 0.0;
            timezone = 0;
            ddpsi = 0.0;
            ddeps = 0.0;

            MathTimeLibr.convtime(year, mon, day, hr, minute, second, timezone, dut1, dat,
                out ut1, out tut1, out jdut1, out jdut1frac, out utc, out tai,
                out tt, out ttt, out jdtt, out jdttfrac,
                out tdb, out ttdb, out jdtdb, out jdtdbfrac);

            // ---convert the eci state into the various other state formats(classical, equinoctial, etc)
            double[,] cartcov = new double[,] {
            { 100.0, 1.0e-2, 1.0e-2, 1.0e-4, 1.0e-4, 1.0e-4 },
            {1.0e-2, 100.0,  1.0e-2, 1.0e-4,   1.0e-4,   1.0e-4},
            {1.0e-2, 1.0e-2, 100.0,  1.0e-4,   1.0e-4,   1.0e-4},
            {1.0e-4, 1.0e-4, 1.0e-4, 0.0001,   1.0e-6,   1.0e-6},
            {1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   0.0001,   1.0e-6},
            {1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   1.0e-6,   0.0001} };
            cartstate = new double[] { reci[0], reci[1], reci[2], veci[0], veci[1], veci[2] };  // in km


            // test position and velocity going back
            avec = new double[] { 0.0, 0.0, 0.0 };

            AstroLibr.eci_ecef(ref reci, ref veci, MathTimeLib.Edirection.eto, ref recef, ref vecef,
                  iau80arr, jdtt, jdttfrac, jdut1, lod, xp, yp, ddpsi, ddeps);

            strbuild.AppendLine("==================== do the sensitivity tests \n");
            strbuild.AppendLine("1.  Cartesian Covariance \n");
            printcov(cartcov, "ct", 'm', anom, out strout);
            strbuild.AppendLine(strout);

            strbuild.AppendLine("2.  RSW Covariance from Cartesian #1 above  ------------------- \n");
            AstroLibr.covct_rsw(ref cartcov, cartstate, MathTimeLib.Edirection.eto, ref cartcovrsw, out tmct2cl);
            printcov(cartcovrsw, "ct", 'm', anom, out strout);
            strbuild.AppendLine(strout);

            strbuild.AppendLine("2.  NTW Covariance from Cartesian #1 above  ------------------- \n");
            AstroLibr.covct_ntw(ref cartcov, cartstate, MathTimeLib.Edirection.eto, ref cartcovntw, out tmct2cl);
            printcov(cartcovntw, "ct", 'm', anom, out strout);
            strbuild.AppendLine(strout);
            strbuild.AppendLine("\n");
        }


        public void testcovct2ntw()
        {

        }


        // test eci_ecef too
        public void testcovct2clmean()
        {
            int year, mon, day, hr, minute, timezone, dat, terms;
            double second, dut1, ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac;
            double p, a, n, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper;
            double af, ag, chi, psi, meanlonNu, meanlonM;
            double latgc, lon, rtasc, decl, fpa, lod, xp, yp, ddpsi, ddeps, az, magr, magv;
            Int16 fr;
            double[,] tm = new double[,] { { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 } };
            string anom = "meana";  // truea/n, meana/n
            string anomflt = "latlon"; // latlon  radec
            double[] cartstate = new double[6];
            double[] classstate = new double[6];
            double[] eqstate = new double[6];
            double[] fltstate = new double[6];
            double[] recef = new double[3];
            double[] vecef = new double[3];
            double[] avec = new double[3];

            double[,] classcovmeana = new double[6, 6];
            double[,] cartcovmeanarev = new double[6, 6];
            double[,] tmct2cl = new double[6, 6];
            double[,] tmcl2ct = new double[6, 6];
            string strout;

            double[] reci = new double[3] { -605.79221660, -5870.22951108, 3493.05319896 };
            double[] veci = new double[3] { -1.56825429, -3.70234891, -6.47948395 };
            double[] aeci = new double[3] { 0.001, 0.002, 0.003 };

            // StringBuilder strbuild = new StringBuilder();
            // strbuild.Clear();

            EOPSPWLib.iau80Class iau80arr;
            string fileLoc;
            fileLoc = @"D:\Codes\LIBRARY\DataLib\nut80.dat";
            EOPSPWLibr.iau80in(fileLoc, out iau80arr);

            year = 2000;
            mon = 12;
            day = 15;
            hr = 16;
            minute = 58;
            second = 50.208;
            dut1 = 0.10597;
            dat = 32;
            timezone = 0;
            xp = 0.0;
            yp = 0.0;
            lod = 0.0;
            terms = 2;
            timezone = 0;
            ddpsi = 0.0;
            ddeps = 0.0;

            MathTimeLibr.convtime(year, mon, day, hr, minute, second, timezone, dut1, dat,
                out ut1, out tut1, out jdut1, out jdut1frac, out utc, out tai,
                out tt, out ttt, out jdtt, out jdttfrac,
                out tdb, out ttdb, out jdtdb, out jdtdbfrac);

            // ---convert the eci state into the various other state formats(classical, equinoctial, etc)
            double[,] cartcov = new double[,] {
            { 100.0, 1.0e-2, 1.0e-2, 1.0e-4, 1.0e-4, 1.0e-4 },
            {1.0e-2, 100.0,  1.0e-2, 1.0e-4,   1.0e-4,   1.0e-4},
            {1.0e-2, 1.0e-2, 100.0,  1.0e-4,   1.0e-4,   1.0e-4},
            {1.0e-4, 1.0e-4, 1.0e-4, 0.0001,   1.0e-6,   1.0e-6},
            {1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   0.0001,   1.0e-6},
            {1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   1.0e-6,   0.0001} };
            cartstate = new double[] { reci[0], reci[1], reci[2], veci[0], veci[1], veci[2] };  // in km

            // --------convert to a classical orbit state
            AstroLibr.rv2coe(reci, veci,
                out p, out a, out ecc, out incl, out raan, out argp, out nu, out m, out arglat, out truelon, out lonper);
            classstate[0] = a;   // km
            classstate[1] = ecc;
            classstate[2] = incl;
            classstate[3] = raan;
            classstate[4] = argp;
            if (anom.Contains("mean")) // meann or meana
                classstate[5] = m;
            else  // truea or truen
                classstate[5] = nu;

            // -------- convert to an equinoctial orbit state
            AstroLibr.rv2eq(reci, veci, out a, out n, out af, out ag, out chi, out psi, out meanlonM, out meanlonNu, out fr);
            if (anom.Equals("meana") || anom.Equals("truea"))
                eqstate[0] = a;  // km
            else // meann or truen
                eqstate[0] = n;
            eqstate[1] = af;
            eqstate[2] = ag;
            eqstate[3] = chi;
            eqstate[4] = psi;
            if (anom.Contains("mean")) //  meana or meann
                eqstate[5] = meanlonM;
            else // truea or truen
                eqstate[5] = meanlonNu;

            // --------convert to a flight orbit state
            AstroLibr.rv2flt(reci, veci, jdtt, jdttfrac, jdut1, lod, xp, yp, terms, ddpsi, ddeps,
                iau80arr, out lon, out latgc, out rtasc, out decl, out fpa, out az, out magr, out magv);
            if (anomflt.Equals("radec"))
            {
                fltstate[0] = rtasc;
                fltstate[1] = decl;
            }
            else
            if (anomflt.Equals("latlon"))
            {
                fltstate[0] = lon;
                fltstate[1] = latgc;
            }
            fltstate[2] = fpa;
            fltstate[3] = az;
            fltstate[4] = magr;  // km
            fltstate[5] = magv;

            // test position and velocity going back
            avec = new double[] { 0.0, 0.0, 0.0 };

            AstroLibr.eci_ecef(ref reci, ref veci, MathTimeLib.Edirection.eto, ref recef, ref vecef,
                iau80arr, jdtt, jdttfrac, jdut1, lod, xp, yp, ddpsi, ddeps);
            //vx = magv* ( -cos(lon)*sin(latgc)*cos(az)*cos(fpa) - sin(lon)*sin(az)*cos(fpa) + cos(lon)*cos(latgc)*sin(fpa) ); 
            //vy = magv* ( -sin(lon)*sin(latgc)*cos(az)*cos(fpa) + cos(lon)*sin(az)*cos(fpa) + sin(lon)*cos(latgc)*sin(fpa) );  
            //vz = magv* (sin(latgc) * sin(fpa) + cos(latgc)*cos(az)*cos(fpa) );
            //// correct:
            //ve1 = magv* ( -cos(rtasc)*sin(decl)*cos(az)*cos(fpa) - sin(rtasc)*sin(az)*cos(fpa) + cos(rtasc)*cos(decl)*sin(fpa) ); // m/s
            // ve2 = magv* (-sin(rtasc) * sin(decl) * cos(az) * cos(fpa) + cos(rtasc) * sin(az) * cos(fpa) + sin(rtasc) * cos(decl) * sin(fpa));
            //ve3 = magv* (sin(decl) * sin(fpa) + cos(decl)*cos(az)*cos(fpa) );

            strbuild.AppendLine("==================== do the sensitivity tests \n");

            strbuild.AppendLine("1.  Cartesian Covariance \n");
            printcov(cartcov, "ct", 'm', anom, out strout);
            strbuild.AppendLine(strout);

            strbuild.AppendLine("2.  Classical Covariance from Cartesian #1 above (" + anom + ") ------------------- \n");

            AstroLibr.covct2cl(cartcov, cartstate, anom, out classcovmeana, out tmct2cl);
            printcov(classcovmeana, "cl", 'm', anom, out strout);
            strbuild.AppendLine(strout);

            strbuild.AppendLine("  Cartesian Covariance from Classical #2 above \n");
            AstroLibr.covcl2ct(classcovmeana, classstate, anom, out cartcovmeanarev, out tmcl2ct);
            printcov(cartcovmeanarev, "ct", 'm', anom, out strout);
            strbuild.AppendLine(strout);
            strbuild.AppendLine("\n");

            printdiff(" cartcov - cartcovmeanarev \n", cartcov, cartcovmeanarev, out strout);
            strbuild.AppendLine(strout);

            double[,] ecefcartcov = new double[6, 6];

            //AstroLibr.coveci_ecef(ref cartcov, cartstate, MathTimeLib.Edirection.eto,  ref ecefcartcov, out tm, iau80arr,
            //            ttt, jdut1, lod, xp, yp, 2, ddpsi, ddeps, AstroLib.EOpt.e80);
            //printcov(cartcovmeanarev, "ct", 'm', anom, out strout);
            //strbuild.AppendLine(strout);
            //strbuild.AppendLine("\n");

        }  // testcovct2clmean


        public void testcovct2cltrue()
        {
            int year, mon, day, hr, minute, timezone, dat, terms;
            double second, dut1, ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac;
            double p, a, n, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper;
            double af, ag, chi, psi, meanlonNu, meanlonM;
            double latgc, lon, rtasc, decl, fpa, lod, xp, yp, ddpsi, ddeps, az, magr, magv;
            Int16 fr;
            double[,] tm = new double[,] { { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 } };
            string anom = "truea";  // truea/n, meana/n
            string anomflt = "latlon"; // latlon  radec
            double[] cartstate = new double[6];
            double[] classstate = new double[6];
            double[] eqstate = new double[6];
            double[] fltstate = new double[6];
            double[] recef = new double[3];
            double[] vecef = new double[3];
            double[] avec = new double[3];

            double[,] classcovtruea = new double[6, 6];
            double[,] cartcovtruearev = new double[6, 6];
            double[,] tmct2cl = new double[6, 6];
            double[,] tmcl2ct = new double[6, 6];
            string strout;

            double[] reci = new double[3] { -605.79221660, -5870.22951108, 3493.05319896 };
            double[] veci = new double[3] { -1.56825429, -3.70234891, -6.47948395 };
            double[] aeci = new double[3] { 0.001, 0.002, 0.003 };

            //StringBuilder strbuild = new StringBuilder();
            //strbuild.Clear();

            EOPSPWLib.iau80Class iau80arr;
            string fileLoc;
            fileLoc = @"D:\Codes\LIBRARY\DataLib\nut80.dat";
            EOPSPWLibr.iau80in(fileLoc, out iau80arr);

            year = 2000;
            mon = 12;
            day = 15;
            hr = 16;
            minute = 58;
            second = 50.208;
            dut1 = 0.10597;
            dat = 32;
            timezone = 0;
            xp = 0.0;
            yp = 0.0;
            lod = 0.0;
            terms = 2;
            timezone = 0;
            ddpsi = 0.0;
            ddeps = 0.0;

            MathTimeLibr.convtime(year, mon, day, hr, minute, second, timezone, dut1, dat,
                out ut1, out tut1, out jdut1, out jdut1frac, out utc, out tai,
                out tt, out ttt, out jdtt, out jdttfrac,
                out tdb, out ttdb, out jdtdb, out jdtdbfrac);

            // ---convert the eci state into the various other state formats(classical, equinoctial, etc)
            double[,] cartcov = new double[,] {
            { 100.0, 1.0e-2, 1.0e-2, 1.0e-4, 1.0e-4, 1.0e-4 },
            {1.0e-2, 100.0,  1.0e-2, 1.0e-4,   1.0e-4,   1.0e-4},
            {1.0e-2, 1.0e-2, 100.0,  1.0e-4,   1.0e-4,   1.0e-4},
            {1.0e-4, 1.0e-4, 1.0e-4, 0.0001,   1.0e-6,   1.0e-6},
            {1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   0.0001,   1.0e-6},
            {1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   1.0e-6,   0.0001} };
            cartstate = new double[] { reci[0], reci[1], reci[2], veci[0], veci[1], veci[2] };  // in km

            // --------convert to a classical orbit state
            AstroLibr.rv2coe(reci, veci,
                out p, out a, out ecc, out incl, out raan, out argp, out nu, out m, out arglat, out truelon, out lonper);
            classstate[0] = a;  // in km
            classstate[1] = ecc;
            classstate[2] = incl;
            classstate[3] = raan;
            classstate[4] = argp;
            if (anom.Contains("mean")) // meann or meana
                classstate[5] = m;
            else  // truea or truen
                classstate[5] = nu;

            // -------- convert to an equinoctial orbit state
            AstroLibr.rv2eq(reci, veci, out a, out n, out af, out ag, out chi, out psi, out meanlonM, out meanlonNu, out fr);
            if (anom.Equals("meana") || anom.Equals("truea"))
                eqstate[0] = a;  // km
            else // meann or truen
                eqstate[0] = n;
            eqstate[1] = af;
            eqstate[2] = ag;
            eqstate[3] = chi;
            eqstate[4] = psi;
            if (anom.Contains("mean")) //  meana or meann
                eqstate[5] = meanlonM;
            else // truea or truen
                eqstate[5] = meanlonNu;

            // --------convert to a flight orbit state
            AstroLibr.rv2flt(reci, veci, jdtt, jdttfrac, jdut1, lod, xp, yp, terms, ddpsi, ddeps,
                iau80arr, out lon, out latgc, out rtasc, out decl, out fpa, out az, out magr, out magv);
            if (anomflt.Equals("radec"))
            {
                fltstate[0] = rtasc;
                fltstate[1] = decl;
            }
            else
            if (anomflt.Equals("latlon"))
            {
                fltstate[0] = lon;
                fltstate[1] = latgc;
            }
            fltstate[2] = fpa;
            fltstate[3] = az;
            fltstate[4] = magr;  // in km
            fltstate[5] = magv;

            // test position and velocity going back
            avec = new double[] { 0.0, 0.0, 0.0 };

            AstroLibr.eci_ecef(ref reci, ref veci, MathTimeLib.Edirection.eto, ref recef, ref vecef,
                iau80arr, jdtt, jdttfrac, jdut1, lod, xp, yp, ddpsi, ddeps);
            //vx = magv* ( -cos(lon)*sin(latgc)*cos(az)*cos(fpa) - sin(lon)*sin(az)*cos(fpa) + cos(lon)*cos(latgc)*sin(fpa) ); 
            //vy = magv* ( -sin(lon)*sin(latgc)*cos(az)*cos(fpa) + cos(lon)*sin(az)*cos(fpa) + sin(lon)*cos(latgc)*sin(fpa) );  
            //vz = magv* (sin(latgc) * sin(fpa) + cos(latgc)*cos(az)*cos(fpa) );
            //// correct:
            //ve1 = magv* ( -cos(rtasc)*sin(decl)*cos(az)*cos(fpa) - sin(rtasc)*sin(az)*cos(fpa) + cos(rtasc)*cos(decl)*sin(fpa) ); // m/s
            // ve2 = magv* (-sin(rtasc) * sin(decl) * cos(az) * cos(fpa) + cos(rtasc) * sin(az) * cos(fpa) + sin(rtasc) * cos(decl) * sin(fpa));
            //ve3 = magv* (sin(decl) * sin(fpa) + cos(decl)*cos(az)*cos(fpa) );

            strbuild.AppendLine("==================== do the sensitivity tests \n");

            strbuild.AppendLine("1.  Cartesian Covariance \n");
            printcov(cartcov, "ct", 'm', anom, out strout);
            strbuild.AppendLine(strout);

            strbuild.AppendLine("2.  Classical Covariance from Cartesian #1 above (" + anom + ") ------------------- \n");

            AstroLibr.covct2cl(cartcov, cartstate, anom, out classcovtruea, out tmct2cl);
            printcov(classcovtruea, "cl", 'm', anom, out strout);
            strbuild.AppendLine(strout);

            strbuild.AppendLine("  Cartesian Covariance from Classical #2 above \n");
            AstroLibr.covcl2ct(classcovtruea, classstate, anom, out cartcovtruearev, out tmcl2ct);
            printcov(cartcovtruearev, "ct", 'm', anom, out strout);
            strbuild.AppendLine(strout);
            strbuild.AppendLine("\n");

            printdiff(" cartcov - cartcovtruearev \n", cartcov, cartcovtruearev, out strout);
            strbuild.AppendLine(strout);
        }  // testcovct2cltrue



        public void testcovcl2eq(string anom)
        {
            int year, mon, day, hr, minute, timezone, dat, terms;
            double second, dut1, ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac;
            double p, a, n, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper;
            double af, ag, chi, psi, meanlonNu, meanlonM;
            double latgc, lon, rtasc, decl, fpa, lod, xp, yp, ddpsi, ddeps, az, magr, magv;
            Int16 fr;
            double[,] tm = new double[,] { { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 } };
            string anomflt = "latlon"; // latlon  radec
            double[] cartstate = new double[6];
            double[] classstate = new double[6];
            double[] eqstate = new double[6];
            double[] fltstate = new double[6];
            double[] recef = new double[3];
            double[] vecef = new double[3];
            double[] avec = new double[3];

            double[,] classcovmeana = new double[6, 6];
            double[,] cartcovmeanarev = new double[6, 6];
            double[,] eqcovmeana = new double[6, 6];
            double[,] tmct2cl = new double[6, 6];
            double[,] tmcl2ct = new double[6, 6];
            double[,] tmcl2eq = new double[6, 6];
            double[,] tmeq2cl = new double[6, 6];
            string strout;

            double[] reci = new double[3] { -605.79221660, -5870.22951108, 3493.05319896 };
            double[] veci = new double[3] { -1.56825429, -3.70234891, -6.47948395 };
            double[] aeci = new double[3] { 0.001, 0.002, 0.003 };

            // StringBuilder strbuild = new StringBuilder();
            // strbuild.Clear();

            EOPSPWLib.iau80Class iau80arr;
            string fileLoc;
            fileLoc = @"D:\Codes\LIBRARY\DataLib\nut80.dat";
            EOPSPWLibr.iau80in(fileLoc, out iau80arr);

            year = 2000;
            mon = 12;
            day = 15;
            hr = 16;
            minute = 58;
            second = 50.208;
            dut1 = 0.10597;
            dat = 32;
            timezone = 0;
            xp = 0.0;
            yp = 0.0;
            lod = 0.0;
            terms = 2;
            timezone = 0;
            ddpsi = 0.0;
            ddeps = 0.0;

            MathTimeLibr.convtime(year, mon, day, hr, minute, second, timezone, dut1, dat,
                out ut1, out tut1, out jdut1, out jdut1frac, out utc, out tai,
                out tt, out ttt, out jdtt, out jdttfrac,
                out tdb, out ttdb, out jdtdb, out jdtdbfrac);

            // ---convert the eci state into the various other state formats(classical, equinoctial, etc)
            double[,] cartcov = new double[,] {
            { 100.0, 1.0e-2, 1.0e-2, 1.0e-4, 1.0e-4, 1.0e-4 },
            {1.0e-2, 100.0,  1.0e-2, 1.0e-4,   1.0e-4,   1.0e-4},
            {1.0e-2, 1.0e-2, 100.0,  1.0e-4,   1.0e-4,   1.0e-4},
            {1.0e-4, 1.0e-4, 1.0e-4, 0.0001,   1.0e-6,   1.0e-6},
            {1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   0.0001,   1.0e-6},
            {1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   1.0e-6,   0.0001} };
            cartstate = new double[] { reci[0], reci[1], reci[2], veci[0], veci[1], veci[2] };  // in km

            // --------convert to a classical orbit state
            AstroLibr.rv2coe(reci, veci,
                out p, out a, out ecc, out incl, out raan, out argp, out nu, out m, out arglat, out truelon, out lonper);
            classstate[0] = a;   // km
            classstate[1] = ecc;
            classstate[2] = incl;
            classstate[3] = raan;
            classstate[4] = argp;
            if (anom.Contains("mean")) // meann or meana
                classstate[5] = m;
            else  // truea or truen
                classstate[5] = nu;

            // -------- convert to an equinoctial orbit state
            AstroLibr.rv2eq(reci, veci, out a, out n, out af, out ag, out chi, out psi, out meanlonM, out meanlonNu, out fr);
            if (anom.Equals("meana") || anom.Equals("truea"))
                eqstate[0] = a;  // km
            else // meann or truen
                eqstate[0] = n;
            eqstate[1] = af;
            eqstate[2] = ag;
            eqstate[3] = chi;
            eqstate[4] = psi;
            if (anom.Contains("mean")) //  meana or meann
                eqstate[5] = meanlonM;
            else // truea or truen
                eqstate[5] = meanlonNu;

            // --------convert to a flight orbit state
            AstroLibr.rv2flt(reci, veci, jdtt, jdttfrac, jdut1, lod, xp, yp, terms, ddpsi, ddeps,
                iau80arr, out lon, out latgc, out rtasc, out decl, out fpa, out az, out magr, out magv);
            if (anomflt.Equals("radec"))
            {
                fltstate[0] = rtasc;
                fltstate[1] = decl;
            }
            else
            if (anomflt.Equals("latlon"))
            {
                fltstate[0] = lon;
                fltstate[1] = latgc;
            }
            fltstate[2] = fpa;
            fltstate[3] = az;
            fltstate[4] = magr;  // km
            fltstate[5] = magv;

            // test position and velocity going back
            avec = new double[] { 0.0, 0.0, 0.0 };

            AstroLibr.eci_ecef(ref reci, ref veci, MathTimeLib.Edirection.eto, ref recef, ref vecef,
                  iau80arr, jdtt, jdttfrac, jdut1, lod, xp, yp, ddpsi, ddeps);
            //vx = magv* ( -cos(lon)*sin(latgc)*cos(az)*cos(fpa) - sin(lon)*sin(az)*cos(fpa) + cos(lon)*cos(latgc)*sin(fpa) ); 
            //vy = magv* ( -sin(lon)*sin(latgc)*cos(az)*cos(fpa) + cos(lon)*sin(az)*cos(fpa) + sin(lon)*cos(latgc)*sin(fpa) );  
            //vz = magv* (sin(latgc) * sin(fpa) + cos(latgc)*cos(az)*cos(fpa) );
            //// correct:
            //ve1 = magv* ( -cos(rtasc)*sin(decl)*cos(az)*cos(fpa) - sin(rtasc)*sin(az)*cos(fpa) + cos(rtasc)*cos(decl)*sin(fpa) ); // m/s
            // ve2 = magv* (-sin(rtasc) * sin(decl) * cos(az) * cos(fpa) + cos(rtasc) * sin(az) * cos(fpa) + sin(rtasc) * cos(decl) * sin(fpa));
            //ve3 = magv* (sin(decl) * sin(fpa) + cos(decl)*cos(az)*cos(fpa) );

            strbuild.AppendLine("==================== do the sensitivity tests \n");

            strbuild.AppendLine("1.  Cartesian Covariance \n");
            printcov(cartcov, "ct", 'm', anom, out strout);
            strbuild.AppendLine(strout);

            strbuild.AppendLine("3.  Equinoctial Covariance from Classical (Cartesian) #1 above (" + anom + ") ------------------- \n");
            AstroLibr.covct2cl(cartcov, cartstate, anom, out classcovmeana, out tmct2cl);
            AstroLibr.covcl2eq(classcovmeana, classstate, anom, anom, fr, out eqcovmeana, out tmcl2eq);

            printcov(eqcovmeana, "eq", 'm', anom, out strout);
            strbuild.AppendLine(strout);

            strbuild.AppendLine("  Cartesian Covariance from Classical #3 above \n");
            AstroLibr.coveq2cl(eqcovmeana, eqstate, anom, anom, fr, out classcovmeana, out tmeq2cl);
            printcov(classcovmeana, "cl", 'm', anom, out strout);
            strbuild.AppendLine(strout);

            AstroLibr.covcl2ct(classcovmeana, classstate, anom, out cartcovmeanarev, out tmcl2ct);
            printcov(cartcovmeanarev, "ct", 'm', anom, out strout);
            strbuild.AppendLine(strout);
            strbuild.AppendLine("\n");

            printdiff(" cartcov - cartcov" + anom + "rev \n", cartcov, cartcovmeanarev, out strout);
            strbuild.AppendLine(strout);
        }  // testcovcl2eq

        public void testcovct2eq(string anom)
        {
            int year, mon, day, hr, minute, timezone, dat, terms;
            double second, dut1, ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac;
            double p, a, n, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper;
            double af, ag, chi, psi, meanlonNu, meanlonM;
            double latgc, lon, rtasc, decl, fpa, lod, xp, yp, ddpsi, ddeps, az, magr, magv;
            Int16 fr;
            double[,] tm = new double[,] { { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 } };
            string anomflt = "latlon"; // latlon  radec
            double[] cartstate = new double[6];
            double[] classstate = new double[6];
            double[] eqstate = new double[6];
            double[] fltstate = new double[6];
            double[] recef = new double[3];
            double[] vecef = new double[3];
            double[] avec = new double[3];

            double[,] classcovmeana = new double[6, 6];
            double[,] cartcovmeanarev = new double[6, 6];
            double[,] eqcovmeana = new double[6, 6];
            double[,] tmct2cl = new double[6, 6];
            double[,] tmcl2ct = new double[6, 6];
            double[,] tmct2eq = new double[6, 6];
            double[,] tmeq2ct = new double[6, 6];
            string strout;

            double[] reci = new double[3] { -605.79221660, -5870.22951108, 3493.05319896 };
            double[] veci = new double[3] { -1.56825429, -3.70234891, -6.47948395 };
            double[] aeci = new double[3] { 0.001, 0.002, 0.003 };

            // StringBuilder strbuild = new StringBuilder();
            // strbuild.Clear();

            EOPSPWLib.iau80Class iau80arr;
            string fileLoc;
            fileLoc = @"D:\Codes\LIBRARY\DataLib\nut80.dat";
            EOPSPWLibr.iau80in(fileLoc, out iau80arr);

            year = 2000;
            mon = 12;
            day = 15;
            hr = 16;
            minute = 58;
            second = 50.208;
            dut1 = 0.10597;
            dat = 32;
            timezone = 0;
            xp = 0.0;
            yp = 0.0;
            lod = 0.0;
            terms = 2;
            timezone = 0;
            ddpsi = 0.0;
            ddeps = 0.0;

            MathTimeLibr.convtime(year, mon, day, hr, minute, second, timezone, dut1, dat,
                out ut1, out tut1, out jdut1, out jdut1frac, out utc, out tai,
                out tt, out ttt, out jdtt, out jdttfrac,
                out tdb, out ttdb, out jdtdb, out jdtdbfrac);

            // ---convert the eci state into the various other state formats(classical, equinoctial, etc)
            double[,] cartcov = new double[,] {
            { 100.0, 1.0e-2, 1.0e-2, 1.0e-4, 1.0e-4, 1.0e-4 },
            {1.0e-2, 100.0,  1.0e-2, 1.0e-4,   1.0e-4,   1.0e-4},
            {1.0e-2, 1.0e-2, 100.0,  1.0e-4,   1.0e-4,   1.0e-4},
            {1.0e-4, 1.0e-4, 1.0e-4, 0.0001,   1.0e-6,   1.0e-6},
            {1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   0.0001,   1.0e-6},
            {1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   1.0e-6,   0.0001} };
            cartstate = new double[] { reci[0], reci[1], reci[2], veci[0], veci[1], veci[2] };  // in km

            // --------convert to a classical orbit state
            AstroLibr.rv2coe(reci, veci,
                out p, out a, out ecc, out incl, out raan, out argp, out nu, out m, out arglat, out truelon, out lonper);
            classstate[0] = a;   // km
            classstate[1] = ecc;
            classstate[2] = incl;
            classstate[3] = raan;
            classstate[4] = argp;
            if (anom.Contains("mean")) // meann or meana
                classstate[5] = m;
            else  // truea or truen
                classstate[5] = nu;

            // -------- convert to an equinoctial orbit state
            AstroLibr.rv2eq(reci, veci, out a, out n, out af, out ag, out chi, out psi, out meanlonM, out meanlonNu, out fr);
            if (anom.Equals("meana") || anom.Equals("truea"))
                eqstate[0] = a;  // km
            else // meann or truen
                eqstate[0] = n;
            eqstate[1] = af;
            eqstate[2] = ag;
            eqstate[3] = chi;
            eqstate[4] = psi;
            if (anom.Contains("mean")) //  meana or meann
                eqstate[5] = meanlonM;
            else // truea or truen
                eqstate[5] = meanlonNu;

            // --------convert to a flight orbit state
            AstroLibr.rv2flt(reci, veci, jdtt, jdttfrac, jdut1, lod, xp, yp, terms, ddpsi, ddeps,
                iau80arr, out lon, out latgc, out rtasc, out decl, out fpa, out az, out magr, out magv);
            if (anomflt.Equals("radec"))
            {
                fltstate[0] = rtasc;
                fltstate[1] = decl;
            }
            else
            if (anomflt.Equals("latlon"))
            {
                fltstate[0] = lon;
                fltstate[1] = latgc;
            }
            fltstate[2] = fpa;
            fltstate[3] = az;
            fltstate[4] = magr;  // km
            fltstate[5] = magv;

            // test position and velocity going back
            avec = new double[] { 0.0, 0.0, 0.0 };

            AstroLibr.eci_ecef(ref reci, ref veci, MathTimeLib.Edirection.eto, ref recef, ref vecef,
                  iau80arr, jdtt, jdttfrac, jdut1, lod, xp, yp, ddpsi, ddeps);
            //vx = magv* ( -cos(lon)*sin(latgc)*cos(az)*cos(fpa) - sin(lon)*sin(az)*cos(fpa) + cos(lon)*cos(latgc)*sin(fpa) ); 
            //vy = magv* ( -sin(lon)*sin(latgc)*cos(az)*cos(fpa) + cos(lon)*sin(az)*cos(fpa) + sin(lon)*cos(latgc)*sin(fpa) );  
            //vz = magv* (sin(latgc) * sin(fpa) + cos(latgc)*cos(az)*cos(fpa) );
            //// correct:
            //ve1 = magv* ( -cos(rtasc)*sin(decl)*cos(az)*cos(fpa) - sin(rtasc)*sin(az)*cos(fpa) + cos(rtasc)*cos(decl)*sin(fpa) ); // m/s
            // ve2 = magv* (-sin(rtasc) * sin(decl) * cos(az) * cos(fpa) + cos(rtasc) * sin(az) * cos(fpa) + sin(rtasc) * cos(decl) * sin(fpa));
            //ve3 = magv* (sin(decl) * sin(fpa) + cos(decl)*cos(az)*cos(fpa) );

            strbuild.AppendLine("==================== do the sensitivity tests \n");

            strbuild.AppendLine("1.  Cartesian Covariance \n");
            printcov(cartcov, "ct", 'm', anom, out strout);
            strbuild.AppendLine(strout);

            strbuild.AppendLine("3.  Equinoctial Covariance from Cartesian #1 above (" + anom + ") ------------------- \n");
            AstroLibr.covct2eq(cartcov, cartstate, anom, fr, out eqcovmeana, out tmct2eq);

            printcov(eqcovmeana, "eq", 'm', anom, out strout);
            strbuild.AppendLine(strout);

            strbuild.AppendLine("  Cartesian Covariance from Classical #3 above \n");
            AstroLibr.coveq2ct(eqcovmeana, eqstate, anom, fr, out cartcovmeanarev, out tmeq2ct);

            printcov(cartcovmeanarev, "ct", 'm', anom, out strout);
            strbuild.AppendLine(strout);
            strbuild.AppendLine("\n");

            printdiff(" cartcov - cartcov" + anom + "rev \n", cartcov, cartcovmeanarev, out strout);
            strbuild.AppendLine(strout);
        }  // testcoveq2clmeann


        public void testcovct2fl(string anomflt)
        {
            int year, mon, day, hr, minute, timezone, dat, terms;
            double second, dut1, ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac;
            double p, a, n, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper;
            double af, ag, chi, psi, meanlonNu, meanlonM;
            double latgc, lon, rtasc, decl, fpa, lod, xp, yp, ddpsi, ddeps,az, magr, magv;
            Int16 fr;
            double[,] tm = new double[,] { { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 } };
            string anom = "meann";
            double[] cartstate = new double[6];
            double[] classstate = new double[6];
            double[] eqstate = new double[6];
            double[] fltstate = new double[6];
            double[] recef = new double[3];
            double[] vecef = new double[3];
            double[] avec = new double[3];

            double[,] classcovmeana = new double[6, 6];
            double[,] cartcovmeanarev = new double[6, 6];
            double[,] fltcovmeana = new double[6, 6];
            double[,] tmct2cl = new double[6, 6];
            double[,] tmcl2ct = new double[6, 6];
            double[,] tmct2fl = new double[6, 6];
            double[,] tmfl2ct = new double[6, 6];
            string strout;

            double[] reci = new double[3] { -605.79221660, -5870.22951108, 3493.05319896 };
            double[] veci = new double[3] { -1.56825429, -3.70234891, -6.47948395 };
            double[] aeci = new double[3] { 0.001, 0.002, 0.003 };

            // StringBuilder strbuild = new StringBuilder();
            // strbuild.Clear();

            EOPSPWLib.iau80Class iau80arr;
            string fileLoc;
            fileLoc = @"D:\Codes\LIBRARY\DataLib\nut80.dat";
            EOPSPWLibr.iau80in(fileLoc, out iau80arr);

            year = 2000;
            mon = 12;
            day = 15;
            hr = 16;
            minute = 58;
            second = 50.208;
            dut1 = 0.10597;
            dat = 32;
            timezone = 0;
            xp = 0.0;
            yp = 0.0;
            lod = 0.0;
            terms = 2;
            timezone = 0;
            ddpsi = 0.0;
            ddeps = 0.0;

            MathTimeLibr.convtime(year, mon, day, hr, minute, second, timezone, dut1, dat,
                out ut1, out tut1, out jdut1, out jdut1frac, out utc, out tai,
                out tt, out ttt, out jdtt, out jdttfrac,
                out tdb, out ttdb, out jdtdb, out jdtdbfrac);

            // ---convert the eci state into the various other state formats(classical, equinoctial, etc)
            double[,] cartcov = new double[,] {
            { 100.0, 1.0e-2, 1.0e-2, 1.0e-4, 1.0e-4, 1.0e-4 },
            {1.0e-2, 100.0,  1.0e-2, 1.0e-4,   1.0e-4,   1.0e-4},
            {1.0e-2, 1.0e-2, 100.0,  1.0e-4,   1.0e-4,   1.0e-4},
            {1.0e-4, 1.0e-4, 1.0e-4, 0.0001,   1.0e-6,   1.0e-6},
            {1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   0.0001,   1.0e-6},
            {1.0e-4, 1.0e-4, 1.0e-4, 1.0e-6,   1.0e-6,   0.0001} };
            cartstate = new double[] { reci[0], reci[1], reci[2], veci[0], veci[1], veci[2] };  // in km

            // --------convert to a classical orbit state
            AstroLibr.rv2coe(reci, veci,
                out p, out a, out ecc, out incl, out raan, out argp, out nu, out m, out arglat, out truelon, out lonper);
            classstate[0] = a;   // km
            classstate[1] = ecc;
            classstate[2] = incl;
            classstate[3] = raan;
            classstate[4] = argp;
            if (anom.Contains("mean")) // meann or meana
                classstate[5] = m;
            else  // truea or truen
                classstate[5] = nu;

            // -------- convert to an equinoctial orbit state
            AstroLibr.rv2eq(reci, veci, out a, out n, out af, out ag, out chi, out psi, out meanlonM, out meanlonNu, out fr);
            if (anom.Equals("meana") || anom.Equals("truea"))
                eqstate[0] = a;  // km
            else // meann or truen
                eqstate[0] = n;
            eqstate[1] = af;
            eqstate[2] = ag;
            eqstate[3] = chi;
            eqstate[4] = psi;
            if (anom.Contains("mean")) //  meana or meann
                eqstate[5] = meanlonM;
            else // truea or truen
                eqstate[5] = meanlonNu;

            // --------convert to a flight orbit state
            AstroLibr.rv2flt(reci, veci, jdtt, jdttfrac, jdut1, lod, xp, yp, terms, ddpsi, ddeps,
                iau80arr, out lon, out latgc, out rtasc, out decl, out fpa, out az, out magr, out magv);
            if (anomflt.Equals("radec"))
            {
                fltstate[0] = rtasc;
                fltstate[1] = decl;
            }
            else
            if (anomflt.Equals("latlon"))
            {
                fltstate[0] = lon;
                fltstate[1] = latgc;
            }
            fltstate[2] = fpa;
            fltstate[3] = az;
            fltstate[4] = magr;  // km
            fltstate[5] = magv;

            // test position and velocity going back
            avec = new double[] { 0.0, 0.0, 0.0 };

            AstroLibr.eci_ecef(ref reci, ref veci, MathTimeLib.Edirection.eto, ref recef, ref vecef,
                  iau80arr, jdtt, jdttfrac, jdut1, lod, xp, yp, ddpsi, ddeps);
            //vx = magv* ( -cos(lon)*sin(latgc)*cos(az)*cos(fpa) - sin(lon)*sin(az)*cos(fpa) + cos(lon)*cos(latgc)*sin(fpa) ); 
            //vy = magv* ( -sin(lon)*sin(latgc)*cos(az)*cos(fpa) + cos(lon)*sin(az)*cos(fpa) + sin(lon)*cos(latgc)*sin(fpa) );  
            //vz = magv* (sin(latgc) * sin(fpa) + cos(latgc)*cos(az)*cos(fpa) );
            //// correct:
            //ve1 = magv* ( -cos(rtasc)*sin(decl)*cos(az)*cos(fpa) - sin(rtasc)*sin(az)*cos(fpa) + cos(rtasc)*cos(decl)*sin(fpa) ); // m/s
            // ve2 = magv* (-sin(rtasc) * sin(decl) * cos(az) * cos(fpa) + cos(rtasc) * sin(az) * cos(fpa) + sin(rtasc) * cos(decl) * sin(fpa));
            //ve3 = magv* (sin(decl) * sin(fpa) + cos(decl)*cos(az)*cos(fpa) );

            strbuild.AppendLine("==================== do the sensitivity tests \n");

            strbuild.AppendLine("1.  Cartesian Covariance \n");
            printcov(cartcov, "ct", 'm', anomflt, out strout);
            strbuild.AppendLine(strout);

            strbuild.AppendLine("7.  Flight Covariance from Cartesian #1 above (" + anomflt + ") ------------------- \n");
            AstroLibr.covct2fl(cartcov, cartstate, anomflt, jdtt, jdttfrac, jdut1, lod, xp, yp, 2, ddpsi, ddeps,
                iau80arr, out fltcovmeana, out tmct2fl);

            if (anomflt.Equals("latlon"))
                printcov(fltcovmeana, "fl", 'm', anomflt, out strout);
            else
                printcov(fltcovmeana, "sp", 'm', anomflt, out strout);

            strbuild.AppendLine(strout);

            strbuild.AppendLine("  Cartesian Covariance from Flight #7 above \n");
            AstroLibr.covfl2ct(fltcovmeana, fltstate, anomflt, jdtt, jdttfrac, jdut1,
                lod, xp, yp, 2, ddpsi, ddeps, iau80arr, out cartcovmeanarev, out tmfl2ct);

            printcov(cartcovmeanarev, "ct", 'm', anomflt, out strout);
            strbuild.AppendLine(strout);
            strbuild.AppendLine("\n");

            printdiff(" cartcov - cartcov" + anomflt + "rev \n", cartcov, cartcovmeanarev, out strout);
            strbuild.AppendLine(strout);
        }  // testcovct2fl


        public void testlegacc()
        {
            Int32 orderSize = 1500;
            double[] aPertGot = new double[3];
            double[] norm1 = new double[orderSize + 2];
            double[] norm2 = new double[orderSize + 2];
            double[] norm10 = new double[orderSize + 2];
            double[] norm11 = new double[orderSize + 2];
            double[] normn10 = new double[orderSize + 2];
            double[,] normn1 = new double[orderSize + 2, orderSize + 2];
            double[,] norm1m = new double[orderSize + 2, orderSize + 2];
            double[,] norm2m = new double[orderSize + 2, orderSize + 2];
            double[,,] normarr = new double[orderSize + 2, orderSize + 2, 7];

            double latgc, lon;
            Int32  order, orderxx, startKtr;
            double[,] LegGottN; // Gottlieb

            StringBuilder strbuildall = new StringBuilder();
            AstroLib.gravityConst gravData;

            order = orderSize;
            double rad = 180.0 / Math.PI;              // deg to rad

            string fname = "D:/Dataorig/Gravity/EGM2008_to2190_TideFree.txt";
            startKtr = 0;
            char normalized = 'y';  // if file has normalized coefficients
            AstroLibr.readGravityField(fname, normalized, startKtr, out orderxx, out gravData);

            // test case
            double[] recef = new double[] { -2110.2895393, -5511.9160245, 3491.9133986 };
            latgc = 30.610308418 / rad;
            lon = -110.949820402 / rad;

            // Gottlieb acceleration
            strbuildall.AppendLine("Gottlieb acceleration ");
            
            // this part can be done one time
            AstroLibr.LegPolyGottN(order, out norm1, out norm2, out norm11, out normn10, out norm1m,
                out norm2m, out normn1);

            AstroLibr.FullGeopGott(recef, order, gravData, norm1, norm2, norm11, normn10,
                   norm1m, norm2m, normn1, out LegGottN, out aPertGot);

            strbuildall.AppendLine("0 " + LegGottN[0, 0].ToString());
            strbuildall.AppendLine("1 " + LegGottN[1, 0].ToString() + " " + LegGottN[1, 1].ToString());
            strbuildall.AppendLine("2 " + LegGottN[2, 0].ToString() + " " + LegGottN[2, 1].ToString()
                + " " + LegGottN[2, 2].ToString());
            strbuildall.AppendLine("3 " + LegGottN[3, 0].ToString() + " " + LegGottN[3, 1].ToString()
                + " " + LegGottN[3, 2].ToString() + " " + LegGottN[3, 3].ToString());
            strbuildall.AppendLine("acceleration " + order + " " + aPertGot[0] + " " + aPertGot[1]
                + " " + aPertGot[2]);

            string directory = @"D:\Codes\LIBRARY\cs\TestAll\";
            File.WriteAllText(directory + "testdan.out", strbuildall.ToString());

        }  // testlegacc


        public void runtests(int testnum)
        {
            Int32 opt, optstart, optstop;

            if (testnum == -10)
            {
                optstart = 1;
                optstop = 103; // 103
            }
            else
            {
                optstart = testnum;
                optstop = testnum;
            }

            string directory = @"D:\Codes\LIBRARY\cs\TestAll\";

            for (opt = optstart; opt <= optstop; opt++)  //102
            {
                strbuild.AppendLine("\n\n=================================== Case " + opt.ToString() + " =======================================");
                this.opsStatus.Text = "Status:  on case " + opt.ToString();
                Refresh();
                switch (opt)
                {
                    case 1:
                        testvecouter();
                        break;
                    case 2:
                        testmatadd();
                        break;
                    case 3:
                        testmatsub();
                        break;
                    case 4:
                        testmatmult();
                        break;
                    case 5:
                        testmattrans();
                        break;
                    case 6:
                        testmattransx();
                        break;
                    case 7:
                        testmatinverse();
                        break;
                    case 8:
                        testdeterminant();
                        break;
                    case 9:
                        testcholesky();
                        break;
                    case 10:
                        testposvelcov2pts();
                        break;
                    case 11:
                        testposcov2pts();
                        break;
                    case 12:
                        testremakecovpv();
                        break;
                    case 13:
                        testremakecovp();
                        break;
                    case 14:
                        testmatequal();
                        break;
                    case 15:
                        testmatscale();
                        break;
                    case 16:
                        testnorm();
                        break;
                    case 17:
                        testmag();
                        break;
                    case 18:
                        testcross();
                        break;
                    case 19:
                        testdot();
                        break;
                    case 20:
                        testangle();
                        break;
                    case 21:
                        testasinh();
                        break;
                    case 22:
                        testcot();
                        break;
                    case 23:
                        testacosh();
                        break;
                    case 24:
                        testaddvec();
                        break;
                    case 25:
                        testPercentile();
                        break;
                    case 26:
                        testrot1();
                        break;
                    case 27:
                        testrot2();
                        break;
                    case 28:
                        testrot3();
                        break;
                    case 29:
                        testfactorial();
                        break;
                    case 30:
                        testcubicspl();
                        break;
                    case 31:
                        testcubic();
                        break;
                    case 32:
                        testcubicinterp();
                        break;
                    case 33:
                        testquadratic();
                        break;
                    case 34:
                        testconvertMonth();
                        break;
                    case 35:
                        testjday();
                        break;
                    case 36:
                        testdays2mdhms();
                        break;
                    case 37:
                        testinvjday();
                        break;
                    case 38:
                        // tests eop, spw, and fk5 iau80
                        testiau80in();
                        break;
                    case 39:
                        testfundarg();
                        break;
                    case 40:
                        testprecess();
                        break;
                    case 41:
                        testnutation();
                        break;
                    case 42:
                        testnutationqmod();
                        break;
                    case 43:
                        testsidereal();
                        break;
                    case 44:
                        testpolarm();
                        break;
                    case 45:
                        testgstime();
                        break;
                    case 46:
                        testlstime();
                        break;
                    case 47:
                        testhms_sec();
                        break;
                    case 48:
                        testhms_ut();
                        break;
                    case 49:
                        testhms_rad();
                        break;
                    case 50:
                        testdms_rad();
                        break;
                    case 51:
                        testeci_ecef();
                        break;
                    case 52:
                        testtod2ecef();
                        break;
                    case 53:
                        testteme_ecef();
                        break;
                    case 54:
                        testteme_eci();
                        break;
                    case 55:
                        testqmod2ecef();
                        break;
                    case 56:
                        testcsm2efg();
                        break;
                    case 57:
                        testrv_elatlon();
                        break;
                    case 58:
                        testrv2radec();
                        break;
                    case 59:
                        testrv_razel();
                        break;
                    case 60:
                        testrv_tradec();
                        break;
                    case 61:
                        testrvsez_razel();
                        break;
                    case 62:
                        testrv2rsw();
                        break;
                    case 63:
                        testrv2pqw();
                        break;
                    case 64:
                        testrv2coe();
                        break;
                    case 65:
                        testcoe2rv();
                        break;
                    case 66:
                        testlon2nu();
                        break;
                    case 67:
                        testnewtonmx();
                        break;
                    case 68:
                        testnewtonm();
                        break;
                    case 69:
                        testnewtonnu();
                        break;
                    case 70:
                        testfindc2c3();
                        break;
                    case 71:
                        testfindfandg();
                        break;
                    case 72:
                        testcheckhitearth();
                        testcheckhitearthc();
                        break;
                    case 73:
                        testgibbs();
                        testhgibbs();
                        break;
                    case 74:
                        // in sln directory, testall-Angles.out
                        testangles();
                        // D:\faabook\current\excel\testgeo.out for ch9 plot
                        testgeo();
                        break;
                    case 75:
                        testlambertumins();
                        testlambertminT();
                        break;
                    case 76:
                        testlambhodograph();
                        break;
                    case 77:
                        testlambertbattin();
                        break;
                    case 78:
                        testeq2rv();
                        break;
                    case 79:
                        testrv2eq();
                        break;
                    case 80:
                        string directoryx = @"d:\codes\library\matlab\";
                        // book example, simple
                        testlambertuniv();
                        strbuild.AppendLine("lambert envelope test case results written to " + directoryx + "testall.out ");

                        // old approach? yes...
                        testAllMoving();
                        strbuild.AppendLine("lambert all moving test case results written to " + directoryx + "tlambertAllx.out ");
                        strbuild.AppendLine("lambert all moving test case results written to " + directoryx + "tlamb3dx.out ");

                        // envelope testing
                        testAll();
                        strbuild.AppendLine("lambert envelope test case results written to " + directoryx + "testall.out ");

                        // known problem cases testall-lambertknown.out in sln directory
                        testknowncases();
                        strbuild.AppendLine("lambert known test case results written to " + directory + "testall-lambertknown.out ");
                        break;
                    case 81:
                        testradecgeo2azel();
                        break;
                    case 82:
                        testijk2ll();
                        break;
                    case 83:
                        testgd2gc();
                        break;
                    case 84:
                        testsite();
                        break;
                    case 85:
                        testsun();
                        break;
                    case 86:
                        testmoon();
                        break;
                    case 87:
                        testkepler();
                        break;
                    case 88:
                        // mean
                        testcovct2clmean();
                        break;
                    case 89:
                        // true
                        testcovct2cltrue();
                        break;
                    case 90:
                        testhill();
                        break;
                    case 91:
                        break;
                    case 92:
                        testcovct2rsw();
                        testcovct2ntw();
                        testcovcl2eq("truea");
                        testcovcl2eq("truen");
                        testcovcl2eq("meana");
                        testcovcl2eq("meann");
                        testcovct2eq("truea");
                        testcovct2eq("truen");
                        testcovct2eq("meana");
                        testcovct2eq("meann");
                        testcovct2fl("latlon");
                        testcovct2fl("radec");
                        break;
                    case 93:
                        break;
                    case 94:
                        break;
                    case 95:
                        break;
                    case 96:
                        testcovct2rsw();
                        break;
                    case 97:
                        testcovct2ntw();
                        break;
                    case 98:
                        testsunmoonjpl();
                        break;
                    case 99:
                        testkp2ap();
                        break;
                    case 100:
                        testproporbit();
                        string directoryy = @"d:\codes\library\matlab\";
                        strbuild.AppendLine("testproporbit (legendre) results written to " + directoryy + "legpoly.out ");
                        strbuild.AppendLine("testproporbit (legendre) results written to " + directoryy + "legendreAcc.out ");
                        break;
                    case 101:
                        //testsemianaly();
                        break;
                    case 102:
                        testazel2radec();
                        break;
                    case 103:
                        testlegacc();
                        break;
                }

            } // for

            // write data out
            File.WriteAllText(directory + "testall.out", strbuild.ToString());
        }  // runtests



        // -----------------------------------------------------------------
        //
        // this file tests all the various functions.
        //
        // companion code for
        // fundamentals of astrodyanmics and applications
        // 2019
        // by david vallado
        //
        // (w)719-573-2600, email dvallado@agi.com, davallado@gmail.com
        //
        // *****************************************************************
        //
        // current :
        // 24 jan 19  david vallado
        // original baseline
        //
        // *****************************************************************

        private void button4_Click_1(object sender, EventArgs e)
        {
            double rad = 180.0 / Math.PI;
            double rad2 = rad * rad;


            this.opsStatus.Text = "Status: Find Test Solutions ";
            Refresh();

            // -----------------------------------------------------------------
            // 74 angles-only tests
            // 80 lambert tests
            // 94 covariance tests
            // 101 semianalytical (not built yet)
            //

            runtests(-10);  // run all

            this.opsStatus.Text = "Done";
            Refresh();
            int pauseTime = 1000;
            System.Threading.Thread.Sleep(pauseTime);
            this.opsStatus.Text = "Status: ";

        }


        // setup Lambert test page
        void button1_Click(object sender, EventArgs e)
        {

            frmLambert form = new frmLambert();
            form.Show();

        }

        private void button2_Click(object sender, EventArgs e)
        {
            this.opsStatus.Text = "Status: Running test ";
            Refresh();

            runtests(Convert.ToInt32(this.testNum.Text));  // run specific

            this.opsStatus.Text = "Status: Done ";
            Refresh();
        }

        // setup Gauss test page
        private void button3_Click(object sender, EventArgs e)
        {

            frmGauss form = new frmGauss();
            form.Show();

        }

        private void time2_TextChanged(object sender, EventArgs e)
        {
            double jd, jdf;
            int year, month, day, h, m, dayofyr;
            double s;
            
            MathTimeLibr.STKtime2JD(this.time2.Text, out jd, out jdf);
            MathTimeLibr.invjday(jd, jdf, out year, out month, out day, out h, out m, out s);
            MathTimeLibr.findDays(year, month, day, h, m, s, out dayofyr);
            this.doy.Text = dayofyr.ToString();
        }

        //private void textChangedEventHandler(object sender, TextChangedEventArgs args)
        //{
        //    // Omitted Code: Insert code that does something whenever
        //    // the text changes...
        //} // end textChangedEventHandler

        private void rtasc_TextChanged(object sender, EventArgs e)
        {
            int h, m;
            double hms, s, rad;
            rad = 180.0 / Math.PI;

            h = m = 0;
            s = 0.0;

            hms = Convert.ToDouble(this.rtasc.Text) / rad;  // get in rads first, not hrs
            MathTimeLibr.hms_rad(ref h, ref m, ref s, MathTimeLib.Edirection.efrom, ref hms);

            this.hms.Text = h + " " + m + " " + s;
        }

        private void decl_TextChanged(object sender, EventArgs e)
        {
            int d, m;
            double dms, s, rad;
            rad = 180.0 / Math.PI;

            d = m = 0;
            s = 0.0;

            dms = Convert.ToDouble(this.decl.Text) / rad;
            MathTimeLibr.dms_rad(ref d, ref m, ref s, MathTimeLib.Edirection.efrom, ref dms);

            this.dms.Text = d + " " + m + " " + s;
        }

        private void hms_TextChanged(object sender, EventArgs e)
        {
            int h, m;
            double hms, s;

            hms = 0.0;
            string line;
            line = this.hms.Text;
            line.Replace(@"\s+", " ");
            string[] linesplt = line.Split(' ');
            h = Convert.ToInt32(linesplt[0]);
            m = Convert.ToInt32(linesplt[1]);
            s = Convert.ToDouble(linesplt[2]);

            MathTimeLibr.hms_rad(ref h, ref m, ref s, MathTimeLib.Edirection.eto, ref hms);

            this.rtasc.Text = hms.ToString();
        }

        private void dms_TextChanged(object sender, EventArgs e)
        {
            int d, m;
            double dms, s;

            dms = 0.0;
            string line;
            line = this.dms.Text;
            line.Replace(@"\s+", " ");
            string[] linesplt = line.Split(' ');
            d = Convert.ToInt32(linesplt[0]);
            m = Convert.ToInt32(linesplt[1]);
            s = Convert.ToDouble(linesplt[2]);

            MathTimeLibr.dms_rad(ref d, ref m, ref s, MathTimeLib.Edirection.eto, ref dms);

            this.decl.Text = dms.ToString();
        }
    }
}
