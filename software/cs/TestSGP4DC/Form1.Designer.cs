﻿namespace TestSGP4DC
{
    partial class Form1
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.button4 = new System.Windows.Forms.Button();
            this.label2 = new System.Windows.Forms.Label();
            this.label1 = new System.Windows.Forms.Label();
            this.label12 = new System.Windows.Forms.Label();
            this.OFLoc = new System.Windows.Forms.TextBox();
            this.label3 = new System.Windows.Forms.Label();
            this.label4 = new System.Windows.Forms.Label();
            this.label5 = new System.Windows.Forms.Label();
            this.label6 = new System.Windows.Forms.Label();
            this.obs2skip = new System.Windows.Forms.TextBox();
            this.label7 = new System.Windows.Forms.Label();
            this.obs2read = new System.Windows.Forms.TextBox();
            this.cbTypeODRun = new System.Windows.Forms.ComboBox();
            this.cbmatinv = new System.Windows.Forms.ComboBox();
            this.cbproptype = new System.Windows.Forms.ComboBox();
            this.SuspendLayout();
            // 
            // button4
            // 
            this.button4.BackColor = System.Drawing.Color.FromArgb(((int)(((byte)(255)))), ((int)(((byte)(192)))), ((int)(((byte)(255)))));
            this.button4.Location = new System.Drawing.Point(514, 647);
            this.button4.Margin = new System.Windows.Forms.Padding(7, 6, 7, 6);
            this.button4.Name = "button4";
            this.button4.Size = new System.Drawing.Size(307, 44);
            this.button4.TabIndex = 294;
            this.button4.Text = "Run Tests";
            this.button4.UseVisualStyleBackColor = false;
            this.button4.Click += new System.EventHandler(this.button4_Click);
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Font = new System.Drawing.Font("Cambria", 12F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.label2.Location = new System.Drawing.Point(595, 70);
            this.label2.Margin = new System.Windows.Forms.Padding(4, 0, 4, 0);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(273, 37);
            this.label2.TabIndex = 293;
            this.label2.Text = "Test All Summary";
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Font = new System.Drawing.Font("Microsoft Sans Serif", 6F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.label1.Location = new System.Drawing.Point(157, 120);
            this.label1.Margin = new System.Windows.Forms.Padding(4, 0, 4, 0);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(243, 20);
            this.label1.TabIndex = 292;
            this.label1.Text = "Program to generate debris fields";
            // 
            // label12
            // 
            this.label12.AutoSize = true;
            this.label12.Location = new System.Drawing.Point(178, 170);
            this.label12.Margin = new System.Windows.Forms.Padding(4, 0, 4, 0);
            this.label12.Name = "label12";
            this.label12.Size = new System.Drawing.Size(191, 25);
            this.label12.TabIndex = 291;
            this.label12.Text = "Output file location";
            // 
            // OFLoc
            // 
            this.OFLoc.Location = new System.Drawing.Point(183, 200);
            this.OFLoc.Margin = new System.Windows.Forms.Padding(4);
            this.OFLoc.Name = "OFLoc";
            this.OFLoc.Size = new System.Drawing.Size(1130, 31);
            this.OFLoc.TabIndex = 290;
            this.OFLoc.Text = "D:\\POE\\TestSLR\\Automation\\CurrentScenarioLEOGrace\\";
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(473, 298);
            this.label3.Margin = new System.Windows.Forms.Padding(4, 0, 4, 0);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(121, 25);
            this.label3.TabIndex = 296;
            this.label3.Text = "Type of run";
            // 
            // label4
            // 
            this.label4.AutoSize = true;
            this.label4.Location = new System.Drawing.Point(473, 366);
            this.label4.Margin = new System.Windows.Forms.Padding(4, 0, 4, 0);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(328, 50);
            this.label4.TabIndex = 298;
            this.label4.Text = "Matrix Inversion";
            // 
            // label5
            // 
            this.label5.AutoSize = true;
            this.label5.Location = new System.Drawing.Point(473, 425);
            this.label5.Margin = new System.Windows.Forms.Padding(4, 0, 4, 0);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(350, 50);
            this.label5.TabIndex = 300;
            this.label5.Text = "Propagation type";
            // 
            // label6
            // 
            this.label6.AutoSize = true;
            this.label6.Location = new System.Drawing.Point(473, 487);
            this.label6.Margin = new System.Windows.Forms.Padding(4, 0, 4, 0);
            this.label6.Name = "label6";
            this.label6.Size = new System.Drawing.Size(240, 50);
            this.label6.TabIndex = 302;
            this.label6.Text = "Obs to skip";
            // 
            // obs2skip
            // 
            this.obs2skip.Location = new System.Drawing.Point(183, 481);
            this.obs2skip.Margin = new System.Windows.Forms.Padding(4);
            this.obs2skip.Name = "obs2skip";
            this.obs2skip.Size = new System.Drawing.Size(92, 31);
            this.obs2skip.TabIndex = 301;
            this.obs2skip.Text = "1";
            // 
            // label7
            // 
            this.label7.AutoSize = true;
            this.label7.Location = new System.Drawing.Point(473, 554);
            this.label7.Margin = new System.Windows.Forms.Padding(4, 0, 4, 0);
            this.label7.Name = "label7";
            this.label7.Size = new System.Drawing.Size(348, 50);
            this.label7.TabIndex = 304;
            this.label7.Text = "Total obs to read";
            // 
            // obs2read
            // 
            this.obs2read.Location = new System.Drawing.Point(183, 548);
            this.obs2read.Margin = new System.Windows.Forms.Padding(4);
            this.obs2read.Name = "obs2read";
            this.obs2read.Size = new System.Drawing.Size(92, 31);
            this.obs2read.TabIndex = 303;
            this.obs2read.Text = "20";
            // 
            // cbTypeODRun
            // 
            this.cbTypeODRun.FormattingEnabled = true;
            this.cbTypeODRun.Items.AddRange(new object[] {
            "Range - az - el",
            "TLE",
            "Ephemeris"});
            this.cbTypeODRun.Location = new System.Drawing.Point(183, 290);
            this.cbTypeODRun.Margin = new System.Windows.Forms.Padding(7, 6, 7, 6);
            this.cbTypeODRun.Name = "cbTypeODRun";
            this.cbTypeODRun.Size = new System.Drawing.Size(238, 33);
            this.cbTypeODRun.TabIndex = 305;
            this.cbTypeODRun.Text = "Select...";
            // 
            // cbmatinv
            // 
            this.cbmatinv.FormattingEnabled = true;
            this.cbmatinv.Items.AddRange(new object[] {
            "SVD",
            "LuBkSub"});
            this.cbmatinv.Location = new System.Drawing.Point(183, 353);
            this.cbmatinv.Margin = new System.Windows.Forms.Padding(7, 6, 7, 6);
            this.cbmatinv.Name = "cbmatinv";
            this.cbmatinv.Size = new System.Drawing.Size(238, 33);
            this.cbmatinv.TabIndex = 306;
            this.cbmatinv.Text = "Select...";
            // 
            // cbproptype
            // 
            this.cbproptype.FormattingEnabled = true;
            this.cbproptype.Items.AddRange(new object[] {
            "Two-body",
            "SGP4",
            "Numerical",
            "Semianalytical"});
            this.cbproptype.Location = new System.Drawing.Point(183, 422);
            this.cbproptype.Margin = new System.Windows.Forms.Padding(7, 6, 7, 6);
            this.cbproptype.Name = "cbproptype";
            this.cbproptype.Size = new System.Drawing.Size(238, 33);
            this.cbproptype.TabIndex = 307;
            this.cbproptype.Text = "Select...";
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(12F, 25F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(1408, 835);
            this.Controls.Add(this.cbproptype);
            this.Controls.Add(this.cbmatinv);
            this.Controls.Add(this.cbTypeODRun);
            this.Controls.Add(this.label7);
            this.Controls.Add(this.obs2read);
            this.Controls.Add(this.label6);
            this.Controls.Add(this.obs2skip);
            this.Controls.Add(this.label5);
            this.Controls.Add(this.label4);
            this.Controls.Add(this.label3);
            this.Controls.Add(this.button4);
            this.Controls.Add(this.label2);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.label12);
            this.Controls.Add(this.OFLoc);
            this.Name = "Form1";
            this.Text = "Form1";
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Button button4;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.Label label12;
        private System.Windows.Forms.TextBox OFLoc;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.Label label4;
        private System.Windows.Forms.Label label5;
        private System.Windows.Forms.Label label6;
        private System.Windows.Forms.TextBox obs2skip;
        private System.Windows.Forms.Label label7;
        private System.Windows.Forms.TextBox obs2read;
        private System.Windows.Forms.ComboBox cbTypeODRun;
        private System.Windows.Forms.ComboBox cbmatinv;
        private System.Windows.Forms.ComboBox cbproptype;
    }
}

