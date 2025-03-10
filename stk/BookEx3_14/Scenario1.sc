stk.v.12.0
WrittenBy    STK_v12.0.0
BEGIN Scenario
    Name		 Scenario1

    BEGIN Epoch

        Epoch		 6 Apr 2004 07:51:28.386000000
        SmartEpoch		
        BEGIN EVENT
            Epoch		 6 Apr 2004 07:51:28.386000000
            EpochState		 Explicit
        END EVENT


    END Epoch

    BEGIN Interval

        Start		 6 Apr 2004 07:51:28.386009000
        Stop		 7 Apr 2004 19:00:00.000000000
        SmartInterval		
        BEGIN EVENTINTERVAL
            BEGIN Interval
                Start		 6 Apr 2004 07:51:28.386009000
                Stop		 7 Apr 2004 19:00:00.000000000
            END Interval
            IntervalState		 Explicit
        END EVENTINTERVAL

        EpochUsesAnalStart		 No
        AnimStartUsesAnalStart		 Yes
        AnimStopUsesAnalStop		 Yes

    END Interval

    BEGIN EOPFile

        InheritEOPSource		 No
        EOPFilename		 EOP-All-v1.1.txt

    END EOPFile

    BEGIN GlobalPrefs
        SatelliteNoOrbWarning		 No
        MissilePerigeeWarning		 No
        MissileStopTimeWarning		 No
        AircraftWGS84Warning		 Always
    END GlobalPrefs

    BEGIN CentralBody

        PrimaryBody		 Earth

    END CentralBody

    BEGIN CentralBodyTerrain

        BEGIN CentralBody
            Name		 Earth
            UseTerrainCache		 Yes
            TotalCacheSize		 402653184

            BEGIN StreamingTerrain
                UseCurrentStreamingTerrainServer		 No
                StreamingTerrainServerName		 assets.agi.com/stk-terrain/
                StreamingTerrainAzimuthElevationMaskEnabled		 Yes
                StreamingTerrainObscurationEnabled		 Yes
                StreamingTerrainCoverageGridObscurationEnabled		 Yes
            END StreamingTerrain
        END CentralBody

    END CentralBodyTerrain

    BEGIN StarCollection

        Name		 Hipparcos 2 Mag 6

    END StarCollection

    BEGIN ScenarioLicenses
        Module		 AMMProv12.0
        Module		 AMMv12.0
        Module		 ASTGv12.0
        Module		 AnalysisWBv12.0
        Module		 CATv12.0
        Module		 COVv12.0
        Module		 Commv12.0
        Module		 DISv12.0
        Module		 EOIRv12.0
        Module		 RT3Clientv12.0
        Module		 Radarv12.0
        Module		 SEETv12.0
        Module		 STKIntegrationv12.0
        Module		 STKParallelComputingv12.0
        Module		 STKProfessionalv12.0
        Module		 STKTIMv12.0
        Module		 STKv12.0
        Module		 SatProv12.0
        Module		 TIREMv12.0
        Module		 UPropv12.0
    END ScenarioLicenses

    BEGIN Extensions

        BEGIN ClsApp
            RangeConstraint		 5000
            ApoPeriPad		 30000
            OrbitPathPad		 30000
            TimeDistPad		 30000
            OutOfDate		 2592000
            MaxApoPeriStep		 900
            ApoPeriAngle		 0.785
            UseApogeePerigeeFilter		 Yes
            UsePathFilter		 No
            UseTimeFilter		 No
            UseOutOfDate		 Yes
            CreateSats		 No
            MaxSatsToCreate		 500
            UseModelScale		 No
            ModelScale		 4
            UseCrossRefDb		 Yes
            CollisionDB		 stkAllTLE.tce
            CollisionCrossRefDB		 stkAllTLE.sd
            ShowLine		 Yes
            AnimHighlight		 Yes
            StaticHighlight		 Yes
            UseLaunchWindow		 No
            LaunchWindowUseEntireTraj		 Yes
            LaunchWindowTrajMETStart		 0
            LaunchWindowTrajMETStop		 900
            LaunchWindowStart		 0
            LaunchWindowStop		 0
            LaunchMETOffset		 0
            LaunchWindowUseSecEphem		 No 
            LaunchWindowUseScenFolderForSecEphem		 Yes
            LaunchWindowUsePrimEphem		 No 
            LaunchWindowUseScenFolderForPrimEphem		 Yes
            LaunchWindowIntervalPtr		
            BEGIN EVENTINTERVAL
                BEGIN Interval
                    Start		 6 Apr 2004 07:51:28.386009000
                    Stop		 7 Apr 2004 19:00:00.000000000
                END Interval
                IntervalState		 Explicit
            END EVENTINTERVAL

            LaunchWindowUsePrimMTO		 No 
            GroupLaunches		 No 
            LWTimeConvergence		 0.001
            LWRelValueConvergence		 1e-08
            LWTSRTimeConvergence		 0.0001
            LWTSRRelValueConvergence		 1e-10
            LaunchWindowStep		 300
            MaxTSRStep		  1.8000000000000000e+02
            MaxTSRRelMotion		  2.0000000000000000e+01
            UseLaunchArea		 No 
            LaunchAreaOrientation		 North
            LaunchAreaAzimuth		 0
            LaunchAreaXLimits		 -10000 10000
            LaunchAreaYLimits		 -10000 10000
            LaunchAreaNumXIntrPnts		 1
            LaunchAreaNumYIntrPnts		 1
            LaunchAreaAltReference		 Ellipsoid
            TargetSameStop		 No 
            SkipSurfaceMetric		 No 
            LWAreaTSRRelValueConvergence		 1e-10
            AreaLaunchWindowStep		 300
            AreaMaxTSRStep		  3.0000000000000000e+01
            AreaMaxTSRRelMotion		 1
            ShowLaunchArea		 No 
            ShowBlackoutTracks		 No 
            ShowClearedTracks		 No 
            UseObjectForClearedColor		 No 
            BlackoutColor		 #ff0000
            ClearedColor		 #ffffff
            ShowTracksSegments		 No 
            ShowMinRangeTracks		 No 
            MinRangeTrackTimeStep		 0.5
            UsePrimStepForTracks		 Yes
            GfxTracksTimeStep		 30
            GfxAreaNumXIntrPnts		 1
            GfxAreaNumYIntrPnts		 1
            CreateLaunchMTO		 No 
            CovarianceSigmaScale		 3
            CovarianceMode		 None 
        END ClsApp

        BEGIN Units
            DistanceUnit		 Kilometers
            TimeUnit		 Seconds
            DateFormat		 GregorianUTC
            AngleUnit		 Degrees
            MassUnit		 Kilograms
            PowerUnit		 dBW
            FrequencyUnit		 Gigahertz
            SmallDistanceUnit		 Meters
            LatitudeUnit		 Degrees
            LongitudeUnit		 Degrees
            DurationUnit		 Hr:Min:Sec
            Temperature		 Kelvin
            SmallTimeUnit		 Seconds
            RatioUnit		 Decibel
            RcsUnit		 Decibel
            DopplerVelocityUnit		 MetersperSecond
            SARTimeResProdUnit		 Meter-Second
            ForceUnit		 Newtons
            PressureUnit		 Pascals
            SpecificImpulseUnit		 Seconds
            PRFUnit		 Kilohertz
            BandwidthUnit		 Megahertz
            SmallVelocityUnit		 CentimetersperSecond
            Percent		 Percentage
            AviatorDistanceUnit		 NauticalMiles
            AviatorTimeUnit		 Hours
            AviatorAltitudeUnit		 Feet
            AviatorFuelQuantityUnit		 Pounds
            AviatorRunwayLengthUnit		 Kilofeet
            AviatorBearingAngleUnit		 Degrees
            AviatorAngleOfAttackUnit		 Degrees
            AviatorAttitudeAngleUnit		 Degrees
            AviatorGUnit		 StandardSeaLevelG
            SolidAngle		 Steradians
            AviatorTSFCUnit		 TSFCLbmHrLbf
            AviatorPSFCUnit		 PSFCLbmHrHp
            AviatorForceUnit		 Pounds
            AviatorPowerUnit		 Horsepower
            SpectralBandwidthUnit		 Hertz
            AviatorAltTimeUnit		 Minutes
            AviatorSmallTimeUnit		 Seconds
            AviatorEnergyUnit		 kilowatt-hours
            BitsUnit		 MegaBits
            MagneticFieldUnit		 nanoTesla
        END Units

        BEGIN ReportUnits
            DistanceUnit		 Kilometers
            TimeUnit		 Seconds
            DateFormat		 GregorianUTC
            AngleUnit		 Degrees
            MassUnit		 Kilograms
            PowerUnit		 dBW
            FrequencyUnit		 Gigahertz
            SmallDistanceUnit		 Meters
            LatitudeUnit		 Degrees
            LongitudeUnit		 Degrees
            DurationUnit		 Hr:Min:Sec
            Temperature		 Kelvin
            SmallTimeUnit		 Seconds
            RatioUnit		 Decibel
            RcsUnit		 Decibel
            DopplerVelocityUnit		 MetersperSecond
            SARTimeResProdUnit		 Meter-Second
            ForceUnit		 Newtons
            PressureUnit		 Pascals
            SpecificImpulseUnit		 Seconds
            PRFUnit		 Kilohertz
            BandwidthUnit		 Megahertz
            SmallVelocityUnit		 CentimetersperSecond
            Percent		 Percentage
            AviatorDistanceUnit		 NauticalMiles
            AviatorTimeUnit		 Hours
            AviatorAltitudeUnit		 Feet
            AviatorFuelQuantityUnit		 Pounds
            AviatorRunwayLengthUnit		 Kilofeet
            AviatorBearingAngleUnit		 Degrees
            AviatorAngleOfAttackUnit		 Degrees
            AviatorAttitudeAngleUnit		 Degrees
            AviatorGUnit		 StandardSeaLevelG
            SolidAngle		 Steradians
            AviatorTSFCUnit		 TSFCLbmHrLbf
            AviatorPSFCUnit		 PSFCLbmHrHp
            AviatorForceUnit		 Pounds
            AviatorPowerUnit		 Horsepower
            SpectralBandwidthUnit		 Hertz
            AviatorAltTimeUnit		 Minutes
            AviatorSmallTimeUnit		 Seconds
            AviatorEnergyUnit		 kilowatt-hours
            BitsUnit		 MegaBits
            MagneticFieldUnit		 nanoTesla
        END ReportUnits

        BEGIN ConnectReportUnits
            DistanceUnit		 Kilometers
            TimeUnit		 Seconds
            DateFormat		 GregorianUTC
            AngleUnit		 Degrees
            MassUnit		 Kilograms
            PowerUnit		 dBW
            FrequencyUnit		 Gigahertz
            SmallDistanceUnit		 Meters
            LatitudeUnit		 Degrees
            LongitudeUnit		 Degrees
            DurationUnit		 Hr:Min:Sec
            Temperature		 Kelvin
            SmallTimeUnit		 Seconds
            RatioUnit		 Decibel
            RcsUnit		 Decibel
            DopplerVelocityUnit		 MetersperSecond
            SARTimeResProdUnit		 Meter-Second
            ForceUnit		 Newtons
            PressureUnit		 Pascals
            SpecificImpulseUnit		 Seconds
            PRFUnit		 Kilohertz
            BandwidthUnit		 Megahertz
            SmallVelocityUnit		 CentimetersperSecond
            Percent		 Percentage
            AviatorDistanceUnit		 NauticalMiles
            AviatorTimeUnit		 Hours
            AviatorAltitudeUnit		 Feet
            AviatorFuelQuantityUnit		 Pounds
            AviatorRunwayLengthUnit		 Kilofeet
            AviatorBearingAngleUnit		 Degrees
            AviatorAngleOfAttackUnit		 Degrees
            AviatorAttitudeAngleUnit		 Degrees
            AviatorGUnit		 StandardSeaLevelG
            SolidAngle		 Steradians
            AviatorTSFCUnit		 TSFCLbmHrLbf
            AviatorPSFCUnit		 PSFCLbmHrHp
            AviatorForceUnit		 Pounds
            AviatorPowerUnit		 Horsepower
            SpectralBandwidthUnit		 Hertz
            AviatorAltTimeUnit		 Minutes
            AviatorSmallTimeUnit		 Seconds
            AviatorEnergyUnit		 kilowatt-hours
            BitsUnit		 MegaBits
            MagneticFieldUnit		 nanoTesla
        END ConnectReportUnits

        BEGIN ReportFavorites
            BEGIN Class
                Name		 Satellite
                BEGIN Favorite
                    Type		 Report
                    BaseDir		 User
                    Style		 TEMEofDate Position Velocity
                END Favorite
                BEGIN Favorite
                    Type		 Report
                    BaseDir		 User
                    Style		 TEMEofDate_Position_Velocity
                END Favorite
                BEGIN Favorite
                    Type		 Report
                    BaseDir		 User
                    Style		 J2000_ECI_Position_Velocity
                END Favorite
                BEGIN Favorite
                    Type		 Report
                    BaseDir		 Install
                    Style		 Fixed Position Velocity
                END Favorite
                BEGIN Favorite
                    Type		 Report
                    BaseDir		 User
                    Style		 ICRF_ECI_Position_Velocity_Accel
                END Favorite
                BEGIN Favorite
                    Type		 Report
                    BaseDir		 Install
                    Style		 B1950 Position Velocity
                END Favorite
                BEGIN Favorite
                    Type		 Report
                    BaseDir		 Install
                    Style		 MOD Position Velocity
                END Favorite
                BEGIN Favorite
                    Type		 Report
                    BaseDir		 Install
                    Style		 TOD Position Velocity
                END Favorite
            END Class
        END ReportFavorites

        BEGIN ADFFileData
        END ADFFileData

        BEGIN GenDb

            BEGIN Database
                DbType		 Satellite
                DefDb		 stkAllTLE.sd
                UseMyDb		 Off
                MaxMatches		 2000
                Use4SOC		 On

                BEGIN FieldDefaults

                END FieldDefaults

            END Database

            BEGIN Database
                DbType		 City
                DefDb		 stkCityDb.cd
                UseMyDb		 Off
                MaxMatches		 2000
                Use4SOC		 On

                BEGIN FieldDefaults

                END FieldDefaults

            END Database

            BEGIN Database
                DbType		 Facility
                DefDb		 stkFacility.fd
                UseMyDb		 Off
                MaxMatches		 2000
                Use4SOC		 On

                BEGIN FieldDefaults

                END FieldDefaults

            END Database
        END GenDb

        BEGIN SOCDb
            BEGIN Defaults
            END Defaults
        END SOCDb

        BEGIN Msgp4Ext
        END Msgp4Ext

        BEGIN FileLocations
        END FileLocations

        BEGIN Author
            Optimize		 No
            UseBasicGlobe		 No
            SaveEphemeris		 Yes
            SaveScenFolder		 No
            BEGIN ExternalFileTypes
                BEGIN Type
                    FileType		 Calculation Scalar
                    Include		 Yes
                END Type
                BEGIN Type
                    FileType		 Celestial Image
                    Include		 Yes
                END Type
                BEGIN Type
                    FileType		 Cloud
                    Include		 Yes
                END Type
                BEGIN Type
                    FileType		 EOP
                    Include		 Yes
                END Type
                BEGIN Type
                    FileType		 External Vector Data
                    Include		 Yes
                END Type
                BEGIN Type
                    FileType		 Globe
                    Include		 Yes
                END Type
                BEGIN Type
                    FileType		 Globe Data
                    Include		 Yes
                END Type
                BEGIN Type
                    FileType		 Map
                    Include		 Yes
                END Type
                BEGIN Type
                    FileType		 Map Image
                    Include		 Yes
                END Type
                BEGIN Type
                    FileType		 Marker/Label
                    Include		 Yes
                END Type
                BEGIN Type
                    FileType		 Model
                    Include		 Yes
                END Type
                BEGIN Type
                    FileType		 Object Break-up File
                    Include		 No
                END Type
                BEGIN Type
                    FileType		 Planetary Ephemeris
                    Include		 No
                END Type
                BEGIN Type
                    FileType		 Report Style Script
                    Include		 Yes
                END Type
                BEGIN Type
                    FileType		 Report/Graph Style
                    Include		 Yes
                END Type
                BEGIN Type
                    FileType		 Scalar Calculation File
                    Include		 Yes
                END Type
                BEGIN Type
                    FileType		 Terrain
                    Include		 Yes
                END Type
                BEGIN Type
                    FileType		 Volume Grid Intervals File
                    Include		 Yes
                END Type
                BEGIN Type
                    FileType		 Volumetric File
                    Include		 Yes
                END Type
                BEGIN Type
                    FileType		 WTM
                    Include		 Yes
                END Type
            END ExternalFileTypes
            ReadOnly		 No
            ViewerPassword		 No
            STKPassword		 No
            ExcludeInstallFiles		 No
            BEGIN ExternalFileList
            END ExternalFileList
        END Author

        BEGIN ExportDataFile
            FileType		 Ephemeris
            IntervalType		 Ephemeris
            TimePeriodStart		 0
            TimePeriodStop		 0
            StepType		 Ephemeris
            StepSize		 60
            EphemType		 STK
            UseVehicleCentralBody		 Yes
            CentralBody		 Earth
            SatelliteID		 -200000
            CoordSys		 ICRF
            NonSatCoordSys		 Fixed
            InterpolateBoundaries		 Yes
            EphemFormat		 Current
            InterpType		 9
            InterpOrder		 5
            AttCoordSys		 Fixed
            Quaternions		 0
            ExportCovar		 Position
            AttitudeFormat		 Current
            TimePrecision		 6
            CCSDSDateFormat		 YMD
            CCSDSEphFormat		 SciNotation
            CCSDSTimeSystem		 UTC
            CCSDSRefFrame		 ICRF
            UseSatCenterAndFrame		 No
            IncludeCovariance		 No
            IncludeAcceleration		 No
            CCSDSFileFormat		 KVN
        END ExportDataFile

        BEGIN Desc
        END Desc

        BEGIN RfEnv
<?xml version = "1.0" standalone = "yes"?>
<SCOPE>
    <VAR name = "PropagationChannel">
        <SCOPE>
            <VAR name = "UseITU618Section2p5">
                <BOOL>false</BOOL>
            </VAR>
            <VAR name = "UseCloudFogModel">
                <BOOL>false</BOOL>
            </VAR>
            <VAR name = "CloudFogModel">
                <VAR name = "ITU-R_P840-7">
                    <SCOPE Class = "CloudFogLossModel">
                        <VAR name = "Version">
                            <STRING>&quot;1.0.0 a&quot;</STRING>
                        </VAR>
                        <VAR name = "ComponentName">
                            <STRING>&quot;ITU-R_P840-7&quot;</STRING>
                        </VAR>
                        <VAR name = "Type">
                            <STRING>&quot;ITU-R P840-7&quot;</STRING>
                        </VAR>
                        <VAR name = "LiquidWaterDensityValueChoice">
                            <STRING>&quot;Liquid Water Content Density Value&quot;</STRING>
                        </VAR>
                        <VAR name = "CloudCeiling">
                            <QUANTITY Dimension = "DistanceUnit" Unit = "m">
                                <REAL>3000</REAL>
                            </QUANTITY>
                        </VAR>
                        <VAR name = "CloudLayerThickness">
                            <QUANTITY Dimension = "DistanceUnit" Unit = "m">
                                <REAL>500</REAL>
                            </QUANTITY>
                        </VAR>
                        <VAR name = "CloudTemp">
                            <QUANTITY Dimension = "Temperature" Unit = "K">
                                <REAL>273.15</REAL>
                            </QUANTITY>
                        </VAR>
                        <VAR name = "CloudLiqWaterDensity">
                            <QUANTITY Dimension = "SmallDensity" Unit = "kg*m^-3">
                                <REAL>0.0001</REAL>
                            </QUANTITY>
                        </VAR>
                        <VAR name = "AnnualAveragePercentValue">
                            <QUANTITY Dimension = "Percent" Unit = "unitValue">
                                <REAL>0.01</REAL>
                            </QUANTITY>
                        </VAR>
                        <VAR name = "MonthlyAveragePercentValue">
                            <QUANTITY Dimension = "Percent" Unit = "unitValue">
                                <REAL>0.01</REAL>
                            </QUANTITY>
                        </VAR>
                        <VAR name = "LiqWaterAverageDataMonth">
                            <INT>1</INT>
                        </VAR>
                        <VAR name = "UseRainHeightAsCloudThickness">
                            <BOOL>false</BOOL>
                        </VAR>
                    </SCOPE>
                </VAR>
            </VAR>
            <VAR name = "UseTropoScintModel">
                <BOOL>false</BOOL>
            </VAR>
            <VAR name = "TropoScintModel">
                <VAR name = "ITU-R_P618-8">
                    <SCOPE Class = "TropoScintLossModel">
                        <VAR name = "Version">
                            <STRING>&quot;1.0.0 a&quot;</STRING>
                        </VAR>
                        <VAR name = "ComponentName">
                            <STRING>&quot;ITU-R_P618-8&quot;</STRING>
                        </VAR>
                        <VAR name = "Type">
                            <STRING>&quot;ITU-R P618-8&quot;</STRING>
                        </VAR>
                        <VAR name = "Category">
                            <STRING>&quot;Previous Versions&quot;</STRING>
                        </VAR>
                        <VAR name = "ComputeDeepFade">
                            <BOOL>false</BOOL>
                        </VAR>
                        <VAR name = "FadeOutage">
                            <QUANTITY Dimension = "Percent" Unit = "unitValue">
                                <REAL>0.001</REAL>
                            </QUANTITY>
                        </VAR>
                        <VAR name = "PercentTimeRefracGrad">
                            <QUANTITY Dimension = "Percent" Unit = "unitValue">
                                <REAL>0.1</REAL>
                            </QUANTITY>
                        </VAR>
                        <VAR name = "SurfaceTemperature">
                            <QUANTITY Dimension = "Temperature" Unit = "K">
                                <REAL>273.15</REAL>
                            </QUANTITY>
                        </VAR>
                    </SCOPE>
                </VAR>
            </VAR>
            <VAR name = "UseIonoFadingModel">
                <BOOL>false</BOOL>
            </VAR>
            <VAR name = "IonoFadingModel">
                <VAR name = "ITU-R_P531-13">
                    <SCOPE Class = "IonoFadingLossModel">
                        <VAR name = "Version">
                            <STRING>&quot;1.0.0 a&quot;</STRING>
                        </VAR>
                        <VAR name = "ComponentName">
                            <STRING>&quot;ITU-R_P531-13&quot;</STRING>
                        </VAR>
                        <VAR name = "Type">
                            <STRING>&quot;ITU-R P531-13&quot;</STRING>
                        </VAR>
                        <VAR name = "UseAlternateAPFile">
                            <BOOL>false</BOOL>
                        </VAR>
                        <VAR name = "AlternateAPDataFile">
                            <STRING>
                                <PROP name = "FullName">
                                    <STRING>&quot;&quot;</STRING>
                                </PROP>&quot;&quot;</STRING>
                        </VAR>
                    </SCOPE>
                </VAR>
            </VAR>
            <VAR name = "UseRainModel">
                <BOOL>false</BOOL>
            </VAR>
            <VAR name = "RainModel">
                <VAR name = "ITU-R_P618-12">
                    <SCOPE Class = "RainLossModel">
                        <VAR name = "Version">
                            <STRING>&quot;1.0.0 a&quot;</STRING>
                        </VAR>
                        <VAR name = "ComponentName">
                            <STRING>&quot;ITU-R_P618-12&quot;</STRING>
                        </VAR>
                        <VAR name = "Type">
                            <STRING>&quot;ITU-R P618-12&quot;</STRING>
                        </VAR>
                        <VAR name = "SurfaceTemperature">
                            <QUANTITY Dimension = "Temperature" Unit = "K">
                                <REAL>273.15</REAL>
                            </QUANTITY>
                        </VAR>
                        <VAR name = "EnableDepolarizationLoss">
                            <BOOL>false</BOOL>
                        </VAR>
                    </SCOPE>
                </VAR>
            </VAR>
            <VAR name = "UseAtmosAbsorptionModel">
                <BOOL>false</BOOL>
            </VAR>
            <VAR name = "AtmosAbsorptionModel">
                <VAR name = "Simple_Satcom">
                    <SCOPE Class = "AtmosphericAbsorptionModel">
                        <VAR name = "Version">
                            <STRING>&quot;1.0.1 a&quot;</STRING>
                        </VAR>
                        <VAR name = "ComponentName">
                            <STRING>&quot;Simple_Satcom&quot;</STRING>
                        </VAR>
                        <VAR name = "Type">
                            <STRING>&quot;Simple Satcom&quot;</STRING>
                        </VAR>
                        <VAR name = "SurfaceTemperature">
                            <QUANTITY Dimension = "Temperature" Unit = "K">
                                <REAL>293.15</REAL>
                            </QUANTITY>
                        </VAR>
                        <VAR name = "WaterVaporConcentration">
                            <QUANTITY Dimension = "Density" Unit = "g*m^-3">
                                <REAL>7.5</REAL>
                            </QUANTITY>
                        </VAR>
                    </SCOPE>
                </VAR>
            </VAR>
            <VAR name = "UseUrbanTerresPropLossModel">
                <BOOL>false</BOOL>
            </VAR>
            <VAR name = "UrbanTerresPropLossModel">
                <VAR name = "Two_Ray">
                    <SCOPE Class = "UrbanTerrestrialPropagationLossModel">
                        <VAR name = "Version">
                            <STRING>&quot;1.0.0 a&quot;</STRING>
                        </VAR>
                        <VAR name = "ComponentName">
                            <STRING>&quot;Two_Ray&quot;</STRING>
                        </VAR>
                        <VAR name = "Type">
                            <STRING>&quot;Two Ray&quot;</STRING>
                        </VAR>
                        <VAR name = "SurfaceTemperature">
                            <QUANTITY Dimension = "Temperature" Unit = "K">
                                <REAL>273.15</REAL>
                            </QUANTITY>
                        </VAR>
                        <VAR name = "LossFactor">
                            <REAL>1</REAL>
                        </VAR>
                    </SCOPE>
                </VAR>
            </VAR>
            <VAR name = "UseCustomA">
                <BOOL>false</BOOL>
            </VAR>
            <VAR name = "UseCustomB">
                <BOOL>false</BOOL>
            </VAR>
            <VAR name = "UseCustomC">
                <BOOL>false</BOOL>
            </VAR>
        </SCOPE>
    </VAR>
    <VAR name = "EarthTemperature">
        <QUANTITY Dimension = "Temperature" Unit = "K">
            <REAL>290</REAL>
        </QUANTITY>
    </VAR>
    <VAR name = "RainOutagePercent">
        <PROP name = "FormatString">
            <STRING>&quot;%#6.3f&quot;</STRING>
        </PROP>
        <REAL>0.1</REAL>
    </VAR>
    <VAR name = "ActiveCommSystem">
        <LINKTOOBJ>
            <STRING>&quot;None&quot;</STRING>
        </LINKTOOBJ>
    </VAR>
</SCOPE>        END RfEnv

        BEGIN CommRad
        END CommRad

        BEGIN RadarCrossSection
<?xml version = "1.0" standalone = "yes"?>
<SCOPE>
    <VAR name = "Model">
        <VAR name = "Radar_Cross_Section">
            <SCOPE Class = "RCS">
                <VAR name = "Version">
                    <STRING>&quot;1.0.0 a&quot;</STRING>
                </VAR>
                <VAR name = "ComponentName">
                    <STRING>&quot;Radar_Cross_Section&quot;</STRING>
                </VAR>
                <VAR name = "Type">
                    <STRING>&quot;Radar Cross Section&quot;</STRING>
                </VAR>
                <VAR name = "FrequencyBandList">
                    <LIST>
                        <SCOPE>
                            <VAR name = "MinFrequency">
                                <QUANTITY Dimension = "BandwidthUnit" Unit = "Hz">
                                    <REAL>2997920</REAL>
                                </QUANTITY>
                            </VAR>
                            <VAR name = "ComputeTypeStrategy">
                                <VAR name = "Constant Value">
                                    <SCOPE Class = "RCS Compute Strategy">
                                        <VAR name = "ConstantValue">
                                            <QUANTITY Dimension = "RcsUnit" Unit = "sqm">
                                                <REAL>1</REAL>
                                            </QUANTITY>
                                        </VAR>
                                        <VAR name = "Type">
                                            <STRING>&quot;Constant Value&quot;</STRING>
                                        </VAR>
                                        <VAR name = "ComponentName">
                                            <STRING>&quot;Constant Value&quot;</STRING>
                                        </VAR>
                                    </SCOPE>
                                </VAR>
                            </VAR>
                            <VAR name = "SwerlingCase">
                                <STRING>&quot;0&quot;</STRING>
                            </VAR>
                        </SCOPE>
                        <SCOPE>
                            <VAR name = "MinFrequency">
                                <QUANTITY Dimension = "BandwidthUnit" Unit = "Hz">
                                    <REAL>3000000</REAL>
                                </QUANTITY>
                            </VAR>
                            <VAR name = "ComputeTypeStrategy">
                                <VAR name = "Constant Value">
                                    <SCOPE Class = "RCS Compute Strategy">
                                        <VAR name = "ConstantValue">
                                            <QUANTITY Dimension = "RcsUnit" Unit = "sqm">
                                                <REAL>1</REAL>
                                            </QUANTITY>
                                        </VAR>
                                        <VAR name = "Type">
                                            <STRING>&quot;Constant Value&quot;</STRING>
                                        </VAR>
                                        <VAR name = "ComponentName">
                                            <STRING>&quot;Constant Value&quot;</STRING>
                                        </VAR>
                                    </SCOPE>
                                </VAR>
                            </VAR>
                            <VAR name = "SwerlingCase">
                                <STRING>&quot;0&quot;</STRING>
                            </VAR>
                        </SCOPE>
                    </LIST>
                </VAR>
            </SCOPE>
        </VAR>
    </VAR>
</SCOPE>        END RadarCrossSection

        BEGIN RadarClutter
<?xml version = "1.0" standalone = "yes"?>
<SCOPE>
    <VAR name = "ClutterMap">
        <VAR name = "Constant Coefficient">
            <SCOPE Class = "Clutter Map">
                <VAR name = "ClutterCoefficient">
                    <QUANTITY Dimension = "RatioUnit" Unit = "units">
                        <REAL>1</REAL>
                    </QUANTITY>
                </VAR>
                <VAR name = "Type">
                    <STRING>&quot;Constant Coefficient&quot;</STRING>
                </VAR>
                <VAR name = "ComponentName">
                    <STRING>&quot;Constant Coefficient&quot;</STRING>
                </VAR>
            </SCOPE>
        </VAR>
    </VAR>
</SCOPE>        END RadarClutter

        BEGIN Gator
        END Gator

        BEGIN Crdn
        END Crdn

        BEGIN SpiceExt
        END SpiceExt

        BEGIN FlightScenExt
        END FlightScenExt

        BEGIN Graphics

            BEGIN Animation

                StartTime		 6 Apr 2004 07:51:28.386009000
                EndTime		 7 Apr 2004 19:00:00.000000000
                CurrentTime		 6 Apr 2004 07:51:28.386009008
                Direction		 Forward
                UpdateDelta		 60
                RefreshDelta		 0.010000
                XRealTimeMult		 1
                RealTimeOffset		 0
                XRtStartFromPause		                Yes		

            END Animation


            BEGIN DisplayFlags
                ShowLabels		 On
                ShowPassLabel		 Off
                ShowElsetNum		 Off
                ShowGndTracks		 On
                ShowGndMarkers		 On
                ShowOrbitMarkers		 On
                ShowPlanetOrbits		 Off
                ShowPlanetCBIPos		 On
                ShowPlanetCBILabel		 On
                ShowPlanetGndPos		 On
                ShowPlanetGndLabel		 On
                ShowSensors		 On
                ShowWayptMarkers		 Off
                ShowWayptTurnMarkers		 Off
                ShowOrbits		 On
                ShowDtedRegions		 Off
                ShowAreaTgtCentroids		 On
                ShowToolBar		 On
                ShowStatusBar		 On
                ShowScrollBars		 On
                AllowAnimUpdate		 On
                AccShowLine		 On
                AccAnimHigh		 On
                AccStatHigh		 On
                ShowPrintButton		 On
                ShowAnimButtons		 On
                ShowAnimModeButtons		 On
                ShowZoomMsrButtons		 On
                ShowMapCbButton		 Off
            END DisplayFlags

            BEGIN WinFonts

                System
                MS Sans Serif,16,0,0
                MS Sans Serif,20,0,0

            END WinFonts

            BEGIN MapData

                BEGIN TerrainConverterData
                    NorthLat		  0.0000000000000000e+00
                    EastLon		  0.0000000000000000e+00
                    SouthLat		  0.0000000000000000e+00
                    WestLon		  0.0000000000000000e+00
                    ColorByRGB		 No
                    AltsFromMSL		 No
                    UseColorRamp		 Yes
                    UseRegionMinMax		 Yes
                    SizeSameAsSrc		 Yes
                    MinAltHSV		  0.0000000000000000e+00  6.9999999999999996e-01  8.0000000000000004e-01  4.0000000000000002e-01
                    MaxAltHSV		  1.0000000000000000e+06  0.0000000000000000e+00  2.0000000000000001e-01  1.0000000000000000e+00
                    SmoothColors		 Yes
                    CreateChunkTrn		 No
                    OutputFormat		 PDTTX
                END TerrainConverterData

                DisableDefKbdActions		 Off
                TextShadowStyle		 None
                TextShadowColor		 #000000
                BingLevelOfDetailScale		 2
                BEGIN Map
                    MapNum		 1
                    TrackingMode		 LatLon
                    PickEnabled		 On
                    PanEnabled		 On

                    BEGIN MapAttributes
                        PrimaryBody		 Earth
                        SecondaryBody		 Sun
                        CenterLatitude		 0
                        CenterLongitude		 0
                        ProjectionAltitude		 63621860
                        FieldOfView		 35
                        OrthoDisplayDistance		 20000000
                        TransformTrajectory		 On
                        EquatorialRadius		 6378137
                        BackgroundColor		 #000000
                        LatLonLines		 On
                        LatSpacing		 30
                        LonSpacing		 30
                        LatLonLineColor		 #999999
                        LatLonLineStyle		 2
                        ShowOrthoDistGrid		 Off
                        OrthoGridXSpacing		 5
                        OrthoGridYSpacing		 5
                        OrthoGridColor		 #ffffff
                        ShowImageExtents		 Off
                        ImageExtentLineColor		 #ffffff
                        ImageExtentLineStyle		 0
                        ImageExtentLineWidth		 1
                        ShowImageNames		 Off
                        ImageNameFont		 0
                        Projection		 EquidistantCylindrical
                        Resolution		 VeryLow
                        CoordinateSys		 ECF
                        UseBackgroundImage		 On
                        UseBingForBackground		 Off
                        BingType		 Aerial
                        BingLogoHorizAlign		 Right
                        BingLogoVertAlign		 Bottom
                        BackgroundImageFile		 Basic.bmp
                        UseNightLights		 Off
                        NightLightsFactor		 3.5
                        UseCloudsFile		 Off
                        BEGIN ZoomLocations
                            BEGIN ZoomLocation
                                CenterLat		 0
                                CenterLon		 0
                                ZoomWidth		 359.999998
                                ZoomHeight		 180
                            END ZoomLocation
                        END ZoomLocations
                        UseVarAspectRatio		 No
                        SwapMapResolution		 Yes
                        NoneToVLowSwapDist		 2000000
                        VLowToLowSwapDist		 20000
                        LowToMediumSwapDist		 10000
                        MediumToHighSwapDist		 5000
                        HighToVHighSwapDist		 1000
                        VHighToSHighSwapDist		 100
                        BEGIN Axes
                            DisplayAxes		 no
                            CoordSys		 CBI
                            2aryCB		 Sun
                            Display+x		 yes
                            Label+x		 yes
                            Color+x		 #ffffff
                            Scale+x		 3
                            Display-x		 yes
                            Label-x		 yes
                            Color-x		 #ffffff
                            Scale-x		 3
                            Display+y		 yes
                            Label+y		 yes
                            Color+y		 #ffffff
                            Scale+y		 3
                            Display-y		 yes
                            Label-y		 yes
                            Color-y		 #ffffff
                            Scale-y		 3
                            Display+z		 yes
                            Label+z		 yes
                            Color+z		 #ffffff
                            Scale+z		 3
                            Display-z		 yes
                            Label-z		 yes
                            Color-z		 #ffffff
                            Scale-z		 3
                        END Axes

                    END MapAttributes

                    BEGIN MapList
                        BEGIN Detail
                            Alias		 RWDB2_Coastlines
                            Show		 Yes
                            Color		 #8fbc8f
                        END Detail
                        BEGIN Detail
                            Alias		 RWDB2_International_Borders
                            Show		 Yes
                            Color		 #8fbc8f
                        END Detail
                        BEGIN Detail
                            Alias		 RWDB2_Islands
                            Show		 Yes
                            Color		 #8fbc8f
                        END Detail
                        BEGIN Detail
                            Alias		 RWDB2_Lakes
                            Show		 No
                            Color		 #87cefa
                        END Detail
                        BEGIN Detail
                            Alias		 RWDB2_Provincial_Borders
                            Show		 No
                            Color		 #8fbc8f
                        END Detail
                        BEGIN Detail
                            Alias		 RWDB2_Rivers
                            Show		 No
                            Color		 #87cefa
                        END Detail
                    END MapList


                    BEGIN MapAnnotations
                    END MapAnnotations

                    BEGIN DisplayFlags
                        ShowLabels		 On
                        ShowPassLabel		 Off
                        ShowElsetNum		 Off
                        ShowGndTracks		 On
                        ShowGndMarkers		 On
                        ShowOrbitMarkers		 On
                        ShowPlanetOrbits		 Off
                        ShowPlanetCBIPos		 On
                        ShowPlanetCBILabel		 On
                        ShowPlanetGndPos		 On
                        ShowPlanetGndLabel		 On
                        ShowSensors		 On
                        ShowWayptMarkers		 Off
                        ShowWayptTurnMarkers		 Off
                        ShowOrbits		 On
                        ShowDtedRegions		 Off
                        ShowAreaTgtCentroids		 On
                        ShowToolBar		 On
                        ShowStatusBar		 On
                        ShowScrollBars		 On
                        AllowAnimUpdate		 Off
                        AccShowLine		 On
                        AccAnimHigh		 On
                        AccStatHigh		 On
                        ShowPrintButton		 On
                        ShowAnimButtons		 On
                        ShowAnimModeButtons		 On
                        ShowZoomMsrButtons		 On
                        ShowMapCbButton		 Off
                    END DisplayFlags

                    BEGIN RecordMovie
                        OutputFormat		 WMV
                        SdfSelected		 No
                        Directory		 C:\Users\dvallado\Documents\STK 9\Scenario1
                        BaseName		 Frame
                        Digits		 4
                        Frame		 0
                        LastAnimTime		 0
                        OutputMode		 Normal
                        HiResAssembly		 Assemble
                        HRWidth		 6000
                        HRHeight		 4500
                        HRDPI		 600
                        UseSnapInterval		 No
                        SnapInterval		 0
                        WmvCodec		 "Windows Media Video 9"
                        Framerate		 30
                        Bitrate		 3000000
                    END RecordMovie


                    BEGIN TimeDisplay
                        Show		 0
                        TextColor		 #ffffff
                        TextTranslucency		 0
                        ShowBackground		 0
                        BackColor		 #4d4d4d
                        BackTranslucency		 0.4
                        XPosition		 20
                        YPosition		 -20
                    END TimeDisplay

                    BEGIN LightingData
                        DisplayAltitude		 0
                        SubsolarPoint		 Off
                        SubsolarPointColor		 #ffff00
                        SubsolarPointMarkerStyle		 2

                        ShowUmbraLine		 Off
                        UmbraLineColor		 #0000ff
                        UmbraLineStyle		 0
                        UmbraLineWidth		 2
                        FillUmbra		 Off
                        UmbraFillColor		 #0000ff
                        ShowSunlightLine		 Off
                        SunlightLineColor		 #ffff00
                        SunlightLineStyle		 0
                        SunlightLineWidth		 2
                        FillSunlight		 Off
                        SunlightFillColor		 #ffff00
                        SunlightMinOpacity		 0
                        SunlightMaxOpacity		 0.2
                        UmbraMaxOpacity		 0.7
                        UmbraMinOpacity		 0.4
                    END LightingData
                END Map

                BEGIN MapStyles

                    UseStyleTime		 No

                    BEGIN Style
                        Name		 No_Map_Bckgrnd
                        Time		 241873713.613991
                        UpdateDelta		 60

                        BEGIN MapAttributes
                            PrimaryBody		 Earth
                            SecondaryBody		 Sun
                            CenterLatitude		 0
                            CenterLongitude		 0
                            ProjectionAltitude		 63621860
                            FieldOfView		 35
                            OrthoDisplayDistance		 20000000
                            TransformTrajectory		 On
                            EquatorialRadius		 6378137
                            BackgroundColor		 #000000
                            LatLonLines		 On
                            LatSpacing		 30
                            LonSpacing		 30
                            LatLonLineColor		 #999999
                            LatLonLineStyle		 2
                            ShowOrthoDistGrid		 Off
                            OrthoGridXSpacing		 5
                            OrthoGridYSpacing		 5
                            OrthoGridColor		 #ffffff
                            ShowImageExtents		 Off
                            ImageExtentLineColor		 #ffffff
                            ImageExtentLineStyle		 0
                            ImageExtentLineWidth		 1
                            ShowImageNames		 Off
                            ImageNameFont		 0
                            Projection		 EquidistantCylindrical
                            Resolution		 VeryLow
                            CoordinateSys		 ECF
                            UseBackgroundImage		 Off
                            UseBingForBackground		 Off
                            BingType		 Aerial
                            BingLogoHorizAlign		 Right
                            BingLogoVertAlign		 Bottom
                            UseNightLights		 Off
                            NightLightsFactor		 3.5
                            UseCloudsFile		 Off
                            BEGIN ZoomLocations
                                BEGIN ZoomLocation
                                    CenterLat		 0
                                    CenterLon		 0
                                    ZoomWidth		 359.999998
                                    ZoomHeight		 180
                                END ZoomLocation
                            END ZoomLocations
                            UseVarAspectRatio		 No
                            SwapMapResolution		 Yes
                            NoneToVLowSwapDist		 2000000
                            VLowToLowSwapDist		 20000
                            LowToMediumSwapDist		 10000
                            MediumToHighSwapDist		 5000
                            HighToVHighSwapDist		 1000
                            VHighToSHighSwapDist		 100
                            BEGIN Axes
                                DisplayAxes		 no
                                CoordSys		 CBI
                                2aryCB		 Sun
                                Display+x		 yes
                                Label+x		 yes
                                Color+x		 #ffffff
                                Scale+x		 3
                                Display-x		 yes
                                Label-x		 yes
                                Color-x		 #ffffff
                                Scale-x		 3
                                Display+y		 yes
                                Label+y		 yes
                                Color+y		 #ffffff
                                Scale+y		 3
                                Display-y		 yes
                                Label-y		 yes
                                Color-y		 #ffffff
                                Scale-y		 3
                                Display+z		 yes
                                Label+z		 yes
                                Color+z		 #ffffff
                                Scale+z		 3
                                Display-z		 yes
                                Label-z		 yes
                                Color-z		 #ffffff
                                Scale-z		 3
                            END Axes

                        END MapAttributes

                        BEGIN MapList
                            BEGIN Detail
                                Alias		 RWDB2_Coastlines
                                Show		 Yes
                                Color		 #00ff00
                            END Detail
                            BEGIN Detail
                                Alias		 RWDB2_International_Borders
                                Show		 Yes
                                Color		 #00ff00
                            END Detail
                            BEGIN Detail
                                Alias		 RWDB2_Islands
                                Show		 Yes
                                Color		 #00ff00
                            END Detail
                            BEGIN Detail
                                Alias		 RWDB2_Lakes
                                Show		 Yes
                                Color		 #87cefa
                            END Detail
                            BEGIN Detail
                                Alias		 RWDB2_Provincial_Borders
                                Show		 Yes
                                Color		 #00ff00
                            END Detail
                            BEGIN Detail
                                Alias		 RWDB2_Rivers
                                Show		 No
                                Color		 #00ff00
                            END Detail
                        END MapList


                        BEGIN MapAnnotations
                        END MapAnnotations

                        BEGIN RecordMovie
                            OutputFormat		 WMV
                            SdfSelected		 No
                            Directory		 C:\Users\dvallado\Documents\STK 9\Scenario1
                            BaseName		 Frame
                            Digits		 4
                            Frame		 0
                            LastAnimTime		 0
                            OutputMode		 Normal
                            HiResAssembly		 Assemble
                            HRWidth		 6000
                            HRHeight		 4500
                            HRDPI		 600
                            UseSnapInterval		 No
                            SnapInterval		 0
                            WmvCodec		 "Windows Media Video 9"
                            Framerate		 30
                            Bitrate		 3000000
                        END RecordMovie


                        BEGIN TimeDisplay
                            Show		 0
                            TextColor		 #00ffff
                            TextTranslucency		 0
                            ShowBackground		 0
                            BackColor		 #000000
                            BackTranslucency		 0
                            XPosition		 20
                            YPosition		 -20
                        END TimeDisplay

                        BEGIN LightingData
                            DisplayAltitude		 0
                            SubsolarPoint		 Off
                            SubsolarPointColor		 #ffff00
                            SubsolarPointMarkerStyle		 2

                            ShowUmbraLine		 Off
                            UmbraLineColor		 #ffff00
                            UmbraLineStyle		 0
                            UmbraLineWidth		 1
                            FillUmbra		 Off
                            UmbraFillColor		 #000000
                            ShowSunlightLine		 Off
                            SunlightLineColor		 #ffff00
                            SunlightLineStyle		 0
                            SunlightLineWidth		 1
                            FillSunlight		 Off
                            SunlightFillColor		 #ffff00
                            SunlightMinOpacity		 0.1
                            SunlightMaxOpacity		 0.5
                            UmbraMaxOpacity		 0.5
                            UmbraMinOpacity		 0.2
                        END LightingData

                        ShowDtedRegions		 Off

                    END Style

                    BEGIN Style
                        Name		 Basic_Map_Bckgrnd
                        Time		 241873713.613991
                        UpdateDelta		 60

                        BEGIN MapAttributes
                            PrimaryBody		 Earth
                            SecondaryBody		 Sun
                            CenterLatitude		 0
                            CenterLongitude		 0
                            ProjectionAltitude		 63621860
                            FieldOfView		 35
                            OrthoDisplayDistance		 20000000
                            TransformTrajectory		 On
                            EquatorialRadius		 6378137
                            BackgroundColor		 #000000
                            LatLonLines		 On
                            LatSpacing		 30
                            LonSpacing		 30
                            LatLonLineColor		 #999999
                            LatLonLineStyle		 2
                            ShowOrthoDistGrid		 Off
                            OrthoGridXSpacing		 5
                            OrthoGridYSpacing		 5
                            OrthoGridColor		 #ffffff
                            ShowImageExtents		 Off
                            ImageExtentLineColor		 #ffffff
                            ImageExtentLineStyle		 0
                            ImageExtentLineWidth		 1
                            ShowImageNames		 Off
                            ImageNameFont		 0
                            Projection		 EquidistantCylindrical
                            Resolution		 VeryLow
                            CoordinateSys		 ECF
                            UseBackgroundImage		 On
                            UseBingForBackground		 Off
                            BingType		 Aerial
                            BingLogoHorizAlign		 Right
                            BingLogoVertAlign		 Bottom
                            BackgroundImageFile		 Basic.bmp
                            UseNightLights		 Off
                            NightLightsFactor		 3.5
                            UseCloudsFile		 Off
                            BEGIN ZoomLocations
                                BEGIN ZoomLocation
                                    CenterLat		 0
                                    CenterLon		 0
                                    ZoomWidth		 359.999998
                                    ZoomHeight		 180
                                END ZoomLocation
                            END ZoomLocations
                            UseVarAspectRatio		 No
                            SwapMapResolution		 Yes
                            NoneToVLowSwapDist		 2000000
                            VLowToLowSwapDist		 20000
                            LowToMediumSwapDist		 10000
                            MediumToHighSwapDist		 5000
                            HighToVHighSwapDist		 1000
                            VHighToSHighSwapDist		 100
                            BEGIN Axes
                                DisplayAxes		 no
                                CoordSys		 CBI
                                2aryCB		 Sun
                                Display+x		 yes
                                Label+x		 yes
                                Color+x		 #ffffff
                                Scale+x		 3
                                Display-x		 yes
                                Label-x		 yes
                                Color-x		 #ffffff
                                Scale-x		 3
                                Display+y		 yes
                                Label+y		 yes
                                Color+y		 #ffffff
                                Scale+y		 3
                                Display-y		 yes
                                Label-y		 yes
                                Color-y		 #ffffff
                                Scale-y		 3
                                Display+z		 yes
                                Label+z		 yes
                                Color+z		 #ffffff
                                Scale+z		 3
                                Display-z		 yes
                                Label-z		 yes
                                Color-z		 #ffffff
                                Scale-z		 3
                            END Axes

                        END MapAttributes

                        BEGIN MapList
                            BEGIN Detail
                                Alias		 RWDB2_Coastlines
                                Show		 Yes
                                Color		 #8fbc8f
                            END Detail
                            BEGIN Detail
                                Alias		 RWDB2_International_Borders
                                Show		 Yes
                                Color		 #8fbc8f
                            END Detail
                            BEGIN Detail
                                Alias		 RWDB2_Islands
                                Show		 Yes
                                Color		 #8fbc8f
                            END Detail
                            BEGIN Detail
                                Alias		 RWDB2_Lakes
                                Show		 No
                                Color		 #00ff00
                            END Detail
                            BEGIN Detail
                                Alias		 RWDB2_Provincial_Borders
                                Show		 No
                                Color		 #00ff00
                            END Detail
                            BEGIN Detail
                                Alias		 RWDB2_Rivers
                                Show		 No
                                Color		 #00ff00
                            END Detail
                        END MapList


                        BEGIN MapAnnotations
                        END MapAnnotations

                        BEGIN RecordMovie
                            OutputFormat		 WMV
                            SdfSelected		 No
                            Directory		 C:\Users\dvallado\Documents\STK 9\Scenario1
                            BaseName		 Frame
                            Digits		 4
                            Frame		 0
                            LastAnimTime		 0
                            OutputMode		 Normal
                            HiResAssembly		 Assemble
                            HRWidth		 6000
                            HRHeight		 4500
                            HRDPI		 600
                            UseSnapInterval		 No
                            SnapInterval		 0
                            WmvCodec		 "Windows Media Video 9"
                            Framerate		 30
                            Bitrate		 3000000
                        END RecordMovie


                        BEGIN TimeDisplay
                            Show		 0
                            TextColor		 #00ffff
                            TextTranslucency		 0
                            ShowBackground		 0
                            BackColor		 #000000
                            BackTranslucency		 0
                            XPosition		 20
                            YPosition		 -20
                        END TimeDisplay

                        BEGIN LightingData
                            DisplayAltitude		 0
                            SubsolarPoint		 Off
                            SubsolarPointColor		 #ffff00
                            SubsolarPointMarkerStyle		 2

                            ShowUmbraLine		 Off
                            UmbraLineColor		 #ffff00
                            UmbraLineStyle		 0
                            UmbraLineWidth		 1
                            FillUmbra		 Off
                            UmbraFillColor		 #000000
                            ShowSunlightLine		 Off
                            SunlightLineColor		 #ffff00
                            SunlightLineStyle		 0
                            SunlightLineWidth		 1
                            FillSunlight		 Off
                            SunlightFillColor		 #ffff00
                            SunlightMinOpacity		 0.1
                            SunlightMaxOpacity		 0.5
                            UmbraMaxOpacity		 0.5
                            UmbraMinOpacity		 0.2
                        END LightingData

                        ShowDtedRegions		 Off

                    END Style

                    BEGIN Style
                        Name		 Orthographic_Projection
                        Time		 241873713.613991
                        UpdateDelta		 60

                        BEGIN MapAttributes
                            PrimaryBody		 Earth
                            SecondaryBody		 Sun
                            CenterLatitude		 0
                            CenterLongitude		 0
                            ProjectionAltitude		 63621860
                            FieldOfView		 35
                            OrthoDisplayDistance		 20000000
                            TransformTrajectory		 On
                            EquatorialRadius		 6378137
                            BackgroundColor		 #000000
                            LatLonLines		 On
                            LatSpacing		 30
                            LonSpacing		 30
                            LatLonLineColor		 #999999
                            LatLonLineStyle		 2
                            ShowOrthoDistGrid		 Off
                            OrthoGridXSpacing		 5
                            OrthoGridYSpacing		 5
                            OrthoGridColor		 #ffffff
                            ShowImageExtents		 Off
                            ImageExtentLineColor		 #ffffff
                            ImageExtentLineStyle		 0
                            ImageExtentLineWidth		 1
                            ShowImageNames		 Off
                            ImageNameFont		 0
                            Projection		 Orthographic
                            Resolution		 VeryLow
                            CoordinateSys		 ECF
                            UseBackgroundImage		 Off
                            UseBingForBackground		 Off
                            BingType		 Aerial
                            BingLogoHorizAlign		 Right
                            BingLogoVertAlign		 Bottom
                            UseNightLights		 Off
                            NightLightsFactor		 3.5
                            UseCloudsFile		 Off
                            BEGIN ZoomLocations
                                BEGIN ZoomLocation
                                    CenterLat		 0
                                    CenterLon		 0
                                    ZoomWidth		 359.99999
                                    ZoomHeight		 180
                                END ZoomLocation
                            END ZoomLocations
                            UseVarAspectRatio		 No
                            SwapMapResolution		 Yes
                            NoneToVLowSwapDist		 2000000
                            VLowToLowSwapDist		 20000
                            LowToMediumSwapDist		 10000
                            MediumToHighSwapDist		 5000
                            HighToVHighSwapDist		 1000
                            VHighToSHighSwapDist		 100
                            BEGIN Axes
                                DisplayAxes		 no
                                CoordSys		 CBI
                                2aryCB		 Sun
                                Display+x		 yes
                                Label+x		 yes
                                Color+x		 #ffffff
                                Scale+x		 3
                                Display-x		 yes
                                Label-x		 yes
                                Color-x		 #ffffff
                                Scale-x		 3
                                Display+y		 yes
                                Label+y		 yes
                                Color+y		 #ffffff
                                Scale+y		 3
                                Display-y		 yes
                                Label-y		 yes
                                Color-y		 #ffffff
                                Scale-y		 3
                                Display+z		 yes
                                Label+z		 yes
                                Color+z		 #ffffff
                                Scale+z		 3
                                Display-z		 yes
                                Label-z		 yes
                                Color-z		 #ffffff
                                Scale-z		 3
                            END Axes

                        END MapAttributes

                        BEGIN MapList
                            BEGIN Detail
                                Alias		 RWDB2_Coastlines
                                Show		 Yes
                                Color		 #00ff00
                            END Detail
                            BEGIN Detail
                                Alias		 RWDB2_International_Borders
                                Show		 Yes
                                Color		 #00ff00
                            END Detail
                            BEGIN Detail
                                Alias		 RWDB2_Islands
                                Show		 Yes
                                Color		 #00ff00
                            END Detail
                            BEGIN Detail
                                Alias		 RWDB2_Lakes
                                Show		 Yes
                                Color		 #87cefa
                            END Detail
                            BEGIN Detail
                                Alias		 RWDB2_Provincial_Borders
                                Show		 Yes
                                Color		 #00ff00
                            END Detail
                            BEGIN Detail
                                Alias		 RWDB2_Rivers
                                Show		 No
                                Color		 #00ff00
                            END Detail
                        END MapList


                        BEGIN MapAnnotations
                        END MapAnnotations

                        BEGIN RecordMovie
                            OutputFormat		 WMV
                            SdfSelected		 No
                            Directory		 C:\Users\dvallado\Documents\STK 9\Scenario1
                            BaseName		 Frame
                            Digits		 4
                            Frame		 0
                            LastAnimTime		 0
                            OutputMode		 Normal
                            HiResAssembly		 Assemble
                            HRWidth		 6000
                            HRHeight		 4500
                            HRDPI		 600
                            UseSnapInterval		 No
                            SnapInterval		 0
                            WmvCodec		 "Windows Media Video 9"
                            Framerate		 30
                            Bitrate		 3000000
                        END RecordMovie


                        BEGIN TimeDisplay
                            Show		 0
                            TextColor		 #00ffff
                            TextTranslucency		 0
                            ShowBackground		 0
                            BackColor		 #000000
                            BackTranslucency		 0
                            XPosition		 20
                            YPosition		 -20
                        END TimeDisplay

                        BEGIN LightingData
                            DisplayAltitude		 0
                            SubsolarPoint		 Off
                            SubsolarPointColor		 #ffff00
                            SubsolarPointMarkerStyle		 2

                            ShowUmbraLine		 Off
                            UmbraLineColor		 #ffff00
                            UmbraLineStyle		 0
                            UmbraLineWidth		 1
                            FillUmbra		 Off
                            UmbraFillColor		 #000000
                            ShowSunlightLine		 Off
                            SunlightLineColor		 #ffff00
                            SunlightLineStyle		 0
                            SunlightLineWidth		 1
                            FillSunlight		 Off
                            SunlightFillColor		 #ffff00
                            SunlightMinOpacity		 0.1
                            SunlightMaxOpacity		 0.5
                            UmbraMaxOpacity		 0.5
                            UmbraMinOpacity		 0.2
                        END LightingData

                        ShowDtedRegions		 Off

                    END Style

                    BEGIN Style
                        Name		 Zoom_North_America
                        Time		 241873713.613991
                        UpdateDelta		 60

                        BEGIN MapAttributes
                            PrimaryBody		 Earth
                            SecondaryBody		 Sun
                            CenterLatitude		 0
                            CenterLongitude		 0
                            ProjectionAltitude		 63621860
                            FieldOfView		 35
                            OrthoDisplayDistance		 20000000
                            TransformTrajectory		 On
                            EquatorialRadius		 6378137
                            BackgroundColor		 #000000
                            LatLonLines		 On
                            LatSpacing		 30
                            LonSpacing		 30
                            LatLonLineColor		 #999999
                            LatLonLineStyle		 2
                            ShowOrthoDistGrid		 Off
                            OrthoGridXSpacing		 5
                            OrthoGridYSpacing		 5
                            OrthoGridColor		 #ffffff
                            ShowImageExtents		 Off
                            ImageExtentLineColor		 #ffffff
                            ImageExtentLineStyle		 0
                            ImageExtentLineWidth		 1
                            ShowImageNames		 Off
                            ImageNameFont		 0
                            Projection		 EquidistantCylindrical
                            Resolution		 Low
                            CoordinateSys		 ECF
                            UseBackgroundImage		 On
                            UseBingForBackground		 Off
                            BingType		 Aerial
                            BingLogoHorizAlign		 Right
                            BingLogoVertAlign		 Bottom
                            BackgroundImageFile		 Basic.bmp
                            UseNightLights		 Off
                            NightLightsFactor		 3.5
                            UseCloudsFile		 Off
                            BEGIN ZoomLocations
                                BEGIN ZoomLocation
                                    CenterLat		 0
                                    CenterLon		 0
                                    ZoomWidth		 359.999998
                                    ZoomHeight		 180
                                END ZoomLocation
                                BEGIN ZoomLocation
                                    CenterLat		 48.235294
                                    CenterLon		 -106.229507
                                    ZoomWidth		 129.836065
                                    ZoomHeight		 64.918032
                                END ZoomLocation
                            END ZoomLocations
                            UseVarAspectRatio		 No
                            SwapMapResolution		 Yes
                            NoneToVLowSwapDist		 2000000
                            VLowToLowSwapDist		 20000
                            LowToMediumSwapDist		 10000
                            MediumToHighSwapDist		 5000
                            HighToVHighSwapDist		 1000
                            VHighToSHighSwapDist		 100
                            BEGIN Axes
                                DisplayAxes		 no
                                CoordSys		 CBI
                                2aryCB		 Sun
                                Display+x		 yes
                                Label+x		 yes
                                Color+x		 #ffffff
                                Scale+x		 3
                                Display-x		 yes
                                Label-x		 yes
                                Color-x		 #ffffff
                                Scale-x		 3
                                Display+y		 yes
                                Label+y		 yes
                                Color+y		 #ffffff
                                Scale+y		 3
                                Display-y		 yes
                                Label-y		 yes
                                Color-y		 #ffffff
                                Scale-y		 3
                                Display+z		 yes
                                Label+z		 yes
                                Color+z		 #ffffff
                                Scale+z		 3
                                Display-z		 yes
                                Label-z		 yes
                                Color-z		 #ffffff
                                Scale-z		 3
                            END Axes

                        END MapAttributes

                        BEGIN MapList
                            BEGIN Detail
                                Alias		 RWDB2_Coastlines
                                Show		 Yes
                                Color		 #00ff00
                            END Detail
                            BEGIN Detail
                                Alias		 RWDB2_International_Borders
                                Show		 Yes
                                Color		 #00ff00
                            END Detail
                            BEGIN Detail
                                Alias		 RWDB2_Islands
                                Show		 Yes
                                Color		 #00ff00
                            END Detail
                            BEGIN Detail
                                Alias		 RWDB2_Lakes
                                Show		 Yes
                                Color		 #87cefa
                            END Detail
                            BEGIN Detail
                                Alias		 RWDB2_Provincial_Borders
                                Show		 Yes
                                Color		 #00ff00
                            END Detail
                            BEGIN Detail
                                Alias		 RWDB2_Rivers
                                Show		 No
                                Color		 #00ff00
                            END Detail
                        END MapList


                        BEGIN MapAnnotations
                        END MapAnnotations

                        BEGIN RecordMovie
                            OutputFormat		 WMV
                            SdfSelected		 No
                            Directory		 C:\Users\dvallado\Documents\STK 9\Scenario1
                            BaseName		 Frame
                            Digits		 4
                            Frame		 0
                            LastAnimTime		 0
                            OutputMode		 Normal
                            HiResAssembly		 Assemble
                            HRWidth		 6000
                            HRHeight		 4500
                            HRDPI		 600
                            UseSnapInterval		 No
                            SnapInterval		 0
                            WmvCodec		 "Windows Media Video 9"
                            Framerate		 30
                            Bitrate		 3000000
                        END RecordMovie


                        BEGIN TimeDisplay
                            Show		 0
                            TextColor		 #00ffff
                            TextTranslucency		 0
                            ShowBackground		 0
                            BackColor		 #000000
                            BackTranslucency		 0
                            XPosition		 20
                            YPosition		 -20
                        END TimeDisplay

                        BEGIN LightingData
                            DisplayAltitude		 0
                            SubsolarPoint		 Off
                            SubsolarPointColor		 #ffff00
                            SubsolarPointMarkerStyle		 2

                            ShowUmbraLine		 Off
                            UmbraLineColor		 #ffff00
                            UmbraLineStyle		 0
                            UmbraLineWidth		 1
                            FillUmbra		 Off
                            UmbraFillColor		 #000000
                            ShowSunlightLine		 Off
                            SunlightLineColor		 #ffff00
                            SunlightLineStyle		 0
                            SunlightLineWidth		 1
                            FillSunlight		 Off
                            SunlightFillColor		 #ffff00
                            SunlightMinOpacity		 0.1
                            SunlightMaxOpacity		 0.5
                            UmbraMaxOpacity		 0.5
                            UmbraMinOpacity		 0.2
                        END LightingData

                        ShowDtedRegions		 Off

                    END Style

                END MapStyles

            END MapData

            BEGIN GfxClassPref

            END GfxClassPref


            BEGIN ConnectGraphicsOptions

                AsyncPickReturnUnique		 OFF

            END ConnectGraphicsOptions

        END Graphics

        BEGIN Overlays
        END Overlays

        BEGIN VO
        END VO

        BEGIN DIS

            BEGIN General

                Verbose		 Off
                Processing		 Off
                Statistics		 Off
                ExerciseID		 -1
                ForceID		 -1

            END General


            BEGIN Output

                Version		 5
                ExerciseID		 1
                forceID		 1
                HeartbeatTimer		 5
                DistanceThresh		 1
                OrientThresh		 3

            END Output


            BEGIN Time

                Mode		 rtPDUTimestamp

            END Time


            BEGIN PDUInfo


            END PDUInfo


            BEGIN Parameters

                ParmData		 COLORFRIENDLY blue
                ParmData		 COLORNEUTRAL white
                ParmData		 COLOROPFORCE red
                ParmData		 MAXDRELSETS 1000

            END Parameters


            BEGIN Network

                NetIF		 Default
                Mode		 Broadcast
                McastIP		 224.0.0.1
                Port		 3000
                rChannelBufferSize		 65000
                ReadBufferSize		 1500
                QueuePollPeriod		 20
                MaxRcvQueueEntries		 1000
                MaxRcvIOThreads		 4
                sChannelBufferSize		 65000

            END Network


            BEGIN EntityTypeDef


# order: kind:domain:country:catagory:subCatagory:specific:xtra ( -1 = * )


            END EntityTypeDef


            BEGIN EntityFilter
                Include		 *:*:*
            END EntityFilter

        END DIS

    END Extensions

    BEGIN SubObjects

        Class Satellite

            Satellite		

        END Class

    END SubObjects

    BEGIN References
        Instance *
            *		
        END Instance
        Instance Satellite/Satellite
            Satellite/Satellite		
        END Instance
    END References

END Scenario
