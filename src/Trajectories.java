public class Trajectories {

    public void OptimalSpanTimeTrajectoryLambert(int departYear)
    {
        depart_and_arrive_data depart = new depart_and_arrive_data();
        depart_and_arrive_data arrival = new depart_and_arrive_data();

        depart.planet_id = 3;
        arrival.planet_id = 4;

        double valorisationScalar = Double.MAX_VALUE;
        double tmp_valorisationScalar = 0;

        int selectedIteration = 0;

        interplanetary_data tmp;

        for(int i = 0 ; i < 40; i++)
        {
            depart.year = departYear;
            depart.month = (i) % 12+1;
            depart.day = 15;

            if(i % 12 + 9 <= 12) {
                arrival.year = departYear;
                arrival.month = depart.month + 8;
            }
            else{
                arrival.year = departYear + 1 ;
                arrival.month = depart.month - 4;
            }
            arrival.day = 15;

            tmp = interplanet(depart, arrival);

            tmp_valorisationScalar = Math.sqrt((tmp.trajectory.V1.x - tmp.planet1.Vp.x) * (tmp.trajectory.V1.x - tmp.planet1.Vp.x) + (tmp.trajectory.V1.y - tmp.planet1.Vp.y) * (tmp.trajectory.V1.y - tmp.planet1.Vp.y) + (tmp.trajectory.V1.z - tmp.planet1.Vp.z) * (tmp.trajectory.V1.z - tmp.planet1.Vp.z));

            if(tmp_valorisationScalar < valorisationScalar){
                valorisationScalar = tmp_valorisationScalar;
                selectedIteration = i;
            }
        }

        boolean success = false;
        if(selectedIteration == 0)
        {
            if(valorisationScalar < 10){
                success = true;
            }
        }
        else {
            success = true;
        }

        if(success) {
            System.out.println(valorisationScalar);

            double DateDepartJulian = J0(departYear, 1, 15, depart.hour, depart.minute, depart.second) + selectedIteration * 30;
            //double DateArrivalJulian = J0(arrival.year, arrival.month, arrival.day, arrival.hour, arrival.minute, arrival.second);

            valorisationScalar = Double.MAX_VALUE;

            int kk = 0;
            int hh = 0;

            for (int k = -25; k < 25; k++) // Julian date are required
            {
                for(int h = -25; h < 25; h++)
                {
                    int[] fff = fromJulian(DateDepartJulian + k);
                    int[] ggg = fromJulian(DateDepartJulian + 240 + h);

                    depart.year = fff[0];
                    depart.month = fff[1];
                    depart.day = fff[2];

                    arrival.month = ggg[1];
                    arrival.day = ggg[2];
                    arrival.year = ggg[0];

                    tmp = interplanet(depart, arrival);

                    valorisationScalar = Double.MAX_VALUE;

                    tmp_valorisationScalar = Math.sqrt((tmp.trajectory.V1.x - tmp.planet1.Vp.x) * (tmp.trajectory.V1.x - tmp.planet1.Vp.x) + (tmp.trajectory.V1.y - tmp.planet1.Vp.y) * (tmp.trajectory.V1.y - tmp.planet1.Vp.y) + (tmp.trajectory.V1.z - tmp.planet1.Vp.z) * (tmp.trajectory.V1.z - tmp.planet1.Vp.z));

                    if(tmp_valorisationScalar < valorisationScalar){
                        valorisationScalar = tmp_valorisationScalar;
                        kk = k;
                        hh = h;
                    }
                }
            }

            System.out.println("Good!!! " + String.format("%.2f", (double)valorisationScalar) + "km/s Year: " + fromJulian(DateDepartJulian + kk)[0] + " Month: " + fromJulian(DateDepartJulian + kk)[1] + " Day: " + fromJulian(DateDepartJulian + kk)[2] + " Duration: " + (240 + hh + kk));
        }
        else{
            System.out.println("Bad!!! " + valorisationScalar);
        }
    }

    private static int[] fromJulian(double injulian) {
        int JGREG= 15 + 31*(10+12*1582);
        double HALFSECOND = 0.5;

        int jalpha,ja,jb,jc,jd,je,year,month,day;
        double julian = injulian + HALFSECOND / 86400.0;
        ja = (int) injulian;
        if (ja>= JGREG) {
            jalpha = (int) (((ja - 1867216) - 0.25) / 36524.25);
            ja = ja + 1 + jalpha - jalpha / 4;
        }

        jb = ja + 1524;
        jc = (int) (6680.0 + ((jb - 2439870) - 122.1) / 365.25);
        jd = 365 * jc + jc / 4;
        je = (int) ((jb - jd) / 30.6001);
        day = jb - jd - (int) (30.6001 * je);
        month = je - 1;
        if (month > 12) month = month - 12;
        year = jc - 4715;
        if (month > 2) year--;
        if (year <= 0) year--;

        return new int[] {year, month, day};
    }
    public interplanetary_data Lambert()
    {
        depart_and_arrive_data departure = new depart_and_arrive_data();

        departure.planet_id = 3;
        departure.year = 2035;
        departure.month = 5;
        departure.day = 13;
        departure.hour = 10;
        departure.minute = 15;
        departure.second = 30;

        depart_and_arrive_data arrival = new depart_and_arrive_data();

        arrival.planet_id = 4;
        arrival.year = 2036;
        arrival.month = 1;
        arrival.day = 8;
        arrival.hour = 10;
        arrival.minute = 15;
        arrival.second = 30;

        return interplanet(departure, arrival);
    }

    private interplanetary_data interplanet(depart_and_arrive_data depart, depart_and_arrive_data arrive)
    {

        //double mu = 1.327124E11; ////////////////////////////////////////// insert mu ///////////////////////////

        int planet_id = depart.planet_id;
        double year = depart.year;
        double month = depart.month;
        double day = depart.day;
        double hour = depart.hour;
        double minute = depart.minute;
        double second = depart.second;

        state_planet_epoch temp = planet_elements_and_sv(planet_id, year, month, day, hour, minute, second);

        planet_id = arrive.planet_id;
        year = arrive.year;
        month = arrive.month;
        day = arrive.day;
        hour = arrive.hour;
        minute = arrive.minute;
        second = arrive.second;

        state_planet_epoch temp2 = planet_elements_and_sv(planet_id, year, month, day, hour, minute, second);

        double tof = (temp2.jd - temp.jd) * 24 * 3600;

        combine_state_vector lam = lambert(temp.r, temp2.r, tof, false, true);

        interplanetary_data ret = new interplanetary_data();

        ret.planet1.coe = temp.coe;
        ret.planet1.Rp = temp.r;
        ret.planet1.Vp = temp.v;
        ret.planet1.jd = temp.jd;

        ret.planet2.Rp = temp2.r;
        ret.planet2.Vp = temp2.v;
        ret.planet2.jd = temp2.jd;

        ret.trajectory.V1 = lam.R0;
        ret.trajectory.V2 = lam.V0;

        return ret;
    }
    private state_planet_epoch planet_elements_and_sv(int planet_id, double year, double month, double day, double hour, double minute, double second)
    {
        double mu = 1.327124E11; ///////// Sun parameter
        double deg = Math.PI / 180; /// to translate degree to radian
        double j0 = J0(year, month, day, hour, minute, second);
        combine_orbital_elements rrr;
        rrr = planetary_elements(planet_id);
        double t0 = (j0 - 2451545) / 36525; /// time after 2000 january first, divide by a year, divide by a centenary
        double a = rrr.J2000.a + rrr.rates.a * t0;
        double e = rrr.J2000.e + rrr.rates.e * t0;
        double h = Math.sqrt(mu * a * (1 - Math.pow(e, 2)));
        double incl = rrr.J2000.i + rrr.rates.i * t0;
        double RA = zero_to_360(rrr.J2000.RA + rrr.rates.RA * t0);
        double w_hat = zero_to_360(rrr.J2000.w_hat + rrr.rates.w_hat * t0);
        double L = zero_to_360(rrr.J2000.L + rrr.rates.L * t0);
        double w = zero_to_360(w_hat - RA);
        double M = zero_to_360((L - w_hat));
        double E = kepler_E(e, M * deg);
        double TA = zero_to_360(2 * Math.atan(Math.sqrt((1 - e) / (1 + e)) * Math.tan(E / 2)) / deg);

        heliocentric_elements coe = new heliocentric_elements();

        coe.h = h;
        coe.e = e;
        coe.RA = RA;
        coe.incl = incl;
        coe.w = w;
        coe.TA = TA;
        coe.a = a;
        coe.w_hat = w_hat;
        coe.L = L;
        coe.M = M;
        coe.E = E / deg;
        orbital_elements oe = new orbital_elements();
        oe.h = h;
        oe.e = e;
        oe.RA = RA * deg;
        oe.incl = incl * deg;
        oe.w = w * deg;
        oe.TA = TA * deg;
        combine_state_vector csv = sv_from_coe(oe);
        state_planet_epoch ret = new state_planet_epoch();
        ret.coe = coe;
        ret.r = csv.R0;
        ret.v = csv.V0;
        ret.jd = j0;
        return ret;
    }

    private double J0(double year, double month, double day, double hour, double minute, double second)
    {
        return 367 * year - Math.floor(7 * (year + Math.floor((month + 9) / 12)) / 4) + Math.floor(275 * month / 9) + day + 1721013.5 + (hour + minute / 60 + second / 3600) / 24;
    }

    private combine_orbital_elements planetary_elements(int planet_id)
    {

        double[][] J2000_elements = new double[9][6];

        J2000_elements[0] = new double[]{ 0.38709893, 0.20563069, 7.00487, 48.33167, 77.45645, 252.25084 };
        J2000_elements[1] = new double[]{ 0.72333199, 0.00677323, 3.39471, 76.68069, 131.53298, 181.97973 };
        J2000_elements[2] = new double[]{ 1.00000011, 0.01671022, 0.00005, - 11.26064, 102.94719, 100.46435 };
        J2000_elements[3] = new double[]{ 1.52366231, 0.09341233, 1.85061, 49.57854, 336.04084, 355.45332};
        J2000_elements[4] = new double[]{ 5.20336301, 0.04839266, 1.30530, 100.55615, 14.75385, 34.40438 };
        J2000_elements[5] = new double[]{ 9.53707032, 0.05415060, 2.48446, 113.71504, 92.43194, 49.94432 };
        J2000_elements[6] = new double[]{ 19.19126393, 0.04716771, 0.76986, 74.22988, 170.96424, 313.23218 };
        J2000_elements[7] = new double[]{ 30.06896348, 0.00858587, 1.76917, 131.72169, 44.97135, 304.88003 };
        J2000_elements[8] = new double[]{ 39.48168677, 0.24880766, 17.14175, 110.30347, 224.06676, 238.92881 };

        double[][] cent_rates = new double[9][6];

        cent_rates[0] = new double[]{ 0.00000066, 0.00002527, - 23.51, - 446.30, 573.57, 538101628.29 };
        cent_rates[1] = new double[]{ 0.00000092, - 0.00004938, - 2.86, - 996.89, - 108.80, 210664136.06 };
        cent_rates[2] = new double[]{ -0.00000005, - 0.00003804, - 46.94, - 18228.25, 1198.28, 129597740.63 };
        cent_rates[3] = new double[]{ -0.00007221, 0.00011902, - 25.47, - 1020.19, 1560.78, 68905103.78 };
        cent_rates[4] = new double[]{ 0.00060737, - 0.00012880, - 4.15, 1217.17, 839.93, 10925078.35 };
        cent_rates[5] = new double[]{ -0.00301530, - 0.00036762, 6.11, - 1591.05, - 1948.89, 4401052.95 };
        cent_rates[6] = new double[]{ 0.00152025, - 0.00019150, - 2.09, - 1681.4, 1312.56, 1542547.79 };
        cent_rates[7] = new double[]{ -0.00125196, 0.00002514, - 3.64, - 151.25, - 844.43, 786449.21 };
        cent_rates[8] = new double[]{ -0.00076912, 0.00006465, 11.07, - 37.33, - 132.25, 522747.90 };


        combine_orbital_elements ret = new combine_orbital_elements();

        ret.J2000.a = J2000_elements[planet_id - 1][0];
        ret.J2000.e = J2000_elements[planet_id - 1][1];
        ret.J2000.i = J2000_elements[planet_id - 1][2];
        ret.J2000.RA = J2000_elements[planet_id - 1][3];
        ret.J2000.w_hat = J2000_elements[planet_id - 1][4];
        ret.J2000.L = J2000_elements[planet_id - 1][5];

        ret.rates.a = cent_rates[planet_id - 1][0];
        ret.rates.e = cent_rates[planet_id - 1][1];
        ret.rates.i = cent_rates[planet_id - 1][2];
        ret.rates.RA = cent_rates[planet_id - 1][3];
        ret.rates.w_hat = cent_rates[planet_id - 1][4];
        ret.rates.L = cent_rates[planet_id - 1][5];

        double AU = 149597871;
        ret.J2000.a = ret.J2000.a * AU;
        ret.rates.a = ret.rates.a * AU;

        ret.rates.i = ret.rates.i / 3600.0;
        ret.rates.L = ret.rates.L / 3600.0;
        ret.rates.RA = ret.rates.RA / 3600.0;
        ret.rates.w_hat = ret.rates.w_hat / 3600.0;

        return ret;
    }

    private double zero_to_360(double elem)
    {

        if (elem >= 360) elem = elem - Math.floor(elem / 360) * 360;
        else if (elem < 0)
        {
            //MessageBox.Show("floor");
            elem = elem + (Math.ceil(Math.abs(elem) / 360)) * 360; //// due to floor definition between matlab and c#
        }

        //MessageBox.Show("new angle: " + elem);

        return elem;
    }

    public double kepler_E(double e, double M)
    {
        double E;
        double error = 1E-11;

        if (M < Math.PI) E = M + e / 2;
        else E = M - e / 2;

        double ratio = 1;
        while (Math.abs(ratio) > error)
        {
            ratio = (E - e * Math.sin(E) - M) / (1 - e * Math.cos(E));
            E = E - ratio;
            //MessageBox.Show("@@@: " + E.ToString());
        }

        return E;
    }

    public combine_state_vector sv_from_coe(orbital_elements coe)
    {
        double mu = 1.327124E11; /////////////////////////////////////////// insert mu  ////////////////////////////////////////

        state_vector rp = new state_vector();
        rp.x = (Math.pow(coe.h, 2) / mu) * (1 / (1 + coe.e * Math.cos(coe.TA))) * Math.cos(coe.TA);
        rp.y = (Math.pow(coe.h, 2) / mu) * (1 / (1 + coe.e * Math.cos(coe.TA))) * Math.sin(coe.TA);
        rp.z = 0.0;

        state_vector vp = new state_vector();
        vp.x = (mu / coe.h) * (-Math.sin(coe.TA));
        vp.y = (mu / coe.h) * (coe.e + Math.cos(coe.TA));
        vp.z = 0.0;

        matrix_3X3 R3_W = new matrix_3X3();
        R3_W.m11 = Math.cos(coe.RA);
        R3_W.m12 = Math.sin(coe.RA);
        R3_W.m13 = 0.0;
        R3_W.m21 = -Math.sin(coe.RA);
        R3_W.m22 = Math.cos(coe.RA);
        R3_W.m23 = 0.0;
        R3_W.m31 = 0.0;
        R3_W.m32 = 0.0;
        R3_W.m33 = 1.0;

        matrix_3X3 R1_i = new matrix_3X3();
        R1_i.m11 = 1.0;
        R1_i.m12 = 0.0;
        R1_i.m13 = 0.0;
        R1_i.m21 = 0.0;
        R1_i.m22 = Math.cos(coe.incl);
        R1_i.m23 = Math.sin(coe.incl);
        R1_i.m31 = 0.0;
        R1_i.m32 = -Math.sin(coe.incl);
        R1_i.m33 = Math.cos(coe.incl);

        matrix_3X3 R3_w = new matrix_3X3();
        R3_w.m11 = Math.cos(coe.w);
        R3_w.m12 = Math.sin(coe.w);
        R3_w.m13 = 0.0;
        R3_w.m21 = -Math.sin(coe.w);
        R3_w.m22 = Math.cos(coe.w);
        R3_w.m23 = 0.0;
        R3_w.m31 = 0.0;
        R3_w.m32 = 0.0;
        R3_w.m33 = 1.0;

        matrix_3X3 m1 = matrix_transpose(R3_W);
        matrix_3X3 m2 = matrix_transpose(R1_i);
        matrix_3X3 m3 = matrix_transpose(R3_w);


        matrix_3X3 Q_pX = matrix_product(m1, m2);
        Q_pX = matrix_product(Q_pX, m3);          ///// to compute without quaternion: uncomment this, and comment the quaternion calculation


        combine_state_vector ret = new combine_state_vector();

        ret.R0 = matrix_vect_product(Q_pX, rp);
        ret.V0 = matrix_vect_product(Q_pX, vp);


        ret.R0 = rp;
        ret.V0 = vp;


        //*************************************************** Quaternions calculation modification to avoid Gimbal ******************************************************************************



        quaternion A = new quaternion(0, 0, Math.sin(coe.w / 2.0), Math.cos(coe.w / 2.0));

        state_vector V = rotate_vector_by_quaternion(new state_vector(ret.R0.x, ret.R0.y, ret.R0.z), A);

        ret.R0.x = V.x;
        ret.R0.y = V.y;
        ret.R0.z = V.z;

        V = rotate_vector_by_quaternion(new state_vector(ret.V0.x, ret.V0.y, ret.V0.z), A);

        ret.V0.x = V.x;
        ret.V0.y = V.y;
        ret.V0.z = V.z;


        A = new quaternion(Math.sin(coe.incl/ 2.0),0, 0,  Math.cos(coe.incl / 2.0));

        V = rotate_vector_by_quaternion(new state_vector(ret.R0.x, ret.R0.y, ret.R0.z), A);

        ret.R0.x = V.x;
        ret.R0.y = V.y;
        ret.R0.z = V.z;

        V = rotate_vector_by_quaternion(new state_vector(ret.V0.x, ret.V0.y, ret.V0.z), A);

        ret.V0.x = V.x;
        ret.V0.y = V.y;
        ret.V0.z = V.z;


        A = new quaternion(0, 0, Math.sin(coe.RA / 2.0), Math.cos(coe.RA / 2.0));

        V = rotate_vector_by_quaternion(new state_vector(ret.R0.x, ret.R0.y, ret.R0.z), A);

        ret.R0.x = V.x;
        ret.R0.y = V.y;
        ret.R0.z = V.z;

        V = rotate_vector_by_quaternion(new state_vector(ret.V0.x, ret.V0.y, ret.V0.z), A);

        ret.V0.x = V.x;
        ret.V0.y = V.y;
        ret.V0.z = V.z;


            /*
            matrix_3X3 tmp = m3;
            MessageBox.Show(coe.RA.ToString());
            //m3 = m1;

            m3 = matrix_transpose(m3);
            Matrix3x3 M = new Matrix3x3((float)m3.m11, (float)m3.m12, (float)m3.m13, (float)m3.m21, (float)m3.m22, (float)m3.m23, (float)m3.m31, (float)m3.m32, (float)m3.m33);
            M = Matrix3x3.Orthogonalize(M);

            if (Math.Abs(M.Determinant() - 1) > 0.001) MessageBox.Show("Interplanetary: SV-From_COE: quaternion calculus: Matrix deternimant <> 1: " + M.Determinant().ToString());

            Quaternion Q = Mat_To_Quat(M);

            Vector3 V = Vector3.Transform(new Vector3((float)ret.R0.x, (float)ret.R0.y, (float)ret.R0.z), Q);

            ret.R0.x = V.X;
            ret.R0.y = V.Y;
            ret.R0.z = V.Z;

            V = Vector3.Transform(new Vector3((float)ret.V0.x, (float)ret.V0.y, (float)ret.V0.z), Q);

            ret.V0.x = V.X;
            ret.V0.y = V.Y;
            ret.V0.z = V.Z;


            //---------------------
            //MessageBox.Show(coe.incl.ToString());

            m2 = matrix_transpose(m2);
            M = new Matrix3x3((float)m2.m11, (float)m2.m12, (float)m2.m13, (float)m2.m21, (float)m2.m22, (float)m2.m23, (float)m2.m31, (float)m2.m32, (float)m2.m33);
            M = Matrix3x3.Orthogonalize(M);

            if (Math.Abs(M.Determinant() - 1) > 0.001) MessageBox.Show("Interplanetary: SV-From_COE: quaternion calculus: Matrix deternimant <> 1: " + M.Determinant().ToString());

            Q = Mat_To_Quat(M);

            V = Vector3.Transform(new Vector3((float)ret.R0.x, (float)ret.R0.y, (float)ret.R0.z), Q);

            ret.R0.x = V.X;
            ret.R0.y = V.Y;
            ret.R0.z = V.Z;

            V = Vector3.Transform(new Vector3((float)ret.V0.x, (float)ret.V0.y, (float)ret.V0.z), Q);

            ret.V0.x = V.X;
            ret.V0.y = V.Y;
            ret.V0.z = V.Z;

            //------------------------
            //MessageBox.Show(coe.w.ToString());
            //m1 = tmp;

            m1 = matrix_transpose(m1);
             M = new Matrix3x3((float)m1.m11, (float)m1.m12, (float)m1.m13, (float)m1.m21, (float)m1.m22, (float)m1.m23, (float)m1.m31, (float)m1.m32, (float)m1.m33);
            M = Matrix3x3.Orthogonalize(M);

            if (Math.Abs(M.Determinant() - 1) > 0.001) MessageBox.Show("Interplanetary: SV-From_COE: quaternion calculus: Matrix deternimant <> 1: " + M.Determinant().ToString());

             Q = Mat_To_Quat(M);

             V = Vector3.Transform(new Vector3((float)ret.R0.x, (float)ret.R0.y, (float)ret.R0.z), Q);

            ret.R0.x = V.X;
            ret.R0.y = V.Y;
            ret.R0.z = V.Z;

            V = Vector3.Transform(new Vector3((float)ret.V0.x, (float)ret.V0.y, (float)ret.V0.z), Q);

            ret.V0.x = V.X;
            ret.V0.y = V.Y;
            ret.V0.z = V.Z;

    */
        //*************************************************** End of Quaternions calculation modification ******************************************************************************


        return ret;
    }

    private matrix_3X3 matrix_transpose(matrix_3X3 A)
    {
        matrix_3X3 ret = new matrix_3X3();

        ret.m11 = A.m11;
        ret.m21 = A.m12;
        ret.m31 = A.m13;

        ret.m12 = A.m21;
        ret.m22 = A.m22;
        ret.m32 = A.m23;

        ret.m13 = A.m31;
        ret.m23 = A.m32;
        ret.m33 = A.m33;

        return ret;
    }

    private matrix_3X3 matrix_product(matrix_3X3 A, matrix_3X3 B)
    {
        matrix_3X3 C = new matrix_3X3();

        C.m11 = A.m11 * B.m11 + A.m12 * B.m21 + A.m13 * B.m31;
        C.m12 = A.m11 * B.m12 + A.m12 * B.m22 + A.m13 * B.m32;
        C.m13 = A.m11 * B.m13 + A.m12 * B.m23 + A.m13 * B.m33;

        C.m21 = A.m21 * B.m11 + A.m22 * B.m21 + A.m23 * B.m31;
        C.m22 = A.m21 * B.m12 + A.m22 * B.m22 + A.m23 * B.m32;
        C.m23 = A.m21 * B.m13 + A.m22 * B.m23 + A.m23 * B.m33;

        C.m31 = A.m31 * B.m11 + A.m32 * B.m21 + A.m33 * B.m31;
        C.m32 = A.m31 * B.m12 + A.m32 * B.m22 + A.m33 * B.m32;
        C.m33 = A.m31 * B.m13 + A.m32 * B.m23 + A.m33 * B.m33;

        return C;
    }

    private state_vector matrix_vect_product(matrix_3X3 A, state_vector B)
    {
        state_vector ret = new state_vector();

        ret.x = A.m11 * B.x + A.m12 * B.y + A.m13 * B.z;
        ret.y = A.m21 * B.x + A.m22 * B.y + A.m23 * B.z;
        ret.z = A.m31 * B.x + A.m32 * B.y + A.m33 * B.z;

        return ret;
    }

    public combine_state_vector lambert(state_vector R1, state_vector R2, double t, boolean retrograde, boolean warning)
    {

        double mu = 1.327124E11; ///////////////cvdc

        double r1 = Math.sqrt(Math.pow(R1.x, 2) + Math.pow(R1.y, 2) + Math.pow(R1.z, 2));
        double r2 = Math.sqrt(Math.pow(R2.x, 2) + Math.pow(R2.y, 2) + Math.pow(R2.z, 2));

        state_vector c12 = new state_vector();
        c12.x = R1.y*R2.z - R1.z*R2.y;
        c12.y = -R1.x*R2.z + R1.z*R2.x;
        c12.z = R1.x*R2.y - R1.y*R2.x;

        double theta = Math.acos((R1.x*R2.x + R1.y*R2.y + R1.z*R2.z) / r1 / r2);

        //MessageBox.Show(theta.ToString());



        if (retrograde == false)
        {
            if (c12.z <= 0) theta = 2 * Math.PI - theta;
            //MessageBox.Show("allo");
        }
        else if (retrograde == true)
        {
            if (c12.z >= 0) theta = 2 * Math.PI - theta;
        }

        //if (Math.Abs(Math.Abs(theta) - 3.14159) < 0.1) theta = Math.Abs(theta) / theta * 3.14159 + Math.Abs(Math.Abs(theta) - 3.14159) / (theta - 3.14159) * Math.Abs(Math.Abs(theta) - 3.14159)*1.1;

        //MessageBox.Show(theta.ToString());


        double A = Math.sin(theta) * Math.sqrt(r1 * r2 / (1 - Math.cos(theta)));



        double z = -100.0f;


        double cntr = F(z, t, mu, A, r1, r2);

        boolean usedata = true;

        if (cntr <= 0 || cntr > 0) ;            // for detecting invalid number
        else
        {
            cntr = -1e6;
            //MessageBox.Show("z too high");                // for detecting invalide number
            usedata = false;
        }
        //else cntr=-1;

        //int test = 0;

        while ( cntr < 0)
        {
            //test++;
            z = z + 0.1f;
            cntr = F(z, t, mu, A, r1, r2);
            if (cntr <= 0 || cntr > 0) ;            // for detecting invalid number
            else
            {
                cntr = -1e6;
                //MessageBox.Show("Invalid: " + z);
                usedata = false;                           ///////////////////// need to verify that fact
            }
        }

        //
        //std::cout << "@@@@@@@@@@@@@@@@@@@  " << z << "  @@@@@@@@@@@@@@@@@@@" << std::endl;

        double tol = 1E-11;
        int nmax = 50000;

        double fact = 0.999;

        double ratio = 1.0;
        int n = 0;
        while (Math.abs(ratio) > tol && n <= nmax)
        {
            //MessageBox.Show(ratio.ToString());
            n = n + 1;
            ratio = F(z, t, mu, A, r1, r2) / dFdz(z,A,r1,r2);
            z = z - (ratio);
                /*if (F(z * fact, t, mu, A, r1, r2) > 0) z *= fact;
                else
                {
                    fact = 1.0 - (1.0 - fact) / 10.0;
                    if (fact == 1) break;
                }*/
            //
        }
//MessageBox.Show(F(z, t, mu, A, r1, r2).ToString() + "::" + z);


        combine_state_vector ret = new combine_state_vector();


        if (warning && n >= nmax)
        {
            //MessageBox.Show("Lambert: nombre d'iteration depassé z = " + z.ToString());
        }
        else if((!warning && n >= nmax) || !usedata)
        {
            //MessageBox.Show("WarningOff et Lambert: nombre d'iteration depassé z = " + z.ToString());
            ret.R0.x = Double.MAX_VALUE;
            ret.R0.y = Double.MAX_VALUE;
            ret.R0.z = Double.MAX_VALUE;

            ret.V0.x = Double.MAX_VALUE;
            ret.V0.y = Double.MAX_VALUE;
            ret.V0.z = Double.MAX_VALUE;
        }
        else
        {
            double f = 1.0 - y(z, r1, r2, A) / r1;
            double g = A * Math.sqrt(y(z, r1, r2, A) / mu);
            double gdot = 1.0 - y(z, r1, r2, A) / r2;

            //MessageBox.Show(f + " :: " + g);

            ret.R0.x = 1.0 / g * (R2.x - f * R1.x);
            ret.R0.y = 1.0 / g * (R2.y - f * R1.y);
            ret.R0.z = 1.0 / g * (R2.z - f * R1.z);

            ret.V0.x = 1.0 / g * (gdot * R2.x - R1.x);
            ret.V0.y = 1.0 / g * (gdot * R2.y - R1.y);
            ret.V0.z = 1.0 / g * (gdot * R2.z - R1.z);
        }

        return ret;

    }

    private double F(double z, double t, double mu, double A, double r1, double r2)
    {
        //std::cout << z << " " << t << " " << mu << " " << A << " " << r1 << " " << r2 << "*****************" << std::endl;
        double ret;
        ret = Math.pow((y(z, r1, r2, A) / C(z)), 1.5) * S(z) + A * Math.sqrt(y(z, r1, r2, A)) - Math.sqrt(mu) * t;
        //std::cout << ret << ".....................";
        return ret;
    }

    private double y(double z, double r1, double r2, double A)
    {
        double ret;
        ret = r1 + r2 + A * (z * S(z) - 1.0) / Math.sqrt(C(z));
        return ret;
    }

    private double S(double z)
    {
        double ret = stumpS(z);
        return ret;
    }

    private double C(double z)
    {
        double ret = stumpC(z);
        return ret;
    }

    private double dFdz(double z, double A, double r1, double r2)
    {
        double ret;

        if (z == 0)
        {
            ret = Math.sqrt(2.0) / 40.0 * Math.pow(y(0.0, r1, r2, A), 1.5) + A / 8.0 * Math.sqrt(y(0.0, r1, r2, A)) + A * Math.sqrt(1.0 / 2.0 / y(0.0, r1, r2, A));
        }
        else ret = Math.pow((y(z, r1, r2, A) / C(z)), 1.5) * (1.0 / 2.0 / z * (C(z) - 3.0 * S(z) / 2.0 / C(z)) + 3.0 * Math.pow(S(z), 2) / 4.0 / C(z)) + A / 8.0 * (3.0 * S(z) / C(z) * Math.sqrt(y(z, r1, r2, A)) + A * Math.sqrt(C(z) / y(z, r1, r2, A)));

        return ret;
    }

    //--------------------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------------------
    private double stumpS(double z)
    {
        double s = 0;

        if (z > 0) s = (Math.sqrt(z) - Math.sin(Math.sqrt(z))) / Math.pow((Math.sqrt(z)), 3);
        if (z < 0) s = (Math.sinh(Math.sqrt(-z)) - Math.sqrt(-z)) / Math.pow(Math.sqrt(-z), 3);
        if (z == 0) s = 1 / 6;

        return s;
    }

    private double stumpC(double z)
    {
        double c = 0;

        if (z > 0) c = (1 - Math.cos(Math.sqrt(z))) / z;
        if (z < 0) c = (Math.cosh(Math.sqrt(-z)) - 1) / (-z);
        if (z == 0) c = 1 / 2;

        return c;
    }

    private state_vector rotate_vector_by_quaternion(state_vector v, quaternion q)
    {
        // Extract the vector part of the quaternion
        state_vector u = new state_vector();

        u.x = q.x;
        u.y = q.y;
        u.z = q.z;

        // Extract the scalar part of the quaternion
        double s = q.w;

        double x;
        double y;
        double z;

        double dotuv = (u.x * v.x + u.y * v.y + u.z * v.z);
        double dotuu = (u.x * u.x + u.y * u.y + u.z * u.z);

        state_vector cross = new state_vector();
        cross.x = (u.y*v.z - u.z * v.y);
        cross.y = -(u.x*v.z - u.z * v.x);
        cross.z = (u.x*v.y - u.y * v.x);

        state_vector ret = new state_vector();

        ret.x = 2 * dotuv * u.x + (s * s - dotuu) * v.x + 2 * s * cross.x;
        ret.y = 2 * dotuv * u.y + (s * s - dotuu) * v.y + 2 * s * cross.y;
        ret.z = 2 * dotuv * u.z + (s * s - dotuu) * v.z + 2 * s * cross.z;

        return ret;
    }
}