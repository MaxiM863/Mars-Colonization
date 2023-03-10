class combine_orbital_elements
{
    public J2000_elements J2000;
    public J2000_elements rates;
    public combine_orbital_elements()
    {
        J2000 = new J2000_elements();
        rates = new J2000_elements();
    }
}

class combine_state_vector
{
    public state_vector R0;
    public state_vector V0;
    public combine_state_vector()
    {
        R0 = new state_vector();
        V0 = new state_vector();
    }

}

class depart_and_arrive_data
{
    public int planet_id;
    public double year, month, day, hour, minute, second;
}

class DepartureAndArrival{
    int planet_id;
    double year;
    double month;
    double day;
    double hour;
    double minute;
    double second;
}

class heliocentric_elements
{
    public double h, e, RA, incl, w, TA, a, w_hat, L, M, E;
    public orbital_elements Get_Coe()
    {
        double deg = Math.PI / 180;

        orbital_elements oe = new orbital_elements();
        oe.a = a;
        oe.e = e;
        oe.h = h;
        oe.incl = incl*deg;
        oe.RA = RA*deg;
        oe.TA = TA*deg;
        oe.w = w*deg;

        return oe;
    }
}

class interplanetary_data
{
    public planet_data planet1;
    public planet_data planet2;
    public trajectory_data trajectory;
    public interplanetary_data()
    {
        planet1 = new planet_data();
        planet2 = new planet_data();
        trajectory = new trajectory_data();
    }
}

class J2000_elements
{
    public double a;
    public double e;
    public double i;
    public double RA;
    public double w_hat;
    public double L;
}

class matrix_3X3
{
    public double m11;
    public double m12;
    public double m13;
    public double m21;
    public double m22;
    public double m23;
    public double m31;
    public double m32;
    public double m33;
}

class orbital_elements
{
    public double h;
    public double e;
    public double RA;
    public double incl;
    public double w;
    public double TA;
    public double a;
}

class planet_data
{
    public heliocentric_elements coe;
    public state_vector Rp;
    public state_vector Vp;
    public double jd;
    public planet_data()
    {
        Rp = new state_vector();
        Vp = new state_vector();
        coe = new heliocentric_elements();
    }

}

class state_planet_epoch
{
    public heliocentric_elements coe;
    public state_vector r;
    public state_vector v;
    public double jd;
    public state_planet_epoch()
    {
        coe = new heliocentric_elements();
        r = new state_vector();
        v = new state_vector();
    }
}

class state_vector
{
    public double x;
    public double y;
    public double z;

    public state_vector(double x, double y, double z){
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public state_vector(){}
}

class trajectory_data
{
    public state_vector V1;
    public state_vector V2;
    public trajectory_data()
    {
        V1 = new state_vector();
        V2 = new state_vector();
    }
}

class quaternion
{
    public double w;
    public double x;
    public double y;
    public double z;

    public quaternion(double x, double y, double z, double w){
        this.w = w;
        this.x = x;
        this.y = y;
        this.z = z;
    }
}