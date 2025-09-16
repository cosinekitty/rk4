#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>
#include "rk4_integrator.hpp"
#include "rk4_simulator.hpp"

using test_func_t = int (*)();

struct Test
{
    const char *name;
    test_func_t func;
};


static int Logarithm();
static int Pendulum();
static int SolarSystem();
static int Catenary();


static Test TestList[] =
{
    { "log", Logarithm },
    { "pendulum", Pendulum },
    { "solsys", SolarSystem },
    { "catenary", Catenary }
};


constexpr int NumTests = static_cast<int>(sizeof(TestList) / sizeof(TestList[0]));
bool Verbose = false;


int main(int argc, const char *argv[])
{
    if (argc>1 && !strcmp(argv[1], "-v"))
    {
        Verbose = true;
        --argc;
        ++argv;
    }

    if (argc == 2)
    {
        const char *verb = argv[1];
        if (!strcmp(verb, "all"))
        {
            for (int i = 0; i < NumTests; ++i)
            {
                printf("Running: %s\n", TestList[i].name);
                if (TestList[i].func())
                {
                    printf("FAIL: %s\n", TestList[i].name);
                    return 1;
                }
                printf("\n");
            }
        }
        else
        {
            Test* test = nullptr;
            for (int i = 0; i < NumTests; ++i)
            {
                if (!strcmp(verb, TestList[i].name))
                {
                    test = &TestList[i];
                    break;
                }
            }
            if (test == nullptr)
            {
                printf("FAIL: Unknown test name '%s'\n", verb);
                return 1;
            }
            if (test->func())
            {
                printf("FAIL: %s\n", test->name);
                return 1;
            }
        }
        printf("SUCCESS\n");
        return 0;
    }
    else
    {
        printf("unittest: Invalid command line arguments\n");
        return 1;
    }
}


static int Logarithm()
{
    using state_t = RungeKutta::StateVector<2, double, double>;

    // Use the numerical integrator to estimate ln(2).
    // Let t be the independent variable that ranges [1, 2].
    // Let v = 1/t be the "velocity" or slope.
    // Let x = integral(v) be the "position" or integral of v over the range t=1 to t=2.
    // The state necessarily includes t because it is needed to calculate v.

    auto deriv = [](const state_t& state) -> state_t
    {
        state_t m;
        m.coord[0] = 1;                     // dt/dt = 1
        m.coord[1] = 1 / state.coord[0];    // v = dx/dt = 1/t
        return m;
    };

    using integ_t = RungeKutta::Integrator<double, state_t, decltype(deriv)>;
    integ_t integ(deriv);

    integ.state.coord[1] = 0;
    const int nSteps = 200;
    const double dt = 1.0 / nSteps;
    integ.state.coord[0] = 1;
    for (int n = 0; n < nSteps; ++n)
        integ.step(dt);

    const double x2 = integ.state.coord[0];
    const double xdiff = std::abs(2 - x2);
    printf("x2 = %0.16f, xdiff = %g\n", x2, xdiff);
    if (xdiff > 2.14e-14)
    {
        printf("FAIL(Logarithm): EXCESSIVE xdiff\n");
        return 1;
    }

    const double correct = std::log(2.0);
    const double diff = std::abs(correct - integ.state.coord[1]);
    const double tolerance = 1.3e-12;
    printf("Integral = %0.16lf\n", integ.state.coord[1]);
    printf("Correct  = %0.16lf\n", correct);
    printf("Diff     = %g\n", diff);

    if (diff > tolerance)
    {
        printf("Logarithm: EXCESSIVE ERROR\n");
        return 1;
    }

    printf("Logarithm: PASS\n");
    return 0;
}


static int Pendulum()
{
    // Goal: model a 2D pendulum of fixed length L, with mass m concentrated at a point on the end.
    // Represent the state in angular terms.
    // The displacement angle away from the vertical [radians] is called theta.
    // The derivative of theta with respect to time, i.e. the angular speed, is called omega.
    // The angular acceleration (the second time derivative) is called alpha.
    // The kinematic model is based on analysis in the tangential direction only,
    // because we assume the pendulum has a fixed length L. Because L is constant,
    // movement in the radial direction is not possible.
    // Therefore, we can exclude radial variables from the model without missing anything.
    // In terms of the RK4 integrator, we consider the state to be the current angular
    // position and angular speed as a pair: [theta, omega].
    // The derivative of the state gives another pair: [omega, alpha].
    // Thus the state type needs 2 numbers.
    using state_t = RungeKutta::StateVector<2, double, double>;

    static constexpr double g = 9.8;    // gravitational acceleration [m/s^2]
    static constexpr double L = 1;      // length of the pendulum [m]

    auto deriv = [](const state_t& state) -> state_t
    {
        // Calculate the derivative of [theta, omega] to obtain [omega, alpha].
        state_t slope;
        slope.coord[0] = state.coord[1];    // copy angular speed (omega) [rad/s]
        slope.coord[1] = (-g/L) * std::sin(state.coord[0]);    // calculate angular acceleration alpha [rad/s^2]
        return slope;
    };

    using integ_t = RungeKutta::Integrator<double, state_t, decltype(deriv)>;

    auto energy_over_mass = [](const state_t& state) -> double
    {
        const double theta = state.coord[0];
        const double v = L * state.coord[1];
        return L*g*(1 - std::cos(theta)) + v*v/2;
    };

    // We test the pendulum model by running the simulation
    // and measuring the following properties:
    //
    // 1. Time between consecutive axis-crossings (both left and right) should be identical.
    // 2. In general, because this is a frictionless pendulum, kinetic+potential energy must be conserved.
    // 3. The maximum and minimum amplitude should be identical.
    // 4. When run backwards in time, the results should be the same. (Pass in dt < 0)

    // Start the pendulum displaced away from the vertical, but at rest (speed = 0).
    static constexpr double A = 5.0 * (M_PI/180);   // initial angle [rad]
    integ_t integ(deriv);
    integ.state.coord[0] = A;
    integ.state.coord[1] = 0.0;

    const double dt = 0.01;
    int zeroCrossingCount = 0;
    double periodTimeSum = 0;
    double minCrossingTime = 0;
    double maxCrossingTime = 0;
    double prevCrossingTime = -1;
    const int zeroCrossingLimit = 100;
    double t = 0;
    double minEnergyOverMass = energy_over_mass(integ.state);
    double maxEnergyOverMass = minEnergyOverMass;
    while (zeroCrossingCount < zeroCrossingLimit && t < 600.0)
    {
        if (Verbose)
            printf("Pendulum: t=%0.2lf, theta=%0.6lf, omega=%0.6lf\n", t, integ.state.coord[0], integ.state.coord[1]);

        const double prevAngle = integ.state.coord[0];
        integ.step(dt);
        const double angle = integ.state.coord[0];

        const double em = energy_over_mass(integ.state);
        minEnergyOverMass = std::min(minEnergyOverMass, em);
        maxEnergyOverMass = std::max(maxEnergyOverMass, em);

        if (angle * prevAngle <= 0)
        {
            assert(angle != prevAngle);
            const double rho = prevAngle/(prevAngle - angle);
            assert(rho>=0 && rho<=1);
            const double crossingTime = t + rho*dt;
            if (prevCrossingTime >= 0)
            {
                ++zeroCrossingCount;
                double elapsed = crossingTime - prevCrossingTime;
                periodTimeSum += elapsed;
                if (maxCrossingTime == 0)
                {
                    minCrossingTime = elapsed;
                    maxCrossingTime = elapsed;
                }
                else
                {
                    minCrossingTime = std::min(minCrossingTime, elapsed);
                    maxCrossingTime = std::max(maxCrossingTime, elapsed);
                }
                if (Verbose)
                    printf("Pendulum: zero crossing %d at time %0.6lf\n", zeroCrossingCount, crossingTime);
            }
            prevCrossingTime = crossingTime;
        }
        t += dt;
    }

    if (zeroCrossingCount != zeroCrossingLimit)
    {
        printf("FAIL(Pendulum): expected %d zero-crossings, found %d\n", zeroCrossingLimit, zeroCrossingCount);
        return 1;
    }

    const double emRange = std::abs(maxEnergyOverMass - minEnergyOverMass);
    printf("          energy/mass: min=%0.16f, max=%0.16f, range=%g\n", minEnergyOverMass, maxEnergyOverMass, emRange);

    static constexpr double largeAngleCorrection = 1 + A*A/16 + A*A*A*A*(11.0/3072);
    printf("          large angle correction factor = %0.16lf\n", largeAngleCorrection);
    const double meanCrossingTime = periodTimeSum / zeroCrossingCount;
    const double expectedCrossingTime = M_PI * std::sqrt(L/g) * largeAngleCorrection;       // half-trip time
    const double diff = std::abs(expectedCrossingTime - meanCrossingTime);
    const double spread = std::abs(maxCrossingTime - minCrossingTime);
    printf("          %d zero-crossings\n", zeroCrossingCount);
    printf("          mean interval = %0.16f seconds\n", meanCrossingTime);
    printf("          expected      = %0.16f seconds\n", expectedCrossingTime);
    printf("          diff          = %g\n", diff);
    if (Verbose)
    {
        printf("          min           = %0.16f seconds\n", minCrossingTime);
        printf("          max           = %0.16f seconds\n", maxCrossingTime);
    }
    printf("          spread        = %g\n", spread);
    if (diff > 7.08e-09)
    {
        printf("FAIL(Pendulum): EXCESSIVE diff\n");
        return 1;
    }
    if (spread > 5.89e-07)
    {
        printf("FAIL(Pendulum): EXCESSIVE spread\n");
        return 1;
    }
    if (emRange > 4.91e-09)
    {
        printf("FAIL(Pendulum): EXCESSIVE energy/mass range\n");
        return 1;
    }
    printf("Pendulum: PASS\n");
    return 0;
}


//-------------------------------------------------------------------------------------

struct vec_t
{
    double x{};
    double y{};
    double z{};

    void operator += (const vec_t& other)
    {
        x += other.x;
        y += other.y;
        z += other.z;
    }

    void operator -= (const vec_t& other)
    {
        x -= other.x;
        y -= other.y;
        z -= other.z;
    }

    double magSquared() const
    {
        return x*x + y*y + z*z;
    }

    double mag() const
    {
        return std::sqrt(magSquared());
    }
};


inline vec_t operator * (double s, const vec_t& v)
{
    vec_t p;
    p.x = s * v.x;
    p.y = s * v.y;
    p.z = s * v.z;
    return p;
}


inline vec_t operator + (const vec_t& a, const vec_t& b)
{
    vec_t c;
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    c.z = a.z + b.z;
    return c;
}


inline vec_t operator - (const vec_t& a, const vec_t& b)
{
    vec_t c;
    c.x = a.x - b.x;
    c.y = a.y - b.y;
    c.z = a.z - b.z;
    return c;
}



struct body_state_t
{
    vec_t   pos;
    vec_t   vel;

    void operator += (const body_state_t& other)
    {
        pos += other.pos;
        vel += other.vel;
    }
};


inline body_state_t operator * (double s, const body_state_t& v)
{
    body_state_t b;
    b.pos = s * v.pos;
    b.vel = s * v.vel;
    return b;
}


inline body_state_t operator + (const body_state_t& a, const body_state_t& b)
{
    body_state_t c;
    c.pos = a.pos + b.pos;
    c.vel = a.vel + b.vel;
    return c;
}

constexpr unsigned nbodies = 5;

using system_state_t = RungeKutta::StateVector<nbodies, body_state_t, double>;


static void InitSolarSystem(system_state_t& ss)
{
    // Initialize the solar system with barycentric states.
    // See input/solsys/*.txt for reference. These files were generated by JPL Horizons.

    // Sun
    // 2451544.500000000, A.D. 2000-Jan-01 00:00:00.0000, -7.139867342351965E-03, -2.792293024590130E-03,  2.063456515496323E-04,  5.374261984736907E-06, -7.410968933449047E-06, -9.423856739196385E-08,  4.429381450148228E-05,  7.669236241600089E-03, -2.307578738849755E-06,
    ss.coord[0].pos.x = -7.139867342351965E-03;
    ss.coord[0].pos.y = -2.792293024590130E-03;
    ss.coord[0].pos.z =  2.063456515496323E-04;
    ss.coord[0].vel.x =  5.374261984736907E-06;
    ss.coord[0].vel.y = -7.410968933449047E-06;
    ss.coord[0].vel.z = -9.423856739196385E-08;

    // Jupiter
    // 2451544.500000000, A.D. 2000-Jan-01 00:00:00.0000,  3.996320621351185E+00,  2.932560894863292E+00, -1.016166987472685E-01, -4.558376533394469E-03,  6.439863253809189E-03,  7.537585811287700E-05,  2.863446881048567E-02,  4.957904584013533E+00,  1.333017909147770E-04,
    ss.coord[1].pos.x =  3.996320621351185E+00;
    ss.coord[1].pos.y =  2.932560894863292E+00;
    ss.coord[1].pos.z = -1.016166987472685E-01;
    ss.coord[1].vel.x = -4.558376533394469E-03;
    ss.coord[1].vel.y =  6.439863253809189E-03;
    ss.coord[1].vel.z =  7.537585811287700E-05;

    // Saturn
    // 2451544.500000000, A.D. 2000-Jan-01 00:00:00.0000,  6.401416168163572E+00,  6.565250459597368E+00, -3.689209424165721E-01, -4.285166236914331E-03,  3.884579924219014E-03,  1.025155157793841E-04,  5.300174802621731E-02,  9.176968193092039E+00, -2.142006356308177E-04,
    ss.coord[2].pos.x =  6.401416168163572E+00;
    ss.coord[2].pos.y =  6.565250459597368E+00;
    ss.coord[2].pos.z = -3.689209424165721E-01;
    ss.coord[2].vel.x = -4.285166236914331E-03;
    ss.coord[2].vel.y =  3.884579924219014E-03;
    ss.coord[2].vel.z =  1.025155157793841E-04;

    // Uranus
    // 2451544.500000000, A.D. 2000-Jan-01 00:00:00.0000,  1.442337622913901E+01, -1.373845196763103E+01, -2.379230398823268E-01,  2.683839889869174E-03,  2.665015934138373E-03, -2.484537350873912E-05,  1.150525366451675E-01,  1.992072919566711E+01,  1.055558926764745E-04,
    ss.coord[3].pos.x =  1.442337622913901E+01;
    ss.coord[3].pos.y = -1.373845196763103E+01;
    ss.coord[3].pos.z = -2.379230398823268E-01;
    ss.coord[3].vel.x =  2.683839889869174E-03;
    ss.coord[3].vel.y =  2.665015934138373E-03;
    ss.coord[3].vel.z = -2.484537350873912E-05;

    // Neptune
    // 2451544.500000000, A.D. 2000-Jan-01 00:00:00.0000,  1.680361785911107E+01, -2.499544357157440E+01,  1.274769173927470E-01,  2.584591104377939E-03,  1.768944210379373E-03, -9.629428153854004E-05,  1.739524398661104E-01,  3.011893130340558E+01, -2.647302594893257E-05,
    ss.coord[4].pos.x =  1.680361785911107E+01;
    ss.coord[4].pos.y = -2.499544357157440E+01;
    ss.coord[4].pos.z =  1.274769173927470E-01;
    ss.coord[4].vel.x =  2.584591104377939E-03;
    ss.coord[4].vel.y =  1.768944210379373E-03;
    ss.coord[4].vel.z = -9.629428153854004E-05;
}


static void PrintSolarSystem(const system_state_t& ss)
{
    const char *name[] = { "Sun", "Jupiter", "Saturn", "Uranus", "Neptune" };
    for (unsigned i = 0; i < nbodies; ++i)
    {
        printf("%-10s  rx=%20.16f ry=%20.16f rz=%20.16f  vx=%20.16f vy=%20.16f vz=%20.16f\n",
            name[i],
            ss.coord[i].pos.x,
            ss.coord[i].pos.y,
            ss.coord[i].pos.z,
            ss.coord[i].vel.x,
            ss.coord[i].vel.y,
            ss.coord[i].vel.z
        );
    }
    printf("\n");
}


static double VecDiff(const vec_t& vec, double x, double y, double z)
{
    double dx = vec.x - x;
    double dy = vec.y - y;
    double dz = vec.z - z;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}


static int VerifyBody(
    const char *name,
    const body_state_t& body,
    double rtol, double vtol,
    double rx, double ry, double rz,
    double vx, double vy, double vz)
{
    double posdiff = VecDiff(body.pos, rx, ry, rz);
    double veldiff = VecDiff(body.vel, vx, vy, vz);
    printf("%-10s:  posdiff=%8.6f, veldiff=%g\n", name, posdiff, veldiff);
    int rc = 0;
    if (posdiff > rtol)
    {
        printf("EXCESSIVE POSITION ERROR\n");
        rc = 1;
    }
    if (veldiff > vtol)
    {
        printf("EXCESSIVE VELOCITY ERROR\n");
        rc = 1;
    }
    return rc;
}


static int SolarSystem()
{
    /*
        https://web.archive.org/web/20120220062549/http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf

        Page 10 in the above document describes the constants used in the DE405 ephemeris.
        The following are GM values (gravity constant * mass) in [au^3 / day^2].
        This side-steps issues of not knowing the exact values of G and masses M[i];
        the products GM[i] are known extremely accurately.
    */
    static const double GM[] =
    {
        0.2959122082855911e-03,     // [0] Sun
        0.2825345909524226e-06,     // [1] Jupiter
        0.8459715185680659e-07,     // [2] Saturn
        0.1292024916781969e-07,     // [3] Uranus
        0.1524358900784276e-07      // [4] Neptune
    };

    auto deriv = [](const system_state_t& state) -> system_state_t
    {
        // Calculate sum of accelerations on each body caused by the other bodies.
        vec_t acc[nbodies]{};

        for (unsigned a = 0; a+1 < nbodies; ++a)
        {
            for (unsigned b = a+1; b < nbodies; ++b)
            {
                // Compute the vector from position a to position b.
                vec_t dr = state.coord[b].pos - state.coord[a].pos;
                double r2 = dr.magSquared();
                double pull = 1 / (r2 * std::sqrt(r2));
                acc[a] += (pull * GM[b]) * dr;
                acc[b] -= (pull * GM[a]) * dr;
            }
        }

        // This is the confusing thing about this approach.
        // We have to put dr/dt in pos and dv/dt in vel.

        system_state_t deriv;
        for (unsigned a = 0; a < nbodies; ++a)
        {
            deriv.coord[a].pos = state.coord[a].vel;
            deriv.coord[a].vel = acc[a];
        }
        return deriv;
    };

    using integ_t = RungeKutta::Integrator<double, system_state_t, decltype(deriv)>;
    integ_t integ(deriv);
    InitSolarSystem(integ.state);

    static constexpr double dt = 1.0;           // days
    static constexpr unsigned oversample = 2;
    static constexpr unsigned ndays = 36520;    // from year 2000 to year 2100
    static constexpr unsigned nsteps = oversample * ndays;

    for (unsigned i = 0; i < nsteps; ++i)
        integ.step(dt/oversample);

    PrintSolarSystem(integ.state);

    if (VerifyBody(
        "Sun",
        integ.state.coord[0],
        0.0027, 1.4e-07,
        8.413897447213910E-03,  9.903789233113485E-04, -2.534903517530939E-04,
        -1.318249797153773E-06,  8.245353286641178E-06, -2.091998687543678E-08
    )) return 1;

    if (VerifyBody(
        "Jupiter",
        integ.state.coord[1],
        0.067, 9.0e-05,
        -5.371058248261468E+00, -8.719764486851470E-01,  1.236123797189617E-01,
        1.116401240893049E-03, -7.098409700584834E-03,  4.923287714072183E-06
    )) return 1;

    if (VerifyBody(
        "Saturn",
        integ.state.coord[2],
        0.011, 7.7e-06,
        -9.152012613180709E+00, -3.051457308928910E+00,  4.178207655560304E-01,
        1.458944935312708E-03, -5.302701887770315E-03,  3.316613855544254E-05
    )) return 1;

    if (VerifyBody(
        "Uranus",
        integ.state.coord[3],
        0.0031, 1.1e-07,
        1.888005010509793E+01,  6.532300850250829E+00, -2.200617931964084E-01,
        -1.315257453228866E-03,  3.534033021483668E-03,  3.020648507746295E-05
    )) return 1;

    if (VerifyBody(
        "Neptune",
        integ.state.coord[4],
        0.047, 5.0e-06,
        -2.904645113602549E+01,  8.251202310543533E+00,  4.995306486561001E-01,
        -8.755137439452106E-04, -2.999579728759390E-03,  8.225104197814640E-05
    )) return 1;

    printf("SolarSystem: PASS\n");
    return 0;
}


//----------------------------------------------------------------------------------------------------
// Test RungeKutta::Simulator, which is designed for large or variable-sized states.


using catenary_state_t = std::vector<body_state_t>;
static constexpr double CatSpan = 10.0;        // meters separating the two anchors
static constexpr double CatLink = 0.25;        // rest length of spring connecting each link


struct CatenaryDeriv
{
    vec_t anchor1 { 0.0, 0.0, 0.0 };
    vec_t anchor2 { CatSpan, 0.0, 0.0 };

    void operator() (catenary_state_t& slope, const catenary_state_t& state)
    {
        // The catenary state includes a chain of mobile particles.
        // There are spring forces between consecutive particles,
        // and a single gravity acceleration that acts on all of them.

        const std::size_t n = state.size();
        assert(n == slope.size());
        for (std::size_t i=0; i < n; ++i)
        {
            vec_t pos1 = (i == 0)   ? anchor1 : state.at(i-1).pos;
            vec_t pos2 = (i == n-1) ? anchor2 : state.at(i+1).pos;
            (void)pos1;
            (void)pos2;
            // spring force
            // gravity
        }
    }
};


static int Catenary()
{
    CatenaryDeriv deriv;
    RungeKutta::ListAdd<body_state_t> add;
    RungeKutta::ListMul<body_state_t, double> mul;
    using catenary_sim_t = RungeKutta::Simulator<double, catenary_state_t, decltype(deriv), decltype(add), decltype(mul)>;
    catenary_sim_t sim(deriv, add, mul);

    const std::size_t nparticles = 50;
    sim.resize(nparticles);

    // Set up initial state: particles at locations around the line from anchor1 to anchor2.
    for (std::size_t i = 0; i < nparticles; ++i)
    {
        body_state_t& s = sim.state.at(i);

        // Linear interpolation between the two endpoints.
        // But there are 2+nparticles including the anchors conceptually at indexes [-1] and [nparticles].
        double frac = (i + 1.0) / (nparticles + 1);
        s.pos = deriv.anchor1 + frac*(deriv.anchor2 - deriv.anchor1);
        s.vel = vec_t{};
    }

    const double dt = 0.1;
    sim.step(dt);

    printf("Catenary: PASS\n");
    return 0;
}
