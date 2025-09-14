#include <cstdio>
#include <cstring>
#include <cmath>
#include "rk4_integrator.hpp"

using test_func_t = int (*)();

struct Test
{
    const char *name;
    test_func_t func;
};


static int Logarithm();
static int Pendulum();


static Test TestList[] =
{
    { "log", Logarithm },
    { "pendulum", Pendulum }
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
    using state_t = CosineKitty::StateVector<2, double>;

    // Use the numerical integrator to estimate ln(2).
    // Let t be the independent variable that ranges [1, 2].
    // Let v = 1/t be the "velocity" or slope.
    // Let x = integral(v) be the "position" or integral of v over the range t=1 to t=2.
    // The state necessarily includes t because it is needed to calculate v.

    auto deriv = [](const state_t& state)
    {
        state_t m;
        m.coord[0] = 1;                     // dt/dt = 1
        m.coord[1] = 1 / state.coord[0];    // v = dx/dt = 1/t
        return m;
    };

    CosineKitty::Integrator<double, state_t> integ(deriv);

    integ.state.coord[1] = 0;
    const int nSteps = 200;
    const double dt = 1.0 / nSteps;
    for (int n = 0; n < nSteps; ++n)
    {
        integ.state.coord[0] = 1 + dt*n;
        integ.step(dt);
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
    using state_t = CosineKitty::StateVector<2, double>;
    using integ_t = CosineKitty::Integrator<double, state_t>;

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
    while (zeroCrossingCount < zeroCrossingLimit && t < 600.0)
    {
        if (Verbose)
            printf("Pendulum: t=%0.2lf, theta=%0.6lf, omega=%0.6lf\n", t, integ.state.coord[0], integ.state.coord[1]);

        const double prevAngle = integ.state.coord[0];
        integ.step(dt);
        const double angle = integ.state.coord[0];
        if (angle * prevAngle <= 0)
        {
            ++zeroCrossingCount;
            double crossingTime = t; // FIXFIXFIX: Interpolate how far between samples the zero-crossing occurred.
            if (prevCrossingTime >= 0)
            {
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

    static constexpr double largeAngleCorrection = 1 + A*A/16 + A*A*A*A*(11.0/3072);
    printf("          large angle correction factor = %0.6lf\n", largeAngleCorrection);
    const double meanCrossingTime = periodTimeSum / zeroCrossingCount;
    const double expectedCrossingTime = M_PI * std::sqrt(L/g) * largeAngleCorrection;       // half-trip time
    const double diff = std::abs(expectedCrossingTime - meanCrossingTime);
    printf("          %d zero-crossings\n", zeroCrossingCount);
    printf("          mean interval = %0.6lf seconds\n", meanCrossingTime);
    printf("          expected      = %0.6lf seconds\n", expectedCrossingTime);
    printf("          diff          = %g\n", diff);
    printf("          min           = %0.6lf seconds\n", minCrossingTime);
    printf("          max           = %0.6lf seconds\n", maxCrossingTime);
    printf("          spread        = %g\n", maxCrossingTime - minCrossingTime);
    printf("Pendulum: PASS\n");
    return 0;
}
