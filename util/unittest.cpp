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


static Test TestList[] =
{
    { "log", Logarithm }
};


constexpr int NumTests = static_cast<int>(sizeof(TestList) / sizeof(TestList[0]));


int main(int argc, const char *argv[])
{
    if (argc == 2)
    {
        int rc;
        const char *verb = argv[1];
        if (!strcmp(verb, "all"))
        {
            for (int i = 0; i < NumTests; ++i)
            {
                printf("Running: %s\n", TestList[i].name);
                rc = TestList[i].func();
                if (rc)
                {
                    printf("FAIL: %s returned %d\n", TestList[i].name, rc);
                    return 1;
                }
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
            rc = test->func();
            if (rc)
            {
                printf("FAIL: %s returned %d\n", test->name, rc);
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
    // Use the numerical integrator to estimate ln(2).
    // Let t be the independent variable that ranges [1, 2].
    // Let v = 1/t be the "velocity" or slope.
    // Let x = integral(v) be the "position" or integral of v over the range t=1 to t=2.
    CosineKitty::Integrator<double, double> integ([](double t, const double& x){return 1/t;});

    integ.state = 0;
    const int nSteps = 200;
    const double dt = 1.0 / nSteps;
    for (int n = 0; n < nSteps; ++n)
        integ.step(1 + dt*n, dt);

    const double correct = std::log(2.0);
    const double diff = std::abs(correct - integ.state);
    const double tolerance = 1.3e-12;
    printf("Integral = %0.16lf\n", integ.state);
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
