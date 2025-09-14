#pragma once
#include <functional>

namespace CosineKitty
{
    template <typename real_t, typename state_t>
    class Integrator
    {
    public:
        using deriv_func_t = std::function<state_t(const state_t&)>;
        deriv_func_t deriv;
        state_t state{};

        explicit Integrator(deriv_func_t _deriv)
            : deriv(_deriv)
            {}

        void step(real_t dt)
        {
            state_t k1 = deriv(state);
            state_t k2 = deriv(state + (dt/2)*k1);
            state_t k3 = deriv(state + (dt/2)*k2);
            state_t k4 = deriv(state + dt*k3);
            state += (dt/6)*(k1 + 2*(k2+k3) + k4);
        }
    };

    // `StateVector` is handy for use as a `state_t` with Integrator,
    // but it is NOT REQUIRED. You can use `double` if your model has a 1D state.
    // But most of the time you can use StateVector to express any compound state as
    // an ordered list of numbers.
    // For example, suppose 4 particles in 3D space have position and velocity (two 3D vectors):
    // dimensions = 4*3*2 = 24
    // The aggregate system state is a vector in a 24-dimensional space:
    // [
    //      rx0, ry0 rz0,  rx1, ry1, rz1,  rx2, ry2, rz2,  rx3, ry3, rz3,
    //      vx0, vy0 vz0,  vx1, vy1, vz1,  vx2, vy2, vz2,  vx3, vy3, vz3
    // ]
    // This is one example, and you don't have to follow it as a convention.
    // It is up to the caller to understand what the array index values mean,
    // and it is an arbitrary mapping for each use case.
    // This code tries to remain as unopinionated as possible in order
    // to remain flexibly reusable for a variety of physical systems.

    template <unsigned dimensions, typename coord_t, typename real_t>
    struct StateVector
    {
        coord_t coord[dimensions]{};

        void operator += (const StateVector& other)
        {
            for (unsigned i = 0; i < dimensions; ++i)
                coord[i] += other.coord[i];
        }

        friend StateVector operator * (real_t scalar, const StateVector& vec)
        {
            StateVector product;
            for (unsigned i = 0; i < dimensions; ++i)
                product.coord[i] = scalar * vec.coord[i];
            return product;
        }

        friend StateVector operator + (const StateVector& a, const StateVector& b)
        {
            StateVector sum;
            for (unsigned i = 0; i < dimensions; ++i)
                sum.coord[i] = a.coord[i] + b.coord[i];
            return sum;
        }
    };
}
