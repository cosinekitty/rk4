#pragma once
#include <functional>

namespace CosineKitty
{
    template <typename real_t, typename state_t>
    class Integrator
    {
    public:
        using deriv_func_t = std::function<state_t(real_t, const state_t&)>;
        deriv_func_t deriv;
        state_t state{};

        explicit Integrator(deriv_func_t _deriv)
            : deriv(_deriv)
            {}

        void step(real_t t, real_t dt)
        {
            state_t k1 = deriv(t, state);
            state_t k2 = deriv(t + dt/2, state + (dt/2)*k1);
            state_t k3 = deriv(t + dt/2, state + (dt/2)*k2);
            state_t k4 = deriv(t + dt, state + dt*k3);
            state += (dt/6)*(k1 + 2*(k2+k3) + k4);
        }
    };
}
