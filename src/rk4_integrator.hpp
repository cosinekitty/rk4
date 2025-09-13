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
        }
    };
}
