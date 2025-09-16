#pragma once

namespace RungeKutta
{
    // Simulator is like Integrator, only designed for large or
    // variable-sized states. No return-by-value; instead, the slope
    // vectors k1..k4 are members, and they are updated in place as needed.
    // This makes Simulator less functional-oriented, but eliminates
    // constantly having to construct state_t objects.

    template <typename real_t, typename state_t, typename deriv_proc_t>
    class Simulator
    {
    private:
        state_t k1, k2, k3, k4;     // slope vectors
        state_t w;                  // work register
        deriv_proc_t deriv;

        static void extrapolate(state_t& target, const state_t& origin, const state_t& slope, real_t dt)
        {
            // target = origin + slope*dt
            mul(target, slope, dt);
            add(target, origin);
        }

    public:
        state_t state{};

        explicit ProcIntegrator(deriv_proc_t _deriv)
            : deriv(_deriv)
            {}

        void step(real_t dt)
        {
            deriv(k1, state);

            mul(w, k1, dt/2);
            add(w, w, state);
            deriv(k2, w);

            mul(w, k2, dt/2);
            add(w, w, state);
            deriv(k3, w);

            mul(w, k3, dt);
            add(w, w, state);
            deriv(k4, w);

            add(w, k2, k3);
            mul(w, w, 2);
            add(w, w, k1);
            add(w, w, k4);
            mul(w, w, dt/6);
            add(state, state, w);
        }
    };
}
