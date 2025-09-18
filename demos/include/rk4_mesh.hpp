#pragma once
#include "rk4_simulator.hpp"

namespace RungeKutta
{
    struct MeshVector
    {
        double x{};
        double y{};
        double z{};

        explicit MeshVector() {}

        explicit MeshVector(double _x, double _y, double _z)
            : x(_x)
            , y(_y)
            , z(_z)
            {}

        friend MeshVector operator * (double factor, const MeshVector& vec)
        {
            return MeshVector(factor*vec.x, factor*vec.y, factor*vec.z);
        }

        friend MeshVector operator + (const MeshVector& a, const MeshVector& b)
        {
            return MeshVector(a.x + b.x, a.y + b.y, a.z + b.z);
        }
    };


    struct MeshParticle
    {
        MeshVector pos;
        MeshVector vel;

        explicit MeshParticle() {}

        explicit MeshParticle(const MeshVector& _pos, const MeshVector& _vel)
            : pos(_pos)
            , vel(_vel)
            {}

        friend MeshParticle operator * (double factor, const MeshParticle& state)
        {
            return MeshParticle(factor * state.pos, factor * state.vel);
        }

        friend MeshParticle operator + (const MeshParticle& a, const MeshParticle& b)
        {
            return MeshParticle(a.pos + b.pos, a.vel + b.vel);
        }
    };


    using mesh_list_t = std::vector<MeshParticle>;


    class MeshDeriv
    {
    public:
        void operator() (mesh_list_t& slope, const mesh_list_t& state)
        {
            const std::size_t n = slope.size();
            if (n == state.size())
            {
                for (std::size_t i = 0; i < n; ++i)
                {
                    slope[i].pos = state[i].vel;
                    slope[i].vel = MeshVector(0, 0, 0);
                }
            }
        }
    };


    using mesh_base_t = ListSimulator<double, MeshParticle, MeshDeriv>;


    class MeshSimulator : public mesh_base_t
    {
    private:

    public:
        explicit MeshSimulator(std::size_t particleCount)
            : mesh_base_t(MeshDeriv(), particleCount)
            {}
    };


    class RibbonSimulator : public MeshSimulator
    {
    public:
        static constexpr std::size_t MobileLength = 13;
        static constexpr std::size_t RibbonLength = 2 + MobileLength;     // anchors on both ends of the ribbon
        static constexpr std::size_t RibbonWidth = 3;
        static constexpr std::size_t ParticleCount = RibbonWidth * RibbonLength;
        static constexpr std::size_t MobileCount   = RibbonWidth * MobileLength;

        explicit RibbonSimulator()
            : MeshSimulator(ParticleCount)
        {
            makeRibbon();
        }

    private:
        void makeRibbon()
        {
            // Consider this example (smaller than the real ribbon):
            //
            //     A----M----M----M----M----A
            //          |    |    |    |
            //     A----M----M----M----M----A
            //          |    |    |    |
            //     A----M----M----M----M----A
            //
            // A = anchor
            // M = mobile particle
        }
    };
}
