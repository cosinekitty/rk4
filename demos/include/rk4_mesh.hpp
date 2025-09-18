#pragma once
#include "rk4_simulator.hpp"

namespace RungeKutta
{
    struct MeshVector
    {
        double x{};
        double y{};
        double z{};
    };


    struct MeshParticle
    {
        MeshVector pos;
        MeshVector vel;
    };


    using mesh_list_t = std::vector<MeshParticle>;


    class MeshDeriv
    {
    public:
        void operator() (const mesh_list_t& slope, const mesh_list_t& state)
        {

        }
    };


    using mesh_base_t = ListSimulator<double, MeshVector, MeshDeriv>;


    class MeshSimulator : mesh_base_t
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

        RibbonSimulator()
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
