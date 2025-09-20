#pragma once
#include "rk4_simulator.hpp"
#include <cassert>
#include <stdexcept>

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
    private:
        const std::size_t mobileCount;
        const std::size_t anchorCount;

    public:
        explicit MeshDeriv(std::size_t _mobileCount, std::size_t _anchorCount)
            : mobileCount(_mobileCount)
            , anchorCount(_anchorCount)
            {}

        void operator() (mesh_list_t& slope, const mesh_list_t& state)
        {
            const std::size_t numParticles = mobileCount + anchorCount;
            assert(numParticles == slope.size());
            assert(numParticles == state.size());

            for (std::size_t i = 0; i < numParticles; ++i)
            {
                slope[i].pos = state[i].vel;
                slope[i].vel = MeshVector(0, 0, 0);
            }
        }
    };


    using mesh_base_t = ListSimulator<double, MeshParticle, MeshDeriv>;


    class MeshSimulator : public mesh_base_t
    {
    public:
        explicit MeshSimulator(std::size_t mobileCount, std::size_t anchorCount)
            : mesh_base_t(MeshDeriv(mobileCount, anchorCount), mobileCount + anchorCount)
            {}
    };


    class RibbonSimulator : public MeshSimulator
    {
    public:
        static constexpr double HorizontalSpacing = 0.01;  // distance between columns [m]
        static constexpr double VerticalSpacing = 0.01;    // distance between rows [m]
        static constexpr std::size_t MobileColumns = 13;
        static constexpr std::size_t RibbonColumns = 2 + MobileColumns;     // anchors on both ends of the ribbon
        static constexpr std::size_t RibbonRows = 3;
        static constexpr std::size_t ParticleCount = RibbonRows * RibbonColumns;
        static constexpr std::size_t MobileCount   = RibbonRows * MobileColumns;
        static constexpr std::size_t AnchorCount   = 2*RibbonRows;

        explicit RibbonSimulator()
            : MeshSimulator(MobileCount, AnchorCount)
        {
            makeRibbon();
        }

    private:
        static std::size_t index(std::size_t col, std::size_t row)
        {
            // Diagram:
            //
            // col: 0    1    2    3    4    5
            //
            //      A----M----M----M----M----A   row: 2
            //           |    |    |    |
            //           |    |    |    |
            //      A----M----M----M----M----A   row: 1
            //           |    |    |    |
            //           |    |    |    |
            //      A----M----M----M----M----A   row: 0
            //
            //
            // A = anchor
            // M = mobile particle

            if (col >= RibbonColumns)
                throw std::out_of_range("Invalid mesh column: " + std::to_string(col) + " >= " + std::to_string(RibbonColumns));

            if (row >= RibbonRows)
                throw std::out_of_range("Invalid mesh row: " + std::to_string(row) + " >= " + std::to_string(RibbonRows));

            // Cooperate with MeshDeriv to make things run faster and to make the code simpler.
            // Put all mobile balls at the front of the list, and all anchor balls at the back.
            // Anything that wants to iterate over either or both has a very simple contiguous range.

            int x;
            if (col == 0)
            {
                // The column of left anchors goes after mobile balls and right anchors in the array.
                x = RibbonColumns - 1;
            }
            else
            {
                // Everything else shifts left by one column.
                // This puts the mobile balls at the front.
                x = col - 1;
            }

            return row + RibbonRows*x;
        }

        MeshParticle& particle(std::size_t col, std::size_t row)
        {
            return state.at(index(col, row));
        }

        void makeRibbon()
        {
            for (std::size_t col = 0; col < RibbonColumns; ++col)
            {
                for (std::size_t row = 0; row < RibbonRows; ++row)
                {
                    MeshParticle&p = particle(col, row);
                    p.pos = MeshVector(col * HorizontalSpacing, row * VerticalSpacing, 0);
                    p.vel = MeshVector(0, 0, 0);
                }
            }
        }
    };
}
