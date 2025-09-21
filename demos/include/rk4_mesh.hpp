#pragma once
#include "rk4_simulator.hpp"
#include <cassert>
#include <stdexcept>

namespace RungeKutta
{
    using uint = std::size_t;

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
        const uint mobileCount;
        const uint anchorCount;

    public:
        explicit MeshDeriv(uint _mobileCount, uint _anchorCount)
            : mobileCount(_mobileCount)
            , anchorCount(_anchorCount)
            {}

        void operator() (mesh_list_t& slope, const mesh_list_t& state)
        {
            const uint numParticles = mobileCount + anchorCount;
            assert(numParticles == slope.size());
            assert(numParticles == state.size());

            for (uint i = 0; i < numParticles; ++i)
            {
                slope[i].pos = state[i].vel;
                slope[i].vel = MeshVector(0, 0, 0);
            }
        }
    };


    using mesh_base_t = ListSimulator<double, MeshParticle, MeshDeriv>;


    struct MeshSpring
    {
        const uint ia;
        const uint ib;

        explicit MeshSpring()
            : ia(-1)
            , ib(-1)
            {}

        explicit MeshSpring(uint aIndex, uint bIndex)
            : ia(aIndex)
            , ib(bIndex)
            {}
    };


    using spring_list_t = std::vector<MeshSpring>;


    class MeshSimulator : public mesh_base_t
    {
    public:
        spring_list_t springs;

        explicit MeshSimulator(uint mobileCount, uint anchorCount)
            : mesh_base_t(MeshDeriv(mobileCount, anchorCount), mobileCount + anchorCount)
            {}
    };


    class RibbonSimulator : public MeshSimulator
    {
    public:
        static constexpr double HorizontalSpacing = 0.01;  // distance between columns [m]
        static constexpr double VerticalSpacing = 0.01;    // distance between rows [m]
        static constexpr uint MobileColumns = 13;
        static constexpr uint RibbonColumns = 2 + MobileColumns;     // anchors on both ends of the ribbon
        static constexpr uint RibbonRows = 3;
        static constexpr uint ParticleCount = RibbonRows * RibbonColumns;
        static constexpr uint MobileCount   = RibbonRows * MobileColumns;
        static constexpr uint AnchorCount   = 2*RibbonRows;

        explicit RibbonSimulator()
            : MeshSimulator(MobileCount, AnchorCount)
        {
            makeRibbon();
        }

        static constexpr bool valid(uint col, uint row)
        {
            return (col < RibbonColumns) && (row < RibbonRows);
        }

        static constexpr bool isAnchor(uint col, uint row)
        {
            return (col == 0) || (col+1 == RibbonColumns);
        }

        static constexpr bool isMobile(uint col, uint row)
        {
            return !isAnchor(col, row);
        }

    private:

        static constexpr uint index(uint col, uint row)
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

            // The column of left anchors goes after mobile balls and right anchors in the array.
            // Everything else shifts left by one column.
            // This puts all the mobile balls at the front as a contiguous block.
            int x = (col==0) ? (RibbonColumns-1) : (col-1);
            return row + RibbonRows*x;
        }

        MeshParticle& particle(uint col, uint row)
        {
            return state.at(index(col, row));
        }

        void addSpring(uint acol, uint arow, uint bcol, uint brow)
        {
            if (valid(acol,arow) && valid(bcol,brow))
            {
                // Never connect two anchors with a spring.
                // Such a spring has no purpose or effect; it just wastes memory and CPU.
                if (isMobile(acol,arow) || isMobile(bcol,brow))
                {
                    int aindex = index(acol, arow);
                    int bindex = index(bcol, brow);
                    springs.push_back(MeshSpring(aindex, bindex));
                }
            }
        }

        void makeRibbon()
        {
            // Intialize particle positions and velocies.
            for (uint col = 0; col < RibbonColumns; ++col)
            {
                for (uint row = 0; row < RibbonRows; ++row)
                {
                    MeshParticle& p = particle(col, row);
                    p.pos = MeshVector(col * HorizontalSpacing, row * VerticalSpacing, 0);
                    p.vel = MeshVector(0, 0, 0);

                    addSpring(col-1, row, col, row);
                    addSpring(col, row-1, col, row);
                }
            }
        }
    };
}
