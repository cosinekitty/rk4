#include <cstdio>
#include <cmath>
#include <raylib.h>
#include "rk4_mesh.hpp"


struct RenderContext
{
    const int screenWidth  = 1200;
    const int screenHeight =  600;
    float zoom = 7000;     // pixels per meter
    float xCenter{};
    float yCenter{};
    float zCenter{};
    float xRotation = 60.0;      // rotation around the x-axis, in degrees
    const RungeKutta::RibbonSimulator& ribbon;

    static constexpr float dcos(float deg)
    {
        return std::cos(static_cast<float>(M_PI/180) * deg);
    }

    static constexpr float dsin(float deg)
    {
        return std::sin(static_cast<float>(M_PI/180) * deg);
    }

    explicit RenderContext(const RungeKutta::RibbonSimulator& _ribbon)
        : ribbon(_ribbon)
    {
        using namespace RungeKutta;

        const MeshParticle& p0 = ribbon.state.at(0);
        double xmin = p0.pos.x;
        double xmax = p0.pos.x;
        double ymin = p0.pos.y;
        double ymax = p0.pos.y;
        double zmin = p0.pos.z;
        double zmax = p0.pos.z;
        for (const MeshParticle& p : ribbon.state)
        {
            xmin = std::min(xmin, p.pos.x);
            xmax = std::max(xmax, p.pos.x);
            ymin = std::min(ymin, p.pos.y);
            ymax = std::max(ymax, p.pos.y);
            zmin = std::min(zmin, p.pos.z);
            zmax = std::min(zmax, p.pos.z);
        }
        xCenter = (xmin + xmax) / 2;
        yCenter = (ymin + ymax) / 2;
        zCenter = (zmin + zmax) / 2;
    }

    void getScreenCoords(int& hor, int& ver, const RungeKutta::MeshParticle& p) const
    {
        const float c = dcos(xRotation);
        const float s = dsin(xRotation);
        const float dx = p.pos.x - xCenter;
        const float dy = p.pos.y - yCenter;
        const float dz = p.pos.z - zCenter;
        const float x = dx;
        const float y = c*dy + s*dz;
        const float z = c*dz - s*dy;
        const float denom = 0.5;
        const float pers = (denom + z) / denom;
        hor = (screenWidth /2) + static_cast<int>(std::round(zoom * pers * x));
        ver = (screenHeight/2) - static_cast<int>(std::round(zoom * pers * y));
    }

    int scale(float r) const
    {
        return static_cast<int>(std::round(zoom * r));
    }

    void draw()
    {
        using namespace RungeKutta;

        const int ballRadius = scale(0.0005);
        for (const MeshParticle& p : ribbon.state)
        {
            int hor, ver;
            getScreenCoords(hor, ver, p);
            DrawCircleGradient(hor, ver, ballRadius, GREEN, DARKGREEN);
        }

        int ah, av, bh, bv;
        for (const MeshSpring& s : ribbon.springs)
        {
            const MeshParticle& a = ribbon.state.at(s.ia);
            const MeshParticle& b = ribbon.state.at(s.ib);
            getScreenCoords(ah, av, a);
            getScreenCoords(bh, bv, b);
            DrawLine(ah, av, bh, bv, YELLOW);
        }
    }
};


int main(int argc, const char *argv[])
{
    // This program does NOT generate audio,
    // but it *is* a feasibility study for audio generation.
    const double sampleRate = 48000;
    const double dt = 1 / sampleRate;
    const double degreesPerAnimationFrame = 0.0;

    RungeKutta::RibbonSimulator ribbon;
    //ribbon.particle(2, 0).pos.z = 0.007;
    ribbon.particle(2, 0).vel.z = 0.007;
    ribbon.deriv.gravity.z = -9.8;
    RenderContext render(ribbon);
    ribbon.decayHalfLife = 0.25;

    InitWindow(render.screenWidth, render.screenHeight, "Ribbon Mesh");
    SetTargetFPS(200);
    while (!WindowShouldClose())
    {
        BeginDrawing();
        ClearBackground(BLACK);
        render.draw();
        EndDrawing();
        render.xRotation += degreesPerAnimationFrame;
        ribbon.update(dt);
    }
    CloseWindow();

    return 0;
}
