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
    const RungeKutta::RibbonSimulator& ribbon;

    explicit RenderContext(const RungeKutta::RibbonSimulator& _ribbon)
        : ribbon(_ribbon)
    {
        using namespace RungeKutta;

        const MeshParticle& p0 = ribbon.state.at(0);
        double xmin = p0.pos.x;
        double xmax = p0.pos.x;
        double ymin = p0.pos.y;
        double ymax = p0.pos.y;
        for (const MeshParticle& p : ribbon.state)
        {
            xmin = std::min(xmin, p.pos.x);
            xmax = std::max(xmax, p.pos.x);
            ymin = std::min(ymin, p.pos.y);
            ymax = std::max(ymax, p.pos.y);
        }

        xCenter = (xmin + xmax) / 2;
        yCenter = (ymin + ymax) / 2;
    }

    int xScreen(float x) const
    {
        return (screenWidth/2) + static_cast<int>(std::round(zoom * (x - xCenter)));
    }

    int yScreen(float y) const
    {
        return (screenHeight/2) - static_cast<int>(std::round(zoom * (y - yCenter)));
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
            int x = xScreen(p.pos.x);
            int y = yScreen(p.pos.y);
            //printf("(%d, %d)\n", x, y);
            DrawCircleGradient(x, y, ballRadius, GREEN, DARKGREEN);
        }
    }
};


int main(int argc, const char *argv[])
{
    // This program does NOT generate audio,
    // but it is a feasibility study for audio generation.
    const double sampleRate = 48000;
    const double dt = 1 / sampleRate;

    RungeKutta::RibbonSimulator ribbon;
    RenderContext render(ribbon);
    InitWindow(render.screenWidth, render.screenHeight, "Ribbon Mesh");
    SetTargetFPS(80);
    while (!WindowShouldClose())
    {
        BeginDrawing();
        ClearBackground(BLACK);
        render.draw();
        EndDrawing();
        ribbon.step(dt);
    }
    CloseWindow();

    return 0;
}
