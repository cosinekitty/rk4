#include <cstdio>
#include <raylib.h>
#include "rk4_mesh.hpp"


struct RenderContext
{
    const int screenWidth  = 1200;
    const int screenHeight =  600;
    float zoom = 8000;     // pixels per meter
    float xCenter = 0.05;
    float yCenter = 0.00;

    void draw(RungeKutta::RibbonSimulator& ribbon)
    {
        
    }
};


int main(int argc, const char *argv[])
{
    // This program does NOT generate audio,
    // but it is a feasibility study for audio generation.
    const double sampleRate = 48000;
    const double dt = 1 / sampleRate;

    RungeKutta::RibbonSimulator ribbon;
    RenderContext render;
    InitWindow(render.screenWidth, render.screenHeight, "Ribbon Mesh");
    SetTargetFPS(240);
    while (!WindowShouldClose())
    {
        BeginDrawing();
        ClearBackground(BLACK);
        render.draw(ribbon);
        EndDrawing();
        ribbon.step(dt);
    }
    CloseWindow();

    return 0;
}
