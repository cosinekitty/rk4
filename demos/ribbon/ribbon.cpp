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
};


int main(int argc, const char *argv[])
{
    RungeKutta::RibbonSimulator ribbon;
    RenderContext render;
    InitWindow(render.screenWidth, render.screenHeight, "Ribbon Mesh");
    SetTargetFPS(240);
    while (!WindowShouldClose())
    {
        BeginDrawing();
        ClearBackground(BLACK);
        // render a frame
        EndDrawing();
        // update the simulation
    }
    CloseWindow();

    return 0;
}
