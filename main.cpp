#include "sdl3d.h"

using namespace SDL3D;

int main(int argc, char** argv)
{
    SDL_Init(SDL_INIT_EVERYTHING);

    screen = SDL_SetVideoMode(640,480,24,SDL_SWSURFACE);

    st.wireframe = true;
    st.textured = false;

    triangle tri[2];
    tri[0].p[0].x = -0.1;
    tri[0].p[0].y = -0.1;
    tri[0].p[0].z = 0;
    tri[0].p[1].x = 0.1;
    tri[0].p[1].y = 0;
    tri[0].p[1].z = 0;
    tri[0].p[2].x = 0;
    tri[0].p[2].y = 0.1;
    tri[0].p[2].z = 0;
    tri[0].r = 255;
    tri[0].g = 255;
    tri[0].b = 255;
    tri[0].uv[0][0] = -1;
    tri[0].uv[1][0] = -1;
    tri[0].uv[0][1] = 1;
    tri[0].uv[1][1] = 0;
    tri[0].uv[0][2] = 0;
    tri[0].uv[1][2] = 1;

    tri[0].texture = SDL_LoadBMP("tex.bmp");

    /*tri[1].p[0].x = -0.1;
    tri[1].p[0].y = -0.1;
    tri[1].p[0].z = 1;
    tri[1].p[1].x = 0.1;
    tri[1].p[1].y = -0.1;
    tri[1].p[1].z = 1;
    tri[1].p[2].x = 0;
    tri[1].p[2].y = 0.1;
    tri[1].p[2].z = 1;
    tri[1].r = 255;
    tri[1].g = 0;
    tri[1].b = 0;*/

    tri[0].perspective(-1,0,0);
    //tri[1].perspective(-1,0,0);

    for(int i = 0;i<100;i++)
    {
        SDL_FillRect(screen,&screen->clip_rect,SDL_MapRGB(screen->format,0,0,0));

        SDL_LockSurface(screen);
        SDL_LockSurface(tri[0].texture);
        //tri[1].render();
        //tri[1].rotz(0.05);
        tri[0].render();
        tri[0].rotz(0.05);

        tri[0].roty(0.1);
        //tri[1].roty(0.1);

        SDL_UnlockSurface(screen);
        SDL_UnlockSurface(tri[0].texture);

        SDL_Flip(screen);

        SDL_Delay(10);
    }

    SDL_FreeSurface(tri[0].texture);

    SDL_Quit();

    return 0;
}
