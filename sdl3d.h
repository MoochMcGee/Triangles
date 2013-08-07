#ifndef SDL3D_H
#define SDL3D_H

#include <cmath>
#include <SDL.h>

SDL_Surface* screen;

namespace SDL3D
{
    struct state
    {
        bool wireframe;
        bool textured;
    } st;
	struct Vector
	{
		float x,y,z,hv;
		Vector()
		{
		    hv = 1;
		}
		float length()
		{
			return sqrt(x*x+y*y+z*z);
		}
	};
	Vector vecadd(Vector a, Vector b)
	{
		Vector c;
		c.x = a.x + b.x;
		c.y = a.y + b.y;
		c.z = a.z + b.z;
		return c;
	}
	Vector vecsub(Vector a, Vector b)
	{
		Vector c;
		c.x = a.x - b.x;
		c.y = a.y - b.y;
		c.z = a.z - b.z;
		return c;
	}
	Vector scamul(Vector a, float b)
	{
		Vector c;
		c.x = a.x * b;
		c.y = a.y * b;
		c.z = a.z * b;
		return c;
	}
	Vector scadiv(Vector a, float b)
	{
		Vector c;
		c.x = a.x / b;
		c.y = a.y / b;
		c.z = a.z / b;
		return c;
	}
	Vector unit(Vector a)
	{
		Vector b;
		double c,d;
		c = a.x * a.x + a.y * a.y + a.z * a.z;

		long i = *(long*)&c;		//Evil, evil hack
		i = 0x5f375a86 - (i >> 1);	//WAT.
		d = *(float*)&i;		//Evil, evil hack part deux
		d = d * (1.5f - (c * 0.5f * d * d));

		b.x = a.x / d;
		b.y = a.y / d;
		b.z = a.z / d;

		return b;
	}
	double vecdot(Vector a, Vector b)
	{
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}
	double vecproj(Vector a, Vector b)
	{
		return vecdot(a,b)/b.length();
	}
	Vector veccross(Vector a, Vector b)
	{
		Vector c;
		c.x = a.y * b.z - a.z * b.y;
		c.y = a.z * b.x - a.x * b.z;
		c.z = a.x * b.y - a.y * b.x;
		return c;
	}
	struct matrix
	{
	    float e[4][4];
	};
	Vector matrixmul(matrix a, Vector b)
	{
	    Vector c;
	    c.x = a.e[0][0] * b.x + a.e[0][1] * b.y + a.e[0][2] * b.z + a.e[0][3] * b.hv;
	    c.y = a.e[1][0] * b.x + a.e[1][1] * b.y + a.e[1][2] * b.z + a.e[1][3] * b.hv;
	    c.z = a.e[2][0] * b.x + a.e[2][1] * b.y + a.e[2][2] * b.z + a.e[2][3] * b.hv;
	    c.hv = a.e[3][0] * b.x + a.e[3][1] * b.y + a.e[3][2] * b.z + a.e[3][3] * b.hv;
	    return c;
	}
	void putpix(int x, int y, int r, int g, int b)
	{
	    unsigned char* p = (unsigned char*)screen->pixels;
	    if(x < 0 || x>=screen->w) return;
	    if(y < 0 || y>=screen->h) return;
	    p[(((y*screen->w) + x) * 3) + 0] = b;
	    p[(((y*screen->w) + x) * 3) + 1] = g;
	    p[(((y*screen->w) + x) * 3) + 2] = r;
	}
	void renderline(int x0, int y0, int x1, int y1, int r, int g, int b)
	{
	    int x = x0;
	    int y = y0;
	    int dx = abs(x1-x0);
	    int dy = abs(y1-y0);
	    int sx = (x0 < x1) ? 1 : -1;
	    int sy = (y0 < y1) ? 1 : -1;
	    int err = dx-dy;
	    for(;;)
        {
            putpix(x,y,r,g,b);
            if(x == x1 && y == y1) break;
            int e2 = err << 1;
            if(e2 >= -dy)
            {
                err -= dy;
                x += sx;
            }
            if(e2 < dx)
            {
                err += dx;
                y += sy;
            }
        }
	}
	void floodfill(int x, int y, int r, int g, int b)
	{
        {
            unsigned char* p = (unsigned char*)screen->pixels;
            if(x < 0 || x>=screen->w) return;
            if(y < 0 || y>=screen->h) return;
            unsigned char b0 = p[(((y*screen->w) + x) * 3) + 0];
            unsigned char g0 = p[(((y*screen->w) + x) * 3) + 1];
            unsigned char r0 = p[(((y*screen->w) + x) * 3) + 2];

            unsigned int col = (b << 16) | (g << 8) | r;
            unsigned int col0 = (b0 << 16) | (g0 << 8) | r0;
            if(col0!=col)
            {
                putpix(x,y,r,g,b);
                floodfill(x-1,y,r,g,b);
                floodfill(x+1,y,r,g,b);
                floodfill(x,y-1,r,g,b);
                floodfill(x,y+1,r,g,b);
            }
        }
	}
	void rendertexturedtri(int x0, int y0, int x1, int y1, int x2, int y2, int u0, int v0, int u1, int v1, int u2, int v2, SDL_Surface* texture)
	{
	    if(y0 != y1 && y0 != y2 && texture != nullptr)
        {
                float startx = x0;
                float endx = x0 + ((x1-x0)/(y1-y0));

                if(startx < 0) startx = 0;
                if(endx > screen->w) endx = screen->w;

                int u = u0;
                int v = v0;

                int r, g, b;

                for(int i = y0;i<=y1;i++)
                {
                    for(int j = startx;j<=endx;j++)
                    {
                        uint8_t* p = (uint8_t*)texture->pixels;
                        //r = p[((u + (v*texture->w))*3)+0];
                        //g = p[((u + (v*texture->w))*3)+1];
                        //b = p[((u + (v*texture->w))*3)+2];
                        r = 255;
                        g = 0;
                        b = 0;
                        putpix(j,i,r,g,b);
                        u++;
                    }
                    startx -= ((x2-x0)/(y2-y0));
                    endx += ((x1-x0)/(y1-y0));
                }

                startx = x2;
                endx = x2;

                for(int i = y2;i>y1;i--)
                {
                    for(int j = startx;j<=endx;j++)
                    {
                        uint8_t* p = (uint8_t*)texture->pixels;
                        //b = p[((u + (v*texture->w))*3)+0];
                        //g = p[((u + (v*texture->w))*3)+1];
                        //r = p[((u + (v*texture->w))*3)+2];
                        r = 255;
                        g = 255;
                        b = 255;
                        putpix(j,i,r,g,b);
                        u++;
                    }
                    startx += ((x2-x0)/(y2-y0));
                    endx += ((x2-x1)/(y2-y1));
                }
        }
	}

	struct triangle
	{
	    Vector p[3];
	    int r,g,b;
	    float uv[2][3];
	    SDL_Surface* texture;
	    void translate(float x, float y, float z)
	    {
	        matrix tm;
	        float tme[4][4] = {{1,0,0,x},{0,1,0,y},{0,0,1,z},{0,0,0,1}};
	        for(int i = 0;i<4;i++)
            {
                for(int j = 0;j<4;j++)
                {
                    tm.e[i][j] = tme[i][j];
                }
            }
	        p[0] = matrixmul(tm,p[0]);
	        p[1] = matrixmul(tm,p[1]);
	        p[2] = matrixmul(tm,p[2]);
	    }
	    void rotx(float angle)
	    {
	        matrix tm;
	        float tme[4][4] = {{1,0,0,0},{0,cos(angle),-sin(angle),0},{0,sin(angle),cos(angle),0},{0,0,0,1}};
	        for(int i = 0;i<4;i++)
            {
                for(int j = 0;j<4;j++)
                {
                    tm.e[i][j] = tme[i][j];
                }
            }
	        p[0] = matrixmul(tm,p[0]);
	        p[1] = matrixmul(tm,p[1]);
	        p[2] = matrixmul(tm,p[2]);
	    }
	    void roty(float angle)
	    {
	        matrix tm;
	        float tme[4][4] = {{cos(angle),0,sin(angle),0},{0,1,0,0},{-sin(angle),0,cos(angle),0},{0,0,0,1}};
	        for(int i = 0;i<4;i++)
            {
                for(int j = 0;j<4;j++)
                {
                    tm.e[i][j] = tme[i][j];
                }
            }
	        p[0] = matrixmul(tm,p[0]);
	        p[1] = matrixmul(tm,p[1]);
	        p[2] = matrixmul(tm,p[2]);
	    }
	    void rotz(float angle)
	    {
	        matrix tm;
	        float tme[4][4] = {{cos(angle),-sin(angle),0,0},{sin(angle),cos(angle),0,0},{0,0,1,0},{0,0,0,1}};
	        for(int i = 0;i<4;i++)
            {
                for(int j = 0;j<4;j++)
                {
                    tm.e[i][j] = tme[i][j];
                }
            }
	        p[0] = matrixmul(tm,p[0]);
	        p[1] = matrixmul(tm,p[1]);
	        p[2] = matrixmul(tm,p[2]);
	    }
	    void scale(float x, float y, float z)
	    {
	        matrix tm;
	        float tme[4][4] = {{x,0,0,0},{0,y,0,0},{0,0,z,0},{0,0,0,1}};
	        for(int i = 0;i<4;i++)
            {
                for(int j = 0;j<4;j++)
                {
                    tm.e[i][j] = tme[i][j];
                }
            }
	        p[0] = matrixmul(tm,p[0]);
	        p[1] = matrixmul(tm,p[1]);
	        p[2] = matrixmul(tm,p[2]);
	    }
	    void perspective(float near, float ex, float ey)
	    {
	        matrix tm;
	        float tme[4][4] = {{1,0,-(ex/near),0},{0,1,-(ey/near),0},{0,0,1,0},{0,0,1/near,0}};
	        for(int i = 0;i<4;i++)
            {
                for(int j = 0;j<4;j++)
                {
                    tm.e[i][j] = tme[i][j];
                }
            }
	        p[0] = matrixmul(tm,p[0]);
	        p[1] = matrixmul(tm,p[1]);
	        p[2] = matrixmul(tm,p[2]);
	        scadiv(p[0],p[0].hv);
	        scadiv(p[1],p[1].hv);
	        scadiv(p[2],p[2].hv);
	    }
	    void render()
	    {
            float x0 = ((p[0].x + 1)/2 * (screen->w-1));
            float y0 = ((p[0].y + 1)/2 * (screen->h-1));
            float x1 = ((p[1].x + 1)/2 * (screen->w-1));
            float y1 = ((p[1].y + 1)/2 * (screen->h-1));
            float x2 = ((p[2].x + 1)/2 * (screen->w-1));
            float y2 = ((p[2].y + 1)/2 * (screen->h-1));

            for(int i = 0;i<2;i++)
            {
                for(int j = 0;j<3;j++)
                {
                    if(uv[i][j] < -1) uv[i][j] = -1;
                    if(uv[i][j] > 1) uv[i][j] = 1;
                }
            }

            float u0 = ((uv[0][0] + 1)/2) * (texture->w-1);
            float v0 = ((uv[1][0] + 1)/2) * (texture->h-1);
            float u1 = ((uv[0][1] + 1)/2) * (texture->w-1);
            float v1 = ((uv[1][1] + 1)/2) * (texture->h-1);
            float u2 = ((uv[0][2] + 1)/2) * (texture->w-1);
            float v2 = ((uv[1][2] + 1)/2) * (texture->h-1);

            if(!st.textured)
            {
                renderline(x0,y0,x1,y1,r,g,b);
                renderline(x0,y0,x2,y2,r,g,b);
                renderline(x1,y1,x2,y2,r,g,b);
                if(!st.wireframe)
                {
                    float avgx = (x0 + x1 + x2) / 3;
                    float avgy = (y0 + y1 + y2) / 3;

                    floodfill(avgx,avgy,r,g,b);
                }
            }
            else
            {
                rendertexturedtri(x0,y0,x1,y1,x2,y2,u0,v0,u1,v1,u2,v2,texture);
            }
	    }
	};
};

#endif //SDL3D_H
