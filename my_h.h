#ifndef MY_H
#define MY_H
# include <mlx.h>
# include <stdio.h>
# include <stdlib.h>
# include <unistd.h>
# include <math.h>
# define RGB(r, g, b)(256 * 256 * (int)(r) + 256 * (int)(g) + (int)(b))
# define WIDTH 800
# define HEIGHT 600

typedef struct                                                                  s_img
{
  double cos;
  void                            *mlx;
  void                            *win;
  unsigned long           img_color;
  char                            *data;
  void                            *img_ptr;
  int                                     sizeline;
  int                                     bpp;
  int                                     x;
  int                                     y;
  int                                     endian;
  void                            *mlx_ptr;
}                                                                                               t_img;


typedef struct                                                                  s_struct
{
  int                                     fracnum;
  void                            *mlx;
  void                            *win;
  t_img                           *img;
  int color;
  double decx;
  double decy;
  double decz;
}                                                                                               t_struct;

typedef struct s_vec3d
{
  double x;
  double y;
  double z;
  double v;
}		t_vec3d;

typedef struct s_ray
{
  t_vec3d *o;
  t_vec3d *d;
  double vx;
  double vy;
  double vz;
  double xeye;
  double yeye;
  double zeye;
}		t_ray;

typedef struct s_sphere
{
  t_vec3d *pos;
  double radius;
}		t_sphere;

typedef struct s_light
{
  double x;
  double y;
  double z;
  double dx;
  double dy;
  double dz;
  double ox;
  double oy;
  double oz;
}		t_light;
  

double mypow(int a, int n);
void my_pixel_put_to_image(t_img *myimg, int x, int y, int color);
#endif
