#include "my_h.h"
#include <math.h>
#include <stdio.h>


void parcour_all(t_struct *mystruct);

void vector_normalize(t_vec3d *v)
{
  double id;
  id =  1/sqrt((v->x * v->x) + (v->y * v->y) + (v->z * v->z));
  v->x = v->x * id;
  v->y = v->y *id;
  v->z = v->z *id;
}

double vector_dot(t_vec3d *a, t_vec3d *b)
{
  double res;

  res = (a->x * b->x) + (a->y * b->y) + (a->z * b->z);

  //  return(a->x * b->x + a->y * b->y + a->z * b->z);
  return (res);
}

t_vec3d *vector_copy(t_vec3d *a)//bad type?
{
  t_vec3d *v;
  v = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  v->x = a->x;
  v->y = a->y;
  v->z = a->z;
  return (v);
}

t_vec3d *vector_sub(t_vec3d *a, t_vec3d *b)
{
  t_vec3d *v;
  v = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  v->x = a->x - b->x;
  v->y = a->y - b->y;
  v->z = a->z - b->z;
  return (v);
}

t_vec3d *vector_mul(t_vec3d *a, t_vec3d *b)
{
  t_vec3d *v;
  v = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  v->x = a->x * b->x;
  v->y = a->y * b->y;
  v->z = a->z * b->z;
  return (v);
}

int             my_pixel_put_to_image2(t_img *img, int x, int y, int color)
{
  int           a;
  int           b;
  unsigned char *col;

  if (x >= WIDTH || y >= HEIGHT || x < 0 || y < 0)
    return (0);
  img->data = mlx_get_data_addr(img->img_ptr,
                                  &img->bpp, &img->sizeline, &img->endian);

  col = (unsigned char *)&color;
  a = img->bpp / 8;
  b = img->sizeline;
  if (img->cos < 0)
    img->cos = 0;
  img->data[(x * a) + 0 + (y * b)] = col[0] * img->cos;
  img->data[(x * a) + 1 + (y * b)] = col[1] * img->cos;
  img->data[(x * a) + 2 + (y * b)] = col[2] * img->cos;
  img->data[(x * a) + 3 + (y * b)] = col[3] * img->cos;
  return (0);
}



void            my_pixel_put_to_image(t_img *myimg, int x, int y, int color)
{

  if(x >= WIDTH || x <= 0 || y >= HEIGHT || y <= 0)
    {
      return;
    } 
  myimg->data = mlx_get_data_addr(myimg->img_ptr,
				  &myimg->bpp, &myimg->sizeline, &myimg->endian);
  myimg->data[y * myimg->sizeline + x * myimg->bpp / 8] = color % 256;
  color /= 256;
  myimg->data[y * myimg->sizeline + x * myimg->bpp / 8 + 1] = color % 256;
  color /= 256;
  myimg->data[y * myimg->sizeline + x * myimg->bpp / 8 + 2] = color % 256;
}

int                     event_mlx(int keycode, t_struct *mystruct)
{
  if (keycode == 65307)
    exit(1);
  if(keycode == 65362)
    {
      //haut
      mystruct->decz += 50;
    }
  if(keycode == 65364)
    {
      //bas
      mystruct->decz -= 50;
    }
  if(keycode == 65361)
    {
      //gauche
      mystruct->decy -= 50;
    }
  if(keycode == 65362)
    {
      //droite
      mystruct->decy += 50;
    }

  mlx_destroy_image(mystruct->mlx, mystruct->img->img_ptr);
  mystruct->img->img_ptr = mlx_new_image(mystruct->mlx, WIDTH, HEIGHT);
  parcour_all(mystruct);
  mlx_put_image_to_window(mystruct->mlx, mystruct->win, mystruct->img->img_ptr, 0, 0);
  return (0);
}

double mypow(int a, int n)
{
  double res;
  int i;

  i = 1;
  res = a;
  while(i < n)
    {
      res = res * a;
      i++;
    }
  return (res);
}

int plan_int(t_ray *ray, t_vec3d *plan, t_light *light)
{
  if(ray->d->z == 0)
    return (0);
  plan->v = (-ray->o->z) / ray->d->z;
  light->ox = 0;
  light->oy = 0;
  light->oz = 100;
  if(plan->v > 0)
    return (2);
  else
    return (1);
}

void sphere_int(t_sphere *sphere, t_ray *ray)
{
  double a;
  double b;
  double c;
  double det;
  double myv;

  a = (sphere->pos->x * sphere->pos->x + sphere->pos->y * sphere->pos->y + sphere->pos->z *sphere->pos->z);
  b = (2 * sphere->pos->x * ray->o->x + 2 * sphere->pos->y * ray->o->y + 2 * sphere->pos->z * ray->o->z);
  c = (ray->o->x * ray->o->x + ray->o->y * ray->o->y + ray->o->z * ray->o->z - sphere->radius * sphere->radius);
  det = (b * b) - (4 * a * c);
  sphere->pos->v = (-b - sqrt(det)) / (2 * a);
  myv = (-b + sqrt(det)) / (2 * a);
  if(det < 0)
    sphere->pos->v = -1;
  else if (sphere->pos->v > myv && myv > 0)
    sphere->pos->v = myv;
}


int light_sphere(t_sphere *sphere, t_light *light, t_ray *ray, t_struct *mystruct)
{
  double v1;
  double n1;
  double n2;
  double tempx;
  double tempy;
  double tempz;

  tempx = ray->o->x + (sphere->pos->x * sphere->pos->v);
  tempy = ray->o->y + (sphere->pos->y * sphere->pos->v);
  tempz = ray->o->z + (sphere->pos->z * sphere->pos->v);
  light->ox = tempx;
  light->oy = tempy;
  light->oz = tempz;
  light->dx = light->x - tempx;
  light->dy = light->y - tempy;
  light->dz = light->z - tempz;
  v1 = (light->ox * light->dx + light->oy * light->dy + light->oz * light->dz);
  n1 = (mypow(light->ox,2) + mypow(light->oy,2) + mypow(light->oz, 2));
  n2 = (mypow(light->dx,2) + mypow(light->dy,2) + mypow(light->dz, 2));
  mystruct->img->cos = v1 / (sqrt(n1) * sqrt(n2));
  return (0xFF0000);
}

int light_plan(t_vec3d *plan, t_ray *ray, t_light *light, t_struct *mystruct)
{
  double v1;
  double n1;
  double n2;
  double tempx;
  double tempy;
  double tempz;

  tempx = ray->o->x + ray->d->x * plan->v;
  tempy = ray->o->y + ray->d->y * plan->v;
  tempz = ray->o->z + ray->d->z * plan->v;
  light->dx = light->x - tempx;
  light->dy = light->y - tempy;
  light->dz = light->z - tempz;
  v1 = (light->ox * light->dx + light->oy * light->dy + light->oz * light->dz);
  n1 = (mypow(light->ox,2) + mypow(light->oy,2) + mypow(light->oz, 2));
  n2 = (mypow(light->dx,2) + mypow(light->dy,2) + mypow(light->dz, 2));
  mystruct->img->cos = v1 / (sqrt(n1) * sqrt(n2));
  return (0x00FF00);
}

void give_value(int x, int y, t_sphere *sphere, t_ray *ray, t_light *light, t_struct *mystruct)
{
  light->x = -3000;
  light->y = 0;
  light->z = 2000;
  ray->o->x = -3000;
  ray->o->y = 0; 
  ray->o->z = 500;   
  ray->d->x = 0; 
  ray->d->y = (WIDTH / 2) - x; 
  ray->d->z = (HEIGHT / 2) - y;
  sphere->pos->x = ray->d->x + mystruct->decx;
  sphere->pos->y = ray->d->y + mystruct->decy;
  sphere->pos->z = ray->d->z + mystruct->decz;
  sphere->radius = 800;
}


int lets_choice(t_sphere *sphere, t_vec3d *plan, t_light *light, t_ray *ray, t_struct *mystruct)
{
  double v;
  double mytab[2];
  int fptr[2];
  int           i;

  if(plan->v < sphere->pos->v)
    {
      v = plan->v;
      i = 0;
      if(v <= 0)
	{
	  v = sphere->pos->v;
	  i = 1;
	}
    }
  else
    {
      v = sphere->pos->v;
      i = 1;
      if(v <= 0)
	{
	  v = plan->v;
	  i = 0;
	}
    }
  if(v <= 0)
    return (0);
  if(i == 0)
    return (light_plan(plan, ray, light, mystruct));
  else if(i == 1)
    return (light_sphere(sphere, light, ray, mystruct));
  return (0);
}

int define_color(int x, int y, t_struct *mystruct)
{
  int color;
  t_sphere *sphere;
  t_vec3d *plan;
  t_ray *ray;
  t_light *light;

  light = (t_light *)malloc(sizeof(t_light) * 1);
  ray = (t_ray *)malloc(sizeof(t_ray) * 1);
  ray->o = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  ray->d = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  plan = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  sphere = (t_sphere *)malloc(sizeof(t_sphere) * 1);
  sphere->pos = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  color = 0;
  give_value(x, y,sphere, ray, light, mystruct);
  plan_int(ray, plan, light);
  sphere_int(sphere, ray);
  color = lets_choice(sphere, plan,light, ray, mystruct);
  return (color);
}

void parcour_all(t_struct *mystruct)
{
  int i;
  int j;
  int mypixelcolor;

  mypixelcolor = 0;
  j = 0;
  while(j < HEIGHT)
    {
      i = 0;
      while(i < WIDTH)
	{
	  mypixelcolor = define_color(i, j, mystruct);
	  if(mypixelcolor != 0)
	    {
	      my_pixel_put_to_image2(mystruct->img, i, j, mypixelcolor);
	    }
	  i++;
	}
      j++;
    }
}


int                     main(int argc, char **argv)
{
  t_struct *mystruct;
  t_img *myimg;

  //create win and img
  mystruct = (t_struct *)malloc(sizeof(t_struct));
  myimg = (t_img *)malloc(sizeof(t_img));
  mystruct->mlx = mlx_init();
  mystruct->win = mlx_new_window(mystruct->mlx, WIDTH, HEIGHT, "RayTracer V1");
  myimg->img_ptr = mlx_new_image(mystruct->mlx, WIDTH, HEIGHT);
  mystruct->img = myimg;
  //////////////
  mystruct->decx = 1000;
  mystruct->decy = 0;
  mystruct->decz = -300;
  parcour_all(mystruct);
  mlx_put_image_to_window(mystruct->mlx, mystruct->win, mystruct->img->img_ptr, 0, 0);
  mlx_key_hook(mystruct->win, event_mlx, mystruct);
  mlx_loop(mystruct->mlx);
  return (0);
}
