#include "my_h.h"
#include <math.h>
#include <stdio.h>


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

int intersection_sphere(t_sphere *sphere, t_ray *ray, double *coef)
{
  t_vec3d *dist;
  dist = vector_copy(sphere->pos);
  dist = vector_sub(dist, ray->o);
  double B = vector_dot(ray->d, dist);
  double D = (B * B) - vector_dot(dist, dist) + (sphere->radius * sphere->radius);
  if(D < 0)
    {
      return (-1);
    }
  float T0 = B - sqrt(D);//bad type?
  float T1 = B + sqrt(D);
  int result = 0;

  if(T0 > 0.1 && T0 < *coef)
    {
      *coef = T0;
      return (0);
    }
  if(T1 > 0.1 && T1 < *coef)
    {
      *coef = T1;
      return (0);
    }
  return (-1);
}

double inter_sphere(t_ray *ray, t_sphere *sphere)
{
  double a;
  double b;
  double c;
  a = (ray->d->x * ray->d->x) + (ray->d->y * ray->d->y) + (ray->d->z * ray->d->z);
  b = 2 * (ray->d->x * (ray->o->x - sphere->pos->x) + ray->d->y * (ray->o->y - sphere->pos->y) + ray->d->z * (ray->o->z - sphere->pos->z));
  c = (((ray->o->x - sphere->pos->x) * (ray->o->x - sphere->pos->x)) + ((ray->o->y - sphere->pos->y) * (ray->o->y - sphere->pos->y)) + ((ray->o->z - sphere->pos->z) * (ray->o->z - sphere->pos->z))) - (sphere->radius * sphere->radius);
  double det;
  det = (b * b) - (4 * a * c);
  if(det < 0)
    {
      return (-1);
    }
  double t1;
  double t2;
  t1 = (-b + sqrt(det)) / (2*a);
  t2 =  (-b - sqrt(det)) / (2*a);
  if(det == 0)
    return (t1);
  if(t1 < t2)
    {
      return (t1);
    }
  else
    return (t2);
}

double inter_p(t_ray *ray)
{
  double a = -10;
  double b = -10;
  double c = -10;

  if(a * ray->o->x + b * ray->o->y + c * ray->o->z == 0)
    {
      return (1);
    }
  return (-1);
}


double inter_cylindre(t_ray *ray, t_sphere *sphere)
{
  double a;
  double b;
  double c;
  a = (ray->o->x * ray->o->x) + (ray->o->y * ray->o->y);
  b = (2 * (ray->o->x * sphere->pos->x)) + (2 * (ray->o->y * sphere->pos->y));
  c = (sphere->pos->x * sphere->pos->x) + (sphere->pos->y * sphere->pos->y) - (sphere->radius * sphere->radius);
  double det;
  det = (b * b) - (4 * a * c);
  if(det < 0)
    {
      return (-1);
    }
  double t1;
  double t2;
  t1 = (-b + sqrt(det)) / (2*a);
  t2 =  (-b - sqrt(det)) / (2*a);
  if(det == 0)
    return (t1);
  if(t1 < t2)
    {
      return (t1);
    }
  else
  return (t2);
  return (2);
}

int ft_plane(t_ray *ray, t_vec3d *plane)
{
  return (0);
}
//x0 y0 z0 a point on the plane
/* plane
nx(x-x0)+ny(y-y0)+nz(z - z0) = 0
nx ny nz normal vector of the plane
*/

//Ray
//x = xs + t * xd
/*y = ys + t * yd
z = zs + t *zd
xs ys zs sart point
xd yd zd direction vector
*/
//nx(xs + t * xd-x0) + ny(ys + t *yd - y0)+nz(zs + t * zd - z0) = 0
//nx * xs + t *xd *nx - x0 * nx + ny * ys + t * yd * ny + y0+ny + nz * zs + t * zd * nz - nz - n * z0 = 0



void do_all(t_struct *mystruct, double z1, double z2, double x1, double y1)
{
  t_sphere *sphere;
  sphere = (t_sphere *)malloc(sizeof(t_sphere) * 3);
  sphere->pos = (t_vec3d *)malloc(sizeof(t_vec3d) * 4);
  sphere->pos->x = 0;
  sphere->pos->y = 0;
  sphere->pos->z = z1;
  sphere->radius = 20;

  t_sphere *sphere2;
  sphere2 = (t_sphere *)malloc(sizeof(t_sphere) * 3);
  sphere2->pos = (t_vec3d *)malloc(sizeof(t_vec3d) * 4);
  sphere2->pos->x = 1;
  sphere2->pos->y = 0;
  sphere2->pos->z = 200;
  sphere2->radius = 20;

  t_sphere *sphere3;
  sphere3 = (t_sphere *)malloc(sizeof(t_sphere) * 3);
  sphere3->pos = (t_vec3d *)malloc(sizeof(t_vec3d) * 4);
  sphere3->pos->x = 80;
  sphere3->pos->y = 0;
  sphere3->pos->z = z2;
  sphere3->radius = 20;

  t_sphere *sphere4;
  sphere4 = (t_sphere *)malloc(sizeof(t_sphere) * 3);
  sphere4->pos = (t_vec3d *)malloc(sizeof(t_vec3d) * 4);
  sphere4->pos->x = 20;
  sphere4->pos->y = 10;
  sphere4->pos->z = 800;
  sphere4->radius = 1;


  t_vec3d *plane;
  plane = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  plane->x = 0;
  plane->y = 0;
  plane->z = 1;

  t_vec3d *a;
  a = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  a->x = x1;
  a->y = y1;
  a->z = -1;
  
  t_vec3d *b;
  b = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  b->x = 0;
  b->y = 0;
  b->z = 0;

  t_vec3d *raydir;
  raydir = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  
  t_ray *rayon;
  rayon = (t_ray *)malloc(sizeof(t_ray) * 1);
  rayon->o = (t_vec3d *)malloc(sizeof(t_vec3d) *1);  
  rayon->d = (t_vec3d *)malloc(sizeof(t_vec3d) *1);
  rayon->o->x = a->x;
  rayon->o->y = a->y;
  rayon->o->z = a->z;
	  
  int debug;
  debug = 0;

  t_vec3d *Lumiere;
  Lumiere = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  Lumiere->x = 0;
  Lumiere->y = 0;//-500;
    Lumiere->z = -1000;
  
  t_vec3d *LumiereCouleur;
  LumiereCouleur = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
  LumiereCouleur->x = 0.5;
  LumiereCouleur->y = 0.5;
  LumiereCouleur->z = 0.5;


for(int pixel_y = 0; pixel_y < HEIGHT - 1; pixel_y++)
    {
      for(int pixel_x = 0; pixel_x < WIDTH - 1; pixel_x++)
	{
	  b->x = pixel_x - (WIDTH / 2);
	  b->y = pixel_y - (HEIGHT / 2);
	  b->z = -(WIDTH / (2 * tan(30 / 2)));
	  vector_normalize(b);

	  raydir->x = b->x - a->x;
	  raydir->y = b->y - a->y;
	  raydir->z = b->z - a->z;
	  vector_normalize(raydir);

	  
	  rayon->d->x = raydir->x;
	  rayon->d->y = raydir->y;
	  rayon->d->z = raydir->z;
      	  
	  double coef = 20000;
	  //int inter = intersection_sphere(sphere,rayon, &coef);
	  int inter2 = intersection_sphere(sphere,rayon, &coef);
	  if(inter2 != -1)
	  {
	      debug++;
	      t_vec3d *Impact;
	      Impact = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
	      Impact->x = rayon->o->x + rayon->d->x * coef;
	      Impact->y = rayon->o->y + rayon->d->x * coef;
	      Impact->z = rayon->o->z + rayon->d->z * coef;
    
	       t_vec3d *Direction;
	      t_vec3d *Distance;
	      Direction = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
	      Distance = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
	      Distance->x = Impact->x - sphere->pos->x;
	      Distance->y = Impact->y - sphere->pos->y;
	      Distance->z = Impact->z - sphere->pos->z;
	      vector_normalize(Distance);
	      double Norme;
	      Norme = vector_dot(Distance,Distance);
	      float Temp;
		Temp = Norme * Norme;
		Temp = 1 / sqrt(Temp);
		Norme = Temp * Norme;
		if(Norme != 0)
		  {
		    t_vec3d *RayLightDist;
		    RayLightDist = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
		    RayLightDist->x = Lumiere->x - Impact->x;
		    RayLightDist->y = Lumiere->y - Impact->y;
		    RayLightDist->z = Lumiere->z - Impact->z;
		    double NormeRayDist = vector_dot(RayLightDist,RayLightDist);
      
		    t_vec3d *RayLightDir;
		    RayLightDir = (t_vec3d *)malloc(sizeof(t_vec3d) * 1);
		    if (NormeRayDist != 0)
		      {
			RayLightDir->x = RayLightDist->x / Norme;
			RayLightDir->y = RayLightDist->y / Norme;
			RayLightDir->z = RayLightDist->z / Norme;
		      }		   
		    t_vec3d *IntensityLight;
		    IntensityLight = vector_mul(RayLightDir,Distance);
		    float lambert = (IntensityLight->x * Norme) + (IntensityLight->y * Norme)+(IntensityLight->z * Norme);
        
		    double Red    = (LumiereCouleur->x * 0.5  * lambert);
		    double Green  = (LumiereCouleur->y * 0.5  * lambert);
		    double Blue   = (LumiereCouleur->z * 0.5  * lambert) ;
		    if(Red < 0)
		      Red = 0;
		    if(Green < 0)
		      Green = 0;
		    if(Blue < 0)
		      Blue = 0;
		    if(Red > 255)
		      Red = 255;
		    if(Blue > 255)
		      Blue = 255;
		    if(Green > 255)
		      Green = 255;
		    my_pixel_put_to_image(mystruct->img, pixel_x, pixel_y, RGB(Red,Green,Blue));		    

		  }
	    }
	}
    }
/* if(debug == 0)
   {
     printf("debug work\n");
     do_all(mystruct,z1,z2);
     }*/
}

void line(t_img *myimg, float xi, float yi, float xf, float yf, int color)
{
  int dx,dy,i,xinc,yinc,cumul,x,y ;
  x = xi ;
  y = yi ;
  dx = xf - xi ;
  dy = yf - yi ;
  xinc = ( dx > 0 ) ? 1 : -1 ;
  yinc = ( dy > 0 ) ? 1 : -1 ;
  dx = abs(dx) ;
  dy = abs(dy) ;
  //  if(ok_draw(yi,yf))
  my_pixel_put_to_image(myimg, x, y, color);
  if ( dx > dy ) {
    cumul = dx / 2 ;
    for ( i = 1 ; i <= dx ; i++ ) {
      x += xinc ;
      cumul += dy ;
      if ( cumul >= dx ) {
        cumul -= dx ;
        y += yinc ; }
      // if(ok_draw(yi,yf))
      my_pixel_put_to_image(myimg, x, y, color); } }
  else {
    cumul = dy / 2 ;
    for ( i = 1 ; i <= dy ; i++ ) {
      y += yinc ;
      cumul += dx ;
      if ( cumul >= dy ) {
        cumul -= dy ;
        x += xinc ; }
      // if(ok_draw(yi,yf))
      my_pixel_put_to_image(myimg, x, y, color); } }
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
  static double z1;
  static double z2;
  static double x1;
  static double y1;
  static double z;

  if(!z)
    z = -1;
  if(!x1)
    x1 = 0;
  if(!y1)
    y1 = 0;
  if(!z1)
    z1 = 400;
  if(!z2)
    z2 = 400;

  if (keycode == 53)
    exit(1);
  if(keycode == 123)
    {
      z1 -= 20;
      mlx_destroy_image(mystruct->mlx, mystruct->img->img_ptr);
      mystruct->img->img_ptr = mlx_new_image(mystruct->mlx, WIDTH, HEIGHT);
      do_all(mystruct, z1, z1, x1, y1);
      mlx_put_image_to_window(mystruct->mlx, mystruct->win, mystruct->img->img_ptr, 0, 0);
    }
  if(keycode == 124)
    {
      z1 += 20;
      mlx_destroy_image(mystruct->mlx, mystruct->img->img_ptr);
      mystruct->img->img_ptr = mlx_new_image(mystruct->mlx, WIDTH, HEIGHT);
      do_all(mystruct, z1, z1, x1, y1);
      mlx_put_image_to_window(mystruct->mlx, mystruct->win, mystruct->img->img_ptr, 0, 0);
    }
  if(keycode == 0)
    {
      //gauche
      x1 += 0.03;
      mlx_destroy_image(mystruct->mlx, mystruct->img->img_ptr);
      mystruct->img->img_ptr = mlx_new_image(mystruct->mlx, WIDTH, HEIGHT);
      do_all(mystruct, z1, z1, x1, y1);
      mlx_put_image_to_window(mystruct->mlx, mystruct->win, mystruct->img->img_ptr, 0, 0);
    }
  if(keycode == 2)
    {
      //droite
      x1 -= 0.03;
      mlx_destroy_image(mystruct->mlx, mystruct->img->img_ptr);
      mystruct->img->img_ptr = mlx_new_image(mystruct->mlx, WIDTH, HEIGHT);
      do_all(mystruct, z1, z1, x1, y1);
      mlx_put_image_to_window(mystruct->mlx, mystruct->win, mystruct->img->img_ptr, 0, 0);
    }
  if(keycode == 13)
    {
      //haut
      y1 += 0.03;
      mlx_destroy_image(mystruct->mlx, mystruct->img->img_ptr);
      mystruct->img->img_ptr = mlx_new_image(mystruct->mlx, WIDTH, HEIGHT);
      do_all(mystruct, z1, z1, x1, y1);
      mlx_put_image_to_window(mystruct->mlx, mystruct->win, mystruct->img->img_ptr, 0, 0);
    }
  if(keycode == 1)
    {
      y1 -= 0.03;
      //bas
      mlx_destroy_image(mystruct->mlx, mystruct->img->img_ptr);
      mystruct->img->img_ptr = mlx_new_image(mystruct->mlx, WIDTH, HEIGHT);
      do_all(mystruct, z1, z1, x1, y1);
      mlx_put_image_to_window(mystruct->mlx, mystruct->win, mystruct->img->img_ptr, 0, 0);
    }
  return (0);
}



void DrawFilledCircle(t_struct *mystruct, int x0, int y0, int z0, int radius, int color)
{
  x0 += z0;
  y0 -= z0;
  int x = radius;
  int y = 0;
  int xChange = 1 - (radius << 1);
  int yChange = 0;
  int radiusError = 0;

  while (x >= y)
    {
      for (int i = x0 - x; i <= x0 + x; i++)
        {
	  my_pixel_put_to_image(mystruct->img, i, y0 + y, color / 3);
	  my_pixel_put_to_image(mystruct->img, i, y0 - y, color / 3);
	}
      for (int i = x0 - y; i <= x0 + y; i++)
        {
	    my_pixel_put_to_image(mystruct->img, i, y0 + x, color / 3);
	  	  my_pixel_put_to_image(mystruct->img, i, y0 - x, color / 3);
	}
      y++;
      radiusError += yChange;
      yChange += 2;
      if (((radiusError << 1) + xChange) > 0)
        {
	  x--;
	  radiusError += xChange;
	  xChange += 2;
	  }
	}
}

void draw_floor(t_struct *mystruct)
{
  int i;

  i = 300;

  while(i < HEIGHT)
    {
      line(mystruct->img, 0, i, WIDTH, i, 0x00FFFF);
      i++;
    }
}

void draw_center(t_struct *mystruct)
{
  line(mystruct->img, WIDTH / 2,  0, WIDTH / 2, HEIGHT, 0xFF00);
  line(mystruct->img, 0,  HEIGHT / 2, WIDTH, HEIGHT / 2, 0xFF00);
  line(mystruct->img, WIDTH / 2,  HEIGHT / 2, WIDTH, 0, 0xFF00);
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

void draw_spheres(t_struct *mystruct, int color, int xc, int yc,int zc, int r)
{
  int x;
  int y;
  int z;
  int xp;
  int yp;
  int i;

  i = 0;
  xp = 0;
  yp = 0;
  z = 0;
  y = 0;
  while(y < HEIGHT)
    {
      x = 0;
      while(x < WIDTH)
	{
	  while(z < 50)
	    {
	      if(mypow((x-xc),2) + mypow((y-yc),2) + mypow((z-zc),2) <= mypow(r,2))
		{
		  i++;
		  my_pixel_put_to_image(mystruct->img, x, y, color);
		  if(i >= 2)
		    {
		      //line(mystruct->img, xp, yp, x, y, color / 3);
		    }
		    xp = x;
		  yp = y;
		}
	      z++;
	    }
	  z = 0;
	  x++;
	}
    y++;
    }
}

void draw_cylindre(t_struct *mystruct, int color, int xc, int yc,int zc, int r)
{
  int x;
  int y;
  int z;

  z = 0;
  y = 0;
  while(y < HEIGHT)
    {
      x = 0;
      while(x < WIDTH)
        {
          if(mypow((x-xc),2) + mypow((y-yc),2) <= mypow(r,2))
            {
              my_pixel_put_to_image(mystruct->img, x, y, color);
            }
          x++;
        }
      y++;
    }
}

void draw_cone(t_struct *mystruct, int color, int xc, int yc,int h, int r)
{
  int x;
  int y;
  int z;

  z = 0;
  y = 0;
  while(y < HEIGHT)
    {
      x = 0;
      while(x < WIDTH)
        {
          if(mypow((x-xc),2) + mypow((y-yc),2) <= mypow(h-z,2) * mypow(r,2)/mypow(h,2))
            {
              my_pixel_put_to_image(mystruct->img, x, y, color);
            }
          x++;
        }
      y++;
    }
}

void draw_plan(t_struct *mystruct, int color, int xc, int yc,int zc, int d)
{
  int x;
  int y;
  int z;

  z = 0;
  y = 0;
  while(y < HEIGHT)
    {
      x = 0;
      while(x < WIDTH)
        {
	  if((x - xc) + (y - yc) + (z - zc) + d >= 0)
            {
              my_pixel_put_to_image(mystruct->img, x, y, color);
            }
          x++;
        }
      y++;
    }
}

int                     main(int argc, char **argv)
{
  t_struct *mystruct;
  t_img *myimg;

  mystruct = (t_struct *)malloc(sizeof(t_struct));
  myimg = (t_img *)malloc(sizeof(t_img));
  mystruct->mlx = mlx_init();
  mystruct->win = mlx_new_window(mystruct->mlx, WIDTH, HEIGHT, "RayTracer V1");
  myimg->img_ptr = mlx_new_image(mystruct->mlx, WIDTH, HEIGHT);
  mystruct->img = myimg;

  //  draw_floor(mystruct);
  //DrawFilledCircle(mystruct, WIDTH / 2, 410, -100, 50, 0xFF0000);
  //draw_center(mystruct);
  //   mystruct->color = 0x00FF00;
   //draw_plan(mystruct, 0x0000FF / 3, 0, 400, 0, 100);
  //draw_spheres(mystruct,0x0000FF, 300, 200, 80, 100);
  //draw_spheres(mystruct, 0xFF0000, 600, 200, 40, 100);
  //draw_spheres(mystruct, 0xFFFFFF, 100, 400, 40, 50);
  //  draw_cylindre(mystruct, 0x00FF00, 400, 300, 0, 30);
  //draw_cone(mystruct, 0xFFFF00, WIDTH / 2 - 100, HEIGHT / 2 - 20, 50, 50);
//  DrawFilledCircle(mystruct, (WIDTH / 2) + 20, (HEIGHT / 2) - 20, 100, 0xFF0000 / 2);  
  // my_pixel_put_to_image(mystruct->img, 200, 200, 0x00FFFF);
  do_all(mystruct, 400, 400, 0, 0);  
  mlx_put_image_to_window(mystruct->mlx, mystruct->win, mystruct->img->img_ptr, 0, 0);
  mlx_key_hook(mystruct->win, event_mlx, mystruct);
  mlx_loop(mystruct->mlx);
  return (0);
}
