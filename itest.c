#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include <glib.h>

#include "hefsi.h"

gchar *progname ;

gint surface_sphere(gdouble u, gdouble v, gdouble *x, gpointer data)

{
  gdouble *x0 = data, r ;

  r = x0[3] ;
  
  x[0] = r*sin(M_PI*u)*cos(M_PI*v) + x0[0] ;
  x[1] = r*sin(M_PI*u)*sin(M_PI*v) + x0[1] ;
  x[2] = r*cos(M_PI*u)             + x0[2] ;
  
  return 0 ;
}

gint surface_sphere_d(gdouble u, gdouble v, gdouble *x,
		      gdouble *xu, gdouble *xv, gpointer data)

{
  gdouble *x0 = data, r ;

  r = x0[3] ;
  
  x[0] = r*sin(M_PI*u)*cos(M_PI*v) + x0[0] ;
  x[1] = r*sin(M_PI*u)*sin(M_PI*v) + x0[1] ;
  x[2] = r*cos(M_PI*u)             + x0[2] ;

  xu[0] =  M_PI*r*cos(M_PI*u)*cos(M_PI*v) ;
  xu[1] =  M_PI*r*cos(M_PI*u)*sin(M_PI*v) ;
  xu[2] = -M_PI*r*sin(M_PI*u) ;

  xv[0] = -M_PI*r*sin(M_PI*u)*sin(M_PI*v) ;
  xv[1] =  M_PI*r*sin(M_PI*u)*cos(M_PI*v) ;
  xv[2] =  0 ;
  
  return 0 ;
}

gint surface_cylinder(gdouble u, gdouble v, gdouble *x, gpointer data)

{
  gdouble *x0 = data, r, len ;

  r = x0[3] ;
  len = x0[4] ;
  
  x[0] = r*cos(M_PI*v) + x0[0] ;
  x[1] = r*sin(M_PI*v) + x0[1] ;
  x[2] = len*u         + x0[2] ;
  
  return 0 ;
}

gint surface_cylinder_d(gdouble u, gdouble v, gdouble *x,
			gdouble *xu, gdouble *xv, gpointer data)

{
  gdouble *x0 = data, r, len ;

  r = x0[3] ;
  len = x0[4] ;
  
  x[0] = r*cos(M_PI*v) + x0[0] ;
  x[1] = r*sin(M_PI*v) + x0[1] ;
  x[2] = len*u         + x0[2] ;

  xu[0] = 0.0 ;
  xu[1] = 0.0 ;
  xu[2] = len ;
  
  xv[0] = -M_PI*r*sin(M_PI*v) ;
  xv[1] =  M_PI*r*cos(M_PI*v) ;
  xv[2] =  0.0 ;
  
  return 0 ;
}

static void print_help_message(FILE *f, gdouble scale)

{
  fprintf(f, "%s: usage\n\n", progname) ;
  fprintf(f,
	  "  %s [opts] > curve.dat\n\n", progname) ;
  fprintf(f,
	  "Options:\n\n"
	  "  -h print this message and exit\n"
	  "  -d use derivative-free surface evaluators\n"
	  "  -s # scaling factor for element "
	  "bounding box sizes (%lg)\n", scale) ;
  
  return ;
}

gint main(gint argc, gchar **argv)

{
  hefsi_surface_t *s1, *s2 ;
  gpointer func1, func2 ;
  hefsi_segment_t *seg ;
  hefsi_workspace_t *w ;
  GSList *il ;
  gdouble umin, umax, vmin, vmax, x1[5], x2[5], tol, scale, emax ;
  gdouble x[3], uv1[2], uv2[2] ;
  gint dmin, dmax, i, j ;
  FILE *output ;
  gboolean derivative ;
  gchar ch ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;
  scale = 1.0625 ;
  umin = 0 ; umax = 1 ;
  vmin = -1 ; vmax = 1 ;
  x1[0] =  0.3 ; x1[1] = 0.1 ; x1[2] = 1.5 ; x1[3] = 0.9 ; x1[4] = 2.5 ;
  x2[0] = -0.1 ; x2[1] = 0.1 ; x2[2] = 0.1 ; x2[3] = 0.401 ; x2[4] = 2.5 ;

  derivative = TRUE ;

  while ( (ch = getopt(argc, argv, "hds:")) != EOF ) {
    switch (ch) {
    default: return 1 ; break ;
    case 'h': print_help_message(stderr, scale) ; return 0 ; break ;
    case 'd': derivative = FALSE ; break ;
    case 's': scale = atof(optarg) ; break ;
    }
  }
  
  if ( derivative ) {
    func1 = surface_sphere_d ; func2 = surface_cylinder_d ;
  } else {
    func1 = surface_sphere ; func2 = surface_cylinder ;
  }
  s1 = hefsi_surface_new(func1, x1, derivative, umin, umax, vmin, vmax) ;
  s2 = hefsi_surface_new(func2, x2, derivative, umin, umax, vmin, vmax) ;

  dmin = 5 ; dmax = 8 ; tol = 1e-9 ;
  hefsi_surface_initialize(s1, dmin, dmax) ;
  hefsi_set_bounding_boxes(s1, scale) ;
  fprintf(stderr, "%s: surface 1 initialized\n", progname) ;
  hefsi_surface_initialize(s2, dmin, dmax) ;
  hefsi_set_bounding_boxes(s2, scale) ;
  fprintf(stderr, "%s: surface 2 initialized\n", progname) ;
  
  w = hefsi_workspace_new() ;

  hefsi_surface_intersections(s2, s1, tol, w) ;
  /*do it twice to check that the workspace is properly reinitialized
    when it is reused*/
  hefsi_surface_intersections(s1, s2, tol, w) ;
  
  fprintf(stderr, "%d curves; %d segments\n",
	  hefsi_workspace_curve_number(w),
	  hefsi_workspace_segment_number(w)) ;

  i = hefsi_surface_intersection_check(s1, s2, w, tol, &emax) ;

  fprintf(stderr, "%d segments fail intersection check; emax=%lg\n",
	  i, emax) ;
  
  for ( i = 0 ; i < hefsi_workspace_curve_number(w) ; i ++ ) {
    for ( il = hefsi_workspace_curve(w,i) ; il != NULL ; il = il->next ) {
      j = GPOINTER_TO_INT(il->data) ;
      seg = hefsi_workspace_segment(w,j) ;
      x[0] = 0.5*(seg->x1[0]+seg->x2[0]) ;
      x[1] = 0.5*(seg->x1[1]+seg->x2[1]) ;
      x[2] = 0.5*(seg->x1[2]+seg->x2[2]) ;
      uv1[0] = 0.5*(seg->x1[3]+seg->x2[3]) ;
      uv1[1] = 0.5*(seg->x1[4]+seg->x2[4]) ;
      uv2[0] = 0.5*(seg->x1[5]+seg->x2[5]) ;
      uv2[1] = 0.5*(seg->x1[6]+seg->x2[6]) ;
      hefsi_refine_point(s1, uv1, s2, uv2, x, tol, 8) ;
      fprintf(stdout, "%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %d\n",
	      x[0], x[1], x[2], uv1[0], uv1[1], uv2[0], uv2[1], i) ;
    }
  }
  
  output = fopen("surface1.dat", "w") ;
  hefsi_write_surface(output, s1) ;
  fclose(output) ;

  output = fopen("surface2.dat", "w") ;
  hefsi_write_surface(output, s2) ;
  fclose(output) ;
  
  output = fopen("boxes1.dat", "w") ;
  hefsi_write_bounding_boxes(output, s1) ;
  fclose(output) ;

  output = fopen("boxes2.dat", "w") ;
  hefsi_write_bounding_boxes(output, s2) ;
  fclose(output) ;
  
  return 0 ;
}
