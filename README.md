* Calculation of intersections between parametric surfaces

HEFSI calculates surface intersections based on the method of
Houghton, Emnett, Factor, and Sabharwal, "Implementation of a
divide-and-conquer method for intersection of parametric surfaces",
Computer Aided Geometric Design Volume 2, Issues 1-3, September 1985,
173-183, https://doi.org/10.1016/0167-8396(85)90022-6.

HEFSI includes Tomas Moller's triangle-triangle intersection code from
A Fast Triangle-Triangle Intersection Test, Journal of Graphics Tools,
2(2), 1997. The code and paper are available from Professor Moller's
page:

https://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/

* Installation

HEFSI uses the glib-2.0 library for a number of utility functions and
should compile on any standard Linux system, using the supplied
Makefile. If you want to use HEFSI in your own code (this is how it is
intended to be used) you need to add hefsi.c and the header hefsi.h to
your source code.

To run the sample calculation, for the intersection of a sphere and a
cylinder, run

> ./itest > curve.dat

This will output two files of the triangulations of the surfaces,
which can be read in Octave using the readtri function supplied, and
visualized using Octave's built-in trisurf function:

> [x1,t1]=readtri("surface1.dat") ;
> trisurf(t1(:,1:3)+1, x1(:,1), x1(:,2), x1(:,3))
> [x2,t2]=readtri("surface2.dat") ;
> trisurf(t2(:,1:3)+1, x2(:,1), x2(:,2), x2(:,3))

The intersection curve is written as a connected sequence of
intersection points, one per line:

> xi yi zi u1i v1i u2i v2i inter

where (xi,yi,zi) is the ith intersection point, (u1i,v1i) are the
parametric coordinates of the point on surface 1, (u2i,v2i) are the
parametric coordinates of the point on surface 2, and inter is the
number of the intersection curve (there may be more than one curve of
intersection between the surfaces).

There are some command line options for the test code, which are
listed by the help option

>  ./itest -h


