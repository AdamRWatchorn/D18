/*
  CSC D18 - Assignment 1 - 2D light propagation

  This is the place where you will be doing most of your work for solving this
  assignment. Before you start working here, you shold have become familiar
  with the data structures and functions available to you from light2D.h, and
  with the way a scene is built in buildScene.c

  Go over each of the functions below, and implement the different components
  of the solution in the sections marked

  /************************
  / TO DO:
  ************************ /

  Do not add or modify code outside these sections.

  Details about what needs to be implemented are described in the comments, as
  well as in the assignment handout. You must read both carefully. 

  Starter by: F.J. Estrada, Aug. 2017
*/

/****************************************************************************
 * Uncomment the #define below to enable debug code, add whatever you need
 * to help you debug your program between #ifdef - #endif blocks
 * ************************************************************************/
#define __DEBUG_MODE

/*****************************************************************************
* COMPLETE THIS TEXT BOX:
*
* 1) Student Name:		
* 2) Student Name:		
*
* 1) Student number:
* 2) Student number:
* 
* 1) UtorID
* 2) UtorID
* 
* We hereby certify that the work contained here is our own
*
* ____________________             _____________________
* (sign with your name)            (sign with your name)
********************************************************************************/

struct ray2D makeLightSourceRay(void)
{
 /*
   This function should return a light ray that has its origin at the light
   source, and whose direction depends on the type of light source.

   For point light sources (which emit light in all directions) the direction
    has to be chosen randomly and uniformly over a unit circle (i.e. any
    direction is equally likely)

   For a laser type light source, the direction of the ray is exactly the same
    as the direction of the lightsource.
    
   Set the colour of the ray to the same colour as the light source, and
    set the inside_outside flag to 0 (we assume light sources are 
    outside objects)

   In either case, the direction vector *must be unit length*.
*/
 
 /************************************************************************
 *  TO DO: Complete this function so that we can sample rays from the
 *         lightsource for the light propagation process.
 ************************************************************************/
 
 struct ray2D ray;

 // This creates a dummy ray (which won't go anywhere since the direction
 // vector d=[0 0]!. But it's here so you see which data values you have to
 // provide values for given the light source position, and type.  
 // ** REPLACE THE CODE BELOW with code that provides a valid ray that
 //    is consistent with the lightsource.
 
// ray.p.px=0;			// Ray's origin
// ray.p.py=0;
 ray.p = lightsource.l.p;
// ray.d.px=1.0;			// Ray's direction
// ray.d.py=-0.8;
 ray.d = lightsource.l.d;
 normalize(&ray.d);
 ray.inside_out=0;		// Initially 0 since the ray starts outside an object
 ray.monochromatic=0;		// Initially 0 since the ray is white (from lightsource)
// ray.R=1.0;			// Ray colour in RGB must be the same as the lightsource
 ray.R = lightsource.R;
// ray.G=1.0;
 ray.G = lightsource.G;
// ray.B=1.0;
 ray.B = lightsource.B;
 
 return(ray);
}

void angle(double *theta, struct point2D *vector) {

// probably won't work. may need to rethink + signs

   if(vector->px > 0 && vector->py > 0) {
      *theta = atan(vector->py/vector->px);
//      fprintf(stderr,"vector=(%f, %f)\n",vector->px, vector->py);
   }
   else if(vector->px < 0 && vector->py > 0) {
      *theta = PI + atan(vector->py/vector->px);
   }
   else if(vector->px < 0 && vector->py < 0) {
      *theta = PI + atan(vector->py/vector->px);
   }
   else if(vector->px > 0 && vector->py < 0) {
      *theta = (2*PI) + atan(vector->py/vector->px);
   }
   else if(vector->px > 0 && vector->py == 0) {
      *theta = atan(vector->py/vector->px);
   }
   else if(vector->px == 0 && vector->py > 0) {
      *theta = (PI/2);
   }
   else if(vector->px < 0 && vector->py == 0) {
      *theta = PI;
   }
   else if(vector->px == 0 && vector->py < 0) {
      *theta = ((3*PI)/2);
   }

}

void rotation(double *theta, struct point2D *vector) {

   double px = vector->px, py = vector->py;

   vector->px = px*cos(*theta) - py*sin(*theta);

   vector->py = px*sin(*theta) + py*cos(*theta);

}

void invertDirection(struct point2D *vector) {

   vector->px = -vector->px;
   vector->py = -vector->py;
}

void propagateRay(struct ray2D *ray, int depth)
{
 /*
   This function carries out the light propagation process. It is provided with access
   to a ray data structure, and must perform the following steps (in order!):
   
   - Check if maximum recursion depth has been reached (in which case, it just returns)
   - Find the *closest* intersection between the ray and objects in the scene. This
     means you have to check against the 4 walls, and any circles added in buildScene,
     determine which intersection is closest, and obtain the intersection point, the
     normal at the intersection, and the lambda value at which the intersection happens.
   - Renders the ray onto the image from its starting point all the way up to the 
     intersection point.
   - At the intersection, use the material properties to determine how the propagation
     process proceeds:
         * For mirror materials, compute the mirror reflection direction, create a ray
           along that direction whose origin is the intersection point, and propagate
           that ray
         * For scattering materials, choose a random direction within +- 90 degrees of
           the normal, create a ray with that direction with origin at the intersection
           point, and propagate that ray
         * For refracting materials you will need to propagate two rays - one in the
           mirror reflection direction (just like for reflecting materials), and
           another in the refraction direction. Propagate both rays one after the other!
           
   NOTE: You should only care about intersections for which lambda is POSITIVE (in front
         of the ray), and greater than zero (e.g. if the ray is reflected from some
         object, you do not care about the intersection with the object itself which will
         have a lambda of *very close* to zero)
    
   In every case, make sure the ray's direction vector has unit length. You will need to
   complete other functions as part of your work here.
*/
  
 /*********************************************************************************
  * TO DO: Complete this function to implement light propagation!
  ********************************************************************************/
 
 // Define your local variables here
 struct point2D p1wall, norm, diff, intersection, closestNorm, circInter;
 struct wall2D closestWall;
 double lambda, prevlambda = INFINITY, theta, phi, r_dx;
 int inter_mat_type, circ_mat_type;

 if (depth>=max_depth) return;	 	// Leave this be, it makes sure you don't
					// recurse forever
 

 // Step 1 - Find *closest* intersection with the 4 walls (the written part of A1
 //          should help you figure out how to do that.

 // How many walls can the ray intersect? how many walls can the ray intersect in the
 // forward direction?

 // normal of lines: n = [dy, -dx], then normalize

 //lambda = [(p1wall - p0ray) dot n]/ [dray dot n]

 for (int i = 0; i < 4; i++) {

   // origin + d
   p1wall.px = walls[i].w.p.px + walls[i].w.d.px;
   p1wall.py = walls[i].w.p.py + walls[i].w.d.py;

//   fprintf(stderr,"wallpoint=(%f,%f)\n",p1wall.px,p1wall.py);

   diff.px = p1wall.px - ray->p.px;
   diff.py = p1wall.py - ray->p.py;

//   fprintf(stderr,"diff=(%f,%f)\n",diff.px,diff.py);

   norm.px = -walls[i].w.d.py;
   norm.py = walls[i].w.d.px;

   normalize(&norm);

   /// points norm inwards
   //invertDirection(&norm);

//   fprintf(stderr,"norm=(%f,%f)\n",norm.px,norm.py);

   lambda = dot(&diff, &norm) / dot(&ray->d, &norm);

   fprintf(stderr,"lambda=(%f)\n",lambda);

   if (lambda > TOL) {

      if (lambda < prevlambda) {
         //closestWall = walls[i];
         intersection.px = ray->p.px + lambda*ray->d.px;
         intersection.py = ray->p.py + lambda*ray->d.py;
         inter_mat_type = walls[i].material_type;
         closestNorm.px = norm.px;
         closestNorm.py = norm.py;
         prevlambda = lambda;
      }      
   }


 }

 // Step 2 - Check for intersection against objects in the object array - you must
 //          complete the intersectRay() function, call it, and obtain the closest
 //          intersection (in the forward ray direction) with objects in the scene.
 //          Note that you must provide variables for intersectRay() to return
 //          the point of intersection, normal at intersection, lambda, material type,
 //          and refraction index for the closest object hit by the ray.


 intersectRay(ray, &circInter, &norm, &lambda, &circ_mat_type, &r_dx);

 
 // Step 3 - Check whether the closest intersection with objects is closer than the
 //          closest intersection with a wall. Choose whichever is closer.

   fprintf(stderr,"circle lambda=(%f)\n",lambda);

   if (lambda > TOL) {
      if (lambda < prevlambda) {
         intersection.px = circInter.px;
         intersection.py = circInter.py;
         inter_mat_type = circ_mat_type;
         closestNorm.px = norm.px;
         closestNorm.py = norm.py;
         prevlambda = lambda;
      }
   }

 // Step 4 - Render the ray onto the image. Use renderRay(). Provide renderRay() with
 //          the origin of the ray, and the intersection point (it will then draw a
 //          ray from the origin to the intersection). You also need to provide the
 //          ray's colour.

//   fprintf(stderr,"intersection=(%f,%f)\n",intersection.px,intersection.py);

 renderRay(&ray->p,&intersection,ray->R,ray->G,ray->B);


 // Step 5 - Decide how to handle the ray's bounce at the intersection. You will have
 //          to provide code for 3 cases:
 //          If material type = 0, you have a mirror-reflecting object. 
 //                                Create a ray in the mirror reflection direction,
 //                                with the same colour as the incoming ray, and
 //                                with origin at the intersection point.
 //                                Then call propagateRay() recursively to trace it.
 //          if material type = 1, you have a scattering surface. 
 //                                Choose a random direction within +- 90 degrees 
 //                                from the normal at the intersection. Create a
 //                                ray in this direction, with the same colour as
 //                                the incoming ray, and origin at the intersection,
 //                                then call propagateRay() recursively to trace it.
 //          if material type = 2, you have a refracting (transparent) material.
 // 				   Here you need to process two rays:
 //                                * First, determine how much of the incoming light is
 //                                  reflected and how much is transmitted, using 
 //				     Schlick's approximation:
 // 					 R0 = ((n1-n2)/(n1+n2))^2   
 // 					 R(theta)=R0+((1-R0)*(1-cos(theta))^5)
 //				     If the ray is travelling from air to the inside
 //                                  of an object, n1=1, n2=object's index of refraction.
 //                                  If the ray is travelling from inside an object
 //                                  back onto air, n1=object's index of refraction, n2=1
 //				     And 'theta' is the angle between the normal and the
 // 				     ray direction.
 //				     R(theta) gives the amount Rs of reflected light, 
 //				     1.0-R(theta) gives the amount Rt of transmitted light.
 //                                * Now, make a ray in the mirror-reflection direction
 //				     (same as for material type 0), with the same colour
 //				     as the incoming ray, but with intensity modulated
 //				     by Rs. (e.g. if the incoming's colour is R,G,B,
 //                                  the reflected ray's colour will be R*Rs, G*Rs, B*Rs)
 //				     trace this ray.
 //				   * Make a ray in the refracted-ray direction. The 
 //				     angle for the transmitted ray is given by Snell's law
 //				     n1*sin(theta1) = n2*sin(theta2). The colour of the
 //				     transmitted ray is the same as the incoming ray but
 //			             modulated by Rt. Trace this ray.
 //	That's it! you're done!

   fprintf(stderr,"chosen lambda=(%f)\n",prevlambda);

 // norm has only a point, no direction
 // functions require double pointer and point2D pointer
 if (inter_mat_type == 0) {

   fprintf(stderr,"norm=(%f,%f)\n",closestNorm.px,closestNorm.py);

   // step 1
   angle(&theta, &closestNorm);

//   fprintf(stderr,"theta=(%f)\n",theta);

   theta = -theta;

   // Step 2
   rotation(&theta, &closestNorm);

   fprintf(stderr,"step 2=(%f,%f)\n",closestNorm.px,closestNorm.py);

   // Step 3
   rotation(&theta, &ray->d);

   fprintf(stderr,"step 3=(%f,%f)\n",ray->d.px,ray->d.py);

   // Step 4
   angle(&phi, &ray->d);

   fprintf(stderr,"step 4=(%f)\n",phi);

   phi = -phi;

   // Step 5
   rotation(&phi, &closestNorm);

   fprintf(stderr,"step 5=(%f,%f)\n",closestNorm.px,closestNorm.py);

   // Step 6
   theta = -theta;
   rotation(&theta, &closestNorm);

   fprintf(stderr,"step 6=(%f,%f)\n",closestNorm.px,closestNorm.py);

   // Step 7
   invertDirection(&closestNorm);

   //Step 8
   ray->d.px = closestNorm.px;
   ray->d.py = closestNorm.py;

   normalize(&ray->d);

   fprintf(stderr,"rayD step 7=(%f,%f)\n",ray->d.px,ray->d.py);

   ray->p.px = intersection.px;
   ray->p.py = intersection.py;

   fprintf(stderr,"rayP=(%f,%f)\n",ray->p.px,ray->p.py);


 }
 else if(inter_mat_type == 1) {


 }
 else if(inter_mat_type == 2) {


 }

 depth += 1;

 propagateRay(ray, depth);
   
}

void intersectRay(struct ray2D *ray, struct point2D *p, struct point2D *n, double *lambda, int *type, double *r_idx)
{
 /*
  This function checks for intersection between the ray and any objects in the objects 
  array. The objects are circles, so we are in fact solving for the intersection
  between a ray and a circle.
  
  For a unit circle centered at the origin, we would have the equation
  
  x^2 + y^2 = 1
  
  Using vector notation, with C=[x y]', we get
  
  ||C||^2 = 1
  
  A point on the ray is given by p + lambda*d
  
  Substituting in the equation for the circle we have 
  
  (p + lambda*d)(p + lambda*d) - 1 = 0
  
  If we expand the product above (here the product of two vectors is a DOT product), 
  we can form a quadratic equation
  
  A*lambda^2 + B*lambda + C = 0
  
  Which as you know, has a very simple solution. 
  
  Your task is to 
  * Figure out A, B, and C, considering that your circles don't necessarily have r=1
  * Figure out how to deal with the fact that circles in the scene are likely
    *not* centered at the origin
    
  Then implement the code that will find the value(s) of lambda at the intersection(s).
  
  Note that you can have no intersections, 1 intersection, or 2 intersections
  
  This function *must* find the closest intersection (if any) and update the value
  of lambda, the intersection point p, the normal n at the intersection, 
  the corresponding object's material type (needed outside here to figure out how
  to handle the light's bouncing off this object), and the index of refraction for
  the object (needed if this is a transparent object). 
  
  You must make sure that n is a unit-length vector.
 */
 
 /**********************************************************************************
  * TO DO: Complete this function to find the closest intersection between the
  *        ray and any objects in the scene, as well as the values at the
  *        intersection that will be needed to determine how to bounce/refract the
  *	   ray.
  * *******************************************************************************/


   struct point2D norm1, norm2;
   double A, B, C, Constant, conB, conX, conY, h, k, radius, lambda1, lambda2, prevlambda = INFINITY;

   for (int i = 0; i < MAX_OBJECTS; i++) {

      if(objects[i].r > 0) {

         h = objects[i].c.px;
         k = objects[i].c.py;
         radius = objects[i].r;

         A = dot(&ray->d, &ray->d);

//         fprintf(stderr,"A=(%f)\n",A);

         conB = (((2 * h) * ray->d.px) + ((2 * k) * ray->d.py));

         B = ((2 * dot(&ray->p, &ray->d)) - conB);

//         fprintf(stderr,"B=(%f)\n",B);

         conX = ((h * h) - ((2 * h) * ray->p.px));

         conY = ((k * k) - ((2 * k) * ray->p.py));

         Constant = (conX + conY);

         C = ((dot(&ray->p, &ray->p) - (radius * radius)) + Constant);

//         fprintf(stderr,"C=(%f)\n",C);

//         fprintf(stderr,"h=(%f)\n",h);

//         fprintf(stderr,"h^2=(%f)\n",(h*h));

//         fprintf(stderr,"k=(%f)\n",k);

//         fprintf(stderr,"k^2=(%f)\n",(k*k));

//         fprintf(stderr,"um?=(%f)\n",((B * B) - ((4 * A) * C)));

         lambda1 = ((-B + sqrt((B * B) - ((4 * A) * C))) / (2 * A));

         lambda2 = ((-B - sqrt((B * B) - ((4 * A) * C))) / (2 * A));

         fprintf(stderr,"lambda1=(%f)\n",lambda1);
         fprintf(stderr,"lambda2=(%f)\n",lambda2);

// + k and h maybe need to be removed

//         fprintf(stderr,"norm1=(%f,%f)\n",norm1.px,norm1.py);

       // norm2.px = ray->p.px + (lambda2 * ray->d.px);
       //  norm2.py = ray->p.py + (lambda2 * ray->d.py);
       //  invertDirection(&norm2);
       //  normalize(&norm2);

//         fprintf(stderr,"norm2=(%f,%f)\n",norm2.px,norm2.py);

         if (lambda1 > TOL && lambda1 < lambda2) {

            if (lambda1 < prevlambda) {
               p->px = ray->p.px + lambda2*ray->d.px;
               p->py = ray->p.py + lambda2*ray->d.py;
               type = &objects[i].material_type;

               norm1.px = (2 * (p->px - h));
               norm1.py = (2 * (p->py - k));
               invertDirection(&norm1);
               normalize(&norm1);

               fprintf(stderr,"norm1=(%f,%f)\n",norm1.px,norm1.py);


               n->px = norm1.px;
               n->py = norm1.py;
               *lambda = lambda1;
               r_idx = &objects[i].r_idx;
            }      
         }
         else if (lambda2 > TOL && lambda2 < lambda1) {
            if (lambda2 < prevlambda) {
               p->px = ray->p.px + lambda2*ray->d.px;
               p->py = ray->p.py + lambda2*ray->d.py;
               type = &objects[i].material_type;

               norm2.px = (2 * (p->px - h));
               norm2.py = (2 * (p->py - k));
               invertDirection(&norm2);
               normalize(&norm2);

               fprintf(stderr,"norm2=(%f,%f)\n",norm2.px,norm2.py);

               n->px = norm2.px;
               n->py = norm2.py;
//               lambda = &lambda2;
               *lambda = lambda2;
               r_idx = &objects[i].r_idx;

               fprintf(stderr,"chosen lambda=(%f)\n",*lambda);
            }

         }

      }
   }   
}
