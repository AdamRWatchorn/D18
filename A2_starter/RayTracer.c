/*
  CSC D18 - RayTracer code.

  Written Dec. 9 2010 - Jan 20, 2011 by F. J. Estrada
  Freely distributable for adacemic purposes only.

  Uses Tom F. El-Maraghi's code for computing inverse
  matrices. You will need to compile together with
  svdDynamic.c

  You need to understand the code provided in
  this file, the corresponding header file, and the
  utils.c and utils.h files. Do not worry about
  svdDynamic.c, we need it only to compute
  inverse matrices.

  You only need to modify or add code in sections
  clearly marked "TO DO" - remember to check what
  functionality is actually needed for the corresponding
  assignment!

  Last updated: Aug. 2017   - F.J.E.
*/

/*****************************************************************************
* COMPLETE THIS TEXT BOX:
*
* 1) Student Name:	Adam Watchorn
* 2) Student Name:	N/A	
*
* 1) Student number:    1000867788
* 2) Student number:    N/A
* 
* 1) UtorID             watchor3
* 2) UtorID             N/A
* 
* We hereby certify that the work contained here is our own
*
* Adam Watchorn                    N/A
* (sign with your name)            (sign with your name)
********************************************************************************/

#include "utils.h"	// <-- This includes RayTracer.h

// A couple of global structures and data: An object list, a light list, and the
// maximum recursion depth
struct object3D *object_list;
struct pointLS *light_list;
struct textureNode *texture_list;
int MAX_DEPTH;

void buildScene(void)
{
#include "buildscene.c"		// <-- Import the scene definition! 
}

void rtShade(struct object3D *obj, struct point3D *p, struct point3D *n, struct ray3D *ray, int depth, double a, double b, struct colourRGB *col)
{
 // This function implements the shading model as described in lecture. It takes
 // - A pointer to the first object intersected by the ray (to get the colour properties)
 // - The coordinates of the intersection point (in world coordinates)
 // - The normal at the point
 // - The ray (needed to determine the reflection direction to use for the global component, as well as for
 //   the Phong specular component)
 // - The current racursion depth
 // - The (a,b) texture coordinates (meaningless unless texture is enabled)
 //
 // Returns:
 // - The colour for this ray (using the col pointer)
 //

 struct colourRGB tmp_col;	// Accumulator for colour components
 double R,G,B;			// Colour for the object in R G and B

 // This will hold the colour as we process all the components of
 // the Phong illumination model
 tmp_col.R=0;
 tmp_col.G=0;
 tmp_col.B=0;

 if (obj->texImg==NULL)		// Not textured, use object colour
 {
  R=obj->col.R;
  G=obj->col.G;
  B=obj->col.B;
 }
 else
 {
  // Get object colour from the texture given the texture coordinates (a,b), and the texturing function
  // for the object. Note that we will use textures also for Photon Mapping.
  obj->textureMap(obj->texImg,a,b,&R,&G,&B);
 }

 //////////////////////////////////////////////////////////////
 // TO DO: Implement this function. Refer to the notes for
 // details about the shading model.
 //////////////////////////////////////////////////////////////

 struct pointLS *light_source = light_list;
 struct point3D *s, *c, *m, *p_sh, *n_sh, *n_b;
 struct object3D *obj_sh;
 double *lambda_sh, *a_sh, *b_sh, spec, max_spec, m_a = 1.0;
 int shadowDepth = MAX_DEPTH - 1, i;
 struct ray3D *ray_sh;
 struct colourRGB col_sh;

 s = newPoint(light_source->p0.px, light_source->p0.py, light_source->p0.pz);
 c = newPoint(ray->d.px, ray->d.py, ray->d.pz);
 n_b = newPoint(-n->px, -n->py, -n->pz);

 c->px = -c->px;
 c->py = -c->py;
 c->pz = -c->pz;

 normalize(c); 

 // Be sure to update 'col' with the final colour computed here!

 // Updating ambient component of light
 tmp_col.R += obj->alb.ra*R;
 tmp_col.G += obj->alb.ra*G;
 tmp_col.B += obj->alb.ra*B;

 subVectors(p, s);

 // Need to normalize s below as it is operated on above
 normalize(s);

 // Below checks if both sides can be illuminated and
 // Updates diffuse component of light
 if(obj->frontAndBack == 1) {
  tmp_col.R += obj->alb.rd * R * light_source->col.R * max(0,dot(n_b,s));
  tmp_col.G += obj->alb.rd * G * light_source->col.G * max(0,dot(n_b,s));
  tmp_col.B += obj->alb.rd * B * light_source->col.B * max(0,dot(n_b,s));
 } else {
  tmp_col.R += obj->alb.rd * R * light_source->col.R * max(0,dot(n,s));
  tmp_col.G += obj->alb.rd * G * light_source->col.G * max(0,dot(n,s));
  tmp_col.B += obj->alb.rd * B * light_source->col.B * max(0,dot(n,s));
 }

 // Setup for calculating specular component of light
 spec = (2 * dot(s, n));
 m = newPoint(spec * n->px, spec * n->py, spec * n->pz);
 subVectors(s, m);
 normalize(m);

 max_spec = max(0,dot(m,c));

 // Applies the correct exponent alpha
 for(i = 0; i < obj->shinyness; i++) {
  m_a = m_a * max_spec;
 }


 tmp_col.R += obj->alb.rs * R * light_source->col.R * m_a;
 tmp_col.G += obj->alb.rs * G * light_source->col.G * m_a;
 tmp_col.B += obj->alb.rs * B * light_source->col.B * m_a;

 // Below is for testing specular
// tmp_col.R += obj->alb.rs * m_a;
// tmp_col.G += obj->alb.rs * m_a;
// tmp_col.B += obj->alb.rs * m_a;


 // Make the returned color equal the calculated color
 col->R = tmp_col.R;
 col->G = tmp_col.G;
 col->B = tmp_col.B;

 // Free up allocated memory
 free(s);
 free(c);
 free(n_b);
 free(m);

 return;

}

void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p, struct point3D *n, double *a, double *b)
{
 // Find the closest intersection between the ray and any objects in the scene.
 // Inputs:
 //   *ray    -  A pointer to the ray being traced
 //   *Os     -  'Object source' is a pointer toward the object from which the ray originates. It is used for reflected or refracted rays
 //              so that you can check for and ignore self-intersections as needed. It is NULL for rays originating at the center of
 //              projection
 // Outputs:
 //   *lambda -  A pointer toward a double variable 'lambda' used to return the lambda at the intersection point
 //   **obj   -  A pointer toward an (object3D *) variable so you can return a pointer to the object that has the closest intersection with
 //              this ray (this is required so you can do the shading)
 //   *p      -  A pointer to a 3D point structure so you can store the coordinates of the intersection point
 //   *n      -  A pointer to a 3D pointnewR structure so you can return the normal at the intersection point
 //   *a, *b  -  Pointers toward double variables so you can return the texture coordinates a,b at the intersection point

 /////////////////////////////////////////////////////////////
 // TO DO: Implement this function. See the notes for
 // reference of what to do in here
 /////////////////////////////////////////////////////////////

    struct object3D *obj_search = object_list, *closest_obj = NULL;
    struct point3D *closest_p, *closest_n;
    double closestLambda = INFINITY;

    closest_p = newPoint(0.0,  0.0, 0.0);
    closest_n = newPoint(0.0, 0.0, 0.0);

    // Search entire object list (just removed &s)
    while(obj_search != NULL) {

        // Only perform operations if object isn't the source
        if(Os != obj_search) {

            // Check for intersections
            obj_search->intersect(obj_search, ray, lambda, p, n, a, b);

            // if sets closest lambda value thus far
            if((*lambda > 0.0) && (*lambda < closestLambda)) {
                closestLambda = *lambda;
                closest_obj = obj_search;

                closest_p->px = p->px;
                closest_p->py = p->py;
                closest_p->pz = p->pz;

                closest_n->px = n->px;
                closest_n->py = n->py;
                closest_n->pz = n->pz;

            }
        }

        obj_search = obj_search->next;
    }

    *lambda = closestLambda;
    *obj = closest_obj;

    p->px = closest_p->px;
    p->py = closest_p->py;
    p->pz = closest_p->pz;

    n->px = closest_n->px;
    n->py = closest_n->py;
    n->pz = closest_n->pz;

    free(closest_p);
    free(closest_n);

}

void rayTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os)
{
 // Trace one ray through the scene.
 //
 // Parameters:
 //   *ray   -  A pointer to the ray being traced
 //   depth  -  Current recursion depth for recursive raytracing
 //   *col   - Pointer to an RGB colour structure so you can return the object colour
 //            at the intersection point of this ray with the closest scene object.
 //   *Os    - 'Object source' is a pointer to the object from which the ray 
 //            originates so you can discard self-intersections due to numerical
 //            errors. NULL for rays originating from the center of projection. 
 
 double lambda;		// Lambda at intersection
 double a,b;		// Texture coordinates
 struct object3D *obj;	// Pointer to object at intersection
 struct point3D p;	// Intersection point
 struct point3D n;	// Normal at intersection
 struct colourRGB I;	// Colour returned by shading function

 struct ray3D *ray_sh;
 double lambda_sh;
 double a_sh,b_sh;
 struct object3D *obj_sh;
 struct point3D p_sh;
 struct point3D n_sh;
 struct colourRGB *I_sh;
 struct point3D *ds;
 struct pointLS *light_source = light_list;
 int depth_sh = MAX_DEPTH;

 struct point3D *c, *dg;
 struct ray3D *ray_dg;
 double dx, dy, dz;

 if (depth>MAX_DEPTH)	// Max recursion depth reached. Return invalid colour.
 {

  col->R=-1;
  col->G=-1;
  col->B=-1;
  return;
 }

 ///////////////////////////////////////////////////////
 // TO DO: Complete this function. Refer to the notes
 // if you are unsure what to do here.
 ///////////////////////////////////////////////////////

     depth += 1;

    // Finds closest intersection
    findFirstHit(ray, &lambda, Os, &obj, &p, &n, &a, &b);

    // If there is no intersection, set pixel colour to that of the background
    if(obj == NULL) {
        col->R = 0;
        col->G = 0;
        col->B = 0;
        return;
    } else {
        // Shade component primarily focused on ambient, diffuse, and specular
        rtShade(obj, &p, &n, ray, depth, a, b, col);
    }

    // Builds a ray pointed in the direction of the light source
    ds = newPoint(light_source->p0.px, light_source->p0.py, light_source->p0.pz);
    subVectors(&p, ds);
    ray_sh = newRay(&p, ds);


    findFirstHit(ray_sh, &lambda_sh, obj, &obj_sh, &p_sh, &n_sh, &a_sh, &b_sh);
    

    // checks if there is an intersection in the shadow direction
    if(obj_sh != NULL) {
        // checks if the intersection is before the light source
        if(lambda_sh < 1 && lambda_sh > 0) {
            col->R = obj->alb.ra*obj->col.R;
            col->G = obj->alb.ra*obj->col.G;
            col->B = obj->alb.ra*obj->col.B;
        }
    }

    // Below code checks for correct transformed normal direction vectors
    // col->R = (n.px + 1)/2;
    // col->G = (n.py + 1)/2;
    // col->B = (n.pz + 1)/2;

    c = newPoint(ray->d.px, ray->d.py, ray->d.pz);

    c->px = -c->px;
    c->py = -c->py;
    c->pz = -c->pz;

    normalize(c);

    dx = ((2 * dot(c,&n)) * n.px);
    dy = ((2 * dot(c,&n)) * n.py);
    dz = ((2 * dot(c,&n)) * n.pz);

    dg = newPoint(dx, dy, dz);
    subVectors(c, dg);

    ray_dg = newRay(&p, dg);

    rayTrace(ray_dg, depth, &I, obj);

    col->R += obj->alb.rg * I.R;
    col->G += obj->alb.rg * I.G;
    col->B += obj->alb.rg * I.B;

    if(col->R > 1) {
        col->R = 1;
    }

    if(col->G > 1) {
        col->G = 1;
    }

    if(col->B > 1) {
        col->B = 1;
    }

    free(ds);
    free(ray_sh);
    free(c);
    free(dg);
    free(ray_dg);

}

int main(int argc, char *argv[])
{
 // Main function for the raytracer. Parses input parameters,
 // sets up the initial blank image, and calls the functions
 // that set up the scene and do the raytracing.
 struct image *im;	// Will hold the raytraced image
 struct view *cam;	// Camera and view for this scene
 int sx;		// Size of the raytraced image
 int antialiasing;	// Flag to determine whether antialiaing is enabled or disabled
 char output_name[1024];	// Name of the output file for the raytraced .ppm image
 struct point3D e;		// Camera view parameters 'e', 'g', and 'up'
 struct point3D g;
 struct point3D up;
 double du, dv;			// Increase along u and v directions for pixel coordinates
 struct point3D pc,d;		// Point structures to keep the coordinates of a pixel and
				// the direction or a ray
 struct ray3D *ray;		// Structure to keep the ray from e to a pixel
 struct colourRGB col;		// Return colour for raytraced pixels
 struct colourRGB background;   // Background colour
 int i,j;			// Counters for pixel coordinates
 unsigned char *rgbIm;

 if (argc<5)
 {
  fprintf(stderr,"RayTracer: Can not parse input parameters\n");
  fprintf(stderr,"USAGE: RayTracer size rec_depth antialias output_name\n");
  fprintf(stderr,"   size = Image size (both along x and y)\n");
  fprintf(stderr,"   rec_depth = Recursion depth\n");
  fprintf(stderr,"   antialias = A single digit, 0 disables antialiasing. Anything else enables antialiasing\n");
  fprintf(stderr,"   output_name = Name of the output file, e.g. MyRender.ppm\n");
  exit(0);
 }
 sx=atoi(argv[1]);
 MAX_DEPTH=atoi(argv[2]);
 if (atoi(argv[3])==0) antialiasing=0; else antialiasing=1;
 strcpy(&output_name[0],argv[4]);

 fprintf(stderr,"Rendering image at %d x %d\n",sx,sx);
 fprintf(stderr,"Recursion depth = %d\n",MAX_DEPTH);
 if (!antialiasing) fprintf(stderr,"Antialising is off\n");
 else fprintf(stderr,"Antialising is on\n");
 fprintf(stderr,"Output file name: %s\n",output_name);

 object_list=NULL;
 light_list=NULL;
 texture_list=NULL;

 // Allocate memory for the new image
 im=newImage(sx, sx);
 if (!im)
 {
  fprintf(stderr,"Unable to allocate memory for raytraced image\n");
  exit(0);
 }
 else rgbIm=(unsigned char *)im->rgbdata;

 ///////////////////////////////////////////////////
 // TO DO: You will need to implement several of the
 //        functions below. For Assignment 2, you can use
 //        the simple scene already provided. But
 //        for Assignment 3 you need to create your own
 //        *interesting* scene.
 ///////////////////////////////////////////////////
 buildScene();		// Create a scene. This defines all the
			// objects in the world of the raytracer

 //////////////////////////////////////////
 // TO DO: For Assignment 2 you can use the setup
 //        already provided here. For Assignment 3
 //        you may want to move the camera
 //        and change the view parameters
 //        to suit your scene.
 //////////////////////////////////////////

 // Mind the homogeneous coordinate w of all vectors below. DO NOT
 // forget to set it to 1, or you'll get junk out of the
 // geometric transformations later on.

 // Camera center is at (0,0,-1)
 e.px=0;
 e.py=0;
 e.pz=-1;
 e.pw=1;

 // To define the gaze vector, we choose a point 'pc' in the scene that
 // the camera is looking at, and do the vector subtraction pc-e.
 // Here we set up the camera to be looking at the origin.
 g.px=0-e.px;
 g.py=0-e.py;
 g.pz=0-e.pz;
 g.pw=1;
 // In this case, the camera is looking along the world Z axis, so
 // vector w should end up being [0, 0, -1]

 // Define the 'up' vector to be the Y axis
 up.px=0;
 up.py=1;
 up.pz=0;
 up.pw=1;

 // Set up view with given the above vectors, a 4x4 window,
 // and a focal length of -1 (why? where is the image plane?)
 // Note that the top-left corner of the window is at (-2, 2)
 // in camera coordinates.
 cam=setupView(&e, &g, &up, -1, -2, 2, 4);

 if (cam==NULL)
 {
  fprintf(stderr,"Unable to set up the view and camera parameters. Our of memory!\n");
  cleanup(object_list,light_list, texture_list);
  deleteImage(im);
  exit(0);
 }

 // Set up background colour here
 background.R=0;
 background.G=0;
 background.B=0;

 // Do the raytracing
 //////////////////////////////////////////////////////
 // TO DO: You will need code here to do the raytracing
 //        for each pixel in the image. Refer to the
 //        lecture notes, in particular, to the
 //        raytracing pseudocode, for details on what
 //        to do here. Make sure you undersand the
 //        overall procedure of raytracing for a single
 //        pixel.
 //////////////////////////////////////////////////////
 du=cam->wsize/(sx-1);		// du and dv. In the notes in terms of wl and wr, wt and wb,
 dv=-cam->wsize/(sx-1);		// here we use wl, wt, and wsize. du=dv since the image is
				// and dv is negative since y increases downward in pixel
				// coordinates and upward in camera coordinates.

 fprintf(stderr,"View parameters:\n");
 fprintf(stderr,"Left=%f, Top=%f, Width=%f, f=%f\n",cam->wl,cam->wt,cam->wsize,cam->f);
 fprintf(stderr,"Camera to world conversion matrix (make sure it makes sense!):\n");
 printmatrix(cam->C2W);
 fprintf(stderr,"World to camera conversion matrix:\n");
 printmatrix(cam->W2C);
 fprintf(stderr,"\n");

 fprintf(stderr,"Rendering row: ");
 for (j=0;j<sx;j++)		// For each of the pixels in the image
 {
  fprintf(stderr,"%d/%d, ",j,sx);
  for (i=0;i<sx;i++)
  {
    ///////////////////////////////////////////////////////////////////
    // TO DO - complete the code that should be in this loop to do the
    //         raytracing!
    ///////////////////////////////////////////////////////////////////

    // Set up (x,y) coordinates for pixel in camera coordinates
    double pixel_x = (cam->wl + (i * du));
    double pixel_y = (cam->wt + (j * dv));
    int depth = 0;
    struct point3D *pc_hold, *d_hold;

    // Create new points at the above coordinates
    pc_hold = newPoint(pixel_x, pixel_y, cam->f);
    d_hold = newPoint(pixel_x, pixel_y, cam->f);

    pc = *pc_hold;
    d = *d_hold;

    // Convert camera coordinates to world coordinates
    matVecMult(cam->C2W, &pc);
    matVecMult(cam->C2W, &d);

    // Translates d to align point with its true location
    subVectors(&(cam->e),&d);

    // Create ray that leaves camera and intersects the plane at pc
    ray = newRay(&pc, &d);

    // Trace the ray
    rayTrace(ray, depth, &col, NULL);

    *(rgbIm+((i+(j*sx))*3)+0)=(unsigned char)(255*col.R);
    *(rgbIm+((i+(j*sx))*3)+1)=(unsigned char)(255*col.G);
    *(rgbIm+((i+(j*sx))*3)+2)=(unsigned char)(255*col.B);

    free(pc_hold);
    free(d_hold);
    free(ray);

  } // end for i
 } // end for j


 fprintf(stderr,"\nDone!\n");

 // Output rendered image
 imageOutput(im,output_name);

 // Exit section. Clean up and return.
 cleanup(object_list,light_list,texture_list);		// Object, light, and texture lists
 deleteImage(im);					// Rendered image
 free(cam);						// camera view
 exit(0);
}

