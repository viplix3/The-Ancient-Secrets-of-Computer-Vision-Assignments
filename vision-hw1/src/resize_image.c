#include <math.h>
#include "image.h"

float nn_interpolate(image im, float x, float y, int c)
{
    // TODO Fill in
   	/* Logic for nearest neighbour interpolation
   		round() C function of math library will return the closest integer to the given float number
   	*/
    return get_pixel(im, round(x), round(y), c);
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix that first line)
    float w_scaling_factor, w_displacement_factor, h_scaling_factor, h_displacement_factor;
    int i, j, k;

    image resized_image = make_image(w, h, im.c);

    /* Logic behind new image to original image pixel coordinate mapping
    	Step 1.
    	If original image is of size 4x4x3 and resized image is tobe of size 7x7x3, then:
    	7 pixels of new image = 4 pixels of original image
    	1 pixel of new image = 4/7 pixel of original image
    	Step 2.
		Aprt from this, we will also have to take into account the fact the the pixel actually at (-0.5, -0.5) instead of (0, 0)
    	This logic will be used to map pixel coordinate of new image to old image, and then will be interpolated using nn_interpolate function to get exact pixel values of new image.
    */
    w_scaling_factor = (float)im.w/w;
    h_scaling_factor = (float)im.h/h;
    w_displacement_factor = -0.5 + 0.5*w_scaling_factor;
    h_displacement_factor = -0.5 + 0.5*h_scaling_factor;

    for(i=0; i<w; i++){
    	for(j=0; j<h; j++){
    		for(k=0; k<im.c; k++){
    			set_pixel(resized_image, i, j, k, nn_interpolate(im, (i*w_scaling_factor) + w_displacement_factor, (j*h_scaling_factor) + h_displacement_factor, k));
    		}
    	}
    }

    return resized_image;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // TODO
    /* Logic for bilinear interpolation
    	Step 1.
    	Get the floor() and ceil() of given float pixel coordinate of resized image.
    	Step 2.
    	Calculate the distance of floor() from given float val and ceil() from given float val
    	Step 3.
    	Value will be interpolated using 2 intermediate points.
    */
	int floor_x, floor_y, ceil_x, ceil_y;
	float dx1, dx2, dy1, dy2, inter_pixel_1, inter_pixel_2, pixel_value;
	
	floor_x = floor(x);
	floor_y = floor(y);
	ceil_x = ceil(x);
	ceil_y = ceil(y);

	dx1 = (x - floor_x);
	dx2 = (ceil_x - x);
	dy1 = (y - floor_y);
	dy2 = (ceil_y - y);

	inter_pixel_1 = (dx1 * get_pixel(im, ceil_x, ceil_y, c)) + (dx2 * get_pixel(im, floor_x, ceil_y, c));
	inter_pixel_2 = (dx1 * get_pixel(im, ceil_x, floor_y, c)) + (dx2 * get_pixel(im, floor_x, floor_y, c));
	pixel_value = (dy2 * inter_pixel_2) + (dy1 * inter_pixel_1);

    return pixel_value;
}

image bilinear_resize(image im, int w, int h)
{
    // TODO
    float w_scaling_factor, w_displacement_factor, h_scaling_factor, h_displacement_factor;
    int i, j, k;

    image resized_image = make_image(w, h, im.c);

    /* Logic behind new image to original image pixel coordinate mapping
    	Step 1.
    	If original image is of size 4x4x3 and resized image is tobe of size 7x7x3, then:
    	7 pixels of new image = 4 pixels of original image
    	1 pixel of new image = 4/7 pixel of original image
    	Step 2.
		Aprt from this, we will also have to take into account the fact the the pixel actually at (-0.5, -0.5) instead of (0, 0)
    	This logic will be used to map pixel coordinate of new image to old image, and then will be interpolated using bilinear_interpolate function to get exact pixel values of new image.
    */
    w_scaling_factor = (float)im.w/w;
    h_scaling_factor = (float)im.h/h;
    w_displacement_factor = -0.5 + 0.5*w_scaling_factor;
    h_displacement_factor = -0.5 + 0.5*h_scaling_factor;

    for(i=0; i<w; i++){
    	for(j=0; j<h; j++){
    		for(k=0; k<im.c; k++){
    			set_pixel(resized_image, i, j, k, bilinear_interpolate(im, (i*w_scaling_factor) + w_displacement_factor, (j*h_scaling_factor) + h_displacement_factor, k));
    		}
    	}
    }

    return resized_image;
}

