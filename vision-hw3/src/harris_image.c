#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"
#include <time.h>

// Frees an array of descriptors.
// descriptor *d: the array.
// int n: number of elements in array.
void free_descriptors(descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        free(d[i].data);
    }
    free(d);
}

// Create a feature descriptor for an index in an image.
// image im: source image.
// int i: index in image for the pixel we want to describe.
// returns: descriptor for that index.
descriptor describe_index(image im, int i)
{
    int w = 5;
    descriptor d;
    d.p.x = i%im.w;
    d.p.y = i/im.w;
    d.data = calloc(w*w*im.c, sizeof(float));
    d.n = w*w*im.c;
    int c, dx, dy;
    int count = 0;
    // If you want you can experiment with other descriptors
    // This subtracts the central value from neighbors
    // to compensate some for exposure/lighting changes.
    for(c = 0; c < im.c; ++c){
        float cval = im.data[c*im.w*im.h + i];
        for(dx = -w/2; dx < (w+1)/2; ++dx){
            for(dy = -w/2; dy < (w+1)/2; ++dy){
                float val = get_pixel(im, i%im.w+dx, i/im.w+dy, c);
                d.data[count++] = cval - val;
            }
        }
    }
    return d;
}

// Marks the spot of a point in an image.
// image im: image to mark.
// ponit p: spot to mark in the image.
void mark_spot(image im, point p)
{
    int x = p.x;
    int y = p.y;
    int i;
    for(i = -9; i < 10; ++i){
        set_pixel(im, x+i, y, 0, 1);
        set_pixel(im, x, y+i, 0, 1);
        set_pixel(im, x+i, y, 1, 0);
        set_pixel(im, x, y+i, 1, 0);
        set_pixel(im, x+i, y, 2, 1);
        set_pixel(im, x, y+i, 2, 1);
    }
}

// Marks corners denoted by an array of descriptors.
// image im: image to mark.
// descriptor *d: corners in the image.
// int n: number of descriptors to mark.
void mark_corners(image im, descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        mark_spot(im, d[i].p);
    }
}


float compute_1d_gaussian(int x, float sigma)
{
    float gaussian = (1/(TWOPI*pow(sigma, 2))) * (exp(-(pow(x, 2))/(2*pow(sigma, 2))));
    return gaussian;
}


// Creates a 1d Gaussian filter.
// float sigma: standard deviation of Gaussian.
// returns: single row image of the filter.
image make_1d_gaussian(float sigma)
{
    // TODO: optional, make separable 1d Gaussian.

    int filter_size = 6*sigma;
    filter_size = (filter_size % 2 == 0) ? (filter_size+1) : (filter_size);
    image gaussian_1d_filter = make_image(filter_size, 1, 1);
    int j;

    for(j=0; j<gaussian_1d_filter.w; j++){
        set_pixel(gaussian_1d_filter, j, 0, 0, compute_1d_gaussian((gaussian_1d_filter.w/2)-j, sigma));
    }

    l1_normalize(gaussian_1d_filter);

    return gaussian_1d_filter;
}

image transpose_image(image im)
{
    image transposed_im = make_image(im.h, im.w, im.c);
    int i, j, k;

    for(i=0; i<im.h; i++){
        for(j=0; j<im.w; j++){
            for(k=0; k<im.c; k++)
                set_pixel(transposed_im, i, j, k, get_pixel(im, j, i, k));
        }
    }

    return transposed_im;
}

// Smooths an image using separable Gaussian filter.
// image im: image to smooth.
// float sigma: std dev. for Gaussian.
// returns: smoothed image.
image smooth_image(image im, float sigma)
{
    if(0){
        image g = make_gaussian_filter(sigma);
        image s = convolve_image(im, g, 1);
        free_image(g);
        return s;
    } else {
        // TODO: optional, use two convolutions with 1d gaussian filter.
        // If you implement, disable the above if check.
        image g_nx1 = make_1d_gaussian(sigma);
        image g_1xn = transpose_image(g_nx1);

        image s = convolve_image(im, g_1xn, 1);
        s = convolve_image(s, g_nx1, 1);
        free_image(g_nx1);
        free_image(g_1xn);
        return s;
    }
}

// Calculate the structure matrix of an image.
// image im: the input image.
// float sigma: std dev. to use for weighted sum.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          third channel is IxIy.
image structure_matrix(image im, float sigma)
{
    image S = make_image(im.w, im.h, 3);
    // TODO: calculate structure matrix for im.

    image gx_filter = make_gx_filter(); // Making sobel filter to give x-directional gradient
    image gy_filter = make_gy_filter(); // Making sobel filter to give y-directional gradient
    // image gaussian_filter = make_gaussian_filter(sigma);
    int i, j, k;
    float pixel_value;

    image Ix = convolve_image(im, gx_filter, 0); // Convolving input image with x-directional sobel filter to get Ix (derivative along x-axis)
    image Iy = convolve_image(im, gy_filter, 0); // Convolving input image with y-directional sobel filter to get Iy (derivative along y-axis)

    // Setting pixel values of S such that 1st channel is Ix^2, 2nd channel is Iy^2 and 3rd channel is IxIy
    for(i=0; i<S.w; i++){
        for(j=0; j<S.h; j++){
            for(k=0; k<S.c; k++){
                if(k==0){
                    pixel_value = pow(get_pixel(Ix, i, j, 0), 2);
                    set_pixel(S, i, j, k, pixel_value);
                }
                else if(k==1){
                    pixel_value = pow(get_pixel(Iy, i, j, 0), 2);
                    set_pixel(S, i, j, k, pixel_value);
                }
                else{
                    pixel_value = get_pixel(Ix, i, j, 0) * get_pixel(Iy, i, j, 0);
                    set_pixel(S, i, j, k, pixel_value);
                }
            }
        }
    }

    // Calculating structure matrix S using weighted sun (with the help of gaussian filter having provided sigma)
    // S = convolve_image(S, gaussian_filter, 1);
    S = smooth_image(S, sigma);

    // free_image(gaussian_filter);
    free_image(gx_filter);
    free_image(gy_filter);
    free_image(Ix);
    free_image(Iy);

    return S;
}

// Estimate the cornerness of each pixel given a structure matrix S.
// image S: structure matrix for an image.
// returns: a response map of cornerness calculations.
image cornerness_response(image S)
{
    image R = make_image(S.w, S.h, 1);
    // TODO: fill in R, "cornerness" for each pixel using the structure matrix.
    // We'll use formulation det(S) - alpha * trace(S)^2, alpha = .06.

    int i, j;
    float alpha = 0.06;

    for(i=0; i<R.w; i++){
        for(j=0; j<R.h; j++){
            float Ixx = get_pixel(S, i, j, 0);
            float Iyy = get_pixel(S, i, j, 1);
            float Ixy = get_pixel(S, i, j, 2);
            float cornerness = ((Ixx * Iyy) - pow(Ixy, 2)) - (alpha * pow((Ixx + Iyy), 2)); // Used normal determinanat and trace formulas
            set_pixel(R, i, j, 0, cornerness);
        }
    }

    return R;
}

// Perform non-max supression on an image of feature responses.
// image im: 1-channel image of feature responses.
// int w: distance to look for larger responses.
// returns: image with only local-maxima responses within w pixels.
image nms_image(image im, int w)
{
    image r = copy_image(im);
    // TODO: perform NMS on the response map.
    // for every pixel in the image:
    //     for neighbors within w:
    //         if neighbor response greater than pixel response:
    //             set response to be very low (I use -999999 [why not 0??])

    int i, j, k, l;
    for(i=0; i<r.w; i++){
        for(j=0; j<r.h; j++){
            float current_pixel_response = get_pixel(im, i, j, 0);
            
            for(k=-w; k<w+1; k++){
                for(l=-w; l<w+1; l++){
                    float comparator_pixel_response = get_pixel(im, i+k, j+l, 0);
                    if(comparator_pixel_response > current_pixel_response){
                        set_pixel(r, i, j, 0, -999999);
                        goto done_with_current_pixel;
                    }
                }
            }
            done_with_current_pixel:;
        }
    }

    return r;
}

// Perform harris corner detection and extract features from the corners.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
// int *n: pointer to number of corners detected, should fill in.
// returns: array of descriptors of the corners in the image.
descriptor *harris_corner_detector(image im, float sigma, float thresh, int nms, int *n)
{
    // Calculate structure matrix
    image S = structure_matrix(im, sigma);

    // Estimate cornerness
    image R = cornerness_response(S);

    // Run NMS on the responses
    image Rnms = nms_image(R, nms);


    //TODO: count number of responses over threshold
    int count = 0; // change this
    int i, j;

    for(i=0; i<Rnms.w; i++){
        for(j=0; j<Rnms.h; j++){
            float pixel_response = get_pixel(Rnms, i, j, 0);
            if(pixel_response > thresh)
                count++;
        }
    }
    
    *n = count; // <- set *n equal to number of corners in image.
    descriptor *d = calloc(count, sizeof(descriptor));
    //TODO: fill in array *d with descriptors of corners, use describe_index.
    for (i = 0; i < Rnms.h*Rnms.w; ++i){
      if (Rnms.data[i] > thresh) {
        *d++ = describe_index(im, i);
      }
    }
    d -= count;

    free_image(S);
    free_image(R);
    free_image(Rnms);
    return d;
}

// Find and draw corners on an image.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
void detect_and_draw_corners(image im, float sigma, float thresh, int nms)
{
    int n = 0;
    descriptor *d = harris_corner_detector(im, sigma, thresh, nms, &n);
    mark_corners(im, d, n);
}
