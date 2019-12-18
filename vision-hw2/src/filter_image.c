#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

void l1_normalize(image im)
{
    // TODO
    float pixel_sum = 0;
    int i, j, k;
    for(i=0; i<im.w; i++){
        for(j=0; j<im.h; j++){
            for(k=0; k<im.c; k++){
                    pixel_sum += get_pixel(im, i, j, k);
            }
        }
    }

    for(i=0; i<im.c; i++)
        scale_image(im, i, 1.0/pixel_sum);
}

image make_box_filter(int w)
{
    // TODO
    image box_filter = make_image(w, w, 1);
    int i, j, k;

    for(i=0; i<box_filter.w; i++){
        for(j=0; j<box_filter.h; j++){
            for(k=0; k<box_filter.c; k++)
                set_pixel(box_filter, i, j, k, 1.0);
        }
    }

    l1_normalize(box_filter);

    return box_filter;
}

image convolve_image(image im, image filter, int preserve)
{
    // TODO
    assert(filter.c == im.c || filter.c == 1); //Assertion check to make sure filter is of expected channels

    // Conditional check for number of channels between image and filters
    int out_channels = (preserve) ? (im.c) : (1);

    // This will hold the convolved image
    image convolved_image = make_image(im.w, im.h, out_channels);

    // Convolution loop
    int filter_i, filter_j, image_i, image_j, channels;
    float pixel_value;

    for(image_j=0; image_j<im.h; image_j++)
    {
        for(image_i=0; image_i<im.w; image_i++)
        {

            for(filter_j=0; filter_j<filter.h; filter_j++){
                for(filter_i=0; filter_i<filter.w; filter_i++){
                        for(channels=0; channels<im.c; channels++){

                            int pixel_w, pixel_h;
                            pixel_w = image_i+filter_i-(filter.w/2);
                            pixel_h = image_j+filter_j-(filter.h/2);

                            pixel_value = get_pixel(convolved_image, image_i, image_j, channels);

                            if(filter.c == 1)
                                pixel_value += get_pixel(filter, filter_i, filter_j, 0) * get_pixel(im, pixel_w, pixel_h, channels);
                            
                            else
                                pixel_value += get_pixel(filter, filter_i, filter_j, channels) * get_pixel(im, pixel_w, pixel_h, channels);

                            if(preserve)
                                set_pixel(convolved_image, image_i, image_j, channels, pixel_value);
                            else
                                set_pixel(convolved_image, image_i, image_j, 0, pixel_value);
                    }
                }
            }
        }

    }

    return convolved_image;
}

image make_highpass_filter()
{
    // TODO
    image high_pass_filter = make_image(3, 3, 1);
    set_pixel(high_pass_filter, 1, 0, 0, -1);
    set_pixel(high_pass_filter, 0, 1, 0, -1);
    set_pixel(high_pass_filter, 2, 1, 0, -1);
    set_pixel(high_pass_filter, 1, 2, 0, -1);
    set_pixel(high_pass_filter, 1, 1, 0, 4);

    return high_pass_filter;
}

image make_sharpen_filter()
{
    // TODO
    image sharpen_filter = make_image(3, 3, 1);
    set_pixel(sharpen_filter, 1, 0, 0, -1);
    set_pixel(sharpen_filter, 0, 1, 0, -1);
    set_pixel(sharpen_filter, 2, 1, 0, -1);
    set_pixel(sharpen_filter, 1, 2, 0, -1);
    set_pixel(sharpen_filter, 1, 1, 0, 5);

    return sharpen_filter;
}

image make_emboss_filter()
{
    // TODO
    image emboss_filter = make_image(3, 3, 1);
    set_pixel(emboss_filter, 0, 0, 0, -2);
    set_pixel(emboss_filter, 1, 0, 0, -1);
    set_pixel(emboss_filter, 0, 1, 0, -1);
    set_pixel(emboss_filter, 1, 1, 0, 1);
    set_pixel(emboss_filter, 2, 1, 0, 1);
    set_pixel(emboss_filter, 1, 2, 0, 1);
    set_pixel(emboss_filter, 2, 2, 0, 2);

    return emboss_filter;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: Sharpen and Emboss filters will require the preserve to be set True as we are expecting coloured output from these filters.

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: Yes, high_pass filter will require some post-processing to get rid of false activations of edges.
float compute_gaussian(int x, int y, float sigma)
{
    float gaussian = (1/(TWOPI*pow(sigma, 2))) * (exp(-(pow(x, 2) + pow(y, 2))/(2*pow(sigma, 2))));
    return gaussian;
}

image make_gaussian_filter(float sigma)
{
    // TODO
    int filter_size = 6*sigma;
    filter_size = (filter_size % 2 == 0) ? (filter_size+1) : (filter_size);
    image gaussian_filter = make_image(filter_size, filter_size, 1);
    int i, j;

    for(i=0; i<gaussian_filter.h; i++){
        for(j=0; j<gaussian_filter.w; j++){
            set_pixel(gaussian_filter, j, i, 0, compute_gaussian((gaussian_filter.w/2)-j, (gaussian_filter.h/2)-i, sigma));
        }
    }


    l1_normalize(gaussian_filter);

    return gaussian_filter;
}

image add_image(image a, image b)
{
    // TODO
    assert(a.w == b.w);
    assert(a.h == b.h);

    image added_image = make_image(a.w, a.h, a.c);
    int i, j, k;

    for(i=0; i<a.w; i++){
        for(j=0; j<a.h; j++){
            for(k=0; k<a.c; k++)
                set_pixel(added_image, i, j, k, get_pixel(a, i, j, k) + get_pixel(b, i, j, k));
        }
    }
    return added_image;
}

image sub_image(image a, image b)
{
    // TODO
    assert(a.w == b.w);
    assert(a.h == b.h);

    image subtracted_image = make_image(a.w, a.h, a.c);
    int i, j, k;

    for(i=0; i<a.w; i++){
        for(j=0; j<a.h; j++){
            for(k=0; k<a.c; k++)
                set_pixel(subtracted_image, i, j, k, get_pixel(a, i, j, k) - get_pixel(b, i, j, k));
        }
    }
    return subtracted_image;
}

image make_gx_filter()
{
    image gx_filter = make_image(3, 3, 1);

    set_pixel(gx_filter, 0, 0, 0, -1);
    set_pixel(gx_filter, 0, 1, 0, -2);
    set_pixel(gx_filter, 0, 2, 0, -1);
    
    set_pixel(gx_filter, 2, 0, 0, 1);
    set_pixel(gx_filter, 2, 1, 0, 2);
    set_pixel(gx_filter, 2, 2, 0, 1);

    return gx_filter;
}

image make_gy_filter()
{
    // TODO
    image gy_filter = make_image(3, 3, 1);

    set_pixel(gy_filter, 0, 0, 0, -1);
    set_pixel(gy_filter, 1, 0, 0, -2);
    set_pixel(gy_filter, 2, 0, 0, -1);
    
    set_pixel(gy_filter, 0, 2, 0, 1);
    set_pixel(gy_filter, 1, 2, 0, 2);
    set_pixel(gy_filter, 2, 2, 0, 1);
    
    return gy_filter;
}

void feature_normalize(image im)
{
    int p;
    float min = INFINITY;
    float max = -1.0;

    for(p = 0; p < im.w * im.h * im.c; p++) {
       float pix = im.data[p];
       if(pix < min) {
           min = pix;
       } 

       if(pix > max) {
           max = pix;
       }
    }

    int channel, row, col;
    for(channel = 0; channel < im.c; channel++){
        for(row = 0; row < im.h; row++) {
            for(col = 0; col < im.w; col++) {
                if(max - min == 0 ) {
                    set_pixel(im, col, row, channel, 0.0);
                } else {
                    float x = get_pixel(im, col, row, channel);
                    float val = (x - min)/(max - min);
                    set_pixel(im, col, row, channel, val);
                }
            }     
        }
    }
    
}

image *sobel_image(image im)
{
    // TODO
    image gx_filter = make_gx_filter();
    image gy_filter = make_gy_filter();
    image *result_image = calloc(2, sizeof(image));

    image GX = convolve_image(im, gx_filter, 0);
    image GY = convolve_image(im, gy_filter, 0);

    result_image[0] = make_image(GX.w, GX.h, GX.c);
    result_image[1] = make_image(GX.w, GX.h, GX.c);

    int i, j, k;

    
    for(i=0; i<GX.w; i++){
        for(j=0; j<GX.h; j++){
            for(k=0; k<GX.c; k++){
                set_pixel(result_image[0], i, j, k, sqrt(pow(get_pixel(GX, i, j, k), 2) + pow(get_pixel(GY, i, j, k), 2)));
                set_pixel(result_image[1], i, j, k, atan2f(get_pixel(GY, i, j, k), get_pixel(GX, i, j, k)));
            }
        }
    }

    free_image(gx_filter);
    free_image(gy_filter);
    free_image(GX);
    free_image(GY);
    return result_image;
}

image colorize_sobel(image im)
{
    
    image *pimage = sobel_image(im);
    image magnitude = pimage[0];
    feature_normalize(magnitude);
    clamp_image(magnitude);
    image direction = pimage[1];
    image color = make_image(im.w, im.h, 3);
    int i;

    for (i=0; i<im.h * im.w; i++) {
        // hue
        color.data[0 * im.w * im.h + i] = direction.data[i];
        // saturation
        color.data[1 * im.w * im.h + i] = 10*magnitude.data[i];
        // value
        color.data[2 * im.w * im.h + i] = 0.2 - magnitude.data[i];
    }
    free_image(magnitude);
    free_image(direction);
    hsv_to_rgb(color);
    return color;
}
