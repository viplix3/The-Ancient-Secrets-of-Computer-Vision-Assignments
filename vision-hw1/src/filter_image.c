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
// Answer: TODO

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: TODO

image make_gaussian_filter(float sigma)
{
    // TODO
    return make_image(1,1,1);
}

image add_image(image a, image b)
{
    // TODO
    return make_image(1,1,1);
}

image sub_image(image a, image b)
{
    // TODO
    return make_image(1,1,1);
}

image make_gx_filter()
{
    // TODO
    return make_image(1,1,1);
}

image make_gy_filter()
{
    // TODO
    return make_image(1,1,1);
}

void feature_normalize(image im)
{
    // TODO
}

image *sobel_image(image im)
{
    // TODO
    return calloc(2, sizeof(image));
}

image colorize_sobel(image im)
{
    // TODO
    return make_image(1,1,1);
}
