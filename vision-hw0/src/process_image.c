#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

float get_pixel(image im, int x, int y, int c)
{
    // TODO Fill this in
    // x = column coordinate, y = row coordinate, c = channel coordinate
    int columns= im.w;
    int rows = im.h;
    int channels = im.c;
    x = ((x < columns) ? ( (x <= 0) ? 0 : x) : columns-1);
    y = ((y < rows) ? ( (y <= 0) ? 0 : y) : rows-1);
    c = ((c < channels) ? ( (c <= 0) ? 0 : c) : channels-1);

    float pixel_value = im.data[c*columns*rows + y*columns + x];

    return pixel_value;
}

void set_pixel(image im, int x, int y, int c, float v)
{
    // TODO Fill this in
    int columns= im.w;
    int rows = im.h;
    int channels = im.c;
    
    int columns_assertion = (x < columns) ? ( (x < 0) ? 0 : 1) : 0;
    int rows_assertion = (y < rows) ? ( (y < 0) ? 0 : 1) : 0;
    int channels_assertion = (c < channels) ? ( (c < 0) ? 0 : 1) : 0;

    if (!columns_assertion)
    {
        printf("Column index is out of bounds");
        return;
    }

    if (!rows_assertion)
    {
        printf("Row index is out of bounds");
        return;
    }

    if (!channels_assertion)
    {
        printf("Channel index is out of bounds");
        return;
    }

    im.data[c*columns*rows + y*columns + x] = v;

}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    // TODO Fill this in
    int i;

    for(i=0; i<im.w*im.h*im.c; i++){
        copy.data[i] = im.data[i];
    }

    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    // TODO Fill this in
    
    int i;
    float weight_r = 0.299;
    float weight_g = 0.587;
    float weight_b = 0.114;

    for(i=0; i<gray.w*gray.h*gray.c; i++){
        gray.data[i] = weight_r*im.data[i] + weight_g*im.data[i+(im.w*im.h)] + weight_b*im.data[i+(2*im.w*im.h)];
    }
    
    return gray;
}

void shift_image(image im, int c, float v)
{
    // TODO Fill this in
    int i;
    for(i=0; i<im.w*im.h; i++){
        im.data[c*im.w*im.h + i] = v + im.data[c*im.w*im.h + i];
    }
}

void clamp_image(image im)
{
    // TODO Fill this in
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    // TODO Fill this in
}

void hsv_to_rgb(image im)
{
    // TODO Fill this in
}
