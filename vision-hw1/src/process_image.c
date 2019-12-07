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
		printf("Column index is out of bounds\n");
		return;
	}

	if (!rows_assertion)
	{
		printf("Row index is out of bounds\n");
		return;
	}

	if (!channels_assertion)
	{
		printf("Channel index is out of bounds\n");
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

void scale_image(image im, int c, float v)
{
    // TODO Fill this in
    for (int h=0; h<im.h; h++){
        for (int w=0; w<im.w; w++){
            float value = get_pixel(im, w, h, c);
            set_pixel(im, w, h, c, value*v);
        }
    }
}

void clamp_image(image im)
{
	// TODO Fill this in
	int i;
	for(i=0; i<im.w*im.h*im.c; i++){
		im.data[i] = (im.data[i] > 0.) ? ( (im.data[i] > 1.) ? 1. : im.data[i]) : 0.;
	}
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
	int i;
	float hue, saturation, value, min, max, hue_prime, r, g, b;

	for(i=0; i<im.w*im.h; i++){
		r = im.data[i];
		g = im.data[i+(im.w*im.h)];
		b = im.data[i+(2*im.w*im.h)];
		max = three_way_max(r, g, b);
		min = three_way_min(r, g, b);

		value = max;

		float C = value - min;
		saturation = (C == 0) ? 0 : ((value == 0) ? 0 : C / value);

		hue_prime = (C == 0) ? 0 : ((value == r) ? (g - b) / C : ((value == g) ? ((b - r) / C) + 2.0 : ((r - g) / C) + 4.0));
		hue = (hue_prime < 0) ? (hue_prime / 6.0) + 1 : hue_prime / 6.0;

		im.data[i] = hue;
		im.data[i+(im.w*im.h)] = saturation;
		im.data[i+(2*im.w*im.h)] = value;

	}
}

void hsv_to_rgb(image im)
{
	// TODO Fill this in
	int i;
	float r, g, b, hue, saturation, value, hue_prime, C, max, min;

	for(i=0; i<im.w*im.h; i++){
		hue = im.data[i]; // max value
		saturation = im.data[i+(im.w*im.h)];
		value = im.data[i+(2*im.w*im.h)];

		max = value;
		C = saturation * value; // max - min
		min = max - C;

		hue_prime = hue * 6.;
		float X = (1 - fabs(fmod(hue_prime, 2) - 1)); // Calulating hue values in range [0-1)

		// In this part of HSV hexagon, R component will have the maximum value, B will have the minimum, and G will have value somehwere between R and B
		if ( (hue_prime >= 0.) && (hue_prime < 1.) )
		{
			r = max;
			b = min;
			g = (C * X) + b;

			im.data[i] = r;
			im.data[i+(im.w*im.h)] = g;
			im.data[i+(2*im.w*im.h)] = b;
		}
		else if ( (hue_prime >= 1.) && (hue_prime < 2.) ) // In this part of HSV hexagon, G component will have the maximum value, B will have the minimum, and R will have value somehwere between G and B
		{
			g = max;
			b = min;
			r = (C * X) + b;

			im.data[i] = r;
			im.data[i+(im.w*im.h)] = g;
			im.data[i+(2*im.w*im.h)] = b;
		}
		else if ( (hue_prime >= 2.) && (hue_prime < 3.) ) // Same logic exteded
		{
			g = max;
			r = min;
			b = (C* X) + r;

			im.data[i] = r;
			im.data[i+(im.w*im.h)] = g;
			im.data[i+(2*im.w*im.h)] = b;
		}
		else if ( (hue_prime >= 3.) && (hue_prime < 4.) )
		{
			b = max;
			r = min;
			g = (C* X) + r;

			im.data[i] = r;
			im.data[i+(im.w*im.h)] = g;
			im.data[i+(2*im.w*im.h)] = b;
		}
		else if ( (hue_prime >= 4.) && (hue_prime < 5.) )
		{
			b = max;
			g = min;
			r = (C* X) + g;

			im.data[i] = r;
			im.data[i+(im.w*im.h)] = g;
			im.data[i+(2*im.w*im.h)] = b;
		}
		else if ( (hue_prime >= 5.) && (hue_prime < 6.) )
		{
			r = max;
			g = min;
			b = (C*X) + g;

			im.data[i] = r;
			im.data[i+(im.w*im.h)] = g;
			im.data[i+(2*im.w*im.h)] = b;
		}
		else{
			im.data[i] = 0;
			im.data[i+(im.w*im.h)] = 0;
			im.data[i+(2*im.w*im.h)] = 0;	
		}
	}
}