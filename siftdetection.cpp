

#include "siftdetection.h"





#define zooming_func bilinear_zooming

#define max(a,b) (((a)>(b))?(a):(b))



static unsigned int * byte2int (BYTE *image_b, int width, int height)
{
    unsigned int *image_i = (unsigned int *) calloc(width * height, sizeof(unsigned int));
    int i, j;
         
    for( i=0; i <= height-1; i++)
    {
        for( j=0; j <= width-1; j++)
        {
            *(image_i + i*width + j) = *(image_b + i*width + j);
        }
    }

    return image_i;

}

/////////////////////////////////////////////////////////////////

/*  Normalization function to normalize the image to the range of [0, 1]

  for the following processing */

/*   float *normalize(int *image, int width, int height)
    input:  int *image, int width, int height (640 * 480, width = 640, height = 480) 
*/
/////////////////////////////////////////////////////////////////
float *normalize(unsigned int *image, int width, int height)
{
    unsigned int max = *image, min = *image;

    float range;

    unsigned int i,j;

    float *nor_image = (float *)calloc(width * height, sizeof(float *));

    for ( i = 0; i <= height - 1; i++)
    {
        for ( j = 0; j <= width - 1; j++)
        {
            if (*(image + i*width+j) > max)
                max = *(image + i*width+j);
            
            if (*(image + i*width+j) < min)
                min = *(image + i*width+j);

        }
    }

    range = (float)(max - min);


    for ( i = 0; i <= height - 1; i++)
    {
        for (  j = 0; j <= width - 1; j++)
        {

            *(nor_image + i*width+j) = (float)(*(image + i*width+j) - min)/range;    
            
            
        }

    }

    return nor_image;

}



float *filtering(float *image, int width, int height, float sigma, int length)
{
    float sum;
        
    int    i,j,a,b;
    
    float norm = 0.0f;

    //float max = 1.0f;

    //float *temp = image;

    float *filtered_image = (float *)calloc(width*height, sizeof(float *));


    float **gaussianKernel = gaussian_coeff(sigma, length);

    for (i = 0; i <= length-1; i++)
        for (j = 0; j <= length-1; j++)
            norm = norm + gaussianKernel[i][j];


    for(j=2; j < height-2; j++)
    {
        for(i=2; i < width-2; i++)
        {

            sum = 0;
            for(a=-2; a<3; a++){
                for(b=-2; b<3; b++){
                    sum = sum + *(image+(j+a)*width+(i+b)) * gaussianKernel[a+2][b+2];
                }
            }
            //sum = sum/norm;
        
        
            *(filtered_image + j*width + i) = sum;
        }  

    }  
    

    for (i = 0; i <= length-1; i++)
        free(gaussianKernel[i]);
    free(gaussianKernel);

    fix_edges(filtered_image, 2, width, height);

    return filtered_image;
}


 
void fix_edges(float *image, int edge_width, int width, int height)
{
    int i, j;


    for(i = edge_width; i > 0; i--)
    {
        *(image + (i-1)*width + i-1) = *(image + i*width + i);
    
        *(image + (i-1)*width + width - i) = *(image + i*width + width - i -1 );

        *(image + (height-i)*width + i-1) = *(image + (height-i-1)*width + i );

        *(image + (height-i)*width + width - i) = *(image + (height-i-1)*width + width - i -1 );

    
    }  


    for( j = edge_width; j < height-edge_width; j++ )
    {
        for (i = edge_width; i > 0; i-- )
        {
            *(image + j*width + i-1) = *(image + j*width + i);

            *(image + j*width + width - i) = *(image + j*width + width - i - 1);

        }

    }

    for( j = edge_width; j >0; j-- )
    {
        for (i = edge_width; i < width-edge_width; i++ )
        {
            *(image + (j-1)*width + i) = *(image + j*width + i);

            *(image + (height-j)*width + i) = *(image + (height-j-1)*width + i);

        }

    }
} 

/*
    Generating the Gaussian coefficient matrix accoding to the specified standard 
    deviation (float sigma), (int length) is the length of the coefficients, which is fixed
    at 5. 

  */

float **gaussian_coeff(float sigma, int length )
{
    double pi = 3.1415926; 


    
    float **coefficients_2d = (float **) calloc (length, sizeof(float**));
    for (int i = 0; i < length; i++)
    {
        coefficients_2d[i] = (float *) calloc (length, sizeof(float*));
    }

    float *coefficients_1d = (float *) calloc (length, sizeof(float*));

    for ( i = 0; i < length; i++)
    {

        int temp = (int)(i - ceil((length+1)/2) + 1);
        
        *(coefficients_1d + i) =  (float)( exp(-pow(temp, 2)/(2*pow(sigma,2)))/(sigma*sqrt(2*pi)));
    }

    for ( i = 0; i < length; i++)
    {
        for (int j = 0; j < length; j++)
        {
            
            coefficients_2d[i][j] = (float)((*(coefficients_1d + i))*(*(coefficients_1d + j)));
            
        }
    }
    
    
    free(coefficients_1d);
    return coefficients_2d;

}




float bspline_func(const float x)
{
    if (x>2.0f) return 0.0f;

    float a, b, c, d;
    float xm1 = x - 1.0f; 
    float xp1 = x + 1.0f;
    float xp2 = x + 2.0f;

    if ((xp2) <= 0.0f) a = 0.0f; else a = xp2*xp2*xp2; 
    if ((xp1) <= 0.0f) b = 0.0f; else b = xp1*xp1*xp1;
    if (x <= 0) c = 0.0f; else c = x*x*x;  
    if ((xm1) <= 0.0f) d = 0.0f; else d = xm1*xm1*xm1;

    return (0.16666666666666666667f * (a - (4.0f * b) + (6.0f * c) - (4.0f * d)));


}


float *bicubic_zooming(float *image, int width, int height, int new_width, int new_height)
{
    float xScale, yScale;
    
    xScale = (float)width  / (float)new_width;
    yScale = (float)height / (float)new_height;

    float dx, dy;
    float f_x, f_y;
    int i_x, i_y;
    int xx, yy;
    float r1, r2;

    float *newimage = (float *)calloc(new_width*new_height, sizeof(float*));

    float max = 1.0f;

    float temp = 0.0f;
    

    for(int y=0; y<new_height; y++)
    {
        f_y = (float) y * yScale - 0.5f;
        i_y = (int)f_y;
        dy   = f_y - i_y;
        
        for(int x=0; x<new_width; x++)
        {
            f_x = (float) x * xScale - 0.5f;
            i_x = (int) f_x;
            dx   = f_x - i_x;

            for(int m=-1; m<3; m++) 
            {
                r1 = bspline_func(dy - (float)m);
                yy = i_y+m;
                if (yy < 0) yy = 0;
                if (yy >= height) yy = height-1;
                
                for(int n=-1; n<3; n++) 
                {
                    r2 = r1 * bspline_func((float)n - dx);
                    xx = i_x+n;
                    if (xx < 0) xx = 0;
                    if (xx >= width) xx = width-1;
        
                    //*(newimage + y*new_width + x) += *(image + yy * width + xx)*r2; 
                    temp += *(image + yy * width + xx)*r2;
                                        
                }
            }

            //if(temp < 0)   temp = 0.01f;
            //if(temp > max) temp = max;
        
            *(newimage + y*new_width + x) = temp;

            temp = 0.0f;

        }

    }
    return newimage;

}


float *bilinear_zooming(float *image, int width, int height, int new_width, int new_height)
{
    float xScale, yScale, fX, fY;
 
    xScale = (float)width  / (float)new_width;
    yScale = (float)height / (float)new_height;
 
 
    int ifX, ifY, ifX1, ifY1;
    float dx, dy;
 
    int xmax = width  - 1;
    int ymax = height - 1;
 
    float *newimage = (float *)malloc(new_width*new_height*sizeof(float));
 
    for(int y = 0; y < new_height; y++)
    {
            fY = (float)y * yScale;
            ifY = (int)fY;
            ifY1 = min_bi(ymax, ifY+1);
            dy = fY - ifY;
   
            
            for(int x = 0; x < new_width; x++)
            {
                    fX = (float)x * xScale;
                    ifX = (int)fX;
                    ifX1 = min_bi(xmax, ifX+1);
                    dx = fX - ifX;
    
     
                    *(newimage + y*new_width + x) = (*(image + ifY * width + ifX)*(1-dx)*(1-dy) 
                                + *(image + ifY1 * width + ifX)*dy*(1-dx)
                                + *(image + ifY * width+ ifX1)*dx*(1-dy)
                                + *(image + ifY1 * width + ifX1)*dx*dy); 
                    //*(newimage + y*new_width + x) = *(image + ifY * width + ifX);
                   
    
            }  
            
    }
     
    return  newimage;
     
}


int min_bi(int a, int b)
{
    return a<=b?a:b;

}
//const int octaves = 4;
//const int intervals = 2;


float ***sift_extraction(BYTE *image_b, int width,int height, 
                     float contrast_threshold, float curvature_threshold)
{

    //////////////////////////////////////
    // test part
    orientation_calc_test(image_b, width, height);



    //////////////////////////////////////












    int n = 0;
    int m = 0;
    float sub_sample = 0.5;

    /* absolute sigma of the gaussian filtering applying to each image */
    float ***abs_sigma = (float ***)calloc(octaves, sizeof(float ***));
    for(n = 0; n <= octaves-1; n++)
        abs_sigma[n] = (float **)calloc((intervals+3), sizeof(float **));

    for(n = 0; n <= octaves-1; n++)
        for( m = 0; m <= (intervals+2); m++)
            abs_sigma[n][m] = (float *)calloc(1, sizeof(float *));
    /////////////////////////////////////////////////////////////////////



    /* define the data structure for the pyrmiad of the Gaussian and difference of Gaussian */

    
    float ***gaussain_pyrmiad    = (float ***)calloc(octaves, sizeof(float ***));;
    
    for( m = 0; m <= octaves-1; m++)
        gaussain_pyrmiad[m] = (float **)calloc((intervals+3), sizeof(float **));

    //float *differenceofgaussian_pyrmiad[octaves][intervals+2];

    float ***differenceofgaussian_pyrmiad = (float ***)calloc(octaves, sizeof(float ***));

    //for(int m = 0; m <= intervals+1; m++)
    for( m = 0; m <= octaves-1; m++)
        differenceofgaussian_pyrmiad[m] = (float **)calloc((intervals+2), sizeof(float **));


    /* convert the image in BYTE to int format  */
    unsigned int *image = byte2int(image_b, width, height);



    /*  normalized the image to be with range of [0 1] */
    float *nor_image = normalize(image, width, height); 


    /*  release memory */
    free(image);

    /* Prefiltering the image by Gaussian filter with standard deviation of 0.5.  */
    /* Antialias procedure needed for the scaling operation                       */
    
    float antialias_sigma = 0.5;

    int filter_length = 5;

    float *filtered_image = filtering(nor_image, width, height, antialias_sigma, filter_length);

    free(nor_image);

    /* Scaling the prefiltered image by factor of 2. */

    int new_width = width*2;

    int new_height = height*2;

    float *pre_image = zooming_func(filtered_image, width, height, new_width, new_height);


    free(filtered_image);

    
    /* The actual Gaussian filter to be applied to the upsampled image is of standard deviation sqrt(2) */
    /* Since the image has been prefiltered to avoid the alias caused by the upsampling, only the Gaussian */
    /* filtering with the difference of stanard deviation need to be applied to the upsampled image  */

    float preblur_sigma = sqrt(pow(sqrt(2), 2) - pow((2*antialias_sigma), 2));

    //float preblur_sigma = (float)sqrt(pow(1.6, 2) - pow((2*antialias_sigma), 2));


    /* gaussain_pyrmiad[0][0] is the first filtered image of the first octave, upsampled by 2 and filtered by sqrt(20 */
    gaussain_pyrmiad[0][0] = filtering(pre_image, new_width, new_height, preblur_sigma, filter_length);

    free(pre_image);

    
    float initial_sigma = (float)sqrt( pow((2*antialias_sigma), 2) + pow(preblur_sigma, 2) );

    *abs_sigma[0][0] = initial_sigma*sub_sample;
    
    float sigma, sigma_difference;

    
    int i, j, k;
    
    for (i = 0; i <= octaves-1; i++)
    {

        sigma = initial_sigma;

        /* allocate the memeory space for the difference of Gaussian data array */
        for (k = 0; k <= intervals + 1; k++)

            differenceofgaussian_pyrmiad[i][k] = (float *)calloc(new_width*new_height, sizeof(float *));

        for (j = 1; j<= intervals + 2; j++)
        {
            
            sigma_difference = sigma;

            sigma = 1.4142*sigma;

            *abs_sigma[i][j] = sigma*sub_sample;

            /* calculate the successively Gaussian smoothed image */
            gaussain_pyrmiad[i][j] = filtering(gaussain_pyrmiad[i][j-1], 
                                    new_width, new_height, sigma_difference, filter_length);
            
            /* calculated the Difference of Gaussian image*/
            for (int column = 0; column <= new_height-1; column++)
            {
                for (int row = 0; row <= new_width-1; row++)
                {
                                
                    *(differenceofgaussian_pyrmiad[i][j-1]+column*new_width+row) =  
                        *(gaussain_pyrmiad[i][j]+column*new_width+row) - *(gaussain_pyrmiad[i][j-1]+column*new_width+row);
                }
            }

        }

        /* downsampling the third Gaussian filtered image gaussain_pyrmiad[i][intervals] by 2 to get */
        /* the first Gaussian filtered image in the next octave gaussain_pyrmiad[i+1][0] */
        /* then go back to generate the remaing 4 filtered images and the difference of Gaussian images */
        if(i < octaves-1)
        {
            new_width = new_width/2;

            new_height = new_height/2;

            gaussain_pyrmiad[i+1][0] = zooming_func(gaussain_pyrmiad[i][intervals], new_width*2, new_height*2, 
                                                                                    new_width, new_height); 
            *abs_sigma[i+1][0] = *abs_sigma[i][intervals];

            sub_sample = sub_sample*2;

        }

    }

    
    ///////////////////////////////////////////////////////////
    
    Stash ***pos = keypoint_location(differenceofgaussian_pyrmiad, width*2, height*2);

    ///////////////////////////////////////////////////////////


    ///////////////////////////////////////////////////////////

    Stash *orientation_pos = orientation_assignment(gaussain_pyrmiad, pos, abs_sigma, width*2, height*2);

    ///////////////////////////////////////////////////////////


    ///////////////////////////////////////////////////////////

    local_descriptor_calculation(gaussain_pyrmiad, orientation_pos, width*2, height*2);

    ///////////////////////////////////////////////////////////

    SiftFindMatches();

    ///////////////////////////////////////////////////////////
    

    ///////////////////////////////////////////////////////////


    for (i = 0; i <= octaves-1; i++)
    {
        for (int j = 0; j <= intervals + 2; j++)
        {
            free(gaussain_pyrmiad[i][j]);
            
            float iii = *abs_sigma[i][j];
            
            free(abs_sigma[i][j]);
        }
    }


    for (i = 0; i <= octaves - 1; i++) 
    {
        for (int j = 0; j <= intervals + 2; j++)
        {
            cleanup(pos[i][j]);
            free(pos[i][j]);
        }
    }

    for (i = 0; i <= octaves - 1; i++)
    {
        free(pos[i]);
        free(gaussain_pyrmiad[i]);
        free(abs_sigma[i]);
        

    }

    free(pos);
    free(abs_sigma);
    free(gaussain_pyrmiad);

    cleanup(orientation_pos);


    return differenceofgaussian_pyrmiad;

}



 
/*

float *bilinear_zooming(float *image, int width, int height, int new_width, int new_height)
{
    float xScale, yScale, fX, fY;
    
    xScale = (float)width  / (float)new_width;
    yScale = (float)height / (float)new_height;
    
    

    int ifX, ifY, ifX1, ifY1;
    float dx, dy;
    
    int xmax = width  - 1;
    int ymax = height - 1;
    
    float *newimage = (float *)calloc(new_width*new_height, sizeof(float *));
    
    for(int y = 0; y < new_height; y++)
    {
            fY = y * yScale - 0.5f;
            ifY = (int)fY;
            ifY1 = min_bi(ymax, ifY+1);
            dy = fY - ifY;
            
            for(int x = 0; x < new_width; x++)
            {
                    fX = x * xScale - 0.5f;
                    ifX = (int)fX;
                    ifX1 = min_bi(xmax, ifX+1);
                    dx = fX - ifX;
                
                    
                    *(newimage + y*new_width + x) = (*(image + ifY * width + ifX)*(1-dx)*(1-dy) 
                                               + *(image + ifY1 * width + ifX)*dy*(1-dx)
                                               + *(image + ifY * width + ifX1)*dx*(1-dy)
                                               + *(image + ifY1 * width + ifX1)*dx*dy); 
                   
            }  
            
    }
     
    return  newimage;
     
}

*/





Stash *** keypoint_location(float ***differenceofgaussian, int width, int height, 
                                float contrast_threshold, float curvature_threshold)
{

    int i, j, x, y, m, n, k;

    ////////////////////////////////////////////////////
    Stash ***pos;
    
    pos = (Stash ***)calloc(octaves, sizeof(Stash ***));
    
    for (i = 0; i <= octaves - 1; i++)
        pos[i] = (Stash **)calloc(intervals+3, sizeof(Stash **));
       
    for (i = 0; i <= octaves - 1; i++) 
    {
        for (int j = 0; j <= intervals + 2; j++)
        {
            pos[i][j] = (Stash *)malloc(sizeof(Stash));
            init(pos[i][j]);
        }
    } 

    ////////////////////////////////////////////////////

    

    float max, min;

    Stash raw_keypoints, contrast_keypoints, curve_keypoints, final_keypoints; 

    init(&raw_keypoints);
    init(&contrast_keypoints);
    init(&curve_keypoints);
    init(&final_keypoints);

    // The second order derivative kernel to compute the Hessian matrix 
    float xx[3] = { 1, -2, 1 };
    float yy[3] = { 1, -2, 1 };
    float xy[3][3] = { {1/4, 0, -1/4}, {0, 0, 0}, {-1/4, 0, 1/4} };

    float Tr_H = 0;
    float Det_H = 0;
    float curvature_ratio = 0;



    // because the length of the Gaussian filter is 5, the edge is 2 due to the property of 2D convolution. 
    int edge = 2; 

    float extrema_matrix[3][3][3];

    float Dxx = 0, Dyy = 0, Dxy = 0;

    float scale_ratio[4] = {0.5, 1.0, 2.0, 4.0};


    for (i = 0; i <= octaves-1; i++)
    {
        for (j = 1; j<= intervals; j++)
        {

            

            

            for (y = edge+1; y <= height-edge-1; y++)
            {
                for (x = edge+1; x <= width-edge-1; x++) 
                {


                    for (m = 0; m <= 2; m++)
                    {
                        for (n = 0; n <= 2; n++)
                        {
                            extrema_matrix[m][n][0] = *(differenceofgaussian[i][j-1]+  x + (n-1)  +   (y+(m-1))*width);
                            extrema_matrix[m][n][1] = *(differenceofgaussian[i][j]+ x + (n-1)  +   (y+(m-1))*width);
                            extrema_matrix[m][n][2] = *(differenceofgaussian[i][j+1]+ x + (n-1)  +   (y+(m-1))*width);
                        }
                    }

                    extrema_matrix[1][1][1] = extrema_matrix[0][0][0];

                    max = findmax_inmatrix(extrema_matrix) ;
                    min = findmin_inmatrix(extrema_matrix) ;
                
                    

                    if ( (*(differenceofgaussian[i][j]+x+y*width) <1.005*min ) || (*(differenceofgaussian[i][j]+x+y*width) >1.005*max))
                    {

                        
                        add_position(&raw_keypoints, (int)(x*scale_ratio[i]), (int)(y*scale_ratio[i]));

                        if ( fabs(*(differenceofgaussian[i][j]+x+y*width)) >= (contrast_threshold/scale_ratio[i]) )
                        {
                            add_position(&contrast_keypoints, (int)(x*scale_ratio[i]), (int)(y*scale_ratio[i]));

                                
                            for (k = 0; k <= 2; k++)
                            {
                                Dxx += xx[k] * (*(differenceofgaussian[i][j]+x+y*width+k-1));
                
                                Dyy += yy[k] * (*(differenceofgaussian[i][j]+x+y*width+(k-1)*width));
                
                                Dxy += xy[k][0]* (*(differenceofgaussian[i][j]+x+y*width+(k-1)*width-1))
                                            + xy[k][1]* (*(differenceofgaussian[i][j]+x+y*width+(k-1)*width))
                                                + xy[k][2]* (*(differenceofgaussian[i][j]+x+y*width+(k-1)*width+1)); 
                
                            }


                            Tr_H = Dxx + Dyy;
                            Det_H = Dxx*Dyy - Dxy*Dxy;

                            curvature_ratio = (Tr_H*Tr_H)/Det_H;

                            if ((Det_H >= 0) && (curvature_ratio < curvature_threshold))
                            {
                
                                add_position(&curve_keypoints, (int)(x*scale_ratio[i]), (int)(y*scale_ratio[i]));

                                add_position(&final_keypoints, (int)(x*scale_ratio[i]), (int)(y*scale_ratio[i]));

                                /////////////////////////
                    

                                add_position(pos[i][j], x, y);
                                
                                /////////////////////////

                            }

                        }
                    }

                    Dxx = Dyy = Dxy = 0;
                }
            }

        }
    
        width = width/2;
        height = height/2;

    }

    record_position(&final_keypoints);

    cleanup(&raw_keypoints);
    cleanup(&contrast_keypoints);
    cleanup(&curve_keypoints);
    cleanup(&final_keypoints);

    return pos;


}

float findmax_inmatrix(float matrix[3][3][3])
{
    float mx = matrix[0][0][0];

    for(int i = 0; i <= 2; i++)
    {
        if ( matrix[i][0][0] > mx )
            mx = matrix[i][0][0];

        if ( matrix[i][0][1] > mx )
            mx = matrix[i][0][1];
        
        if ( matrix[i][0][2] > mx )
            mx = matrix[i][0][2];
        
        if ( matrix[i][1][0] > mx )
            mx = matrix[i][1][0];
        
        if ( matrix[i][1][1] > mx )
            mx = matrix[i][1][1];
        
        if ( matrix[i][1][2] > mx )
            mx = matrix[i][1][2];
        
        if ( matrix[i][2][0] > mx )
            mx = matrix[i][2][0];
        
        if ( matrix[i][2][1] > mx )
            mx = matrix[i][2][1];
        
        if ( matrix[i][2][2] > mx )
            mx = matrix[i][2][2];
    }

    return mx;
        


}

float findmin_inmatrix(float matrix[3][3][3])
{
    float mi =  matrix[0][0][0];

    for(int i = 0; i <= 2; i++)
    {
        if ( mi > matrix[i][0][0] )
            mi = matrix[i][0][0];

        if ( mi > matrix[i][0][1] )
            mi = matrix[i][0][1];
        
        if ( mi > matrix[i][0][2] )
            mi = matrix[i][0][2];
        
        if ( mi > matrix[i][1][0] )
            mi = matrix[i][1][0];
        
        if ( mi > matrix[i][1][1] )
            mi = matrix[i][1][1];
        
        if ( mi > matrix[i][1][2] )
            mi = matrix[i][1][2];
        
        if ( mi > matrix[i][2][0] )
            mi = matrix[i][2][0];
        
        if ( mi > matrix[i][2][1] )
            mi = matrix[i][2][1];
        
        if ( mi > matrix[i][2][2] )
            mi = matrix[i][2][2];
    }

    return mi;


}



#define PI 3.14159265


Stash * orientation_assignment(float ***gaussain_pyrmiad, Stash ***pos, float ***abs_sigma, int width, int height)
{
    int i,j;

    float sub_sample[] = {0.5, 1.0, 2.0, 3.0 };


    /* define the data structure to store all the orientation information */

    Stash *orientation_keypoints = (Stash *)malloc(sizeof(Stash));

    init(orientation_keypoints);

    int zero_padding = 2;

      float *magnitude_pyrmiad[octaves][intervals+3];
    
    float *gradient_pyrmiad[octaves][intervals+3];
    
    float diff_x, diff_y;

    int ori_width = width; int ori_height = height;
    
    for ( i = 0; i <= octaves-1; i++)
    {
        for ( j = 1; j<= intervals; j++)
        {

        
            /*   */
            magnitude_pyrmiad[i][j] = (float *)calloc((width+2*zero_padding)*(height+2*zero_padding), 
                                                                                        sizeof(float *));
            gradient_pyrmiad[i][j] = (float *)calloc((width+2*zero_padding)*(height+2*zero_padding), 
                                                                                        sizeof(float *));

            for (int jj = 1; jj <= height-2; jj++)
            {
                for(int ii = 1; ii <= width-2; ii++)
                {

                    /*  compute the first order derivative to get the gradient and magnitude */

                    diff_x = *(gaussain_pyrmiad[i][j] + jj*width + ii-1)
                                    - *(gaussain_pyrmiad[i][j] + jj*width + ii+1);

                    diff_y = *(gaussain_pyrmiad[i][j] + (jj-1)*width + ii)
                                    - *(gaussain_pyrmiad[i][j] + (jj+1)*width + ii);

                    diff_x = diff_x/2;

                    diff_y = diff_y/2;


                    *(magnitude_pyrmiad[i][j] + (jj+zero_padding)*(width+2*zero_padding) + ii + zero_padding)
                                    = (float)sqrt(diff_x*diff_x + diff_y*diff_y);

                    float temp  = (float)atan2(diff_y, diff_x);
                    float temp1 = (float)temp*180/PI;
                    
                    *(gradient_pyrmiad[i][j] + (jj+zero_padding)*(width+2*zero_padding) + ii + zero_padding)
                                    = (temp1 == 180)?(-temp):temp;

                }

            }

        }

        width  = width/2;
        height = height/2;

    }    

    width = ori_width; height = ori_height;
    

    /* Generate the orientation histogram bins, 36 bins from [-PI, PI)  */
    const int num_bins = 36;
    float hist_step = 2*PI/num_bins;
    float hist_orient[num_bins];
    
    for(i = 0; i <= num_bins-1; i++)
    {

        hist_orient[i] = -PI + i*hist_step;
        
    }



    
    float **gaussianKernel;
    float weight[5][5];
    float diff[5][5];
    float temp[5];
      
    for ( i = 0; i <= octaves-1; i++)
    {
        for ( j = 1; j<= intervals; j++)
        {

            gaussianKernel = gaussian_coeff(1.5*(*(abs_sigma[i][j]))/sub_sample[i], 5);
            
            int num_keypoints = count_position(pos[i][j])/2;

            for (int ii = 0; ii <= num_keypoints-1; ii++)
            {
                int x = pos[i][j]->position[2*ii] ;

                int y = pos[i][j]->position[2*ii+1];

                //x += zero_padding;

            
                //y += zero_padding;

                float *orient_hist = (float *)calloc(num_bins, sizeof(float));

                for (int jj = 0; jj < 5; jj++)
                {
                    weight[jj][0] = (*(magnitude_pyrmiad[i][j] + (y+jj)*(width+2*zero_padding) + x    ))*gaussianKernel[jj][0];
                    weight[jj][1] = (*(magnitude_pyrmiad[i][j] + (y+jj)*(width+2*zero_padding) + x + 1))*gaussianKernel[jj][1];
                    weight[jj][2] = (*(magnitude_pyrmiad[i][j] + (y+jj)*(width+2*zero_padding) + x + 2))*gaussianKernel[jj][2];
                    weight[jj][3] = (*(magnitude_pyrmiad[i][j] + (y+jj)*(width+2*zero_padding) + x + 3))*gaussianKernel[jj][3];
                    weight[jj][4] = (*(magnitude_pyrmiad[i][j] + (y+jj)*(width+2*zero_padding) + x + 4))*gaussianKernel[jj][4];


                    diff[jj][0] = *(gradient_pyrmiad[i][j] + (y+jj)*(width+2*zero_padding) + x     );
                    diff[jj][1] = *(gradient_pyrmiad[i][j] + (y+jj)*(width+2*zero_padding) + x + 1 );
                    diff[jj][2] = *(gradient_pyrmiad[i][j] + (y+jj)*(width+2*zero_padding) + x + 2 );
                    diff[jj][3] = *(gradient_pyrmiad[i][j] + (y+jj)*(width+2*zero_padding) + x + 3 );
                    diff[jj][4] = *(gradient_pyrmiad[i][j] + (y+jj)*(width+2*zero_padding) + x + 4 );
                    
                }

                for (int bin = 0; bin <= num_bins-1; bin++)
                {
                    for (int jj = 0; jj < 5; jj++)
                    {
                        temp[0] = fmod((diff[jj][0]-hist_orient[bin] + PI), 2*PI) - PI;
                        temp[1] = fmod((diff[jj][1]-hist_orient[bin] + PI), 2*PI) - PI;
                        temp[2] = fmod((diff[jj][2]-hist_orient[bin] + PI), 2*PI) - PI;
                        temp[3] = fmod((diff[jj][3]-hist_orient[bin] + PI), 2*PI) - PI;
                        temp[4] = fmod((diff[jj][4]-hist_orient[bin] + PI), 2*PI) - PI;
/*
                        temp[0] = (fabs(temp[0]) <= hist_step)?1:0;
                        temp[1] = (fabs(temp[1]) <= hist_step)?1:0;
                        temp[2] = (fabs(temp[2]) <= hist_step)?1:0;
                        temp[3] = (fabs(temp[3]) <= hist_step)?1:0;
                        temp[4] = (fabs(temp[4]) <= hist_step)?1:0;
*/
                        temp[0] = max((1- fabs(temp[0])/hist_step), 0);
                        temp[1] = max((1- fabs(temp[1])/hist_step), 0);
                        temp[2] = max((1- fabs(temp[2])/hist_step), 0);
                        temp[3] = max((1- fabs(temp[3])/hist_step), 0);
                        temp[4] = max((1- fabs(temp[4])/hist_step), 0);


                        temp[0] = temp[0]*weight[jj][0];
                        temp[1] = temp[1]*weight[jj][1];
                        temp[2] = temp[2]*weight[jj][2];
                        temp[3] = temp[3]*weight[jj][3];
                        temp[4] = temp[4]*weight[jj][4];

                        orient_hist[bin] += temp[0] + temp[1] + temp[2] + temp[3] + temp[4];

            
                       }

                    
                    
                }

                /*  suppressing the nonmaxima  */

                float *peaks = (float *)malloc(num_bins*sizeof(float *));

                for ( bin = 1; bin <= num_bins-2; bin++)
                {
                    if ((orient_hist[bin] < orient_hist[bin-1]) || (orient_hist[bin] < orient_hist[bin+1]))
                    {
                        peaks[bin] = 0;
                    }
                    else
                    {
                        peaks[bin] = orient_hist[bin];

                    }
                }

                if ((orient_hist[0] < orient_hist[1]) || (orient_hist[0] < orient_hist[num_bins-1]))
                    peaks[0] = 0;
                else
                    peaks[0] = orient_hist[0];


                if ((orient_hist[num_bins-1] < orient_hist[num_bins-2]) || (orient_hist[num_bins-1] < orient_hist[0]))
                    peaks[num_bins-1] = 0;
                else
                    peaks[num_bins-1] = orient_hist[num_bins-1];

                int peak_pos = finding_peak(peaks, num_bins);

                float max = peaks[peak_pos];

                float peak_value = max;

                
                while (peak_value > 0.8*max)
                {


                    float a = hist_orient[peak_pos]-hist_step;
                    float b = hist_orient[peak_pos];
                    float c = hist_orient[peak_pos]+hist_step;


                    bin = fmod( peak_pos -1 + num_bins - 1, num_bins ) + 1;
                    float fa = orient_hist[bin];

                    bin = fmod( peak_pos + num_bins - 1, num_bins ) + 1;
                    float fb = orient_hist[bin];

                    bin = fmod( peak_pos + 1 + num_bins - 1, num_bins ) + 1;
                    float fc = orient_hist[bin];
    
                    float max_orient = parabolic_interp(a, fa, b, fb, c, fc);
                    
                    while( max_orient < -PI )
                    {
                        max_orient = max_orient + 2*PI;
                    }

                    while( max_orient >= PI )
                    {
                        max_orient = max_orient - 2*PI;
                    }  

                    //float max_orient = hist_orient[peak_pos];

                    //add_orientation(&orientation_keypoints, (x-zero_padding), (y-zero_padding), (max_orient*180/PI));
                    add_orientation(orientation_keypoints, x, y, (max_orient*180/PI));

                    add_orientation(orientation_keypoints, i, j, *abs_sigma[i][j]);

                    peaks[peak_pos] = 0;
                    peak_pos = finding_peak(peaks, num_bins);
                    peak_value = peaks[peak_pos];

                }

                free(peaks);
                free(orient_hist);
                

            }

            for (int iii = 0; iii <= 4; iii++)
                free(gaussianKernel[iii]);
            free(gaussianKernel);


        }

        width = width/2;
        height = height/2;
    
    }

    record_orientation(orientation_keypoints);

    return orientation_keypoints;

    //cleanup(&orientation_keypoints);

}

int finding_peak(float *peaks, int num_bins)
{
    float max = peaks[0];

    int pos = 0;


    for(int i  = 1; i <= num_bins-1; i++)
    {
        if (peaks[i] >= max)
        {
            pos = i;
            max = peaks[i];

        }

    }

    return pos;
    
}


// parabolic interpolation
float parabolic_interp(float a, float fa, float b, float fb, float c, float fc)
{
       
    float p = (b - a) * (fb - fc);
    float q = (b - c) * (fb - fa);
    float x = (b - c) * q - (b - a) * p;
    float y = (float)2.0 * (p - q);
       
    float result = b + x/y;
       
    return result;
}




void local_descriptor_calculation(float ***gaussain_pyrmiad, Stash *orientation_pos, int width, int height)
{
    typedef struct ftag{  
           float x;
           float y;
         
           }fCoords; 


    int i,j;

    int ii, jj;

    float orient_bin_step = PI/4.0;
    
    float orient_angles[8] = {-PI, -PI+orient_bin_step, -PI+2*orient_bin_step, -PI+3*orient_bin_step, 0, 
                                   orient_bin_step, orient_bin_step*2, orient_bin_step*3};
                                       
            
    float center_cor[4] = {-6.0, -2.0, 2.0, 6.0};       
    
    int grid_space = 4;
    
    int feature_weight_window = 8;
            
    fCoords center_coordinations[4][4];
    
    
/*    
    for (i = 0; i <=3; i++)
    {
        
        for (j = 0; j <=3; j++)
        {
            center_coordinations[i][j].x = center_cor[i];
            center_coordinations[i][j].y = center_cor[j];
       
        }        
    }
*/    
    float sample_cor[16] = {-7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5};
    
    fCoords sample_coordinations[16][16];
/*    
    for (i = 0; i <=15; i++)
    {
        for (j = 0; j <=15; j++)
        {
            sample_coordinations[i][j].x = sample_cor[i];
            sample_coordinations[i][j].y = sample_cor[j];
          
        }
    }
*/

    int i_cor, j_cor;

    int inn_width;

    int inn_height;


    //////////////////////////////////////////////////////////////

    int num_samples = 256;
       
    int num_feature_point = count_orientation(orientation_pos)/2;



    
    FILE * pFile;
    pFile = fopen ("local_descriptor.key","w");

    
    fprintf(pFile, "%i", num_feature_point);

    fprintf(pFile, "\r\n"); 
    
    for(i = 0; i <= num_feature_point-1; i ++)
    {
        float feature_descriptor[128] = {0.0f}; 
        float x_weight[4][4], y_weight[4][4];
        float position_weight[4][4];
        float orientation_weight[8]; 
        float gaussian_weight;


        /*  For each key point, the local descriptor is calculated from a 16*16 block centered at the key point position  */
        /*  First the 16*16 block will rotated to align with the orientation of the key point to achieve the rotation invariance */
          
        ////////////////////////////////////////////////////////////////
        
        
        int x = orientation_pos->position[4*i];
        
        int y = orientation_pos->position[4*i+1];
                
        float feature_orentation = orientation_pos->orientation[2*i]*PI/180;

        int octave = orientation_pos->position[4*i+2];
    
        int interval = orientation_pos->position[4*i+3];

        float scale = orientation_pos->orientation[2*i+1];

    

        /////////////////////////////////////////////////////////////////


        inn_width = width/pow(2, octave);
        inn_height = height/pow(2, octave);

        for (i_cor = 0; i_cor <=3; i_cor++)
        {
        
            for (j_cor = 0; j_cor <=3; j_cor++)
            {
                center_coordinations[i_cor][j_cor].x = center_cor[i_cor];
                center_coordinations[i_cor][j_cor].y = center_cor[j_cor];
       
            }        
        }

        for (i_cor = 0; i_cor <=15; i_cor++)
        {
            for (j_cor = 0; j_cor <=15; j_cor++)
            {
                sample_coordinations[i_cor][j_cor].x = sample_cor[i_cor];
                sample_coordinations[i_cor][j_cor].y = sample_cor[j_cor];
          
            }
        }

    
        float rotation_matrix[2][2];
        rotation_matrix[0][0] = cos(feature_orentation); rotation_matrix[0][1] = -sin(feature_orentation); 
        rotation_matrix[1][0] = sin(feature_orentation); rotation_matrix[1][1] =  cos(feature_orentation);
        
          
        float temp1, temp2;
        
        int iii, jjj;
    
        for (iii = 0; iii <=3; iii++)
        {
            for (jjj = 0; jjj <=3; jjj++)
            {
                temp1 = rotation_matrix[0][0]*center_coordinations[iii][jjj].x 
                                         + rotation_matrix[0][1]*center_coordinations[iii][jjj].y;

                
            
                temp2 = rotation_matrix[1][0]*center_coordinations[iii][jjj].x 
                                         + rotation_matrix[1][1]*center_coordinations[iii][jjj].y;

                
                                         
                center_coordinations[iii][jjj].x = temp1;
            
                center_coordinations[iii][jjj].y = temp2;
             }
        }
    
        for (iii = 0; iii <=15; iii++)
        {
            for (jjj = 0; jjj <=15; jjj++)
            {
                temp1 = rotation_matrix[0][0]*sample_coordinations[iii][jjj].x 
                                         + rotation_matrix[0][1]*sample_coordinations[iii][jjj].y;
                
            
                temp2 = rotation_matrix[1][0]*sample_coordinations[iii][jjj].x 
                                         + rotation_matrix[1][1]*sample_coordinations[iii][jjj].y;
            
                                         
                sample_coordinations[iii][jjj].x = temp1;
            
                sample_coordinations[iii][jjj].y = temp2;
            }
        }
    
        for (iii = 0; iii <=3; iii++)
        {
            for (jjj = 0; jjj <=3; jjj++)
            {

                temp1 = center_coordinations[iii][jjj].x+x;

                if (temp1 < 0) temp1 = 0;
                if (temp1 > inn_width-1) temp1 = inn_width-1;

                center_coordinations[iii][jjj].x = temp1; 

                temp2 = center_coordinations[iii][jjj].y+y;

                if (temp2 < 0) temp2 = 0;
                if (temp2 > inn_height-1) temp2 = inn_height-1;
                                
                center_coordinations[iii][jjj].y = temp2;
                
            }
        }
  
    
        for(iii = 0; iii <= 15; iii++)
        {        
            for(jjj = 0; jjj <= 15; jjj++ )
            {
                temp1 = sample_coordinations[iii][jjj].x+x;

                if (temp1 < 0) temp1 = 0;
                if (temp1 > inn_width-1) temp1 = inn_width-1;

                sample_coordinations[iii][jjj].x = temp1; 
                
                temp2 = sample_coordinations[iii][jjj].y+y;

                if (temp2 < 0) temp2 = 0;
                if (temp2 > inn_height-1) temp2 = inn_height-1;
                                
                sample_coordinations[iii][jjj].y = temp2;
                            
                float samples[3][3];


                
                float x_temp;// = sample_coordinations[iii][jjj].x;
                float y_temp;// = sample_coordinations[iii][jjj].y-1;
                float xy_temp1, xy_temp2;

                x_temp = -1;
                y_temp = 0;

                xy_temp1 = rotation_matrix[0][0]*x_temp 
                                         + rotation_matrix[0][1]*y_temp;
                        
            
                xy_temp2 = rotation_matrix[1][0]*x_temp
                                         + rotation_matrix[1][1]*y_temp;

                x_temp = xy_temp1 + sample_coordinations[iii][jjj].x;
                y_temp = xy_temp2 + sample_coordinations[iii][jjj].y;

                samples[1][0] = interp2(gaussain_pyrmiad[octave][interval], inn_width, inn_height, x_temp, y_temp);
              
                //samples[1][0] = interp2(gaussain_pyrmiad[octave][interval], inn_width, inn_height, sample_coordinations[iii][jjj].x, sample_coordinations[iii][jjj].y-1);
              

                x_temp = 1;
                y_temp = 0;

                xy_temp1 = rotation_matrix[0][0]*x_temp 
                                         + rotation_matrix[0][1]*y_temp;
                            
            
                xy_temp2 = rotation_matrix[1][0]*x_temp
                                         + rotation_matrix[1][1]*y_temp;

                x_temp = xy_temp1 + sample_coordinations[iii][jjj].x;
                y_temp = xy_temp2 + sample_coordinations[iii][jjj].y;


                samples[1][2] = interp2(gaussain_pyrmiad[octave][interval], inn_width, inn_height,x_temp, y_temp);
            

                //samples[1][2] = interp2(gaussain_pyrmiad[octave][interval], inn_width, inn_height, sample_coordinations[iii][jjj].x, sample_coordinations[iii][jjj].y+1);
              
                x_temp = 0;
                y_temp = 1;

                xy_temp1 = rotation_matrix[0][0]*x_temp 
                                         + rotation_matrix[0][1]*y_temp;
                            
            
                xy_temp2 = rotation_matrix[1][0]*x_temp
                                         + rotation_matrix[1][1]*y_temp;

                x_temp = xy_temp1 + sample_coordinations[iii][jjj].x;
                y_temp = xy_temp2 + sample_coordinations[iii][jjj].y;

                samples[2][1] = interp2(gaussain_pyrmiad[octave][interval], inn_width, inn_height, x_temp, y_temp);

                //samples[2][1] = interp2(gaussain_pyrmiad[octave][interval], inn_width, inn_height, sample_coordinations[iii][jjj].x+1, sample_coordinations[iii][jjj].y);
              
                x_temp = 0;
                y_temp = -1;

                xy_temp1 = rotation_matrix[0][0]*x_temp 
                                         + rotation_matrix[0][1]*y_temp;
            
                            
                xy_temp2 = rotation_matrix[1][0]*x_temp
                                         + rotation_matrix[1][1]*y_temp;

                x_temp = xy_temp1 + sample_coordinations[iii][jjj].x;
                y_temp = xy_temp2 + sample_coordinations[iii][jjj].y;

                samples[0][1] = interp2(gaussain_pyrmiad[octave][interval], inn_width, inn_height, x_temp, y_temp);

                //samples[0][1] = interp2(gaussain_pyrmiad[octave][interval], inn_width, inn_height, sample_coordinations[iii][jjj].x-1, sample_coordinations[iii][jjj].y);
              
                float diff_x = 0.5*(samples[1][2] - samples[1][0]);
                float diff_y = 0.5*(samples[2][1] - samples[0][1]);
              
                float mag_sample = 1;//sqrt(diff_x*diff_x + diff_y*diff_y);
                float grad_sample = (float)atan2(diff_y, diff_x);
           
                temp1 = (float)grad_sample*180/PI;
                    
                grad_sample = (temp1 == 180)?(-grad_sample):grad_sample;
                                  
                for (ii = 0; ii <= 3; ii++)
                {
                    for (jj = 0; jj <= 3; jj++)
                    {
                        x_weight[ii][jj] = max(1-fabs( sample_coordinations[iii][jjj].x - center_coordinations[ii][jj].x)/grid_space, 0);
                        y_weight[ii][jj] = max(1-fabs( sample_coordinations[iii][jjj].y - center_coordinations[ii][jj].y)/grid_space, 0);
                        position_weight[ii][jj] = x_weight[ii][jj]*y_weight[ii][jj]; 
                    }
                }
              
                for(ii = 0; ii <= 7; ii++)
                {
                    //float difference = fmod(grad_sample - feature_orentation - orient_angles[ii] + PI, 2*PI) - PI;
                    float difference = fmod(grad_sample - orient_angles[ii] + PI, 2*PI) - PI;
                    orientation_weight[ii] = max(1 - fabs(difference)/orient_bin_step, 0);
                   
                }
              
                gaussian_weight = exp(-(pow((sample_coordinations[iii][jjj].x - x),2) + pow((sample_coordinations[iii][jjj].y-y), 2))/pow(2*feature_weight_window, 2))/pow(2*PI*feature_weight_window, 2);
              
                for(ii = 0; ii <= 7; ii++)
                {
                    feature_descriptor[ii + 8*0]  += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[0][0];
                    feature_descriptor[ii + 8*1]  += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[0][1];
                    feature_descriptor[ii + 8*2]  += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[0][2];
                    feature_descriptor[ii + 8*3]  += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[0][3];   
                    feature_descriptor[ii + 8*4]  += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[1][0];
                    feature_descriptor[ii + 8*5]  += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[1][1];
                    feature_descriptor[ii + 8*6]  += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[1][2];
                    feature_descriptor[ii + 8*7]  += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[1][3];
                    feature_descriptor[ii + 8*8]  += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[2][0];
                    feature_descriptor[ii + 8*9]  += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[2][1];
                    feature_descriptor[ii + 8*10] += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[2][2];
                    feature_descriptor[ii + 8*11] += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[2][3];
                    feature_descriptor[ii + 8*12] += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[3][0];
                    feature_descriptor[ii + 8*13] += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[3][1];
                    feature_descriptor[ii + 8*14] += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[3][2];
                    feature_descriptor[ii + 8*15] += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[3][3];
                  
                }
            }
        }       
         
      
      
    
        float feature_descriptor_norm = 0.0f;
      
        for(iii = 0; iii <= 127; iii++)
        {
            feature_descriptor_norm += feature_descriptor[iii]*feature_descriptor[iii];
        } 
        feature_descriptor_norm = sqrt(feature_descriptor_norm);
        for(iii = 0; iii <= 127; iii++)
        {
            feature_descriptor[iii] = feature_descriptor[iii]/feature_descriptor_norm;
            if (feature_descriptor[iii] > 0.2)
                feature_descriptor[iii] = 0.2;       
        } 
        feature_descriptor_norm = 0;
        for(iii = 0; iii <= 127; iii++)
        {
            feature_descriptor_norm += feature_descriptor[iii]*feature_descriptor[iii];
        }     
        for(iii = 0; iii <= 127; iii++)
        {
            feature_descriptor[iii] = feature_descriptor[iii]/feature_descriptor_norm;
        }          
    
    
    

        int temp111 = x/pow(2, 1-octave);
        int temp222 = y/pow(2, 1-octave);

        //fprintf(pFile, "%i %i %f %f", x/pow(2, 1-octave), y/pow(2, 1-octave), feature_orentation, scale);
        fprintf(pFile, "%i %i %f %f", temp111, temp222, feature_orentation, scale);
        fprintf(pFile, "\r\n"); 
        for(iii = 0; iii <= 127; iii++)
        {
            fprintf(pFile, "%f ", feature_descriptor[iii]);
        } 
        fprintf(pFile, "\r\n"); 
    }

    fclose (pFile);
    /////////////////////////





}



float interp2(float *image, int width, int height, float fx, float fy)
{
    int ix = floor(fx);
    int iy = floor(fy);
    
    if ((ix < 0)||(iy < 0)||(ix >= width-1) || (iy >= height-1))
        return 0.0f;
       
    int ix1 = ix + 1;
    int iy1 = iy + 1;
    
    float dx = fx - ix;
    float dy = fy - iy;
    
    float result = (*(image+iy*width+ix))*(1-dx)*(1-dy) 
                       +   (*(image+iy*width+ix1))*dx*(1-dy) 
                           + (*(image+iy1*width+ix))*(1-dx)*dy 
                               + (*(image+iy1*width+ix1))*dx*dy;   
                               
    return result;
      
}


//////////////////////////////////////////////////////////////////////

Keypoint ReadKeys(FILE *fp)
{
    int i, j, num_keypoints; 
    Keypoint k, keys = NULL;

    int length_descriptor = 128;

    float val;
    
    fscanf(fp, "%i", &num_keypoints);


    for (i = 0; i < num_keypoints; i++) 
    {
      // Allocate memory for the keypoint. 
      k = (Keypoint) malloc(sizeof(struct KeypointSt));
      k->next = keys;
      keys = k;
      k->feature_descriptor = (float *)malloc(length_descriptor * sizeof(float));


      if (fscanf(fp, "%i %i %f %f", &(k->x), &(k->y), &(k->feature_orentation),     &(k->scale)) != 4)
          FatalError("Invalid keypoint file format.");

      for (j = 0; j < length_descriptor; j++) 
      {
          if (fscanf(fp, "%f", &val) != 1 || val < 0 || val > 255)
             FatalError("Invalid keypoint file value.");
          k->feature_descriptor[j] = val;
      }
    }

    return keys;
}

int MatchCalculation(float *keypoint, Keypoint klist)
{

    int i,j;

    const int length_descriptor = 128;

    float feature_descriptor[length_descriptor];

    for (j = 0; j < length_descriptor; j++) 
    {

        feature_descriptor[j] = keypoint[j];

    }

    
    float min, second_to_min;

    Keypoint head = klist;

    int num = 0;

    if (klist != NULL)
    {
        do    {
            num++;
            klist = klist->next;
        }while (klist != NULL);
    }


    if ( num >0 )
    {
        float *result = (float *)calloc(num, sizeof(float));

        

        for (i = 0; i < num; i++)
        {
        

            for (j = 0; j < length_descriptor; j++) 
            {
                result[i] += pow((feature_descriptor[j] - head->feature_descriptor[j]), 2);
            }

            head = head->next;

        }

        int position = 0;

        min = result[0];


        for (i = 0; i < num; i++)
        {

            if (result[i] < min)
            {
                position = i;
                min = result[i];
            }

        }

        result[position] = 100000;

        second_to_min = result[0];

        for (i = 0; i < num; i++)
        {

            if (result[i] < second_to_min)
            {
                
                second_to_min = result[i];
            }

        }

        if ( min < 0.36*second_to_min )
        {
            return position;
        }


    }

    return -1;
}


void FatalError(char *fmt, ...)
{
    va_list args;

    va_start(args, fmt);
    fprintf(stderr, "Error: ");
    vfprintf(stderr, fmt, args);
    fprintf(stderr,"\n");
    va_end(args);
    exit(1);
}




void SiftDrawLine(BYTE *image, int width, int height, int x1, int y1, int x2, int y2)
{
    int i, dx, dy, temp, x_incr, y_incr;

    if (x1 == x2 && y1 == y2)  // Line of zero length. 
      return;

    // Is line more horizontal than vertical? 
    if (abs(x2 - x1) < abs(y2 - y1)) 
    {

        // Put points in increasing order by column. 
        if (y1 > y2) 
        {
            temp = x1; x1 = x2; x2 = temp;
            temp = y1; y1 = y2; y2 = temp;
        }
        dx = x2 - x1;
        dy = y2 - y1;
        for (i = y1; i <= y2; i++)
        {
            x_incr = x1 + (i - y1) * dx / dy;
            
            *(image + i*width + x_incr) = 255;
        }

    } 
    else 
    {

        if (x1 > x2) 
        {
            temp = x1; x1 = x2; x2 = temp;
            temp = y1; y1 = y2; y2 = temp;
        }
        dx = x2 - x1;
        dy = y2 - y1;
        for (i = x1; i <= x2; i++)
        {
            y_incr = y1 + (i - x1) * dy / dx;
            *(image + y_incr*width + i) = 255;
        }
    }
}



void SiftFindMatches()
{
    FILE * pMainFile;
    pMainFile = fopen ("local_descriptor.key","r");

    FILE * pMatchFile;
    pMatchFile = fopen("match_result.key","w");

    int position;

    Keypoint k_main = NULL;
    Keypoint k_loc = NULL;
    Keypoint k_template = NULL;

    FILE * pTemplateFile;
    pTemplateFile = fopen ("template_descriptor.key","r");


    if ( (pMainFile != NULL) && (pTemplateFile != NULL) )
    {
        k_main = ReadKeys(pMainFile);
    
        k_loc = k_main;

        k_template = ReadKeys(pTemplateFile);
    
    

        while( (k_template != NULL) && (k_main != NULL))
        {
            position = MatchCalculation(k_template->feature_descriptor, k_main);

            if (position != -1)
            {
                while (position != 0)
                {
                    k_loc = k_loc->next;
                    position --;
                }

                fprintf(pMatchFile, "%i %i", k_loc->x, k_loc->y);

                fprintf(pMatchFile, "\r\n"); 
            
            }

            k_loc = k_main;

            k_template = k_template->next;

        }

    }
    if (pMainFile != NULL)
    {
        fclose(pMainFile);
    }
    if (pTemplateFile != NULL)
    {
        fclose(pTemplateFile);
    }    

    fclose(pMatchFile);

}


/////////////////////////////////////////////////////////////////////////
//               Test
/////////////////////////////////////////////////////////////////////////
float orientation_calc_test(BYTE *image, int width, int height)
{

    unsigned int *image_i = byte2int(image, width, height);



    float *nor_image = normalize(image_i, width, height); 

    float *filtered_image = filtering(nor_image, width, height, 3.5, 5);

    free(nor_image);


    free(image_i);


    const int zero_padding = 2;

      float *magnitude_pyrmiad;
    
    float *gradient_pyrmiad;
    
    float diff_x, diff_y, diff_x_1, diff_x_2, diff_x_3, diff_y_1, diff_y_2, diff_y_3;
    
    magnitude_pyrmiad = (float *)calloc((width+2*zero_padding)*(height+2*zero_padding), 
                                                                                        sizeof(float *));
    gradient_pyrmiad = (float *)calloc((width+2*zero_padding)*(height+2*zero_padding), 
                                                                                        sizeof(float *));
    
    
    for (int jj = 1; jj <= height-2; jj++)
    {
        for(int ii = 1; ii <= width-2; ii++)
        {
            
            
            diff_x_1 = *(filtered_image + (jj-1)*width + ii+1)
                                    - *(filtered_image + (jj-1)*width + ii-1);

            diff_x_2 = *(filtered_image + jj*width + ii+1)
                                    - *(filtered_image + jj*width + ii-1);

            diff_x_3 = *(filtered_image + (jj+1)*width + ii+1)
                                    - *(filtered_image + (jj+1)*width + ii-1);

            diff_y_1 = *(filtered_image + (jj+1)*width + ii-1)
                                    - *(filtered_image + (jj-1)*width + ii-1);

            diff_y_2 = *(filtered_image + (jj+1)*width + ii)
                                    - *(filtered_image + (jj-1)*width + ii);

            diff_y_3 = *(filtered_image + (jj+1)*width + ii+1)
                                    - *(filtered_image + (jj-1)*width + ii+1);

            diff_x = (0*diff_x_1 + 1*diff_x_2 + 0*diff_x_3)/2;

            diff_y = (0*diff_y_1 + 1*diff_y_2 + 0*diff_y_3)/2;

                    

            *(magnitude_pyrmiad + (jj+zero_padding)*(width+2*zero_padding) + ii + zero_padding)
                                    = (float)sqrt(diff_x*diff_x + diff_y*diff_y);

            float temp  = (float)atan2(diff_y, diff_x);
            float temp1 = (float)temp*180/PI;
                    
            *(gradient_pyrmiad + (jj+zero_padding)*(width+2*zero_padding) + ii + zero_padding) 
                                    = (temp1 == 180)?(-temp):temp;

        }

    }

//    diff_x = *(filtered_image + 256*width + 257)
//                                    - *(filtered_image + 256*width + 255);

//    diff_y = *(filtered_image + 257*width + 256)
//                                    - *(filtered_image + 255*width + 256);

//    float test = atan2(diff_y, diff_x);

    
    const int num_bins = 360;
    float hist_step = 2*PI/num_bins;
    float hist_orient[num_bins];
    
    for(int i = 0; i <= num_bins-1; i++)
    {

        hist_orient[i] = -PI + i*hist_step;
        
    }
    
    float **gaussianKernel;
    float weight[9][9]={0};
    float diff[9][9]={0};
    float temp[9]={0};
      

    gaussianKernel = gaussian_coeff(5.25, 9);

    float norm = 1;//gaussianKernel[4][4]; 

    //for (i = 0; i <= 8; i++)
        //for (int j = 0; j <= 8; j++)
            //norm = (gaussianKernel[i][j]>norm)?gaussianKernel[i][j]:norm;
            //norm = norm + gaussianKernel[i][j];
            
    int x = 300;

    int y = 300;

    float *orient_hist = (float *)calloc(num_bins, sizeof(float));

    //r = 5;


    for ( jj = 0; jj < 9; jj++)
    //for ( jj = 0; jj < 9; jj++)
    {
        weight[jj][0] = (*(magnitude_pyrmiad + (y+jj-2)*(width+2*zero_padding) + x - 2))*gaussianKernel[jj][0]/norm;
        weight[jj][1] = (*(magnitude_pyrmiad + (y+jj-2)*(width+2*zero_padding) + x - 1))*gaussianKernel[jj][1]/norm;
        weight[jj][2] = (*(magnitude_pyrmiad + (y+jj-2)*(width+2*zero_padding) + x + 0))*gaussianKernel[jj][2]/norm;
        weight[jj][3] = (*(magnitude_pyrmiad + (y+jj-2)*(width+2*zero_padding) + x + 1))*gaussianKernel[jj][3]/norm;
        weight[jj][4] = (*(magnitude_pyrmiad + (y+jj-2)*(width+2*zero_padding) + x + 2))*gaussianKernel[jj][4]/norm;
        weight[jj][5] = (*(magnitude_pyrmiad + (y+jj-2)*(width+2*zero_padding) + x + 3))*gaussianKernel[jj][5]/norm;
        weight[jj][6] = (*(magnitude_pyrmiad + (y+jj-2)*(width+2*zero_padding) + x + 4))*gaussianKernel[jj][6]/norm;
        weight[jj][7] = (*(magnitude_pyrmiad + (y+jj-2)*(width+2*zero_padding) + x + 5))*gaussianKernel[jj][7]/norm;
        weight[jj][8] = (*(magnitude_pyrmiad + (y+jj-2)*(width+2*zero_padding) + x + 6))*gaussianKernel[jj][8]/norm;
    

        diff[jj][0] = *(gradient_pyrmiad + (y+jj-2)*(width+2*zero_padding) + x - 2 );
        diff[jj][1] = *(gradient_pyrmiad + (y+jj-2)*(width+2*zero_padding) + x - 1 );
        diff[jj][2] = *(gradient_pyrmiad + (y+jj-2)*(width+2*zero_padding) + x + 0 );
        diff[jj][3] = *(gradient_pyrmiad + (y+jj-2)*(width+2*zero_padding) + x + 1 );
        diff[jj][4] = *(gradient_pyrmiad + (y+jj-2)*(width+2*zero_padding) + x + 2 );
        diff[jj][5] = *(gradient_pyrmiad + (y+jj-1)*(width+2*zero_padding) + x + 3 );
        diff[jj][6] = *(gradient_pyrmiad + (y+jj-2)*(width+2*zero_padding) + x + 4 );
        diff[jj][7] = *(gradient_pyrmiad + (y+jj-2)*(width+2*zero_padding) + x + 5 );
        diff[jj][8] = *(gradient_pyrmiad + (y+jj-2)*(width+2*zero_padding) + x + 6 );
        
                    
    }

    weight[0][0] = weight[0][1] = weight[1][0] = weight[1][1] = 0;
    weight[8][8] = weight[8][7] = weight[7][8] = weight[7][7] = 0;
    weight[0][8] = weight[0][7] = weight[1][8] = weight[1][7] = 0;
    weight[8][0] = weight[7][0] = weight[8][1] = weight[7][1] = 0;

    for (int bin = 0; bin <= num_bins-1; bin++)
    {
        for (int jj = 0; jj < 9; jj++)
        {
            float temporary;

            temp[0] = fmod((diff[jj][0]-hist_orient[bin] + PI), 2*PI) - PI;
            temp[1] = fmod((diff[jj][1]-hist_orient[bin] + PI), 2*PI) - PI;
            temp[2] = fmod((diff[jj][2]-hist_orient[bin] + PI), 2*PI) - PI;
            temp[3] = fmod((diff[jj][3]-hist_orient[bin] + PI), 2*PI) - PI;
            temp[4] = fmod((diff[jj][4]-hist_orient[bin] + PI), 2*PI) - PI;
            temp[5] = fmod((diff[jj][5]-hist_orient[bin] + PI), 2*PI) - PI;
            temp[6] = fmod((diff[jj][6]-hist_orient[bin] + PI), 2*PI) - PI;
            temp[7] = fmod((diff[jj][7]-hist_orient[bin] + PI), 2*PI) - PI;
            temp[8] = fmod((diff[jj][8]-hist_orient[bin] + PI), 2*PI) - PI;
            
/*
            temp[0] = (fabs(temp[0]) <= hist_step)?1:0;
            temp[1] = (fabs(temp[1]) <= hist_step)?1:0;
            temp[2] = (fabs(temp[2]) <= hist_step)?1:0;
            temp[3] = (fabs(temp[3]) <= hist_step)?1:0;
            temp[4] = (fabs(temp[4]) <= hist_step)?1:0;
*/
/*
            temp[0] = diff[jj][0]-hist_orient[bin] ;
            temp[1] = diff[jj][1]-hist_orient[bin] ;
            temp[2] = diff[jj][2]-hist_orient[bin] ;
            temp[3] = diff[jj][3]-hist_orient[bin] ;
            temp[4] = diff[jj][4]-hist_orient[bin] ;
            temp[5] = diff[jj][5]-hist_orient[bin] ;
            temp[6] = diff[jj][6]-hist_orient[bin] ;
            temp[7] = diff[jj][7]-hist_orient[bin] ;
            temp[8] = diff[jj][8]-hist_orient[bin] ;
*/


            temporary = max((1- fabs(temp[0])/hist_step), 0);
            temp[0] = temporary;
            temporary = max((1- fabs(temp[1])/hist_step), 0);
            temp[1] = temporary;
            temporary = max((1- fabs(temp[2])/hist_step), 0);
            temp[2] = temporary;
            temporary = max((1- fabs(temp[3])/hist_step), 0);
            temp[3] = temporary;
            temporary = max((1- fabs(temp[4])/hist_step), 0);
            temp[4] = temporary;
            temporary = max((1- fabs(temp[5])/hist_step), 0);
            temp[5] = temporary;
            temporary = max((1- fabs(temp[6])/hist_step), 0);
            temp[6] = temporary;
            temporary = max((1- fabs(temp[7])/hist_step), 0);
            temp[7] = temporary;
            temporary = max((1- fabs(temp[8])/hist_step), 0);
            temp[8] = temporary;


            temp[0] = temp[0]*fabs(weight[jj][0]);
            temp[1] = temp[1]*fabs(weight[jj][1]);
            temp[2] = temp[2]*fabs(weight[jj][2]);
            temp[3] = temp[3]*fabs(weight[jj][3]);
            temp[4] = temp[4]*fabs(weight[jj][4]);
            temp[5] = temp[5]*fabs(weight[jj][5]);
            temp[6] = temp[6]*fabs(weight[jj][6]);
            temp[7] = temp[7]*fabs(weight[jj][7]);
            temp[8] = temp[8]*fabs(weight[jj][8]);

            orient_hist[bin] += temp[0] + temp[1] + temp[2] + temp[3] + temp[4] + temp[5] + temp[6] + temp[7] + temp[8] ;

            
        }
                
    }

    /*  suppressing the nonmaxima  */

    float *peaks = (float *)malloc(num_bins*sizeof(float *));

    for ( bin = 1; bin <= num_bins-2; bin++)
    {
        if ((orient_hist[bin] < orient_hist[bin-1]) || (orient_hist[bin] < orient_hist[bin+1]))
        {
            peaks[bin] = 0;
        }
        else
        {
            peaks[bin] = orient_hist[bin];

        }
    }

    if ((orient_hist[0] < orient_hist[1]) || (orient_hist[0] < orient_hist[num_bins-1]))
        peaks[0] = 0;
    else
        peaks[0] = orient_hist[0];


    if ((orient_hist[num_bins-1] < orient_hist[num_bins-2]) || (orient_hist[num_bins-1] < orient_hist[0]))
        peaks[num_bins-1] = 0;
    else
        peaks[num_bins-1] = orient_hist[num_bins-1];

    int peak_pos = finding_peak(peaks, num_bins);

    float max = peaks[peak_pos];

    float peak_value = max;

    float max_orient[10]={0};
    
    int index = 0;
    
    while (peak_value > 0.8*max)
    {


        float a = hist_orient[peak_pos]-hist_step;
        float b = hist_orient[peak_pos];
        float c = hist_orient[peak_pos]+hist_step;


        bin = fmod( peak_pos -1 + num_bins - 1, num_bins ) + 1;
        float fa = orient_hist[bin];

        bin = fmod( peak_pos + num_bins - 1, num_bins ) + 1;
        float fb = orient_hist[bin];

        bin = fmod( peak_pos + 1 + num_bins - 1, num_bins ) + 1;
        float fc = orient_hist[bin];
    
        max_orient[index] = parabolic_interp(a, fa, b, fb, c, fc);
                    
        while( max_orient[index] < -PI )
        {
            max_orient[index] = max_orient[index] + 2*PI;
        }

        while( max_orient[index] >= PI )
        {
            max_orient[index] = max_orient[index] - 2*PI;
        }  

        //max_orient[index] = hist_orient[peak_pos];
    
        peaks[peak_pos] = 0;
        peak_pos = finding_peak(peaks, num_bins);
        peak_value = peaks[peak_pos];

        index++;

    }

    free(peaks);
    free(orient_hist);

    for (int iii = 0; iii <= 4; iii++)
        free(gaussianKernel[iii]);
    free(gaussianKernel);

    free(magnitude_pyrmiad);
    
    free(gradient_pyrmiad);

    //free(nor_image);


    float degree[10]={0};
    
    degree[0] = max_orient[0]/PI*180;
    degree[1] = max_orient[1]/PI*180;
    degree[2] = max_orient[2]/PI*180;

    return max_orient[0];

}

void local_descriptor_calc(BYTE *image, int width, int height)
{

    float feature_orentation =  orientation_calc_test(image, width, height);

    unsigned int *image_i = byte2int(image, width, height);


    float *nor_image = normalize(image_i, width, height); 


    free(image_i);


    typedef struct ftag{  
           float x;
           float y;
         
           }fCoords; 


    int i,j;

    int ii, jj;

    int i_cor, j_cor;

    int grid_space = 4;
    
    int feature_weight_window = 8;

    float orient_bin_step = PI/4.0;
    
    float orient_angles[8] = {-PI, -PI+orient_bin_step, -PI+2*orient_bin_step, -PI+3*orient_bin_step, 0, 
                                   orient_bin_step, orient_bin_step*2, orient_bin_step*3};
                                       
            
    float center_cor[4] = {-6.0, -2.0, 2.0, 6.0};       
    
               
    fCoords center_coordinations[4][4];
    
    
    float sample_cor[16] = {-7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5};
    
    fCoords sample_coordinations[16][16];

    
    for (i_cor = 0; i_cor <=3; i_cor++)
    {
        
        for (j_cor = 0; j_cor <=3; j_cor++)
        {
            center_coordinations[i_cor][j_cor].x = center_cor[i_cor];
            center_coordinations[i_cor][j_cor].y = center_cor[j_cor];
       
        }        
    }

    for (i_cor = 0; i_cor <=15; i_cor++)
    {
        for (j_cor = 0; j_cor <=15; j_cor++)
        {
            sample_coordinations[i_cor][j_cor].x = sample_cor[i_cor];
            sample_coordinations[i_cor][j_cor].y = sample_cor[j_cor];
          
        }
    }

    float rotation_matrix[2][2];
    rotation_matrix[0][0] = cos(feature_orentation); rotation_matrix[0][1] = -sin(feature_orentation); 
    rotation_matrix[1][0] = sin(feature_orentation); rotation_matrix[1][1] =  cos(feature_orentation);
        
          
    float temp1, temp2;
        
    int iii, jjj;
    
    for (iii = 0; iii <=3; iii++)
    {
        for (jjj = 0; jjj <=3; jjj++)
        {
            temp1 = rotation_matrix[0][0]*center_coordinations[iii][jjj].x 
                                         + rotation_matrix[0][1]*center_coordinations[iii][jjj].y;

                
            
            temp2 = rotation_matrix[1][0]*center_coordinations[iii][jjj].x 
                                         + rotation_matrix[1][1]*center_coordinations[iii][jjj].y;

                
                                         
            center_coordinations[iii][jjj].x = temp1;
            
            center_coordinations[iii][jjj].y = temp2;
        }
    }
    
    for (iii = 0; iii <=15; iii++)
    {
        for (jjj = 0; jjj <=15; jjj++)
        {
            temp1 = rotation_matrix[0][0]*sample_coordinations[iii][jjj].x 
                                         + rotation_matrix[0][1]*sample_coordinations[iii][jjj].y;
                
            
            temp2 = rotation_matrix[1][0]*sample_coordinations[iii][jjj].x 
                                         + rotation_matrix[1][1]*sample_coordinations[iii][jjj].y;
            
                                         
            sample_coordinations[iii][jjj].x = temp1;
            
            sample_coordinations[iii][jjj].y = temp2;
        }
    }



    int x = 300;

    int y = 300;

    float feature_descriptor[128] = {0.0f}; 
    float x_weight[4][4], y_weight[4][4];
    float position_weight[4][4];
    float orientation_weight[8]; 
    float gaussian_weight;

    for (iii = 0; iii <=3; iii++)
    {
        for (jjj = 0; jjj <=3; jjj++)
        {

            temp1 = center_coordinations[iii][jjj].x+x;

            if (temp1 < 0) temp1 = 0;
            if (temp1 > width-1) temp1 = width-1;

            center_coordinations[iii][jjj].x = temp1; 

            temp2 = center_coordinations[iii][jjj].y+y;

            if (temp2 < 0) temp2 = 0;
            if (temp2 > height-1) temp2 = height-1;
                                
            center_coordinations[iii][jjj].y = temp2;
                
        }
    }
  
    for(iii = 0; iii <= 15; iii++)
    {        
        for(jjj = 0; jjj <= 15; jjj++ )
        {

            
            temp1 = sample_coordinations[iii][jjj].x+x;

            if (temp1 < 0) temp1 = 0;
            if (temp1 > width-1) temp1 = width-1;

            sample_coordinations[iii][jjj].x = temp1; 
                
            temp2 = sample_coordinations[iii][jjj].y+y;

            if (temp2 < 0) temp2 = 0;
            if (temp2 > height-1) temp2 = height-1;
                                
            sample_coordinations[iii][jjj].y = temp2;
                            
            float samples[3][3];
              

            float x_temp;// = sample_coordinations[iii][jjj].x;
            float y_temp;// = sample_coordinations[iii][jjj].y-1;
            float xy_temp1, xy_temp2;

            x_temp = -1;
            y_temp = 0;

            xy_temp1 = rotation_matrix[0][0]*x_temp 
                                         + rotation_matrix[0][1]*y_temp;
                        
            
            xy_temp2 = rotation_matrix[1][0]*x_temp
                                         + rotation_matrix[1][1]*y_temp;

            x_temp = xy_temp1 + sample_coordinations[iii][jjj].x;
            y_temp = xy_temp2 + sample_coordinations[iii][jjj].y;
            
            
            //samples[1][0] = interp2(nor_image, width, height, sample_coordinations[iii][jjj].x, sample_coordinations[iii][jjj].y-1);
            samples[1][0] = interp2(nor_image, width, height, x_temp, y_temp);

            x_temp = 1;
            y_temp = 0;

            xy_temp1 = rotation_matrix[0][0]*x_temp 
                                         + rotation_matrix[0][1]*y_temp;
            
                
            
            xy_temp2 = rotation_matrix[1][0]*x_temp
                                         + rotation_matrix[1][1]*y_temp;

            x_temp = xy_temp1 + sample_coordinations[iii][jjj].x;
            y_temp = xy_temp2 + sample_coordinations[iii][jjj].y;
              
            //samples[1][2] = interp2(nor_image, width, height, sample_coordinations[iii][jjj].x, sample_coordinations[iii][jjj].y+1);
            samples[1][2] = interp2(nor_image, width, height,x_temp, y_temp);
            
            x_temp = 0;
            y_temp = 1;

            xy_temp1 = rotation_matrix[0][0]*x_temp 
                                         + rotation_matrix[0][1]*y_temp;
            
                
            
            xy_temp2 = rotation_matrix[1][0]*x_temp
                                         + rotation_matrix[1][1]*y_temp;

            x_temp = xy_temp1 + sample_coordinations[iii][jjj].x;
            y_temp = xy_temp2 + sample_coordinations[iii][jjj].y;

            samples[2][1] = interp2(nor_image, width, height, x_temp, y_temp);
            //samples[2][1] = interp2(nor_image, width, height, sample_coordinations[iii][jjj].x+1, sample_coordinations[iii][jjj].y);
            
            x_temp = 0;
            y_temp = -1;

            xy_temp1 = rotation_matrix[0][0]*x_temp 
                                         + rotation_matrix[0][1]*y_temp;
            
                
            
            xy_temp2 = rotation_matrix[1][0]*x_temp
                                         + rotation_matrix[1][1]*y_temp;

            x_temp = xy_temp1 + sample_coordinations[iii][jjj].x;
            y_temp = xy_temp2 + sample_coordinations[iii][jjj].y;

            samples[0][1] = interp2(nor_image, width, height, x_temp, y_temp);

            //samples[0][1] = interp2(nor_image, width, height, sample_coordinations[iii][jjj].x-1, sample_coordinations[iii][jjj].y);
              
            float diff_x = 0.5*(samples[1][2] - samples[1][0]);
            float diff_y = 0.5*(samples[2][1] - samples[0][1]);
              
            float mag_sample = sqrt(diff_x*diff_x + diff_y*diff_y);
            float grad_sample = (float)atan2(diff_y, diff_x);
           
            temp1 = (float)grad_sample*180/PI;
                    
            grad_sample = (temp1 == 180)?(-grad_sample):grad_sample;
                                  
            for (ii = 0; ii <= 3; ii++)
            {
                for (jj = 0; jj <= 3; jj++)
                {
                    x_weight[ii][jj] = max(1-fabs( sample_coordinations[iii][jjj].x - center_coordinations[ii][jj].x)/grid_space, 0);
                    y_weight[ii][jj] = max(1-fabs( sample_coordinations[iii][jjj].y - center_coordinations[ii][jj].y)/grid_space, 0);
                    position_weight[ii][jj] = x_weight[ii][jj]*y_weight[ii][jj]; 
                }
            }
              
            for(ii = 0; ii <= 7; ii++)
            {
                //float difference = fmod(grad_sample - feature_orentation - orient_angles[ii] + PI, 2*PI) - PI;
                float difference = fmod(grad_sample - orient_angles[ii] + PI, 2*PI) - PI;
                orientation_weight[ii] = max(1 - fabs(difference)/orient_bin_step, 0);
                   
            }
              
            gaussian_weight = exp(-(pow((sample_coordinations[iii][jjj].x - x),2) + pow((sample_coordinations[iii][jjj].y-y), 2))/pow(2*feature_weight_window, 2))/pow(2*PI*feature_weight_window, 2);
              
            for(ii = 0; ii <= 7; ii++)
            {
                feature_descriptor[ii + 8*0]  += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[0][0];
                feature_descriptor[ii + 8*1]  += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[0][1];
                feature_descriptor[ii + 8*2]  += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[0][2];
                feature_descriptor[ii + 8*3]  += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[0][3];   
                feature_descriptor[ii + 8*4]  += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[1][0];
                feature_descriptor[ii + 8*5]  += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[1][1];
                feature_descriptor[ii + 8*6]  += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[1][2];
                feature_descriptor[ii + 8*7]  += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[1][3];
                feature_descriptor[ii + 8*8]  += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[2][0];
                feature_descriptor[ii + 8*9]  += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[2][1];
                feature_descriptor[ii + 8*10] += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[2][2];
                feature_descriptor[ii + 8*11] += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[2][3];
                feature_descriptor[ii + 8*12] += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[3][0];
                feature_descriptor[ii + 8*13] += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[3][1];
                feature_descriptor[ii + 8*14] += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[3][2];
                feature_descriptor[ii + 8*15] += orientation_weight[ii]*gaussian_weight*mag_sample*position_weight[3][3];
                  
            }
        }
    }       
         
      
      
    
    float feature_descriptor_norm = 0.0f;
      
    for(iii = 0; iii <= 127; iii++)
    {
        feature_descriptor_norm += feature_descriptor[iii]*feature_descriptor[iii];
    } 
    feature_descriptor_norm = sqrt(feature_descriptor_norm);
    for(iii = 0; iii <= 127; iii++)
    {
        feature_descriptor[iii] = feature_descriptor[iii]/feature_descriptor_norm;
        if (feature_descriptor[iii] > 0.2)
            feature_descriptor[iii] = 0.2;       
    } 
    feature_descriptor_norm = 0;
    for(iii = 0; iii <= 127; iii++)
    {
        feature_descriptor_norm += feature_descriptor[iii]*feature_descriptor[iii];
    }     
    for(iii = 0; iii <= 127; iii++)
    {
        feature_descriptor[iii] = feature_descriptor[iii]/feature_descriptor_norm;
    }  
    
    FILE * pFile;
    pFile = fopen ("test_descriptor.key","w");

    for(iii = 0; iii <= 127; iii++)
    {
        fprintf(pFile, "%f ", feature_descriptor[iii]);
    } 
    fprintf(pFile, "\r\n"); 

    fclose(pFile);


}











