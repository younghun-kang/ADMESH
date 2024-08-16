//
//  SpacePolyPoints.c
//
//
//  Created by Dustin West on 11/26/13.
//
//
/*-----------------------------------
 * MEX file
 *
 * -----------------------------------*/

/*************************************************************************
 * %  To compile
 *  mex SpacePolyPoints.c
 **************************************************************************/

/* Author:                                                              */
/* Dustin West      email: west.425@osu.edu                             */
#include "mex.h"
#include <math.h>
#include <matrix.h>
#include <limits.h>

#define EPS 2.2204e-16

#define MAX(a,b)    ((a) > (b) ? (a) : (b))
#define MIN(a,b)    ((a) > (b) ? (b) : (a))

void mexFunction(
        int nOUT, mxArray *pOUT[],
        int nINP, const mxArray *pINP[])
{
    int nPM,nPN;
    int nP;
    int nSM,nSN;
    int nS;
    int nSMP,nSNP;
    int nSP;
    int i, j, ind0, ind1, flag;
    double *polyX, *polyY, *s, *s_spaced;
    double *px, *py;
    
    double pt0[2], pt1[2], w0, w1;
    double inf;
    
    
    /* unpack inputs */
    // (x,y) coordinates to polygon
    nPM     = mxGetM(pINP[0]);
    nPN     = mxGetN(pINP[0]);
    nP      = nPM*nPN;
    polyX   = mxGetPr(pINP[0]);
    polyY   = mxGetPr(pINP[1]);
    
    // Cumulative distance vector
    nSM     = mxGetM(pINP[2]);
    nSN     = mxGetN(pINP[2]);
    nS      = nSM*nSN;
    s       = mxGetPr(pINP[2]);
    
    // Evenly distributed vector
    nSMP     = mxGetM(pINP[3]);
    nSNP     = mxGetN(pINP[3]);
    nSP      = nSMP*nSNP;
    s_spaced = mxGetPr(pINP[3]);
    
    // Assign infinity vector
    inf = mxGetInf();
    
    /* pack outputs */
    
    pOUT[0]      = mxCreateDoubleMatrix(nSMP, nSNP, mxREAL);
    px           = mxGetPr(pOUT[0]);
    pOUT[1]      = mxCreateDoubleMatrix(nSMP, nSNP, mxREAL);
    py           = mxGetPr(pOUT[1]);
    
    flag = 0;
    
    /* Perform Calculations */
    for(i=0; i<nSP; i++)
    {
        //Determine the index of surrounding points
        for(j=0; j<nS; j++)
        {
            // Find the last logically true comparison
            if( s[j] <= s_spaced[i] )
            {
                ind0 = j;
            }
            // Find the first logically true comparison
            if( s[j] >= s_spaced[i] && flag == 0 )
            {
                ind1 = j;
                flag = 1;
            }
            
        } // End j loop
        
        // Reset ind1 equal to 0
        flag = 0;
        
        //mexPrintf("\n ind0 is %i\n", ind0+1);
        //mexPrintf("\n ind1 is %i\n", ind1+1);
        
        // For beginning & ending points
        if( ind0 == ind1 )
        {
            px[i] = polyX[ind0];
            py[i] = polyY[ind0];
        }
        else
        {
            // Store surrounding points
            pt0[0] = polyX[ind0];
            pt0[1] = polyY[ind0];
            pt1[0] = polyX[ind1];
            pt1[1] = polyY[ind1];
            
            // Determine the weights associated to each surrounding point
            w0 = s_spaced[i] - s[ind0];
            w1 = s[ind1] - s_spaced[i];
            
            // Compute weighted average with neighbor positions
            px[i] = (pt0[0] * w1 + pt1[0] * w0) / (w0 + w1);
            py[i] = (pt0[1] * w1 + pt1[1] * w0) / (w0 + w1);
            
        }
        

    } // End i loop
    
}