//
//  Point2PolyDistance.c
//
//
//  Created by Dustin West on 11/25/13.
//
//
/*-----------------------------------
 * MEX file
 *
 * S = Point2PolyDistance(px,py,cx,cy);
 *
 * input: px - nPM x nPN (Matrix) X-coordinates of the background grid
 *        py - nPM x nPN (Matrix) Y-coordinates of the background grid
 *        cx - nC  x 1   (Vector) X-coordinates of the polygon
 *        cy - nC  x 1   (Vector) Y-coordinates of the polygon
 * output:
 *        S            - Distnace Function, nPM x nPN
 *
 * -----------------------------------*/

/*************************************************************************
 * %  To compile
 *  mex Point2EdgeDistance.c
 **************************************************************************/

/* Author:                                                              */
/* Dustin West      email: west.425@osu.edu                             */

// Include Libraries
#include "mex.h"
#include <math.h>
#include <matrix.h>
#include <limits.h>

// Define variables & basic functions
#define eps 1.0e-10
#define MAX(a,b)    ((a) > (b) ? (a) : (b))
#define MIN(a,b)    ((a) > (b) ? (b) : (a))
#define abs(a)      ((a) < eps ? (-a) : (a))

void mexFunction(int nOUT, mxArray *pOUT[], int nIN, const mxArray *pIN[])
{
    //--------------------------------------------------------------------------
    // MEX Gateway
    //--------------------------------------------------------------------------
    int nX,nY,nU,nV;
    int nP, nC;
    int k, j;
    int ind;
    double *S, *X, *Y, *u, *v;
    double L2, dmin, t, d, up, vp, du, dv, duX, dvY;
    double inf;    
    
    //--------------------------------------------------------------------------
    // Unpack Inputs 
    //--------------------------------------------------------------------------
    
    // Background grid
    nX = mxGetM(pIN[0]); nY = mxGetN(pIN[0]); 
    nP = nX*nY;
    X = mxGetPr(pIN[0]); Y = mxGetPr(pIN[1]);
    
    // Polygon points
    nU = mxGetM(pIN[2]); nV = mxGetN(pIN[2]); 
    nC = nU*nV;
    u = mxGetPr(pIN[2]);
    v = mxGetPr(pIN[3]);
    
    // Get infinity variable
    inf = mxGetInf();
    
    //--------------------------------------------------------------------------
    // Pack Outputs
    //--------------------------------------------------------------------------
    pOUT[0]      = mxCreateDoubleMatrix(nX, nY, mxREAL);
    S            = mxGetPr(pOUT[0]);
        
    
    //--------------------------------------------------------------------------
    // Enter main loop. Loop over each background grid point.
    //--------------------------------------------------------------------------
    // Loop over grid and compute nearest distance
    for (j = 0; j<nP; j++)
    {
        
        // initialize minimum distance to infinity
        dmin = inf;
        
        //---------------------------------------------------
        // Compute the distance from (u,v) to (X,Y)
        //---------------------------------------------------

        for (k = 0; k<nC-1; k++)
        {
            
            // Compute squared length of edge
            du = u[k] - u[k+1]; dv = v[k] - v[k+1];
            L2 = du*du + dv*dv;
            
            // Consider the line extending the segment, parameterized as 
            // v + t (w - v).
            // We find projection of point p onto the line.
            // It falls where t = [(p-v) . (w-v)] / |w-v|^2
            t = ((X[j] - u[k+1])*du + (Y[j] - v[k+1])*dv)/L2;
            
            if ( t < eps ) // Beyond the 'k+1' end of the segment
            {
                duX = u[k+1]-X[j]; dvY = v[k+1]-Y[j];
                d = duX*duX + dvY*dvY;
            }
            else if ( t > 1.00 ) // Beyond the 'k' end of the segment
            {
                duX = u[k]-X[j]; dvY = v[k]-Y[j];
                d = duX*duX + dvY*dvY;
            }
            else // Projection falls on the segment
            {
                up = u[k+1] + t*du;
                vp = v[k+1] + t*dv;
                
                duX = up-X[j]; dvY = vp-Y[j];
                d = duX*duX + dvY*dvY;
            }
            
            // Save minimum value
            dmin = MIN(dmin,d);

        }
        
        // Store the minimum distance (squared here to save computation time)
        S[j] = sqrt(dmin);

        
    } // End j Loop
    
    //**************************************************************************

}