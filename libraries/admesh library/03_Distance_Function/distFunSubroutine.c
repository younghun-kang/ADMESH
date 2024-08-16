//
//  distFunSubroutine.c
//
//
//  Created by Dustin West on 05/12/14.
//
//
/*-----------------------------------
 * MEX file
 *
 * S = distFunSubroutine(X,Y,u,v,id);
 *
 * input: X     - nPM x nPN (Matrix) X-coordinates of the background grid
 *        Y     - nPM x nPN (Matrix) Y-coordinates of the background grid
 *        u     - nC  x 1   (Vector) X-coordinates of the polygon
 *        v     - nC  x 1   (Vector) Y-coordinates of the polygon
 *        id    - nPM x nPN (Vector) u & v indexing vector
 * output:
 *        S     - Distance Function, nPM x nPN
 *
 * -----------------------------------*/

/*************************************************************************
 * %  To compile
 *  mex distFunSubroutine.c
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
    int k, n0, n1, n2;
    double *S, *X, *Y, *u, *v, *id;
    double L2, t, d, up, vp, du, dv, duX, dvY;
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
    
    // Get id vector
    id = mxGetPr(pIN[4]);
    
    // Get infinity variable
    inf = mxGetInf();
    
    //--------------------------------------------------------------------------
    // Pack Outputs
    //--------------------------------------------------------------------------
    pOUT[0]      = mxCreateDoubleMatrix(nX, nY, mxREAL);
    S            = mxGetPr(pOUT[0]);
    
    //--------------------------------------------------------------------------
    // Begin computational routine
    //--------------------------------------------------------------------------
    
    // Convert ID to 0 indexing order. Check ID, if ID == nP-1, ID = 0
    for(k = 0; k < nP; k++)
    {
        S[k] = inf;
        
        id[k] = id[k] - 1;
        
        if(id[k] == nP-1)
        {
            id[k] = 0;
        }
    }
    
    // Loop over grid and compute nearest distance
    for(k = 0; k < nP; k++)
    {
        // Looking at 2 edge segments
        //
        // o-------------o-------------o
        // ID-1          ID           ID+1  (Index representation)
        // n0            n1            n2   (Variable representation)
        //
        
        // Check ID & ID+1 segment first
        n1 = (int)id[k];
        n2 = (int)id[k]+1;
        
        // Compute squared length of edge
        du = u[n1] - u[n2];
        dv = v[n1] - v[n2];
        L2 = du*du + dv*dv;
        
        // Consider the line extending the segment, parameterized as
        // v + t (w - v).
        // We find projection of point p onto the line.
        // It falls where t = [(p-v) . (w-v)] / |w-v|^2
        t = ((X[k] - u[n2])*du + (Y[k] - v[n2])*dv)/L2;
        
        if( (t < eps) ) // t < eps Projection falls on'ID+1' end of the segment
        {
            // No computation needed. Just return inf
            continue;
        }
        else if(t > 1.00)
        {
            
            // Check ID-1 & ID segment
            n1 = (int)id[k];
            n0 = (int)id[k]-1;
            
            // If n0 is less than 0 then we need to look at the second to last
            // point in [px,py].
            if(n0 < 0)
            {
                n0 = nC-2;
            }
            
            // Compute squared length of edge
            du = u[n0] - u[n1];
            dv = v[n0] - v[n1];
            L2 = du*du + dv*dv;
            
            // Consider the line extending the segment, parameterized as
            // v + t (w - v).
            // We find projection of point p onto the line.
            // It falls where t = [(p-v) . (w-v)] / |w-v|^2
            t = ((X[k] - u[n1])*du + (Y[k] - v[n1])*dv)/L2;
            
            if( (t > eps) && (t < 1))
            {
                up = u[n1] + t*du;
                vp = v[n1] + t*dv;
                duX = up-X[k]; dvY = vp-Y[k];
                S[k] = sqrt(duX*duX + dvY*dvY);
            }
        }
        else
        {
            up = u[n2] + t*du;
            vp = v[n2] + t*dv;
            duX = up-X[k]; dvY = vp-Y[k];
            S[k] = sqrt(duX*duX + dvY*dvY);
        }
        
    }
    
    //**************************************************************************
    
}