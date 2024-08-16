/*-----------------------------------
 * MeshSizeIterativeSolver
 * MEX file
 *
 * h = MeshSizeIter(hn,h0,g,delta);
 *
 * input: 
 * output:
 *
 * -----------------------------------*/

// Include Libraries
#include "mex.h"
#include <math.h>
#include <matrix.h>
#include <limits.h>

// Define variables and functions
#define eps 2.2204e-16
#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))


// Mex Function
void mexFunction(int nOUT, mxArray *pOUT[], int nIN, const mxArray *pIN[])
{
        
    // Define variable types
    mwSize LY,LX,N,L;
    int i, j, k, ki, *r, *c;
    double *h, *h0, delta, g, hmax, hmin;
    double tol, R, xfordiff, xbackdiff, yfordiff, ybackdiff, dt, Delta;
    double inf,hn;
            
    // Get sizes
    LY      = mxGetM(pIN[0]);
    LX      = mxGetN(pIN[0]);
    N       = LY*LX;
    L       = mxGetM(pIN[1])*mxGetN(pIN[1]); 
    
    // Read in inputs
    h0      = mxGetPr(pIN[0]);              // Initial mesh size
    r       = mxGetPr(pIN[1]);              // Row index (0-based)
    c       = mxGetPr(pIN[2]);              // Column index (0-based)
    hmax    = mxGetScalar(pIN[3]);          // Maximum element size
    hmin    = mxGetScalar(pIN[4]);          // Minimum element size
    g       = mxGetScalar(pIN[5]);          // Gradient ratio
    delta   = mxGetScalar(pIN[6]);          // Grid Spacing
    
    // Get infinity variable
    inf = mxGetInf();

    // Prepare output
    pOUT[0] = mxCreateDoubleMatrix(LY, LX, mxREAL);
    h       = mxGetPr(pOUT[0]);
    
    // Define variables
    dt  = delta/2; // Time step
    tol = (10e-6); // Tolerance for convergence
    R   = 0;       // Residual    
    
    // While residual is less than specified tolerance
     while(1)
     {
         
         R = 0; // Initialize residual
         
         for ( i = 0; i<L; i++ )
         {
             
             k = LY * r[i] + c[i]; // 0-based linear indexing
             
//              // Compute upwind differences
//              ki = LY * (r[i]+1) + c[i]; // 0-based
//              xfordiff    = min((h0[ki] - h0[k])/delta,0);
//              xfordiff    = xfordiff*xfordiff;
//              
//              ki = LY * (r[i]-1) + c[i]; // 0-based
//              xbackdiff   = max((h0[k] - h0[ki])/delta,0);
//              xbackdiff   = xbackdiff*xbackdiff;
//              
//              ki = LY * (r[i]) + c[i]+1; // 0-based
//              yfordiff    = min((h0[ki] - h0[k])/delta,0);
//              yfordiff    = yfordiff*yfordiff;
//              
//              ki = LY * (r[i]) + c[i]-1; // 0-based
//              ybackdiff   = max((h0[k] - h0[ki])/delta,0);
//              ybackdiff   = ybackdiff*ybackdiff;
//              
//              // Compute Delta
//              Delta = sqrt(xfordiff + xbackdiff + yfordiff + ybackdiff);
//              
//              // Compute next time step
//              hn = h0[k] + dt*(min(Delta,g) - Delta);
//              
//              // Compute max Residual
//              R = fabs( (hn-h0[k]) ) + R;
             
             // Next time-step
             h0[k] = hn;
 
         }
         
         // Check tolerance
         if ( R <= tol)
         {
             break;
         }

         
         
//          for (i = 1; i<LX-1; i++) // Outer loop
//          {
//              for (j = 1; j<LY-1; j++) // Inner loop
//              {
//                   
//                  k = LY * i + j; // 0-based linear indexing
//                  
//                  //mexPrintf("\n k=%i",k);
//                  
//                  if( D[k] > delta)
//                  {
//                      continue;
//                  }
//                                  
//                  // Compute upwind differences
//                  ki = LY * (i+1) + j; // 0-based
//                  xfordiff    = min((h0[ki] - h0[k])/delta,0);
//                  xfordiff    = xfordiff*xfordiff;
//                  
//                  ki = LY * (i-1) + j; // 0-based
//                  xbackdiff   = max((h0[k] - h0[ki])/delta,0);
//                  xbackdiff   = xbackdiff*xbackdiff;
//                  
//                  ki = LY * (i) + j+1; // 0-based
//                  yfordiff    = min((h0[ki] - h0[k])/delta,0);
//                  yfordiff    = yfordiff*yfordiff;
//                  
//                  ki = LY * (i) + j-1; // 0-based
//                  ybackdiff   = max((h0[k] - h0[ki])/delta,0);
//                  ybackdiff   = ybackdiff*ybackdiff;
//                  
//                  // Compute Delta
//                  Delta = sqrt(xfordiff + xbackdiff + yfordiff + ybackdiff);
//                  
//                  // Compute next time step
//                  hn = h0[k] + deltat*(min(Delta,g) - Delta);
//                  
//                  // Compute max Residual
//                  R = fabs( (hn-h0[k]) ) + R;
//                  
//                  // Next time-step
//                  h0[k] = hn;
//                  
//              } // j loop
//          } // i loop
//          
//          // Check tolerance
//          if ( R <= tol)
//          {
//              break;
//          }
//          
     } // while loop
    
    for (i = 0; i<LX; i++)
    {
        for (j = 0; j<LY; j++)
        {
            k = LY * i + j; // 0-based

            h[k] = h0[k];
        }
    }
        
}
