#include <stdio.h>

double getdeltaT(int i, int poly_deg, double* knots)
{
    return knots[i+poly_deg+1] - knots[i+1];
}

double getKnots(double start_knot, double* knot_steps, int poly_deg, int index, int step_length)
{
    if (index <= poly_deg)
    {
        return start_knot;
    }
    else if (index <= poly_deg+step_length)
    {
        double knot = start_knot;
        for (int i = 0; i < index-poly_deg-1; i++)
        {
            knot += knot_steps[i];
        }
        return knot;
    }
    else
    {
        double knot = start_knot;
        for (int i = 0; i < step_length; i++)
        {
            knot += knot_steps[i];
        }
        return knot;
    }
    
}

double controlPointDerivative(int i, int derivative, int poly_deg, int n, double start_knot, double* knot_steps, double* control_points, int step_length)
{
    if (derivative == 0)
    {
        return control_points[i];
    }
    else
    {
        double temp1 = getKnots(start_knot, knot_steps, poly_deg, i+poly_deg+2, step_length);
        double temp2 = getKnots(start_knot, knot_steps, poly_deg, i+derivative+1, step_length);

        double denom = temp1 - temp2;

        return ((poly_deg - derivative + 1)/denom) * (controlPointDerivative(i+1, derivative-1, poly_deg, n, start_knot, knot_steps, control_points, step_length) - controlPointDerivative(i, derivative-1, poly_deg, n, start_knot, knot_steps, control_points, step_length));
    }
    
}

int main()
{
    int d = 1;
    int k = 5;
    int n = 11;
    int knot1 = 0;

    double control[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    double knot_steps[7] = {3, 3, 3, 3, 3, 3, 3};

    int opto_vars = 19;

    int step_length = 7;

    double knots[step_length+(2*k)+1];
    for (int i = 0; i < k+1; i++)
    {
        knots[i] = 0;
    }
    for (int i = k+1; i < step_length+k+1; i++)
    {
        knots[i] = knots[i-1] + knot_steps[i-k-1];
    }
    for (int i = step_length+k+1; i < step_length+(2*k)+1; i++)
    {
        knots[i] = knots[i-1];
    }

    double knots_v[step_length+(2*k)-1];
    for (int i = 0; i < step_length+(2*k)-1; i++)
    {
        knots_v[i] = knots[i+1];
    }

    double knots_a[step_length+(2*k)-3]; // length is wrong
    for (int i = 0; i < step_length+(2*k)-3; i++)
    {
        knots_a[i] = knots_v[i+1];
    }
    
    int n_c = 12;

    double d1_control[n_c-d];
    double d1_dtmethod[n_c-d];
    double d2_control[n_c-d-1];
    double d2_dtmethod[n_c-d-1];
    double d3_control[n_c-d-2];
    double d3_dtmethod[n_c-d-2];

    int m = 27;

    double result[m];

    for (int i = 0; i < n; i++)
    {
        d1_control[i] = controlPointDerivative(i, d, k, n, 0, knot_steps, control, step_length);
        d1_dtmethod[i] = (control[i+1]-control[i])*k/getdeltaT(i, k, knots);
        result[i] = (d1_control[i]*d1_control[i]) - 4; // constraint
    }

    for (int i = 0; i < n_c-d-1; i++)
    {
        d2_control[i] = controlPointDerivative(i, d, k-1, n_c-d-1, 0, knot_steps, d1_control, step_length);
        d2_dtmethod[i] = (d1_dtmethod[i+1]-d1_dtmethod[i])*(k-1)/getdeltaT(i, k-1, knots_v);
        result[i+n] = (d2_control[i]*d2_control[i]) - 9;
    }

    for (int i = 0; i < n_c-d-2; i++)
    {
        d3_control[i] = controlPointDerivative(i, d, k-2, n_c-d-2, 0, knot_steps, d2_control, step_length);
        d3_dtmethod[i] = (d2_dtmethod[i+1]-d2_dtmethod[i])*(k-2)/getdeltaT(i, k-2, knots_a);
        result[n+n_c-d-1] = (d3_control[i]*d3_control[i]) - 16;
    }


    // This gradient and control point calculation both work. Note: getdeltaT has changed. Implement this with SLSQP
    
    double grad[n*opto_vars];
    int start = 0;
    //printf("i | j |     dt    |   grad   | \n");
    for (int i = 0; i < n; i++) // go through the first n constraints
    {
        for(int j = 0; j < opto_vars; j++)  // differentiate wrt to each of the opto variables
        {
            double dt = getdeltaT(i,k,knots);
            if (j < step_length)  // if we are taking wrt to a knot step then use this eq. BUT also need to make sure that knot step is used
            {
                if (j <= i && j >= i-k+1)   // chect that that knot step is used
                {
                    grad[i*opto_vars+j] = -2*d1_dtmethod[i]*k*(control[i+1]-control[i])/(dt*dt);
                }
                else
                {
                    grad[i*opto_vars+j] = 0;
                }
            }
            else    // differentiate wrt to the control points
            {
                if (j-step_length == i)
                {
                    grad[i*opto_vars+j] = -2*d1_dtmethod[i]*k/dt;
                }
                else if(j-step_length == i+1)
                {
                    grad[i*opto_vars+j] = 2*d1_dtmethod[i]*k/dt;
                }
                else
                {
                    grad[i*opto_vars+j] = 0;
                }
            }
            //printf("%d | %d | %f | %f \n", i, j, dt, grad[i*opto_vars+j]);
            printf("%f ", grad[i*opto_vars + j]);
        }
        printf("\n");
    }

    printf("\n");
    int k2 = k - 1;
    double grad2[(n-1)*opto_vars];
//    start = step_length*opto_vars;
    for (int i = 0; i < n-1; i++) // go through the second n-1 constraints
    {
        for(int j = 0; j < opto_vars; j++)  // differentiate wrt to each of the opto variables
        {
            double dt = getdeltaT(i, k2, knots_v);
            if (j < step_length)  // if we are taking wrt to a knot step then use this eq. BUT also need to make sure that knot step is used
            {
                if (j <= i && j >= i-k2+1)   // chect that that knot step is used
                {
                    grad2[(i*opto_vars)+j] = -2*d2_dtmethod[i]*k2*(d1_dtmethod[i+1]-d1_dtmethod[i])/(dt*dt);
                }
                else
                {
                    grad2[i*opto_vars+j] = 0; 
                }
                
            }
            else    // differentiate wrt to the control points
            {
                if (j-step_length == i)
                {
                    grad2[(i*opto_vars)+j] = -2*d2_dtmethod[i]*k2/dt;
                }
                else if (j-step_length == i+1)
                {
                    grad2[(i*opto_vars)+j] = 2*d2_dtmethod[i]*k2/dt;
                }
                
                else
                {
                    grad2[i*opto_vars+j] = 0;
                }
                
            }
            printf("%f ", grad2[(i*opto_vars)+j]);
        }
        printf("\n");
    }

    printf("\n");
    int k3 = k - 2;
    double grad3[(n-2)*opto_vars];
    for (int i = 0; i < n-2; i++) // go through the third n-2 constraints
    {
        for(int j = 0; j < opto_vars; j++)  // differentiate wrt to each of the opto variables
        {
            double dt = getdeltaT(i,k3,knots_a);
            if (j < step_length)  // if we are taking wrt to a knot step then use this eq. BUT also need to make sure that knot step is used
            {
                if (j <= i && j >= i-k3+1)   // chect that that knot step is used
                {
                    grad3[(i*opto_vars)+j] = -2*d3_dtmethod[i]*k3*(d2_dtmethod[i+1]-d2_dtmethod[i])/(dt * dt);
                }
                else
                {
                    grad3[i*opto_vars+j] = 0; 
                }
                
            }
            else    // differentiate wrt to the control points
            {
                if (j-step_length == i)
                {
                    grad3[(i*opto_vars)+j] = -2*d3_dtmethod[i]*k3/dt;
                }
                else if(j-step_length == i+1)
                {
                    grad3[(i*opto_vars)+j] = 2*d3_dtmethod[i]*k3/dt;
                }
                else
                {
                    grad3[i*opto_vars+j] = 0;
                }
                
            }
            printf("%f ", grad3[(i*opto_vars)+j]);
        }
        printf("\n");
    }

    printf("\n");
    // Concatenate grad, grad2 and grad3 to get the gradient vector that NLOPT needs

    printf("d1: \n");
    for (int i = 0; i < n_c-d; i++)
    {
        printf("%f ", d1_control[i]);
    }
    printf("\n");

    printf("d1_dtmethod: \n");
    for (int i = 0; i < n_c-d; i++)
    {
        printf("%f ", d1_dtmethod[i]);
    }
    printf("\n");
    printf("\n");

    printf("d2: \n");
    for (int i = 0; i < n_c-d-1; i++)
    {
        printf("%f ", d2_control[i]);
    }
    printf("\n");
    printf("d2_dtmethod: \n");
    for (int i = 0; i < n_c-d-1; i++)
    {
        printf("%f ", d2_dtmethod[i]);
    }
    printf("\n");
    printf("\n");

    printf("d3: \n");
    for (int i = 0; i < n_c-d-2; i++)
    {
        printf("%f ", d3_control[i]);
    }
    printf("\n");
    printf("d3_dtmethod: \n");
    for (int i = 0; i < n_c-d-2; i++)
    {
        printf("%f ", d3_dtmethod[i]);
    }
    printf("\n");
}