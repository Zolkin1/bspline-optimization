/*
This will use the nonlinear optimization tools in NLOPT to optimize a BSpline based trajectory and time the algorithm.
*/

#include <nlopt.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "optimizeSpline.h"

int iter_count = 0;
int d_con_count = 0;
int c_con_count = 0;

double prevX[19];

int main()
{
    bspline spline_dat;
    spline_dat.start_knot = 0;
    spline_dat.poly_deg = 5;
    spline_dat.num_control_points = 12;
    spline_dat.step_length = 7;
    spline_dat.v_max = 2;
    spline_dat.a_max = 3;
    spline_dat.j_max = 11;
    spline_dat.dist = 15;

    int n = spline_dat.step_length + spline_dat.num_control_points; // dimensionality of the problem

    unsigned int m = 27;// + 11; // double check this. it is the number of constraints
    double tol[m];
    for (int i = 0; i < m; i++)
    {
        tol[i] = 1e-8;
    }

    double tol2[4] = {1e-8, 1e-8, 1e-8, 1e-8};
    int dist = spline_dat.dist;
    double lb[n];      // lower boungs on the optimization params. maybe needs to be const
    double ub[n];
    // COBYLA does NOT work in some initial conditions 
    double x[19] = {.4, 5, 5, 5, 5, 5, .4, 0, 0, 0, 3*dist/10, 4*dist/10, 5*dist/10, 6*dist/10, 7*dist/10, 8*dist/10, dist, dist, dist};
    for (int i = 0; i < 19; i++)
    {
        prevX[i] = x[i];
    }

    for (int i = 0; i < spline_dat.step_length; i++)
    {
        lb[i] = 1e-4;  // knot step lower bounds
        ub[i] = 200; // spam upper value
        //x[i] = 1;   // inital guess
    }
    for (int i = spline_dat.step_length; i < spline_dat.num_control_points+spline_dat.step_length; i++)
    {
        lb[i] = 0;   // could probably adjust this
        ub[i] = dist; // spam upper value
        //x[i] = spline_dat.dist/2;  // intial guess
    }
    x[spline_dat.step_length] = 0;  // adjust the first inital guess
    ub[spline_dat.step_length] = 0; // clamping the start
    lb[spline_dat.step_length] = 0; // clamping the start
    x[n-1] = spline_dat.dist; // adjust the last initial guess
    ub[n-1] = spline_dat.dist; // clamping the end
    lb[n-1] = spline_dat.dist; // clamping the end

    nlopt_opt opt;                  // optimization "object"

    opt = nlopt_create(NLOPT_LD_SLSQP, n); // select the algo and dimensionality
    int res;
    res = nlopt_set_lower_bounds(opt, lb);    // set the lower bounds
    printf("lower bounds result: %d \n", res);

    res = nlopt_set_upper_bounds(opt, ub);     // set the upper bounds
    printf("upper bounds result: %d \n", res);
    //res = nlopt_set_vector_storage(opt, m);
    //printf("set vector storage: %d\n", res );

    res = nlopt_set_min_objective(opt, objectiveFunc, &spline_dat); // set the objective function
    printf("set objective result: %d \n", res);

    res = nlopt_add_inequality_mconstraint(opt, m, derivConstraints, &spline_dat, tol);     // add the first inequality constraint
    printf("Inequality result: %d \n", res);

    res = nlopt_add_equality_mconstraint(opt, 4, clampConstraints, &spline_dat, tol2); // using m and tol I get a -2 here. Using 4 and tol2 I get a -4 during the opto
    printf("Equality result: %d\n", res);

    res = nlopt_set_xtol_rel(opt, 1e-4);      // set the stopping param as a relative dx value
    printf("xtol result: %d\n", res);

    res = nlopt_set_maxeval(opt, 100);
    printf("max eval result: %d\n", res);

    /*nlopt_opt local_opt = nlopt_create(NLOPT_LD_LBFGS, n);
    nlopt_set_xtol_rel(local_opt, 1e-4);
    nlopt_set_maxeval(local_opt, 200);
    res = nlopt_set_local_optimizer(opt, local_opt);*/


    double minf;                        // will hold the min value
    clock_t begin = clock();

    int code = nlopt_optimize(opt,x,&minf);     // optimize
    
    clock_t end = clock();

    double total_time = ((double)(end-begin))/ CLOCKS_PER_SEC;
    printf("Optimization execution time: %f seconds.\n", total_time);

    if (code < 0)
    {
        printf("nlopt failed!\n");
    }
    else
    {
    }

    printf("Return Code: %d \n", code);
    printf("found minimum after %d evaluations\n", iter_count);

    for (int i = 0; i < n; i++)
    {
        printf("%f ", x[i]);
    }
    printf("\n");

    printf("min val: %f \n", minf);

    int step_length = spline_dat.step_length; 
    double knot_steps[step_length];
    for (int i = 0; i < step_length; i++)
    {
        knot_steps[i] = x[i]; // get the knot steps from the param array
    }

    int control_length = spline_dat.num_control_points;
    double control_points[control_length];
    for (int i = 0; i < control_length; i++)
    {
        control_points[i] = x[i+spline_dat.step_length];   // get the control points from the param aray
    }

    int d = 1;
    double d1_control[spline_dat.num_control_points-d];
    double d2_control[spline_dat.num_control_points-d-1];
    double d3_control[spline_dat.num_control_points-d-2];
    for (int i = 0; i < n; i++)
    {
        d1_control[i] = controlPointDerivative(i, d, spline_dat.poly_deg, spline_dat.num_control_points-d, 0, knot_steps, control_points, step_length);
    }

    for (int i = 0; i < spline_dat.num_control_points-d-1; i++)
    {
        d2_control[i] = controlPointDerivative(i, d, spline_dat.poly_deg-1, spline_dat.num_control_points-d-1, 0, knot_steps, d1_control, step_length);
    }

    for (int i = 0; i < spline_dat.num_control_points-d-2; i++)
    {
        d3_control[i] = controlPointDerivative(i, d, spline_dat.poly_deg-2, spline_dat.num_control_points-d-2, 0, knot_steps, d2_control, step_length);
    }

    printf("d1: \n");
    for (int i = 0; i < spline_dat.num_control_points-d; i++)
    {
        printf("%f ", d1_control[i]);
    }
    printf("\n");

    printf("d2: \n");
    for (int i = 0; i < spline_dat.num_control_points-d-1; i++)
    {
        printf("%f ", d2_control[i]);
    }
    printf("\n");

    printf("d3: \n");
    for (int i = 0; i < spline_dat.num_control_points-d-2; i++)
    {
        printf("%f ", d3_control[i]);
    }
    printf("\n");    

    printf("Deriv constraint called %d times\n", d_con_count);
    printf("Clamp constraint called %d times\n", c_con_count);
    nlopt_destroy(opt);

    return 0;
}


// TODO: Try numeric gradient calcs
// TODO: Verify the gradient calcs
void derivConstraints(unsigned m, double* result, unsigned n, const double* x, double* grad, void* data)
{
    
    d_con_count++;
    bspline* spline_data = (bspline*) data;

    int step_length = spline_data->step_length; 
    double knot_steps[step_length];
    for (int i = 0; i < step_length; i++)
    {
        knot_steps[i] = x[i]; // get the knot steps from the param array
    }

    int control_length = spline_data->num_control_points;
    double control_points[control_length];
    for (int i = 0; i < control_length; i++)
    {
        control_points[i] = x[i+spline_data->step_length];   // get the control points from the param aray
    }
    
    double knots[step_length+(2*spline_data->poly_deg)+1];
    for (int i = 0; i < spline_data->poly_deg+1; i++)
    {
        knots[i] = 0;
    }
    for (int i = spline_data->poly_deg+1; i < step_length+spline_data->poly_deg+1; i++)
    {
        knots[i] = knots[i-1] + knot_steps[i-spline_data->poly_deg-1];
    }
    for (int i = step_length+spline_data->poly_deg+1; i < step_length+(2*spline_data->poly_deg)+1; i++)
    {
        knots[i] = knots[i-1];
    }

    double knots_v[step_length+(2*spline_data->poly_deg)-1];
    for (int i = 0; i < step_length+(2*spline_data->poly_deg)-1; i++)
    {
        knots_v[i] = knots[i+1];
    }

    double knots_a[step_length+(2*spline_data->poly_deg)-3]; // length is wrong
    for (int i = 0; i < step_length+(2*spline_data->poly_deg)-3; i++)
    {
        knots_a[i] = knots[i+2];
    }

    int n_c = spline_data->num_control_points; //- 1;
    int d = 1;

    double d1_control[n_c-d];
    double d2_control[n_c-d-1];
    double d3_control[n_c-d-2];
    for (int i = 0; i < n_c-d; i++)
    {
        // TODO: Fix nan issue
        d1_control[i] = controlPointDerivative(i, d, spline_data->poly_deg, n_c-d, 0, knot_steps, control_points, step_length); // this is giving me a nan at some point at i = 0
    }

    for (int i = 0; i < n_c-d-1; i++)
    {
        d2_control[i] = controlPointDerivative(i, d, spline_data->poly_deg-1, n_c-d-1, 0, knot_steps, d1_control, step_length);
    }

    for (int i = 0; i < n_c-d-2; i++)
    {
        d3_control[i] = controlPointDerivative(i, d, spline_data->poly_deg-2, n_c-d-2, 0, knot_steps, d2_control, step_length);
    }

    int temp = m;//-n_c+d;
    for (int i = 0; i < temp; i++)
    {
        if (i < n_c-d)
        {
            result[i] = (d1_control[i]*d1_control[i]) - (spline_data->v_max*spline_data->v_max);
        }
        else if(i-(n_c-d)<n_c-d-1)
        {
            result[i] = (d2_control[i-(n_c-d)]*d2_control[i-(n_c-d)]) - (spline_data->a_max*spline_data->a_max);
        }
        else
        {
            result[i] = (d3_control[i-2*(n_c)-(2*d)-1]*d3_control[i-2*(n_c)-(2*d)-1]) - (spline_data->j_max*spline_data->j_max);
        }
        if (isnan(result[i]))
        {
            printf("NaN: result[%d] \n", i);
        } 
    }

    /*for (int i = 0; i < n_c-d; i++)
    {
        result[temp+i] = -d1_control[i];
        if (isnan(result[temp+i]))
        {
            printf("NaN: result[%d] \n", temp+i);
        }
    }*/

    // numeric calc
    double dx[n];
    double control_f[n];
    double control_b[n];
    for (int i = 0; i < n; i++)
    {
        dx[i] = x[i] - prevX[i];
        control_f[i] = control_points[i] + dx[i];
        control_b[i] = control_points[i] - dx[i];
    }

    double d1_control_f[n_c-d];
    double d2_control_f[n_c-d-1];
    double d3_control_f[n_c-d-2];
    double d1_control_b[n_c-d];
    double d2_control_b[n_c-d-1];
    double d3_control_b[n_c-d-2];
    for (int i = 0; i < n_c-d; i++)
    {
        d1_control_f[i] = controlPointDerivative(i, d, spline_data->poly_deg, n_c-d, 0, knot_steps, control_f, step_length);
        d1_control_b[i] = controlPointDerivative(i, d, spline_data->poly_deg, n_c-d, 0, knot_steps, control_b, step_length);
    }

    for (int i = 0; i < n_c-d-1; i++)
    {
        d2_control_f[i] = controlPointDerivative(i, d, spline_data->poly_deg-1, n_c-d-1, 0, knot_steps, d1_control_f, step_length);
        d2_control_b[i] = controlPointDerivative(i, d, spline_data->poly_deg-1, n_c-d-1, 0, knot_steps, d1_control_b, step_length);
    }

    for (int i = 0; i < n_c-d-2; i++)
    {
        d3_control_f[i] = controlPointDerivative(i, d, spline_data->poly_deg-2, n_c-d-2, 0, knot_steps, d2_control_f, step_length);
        d3_control_b[i] = controlPointDerivative(i, d, spline_data->poly_deg-2, n_c-d-2, 0, knot_steps, d2_control_b, step_length);
    }   

    // un commenting this seems to increase the number of calls made in COBYLA
    if(grad)
    {
        int mod = 0;
        int offset = 0;
        for (int i = 0; i < m; i++)
        {
            double deltaT = 0;
            if (i < n_c-d)
            {
                mod = 0;
                offset = 0;
                deltaT = getdeltaT(i, spline_data->poly_deg-mod, knots);
                //double temp1 = getKnots(0, knot_steps, spline_data->poly_deg, i+spline_data->poly_deg+2, step_length);
                //double temp2 = getKnots(0, knot_steps, spline_data->poly_deg, i+1+1, step_length);

                //deltaT = temp1 - temp2;
            }
            else if(i < 2*(n_c)-(2*d)-2)
            {
                mod = 1;
                offset = n_c-d;
                deltaT = getdeltaT(i-offset, spline_data->poly_deg-mod, knots_v);

                //double temp1 = getKnots(0, knot_steps, spline_data->poly_deg-1, i+spline_data->poly_deg+1-offset, step_length);
                //double temp2 = getKnots(0, knot_steps, spline_data->poly_deg-1, i+1+1-offset, step_length);

                //deltaT = temp1 - temp2;
            }
            else
            {
                mod = 2;
                offset = 2*(n_c)-(2*d)-2;
                deltaT = getdeltaT(i-offset, spline_data->poly_deg-mod, knots_a);

                //double temp1 = getKnots(0, knot_steps, spline_data->poly_deg-2, i+spline_data->poly_deg-offset, step_length);
                //double temp2 = getKnots(0, knot_steps, spline_data->poly_deg-2, i+1+1-offset, step_length);

                //deltaT = temp1 - temp2;
            }
            for (int j = 0; j < n; j++)
            {
                if (j < spline_data->step_length)
                {
                    // Calculate the derivative of the constraint wrt the knot SPACE at x[j]
                    
                    switch(i-offset)
                    {
                        case 0:
                            if (j == 0)
                            {
                                grad[i*n + j] = -(spline_data->poly_deg-mod)/(deltaT* deltaT);   
                            }
                            else
                            {
                                grad[i*n + j] = 0;
                            }
                            break;

                        case 1:
                           if (j <= 1)
                            {
                                grad[i*n + j] = -(spline_data->poly_deg-mod)/(deltaT* deltaT);   
                            }
                            else
                            {
                                grad[i*n + j] = 0;
                            }
                            break; 

                        case 2:
                           if (j <= 2)
                            {
                                grad[i*n + j] = -(spline_data->poly_deg-mod)/(deltaT* deltaT);   
                            }
                            else
                            {
                                grad[i*n + j] = 0;
                            }
                            break; 

                        case 3:
                           if (j <= 3)
                            {
                                grad[i*n + j] = -(spline_data->poly_deg-mod)/(deltaT* deltaT);   
                            }
                            else
                            {
                                grad[i*n + j] = 0;
                            }
                            break; 
                        case 4:
                           if (j <= 4)
                            {
                                grad[i*n + j] = -(spline_data->poly_deg-mod)/(deltaT* deltaT);   
                            }
                            else
                            {
                                grad[i*n + j] = 0;
                            }
                            break; 
                        case 5:
                           if (j <= 5 && j > 0)
                            {
                                grad[i*n + j] = -(spline_data->poly_deg-mod)/(deltaT* deltaT);   
                            }
                            else
                            {
                                grad[i*n + j] = 0;
                            }
                            break;
                        case 6:
                           if (j <= 6 && j > 1)
                            {
                                grad[i*n + j] = -(spline_data->poly_deg-mod)/(deltaT* deltaT);   
                            }
                            else
                            {
                                grad[i*n + j] = 0;
                            }
                            break; 
                        case 7:
                           if (j <= 7 && j > 2)
                            {
                                grad[i*n + j] = -(spline_data->poly_deg-mod)/(deltaT* deltaT);   
                            }
                            else
                            {
                                grad[i*n + j] = 0;
                            }
                            break; 
                        case 8:
                           if (j <= 7 && j > 3)
                            {
                                grad[i*n + j] = -(spline_data->poly_deg-mod)/(deltaT* deltaT);   
                            }
                            else
                            {
                                grad[i*n + j] = 0;
                            }
                            break; 
                        case 9:
                           if (j <= 7 && j > 4)
                            {
                                grad[i*n + j] = -(spline_data->poly_deg-mod)/(deltaT* deltaT);   
                            }
                            else
                            {
                                grad[i*n + j] = 0;
                            }
                            break; 
                        case 10:
                        {
                            if (j <= 7 && j > 5)
                            {
                                grad[i*n + j] = -(spline_data->poly_deg-mod)/(deltaT* deltaT);   
                            }
                            else
                            {
                                grad[i*n + j] = 0;
                            }
                            break;
                        } 
                        case 11:
                        {
                            if (j <= 7 && j > 6)
                            {
                                grad[i*n + j] = -(spline_data->poly_deg-mod)/(deltaT* deltaT);   
                            }
                            else
                            {
                                grad[i*n + j] = 0;
                            }
                            break; 
                        }
                    }
                    if (grad[i*n+j] != 0)
                    {
                        grad[i*n + j] = 2*result[i]*grad[i*n + j]*(control_points[i+1-offset]+control_points[i-offset]);
                    }

                }
                else
                {
                    // Calculate the derivatie of the constraint wrt the control point at x[j]
                    if (j-spline_data->step_length == i-offset || j-spline_data->step_length == i+1-offset)
                    {
                        grad[i*n + j] = 2*result[i]*(spline_data->poly_deg-mod)/(deltaT);

                        int step_len = spline_data->step_length; // is this what I want to be using? Should it be spline_data->poly_deg?

                        if (j-spline_data->step_length < n_c-d && dx[j-step_len] != 0)
                        {
                            // Check the control_f and control_b arrays
                            // maybe check by hand
                            //grad[i*n+j] = (control_f[j-step_len] - control_b[j-step_len])/(dx[j-step_len]*2);
                            int xxxx = 0;
                        }                       
                    }
                    else
                    {
                        grad[i*n + j] = 0;
                    }
                }
                if (isinf(grad[i*n+j]))
                {
                    printf("inf deriv grad: grad[%d], deriv constraint counter: %d \n", i*n+j, d_con_count);
                    grad[i*n+j] = 1;
                }
                if (isnan(grad[i*n+j]))
                {
                    printf("NaN: Deriv grad[%d] \n", i*n+j);
                }
                if (grad[i*n+j] < 1e-10 && grad[i*n+j] > 0)
                {
                    printf("Very small deriv grad: grad[%d], deriv constraint counter: %d \n", i*n+j, d_con_count);
                    // grad[i*n+j] = 0; // gets here at i = 22, j = 0
                }
                // grad[i*n+j] = 1;
            }
        }
    }
}

void clampConstraints(unsigned m, double* result, unsigned n, const double* x, double* grad, void* data)
{
    // the problem relating to the -4 ret is in this function
    c_con_count++;
    //printf("here");
    /*for (int i = 0; i < m; i++)
    {
        result[i] = 1;
    }
    for (int i = 0; i < m*n; i++)
    {
        grad[i] = 1;
    }*/
    
    bspline* spline_data = (bspline*) data;

    int step_length = spline_data->step_length; 
    double knot_steps[step_length];
    for (int i = 0; i < step_length; i++)
    {
        knot_steps[i] = x[i]; // get the knot steps from the param array
    }

    int control_length = spline_data->num_control_points;
    double control_points[control_length];
    for (int i = 0; i < control_length; i++)
    {
        control_points[i] = x[i+spline_data->step_length];   // get the control points from the param aray
    }
    
    double knots[step_length+(2*spline_data->poly_deg)+1];
    for (int i = 0; i < spline_data->poly_deg+1; i++)
    {
        knots[i] = 0;
    }
    for (int i = spline_data->poly_deg+1; i < step_length+spline_data->poly_deg+1; i++)
    {
        knots[i] = knots[i-1] + knot_steps[i-spline_data->poly_deg-1];
    }
    for (int i = step_length+spline_data->poly_deg+1; i < step_length+(2*spline_data->poly_deg)+1; i++)
    {
        knots[i] = knots[i-1];
    }

    double knots_v[step_length+(2*spline_data->poly_deg)-1];
    for (int i = 0; i < step_length+(2*spline_data->poly_deg)-1; i++)
    {
        knots_v[i] = knots[i+1];
    }

    double knots_a[step_length+(2*spline_data->poly_deg)-3]; // length is wrong
    for (int i = 0; i < step_length+(2*spline_data->poly_deg)-3; i++)
    {
        knots_a[i] = knots[i+2];
    }

    int n_c = spline_data->num_control_points;
    int d = 1;

    double d1_control[n_c-d];
    double d2_control[n_c-d-1];
    double d3_control[n_c-d-2];
    for (int i = 0; i < n_c-d; i++)
    {
        d1_control[i] = controlPointDerivative(i, d, spline_data->poly_deg, n_c-d, 0, knot_steps, control_points, step_length);
    }

    for (int i = 0; i < n_c-d-1; i++)
    {
        d2_control[i] = controlPointDerivative(i, d, spline_data->poly_deg-1, n_c-d-1, 0, knot_steps, d1_control, step_length);
    }

    for (int i = 0; i < n_c-d-2; i++)
    {
        d3_control[i] = controlPointDerivative(i, d, spline_data->poly_deg-2, n_c-d-2, 0, knot_steps, d2_control, step_length);
    }

    result[0] = d1_control[0];
    result[1] = d1_control[n_c-d-1];

    result[2] = d2_control[0];
    result[3] = d2_control[9];
    for (int i = 0; i < 4; i++)
    {
        if (isnan(result[i]))
        {
            printf("NaN: result[%d] \n", i);
        }
    }

    if (grad)
    {
        int mod = 0;
        int offset = 0;
        int pos = 0;
        for (int i = 0; i < m; i++)
        {
            double deltaT = 0;
            if (i == 0)// || i == 1)
            {
                mod = 1;
                offset = 0;
                pos = 0;
                deltaT = getdeltaT(0, spline_data->poly_deg, knots_v);
            }
            else if(i == 1)
            {
                mod = 1;
                offset = 0;
                pos = n_c-d-1;
                deltaT = getdeltaT(pos, spline_data->poly_deg, knots_v);
            }
            else if(i == 2)
            {
                mod = 2;
                pos = 0;
                deltaT = getdeltaT(0, spline_data->poly_deg-1, knots_a);
            }
            else
            {
                mod = 2;
                pos = 8;
                deltaT = getdeltaT(pos, spline_data->poly_deg-1, knots_a);
            }
            
            for (int j = 0; j < n; j++)
            {
                if (j < spline_data->step_length)
                {
                    // Calculate the derivative of the constraint wrt the knot SPACE at x[j]
                    switch (i)
                    {
                    case 0:
                        if (j == 0)
                        {
                            grad[i*n + j] = -(spline_data->poly_deg-mod)/(deltaT* deltaT);   
                        }
                        else
                        {
                            grad[i*n + j] = 0;
                        }
                        break;
                    case 1:
                        if (j <= 7 && j > 6)
                        {
                            grad[i*n + j] = -(spline_data->poly_deg-mod)/(deltaT* deltaT);   
                        }
                        else
                        {
                            grad[i*n + j] = 0;
                        }
                        break; 
                    case 2:
                        if (j == 0)
                        {
                            grad[i*n + j] = -(spline_data->poly_deg-mod)/(deltaT* deltaT);   
                        }
                        else
                        {
                            grad[i*n + j] = 0;
                        }
                        break;
                    case 3:
                         if (j <= 7 && j > 6)
                        {
                            grad[i*n + j] = -(spline_data->poly_deg-mod)/(deltaT* deltaT);   
                        }
                        else
                        {
                            grad[i*n + j] = 0;
                        }
                        break; 
                    default:
                        break;
                    }
                    if (grad[i*n+j] != 0) // don't need this if statement
                    {
                        grad[i*n + j] = grad[i*n + j]*(control_points[i+1]+control_points[i]);
                    }
                }
                else
                {
                    // Calculate the derivatie of the constraint wrt the control point at x[j]
                    if (j-spline_data->step_length == pos || j-spline_data->step_length == pos+1)
                    {
                        grad[i*n + j] = (spline_data->poly_deg-mod)/(deltaT);
                    }
                    else
                    {
                        grad[i*n + j] = 0;
                    }
                }
                //grad[i*n+j] = 1;
                if (grad[i*n+j] < 1e-10 && grad[i*n+j] > 0)
                {
                    printf("Very small number clamp grad[%d]\n", i*n+j);
                    //grad[i*n+j] = 0;
                }
                if (isnan(grad[i*n+j]))
                {
                    printf("NaN: Clamp grad[%d] \n", i*n+j);
                }
                if (isinf(grad[i*n+j]))
                {
                    printf("Inf: Clamp grad[%d]\n", i*n+j);
                }
            }
        }
    }
    int xxx = 0;
    // RIGHT NOW EVERYTHING IS A 0 - WHY??
}

double objectiveFunc(unsigned n, const double* x, double* grad, void* data)
{
    iter_count++;

    bspline* spline_data = (bspline*) data;

    double sum = 0;
    for (int i = 0; i < spline_data->step_length; i++)
    {
        sum += x[i];
    }

    if (grad)
    {
        for (int i = 0; i < n; i++)
        {
            if (i < spline_data->step_length)
            {
                grad[i] = 1;
            }
            else
            {
                grad[i] = 0;
            }
            
        }
    }

    return sum;

}

double getKnots(double start_knot, double* knot_steps, int poly_deg, int index, int step_length)
{
    if (index <= poly_deg)
    {
        return start_knot;
    }
    else if (index <= poly_deg+step_length) // removed the equals
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

double getdeltaT(int i, int poly_deg, double* knots)
{
    return knots[i+poly_deg+1] - knots[i+1];
}