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
    spline_dat.dist = 10;

    int n = spline_dat.step_length + spline_dat.num_control_points; // dimensionality of the problem

    unsigned int m = 27; // double check this. it is the number of constraints
    double tol[m];
    for (int i = 0; i < m; i++)
    {
        tol[i] = 1e-8;
    }

    double lb[n];      // lower boungs on the optimization params. maybe needs to be const
    double ub[n];
    double x[19] = {3, 3, 3, 3, 3, 3, 3, 0, 0.5, 1, 1.5, 2.75, 3.75, 4.25, 5.25, 6.75, 8.5, 9, 10};
       // starting guess
    // TODO: Get a better initial guess

    for (int i = 0; i < spline_dat.step_length; i++)
    {
        lb[i] = 0;  // knot step lower bounds
        ub[i] = 200; // spam upper value
        //x[i] = 1;   // inital guess
    }
    for (int i = spline_dat.step_length; i < spline_dat.num_control_points+spline_dat.step_length; i++)
    {
        lb[i] = -200;   // could probably adjust this
        ub[i] = 200; // spam upper value
        //x[i] = spline_dat.dist/2;  // intial guess
    }
    x[spline_dat.step_length] = 0;  // adjust the first inital guess
    ub[spline_dat.step_length] = 0; // clamping the start
    lb[spline_dat.step_length] = 0; // clamping the start
    x[n-1] = spline_dat.dist; // adjust the last initial guess
    ub[n-1] = spline_dat.dist; // clamping the end
    lb[n-1] = spline_dat.dist; // clamping the end

    nlopt_opt opt;                  // optimization "object"

    opt = nlopt_create(NLOPT_LN_COBYLA, n);//NLOPT_LD_SLSQP, n); // select the algo and dimensionality
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

    res = nlopt_add_equality_mconstraint(opt, 4, clampConstraints, &spline_dat, tol);
    printf("Equality result: %d\n", res);

    res = nlopt_set_xtol_rel(opt, 1e-4);      // set the stopping param as a relative dx value
    printf("xtol result: %d\n", res);

    res = nlopt_set_maxeval(opt, 20000);
    printf("max eval result: %d\n", res);

/*    
    nlopt_opt local_opt = nlopt_create(NLOPT_LD_LBFGS, n);
    nlopt_set_xtol_rel(local_opt, 1e-4);
    nlopt_set_maxeval(local_opt, 200);
    res = nlopt_set_local_optimizer(opt, local_opt);
*/

    double minf;                        // will hold the min value
    clock_t begin = clock();

    int code = nlopt_optimize(opt,x,&minf);     // optimize
    
    clock_t end = clock();

    double total_time = ((double)(end-begin))/ CLOCKS_PER_SEC;
    printf("Optimization execution time: %f seconds.\n", total_time);

    if (code < 0)
    {
        printf("nlopt failed!\n");
        printf("Error Code: %d \n", code);
    }
    else
    {
    }

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


// TODO: Verify control point calcs wuth matlab
// TODO: Try numeric gradient calcs
// TODO: Verify the gradient calcs
// TODO: Why is this only getting called once?
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
    /*
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
    for (int i = 0; i < step_length+(2*spline_data->poly_deg); i++)
    {
        knots_v[i] = knots[i+1];
    }

    double knots_a[step_length+(2*spline_data->poly_deg)-2];
    for (int i = 0; i < step_length+(2*spline_data->poly_deg)-1; i++)
    {
        knots_a[i] = knots[i+2];
    }*/

    int n_c = spline_data->num_control_points; //- 1;
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
   /* double deltaT1 = 0;
    for (int i = 0; i < n_c-d; i++)
    {
        deltaT1 = getdeltaT(i, spline_data->poly_deg, knots);
        d1_control[i] = (control_points[i+1] - control_points[i])*spline_data->poly_deg/(deltaT1);
    }

    double deltaT2 = 0;
    for (int i = 0; i < n_c-d-1; i++)
    {
        deltaT2 = getdeltaT(i, spline_data->poly_deg-1, knots_v);
        d2_control[i] = (d1_control[i+1] - d1_control[i])*(spline_data->poly_deg-1)/(deltaT2);
    }

    double deltaT3 = 0;
    for (int i = 0; i < n_c-d-2; i++)
    {
        deltaT3 = getdeltaT(i, spline_data->poly_deg-2, knots_a);

        d3_control[i] = (d2_control[i+1] - d2_control[i])*(spline_data->poly_deg-2)/(deltaT3);
    }*/

    for (int i = 0; i < m; i++)
    {
        // TODO: Check the results are being indexed correctly and are correct
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
    }

/*
    if(grad)
    {
        int mod = 0;
        int offset = 0;
        for (int i = 0; i < m; i++)
        {
            double deltaT = 0;
            if (i < n_c-d)
            {
                deltaT = getdeltaT(i, spline_data->poly_deg, knots);
            }
            else if(i < 2*(n_c)-(2*d)-1)
            {
                mod = 1;
                offset = n_c-d;
                deltaT = getdeltaT(i-offset, spline_data->poly_deg-mod, knots_v);
            }
            else
            {
                mod = 2;
                offset = 2*(n_c)-(2*d)-1;
                deltaT = getdeltaT(i-offset, spline_data->poly_deg-mod, knots_a);
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
                        grad[i*n + j] = grad[i*n + j]*(control_points[i+1]+control_points[i]);
                    }

                }
                else
                {
                    // Calculate the derivatie of the constraint wrt the control point at x[j]
                    if (j-spline_data->step_length == i || j-spline_data->step_length == i+1)
                    {
                        grad[i*n + j] = (spline_data->poly_deg-mod)/(deltaT);
                    }
                    else
                    {
                        grad[i*n + j] = 0;
                    }
                }
            }
        }
    }*/
}

void clampConstraints(unsigned m, double* result, unsigned n, const double* x, double* grad, void* data)
{
    c_con_count++;
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
    for (int i = 0; i < step_length+(2*spline_data->poly_deg); i++)
    {
        knots_v[i] = knots[i+1];
    }

    double knots_a[step_length+(2*spline_data->poly_deg)-2];
    for (int i = 0; i < step_length+(2*spline_data->poly_deg)-1; i++)
    {
        knots_a[i] = knots[i+2];
    }

    // TODO: Add in the gradient of the contraint (which should be doable but using the derivatives I already need to compute)

    int n_c = spline_data->num_control_points - 1;
    int d = 1;

 //   result[0] = x[step_length];
  //  result[1] = x[control_length + step_length] - spline_data->dist;
    
    double control = controlPointDerivative(0, d, spline_data->poly_deg, n_c-d, spline_data->start_knot, knot_steps, control_points, step_length);
    result[0] = control;
    control = controlPointDerivative(n_c-d, d, spline_data->poly_deg, n_c-d, spline_data->start_knot, knot_steps, control_points, step_length);
    result[1] = control;

    control = controlPointDerivative(0, d, spline_data->poly_deg-1, n_c-d-1, spline_data->start_knot, knot_steps, control_points, step_length);
    result[2] = control;
    control = controlPointDerivative(n_c-d-1, d, spline_data->poly_deg-1, n_c-d-1, spline_data->start_knot, knot_steps, control_points, step_length);
    result[3] = control;
    /*
    double deltaT1;

    deltaT1 = getdeltaT(0, spline_data->poly_deg-1, knots_v);
    result[0] = (control_points[1] - control_points[0])*spline_data->poly_deg/(deltaT1);
    deltaT1 = getdeltaT(1, spline_data->poly_deg-1, knots_v);
    double c1 = (control_points[2] - control_points[1])*spline_data->poly_deg/(deltaT1);

    deltaT1 = getdeltaT(n_c-d, spline_data->poly_deg-1, knots_v);
    result[1] = (control_points[n_c-d+1] - control_points[n_c-d])*spline_data->poly_deg/(deltaT1);
    deltaT1 = getdeltaT(1, spline_data->poly_deg-1, knots_v);
    double c2 = (control_points[n_c-d] - control_points[n_c-d-1])*spline_data->poly_deg/(deltaT1);

    deltaT1 = getdeltaT(0, spline_data->poly_deg-2, knots_a);
    result[2] = (c1 - result[0])*(spline_data->poly_deg-1)/(deltaT1);

    deltaT1 = getdeltaT(0, spline_data->poly_deg-2, knots_a);
    result[3] = (result[1] - c2)*(spline_data->poly_deg-1)/(deltaT1);

    if (grad)
    {
        int mod = 0;
        int offset = 0;
        for (int i = 0; i < m; i++)
        {
            double deltaT = 0;
            if (i == 0 || i == 1)
            {
                mod = 1;
                deltaT = getdeltaT(i, spline_data->poly_deg-mod, knots_v);
            }
            else
            {
                mod = 1;
                offset = 2;
                deltaT = getdeltaT(i-offset, spline_data->poly_deg-mod, knots_a);
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
                           if (j <= 7 && j > 5)
                            {
                                grad[i*n + j] = -(spline_data->poly_deg-mod)/(deltaT* deltaT);   
                            }
                            else
                            {
                                grad[i*n + j] = 0;
                            }
                            break; 
                        case 11:
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
                    if (grad[i*n+j] != 0)
                    {
                        grad[i*n + j] = grad[i*n + j]*(control_points[i+1]+control_points[i]);
                    }
                }
                else
                {
                    // Calculate the derivatie of the constraint wrt the control point at x[j]
                    if (j-spline_data->step_length == i || j-spline_data->step_length == i+1)
                    {
                        grad[i*n + j] = (spline_data->poly_deg-mod)/(deltaT);
                    }
                    else
                    {
                        grad[i*n + j] = 0;
                    }
                }
            }
        }
    }
*/
/*
    control = controlPointDerivative(0, d, spline_data->poly_deg, n, spline_data->start_knot, knot_steps, control_points, step_length);
    result[4] = control;
    control = controlPointDerivative(n-d-1, d, spline_data->poly_deg, n, spline_data->start_knot, knot_steps, control_points, step_length);
    result[5] = control;
*/

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

    // TODO: Add in gradient

    return sum;

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

// might need some modification to work on a 0 indexed system
double controlPointDerivative(int i, int derivative, int poly_deg, int n, double start_knot, double* knot_steps, double* control_points, int step_length)
{
    if (derivative == 0)
    {
        return control_points[i];
    }
    else
    {
        double temp1 = getKnots(start_knot, knot_steps, poly_deg, i+poly_deg+2, step_length); // might need to drop the +1
        double temp2 = getKnots(start_knot, knot_steps, poly_deg, i+derivative+1, step_length);

        double denom = temp1 - temp2;

        return ((poly_deg - derivative + 1)/denom) * (controlPointDerivative(i+1, derivative-1, poly_deg, n, start_knot, knot_steps, control_points, step_length) - controlPointDerivative(i, derivative-1, poly_deg, n, start_knot, knot_steps, control_points, step_length));
    }
    
}

double getdeltaT(int i, int poly_deg, double* knots)
{
    return knots[i+poly_deg+1] - knots[i+1];
}