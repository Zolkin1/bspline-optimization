/*
This will use the nonlinear optimization tools in NLOPT to optimize a BSpline based trajectory and time the algorithm.
*/

/*
Notes:
- The initial condition can only work for >= distances. i.e. the one below only works for >= 15m distances

Overall todos:
- Make a few different intial condition (does one really not aggressive one work?) and test their convergence speeds
- Debug the issue with the initial acceleration not being constrained
- Test in simulation with a non-clamped acceleration
- Run on a PI and/or the STM32
*/

#include <nlopt.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>

#include "optimizeSpline.h"

int iter_count = 0;
int d_con_count = 0;
int c_con_count = 0;
int c1_con_count = 0;
int c2_con_count = 0;
int c3_con_count = 0;
int c4_con_count = 0;

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
    spline_dat.dist = 18;

    int n = spline_dat.step_length + spline_dat.num_control_points; // dimensionality of the problem

    unsigned m = 30; // it is the number of constraints
    double tol[m];
    double tol2[4] = {1e-8, 1e-8, 1e-8, 1e-8};

    for (int i = 0; i < m; i++)
    {
        tol[i] = 1e-8;
        //tol2[i] = 1e-8;
    }

    int dist = spline_dat.dist;
    double lb[n];      // lower boungs on the optimization params. maybe needs to be const
    double ub[n];
    
    // initial conditions can only be used on larger distances. So this condition is only valid for >= 15m
    double x[19] = {.4, 5, 5, 5, 5, 5, .4, 0, 0, 0, 3*dist/10, 4*dist/10, 5*dist/10, 6*dist/10, 7*dist/10, 8*dist/10, dist, dist, dist};

    for (int i = 0; i < 19; i++)
    {
        prevX[i] = x[i];
    }

    for (int i = 0; i < spline_dat.step_length; i++)
    {
        lb[i] = 1e-6;  // knot step lower bounds
        ub[i] = 20000; // spam upper value
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

    // putting in the upper bounds seems to prevent it from returing the correct answer even though it finds the correct answer
    res = nlopt_set_upper_bounds(opt, ub);     // set the upper bounds
    printf("upper bounds result: %d \n", res);

    res = nlopt_set_min_objective(opt, objectiveFunc, &spline_dat); // set the objective function
    printf("set objective result: %d \n", res);

    res = nlopt_add_inequality_mconstraint(opt, m, derivConstraints, &spline_dat, tol);     // add the first inequality constraint
    printf("Inequality result: %d \n", res);

    //res = nlopt_add_equality_mconstraint(opt, 4, clampConstraints, &spline_dat, tol2); // using m and tol I get a -2 here. Using 4 and tol2 I get a -4 during the opto
    //printf("Equality result: %d\n", res);

    res = nlopt_add_equality_constraint(opt, c1_constraint, &spline_dat, 1e-8);
    res = nlopt_add_equality_constraint(opt, c2_constraint, &spline_dat, 1e-8);
    //res = nlopt_add_equality_constraint(opt, c3_constraint, &spline_dat, 1e-8);   // not allowing to get to max acc for some reason
    res = nlopt_add_equality_constraint(opt, c4_constraint, &spline_dat, 1e-8);


    res = nlopt_set_xtol_rel(opt, 1e-4);      // set the stopping param as a relative dx value
    printf("xtol result: %d\n", res);

    //res = nlopt_set_maxeval(opt, 10);
    //printf("max eval result: %d\n", res);

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

    printf("Deriv constraint called %d times\n", d_con_count);
    printf("Clamp constraint called %d times\n", c_con_count);
    printf("C1 constraint called %d times \n", c1_con_count);
    printf("C2 constraint called %d times \n", c2_con_count);
    printf("C3 constraint called %d times \n", c3_con_count);
    printf("C4 constraint called %d times \n", c4_con_count);

    printf("real final x: \n");
    for (int i = 0; i < n; i++)
    {
        printf("%f ", spline_dat.x_l[i]);
    }
    printf("\n");
    nlopt_destroy(opt);

    return 0;
}

void derivConstraints(unsigned m, double* result, unsigned n, const double* x, double* grad, void* data)
{
    d_con_count++;
    bspline* spline_data = (bspline*) data;

    for (int i = 0; i < n; i++)
    {
        spline_data->x_l[i] = x[i];
    }
    int n_c = spline_data->num_control_points;
    int step_length = spline_data->step_length;
    int k = spline_data->poly_deg;
    int d = 1;
    
    //printf("Knot Steps: \n");
    double knot_steps[step_length];
    for (int i = 0; i < step_length; i++)
    {
        knot_steps[i] = x[i];
        //printf("%f ", x[i]);
    }

    //printf("\nControl Points: \n");
    double control[n_c];
    for (int i = 0; i < n_c; i++)
    {
        control[i] = x[i+step_length];
        //printf("%f ", x[i+step_length]);
    }
    //printf("\n \n");
    // Recover the knots from the knot steps
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

    // Calculate the knots for the derivatives. They contain the same knot steps but less multiplicity on the ends
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

    double d1_control[n_c-d];
    double d2_control[n_c-d-1];
    double d3_control[n_c-d-2];

    for (int i = 0; i < n_c-d; i++)
    {
        d1_control[i] = (control[i+1]-control[i])*k/getdeltaT(i, k, knots);
        result[i] = (d1_control[i] * d1_control[i]) - (spline_data->v_max * spline_data->v_max); // constraint
    }

    for (int i = 0; i < n_c-d-1; i++)
    {
        d2_control[i] = (d1_control[i+1]-d1_control[i])*(k-1)/getdeltaT(i, k-1, knots_v);
        result[i+n_c-d] = (d2_control[i] * d2_control[i]) - (spline_data->a_max * spline_data->a_max);
    }

    for (int i = 0; i < n_c-d-2; i++)
    {
        d3_control[i] = (d2_control[i+1]-d2_control[i])*(k-2)/getdeltaT(i, k-2, knots_a);
        result[i + n_c-d + n_c-d-1] = (d3_control[i] * d3_control[i]) - (spline_data->j_max * spline_data->j_max);
    }

    /*printf("Result: \n");
    for (int i = 0; i < m; i++)
    {
        printf("%f ", result[i]);
    }
    printf("\n \n");*/
    
    // This is not used eevery time this constraint is called
    if (grad)
    {
        int opto_vars = n;
        double grad1[(n_c-d)*opto_vars];
        int start = 0;
        //printf("i | j |     dt    |   grad   | \n");
        for (int i = 0; i < n_c-d; i++) // go through the first n constraints
        {
            for(int j = 0; j < opto_vars; j++)  // differentiate wrt to each of the opto variables
            {
                double dt = getdeltaT(i,k,knots);
                if (j < step_length)  // if we are taking wrt to a knot step then use this eq. BUT also need to make sure that knot step is used
                {
                    if (j <= i && j >= i-k+1)   // chect that that knot step is used
                    {
                        grad1[i*opto_vars+j] = -2*d1_control[i]*k*(control[i+1]-control[i])/(dt*dt);
                    }
                    else
                    {
                        grad1[i*opto_vars+j] = 0;
                    }
                }
                else    // differentiate wrt to the control points
                {
                    if (j-step_length == i)
                    {
                        grad1[i*opto_vars+j] = -2*d1_control[i]*k/dt;
                    }
                    else if(j-step_length == i+1)
                    {
                        grad1[i*opto_vars+j] = 2*d1_control[i+1]*k/dt;
                    }
                    else
                    {
                        grad1[i*opto_vars+j] = 0;
                    }
                }
                //printf("%d | %d | %f | %f \n", i, j, dt, grad[i*opto_vars+j]);
                //printf("%f ", grad1[i*opto_vars + j]);
            }
            //printf("\n");
        }

        //printf("\n");
        int k2 = k - 1;
        double grad2[(n_c-d-1)*opto_vars];
    //    start = step_length*opto_vars;
        for (int i = 0; i < n_c-d-1; i++) // go through the second n-1 constraints
        {
            for(int j = 0; j < opto_vars; j++)  // differentiate wrt to each of the opto variables
            {
                double dt = getdeltaT(i, k2, knots_v);
                if (j < step_length)  // if we are taking wrt to a knot step then use this eq. BUT also need to make sure that knot step is used
                {
                    if (j <= i && j >= i-k2+1)   // chect that that knot step is used
                    {
                        grad2[(i*opto_vars)+j] = -2*d2_control[i]*k2*(d1_control[i+1]-d1_control[i])/(dt*dt);
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
                        grad2[(i*opto_vars)+j] = -2*d2_control[i]*k2/dt;
                    }
                    else if (j-step_length == i+1)
                    {
                        grad2[(i*opto_vars)+j] = 2*d2_control[i+1]*k2/dt;
                    }
                    
                    else
                    {
                        grad2[i*opto_vars+j] = 0;
                    }
                    
                }
                //printf("%f ", grad2[(i*opto_vars)+j]);
            }
            //printf("\n");
        }

        //printf("\n");
        int k3 = k - 2;
        double grad3[(n_c-d-2)*opto_vars];
        for (int i = 0; i < n_c-d-2; i++) // go through the third n-2 constraints
        {
            for(int j = 0; j < opto_vars; j++)  // differentiate wrt to each of the opto variables
            {
                double dt = getdeltaT(i,k3,knots_a);
                if (j < step_length)  // if we are taking wrt to a knot step then use this eq. BUT also need to make sure that knot step is used
                {
                    if (j <= i && j >= i-k3+1)   // chect that that knot step is used
                    {
                        grad3[(i*opto_vars)+j] = -2*d3_control[i]*k3*(d2_control[i+1]-d2_control[i])/(dt * dt);
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
                        grad3[(i*opto_vars)+j] = -2*d2_control[i]*k3/dt;
                    }
                    else if(j-step_length == i+1)
                    {
                        grad3[(i*opto_vars)+j] = 2*d2_control[i+1]*k3/dt;
                    }
                    else
                    {
                        grad3[i*opto_vars+j] = 0;
                    }
                    
                }
                //printf("%f ", grad3[(i*opto_vars)+j]);
            }
            //printf("\n");
        }

        // Check to be sure that this is being concatenated in the correct order
        // make sure this is correct
        for (int i = 0; i < n*m; i++)
        {
            if (i < (n_c-d)*opto_vars)
            {
                grad[i] = grad1[i]; 
            }
            else if (i < ((n_c-d-1)*opto_vars) + (n_c-d)*opto_vars)
            {
                grad[i] = grad2[i-((n_c-d)*opto_vars)];
            }
            else
            {
                grad[i] = grad3[i - ((n_c-d-1)*opto_vars) - (n_c-d)*opto_vars];
            }
        }
    }
    
    /*printf("\n");
    printf("d1: \n");
    for (int i = 0; i < n_c-d; i++)
    {
        printf("%f ", d1_control[i]);
    }

    printf("\n");
    printf("\n");

    printf("d2: \n");
    for (int i = 0; i < n_c-d-1; i++)
    {
        printf("%f ", d2_control[i]);
    }
    printf("\n");
    printf("\n");

    printf("d3: \n");
    for (int i = 0; i < n_c-d-2; i++)
    {
        printf("%f ", d3_control[i]);
    }
    printf("\n");
    printf("\n"); */
}

void clampConstraints(unsigned m, double* result, unsigned n, const double* x, double* grad, void* data)
{
    c_con_count++;
    bspline* spline_data = (bspline*) data;

    int n_c = spline_data->num_control_points;
    int step_length = spline_data->step_length;
    int k = spline_data->poly_deg;
    int d = 1;
    
    printf("Knot Steps Clamp: \n");
    double knot_steps[step_length];
    for (int i = 0; i < step_length; i++)
    {
        knot_steps[i] = x[i];
        printf("%f ", x[i]);
    }

    printf("\nControl Points Clamp: \n");
    double control[n_c];
    for (int i = 0; i < n_c; i++)
    {
        control[i] = x[i+step_length];
        printf("%f ", x[i+step_length]);
    }
    printf("\n \n");
    // Recover the knots from the knot steps
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

    // Calculate the knots for the derivatives. They contain the same knot steps but less multiplicity on the ends
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

    double d1_control[n_c-d];
    double d2_control[n_c-d-1];
    double d3_control[n_c-d-2];

    for (int i = 0; i < n_c-d; i++)
    {
        d1_control[i] = (control[i+1]-control[i])*k/getdeltaT(i, k, knots);
        //result[i] = (d1_control[i] * d1_control[i]) - (spline_data->v_max * spline_data->v_max); // constraint
    }

    for (int i = 0; i < n_c-d-1; i++)
    {
        d2_control[i] = (d1_control[i+1]-d1_control[i])*(k-1)/getdeltaT(i, k-1, knots_v);
        //result[i+n_c-d] = (d2_control[i] * d2_control[i]) - (spline_data->a_max * spline_data->a_max);
    }

    for (int i = 0; i < n_c-d-2; i++)
    {
        d3_control[i] = (d2_control[i+1]-d2_control[i])*(k-2)/getdeltaT(i, k-2, knots_a);
        //result[i + n_c-d + n_c-d-1] = (d3_control[i] * d3_control[i]) - (spline_data->j_max * spline_data->j_max);
    }

    result[0] = d1_control[0];
    result[1] = d1_control[n_c-d-1];

    result[2] = d2_control[0];
    result[3] = d2_control[n_c-d-2];

    /*for (int i = 0; i < m; i++)
    {
        result[i] = 0;
    }*/

    printf("Clamp Result: \n");
    for (int i = 0; i < m; i++)
    {
        printf("%f ", result[i]);
    }
    printf("\n \n");

    if (grad)
    {
        int m2 = 4;
        int opto_vars = n;
        int i = 1;
        for (int h = 0; h < 2; h++)
        {
            if (h == 0)
            {
                i = 0;
            }
            else
            {
                i = n_c-d-1;
            }
            
            for(int j = 0; j < opto_vars; j++)  // differentiate wrt to each of the opto variables
            {
                double dt = getdeltaT(i,k,knots);
                if (j < step_length)  // if we are taking wrt to a knot step then use this eq. BUT also need to make sure that knot step is used
                {
                    if (j <= i && j >= i-k+1)   // chect that that knot step is used
                    {
                        grad[h*opto_vars + j] = -k*(control[i+1] - control[i])/(dt*dt);
                    }
                    else
                    {
                        grad[h*opto_vars + j] = 0;
                    }
                }
                else    // differentiate wrt to the control points
                {
                    if (j-step_length == i)
                    {
                        grad[h*opto_vars+j] = -k/dt;
                    }
                    else if(j-step_length == i+1)
                    {
                        grad[h*opto_vars+j] = k/dt;
                    }
                    else
                    {
                        grad[h*opto_vars+j] = 0;
                    }
                }
                //printf("%d | %d | %f | %f \n", i, j, dt, grad[i*opto_vars+j]);
                //printf("%f ", grad[h*opto_vars + j]);
            }
            //printf("\n");
        }
        int k2 = k - 1;
        for (int h = 2; h < 4; h++)
        {
            if (h == 2)
            {
                i = 1;
            }
            else
            {
                i = n_c-d-2;
            }
            
            for(int j = 0; j < opto_vars; j++)  // differentiate wrt to each of the opto variables
            {
                double dt = getdeltaT(i, k2, knots_v);
                if (j < step_length)  // if we are taking wrt to a knot step then use this eq. BUT also need to make sure that knot step is used
                {
                    if (j <= i && j >= i-k2+1)   // chect that that knot step is used
                    {
                        grad[(h*opto_vars)+j] = -k2*(d1_control[i+1] - d1_control[i])/(dt*dt);
                    }
                    else
                    {
                        grad[h*opto_vars+j] = 0; 
                    }
                }
                else    // differentiate wrt to the control points
                {
                    if (j-step_length == i)
                    {
                        grad[(h*opto_vars)+j] = -k2/dt;
                    }
                    else if (j-step_length == i+1)
                    {
                        grad[(h*opto_vars)+j] = k2/dt;
                    }
                    
                    else
                    {
                        grad[h*opto_vars+j] = 0;
                    }
                    
                }
                printf("%f ", grad[(h*opto_vars)+j]);
            }
            printf("\n");
        }

        /*for (int i = 0; i < m*n; i++)
        {
            grad[i] = 0;
        }*/
    }
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

double getdeltaT(int i, int poly_deg, double* knots)
{
    return knots[i+poly_deg+1] - knots[i+1];
}

double c1_constraint(unsigned n, const double* x, double * grad, void* data)
{
    c1_con_count++;
    bspline* spline_data = (bspline*) data;

    int n_c = spline_data->num_control_points;
    int step_length = spline_data->step_length;
    int k = spline_data->poly_deg;
    int d = 1;
    
    //printf("Knot Steps Clamp: \n");
    double knot_steps[step_length];
    for (int i = 0; i < step_length; i++)
    {
        knot_steps[i] = x[i];
      //  printf("%f ", x[i]);
    }

    //printf("\nControl Points Clamp: \n");
    double control[n_c];
    for (int i = 0; i < n_c; i++)
    {
        control[i] = x[i+step_length];
        //printf("%f ", x[i+step_length]);
    }
    //printf("\n \n");
    // Recover the knots from the knot steps
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

    // Calculate the knots for the derivatives. They contain the same knot steps but less multiplicity on the ends
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

    double d1_control[n_c-d];
    double d2_control[n_c-d-1];
    double d3_control[n_c-d-2];

    for (int i = 0; i < n_c-d; i++)
    {
        d1_control[i] = (control[i+1]-control[i])*k/getdeltaT(i, k, knots);
        //result[i] = (d1_control[i] * d1_control[i]) - (spline_data->v_max * spline_data->v_max); // constraint
    }

    if (grad)
    {
        int m2 = 4;
        int opto_vars = n;
        int i = 0;
        int h = 0;
    
        for(int j = 0; j < opto_vars; j++)  // differentiate wrt to each of the opto variables
        {
            double dt = getdeltaT(i,k,knots);
            if (j < step_length)  // if we are taking wrt to a knot step then use this eq. BUT also need to make sure that knot step is used
            {
                if (j <= i && j >= i-k+1)   // chect that that knot step is used
                {
                    grad[h*opto_vars + j] = -k*(control[i+1] - control[i])/(dt*dt);
                }
                else
                {
                    grad[h*opto_vars + j] = 0;
                }
            }
            else    // differentiate wrt to the control points
            {
                if (j-step_length == i)
                {
                    grad[h*opto_vars+j] = -k/dt;
                }
                else if(j-step_length == i+1)
                {
                    grad[h*opto_vars+j] = k/dt;
                }
                else
                {
                    grad[h*opto_vars+j] = 0;
                }
            }
            //printf("%d | %d | %f | %f \n", i, j, dt, grad[i*opto_vars+j]);
            //printf("%f ", grad[h*opto_vars + j]);
        }
    }
    return d1_control[0];
}

double c2_constraint(unsigned n, const double* x, double * grad, void* data)
{
    c2_con_count++;
    bspline* spline_data = (bspline*) data;

    int n_c = spline_data->num_control_points;
    int step_length = spline_data->step_length;
    int k = spline_data->poly_deg;
    int d = 1;
    
    //printf("Knot Steps Clamp: \n");
    double knot_steps[step_length];
    for (int i = 0; i < step_length; i++)
    {
        knot_steps[i] = x[i];
        //printf("%f ", x[i]);
    }

    //printf("\nControl Points Clamp: \n");
    double control[n_c];
    for (int i = 0; i < n_c; i++)
    {
        control[i] = x[i+step_length];
        //printf("%f ", x[i+step_length]);
    }
    //printf("\n \n");
    // Recover the knots from the knot steps
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

    // Calculate the knots for the derivatives. They contain the same knot steps but less multiplicity on the ends
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

    double d1_control[n_c-d];
    double d2_control[n_c-d-1];
    double d3_control[n_c-d-2];

    for (int i = 0; i < n_c-d; i++)
    {
        d1_control[i] = (control[i+1]-control[i])*k/getdeltaT(i, k, knots);
        //result[i] = (d1_control[i] * d1_control[i]) - (spline_data->v_max * spline_data->v_max); // constraint
    }

    if (grad)
    {
        int m2 = 4;
        int opto_vars = n;
        int i = n_c-d-1;
        int h = 0;
    
        for(int j = 0; j < opto_vars; j++)  // differentiate wrt to each of the opto variables
        {
            double dt = getdeltaT(i,k,knots);
            if (j < step_length)  // if we are taking wrt to a knot step then use this eq. BUT also need to make sure that knot step is used
            {
                if (j <= i && j >= i-k+1)   // chect that that knot step is used
                {
                    grad[h*opto_vars + j] = -k*(control[i+1] - control[i])/(dt*dt);
                }
                else
                {
                    grad[h*opto_vars + j] = 0;
                }
            }
            else    // differentiate wrt to the control points
            {
                if (j-step_length == i)
                {
                    grad[h*opto_vars+j] = -k/dt;
                }
                else if(j-step_length == i+1)
                {
                    grad[h*opto_vars+j] = k/dt;
                }
                else
                {
                    grad[h*opto_vars+j] = 0;
                }
            }
            //printf("%d | %d | %f | %f \n", i, j, dt, grad[i*opto_vars+j]);
            //printf("%f ", grad[h*opto_vars + j]);
        }
    }
    return d1_control[n_c-d-1];
}

double c3_constraint(unsigned n, const double* x, double * grad, void* data)
{
    c3_con_count++;
    bspline* spline_data = (bspline*) data;

    int n_c = spline_data->num_control_points;
    int step_length = spline_data->step_length;
    int k = spline_data->poly_deg;
    int d = 1;
    
    //printf("Knot Steps Clamp: \n");
    double knot_steps[step_length];
    for (int i = 0; i < step_length; i++)
    {
        knot_steps[i] = x[i];
        //printf("%f ", x[i]);
    }

    //printf("\nControl Points Clamp: \n");
    double control[n_c];
    for (int i = 0; i < n_c; i++)
    {
        control[i] = x[i+step_length];
        //printf("%f ", x[i+step_length]);
    }
    //printf("\n \n");
    // Recover the knots from the knot steps
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

    // Calculate the knots for the derivatives. They contain the same knot steps but less multiplicity on the ends
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

    double d1_control[n_c-d];
    double d2_control[n_c-d-1];

    for (int i = 0; i < n_c-d; i++)
    {
        d1_control[i] = (control[i+1]-control[i])*k/getdeltaT(i, k, knots);
        //result[i] = (d1_control[i] * d1_control[i]) - (spline_data->v_max * spline_data->v_max); // constraint
    }
    for (int i = 0; i < n_c-d-1; i++)
    {
        d2_control[i] = (d1_control[i+1]-d1_control[i])*(k-1)/getdeltaT(i, k-1, knots_v);
        //result[i+n_c-d] = (d2_control[i] * d2_control[i]) - (spline_data->a_max * spline_data->a_max);
    }

    if (grad)
    {
        int m2 = 4;
        int opto_vars = n;
        int i = 1;
        int h = 0;
    
        for(int j = 0; j < opto_vars; j++)  // differentiate wrt to each of the opto variables
        {
            double dt = getdeltaT(i,k,knots);
            if (j < step_length)  // if we are taking wrt to a knot step then use this eq. BUT also need to make sure that knot step is used
            {
                if (j <= i && j >= i-k+1)   // chect that that knot step is used
                {
                    grad[h*opto_vars + j] = -k*(d1_control[i+1] - d1_control[i])/(dt*dt);
                }
                else
                {
                    grad[h*opto_vars + j] = 0;
                }
            }
            else    // differentiate wrt to the control points
            {
                if (j-step_length == i)
                {
                    grad[h*opto_vars+j] = -k/dt;
                }
                else if(j-step_length == i+1)
                {
                    grad[h*opto_vars+j] = k/dt;
                }
                else
                {
                    grad[h*opto_vars+j] = 0;
                }
            }
            //printf("%d | %d | %f | %f \n", i, j, dt, grad[i*opto_vars+j]);
            //printf("%f ", grad[h*opto_vars + j]);
        }
    }
    return d2_control[0];
}

double c4_constraint(unsigned n, const double* x, double * grad, void* data)
{
    c4_con_count++;
    bspline* spline_data = (bspline*) data;

    int n_c = spline_data->num_control_points;
    int step_length = spline_data->step_length;
    int k = spline_data->poly_deg;
    int d = 1;
    
    //printf("Knot Steps Clamp: \n");
    double knot_steps[step_length];
    for (int i = 0; i < step_length; i++)
    {
        knot_steps[i] = x[i];
        //printf("%f ", x[i]);
    }

    //printf("\nControl Points Clamp: \n");
    double control[n_c];
    for (int i = 0; i < n_c; i++)
    {
        control[i] = x[i+step_length];
        //printf("%f ", x[i+step_length]);
    }
    //printf("\n \n");
    // Recover the knots from the knot steps
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

    // Calculate the knots for the derivatives. They contain the same knot steps but less multiplicity on the ends
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

    double d1_control[n_c-d];
    double d2_control[n_c-d-1];

    for (int i = 0; i < n_c-d; i++)
    {
        d1_control[i] = (control[i+1]-control[i])*k/getdeltaT(i, k, knots);
        //result[i] = (d1_control[i] * d1_control[i]) - (spline_data->v_max * spline_data->v_max); // constraint
    }
    for (int i = 0; i < n_c-d-1; i++)
    {
        d2_control[i] = (d1_control[i+1]-d1_control[i])*(k-1)/getdeltaT(i, k-1, knots_v);
        //result[i+n_c-d] = (d2_control[i] * d2_control[i]) - (spline_data->a_max * spline_data->a_max);
    }

    if (grad)
    {
        int m2 = 4;
        int opto_vars = n;
        int i = n_c-d-2;
        int h = 0;
    
        for(int j = 0; j < opto_vars; j++)  // differentiate wrt to each of the opto variables
        {
            double dt = getdeltaT(i,k,knots);
            if (j < step_length)  // if we are taking wrt to a knot step then use this eq. BUT also need to make sure that knot step is used
            {
                if (j <= i && j >= i-k+1)   // chect that that knot step is used
                {
                    grad[h*opto_vars + j] = -k*(control[i+1] - control[i])/(dt*dt);
                }
                else
                {
                    grad[h*opto_vars + j] = 0;
                }
            }
            else    // differentiate wrt to the control points
            {
                if (j-step_length == i)
                {
                    grad[h*opto_vars+j] = -k/dt;
                }
                else if(j-step_length == i+1)
                {
                    grad[h*opto_vars+j] = k/dt;
                }
                else
                {
                    grad[h*opto_vars+j] = 0;
                }
            }
            //printf("%d | %d | %f | %f \n", i, j, dt, grad[i*opto_vars+j]);
            //printf("%f ", grad[h*opto_vars + j]);
        }
    }
    return d1_control[n_c-d-2];
}