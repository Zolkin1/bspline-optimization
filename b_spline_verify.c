#include <stdio.h>

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

int main()
{
    int d = 1;
    int k = 5;
    int n = 11;
    int knot1 = 0;

    double control[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    double knot_steps[7] = {3, 3, 3, 3, 3, 3, 3};

    int step_length = 7;
    
    int n_c = 12;

    double d1_control[n_c-d];
    double d2_control[n_c-d-1];
    double d3_control[n_c-d-2];
    for (int i = 0; i < n; i++)
    {
        d1_control[i] = controlPointDerivative(i, d, k, n, 0, knot_steps, control, step_length);
    }

    for (int i = 0; i < n_c-d-1; i++)
    {
        d2_control[i] = controlPointDerivative(i, d, k-1, n_c-d-1, 0, knot_steps, d1_control, step_length);
    }

    for (int i = 0; i < n_c-d-2; i++)
    {
        d3_control[i] = controlPointDerivative(i, d, k-2, n_c-d-2, 0, knot_steps, d2_control, step_length);
    }

    printf("d1: \n");
    for (int i = 0; i < n_c-d; i++)
    {
        printf("%f ", d1_control[i]);
    }
    printf("\n");

    printf("d2: \n");
    for (int i = 0; i < n_c-d-1; i++)
    {
        printf("%f ", d2_control[i]);
    }
    printf("\n");

    printf("d3: \n");
    for (int i = 0; i < n_c-d-2; i++)
    {
        printf("%f ", d3_control[i]);
    }
    printf("\n");
}