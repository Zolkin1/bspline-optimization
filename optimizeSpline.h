// may want to use a struct for the spline but it might not be worth it
typedef struct
{
    int num_control_points;
    int poly_deg;
    double start_knot;
    int step_length;
    double v_max;
    double a_max;
    double j_max;
    double dist;
    double x_l[27];
} bspline;

void derivConstraints(unsigned m, double* result, unsigned n, const double* x, double* grad, void* data);
void clampConstraints(unsigned m, double* result, unsigned n, const double* x, double* grad, void* data);

double c1_constraint(unsigned n, const double* x, double * grad, void* data);
double c2_constraint(unsigned n, const double* x, double * grad, void* data);
double c3_constraint(unsigned n, const double* x, double * grad, void* data);
double c4_constraint(unsigned n, const double* x, double * grad, void* data);

double objectiveFunc(unsigned n, const double* x, double* grad, void* data);

double getKnots(double start_knot, double* knot_steps, int poly_deg, int index, int step_length);
double controlPointDerivative(int i, int derivative, int poly_deg, int n, double start_knot, double* knot_steps, double* control_points, int step_length);

double getdeltaT(int i, int poly_deg, double* knot_steps);