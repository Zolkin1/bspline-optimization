#include <vector>

class BSpline
{
    private:
    int poly_deg;
    int num_control_points;
    int num_knots;

    int start_knot;
    std::vector<double> knot_spaces;
    
    std::vector<double> getKnots();

    public:
    BSpline getDerivative();
}