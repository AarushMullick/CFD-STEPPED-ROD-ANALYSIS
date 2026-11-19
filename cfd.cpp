#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <climits>

using namespace std;

const double dX = 0.1;
const double dY = 0.1;

const int NX = 31;
const int NY = 21;

const double omega = 1.75;
const double tolerance = 1e-10;
const int maxIter = INT_MAX;

double Bi_values[] = {
    1e-3, 1e-2, 1e-1, 1, 10, 100
};

bool solid[NY][NX];

int number[NY][NX];

vector<pair<int, int>> points;

struct Equation {
    map<int, double> coefficient;
    double rhs;
};

vector<Equation> A;

double theta[NY][NX];

// Directions: East, West, North, South
int di[] = {1, -1, 0, 0};
int dj[] = {0, 0, 1, -1};

//crating stepped geometry

void createGeometry()
{
    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            double X = i * dX;
            double Y = j * dY;

            bool leftSection =
                (X <= 2.0 + 1e-12);

            bool rightSection =
                (X >= 2.0 - 1e-12 &&
                 Y >= 0.5 - 1e-12 &&
                 Y <= 1.5 + 1e-12);

            solid[j][i] = leftSection || rightSection;

            number[j][i] = -1;
        }
    }

    int count = 0;

    for (int j = 0; j < NY; j++)  //Pseudo codes in notes 
    {
        for (int i = 0; i < NX; i++)
        {
            if (solid[j][i])
            {
                number[j][i] = count++;

                points.push_back({j, i});
            }
        }
    }
}


//Finite difference matrix formation
void constructMatrix(double Bi)
{
    int N = points.size();

    A.clear();
    A.resize(N);

    for (int p = 0; p < N; p++)
    {
        int j = points[p].first;
        int i = points[p].second;

        double diagonal = 4.0;

        // RHS of finite difference equation
        A[p].rhs = dX * dX;

        for (int d = 0; d < 4; d++)
        {
            int jj = j + dj[d];
            int ii = i + di[d];

            bool inside =
                (jj >= 0 && jj < NY &&
                 ii >= 0 && ii < NX &&
                 solid[jj][ii]);

            if (inside)
            {
                // Interior neighbouring point
                int q = number[jj][ii];

                A[p].coefficient[q] -= 1.0;
            }

            else
            {
                // Check convective boundaries
                bool convection = false;

                if (d == 1 && i == 0)
                    convection = true;

                if (d == 0 && i == NX - 1)
                    convection = true;

                // Opposite inward grid point
                int jo = j - dj[d];
                int io = i - di[d];

                int q = number[jo][io];

                // Ghost-point / image-point contribution
                A[p].coefficient[q] -= 1.0;

                if (convection)
                {
                    diagonal += 2.0 * Bi * dX;
                }
            }
        }

        A[p].coefficient[p] += diagonal;
    }
}


//matrix solving function 

bool solveMatrix()
{
    int N = points.size();

    vector<double> T(N, 0.0);

    for (int iteration = 0;
         iteration < maxIter;
         iteration++)
    {
        double maxChange = 0.0;

        for (int p = 0; p < N; p++)
        {
            double diagonal =
                A[p].coefficient[p];

            double sum = 0.0;

            for (auto entry : A[p].coefficient)
            {
                int q = entry.first;
                double value = entry.second;

                if (q != p)
                {
                    sum += value * T[q];
                }
            }

            double newValue =
                (A[p].rhs - sum) / diagonal;

            // Successive Over-Relaxation
            double updatedValue =
                (1.0 - omega) * T[p]
                + omega * newValue;

            maxChange = max(
                maxChange,
                fabs(updatedValue - T[p])
            );

            T[p] = updatedValue;
        }

        if (maxChange < tolerance)
        {
            // Store calculated temperatures
            for (int p = 0; p < N; p++)
            {
                int j = points[p].first;
                int i = points[p].second;

                theta[j][i] = T[p];
            }

            cout << "Converged in "
                 << iteration + 1
                 << " iterations\n";

            return true;
        }
    }

    cout << "Maximum iterations reached.\n";

    return false;
}

void calculateResults(double Bi)
{
    double Tmin = INT_MAX;
    double Tmax = INT_MIN;

    for (auto point : points)
    {
        int j = point.first;
        int i = point.second;

        Tmin = min(Tmin, theta[j][i]);
        Tmax = max(Tmax, theta[j][i]);
    }

    cout << "\nBi = " << Bi << endl;

    cout << "Minimum theta = "
         << Tmin << endl;

    cout << "Maximum theta = "
         << Tmax << endl;
}


//file creation to create csv files with data from code

void saveResults(double Bi)
{
    string filename =
        "Temperature_Bi_" + to_string(Bi) + ".csv";

    ofstream file(filename);

    file << "X,Y,Theta\n";

    for (auto point : points)
    {
        int j = point.first;
        int i = point.second;

        file << fixed << setprecision(6)
             << i * dX << ","
             << j * dY << ","
             << theta[j][i] << "\n";
    }

    file.close();
}


//calling functions in int main function
int main()
{
    cout << fixed << setprecision(6);

    createGeometry();

    cout << "Number of solid grid points = "
         << points.size() << endl;

    for (double Bi : Bi_values)
    {
        cout << "\n-----------------------------\n";

        cout << "Solving for Bi = " << Bi << endl;

        // Construct [A]{theta} = {b}
        constructMatrix(Bi);

        // Solve simultaneous equations
        bool converged = solveMatrix();

        if (converged)
        {
            calculateResults(Bi);

            saveResults(Bi);
        }
        else
        {
            cout << "Results not saved for this Bi.\n";
        }
    }

    return 0;
}