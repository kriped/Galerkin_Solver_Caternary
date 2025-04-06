#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <stdio.h>
#include <direct.h>
#include <cmath>

using namespace std;

// Create Galerkin solver class for catenary equation d^2y/dx^2 = sqrt(1+(dy/dx)^2)
class Galerkin_solver_catenary
{
    vector<double> y, dydx, f;
    vector<vector<double>> basisk;
    vector<double> weights;
    double xIndex;
    double spaceStep;
    double gridSize;    
    int N;
    string filename;

public:
    
    Galerkin_solver_catenary(int gridSize, int N): gridSize(gridSize), N(N) //Constructor
    {
        y.resize(gridSize);
        dydx.resize(gridSize);
        f.resize(gridSize);
        basisk = vector<vector<double>>(gridSize, vector<double>(N));
        weights.resize(N);
        spaceStep = 2/double(gridSize);
    };

    void xBegin()
    {
        xIndex = 0;
    }

    void y_set(double Y)
    {
        y[xIndex] = Y;    
    }
    // Update spatial point and test if we are at the end
    
    bool x_next(double &CurrentCoordinate)
    {
        xIndex++;
        if (xIndex > gridSize-1) return false;
        CurrentCoordinate = -1.0 + (xIndex + 0.5)*spaceStep;
        return true;
    }
    void make_basisk(double Y,int k)
    {
        basisk[xIndex][k] = Y;
    }
    // Define boundary condition [-1,1]; y(-1) = y(1) = 0;
    void f_set()
    { 
        if (xIndex == 0 || xIndex == gridSize - 1) {
            // Enforce boundary conditions
            y[xIndex] = 0;
            dydx[xIndex] = 0;
            f[xIndex] = 0;
        } else {
            dydx[xIndex] = (y[xIndex + 1] - y[xIndex]) / spaceStep;
            f[xIndex] = sqrt(1 + pow(dydx[xIndex], 2));
        }
        // Debugging output
        cout << "xIndex: " << xIndex << ", dydx: " << dydx[xIndex] << ", f: " << f[xIndex] << endl;
    }
        

    void set_weights(int k)
    {
        weights[k] = 0;
        for (int i = 0; i < gridSize; i++)
        {
            weights[k] += f[i]*basisk[i][k];
        }
        
        
    } 

    void y_add_basisk(int k) // Sinus functions second derivative is a sin functions and inner product becomes delta function
    { 
        y[xIndex] += weights[k] * basisk[xIndex][k];
    } 
    
    void write_y_to_file(const string& filename)
    {
        ofstream outFile(filename);
        if (!outFile)
        {
            cerr << "Error: Could not open file " << filename << " for writing." << endl;
            return;
        }
        for (const auto& value : y)
        {
            outFile << value << endl;
        }
        outFile.close();
        cout << "Results written to " << filename << endl;
    }
     
};



int main()
{   
    cout << "Solving catenary equation..." << endl;
    int gridSize = 100;
    double spaceStep = 2/double(gridSize);
    double alpha = 1.0; 
    int N = 2;
    cout <<"spaceStep = " << spaceStep << endl;

    double x = -1.0 + 0.5 * spaceStep; // Initialize x with the first spatial coordinate
    Galerkin_solver_catenary mySolver(gridSize,N);
    // initialise problem
    mySolver.xBegin();
    // Define initial guess parabolic function y0 = alpha*x^2-1 and calculate f
    while(mySolver.x_next(x))
    {
       mySolver.y_set(alpha*x*x-1.0); // Corrected to x*x for proper squaring
       mySolver.f_set();
    }
    
    // Add N basis functions     
    for (int i = 1; i <= N; i++) 
    {
        mySolver.xBegin();
        while(mySolver.x_next(x))
        {
            mySolver.make_basisk(sin(i*x),i-1);
        }
        mySolver.set_weights(i-1);
        mySolver.xBegin();
        while(mySolver.x_next(x))
        {
            mySolver.y_add_basisk(i-1);
        }
    }

    // Write results to file
    mySolver.write_y_to_file("results.out");

    return 0;
}












