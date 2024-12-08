#include <bits/stdc++.h>
using namespace std;
#define db double

vector<vector<db>> X(1, vector<db>(166792, 1));
vector<vector<db>> Y;
vector<vector<db>> a(19, vector<db>(1, 0));

// Edit this function Later
//  Function to compute the transpose of a matrix
vector<vector<db>> transposeMatrix(const vector<vector<db>> &matrix)
{
    int rows = matrix.size();
    int cols = matrix[0].size();

    // Initialize the transposed matrix with reversed dimensions
    vector<vector<db>> transposed(cols, vector<db>(rows));

    // Compute the transpose
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            transposed[j][i] = matrix[i][j];
        }
    }

    return transposed;
}

void compute_X_Y()
{
    ifstream RE("../Columns/Region_East.csv");                    // X1
    ifstream RW("../Columns/Region_West.csv");                    // X2
    ifstream RN("../Columns/Region_North.csv");                   // X3
    ifstream RS("../Columns/Region_South.csv");                   // X4
    ifstream STCh("../Columns/Soil_Type_Chalky.csv");             // X5
    ifstream STCl("../Columns/Soil_Type_Clay.csv");               // X6
    ifstream STLo("../Columns/Soil_Type_Loam.csv");               // X7
    ifstream STPe("../Columns/Soil_Type_Peaty.csv");              // X8
    ifstream STSa("../Columns/Soil_Type_Sandy.csv");              // X9
    ifstream STSl("../Columns/Soil_Type_Slit.csv");               // X10
    ifstream DTH("../Columns/Days_to_Harvest.csv");               // X11  // Normalized is not used
    ifstream FU("../Columns/Fertilizer_Used.csv");                // X12
    ifstream IU("../Columns/Irrigation_Used.csv");                // X13
    ifstream Rm("../Columns/Rainfall_mm_Normalized.csv");         // X14
    ifstream TC("../Columns/Temperature_Celsius_Normalized.csv"); // X15
    ifstream WCC("../Columns/Weather_Condition_Cloudy.csv");      // X16
    ifstream WCR("../Columns/Weather_Condition_Rainy.csv");       // X17
    ifstream WCS("../Columns/Weather_Condition_Sunny.csv");       // X18
    ifstream Yield("../Columns/Yield_tons_per_hectare.csv");      // Y
    // Defining X as a 2D vector for all ../columns
    vector<db> X1(166792);
    vector<db> X2(166792);
    vector<db> X3(166792);
    vector<db> X4(166792);
    vector<db> X5(166792);
    vector<db> X6(166792);
    vector<db> X7(166792);
    vector<db> X8(166792);
    vector<db> X9(166792);
    vector<db> X10(166792);
    vector<db> X11(166792);
    vector<db> X12(166792);
    vector<db> X13(166792);
    vector<db> X14(166792);
    vector<db> X15(166792);
    vector<db> X16(166792);
    vector<db> X17(166792);
    vector<db> X18(166792);

    vector<db> y(166792); // For Y values

    // Reading data from each file and storing in vectors
    for (int i = 0; i < 166792; i++)
    {
        RE >> X1[i];
        RW >> X2[i];
        RN >> X3[i];
        RS >> X4[i];
        STCh >> X5[i];
        STCl >> X6[i];
        STLo >> X7[i];
        STPe >> X8[i];
        STSa >> X9[i];
        STSl >> X10[i];
        DTH >> X11[i];
        FU >> X12[i];
        IU >> X13[i];
        Rm >> X14[i];
        TC >> X15[i];
        WCC >> X16[i];
        WCR >> X17[i];
        WCS >> X18[i];
        Yield >> y[i];
    }

    // Pushing each vector into X as columns
    X.push_back(X1);
    X.push_back(X2);
    X.push_back(X3);
    X.push_back(X4);
    X.push_back(X5);
    X.push_back(X6);
    X.push_back(X7);
    X.push_back(X8);
    X.push_back(X9);
    X.push_back(X10);
    X.push_back(X11);
    X.push_back(X12);
    X.push_back(X13);
    X.push_back(X14);
    X.push_back(X15);
    X.push_back(X16);
    X.push_back(X17);
    X.push_back(X18);

    // Storing Y in a vector
    Y.push_back(y);

    X = transposeMatrix(X);
    Y = transposeMatrix(Y);
}

db average(vector<vector<db>> &A)
{
    db avg = 0.0;
    int n = 0;
    int size = A.size();
    for (int i = 0; i < size; i++)
    {
        n++;
        avg += (A[i][0] - avg) / n;
    }
    return avg;
}

// Edit this function later
//  Function to multiply two matrices
vector<vector<db>> multiplyMatrices(const vector<vector<db>> &mat1, const vector<vector<db>> &mat2)
{
    int r1 = mat1.size();
    int c1 = mat1[0].size();
    int r2 = mat2.size();
    int c2 = mat2[0].size();

    // Initialize the resultant matrix with 0
    vector<vector<db>> result(r1, vector<db>(c2, 0));

    // Matrix multiplication logic
    for (int i = 0; i < r1; i++)
    {
        for (int j = 0; j < c2; j++)
        {
            for (int k = 0; k < c1; k++)
            {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }

    return result;
}

// Function to calculate the inverse of the Jacobian matrix using Gaussian elimination
vector<vector<db>> GaussElim(vector<vector<db>> A, vector<vector<db>> B)
{
    int n = A.size();
    vector<vector<double>> X(n, vector<db>(1));

    // Gauss Elimination with partial pivoting
    for (int k = 0; k < n - 1; k++)
    {
        // Finding row with maximum coefficient
        double max = 0;
        int maxrow = k;
        for (int i = k; i < n; i++)
            if (abs(A[i][k]) > max)
            {
                max = abs(A[i][k]);
                maxrow = i;
            }

        // switch with pivot row
        swap(A[k], A[maxrow]);
        swap(B[k], B[maxrow]);

        // Eliminate col if pivot coeff is non zero
        if (A[k][k] != 0)
            for (int i = k + 1; i < n; i++)
            {
                double factor = A[i][k] / A[k][k];
                A[i][k] = 0;
                for (int j = k + 1; j < n; j++)
                    A[i][j] -= A[k][j] * factor;
                B[i][0] -= B[k][0] * factor;
            }
    }

    // Back Substitution
    for (int i = n - 1; i >= 0; i--)
    {
        double sum = B[i][0];
        for (int j = i + 1; j < n; j++)
            sum -= X[j][0] * A[i][j];
        X[i][0] = sum / A[i][i];
    }

    return X;
}

vector<vector<db>> linear_regression(vector<vector<db>> X, vector<vector<db>> Y)
{
    vector<vector<db>> Xt = transposeMatrix(X);
    vector<vector<db>> XtX = multiplyMatrices(Xt, X);
    vector<vector<db>> XtY = multiplyMatrices(Xt, Y);

    return GaussElim(XtX, XtY);
}

// Ypredicted = a[0]*X[0] + a[1]*X[1] + .......

void generateYpred(vector<vector<db>> &Ypredicted)
{
    ofstream Ypredicted_out("../Output/Ypredicted.csv");

    for (auto x : Ypredicted)
    {
        Ypredicted_out << x[0] << "\n";
    }
}

int main()
{
    compute_X_Y();
    // for (int i = 0; i < 10; i++)
    // {
    //     for (int j = 0; j < 19; j++)
    //     {
    //         cout << X[i][j] << " ";
    //     }
    //     cout << "\n";
    // }
    a = linear_regression(X, Y); // Coefficients

    vector<vector<db>> Ypredicted(166792, vector<db>(1, 0));
    vector<db> error2(166792, 0);
    db Yavg = average(Y);
    // cout << Yavg;
    // calculating error square
    for (int i = 0; i < 166792; i++)
    {
        Ypredicted[i][0] = 0;
        for (int j = 0; j < 19; j++)
        {
            Ypredicted[i][0] += a[j][0] * X[i][j];
        }
        db error = Y[i][0] - Ypredicted[i][0];
        error2[i] = error * error;
    }

    // Calculating Sr
    db Sr = 0;
    for (auto x : error2)
    {
        Sr += x;
    }

    // Calculating St
    db St = 0;
    for (auto x : Y)
    {
        db diff = x[0] - Yavg;
        St += diff * diff;
    }

    db r2 = 1 - Sr / St;

    cout << "R square is " << r2 << "\n";

    generateYpred(Ypredicted);

    cout << "Coefficients a0,a1,...,X18 are\n";

    for (int i = 0; i < 19; i++)
    {
        cout << fixed << setprecision(10) << "a" << i << " = " << a[i][0] << "\n";
    }
    cout << "Adjusted R square is " << 1 - (1 - r2) * (166792 - 1) / (166792 - 18 - 1) << "\n";

    // temperature_dependence_();
}