#ifndef VECTOR_MATRIX_H
#define VECTOR_MATRIX_H

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

template <typename T> class vectormath {
    private:
    T data;

    public:
    // make a vector from input data
    // product of a vector * a scalar
    vector<T> vec_scalar(vector<T> v1, double scalar) {
        int length = v1.size();
        vector<T> v_new = vector<T>(length,0.0);
        for(int i = 0; i < length; i++) {
            v_new[i] = v1[i] * scalar;
        }
    return v_new;
    }
    // divide a scalar by a vector
    vector<double> scalar_vec_div(vector<double> f, double S) {
    int size = f.size();
    vector<double> s_vec = vector<double>(size,0.0);
    for(int i = 0; i < size; i++) {
        s_vec[i] = S / f[i];
        }
    return s_vec;
    }
    // dividing a vector by a scalar
    vector<T> vec_scalar_div(vector<T> v1, double scalar) {
        int length = v1.size();
        vector<T> v_new = vector<T>(length,0.0);
        for(int i = 0; i < length; i++) {
            v_new[i] = v1[i] / scalar;
        }
    return v_new;
    }
    // dot product of two vectors
    double dot_product(vector<double> v1, vector<double> v2) {
        double v_dot_sum = 0;
        int length = v1.size();
        for(int i=0; i < length; i++) {
            v_dot_sum += v1[i] * v2[i];
        }
    return v_dot_sum;
    }   
    // addition of two vectors 
    vector<double> v_addition(vector<double> v1, vector<double> v2) {
        int length = v1.size();
        vector<double> v_add = vector<double>(length, 0.0); // initialise as zeros
        for(int i =0; i < length; i++) {
            v_add[i] = v1[i] + v2[i];
        } 
        return v_add;
    }
    // two vectors timesed together
    vector<double> v_times(vector<double> v1, vector<double> v2) {
        int length = v1.size();
        vector<double> v_add = vector<double>(length, 0.0); // initialise as zeros
        for(int i =0; i < length; i++) {
            v_add[i] = v1[i] * v2[i];
        } 
        return v_add;
    }
    // subtraction of two vectors
    vector<double> v_subtraction(vector<double> v1, vector<double> v2) {
        int length = v1.size();
        vector<double> v_sub = vector<double>(length, 0.0); // initialise as zeros
        for(int i =0; i < length; i++) {
            v_sub[i] = v1[i] - v2[i];
        } 
        return v_sub;
    }
    // dyasdic product of two vectors
    vector<vector<double>> tensor(vector<double> v1, vector<double> v2) {
        int length = v1.size();
        int length_1 = v2.size();
        vector<vector<double>> dy_product = vector<vector<double>>(length, vector<double>(length_1, 0.0));
        for(int i = 0; i < length; i++) {
            for(int j = 0 ; j < length_1; j++) {
                dy_product[i][j] = v1[i] * v2[j];
            }
        }
    return dy_product;
    }
    double sum(vector<double> f) {
        int size = f.size();
        double sum = 0;
        for(int i = 0; i < size; i++) {
            sum += f[i]; 
        }
    return sum;
    }
};

template <typename T> class matrixmath {
    private:
    T data;
    public:
    vector<vector<double>> matrix_addition(vector<vector<double>>m1, vector<vector<double>> m2) {
        int length = m1[0].size(); //access number of columns 
        int length_1 = m2.size(); // number of rows
        vector<vector<double>> m_add = vector<vector<double>>(length_1, vector<double>(length, 0.0));
        for(int i = 0; i < length; i++) {
            for(int j = 0; j<length; j++) {
                m_add[i][j] = m1[i][j] + m2[i][j]; //add each position of the matrix together
            }
        }
    return m_add;
    }
    vector<vector<double>> matrix_subtraction(vector<vector<double>>m1, vector<vector<double>> m2) {
        int length = m1[0].size(); //access number of columns 
        int length_1 = m2.size(); // number of rows
        vector<vector<double>> m_sub = vector<vector<double>>(length_1, vector<double>(length, 0.0));
        for(int i = 0; i < length; i++) {
            for(int j = 0; j<length; j++) {
                m_sub[i][j] = m1[i][j] - m2[i][j]; //add each position of the matrix together
            }
        }
    return m_sub;
    }
    // transpose a matrix
    vector<vector<double>> matrix_transpose(vector<vector<double>> m1) {
        int length = m1[0].size(); //access number of columns 
        int length_1 = m1.size(); // number of rows
        vector<vector<double>> m_T = vector<vector<double>>(length, vector<double>(length,0.0)); //swap the rows + cols
        for(int i = 0; i<length; i++) {
            for(int j = 0; j< length_1; j++) {
                m_T[i][j] = m1[j][i];
            }
        }
    return m_T;
    }
    vector<vector<double>> matrix_multiplication(vector<vector<double>>m1, vector<vector<double>> m2) {
        int length = m1[0].size(); //access number of columns 
        int length_1 = m2.size(); // number of rows
        vector<vector<double>> m_mult = vector<vector<double>>(length_1, vector<double>(length, 0.0));
        for(int i = 0; i < length; i++) {
            for(int j = 0; j < length_1; j++) {
                for(int k = 0; k < length; k++) {
                    m_mult[i][j] += m1[i][k] * m2[j][k];
                }
            }
        }
    return m_mult;
    }
    pair<vector<vector<double>>, vector<vector<double>>> LU_decomp(vector<vector<double>> f) {
    int rows = f.size(); //rows
    int cols = f[0].size(); //columns
    vector<vector<double>> U = vector<vector<double>> (rows, vector<double>(cols, 0.0));
    vector<vector<double>> L = vector<vector<double>> (rows, vector<double>(cols, 0.0));
    for(int k = 0; k < rows; k++) {
        int max_row = k;
        double max_val = abs(f[k][k]);
        for(int p = k + 1; p < rows; p++) {
            if(abs(f[p][k]) > max_val) {
                max_val = abs(f[p][k]);
                max_row = p;
            }
        }
        swap(f[k], f[max_row]);
        L[k][k] = 1.0;
        for(int i = k; i < cols; i++) {
            double sum = 0.0;
            for(int j = 0; j < k; j++) {
                sum += L[k][j] * U[j][i];
            }
            U[k][i] = f[k][i] - sum;
        }
        for(int z = k + 1; z < rows; z++) {
            double sum_L = 0.0;
            for(int d = 0; d < k; d++) {
                sum_L += L[z][d] * U[d][k];
            }
            L[z][k] = (f[z][k] - sum_L) / U[k][k]; 
        }
    }
    return make_pair(L, U);
    }

    vector<double> LU_substitution(vector<vector<double>> L, vector<vector<double>> U, vector<double> b) {
    int length = b.size();
    vector<double> y = vector<double>(length, 0.0);
    y[0] = b[0];
    for(int i = 1; i < length; i++) {
        double sum = 0 ;
        for(int j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }
    y[i] = b[i] - sum;
    }
    vector<double> x = vector<double>(length, 0.0);
    for(int ii = length - 1 ; ii >= 0; ii--) {
        double sum_b = 0;
        for(int jj = ii + 1 ; jj < length; jj++) {
            sum_b += U[ii][jj] * x[jj];
        }
    x[ii] = (y[ii] - sum_b) / U[ii][ii];
    }
    return x;
}
};

#endif

/*
int main() {
    vectormath<vector<double>> vm; //calling the class for the vector maths
    cout << vm.dot_product(v3,v2) << endl; //calling the function dot_product int he class
    for (double x : vm.vec_scalar(v1,scalar)) {
    cout << x << " ";
    };
    cout << endl; // for loop that prints the outlet for the new vector using the scalar
    // output code for 2D matrix (3D just add int z for loop)
    vector<vector<double>> result = vm.tensor(v1,v2);
    for (int i = 0; i < result.size(); i++) {
    for (int j = 0; j < result[i].size(); j++) {
        cout << result[i][j] << " ";
    }
    cout << endl;}
    for (int i = 0; i < m3.size(); i++) {
    for (int j = 0; j < m3[i].size(); j++) {
        cout << result[i][j] << " ";
    }
    cout << endl;}
    return 0;
}
*/
