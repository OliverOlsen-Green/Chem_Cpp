#ifndef ODE_functions_H
#define ODE_functions_H

#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include "vector_matrix.h"
#include <fstream>

using namespace std;

template <typename T> class derivatives {
    private:
    T data;
    public:

    double dev(double x, function<double(double)> f) {
        double h = 1e-7;
        double derivation = (f(x+h) - f(x)) / h;
        return derivation;
    }
    vector<double> vector_dev(double x, vector<function<double(double)>>f) {
        double h = 1e-7;
        vector<double> vec_derivation; //initialise to f size
        for(int i=0; i < f.size(); i++) {
            vec_derivation.push_back((f[i](x+h) - f[i](x)) / h); //loop through derivatives of vector
        }
        return vec_derivation;
    }

    double integrals(function<double(double)> f, double upper_lim, double lower_lim) {
        double steps = 1000;
        double step_size = (upper_lim - lower_lim) / steps ;
        double trap_sum = 0;
        for(int i = 1; i< steps; i++) {
            trap_sum += 2 * (f(lower_lim + (i * step_size)));
        }
        double integral = (0.5 * step_size) * (f(lower_lim) + trap_sum + f(upper_lim));
        return integral;
    }
    // function that returns integrals for vectors (multiple functions)
    vector<double> vec_integrals(vector<function<double(double)>> f, double upper_lim, double lower_lim) {
        double steps = 1000;
        double step_size = (upper_lim - lower_lim) / steps ;
        vector<double> vec_trap_sum = vector<double>(f.size(),0.0);
        for(int i = 1; i< steps; i++) {
        for(int j = 0; j< f.size(); j++){
            vec_trap_sum[j] += 2 * (f[j](lower_lim + (i * step_size)));
        }
    }
        vector<double> vec_integral = vector<double>(f.size(),0.0);
        for(int k = 0; k < f.size(); k++) {
            vec_integral[k] = (0.5 * step_size) * (f[k](lower_lim) + vec_trap_sum[k] + f[k](upper_lim));
        }
        return vec_integral;
    }
    vector<vector<double>> euler_ODE(vector<function<double(double)>> f, double y_0, 
    double y_final, double steps, double x) {
        double step_size = (y_final - y_0) / steps;
        vector<double> y(f.size(), y_0);
        vector<vector<double>> all_points;
        all_points.push_back(y);
        for(int i = 1 ; i < steps; i++){
            vector<double> derivs = vector_dev(x, f);
            for(int j = 0; j < f.size(); j++) {
                y[j] = y[j] + (step_size * derivs[j]);
            }
            x += step_size;
            all_points.push_back(y);
        }
        return all_points;
    }
    vector<vector<double>> Runga_Kutta(function<vector<double>(vector<double>)> f, vector<double> y_0, 
    double y_final, double y_i, double steps) {
        double h = (y_final - y_i) / steps;
        double x = y_i;
        vector<double> y = y_0 ;
        vector<vector<double>> all_points;
        all_points.push_back(y); 
        for(int i = 1; i < steps; i++) {
            vectormath<double> vm_d;
            vector<double> k1 = vm_d.vec_scalar(f(y), h);
            vector<double> k2 = vm_d.vec_scalar(f(vm_d.v_addition(y, vm_d.vec_scalar(k1, 0.5))), h);
            vector<double> k3 = vm_d.vec_scalar(f(vm_d.v_addition(y, vm_d.vec_scalar(k2, 0.5))), h);
            vector<double> k4 = vm_d.vec_scalar(f(vm_d.v_addition(y, k3)), h);
            for(int j = 0; j < y.size(); j++) {
                y[j] = y[j] + (1.0 / 6.0) * (k1[j] + (2 * k2[j]) + (2 * k3[j]) + k4[j]);
            }
            x += h;
            all_points.push_back(y);
        }
        return all_points;
    }
    vector<vector<double>> Jacobian(function<vector<double>(vector<double>)> odes, vector<double> c) {
        vector<double> c_J = c;
        int length = c.size();
        double h = 1e-7;
        vector<vector<double>> J = vector<vector<double>> (length, vector<double>(length, 0.0));
        for(int i = 0; i < length; i++) {
            c_J[i] += h;
            vector<double>F_J = odes(c_J);
            vector<double>F_0 = odes(c);
            for(int j =0; j < length; j++) {
                J[j][i] = (F_J[j] - F_0[j]) / h;
            }
            c_J[i] -= h;
        }
        return J;
    }
    double newton_raphson(double x, function<double(double)> f) {
    double error = 1;
    while (error > 0.000001) {
        double x_next = x - f(x) / dev(x,f);
        error = abs(x_next - x);
        x = x_next;
    }
    return x;
    }
    vector<double> Newton_roots(vector<double> x0, function<vector<double>(vector<double>)> F) {
    vectormath<double> vm_d;
        matrixmath<double> mm_d;
        double error = 1;
        double epsillon = 1e-7;
        while(error>epsillon) {
            vector<vector<double>> J = Jacobian(F, x0);
            auto LU = mm_d.LU_decomp(J);
            vector<double> F_x = F(x0); 
            vector<double> dx = mm_d.LU_substitution(LU.first, LU.second, F_x);
            vector<double> x_new = vm_d.v_subtraction(x0, dx);  // x = x - J^-1 * F(x)
            error = 0;
            for(int i = 0; i < dx.size(); i++) error += dx[i] * dx[i];
            error = sqrt(error);
            x0 = x_new;
        };
        return x0;
    }
    vector<double> Newton_ODE_system(vector<double> x0, function<vector<double>(vector<double>)> odes,
double step_size, vector<double> y_old) {
        vectormath<double> vm_d;
        matrixmath<double> mm_d;
        double error = 1;
        double epsillon = 1e-7;
        int iter = 0;
        while(error>epsillon) {
            iter++;
            vector<vector<double>> J = Jacobian(odes, x0);
            for(int i = 0; i < x0.size(); i++) {
                for(int j = 0; j < x0.size(); j++) {
                    J[i][j] = -step_size * J[i][j];
                }
                J[i][i] += 1.0;  // add identity
            };
            auto LU = mm_d.LU_decomp(J);
            int length = x0.size();
            vector<vector<double>> L = LU.first;
            vector<vector<double>> U = LU.second;
            vector<double> dx_dt = odes(x0);
            vector<double> F_x(length);
            for(int i = 0; i < length; i++) {
                F_x[i] = x0[i] - step_size * dx_dt[i] - y_old[i];
            };
            vector<double> x_LU = mm_d.LU_substitution(L, U, F_x);
            double alpha = 0.5;
            vector<double> x_new(length);
            for(int i = 0; i < length; i++) {
                x_new[i] = x0[i] - alpha * x_LU[i];
                if(x_new[i] < 0) x_new[i] = 1e-10;  
            };
            error = 0;
            for(int i = 0; i < x_LU.size(); i++) {
                error += x_LU[i] * x_LU[i];
            }
            error = sqrt(error);
            x0 = x_new;
        }
        cout << "Newton iterations: " << iter << endl;
        return x0;
    }
    vector<vector<double>> Implicit_newton(function<vector<double>(vector<double>)> f, vector<double> y_0, 
    double y_final, double y_i, double steps) {
        double step_size = (y_final - y_i) / steps;
        vector<double> y = y_0;
        vector<vector<double>> all_points;
        all_points.push_back(y);
        for(int i = 1; i < steps; i++) {
            vector<double> y_old = y;
            y = Newton_ODE_system(y, f,step_size,y_old);  // Newton solves for next y implicitly
            all_points.push_back(y);
        }
    return all_points;
    }
    double interpolate(vector<double> x_data, double x, vector<double> y_data) {
        int size = x_data.size();
        int i = 0;
        while(i < size - 2 && x < x_data[i+1]) {
            i++;
        }
        double xL = x_data[i];
        double yL = y_data[i];
        double xR = x_data[i+1];
        double yR = y_data[i+1];
        double dydx = (yR - yL) / (xR - xL);
        double y = yL + dydx;
        return y;
    }
};

#endif

/*
int main() {
    cout << fm.dev(3,f_1) << endl ;
    cout << fm.integrals(f_1,upper_lim,lower_lim) << endl; 
    for (double x : fm.vector_dev(3,f)) {
    cout << x << " ";
    };
    for (double x : fm.vec_integrals(f,upper_lim,lower_lim)) {
    cout << x << " ";
    };
    vector<vector<double>> result_ODE = fm.Runga_Kutta(f,y_0,y_final,steps,x_inital);
    for (int i = 0; i < result_ODE.size(); i++) {
    for (int j = 0; j < result_ODE[i].size(); j++) {
        cout << result_ODE[i][j] << " ";
    }
    cout << endl;}
// after your euler_ODE call
    ofstream file("/Users/oliverolsen-green/Desktop/MyCode/output.csv");
    for (int i = 0; i < result_ODE.size(); i++) {
        for (int j = 0; j < result_ODE[i].size(); j++) {
            file << result_ODE[i][j];
            if (j < result_ODE[i].size() - 1) file << ",";
        }
        file << "\n";
}
file.close();
};
*/
