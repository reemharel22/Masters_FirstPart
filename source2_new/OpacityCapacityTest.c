#define N 200
double rho = 0.05;
double arad = 7.56E-15;
double alpha = 3.5;
double beta = 1.1;
double mu_sio = 0.1;
double s_lambda = 0.75;
double s_f = 8.8E13;
double s_g = 1.0/9175.0;

double getCv1(double t) {
    double abb = beta * s_f * 
        pow(t , beta - 1.0) * pow(rho, -mu_sio);
        return abb;
    
}

double getOpacity1(double t) {
    return (pow(rho,s_lambda) * rho
                          / (s_g * pow(t, alpha)));

}


int main() {
    int i, j;
    double min_T = 100, max_T = 300; // in ev
    double T[N];
    double sigma_a[N];
    double cv[N];
    s_f = s_f / (pow(1160452, beta));
    s_g = s_g / pow(1160452, alpha);
    for (i = 0; i < N; i++) {
        T[i] = min_T + (max_T - min_T) * (i + 1) / (N + 1);
        T[i] = 11604.52 * T[i];
    }
    for (i = 0; i < N; i++) {
        cv[i] = getCv1(T[i]);
        sigma_a[i] = getOpacity1(T[i]);
        printf("Cv: %10e\t Opacity: %10e\n",cv[i],sigma_a[i]);
    }
    return 0;
}

