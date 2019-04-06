#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>
#include <string.h>
int n,m;
double h,k;
double a=0.0, b=2.0, c=0.0, d=1.0;
double f(double x, double y){
    return x*exp(y);
}
double g(double x, double y){
    return x*exp(y);
}
double infinite_norm(double w1[n][m], double w2[n][m])
{   
    double norm = fabs(w1[1][1] - w2[1][1]);
    for (int i = 1; i < n; i++){
        for (int j = 1; j < m; j++){
            if (norm < fabs(w1[i][j] - w2[i][j])){
                norm = fabs(w1[i][j] - w2[i][j]);
            }
        }
    }
    return norm;
}



void gs(double*x, double*y, double lambda, double mu, double tolerance, int max_iter, double w[n][m], double* final_error, int* iter){
    *iter = 0;
    *final_error=1000.0;    
    double norm;
    for(int i=1;i<n;i++){
        for(int j=1; j<m; j++){
            w[i][j]=0.0;
        }
    }
    do{
        double z;
        
        z = (-h*h*f(x[1],y[m-1])+g(a, y[m-1])+lambda*g(x[1], d)+lambda*w[1][m-2]+w[2][m-1])/mu;
        norm = fabs(z-w[1][m-1]);
        w[1][m-1]=z;

        for (int i=1;i<n-1;i++){
            z = (-h*h*f(x[i],y[m-1])+lambda*g(x[i], d)+w[i-1][m-1]+w[i+1][m-1]+lambda*w[i][m-2])/mu;
            if (fabs(w[i][m-1]-z)>norm){
                norm = fabs(w[i][m-1]-z);
            };
            w[i][m-1]=z;
        }
        
        z = (-h*h*f(x[n-1],y[m-1])+g(b, y[m-1])+lambda*g(x[n-1], d)+w[n-2][m-1]+lambda*w[n-1][m-2])/mu;
        if (fabs(w[n-1][m-1]-z)>norm){
            norm = fabs(w[n-1][m-1]-z);
        };
        w[n-1][m-1]=z;

        for (int j=m-2;j>1;j--){            
            z = (-h*h*f(x[1],y[j])+g(a, y[j])+lambda*w[1][j+1]+lambda*w[1][j-1]+w[2][j])/mu;
            if (fabs(w[1][j]-z)>norm){
                norm = fabs(w[1][j]-z);
            };
            w[1][j]=z;

            for(int i=2; i<n-1;i++){
                z = (-h*h*f(x[i],y[j])+w[i-1][j]+lambda*w[i][j+1]+w[i+1][j]+lambda*w[i][j-1])/mu;
                if (fabs(w[i][j]-z)>norm){
                    norm = fabs(w[i][j]-z);
                };
                w[i][j]=z;
            }

            z = (-h*h*f(x[n-1],y[j])+g(b, y[j])+w[n-2][j]+lambda*w[n-1][j+1]+lambda*w[n-1][j-1])/mu;
            if (fabs(w[n-1][j]-z)>norm){
                norm = fabs(w[n-1][j]-z);
            };
            w[n-1][j]=z;
        }

        //AAAAAARIGHT
        z = (-h*h*f(x[1],y[1])+g(a, y[1])+lambda*g(x[1], c)+lambda*w[1][2]+w[2][1])/mu;
        if (fabs(w[1][1]-z)>norm){
            norm = fabs(w[1][1]-z);
        };
        w[1][1]=z;

        for (int i=2;i<n-1;i++){
            z = (-h*h*f(x[i],y[1])+lambda*g(x[i], c)+w[i-1][1]+lambda*w[i][2]+w[i+1][1])/mu;
            if (fabs(w[i][1]-z)>norm){
                norm = fabs(w[i][1]-z);
            };
            w[i][1]=z;
        }
        z = (-h*h*f(x[n-1],y[1])+g(b, y[1])+lambda*g(x[n-1], c)+w[n-2][1]+lambda*w[n-1][2])/mu;
        if (fabs(w[n-1][1]-z)>norm){
            norm = fabs(w[n-1][1]-z);
        };
        w[n-1][1]=z;

        (*iter)++;
    }while(*iter<max_iter && norm>tolerance); //&& *final_error>tolerance
    *final_error = norm;
}

void jacobi(double*x, double*y, double lambda, double mu, double tolerance, int max_iter, double w[n][m], double* final_error, int* iter){
    *iter = 0;
    *final_error=1000.0;    
    double norm;
    double w_next[n][m];
    for(int i=1;i<n;i++){
        for(int j=1; j<m; j++){
            w[i][j]=0.0;
        }
    }
    do{
        double z;
        
        z = (-h*h*f(x[1],y[m-1])+g(a, y[m-1])+lambda*g(x[1], d)+lambda*w[1][m-2]+w[2][m-1])/mu;
        w_next[1][m-1]=z;

        for (int i=1;i<n-1;i++){
            z = (-h*h*f(x[i],y[m-1])+lambda*g(x[i], d)+w[i-1][m-1]+w[i+1][m-1]+lambda*w[i][m-2])/mu;
            w_next[i][m-1]=z;
        }
        
        z = (-h*h*f(x[n-1],y[m-1])+g(b, y[m-1])+lambda*g(x[n-1], d)+w[n-2][m-1]+lambda*w[n-1][m-2])/mu;
        w_next[n-1][m-1]=z;

        for (int j=m-2;j>1;j--){            
            z = (-h*h*f(x[1],y[j])+g(a, y[j])+lambda*w[1][j+1]+lambda*w[1][j-1]+w[2][j])/mu;
            w_next[1][j]=z;

            for(int i=2; i<n-1;i++){
                z = (-h*h*f(x[i],y[j])+w[i-1][j]+lambda*w[i][j+1]+w[i+1][j]+lambda*w[i][j-1])/mu;
                w_next[i][j]=z;
            }

            z = (-h*h*f(x[n-1],y[j])+g(b, y[j])+w[n-2][j]+lambda*w[n-1][j+1]+lambda*w[n-1][j-1])/mu;
            w_next[n-1][j]=z;
        }

        z = (-h*h*f(x[1],y[1])+g(a, y[1])+lambda*g(x[1], c)+lambda*w[1][2]+w[2][1])/mu;
        w_next[1][1]=z;

        for (int i=2;i<n-1;i++){
            z = (-h*h*f(x[i],y[1])+lambda*g(x[i], c)+w[i-1][1]+lambda*w[i][2]+w[i+1][1])/mu;
            w_next[i][1]=z;
        }
        z = (-h*h*f(x[n-1],y[1])+g(b, y[1])+lambda*g(x[n-1], c)+w[n-2][1]+lambda*w[n-1][2])/mu;
        w_next[n-1][1]=z;
        norm = infinite_norm(w, w_next);
        for(int i=1;i<n;i++){
            for(int j=1; j<m; j++){
                w[i][j]=w_next[i][j];
            }
        }
        
        (*iter)++;
    }while(*iter<max_iter && norm>tolerance); //&& *final_error>tolerance
    *final_error = norm;
}

void sor(double*x, double*y, double omega, double lambda, double mu, double tolerance, int max_iter, double w[n][m], double* final_error, int* iter){
    *iter = 0;
    *final_error=1000.0;    
    double norm;
    double w_next[n][m];
    for(int i=1;i<n;i++){
        for(int j=1; j<m; j++){
            w[i][j]=0.0;
        }
    }
    /* for(int i=1;i<n;i++){
        printf("%d,1: %f\n", i,w[i][1]);
        for(int j=1; j<m; j++){
            w_next[i][j]=w[i][j];
        }
    } */
    do{
        double z;
        
        z = (-h*h*f(x[1],y[m-1])+g(a, y[m-1])+lambda*g(x[1], d)+lambda*w_next[1][m-2]+w_next[2][m-1])/mu;
        w_next[1][m-1]=z;

        for (int i=1;i<n-1;i++){
            z = (-h*h*f(x[i],y[m-1])+lambda*g(x[i], d)+w_next[i-1][m-1]+w_next[i+1][m-1]+lambda*w_next[i][m-2])/mu;
            w_next[i][m-1]=z;
        }

        z = (-h*h*f(x[n-1],y[m-1])+g(b, y[m-1])+lambda*g(x[n-1], d)+w_next[n-2][m-1]+lambda*w_next[n-1][m-2])/mu;
        w_next[n-1][m-1]=z;

        for (int j=m-2;j>1;j--){            
            z = (-h*h*f(x[1],y[j])+g(a, y[j])+lambda*w_next[1][j+1]+lambda*w_next[1][j-1]+w_next[2][j])/mu;
            w_next[1][j]=z;

            for(int i=2; i<n-1;i++){
                z = (-h*h*f(x[i],y[j])+w_next[i-1][j]+lambda*w_next[i][j+1]+w_next[i+1][j]+lambda*w_next[i][j-1])/mu;
                w_next[i][j]=z;
            }

            z = (-h*h*f(x[n-1],y[j])+g(b, y[j])+w_next[n-2][j]+lambda*w_next[n-1][j+1]+lambda*w_next[n-1][j-1])/mu;
            w_next[n-1][j]=z;
        }

        z = (-h*h*f(x[1],y[1])+g(a, y[1])+lambda*g(x[1], c)+lambda*w_next[1][2]+w_next[2][1])/mu;
        w_next[1][1]=z;

        for (int i=2;i<n-1;i++){
            z = (-h*h*f(x[i],y[1])+lambda*g(x[i], c)+w_next[i-1][1]+lambda*w_next[i][2]+w_next[i+1][1])/mu;
            w_next[i][1]=z;
        }
        z = (-h*h*f(x[n-1],y[1])+g(b, y[1])+lambda*g(x[n-1], c)+w_next[n-2][1]+lambda*w_next[n-1][2])/mu;
        w_next[n-1][1]=z;
        
        for(int i=1;i<n;i++){
            for(int j=1; j<m; j++){
                w_next[i][j]=(1-omega)*w[i][j]+omega*w_next[i][j];//(1-w)*u_k[i]+w*gs;
            }
        }
        norm = infinite_norm(w, w_next);
        for(int i=1;i<n;i++){
            for(int j=1; j<m; j++){
                w[i][j]=w_next[i][j];
            }
        }
        (*iter)++;
    }while(*iter<max_iter && norm>tolerance); //&& *final_error>tolerance
    *final_error = norm;
}
void fprint_vector(FILE *out,double w[n][m]){
    
}
void save_result(char* name_file, int* iter, double* final_error, double*x, double*y, double w[n][m]){
    FILE *out = fopen(name_file,"w");
    fprintf(out,"\niters:%d\nerror:%.10lf\nsolucio: \n",*iter, *final_error);
    for (int i=1; i<n; i++) {
        for(int j=1; j<m; j++){
            fprintf(out,"   %d %d | %f %f | %f %f \n",i,j,x[i],y[j],w[i][j], f(x[i],y[j]));
        }
    }
    
    fclose(out); 
}

void main(){
    double omega;
    printf( "Entra els valors: n m w (int, int, float [0.0,2.0]) : ");
    scanf("%d %d %lf", &n, &m, &omega);
    h=(b-a)/n;
    k=(d-c)/m;
    double* x = (double*)malloc(sizeof(double) * n);
    double* y = (double*)malloc(sizeof(double) * m);
    //type arrayName [ num_rows ][ num_columns ];
    double w[n][m]; 
    for(int i=1; i<n; i++) x[i]=a+i*h;
    for(int j=1; j<m; j++) y[j]=a+j*k;
    for (int i=1; i<n; i++) {
        for(int j=1; j<m; j++)
            w[i][j] = 0.0;
    }
    double lambda = (h*h)/(k*k);
    double mu = 2*(1+lambda);
    int max_iter = 1000;
    double tolerance = 1e-10;
    //gs(double*x, double*y, double lambda, double mu, double tolerance, int max_iter, double**w, double* final_error, int* iter)
    int* iter = malloc(sizeof(int));
    double* final_error = malloc(sizeof(double));
    jacobi(x,y,lambda,mu,tolerance,max_iter,w,final_error,iter);
    save_result("jacobi-BernatEsquirol.txt",iter,final_error,x,y,w); 
    gs(x,y,lambda,mu,tolerance,max_iter,w,final_error,iter);
    save_result("gauss-seidel-BernatEsquirol.txt",iter,final_error,x,y,w); 
    sor(x,y,omega,lambda,mu,tolerance,max_iter,w,final_error,iter);
    save_result("sor-BernatEsquirol.txt",iter,final_error,x,y,w); 
    /* for (int i=1; i<n; i++) {
        for(int j=1; j<m; j++){
            printf("%d %d | %f %f | %f %f\n",i,j,x[i], y[j], w[i][j], f(x[i],y[j]));             
        }
    } */
    


}