/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: sergei
 *
 * Created on June 3, 2017, 7:12 PM
 */

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

const int m=500;     //pocet kroku
double e;
int a=0;            //krajni left bod
int b=1;            //krajni right bod
double h=double(b)/double(m);
double** F=new double*[m];
double** W=new double*[m];       
double* u=new double[m];        

void printMatrix(double **M){
    for (int i = 0; i < m; i++){
        for (int j = 0; j < m; j++){
            cout << M[i][j] << "|";
        }
        cout << endl;
    }
    cout << endl;
}

void printVector(double *L){
    for (int i = 0; i < m; i++){
            cout << L[i] << "|";
    } 
    cout << endl;
}

void emptyVector( double *L){
    for(int i = 0; i < m; i++){
            L[i]=0;
    }
}

void fillPol(double *Q, double *k){
    for (int i = 1; i < m-1; i++){
            Q[i]=-k[i]*k[i]*k[i]+3.0/2.0*k[i]*k[i]-1.0/2.0*k[i];
    } 
}

void fillPolDer(double *F, double *k){
    for (int i = 1; i < m-1; i++){
            F[i]=-3*k[i]*k[i]+3*k[i]-0.5-2*e*e/h/h;
    } 
}

void vecMat(double *L, double **M, double *u){
    for(int i = 0; i < m; i++){
        for(int j = 0; j < m; j++){
            L[i]+=M[i][j]*u[j];
        }
    }
}

void eqVec(double *u,double *v){
    for(int i = 0; i < m; i++)
        u[i]=v[i];
}

void fillMat(double **M, double *L){
    for (int i = 1; i < m-1; ++i){
        int j=i;
        W[i][j]=L[i];
        W[i][j-1]=e*e/h/h;
        W[i][j+1]=e*e/h/h; 
    }  
}

void solMatNewton(double**M, double*L, double*u, int n){
    for(int k = 1; k < n; k++){
        L[k] = (L[k]-M[k][k-1] * L[k-1]) / (M[k][k]-M[k][k-1] * M[k-1][k]);
        if( k < n )
        {
        M[k][k+1]= M[k][k+1]/ (M[k][k]-M[k][k-1] * M[k-1] [k ]);
        M[k][k] = 1.0;
        M[k][k-1]=0.0;
        }
    }
    for(int k = n - 1; k >= 0; k--)
    L[k] = L[k]-M[k][k+1] * L[k + 1];
    for(int i=0;i<=n;i++)
    u[i]=L[i];
}

double MaxSlozka(double *L){
    double g;
    for (int i=0; i < m; i++){
        if(L[i] > g)
            g=L[i];
    }
    return g;
}

void soltotxt(double* solution, int m){
    double X[m];
    fstream file;
    file.open( "solution.txt", ios::out );
    for(int i = 0; i <= m-1; i++){
        X[i]=(double)i/(double)m;
        file << X[i] <<" "<< solution[i] << endl;
    }
    cout << "Data jsou v txt." << endl;
}





int main(int argc, char** argv) {
/*cout << "Zadejte parametr 'e' z intervalu (0;0.5):" << endl;
    cin >> e;
    while ( e <= 0 || e >= 0.5){
        cout << "Chyba. Zadejte parametr 'e' z intervalu (0;0.5):" << endl;
        cin >> e;   
    } 
    cout << "e = " << e << endl;
 */
    e=0.369;
    for(int i = 0; i < m; i++){
        F[i]=new double[m];
        W[i]=new double[m];
    }
    double *Q=new double[m]; 
    double *p=new double[m];
    double *L=new double[m];
    double *P=new double[m];
    double *v=new double[m];
    double *PS=new double[m];
    double *U=new double[m];
    int it=10000;
    double t=0.000001;
    emptyVector(PS);
    emptyVector(L);
    emptyVector(p);
    emptyVector(P);
    emptyVector(Q);
    emptyVector(v);
    emptyVector(U);
     
    W[0][0]=0;
    W[m-1][m-1]=1;   
    F[0][0]=0;
    F[m-1][m-1]=1;
    Q[0]=0;
    Q[m-1]=-50;
    
    for(int i = 1; i < m-1; i++)
        u[i]=1;
        u[0]=0;
        u[m-1]=1;
    
    for (int i = 1; i < m-1; ++i){
        int j=i;
            F[i][j]=-2*e*e/h/h;
            F[i][j-1]=e*e/h/h;
            F[i][j+1]=e*e/h/h;
    }
        //printMatrix(F);
        //cout << endl;
    
    for (int i=0; i<m; i++)
        U[i]=abs(u[i]-p[i]);
    
    //cout << "Zadejte pocet iteracnich kroku 'it' a presnost 't' : " << endl;
    //cin >> it;
    //cin >> t;
    int i=0;
    while (MaxSlozka(U) > t && i < it){
        fillPolDer(L,u);
        fillMat(W,L);
        //printMatrix(W);
        fillPol(Q,u);
        emptyVector(P);
        vecMat(P,F,u);
        for(int i = 0; i < m; i++)
            PS[i]=-P[i]-Q[i];
        solMatNewton(W,PS,v,m);     
        for (int i=0; i<m; i++){
            p[i]=v[i]+u[i];
            U[i]=abs(v[i]);       
        }
        
        eqVec(u,p);
        i++;
    }
    //cout << "Konec iterace na kroku :" << i << endl;
    //cout << "Reseni : " << endl;
    //printVector(p);
    soltotxt(p,m);
    
    return 0;
}

