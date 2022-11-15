//
// Created by drn on 05/12/19.
// Resolution de l'equation de Laplace en 2D par la methode des differences finis!
//

#include "laplace.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <time.h>

int main (void){

    double L = 1;
    int N = 500;    // Discretidation en abcisse
   
    // double H = 1.25;
    // int M = N*H/L;    // Discretidation en ordonne
    double M = 100;
    double H = M*L/N;


    laplace lplc;

    time_t start, end;

    start = time(NULL);

    laplace_init(&lplc, L, N, H, M);

    laplace_solve(&lplc);

    end = time(NULL);

    printf("Time taken: %f seconds\n", difftime(end, start)); 

    laplace_display(&lplc);
}

void laplace_init (laplace *lplc, double L, int N, double H, int M){
    lplc->L = L;
    lplc->N = N;
    lplc->H = H;
    lplc->M = M;
    lplc->n = (N-1)*(M-1);
    lplc->nb_aretes = (M-1)*(N-2) + (N-1)*(M-2);
    // assert (L/N == H/M)

    lplc->extremite_arete = malloc(lplc->nb_aretes * sizeof(int *));
    for (int i = 0; i < lplc->nb_aretes; i++){
        lplc->extremite_arete[i] = malloc(2 * sizeof(int));
    }
    
    // printf("Test: %d\n", lplc->extremite_arete[lplc->nb_aretes+1][lplc->nb_aretes+1]);
    
    // On verifie que le pas h est le meme dans les deux axes
    assert(L/N == H/M);
    lplc->h = L/N;

    lplc->x = malloc((N+1) * sizeof(double));
    lplc->y = malloc((M+1) * sizeof(double));
    lplc->u = malloc((N+1) * sizeof(double *));
    for (int i = 0; i < N+1; i++){
        lplc->u[i] = malloc((M+1) * sizeof(double));
    }
    lplc->U = malloc(lplc->n * sizeof(double));
    lplc->f = malloc((N+1) * sizeof(double *));       // A noter que fait f est connu sur toute la grille 
    for (int i = 0; i < N+1; i++){
        lplc->f[i] = malloc((M+1) * sizeof(double));
    }
    lplc->F = malloc(lplc->n * sizeof(double));

    // On remplit les differents points de discretisation en abcisse
    for (int i = 0; i < N+1; i++){
        lplc->x[i] = lplc->h * i;
    }

    // On remplit les differents points de discretisation en ordonnee
    for (int j = 0; j < M+1; j++){
        lplc->y[j] = lplc->h * j;
    }

    Skyline *A = &(lplc->A);

    InitSkyline(A, (N-1)*(M-1));        // (N-1)*(M-1) est le nombre d'equations qu'on a dans notre systeme Aline a resoudre

    // Si l'on SetSkyline alors qu'on a pas SwitchOn, alors il y a probleme  
    // for (int i = 0; i<lplc->n; i++){
    //     for (int j = 0; j<lplc->n; j++)
    //         SwitchOn(A, i, j);
    // }
    for (int i = 0; i<lplc->n-1; i++){
            SwitchOn(A, i, i+1);
    }
    // for (int i = 0; i<lplc->n-2; i++){
    //     SwitchOn(A, i, i+2);
    // }
    for (int i = 0; i<lplc->n-N+1; i++){
        SwitchOn(A, i, i+N-1);
    }
    for (int j = 0; j<lplc->n-1; j++){
        SwitchOn(A, j+1, j);
    }
    // for (int j = 0; j<lplc->n-2; j++){
    //     SwitchOn(A, j+2, j);
    // }
    for (int j = 0; j<lplc->n-N+1; j++){
        SwitchOn(A, j+N-1, j);
    }

    AllocateSkyline(A);

    double h = lplc->h;
    double val = 4 / (h*h);      // Valeur a utiliser pour remplir la diagonale

    // On remplit la diagonale principale et les diagonales inf et sup 2
    for(int i=0; i<lplc->n; i++){
        SetSkyline(A, i, i, val);
    }
    // for (int i = 0; i<lplc->n-2; i++){
    //     SetSkyline(A, i, i+2, 0);
    // }
    // for (int j = 0; j<lplc->n-2; j++){
    //     SetSkyline(A, j+2, j, 0);
    // }

    // DisplaySkyline(A);

    // // connection des extremites de chauque arete
    // // Pour les aretes horizontales
    int a = 0;
    int l;
    int k;
    // for (int j = 0; j < M-1; j++){
    //     for (int i = 0; i < N-2; i++){
    //         l = i + j*(N-1);
    //         k = i+1 + j*(N-1);
    //         lplc->extremite_arete[a][0] = l;
    //         lplc->extremite_arete[a][1] = k;
    //         a += 1;
    //     }
    // }
    
    // // Pour les aretes verticales
    // // int a = (M-1)*(N-2);
    // for (int i = 0; i < N-1; i++){
    //     for (int j = 0; j < M-2; j++){
    //         l = i + j*(N-1);
    //         k = i + (j+1)*(N-1);
    //         lplc->extremite_arete[a][0] = l;
    //         lplc->extremite_arete[a][1] = k;
    //         a += 1;
    //     }
    // }

    // VIsulisationd es aretes
    for (int a = 0; a < lplc->nb_aretes; a++){
        // printf("arete %d: %d, %d\n", a,  lplc->extremite_arete[a][0], lplc->extremite_arete[a][1]);
    }
    
    // assert(a== nb_aretes + 1);

    val = -1 / (h*h); 
    // On remplit la matrice A par la methode des extremites
    // for (int a = 0; a < lplc->nb_aretes; a++){
    //     l = lplc->extremite_arete[a][0];
    //     k = lplc->extremite_arete[a][1];
    //     SetSkyline(A, l, k, val);
    //     SetSkyline(A, k, l, val);
    // }
    
    // On remplit la matrice A par la methode des restes
    for (int k = 0; k<lplc->n-N+1; k++){
        SetSkyline(A, k, k+N-1, val);
        SetSkyline(A, k+N-1, k, val);
    }
    for (int k = 0; k < lplc->n-1; k++){
        if (!((k+1)%(N-1)==0))
            SetSkyline(A, k, k+1, val);     // Rempli les diagonales superieures +1
    }
    for (int k = 1; k < lplc->n; k++){
        if (!((k)%(N-1)==0))
            SetSkyline(A, k, k-1, val);     // Rempli les diagonales inferieures -1
    }
    
    
    
    // Ceci display la matrice avec des etoiles a moins qu'on definisse _FULL dans la declaration de DisplayAline
    // DisplaySkyline(A);

    // printf("nb d'aretes = %d\n", lplc->nb_aretes);
    // printf("derniere aretes = %d, origine=%d, fin=%d\n", a, lplc->extremite_arete[a-1][0], lplc->extremite_arete[a-1][1]);
    // printf("element A[%d][%d] = %f\n", 11, 15, val);
    // printf("element A[%d][%d] = %f\n", 11, 15, GetSkyline(A, 11, 15));

    // Facto LU
    FactoLU(A);

    // DisplaySkyline(A);

    // Second membre
    // for (int i = N/4; i <= 3*N/4; i++){
    for (int i = 0; i < N+1; i++){
        // for (int j = M/4; j <= 3*M/4; j++)
        for (int j = 0; j < M+1; j++)
            lplc->f[i][j] = 10;         // A calculer plus tard
    }
    
    // Remplissage du vecteur membre de droite
    for (int i = 0; i < N-1; i++){
        for (int j = 0; j < M-1; j++){
            l = i + j*(N-1);
            lplc->F[l] = lplc->f[i][j];
            // printf("F[%d] = %f\n", l, lplc->F[l]);
        } 
    }   
}

void laplace_solve(laplace *lplc){

    // Resolution du systeme
    SolveSkyline(&lplc->A, lplc->F, lplc->U);
    
    int l;
    // Retrouvons les resultats 2D
    for (int j = 0; j < lplc->M-1; j++){
        for (int i = 0; i < lplc->N-1; i++){
            l = i + j*(lplc->N-1);
            // printf("U[%d] = %f     ", l, lplc->U[l]);
            lplc->u[i+1][j+1] = lplc->U[l];
            // printf("u[%d][%d] = %f\n", i+1, j+1, lplc->U[l]);
        }      
    }
}

void laplace_display(laplace *lplc){
    FILE *plotfile;

    plotfile = fopen("plot.csv", "w");
    fprintf(plotfile, "%s, %s, %s, %s\n", "x", "y", "etat final", "etat initial");
    for (int i = 0; i < lplc->N+1; i++){
        for (int j = 0; j <lplc->M+1; j++){
            if (i == 0 || i == lplc->N){
                fprintf(plotfile, "%f, %f, %f, %f\n", lplc->x[i], lplc->y[j], 0.0, 0.0);
            }else{ 
                if (j == 0 || j == lplc->M)
                    fprintf(plotfile, "%f, %f, %f, %f\n", lplc->x[i], lplc->y[j], 0.0, 0.0);
                else
                    fprintf(plotfile, "%f, %f, %f, %f\n", lplc->x[i], lplc->y[j], lplc->u[i][j], 0.0);
            }
        }
        fprintf(plotfile, "\n");        // Il faut une nouvelle ligne a chaque fois que x change
    }
    
    fclose(plotfile);

    system("gnuplot plotcom");      // plotcom contient les commandes a effectuer pour afficher les graphes issues de plot.dat

}