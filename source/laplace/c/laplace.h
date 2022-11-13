//
// Created by drn on 05/12/19.
// Resolution de l'equation de Laplace en 2D par la methode des differences finis!
//

#ifndef SKYLINE_NEUMANN_H
#define SKYLINE_NEUMANN_H

#include <stdbool.h>
#include "skyline.h"

typedef struct laplace{

    // LES ABCISSES
    double L;   // Longueur du maillage en abcisse   
    int N;      // Nombre de points du maillage en abcisse

    // LES ORDONNEES
    double H;   // Longueur du maillage en ordonees
    int M;      // Nombre de ponts du maillage en ordonnees

    // Pas d'espace en abcisse  == pas d'espace en ordonnee
    double h;

    // Taille du probleme
    int n;
    int nb_aretes;

    // Tableau des origines et des extremites de chaque aretes du malillage
    int ** extremite_arete;

    //  Differents points du maillage en abcisse et en ordonnees
    double  *x;
    double *y;

    // Solution a calculer
    double ** u;
    double *U;
    
    // Vecteur tel que (1/h**2)AU = F
    double ** f;
    double *F;

    // Matrice creuse du systeme (1/h**2)AU = F
    Skyline A;

} laplace;

//Juste les signatures des fonctions a definir dans neumann.c
void laplace_init(laplace *lplc, double L, int N, double H, int M); // Afin d'eviter les copies
void laplace_solve(laplace *lplc);
void laplace_display(laplace *lplc);

#endif //SKYLINE_NEUMANN_H
