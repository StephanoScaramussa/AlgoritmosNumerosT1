#ifndef Algoritmos
#define Algoritmos
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

void imprimeMatriz(int linha, int coluna, float **matriz, char *opcao);
void imprimeVetor(float vetor[], int tam, char *opcao);
void leMatrizes(FILE* arq, int tam, int qtd, float **A, float **B);
void gaussJacobi(int tam, float **matriz, float *B, float acc);
void fatoracaoLU(int tam, float **A, float **L, float **U);
void resolveLU(int tam, float **L, float **U, float *B, float X[]);
#endif