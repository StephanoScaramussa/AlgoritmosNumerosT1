#ifndef Algoritmos
#define Algoritmos
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

void imprimeMatriz(int linha, int coluna, float matriz[linha][coluna], char *opcao);
void imprimeVetor(float vetor[], int tam, char *opcao);
void leMatrizes(FILE* arq, int tam, int qtd, float A[tam][tam], float B[qtd][tam]);
void gaussJacobi(int tam, float matriz[tam][tam], float B[], float acc);
#endif