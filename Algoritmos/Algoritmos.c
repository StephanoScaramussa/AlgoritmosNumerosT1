#include "Algoritmos.h"

void imprimeMatriz(int linha, int coluna, float **matriz, char *opcao){
    printf("matriz %s = \n", opcao);
    for(int i=0; i<linha; i++){
        for(int j=0; j<coluna; j++){
            printf("%3.1f ", matriz[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void imprimeVetor(float vetor[], int tam, char *opcao){
        printf("matriz %s = \n", opcao);
        for(int i=0; i<tam; i++){
            printf("%3.1f ", vetor[i]);
        }
        printf("\n");
        printf("\n");
}

void leMatrizes(FILE* arq, int tam, int qtd, float **A, float **B){
    float digit;
    int linha = 0, coluna = 0;

    // Continua lendo de onde o main parou (a matriz A)
    while(fscanf(arq, "%f", &digit) == 1){
        A[linha][coluna] = digit;
        coluna++;
        if(coluna == tam){
            coluna = 0; 
            linha++;
            if (linha == tam) break;
        }
    }

    linha = 0; coluna = 0;
    
    // Continua lendo a matriz B
    while(fscanf(arq, "%f", &digit) == 1){
        B[linha][coluna] = digit;
        coluna++;
        if(coluna == tam){
            coluna = 0; 
            linha++;
            if (linha == qtd) break; // Cuidado: linha < qtd mudou para linha == qtd
        }
    }
}

// void gaussSeidel(int tam, float matriz[tam][tam], float B[]){
//     // Monta as variaveis iniciais para resolver o problema
//     float variaveis[tam];
//     for(int i=0; i<tam; i++){
//         variaveis[i] = B[i]/(matriz[i][i]);
//     }
//     imprimeVetor(variaveis, tam, "variaveis");
    

//     // Itera as linhas, uma variavel por linha
//     for(int j=0; j<tam; j++){
//         // Aqui ta calculando o valor de uma variavel usando a equação doida lá usando acumulador, pq senao ia ser hard codado
//         int i= tam -1;
//         float soma = B[j];
//         printf("valor de B: %f - ", soma);

//         for(i; i>=0; i--){
//             if(i != j){
//                 printf("%f ", matriz[j][i]);
//                 soma = soma - matriz[j][i]*variaveis[i];
//             }
//         }
//         printf("/%f \n", matriz[j][j]);
//         // Divide pelo valor da variavel naquele momento
//         soma = soma/matriz[j][j];
//         variaveis[j] = soma;
//         }
    
//     printf("\n");
//     }


// Funfando
void gaussJacobi(int tam, float **matriz, float *B, float acc){
    // Monta as variaveis iniciais para resolver o problema
    float *variaveis = (float*)malloc(tam * sizeof(float));
    for(int i=0; i<tam; i++){
        variaveis[i] = B[i]/(matriz[i][i]);
    }

    float *variaveisNovas = (float*)malloc(tam * sizeof(float));
    float erroRelativo;

    do{
    // Itera as linhas, uma variavel por linha
    for(int j=0; j<tam; j++){
        // Aqui ta calculando o valor de uma variavel usando a equação doida lá usando acumulador, pq senao ia ser hard codado
        
        float soma = B[j];

        for(int i= tam -1; i>=0; i--){
            if(i != j){
                soma = soma - matriz[j][i]*variaveis[i];
            }
        }
        // Divide pelo valor da variavel naquele momento
        variaveisNovas[j] = soma/matriz[j][j];
        }
    
    float maiorDif = 0;
    float maiorVarNova = 0;

    for(int i=0; i< tam; i++){
        if(fabs(variaveisNovas[i]) > maiorVarNova){
            maiorVarNova = variaveisNovas[i]; 
        }
        float dif = fabs(variaveisNovas[i] - variaveis[i]);
        if ( dif > maiorDif){
            maiorDif = dif;
        }
    }
    
    erroRelativo = maiorDif/maiorVarNova;

    for(int i=0; i<tam; i++){
            variaveis[i] = variaveisNovas[i];
        }

    } while (erroRelativo > acc);
    imprimeVetor(variaveisNovas, tam, "Respostas");
    free(variaveis);
    free(variaveisNovas);
}

void fatoracaoLU(int tam, float **A, float **L, float **U){
    int i, j, k;
 
    for(i = 0; i < tam; i++)
        for(j = 0; j < tam; j++){
            L[i][j] = 0.0f;
            U[i][j] = 0.0f;
        }
 
    for(i = 0; i < tam; i++){
        for(j = i; j < tam; j++){
            float soma = 0.0f;
            for(k = 0; k < i; k++)
                soma += L[i][k] * U[k][j];
            U[i][j] = A[i][j] - soma;
        }
 
        L[i][i] = 1.0f;
        for(j = i + 1; j < tam; j++){
            float soma = 0.0f;
            for(k = 0; k < i; k++)
                soma += L[j][k] * U[k][i];
            L[j][i] = (A[j][i] - soma) / U[i][i];
        }
    }
}


void resolveLU(int tam, float **L, float **U, float B[], float X[]){
    float *y = (float*)malloc(tam * sizeof(float));
    int i, j;
 
    for(i = 0; i < tam; i++){
        float soma = 0.0f;
        for(j = 0; j < i; j++){
            soma += L[i][j] * y[j];
        }
        y[i] = (B[i] - soma) / L[i][i];
    }

    for(i = tam - 1; i >= 0; i--){
        float soma = 0.0f;
        for(j = i + 1; j < tam; j++){
            soma += U[i][j] * X[j];
        }
        X[i] = (y[i] - soma) / U[i][i];
    }
    
    free(y);
}