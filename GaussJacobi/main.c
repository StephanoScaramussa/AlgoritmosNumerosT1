#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void imprimeMatriz(int linha, int coluna, float matriz[linha][coluna], char *opcao){
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

void leMatrizes(FILE* arq, int tam, int qtd, float A[tam][tam], float B[qtd][tam]){
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
void gaussJacobi(int tam, float matriz[tam][tam], float B[], float acc){
    // Monta as variaveis iniciais para resolver o problema
    float variaveis[tam];
    for(int i=0; i<tam; i++){
        variaveis[i] = B[i]/(matriz[i][i]);
    }

    float variaveisNovas[tam];
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
}

int main(int argc, char **argv){
    if (argc < 2) {
        printf("Uso: %s <nome_do_arquivo>\n", argv[0]);
        return 1;
    }

    char *nomearq = argv[1];
    printf("Lendo arquivo: %s\n", nomearq);
    
    FILE* arq = fopen(nomearq, "r");
    if (arq == NULL) {
        printf("Erro ao abrir o arquivo!\n");
        return 1;
    }

    int qtd, tam;
    float acc, digit;

    // Lê o "cabeçalho" para descobrir o tamanho das matrizes
    if (fscanf(arq, "%f", &digit) == 1) qtd = (int)digit;
    if (fscanf(arq, "%f", &digit) == 1) tam = (int)digit;
    if (fscanf(arq, "%f", &digit) == 1) acc = digit;

    float A[tam][tam];
    float B[qtd][tam];

    leMatrizes(arq, tam, qtd, A, B);
    imprimeMatriz(tam, tam, A, "A");

    fclose(arq);


    // virar funcao
    float C[tam];

    for(int j=0; j<tam; j++){
        C[j]=B[0][j];
    }
    // ate aqui

    gaussJacobi(tam, A, C, acc);


    return 0;
}