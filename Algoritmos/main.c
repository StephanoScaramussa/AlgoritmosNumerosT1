#include "Algoritmos.h"

int main(int argc, char **argv){
    //Arquivo e Inicialização da matriz
    #pragma region 
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
    fclose(arq);
    //imprimeMatriz(tam, tam, A, "A");

    // virar funcao
    float C[tam];

    for(int j=0; j<tam; j++){
        C[j]=B[0][j];
    }
    // ate aqui
    #pragma endregion
    
    // Gauss Seidel
    #pragma region
    clock_t iniSeidel, fimSeidel;
    iniSeidel = clock();
    gaussSeidel(tam, A, C, acc);
    fimSeidel = clock();

    printf("\nGauss Seidel: %5.2f seg.\n", ((double)(iniSeidel - fimSeidel)) / CLOCKS_PER_SEC);
    #pragma endregion 

    // Gauss Jacobi
    #pragma region
    clock_t start, end;
    start = clock();

    gaussJacobi(tam, A, C, acc);

    end = clock();
    printf("\nGauss Jac: %5.2f seg.\n", ((double)(end - start)) / CLOCKS_PER_SEC);
    
    #pragma endregion

    // LU
    #pragma region
    //Fatoração LU feita uma única vez para matriz A
    float L[tam][tam];
    float U[tam][tam];

    clock_t lu_start = clock();
    fatoracaoLU(tam, A, L, U);
    clock_t lu_end = clock();

    //imprimeMatriz(tam, tam, L, "L");
    //imprimeMatriz(tam, tam, U, "U");
    printf("Tempo fatoracao LU (matriz A): %5.6f seg.\n\n", ((double)(lu_end - lu_start)) / CLOCKS_PER_SEC);

    //Calculo de tempo da resolução das matrizes L e U ja feitas na função anterior
    for(int s = 0; s < qtd; s++){
        float C[tam];
        for(int j = 0; j < tam; j++)
            C[j] = B[s][j];

        float X[tam];

        clock_t s_start = clock();
        resolveLU(tam, L, U, C, X);
        clock_t s_end = clock();

        printf("=== Fatoracao LU - Sistema %d ===\n", s + 1);
        imprimeVetor(X, tam, "X (LU)");
        printf("Tempo resolucao sistema %d (LU): %5.6f seg.\n\n",
               s + 1, ((double)(s_end - s_start)) / CLOCKS_PER_SEC);
    }

    #pragma endregion


     return 0;
}