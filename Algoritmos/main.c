#include "Algoritmos.h"


int main(int argc, char **argv){
    // Arquivo e Inicialização da Matriz
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

    if (fscanf(arq, "%f", &digit) == 1) qtd = (int)digit;
    if (fscanf(arq, "%f", &digit) == 1) tam = (int)digit;
    if (fscanf(arq, "%f", &digit) == 1) acc = digit;

    float **A = (float**)malloc(tam * sizeof(float*));
    for(int i = 0; i < tam; i++) A[i] = (float*)malloc(tam * sizeof(float));
    
    float **B = (float**)malloc(qtd * sizeof(float*));
    for(int i = 0; i < qtd; i++) B[i] = (float*)malloc(tam * sizeof(float));

    leMatrizes(arq, tam, qtd, A, B);
    fclose(arq);

    float *C = (float*)malloc(tam * sizeof(float));

    for(int j = 0; j < tam; j++){
        C[j] = B[0][j];
    }
    #pragma endregion
    
    relatorioDesempenho(tam, A, B, qtd, acc);

    for(int i = 0; i < tam; i++) free(A[i]);
    free(A);
    
    for(int i = 0; i < qtd; i++) free(B[i]);
    free(B);

    return 0;
}