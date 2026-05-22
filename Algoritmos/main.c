#include algoritmos.h

int main(int argc, char **argv){

    clock_t start, end;
    start = clock();

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

    end = clock();
    printf("\nTempo de execucao: %5.2f seg.\n", ((double)(end - start)) / CLOCKS_PER_SEC);
    return 0;
}