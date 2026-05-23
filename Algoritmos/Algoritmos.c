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

void gaussSeidel(int tam, float **matriz, float B[], float acc){
    // Monta as variaveis iniciais para resolver o problema
    float variaveis[tam];
    

    for(int i=0; i<tam; i++){
        if(fabs(matriz[i][i]) < 1e-10){
            printf("Diagonal 0 ou proximo\n");
        }
        variaveis[i] = B[i]/(matriz[i][i]);
    }

    float erroRelativo;
    int iteracao = 0;
    int maxIteracoes = 10000;

    do{
    float maiorDif = 0, maiorVarNova = 0;
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
        if(fabs(matriz[j][j]) < 1e-10){
            printf("Divisão por zero na linha %d\n", j);
            return;
        }
        float valorAntigo = variaveis[j];
        variaveis[j] = soma/matriz[j][j];

        // Detecta valores infinitos ou muito grandes
        if(!isfinite(variaveis[j])){
            printf("Valor infinito detectado na iteracao %d, variavel %d\n", iteracao, j);
            return;
        }

        if(fabs(variaveis[j]) > maiorVarNova){
            maiorVarNova = fabs(variaveis[j]); 
        }
        float dif = fabs(variaveis[j] - valorAntigo);
        if ( dif > maiorDif){
            maiorDif = dif;
        }
        }
    
    if(maiorVarNova < 1e-10){
        printf("Nenhuma permuta possivel\n");
        erroRelativo = 0;
    } else {
        erroRelativo = maiorDif/maiorVarNova;
    }
    
    iteracao++;

    } while (erroRelativo > acc && iteracao < maxIteracoes);
    
    if(iteracao >= maxIteracoes){
        printf("Seidel atingiu limite de %d iteracoes e nao converge\n", maxIteracoes);
    }
    
    if(tam<=200){
        imprimeVetor(variaveis, tam, "Respostas");

    }
}

// Funfando 
void gaussJacobi(int tam, float **matriz, float B[], float acc){
    // Monta as variaveis iniciais para resolver o problema
    float *variaveis = (float*)malloc(tam * sizeof(float));
    for(int i=0; i<tam; i++){
        if(fabs(matriz[i][i]) < 1e-10){
            printf("Valor 0 ou proximo\n");
        }
        variaveis[i] = B[i]/(matriz[i][i]);
    }

    float variaveisNovas[tam], erroRelativo;
    int iteracao = 0;
    int maxIteracoes = 10000;

    do{
    float maiorDif = 0, maiorVarNova = 0;
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
        if(fabs(matriz[j][j]) < 1e-10){
            printf("Divisão por zero na linha %d\n", j);
            free(variaveis);
            return;
        }
        variaveisNovas[j] = soma/matriz[j][j];

        if(!isfinite(variaveisNovas[j])){
            printf("Valor infinito na iteracao %d, variavel %d\n", iteracao, j);
            free(variaveis);
            return;
        }

        if(fabs(variaveisNovas[j]) > maiorVarNova){
            maiorVarNova = fabs(variaveisNovas[j]); 
        }
        float dif = fabs(variaveisNovas[j] - variaveis[j]);
        if ( dif > maiorDif){
            maiorDif = dif;
        }
        }
    
    if(maiorVarNova < 1e-10){
        printf("Nenhuma permuta possivel\n");
        erroRelativo = 0;
    } else {
        erroRelativo = maiorDif/maiorVarNova;
    }

    for(int i=0; i<tam; i++){
            variaveis[i] = variaveisNovas[i];
        }
    
    iteracao++;

    } while (erroRelativo > acc && iteracao < maxIteracoes);
    
    if(iteracao >= maxIteracoes){
        printf("Jacobi atingiu limite de %d iteracoes e nao converge\n", maxIteracoes);
    }
    
    if(tam<=200){
        imprimeVetor(variaveisNovas, tam, "Respostas");
    }
    free(variaveis);
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


void resolveLU(int tam, float **L, float **U, float *B, float X[]){
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

void escalonamentoGauss(int tam, float **A, float *B, float *X){
    float **M = (float**)malloc(tam * sizeof(float*));
    for(int i = 0; i < tam; i++){
        M[i] = (float*)malloc(tam * sizeof(float));
        for(int j = 0; j < tam; j++){
            M[i][j] = A[i][j];
        }
    }
    
    float *b = (float*)malloc(tam * sizeof(float));
    for(int i = 0; i < tam; i++){
        b[i] = B[i];
    }
    
    for(int k = 0; k < tam - 1; k++){
        int pivo = k;
        for(int i = k + 1; i < tam; i++){
            if(fabs(M[i][k]) > fabs(M[pivo][k])){
                pivo = i;
            }
        }
        
        if(fabs(M[pivo][k]) < 1e-10){
            printf("Aviso: pivô é zero na coluna %d\n", k);
            int encontrado = 0;
            for(int i = k + 1; i < tam; i++){
                if(fabs(M[i][k]) > 1e-10){
                    pivo = i;
                    encontrado = 1;
                    break;
                }
            }
            if(!encontrado){
                printf("Nenhum pivô não-nulo na coluna %d\n", k);
                for(int i = 0; i < tam; i++) free(M[i]);
                free(M);
                free(b);
                return;
            }
        }

        // Se pivô é zero, permuta linha pra usar um pivô diferente de zero
        
        if(pivo != k){
            float *temp = M[k];
            M[k] = M[pivo];
            M[pivo] = temp;
            
            float temp_b = b[k];
            b[k] = b[pivo];
            b[pivo] = temp_b;
        }
        
        for(int i = k + 1; i < tam; i++){
            float fator = M[i][k] / M[k][k];
            
            for(int j = k; j < tam; j++){
                M[i][j] = M[i][j] - fator * M[k][j];
            }
            b[i] = b[i] - fator * b[k];
        }
    }
    
    for(int i = tam - 1; i >= 0; i--){
        float soma = 0.0f;
        for(int j = i + 1; j < tam; j++){
            soma += M[i][j] * X[j];
        }
        X[i] = (b[i] - soma) / M[i][i];
    }
    
    for(int i = 0; i < tam; i++){
        free(M[i]);
    }
    free(M);
    free(b);
}
void relatorioDesempenho(int tam, float **A, float **B, int qtd, float acc) {
    float *C = (float*)malloc(tam * sizeof(float));
    clock_t start, end;

    // GAUSS-SEIDEL PARA CADA MATRIZ
    printf("\n--------------GAUSS SEIDEL--------------");
    for(int s = 0; s < qtd; s++) {
        for(int j = 0; j < tam; j++) C[j] = B[s][j]; // Pega o sistema B atual

        printf("\n[Sistema %d] ", s + 1);
        start = clock();
        preparaDiagonal(tam,A,C);
        gaussSeidel(tam, A, C, acc);
        end = clock();
        printf("Tempo: %5.6f seg.\n", ((double)(end - start)) / CLOCKS_PER_SEC);
    }

    // GAUSS-JACOBI PARA CADA B
    printf("\n--------------GAUSS JACOBI--------------");
    for(int s = 0; s < qtd; s++) {
        for(int j = 0; j < tam; j++) C[j] = B[s][j]; // Pega o sistema B atual

        printf("\n[Sistema %d] ", s + 1);
        start = clock();
        preparaDiagonal(tam,A,C);
        gaussJacobi(tam, A, C, acc);
        end = clock();
        printf("Tempo: %5.6f seg.\n", ((double)(end - start)) / CLOCKS_PER_SEC);
    }

    //ESCALONAMENTO GAUSS
    printf("\n--------------ESCALONAMENTO GAUSS--------------");
    for(int s = 0; s < qtd; s++) {
        for(int j = 0; j < tam; j++) C[j] = B[s][j]; // Pega o sistema B atual

        float *X = (float*)malloc(tam * sizeof(float));
        printf("\n[Sistema %d] ", s + 1);
        start = clock();
        escalonamentoGauss(tam, A, C, X);
        end = clock();
        if(tam <= 200){
            imprimeVetor(X, tam, "Respostas");
        }
        printf("Tempo: %5.6f seg.\n", ((double)(end - start)) / CLOCKS_PER_SEC);
        free(X);
    }
    
    //FATORAÇÃO LU
    #pragma region
    float **L = (float**)malloc(tam * sizeof(float*));
    float **U = (float**)malloc(tam * sizeof(float*));
    for(int i = 0; i < tam; i++) {
        L[i] = (float*)malloc(tam * sizeof(float));
        U[i] = (float*)malloc(tam * sizeof(float));
    }

    start = clock();
    fatoracaoLU(tam, A, L, U);
    end = clock();
    printf("\n--------------FATORACAO LU (matriz A) = %5.6f seg.--------------\n", ((double)(end - start)) / CLOCKS_PER_SEC);

    // RESOLUÇÃO DO LU PARA CADA MATRIZ B
    for(int s = 0; s < qtd; s++){
        for(int j = 0; j < tam; j++) C[j] = B[s][j];

        float *X = (float*)malloc(tam * sizeof(float));

        start = clock();
        resolveLU(tam, L, U, C, X);
        end = clock();

        printf("[Sistema %d] ", s + 1);
        if(tam <= 200){
            imprimeVetor(X, tam, "Respostas");
        }
        printf("Tempo: %d (LU): %5.6f seg.\n\n",
               s + 1, ((double)(end - start)) / CLOCKS_PER_SEC);
        
        free(X);
    }

    for(int i = 0; i < tam; i++) {
        free(L[i]);
        free(U[i]);
    }
    free(L);
    free(U);
    free(C);
    #pragma endregion    
}

void preparaDiagonal(int tam, float **matriz, float *B) {
    for (int i = 0; i < tam; i++) {
        int linhaMaior = i;
        float maiorValor = fabs(matriz[i][i]);
        
        for (int k = i + 1; k < tam; k++) {
            if (fabs(matriz[k][i]) > maiorValor) {
                maiorValor = fabs(matriz[k][i]);
                linhaMaior = k;
            }
        }
        
        if (linhaMaior != i) {
            float *temp = matriz[i];
            matriz[i] = matriz[linhaMaior];
            matriz[linhaMaior] = temp;
            
            float tempB = B[i];
            B[i] = B[linhaMaior];
            B[linhaMaior] = tempB;
        }
    }
}