import numpy as np

def ler_matriz():
    while True:
        try:
            n = int(input("Digite o número de linhas/colunas da matriz: "))
            if n <= 0:
                print("O número de linhas/colunas deve ser positivo.")
                continue
            break
        except ValueError:
            print("Entrada inválida. Digite um número inteiro.")

    matriz = np.zeros((n, n))
    print(f"Digite os elementos da matriz {n}x{n}, termo a termo:")
    for i in range(n):
        for j in range(n):
            while True:
                try:
                    matriz[i, j] = float(input(f"Elemento [{i+1}][{j+1}]: "))
                    break
                except ValueError:
                    print("Entrada inválida. Digite um número.")
    return matriz

def decomposicao_lu(A):
    n = A.shape[0]
    L = np.eye(n)  # Matriz L inicializada como identidade
    U = A.copy()  # Matriz U inicializada como uma cópia de A

    for k in range(n - 1):
        if U[k, k] == 0:
            raise ValueError("Decomposição LU sem pivotamento não é possível (elemento diagonal zero).")
        for i in range(k + 1, n):
            multiplicador = U[i, k] / U[k, k]
            L[i, k] = multiplicador
            for j in range(k, n):
                U[i, j] -= multiplicador * U[k, j]
    return L, U

def resolver_sistema_triangular_inferior(L, b):
    n = L.shape[0]
    y = np.zeros(n)
    for i in range(n):
        soma = 0
        for j in range(i):
            soma += L[i, j] * y[j]
        y[i] = (b[i] - soma) / L[i, i]
    return y

def resolver_sistema_triangular_superior(U, y):
    n = U.shape[0]
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        soma = 0
        for j in range(i + 1, n):
            soma += U[i, j] * x[j]
        x[i] = (y[i] - soma) / U[i, i]
    return x

def calcular_inversa_lu(A):
    n = A.shape[0]
    identidade = np.eye(n)
    inversa = np.zeros((n, n))

    try:
        L, U = decomposicao_lu(A)
    except ValueError as e:
        print(f"Erro na decomposição LU: {e}")
        return None

    for i in range(n):
        # Resolve Ly = e_i (coluna da identidade)
        y = resolver_sistema_triangular_inferior(L, identidade[:, i])
        # Resolve Ux = y
        x = resolver_sistema_triangular_superior(U, y)
        inversa[:, i] = x
    return inversa

def main():
    print("--- Cálculo da Inversa de uma Matriz por Decomposição LU ---")
    A = ler_matriz()

    if A.shape[0] != A.shape[1]:
        print("\nA matriz digitada não é quadrada. Não é possível calcular a inversa.")
        return

    print("\nMatriz Original:")
    print(A)

    inversa_A = calcular_inversa_lu(A)

    if inversa_A is not None:
        print("\nMatriz Inversa (A^-1):")
        print(inversa_A)

        # Verificar se A * A^-1 é a matriz identidade
        produto = np.dot(A, inversa_A)
        print("\nVerificação: A * A^-1 (deve ser a matriz identidade):")
        # Arredondar para lidar com pequenas imprecisões de ponto flutuante
        print(np.round(produto, decimals=6))

        # Comparar com a matriz identidade
        identidade_esperada = np.eye(A.shape[0])
        if np.allclose(produto, identidade_esperada):
            print("\nO produto A * A^-1 é a matriz identidade. A inversa foi calculada corretamente.")
        else:
            print("\nO produto A * A^-1 NÃO é a matriz identidade com precisão total. Podem haver erros.")

if __name__ == "__main__":
    main()
