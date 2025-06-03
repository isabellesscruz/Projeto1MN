import numpy as np
import copy

def gerar_sistema_hilbert(n):
    H = np.array([[1 / (i + j + 1) for j in range(n)] for i in range(n)], dtype=float)
    b = np.sum(H, axis=1)
    return H.tolist(), b.tolist()

def Resultados(resultado):
    print("Solução do sistema:")
    for i, x in enumerate(resultado):
        print(f"x[{i}] = {x:.10f}")

def safe_norm(vetor):
    return sum(x**2 for x in vetor) ** 0.5

def eliminacaoGauss(matriz, extensao):
    n = len(matriz)
    for i in range(n - 1):
        for j in range(i + 1, n):
            m = matriz[j][i] / matriz[i][i]
            matriz[j][i] = 0
            for k in range(i + 1, n):
                matriz[j][k] -= m * matriz[i][k]
            extensao[j] -= m * extensao[i]
    resultado = [0] * n
    resultado[n - 1] = extensao[n - 1] / matriz[n - 1][n - 1]
    for i in range(n - 2, -1, -1):
        soma = sum(matriz[i][j] * resultado[j] for j in range(i + 1, n))
        resultado[i] = (extensao[i] - soma) / matriz[i][i]
    Resultados(resultado)
    return resultado

def eliminacaoGauss_parcial(matriz, extensao):
    n = len(matriz)
    for i in range(n):
        max_linha = max(range(i, n), key=lambda r: abs(matriz[r][i]))
        if i != max_linha:
            matriz[i], matriz[max_linha] = matriz[max_linha], matriz[i]
            extensao[i], extensao[max_linha] = extensao[max_linha], extensao[i]
        for j in range(i + 1, n):
            m = matriz[j][i] / matriz[i][i]
            matriz[j][i] = 0
            for k in range(i + 1, n):
                matriz[j][k] -= m * matriz[i][k]
            extensao[j] -= m * extensao[i]
    resultado = [0] * n
    resultado[n - 1] = extensao[n - 1] / matriz[n - 1][n - 1]
    for i in range(n - 2, -1, -1):
        soma = sum(matriz[i][j] * resultado[j] for j in range(i + 1, n))
        resultado[i] = (extensao[i] - soma) / matriz[i][i]
    Resultados(resultado)
    return resultado

def eliminacaoGauss_escala(matriz, extensao):
    n = len(matriz)
    escala = [max(abs(x) for x in linha) for linha in matriz]
    for i in range(n):
        razoes = [abs(matriz[r][i]) / escala[r] for r in range(i, n)]
        max_r = razoes.index(max(razoes)) + i
        if i != max_r:
            matriz[i], matriz[max_r] = matriz[max_r], matriz[i]
            extensao[i], extensao[max_r] = extensao[max_r], extensao[i]
            escala[i], escala[max_r] = escala[max_r], escala[i]
        for j in range(i + 1, n):
            m = matriz[j][i] / matriz[i][i]
            matriz[j][i] = 0
            for k in range(i + 1, n):
                matriz[j][k] -= m * matriz[i][k]
            extensao[j] -= m * extensao[i]
    resultado = [0] * n
    resultado[n - 1] = extensao[n - 1] / matriz[n - 1][n - 1]
    for i in range(n - 2, -1, -1):
        soma = sum(matriz[i][j] * resultado[j] for j in range(i + 1, n))
        resultado[i] = (extensao[i] - soma) / matriz[i][i]
    Resultados(resultado)
    return resultado

def eliminacaoGauss_total(matriz, extensao):
    n = len(matriz)
    pos = list(range(n))
    for i in range(n):
        max_val = 0
        max_row, max_col = i, i
        for r in range(i, n):
            for c in range(i, n):
                if abs(matriz[r][c]) > abs(max_val):
                    max_val = matriz[r][c]
                    max_row, max_col = r, c
        if i != max_row:
            matriz[i], matriz[max_row] = matriz[max_row], matriz[i]
            extensao[i], extensao[max_row] = extensao[max_row], extensao[i]
        if i != max_col:
            for linha in matriz:
                linha[i], linha[max_col] = linha[max_col], linha[i]
            pos[i], pos[max_col] = pos[max_col], pos[i]
        for j in range(i + 1, n):
            m = matriz[j][i] / matriz[i][i]
            matriz[j][i] = 0
            for k in range(i + 1, n):
                matriz[j][k] -= m * matriz[i][k]
            extensao[j] -= m * extensao[i]
    resultado = [0] * n
    resultado[n - 1] = extensao[n - 1] / matriz[n - 1][n - 1]
    for i in range(n - 2, -1, -1):
        soma = sum(matriz[i][j] * resultado[j] for j in range(i + 1, n))
        resultado[i] = (extensao[i] - soma) / matriz[i][i]
    resultado_final = [0] * n
    for i, posicao in enumerate(pos):
        resultado_final[posicao] = resultado[i]
    Resultados(resultado_final)
    return resultado_final

def jacobi(matriz, extensao):
    n = len(matriz)
    resultado = [0.0] * n
    x = copy.deepcopy(resultado)
    for _ in range(10000):
        for i in range(n):
            soma = sum(matriz[i][j] * resultado[j] for j in range(n) if j != i)
            try:
                x[i] = (extensao[i] - soma) / matriz[i][i]
            except ZeroDivisionError:
                print("Divisão por zero detectada. Método interrompido.")
                return
        # Verifica divergência (valores muito grandes)
        if any(abs(val) > 1e10 for val in x):
            print("Divergência detectada! Interrompendo o método de Jacobi.")
            return
        if abs(safe_norm(x) - safe_norm(resultado)) < 0.01:
            break
        resultado = copy.deepcopy(x)
    Resultados(resultado)
    return resultado

def gaussSeidel(matriz, extensao):
    n = len(matriz)
    resultado = [0] * n
    x = copy.deepcopy(resultado)
    for _ in range(10000):
        for i in range(n):
            soma = sum(matriz[i][j] * resultado[j] for j in range(n) if j != i)
            x[i] = (extensao[i] - soma) / matriz[i][i]
            resultado[i] = x[i]
        if abs(safe_norm(x) - safe_norm(resultado)) < 0.01:
            break
    Resultados(resultado)
    return resultado

def sobre_relaxacao(matriz, extensao, w=1.1):
    n = len(matriz)
    resultado = [0] * n
    x = copy.deepcopy(resultado)
    for _ in range(10000):
        for i in range(n):
            soma = sum(matriz[i][j] * resultado[j] for j in range(n) if j != i)
            x[i] = (1 - w) * resultado[i] + (w * (extensao[i] - soma) / matriz[i][i])
            resultado[i] = x[i]
        if abs(safe_norm(x) - safe_norm(resultado)) < 0.01:
            break
    Resultados(resultado)
    return resultado

def calcular_erro_verdadeiro(sol_aproximada, matriz, extensao):
    matriz_np = np.array(matriz)
    extensao_np = np.array(extensao)
    sol_exata = np.linalg.solve(matriz_np, extensao_np)
    erro = np.linalg.norm(np.array(sol_aproximada) - sol_exata)
    print(f"Erro verdadeiro (norma da diferença): {erro:.10e}\n")


def resolver_sistema(n, metodo):
    matriz, extensao = gerar_sistema_hilbert(n)
    resultado = None
    if metodo == "simples":
        resultado = eliminacaoGauss(copy.deepcopy(matriz), copy.deepcopy(extensao))
    elif metodo == "parcial":
        resultado = eliminacaoGauss_parcial(copy.deepcopy(matriz), copy.deepcopy(extensao))
    elif metodo == "escala":
        resultado = eliminacaoGauss_escala(copy.deepcopy(matriz), copy.deepcopy(extensao))
    elif metodo == "total":
        resultado = eliminacaoGauss_total(copy.deepcopy(matriz), copy.deepcopy(extensao))
    elif metodo == "jacobi":
        resultado = jacobi(copy.deepcopy(matriz), copy.deepcopy(extensao))
    elif metodo == "seidel":
        resultado = gaussSeidel(copy.deepcopy(matriz), copy.deepcopy(extensao))
    elif metodo == "sor":
        resultado = sobre_relaxacao(copy.deepcopy(matriz), copy.deepcopy(extensao))
    else:
        print("Método inválido.")
        return
    calcular_erro_verdadeiro(resultado, matriz, extensao)

if __name__ == "__main__":
    n = 15
    metodo = "sor"
    resolver_sistema(n, metodo)
