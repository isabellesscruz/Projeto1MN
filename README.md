# Métodos Numéricos para Sistemas Lineares
Este repositório contém implementações em Python de diversos métodos numéricos para resolver sistemas de equações lineares, Ax=b. Ele explora tanto as técnicas diretas de eliminação de Gauss (com e sem pivotamento) quanto métodos iterativos (Jacobi, Gauss-Seidel, SOR), além de uma implementação para cálculo de inversa de matrizes usando decomposição LU. O foco principal é demonstrar o comportamento desses algoritmos, especialmente ao lidar com matrizes mal-condicionadas, como a Matriz de Hilbert.

---

## 🏛️ Estrutura do Repositório

O repositório é composto por dois arquivos principais:

1.  `Metodos_sistemas.py`: Implementa os métodos de Eliminação de Gauss (simples, com pivotamento parcial, com escala e total) e os métodos iterativos (Jacobi, Gauss-Seidel, SOR) para resolver $Ax=b$. Inclui um gerador de Matrizes de Hilbert para testes.
2.  `inversa_lu.py`: Implementa a decomposição LU e a utiliza para calcular a inversa de uma matriz arbitrária, além de resolver sistemas triangulares.

---

## 1. `Metodos_sistemas.py`

Este arquivo é dedicado à exploração de diferentes abordagens para a resolução de sistemas de equações lineares, $Ax=b$. Ele é particularmente útil para entender como a **estabilidade numérica** e a **precisão** dos algoritmos são afetadas por matrizes com características desafiadoras, como a **Matriz de Hilbert**.

### 📦 Funcionalidades

* **Geração de Matrizes de Hilbert:** A função `gerar_sistema_hilbert(n)` cria uma Matriz de Hilbert de ordem `n` e um vetor `b` correspondente, onde a solução exata $x$ é um vetor de uns. Isso permite avaliar a precisão dos métodos de forma objetiva.
* **Métodos de Eliminação de Gauss:**
    * `eliminacaoGauss`: Implementação básica da Eliminação de Gauss (sem pivotamento).
    * `eliminacaoGauss_parcial`: Eliminação de Gauss com **pivotamento parcial**, que seleciona o maior elemento em valor absoluto na coluna para ser o pivô, melhorando a estabilidade.
    * `eliminacaoGauss_escala`: Eliminação de Gauss com **pivotamento parcial com escala**, aprimorando a escolha do pivô ao considerar a magnitude relativa dos elementos nas linhas.
    * `eliminacaoGauss_total`: Eliminação de Gauss com **pivotamento total**, que busca o maior elemento em valor absoluto em toda a submatriz restante, oferecendo a maior estabilidade numérica, mas com maior custo computacional.
* **Métodos Iterativos:**
    * `jacobi`: Implementa o **Método de Jacobi**.
    * `gaussSeidel`: Implementa o **Método de Gauss-Seidel**, que geralmente converge mais rápido que Jacobi.
    * `sobre_relaxacao`: Implementa o **Método de Sobre-Relaxamento (SOR)**, uma aceleração do Gauss-Seidel, que pode otimizar a taxa de convergência com um fator de relaxamento ($\omega$).
* **Avaliação de Erro:** A função `calcular_erro_verdadeiro` utiliza `numpy.linalg.solve` para obter a solução exata e calcula a norma da diferença entre a solução aproximada e a exata, fornecendo uma métrica clara da precisão.

### Como Funciona com a Matriz de Hilbert

A **Matriz de Hilbert** é notoriamente **mal-condicionada**, o que significa que pequenos erros (como os de arredondamento) podem levar a grandes desvios na solução.

* **Métodos Diretos:** Embora as variações da Eliminação de Gauss com pivotamento (parcial, com escala, total) ofereçam a maior estabilidade numérica para esta matriz, sua **precisão final é limitada pela ordem da matriz e pela precisão de ponto flutuante da máquina**. Para ordens maiores (e.g., $n=15$), a solução pode conter **erros significativos**, mesmo com pivotamento total, devido ao condicionamento extremo da Matriz de Hilbert. A Eliminação de Gauss sem pivotamento é **altamente instável** e resulta em soluções incorretas.

* **Métodos Iterativos:** Jacobi, Gauss-Seidel e SOR, geralmente **não convergem para a solução correta** ou produzem **erros significativos e inaceitáveis** em suas aproximações finais quando aplicados à Matriz de Hilbert. Isso ocorre porque a Matriz de Hilbert não satisfaz as condições de dominância diagonal (ou outras propriedades de convergência) exigidas por esses métodos.

### 🚀 Como usar

Para executar e testar os métodos:

1.  Certifique-se de ter `numpy` instalado (`pip install numpy`).
2.  Execute o arquivo: `python metodos_sistemas_lineares.py`
3.  Modifique as variáveis `n` (ordem da matriz) e `metodo` na seção `if __name__ == "__main__":` para testar diferentes configurações.

---

### 🔧 Parâmetros

- `n`: define a ordem da matriz de Hilbert (ex: 3, 5, 8 etc.)
- `metodo`: define qual método será utilizado para resolver o sistema

Métodos disponíveis:

- `"simples"` — Eliminação de Gauss sem pivoteamento
- `"parcial"` — Pivoteamento parcial
- `"escala"` — Pivoteamento por escala
- `"total"` — Pivoteamento total
- `"jacobi"` — Método de Jacobi (iterativo)
- `"seidel"` — Método de Gauss-Seidel (iterativo)
- `"sor"` — Método SOR com `w = 1.1` (iterativo)

---

### 📈 Exemplo de execução

```bash
python main.py
```

Saída esperada:

```
Solução do sistema:
x[0] = ...
x[1] = ...
...
Erro verdadeiro (norma da diferença): 3.1024347711e-07
```

---

### 📌 Observações

- Matrizes de Hilbert são mal-condicionadas, por isso métodos iterativos podem divergir.
- O código possui verificação automática de divergência e interrupção do Jacobi em caso de divisão por zero.
- O erro verdadeiro é calculado em relação à solução exata obtida com `numpy.linalg.solve`.

---


## 2. `inversa_lu.py`

Este arquivo implementa o método da **Decomposição LU (Lower-Upper)** para calcular a inversa de uma matriz. A decomposição LU é uma técnica fundamental em álgebra linear numérica, que fatoriza uma matriz $A$ em um produto de duas matrizes triangulares: uma matriz triangular inferior $L$ e uma matriz triangular superior $U$ ($A=LU$).

### Funcionalidades

* **Leitura de Matriz:** A função `ler_matriz()` permite que o usuário insira os elementos de uma matriz quadrada interativamente.
* **Decomposição LU:** A função `decomposicao_lu(A)` realiza a fatoração de uma matriz `A` em suas componentes $L$ e $U$. **Importante:** Esta implementação **não utiliza pivotamento**, o que significa que pode falhar se encontrar um elemento diagonal zero durante o processo de decomposição.
* **Resolução de Sistemas Triangulares:**
    * `resolver_sistema_triangular_inferior(L, b)`: Resolve sistemas do tipo $Ly=b$ (substituição progressiva).
    * `resolver_sistema_triangular_superior(U, y)`: Resolve sistemas do tipo $Ux=y$ (substituição regressiva).
* **Cálculo da Inversa por LU:** A função `calcular_inversa_lu(A)` utiliza a decomposição LU para encontrar a inversa de $A$. A ideia é resolver $n$ sistemas lineares, onde cada coluna da matriz identidade é usada como vetor `b` e a solução é uma coluna da inversa de $A$.
* **Verificação:** O programa calcula o produto $A \cdot A^{-1}$ para verificar se o resultado é uma matriz identidade, validando a precisão da inversa calculada (considerando a precisão de ponto flutuante).

###  ⚙️ Como Funciona

1.  O usuário insere os elementos de uma matriz quadrada.
2.  A matriz é decomposta em $L$ e $U$. Se um pivô zero for encontrado, a decomposição falha.
3.  Para cada coluna da matriz identidade ($e_i$), dois sistemas triangulares são resolvidos: $Ly = e_i$ (para encontrar $y$) e $Ux = y$ (para encontrar a coluna $i$ da inversa).
4.  As colunas resultantes formam a matriz inversa.
5.  Uma verificação final é realizada multiplicando a matriz original pela sua inversa para confirmar a proximidade com a matriz identidade.

### 🚀 Como usar

Para executar e calcular a inversa de uma matriz:

1.  Certifique-se de ter `numpy` instalado (`pip install numpy`).
2.  Execute o arquivo: `python inversa_lu.py`
3.  Siga as instruções para inserir os elementos da matriz.

### 📌 Observações

- A decomposição LU usada aqui não realiza pivotamento, portanto o algoritmo pode falhar se encontrar um pivô nulo. Nestes casos, uma mensagem de erro será exibida.
- Pequenas diferenças de ponto flutuante são tratadas com np.allclose.

---



## 📄 Licença

Este projeto é de livre uso para fins acadêmicos e educacionais.  
Desenvolvido como parte de um estudo comparativo sobre métodos de resolução de sistemas lineares.




