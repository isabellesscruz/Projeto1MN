# Resolução de Sistemas Lineares — Métodos Diretos e Iterativos

Este projeto tem como objetivo comparar diferentes métodos de resolução de sistemas lineares aplicados a matrizes de Hilbert, conhecidas por seu alto grau de mal condicionamento. O código implementa métodos diretos e iterativos, mede o erro verdadeiro das soluções e facilita testes com diferentes ordens de sistemas.

---

## 📦 Funcionalidades

- Geração de sistemas lineares com matriz de Hilbert (`n x n`)
- Implementação de 7 métodos de resolução:
  - Eliminação de Gauss (simples)
  - Eliminação com pivoteamento parcial
  - Eliminação com pivoteamento por escala
  - Eliminação com pivoteamento total
  - Método de Jacobi
  - Método de Gauss-Seidel
  - Método SOR (sobre-relaxação)
- Cálculo do **erro verdadeiro** comparando a solução obtida com a fornecida por `numpy.linalg.solve`
- Detecção de divergência ou falha nos métodos iterativos

---

## ⚙️ Requisitos

- Python 3.8+
- Biblioteca NumPy

Instale o NumPy com:

```bash
pip install numpy
```

---

## 🚀 Como usar

No final do arquivo Python, o código possui um bloco de execução principal:

```python
if __name__ == "__main__":
    n = 5
    metodo = "seidel"
    resolver_sistema(n, metodo)
```

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

## 📈 Exemplo de execução

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

## 📌 Observações

- Matrizes de Hilbert são mal-condicionadas, por isso métodos iterativos podem divergir.
- O código possui verificação automática de divergência e interrupção do Jacobi em caso de divisão por zero.
- O erro verdadeiro é calculado em relação à solução exata obtida com `numpy.linalg.solve`.

---

## 📄 Licença

Este projeto é de livre uso para fins acadêmicos e educacionais.  
Desenvolvido como parte de um estudo comparativo sobre métodos de resolução de sistemas lineares.
