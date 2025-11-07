# SUMMITVII (C++ + OpenMP, Linux)

**Archivo principal:** `SUMMITVII.cpp`

Implementa:

- **ShearSort**: secuencial, paralelo y paralelo con **transposición**  
- **Búsquedas**: **BSA** (binaria secuencial), **P-BSA** (binaria paralela por subsecuencias) y **Lineal** (baseline para no ordenados)  
- **Métricas**: tiempos (wall-time) y **speedup** automáticos

---

## 1) Compilación (Linux, lo más simple)


Compila en (GCC + OpenMP):

```bash
g++ SUMMITVII.cpp -fopenmp -O2 -o app
```

> Si tu GCC exige estándar explícito, usa:
>
> ```bash
> g++ SUMMITVII.cpp -fopenmp -O2 -std=c++17 -o app
> ```

---

## 2) Ejecución

```bash
./app [n] [Nvec] [threads] [seed]
```

- `n`: tamaño de la **matriz n×n** para ShearSort (default: `512`)
- `Nvec`: tamaño del **vector** para búsquedas (default: `1000000`)
- `threads`: número de hilos OpenMP (default: **máximo disponible**)
- `seed`: semilla RNG para reproducibilidad (default: `42`)

**Ejemplos:**
```bash
# Valores por defecto
./app

# Matriz 1024×1024, vector 2,000,000, 8 hilos, semilla 123
./app 1024 2000000 8 123
```

---

## 3) Qué funcionalidades corre y cómo se muestran en la salida

Al iniciar, el programa imprime la **configuración** usada:

```
Config: n=<...>  Nvec=<...>  threads=<...>  seed=<...>
```

### 3.1 ShearSort (3 variantes)

El programa ejecuta **las tres** variantes y para cada una imprime **tiempo** y validación:

1) **Secuencial**  
   - Rondas: `⌊log2(n)⌋+1`.  
   - Por ronda: filas en “serpiente” (pares ascendente, impares descendente) y columnas ascendente.  
   - **Salida**:
     ```
     [ShearSort] Secuencial    : <tiempo> s  | snake-ok=<true/false>
     ```
   - `snake-ok=true` → confirma que el orden final global (recorriendo en “serpiente”) es correcto.

2) **Paralelo (OpenMP)**  
   - Paraleliza filas y columnas con `#pragma omp parallel for`.  
   - **Salida**:
     ```
     [ShearSort] Paralelo      : <tiempo> s  | snake-ok=<...>  | speedup=<t_seq/t_par>
     ```
   - `speedup` compara contra la versión **secuencial**.

3) **Paralelo con Transposición**  
   - Alterna: ordenar **filas** → **transponer** → ordenar **filas** (equivale a columnas) → **destransponer**.  
   - Ventaja: accesos más contiguos (mejor caché).  
   - **Salida**:
     ```
     [ShearSort] Con transpose : <tiempo> s  | snake-ok=<...>  | speedup(vs seq)=<t_seq/t_transpose>
     ```

### 3.2 Búsquedas (ordenados y no ordenados)

El programa genera vectores **ordenados** y **no ordenados** y corre:

1) **BSA — Búsqueda Binaria (secuencial, ordenado)**  
   - **Salida**:
     ```
     [BSA] Secuencial (ordenado)   : <tiempo> s  | found=1
     ```
   - `found=1` → se encontró el valor (ver sección 4 para cómo se elige).

2) **P-BSA — Búsqueda Binaria Paralela por subsecuencias (ordenado)**  
   - Divide el rango en `P` subrangos (P = `threads`); cada hilo compara **extremos** `[s,e]`.  
   - **Salida**:
     ```
     [P-BSA] Paralelo  (ordenado)  : <tiempo> s  | found=1  | speedup=<t_BSA/t_PBSA>
     ```

3) **Caso NO ordenado** (mostrar límites de BSA/P-BSA)  
   - **BSA** en no ordenado → **no válido** (solo ilustrativo).  
   - **P-BSA** en no ordenado → **no aplica** (se reporta explícitamente).  
   - **Lineal** → baseline correcto para no ordenados.  
   - **Salida**:
     ```
     [BSA] Secuencial (NO ordenado): <tiempo> s  | (resultado NO válido; se muestra solo a modo ilustrativo)
     [P-BSA] Paralelo  (NO orden.) : <tiempo> s  | status=no-aplica
     [Lineal] (NO ordenado, correcto): <tiempo> s  | found=<0/1>
     ```

---

## 4) ¿Qué valor se busca en las búsquedas binarias? (y cómo cambiarlo)

Para que las pruebas sean **reproducibles y válidas**, el programa **elige automáticamente** un valor que **sí está** en el vector ordenado, y también un valor del vector no ordenado para ilustrar los casos.

### 4.1 En el vector **ordenado** (BSA y P-BSA)

En el código (dentro de `main`):

```cpp
auto Aord = random_vector(Nvec, 0, 2'000'000, g, /*sorted*/true);
int target = Aord[rndi(0, (int)Aord.size()-1, g)]; // elige un elemento que SÍ está
```

- `Aord` es un vector **ordenado** de tamaño `Nvec`.  
- `target` se toma **aleatoriamente entre los elementos existentes** de `Aord`, garantizando `found=1` si el algoritmo es correcto.

**Si quieres forzar “no encontrado” en ordenado**, después de esas líneas puedes **sobrescribir** `target` con un valor fuera del rango generado:
```cpp
int target = Aord[rndi(0, (int)Aord.size()-1, g)];
// Fuerza "miss":
target = 2'000'001; // mayor que el máximo 0..2'000'000
// o: target = -1; // menor que el mínimo
```

### 4.2 En el vector **no ordenado**

```cpp
auto Aunord = random_vector(Nvec, 0, 2'000'000, g, /*sorted*/false);
int target2 = Aunord[rndi(0, (int)Aunord.size()-1, g)];
```

- `Aunord` es **no ordenado**; `target2` es un elemento que **está** en él.  
- El programa muestra que **BSA** y **P-BSA** **no aplican** en no ordenados (por diseño de esos algoritmos), y por eso también imprime la **búsqueda lineal** como baseline correcto.

**Si quieres forzar “no encontrado” en no ordenado**, puedes sobrescribir:
```cpp
int target2 = Aunord[rndi(0, (int)Aunord.size()-1, g)];
// Fuerza "miss":
target2 = 2'000'001; // fuera del rango
// o: target2 = -1;
```

> **Resumen**  
> - En **ordenado**, el programa elige un valor existente para validar BSA/P-BSA.  
> - En **no ordenado**, también elige un valor existente, pero se usa para demostrar que BSA/P-BSA **no aplican** y que la **lineal** sí.  
> - Puedes **forzar misses** asignando manualmente un valor fuera del rango generado.

---

## 5) Ejemplo de salida comentado

```
Config: n=512  Nvec=1000000  threads=8  seed=42

[ShearSort] Secuencial    : 1.234 s  | snake-ok=true
[ShearSort] Paralelo      : 0.410 s  | snake-ok=true  | speedup=3.01
[ShearSort] Con transpose : 0.355 s  | snake-ok=true  | speedup(vs seq)=3.48

[BSA] Secuencial (ordenado)   : 0.0021 s  | found=1
[P-BSA] Paralelo  (ordenado)  : 0.0008 s  | found=1  | speedup=2.62
[BSA] Secuencial (NO ordenado): 0.0018 s  | (resultado NO válido; se muestra solo a modo ilustrativo)
[P-BSA] Paralelo  (NO orden.) : 0.0006 s  | status=no-aplica
[Lineal] (NO ordenado, correcto): 0.094 s  | found=1
```

- **ShearSort**: `snake-ok=true` confirma correctitud; compara tiempos y `speedup`.  
- **Búsquedas**: en **ordenado**, compara BSA vs P-BSA (mira `speedup`); en **no ordenado**, solo la **lineal** es válida.

---

## 6) Experimentos sugeridos (opcional)

**Speedup por hilos:**
```bash
./app 1024 2000000 1
./app 1024 2000000 2
./app 1024 2000000 4
./app 1024 2000000 8
./app 1024 2000000 16
```

**Escalado por tamaño:**
- ShearSort: `n = 256, 512, 1024, 2048`
- Búsqueda: `Nvec = 1e6, 2e6, 4e6, 8e6`

**Comparar “Paralelo” vs “Con transpose” (ShearSort)**:
- Para `n` grandes, la transposición suele ser más **cache-friendly** y rendir mejor.





