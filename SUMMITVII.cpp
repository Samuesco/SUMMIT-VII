#include <bits/stdc++.h>
#include <omp.h>
using namespace std;

// ---------- utilidades de tiempo ----------
struct Timer {
    double t0;
    void start() { t0 = omp_get_wtime(); }
    double stop() const { return omp_get_wtime() - t0; }
};

// ---------- helpers RNG ----------
static inline int rndi(int lo, int hi, std::mt19937 &g){ // [lo, hi]
    std::uniform_int_distribution<int> d(lo, hi);
    return d(g);
}

// ---------- chequeos ----------
static bool is_snake_sorted(const vector<vector<int>>& M){
    int n = (int)M.size();
    if(n==0) return true;
    int last = INT_MIN;
    for(int i=0;i<n;++i){
        if(i%2==0){
            // fila par: ascendente
            for(int j=0;j<n;++j){
                if(j>0 && M[i][j-1] > M[i][j]) return false;
                if(i==0 && j==0){ last = M[i][j]; }
                else { if(M[i][j] < last) return false; last = M[i][j]; }
            }
        } else {
            // fila impar: descendente
            for(int j=0;j<n;++j){
                int jj = n-1-j; // recorre como "serpiente" izquierda->derecha global
                if(j>0 && M[i][n-j] > M[i][n-j-1]) return false; // dentro de la fila: desc
                if(M[i][jj] < last) return false;
                last = M[i][jj];
            }
        }
    }
    return true;
}

// ---------- sort helpers ----------
static void sort_row(vector<vector<int>>& M, int row, bool ascending){
    auto &r = M[row];
    if(ascending) std::sort(r.begin(), r.end());
    else          std::sort(r.begin(), r.end(), std::greater<int>());
}

static void sort_col(vector<vector<int>>& M, int col){
    int n = (int)M.size();
    vector<int> c(n);
    for(int i=0;i<n;++i) c[i]=M[i][col];
    std::sort(c.begin(), c.end());
    for(int i=0;i<n;++i) M[i][col]=c[i];
}

static void transpose_square(vector<vector<int>>& M){
    int n=(int)M.size();
    for(int i=0;i<n;++i)
        for(int j=i+1;j<n;++j)
            std::swap(M[i][j], M[j][i]);
}

// ---------- ShearSort ----------
static void shearsort_seq(vector<vector<int>>& M){
    int n=(int)M.size();
    int rounds = (int)floor(log2(n)) + 1; // |log2(n)| + 1
    for(int r=0; r<rounds; ++r){
        // 1) ordenar filas en forma alternada
        for(int i=0;i<n;++i){
            bool asc = (i%2==0);
            sort_row(M, i, asc);
        }
        // 2) ordenar columnas ascendente
        for(int j=0;j<n;++j){
            sort_col(M, j);
        }
    }
}

// ---------- ShearSort (paralelo con OpenMP) ----------
static void shearsort_omp(vector<vector<int>>& M, int threads){
    int n=(int)M.size();
    int rounds = (int)floor(log2(n)) + 1;

    for(int r=0; r<rounds; ++r){
        // 1) filas (cada hilo toma varias filas independientes)
        #pragma omp parallel for num_threads(threads) schedule(static)
        for(int i=0;i<n;++i){
            bool asc = (i%2==0);
            sort_row(M, i, asc);
        }
        // 2) columnas (independientes entre sí)
        #pragma omp parallel for num_threads(threads) schedule(static)
        for(int j=0;j<n;++j){
            sort_col(M, j);
        }
    }
}

// ---------- ShearSort alternativo (solo sort-row + transposición) ----------
// Idea (pág. 2 del PDF): ordenar filas alternadas, transponer, ordenar filas, transponer…
static void shearsort_transpose_omp(vector<vector<int>>& M, int threads){
    int n=(int)M.size();
    int rounds = (int)floor(log2(n)) + 1;

    for(int r=0; r<rounds; ++r){
        // fase 1: ordenar filas en snake
        #pragma omp parallel for num_threads(threads) schedule(static)
        for(int i=0;i<n;++i){
            bool asc = (i%2==0);
            sort_row(M, i, asc);
        }
        // Transponer para que columnas se conviertan en filas
        transpose_square(M);
        // fase 2: ordenar filas todas ascendentes porque eran columnas
        #pragma omp parallel for num_threads(threads) schedule(static)
        for(int i=0;i<n;++i){
            sort_row(M, i, true);
        }
        // Deshacer la transposición
        transpose_square(M);
    }
}

// ---------- Búsqueda Binaria ----------
static int binary_search_seq(const vector<int>& a, int x){
    int L=0, R=(int)a.size()-1;
    while(L<=R){
        int mid = L + ((R-L)>>1);
        if(a[mid]==x) return mid;
        else if(a[mid]<x) L=mid+1;
        else R=mid-1;
    }
    return -1;
}

// ---------- Verificador de orden ----------
static bool is_sorted_asc(const vector<int>& a){
    for(size_t i=1;i<a.size();++i)
        if(a[i-1]>a[i]) return false;
    return true;
}

// ---------- Búsqueda Lineal  ----------
static int linear_search_seq(const vector<int>& a, int x){
    for(size_t i=0;i<a.size();++i) if(a[i]==x) return (int)i;
    return -1;
}

// ---------- P-BSA bUsqueda binaria paralela por subsecuencias ----------
// Si x cae dentro [inicio, fin], esa subsecuencia se convierte en el nuevo rango.
static int pbsa_parallel(const vector<int>& a, int x, int threads){
    if(a.empty()) return -1;
    if(!is_sorted_asc(a)){
        // La P-BSA requiere orden total (como BSA). En no ordenados, no aplica.
        // Para ser explícitos, devolvemos -2 como “no aplica”.
        return -2;
    }
    int L=0, R=(int)a.size()-1;
    int found = -1;

    while(L<=R && found==-1){
        int P = threads;
        int len = R-L+1;
        int chunk = max(1, len / P);

        // Variables compartidas entre hilos para decidir el nuevo rango
        int newL = -1, newR = -1;
        bool anyOwner = false;

        #pragma omp parallel for num_threads(threads) schedule(static)
        for(int p=0; p<P; ++p){
            int s = L + p*chunk;
            int e = (p==P-1) ? R : min(R, s + chunk - 1);
            if(s>R || e<L) continue;

            // Comparo extremos
            if(a[s]==x){ #pragma omp critical { found = s; } }
            else if(a[e]==x){ #pragma omp critical { found = e; } }
            else if(a[s] < x && x < a[e]){
                // Este bloque "posee" a x -> será el siguiente rango
                #pragma omp critical
                {
                    if(!anyOwner){
                        newL = s; newR = e;
                        anyOwner = true;
                    } else {
                        // unimos por seguridad si colindan solapes
                        newL = min(newL, s);
                        newR = max(newR, e);
                    }
                }
            }
        }
        if(found!=-1) return found;
        if(!anyOwner) return -1; // x no está en ningún subrango

        // Reducimos el rango y repetimos
        L = newL; R = newR;
    }
    return found;
}

// ---------- GeneraciOn de datos ----------
static vector<vector<int>> random_matrix(int n, int lo, int hi, std::mt19937& g){
    vector<vector<int>> M(n, vector<int>(n));
    for(int i=0;i<n;++i) for(int j=0;j<n;++j) M[i][j]=rndi(lo,hi,g);
    return M;
}

static vector<int> random_vector(int N, int lo, int hi, std::mt19937& g, bool sorted){
    vector<int> a(N);
    for(int i=0;i<N;++i) a[i]=rndi(lo,hi,g);
    if(sorted) std::sort(a.begin(), a.end());
    return a;
}

// ---------- metricas ----------
int main(int argc, char** argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n      = (argc>1)? atoi(argv[1]) : 512;        // matriz n x n
    int Nvec   = (argc>2)? atoi(argv[2]) : 1000000;    // tamaño vector para busqueda
    int th     = (argc>3)? atoi(argv[3]) : omp_get_max_threads();
    int seed   = (argc>4)? atoi(argv[4]) : 42;

    std::mt19937 g(seed);
    cout << "Config: n="<<n<<"  Nvec="<<Nvec<<"  threads="<<th<<"  seed="<<seed<<"\n\n";

    // ===== ShearSort =====
    auto M0 = random_matrix(n, 0, 1'000'000, g);

    // Secuencial
    auto M1 = M0;
    Timer t; t.start();
    shearsort_seq(M1);
    double t_seq = t.stop();
    cout << "[ShearSort] Secuencial    : " << t_seq << " s  | snake-ok=" << (is_snake_sorted(M1)?"true":"false") << "\n";

    // Paralelo (filas+columnas)
    auto M2 = M0;
    t.start();
    shearsort_omp(M2, th);
    double t_par = t.stop();
    cout << "[ShearSort] Paralelo      : " << t_par << " s  | snake-ok=" << (is_snake_sorted(M2)?"true":"false")
         << "  | speedup=" << (t_seq / t_par) << "\n";

    // Variante con transposición
    auto M3 = M0;
    t.start();
    shearsort_transpose_omp(M3, th);
    double t_tr = t.stop();
    cout << "[ShearSort] Con transpose : " << t_tr << " s  | snake-ok=" << (is_snake_sorted(M3)?"true":"false")
         << "  | speedup(vs seq)=" << (t_seq / t_tr) << "\n\n";

    // ===== Búsqueda =====
    // Vector ORDENADO
    auto Aord = random_vector(Nvec, 0, 2'000'000, g, /*sorted*/true);
    int target = Aord[rndi(0, (int)Aord.size()-1, g)]; // elige un elemento que está

    // BSA secuencial
    t.start();
    int pos_seq = binary_search_seq(Aord, target);
    double tb_seq = t.stop();
    cout << "[BSA] Secuencial (ordenado)   : " << tb_seq << " s  | found="<< (pos_seq!=-1) << "\n";

    // P-BSA paralelo
    t.start();
    int pos_pbsa = pbsa_parallel(Aord, target, th);
    double tb_par = t.stop();
    cout << "[P-BSA] Paralelo  (ordenado)  : " << tb_par << " s  | found="<< (pos_pbsa>=0)
         << "  | speedup=" << (tb_seq / tb_par) << "\n";

    // Vector NO ORDENADO (para mostrar que BSA/P-BSA no aplican)
    auto Aunord = random_vector(Nvec, 0, 2'000'000, g, /*sorted*/false);
    int target2 = Aunord[rndi(0, (int)Aunord.size()-1, g)];

    // BSA sobre no ordenado 
    t.start();
    int bogus = binary_search_seq(Aunord, target2);
    double t_bogus = t.stop();
    cout << "[BSA] Secuencial (NO ordenado): " << t_bogus << " s  | (resultado NO válido; se muestra solo a modo ilustrativo)\n";

    // P-BSA sobre no ordenado 
    t.start();
    int np = pbsa_parallel(Aunord, target2, th);
    double t_np = t.stop();
    cout << "[P-BSA] Paralelo  (NO orden.) : " << t_np << " s  | status="<< (np==-2? "no-aplica":"(indef)") << "\n";

    // Lineal como baseline correcto en no ordenados
    t.start();
    int pos_lin = linear_search_seq(Aunord, target2);
    double t_lin = t.stop();
    cout << "[Lineal] (NO ordenado, correcto): " << t_lin << " s  | found="<< (pos_lin!=-1) << "\n";

    cout << "\nListo. Copia/pega esta salida en tu reporte.\n";
    return 0;
}
