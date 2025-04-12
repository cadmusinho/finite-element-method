#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <cmath>

using namespace std;

#define npc 4

double dn1ksi(double a);
double dn2ksi(double a);
double dn3ksi(double a);
double dn4ksi(double a);
double dn1eta(double a);
double dn2eta(double a);
double dn3eta(double a);
double dn4eta(double a);

double n1(double ksi, double eta);
double n2(double ksi, double eta);
double n3(double ksi, double eta);
double n4(double ksi, double eta);

struct GaussIntegration {
    double* w = new double[sqrt(npc)];
    double* pc = new double[sqrt(npc) + 2];

    void initialize() {
        if (npc == 4) {
            w[0] = 1.0;
            w[1] = 1.0;
            pc[0] = -1.0;
            pc[1] = 1.0;
            pc[2] = -(1.0 / sqrt(3));
            pc[3] = 1.0 / sqrt(3);
        }
        if (npc == 9) {
            w[0] = 5.0 / 9.0;
            w[1] = 8.0 / 9.0;
            w[2] = 5.0 / 9.0;
            pc[0] = -1.0;
            pc[1] = 1.0;
            pc[2] = -sqrt((3.0 / 5.0));
            pc[3] = 0.0;
            pc[4] = sqrt((3.0 / 5.0));
        }
        if (npc == 16) {
            w[0] = (18.0 - sqrt(30.0)) / 36.0;
            w[1] = (18.0 + sqrt(30.0)) / 36.0;
            w[2] = (18.0 + sqrt(30.0)) / 36.0;
            w[3] = (18.0 - sqrt(30.0)) / 36.0;
            pc[0] = -1.0;
            pc[1] = 1.0;
            pc[2] = -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));
            pc[3] = -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));
            pc[4] = sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));
            pc[5] = sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));
        }
    }
    ~GaussIntegration() {
        delete[] w;
        delete[] pc;
    }
};

struct ElemUniv {
    double** ksi;
    double** eta;

    struct Surface {
        double** N;

        void initialize() {
            N = new double* [sqrt(npc)];
            for (int i = 0; i < sqrt(npc); i++) {
                N[i] = new double[4] {0.0};
            }
        }

        void calculateN(int i, GaussIntegration* gi) {
            int z = 0;
            if (npc == 4) z = 2;
            if (npc == 9) z = 3;
            if (npc == 16) z = 4;

            if (i == 0) {
                for (int j = 0; j < z; ++j) {
                    N[j][0] = n1(gi->pc[2 + j], gi->pc[0]);
                    N[j][1] = n2(gi->pc[2 + j], gi->pc[0]);
                    N[j][2] = n3(gi->pc[2 + j], gi->pc[0]);
                    N[j][3] = n4(gi->pc[2 + j], gi->pc[0]);
                }
            }
            else if (i == 1) {
                for (int j = 0; j < z; ++j) {
                    N[j][0] = n1(gi->pc[1], gi->pc[2 + j]);
                    N[j][1] = n2(gi->pc[1], gi->pc[2 + j]);
                    N[j][2] = n3(gi->pc[1], gi->pc[2 + j]);
                    N[j][3] = n4(gi->pc[1], gi->pc[2 + j]);
                }
            }
            else if (i == 2) {
                for (int j = 0; j < z; ++j) {
                    N[j][0] = n1(gi->pc[2 + j], gi->pc[1]);
                    N[j][1] = n2(gi->pc[2 + j], gi->pc[1]);
                    N[j][2] = n3(gi->pc[2 + j], gi->pc[1]);
                    N[j][3] = n4(gi->pc[2 + j], gi->pc[1]);
                }
            }
            else if (i == 3) {
                for (int j = 0; j < z; ++j) {
                    N[j][0] = n1(gi->pc[0], gi->pc[2 + j]);
                    N[j][1] = n2(gi->pc[0], gi->pc[2 + j]);
                    N[j][2] = n3(gi->pc[0], gi->pc[2 + j]);
                    N[j][3] = n4(gi->pc[0], gi->pc[2 + j]);
                }
            }
        }

        void calculateHbc(double* nodex, double* nodey, bool* node2, int* pom, double** tab, int c, int z, double* p, int t, GaussIntegration* gi) {
            double** pom1 = new double* [4];
            double** pom2 = new double* [4];
            double** pom3 = new double* [4];
            double** pom4 = new double* [4];
            double** tab2 = new double* [4];
            double* p1 = new double[4] {0.0};
            double* p2 = new double[4] {0.0};
            double* p3 = new double[4] {0.0};
            double* p4 = new double[4] {0.0};
            for (int i = 0; i < 4; i++) {
                pom1[i] = new double[4] {0.0};
                pom2[i] = new double[4] {0.0};
                pom3[i] = new double[4] {0.0};
                pom4[i] = new double[4] {0.0};
                tab2[i] = new double[4] {0.0};
            }
            double dx = nodex[pom[z]] - nodex[pom[(z + 1) % 4]];
            double dy = nodey[pom[z]] - nodey[pom[(z + 1) % 4]];
            double det = sqrt(dx * dx + dy * dy) / 2.0;

            if (node2[pom[z]] == true && node2[pom[(z + 1) % 4]] == true) {
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        pom1[i][j] += N[0][i] * N[0][j] * gi->w[0] * (double)c * det;
                        pom2[i][j] += N[1][i] * N[1][j] * gi->w[1] * (double)c * det;
                        if (npc == 9) {
                            pom3[i][j] += N[2][i] * N[2][j] * gi->w[2] * (double)c * det;
                        }
                        if (npc == 16) {
                            pom3[i][j] += N[2][i] * N[2][j] * gi->w[2] * (double)c * det;
                            pom4[i][j] += N[3][i] * N[3][j] * gi->w[3] * (double)c * det;
                        }
                    }
                    p1[i] += (double)c * det * N[0][i] * (double)t * gi->w[0];
                    p2[i] += (double)c * det * N[1][i] * (double)t * gi->w[1];
                    if (npc == 9) {
                        p3[i] += (double)c * det * N[2][i] * (double)t * gi->w[2];
                    }
                    if (npc == 16) {
                        p3[i] += (double)c * det * N[2][i] * (double)t * gi->w[2];
                        p4[i] += (double)c * det * N[3][i] * (double)t * gi->w[3];
                    }
                }
            }
            else {}

            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    tab2[i][j] += (pom1[i][j] + pom2[i][j] + pom3[i][j] + pom4[i][j]);
                }
            }

            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    tab[i][j] += tab2[i][j];
                }
                p[i] += p1[i] + p2[i] + p3[i];
            }

            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    tab2[i][j] = 0.0;
                }
            }

            for (int i = 0; i < 4; i++) {
                delete[] pom1[i];
                delete[] pom2[i];
                delete[] pom3[i];
                delete[] pom4[i];
                delete[] tab2[i];
            }
            delete[] pom1;
            delete[] pom2;
            delete[] pom3;
            delete[] pom4;
            delete[] tab2;
            delete[] p1;
            delete[] p2;
            delete[] p3;
            delete[] p4;
        }

    };
    Surface surface[npc];

    ElemUniv() : ksi(nullptr), eta(nullptr) {}

    ~ElemUniv() {
        for (int i = 0; i < npc; i++) {
            delete[] ksi[i];
            delete[] eta[i];
        }
        delete[] ksi;
        delete[] eta;
    }

    void initialize(GaussIntegration* gi) {
        for (int i = 0; i < 4; i++) surface[i].initialize();
        for (int i = 0; i < 4; i++) {
            surface[i].calculateN(i, gi);
        }
        ksi = new double* [npc];
        eta = new double* [npc];
        for (int i = 0; i < npc; i++) {
            ksi[i] = new double[4];
            eta[i] = new double[4];
        }
        double* punkt;
        int size;
        if (npc == 4) {
            static double punkt_4[] = { 1.0 / sqrt(3.0), -1.0 / sqrt(3.0) };
            punkt = punkt_4;
            size = 2;
        }
        else if (npc == 9) {
            static double punkt_9[] = { sqrt(3.0 / 5.0), 0.0, -sqrt(3.0 / 5.0) };
            punkt = punkt_9;
            size = 3;
        }
        else if (npc == 16) {
            static double punkt_16[] = {
                sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)),
                sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)),
                -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)),
                -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0))
            };
            punkt = punkt_16;
            size = 4;
        }
        int index = 0;
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                ksi[index][0] = dn1ksi(punkt[j]);
                ksi[index][1] = dn2ksi(punkt[j]);
                ksi[index][2] = dn3ksi(punkt[j]);
                ksi[index][3] = dn4ksi(punkt[j]);
                eta[index][0] = dn1eta(punkt[i]);
                eta[index][1] = dn2eta(punkt[i]);
                eta[index][2] = dn3eta(punkt[i]);
                eta[index][3] = dn4eta(punkt[i]);
                index++;
            }
        }
    }

    void calculateHbcForSide(double* nodex, double* nodey, bool* node2, int* tab, double** H_bc_temp, int c, double* p, int t, GaussIntegration* gi) {
        for (int z = 0; z < 4; z++) {
            surface[z].calculateHbc(nodex, nodey, node2, tab, H_bc_temp, c, z, p, t, gi);
        }
    }
};

struct Node {
    double x = 0.0;
    double y = 0.0;
    bool bc = false;
};

struct Jakobian {
    double j[2][2] = { 0.0 };
    double j1[2][2] = { 0.0 };
    double detJ;

    void calculateJ(ElemUniv* eu, Node* nodes, int* tab, int pointIndex) {

        for (int k = 0; k < 4; k++) {
            int nodeId = tab[k] - 1;
            this->j[0][0] += eu->ksi[pointIndex][k] * nodes[nodeId].x;
            this->j[0][1] += eu->ksi[pointIndex][k] * nodes[nodeId].y;
            this->j[1][0] += eu->eta[pointIndex][k] * nodes[nodeId].x;
            this->j[1][1] += eu->eta[pointIndex][k] * nodes[nodeId].y;
        }
    }

    void calculateDetJ() {
        this->detJ = (j[0][0] * j[1][1]) - (j[0][1] * j[1][0]);
    }

    void calculateJ1() {
        calculateDetJ();
        if (detJ != 0) {
            j1[0][0] = j[1][1] / detJ;
            j1[0][1] = -j[0][1] / detJ;
            j1[1][0] = -j[1][0] / detJ;
            j1[1][1] = j[0][0] / detJ;
        }
    }
};

struct Element {
    int* values = new int[4];
    Jakobian* jakobian = new Jakobian[npc];
    double** H = new double* [4];
    double** Hbc = new double* [4];
    double** C = new double* [4];
    double* p = new double[4] {0.0};

    void initialize() {
        for (int i = 0; i < 4; i++) {
            H[i] = new double[4] {0.0};
            Hbc[i] = new double[4] {0.0};
            C[i] = new double[4] {0.0};
        }
    }

    ~Element() {
        delete[] values;
        delete[] jakobian;
        for (int i = 0; i < 4; i++) {
            delete[] H[i];
            delete[] Hbc[i];
            delete[] C[i];
        }
        delete[] H;
        delete[] Hbc;
        delete[] C;
        delete[] p;
    }

    void jacobians(ElemUniv* eu, Node* nodes) {
        for (int j = 0; j < npc; j++) {
            jakobian[j].calculateJ(eu, nodes, values, j);
            jakobian[j].calculateJ1();
        }
    }

    void calculateH(ElemUniv* eu, int conductivity) {
        double** H_dx = new double* [npc];
        double** H_dy = new double* [npc];
        double** H_dx_temp = new double* [4];
        double** H_dy_temp = new double* [4];
        double** H_temp = new double* [4];
        double w[4];
        if (npc == 4) {
            w[0] = 1.0;
            w[1] = 1.0;
        }
        if (npc == 9) {
            w[0] = 5.0 / 9.0;
            w[1] = 8.0 / 9.0;
            w[2] = 5.0 / 9.0;
        }
        if (npc == 16) {
            w[0] = (18.0 - sqrt(30.0)) / 36.0;
            w[1] = (18.0 + sqrt(30.0)) / 36.0;
            w[2] = (18.0 + sqrt(30.0)) / 36.0;
            w[3] = (18.0 - sqrt(30.0)) / 36.0;
        }
        for (int i = 0; i < npc; i++) {
            H_dx[i] = new double[4] {0.0};
        }
        for (int i = 0; i < npc; i++) {
            H_dy[i] = new double[4] {0.0};
        }
        for (int i = 0; i < 4; i++) {
            H_dx_temp[i] = new double[4] {0.0};
        }
        for (int i = 0; i < 4; i++) {
            H_dy_temp[i] = new double[4] {0.0};
        }
        for (int i = 0; i < 4; i++) {
            H_temp[i] = new double[4] {0.0};
        }

        for (int i = 0; i < npc; i++) {
            for (int j = 0; j < 4; j++) {
                H_dx[i][j] = jakobian[i].j1[0][0] * eu->ksi[i][j] + jakobian[i].j1[0][1] * eu->eta[i][j];
                H_dy[i][j] = jakobian[i].j1[1][0] * eu->ksi[i][j] + jakobian[i].j1[1][1] * eu->eta[i][j];
            }
        }

        for (int z = 0; z < npc; z++) {
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    H_dx_temp[i][j] = H_dx[z][i] * H_dx[z][j];
                    H_dy_temp[i][j] = H_dy[z][i] * H_dy[z][j];
                }
            }
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    H_temp[i][j] = (H_dx_temp[i][j] + H_dy_temp[i][j]) * conductivity * jakobian[z].detJ;
                }
            }

            int numW = (npc == 4) ? 2 : (npc == 9) ? 3 : (npc == 16) ? 4 : 5;
            int maxZ = (npc == 4) ? 1 : (npc == 9) ? 9 : (npc == 16) ? 15 : 24;
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    int wIndex1 = z / numW;
                    int wIndex2 = z % numW;
                    H_temp[i][j] *= w[wIndex1] * w[wIndex2];
                }
            }

            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    H[i][j] += H_temp[i][j];
                }
            }
        }

        for (int i = 0; i < npc; i++) {
            delete[] H_dx[i];
            delete[] H_dy[i];
        }
        for (int i = 0; i < 4; i++) {
            delete[] H_dx_temp[i];
            delete[] H_dy_temp[i];
            delete[] H_temp[i];
        }
        delete[] H_dx;
        delete[] H_dy;
        delete[] H_dx_temp;
        delete[] H_dy_temp;
        delete[] H_temp;
    }

    void calculateHbcForEachElement(double* nodex, double* nodey, bool* node2, ElemUniv* euu, int c, int t, GaussIntegration* gi) {
        int* pom = new int[4];
        for (int i = 0; i < 4; i++) {
            pom[i] = values[i] - 1;
        }
        euu->calculateHbcForSide(nodex, nodey, node2, pom, Hbc, c, p, t, gi);
        delete[] pom;
    }

    void  calculateCForEachElement(double sh, double density, GaussIntegration* gi) {
        double** C_pom = new double* [npc];
        double** C_pc = new double* [4];

        for (int i = 0; i < npc; i++) {
            C_pom[i] = new double[6] {0.0};
        }
        for (int i = 0; i < 4; i++) {
            C_pc[i] = new double[4] {0.0};
        }

        int zc = 0;
        if (npc == 4) zc = 3;
        if (npc == 9) zc = 4;
        if (npc == 16) zc = 5;

        int index = 0;
        for (int x = 2; x <= zc; ++x) {
            for (int y = 2; y <= zc; ++y) {
                C_pom[index][0] = n1(gi->pc[x], gi->pc[y]);
                C_pom[index][1] = n2(gi->pc[x], gi->pc[y]);
                C_pom[index][2] = n3(gi->pc[x], gi->pc[y]);
                C_pom[index][3] = n4(gi->pc[x], gi->pc[y]);
                C_pom[index][4] = gi->w[x - 2];
                C_pom[index][5] = gi->w[y - 2];
                ++index;
            }
        }

        for (int z = 0; z < npc; z++) {
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    C_pc[i][j] += C_pom[z][i] * C_pom[z][j] * sh * density * jakobian[z].detJ * C_pom[z][4] * C_pom[z][5];
                }
            }
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    C[i][j] += C_pc[i][j];
                }
            }
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    C_pc[i][j] = 0.0;
                }
            }
        }

        for (int i = 0; i < npc; i++) {
            delete[] C_pom[i];
        }
        for (int i = 0; i < 4; i++) {
            delete[] C_pc[i];
        }

        delete[] C_pom;
        delete[] C_pc;
    }
};

struct Grid {
    int nN;
    int nE;
    Element* element;
    Node* node;

    Grid() : element(nullptr), node(nullptr) {}

    ~Grid() {
        delete[] element;
        delete[] node;
    }

    void initialize() {
        element = new Element[nE];
        node = new Node[nN];
    }

    void Jacobian(ElemUniv* eu) {
        for (int i = 0; i < nE; i++) {
            element[i].jacobians(eu, node);
        }
    }

    void NEH(ElemUniv* eu, int c) {
        for (int i = 0; i < nE; i++) {
            element[i].initialize();
            element[i].calculateH(eu, c);
        }
    }

    void Hbc(Node* node, ElemUniv* euu, int c, int t, GaussIntegration* gi) {
        bool* node2 = new bool[nN];
        double* nodex = new double[nN];
        double* nodey = new double[nN];

        for (int i = 0; i < nN; i++) {
            node2[i] = node[i].bc;
            nodex[i] = node[i].x;
            nodey[i] = node[i].y;
        }

        for (int i = 0; i < nE; i++) {
            element[i].calculateHbcForEachElement(nodex, nodey, node2, euu, c, t, gi);
        }

        delete[] node2;
        delete[] nodex;
        delete[] nodey;
    }

    void C(double sh, double density, GaussIntegration* gi) {
        for (int i = 0; i < nE; i++) {
            element[i].calculateCForEachElement(sh, density, gi);
        }
    }
};

struct GlobalData {
    int simulationTime;
    int simulationStepTime;
    int conductivity;
    int alfa;
    int tot;
    int initialTemp;
    int density;
    int specificHeat;
    int nN;
    int nE;
};

struct HGlobal {
    double** values;
    double** Cg;
    int size;
    double* Pg;
    double* t;

    HGlobal(GlobalData* d) {
        size = d->nN;
        values = new double* [size];
        Cg = new double* [size];
        for (int i = 0; i < size; i++) {
            values[i] = new double[size] {0};
            Cg[i] = new double[size] {0};
        }
        Pg = new double[size] {0};
        t = new double[size] {0};
    }

    ~HGlobal() {
        for (int i = 0; i < size; i++) {
            delete[] values[i];
            delete[] Cg[i];
        }
        delete[] values;
        delete[] Cg;
        delete[] Pg;
        delete[] t;
    }

    void simulation(int time, int step, Grid* grid, int startTemp) {
        aggregated(grid);
        double* tzero = new double[size];
        double** pom = new double* [size];
        for (int i = 0; i < size; i++) {
            tzero[i] = (double)startTemp;
            pom[i] = new double[size] {0.0};
        }
        cout << endl;
        for (int iteracja = step; iteracja <= time; iteracja += step) {
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    values[i][j] += Cg[i][j] / iteracja;
                    if (Cg[i][j] != 0.0) {
                        Pg[i] += (Cg[i][j] / iteracja) * tzero[j];
                    }
                }
            }
            calculateT(grid);
            findMinMax(t);
            for (int i = 0; i < size; i++) {
                tzero[i] = t[i];
                for (int j = 0; j < size; j++) {
                    values[i][j] = 0.0;
                }
                Pg[i] = 0.0;
            }
            aggregated(grid);
        }

        cout << "\n\n\n";

        for (int i = 0; i < size; i++) {
            delete[] pom[i];
        }
        delete[] pom;
        delete[] tzero;
    }

    void findMinMax(double* tab) {
        double min = 0, max = 0;
        min = tab[0], max = tab[0];
        for (int i = 1; i < size; i++) {
            if (tab[i] < min) min = tab[i];
            else if (tab[i] > max) max = tab[i];
            else {}
        }
        printf("\n%0.10lf   %0.10lf", min, max);
    }

    void calculateT(Grid* g) {
        double** tab = new double* [g->nN];
        for (int z = 0; z < g->nN; z++) {
            tab[z] = new double[g->nN + 1] {0.0};
        }

        for (int z = 0; z < g->nN; z++) {
            for (int i = 0; i < g->nN; i++) {
                tab[z][i] = values[z][i];
            }
        }
        for (int i = 0; i < g->nN; i++) {
            tab[i][g->nN] = Pg[i];
        }
        for (int z = 0; z < g->nN; z++) {
            for (int i = z + 1; i < g->nN; i++) {
                for (int j = g->nN; j >= 0; j--) {
                    tab[i][j] -= tab[z][j] * (tab[i][z] / tab[z][z]);
                }
            }
        }

        for (int i = g->nN - 1; i >= 0; i--) {
            t[i] = tab[i][g->nN];
            for (int j = i + 1; j < g->nN; j++) {
                t[i] -= tab[i][j] * t[j];
            }
            t[i] /= tab[i][i];
        }

        for (int i = 0; i < g->nN; i++) {
            delete[] tab[i];
        }
        delete[] tab;
    }

    void aggregated(Grid* g) {
        for (int i = 0; i < g->nE; i++) {
            aggregation(&(g->element[i]));
        }
    }

    void aggregation(Element* e) {
        for (int k = 0; k < 4; k++) {
            for (int j = 0; j < 4; j++) {
                values[e->values[k] - 1][e->values[j] - 1] += e->H[k][j] + e->Hbc[k][j];
                Cg[e->values[k] - 1][e->values[j] - 1] += e->C[k][j];
            }
            Pg[e->values[k] - 1] += e->p[k];
        }
    }
};

double dn1ksi(double a) {
    return -0.25 * (1 - a);
}

double dn2ksi(double a) {
    return 0.25 * (1 - a);
}

double dn3ksi(double a) {
    return 0.25 * (1 + a);
}

double dn4ksi(double a) {
    return -0.25 * (1 + a);
}

double dn1eta(double a) {
    return -0.25 * (1 - a);
}

double dn2eta(double a) {
    return -0.25 * (1 + a);
}

double dn3eta(double a) {
    return 0.25 * (1 + a);
}

double dn4eta(double a) {
    return 0.25 * (1 - a);
}

double n1(double ksi, double eta) {
    return 0.25 * (1 - ksi) * (1 - eta);
}
double n2(double ksi, double eta) {
    return 0.25 * (1 + ksi) * (1 - eta);
}
double n3(double ksi, double eta) {
    return 0.25 * (1 + ksi) * (1 + eta);
}
double n4(double ksi, double eta) {
    return 0.25 * (1 - ksi) * (1 + eta);
}

int main() {
    try {
        if (npc != 4 && npc != 9 && npc != 16) {
            throw std::runtime_error("B³¹d: npc ma nieprawid³ow¹ wartoœæ!");
        }
    }
    catch (const exception& e) {
        std::cerr << "Wyst¹pi³ wyj¹tek: " << e.what() << std::endl;
        return 1;
    }
    GlobalData data;
    ElemUniv eu;
    GaussIntegration gi;
    gi.initialize();
    eu.initialize(&gi);
    Grid grid;
    string line, word;
    ifstream myfile("Test3_31_31_kwadrat.txt");
    int* pom = new int[8];
    if (myfile.is_open()) {
        for (int i = 0; i <= 7; i++) {
            getline(myfile, line);
            size_t space = line.find(' ');
            word = line.substr(space + 1);
            pom[i] = stoi(word);
        }
        data.simulationTime = pom[0];
        data.simulationStepTime = pom[1];
        data.conductivity = pom[2];
        data.alfa = pom[3];
        data.tot = pom[4];
        data.initialTemp = pom[5];
        data.density = pom[6];
        data.specificHeat = pom[7];

        for (int i = 0; i <= 1; i++) {
            getline(myfile, line);
            size_t space = line.find(' ');
            size_t secondSpace = line.find(' ', space + 1);
            word = line.substr(secondSpace + 1);
            pom[i] = stoi(word);
        }
        data.nN = pom[0];
        data.nE = pom[1];
        grid.nN = pom[0];
        grid.nE = pom[1];
        grid.initialize();
        getline(myfile, line);

        for (int i = 0; i < grid.nN; i++) {
            getline(myfile, line);
            stringstream ss(line);
            int id;
            char comma;
            ss >> id >> comma >> grid.node[i].x >> comma >> grid.node[i].y;
        }
        getline(myfile, line);

        for (int i = 0; i < grid.nE; i++) {
            getline(myfile, line);
            stringstream ss(line);
            int id;
            char comma;
            ss >> id >> comma >> grid.element[i].values[0] >> comma >> grid.element[i].values[1]
                >> comma >> grid.element[i].values[2] >> comma >> grid.element[i].values[3];
        }
        getline(myfile, line);

        int* pom2 = new int[data.nN];
        for (int i = 0; i < data.nN; i++) {
            pom2[i] = 0;
        }
        int index = 0;
        getline(myfile, line);
        stringstream ss(line);
        string temp;
        while (getline(ss, temp, ',')) {
            temp.erase(0, temp.find_first_not_of(" \t"));
            temp.erase(temp.find_last_not_of(" \t") + 1);
            if (index < data.nN) {
                pom2[index++] = stoi(temp);
            }
        }
        for (int i = 0; i < data.nN; i++) {
            for (int j = 0; j < index; ++j) {
                if (i + 1 == pom2[j]) {
                    grid.node[i].bc = true;
                    break;
                }
            }
        }
        delete[] pom2;

        myfile.close();
    }
    else {
        cout << "Unable to open file" << endl;
    }

    grid.Jacobian(&eu);
    grid.NEH(&eu, data.conductivity);
    grid.Hbc(grid.node, &eu, data.alfa, data.tot, &gi);
    grid.C(data.specificHeat, data.density, &gi);
    HGlobal hg(&data);
    hg.simulation(data.simulationTime, data.simulationStepTime, &grid, data.initialTemp);

    delete[] pom;

    return 0;
}
