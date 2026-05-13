#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <array>
#include <math.h>
#include <iomanip>
#include <fstream>
using namespace std;
// on vas faire une subdivision uniforme de [-a, a] en N+1 points
vector<double> Subdiv(double a, int N){
    //on cree un vecteur de taille N+1 contenant les points de subdivision
    vector<double> xi(N+1);
    for(int i = 0 ; i<=N ; i++){
        // xi[i] = -a + i * pas, avec pas = 2a/N
        xi[i] = -a + i * (2*a) / N;
    }
    return xi;
}
// numero global d'un noeud (i,j) dans la grille
// on applique la formule : s = (N+1)*j + i
int numgb(int N, int M, int i, int j){
    return (N+1)*j + i;
}
// On inverse de numgb : retrouve (i,j) depuis le numero global S
// i = S mod (N+1)  et  j = S / (N+1)
vector<int> invnumgb(int N, int M, int S){
    int i = S % (N+1); // reste de la division euclidienne = indice en x
    int j = S / (N+1); // quotient de la division euclidienne = indice en y
    return {i, j};
}
//On fait la construction du maillage triangulaire
// On vas retourne une matrice TRG de taille K x 3, K = 2*N*M triangles
// Chaque ligne va contenir les 3 indices globaux des sommets
vector<vector<int>> maillageTR(int N, int M) {
    // K = nombre total de triangles
    int K = 2 * N * M;
    // TRG : table de connectivite
    vector<vector<int>> TRG(K);
    // l : indice du triangle courant
    int l = 0;
    // parcours de chaque cellule rectangulaire (i,j)
    for (int j = 0; j < M; j++) {
        for (int i = 0; i < N; i++) {
            // indices globaux des 4 coins de la cellule
            int s1 = numgb(N,M,i,  j  ); // bas-gauche
            int s2 = numgb(N,M,i+1,j  ); // bas-droite
            int s3 = numgb(N,M,i,  j+1); // haut-gauche
            int s4 = numgb(N,M,i+1,j+1); // haut-droite
            // decoupage en 2 triangles selon la parite de i+j
            if ((i + j) % 2 == 0) {
                // diagonale de s2 vers s3 (bas-droite -> haut-gauche)
                TRG[l++] = {s1, s2, s4}; // triangle T-
                TRG[l++] = {s1, s3, s4}; // triangle T+
            } else {
                // diagonale de s1 vers s4 (bas-gauche -> haut-droite)
                TRG[l++] = {s1, s2, s3}; // triangle T-
                TRG[l++] = {s2, s3, s4}; // triangle T+
            }
        }
    }
    return TRG;
}
// Calcule la matrice jacobienne B_T du triangle T
vector<vector<double>> CalcMatBT(vector<double> xs, vector<double> ys) {
    // BT est une matrice 2x2
    vector<vector<double>> BT(2, vector<double>(2));
    // premiere ligne : differences en x
    BT[0] = {xs[1] - xs[0], xs[2] - xs[0]};
    // deuxieme ligne : differences en y
    BT[1] = {ys[1] - ys[0], ys[2] - ys[0]};
    return BT;
}
// Coefficient de diffusion eta(x,y) = 1 - (eta0/4)*(x^2+y^2)
// ici eta0 = 0.5 fixe pour le test
double eta_coef(double x, double y){
    return 1;
}
// Calcule ET[l] = integrale de eta sur le triangle T_l
// par la formule de quadrature des points milieux
// ET[l] = |det(B_T)| * (1/6) * sum eta(points milieux)
vector<double> integ_eta_triang(
    double (*eta_coef)(double, double),
    const vector<vector<int>>& TRG,
    const vector<double>& xs,
    const vector<double>& ys)
{
    // K = nombre de triangles
    int K = TRG.size();
    // ET[l] contiendra l'integrale de eta sur le triangle l
    vector<double> ET(K);

    for (int l = 0; l < K; l++) {
        // On recupere les indices globaux des 3 sommets
        int n0 = TRG[l][0], n1 = TRG[l][1], n2 = TRG[l][2];
        // coordonnees des 3 sommets
        double x0 = xs[n0], y0 = ys[n0];
        double x1 = xs[n1], y1 = ys[n1];
        double x2 = xs[n2], y2 = ys[n2];
        // determinant de B_T = 2 * aire du triangle
        double detBT = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0);
        // 1er point de quadrature : milieu du cote [A0,A1]
        double px1 = x0 + 0.5*(x1-x0);
        double py1 = y0 + 0.5*(y1-y0);
        // 2eme point de quadrature : milieu du cote [A0,A2]
        double px2 = x0 + 0.5*(x2-x0);
        double py2 = y0 + 0.5*(y2-y0);
        // 3eme point de quadrature : milieu du cote [A1,A2]
        double px3 = x0 + 0.5*(x1-x0) + 0.5*(x2-x0);
        double py3 = y0 + 0.5*(y1-y0) + 0.5*(y2-y0);
        // formule des points milieux : poids 1/6 pour chaque point
        double quadrature = (1.0/6.0) * (eta_coef(px1,py1)
                                        + eta_coef(px2,py2)
                                        + eta_coef(px3,py3));
      // changement de variable : multiplication par |det(B_T)|
        ET[l] = abs(detBT) * quadrature;
    }
    return ET;
}
// Calcule la matrice de reaction locale 3x3
// R[i][r] = |det(B_T)| * int_T_hat lambda_hat_i * lambda_hat_r
// = |det(B_T)| * { 1/12 si i==r, 1/24 sinon }
vector<vector<double>> ReacTerm(vector<double> xs, vector<double> ys) {
    // determinant de B_T
    double det = (xs[1]-xs[0])*(ys[2]-ys[0])
               - (xs[2]-xs[0])*(ys[1]-ys[0]);
    // valeur absolue du determinant
    double absdet = abs(det);
    // matrice de reaction 3x3
    vector<vector<double>> R(3);
    // diagonale : absdet/12, hors-diagonale : absdet/24
    R[0] = {absdet/12.0, absdet/24.0, absdet/24.0};
    R[1] = {absdet/24.0, absdet/12.0, absdet/24.0};
    R[2] = {absdet/24.0, absdet/24.0, absdet/12.0};
    return R;
}
// Calcule la matrice de diffusion locale 3x3
// D[i][r] = val * (grad_lambda_i . grad_lambda_r)
// ou val = ET[t] = integrale de eta sur T
vector<vector<double>> DiffTerm(vector<double> xs, vector<double> ys, double val) {
    // coefficients de B_T
    double a = xs[1]-xs[0], b = xs[2]-xs[0];
    double c = ys[1]-ys[0], d = ys[2]-ys[0];
    // determinant de B_T
    double det = a*d - b*c;
    // inverse de B_T : B_T^{-1} = (1/det) * | d  -b |
    //                                       |-c   a |
    double inv[2][2] = {{ d/det, -b/det},
                        {-c/det,  a/det}};
    // gradients des fonctions de base sur T_hat
    // grad_hat_lambda_0 = (-1,-1), grad_hat_lambda_1 = (1,0), grad_hat_lambda_2 = (0,1)
    double gh[3][2] = {{-1,-1},{1,0},{0,1}};
    // gradients sur T via g[i] = B_T^{-T} * grad_hat[i]
    double g[3][2];
    for (int i = 0; i < 3; i++) {
        // composante x du gradient transforme
        g[i][0] = inv[0][0]*gh[i][0] + inv[1][0]*gh[i][1];
        // composante y du gradient transforme
        g[i][1] = inv[0][1]*gh[i][0] + inv[1][1]*gh[i][1];
    }
    // matrice de diffusion 3x3
    vector<vector<double>> D(3, vector<double>(3));
    // D[i][j] = val * produit scalaire des gradients
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            D[i][j] = val * (g[i][0]*g[j][0] + g[i][1]*g[j][1]);
            return D;
}
// Produit matrice-vecteur global A*V
// W[k] = sum_T sum_r V[j] * (diff[i][r] + lambda*reac[i][r])
vector<double> matvec(vector<double> V, int N, int M, double a, double b, double lambda, double (*eta_coef)(double, double))
{
    // table de connectivite
    vector<vector<int>> TRG = maillageTR(N, M);
    // nombre de noeuds et de triangles
    int G = (N+1)*(M+1);
    int K = 2*M*N;
     // subdivisions en x et y
    vector<double> sub_a = Subdiv(a, N);
    vector<double> sub_b = Subdiv(b, M);
    // coordonnees de tous les noeuds
    vector<double> xs(G), ys(G);
    for (int j = 0; j <= M; j++)
        for (int i = 0; i <= N; i++) {
            int s = numgb(N, M, i, j);
            xs[s] = sub_a[i];
            ys[s] = sub_b[j];
        }
        // integrales de eta sur chaque triangle
    vector<double> ET = integ_eta_triang(eta_coef, TRG, xs, ys);
    // vecteur resultat initialise a zero
    vector<double> W(G, 0.0);
    // matrices locales et coordonnees
    vector<vector<double>> coord(2, vector<double>(3));
    vector<vector<double>> diff(3, vector<double>(3));
    vector<vector<double>> reac(3, vector<double>(3));
    int k, j;
    double res, PROD2;
    // boucle sur les triangles
    for (int t = 0; t < K; t++) {
        // indices globaux des 3 sommets
        int n0 = TRG[t][0], n1 = TRG[t][1], n2 = TRG[t][2];
        // coordonnees des sommets
        coord[0] = {xs[n0], xs[n1], xs[n2]};
        coord[1] = {ys[n0], ys[n1], ys[n2]};
        // matrices locales de diffusion et reaction
        diff = DiffTerm(coord[0], coord[1], ET[t]);
        reac = ReacTerm(coord[0], coord[1]);
        // boucle sur les sommets locaux i vers noeud global k
        for (int i = 0; i < 3; i++) {
            // numero global du sommet i
            k = TRG[t][i];
            res = 0.0;
            // boucle sur les sommets locaux r vers noeud global j
            for (int r = 0; r < 3; r++) {
                // numero global du sommet r
                j = TRG[t][r];
                // PROD2 = a^T(w_j, w_k) = diff + lambda*reac
                PROD2 = diff[i][r] + lambda * reac[i][r];
                res = res + V[j] * PROD2;
            } // fin boucle r
            // accumulation W[k] += res
            W[k] = W[k] + res;
        } // fin boucle i
    } // fin boucle t
    return W;
}
// Second membre f(x,y) = cos(m*pi*x)*cos(m*pi*y)
double rhsf_global = 1; // parametre m global
// Calcule le vecteur second membre B de taille G
// B[s] = sum_T |det(BT)| * sum_t w[t]*f(F_T(v[t]))*lambda_hat_j(v[t])
vector<double> scdmembre(double (*rhsf)(double, double), int N,int M, double a, double b)
{
    // table de connectivite
    vector<vector<int>> TRG = maillageTR(N, M);
    // G = nombre total de noeuds
    int G = (N+1)*(M+1);
    // K = nombre total de triangles
    int K = 2*N*M;
    //On initialise B de taille G avec que des zeros
    vector<double> B(G, 0.0);
    // subdivisions en x et y
    vector<double> sub_a = Subdiv(a, N);
    vector<double> sub_b = Subdiv(b, M);
    // coordonnees de tous les noeuds
    vector<double> xs(G), ys(G);
    for (int j = 0; j <= M; j++)
        for (int i = 0; i <= N; i++) {
            int s = numgb(N, M, i, j);
            xs[s] = sub_a[i]; // coordonnee x du noeud s
            ys[s] = sub_b[j]; // coordonnee y du noeud s
        }
    // points de quadrature sur T_hat (points milieux des cotes)
    vector<vector<double>> v = {{0.5,0.0},{0.0,0.5},{0.5,0.5}};
    // poids de quadrature (tous egaux a 1/6)
    vector<double> w = {1.0/6, 1.0/6, 1.0/6};
    // coordonnees locales du triangle courant
    vector<vector<double>> coord(2, vector<double>(3));
    // matrice B_T du triangle courant
    vector<vector<double>> BT;
    // indice global du sommet courant
    int s;
    // contribution du triangle au noeud, valeur de f, valeur de lambda_hat_j
    double contrib, f_val, lambda_j;
    // boucle sur tous les triangles
    for (int i = 0; i < K; i++) {
        // indices globaux des 3 sommets du triangle i
        int n0 = TRG[i][0], n1 = TRG[i][1], n2 = TRG[i][2];
        // coordonnees x des 3 sommets
        coord[0] = {xs[n0], xs[n1], xs[n2]};
        // coordonnees y des 3 sommets
        coord[1] = {ys[n0], ys[n1], ys[n2]};
        // calcul de B_T pour ce triangle
        BT = CalcMatBT(coord[0], coord[1]);
        // |det(B_T)| = 2 * aire du triangle
        double det = abs(BT[0][0]*BT[1][1] - BT[0][1]*BT[1][0]);
        // boucle sur les 3 sommets du triangle
        for (int j = 0; j < 3; j++) {
            // numero global du sommet j
            s = TRG[i][j];
            // initialisation de la contribution de ce triangle a B[s]
            contrib = 0.0;
            // boucle sur les 3 points de quadrature
            for (int t = 0; t < 3; t++) {
                // transformation du point de quadrature : F_T(v[t])
                double xq = coord[0][0] + BT[0][0]*v[t][0] + BT[0][1]*v[t][1];
                double yq = coord[1][0] + BT[1][0]*v[t][0] + BT[1][1]*v[t][1];
                // evaluation de f au point transforme
                f_val = rhsf(xq, yq);
                // evaluation de lambda_hat_j au point de quadrature v[t]
                // depend du sommet local j (pas du point de quadrature t)
                if (j == 0) lambda_j = 1.0 - v[t][0] - v[t][1]; // lambda_hat_0
                else if (j == 1) lambda_j = v[t][0];                   // lambda_hat_1
                else             lambda_j = v[t][1];                   // lambda_hat_2
                // accumulation : poids * f * lambda_hat_j
                contrib += w[t] * f_val * lambda_j;
            }
            // multiplication par |det(B_T)| et accumulation dans B[s]
            B[s] += det * contrib;
        }
    }
    return B;
}
// Norme L2 de v : sqrt(V^t * A * V) avec eta=0, lambda=1
// On utilise seulement ReacTerm (pas de DiffTerm)
double normL2(vector<double> V, int N, int M, double a, double b) {
    // table de connectivite
    vector<vector<int>> TRG = maillageTR(N, M);
    // nombre de noeuds et de triangles
    int G = (N+1)*(M+1);
    int K = 2*M*N;
    // subdivisions
    vector<double> sub_a = Subdiv(a, N);
    vector<double> sub_b = Subdiv(b, M);
    // coordonnees de tous les noeuds
    vector<double> xs(G), ys(G);
    for (int jj = 0; jj <= M; jj++)
        for (int ii = 0; ii <= N; ii++) {
            int s = numgb(N, M, ii, jj);
            xs[s] = sub_a[ii];
            ys[s] = sub_b[jj];
        }
    // vecteur resultat W = A*V initialise a zero
    vector<double> W(G, 0.0);
    // coordonnees locales et matrice de reaction
    vector<vector<double>> coord(2, vector<double>(3));
    vector<vector<double>> reac(3, vector<double>(3));
    int k, j;
    double res, PROD2;
    // boucle sur les triangles (analogue a matvec avec eta=0, lambda=1)
    for (int t = 0; t < K; t++) {
        // indices globaux des 3 sommets
        int n0 = TRG[t][0], n1 = TRG[t][1], n2 = TRG[t][2];
        // coordonnees des sommets
        coord[0] = {xs[n0], xs[n1], xs[n2]};
        coord[1] = {ys[n0], ys[n1], ys[n2]};
        // seulement ReacTerm : eta=0 donc DiffTerm supprime
        reac = ReacTerm(coord[0], coord[1]);
        // boucle sur les sommets locaux i
        for (int i = 0; i < 3; i++) {
            // numero global du sommet i
            k = TRG[t][i];
            res = 0.0;
            // boucle sur les sommets locaux r
            for (int r = 0; r < 3; r++) {
                // numero global du sommet r
                j = TRG[t][r];
                // PROD2 = reac[i][r] car lambda=1 et pas de diffusion
                PROD2 = reac[i][r];
                res = res + V[j] * PROD2;
            } // fin boucle r
            // accumulation dans W[k]
            W[k] = W[k] + res;
        } // fin boucle i
    } // fin boucle t
    // ||v||^2_L2 = V^t * W, on retourne la racine
    double norme2 = 0.0;
    for (int k = 0; k < G; k++)
        norme2 += V[k] * W[k];

    return sqrt(norme2);
}
// Norme H1 semi-normee de v : sqrt(V^t * A * V) avec eta=1, lambda=0
// On utilise seulement DiffTerm (pas de ReacTerm)
double normL2Grad(vector<double> V, int N, int M, double a, double b) {
    // table de connectivite
    vector<vector<int>> TRG = maillageTR(N, M);
    // nombre de noeuds et de triangles
    int G = (N+1)*(M+1);
    int K = 2*M*N;
    // subdivisions
    vector<double> sub_a = Subdiv(a, N);
    vector<double> sub_b = Subdiv(b, M);
    // coordonnees de tous les noeuds
    vector<double> xs(G), ys(G);
    for (int jj = 0; jj <= M; jj++)
        for (int ii = 0; ii <= N; ii++) {
            int s = numgb(N, M, ii, jj);
            xs[s] = sub_a[ii];
            ys[s] = sub_b[jj];
        }
    // ET[t] = |det(BT)|/2 car eta=1 donc int_T eta = aire(T)
    vector<double> ET(K);
    for (int t = 0; t < K; t++) {
        int n0=TRG[t][0], n1=TRG[t][1], n2=TRG[t][2];
        double x0=xs[n0],y0=ys[n0];
        double x1=xs[n1],y1=ys[n1];
        double x2=xs[n2],y2=ys[n2];
        // |det(BT)| = 2*aire => int_T 1 = |det|/2
        ET[t] = abs((x1-x0)*(y2-y0)-(x2-x0)*(y1-y0)) / 2.0;
    }
    // vecteur resultat W = A*V initialise a zero
    vector<double> W(G, 0.0);
    // coordonnees locales et matrice de diffusion
    vector<vector<double>> coord(2, vector<double>(3));
    vector<vector<double>> diff(3, vector<double>(3));
    int k, j;
    double res, PROD2;
    // boucle sur les triangles (analogue a matvec avec eta=1, lambda=0)
    for (int t = 0; t < K; t++) {
        // indices globaux des 3 sommets
        int n0 = TRG[t][0], n1 = TRG[t][1], n2 = TRG[t][2];
        // coordonnees des sommets
        coord[0] = {xs[n0], xs[n1], xs[n2]};
        coord[1] = {ys[n0], ys[n1], ys[n2]};
        // seulement DiffTerm avec ET[t] = aire(T) car eta=1
        diff = DiffTerm(coord[0], coord[1], ET[t]);
        // boucle sur les sommets locaux i
        for (int i = 0; i < 3; i++) {
            // numero global du sommet i
            k = TRG[t][i];
            res = 0.0;
            // boucle sur les sommets locaux r
            for (int r = 0; r < 3; r++) {
                // numero global du sommet r
                j = TRG[t][r];
                // PROD2 = diff[i][r] car lambda=0 (pas de ReacTerm)
                PROD2 = diff[i][r];
                res = res + V[j] * PROD2;
            } // fin boucle r
            // accumulation dans W[k]
            W[k] = W[k] + res;
        } // fin boucle i
    } // fin boucle t
    // ||grad v||^2_L2 = V^t * W, on retourne la racine
    double norme2 = 0.0;
    for (int k = 0; k < G; k++)
        norme2 += V[k] * W[k];

    return sqrt(norme2);
}
// Produit scalaire euclidien de deux vecteurs
double pdt_scal(vector<double> u, vector<double> v) {
    double res = 0.0; // initialise a zero
    for (int i = 0; i < u.size(); i++)
        res += u[i] * v[i]; // somme des produits terme a terme
    return res;
}
// Norme euclidienne d'un vecteur
double norme_vect(vector<double> u) {
    return sqrt(pdt_scal(u, u)); // racine du produit scalaire avec lui-meme
}
// methode du gradient biconjugue stabilise BICGSTAB pour resoudre A*X = B
vector<double> bicgstab(
    vector<double> B,
    int N, int M,
    double a, double b,
    double lambda,
    double (*eta_coef)(double, double))
{
    // tolerance de convergence
    double tol= 1e-6;
    // nombre max d'iterations
    int itermax = 1000;
    // taille du systeme
    int n = B.size();
    // X0 = vecteur nul (solution initiale)
    vector<double> X(n, 0.0);
    // R0 = B - A*X0 = B  (car X0 = 0)
    vector<double> R = B;
    // R0* = R0 (choix arbitraire)
    vector<double> Rx = R;
    // W0 = R0
    vector<double> W = R;
    // vecteurs intermediaires de l'algorithme
    vector<double> AW(n), Sj(n), ASj(n), Rj(n);
    // scalaires de l'algorithme
    double alpha, omega, beta;
    // convergence initiale = norme du residu initial
    double conv = norme_vect(R);
    // compteur d'iterations
    int j = 0;
    // boucle jusqu'a convergence ou nombre max d'iterations
    while (j < itermax && conv > tol) {
        // AW_j = A * W_j
        AW = matvec(W, N, M, a, b, lambda, eta_coef);
        // alpha_j = <R_j, R0*> / <AW_j, R0*>
        alpha = pdt_scal(R, Rx) / pdt_scal(AW, Rx);
        // S_j = R_j - alpha_j * AW_j
        for (int i =0; i < n; i++)
            Sj[i] = R[i] - alpha * AW[i];
            // AS_j = A * S_j
        ASj = matvec(Sj, N, M, a, b, lambda, eta_coef);
        // omega_j = <AS_j, S_j> / <AS_j, AS_j>
        omega = pdt_scal(ASj, Sj) / pdt_scal(ASj, ASj);
        // X_{j+1} = X_j + alpha_j * W_j + omega_j * S_j
        for (int i = 0; i < n; i++)
            X[i] = X[i] + alpha * W[i] + omega * Sj[i];
        // R_{j+1} = S_j - omega_j * AS_j
        for (int i = 0; i < n; i++)
            Rj[i] = Sj[i] - omega * ASj[i];
        // beta_j = (<R_{j+1},R0*>/<R_j,R0*>) * (alpha_j/omega_j)
        beta = (pdt_scal(Rj, Rx) / pdt_scal(R, Rx)) * (alpha / omega);
        // W_{j+1} = R_{j+1} + beta_j * (W_j - omega_j * AW_j)
        for (int i = 0; i < n; i++)
            W[i] = Rj[i] + beta * (W[i] - omega * AW[i]);
        // mise a jour du residu et du compteur
        R = Rj;
        conv = norme_vect(R);
        j++;
    }
    return X;
}
// Constante pi
const double pi = 3.141592653589793;
// Parametre m global (entier naturel de l'enonce)
int m_global = 1;
// Solution exacte u^p(x,y) = cos(m*pi*x)*cos(m*pi*y)/(2m^2*pi^2+1)
double solExa(double x, double y) {
    return cos(m_global * pi * x) * cos(m_global * pi * y)
           / (2.0 * m_global*m_global * pi*pi + 1.0);
}
// Second membre f(x,y) = cos(m*pi*x)*cos(m*pi*y)
double rhsf(double x, double y) {
    return cos(m_global * pi * x) * cos(m_global * pi * y);
}
// Calcule les 3 erreurs relatives entre u^h et u^p
// e0   = ||uh_hat - uh||_L2    / ||uh_hat||_L2
// e1   = ||grad(uh_hat-uh)||_L2 / ||grad(uh_hat)||_L2
// einf = max|uh_hat - uh|       / max|uh_hat|
vector<double> erreurs(
    double (*solExa)(double, double),
    vector<double> uh,   // solution approchee de taille G
    int N, int M,
    double a, double b)
{
    // nombre total de noeuds
    int G = (N+1)*(M+1);
    // subdivisions
    vector<double> sub_a = Subdiv(a, N);
    vector<double> sub_b = Subdiv(b, M);
    // coordonnees de tous les noeuds
    vector<double> xs(G), ys(G);
    for (int j = 0; j <= M; j++)
        for (int i = 0; i <= N; i++) {
            int s = numgb(N, M, i, j);
            xs[s] = sub_a[i];
            ys[s] = sub_b[j];
        }
    // uh_hat = interpolee de u^p aux noeuds du maillage
    vector<double> uh_hat(G);
    for (int k = 0; k < G; k++)
        uh_hat[k] = solExa(xs[k], ys[k]);
    // diff_vec = uh_hat - uh (difference noeud par noeud)
    vector<double> diff_vec(G);
    for (int k = 0; k < G; k++)
        diff_vec[k] = uh_hat[k] - uh[k];
    // e0 = erreur relative en norme L2
    double e0 = normL2(diff_vec, N, M, a, b)
              / normL2(uh_hat,   N, M, a, b);
    // e1 = erreur relative en semi-norme H1 (gradient)
    double e1 = normL2Grad(diff_vec, N, M, a, b)
              / normL2Grad(uh_hat,   N, M, a, b);
    // calcul de l'erreur infinie
    double maxdiff = 0.0; // max de |uh_hat - uh|
    double maxw    = 0.0; // max de |uh_hat|

    for (int k = 0; k < G; k++) {
        // mise a jour du max de l'erreur absolue
        if (abs(diff_vec[k]) > maxdiff)
            maxdiff = abs(diff_vec[k]);
        // mise a jour du max de la solution exacte
        if (abs(uh_hat[k]) > maxw)
            maxw = abs(uh_hat[k]);
    }
    // einf = erreur relative en norme infinie
    double einf = maxdiff / maxw;
    // retourne les 3 erreurs
    return {e0, e1, einf};
}
int main() {

    double a      = 1.0;
    double b      = 1.0;
    double lambda = 1.0;
    m_global      = 1;    // m fixe a 1

    // En-tete du tableau
    cout << setw(5)  << "N"
         << setw(12) << "h"
         << setw(12) << "e0"
         << setw(12) << "e1"
         << setw(12) << "einf" << "\n";
    cout << string(53, '-') << "\n";

    // Boucle sur N croissant (h decroissant)
    for (int N : {5, 10, 20, 40, 80}) {

        int M = N;
        double h = 2.0 * a / N;  // pas du maillage

        // Second membre B
        vector<double> B = scdmembre(rhsf, N, M, a, b);

        // Resolution par BICGSTAB
        // eta_coef doit renvoyer 1.0 car eta0 = 0
        vector<double> uh = bicgstab(B, N, M, a, b,
                                     lambda, eta_coef);

        // Calcul des erreurs
        vector<double> err = erreurs(solExa, uh, N, M, a, b);

        cout << setw(5)  << N
             << setw(12) << scientific
             << setprecision(3) << h
             << setw(12) << err[0]
             << setw(12) << err[1]
             << setw(12) << err[2] << "\n";
    }

    return 0;
}