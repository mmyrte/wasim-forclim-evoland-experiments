#!/bin/sh
# Regenerate reference.cc from WaSiM's sources: copies SOILtabletype::LookUpTableGeneration, h_calc,
# k_calc and get_theta verbatim from src/defltuzr.cc into a standalone driver.
# usage: extract_reference.sh WASIM_SRC_DIR > reference.cc
set -e
SRC="$1/defltuzr.cc"
fn() { awk -v sig="$1" 'index($0, sig) == 1 {on = 1} on {print} on && /^}/ {exit}' "$SRC"; }
cat <<'HEAD'
// Reference driver for the GT4Py port of WaSiM's soil hydraulic functions (gt4py/soil/van_genuchten.py).
// The function bodies are copied verbatim from wasim-eth src/defltuzr.cc by extract_reference.sh.
// input (stdin): lines "thetas thetar alpha n tau" (one soil horizon each), a line "END", then
// query lines "layer kind x [SE thetarex thetasex]" with kind h (h_calc of theta), k (k_calc) or t (get_theta of h).
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>
using namespace std;
struct SOILtabletype {
	int anz_k = 30, anz_h = 30;
	double **theta_k, **k_tab, **theta_h, **h_tab, *thetas, *thetar, *alpha, *n, *m, *tau;
	void LookUpTableGeneration(int layer);
	double k_calc(int layer, double theta, double SE, double thetarex = -1, double thetasex = -1);
	double h_calc(int layer, double theta, double thetarex = -1, double thetasex = -1);
	double get_theta(int layer, double h, double thetarex = -1, double thetasex = -1);
};
HEAD
fn "void SOILtabletype::LookUpTableGeneration(int layer)"
fn "double SOILtabletype::h_calc("
fn "double SOILtabletype::k_calc("
fn "double SOILtabletype::get_theta("
cat <<'MAIN'
int main()
{
	std::vector<double> p[5];
	char line[512];
	while (fgets(line, sizeof line, stdin) && strncmp(line, "END", 3) != 0) {
		double v[5];
		if (sscanf(line, "%lf %lf %lf %lf %lf", &v[0], &v[1], &v[2], &v[3], &v[4]) == 5)
			for (int i = 0; i < 5; i++) p[i].push_back(v[i]);
	}
	int L = p[0].size();
	SOILtabletype s;
	s.theta_k = new double*[L + 1]; s.k_tab = new double*[L + 1]; s.theta_h = new double*[L + 1]; s.h_tab = new double*[L + 1];
	s.thetas = new double[L + 1]; s.thetar = new double[L + 1]; s.alpha = new double[L + 1]; s.n = new double[L + 1]; s.m = new double[L + 1]; s.tau = new double[L + 1];
	for (int l = 1; l <= L; l++) {
		s.thetas[l] = p[0][l - 1]; s.thetar[l] = p[1][l - 1]; s.alpha[l] = p[2][l - 1]; s.n[l] = p[3][l - 1]; s.tau[l] = p[4][l - 1];
		s.m[l] = 1.0 - 1.0 / s.n[l]; // as in SOILtabletype::initialize
		s.LookUpTableGeneration(l);
		s.h_tab[l][0] = 0; // index 0 is never initialized in WaSiM; get_theta reads it for h <= 0
	}
	int layer; char kind; double x, SE, tr, ts;
	while (fgets(line, sizeof line, stdin)) {
		int got = sscanf(line, "%d %c %lf %lf %lf %lf", &layer, &kind, &x, &SE, &tr, &ts);
		if (got < 3) continue;
		if (got < 6) { tr = -1; ts = -1; }
		if (got < 4) SE = 1.0;
		double r = kind == 'h' ? s.h_calc(layer, x, tr, ts) : kind == 'k' ? s.k_calc(layer, x, SE, tr, ts) : s.get_theta(layer, x, tr, ts);
		printf("%.17g\n", r);
	}
	return 0;
}
MAIN
