#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
//---------------------------------------------------------------------
// program SP
//---------------------------------------------------------------------
//----------
//  Class S:
//----------
//----------
//  Class W:
//----------
//----------
//  Class A:
//----------
//----------
//  Class B:
//----------
//----------
//  Class C:
//----------
//----------
//  Class D:
//----------
//----------
//  Class E:
//----------

struct anon_NAS_SP_c_78 {
   double real;
   double imag;
};

typedef struct anon_NAS_SP_c_78 dcomplex;
/*common /global/*/

int grid_points[3];
int nx2;
int ny2;
int nz2;
/*common /constants/*/

double tx1;
double tx2;
double tx3;
double ty1;
double ty2;
double ty3;
double tz1;
double tz2;
double tz3;
double dx1;
double dx2;
double dx3;
double dx4;
double dx5;
double dy1;
double dy2;
double dy3;
double dy4;
double dy5;
double dz1;
double dz2;
double dz3;
double dz4;
double dz5;
double dssp;
double dt;
double ce[5][13];
double dxmax;
double dymax;
double dzmax;
double xxcon1;
double xxcon2;
double xxcon3;
double xxcon4;
double xxcon5;
double dx1tx1;
double dx2tx1;
double dx3tx1;
double dx4tx1;
double dx5tx1;
double yycon1;
double yycon2;
double yycon3;
double yycon4;
double yycon5;
double dy1ty1;
double dy2ty1;
double dy3ty1;
double dy4ty1;
double dy5ty1;
double zzcon1;
double zzcon2;
double zzcon3;
double zzcon4;
double zzcon5;
double dz1tz1;
double dz2tz1;
double dz3tz1;
double dz4tz1;
double dz5tz1;
double dnxm1;
double dnym1;
double dnzm1;
double c1c2;
double c1c5;
double c3c4;
double c1345;
double conz1;
double c1;
double c2;
double c3;
double c4;
double c5;
double c4dssp;
double c5dssp;
double dtdssp;
double dttx1;
double bt;
double dttx2;
double dtty1;
double dtty2;
double dttz1;
double dttz2;
double c2dttx1;
double c2dtty1;
double c2dttz1;
double comz1;
double comz4;
double comz5;
double comz6;
double c3c4tx3;
double c3c4ty3;
double c3c4tz3;
double c2iv;
double con43;
double con16;
//---------------------------------------------------------------------

// To improve cache performance, grid dimensions padded by 1

// for even number sizes only

//---------------------------------------------------------------------

/*common /fields/*/

double u[36][37][37][5];
double us[36][37][37];
double vs[36][37][37];
double ws[36][37][37];
double qs[36][37][37];
double rho_i[36][37][37];
double speed[36][37][37];
double square[36][37][37];
double rhs[36][37][37][5];
double forcing[36][37][37][5];
//-----------------------------------------------------------------------
// Timer constants
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
int clava_dcg_global[ 42 ] = {0};
void initialize();
void lhsinit(int ni, int nj, double lhs[37][37][5], double lhsp[37][37][5], double lhsm[37][37][5]);
void lhsinitj(int nj, int ni, double lhs[37][37][5], double lhsp[37][37][5], double lhsm[37][37][5]);
void exact_solution(double xi, double eta, double zeta, double dtemp[5]);
void exact_rhs();
void set_constants();
void adi();
void compute_rhs();
void x_solve();
void ninvr();
void y_solve();
void pinvr();
void z_solve();
void tzetar();
void add();
void txinvr();
void error_norm(double rms[5]);
void rhs_norm(double rms[5]);
void verify(int no_time_steps, char *Class, int *verified);
void print_results(char *name, char class, int n1, int n2, int n3, int niter, double t, double mops, char *optype, int verified);
double start[64];
double elapsed[64];
double elapsed_time();
void timer_clear(int n);
void timer_start(int n);
void timer_stop(int n);
double timer_read(int n);
void wtime(double *t);
void clava_call_graph() {
   FILE *log_file_367 = fopen("/home/specs/jbispo/repos/clava-examples/2019-C_Stress_Test/output/NAS_SP.dot", "w+");
   if (log_file_367 == NULL)
   {
       printf("Error opening file /home/specs/jbispo/repos/clava-examples/2019-C_Stress_Test/output/NAS_SP.dot\n");
       exit(1);
   } 
   fprintf(log_file_367, "digraph dynamic_call_graph {\n\n");
   if(clava_dcg_global[0] != 0) {
      fprintf(log_file_367, "	main -> printf [label=\"%d\"];\n", clava_dcg_global[0]);
   }
   if(clava_dcg_global[1] != 0) {
      fprintf(log_file_367, "	main -> set_constants [label=\"%d\"];\n", clava_dcg_global[1]);
   }
   if(clava_dcg_global[2] != 0) {
      fprintf(log_file_367, "	main -> timer_clear [label=\"%d\"];\n", clava_dcg_global[2]);
   }
   if(clava_dcg_global[3] != 0) {
      fprintf(log_file_367, "	main -> exact_rhs [label=\"%d\"];\n", clava_dcg_global[3]);
   }
   if(clava_dcg_global[4] != 0) {
      fprintf(log_file_367, "	main -> initialize [label=\"%d\"];\n", clava_dcg_global[4]);
   }
   if(clava_dcg_global[5] != 0) {
      fprintf(log_file_367, "	main -> adi [label=\"%d\"];\n", clava_dcg_global[5]);
   }
   if(clava_dcg_global[6] != 0) {
      fprintf(log_file_367, "	main -> timer_start [label=\"%d\"];\n", clava_dcg_global[6]);
   }
   if(clava_dcg_global[7] != 0) {
      fprintf(log_file_367, "	main -> timer_stop [label=\"%d\"];\n", clava_dcg_global[7]);
   }
   if(clava_dcg_global[8] != 0) {
      fprintf(log_file_367, "	main -> timer_read [label=\"%d\"];\n", clava_dcg_global[8]);
   }
   if(clava_dcg_global[9] != 0) {
      fprintf(log_file_367, "	main -> verify [label=\"%d\"];\n", clava_dcg_global[9]);
   }
   if(clava_dcg_global[10] != 0) {
      fprintf(log_file_367, "	main -> print_results [label=\"%d\"];\n", clava_dcg_global[10]);
   }
   if(clava_dcg_global[11] != 0) {
      fprintf(log_file_367, "	error_norm -> exact_solution [label=\"%d\"];\n", clava_dcg_global[11]);
   }
   if(clava_dcg_global[12] != 0) {
      fprintf(log_file_367, "	error_norm -> sqrt [label=\"%d\"];\n", clava_dcg_global[12]);
   }
   if(clava_dcg_global[13] != 0) {
      fprintf(log_file_367, "	rhs_norm -> sqrt [label=\"%d\"];\n", clava_dcg_global[13]);
   }
   if(clava_dcg_global[14] != 0) {
      fprintf(log_file_367, "	exact_rhs -> exact_solution [label=\"%d\"];\n", clava_dcg_global[14]);
   }
   if(clava_dcg_global[15] != 0) {
      fprintf(log_file_367, "	adi -> compute_rhs [label=\"%d\"];\n", clava_dcg_global[15]);
   }
   if(clava_dcg_global[16] != 0) {
      fprintf(log_file_367, "	adi -> txinvr [label=\"%d\"];\n", clava_dcg_global[16]);
   }
   if(clava_dcg_global[17] != 0) {
      fprintf(log_file_367, "	adi -> x_solve [label=\"%d\"];\n", clava_dcg_global[17]);
   }
   if(clava_dcg_global[18] != 0) {
      fprintf(log_file_367, "	adi -> y_solve [label=\"%d\"];\n", clava_dcg_global[18]);
   }
   if(clava_dcg_global[19] != 0) {
      fprintf(log_file_367, "	adi -> z_solve [label=\"%d\"];\n", clava_dcg_global[19]);
   }
   if(clava_dcg_global[20] != 0) {
      fprintf(log_file_367, "	adi -> add [label=\"%d\"];\n", clava_dcg_global[20]);
   }
   if(clava_dcg_global[21] != 0) {
      fprintf(log_file_367, "	initialize -> exact_solution [label=\"%d\"];\n", clava_dcg_global[21]);
   }
   if(clava_dcg_global[22] != 0) {
      fprintf(log_file_367, "	compute_rhs -> sqrt [label=\"%d\"];\n", clava_dcg_global[22]);
   }
   if(clava_dcg_global[23] != 0) {
      fprintf(log_file_367, "	verify -> error_norm [label=\"%d\"];\n", clava_dcg_global[23]);
   }
   if(clava_dcg_global[24] != 0) {
      fprintf(log_file_367, "	verify -> compute_rhs [label=\"%d\"];\n", clava_dcg_global[24]);
   }
   if(clava_dcg_global[25] != 0) {
      fprintf(log_file_367, "	verify -> rhs_norm [label=\"%d\"];\n", clava_dcg_global[25]);
   }
   if(clava_dcg_global[26] != 0) {
      fprintf(log_file_367, "	verify -> fabs [label=\"%d\"];\n", clava_dcg_global[26]);
   }
   if(clava_dcg_global[27] != 0) {
      fprintf(log_file_367, "	verify -> printf [label=\"%d\"];\n", clava_dcg_global[27]);
   }
   if(clava_dcg_global[28] != 0) {
      fprintf(log_file_367, "	x_solve -> lhsinit [label=\"%d\"];\n", clava_dcg_global[28]);
   }
   if(clava_dcg_global[29] != 0) {
      fprintf(log_file_367, "	x_solve -> ninvr [label=\"%d\"];\n", clava_dcg_global[29]);
   }
   if(clava_dcg_global[30] != 0) {
      fprintf(log_file_367, "	y_solve -> lhsinitj [label=\"%d\"];\n", clava_dcg_global[30]);
   }
   if(clava_dcg_global[31] != 0) {
      fprintf(log_file_367, "	y_solve -> pinvr [label=\"%d\"];\n", clava_dcg_global[31]);
   }
   if(clava_dcg_global[32] != 0) {
      fprintf(log_file_367, "	z_solve -> lhsinitj [label=\"%d\"];\n", clava_dcg_global[32]);
   }
   if(clava_dcg_global[33] != 0) {
      fprintf(log_file_367, "	z_solve -> tzetar [label=\"%d\"];\n", clava_dcg_global[33]);
   }
   if(clava_dcg_global[34] != 0) {
      fprintf(log_file_367, "	set_constants -> sqrt [label=\"%d\"];\n", clava_dcg_global[34]);
   }
   if(clava_dcg_global[35] != 0) {
      fprintf(log_file_367, "	print_results -> printf [label=\"%d\"];\n", clava_dcg_global[35]);
   }
   if(clava_dcg_global[36] != 0) {
      fprintf(log_file_367, "	print_results -> sprintf [label=\"%d\"];\n", clava_dcg_global[36]);
   }
   if(clava_dcg_global[37] != 0) {
      fprintf(log_file_367, "	print_results -> pow [label=\"%d\"];\n", clava_dcg_global[37]);
   }
   if(clava_dcg_global[38] != 0) {
      fprintf(log_file_367, "	wtime -> gettimeofday [label=\"%d\"];\n", clava_dcg_global[38]);
   }
   if(clava_dcg_global[39] != 0) {
      fprintf(log_file_367, "	elapsed_time -> wtime [label=\"%d\"];\n", clava_dcg_global[39]);
   }
   if(clava_dcg_global[40] != 0) {
      fprintf(log_file_367, "	timer_start -> elapsed_time [label=\"%d\"];\n", clava_dcg_global[40]);
   }
   if(clava_dcg_global[41] != 0) {
      fprintf(log_file_367, "	timer_stop -> elapsed_time [label=\"%d\"];\n", clava_dcg_global[41]);
   }
   fprintf(log_file_367, "}\n");
   fclose(log_file_367);
}

int main(int argc, char *argv[]) {
   atexit(clava_call_graph);
   int i, niter, step, n3;
   double mflops;
   double t;
   double tmax;
   double trecs[16];
   int verified;
   char Class;
   char *t_names[16];
   clava_dcg_global[ 0 ]++;
   printf("\n\n NAS Parallel Benchmarks (NPB3.3-SER-C) - SP Benchmark\n\n");
   niter = 400;
   dt = 0.0015;
   grid_points[0] = 36;
   grid_points[1] = 36;
   grid_points[2] = 36;
   clava_dcg_global[ 0 ]++;
   printf(" Size: %4dx%4dx%4d\n", grid_points[0], grid_points[1], grid_points[2]);
   clava_dcg_global[ 0 ]++;
   printf(" Iterations: %4d    dt: %10.6f\n", niter, dt);
   clava_dcg_global[ 0 ]++;
   printf("\n");
   if((grid_points[0] > 36) || (grid_points[1] > 36) || (grid_points[2] > 36)) {
      clava_dcg_global[ 0 ]++;
      printf(" %d, %d, %d\n", grid_points[0], grid_points[1], grid_points[2]);
      clava_dcg_global[ 0 ]++;
      printf(" Problem size too big for compiled array sizes\n");
      
      return 0;
   }
   nx2 = grid_points[0] - 2;
   ny2 = grid_points[1] - 2;
   nz2 = grid_points[2] - 2;
   clava_dcg_global[ 1 ]++;
   set_constants();
   for(i = 1; i <= 15; i++) {
      clava_dcg_global[ 2 ]++;
      timer_clear(i);
   }
   clava_dcg_global[ 3 ]++;
   exact_rhs();
   clava_dcg_global[ 4 ]++;
   initialize();
   //---------------------------------------------------------------------
   // do one time step to touch all code, and reinitialize
   //---------------------------------------------------------------------
   clava_dcg_global[ 5 ]++;
   adi();
   clava_dcg_global[ 4 ]++;
   initialize();
   for(i = 1; i <= 15; i++) {
      clava_dcg_global[ 2 ]++;
      timer_clear(i);
   }
   clava_dcg_global[ 6 ]++;
   timer_start(1);
   for(step = 1; step <= niter; step++) {
      if((step % 20) == 0 || step == 1) {
         clava_dcg_global[ 0 ]++;
         printf(" Time step %4d\n", step);
      }
      clava_dcg_global[ 5 ]++;
      adi();
   }
   clava_dcg_global[ 7 ]++;
   timer_stop(1);
   clava_dcg_global[ 8 ]++;
   tmax = timer_read(1);
   clava_dcg_global[ 9 ]++;
   verify(niter, &Class, &verified);
   if(tmax != 0.0) {
      n3 = grid_points[0] * grid_points[1] * grid_points[2];
      t = (grid_points[0] + grid_points[1] + grid_points[2]) / 3.0;
      mflops = (881.174 * (double) n3 - 4683.91 * (t * t) + 11484.5 * t - 19272.4) * (double) niter / (tmax * 1000000.0);
   }
   else {
      mflops = 0.0;
   }
   clava_dcg_global[ 10 ]++;
   print_results("SP", Class, grid_points[0], grid_points[1], grid_points[2], niter, tmax, mflops, "          floating point", verified);
   int exitValue = verified ? 0 : 1;
   
   return exitValue;
}

//---------------------------------------------------------------------
// this function computes the norm of the difference between the
// computed solution and the exact solution
//---------------------------------------------------------------------
void error_norm(double rms[5]) {
   int i, j, k, m, d;
   double xi;
   double eta;
   double zeta;
   double u_exact[5];
   double add;
   for(m = 0; m < 5; m++) {
      rms[m] = 0.0;
   }
   for(k = 0; k <= grid_points[2] - 1; k++) {
      zeta = (double) k * dnzm1;
      for(j = 0; j <= grid_points[1] - 1; j++) {
         eta = (double) j * dnym1;
         for(i = 0; i <= grid_points[0] - 1; i++) {
            xi = (double) i * dnxm1;
            clava_dcg_global[ 11 ]++;
            exact_solution(xi, eta, zeta, u_exact);
            for(m = 0; m < 5; m++) {
               add = u[k][j][i][m] - u_exact[m];
               rms[m] = rms[m] + add * add;
            }
         }
      }
   }
   for(m = 0; m < 5; m++) {
      for(d = 0; d < 3; d++) {
         rms[m] = rms[m] / (double) (grid_points[d] - 2);
      }
      clava_dcg_global[ 12 ]++;
      rms[m] = sqrt(rms[m]);
   }
}

void rhs_norm(double rms[5]) {
   int i, j, k, d, m;
   double add;
   for(m = 0; m < 5; m++) {
      rms[m] = 0.0;
   }
   for(k = 1; k <= nz2; k++) {
      for(j = 1; j <= ny2; j++) {
         for(i = 1; i <= nx2; i++) {
            for(m = 0; m < 5; m++) {
               add = rhs[k][j][i][m];
               rms[m] = rms[m] + add * add;
            }
         }
      }
   }
   for(m = 0; m < 5; m++) {
      for(d = 0; d < 3; d++) {
         rms[m] = rms[m] / (double) (grid_points[d] - 2);
      }
      clava_dcg_global[ 13 ]++;
      rms[m] = sqrt(rms[m]);
   }
}

//---------------------------------------------------------------------
// compute the right hand side based on exact solution
//---------------------------------------------------------------------
void exact_rhs() {
   double dtemp[5];
   double xi;
   double eta;
   double zeta;
   double dtpp;
   int m, i, j, k, ip1, im1, jp1, jm1, km1, kp1;
   double ue[36][5];
   double buf[36][5];
   double q[36];
   double cuf[36];
   //---------------------------------------------------------------------
   // initialize
   //---------------------------------------------------------------------
   for(k = 0; k <= grid_points[2] - 1; k++) {
      for(j = 0; j <= grid_points[1] - 1; j++) {
         for(i = 0; i <= grid_points[0] - 1; i++) {
            for(m = 0; m < 5; m++) {
               forcing[k][j][i][m] = 0.0;
            }
         }
      }
   }
   //---------------------------------------------------------------------
   // xi-direction flux differences
   //---------------------------------------------------------------------
   for(k = 1; k <= grid_points[2] - 2; k++) {
      zeta = (double) k * dnzm1;
      for(j = 1; j <= grid_points[1] - 2; j++) {
         eta = (double) j * dnym1;
         for(i = 0; i <= grid_points[0] - 1; i++) {
            xi = (double) i * dnxm1;
            clava_dcg_global[ 14 ]++;
            exact_solution(xi, eta, zeta, dtemp);
            for(m = 0; m < 5; m++) {
               ue[i][m] = dtemp[m];
            }
            dtpp = 1.0 / dtemp[0];
            for(m = 1; m < 5; m++) {
               buf[i][m] = dtpp * dtemp[m];
            }
            cuf[i] = buf[i][1] * buf[i][1];
            buf[i][0] = cuf[i] + buf[i][2] * buf[i][2] + buf[i][3] * buf[i][3];
            q[i] = 0.5 * (buf[i][1] * ue[i][1] + buf[i][2] * ue[i][2] + buf[i][3] * ue[i][3]);
         }
         for(i = 1; i <= grid_points[0] - 2; i++) {
            im1 = i - 1;
            ip1 = i + 1;
            forcing[k][j][i][0] = forcing[k][j][i][0] - tx2 * (ue[ip1][1] - ue[im1][1]) + dx1tx1 * (ue[ip1][0] - 2.0 * ue[i][0] + ue[im1][0]);
            forcing[k][j][i][1] = forcing[k][j][i][1] - tx2 * ((ue[ip1][1] * buf[ip1][1] + c2 * (ue[ip1][4] - q[ip1])) - (ue[im1][1] * buf[im1][1] + c2 * (ue[im1][4] - q[im1]))) + xxcon1 * (buf[ip1][1] - 2.0 * buf[i][1] + buf[im1][1]) + dx2tx1 * (ue[ip1][1] - 2.0 * ue[i][1] + ue[im1][1]);
            forcing[k][j][i][2] = forcing[k][j][i][2] - tx2 * (ue[ip1][2] * buf[ip1][1] - ue[im1][2] * buf[im1][1]) + xxcon2 * (buf[ip1][2] - 2.0 * buf[i][2] + buf[im1][2]) + dx3tx1 * (ue[ip1][2] - 2.0 * ue[i][2] + ue[im1][2]);
            forcing[k][j][i][3] = forcing[k][j][i][3] - tx2 * (ue[ip1][3] * buf[ip1][1] - ue[im1][3] * buf[im1][1]) + xxcon2 * (buf[ip1][3] - 2.0 * buf[i][3] + buf[im1][3]) + dx4tx1 * (ue[ip1][3] - 2.0 * ue[i][3] + ue[im1][3]);
            forcing[k][j][i][4] = forcing[k][j][i][4] - tx2 * (buf[ip1][1] * (c1 * ue[ip1][4] - c2 * q[ip1]) - buf[im1][1] * (c1 * ue[im1][4] - c2 * q[im1])) + 0.5 * xxcon3 * (buf[ip1][0] - 2.0 * buf[i][0] + buf[im1][0]) + xxcon4 * (cuf[ip1] - 2.0 * cuf[i] + cuf[im1]) + xxcon5 * (buf[ip1][4] - 2.0 * buf[i][4] + buf[im1][4]) + dx5tx1 * (ue[ip1][4] - 2.0 * ue[i][4] + ue[im1][4]);
         }
         //---------------------------------------------------------------------
         // Fourth-order dissipation
         //---------------------------------------------------------------------
         for(m = 0; m < 5; m++) {
            i = 1;
            forcing[k][j][i][m] = forcing[k][j][i][m] - dssp * (5.0 * ue[i][m] - 4.0 * ue[i + 1][m] + ue[i + 2][m]);
            i = 2;
            forcing[k][j][i][m] = forcing[k][j][i][m] - dssp * (-4.0 * ue[i - 1][m] + 6.0 * ue[i][m] - 4.0 * ue[i + 1][m] + ue[i + 2][m]);
         }
         for(i = 3; i <= grid_points[0] - 4; i++) {
            for(m = 0; m < 5; m++) {
               forcing[k][j][i][m] = forcing[k][j][i][m] - dssp * (ue[i - 2][m] - 4.0 * ue[i - 1][m] + 6.0 * ue[i][m] - 4.0 * ue[i + 1][m] + ue[i + 2][m]);
            }
         }
         for(m = 0; m < 5; m++) {
            i = grid_points[0] - 3;
            forcing[k][j][i][m] = forcing[k][j][i][m] - dssp * (ue[i - 2][m] - 4.0 * ue[i - 1][m] + 6.0 * ue[i][m] - 4.0 * ue[i + 1][m]);
            i = grid_points[0] - 2;
            forcing[k][j][i][m] = forcing[k][j][i][m] - dssp * (ue[i - 2][m] - 4.0 * ue[i - 1][m] + 5.0 * ue[i][m]);
         }
      }
   }
   //---------------------------------------------------------------------
   // eta-direction flux differences
   //---------------------------------------------------------------------
   for(k = 1; k <= grid_points[2] - 2; k++) {
      zeta = (double) k * dnzm1;
      for(i = 1; i <= grid_points[0] - 2; i++) {
         xi = (double) i * dnxm1;
         for(j = 0; j <= grid_points[1] - 1; j++) {
            eta = (double) j * dnym1;
            clava_dcg_global[ 14 ]++;
            exact_solution(xi, eta, zeta, dtemp);
            for(m = 0; m < 5; m++) {
               ue[j][m] = dtemp[m];
            }
            dtpp = 1.0 / dtemp[0];
            for(m = 1; m < 5; m++) {
               buf[j][m] = dtpp * dtemp[m];
            }
            cuf[j] = buf[j][2] * buf[j][2];
            buf[j][0] = cuf[j] + buf[j][1] * buf[j][1] + buf[j][3] * buf[j][3];
            q[j] = 0.5 * (buf[j][1] * ue[j][1] + buf[j][2] * ue[j][2] + buf[j][3] * ue[j][3]);
         }
         for(j = 1; j <= grid_points[1] - 2; j++) {
            jm1 = j - 1;
            jp1 = j + 1;
            forcing[k][j][i][0] = forcing[k][j][i][0] - ty2 * (ue[jp1][2] - ue[jm1][2]) + dy1ty1 * (ue[jp1][0] - 2.0 * ue[j][0] + ue[jm1][0]);
            forcing[k][j][i][1] = forcing[k][j][i][1] - ty2 * (ue[jp1][1] * buf[jp1][2] - ue[jm1][1] * buf[jm1][2]) + yycon2 * (buf[jp1][1] - 2.0 * buf[j][1] + buf[jm1][1]) + dy2ty1 * (ue[jp1][1] - 2.0 * ue[j][1] + ue[jm1][1]);
            forcing[k][j][i][2] = forcing[k][j][i][2] - ty2 * ((ue[jp1][2] * buf[jp1][2] + c2 * (ue[jp1][4] - q[jp1])) - (ue[jm1][2] * buf[jm1][2] + c2 * (ue[jm1][4] - q[jm1]))) + yycon1 * (buf[jp1][2] - 2.0 * buf[j][2] + buf[jm1][2]) + dy3ty1 * (ue[jp1][2] - 2.0 * ue[j][2] + ue[jm1][2]);
            forcing[k][j][i][3] = forcing[k][j][i][3] - ty2 * (ue[jp1][3] * buf[jp1][2] - ue[jm1][3] * buf[jm1][2]) + yycon2 * (buf[jp1][3] - 2.0 * buf[j][3] + buf[jm1][3]) + dy4ty1 * (ue[jp1][3] - 2.0 * ue[j][3] + ue[jm1][3]);
            forcing[k][j][i][4] = forcing[k][j][i][4] - ty2 * (buf[jp1][2] * (c1 * ue[jp1][4] - c2 * q[jp1]) - buf[jm1][2] * (c1 * ue[jm1][4] - c2 * q[jm1])) + 0.5 * yycon3 * (buf[jp1][0] - 2.0 * buf[j][0] + buf[jm1][0]) + yycon4 * (cuf[jp1] - 2.0 * cuf[j] + cuf[jm1]) + yycon5 * (buf[jp1][4] - 2.0 * buf[j][4] + buf[jm1][4]) + dy5ty1 * (ue[jp1][4] - 2.0 * ue[j][4] + ue[jm1][4]);
         }
         //---------------------------------------------------------------------
         // Fourth-order dissipation
         //---------------------------------------------------------------------
         for(m = 0; m < 5; m++) {
            j = 1;
            forcing[k][j][i][m] = forcing[k][j][i][m] - dssp * (5.0 * ue[j][m] - 4.0 * ue[j + 1][m] + ue[j + 2][m]);
            j = 2;
            forcing[k][j][i][m] = forcing[k][j][i][m] - dssp * (-4.0 * ue[j - 1][m] + 6.0 * ue[j][m] - 4.0 * ue[j + 1][m] + ue[j + 2][m]);
         }
         for(j = 3; j <= grid_points[1] - 4; j++) {
            for(m = 0; m < 5; m++) {
               forcing[k][j][i][m] = forcing[k][j][i][m] - dssp * (ue[j - 2][m] - 4.0 * ue[j - 1][m] + 6.0 * ue[j][m] - 4.0 * ue[j + 1][m] + ue[j + 2][m]);
            }
         }
         for(m = 0; m < 5; m++) {
            j = grid_points[1] - 3;
            forcing[k][j][i][m] = forcing[k][j][i][m] - dssp * (ue[j - 2][m] - 4.0 * ue[j - 1][m] + 6.0 * ue[j][m] - 4.0 * ue[j + 1][m]);
            j = grid_points[1] - 2;
            forcing[k][j][i][m] = forcing[k][j][i][m] - dssp * (ue[j - 2][m] - 4.0 * ue[j - 1][m] + 5.0 * ue[j][m]);
         }
      }
   }
   //---------------------------------------------------------------------
   // zeta-direction flux differences
   //---------------------------------------------------------------------
   for(j = 1; j <= grid_points[1] - 2; j++) {
      eta = (double) j * dnym1;
      for(i = 1; i <= grid_points[0] - 2; i++) {
         xi = (double) i * dnxm1;
         for(k = 0; k <= grid_points[2] - 1; k++) {
            zeta = (double) k * dnzm1;
            clava_dcg_global[ 14 ]++;
            exact_solution(xi, eta, zeta, dtemp);
            for(m = 0; m < 5; m++) {
               ue[k][m] = dtemp[m];
            }
            dtpp = 1.0 / dtemp[0];
            for(m = 1; m < 5; m++) {
               buf[k][m] = dtpp * dtemp[m];
            }
            cuf[k] = buf[k][3] * buf[k][3];
            buf[k][0] = cuf[k] + buf[k][1] * buf[k][1] + buf[k][2] * buf[k][2];
            q[k] = 0.5 * (buf[k][1] * ue[k][1] + buf[k][2] * ue[k][2] + buf[k][3] * ue[k][3]);
         }
         for(k = 1; k <= grid_points[2] - 2; k++) {
            km1 = k - 1;
            kp1 = k + 1;
            forcing[k][j][i][0] = forcing[k][j][i][0] - tz2 * (ue[kp1][3] - ue[km1][3]) + dz1tz1 * (ue[kp1][0] - 2.0 * ue[k][0] + ue[km1][0]);
            forcing[k][j][i][1] = forcing[k][j][i][1] - tz2 * (ue[kp1][1] * buf[kp1][3] - ue[km1][1] * buf[km1][3]) + zzcon2 * (buf[kp1][1] - 2.0 * buf[k][1] + buf[km1][1]) + dz2tz1 * (ue[kp1][1] - 2.0 * ue[k][1] + ue[km1][1]);
            forcing[k][j][i][2] = forcing[k][j][i][2] - tz2 * (ue[kp1][2] * buf[kp1][3] - ue[km1][2] * buf[km1][3]) + zzcon2 * (buf[kp1][2] - 2.0 * buf[k][2] + buf[km1][2]) + dz3tz1 * (ue[kp1][2] - 2.0 * ue[k][2] + ue[km1][2]);
            forcing[k][j][i][3] = forcing[k][j][i][3] - tz2 * ((ue[kp1][3] * buf[kp1][3] + c2 * (ue[kp1][4] - q[kp1])) - (ue[km1][3] * buf[km1][3] + c2 * (ue[km1][4] - q[km1]))) + zzcon1 * (buf[kp1][3] - 2.0 * buf[k][3] + buf[km1][3]) + dz4tz1 * (ue[kp1][3] - 2.0 * ue[k][3] + ue[km1][3]);
            forcing[k][j][i][4] = forcing[k][j][i][4] - tz2 * (buf[kp1][3] * (c1 * ue[kp1][4] - c2 * q[kp1]) - buf[km1][3] * (c1 * ue[km1][4] - c2 * q[km1])) + 0.5 * zzcon3 * (buf[kp1][0] - 2.0 * buf[k][0] + buf[km1][0]) + zzcon4 * (cuf[kp1] - 2.0 * cuf[k] + cuf[km1]) + zzcon5 * (buf[kp1][4] - 2.0 * buf[k][4] + buf[km1][4]) + dz5tz1 * (ue[kp1][4] - 2.0 * ue[k][4] + ue[km1][4]);
         }
         //---------------------------------------------------------------------
         // Fourth-order dissipation
         //---------------------------------------------------------------------
         for(m = 0; m < 5; m++) {
            k = 1;
            forcing[k][j][i][m] = forcing[k][j][i][m] - dssp * (5.0 * ue[k][m] - 4.0 * ue[k + 1][m] + ue[k + 2][m]);
            k = 2;
            forcing[k][j][i][m] = forcing[k][j][i][m] - dssp * (-4.0 * ue[k - 1][m] + 6.0 * ue[k][m] - 4.0 * ue[k + 1][m] + ue[k + 2][m]);
         }
         for(k = 3; k <= grid_points[2] - 4; k++) {
            for(m = 0; m < 5; m++) {
               forcing[k][j][i][m] = forcing[k][j][i][m] - dssp * (ue[k - 2][m] - 4.0 * ue[k - 1][m] + 6.0 * ue[k][m] - 4.0 * ue[k + 1][m] + ue[k + 2][m]);
            }
         }
         for(m = 0; m < 5; m++) {
            k = grid_points[2] - 3;
            forcing[k][j][i][m] = forcing[k][j][i][m] - dssp * (ue[k - 2][m] - 4.0 * ue[k - 1][m] + 6.0 * ue[k][m] - 4.0 * ue[k + 1][m]);
            k = grid_points[2] - 2;
            forcing[k][j][i][m] = forcing[k][j][i][m] - dssp * (ue[k - 2][m] - 4.0 * ue[k - 1][m] + 5.0 * ue[k][m]);
         }
      }
   }
   //---------------------------------------------------------------------
   // now change the sign of the forcing function,
   //---------------------------------------------------------------------
   for(k = 1; k <= grid_points[2] - 2; k++) {
      for(j = 1; j <= grid_points[1] - 2; j++) {
         for(i = 1; i <= grid_points[0] - 2; i++) {
            for(m = 0; m < 5; m++) {
               forcing[k][j][i][m] = -1.0 * forcing[k][j][i][m];
            }
         }
      }
   }
}

//---------------------------------------------------------------------
// this function returns the exact solution at point xi, eta, zeta
//---------------------------------------------------------------------
void exact_solution(double xi, double eta, double zeta, double dtemp[5]) {
   int m;
   for(m = 0; m < 5; m++) {
      dtemp[m] = ce[m][0] + xi * (ce[m][1] + xi * (ce[m][4] + xi * (ce[m][7] + xi * ce[m][10]))) + eta * (ce[m][2] + eta * (ce[m][5] + eta * (ce[m][8] + eta * ce[m][11]))) + zeta * (ce[m][3] + zeta * (ce[m][6] + zeta * (ce[m][9] + zeta * ce[m][12])));
   }
}

void adi() {
   clava_dcg_global[ 15 ]++;
   compute_rhs();
   clava_dcg_global[ 16 ]++;
   txinvr();
   clava_dcg_global[ 17 ]++;
   x_solve();
   clava_dcg_global[ 18 ]++;
   y_solve();
   clava_dcg_global[ 19 ]++;
   z_solve();
   clava_dcg_global[ 20 ]++;
   add();
}

//---------------------------------------------------------------------
// addition of update to the vector u
//---------------------------------------------------------------------
void add() {
   int i, j, k, m;
   for(k = 1; k <= nz2; k++) {
      for(j = 1; j <= ny2; j++) {
         for(i = 1; i <= nx2; i++) {
            for(m = 0; m < 5; m++) {
               u[k][j][i][m] = u[k][j][i][m] + rhs[k][j][i][m];
            }
         }
      }
   }
}

//---------------------------------------------------------------------
// This subroutine initializes the field variable u using
// tri-linear transfinite interpolation of the boundary values
//---------------------------------------------------------------------
void initialize() {
   int i, j, k, m, ix, iy, iz;
   double xi;
   double eta;
   double zeta;
   double Pface[2][3][5];
   double Pxi;
   double Peta;
   double Pzeta;
   double temp[5];
   //---------------------------------------------------------------------
   //  Later (in compute_rhs) we compute 1/u for every element. A few of
   //  the corner elements are not used, but it convenient (and faster)
   //  to compute the whole thing with a simple loop. Make sure those
   //  values are nonzero by initializing the whole thing here.
   //---------------------------------------------------------------------
   for(k = 0; k <= grid_points[2] - 1; k++) {
      for(j = 0; j <= grid_points[1] - 1; j++) {
         for(i = 0; i <= grid_points[0] - 1; i++) {
            u[k][j][i][0] = 1.0;
            u[k][j][i][1] = 0.0;
            u[k][j][i][2] = 0.0;
            u[k][j][i][3] = 0.0;
            u[k][j][i][4] = 1.0;
         }
      }
   }
   //---------------------------------------------------------------------
   // first store the "interpolated" values everywhere on the grid
   //---------------------------------------------------------------------
   for(k = 0; k <= grid_points[2] - 1; k++) {
      zeta = (double) k * dnzm1;
      for(j = 0; j <= grid_points[1] - 1; j++) {
         eta = (double) j * dnym1;
         for(i = 0; i <= grid_points[0] - 1; i++) {
            xi = (double) i * dnxm1;
            for(ix = 0; ix < 2; ix++) {
               Pxi = (double) ix;
               clava_dcg_global[ 21 ]++;
               exact_solution(Pxi, eta, zeta, &Pface[ix][0][0]);
            }
            for(iy = 0; iy < 2; iy++) {
               Peta = (double) iy;
               clava_dcg_global[ 21 ]++;
               exact_solution(xi, Peta, zeta, &Pface[iy][1][0]);
            }
            for(iz = 0; iz < 2; iz++) {
               Pzeta = (double) iz;
               clava_dcg_global[ 21 ]++;
               exact_solution(xi, eta, Pzeta, &Pface[iz][2][0]);
            }
            for(m = 0; m < 5; m++) {
               Pxi = xi * Pface[1][0][m] + (1.0 - xi) * Pface[0][0][m];
               Peta = eta * Pface[1][1][m] + (1.0 - eta) * Pface[0][1][m];
               Pzeta = zeta * Pface[1][2][m] + (1.0 - zeta) * Pface[0][2][m];
               u[k][j][i][m] = Pxi + Peta + Pzeta - Pxi * Peta - Pxi * Pzeta - Peta * Pzeta + Pxi * Peta * Pzeta;
            }
         }
      }
   }
   //---------------------------------------------------------------------
   // now store the exact values on the boundaries
   //---------------------------------------------------------------------
   //---------------------------------------------------------------------
   // west face
   //---------------------------------------------------------------------
   xi = 0.0;
   i = 0;
   for(k = 0; k <= grid_points[2] - 1; k++) {
      zeta = (double) k * dnzm1;
      for(j = 0; j <= grid_points[1] - 1; j++) {
         eta = (double) j * dnym1;
         clava_dcg_global[ 21 ]++;
         exact_solution(xi, eta, zeta, temp);
         for(m = 0; m < 5; m++) {
            u[k][j][i][m] = temp[m];
         }
      }
   }
   //---------------------------------------------------------------------
   // east face
   //---------------------------------------------------------------------
   xi = 1.0;
   i = grid_points[0] - 1;
   for(k = 0; k <= grid_points[2] - 1; k++) {
      zeta = (double) k * dnzm1;
      for(j = 0; j <= grid_points[1] - 1; j++) {
         eta = (double) j * dnym1;
         clava_dcg_global[ 21 ]++;
         exact_solution(xi, eta, zeta, temp);
         for(m = 0; m < 5; m++) {
            u[k][j][i][m] = temp[m];
         }
      }
   }
   //---------------------------------------------------------------------
   // south face
   //---------------------------------------------------------------------
   eta = 0.0;
   j = 0;
   for(k = 0; k <= grid_points[2] - 1; k++) {
      zeta = (double) k * dnzm1;
      for(i = 0; i <= grid_points[0] - 1; i++) {
         xi = (double) i * dnxm1;
         clava_dcg_global[ 21 ]++;
         exact_solution(xi, eta, zeta, temp);
         for(m = 0; m < 5; m++) {
            u[k][j][i][m] = temp[m];
         }
      }
   }
   //---------------------------------------------------------------------
   // north face
   //---------------------------------------------------------------------
   eta = 1.0;
   j = grid_points[1] - 1;
   for(k = 0; k <= grid_points[2] - 1; k++) {
      zeta = (double) k * dnzm1;
      for(i = 0; i <= grid_points[0] - 1; i++) {
         xi = (double) i * dnxm1;
         clava_dcg_global[ 21 ]++;
         exact_solution(xi, eta, zeta, temp);
         for(m = 0; m < 5; m++) {
            u[k][j][i][m] = temp[m];
         }
      }
   }
   //---------------------------------------------------------------------
   // bottom face
   //---------------------------------------------------------------------
   zeta = 0.0;
   k = 0;
   for(j = 0; j <= grid_points[1] - 1; j++) {
      eta = (double) j * dnym1;
      for(i = 0; i <= grid_points[0] - 1; i++) {
         xi = (double) i * dnxm1;
         clava_dcg_global[ 21 ]++;
         exact_solution(xi, eta, zeta, temp);
         for(m = 0; m < 5; m++) {
            u[k][j][i][m] = temp[m];
         }
      }
   }
   //---------------------------------------------------------------------
   // top face
   //---------------------------------------------------------------------
   zeta = 1.0;
   k = grid_points[2] - 1;
   for(j = 0; j <= grid_points[1] - 1; j++) {
      eta = (double) j * dnym1;
      for(i = 0; i <= grid_points[0] - 1; i++) {
         xi = (double) i * dnxm1;
         clava_dcg_global[ 21 ]++;
         exact_solution(xi, eta, zeta, temp);
         for(m = 0; m < 5; m++) {
            u[k][j][i][m] = temp[m];
         }
      }
   }
}

void lhsinit(int ni, int nj, double lhs[37][37][5], double lhsp[37][37][5], double lhsm[37][37][5]) {
   int j, m;
   //---------------------------------------------------------------------
   // zap the whole left hand side for starters
   // set all diagonal values to 1. This is overkill, but convenient
   //---------------------------------------------------------------------
   for(j = 1; j <= nj; j++) {
      for(m = 0; m < 5; m++) {
         lhs[j][0][m] = 0.0;
         lhsp[j][0][m] = 0.0;
         lhsm[j][0][m] = 0.0;
         lhs[j][ni][m] = 0.0;
         lhsp[j][ni][m] = 0.0;
         lhsm[j][ni][m] = 0.0;
      }
      lhs[j][0][2] = 1.0;
      lhsp[j][0][2] = 1.0;
      lhsm[j][0][2] = 1.0;
      lhs[j][ni][2] = 1.0;
      lhsp[j][ni][2] = 1.0;
      lhsm[j][ni][2] = 1.0;
   }
}

void lhsinitj(int nj, int ni, double lhs[37][37][5], double lhsp[37][37][5], double lhsm[37][37][5]) {
   int i, m;
   //---------------------------------------------------------------------
   // zap the whole left hand side for starters
   // set all diagonal values to 1. This is overkill, but convenient
   //---------------------------------------------------------------------
   for(i = 1; i <= ni; i++) {
      for(m = 0; m < 5; m++) {
         lhs[0][i][m] = 0.0;
         lhsp[0][i][m] = 0.0;
         lhsm[0][i][m] = 0.0;
         lhs[nj][i][m] = 0.0;
         lhsp[nj][i][m] = 0.0;
         lhsm[nj][i][m] = 0.0;
      }
      lhs[0][i][2] = 1.0;
      lhsp[0][i][2] = 1.0;
      lhsm[0][i][2] = 1.0;
      lhs[nj][i][2] = 1.0;
      lhsp[nj][i][2] = 1.0;
      lhsm[nj][i][2] = 1.0;
   }
}

//---------------------------------------------------------------------
// block-diagonal matrix-vector multiplication
//---------------------------------------------------------------------
void ninvr() {
   int i, j, k;
   double r1, r2, r3, r4, r5, t1, t2;
   for(k = 1; k <= nz2; k++) {
      for(j = 1; j <= ny2; j++) {
         for(i = 1; i <= nx2; i++) {
            r1 = rhs[k][j][i][0];
            r2 = rhs[k][j][i][1];
            r3 = rhs[k][j][i][2];
            r4 = rhs[k][j][i][3];
            r5 = rhs[k][j][i][4];
            t1 = bt * r3;
            t2 = 0.5 * (r4 + r5);
            rhs[k][j][i][0] = -r2;
            rhs[k][j][i][1] = r1;
            rhs[k][j][i][2] = bt * (r4 - r5);
            rhs[k][j][i][3] = -t1 + t2;
            rhs[k][j][i][4] = t1 + t2;
         }
      }
   }
}

//---------------------------------------------------------------------
// block-diagonal matrix-vector multiplication
//---------------------------------------------------------------------
void pinvr() {
   int i, j, k;
   double r1, r2, r3, r4, r5, t1, t2;
   for(k = 1; k <= nz2; k++) {
      for(j = 1; j <= ny2; j++) {
         for(i = 1; i <= nx2; i++) {
            r1 = rhs[k][j][i][0];
            r2 = rhs[k][j][i][1];
            r3 = rhs[k][j][i][2];
            r4 = rhs[k][j][i][3];
            r5 = rhs[k][j][i][4];
            t1 = bt * r1;
            t2 = 0.5 * (r4 + r5);
            rhs[k][j][i][0] = bt * (r4 - r5);
            rhs[k][j][i][1] = -r3;
            rhs[k][j][i][2] = r2;
            rhs[k][j][i][3] = -t1 + t2;
            rhs[k][j][i][4] = t1 + t2;
         }
      }
   }
}

void compute_rhs() {
   int i, j, k, m;
   double aux, rho_inv, uijk, up1, um1, vijk, vp1, vm1, wijk, wp1, wm1;
   //---------------------------------------------------------------------
   // compute the reciprocal of density, and the kinetic energy,
   // and the speed of sound.
   //---------------------------------------------------------------------
   for(k = 0; k <= grid_points[2] - 1; k++) {
      for(j = 0; j <= grid_points[1] - 1; j++) {
         for(i = 0; i <= grid_points[0] - 1; i++) {
            rho_inv = 1.0 / u[k][j][i][0];
            rho_i[k][j][i] = rho_inv;
            us[k][j][i] = u[k][j][i][1] * rho_inv;
            vs[k][j][i] = u[k][j][i][2] * rho_inv;
            ws[k][j][i] = u[k][j][i][3] * rho_inv;
            square[k][j][i] = 0.5 * (u[k][j][i][1] * u[k][j][i][1] + u[k][j][i][2] * u[k][j][i][2] + u[k][j][i][3] * u[k][j][i][3]) * rho_inv;
            qs[k][j][i] = square[k][j][i] * rho_inv;
            //-------------------------------------------------------------------
            // (don't need speed and ainx until the lhs computation)
            //-------------------------------------------------------------------
            aux = c1c2 * rho_inv * (u[k][j][i][4] - square[k][j][i]);
            clava_dcg_global[ 22 ]++;
            speed[k][j][i] = sqrt(aux);
         }
      }
   }
   //---------------------------------------------------------------------
   // copy the exact forcing term to the right hand side;  because
   // this forcing term is known, we can store it on the whole grid
   // including the boundary
   //---------------------------------------------------------------------
   for(k = 0; k <= grid_points[2] - 1; k++) {
      for(j = 0; j <= grid_points[1] - 1; j++) {
         for(i = 0; i <= grid_points[0] - 1; i++) {
            for(m = 0; m < 5; m++) {
               rhs[k][j][i][m] = forcing[k][j][i][m];
            }
         }
      }
   }
   //---------------------------------------------------------------------
   // compute xi-direction fluxes
   //---------------------------------------------------------------------
   for(k = 1; k <= nz2; k++) {
      for(j = 1; j <= ny2; j++) {
         for(i = 1; i <= nx2; i++) {
            uijk = us[k][j][i];
            up1 = us[k][j][i + 1];
            um1 = us[k][j][i - 1];
            rhs[k][j][i][0] = rhs[k][j][i][0] + dx1tx1 * (u[k][j][i + 1][0] - 2.0 * u[k][j][i][0] + u[k][j][i - 1][0]) - tx2 * (u[k][j][i + 1][1] - u[k][j][i - 1][1]);
            rhs[k][j][i][1] = rhs[k][j][i][1] + dx2tx1 * (u[k][j][i + 1][1] - 2.0 * u[k][j][i][1] + u[k][j][i - 1][1]) + xxcon2 * con43 * (up1 - 2.0 * uijk + um1) - tx2 * (u[k][j][i + 1][1] * up1 - u[k][j][i - 1][1] * um1 + (u[k][j][i + 1][4] - square[k][j][i + 1] - u[k][j][i - 1][4] + square[k][j][i - 1]) * c2);
            rhs[k][j][i][2] = rhs[k][j][i][2] + dx3tx1 * (u[k][j][i + 1][2] - 2.0 * u[k][j][i][2] + u[k][j][i - 1][2]) + xxcon2 * (vs[k][j][i + 1] - 2.0 * vs[k][j][i] + vs[k][j][i - 1]) - tx2 * (u[k][j][i + 1][2] * up1 - u[k][j][i - 1][2] * um1);
            rhs[k][j][i][3] = rhs[k][j][i][3] + dx4tx1 * (u[k][j][i + 1][3] - 2.0 * u[k][j][i][3] + u[k][j][i - 1][3]) + xxcon2 * (ws[k][j][i + 1] - 2.0 * ws[k][j][i] + ws[k][j][i - 1]) - tx2 * (u[k][j][i + 1][3] * up1 - u[k][j][i - 1][3] * um1);
            rhs[k][j][i][4] = rhs[k][j][i][4] + dx5tx1 * (u[k][j][i + 1][4] - 2.0 * u[k][j][i][4] + u[k][j][i - 1][4]) + xxcon3 * (qs[k][j][i + 1] - 2.0 * qs[k][j][i] + qs[k][j][i - 1]) + xxcon4 * (up1 * up1 - 2.0 * uijk * uijk + um1 * um1) + xxcon5 * (u[k][j][i + 1][4] * rho_i[k][j][i + 1] - 2.0 * u[k][j][i][4] * rho_i[k][j][i] + u[k][j][i - 1][4] * rho_i[k][j][i - 1]) - tx2 * ((c1 * u[k][j][i + 1][4] - c2 * square[k][j][i + 1]) * up1 - (c1 * u[k][j][i - 1][4] - c2 * square[k][j][i - 1]) * um1);
         }
      }
      //---------------------------------------------------------------------
      // add fourth order xi-direction dissipation
      //---------------------------------------------------------------------
      for(j = 1; j <= ny2; j++) {
         i = 1;
         for(m = 0; m < 5; m++) {
            rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (5.0 * u[k][j][i][m] - 4.0 * u[k][j][i + 1][m] + u[k][j][i + 2][m]);
         }
         i = 2;
         for(m = 0; m < 5; m++) {
            rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (-4.0 * u[k][j][i - 1][m] + 6.0 * u[k][j][i][m] - 4.0 * u[k][j][i + 1][m] + u[k][j][i + 2][m]);
         }
      }
      for(j = 1; j <= ny2; j++) {
         for(i = 3; i <= nx2 - 2; i++) {
            for(m = 0; m < 5; m++) {
               rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (u[k][j][i - 2][m] - 4.0 * u[k][j][i - 1][m] + 6.0 * u[k][j][i][m] - 4.0 * u[k][j][i + 1][m] + u[k][j][i + 2][m]);
            }
         }
      }
      for(j = 1; j <= ny2; j++) {
         i = nx2 - 1;
         for(m = 0; m < 5; m++) {
            rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (u[k][j][i - 2][m] - 4.0 * u[k][j][i - 1][m] + 6.0 * u[k][j][i][m] - 4.0 * u[k][j][i + 1][m]);
         }
         i = nx2;
         for(m = 0; m < 5; m++) {
            rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (u[k][j][i - 2][m] - 4.0 * u[k][j][i - 1][m] + 5.0 * u[k][j][i][m]);
         }
      }
   }
   //---------------------------------------------------------------------
   // compute eta-direction fluxes
   //---------------------------------------------------------------------
   for(k = 1; k <= nz2; k++) {
      for(j = 1; j <= ny2; j++) {
         for(i = 1; i <= nx2; i++) {
            vijk = vs[k][j][i];
            vp1 = vs[k][j + 1][i];
            vm1 = vs[k][j - 1][i];
            rhs[k][j][i][0] = rhs[k][j][i][0] + dy1ty1 * (u[k][j + 1][i][0] - 2.0 * u[k][j][i][0] + u[k][j - 1][i][0]) - ty2 * (u[k][j + 1][i][2] - u[k][j - 1][i][2]);
            rhs[k][j][i][1] = rhs[k][j][i][1] + dy2ty1 * (u[k][j + 1][i][1] - 2.0 * u[k][j][i][1] + u[k][j - 1][i][1]) + yycon2 * (us[k][j + 1][i] - 2.0 * us[k][j][i] + us[k][j - 1][i]) - ty2 * (u[k][j + 1][i][1] * vp1 - u[k][j - 1][i][1] * vm1);
            rhs[k][j][i][2] = rhs[k][j][i][2] + dy3ty1 * (u[k][j + 1][i][2] - 2.0 * u[k][j][i][2] + u[k][j - 1][i][2]) + yycon2 * con43 * (vp1 - 2.0 * vijk + vm1) - ty2 * (u[k][j + 1][i][2] * vp1 - u[k][j - 1][i][2] * vm1 + (u[k][j + 1][i][4] - square[k][j + 1][i] - u[k][j - 1][i][4] + square[k][j - 1][i]) * c2);
            rhs[k][j][i][3] = rhs[k][j][i][3] + dy4ty1 * (u[k][j + 1][i][3] - 2.0 * u[k][j][i][3] + u[k][j - 1][i][3]) + yycon2 * (ws[k][j + 1][i] - 2.0 * ws[k][j][i] + ws[k][j - 1][i]) - ty2 * (u[k][j + 1][i][3] * vp1 - u[k][j - 1][i][3] * vm1);
            rhs[k][j][i][4] = rhs[k][j][i][4] + dy5ty1 * (u[k][j + 1][i][4] - 2.0 * u[k][j][i][4] + u[k][j - 1][i][4]) + yycon3 * (qs[k][j + 1][i] - 2.0 * qs[k][j][i] + qs[k][j - 1][i]) + yycon4 * (vp1 * vp1 - 2.0 * vijk * vijk + vm1 * vm1) + yycon5 * (u[k][j + 1][i][4] * rho_i[k][j + 1][i] - 2.0 * u[k][j][i][4] * rho_i[k][j][i] + u[k][j - 1][i][4] * rho_i[k][j - 1][i]) - ty2 * ((c1 * u[k][j + 1][i][4] - c2 * square[k][j + 1][i]) * vp1 - (c1 * u[k][j - 1][i][4] - c2 * square[k][j - 1][i]) * vm1);
         }
      }
      //---------------------------------------------------------------------
      // add fourth order eta-direction dissipation
      //---------------------------------------------------------------------
      j = 1;
      for(i = 1; i <= nx2; i++) {
         for(m = 0; m < 5; m++) {
            rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (5.0 * u[k][j][i][m] - 4.0 * u[k][j + 1][i][m] + u[k][j + 2][i][m]);
         }
      }
      j = 2;
      for(i = 1; i <= nx2; i++) {
         for(m = 0; m < 5; m++) {
            rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (-4.0 * u[k][j - 1][i][m] + 6.0 * u[k][j][i][m] - 4.0 * u[k][j + 1][i][m] + u[k][j + 2][i][m]);
         }
      }
      for(j = 3; j <= ny2 - 2; j++) {
         for(i = 1; i <= nx2; i++) {
            for(m = 0; m < 5; m++) {
               rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (u[k][j - 2][i][m] - 4.0 * u[k][j - 1][i][m] + 6.0 * u[k][j][i][m] - 4.0 * u[k][j + 1][i][m] + u[k][j + 2][i][m]);
            }
         }
      }
      j = ny2 - 1;
      for(i = 1; i <= nx2; i++) {
         for(m = 0; m < 5; m++) {
            rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (u[k][j - 2][i][m] - 4.0 * u[k][j - 1][i][m] + 6.0 * u[k][j][i][m] - 4.0 * u[k][j + 1][i][m]);
         }
      }
      j = ny2;
      for(i = 1; i <= nx2; i++) {
         for(m = 0; m < 5; m++) {
            rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (u[k][j - 2][i][m] - 4.0 * u[k][j - 1][i][m] + 5.0 * u[k][j][i][m]);
         }
      }
   }
   //---------------------------------------------------------------------
   // compute zeta-direction fluxes
   //---------------------------------------------------------------------
   for(k = 1; k <= nz2; k++) {
      for(j = 1; j <= ny2; j++) {
         for(i = 1; i <= nx2; i++) {
            wijk = ws[k][j][i];
            wp1 = ws[k + 1][j][i];
            wm1 = ws[k - 1][j][i];
            rhs[k][j][i][0] = rhs[k][j][i][0] + dz1tz1 * (u[k + 1][j][i][0] - 2.0 * u[k][j][i][0] + u[k - 1][j][i][0]) - tz2 * (u[k + 1][j][i][3] - u[k - 1][j][i][3]);
            rhs[k][j][i][1] = rhs[k][j][i][1] + dz2tz1 * (u[k + 1][j][i][1] - 2.0 * u[k][j][i][1] + u[k - 1][j][i][1]) + zzcon2 * (us[k + 1][j][i] - 2.0 * us[k][j][i] + us[k - 1][j][i]) - tz2 * (u[k + 1][j][i][1] * wp1 - u[k - 1][j][i][1] * wm1);
            rhs[k][j][i][2] = rhs[k][j][i][2] + dz3tz1 * (u[k + 1][j][i][2] - 2.0 * u[k][j][i][2] + u[k - 1][j][i][2]) + zzcon2 * (vs[k + 1][j][i] - 2.0 * vs[k][j][i] + vs[k - 1][j][i]) - tz2 * (u[k + 1][j][i][2] * wp1 - u[k - 1][j][i][2] * wm1);
            rhs[k][j][i][3] = rhs[k][j][i][3] + dz4tz1 * (u[k + 1][j][i][3] - 2.0 * u[k][j][i][3] + u[k - 1][j][i][3]) + zzcon2 * con43 * (wp1 - 2.0 * wijk + wm1) - tz2 * (u[k + 1][j][i][3] * wp1 - u[k - 1][j][i][3] * wm1 + (u[k + 1][j][i][4] - square[k + 1][j][i] - u[k - 1][j][i][4] + square[k - 1][j][i]) * c2);
            rhs[k][j][i][4] = rhs[k][j][i][4] + dz5tz1 * (u[k + 1][j][i][4] - 2.0 * u[k][j][i][4] + u[k - 1][j][i][4]) + zzcon3 * (qs[k + 1][j][i] - 2.0 * qs[k][j][i] + qs[k - 1][j][i]) + zzcon4 * (wp1 * wp1 - 2.0 * wijk * wijk + wm1 * wm1) + zzcon5 * (u[k + 1][j][i][4] * rho_i[k + 1][j][i] - 2.0 * u[k][j][i][4] * rho_i[k][j][i] + u[k - 1][j][i][4] * rho_i[k - 1][j][i]) - tz2 * ((c1 * u[k + 1][j][i][4] - c2 * square[k + 1][j][i]) * wp1 - (c1 * u[k - 1][j][i][4] - c2 * square[k - 1][j][i]) * wm1);
         }
      }
   }
   //---------------------------------------------------------------------
   // add fourth order zeta-direction dissipation
   //---------------------------------------------------------------------
   k = 1;
   for(j = 1; j <= ny2; j++) {
      for(i = 1; i <= nx2; i++) {
         for(m = 0; m < 5; m++) {
            rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (5.0 * u[k][j][i][m] - 4.0 * u[k + 1][j][i][m] + u[k + 2][j][i][m]);
         }
      }
   }
   k = 2;
   for(j = 1; j <= ny2; j++) {
      for(i = 1; i <= nx2; i++) {
         for(m = 0; m < 5; m++) {
            rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (-4.0 * u[k - 1][j][i][m] + 6.0 * u[k][j][i][m] - 4.0 * u[k + 1][j][i][m] + u[k + 2][j][i][m]);
         }
      }
   }
   for(k = 3; k <= nz2 - 2; k++) {
      for(j = 1; j <= ny2; j++) {
         for(i = 1; i <= nx2; i++) {
            for(m = 0; m < 5; m++) {
               rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (u[k - 2][j][i][m] - 4.0 * u[k - 1][j][i][m] + 6.0 * u[k][j][i][m] - 4.0 * u[k + 1][j][i][m] + u[k + 2][j][i][m]);
            }
         }
      }
   }
   k = nz2 - 1;
   for(j = 1; j <= ny2; j++) {
      for(i = 1; i <= nx2; i++) {
         for(m = 0; m < 5; m++) {
            rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (u[k - 2][j][i][m] - 4.0 * u[k - 1][j][i][m] + 6.0 * u[k][j][i][m] - 4.0 * u[k + 1][j][i][m]);
         }
      }
   }
   k = nz2;
   for(j = 1; j <= ny2; j++) {
      for(i = 1; i <= nx2; i++) {
         for(m = 0; m < 5; m++) {
            rhs[k][j][i][m] = rhs[k][j][i][m] - dssp * (u[k - 2][j][i][m] - 4.0 * u[k - 1][j][i][m] + 5.0 * u[k][j][i][m]);
         }
      }
   }
   for(k = 1; k <= nz2; k++) {
      for(j = 1; j <= ny2; j++) {
         for(i = 1; i <= nx2; i++) {
            for(m = 0; m < 5; m++) {
               rhs[k][j][i][m] = rhs[k][j][i][m] * dt;
            }
         }
      }
   }
}

//---------------------------------------------------------------------
// block-diagonal matrix-vector multiplication
//---------------------------------------------------------------------
void txinvr() {
   int i, j, k;
   double t1, t2, t3, ac, ru1, uu, vv, ww, r1, r2, r3, r4, r5, ac2inv;
   for(k = 1; k <= nz2; k++) {
      for(j = 1; j <= ny2; j++) {
         for(i = 1; i <= nx2; i++) {
            ru1 = rho_i[k][j][i];
            uu = us[k][j][i];
            vv = vs[k][j][i];
            ww = ws[k][j][i];
            ac = speed[k][j][i];
            ac2inv = ac * ac;
            r1 = rhs[k][j][i][0];
            r2 = rhs[k][j][i][1];
            r3 = rhs[k][j][i][2];
            r4 = rhs[k][j][i][3];
            r5 = rhs[k][j][i][4];
            t1 = c2 / ac2inv * (qs[k][j][i] * r1 - uu * r2 - vv * r3 - ww * r4 + r5);
            t2 = bt * ru1 * (uu * r1 - r2);
            t3 = (bt * ru1 * ac) * t1;
            rhs[k][j][i][0] = r1 - t1;
            rhs[k][j][i][1] = -ru1 * (ww * r1 - r4);
            rhs[k][j][i][2] = ru1 * (vv * r1 - r3);
            rhs[k][j][i][3] = -t2 + t3;
            rhs[k][j][i][4] = t2 + t3;
         }
      }
   }
}

//---------------------------------------------------------------------
// block-diagonal matrix-vector multiplication
//---------------------------------------------------------------------
void tzetar() {
   int i, j, k;
   double t1, t2, t3, ac, xvel, yvel, zvel, r1, r2, r3, r4, r5;
   double btuz, ac2u, uzik1;
   for(k = 1; k <= nz2; k++) {
      for(j = 1; j <= ny2; j++) {
         for(i = 1; i <= nx2; i++) {
            xvel = us[k][j][i];
            yvel = vs[k][j][i];
            zvel = ws[k][j][i];
            ac = speed[k][j][i];
            ac2u = ac * ac;
            r1 = rhs[k][j][i][0];
            r2 = rhs[k][j][i][1];
            r3 = rhs[k][j][i][2];
            r4 = rhs[k][j][i][3];
            r5 = rhs[k][j][i][4];
            uzik1 = u[k][j][i][0];
            btuz = bt * uzik1;
            t1 = btuz / ac * (r4 + r5);
            t2 = r3 + t1;
            t3 = btuz * (r4 - r5);
            rhs[k][j][i][0] = t2;
            rhs[k][j][i][1] = -uzik1 * r2 + xvel * t2;
            rhs[k][j][i][2] = uzik1 * r1 + yvel * t2;
            rhs[k][j][i][3] = zvel * t2 + t3;
            rhs[k][j][i][4] = uzik1 * (-xvel * r2 + yvel * r1) + qs[k][j][i] * t2 + c2iv * ac2u * t1 + zvel * t3;
         }
      }
   }
}

//---------------------------------------------------------------------
// verification routine
//---------------------------------------------------------------------
void verify(int no_time_steps, char *Class, int *verified) {
   double xcrref[5];
   double xceref[5];
   double xcrdif[5];
   double xcedif[5];
   double epsilon;
   double xce[5];
   double xcr[5];
   double dtref = 0.0;
   int m;
   //---------------------------------------------------------------------
   // tolerance level
   //---------------------------------------------------------------------
   epsilon = 1.0e-08;
   //---------------------------------------------------------------------
   // compute the error norm and the residual norm, and exit if not printing
   //---------------------------------------------------------------------
   clava_dcg_global[ 23 ]++;
   error_norm(xce);
   clava_dcg_global[ 24 ]++;
   compute_rhs();
   clava_dcg_global[ 25 ]++;
   rhs_norm(xcr);
   for(m = 0; m < 5; m++) {
      xcr[m] = xcr[m] / dt;
   }
   *Class = 'U';
   *verified = 1;
   for(m = 0; m < 5; m++) {
      xcrref[m] = 1.0;
      xceref[m] = 1.0;
   }
   //---------------------------------------------------------------------
   // reference data for 12X12X12 grids after 100 time steps,
   // with DT = 1.50e-02
   //---------------------------------------------------------------------
   if((grid_points[0] == 12) && (grid_points[1] == 12) && (grid_points[2] == 12) && (no_time_steps == 100)) {
      *Class = 'S';
      dtref = 1.5e-2;
      //---------------------------------------------------------------------
      // Reference values of RMS-norms of residual.
      //---------------------------------------------------------------------
      xcrref[0] = 2.7470315451339479e-02;
      xcrref[1] = 1.0360746705285417e-02;
      xcrref[2] = 1.6235745065095532e-02;
      xcrref[3] = 1.5840557224455615e-02;
      xcrref[4] = 3.4849040609362460e-02;
      //---------------------------------------------------------------------
      // Reference values of RMS-norms of solution error.
      //---------------------------------------------------------------------
      xceref[0] = 2.7289258557377227e-05;
      xceref[1] = 1.0364446640837285e-05;
      xceref[2] = 1.6154798287166471e-05;
      xceref[3] = 1.5750704994480102e-05;
      xceref[4] = 3.4177666183390531e-05;
      //---------------------------------------------------------------------
      // reference data for 36X36X36 grids after 400 time steps,
      // with DT = 1.5e-03
      //---------------------------------------------------------------------
   }
   else if((grid_points[0] == 36) && (grid_points[1] == 36) && (grid_points[2] == 36) && (no_time_steps == 400)) {
      *Class = 'W';
      dtref = 1.5e-3;
      //---------------------------------------------------------------------
      // Reference values of RMS-norms of residual.
      //---------------------------------------------------------------------
      xcrref[0] = 0.1893253733584e-02;
      xcrref[1] = 0.1717075447775e-03;
      xcrref[2] = 0.2778153350936e-03;
      xcrref[3] = 0.2887475409984e-03;
      xcrref[4] = 0.3143611161242e-02;
      //---------------------------------------------------------------------
      // Reference values of RMS-norms of solution error.
      //---------------------------------------------------------------------
      xceref[0] = 0.7542088599534e-04;
      xceref[1] = 0.6512852253086e-05;
      xceref[2] = 0.1049092285688e-04;
      xceref[3] = 0.1128838671535e-04;
      xceref[4] = 0.1212845639773e-03;
      //---------------------------------------------------------------------
      // reference data for 64X64X64 grids after 400 time steps,
      // with DT = 1.5e-03
      //---------------------------------------------------------------------
   }
   else if((grid_points[0] == 64) && (grid_points[1] == 64) && (grid_points[2] == 64) && (no_time_steps == 400)) {
      *Class = 'A';
      dtref = 1.5e-3;
      //---------------------------------------------------------------------
      // Reference values of RMS-norms of residual.
      //---------------------------------------------------------------------
      xcrref[0] = 2.4799822399300195;
      xcrref[1] = 1.1276337964368832;
      xcrref[2] = 1.5028977888770491;
      xcrref[3] = 1.4217816211695179;
      xcrref[4] = 2.1292113035138280;
      //---------------------------------------------------------------------
      // Reference values of RMS-norms of solution error.
      //---------------------------------------------------------------------
      xceref[0] = 1.0900140297820550e-04;
      xceref[1] = 3.7343951769282091e-05;
      xceref[2] = 5.0092785406541633e-05;
      xceref[3] = 4.7671093939528255e-05;
      xceref[4] = 1.3621613399213001e-04;
      //---------------------------------------------------------------------
      // reference data for 102X102X102 grids after 400 time steps,
      // with DT = 1.0e-03
      //---------------------------------------------------------------------
   }
   else if((grid_points[0] == 102) && (grid_points[1] == 102) && (grid_points[2] == 102) && (no_time_steps == 400)) {
      *Class = 'B';
      dtref = 1.0e-3;
      //---------------------------------------------------------------------
      // Reference values of RMS-norms of residual.
      //---------------------------------------------------------------------
      xcrref[0] = 0.6903293579998e+02;
      xcrref[1] = 0.3095134488084e+02;
      xcrref[2] = 0.4103336647017e+02;
      xcrref[3] = 0.3864769009604e+02;
      xcrref[4] = 0.5643482272596e+02;
      //---------------------------------------------------------------------
      // Reference values of RMS-norms of solution error.
      //---------------------------------------------------------------------
      xceref[0] = 0.9810006190188e-02;
      xceref[1] = 0.1022827905670e-02;
      xceref[2] = 0.1720597911692e-02;
      xceref[3] = 0.1694479428231e-02;
      xceref[4] = 0.1847456263981e-01;
      //---------------------------------------------------------------------
      // reference data for 162X162X162 grids after 400 time steps,
      // with DT = 0.67e-03
      //---------------------------------------------------------------------
   }
   else if((grid_points[0] == 162) && (grid_points[1] == 162) && (grid_points[2] == 162) && (no_time_steps == 400)) {
      *Class = 'C';
      dtref = 0.67e-3;
      //---------------------------------------------------------------------
      // Reference values of RMS-norms of residual.
      //---------------------------------------------------------------------
      xcrref[0] = 0.5881691581829e+03;
      xcrref[1] = 0.2454417603569e+03;
      xcrref[2] = 0.3293829191851e+03;
      xcrref[3] = 0.3081924971891e+03;
      xcrref[4] = 0.4597223799176e+03;
      //---------------------------------------------------------------------
      // Reference values of RMS-norms of solution error.
      //---------------------------------------------------------------------
      xceref[0] = 0.2598120500183e+00;
      xceref[1] = 0.2590888922315e-01;
      xceref[2] = 0.5132886416320e-01;
      xceref[3] = 0.4806073419454e-01;
      xceref[4] = 0.5483377491301e+00;
      //---------------------------------------------------------------------
      // reference data for 408X408X408 grids after 500 time steps,
      // with DT = 0.3e-03
      //---------------------------------------------------------------------
   }
   else if((grid_points[0] == 408) && (grid_points[1] == 408) && (grid_points[2] == 408) && (no_time_steps == 500)) {
      *Class = 'D';
      dtref = 0.30e-3;
      //---------------------------------------------------------------------
      // Reference values of RMS-norms of residual.
      //---------------------------------------------------------------------
      xcrref[0] = 0.1044696216887e+05;
      xcrref[1] = 0.3204427762578e+04;
      xcrref[2] = 0.4648680733032e+04;
      xcrref[3] = 0.4238923283697e+04;
      xcrref[4] = 0.7588412036136e+04;
      //---------------------------------------------------------------------
      // Reference values of RMS-norms of solution error.
      //---------------------------------------------------------------------
      xceref[0] = 0.5089471423669e+01;
      xceref[1] = 0.5323514855894e+00;
      xceref[2] = 0.1187051008971e+01;
      xceref[3] = 0.1083734951938e+01;
      xceref[4] = 0.1164108338568e+02;
      //---------------------------------------------------------------------
      // reference data for 1020X1020X1020 grids after 500 time steps,
      // with DT = 0.1e-03
      //---------------------------------------------------------------------
   }
   else if((grid_points[0] == 1020) && (grid_points[1] == 1020) && (grid_points[2] == 1020) && (no_time_steps == 500)) {
      *Class = 'E';
      dtref = 0.10e-3;
      //---------------------------------------------------------------------
      // Reference values of RMS-norms of residual.
      //---------------------------------------------------------------------
      xcrref[0] = 0.6255387422609e+05;
      xcrref[1] = 0.1495317020012e+05;
      xcrref[2] = 0.2347595750586e+05;
      xcrref[3] = 0.2091099783534e+05;
      xcrref[4] = 0.4770412841218e+05;
      //---------------------------------------------------------------------
      // Reference values of RMS-norms of solution error.
      //---------------------------------------------------------------------
      xceref[0] = 0.6742735164909e+02;
      xceref[1] = 0.5390656036938e+01;
      xceref[2] = 0.1680647196477e+02;
      xceref[3] = 0.1536963126457e+02;
      xceref[4] = 0.1575330146156e+03;
   }
   else {
      *verified = 0;
   }
   //---------------------------------------------------------------------
   // verification test for residuals if gridsize is one of
   // the defined grid sizes above (class .ne. 'U')
   //---------------------------------------------------------------------
   //---------------------------------------------------------------------
   // Compute the difference of solution values and the known reference values.
   //---------------------------------------------------------------------
   for(m = 0; m < 5; m++) {
      clava_dcg_global[ 26 ]++;
      xcrdif[m] = fabs((xcr[m] - xcrref[m]) / xcrref[m]);
      clava_dcg_global[ 26 ]++;
      xcedif[m] = fabs((xce[m] - xceref[m]) / xceref[m]);
   }
   //---------------------------------------------------------------------
   // Output the comparison of computed results to known cases.
   //---------------------------------------------------------------------
   if(*Class != 'U') {
      clava_dcg_global[ 27 ]++;
      printf(" Verification being performed for class %c\n", *Class);
      clava_dcg_global[ 27 ]++;
      printf(" accuracy setting for epsilon = %20.13E\n", epsilon);
      clava_dcg_global[ 26 ]++;
      *verified = (fabs(dt - dtref) <= epsilon);
      if(!(*verified)) {
         *Class = 'U';
         clava_dcg_global[ 27 ]++;
         printf(" DT does not match the reference value of %15.8E\n", dtref);
      }
   }
   else {
      clava_dcg_global[ 27 ]++;
      printf(" Unknown class\n");
   }
   if(*Class != 'U') {
      clava_dcg_global[ 27 ]++;
      printf(" Comparison of RMS-norms of residual\n");
   }
   else {
      clava_dcg_global[ 27 ]++;
      printf(" RMS-norms of residual\n");
   }
   for(m = 0; m < 5; m++) {
      if(*Class == 'U') {
         clava_dcg_global[ 27 ]++;
         printf("          %2d%20.13E\n", m + 1, xcr[m]);
      }
      else if(xcrdif[m] <= epsilon) {
         clava_dcg_global[ 27 ]++;
         printf("          %2d%20.13E%20.13E%20.13E\n", m + 1, xcr[m], xcrref[m], xcrdif[m]);
      }
      else {
         *verified = 0;
         clava_dcg_global[ 27 ]++;
         printf(" FAILURE: %2d%20.13E%20.13E%20.13E\n", m + 1, xcr[m], xcrref[m], xcrdif[m]);
      }
   }
   if(*Class != 'U') {
      clava_dcg_global[ 27 ]++;
      printf(" Comparison of RMS-norms of solution error\n");
   }
   else {
      clava_dcg_global[ 27 ]++;
      printf(" RMS-norms of solution error\n");
   }
   for(m = 0; m < 5; m++) {
      if(*Class == 'U') {
         clava_dcg_global[ 27 ]++;
         printf("          %2d%20.13E\n", m + 1, xce[m]);
      }
      else if(xcedif[m] <= epsilon) {
         clava_dcg_global[ 27 ]++;
         printf("          %2d%20.13E%20.13E%20.13E\n", m + 1, xce[m], xceref[m], xcedif[m]);
      }
      else {
         *verified = 0;
         clava_dcg_global[ 27 ]++;
         printf(" FAILURE: %2d%20.13E%20.13E%20.13E\n", m + 1, xce[m], xceref[m], xcedif[m]);
      }
   }
   if(*Class == 'U') {
      clava_dcg_global[ 27 ]++;
      printf(" No reference values provided\n");
      clava_dcg_global[ 27 ]++;
      printf(" No verification performed\n");
   }
   else if(*verified) {
      clava_dcg_global[ 27 ]++;
      printf(" Verification Successful\n");
   }
   else {
      clava_dcg_global[ 27 ]++;
      printf(" Verification failed\n");
   }
}

//---------------------------------------------------------------------
// this function performs the solution of the approximate factorization
// step in the x-direction for all five matrix components
// simultaneously. The Thomas algorithm is employed to solve the
// systems for the x-lines. Boundary conditions are non-periodic
//---------------------------------------------------------------------
void x_solve() {
   int i, j, k, i1, i2, m;
   double ru1, fac1, fac2;
   double rhon[36];
   double cv[36];
   double lhs[37][37][5];
   double lhsp[37][37][5];
   double lhsm[37][37][5];
   for(k = 1; k <= nz2; k++) {
      clava_dcg_global[ 28 ]++;
      lhsinit(nx2 + 1, ny2, lhs, lhsp, lhsm);
      //---------------------------------------------------------------------
      // Computes the left hand side for the three x-factors
      //---------------------------------------------------------------------
      //---------------------------------------------------------------------
      // first fill the lhs for the u-eigenvalue
      //---------------------------------------------------------------------
      for(j = 1; j <= ny2; j++) {
         for(i = 0; i <= grid_points[0] - 1; i++) {
            ru1 = c3c4 * rho_i[k][j][i];
            cv[i] = us[k][j][i];
            rhon[i] = ((((dx2 + con43 * ru1) > (dx5 + c1c5 * ru1) ? (dx2 + con43 * ru1) : (dx5 + c1c5 * ru1))) > (((dxmax + ru1) > (dx1) ? (dxmax + ru1) : (dx1))) ? (((dx2 + con43 * ru1) > (dx5 + c1c5 * ru1) ? (dx2 + con43 * ru1) : (dx5 + c1c5 * ru1))) : (((dxmax + ru1) > (dx1) ? (dxmax + ru1) : (dx1))));
         }
         for(i = 1; i <= nx2; i++) {
            lhs[j][i][0] = 0.0;
            lhs[j][i][1] = -dttx2 * cv[i - 1] - dttx1 * rhon[i - 1];
            lhs[j][i][2] = 1.0 + c2dttx1 * rhon[i];
            lhs[j][i][3] = dttx2 * cv[i + 1] - dttx1 * rhon[i + 1];
            lhs[j][i][4] = 0.0;
         }
      }
      //---------------------------------------------------------------------
      // add fourth order dissipation
      //---------------------------------------------------------------------
      for(j = 1; j <= ny2; j++) {
         i = 1;
         lhs[j][i][2] = lhs[j][i][2] + comz5;
         lhs[j][i][3] = lhs[j][i][3] - comz4;
         lhs[j][i][4] = lhs[j][i][4] + comz1;
         lhs[j][i + 1][1] = lhs[j][i + 1][1] - comz4;
         lhs[j][i + 1][2] = lhs[j][i + 1][2] + comz6;
         lhs[j][i + 1][3] = lhs[j][i + 1][3] - comz4;
         lhs[j][i + 1][4] = lhs[j][i + 1][4] + comz1;
      }
      for(j = 1; j <= ny2; j++) {
         for(i = 3; i <= grid_points[0] - 4; i++) {
            lhs[j][i][0] = lhs[j][i][0] + comz1;
            lhs[j][i][1] = lhs[j][i][1] - comz4;
            lhs[j][i][2] = lhs[j][i][2] + comz6;
            lhs[j][i][3] = lhs[j][i][3] - comz4;
            lhs[j][i][4] = lhs[j][i][4] + comz1;
         }
      }
      for(j = 1; j <= ny2; j++) {
         i = grid_points[0] - 3;
         lhs[j][i][0] = lhs[j][i][0] + comz1;
         lhs[j][i][1] = lhs[j][i][1] - comz4;
         lhs[j][i][2] = lhs[j][i][2] + comz6;
         lhs[j][i][3] = lhs[j][i][3] - comz4;
         lhs[j][i + 1][0] = lhs[j][i + 1][0] + comz1;
         lhs[j][i + 1][1] = lhs[j][i + 1][1] - comz4;
         lhs[j][i + 1][2] = lhs[j][i + 1][2] + comz5;
      }
      //---------------------------------------------------------------------
      // subsequently, fill the other factors (u+c), (u-c) by adding to
      // the first
      //---------------------------------------------------------------------
      for(j = 1; j <= ny2; j++) {
         for(i = 1; i <= nx2; i++) {
            lhsp[j][i][0] = lhs[j][i][0];
            lhsp[j][i][1] = lhs[j][i][1] - dttx2 * speed[k][j][i - 1];
            lhsp[j][i][2] = lhs[j][i][2];
            lhsp[j][i][3] = lhs[j][i][3] + dttx2 * speed[k][j][i + 1];
            lhsp[j][i][4] = lhs[j][i][4];
            lhsm[j][i][0] = lhs[j][i][0];
            lhsm[j][i][1] = lhs[j][i][1] + dttx2 * speed[k][j][i - 1];
            lhsm[j][i][2] = lhs[j][i][2];
            lhsm[j][i][3] = lhs[j][i][3] - dttx2 * speed[k][j][i + 1];
            lhsm[j][i][4] = lhs[j][i][4];
         }
      }
      //---------------------------------------------------------------------
      // FORWARD ELIMINATION
      //---------------------------------------------------------------------
      //---------------------------------------------------------------------
      // perform the Thomas algorithm; first, FORWARD ELIMINATION
      //---------------------------------------------------------------------
      for(j = 1; j <= ny2; j++) {
         for(i = 0; i <= grid_points[0] - 3; i++) {
            i1 = i + 1;
            i2 = i + 2;
            fac1 = 1.0 / lhs[j][i][2];
            lhs[j][i][3] = fac1 * lhs[j][i][3];
            lhs[j][i][4] = fac1 * lhs[j][i][4];
            for(m = 0; m < 3; m++) {
               rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
            }
            lhs[j][i1][2] = lhs[j][i1][2] - lhs[j][i1][1] * lhs[j][i][3];
            lhs[j][i1][3] = lhs[j][i1][3] - lhs[j][i1][1] * lhs[j][i][4];
            for(m = 0; m < 3; m++) {
               rhs[k][j][i1][m] = rhs[k][j][i1][m] - lhs[j][i1][1] * rhs[k][j][i][m];
            }
            lhs[j][i2][1] = lhs[j][i2][1] - lhs[j][i2][0] * lhs[j][i][3];
            lhs[j][i2][2] = lhs[j][i2][2] - lhs[j][i2][0] * lhs[j][i][4];
            for(m = 0; m < 3; m++) {
               rhs[k][j][i2][m] = rhs[k][j][i2][m] - lhs[j][i2][0] * rhs[k][j][i][m];
            }
         }
      }
      //---------------------------------------------------------------------
      // The last two rows in this grid block are a bit different,
      // since they for (not have two more rows available for the
      // elimination of off-diagonal entries
      //---------------------------------------------------------------------
      for(j = 1; j <= ny2; j++) {
         i = grid_points[0] - 2;
         i1 = grid_points[0] - 1;
         fac1 = 1.0 / lhs[j][i][2];
         lhs[j][i][3] = fac1 * lhs[j][i][3];
         lhs[j][i][4] = fac1 * lhs[j][i][4];
         for(m = 0; m < 3; m++) {
            rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
         }
         lhs[j][i1][2] = lhs[j][i1][2] - lhs[j][i1][1] * lhs[j][i][3];
         lhs[j][i1][3] = lhs[j][i1][3] - lhs[j][i1][1] * lhs[j][i][4];
         for(m = 0; m < 3; m++) {
            rhs[k][j][i1][m] = rhs[k][j][i1][m] - lhs[j][i1][1] * rhs[k][j][i][m];
         }
         //---------------------------------------------------------------------
         // scale the last row immediately
         //---------------------------------------------------------------------
         fac2 = 1.0 / lhs[j][i1][2];
         for(m = 0; m < 3; m++) {
            rhs[k][j][i1][m] = fac2 * rhs[k][j][i1][m];
         }
      }
      //---------------------------------------------------------------------
      // for (the u+c and the u-c factors
      //---------------------------------------------------------------------
      for(j = 1; j <= ny2; j++) {
         for(i = 0; i <= grid_points[0] - 3; i++) {
            i1 = i + 1;
            i2 = i + 2;
            m = 3;
            fac1 = 1.0 / lhsp[j][i][2];
            lhsp[j][i][3] = fac1 * lhsp[j][i][3];
            lhsp[j][i][4] = fac1 * lhsp[j][i][4];
            rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
            lhsp[j][i1][2] = lhsp[j][i1][2] - lhsp[j][i1][1] * lhsp[j][i][3];
            lhsp[j][i1][3] = lhsp[j][i1][3] - lhsp[j][i1][1] * lhsp[j][i][4];
            rhs[k][j][i1][m] = rhs[k][j][i1][m] - lhsp[j][i1][1] * rhs[k][j][i][m];
            lhsp[j][i2][1] = lhsp[j][i2][1] - lhsp[j][i2][0] * lhsp[j][i][3];
            lhsp[j][i2][2] = lhsp[j][i2][2] - lhsp[j][i2][0] * lhsp[j][i][4];
            rhs[k][j][i2][m] = rhs[k][j][i2][m] - lhsp[j][i2][0] * rhs[k][j][i][m];
            m = 4;
            fac1 = 1.0 / lhsm[j][i][2];
            lhsm[j][i][3] = fac1 * lhsm[j][i][3];
            lhsm[j][i][4] = fac1 * lhsm[j][i][4];
            rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
            lhsm[j][i1][2] = lhsm[j][i1][2] - lhsm[j][i1][1] * lhsm[j][i][3];
            lhsm[j][i1][3] = lhsm[j][i1][3] - lhsm[j][i1][1] * lhsm[j][i][4];
            rhs[k][j][i1][m] = rhs[k][j][i1][m] - lhsm[j][i1][1] * rhs[k][j][i][m];
            lhsm[j][i2][1] = lhsm[j][i2][1] - lhsm[j][i2][0] * lhsm[j][i][3];
            lhsm[j][i2][2] = lhsm[j][i2][2] - lhsm[j][i2][0] * lhsm[j][i][4];
            rhs[k][j][i2][m] = rhs[k][j][i2][m] - lhsm[j][i2][0] * rhs[k][j][i][m];
         }
      }
      //---------------------------------------------------------------------
      // And again the last two rows separately
      //---------------------------------------------------------------------
      for(j = 1; j <= ny2; j++) {
         i = grid_points[0] - 2;
         i1 = grid_points[0] - 1;
         m = 3;
         fac1 = 1.0 / lhsp[j][i][2];
         lhsp[j][i][3] = fac1 * lhsp[j][i][3];
         lhsp[j][i][4] = fac1 * lhsp[j][i][4];
         rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
         lhsp[j][i1][2] = lhsp[j][i1][2] - lhsp[j][i1][1] * lhsp[j][i][3];
         lhsp[j][i1][3] = lhsp[j][i1][3] - lhsp[j][i1][1] * lhsp[j][i][4];
         rhs[k][j][i1][m] = rhs[k][j][i1][m] - lhsp[j][i1][1] * rhs[k][j][i][m];
         m = 4;
         fac1 = 1.0 / lhsm[j][i][2];
         lhsm[j][i][3] = fac1 * lhsm[j][i][3];
         lhsm[j][i][4] = fac1 * lhsm[j][i][4];
         rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
         lhsm[j][i1][2] = lhsm[j][i1][2] - lhsm[j][i1][1] * lhsm[j][i][3];
         lhsm[j][i1][3] = lhsm[j][i1][3] - lhsm[j][i1][1] * lhsm[j][i][4];
         rhs[k][j][i1][m] = rhs[k][j][i1][m] - lhsm[j][i1][1] * rhs[k][j][i][m];
         //---------------------------------------------------------------------
         // Scale the last row immediately
         //---------------------------------------------------------------------
         rhs[k][j][i1][3] = rhs[k][j][i1][3] / lhsp[j][i1][2];
         rhs[k][j][i1][4] = rhs[k][j][i1][4] / lhsm[j][i1][2];
      }
      //---------------------------------------------------------------------
      // BACKSUBSTITUTION
      //---------------------------------------------------------------------
      for(j = 1; j <= ny2; j++) {
         i = grid_points[0] - 2;
         i1 = grid_points[0] - 1;
         for(m = 0; m < 3; m++) {
            rhs[k][j][i][m] = rhs[k][j][i][m] - lhs[j][i][3] * rhs[k][j][i1][m];
         }
         rhs[k][j][i][3] = rhs[k][j][i][3] - lhsp[j][i][3] * rhs[k][j][i1][3];
         rhs[k][j][i][4] = rhs[k][j][i][4] - lhsm[j][i][3] * rhs[k][j][i1][4];
      }
      //---------------------------------------------------------------------
      // The first three factors
      //---------------------------------------------------------------------
      for(j = 1; j <= ny2; j++) {
         for(i = grid_points[0] - 3; i >= 0; i--) {
            i1 = i + 1;
            i2 = i + 2;
            for(m = 0; m < 3; m++) {
               rhs[k][j][i][m] = rhs[k][j][i][m] - lhs[j][i][3] * rhs[k][j][i1][m] - lhs[j][i][4] * rhs[k][j][i2][m];
            }
            //-------------------------------------------------------------------
            // And the remaining two
            //-------------------------------------------------------------------
            rhs[k][j][i][3] = rhs[k][j][i][3] - lhsp[j][i][3] * rhs[k][j][i1][3] - lhsp[j][i][4] * rhs[k][j][i2][3];
            rhs[k][j][i][4] = rhs[k][j][i][4] - lhsm[j][i][3] * rhs[k][j][i1][4] - lhsm[j][i][4] * rhs[k][j][i2][4];
         }
      }
   }
   //---------------------------------------------------------------------
   // Do the block-diagonal inversion
   //---------------------------------------------------------------------
   clava_dcg_global[ 29 ]++;
   ninvr();
}

//---------------------------------------------------------------------
// this function performs the solution of the approximate factorization
// step in the y-direction for all five matrix components
// simultaneously. The Thomas algorithm is employed to solve the
// systems for the y-lines. Boundary conditions are non-periodic
//---------------------------------------------------------------------
void y_solve() {
   int i, j, k, j1, j2, m;
   double ru1, fac1, fac2;
   double rhoq[36];
   double cv[36];
   double lhs[37][37][5];
   double lhsp[37][37][5];
   double lhsm[37][37][5];
   for(k = 1; k <= grid_points[2] - 2; k++) {
      clava_dcg_global[ 30 ]++;
      lhsinitj(ny2 + 1, nx2, lhs, lhsp, lhsm);
      //---------------------------------------------------------------------
      // Computes the left hand side for the three y-factors
      //---------------------------------------------------------------------
      //---------------------------------------------------------------------
      // first fill the lhs for the u-eigenvalue
      //---------------------------------------------------------------------
      for(i = 1; i <= grid_points[0] - 2; i++) {
         for(j = 0; j <= grid_points[1] - 1; j++) {
            ru1 = c3c4 * rho_i[k][j][i];
            cv[j] = vs[k][j][i];
            rhoq[j] = ((((dy3 + con43 * ru1) > (dy5 + c1c5 * ru1) ? (dy3 + con43 * ru1) : (dy5 + c1c5 * ru1))) > (((dymax + ru1) > (dy1) ? (dymax + ru1) : (dy1))) ? (((dy3 + con43 * ru1) > (dy5 + c1c5 * ru1) ? (dy3 + con43 * ru1) : (dy5 + c1c5 * ru1))) : (((dymax + ru1) > (dy1) ? (dymax + ru1) : (dy1))));
         }
         for(j = 1; j <= grid_points[1] - 2; j++) {
            lhs[j][i][0] = 0.0;
            lhs[j][i][1] = -dtty2 * cv[j - 1] - dtty1 * rhoq[j - 1];
            lhs[j][i][2] = 1.0 + c2dtty1 * rhoq[j];
            lhs[j][i][3] = dtty2 * cv[j + 1] - dtty1 * rhoq[j + 1];
            lhs[j][i][4] = 0.0;
         }
      }
      //---------------------------------------------------------------------
      // add fourth order dissipation
      //---------------------------------------------------------------------
      for(i = 1; i <= grid_points[0] - 2; i++) {
         j = 1;
         lhs[j][i][2] = lhs[j][i][2] + comz5;
         lhs[j][i][3] = lhs[j][i][3] - comz4;
         lhs[j][i][4] = lhs[j][i][4] + comz1;
         lhs[j + 1][i][1] = lhs[j + 1][i][1] - comz4;
         lhs[j + 1][i][2] = lhs[j + 1][i][2] + comz6;
         lhs[j + 1][i][3] = lhs[j + 1][i][3] - comz4;
         lhs[j + 1][i][4] = lhs[j + 1][i][4] + comz1;
      }
      for(j = 3; j <= grid_points[1] - 4; j++) {
         for(i = 1; i <= grid_points[0] - 2; i++) {
            lhs[j][i][0] = lhs[j][i][0] + comz1;
            lhs[j][i][1] = lhs[j][i][1] - comz4;
            lhs[j][i][2] = lhs[j][i][2] + comz6;
            lhs[j][i][3] = lhs[j][i][3] - comz4;
            lhs[j][i][4] = lhs[j][i][4] + comz1;
         }
      }
      for(i = 1; i <= grid_points[0] - 2; i++) {
         j = grid_points[1] - 3;
         lhs[j][i][0] = lhs[j][i][0] + comz1;
         lhs[j][i][1] = lhs[j][i][1] - comz4;
         lhs[j][i][2] = lhs[j][i][2] + comz6;
         lhs[j][i][3] = lhs[j][i][3] - comz4;
         lhs[j + 1][i][0] = lhs[j + 1][i][0] + comz1;
         lhs[j + 1][i][1] = lhs[j + 1][i][1] - comz4;
         lhs[j + 1][i][2] = lhs[j + 1][i][2] + comz5;
      }
      //---------------------------------------------------------------------
      // subsequently, for (the other two factors
      //---------------------------------------------------------------------
      for(j = 1; j <= grid_points[1] - 2; j++) {
         for(i = 1; i <= grid_points[0] - 2; i++) {
            lhsp[j][i][0] = lhs[j][i][0];
            lhsp[j][i][1] = lhs[j][i][1] - dtty2 * speed[k][j - 1][i];
            lhsp[j][i][2] = lhs[j][i][2];
            lhsp[j][i][3] = lhs[j][i][3] + dtty2 * speed[k][j + 1][i];
            lhsp[j][i][4] = lhs[j][i][4];
            lhsm[j][i][0] = lhs[j][i][0];
            lhsm[j][i][1] = lhs[j][i][1] + dtty2 * speed[k][j - 1][i];
            lhsm[j][i][2] = lhs[j][i][2];
            lhsm[j][i][3] = lhs[j][i][3] - dtty2 * speed[k][j + 1][i];
            lhsm[j][i][4] = lhs[j][i][4];
         }
      }
      //---------------------------------------------------------------------
      // FORWARD ELIMINATION
      //---------------------------------------------------------------------
      for(j = 0; j <= grid_points[1] - 3; j++) {
         j1 = j + 1;
         j2 = j + 2;
         for(i = 1; i <= grid_points[0] - 2; i++) {
            fac1 = 1.0 / lhs[j][i][2];
            lhs[j][i][3] = fac1 * lhs[j][i][3];
            lhs[j][i][4] = fac1 * lhs[j][i][4];
            for(m = 0; m < 3; m++) {
               rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
            }
            lhs[j1][i][2] = lhs[j1][i][2] - lhs[j1][i][1] * lhs[j][i][3];
            lhs[j1][i][3] = lhs[j1][i][3] - lhs[j1][i][1] * lhs[j][i][4];
            for(m = 0; m < 3; m++) {
               rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhs[j1][i][1] * rhs[k][j][i][m];
            }
            lhs[j2][i][1] = lhs[j2][i][1] - lhs[j2][i][0] * lhs[j][i][3];
            lhs[j2][i][2] = lhs[j2][i][2] - lhs[j2][i][0] * lhs[j][i][4];
            for(m = 0; m < 3; m++) {
               rhs[k][j2][i][m] = rhs[k][j2][i][m] - lhs[j2][i][0] * rhs[k][j][i][m];
            }
         }
      }
      //---------------------------------------------------------------------
      // The last two rows in this grid block are a bit different,
      // since they for (not have two more rows available for the
      // elimination of off-diagonal entries
      //---------------------------------------------------------------------
      j = grid_points[1] - 2;
      j1 = grid_points[1] - 1;
      for(i = 1; i <= grid_points[0] - 2; i++) {
         fac1 = 1.0 / lhs[j][i][2];
         lhs[j][i][3] = fac1 * lhs[j][i][3];
         lhs[j][i][4] = fac1 * lhs[j][i][4];
         for(m = 0; m < 3; m++) {
            rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
         }
         lhs[j1][i][2] = lhs[j1][i][2] - lhs[j1][i][1] * lhs[j][i][3];
         lhs[j1][i][3] = lhs[j1][i][3] - lhs[j1][i][1] * lhs[j][i][4];
         for(m = 0; m < 3; m++) {
            rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhs[j1][i][1] * rhs[k][j][i][m];
         }
         //---------------------------------------------------------------------
         // scale the last row immediately
         //---------------------------------------------------------------------
         fac2 = 1.0 / lhs[j1][i][2];
         for(m = 0; m < 3; m++) {
            rhs[k][j1][i][m] = fac2 * rhs[k][j1][i][m];
         }
      }
      //---------------------------------------------------------------------
      // for (the u+c and the u-c factors
      //---------------------------------------------------------------------
      for(j = 0; j <= grid_points[1] - 3; j++) {
         j1 = j + 1;
         j2 = j + 2;
         for(i = 1; i <= grid_points[0] - 2; i++) {
            m = 3;
            fac1 = 1.0 / lhsp[j][i][2];
            lhsp[j][i][3] = fac1 * lhsp[j][i][3];
            lhsp[j][i][4] = fac1 * lhsp[j][i][4];
            rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
            lhsp[j1][i][2] = lhsp[j1][i][2] - lhsp[j1][i][1] * lhsp[j][i][3];
            lhsp[j1][i][3] = lhsp[j1][i][3] - lhsp[j1][i][1] * lhsp[j][i][4];
            rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhsp[j1][i][1] * rhs[k][j][i][m];
            lhsp[j2][i][1] = lhsp[j2][i][1] - lhsp[j2][i][0] * lhsp[j][i][3];
            lhsp[j2][i][2] = lhsp[j2][i][2] - lhsp[j2][i][0] * lhsp[j][i][4];
            rhs[k][j2][i][m] = rhs[k][j2][i][m] - lhsp[j2][i][0] * rhs[k][j][i][m];
            m = 4;
            fac1 = 1.0 / lhsm[j][i][2];
            lhsm[j][i][3] = fac1 * lhsm[j][i][3];
            lhsm[j][i][4] = fac1 * lhsm[j][i][4];
            rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
            lhsm[j1][i][2] = lhsm[j1][i][2] - lhsm[j1][i][1] * lhsm[j][i][3];
            lhsm[j1][i][3] = lhsm[j1][i][3] - lhsm[j1][i][1] * lhsm[j][i][4];
            rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhsm[j1][i][1] * rhs[k][j][i][m];
            lhsm[j2][i][1] = lhsm[j2][i][1] - lhsm[j2][i][0] * lhsm[j][i][3];
            lhsm[j2][i][2] = lhsm[j2][i][2] - lhsm[j2][i][0] * lhsm[j][i][4];
            rhs[k][j2][i][m] = rhs[k][j2][i][m] - lhsm[j2][i][0] * rhs[k][j][i][m];
         }
      }
      //---------------------------------------------------------------------
      // And again the last two rows separately
      //---------------------------------------------------------------------
      j = grid_points[1] - 2;
      j1 = grid_points[1] - 1;
      for(i = 1; i <= grid_points[0] - 2; i++) {
         m = 3;
         fac1 = 1.0 / lhsp[j][i][2];
         lhsp[j][i][3] = fac1 * lhsp[j][i][3];
         lhsp[j][i][4] = fac1 * lhsp[j][i][4];
         rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
         lhsp[j1][i][2] = lhsp[j1][i][2] - lhsp[j1][i][1] * lhsp[j][i][3];
         lhsp[j1][i][3] = lhsp[j1][i][3] - lhsp[j1][i][1] * lhsp[j][i][4];
         rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhsp[j1][i][1] * rhs[k][j][i][m];
         m = 4;
         fac1 = 1.0 / lhsm[j][i][2];
         lhsm[j][i][3] = fac1 * lhsm[j][i][3];
         lhsm[j][i][4] = fac1 * lhsm[j][i][4];
         rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
         lhsm[j1][i][2] = lhsm[j1][i][2] - lhsm[j1][i][1] * lhsm[j][i][3];
         lhsm[j1][i][3] = lhsm[j1][i][3] - lhsm[j1][i][1] * lhsm[j][i][4];
         rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhsm[j1][i][1] * rhs[k][j][i][m];
         //---------------------------------------------------------------------
         // Scale the last row immediately
         //---------------------------------------------------------------------
         rhs[k][j1][i][3] = rhs[k][j1][i][3] / lhsp[j1][i][2];
         rhs[k][j1][i][4] = rhs[k][j1][i][4] / lhsm[j1][i][2];
      }
      //---------------------------------------------------------------------
      // BACKSUBSTITUTION
      //---------------------------------------------------------------------
      j = grid_points[1] - 2;
      j1 = grid_points[1] - 1;
      for(i = 1; i <= grid_points[0] - 2; i++) {
         for(m = 0; m < 3; m++) {
            rhs[k][j][i][m] = rhs[k][j][i][m] - lhs[j][i][3] * rhs[k][j1][i][m];
         }
         rhs[k][j][i][3] = rhs[k][j][i][3] - lhsp[j][i][3] * rhs[k][j1][i][3];
         rhs[k][j][i][4] = rhs[k][j][i][4] - lhsm[j][i][3] * rhs[k][j1][i][4];
      }
      //---------------------------------------------------------------------
      // The first three factors
      //---------------------------------------------------------------------
      for(j = grid_points[1] - 3; j >= 0; j--) {
         j1 = j + 1;
         j2 = j + 2;
         for(i = 1; i <= grid_points[0] - 2; i++) {
            for(m = 0; m < 3; m++) {
               rhs[k][j][i][m] = rhs[k][j][i][m] - lhs[j][i][3] * rhs[k][j1][i][m] - lhs[j][i][4] * rhs[k][j2][i][m];
            }
            //-------------------------------------------------------------------
            // And the remaining two
            //-------------------------------------------------------------------
            rhs[k][j][i][3] = rhs[k][j][i][3] - lhsp[j][i][3] * rhs[k][j1][i][3] - lhsp[j][i][4] * rhs[k][j2][i][3];
            rhs[k][j][i][4] = rhs[k][j][i][4] - lhsm[j][i][3] * rhs[k][j1][i][4] - lhsm[j][i][4] * rhs[k][j2][i][4];
         }
      }
   }
   clava_dcg_global[ 31 ]++;
   pinvr();
}

//---------------------------------------------------------------------
// this function performs the solution of the approximate factorization
// step in the z-direction for all five matrix components
// simultaneously. The Thomas algorithm is employed to solve the
// systems for the z-lines. Boundary conditions are non-periodic
//---------------------------------------------------------------------
void z_solve() {
   int i, j, k, k1, k2, m;
   double ru1, fac1, fac2;
   double rhos[36];
   double cv[36];
   double lhs[37][37][5];
   double lhsp[37][37][5];
   double lhsm[37][37][5];
   for(j = 1; j <= ny2; j++) {
      clava_dcg_global[ 32 ]++;
      lhsinitj(nz2 + 1, nx2, lhs, lhsp, lhsm);
      //---------------------------------------------------------------------
      // Computes the left hand side for the three z-factors
      //---------------------------------------------------------------------
      //---------------------------------------------------------------------
      // first fill the lhs for the u-eigenvalue
      //---------------------------------------------------------------------
      for(i = 1; i <= nx2; i++) {
         for(k = 0; k <= nz2 + 1; k++) {
            ru1 = c3c4 * rho_i[k][j][i];
            cv[k] = ws[k][j][i];
            rhos[k] = ((((dz4 + con43 * ru1) > (dz5 + c1c5 * ru1) ? (dz4 + con43 * ru1) : (dz5 + c1c5 * ru1))) > (((dzmax + ru1) > (dz1) ? (dzmax + ru1) : (dz1))) ? (((dz4 + con43 * ru1) > (dz5 + c1c5 * ru1) ? (dz4 + con43 * ru1) : (dz5 + c1c5 * ru1))) : (((dzmax + ru1) > (dz1) ? (dzmax + ru1) : (dz1))));
         }
         for(k = 1; k <= nz2; k++) {
            lhs[k][i][0] = 0.0;
            lhs[k][i][1] = -dttz2 * cv[k - 1] - dttz1 * rhos[k - 1];
            lhs[k][i][2] = 1.0 + c2dttz1 * rhos[k];
            lhs[k][i][3] = dttz2 * cv[k + 1] - dttz1 * rhos[k + 1];
            lhs[k][i][4] = 0.0;
         }
      }
      //---------------------------------------------------------------------
      // add fourth order dissipation
      //---------------------------------------------------------------------
      for(i = 1; i <= nx2; i++) {
         k = 1;
         lhs[k][i][2] = lhs[k][i][2] + comz5;
         lhs[k][i][3] = lhs[k][i][3] - comz4;
         lhs[k][i][4] = lhs[k][i][4] + comz1;
         k = 2;
         lhs[k][i][1] = lhs[k][i][1] - comz4;
         lhs[k][i][2] = lhs[k][i][2] + comz6;
         lhs[k][i][3] = lhs[k][i][3] - comz4;
         lhs[k][i][4] = lhs[k][i][4] + comz1;
      }
      for(k = 3; k <= nz2 - 2; k++) {
         for(i = 1; i <= nx2; i++) {
            lhs[k][i][0] = lhs[k][i][0] + comz1;
            lhs[k][i][1] = lhs[k][i][1] - comz4;
            lhs[k][i][2] = lhs[k][i][2] + comz6;
            lhs[k][i][3] = lhs[k][i][3] - comz4;
            lhs[k][i][4] = lhs[k][i][4] + comz1;
         }
      }
      for(i = 1; i <= nx2; i++) {
         k = nz2 - 1;
         lhs[k][i][0] = lhs[k][i][0] + comz1;
         lhs[k][i][1] = lhs[k][i][1] - comz4;
         lhs[k][i][2] = lhs[k][i][2] + comz6;
         lhs[k][i][3] = lhs[k][i][3] - comz4;
         k = nz2;
         lhs[k][i][0] = lhs[k][i][0] + comz1;
         lhs[k][i][1] = lhs[k][i][1] - comz4;
         lhs[k][i][2] = lhs[k][i][2] + comz5;
      }
      //---------------------------------------------------------------------
      // subsequently, fill the other factors (u+c), (u-c)
      //---------------------------------------------------------------------
      for(k = 1; k <= nz2; k++) {
         for(i = 1; i <= nx2; i++) {
            lhsp[k][i][0] = lhs[k][i][0];
            lhsp[k][i][1] = lhs[k][i][1] - dttz2 * speed[k - 1][j][i];
            lhsp[k][i][2] = lhs[k][i][2];
            lhsp[k][i][3] = lhs[k][i][3] + dttz2 * speed[k + 1][j][i];
            lhsp[k][i][4] = lhs[k][i][4];
            lhsm[k][i][0] = lhs[k][i][0];
            lhsm[k][i][1] = lhs[k][i][1] + dttz2 * speed[k - 1][j][i];
            lhsm[k][i][2] = lhs[k][i][2];
            lhsm[k][i][3] = lhs[k][i][3] - dttz2 * speed[k + 1][j][i];
            lhsm[k][i][4] = lhs[k][i][4];
         }
      }
      //---------------------------------------------------------------------
      // FORWARD ELIMINATION
      //---------------------------------------------------------------------
      for(k = 0; k <= grid_points[2] - 3; k++) {
         k1 = k + 1;
         k2 = k + 2;
         for(i = 1; i <= nx2; i++) {
            fac1 = 1.0 / lhs[k][i][2];
            lhs[k][i][3] = fac1 * lhs[k][i][3];
            lhs[k][i][4] = fac1 * lhs[k][i][4];
            for(m = 0; m < 3; m++) {
               rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
            }
            lhs[k1][i][2] = lhs[k1][i][2] - lhs[k1][i][1] * lhs[k][i][3];
            lhs[k1][i][3] = lhs[k1][i][3] - lhs[k1][i][1] * lhs[k][i][4];
            for(m = 0; m < 3; m++) {
               rhs[k1][j][i][m] = rhs[k1][j][i][m] - lhs[k1][i][1] * rhs[k][j][i][m];
            }
            lhs[k2][i][1] = lhs[k2][i][1] - lhs[k2][i][0] * lhs[k][i][3];
            lhs[k2][i][2] = lhs[k2][i][2] - lhs[k2][i][0] * lhs[k][i][4];
            for(m = 0; m < 3; m++) {
               rhs[k2][j][i][m] = rhs[k2][j][i][m] - lhs[k2][i][0] * rhs[k][j][i][m];
            }
         }
      }
      //---------------------------------------------------------------------
      // The last two rows in this grid block are a bit different,
      // since they for (not have two more rows available for the
      // elimination of off-diagonal entries
      //---------------------------------------------------------------------
      k = grid_points[2] - 2;
      k1 = grid_points[2] - 1;
      for(i = 1; i <= nx2; i++) {
         fac1 = 1.0 / lhs[k][i][2];
         lhs[k][i][3] = fac1 * lhs[k][i][3];
         lhs[k][i][4] = fac1 * lhs[k][i][4];
         for(m = 0; m < 3; m++) {
            rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
         }
         lhs[k1][i][2] = lhs[k1][i][2] - lhs[k1][i][1] * lhs[k][i][3];
         lhs[k1][i][3] = lhs[k1][i][3] - lhs[k1][i][1] * lhs[k][i][4];
         for(m = 0; m < 3; m++) {
            rhs[k1][j][i][m] = rhs[k1][j][i][m] - lhs[k1][i][1] * rhs[k][j][i][m];
         }
         //---------------------------------------------------------------------
         // scale the last row immediately
         //---------------------------------------------------------------------
         fac2 = 1.0 / lhs[k1][i][2];
         for(m = 0; m < 3; m++) {
            rhs[k1][j][i][m] = fac2 * rhs[k1][j][i][m];
         }
      }
      //---------------------------------------------------------------------
      // for (the u+c and the u-c factors
      //---------------------------------------------------------------------
      for(k = 0; k <= grid_points[2] - 3; k++) {
         k1 = k + 1;
         k2 = k + 2;
         for(i = 1; i <= nx2; i++) {
            m = 3;
            fac1 = 1.0 / lhsp[k][i][2];
            lhsp[k][i][3] = fac1 * lhsp[k][i][3];
            lhsp[k][i][4] = fac1 * lhsp[k][i][4];
            rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
            lhsp[k1][i][2] = lhsp[k1][i][2] - lhsp[k1][i][1] * lhsp[k][i][3];
            lhsp[k1][i][3] = lhsp[k1][i][3] - lhsp[k1][i][1] * lhsp[k][i][4];
            rhs[k1][j][i][m] = rhs[k1][j][i][m] - lhsp[k1][i][1] * rhs[k][j][i][m];
            lhsp[k2][i][1] = lhsp[k2][i][1] - lhsp[k2][i][0] * lhsp[k][i][3];
            lhsp[k2][i][2] = lhsp[k2][i][2] - lhsp[k2][i][0] * lhsp[k][i][4];
            rhs[k2][j][i][m] = rhs[k2][j][i][m] - lhsp[k2][i][0] * rhs[k][j][i][m];
            m = 4;
            fac1 = 1.0 / lhsm[k][i][2];
            lhsm[k][i][3] = fac1 * lhsm[k][i][3];
            lhsm[k][i][4] = fac1 * lhsm[k][i][4];
            rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
            lhsm[k1][i][2] = lhsm[k1][i][2] - lhsm[k1][i][1] * lhsm[k][i][3];
            lhsm[k1][i][3] = lhsm[k1][i][3] - lhsm[k1][i][1] * lhsm[k][i][4];
            rhs[k1][j][i][m] = rhs[k1][j][i][m] - lhsm[k1][i][1] * rhs[k][j][i][m];
            lhsm[k2][i][1] = lhsm[k2][i][1] - lhsm[k2][i][0] * lhsm[k][i][3];
            lhsm[k2][i][2] = lhsm[k2][i][2] - lhsm[k2][i][0] * lhsm[k][i][4];
            rhs[k2][j][i][m] = rhs[k2][j][i][m] - lhsm[k2][i][0] * rhs[k][j][i][m];
         }
      }
      //---------------------------------------------------------------------
      // And again the last two rows separately
      //---------------------------------------------------------------------
      k = grid_points[2] - 2;
      k1 = grid_points[2] - 1;
      for(i = 1; i <= nx2; i++) {
         m = 3;
         fac1 = 1.0 / lhsp[k][i][2];
         lhsp[k][i][3] = fac1 * lhsp[k][i][3];
         lhsp[k][i][4] = fac1 * lhsp[k][i][4];
         rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
         lhsp[k1][i][2] = lhsp[k1][i][2] - lhsp[k1][i][1] * lhsp[k][i][3];
         lhsp[k1][i][3] = lhsp[k1][i][3] - lhsp[k1][i][1] * lhsp[k][i][4];
         rhs[k1][j][i][m] = rhs[k1][j][i][m] - lhsp[k1][i][1] * rhs[k][j][i][m];
         m = 4;
         fac1 = 1.0 / lhsm[k][i][2];
         lhsm[k][i][3] = fac1 * lhsm[k][i][3];
         lhsm[k][i][4] = fac1 * lhsm[k][i][4];
         rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
         lhsm[k1][i][2] = lhsm[k1][i][2] - lhsm[k1][i][1] * lhsm[k][i][3];
         lhsm[k1][i][3] = lhsm[k1][i][3] - lhsm[k1][i][1] * lhsm[k][i][4];
         rhs[k1][j][i][m] = rhs[k1][j][i][m] - lhsm[k1][i][1] * rhs[k][j][i][m];
         //---------------------------------------------------------------------
         // Scale the last row immediately (some of this is overkill
         // if this is the last cell)
         //---------------------------------------------------------------------
         rhs[k1][j][i][3] = rhs[k1][j][i][3] / lhsp[k1][i][2];
         rhs[k1][j][i][4] = rhs[k1][j][i][4] / lhsm[k1][i][2];
      }
      //---------------------------------------------------------------------
      // BACKSUBSTITUTION
      //---------------------------------------------------------------------
      k = grid_points[2] - 2;
      k1 = grid_points[2] - 1;
      for(i = 1; i <= nx2; i++) {
         for(m = 0; m < 3; m++) {
            rhs[k][j][i][m] = rhs[k][j][i][m] - lhs[k][i][3] * rhs[k1][j][i][m];
         }
         rhs[k][j][i][3] = rhs[k][j][i][3] - lhsp[k][i][3] * rhs[k1][j][i][3];
         rhs[k][j][i][4] = rhs[k][j][i][4] - lhsm[k][i][3] * rhs[k1][j][i][4];
      }
      //---------------------------------------------------------------------
      // Whether or not this is the last processor, we always have
      // to complete the back-substitution
      //---------------------------------------------------------------------
      //---------------------------------------------------------------------
      // The first three factors
      //---------------------------------------------------------------------
      for(k = grid_points[2] - 3; k >= 0; k--) {
         k1 = k + 1;
         k2 = k + 2;
         for(i = 1; i <= nx2; i++) {
            for(m = 0; m < 3; m++) {
               rhs[k][j][i][m] = rhs[k][j][i][m] - lhs[k][i][3] * rhs[k1][j][i][m] - lhs[k][i][4] * rhs[k2][j][i][m];
            }
            //-------------------------------------------------------------------
            // And the remaining two
            //-------------------------------------------------------------------
            rhs[k][j][i][3] = rhs[k][j][i][3] - lhsp[k][i][3] * rhs[k1][j][i][3] - lhsp[k][i][4] * rhs[k2][j][i][3];
            rhs[k][j][i][4] = rhs[k][j][i][4] - lhsm[k][i][3] * rhs[k1][j][i][4] - lhsm[k][i][4] * rhs[k2][j][i][4];
         }
      }
   }
   clava_dcg_global[ 33 ]++;
   tzetar();
}

void set_constants() {
   ce[0][0] = 2.0;
   ce[0][1] = 0.0;
   ce[0][2] = 0.0;
   ce[0][3] = 4.0;
   ce[0][4] = 5.0;
   ce[0][5] = 3.0;
   ce[0][6] = 0.5;
   ce[0][7] = 0.02;
   ce[0][8] = 0.01;
   ce[0][9] = 0.03;
   ce[0][10] = 0.5;
   ce[0][11] = 0.4;
   ce[0][12] = 0.3;
   ce[1][0] = 1.0;
   ce[1][1] = 0.0;
   ce[1][2] = 0.0;
   ce[1][3] = 0.0;
   ce[1][4] = 1.0;
   ce[1][5] = 2.0;
   ce[1][6] = 3.0;
   ce[1][7] = 0.01;
   ce[1][8] = 0.03;
   ce[1][9] = 0.02;
   ce[1][10] = 0.4;
   ce[1][11] = 0.3;
   ce[1][12] = 0.5;
   ce[2][0] = 2.0;
   ce[2][1] = 2.0;
   ce[2][2] = 0.0;
   ce[2][3] = 0.0;
   ce[2][4] = 0.0;
   ce[2][5] = 2.0;
   ce[2][6] = 3.0;
   ce[2][7] = 0.04;
   ce[2][8] = 0.03;
   ce[2][9] = 0.05;
   ce[2][10] = 0.3;
   ce[2][11] = 0.5;
   ce[2][12] = 0.4;
   ce[3][0] = 2.0;
   ce[3][1] = 2.0;
   ce[3][2] = 0.0;
   ce[3][3] = 0.0;
   ce[3][4] = 0.0;
   ce[3][5] = 2.0;
   ce[3][6] = 3.0;
   ce[3][7] = 0.03;
   ce[3][8] = 0.05;
   ce[3][9] = 0.04;
   ce[3][10] = 0.2;
   ce[3][11] = 0.1;
   ce[3][12] = 0.3;
   ce[4][0] = 5.0;
   ce[4][1] = 4.0;
   ce[4][2] = 3.0;
   ce[4][3] = 2.0;
   ce[4][4] = 0.1;
   ce[4][5] = 0.4;
   ce[4][6] = 0.3;
   ce[4][7] = 0.05;
   ce[4][8] = 0.04;
   ce[4][9] = 0.03;
   ce[4][10] = 0.1;
   ce[4][11] = 0.3;
   ce[4][12] = 0.2;
   c1 = 1.4;
   c2 = 0.4;
   c3 = 0.1;
   c4 = 1.0;
   c5 = 1.4;
   clava_dcg_global[ 34 ]++;
   bt = sqrt(0.5);
   dnxm1 = 1.0 / (double) (grid_points[0] - 1);
   dnym1 = 1.0 / (double) (grid_points[1] - 1);
   dnzm1 = 1.0 / (double) (grid_points[2] - 1);
   c1c2 = c1 * c2;
   c1c5 = c1 * c5;
   c3c4 = c3 * c4;
   c1345 = c1c5 * c3c4;
   conz1 = (1.0 - c1c5);
   tx1 = 1.0 / (dnxm1 * dnxm1);
   tx2 = 1.0 / (2.0 * dnxm1);
   tx3 = 1.0 / dnxm1;
   ty1 = 1.0 / (dnym1 * dnym1);
   ty2 = 1.0 / (2.0 * dnym1);
   ty3 = 1.0 / dnym1;
   tz1 = 1.0 / (dnzm1 * dnzm1);
   tz2 = 1.0 / (2.0 * dnzm1);
   tz3 = 1.0 / dnzm1;
   dx1 = 0.75;
   dx2 = 0.75;
   dx3 = 0.75;
   dx4 = 0.75;
   dx5 = 0.75;
   dy1 = 0.75;
   dy2 = 0.75;
   dy3 = 0.75;
   dy4 = 0.75;
   dy5 = 0.75;
   dz1 = 1.0;
   dz2 = 1.0;
   dz3 = 1.0;
   dz4 = 1.0;
   dz5 = 1.0;
   dxmax = ((dx3) > (dx4) ? (dx3) : (dx4));
   dymax = ((dy2) > (dy4) ? (dy2) : (dy4));
   dzmax = ((dz2) > (dz3) ? (dz2) : (dz3));
   dssp = 0.25 * ((dx1) > (((dy1) > (dz1) ? (dy1) : (dz1))) ? (dx1) : (((dy1) > (dz1) ? (dy1) : (dz1))));
   c4dssp = 4.0 * dssp;
   c5dssp = 5.0 * dssp;
   dttx1 = dt * tx1;
   dttx2 = dt * tx2;
   dtty1 = dt * ty1;
   dtty2 = dt * ty2;
   dttz1 = dt * tz1;
   dttz2 = dt * tz2;
   c2dttx1 = 2.0 * dttx1;
   c2dtty1 = 2.0 * dtty1;
   c2dttz1 = 2.0 * dttz1;
   dtdssp = dt * dssp;
   comz1 = dtdssp;
   comz4 = 4.0 * dtdssp;
   comz5 = 5.0 * dtdssp;
   comz6 = 6.0 * dtdssp;
   c3c4tx3 = c3c4 * tx3;
   c3c4ty3 = c3c4 * ty3;
   c3c4tz3 = c3c4 * tz3;
   dx1tx1 = dx1 * tx1;
   dx2tx1 = dx2 * tx1;
   dx3tx1 = dx3 * tx1;
   dx4tx1 = dx4 * tx1;
   dx5tx1 = dx5 * tx1;
   dy1ty1 = dy1 * ty1;
   dy2ty1 = dy2 * ty1;
   dy3ty1 = dy3 * ty1;
   dy4ty1 = dy4 * ty1;
   dy5ty1 = dy5 * ty1;
   dz1tz1 = dz1 * tz1;
   dz2tz1 = dz2 * tz1;
   dz3tz1 = dz3 * tz1;
   dz4tz1 = dz4 * tz1;
   dz5tz1 = dz5 * tz1;
   c2iv = 2.5;
   con43 = 4.0 / 3.0;
   con16 = 1.0 / 6.0;
   xxcon1 = c3c4tx3 * con43 * tx3;
   xxcon2 = c3c4tx3 * tx3;
   xxcon3 = c3c4tx3 * conz1 * tx3;
   xxcon4 = c3c4tx3 * con16 * tx3;
   xxcon5 = c3c4tx3 * c1c5 * tx3;
   yycon1 = c3c4ty3 * con43 * ty3;
   yycon2 = c3c4ty3 * ty3;
   yycon3 = c3c4ty3 * conz1 * ty3;
   yycon4 = c3c4ty3 * con16 * ty3;
   yycon5 = c3c4ty3 * c1c5 * ty3;
   zzcon1 = c3c4tz3 * con43 * tz3;
   zzcon2 = c3c4tz3 * tz3;
   zzcon3 = c3c4tz3 * conz1 * tz3;
   zzcon4 = c3c4tz3 * con16 * tz3;
   zzcon5 = c3c4tz3 * c1c5 * tz3;
}

void print_results(char *name, char class, int n1, int n2, int n3, int niter, double t, double mops, char *optype, int verified) {
   char size[16];
   int j;
   clava_dcg_global[ 35 ]++;
   printf("\n\n %s Benchmark Completed.\n", name);
   clava_dcg_global[ 35 ]++;
   printf(" Class           =             %12c\n", class);
   // If this is not a grid-based problem (EP, FT, CG), then
   // we only print n1, which contains some measure of the
   // problem size. In that case, n2 and n3 are both zero.
   // Otherwise, we print the grid size n1xn2xn3
   if((n2 == 0) && (n3 == 0)) {
      if((name[0] == 'E') && (name[1] == 'P')) {
         clava_dcg_global[ 36 ]++;
         clava_dcg_global[ 37 ]++;
         sprintf(size, "%15.0lf", pow(2.0, n1));
         j = 14;
         if(size[j] == '.') {
            size[j] = ' ';
            j--;
         }
         size[j + 1] = '\0';
         clava_dcg_global[ 35 ]++;
         printf(" Size            =          %15s\n", size);
      }
      else {
         clava_dcg_global[ 35 ]++;
         printf(" Size            =             %12d\n", n1);
      }
   }
   else {
      clava_dcg_global[ 35 ]++;
      printf(" Size            =           %4dx%4dx%4d\n", n1, n2, n3);
   }
   clava_dcg_global[ 35 ]++;
   printf(" Iterations      =             %12d\n", niter);
   clava_dcg_global[ 35 ]++;
   printf(" Time in seconds =             %12.4lf\n", t);
   clava_dcg_global[ 35 ]++;
   printf(" Mop/s total     =          %15.2lf\n", mops);
   clava_dcg_global[ 35 ]++;
   printf(" Operation type  = %24s\n", optype);
   if(verified) {
      clava_dcg_global[ 35 ]++;
      printf(" Verification    =             %12s\n", "SUCCESSFUL");
   }
   else {
      clava_dcg_global[ 35 ]++;
      printf(" Verification    =             %12s\n", "UNSUCCESSFUL");
   }
}

void wtime(double *t) {
   static int sec = -1;
   struct timeval tv;
   clava_dcg_global[ 38 ]++;
   gettimeofday(&tv, (void *) 0);
   if(sec < 0) sec = tv.tv_sec;
   *t = (tv.tv_sec - sec) + 1.0e-6 * tv.tv_usec;
}

/*****************************************************************/
/******         E  L  A  P  S  E  D  _  T  I  M  E          ******/
/*****************************************************************/
double elapsed_time() {
   double t;
   clava_dcg_global[ 39 ]++;
   wtime(&t);
   
   return (t);
}

/*****************************************************************/
/******            T  I  M  E  R  _  C  L  E  A  R          ******/
/*****************************************************************/
void timer_clear(int n) {
   elapsed[n] = 0.0;
}

/*****************************************************************/
/******            T  I  M  E  R  _  S  T  A  R  T          ******/
/*****************************************************************/
void timer_start(int n) {
   clava_dcg_global[ 40 ]++;
   start[n] = elapsed_time();
}

/*****************************************************************/
/******            T  I  M  E  R  _  S  T  O  P             ******/
/*****************************************************************/
void timer_stop(int n) {
   double t, now;
   clava_dcg_global[ 41 ]++;
   now = elapsed_time();
   t = now - start[n];
   elapsed[n] += t;
}

/*****************************************************************/
/******            T  I  M  E  R  _  R  E  A  D             ******/
/*****************************************************************/
double timer_read(int n) {
   
   return (elapsed[n]);
}
