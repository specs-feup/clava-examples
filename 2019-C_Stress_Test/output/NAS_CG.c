#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
//---------------------------------------------------------------------
// program CG
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

struct anon_NAS_CG_c_75 {
   double real;
   double imag;
};

typedef struct anon_NAS_CG_c_75 dcomplex;
//---------------------------------------------------------------------

/*common / main_int_mem /*/

int colidx[567000];
int rowstr[7001];
int iv[7000];
int arow[7000];
int acol[63000];
/*common / main_flt_mem /*/

double aelt[63000];
double a[567000];
double x[7002];
double z[7002];
double p[7002];
double q[7002];
double r[7002];
/*common / partit_size /*/

int naa;
int nzz;
int firstrow;
int lastrow;
int firstcol;
int lastcol;
/*common /urando/*/

double amult;
double tran;
//---------------------------------------------------------------------
int clava_dcg_global[ 27 ] = {0};
void conj_grad(int colidx[], int rowstr[], double x[], double z[], double a[], double p[], double q[], double r[], double *rnorm);
void makea(int n, int nz, double a[], int colidx[], int rowstr[], int firstrow, int lastrow, int firstcol, int lastcol, int arow[], int acol[][9], double aelt[][9], int iv[]);
void sparse(double a[], int colidx[], int rowstr[], int n, int nz, int nozer, int arow[], int acol[][9], double aelt[][9], int firstrow, int lastrow, int nzloc[], double rcond, double shift);
void sprnvc(int n, int nz, int nn1, double v[], int iv[]);
int icnvrt(double x, int ipwr2);
void vecset(int n, double v[], int iv[], int *nzv, int i, double val);
void print_results(char *name, char class, int n1, int n2, int n3, int niter, double t, double mops, char *optype, int verified);
double randlc(double *x, double a);
void vranlc(int n, double *x, double a, double y[]);
double start[64];
double elapsed[64];
double elapsed_time();
void timer_clear(int n);
void timer_start(int n);
void timer_stop(int n);
double timer_read(int n);
void wtime(double *t);
//---------------------------------------------------------------------
void clava_call_graph() {
   FILE *log_file_361 = fopen("/home/specs/jbispo/repos/clava-examples/2019-C_Stress_Test/output/NAS_CG.dot", "w+");
   if (log_file_361 == NULL)
   {
       printf("Error opening file /home/specs/jbispo/repos/clava-examples/2019-C_Stress_Test/output/NAS_CG.dot\n");
       exit(1);
   } 
   fprintf(log_file_361, "digraph dynamic_call_graph {\n\n");
   if(clava_dcg_global[0] != 0) {
      fprintf(log_file_361, "	main -> timer_clear [label=\"%d\"];\n", clava_dcg_global[0]);
   }
   if(clava_dcg_global[1] != 0) {
      fprintf(log_file_361, "	main -> timer_start [label=\"%d\"];\n", clava_dcg_global[1]);
   }
   if(clava_dcg_global[2] != 0) {
      fprintf(log_file_361, "	main -> printf [label=\"%d\"];\n", clava_dcg_global[2]);
   }
   if(clava_dcg_global[3] != 0) {
      fprintf(log_file_361, "	main -> randlc [label=\"%d\"];\n", clava_dcg_global[3]);
   }
   if(clava_dcg_global[4] != 0) {
      fprintf(log_file_361, "	main -> makea [label=\"%d\"];\n", clava_dcg_global[4]);
   }
   if(clava_dcg_global[5] != 0) {
      fprintf(log_file_361, "	main -> conj_grad [label=\"%d\"];\n", clava_dcg_global[5]);
   }
   if(clava_dcg_global[6] != 0) {
      fprintf(log_file_361, "	main -> sqrt [label=\"%d\"];\n", clava_dcg_global[6]);
   }
   if(clava_dcg_global[7] != 0) {
      fprintf(log_file_361, "	main -> timer_stop [label=\"%d\"];\n", clava_dcg_global[7]);
   }
   if(clava_dcg_global[8] != 0) {
      fprintf(log_file_361, "	main -> timer_read [label=\"%d\"];\n", clava_dcg_global[8]);
   }
   if(clava_dcg_global[9] != 0) {
      fprintf(log_file_361, "	main -> fabs [label=\"%d\"];\n", clava_dcg_global[9]);
   }
   if(clava_dcg_global[10] != 0) {
      fprintf(log_file_361, "	main -> print_results [label=\"%d\"];\n", clava_dcg_global[10]);
   }
   if(clava_dcg_global[11] != 0) {
      fprintf(log_file_361, "	conj_grad -> sqrt [label=\"%d\"];\n", clava_dcg_global[11]);
   }
   if(clava_dcg_global[12] != 0) {
      fprintf(log_file_361, "	makea -> sprnvc [label=\"%d\"];\n", clava_dcg_global[12]);
   }
   if(clava_dcg_global[13] != 0) {
      fprintf(log_file_361, "	makea -> vecset [label=\"%d\"];\n", clava_dcg_global[13]);
   }
   if(clava_dcg_global[14] != 0) {
      fprintf(log_file_361, "	makea -> sparse [label=\"%d\"];\n", clava_dcg_global[14]);
   }
   if(clava_dcg_global[15] != 0) {
      fprintf(log_file_361, "	sparse -> printf [label=\"%d\"];\n", clava_dcg_global[15]);
   }
   if(clava_dcg_global[16] != 0) {
      fprintf(log_file_361, "	sparse -> exit [label=\"%d\"];\n", clava_dcg_global[16]);
   }
   if(clava_dcg_global[17] != 0) {
      fprintf(log_file_361, "	sparse -> pow [label=\"%d\"];\n", clava_dcg_global[17]);
   }
   if(clava_dcg_global[18] != 0) {
      fprintf(log_file_361, "	sprnvc -> randlc [label=\"%d\"];\n", clava_dcg_global[18]);
   }
   if(clava_dcg_global[19] != 0) {
      fprintf(log_file_361, "	sprnvc -> icnvrt [label=\"%d\"];\n", clava_dcg_global[19]);
   }
   if(clava_dcg_global[20] != 0) {
      fprintf(log_file_361, "	wtime -> gettimeofday [label=\"%d\"];\n", clava_dcg_global[20]);
   }
   if(clava_dcg_global[21] != 0) {
      fprintf(log_file_361, "	elapsed_time -> wtime [label=\"%d\"];\n", clava_dcg_global[21]);
   }
   if(clava_dcg_global[22] != 0) {
      fprintf(log_file_361, "	timer_start -> elapsed_time [label=\"%d\"];\n", clava_dcg_global[22]);
   }
   if(clava_dcg_global[23] != 0) {
      fprintf(log_file_361, "	timer_stop -> elapsed_time [label=\"%d\"];\n", clava_dcg_global[23]);
   }
   if(clava_dcg_global[24] != 0) {
      fprintf(log_file_361, "	print_results -> printf [label=\"%d\"];\n", clava_dcg_global[24]);
   }
   if(clava_dcg_global[25] != 0) {
      fprintf(log_file_361, "	print_results -> sprintf [label=\"%d\"];\n", clava_dcg_global[25]);
   }
   if(clava_dcg_global[26] != 0) {
      fprintf(log_file_361, "	print_results -> pow [label=\"%d\"];\n", clava_dcg_global[26]);
   }
   fprintf(log_file_361, "}\n");
   fclose(log_file_361);
}

int main(int argc, char *argv[]) {
   atexit(clava_call_graph);
   int i, j, k, it;
   double zeta;
   double rnorm;
   double norm_temp1, norm_temp2;
   double t, mflops, tmax;
   char Class;
   int verified;
   double zeta_verify_value, epsilon, err;
   char *t_names[3];
   for(i = 0; i < 3; i++) {
      clava_dcg_global[ 0 ]++;
      timer_clear(i);
   }
   clava_dcg_global[ 1 ]++;
   timer_start(0);
   firstrow = 0;
   lastrow = 7000 - 1;
   firstcol = 0;
   lastcol = 7000 - 1;
   if(7000 == 1400 && 8 == 7 && 15 == 15 && 12.0 == 10) {
      Class = 'S';
      zeta_verify_value = 8.5971775078648;
   }
   else if(7000 == 7000 && 8 == 8 && 15 == 15 && 12.0 == 12) {
      Class = 'W';
      zeta_verify_value = 10.362595087124;
   }
   else if(7000 == 14000 && 8 == 11 && 15 == 15 && 12.0 == 20) {
      Class = 'A';
      zeta_verify_value = 17.130235054029;
   }
   else if(7000 == 75000 && 8 == 13 && 15 == 75 && 12.0 == 60) {
      Class = 'B';
      zeta_verify_value = 22.712745482631;
   }
   else if(7000 == 150000 && 8 == 15 && 15 == 75 && 12.0 == 110) {
      Class = 'C';
      zeta_verify_value = 28.973605592845;
   }
   else if(7000 == 1500000 && 8 == 21 && 15 == 100 && 12.0 == 500) {
      Class = 'D';
      zeta_verify_value = 52.514532105794;
   }
   else if(7000 == 9000000 && 8 == 26 && 15 == 100 && 12.0 == 1500) {
      Class = 'E';
      zeta_verify_value = 77.522164599383;
   }
   else {
      Class = 'U';
   }
   clava_dcg_global[ 2 ]++;
   printf("\n\n NAS Parallel Benchmarks (NPB3.3-SER-C) - CG Benchmark\n\n");
   clava_dcg_global[ 2 ]++;
   printf(" Size: %11d\n", 7000);
   clava_dcg_global[ 2 ]++;
   printf(" Iterations: %5d\n", 15);
   clava_dcg_global[ 2 ]++;
   printf("\n");
   naa = 7000;
   nzz = (7000 * (8 + 1) * (8 + 1));
   //---------------------------------------------------------------------
   // Inialize random number generator
   //---------------------------------------------------------------------
   tran = 314159265.0;
   amult = 1220703125.0;
   clava_dcg_global[ 3 ]++;
   zeta = randlc(&tran, amult);
   //---------------------------------------------------------------------
   //
   //---------------------------------------------------------------------
   clava_dcg_global[ 4 ]++;
   makea(naa, nzz, a, colidx, rowstr, firstrow, lastrow, firstcol, lastcol, arow, (int (*)[9]) (void *) acol, (double (*)[9]) (void *) aelt, iv);
   //---------------------------------------------------------------------
   // Note: as a result of the above call to makea:
   //      values of j used in indexing rowstr go from 0 --> lastrow-firstrow
   //      values of colidx which are col indexes go from firstcol --> lastcol
   //      So:
   //      Shift the col index vals from actual (firstcol --> lastcol )
   //      to local, i.e., (0 --> lastcol-firstcol)
   //---------------------------------------------------------------------
   for(j = 0; j < lastrow - firstrow + 1; j++) {
      for(k = rowstr[j]; k < rowstr[j + 1]; k++) {
         colidx[k] = colidx[k] - firstcol;
      }
   }
   //---------------------------------------------------------------------
   // set starting vector to (1, 1, .... 1)
   //---------------------------------------------------------------------
   for(i = 0; i < 7000 + 1; i++) {
      x[i] = 1.0;
   }
   for(j = 0; j < lastcol - firstcol + 1; j++) {
      q[j] = 0.0;
      z[j] = 0.0;
      r[j] = 0.0;
      p[j] = 0.0;
   }
   zeta = 0.0;
   //---------------------------------------------------------------------
   //---->
   // Do one iteration untimed to init all code and data page tables
   //---->                    (then reinit, start timing, to niter its)
   //---------------------------------------------------------------------
   for(it = 1; it <= 1; it++) { // end of do one iteration untimed
      //---------------------------------------------------------------------
      // The call to the conjugate gradient routine:
      //---------------------------------------------------------------------
      clava_dcg_global[ 5 ]++;
      conj_grad(colidx, rowstr, x, z, a, p, q, r, &rnorm);
      //---------------------------------------------------------------------
      // zeta = shift + 1/(x.z)
      // So, first: (x.z)
      // Also, find norm of z
      // So, first: (z.z)
      //---------------------------------------------------------------------
      norm_temp1 = 0.0;
      norm_temp2 = 0.0;
      for(j = 0; j < lastcol - firstcol + 1; j++) {
         norm_temp1 = norm_temp1 + x[j] * z[j];
         norm_temp2 = norm_temp2 + z[j] * z[j];
      }
      clava_dcg_global[ 6 ]++;
      norm_temp2 = 1.0 / sqrt(norm_temp2);
      //---------------------------------------------------------------------
      // Normalize z to obtain x
      //---------------------------------------------------------------------
      for(j = 0; j < lastcol - firstcol + 1; j++) {
         x[j] = norm_temp2 * z[j];
      }
   }
   //---------------------------------------------------------------------
   // set starting vector to (1, 1, .... 1)
   //---------------------------------------------------------------------
   for(i = 0; i < 7000 + 1; i++) {
      x[i] = 1.0;
   }
   zeta = 0.0;
   clava_dcg_global[ 7 ]++;
   timer_stop(0);
   clava_dcg_global[ 2 ]++;
   clava_dcg_global[ 8 ]++;
   printf(" Initialization time = %15.3f seconds\n", timer_read(0));
   clava_dcg_global[ 1 ]++;
   timer_start(1);
   //---------------------------------------------------------------------
   //---->
   // Main Iteration for inverse power method
   //---->
   //---------------------------------------------------------------------
   for(it = 1; it <= 15; it++) { // end of main iter inv pow meth
      //---------------------------------------------------------------------
      // The call to the conjugate gradient routine:
      //---------------------------------------------------------------------
      clava_dcg_global[ 5 ]++;
      conj_grad(colidx, rowstr, x, z, a, p, q, r, &rnorm);
      //---------------------------------------------------------------------
      // zeta = shift + 1/(x.z)
      // So, first: (x.z)
      // Also, find norm of z
      // So, first: (z.z)
      //---------------------------------------------------------------------
      norm_temp1 = 0.0;
      norm_temp2 = 0.0;
      for(j = 0; j < lastcol - firstcol + 1; j++) {
         norm_temp1 = norm_temp1 + x[j] * z[j];
         norm_temp2 = norm_temp2 + z[j] * z[j];
      }
      clava_dcg_global[ 6 ]++;
      norm_temp2 = 1.0 / sqrt(norm_temp2);
      zeta = 12.0 + 1.0 / norm_temp1;
      if(it == 1) {
         clava_dcg_global[ 2 ]++;
         printf("\n   iteration           ||r||                 zeta\n");
      }
      clava_dcg_global[ 2 ]++;
      printf("    %5d       %20.14E%20.13f\n", it, rnorm, zeta);
      //---------------------------------------------------------------------
      // Normalize z to obtain x
      //---------------------------------------------------------------------
      for(j = 0; j < lastcol - firstcol + 1; j++) {
         x[j] = norm_temp2 * z[j];
      }
   }
   clava_dcg_global[ 7 ]++;
   timer_stop(1);
   //---------------------------------------------------------------------
   // End of timed section
   //---------------------------------------------------------------------
   clava_dcg_global[ 8 ]++;
   t = timer_read(1);
   clava_dcg_global[ 2 ]++;
   printf(" Benchmark completed\n");
   epsilon = 1.0e-10;
   if(Class != 'U') {
      clava_dcg_global[ 9 ]++;
      err = fabs(zeta - zeta_verify_value) / zeta_verify_value;
      if(err <= epsilon) {
         verified = 1;
         clava_dcg_global[ 2 ]++;
         printf(" VERIFICATION SUCCESSFUL\n");
         clava_dcg_global[ 2 ]++;
         printf(" Zeta is    %20.13E\n", zeta);
         clava_dcg_global[ 2 ]++;
         printf(" Error is   %20.13E\n", err);
      }
      else {
         verified = 0;
         clava_dcg_global[ 2 ]++;
         printf(" VERIFICATION FAILED\n");
         clava_dcg_global[ 2 ]++;
         printf(" Zeta                %20.13E\n", zeta);
         clava_dcg_global[ 2 ]++;
         printf(" The correct zeta is %20.13E\n", zeta_verify_value);
      }
   }
   else {
      verified = 0;
      clava_dcg_global[ 2 ]++;
      printf(" Problem size unknown\n");
      clava_dcg_global[ 2 ]++;
      printf(" NO VERIFICATION PERFORMED\n");
   }
   if(t != 0.0) {
      mflops = (double) (2 * 15 * 7000) * (3.0 + (double) (8 * (8 + 1)) + 25.0 * (5.0 + (double) (8 * (8 + 1))) + 3.0) / t / 1000000.0;
   }
   else {
      mflops = 0.0;
   }
   clava_dcg_global[ 10 ]++;
   print_results("CG", Class, 7000, 0, 0, 15, t, mflops, "          floating point", verified);
   int exitValue = verified ? 0 : 1;
   
   return exitValue;
}

//---------------------------------------------------------------------
// Floaging point arrays here are named as in NPB1 spec discussion of
// CG algorithm
//---------------------------------------------------------------------
void conj_grad(int colidx[], int rowstr[], double x[], double z[], double a[], double p[], double q[], double r[], double *rnorm) {
   int j, k;
   int cgit, cgitmax = 25;
   double d, sum, rho, rho0, alpha, beta;
   rho = 0.0;
   //---------------------------------------------------------------------
   // Initialize the CG algorithm:
   //---------------------------------------------------------------------
   for(j = 0; j < naa + 1; j++) {
      q[j] = 0.0;
      z[j] = 0.0;
      r[j] = x[j];
      p[j] = r[j];
   }
   //---------------------------------------------------------------------
   // rho = r.r
   // Now, obtain the norm of r: First, sum squares of r elements locally...
   //---------------------------------------------------------------------
   for(j = 0; j < lastcol - firstcol + 1; j++) {
      rho = rho + r[j] * r[j];
   }
   //---------------------------------------------------------------------
   //---->
   // The conj grad iteration loop
   //---->
   //---------------------------------------------------------------------
   for(cgit = 1; cgit <= cgitmax; cgit++) { // end of do cgit=1,cgitmax
      //---------------------------------------------------------------------
      // q = A.p
      // The partition submatrix-vector multiply: use workspace w
      //---------------------------------------------------------------------
      //
      // NOTE: this version of the multiply is actually (slightly: maybe %5)
      //       faster on the sp2 on 16 nodes than is the unrolled-by-2 version
      //       below.   On the Cray t3d, the reverse is 1, i.e., the
      //       unrolled-by-two version is some 10% faster.
      //       The unrolled-by-8 version below is significantly faster
      //       on the Cray t3d - overall speed of code is 1.5 times faster.
      for(j = 0; j < lastrow - firstrow + 1; j++) {
         sum = 0.0;
         for(k = rowstr[j]; k < rowstr[j + 1]; k++) {
            sum = sum + a[k] * p[colidx[k]];
         }
         q[j] = sum;
      }
      /*
      for (j = 0; j < lastrow - firstrow + 1; j++) {
      int i = rowstr[j];
      int iresidue = (rowstr[j+1] - i) % 2;
      double sum1 = 0.0;
      double sum2 = 0.0;
      if (iresidue == 1)
      sum1 = sum1 + a[i]*p[colidx[i]];
      for (k = i + iresidue; k <= rowstr[j+1] - 2; k += 2) {
      sum1 = sum1 + a[k]  *p[colidx[k]];
      sum2 = sum2 + a[k+1]*p[colidx[k+1]];
      }
      q[j] = sum1 + sum2;
      }
      */
      /*
      for (j = 0; j < lastrow - firstrow + 1; j++) {
      int i = rowstr[j];
      int iresidue = (rowstr[j+1] - i) % 8;
      double sum = 0.0;
      for (k = i; k <= i + iresidue - 1; k++) {
      sum = sum + a[k]*p[colidx[k]];
      }
      for (k = i + iresidue; k <= rowstr[j+1] - 8; k += 8) {
      sum = sum + a[k  ]*p[colidx[k  ]]
      + a[k+1]*p[colidx[k+1]]
      + a[k+2]*p[colidx[k+2]]
      + a[k+3]*p[colidx[k+3]]
      + a[k+4]*p[colidx[k+4]]
      + a[k+5]*p[colidx[k+5]]
      + a[k+6]*p[colidx[k+6]]
      + a[k+7]*p[colidx[k+7]];
      }
      q[j] = sum;
      }
      */
      //---------------------------------------------------------------------
      // Obtain p.q
      //---------------------------------------------------------------------
      d = 0.0;
      for(j = 0; j < lastcol - firstcol + 1; j++) {
         d = d + p[j] * q[j];
      }
      //---------------------------------------------------------------------
      // Obtain alpha = rho / (p.q)
      //---------------------------------------------------------------------
      alpha = rho / d;
      //---------------------------------------------------------------------
      // Save a temporary of rho
      //---------------------------------------------------------------------
      rho0 = rho;
      //---------------------------------------------------------------------
      // Obtain z = z + alpha*p
      // and    r = r - alpha*q
      //---------------------------------------------------------------------
      rho = 0.0;
      for(j = 0; j < lastcol - firstcol + 1; j++) {
         z[j] = z[j] + alpha * p[j];
         r[j] = r[j] - alpha * q[j];
      }
      //---------------------------------------------------------------------
      // rho = r.r
      // Now, obtain the norm of r: First, sum squares of r elements locally...
      //---------------------------------------------------------------------
      for(j = 0; j < lastcol - firstcol + 1; j++) {
         rho = rho + r[j] * r[j];
      }
      //---------------------------------------------------------------------
      // Obtain beta:
      //---------------------------------------------------------------------
      beta = rho / rho0;
      //---------------------------------------------------------------------
      // p = r + beta*p
      //---------------------------------------------------------------------
      for(j = 0; j < lastcol - firstcol + 1; j++) {
         p[j] = r[j] + beta * p[j];
      }
   }
   //---------------------------------------------------------------------
   // Compute residual norm explicitly:  ||r|| = ||x - A.z||
   // First, form A.z
   // The partition submatrix-vector multiply
   //---------------------------------------------------------------------
   sum = 0.0;
   for(j = 0; j < lastrow - firstrow + 1; j++) {
      d = 0.0;
      for(k = rowstr[j]; k < rowstr[j + 1]; k++) {
         d = d + a[k] * z[colidx[k]];
      }
      r[j] = d;
   }
   //---------------------------------------------------------------------
   // At this point, r contains A.z
   //---------------------------------------------------------------------
   for(j = 0; j < lastcol - firstcol + 1; j++) {
      d = x[j] - r[j];
      sum = sum + d * d;
   }
   clava_dcg_global[ 11 ]++;
   *rnorm = sqrt(sum);
}

//---------------------------------------------------------------------
// generate the test problem for benchmark 6
// makea generates a sparse matrix with a
// prescribed sparsity distribution
//
// parameter    type        usage
//
// input
//
// n            i           number of cols/rows of matrix
// nz           i           nonzeros as declared array size
// rcond        r*8         condition number
// shift        r*8         main diagonal shift
//
// output
//
// a            r*8         array for nonzeros
// colidx       i           col indices
// rowstr       i           row pointers
//
// workspace
//
// iv, arow, acol i
// aelt           r*8
//---------------------------------------------------------------------
void makea(int n, int nz, double a[], int colidx[], int rowstr[], int firstrow, int lastrow, int firstcol, int lastcol, int arow[], int acol[][9], double aelt[][9], int iv[]) {
   int iouter, ivelt, nzv, nn1;
   int ivc[9];
   double vc[9];
   //---------------------------------------------------------------------
   // nonzer is approximately  (int(sqrt(nnza /n)));
   //---------------------------------------------------------------------
   //---------------------------------------------------------------------
   // nn1 is the smallest power of two not less than n
   //---------------------------------------------------------------------
   nn1 = 1;
   do  {
      nn1 = 2 * nn1;
   }
   while (nn1 < n);
   //---------------------------------------------------------------------
   // Generate nonzero positions and save for the use in sparse.
   //---------------------------------------------------------------------
   for(iouter = 0; iouter < n; iouter++) {
      nzv = 8;
      clava_dcg_global[ 12 ]++;
      sprnvc(n, nzv, nn1, vc, ivc);
      clava_dcg_global[ 13 ]++;
      vecset(n, vc, ivc, &nzv, iouter + 1, 0.5);
      arow[iouter] = nzv;
      for(ivelt = 0; ivelt < nzv; ivelt++) {
         acol[iouter][ivelt] = ivc[ivelt] - 1;
         aelt[iouter][ivelt] = vc[ivelt];
      }
   }
   //---------------------------------------------------------------------
   // ... make the sparse matrix from list of elements with duplicates
   //     (iv is used as  workspace)
   //---------------------------------------------------------------------
   clava_dcg_global[ 14 ]++;
   sparse(a, colidx, rowstr, n, nz, 8, arow, acol, aelt, firstrow, lastrow, iv, 1.0e-1, 12.0);
}

//---------------------------------------------------------------------
// rows range from firstrow to lastrow
// the rowstr pointers are defined for nrows = lastrow-firstrow+1 values
//---------------------------------------------------------------------
void sparse(double a[], int colidx[], int rowstr[], int n, int nz, int nozer, int arow[], int acol[][9], double aelt[][9], int firstrow, int lastrow, int nzloc[], double rcond, double shift) {
   int nrows;
   //---------------------------------------------------
   // generate a sparse matrix from a list of
   // [col, row, element] tri
   //---------------------------------------------------
   int i, j, j1, j2, nza, k, kk, nzrow, jcol;
   double size, scale, ratio, va;
   int cont40;
   //---------------------------------------------------------------------
   // how many rows of result
   //---------------------------------------------------------------------
   nrows = lastrow - firstrow + 1;
   //---------------------------------------------------------------------
   // ...count the number of triples in each row
   //---------------------------------------------------------------------
   for(j = 0; j < nrows + 1; j++) {
      rowstr[j] = 0;
   }
   for(i = 0; i < n; i++) {
      for(nza = 0; nza < arow[i]; nza++) {
         j = acol[i][nza] + 1;
         rowstr[j] = rowstr[j] + arow[i];
      }
   }
   rowstr[0] = 0;
   for(j = 1; j < nrows + 1; j++) {
      rowstr[j] = rowstr[j] + rowstr[j - 1];
   }
   nza = rowstr[nrows] - 1;
   //---------------------------------------------------------------------
   // ... rowstr(j) now is the location of the first nonzero
   //     of row j of a
   //---------------------------------------------------------------------
   if(nza > nz) {
      clava_dcg_global[ 15 ]++;
      printf("Space for matrix elements exceeded in sparse\n");
      clava_dcg_global[ 15 ]++;
      printf("nza, nzmax = %d, %d\n", nza, nz);
      clava_dcg_global[ 16 ]++;
      exit(1);
   }
   //---------------------------------------------------------------------
   // ... preload data pages
   //---------------------------------------------------------------------
   for(j = 0; j < nrows; j++) {
      for(k = rowstr[j]; k < rowstr[j + 1]; k++) {
         a[k] = 0.0;
         colidx[k] = -1;
      }
      nzloc[j] = 0;
   }
   //---------------------------------------------------------------------
   // ... generate actual values by summing duplicates
   //---------------------------------------------------------------------
   size = 1.0;
   clava_dcg_global[ 17 ]++;
   ratio = pow(rcond, (1.0 / (double) (n)));
   for(i = 0; i < n; i++) {
      for(nza = 0; nza < arow[i]; nza++) {
         j = acol[i][nza];
         scale = size * aelt[i][nza];
         for(nzrow = 0; nzrow < arow[i]; nzrow++) {
            jcol = acol[i][nzrow];
            va = aelt[i][nzrow] * scale;
            //--------------------------------------------------------------------
            // ... add the identity * rcond to the generated matrix to bound
            //     the smallest eigenvalue from below by rcond
            //--------------------------------------------------------------------
            if(jcol == j && j == i) {
               va = va + rcond - shift;
            }
            cont40 = 0;
            for(k = rowstr[j]; k < rowstr[j + 1]; k++) {
               if(colidx[k] > jcol) {
                  //----------------------------------------------------------------
                  // ... insert colidx here orderly
                  //----------------------------------------------------------------
                  for(kk = rowstr[j + 1] - 2; kk >= k; kk--) {
                     if(colidx[kk] > -1) {
                        a[kk + 1] = a[kk];
                        colidx[kk + 1] = colidx[kk];
                     }
                  }
                  colidx[k] = jcol;
                  a[k] = 0.0;
                  cont40 = 1;
                  break;
               }
               else if(colidx[k] == -1) {
                  colidx[k] = jcol;
                  cont40 = 1;
                  break;
               }
               else if(colidx[k] == jcol) {
                  //--------------------------------------------------------------
                  // ... mark the duplicated entry
                  //--------------------------------------------------------------
                  nzloc[j] = nzloc[j] + 1;
                  cont40 = 1;
                  break;
               }
            }
            if(cont40 == 0) {
               clava_dcg_global[ 15 ]++;
               printf("internal error in sparse: i=%d\n", i);
               clava_dcg_global[ 16 ]++;
               exit(1);
            }
            a[k] = a[k] + va;
         }
      }
      size = size * ratio;
   }
   //---------------------------------------------------------------------
   // ... remove empty entries and generate final results
   //---------------------------------------------------------------------
   for(j = 1; j < nrows; j++) {
      nzloc[j] = nzloc[j] + nzloc[j - 1];
   }
   for(j = 0; j < nrows; j++) {
      if(j > 0) {
         j1 = rowstr[j] - nzloc[j - 1];
      }
      else {
         j1 = 0;
      }
      j2 = rowstr[j + 1] - nzloc[j];
      nza = rowstr[j];
      for(k = j1; k < j2; k++) {
         a[k] = a[nza];
         colidx[k] = colidx[nza];
         nza = nza + 1;
      }
   }
   for(j = 1; j < nrows + 1; j++) {
      rowstr[j] = rowstr[j] - nzloc[j - 1];
   }
   nza = rowstr[nrows] - 1;
}

//---------------------------------------------------------------------
// generate a sparse n-vector (v, iv)
// having nzv nonzeros
//
// mark(i) is set to 1 if position i is nonzero.
// mark is all zero on entry and is reset to all zero before exit
// this corrects a performance bug found by John G. Lewis, caused by
// reinitialization of mark on every one of the n calls to sprnvc
//---------------------------------------------------------------------
void sprnvc(int n, int nz, int nn1, double v[], int iv[]) {
   int nzv, ii, i;
   double vecelt, vecloc;
   nzv = 0;
   while(nzv < nz) {
      clava_dcg_global[ 18 ]++;
      vecelt = randlc(&tran, amult);
      //---------------------------------------------------------------------
      // generate an integer between 1 and n in a portable manner
      //---------------------------------------------------------------------
      clava_dcg_global[ 18 ]++;
      vecloc = randlc(&tran, amult);
      clava_dcg_global[ 19 ]++;
      i = icnvrt(vecloc, nn1) + 1;
      if(i > n) continue;
      //---------------------------------------------------------------------
      // was this integer generated already?
      //---------------------------------------------------------------------
      int was_gen = 0;
      for(ii = 0; ii < nzv; ii++) {
         if(iv[ii] == i) {
            was_gen = 1;
            break;
         }
      }
      if(was_gen) continue;
      v[nzv] = vecelt;
      iv[nzv] = i;
      nzv = nzv + 1;
   }
}

//---------------------------------------------------------------------
// scale a double precision number x in (0,1) by a power of 2 and chop it
//---------------------------------------------------------------------
int icnvrt(double x, int ipwr2) {
   
   return (int) (ipwr2 * x);
}

//---------------------------------------------------------------------
// set ith element of sparse vector (v, iv) with
// nzv nonzeros to val
//---------------------------------------------------------------------
void vecset(int n, double v[], int iv[], int *nzv, int i, double val) {
   int k;
   int set;
   set = 0;
   for(k = 0; k < *nzv; k++) {
      if(iv[k] == i) {
         v[k] = val;
         set = 1;
      }
   }
   if(set == 0) {
      v[*nzv] = val;
      iv[*nzv] = i;
      *nzv = *nzv + 1;
   }
}

double randlc(double *x, double a) {
   //--------------------------------------------------------------------
   //
   //  This routine returns a uniform pseudorandom double precision number in the
   //  range (0, 1) by using the linear congruential generator
   //
   //  x_{k+1} = a x_k  (mod 2^46)
   //
   //  where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
   //  before repeating.  The argument A is the same as 'a' in the above formula,
   //  and X is the same as x_0.  A and X must be odd double precision integers
   //  in the range (1, 2^46).  The returned value RANDLC is normalized to be
   //  between 0 and 1, i.e. RANDLC = 2^(-46) * x_1.  X is updated to contain
   //  the new seed x_1, so that subsequent calls to RANDLC using the same
   //  arguments will generate a continuous sequence.
   //
   //  This routine should produce the same results on any computer with at least
   //  48 mantissa bits in double precision floating point data.  On 64 bit
   //  systems, double precision should be disabled.
   //
   //  David H. Bailey     October 26, 1990
   //
   //--------------------------------------------------------------------
   // r23 = pow(0.5, 23.0);
   ////  pow(0.5, 23.0) = 1.1920928955078125e-07
   // r46 = r23 * r23;
   // t23 = pow(2.0, 23.0);
   ////  pow(2.0, 23.0) = 8.388608e+06
   // t46 = t23 * t23;
   double const r23 = 1.1920928955078125e-07;
   double const r46 = r23 * r23;
   double const t23 = 8.388608e+06;
   double const t46 = t23 * t23;
   double t1, t2, t3, t4, a1, a2, x1, x2, z;
   double r;
   //--------------------------------------------------------------------
   //  Break A into two parts such that A = 2^23 * A1 + A2.
   //--------------------------------------------------------------------
   t1 = r23 * a;
   a1 = (int) t1;
   a2 = a - t23 * a1;
   //--------------------------------------------------------------------
   //  Break X into two parts such that X = 2^23 * X1 + X2, compute
   //  Z = A1 * X2 + A2 * X1  (mod 2^23), and then
   //  X = 2^23 * Z + A2 * X2  (mod 2^46).
   //--------------------------------------------------------------------
   t1 = r23 * (*x);
   x1 = (int) t1;
   x2 = *x - t23 * x1;
   t1 = a1 * x2 + a2 * x1;
   t2 = (int) (r23 * t1);
   z = t1 - t23 * t2;
   t3 = t23 * z + a2 * x2;
   t4 = (int) (r46 * t3);
   *x = t3 - t46 * t4;
   r = r46 * (*x);
   
   return r;
}

void vranlc(int n, double *x, double a, double y[]) {
   //--------------------------------------------------------------------
   //
   //  This routine generates N uniform pseudorandom double precision numbers in
   //  the range (0, 1) by using the linear congruential generator
   //
   //  x_{k+1} = a x_k  (mod 2^46)
   //
   //  where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
   //  before repeating.  The argument A is the same as 'a' in the above formula,
   //  and X is the same as x_0.  A and X must be odd double precision integers
   //  in the range (1, 2^46).  The N results are placed in Y and are normalized
   //  to be between 0 and 1.  X is updated to contain the new seed, so that
   //  subsequent calls to VRANLC using the same arguments will generate a
   //  continuous sequence.  If N is zero, only initialization is performed, and
   //  the variables X, A and Y are ignored.
   //
   //  This routine is the standard version designed for scalar or RISC systems.
   //  However, it should produce the same results on any single processor
   //  computer with at least 48 mantissa bits in double precision floating point
   //  data.  On 64 bit systems, double precision should be disabled.
   //
   //--------------------------------------------------------------------
   // r23 = pow(0.5, 23.0);
   ////  pow(0.5, 23.0) = 1.1920928955078125e-07
   // r46 = r23 * r23;
   // t23 = pow(2.0, 23.0);
   ////  pow(2.0, 23.0) = 8.388608e+06
   // t46 = t23 * t23;
   double const r23 = 1.1920928955078125e-07;
   double const r46 = r23 * r23;
   double const t23 = 8.388608e+06;
   double const t46 = t23 * t23;
   double t1, t2, t3, t4, a1, a2, x1, x2, z;
   int i;
   //--------------------------------------------------------------------
   //  Break A into two parts such that A = 2^23 * A1 + A2.
   //--------------------------------------------------------------------
   t1 = r23 * a;
   a1 = (int) t1;
   a2 = a - t23 * a1;
   //--------------------------------------------------------------------
   //  Generate N results.   This loop is not vectorizable.
   //--------------------------------------------------------------------
   for(i = 0; i < n; i++) {
      //--------------------------------------------------------------------
      //  Break X into two parts such that X = 2^23 * X1 + X2, compute
      //  Z = A1 * X2 + A2 * X1  (mod 2^23), and then
      //  X = 2^23 * Z + A2 * X2  (mod 2^46).
      //--------------------------------------------------------------------
      t1 = r23 * (*x);
      x1 = (int) t1;
      x2 = *x - t23 * x1;
      t1 = a1 * x2 + a2 * x1;
      t2 = (int) (r23 * t1);
      z = t1 - t23 * t2;
      t3 = t23 * z + a2 * x2;
      t4 = (int) (r46 * t3);
      *x = t3 - t46 * t4;
      y[i] = r46 * (*x);
   }
   
   return;
}

void wtime(double *t) {
   static int sec = -1;
   struct timeval tv;
   clava_dcg_global[ 20 ]++;
   gettimeofday(&tv, (void *) 0);
   if(sec < 0) sec = tv.tv_sec;
   *t = (tv.tv_sec - sec) + 1.0e-6 * tv.tv_usec;
}

/*****************************************************************/
/******         E  L  A  P  S  E  D  _  T  I  M  E          ******/
/*****************************************************************/
double elapsed_time() {
   double t;
   clava_dcg_global[ 21 ]++;
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
   clava_dcg_global[ 22 ]++;
   start[n] = elapsed_time();
}

/*****************************************************************/
/******            T  I  M  E  R  _  S  T  O  P             ******/
/*****************************************************************/
void timer_stop(int n) {
   double t, now;
   clava_dcg_global[ 23 ]++;
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

void print_results(char *name, char class, int n1, int n2, int n3, int niter, double t, double mops, char *optype, int verified) {
   char size[16];
   int j;
   clava_dcg_global[ 24 ]++;
   printf("\n\n %s Benchmark Completed.\n", name);
   clava_dcg_global[ 24 ]++;
   printf(" Class           =             %12c\n", class);
   // If this is not a grid-based problem (EP, FT, CG), then
   // we only print n1, which contains some measure of the
   // problem size. In that case, n2 and n3 are both zero.
   // Otherwise, we print the grid size n1xn2xn3
   if((n2 == 0) && (n3 == 0)) {
      if((name[0] == 'E') && (name[1] == 'P')) {
         clava_dcg_global[ 25 ]++;
         clava_dcg_global[ 26 ]++;
         sprintf(size, "%15.0lf", pow(2.0, n1));
         j = 14;
         if(size[j] == '.') {
            size[j] = ' ';
            j--;
         }
         size[j + 1] = '\0';
         clava_dcg_global[ 24 ]++;
         printf(" Size            =          %15s\n", size);
      }
      else {
         clava_dcg_global[ 24 ]++;
         printf(" Size            =             %12d\n", n1);
      }
   }
   else {
      clava_dcg_global[ 24 ]++;
      printf(" Size            =           %4dx%4dx%4d\n", n1, n2, n3);
   }
   clava_dcg_global[ 24 ]++;
   printf(" Iterations      =             %12d\n", niter);
   clava_dcg_global[ 24 ]++;
   printf(" Time in seconds =             %12.2lf\n", t);
   clava_dcg_global[ 24 ]++;
   printf(" Mop/s total     =          %15.2lf\n", mops);
   clava_dcg_global[ 24 ]++;
   printf(" Operation type  = %24s\n", optype);
   if(verified) {
      clava_dcg_global[ 24 ]++;
      printf(" Verification    =             %12s\n", "SUCCESSFUL");
   }
   else {
      clava_dcg_global[ 24 ]++;
      printf(" Verification    =             %12s\n", "UNSUCCESSFUL");
   }
}
