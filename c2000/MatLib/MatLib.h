#include "stdlib.h"
#include "Math.h"

// https://www.codeproject.com/Articles/5283245/Matrix-Library-in-C

struct Mat_s{
	double* entries;
	int row;
	int col;
};
typedef struct Mat_s Mat;

struct MatList_s{
	Mat* mat;
	struct MatList_s* next;
};
typedef struct MatList_s MatList;

// Matrix initialisations & deletions
void showmat(Mat*);
Mat* newmat(int,int,double);
Mat* from_array(int r,int c,double *d);
void to_array(Mat* M, double *array);
void freemat(Mat* A);
Mat* eye(int n);
Mat* zeros(int r,int c);
Mat* ones(int r,int c);
Mat* randm(int r,int c,double l,double u);
double get(Mat* M,int r,int c);
void set(Mat* M,int r,int c,double d);

// Matrix operations
Mat* scalermultiply(Mat* M,double c);
Mat* sum(Mat* A,Mat* B);
Mat* minus(Mat* A,Mat* B);
Mat* submat(Mat* A,int r1,int r2,int c1,int c2);
void submat2(Mat* A,Mat* B,int r1,int r2,int c1,int c2);
Mat* multiply(Mat* A,Mat* B);
Mat* removerow(Mat* A,int r);
Mat* removecol(Mat* A,int c);
void removerow2(Mat* A,Mat* B,int r);
void removecol2(Mat* A,Mat* B,int c);
Mat* transpose(Mat* A);
double det(Mat* M);
double trace(Mat* A);
Mat* adjoint(Mat* A);
Mat* inverse(Mat* A);
Mat* copyvalue(Mat* A);
Mat* triinverse(Mat* A);
Mat* rowechelon(Mat* A);
Mat* hconcat(Mat* A,Mat* B);
Mat* vconcat(Mat* A,Mat* B);
double norm(Mat* A);
Mat* null(Mat *A);
MatList* lu(Mat* A);
double innermultiply(Mat* a,Mat* b);
MatList* qr(Mat*);




