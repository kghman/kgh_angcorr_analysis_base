
#define MAX_FOLD 4

typedef struct wiko_setup_s{
  double (*orient_function)(int new_spins,double i1,double sigma,int lam);
  int fold;
  double (*f_function1)(int l1,int l2,double i1,double i2,
			int lam1,int lam,int lam2);
  double (*f_function2)(int l1,int l2,double i1,double i2,
			int lam1,int lam,int lam2);
  double (*f_function3)(int l1,int l2,double i1,double i2,
			int lam1,int lam,int lam2);
  double (*f_function4)(int l1,int l2,double i1,double i2,
			int lam1,int lam,int lam2);
} wiko_setup;
