#include <iostream>
#include "../dynare_cpp_driver.hh"
#include <petscdm.h>
#include <petscdmda.h>
#include <petscsnes.h>
#include "main_tao.h"

DynareInfo *preprocessorOutput(void);
PetscErrorCode FormFunction(SNES snes, Vec y, Vec f, void *ctx);
PetscErrorCode FormJacobian(SNES snes, Vec y, Mat J, Mat B, void *ctx);
PetscErrorCode Monitor(SNES, PetscInt, PetscReal, void *);
PetscErrorCode PreCheck(SNESLineSearch, Vec, Vec, PetscBool *, void *);
PetscErrorCode PostCheck(SNESLineSearch, Vec, Vec, Vec, PetscBool *, PetscBool *, void *);
PetscErrorCode PostSetSubKSP(SNESLineSearch, Vec, Vec, Vec, PetscBool *, PetscBool *, void *);
PetscErrorCode MatrixFreePreconditioner(PC, Vec, Vec);

static char  help[] = "Testing";

/*
  User-defined context for monitoring
*/
typedef struct
{
  PetscViewer viewer;
} MonitorCtx;

/*
  User-defined context for checking candidate iterates that are
  determined by line search methods
*/
typedef struct
{
  Vec            last_step;  /* previous iterate */
  PetscReal      tolerance;  /* tolerance for changes between successive iterates */
  AppCtx *user;
} StepCheckCtx;

typedef struct
{
  PetscInt its0; /* num of prevous outer KSP iterations */
} SetSubKSPCtx;

#undef __FUNCT__
#define __FUNCT__ "main"
int
main(int argc, char **argv)
{
  DynareInfo model_info;

  int endo_nbr = model_info.get_endo_nbr();
  int exo_nbr = model_info.get_exo_nbr();
  double *params = model_info.get_params_data();
  // Steady state
  double *steady_state = new double[endo_nbr];
  int info;
  steadystate(NULL, params, steady_state, &info);
  vector<int> NNZDerivatives = model_info.get_NNZDerivatives();
  AppCtx         user;            /* user-defined work context */
  user.exogenous = new double[exo_nbr];
  user.params = params;
  user.steady_state = steady_state;
  user.periods = 4000;
  user.first_col = endo_nbr;
  user.endo_nbr = endo_nbr;
  user.exo_nbr = exo_nbr;
  user.row_ptr = new int[endo_nbr+1];
  user.nnz = NNZDerivatives[0];
  user.col_ptr = new int[NNZDerivatives[0]];
  user.val_ptr = new double[NNZDerivatives[0]];
  user.initial_values = new double[user.endo_nbr];
  user.terminal_values = new double[user.endo_nbr];
  for (int i = 0; i < user.endo_nbr; ++i)
    {
      user.initial_values[i] = user.steady_state[i];
      user.terminal_values[i] = user.steady_state[i];
    }
  /* Initialize PETSc */
  PetscInitialize(&argc, &argv, (char *) 0, help);
  PetscErrorCode ierr;
  /* Get number of processors */
  PetscInt size;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
  /* Set periods a multiple of processor nbr */
  user.periods = size*ceil((double) user.periods/size);
  user.nb_row_x = user.periods + 1;
  user.X = new double[user.nb_row_x*user.exo_nbr];
  for (double *px = user.X; px < user.X+user.nb_row_x*user.exo_nbr; px++)
    *px = 0.0;
  user.X[1] = 0.01;
  user.X[1+user.nb_row_x] = 0.01;
  std::cout << size << " " << user.periods << " " << std::endl;
  PetscInt N = user.periods*user.endo_nbr;
  /* Create DMA */
  DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, N, 1, user.endo_nbr, NULL, &user.da);

  /*     Allocate vector and Jacobian matrix  */ \
  Vec Y, R;
  DMCreateGlobalVector(user.da, &Y);
  VecDuplicate(Y, &R);

  Mat J;
  ierr = MatCreate(PETSC_COMM_WORLD, &J); CHKERRQ(ierr);
  ierr = MatSetSizes(J, PETSC_DECIDE, PETSC_DECIDE, N, N); CHKERRQ(ierr);
  ierr = MatSetFromOptions(J); CHKERRQ(ierr);
  ierr = MatSetUp(J); CHKERRQ(ierr);

  /*
    Get local grid boundaries (for 1-dimensional DMDA):
    xs, xm - starting grid index, width of local grid (no ghost points)
  */
  PetscInt xs, xm, xsg, xmg;
  DMDAGetCorners(user.da, &xs, NULL, NULL, &xm, NULL, NULL);
  std::cout << xs << " " << xm << std::endl;
  DMDAGetGhostCorners(user.da, &xsg, NULL, NULL, &xmg, NULL, NULL);
  std::cout << "ghost " << xsg << " " << xmg << std::endl;
  /*
    Get pointers to vector data
  */
  PetscScalar *YY;
  DMDAVecGetArray(user.da, Y, &YY);

  /*
    Compute local vector entries
  */
  for (int i = xs; i < xs+xm; i += user.endo_nbr)
    for (int j = 0; j < user.endo_nbr; j++)
      YY[i+j] = steady_state[j];
  /*
    Restore vectors
  */
  DMDAVecRestoreArray(user.da, Y, &YY);

  SNES snes;
  ierr = SNESCreate(PETSC_COMM_WORLD, &snes); CHKERRQ(ierr);

  //  SNES snes;
  SNESLineSearch linesearch;          /* SNESLineSearch context */
  MonitorCtx     monP;                 /* monitoring context */
  StepCheckCtx   checkP;               /* step-checking context */
  SetSubKSPCtx   checkP1;
  PetscBool      pre_check, post_check, post_setsubksp; /* flag indicating whether we're checking candidate iterates */
  PetscReal      abstol, rtol, stol, norm;
  PetscInt       its, maxit, maxf;
  PetscBool      flg;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create nonlinear solver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = SNESCreate(PETSC_COMM_WORLD, &snes); CHKERRQ(ierr);

  /*
    Set function evaluation routine and vector.  Whenever the nonlinear
    solver needs to compute the nonlinear function, it will call this
    routine.
    - Note that the final routine argument is the user-defined
    context that provides application-specific data for the
    function evaluation routine.
  */
  ierr = SNESSetFunction(snes, R, FormFunction, &user); CHKERRQ(ierr);

  /*
    Set Jacobian matrix data structure and default Jacobian evaluation
    routine.  Whenever the nonlinear solver needs to compute the
    Jacobian matrix, it will call this routine.
    - Note that the final routine argument is the user-defined
    context that provides application-specific data for the
    Jacobian evaluation routine.
  */
  ierr = SNESSetJacobian(snes, J, J, FormJacobian, &user); CHKERRQ(ierr);

  /*
    Optional allow user provided preconditioner
  */
  ierr = PetscOptionsHasName(NULL, NULL, "-user_precond", &flg); CHKERRQ(ierr);
  if (flg)
    {
      KSP ksp;
      PC  pc;
      ierr = SNESGetKSP(snes, &ksp); CHKERRQ(ierr);
      ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
      ierr = PCSetType(pc, PCSHELL); CHKERRQ(ierr);
      ierr = PCShellSetApply(pc, MatrixFreePreconditioner); CHKERRQ(ierr);
    }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Customize nonlinear solver; set runtime options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
    Set an optional user-defined monitoring routine
  */
  ierr = PetscViewerDrawOpen(PETSC_COMM_WORLD, 0, 0, 0, 0, 400, 400, &monP.viewer); CHKERRQ(ierr);
  ierr = SNESMonitorSet(snes, Monitor, &monP, 0); CHKERRQ(ierr);

  /*
    Set names for some vectors to facilitate monitoring (optional)
  */
  ierr = PetscObjectSetName((PetscObject) Y, "Approximate Solution"); CHKERRQ(ierr);
  //  ierr = PetscObjectSetName((PetscObject)U,"Exact Solution");CHKERRQ(ierr);

  /*
    Set SNES/KSP/KSP/PC runtime options, e.g.,
    -snes_view -snes_monitor -ksp_type <ksp> -pc_type <pc>
  */
  KSP ksp;
  PC pc;
  ierr = SNESGetKSP(snes, &ksp); CHKERRQ(ierr);
  ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, PCLU); CHKERRQ(ierr);
  //  PetscOptionsSetValue(NULL,"-pc_type","lu");
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

  /*
    Set an optional user-defined routine to check the validity of candidate
    iterates that are determined by line search methods
  */
  ierr = SNESGetLineSearch(snes, &linesearch); CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL, NULL, "-post_check_iterates", &post_check); CHKERRQ(ierr);

  if (post_check)
    {
      ierr = PetscPrintf(PETSC_COMM_WORLD, "Activating post step checking routine\n"); CHKERRQ(ierr);
      ierr = SNESLineSearchSetPostCheck(linesearch, PostCheck, &checkP); CHKERRQ(ierr);
      ierr = VecDuplicate(Y, &(checkP.last_step)); CHKERRQ(ierr);

      checkP.tolerance = 1.0;
      checkP.user      = &user;

      ierr = PetscOptionsGetReal(NULL, NULL, "-check_tol", &checkP.tolerance, NULL); CHKERRQ(ierr);
    }

  ierr = PetscOptionsHasName(NULL, NULL, "-post_setsubksp", &post_setsubksp); CHKERRQ(ierr);
  if (post_setsubksp)
    {
      ierr = PetscPrintf(PETSC_COMM_WORLD, "Activating post setsubksp\n"); CHKERRQ(ierr);
      ierr = SNESLineSearchSetPostCheck(linesearch, PostSetSubKSP, &checkP1); CHKERRQ(ierr);
    }

  ierr = PetscOptionsHasName(NULL, NULL, "-pre_check_iterates", &pre_check); CHKERRQ(ierr);
  if (pre_check)
    {
      ierr = PetscPrintf(PETSC_COMM_WORLD, "Activating pre step checking routine\n"); CHKERRQ(ierr);
      ierr = SNESLineSearchSetPreCheck(linesearch, PreCheck, &checkP); CHKERRQ(ierr);
    }

  /*
    Print parameters used for convergence testing (optional) ... just
    to demonstrate this routine; this information is also printed with
    the option -snes_view
  */
  ierr = SNESGetTolerances(snes, &abstol, &rtol, &stol, &maxit, &maxf); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "atol=%g, rtol=%g, stol=%g, maxit=%D, maxf=%D\n", (double) abstol, (double) rtol, (double) stol, maxit, maxf); CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Evaluate initial guess; then solve nonlinear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
    Note: The user should initialize the vector, x, with the initial guess
    for the nonlinear solver prior to calling SNESSolve().  In particular,
    to employ an initial guess of zero, the user should explicitly set
    this vector to zero by calling VecSet().
  */
  ierr = SNESSolve(snes, NULL, Y); CHKERRQ(ierr);
  ierr = SNESGetIterationNumber(snes, &its); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Number of SNES iterations = %D\n", its); CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Check solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
    Free work space.  All PETSc objects should be destroyed when they
    are no longer needed.
  */
  ierr = PetscViewerDestroy(&monP.viewer); CHKERRQ(ierr);
  if (post_check)
    {
      ierr = VecDestroy(&checkP.last_step); CHKERRQ(ierr);
    }
  ierr = SNESDestroy(&snes); CHKERRQ(ierr);
  ierr = DMDestroy(&user.da); CHKERRQ(ierr);

  ierr = MatDestroy(&J); CHKERRQ(ierr);
  ierr = VecDestroy(&Y); CHKERRQ(ierr);
  ierr = VecDestroy(&R); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode
FormFunction(SNES snes, Vec y, Vec f, void *ctx)
{
  AppCtx *user = (AppCtx *) ctx;
  DM             da    = user->da;
  PetscScalar    *yy, *ff;
  PetscInt       M, ys, ym;
  Vec            ylocal;

  DMGetLocalVector(da, &ylocal);
  /*
    Scatter ghost points to local vector, using the 2-step process
    DMGlobalToLocalBegin(), DMGlobalToLocalEnd().
    By placing code between these two statements, computations can
    be done while messages are in transition.
  */
  DMGlobalToLocalBegin(da, y, INSERT_VALUES, ylocal);
  DMGlobalToLocalEnd(da, y, INSERT_VALUES, ylocal);

  /*
    Get pointers to vector data.
    - The vector xlocal includes ghost point; the vectors x and f do
    NOT include ghost points.
    - Using DMDAVecGetArray() allows accessing the values using global ordering
  */
  DMDAVecGetArray(da, ylocal, &yy);
  DMDAVecGetArray(da, f, &ff);

  /*
    Get local grid boundaries (for 1-dimensional DMDA):
    ys, ym  - starting grid index, width of local grid (no ghost points)
  */
  DMDAGetCorners(da, &ys, NULL, NULL, &ym, NULL, NULL);
  DMDAGetInfo(da, NULL, &M, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

  /*
    Set function values for boundary points; define local interior grid point range:
    xsi - starting interior grid index
    xei - ending interior grid index
  */
  if (ys == 0) /* left boundary */
    {
      PetscReal *y1 = new PetscReal[3*user->endo_nbr];
      for (int i = 0; i < user->endo_nbr; ++i)
        y1[i] = user->initial_values[i];
      for (int i = 0; i < 2*user->endo_nbr; ++i)
        y1[i+user->endo_nbr] = yy[i];
      Residuals(y1, user->X, user->nb_row_x, user->params, user->steady_state, 1, ff);
      ys += user->endo_nbr;
      ym -= user->endo_nbr;
    }

  /*
    Compute function over locally owned part of the grid (interior points only)
  */
  while ((ym  >= user->endo_nbr) && (ys + 2*user->endo_nbr <= M))
    {
      int it = ys/user->endo_nbr + 2;
      Residuals(yy+ys-user->endo_nbr, user->X, user->nb_row_x, user->params, user->steady_state, it, ff+ys);
      ys += user->endo_nbr;
      ym -= user->endo_nbr;
    }

  if ((ym >= user->endo_nbr) && (ys + 2*user->endo_nbr >= M))
    {
      int it = ys/user->endo_nbr + 1;
      PetscReal *y1 = new PetscReal[3*user->endo_nbr];
      for (int i = 0; i < 2*user->endo_nbr; ++i)
        y1[i] = yy[ys+i-user->endo_nbr];
      for (int i = 0; i < user->endo_nbr; ++i)
        y1[i+2*user->endo_nbr] = user->terminal_values[i];
      Residuals(y1, user->X, user->nb_row_x, user->params, user->steady_state, it, ff+ys);
    }

  /*
    Restore vectors
  */
  DMDAVecRestoreArray(da, ylocal, &yy);
  DMDAVecRestoreArray(da, f, &ff);
  DMRestoreLocalVector(da, &ylocal);
  return (0);
}

PetscErrorCode
FormJacobian(SNES snes, Vec y, Mat J, Mat B, void *ctx)
{
  AppCtx *user = (AppCtx *) ctx;
  DM             da    = user->da;
  PetscScalar    *yy;
  PetscInt       M, ys, ym, ierr;
  Vec            ylocal;

  DMGetLocalVector(da, &ylocal);
  /*
    Scatter ghost points to local vector, using the 2-step process
    DMGlobalToLocalBegin(), DMGlobalToLocalEnd().
    By placing code between these two statements, computations can
    be done while messages are in transition.
  */
  DMGlobalToLocalBegin(da, y, INSERT_VALUES, ylocal);
  DMGlobalToLocalEnd(da, y, INSERT_VALUES, ylocal);

  /*
    Get pointers to vector data.
    - The vector xlocal includes ghost point; the vectors x and f do
    NOT include ghost points.
    - Using DMDAVecGetArray() allows accessing the values using global ordering
  */
  DMDAVecGetArray(da, ylocal, &yy);

  /*
    Get local grid boundaries (for 1-dimensional DMDA):
    ys, ym  - starting grid index, width of local grid (no ghost points)
  */
  DMDAGetCorners(da, &ys, NULL, NULL, &ym, NULL, NULL);
  DMDAGetInfo(da, NULL, &M, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

  /*
    Set function values for boundary points; define local interior grid point range:
    xsi - starting interior grid index
    xei - ending interior grid index
  */
  int row = ys;
  if (ys == 0) /* left boundary */
    {
      PetscReal *y1 = new PetscReal[3*user->endo_nbr];
      for (int i = 0; i < user->endo_nbr; ++i)
        y1[i] = user->initial_values[i];
      for (int i = 0; i < 2*user->endo_nbr; ++i)
        y1[i+user->endo_nbr] = yy[i];
      FirstDerivatives(y1, user->X, user->nb_row_x, user->params, user->steady_state, 1, NULL, user->row_ptr, user->col_ptr, user->val_ptr);
      for (int *r = user->row_ptr; r < user->row_ptr+user->endo_nbr; r++)
        {
          int first_col = 0;
          int ncol = 0;
          int *pc = user->col_ptr + *r;
          while (*(pc) < user->endo_nbr)
            {
              ++first_col;
              ++pc;
            }
          while (pc < ((user->col_ptr)+*(r+1)))
            {
              if (*pc < 3*(user->endo_nbr))
                {
                  ++ncol;
                  *pc -= user->endo_nbr;
                }
              ++pc;
            }
          ierr = MatSetValues(J, 1, &row, ncol, user->col_ptr + *r + first_col, user->val_ptr + *r  + first_col, INSERT_VALUES);
          CHKERRQ(ierr);
          ++row;
        }
      ys += user->endo_nbr;
      ym -= user->endo_nbr;
    }

  /*
    Compute function over locally owned part of the grid (interior points only)
  */
  while ((ym  >= user->endo_nbr) && (ys + 2*user->endo_nbr <= M))
    {
      int it = ys/user->endo_nbr + 2;
      FirstDerivatives(yy+ys-user->endo_nbr, user->X, user->nb_row_x, user->params, user->steady_state, it, NULL, user->row_ptr, user->col_ptr, user->val_ptr);
      for (int *r = user->row_ptr; r < user->row_ptr+user->endo_nbr; r++)
        {
          int ncol = 0;
          for (int *pc = user->col_ptr + *r; pc < user->col_ptr + *(r+1); ++pc)
            if (*pc < 3*user->endo_nbr)
              {
                *pc += ys - user->endo_nbr;
                ++ncol;
              }
          ierr = MatSetValues(J, 1, &row, ncol, user->col_ptr + *r, user->val_ptr + *r, INSERT_VALUES);
          CHKERRQ(ierr);
          ++row;
        }
      ys += user->endo_nbr;
      ym -= user->endo_nbr;
    }

  if ((ym >= user->endo_nbr) && (ys + 2*user->endo_nbr >= M))
    {
      int it = ys/user->endo_nbr + 1;
      PetscReal *y1 = new PetscReal[3*user->endo_nbr];
      for (int i = 0; i < 2*user->endo_nbr; ++i)
        y1[i] = yy[ys+i-user->endo_nbr];
      for (int i = 0; i < user->endo_nbr; ++i)
        y1[i+2*user->endo_nbr] = user->terminal_values[i];
      FirstDerivatives(y1, user->X, user->nb_row_x, user->params, user->steady_state, it, NULL, user->row_ptr, user->col_ptr, user->val_ptr);
      for (int *r = user->row_ptr; r < user->row_ptr+user->endo_nbr; r++)
        {
          int *pc = user->col_ptr + *r;
          int ncol = 0;
          while ((*pc < 2*user->endo_nbr) && (pc < user->col_ptr + *(r+1)))
            {
              ++ncol;
              *pc += ys - user->endo_nbr;
              ++pc;
            }
          ierr = MatSetValues(J, 1, &row, ncol, user->col_ptr + *r, user->val_ptr + *r, INSERT_VALUES);
          CHKERRQ(ierr);
          ++row;
        }
    }

  /*
    Restore vectors
  */
  ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  DMDAVecRestoreArray(da, ylocal, &yy);
  DMRestoreLocalVector(da, &ylocal);
  ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  return (0);
}

#undef __FUNCT__
#define __FUNCT__ "Monitor"
/*
  Monitor - Optional user-defined monitoring routine that views the
  current iterate with an x-window plot. Set by SNESMonitorSet().

  Input Parameters:
  snes - the SNES context
  its - iteration number
  norm - 2-norm function value (may be estimated)
  ctx - optional user-defined context for private data for the
  monitor routine, as set by SNESMonitorSet()

  Note:
  See the manpage for PetscViewerDrawOpen() for useful runtime options,
  such as -nox to deactivate all x-window output.
*/
PetscErrorCode
Monitor(SNES snes, PetscInt its, PetscReal fnorm, void *ctx)
{
  PetscErrorCode ierr;
  MonitorCtx     *monP = (MonitorCtx *) ctx;
  Vec            x;

  PetscFunctionBeginUser;
  ierr = PetscPrintf(PETSC_COMM_WORLD, "iter = %D,SNES Function norm %g\n", its, (double) fnorm); CHKERRQ(ierr);
  ierr = SNESGetSolution(snes, &x); CHKERRQ(ierr);
  ierr = VecView(x, monP->viewer); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "PreCheck"
/*
  PreCheck - Optional user-defined routine that checks the validity of
  candidate steps of a line search method.  Set by SNESLineSearchSetPreCheck().

  Input Parameters:
  snes - the SNES context
  xcurrent - current solution
  y - search direction and length

  Output Parameters:
  y         - proposed step (search direction and length) (possibly changed)
  changed_y - tells if the step has changed or not
*/
PetscErrorCode
PreCheck(SNESLineSearch linesearch, Vec xcurrent, Vec y, PetscBool *changed_y, void *ctx)
{
  PetscFunctionBeginUser;
  *changed_y = PETSC_FALSE;
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "PostCheck"
/*
  PostCheck - Optional user-defined routine that checks the validity of
  candidate steps of a line search method.  Set by SNESLineSearchSetPostCheck().

  Input Parameters:
  snes - the SNES context
  ctx  - optional user-defined context for private data for the
  monitor routine, as set by SNESLineSearchSetPostCheck()
  xcurrent - current solution
  y - search direction and length
  x    - the new candidate iterate

  Output Parameters:
  y    - proposed step (search direction and length) (possibly changed)
  x    - current iterate (possibly modified)

*/
PetscErrorCode
PostCheck(SNESLineSearch linesearch, Vec xcurrent, Vec y, Vec x, PetscBool  *changed_y, PetscBool  *changed_x, void *ctx)
{
  PetscErrorCode ierr;
  PetscInt       i, iter, xs, xm;
  StepCheckCtx   *check;
  AppCtx *user;
  PetscScalar    *xa, *xa_last, tmp;
  PetscReal      rdiff;
  DM             da;
  SNES           snes;

  PetscFunctionBeginUser;
  *changed_x = PETSC_FALSE;
  *changed_y = PETSC_FALSE;

  ierr  = SNESLineSearchGetSNES(linesearch, &snes); CHKERRQ(ierr);
  check = (StepCheckCtx *) ctx;
  user  = check->user;
  ierr  = SNESGetIterationNumber(snes, &iter); CHKERRQ(ierr);
  ierr  = SNESLineSearchGetPreCheck(linesearch, NULL, (void **) &check); CHKERRQ(ierr);

  /* iteration 1 indicates we are working on the second iteration */
  if (iter > 0)
    {
      da   = user->da;
      ierr = PetscPrintf(PETSC_COMM_WORLD, "Checking candidate step at iteration %D with tolerance %g\n", iter, (double) check->tolerance); CHKERRQ(ierr);

      /* Access local array data */
      ierr = DMDAVecGetArray(da, check->last_step, &xa_last); CHKERRQ(ierr);
      ierr = DMDAVecGetArray(da, x, &xa); CHKERRQ(ierr);
      ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL); CHKERRQ(ierr);

      /*
        If we fail the user-defined check for validity of the candidate iterate,
        then modify the iterate as we like.  (Note that the particular modification
        below is intended simply to demonstrate how to manipulate this data, not
        as a meaningful or appropriate choice.)
      */
      for (i = xs; i < xs+xm; i++)
        {
          if (!PetscAbsScalar(xa[i]))
            rdiff = 2*check->tolerance;
          else
            rdiff = PetscAbsScalar((xa[i] - xa_last[i])/xa[i]);
          if (rdiff > check->tolerance)
            {
              tmp        = xa[i];
              xa[i]      = .5*(xa[i] + xa_last[i]);
              *changed_x = PETSC_TRUE;
              ierr       = PetscPrintf(PETSC_COMM_WORLD, "  Altering entry %D: x=%g, x_last=%g, diff=%g, x_new=%g\n",
                                       i, (double) PetscAbsScalar(tmp), (double) PetscAbsScalar(xa_last[i]), (double) rdiff, (double) PetscAbsScalar(xa[i])); CHKERRQ(ierr);
            }
        }
      ierr = DMDAVecRestoreArray(da, check->last_step, &xa_last); CHKERRQ(ierr);
      ierr = DMDAVecRestoreArray(da, x, &xa); CHKERRQ(ierr);
    }
  ierr = VecCopy(x, check->last_step); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "PostSetSubKSP"
/*
  PostSetSubKSP - Optional user-defined routine that reset SubKSP options when hierarchical bjacobi PC is used
  e.g,
  mpiexec -n 8 ./ex3 -nox -n 10000 -ksp_type fgmres -pc_type bjacobi -pc_bjacobi_blocks 4 -sub_ksp_type gmres -sub_ksp_max_it 3 -post_setsubksp -sub_ksp_rtol 1.e-16
  Set by SNESLineSearchSetPostCheck().

  Input Parameters:
  linesearch - the LineSearch context
  xcurrent - current solution
  y - search direction and length
  x    - the new candidate iterate

  Output Parameters:
  y    - proposed step (search direction and length) (possibly changed)
  x    - current iterate (possibly modified)

*/
PetscErrorCode
PostSetSubKSP(SNESLineSearch linesearch, Vec xcurrent, Vec y, Vec x, PetscBool  *changed_y, PetscBool  *changed_x, void *ctx)
{
  PetscErrorCode ierr;
  SetSubKSPCtx   *check;
  PetscInt       iter, its, sub_its, maxit;
  KSP            ksp, sub_ksp, *sub_ksps;
  PC             pc;
  PetscReal      ksp_ratio;
  SNES           snes;

  PetscFunctionBeginUser;
  ierr    = SNESLineSearchGetSNES(linesearch, &snes); CHKERRQ(ierr);
  check   = (SetSubKSPCtx *) ctx;
  ierr    = SNESGetIterationNumber(snes, &iter); CHKERRQ(ierr);
  ierr    = SNESGetKSP(snes, &ksp); CHKERRQ(ierr);
  ierr    = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
  ierr    = PCBJacobiGetSubKSP(pc, NULL, NULL, &sub_ksps); CHKERRQ(ierr);
  sub_ksp = sub_ksps[0];
  ierr    = KSPGetIterationNumber(ksp, &its); CHKERRQ(ierr);      /* outer KSP iteration number */
  ierr    = KSPGetIterationNumber(sub_ksp, &sub_its); CHKERRQ(ierr); /* inner KSP iteration number */

  if (iter)
    {
      ierr      = PetscPrintf(PETSC_COMM_WORLD, "    ...PostCheck snes iteration %D, ksp_it %d %d, subksp_it %d\n", iter, check->its0, its, sub_its); CHKERRQ(ierr);
      ksp_ratio = ((PetscReal) (its))/check->its0;
      maxit     = (PetscInt) (ksp_ratio*sub_its + 0.5);
      if (maxit < 2)
        maxit = 2;
      ierr = KSPSetTolerances(sub_ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, maxit); CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD, "    ...ksp_ratio %g, new maxit %d\n\n", ksp_ratio, maxit); CHKERRQ(ierr);
    }
  check->its0 = its; /* save current outer KSP iteration number */
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
/*
  MatrixFreePreconditioner - This routine demonstrates the use of a
  user-provided preconditioner.  This code implements just the null
  preconditioner, which of course is not recommended for general use.

  Input Parameters:
  +  pc - preconditioner
  -  x - input vector

  Output Parameter:
  .  y - preconditioned vector
*/
PetscErrorCode
MatrixFreePreconditioner(PC pc, Vec x, Vec y)
{
  PetscErrorCode ierr;
  ierr = VecCopy(x, y); CHKERRQ(ierr);
  return 0;
}
