! @brief Links between the C and Fortran UMFPACK procedures.
! @details partially taken from http://geo.mff.cuni.cz/~lh/Fortran/UMFPACK/README.html
module umfpack_calls
  
  use prec
  use iso_c_binding
  implicit none
  
  ! --------------------------------------------------------------------------------
  ! Data from /usr/include/umfpack.h -- to be adapted to the version of UMFPACK used
  ! --------------------------------------------------------------------------------
  
  ! size of Info and Control arrays
  ! -------------------------------
  integer, parameter :: UMFPACK_INFO    = 90
  integer, parameter :: UMFPACK_CONTROL = 20
  
  ! Version, copyright, and license
  ! -------------------------------
  character(len=*), parameter :: UMFPACK_VERSION   = &
       & "UMFPACK V5.6.2 (Apr 25, 2013)"
  character(len=*), parameter :: UMFPACK_COPYRIGHT = &
       & "UMFPACK:  Copyright (c) 2005-2012 by Timothy A. Davis.  All Rights Reserved."
  character(12),parameter :: UMFPACK_DATE          = &
       & "Apr 25, 2013"
  integer,parameter :: UMFPACK_MAIN_VERSION   = 5
  integer,parameter :: UMFPACK_SUB_VERSION    = 6
  integer,parameter :: UMFPACK_SUBSUB_VERSION = 2
  
  ! contents of Info
  ! ----------------
  enum, bind(c) 
     enumerator :: & 
          ! returned by all routines that use Info:
          UMFPACK_STATUS=0, &                    ! UMFPACK_OK, or other result
          UMFPACK_NROW=1, &                      ! n_row input value
          UMFPACK_NCOL=16, &                     ! n_col input value
          UMFPACK_NZ=2, &                        ! # of entries in A
          ! computed in UMFPACK_*symbolic and UMFPACK_numeric:
          UMFPACK_SIZE_OF_UNIT=3, &              ! sizeof (Unit)
          ! computed in UMFPACK_*symbolic:
          UMFPACK_SIZE_OF_INT=4, &               ! sizeof (int)
          UMFPACK_SIZE_OF_LONG=5, &              ! sizeof (SuiteSparse_long)
          UMFPACK_SIZE_OF_POINTER=6, &           ! sizeof (void *)
          UMFPACK_SIZE_OF_ENTRY=7, &             ! sizeof (Entry), real or complex
          UMFPACK_NDENSE_ROW=8, &                ! number of dense rows
          UMFPACK_NEMPTY_ROW=9, &                ! number of empty rows
          UMFPACK_NDENSE_COL=10, &               ! number of dense rows
          UMFPACK_NEMPTY_COL=11, &               ! number of empty rows
          UMFPACK_SYMBOLIC_DEFRAG=12, &          ! # of memory compactions
          UMFPACK_SYMBOLIC_PEAK_MEMORY=13, &     ! memory used by symbolic analysis
          UMFPACK_SYMBOLIC_SIZE=14, &            ! size of Symbolic object, in Units
          UMFPACK_SYMBOLIC_TIME=15, &            ! time (sec.) for symbolic analysis
          UMFPACK_SYMBOLIC_WALLTIME=17, &        ! wall clock time for sym. analysis
          UMFPACK_STRATEGY_USED=18, &            ! strategy used: sym, unsym
          UMFPACK_ORDERING_USED=19, &            ! ordering used: colamd, amd, given
          UMFPACK_QFIXED=31, &                   ! whether Q is fixed or refined
          UMFPACK_DIAG_PREFERRED=32, &           ! whether diagonal pivoting attempted
          UMFPACK_PATTERN_SYMMETRY=33, &         ! symmetry of pattern of S
          UMFPACK_NZ_A_PLUS_AT=34, &             ! nnz (S+S'), excl. diagonal
          UMFPACK_NZDIAG=35, &                   ! nnz (diag (S))
          ! AMD statistics, computed in UMFPACK_*symbolic:
          UMFPACK_SYMMETRIC_LUNZ=36, &           ! nz in L+U, if AMD ordering used
          UMFPACK_SYMMETRIC_FLOPS=37, &          ! flops for LU, if AMD ordering used
          UMFPACK_SYMMETRIC_NDENSE=38, &         ! # of "dense" rows/cols in S+S'
          UMFPACK_SYMMETRIC_DMAX=39, &           ! max nz in cols of L, for AMD
          ! 51:55 unused
          ! statistics for singleton pruning
          UMFPACK_COL_SINGLETONS=56, &           ! # of column singletons
          UMFPACK_ROW_SINGLETONS=57, &           ! # of row singletons
          UMFPACK_N2=58, &                       ! size of S
          UMFPACK_S_SYMMETRIC=59, &              ! 1 if S square and symmetricly perm.
          ! estimates computed in UMFPACK_*symbolic:
          UMFPACK_NUMERIC_SIZE_ESTIMATE=20, &    ! final size of Numeric->Memory
          UMFPACK_PEAK_MEMORY_ESTIMATE=21, &     ! for symbolic & numeric
          UMFPACK_FLOPS_ESTIMATE=22, &           ! flop count
          UMFPACK_LNZ_ESTIMATE=23, &             ! nz in L, incl. diagonal
          UMFPACK_UNZ_ESTIMATE=24, &             ! nz in U, incl. diagonal
          UMFPACK_VARIABLE_INIT_ESTIMATE=25, &   ! initial size of Numeric->Memory
          UMFPACK_VARIABLE_PEAK_ESTIMATE=26, &   ! peak size of Numeric->Memory
          UMFPACK_VARIABLE_FINAL_ESTIMATE=27, &  ! final size of Numeric->Memory
          UMFPACK_MAX_FRONT_SIZE_ESTIMATE=28, &  ! max frontal matrix size
          UMFPACK_MAX_FRONT_NROWS_ESTIMATE=29, & ! max # rows in any front
          UMFPACK_MAX_FRONT_NCOLS_ESTIMATE=30, & ! max # columns in any front
          ! exact values, (estimates shown above) computed in UMFPACK_numeric:
          UMFPACK_NUMERIC_SIZE=40, &             ! final size of Numeric->Memory
          UMFPACK_PEAK_MEMORY=41, &              ! for symbolic & numeric
          UMFPACK_FLOPS=42, &                    ! flop count
          UMFPACK_LNZ=43, &                      ! nz in L, incl. diagonal
          UMFPACK_UNZ=44, &                      ! nz in U, incl. diagonal
          UMFPACK_VARIABLE_INIT=45, &            ! initial size of Numeric->Memory
          UMFPACK_VARIABLE_PEAK=46, &            ! peak size of Numeric->Memory
          UMFPACK_VARIABLE_FINAL=47, &           ! final size of Numeric->Memory
          UMFPACK_MAX_FRONT_SIZE=48, &           ! max frontal matrix size
          UMFPACK_MAX_FRONT_NROWS=49, &          ! max # rows in any front
          UMFPACK_MAX_FRONT_NCOLS=50, &          ! max # columns in any front
          ! computed in UMFPACK_numeric:
          UMFPACK_NUMERIC_DEFRAG=60, &           ! # of garbage collections
          UMFPACK_NUMERIC_REALLOC=61, &          ! # of memory reallocations
          UMFPACK_NUMERIC_COSTLY_REALLOC=62, &   ! # of costlly memory realloc's
          UMFPACK_COMPRESSED_PATTERN=63, &       ! # of integers in LU pattern
          UMFPACK_LU_ENTRIES=64, &               ! # of reals in LU factors
          UMFPACK_NUMERIC_TIME=65, &             ! numeric factorization time
          UMFPACK_UDIAG_NZ=66, &                 ! nz on diagonal of U
          UMFPACK_RCOND=67, &                    ! est. reciprocal condition #
          UMFPACK_WAS_SCALED=68, &               ! none, max row, or sum row
          UMFPACK_RSMIN=69, &                    ! min (max row) or min (sum row)
          UMFPACK_RSMAX=70, &                    ! max (max row) or max (sum row)
          UMFPACK_UMIN=71, &                     ! min abs diagonal entry of U
          UMFPACK_UMAX=72, &                     ! max abs diagonal entry of U
          UMFPACK_ALLOC_INIT_USED=73, &          ! alloc_init parameter used
          UMFPACK_FORCED_UPDATES=74, &           ! # of forced updates
          UMFPACK_NUMERIC_WALLTIME=75, &         ! numeric wall clock time
          UMFPACK_NOFF_DIAG=76, &                ! number of off-diagonal pivots
          UMFPACK_ALL_LNZ=77, &                  ! nz in L, if no dropped entries
          UMFPACK_ALL_UNZ=78, &                  ! nz in U, if no dropped entries
          UMFPACK_NZDROPPED=79, &                ! # of dropped small entries
          ! computed in UMFPACK_solve:
          UMFPACK_IR_TAKEN=80, &                 ! # of iterative refinement steps taken
          UMFPACK_IR_ATTEMPTED=81, &             ! # of iter. refinement steps attempted
          UMFPACK_OMEGA1=82, &                   ! omega1, sparse backward error estimate
          UMFPACK_OMEGA2=83, &                   ! omega2, sparse backward error estimate
          UMFPACK_SOLVE_FLOPS=84, &              ! flop count for solve
          UMFPACK_SOLVE_TIME=85, &               ! solve time (seconds)
          UMFPACK_SOLVE_WALLTIME=86              ! solve time (wall clock, seconds)
     ! Info(87,88,89) unused
     ! Unused parts of Info may be used in future versions of UMFPACK.
  end enum
  
  ! contents of Control
  ! -------------------
  enum, bind(c) 
     enumerator :: & 
          ! used in all UMFPACK_report_* routines:
          UMFPACK_PRL=0, &                   ! print level
          ! used in UMFPACK_*symbolic only:
          UMFPACK_DENSE_ROW=1, &             ! dense row parameter
          UMFPACK_DENSE_COL=2, &             ! dense col parameter
          UMFPACK_BLOCK_SIZE=4, &            ! BLAS-3 block size
          UMFPACK_STRATEGY=5, &              ! auto, symmetric, or unsym.
          UMFPACK_ORDERING=10, &             ! ordering method to use
          UMFPACK_FIXQ=13, &                 ! -1: no fixQ, 0: default, 1: fixQ
          UMFPACK_AMD_DENSE=14, &            ! for AMD ordering
          UMFPACK_AGGRESSIVE=19, &           ! whether or not to use aggressive
          UMFPACK_SINGLETONS=11, &           ! singleton filter on if true
          ! used in UMFPACK_*numeric only:
          UMFPACK_PIVOT_TOLERANCE=3, &       ! threshold partial pivoting setting
          UMFPACK_ALLOC_INIT=6, &            ! initial allocation ratio
          UMFPACK_SYM_PIVOT_TOLERANCE=15, &  ! threshold, only for diag. entries
          UMFPACK_SCALE=16, &                ! what row scaling to do
          UMFPACK_FRONT_ALLOC_INIT=17, &     ! frontal matrix allocation ratio
          UMFPACK_DROPTOL=18, &              ! drop tolerance for entries in L,U
          ! used in UMFPACK_*solve only:
          UMFPACK_IRSTEP=7, &                ! max # of iterative refinements
          ! compile-time settings - Control(8:11) cannot be changed at run time:
          UMFPACK_COMPILED_WITH_BLAS=8       ! uses the BLAS
  end enum
  ! 9,12: unused

  ! Control(UMFPACK_STRATEGY) is one of the following:
  enum, bind(c) ; enumerator :: & 
       UMFPACK_STRATEGY_AUTO=0, &      ! use sym. or unsym. strategy
       UMFPACK_STRATEGY_UNSYMMETRIC, & ! COLAMD(A), coletree postorder, not prefer diag
       UMFPACK_STRATEGY_OBSOLETE, &    ! 2-by-2 is no longer available
       UMFPACK_STRATEGY_SYMMETRIC      ! AMD(A+A'), no coletree postorder, prefer diagonal
  end enum

  ! Control(UMFPACK_SCALE) is one of the following:
  enum, bind(c) ; enumerator :: & 
       UMFPACK_SCALE_NONE=0, &   ! no scaling
       UMFPACK_SCALE_SUM, &      ! default: divide each row by sum (abs (row))
       UMFPACK_SCALE_MAX         ! divide each row by max (abs (row))
  end enum

  ! Control(UMFPACK_ORDERING) and Info(UMFPACK_ORDERING_USED) are one of:
  enum, bind(c) ; enumerator :: & 
       UMFPACK_ORDERING_CHOLMOD=0, & ! use CHOLMOD (AMD/COLAMD then METIS)
       UMFPACK_ORDERING_AMD, &       ! use AMD/COLAMD
       UMFPACK_ORDERING_GIVEN, &     ! user-provided Qinit
       UMFPACK_ORDERING_METIS, &     ! use METIS
       UMFPACK_ORDERING_BEST, &      ! try many orderings, pick best
       UMFPACK_ORDERING_NONE, &      ! natural ordering
       UMFPACK_ORDERING_USER         ! user-provided function
     ! AMD/COLAMD means: use AMD for symmetric strategy, COLAMD for unsymmetric
  end enum

  ! status codes
  ! ------------
  enum, bind(c) ; enumerator :: & 
       UMFPACK_OK=0, &
       ! status > 0 means a warning, but the method was successful anyway.
       ! A Symbolic or Numeric object was still created.
       UMFPACK_WARNING_singular_matrix=1, &
       ! The following warnings were added in umfpack_*_get_determinant
       UMFPACK_WARNING_determinant_underflow=2, &
       UMFPACK_WARNING_determinant_overflow=3, &
       ! status < 0 means an error, and the method was not successful.
       ! No Symbolic of Numeric object was created.
       UMFPACK_ERROR_out_of_memory=-1, &
       UMFPACK_ERROR_invalid_Numeric_object=-3, &
       UMFPACK_ERROR_invalid_Symbolic_object=-4, &
       UMFPACK_ERROR_argument_missing=-5, &
       UMFPACK_ERROR_n_nonpositive=-6, &
       UMFPACK_ERROR_invalid_matrix=-8, &
       UMFPACK_ERROR_different_pattern=-11, &
       UMFPACK_ERROR_invalid_system=-13, &
       UMFPACK_ERROR_invalid_permutation=-15, &
       UMFPACK_ERROR_internal_error=-911, &       ! yes, call me if you get this      !
       UMFPACK_ERROR_file_IO=-17, &
       UMFPACK_ERROR_ordering_failed=-18
  end enum
  
  ! solve codes
  ! Solve the system ( )x=b, where ( ) is defined below.  "t" refers to the
  ! linear algebraic transpose (complex conjugate if A is complex), or the (')
  ! operator in MATLAB.  "at" refers to the array transpose, or the (.')
  ! operator in MATLAB.
  enum, bind(c) ; enumerator :: & 
       UMFPACK_A=0, &            ! Ax=b   
       UMFPACK_At=1, &           ! A'x=b  
       UMFPACK_Aat=2, &          ! A.'x=b 
       UMFPACK_Pt_L=3, &         ! P'Lx=b 
       UMFPACK_L=4, &            ! Lx=b   
       UMFPACK_Lt_P=5, &         ! L'Px=b 
       UMFPACK_Lat_P=6, &        ! L.'Px=b
       UMFPACK_Lt=7, &           ! L'x=b  
       UMFPACK_Lat=8, &          ! L.'x=b 
       UMFPACK_U_Qt=9, &         ! UQ'x=b 
       UMFPACK_U=10, &           ! Ux=b   
       UMFPACK_Q_Ut=11, &        ! QU'x=b 
       UMFPACK_Q_Uat=12, &       ! QU.'x=b
       UMFPACK_Ut=13, &          ! U'x=b  
       UMFPACK_Uat=14            ! U.'x=b 
  end enum

  ! --------------------------------------------------------------------------------
  ! End of Data from /usr/include/umfpack.h
  ! --------------------------------------------------------------------------------

  interface

     ! int umfpack_di_symbolic(int n_row,int n_col,const int Ap [ ],const int Ai
     ! [ ],const double Ax [ ],
     !
     ! void **Symbolic,const double Control [UMFPACK_CONTROL],double Info
     ! [UMFPACK_INFO]) ;
     integer(c_int) function umfpack_di_symbolic(&
          & n_row,n_col,Ap,Ai,Ax,Symbolic,Control,Info) bind(c)
       import c_int, c_ptr
       integer(c_int), value             :: n_row,n_col
       type(c_ptr),    value, intent(in) :: Ap,Ai
       type(c_ptr),    value, intent(in) :: Ax
       type(c_ptr)                       :: Symbolic
       type(c_ptr),    value             :: Control,Info
     end function umfpack_di_symbolic

     ! int umfpack_di_numeric(const int Ap [ ],const int Ai [ ],const double Ax
     ! [ ],
     !
     ! void *Symbolic,void **Numeric,const double Control
     ! [UMFPACK_CONTROL],double Info [UMFPACK_INFO]) ;
     integer(c_int) function umfpack_di_numeric(&
          & Ap,Ai,Ax,Symbolic,Numeric,Control,Info) bind(c) 
       import c_int, c_ptr
       type(c_ptr), value, intent(in) :: Ap,Ai
       type(c_ptr), value, intent(in) :: Ax
       type(c_ptr), value, intent(in) :: Symbolic
       type(c_ptr)                    :: Numeric
       type(c_ptr), value             :: Control,Info
     end function umfpack_di_numeric
     
     ! int umfpack_di_solve(int sys,const int Ap [ ],const int Ai [ ],const
     ! double Ax [ ],double X [ ],const double B [ ],
     !
     ! void *Numeric,const double Control [UMFPACK_CONTROL],double Info
     ! [UMFPACK_INFO]) ;
     integer(c_int) function umfpack_di_solve(&
          & sys,Ap,Ai,Ax,X,B,Numeric,Control,Info) bind(c)
       import c_int, c_ptr
       integer(c_int), value :: sys
       type(c_ptr), value, intent(in) :: Ap,Ai
       type(c_ptr), value, intent(in) :: Ax
       type(c_ptr), value             :: X
       type(c_ptr), value, intent(in) :: B
       type(c_ptr), value, intent(in) :: Numeric
       type(c_ptr), value             :: Control,Info
     end function umfpack_di_solve
     
     ! void umfpack_di_free_symbolic(void **Symbolic) ;
     subroutine umfpack_di_free_symbolic(Symbolic) bind(c)
       import c_ptr
       type(c_ptr) :: Symbolic
     end subroutine umfpack_di_free_symbolic
     
     ! void umfpack_di_free_numeric(void **Numeric) ;
     subroutine umfpack_di_free_numeric(Numeric) bind(c)
       import c_ptr
       type(c_ptr) :: Numeric
     end subroutine umfpack_di_free_numeric
     
     ! void umfpack_di_defaults(double Control [UMFPACK_CONTROL]) ;
     subroutine umfpack_di_defaults(Control) bind(c)
       import c_ptr
       type(c_ptr),value :: Control
     end subroutine umfpack_di_defaults
     
  end interface
  
contains
  
  integer function f_umfpack_di_symbolic(n_row,n_col,Ap,Ai,Ax,Symbolic,Control,Info)
    integer, intent(in) :: n_row,n_col
    integer, target, intent(in) :: Ap(*), Ai(*)
    real(pr), target, intent(in) :: Ax(*)
    type(c_ptr) :: Symbolic
    real(pr), target, optional :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
    
    integer :: c_n_row,c_n_col
    type(c_ptr) :: c_Ap,c_Ai,c_Ax,c_Symbolic,c_Control,c_Info
    
    c_n_row=n_row
    c_n_col=n_col
    
    c_Ap=c_loc(Ap)
    c_Ai=c_loc(Ai)
    c_Ax=c_loc(Ax)
    c_Symbolic=Symbolic
    c_Control=c_null_ptr
    c_Info=c_null_ptr
    if (present(Control)) c_Control=c_loc(Control)
    if (present(Info)) c_Info=c_loc(Info)

    f_umfpack_di_symbolic = umfpack_di_symbolic(&
         c_n_row,c_n_col,c_Ap,c_Ai,c_Ax,c_Symbolic,c_Control,c_Info)

  end function f_umfpack_di_symbolic

  integer function f_umfpack_di_numeric(Ap,Ai,Ax,Symbolic,Numeric,Control,Info)
    integer, target, intent(in) :: Ap(*),Ai(*)
    real(pr), target, intent(in) :: Ax(*)
    type(c_ptr), intent(in) :: Symbolic
    type(c_ptr) :: Numeric
    real(pr), target, optional :: Control(0:UMFPACK_CONTROL-1),Info(0:UMFPACK_INFO-1)
    
    type(c_ptr) :: c_Ap,c_Ai,c_Ax,c_Symbolic,c_Numeric,c_Control,c_Info
    
    c_Ap=c_loc(Ap)
    c_Ai=c_loc(Ai)
    c_Ax=c_loc(Ax)
    
    c_Symbolic=Symbolic
    c_Numeric=Numeric
    c_Control=c_null_ptr
    c_Info=c_null_ptr
    if (present(Control)) c_Control=c_loc(Control)
    if (present(Info)) c_Info=c_loc(Info)

    f_umfpack_di_numeric = umfpack_di_numeric(&
         c_Ap,c_Ai,c_Ax,c_Symbolic,c_Numeric,c_Control,c_Info)

    Numeric=c_Numeric

  end function f_umfpack_di_numeric

end module umfpack_calls
