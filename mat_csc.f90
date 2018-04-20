!> @brief Matrices creuses au format csc
!> @date Time-stamp: <2017-10-26 03:11:58 yves>
!
!> @details
!> Les indices de colonne et les valeurs des coefficients de la ligne
!> i sont dans col(row(i):row(i+1)-1) et a(row(i):row(i+1)-1).
!>
module mat_csc

  use prec
  use matrices  ! paramètres communs
  use iso_c_binding
  implicit none

  private

  ! Pour la matrice
  !v     [ 3 0 1 0 0 ]
  !v     [ 0 4 0 0 0 ]
  !v     [ 0 7 5 9 0 ]
  !v     [ 0 0 0 0 2 ]
  !v     [ 0 0 0 6 5 ]
  !v col :: 1 2   4   6   8   10
  !v row :: 1 2 3 1 3 3 5 4 5
  !v  a  :: 3 4 7 1 5 9 6 2 5
  !v        -----------------
  !v        1 2 3 4 5 6 7 8 9
  !> @brief la structure de donnée des matrices csc
  type csc
     private
     logical :: c_style = .false. ! True to have array indices starting at 0,
                                  ! false at 1
     integer                             :: nl=0, nc=0, nnz=0
     integer,  dimension(:), pointer :: col => NULL()
     integer,  dimension(:), pointer :: row => NULL()
     real(pr), dimension(:), pointer :: a => NULL()
     ! .true. si les indices de colonne sont par ordre croissant
     logical                             :: sorted = .false.
     logical :: is_factorized = .false. ! whether the LU or Choleski
     type(c_ptr) :: ptr_symbolic = c_null_ptr ! pointers used by the UMFPACK library
     type(c_ptr) :: ptr_numeric = c_null_ptr ! pointers used by the UMFPACK library
  end type csc

  !> @brief libère la mémoire de la matrice
  interface deallocate
     module procedure deallocate_csc
  end interface deallocate

  !> @brief allocation mémoire pour une matrice nulle.
  interface zeros
     module procedure zeros_csc, csc_create_from_cols
     module procedure csc_create_from_nnz_per_col, csc_create_from_cols_list
  end interface zeros

  !> @brief allocation mémoire pour une matrice nulle.
  interface allocate
     module procedure zeros_csc, csc_create_from_cols
     module procedure csc_create_from_nnz_per_col, csc_create_from_cols_list
  end interface allocate

  !> @brief construit une matrice csc à partir d'une autre matrice, s
  !> := a. La matrice s est csc mais la matrice a peut être csc ou
  !> m_list.
  interface assignment(=)
     module procedure csc_assign_from_m_list
  end interface assignment(=)

  !> @brief Sort columns in ascending row index order
  interface sort_cols
     module procedure sort_columns_csc
  end interface sort_cols

  !> @brief construction de la matrice par ajout ou assignation de
  !>        coefficients
!!$  interface mat_add_value
!!$     module procedure construct_entry_csc, add_mat_csc
!!$  end interface mat_add_value

  !> @brief construction de la matrice par ajout ou assignation de
  !>        coefficients diagonaux
!!$  interface mat_add_diag
!!$     module procedure add_diag_csc
!!$  end interface mat_add_diag

  !> @brief transpose une matrice
!!$  interface transpose
!!$     module procedure transpose_csc
!!$  end interface transpose

  !> @brief Matrix-vector product -- overloads the standard *
  !> operator.
  interface operator(*)
     module procedure csc_times_array
  end interface operator(*)

  !> @brief Matrix-vector product
  interface mat_times_vec
     module procedure csc_times_array
  end interface mat_times_vec

  !> @brief Factorize the matrix (UMFPACK)
  interface factorize
     module procedure factorize_csc
  end interface factorize

  !> @brief Solve the factorized system (UMFPACK)
  interface solve
     module procedure solve_csc
  end interface solve

  !> @brief affiche la matrice
  interface print
     module procedure print_csc
  end interface print

  !> @brief Renvoie le nombre de coefficients non-nuls
  interface get_nnz
     module procedure get_nnz_csc
  end interface get_nnz

  public :: csc                          ! structure de données publique
  public :: deallocate                   ! destructeur
  public :: zeros                        ! constructeurs
  public :: allocate, assignment(=)      ! constructeurs
!!$  public :: matrice_poisson2D            ! constructeurs
!!$  public :: transpose                    ! opérateurs unitaires
  public :: get_nl, get_nc, get_nnz      ! interrogations
  public :: print, factorize, solve, sort_cols
  public :: operator(*), mat_times_vec

contains

  ! ************************************************************************
  ! *                                                                      *
  ! *  Destructeur(s)                                                      *
  ! *                                                                      *
  ! ************************************************************************

  subroutine deallocate_csc(s)
    use umfpack_calls

    type(csc), intent(inout) :: s

    s%c_style = .false.
    s%nl = 0
    s%nc = 0
    s%nnz = 0
    if (associated(s%col))  deallocate(s%col)
    if (associated(s%row))  deallocate(s%row)
    if (associated(s%a))    deallocate(s%a)
    nullify(s%col,s%row,s%a)
    s%sorted = .false.
    s%is_factorized = .false.

    ! umfpack_di_free_symbolic : deallocate the Symbolic object and sets the
    ! Symbolic handle to c_null_ptr.
    if (c_associated(s%ptr_symbolic)) call umfpack_di_free_symbolic(s%ptr_symbolic)

    ! umfpack_di_free_numeric : deallocate the Numeric object and sets the
    ! Symbolic handle to c_null_ptr.
    if (c_associated(s%ptr_numeric)) call umfpack_di_free_numeric(s%ptr_numeric)

  end subroutine deallocate_csc

  ! ************************************************************************
  ! *                                                                      *
  ! *  Constructeurs                                                       *
  ! *                                                                      *
  ! ************************************************************************

  !> @brief Allocation initiale pour une matrice de taille n*p nulle
  !>        (sans coefficients) mais avec potentiellement nnz
  !>        coefficients non nuls
  !
  !> @param[in,out] M matrice
  !> @param[in] n,p taille de la matrice
  !> @param[in] nnz taille réservée pour les coefficients
  !
  !> @details alloue de la mémoire une matrice n * p avec nnz
  !>          coefficients non nuls
  subroutine zeros_csc(s,n,p, nnz)

    type(csc), intent(inout) :: s
    integer,   intent(in)    :: n,p, nnz

    integer :: i

    s%c_style = .false.
    s%nl = n
    s%nc = p
    s%nnz = nnz
    allocate(s%col(p+1), s%row(nnz), s%a(nnz))

    do i = 1,nnz
       s%row(i) = 0
       s%a(i) = 0._pr
    end do

    do i = 1, p+1
       s%col(i) = 1
    end do

    s%sorted = .false.
    s%is_factorized = .false.
    s%ptr_symbolic = c_null_ptr
    s%ptr_numeric = c_null_ptr

  end subroutine zeros_csc

  !> @brief Allocation initiale pour une matrice de taille n*p à partir du
  !> tableau nnz(:) du nombre de coefficients non nuls par colonne
  !
  !> @param[in,out] M matrice
  !> @param[in] nnz(:) taille réservée pour les coefficients. La taille de nnz
  !>                   est le nombre de colonne p.
  !> @details On ne connait pas le nombre de lignes et on suppose donc que la
  !>          matrice est carrée
  !> @details alloue de la mémoire une matrice n * p avec nnz coefficients non
  !>          nuls par colonne et prépare la structure de donnée connaissant
  !>          exactement le profil de la matrice
  subroutine csc_create_from_nnz_per_col(s,nnz)

    type(csc)            , intent(inout) :: s
    integer, dimension(:), intent(in)    :: nnz

    integer :: nnz_total, j, k

    s%c_style = .false.
    s%nc = size(nnz)
    s%nl = s%nc
    nnz_total  = sum(nnz(:))
    s%nnz = nnz_total
    allocate(s%col(s%nc+1), s%row(nnz_total), s%a(nnz_total))

    do j = 1, nnz_total
       s%row(j) = 0
       s%a(j)   = 0._pr
    end do

    ! Construction de s%col pour tenir compte du profil donné en entrée
    k = 1
    do j = 1, s%nc
       s%col(j) = k
       k = k + nnz(j)
    end do
    s%col(s%nc+1) = k

    s%sorted = .false.
    s%is_factorized = .false.
    s%ptr_symbolic = c_null_ptr
    s%ptr_numeric = c_null_ptr

  end subroutine csc_create_from_nnz_per_col

  !> @brief Allocation initiale pour une matrice de taille n*p à partir du
  !> tableau cols(:,:) des coefficients non nuls par colonne, mais sans la
  !> valeur des coefficients
  !
  !> @param[in,out] M matrice

  !> @param[in] cols(:,:), où cols(:,j) est la liste des indices des
  !>            coefficients non nuls de la colonne j.
  !
  !> @details On ne connait pas le nombre de lignes et on suppose donc que la
  !>          matrice est carrée
  subroutine csc_create_from_cols(s,cols)

    type(csc),                 intent(inout) :: s
    integer,   dimension(:,:), intent(in)    :: cols

    integer :: nnz, k, i, j, row, p

#ifdef DEBUG
    if (size(shape(cols))/=2) then
       print "('csc_create_from_cols: the shape of cols is not (:,:)')"
       return
    end if
#endif

    s%c_style = .false.
    ! Taille de la matrice
    s%nc = size(cols,2)   ! nombre de colonnes
    s%nl = s%nc           ! matrice supposée carrée
    p = size(cols,1)      ! nombre max d'elt non nul sur une colonne
    ! Decompte du nombre d'elements non nuls au total
    nnz = 0
    do j = 1,s%nc
       nnz = nnz + count(cols(:,j)>0) ! somme des indices positifs de la col j
    end do
    s%nnz = nnz

#ifdef DEBUG
    print "('csc_create_from_cols: nombre d''éléments non nuls : ',I0)", nnz
#endif

    ! Allocation memoire
    allocate(s%col(s%nc+1), s%row(nnz), s%a(nnz))

    do i = 1, nnz
       s%a(i)   = 0._pr ! Les coefficients sont encore inconnus
    end do

    ! Construction de la matrice
    k = 1
    do j = 1,s%nc
       s%col(j) = k
       do i = 1,p
          row = cols(i,j)
          if (row>0) then ! On a un coeff non nul, en ligne row et col j
             s%row(k) = row
             k = k+1
          end if
       end do
    end do
    s%col(s%nc+1) = k

    s%sorted = .false.
    s%is_factorized = .false.
    s%ptr_symbolic = c_null_ptr
    s%ptr_numeric = c_null_ptr

#ifdef DEBUG
    print "('csc_create_from_cols:')"
    print "('nnz : ',I0,' -- nl : ',I0,' -- nc : ',I0)", nnz, s%nl, s%nc
    print "('row : ',10I7)", s%row(:)
    print "('col : ',10I7)", s%col(:)
    print "(' a  : ',10E7.1)", s%a(:)
#endif

  end subroutine csc_create_from_cols


  !> @brief Allocation initiale pour une matrice de taille n*p à partir du
  !> tableau voisins(:) des coefficients non nuls par colonne
  !
  !> @param[in,out] M matrice

  !> @param[in] voisins(:), où voisins(j) est la liste (type liste
  !>            chainée du module list) des coefficients non nuls de
  !>            la colonne j
  !
  !> @details On ne connait pas le nombre de lignes et on suppose donc que la
  !>          matrice est carrée
  subroutine csc_create_from_cols_list(s,voisins)

    use list
    type(csc),                  intent(inout) :: s
    type(i_list), dimension(:), intent(in)    :: voisins

    integer :: nnz, i, k, j

    s%c_style = .false.
    ! Taille de la matrice
    s%nc = size(voisins) ! nombre de colonnes
    s%nl = s%nc          ! matrice carrée
    ! Nombre d'elements non nuls
    nnz = 0
    do i = 1,s%nc
       nnz = nnz + get_length(voisins(i))
    end do
    s%nnz = nnz

#ifdef DEBUG
    print "('csc_create_from_cols_list: nombre d''éléments non nuls : ',I0)", nnz
#endif

    ! Allocation memoire
    allocate(s%col(s%nc+1), s%row(nnz), s%a(nnz))

    do i = 1,nnz
       s%a(i)   = 0._pr
    end do

    ! Construction de la matrice
    k = 1
    do j = 1,s%nc
       s%col(j) = k
       do i = 1,get_length(voisins(i)) ! Nombre de voisins de i
          s%row(k) = get_entry(voisins(j),i)
          k = k+1
       end do
    end do
    s%col(s%nc+1) = k

    s%sorted = .false.
    s%is_factorized = .false.
    s%ptr_symbolic = c_null_ptr
    s%ptr_numeric = c_null_ptr

#ifdef DEBUG
    print "('csc_create_from_cols_list:')"
    print "('nnz : ',I0,' -- nl : ',I0,' -- nc : ',I0)", nnz, s%nl, s%nc
    print "('row : ',10I7)", s%row(:)
    print "('col : ',10I7)", s%col(:)
    print "(' a  : ',10E7.1)", s%a(:)
#endif

  end subroutine csc_create_from_cols_list


  !> @brief Allocation et remplissage pour une matrice csc à partir
  !>        d'une matrice list
  !
  !> @param[out] s matrice csc
  !> @param[in] s_list matrice de type m_list
  subroutine csc_assign_from_m_list(s,s_list)

    use mat_list
    type(csc),    intent(out) :: s
    type(m_list), intent(in)  :: s_list

    integer :: k,i,j
    integer, dimension(:), allocatable :: row
    real(pr), dimension(:), allocatable :: a

    s%c_style = .false.
    ! Dimensions
    s%nl = size(s_list,1)
    s%nc = size(s_list,2)
    s%nnz = get_nnz(s_list)

    ! Allocation memoire
    allocate(s%col(s%nc+1), s%row(s%nnz), s%a(s%nnz))

#ifdef DEBUG
    print "('csc_assign_from_m_list: nombre d''éléments non nuls : ',I0)", nnz
#endif

    ! On construit colone par colonne
    k = 1
    do j = 1,s%nc ! Parcours des colonnes
       s%col(j) = k ! Début des indices de la colonne j dans les tableaux %row et %a
       call get_col(s_list,j, row,a) ! Renvoie la colonnne comme deux tableaux
       do i = 1,size(row)
          s%row(k) = row(i)
          s%a(k)   = a(i)
          k = k+1
       end do
       deallocate(row,a)
    end do
    s%col(s%nc+1) = k

    s%sorted = .false.
    s%is_factorized = .false.
    s%ptr_symbolic = c_null_ptr
    s%ptr_numeric = c_null_ptr

#ifdef DEBUG
    print "('csc_assign_from_m_list:')"
    print "('nnz : ',I0,' -- nl : ',I0,' -- nc : ',I0)", s%nnz, s%nl, s%nc
    print "('row : ',10I7)", s%row(:)
    print "('col : ',10I7)", s%col(:)
    print "(' a  : ',10E7.1)", s%a(:)
#endif

  end subroutine csc_assign_from_m_list


!!$  !> @brief Allocation et remplissage pour une matrice csc à partir
!!$  !>        d'une matrice csc (copie de toute la structure de donnée)
!!$  !
!!$  !> @param[out] s matrice csc
!!$  !> @param[in] a matrice de type m_list
!!$  !
!!$  !> @details alloue de la mémoire une matrice n * p avec nnz
!!$  !>          coefficients non nuls par ligne et prépare la structure
!!$  !>          de donnée connaissant exactement le profil de la matrice
!!$  !>          et aussi les valeurs de coefficients non nuls
!!$  subroutine csc_assign_from_csc(s,a)
!!$
!!$    type(csc), intent(in)    :: a
!!$    type(csc), intent(out)   :: s
!!$
!!$    integer :: nnz,j,k
!!$
!!$    s%c_style = a%c_style
!!$    s%nl = a%nl
!!$    s%nc = a%nc
!!$    s%nnz = a%nnz
!!$    nnz = size(a%row)
!!$    allocate(s%col(s%nc+1), s%row(nnz), s%a(nnz))
!!$
!!$    do j = 1,s%nc+1
!!$       s%col(j) = a%col(j)
!!$    end do
!!$
!!$    do k = 1,nnz
!!$       s%row(k) = a%row(k)
!!$       s%a(k)   = a%a(k)
!!$    end do
!!$
!!$    s%sorted = a%sorted
!!$    s%is_factorized = a%is_factorized
!!$    s%ptr_symbolic = ! We have to make a copy of the structure !!
!!$    s%ptr_numeric = ! We have to make a copy of the structure !!
!!$
!!$  end subroutine csc_assign_from_csc


  ! ************************************************************************
  ! *                                                                      *
  ! *  Méthodes "opérateurs unitaires"                                     *
  ! *                                                                      *
  ! ************************************************************************

    ! @brief Carry out the LU (general band) factorization of the matrix (see
  ! UMFPACK user guide)
  subroutine factorize_csc(A)
    use umfpack_calls

    type(csc), intent(inout) :: A

    integer :: status
    real(pr), dimension(0:UMFPACK_CONTROL-1), target :: control = 0._pr
    real(pr), dimension(0:UMFPACK_INFO-1), target :: info = 0._pr

    ! We need to have indices starting from 0: c_style
    if (.not.A%c_style) then
       A%col(:) = A%col(:) - 1
       A%row(:) = A%row(:) - 1
       A%c_style = .true.
    end if
    ! We have to use *_di_* (double and int) umfpack subroutines

    ! umfpack_di_symbolic : pre-orders the columns of A, returns an opaque
    ! symbolic object
    status = umfpack_di_symbolic(A%nl,A%nc, c_loc(A%col), c_loc(A%row), c_loc(A%a), &
         A%ptr_symbolic, c_loc(control),c_loc(info))
    print "('CSC symbolic factorization (UMFPACK) P*A=L*U done with status ',I0)", status
    ! umfpack_di_numeric : scales and factorizes A into PAQ, PRAQn PR^{-1}AQ,
    ! into LU, depending on the strategy choosen by the method. Returns an
    ! opaque numeric object. This object is used to solve the subsequent linear
    ! systems.
    status = umfpack_di_numeric(c_loc(A%col), c_loc(A%row), c_loc(A%a), &
         A%ptr_symbolic, A%ptr_numeric, c_loc(control),c_loc(info))
    !c_null_ptr, c_null_ptr)
    print "('CSC numeric factorization (UMFPACK) P*A=L*U done with status ',I0)", status

    ! umfpack_di_free_symbolic : deallocate the Symbolic object and sets the
    ! Symbolic handle to c_null_ptr.
    call umfpack_di_free_symbolic(A%ptr_symbolic)

    A%is_factorized = .true.

  end subroutine factorize_csc


  ! @brief Resolution of one linear system Ax=b, after the matrix A has been
  ! factorized.
  !
  ! @details
  subroutine solve_csc(A,x,b)
    use umfpack_calls

    type(csc), intent(in) :: A
    real(pr), dimension(:), target, intent(out)  :: x
    real(pr), dimension(:), target, intent(in)   :: b

    integer :: status
    real(pr), dimension(0:UMFPACK_CONTROL-1), target :: control
    real(pr), dimension(0:UMFPACK_INFO-1), target :: info

    ! sys=UMFPACK_A means that we solve A*x=b
    status = umfpack_di_solve(UMFPACK_A, c_loc(A%col), c_loc(A%row), c_loc(A%a), &
         & c_loc(x), c_loc(b), A%ptr_numeric, c_loc(control), c_loc(info))
    print "('CSC resolution (UMFPACK) done with status ',I0)", status

  end subroutine solve_csc

  !>@brief Sort columns in ascending order of row index
  subroutine sort_columns_csc(s)
    use sort

    type(csc), intent(inout) :: s

    integer :: j,i_deb,i_fin,i, n
    integer, dimension(:), allocatable :: is, new_i
    real(pr), dimension(:), allocatable :: as

    do j = 1, s%nc
       i_deb = s%col(j)
       i_fin = s%col(j+1)-1
       allocate(is(i_fin-i_deb+1), as(i_fin-i_deb+1), new_i(i_fin-i_deb+1))
       is = s%row(i_deb:i_fin)
       as = s%a(i_deb:i_fin)
       call sort_insertion(is,new_i) ! Descending order, and we want ascending
       n = size(is)
       do i = 1,n
          s%row(i_deb-1+new_i(i)) = is(n+1-i)
          s%a(i_deb-1+new_i(i)) = as(n+1-i)
       end do
       deallocate(is,as,new_i)
    end do

    print "('Row index in each column are sorted')"
    s%sorted = .true.

  end subroutine sort_columns_csc

  ! ************************************************************************
  ! *                                                                      *
  ! *  Méthodes "opérateurs binaires"                                     *
  ! *                                                                      *
  ! ************************************************************************

  ! csc_add_lambda_id
  ! csc_times_scalar


  !> @brief Computes the matrix-vector product v=s*u. This is an
  !> _orphaned_ OMP subroutine.

  !> @param[in] s matrix
  !> @param[in] u vector that multiply the matrix
  !> @param[out] v output vector
  function csc_times_array(s,u) result(v)

    type(csc),                 intent(in)  :: s
    real(pr), dimension(s%nc), intent(in)  :: u
    real(pr), dimension(s%nl)              :: v ! := s*u

    ! Les variables locales sont _private_ par defaut.
    integer :: j, kDeb, kFin, i, k
    real(pr) :: y

#ifdef DEBUG
    !$OMP SINGLE
    if (size(shape(u))/=1.or.size(shape(v))/=1.or.&
         size(u,1)/=s%nc.or.size(v,1)/=s%nl) then
       print *,'In axCSC: dimension mismatch'
    end if
    !$OMP END SINGLE
    ! Il y a une barriere implicite en fin de clause SINGLE
#endif

    ! Multiplication
    !$OMP DO SCHEDULE(RUNTIME)
    do i = 1,s%nl
       v(i) = 0._pr
    end do
    !$OMP END DO
    !$OMP DO SCHEDULE(RUNTIME)
    do j = 1,s%nc
       kDeb = s%col(j)
       kFin = s%col(j+1)-1
       y = u(j)
       do k = kDeb,kFin
          i = s%row(k)
          v(i) = v(i) + s%a(k)*y
       end do
    end do
    !$OMP END DO

  end function csc_times_array


  ! ************************************************************************
  ! *                                                                      *
  ! *  Méthodes interrogation de la matrice                                *
  ! *                                                                      *
  ! ************************************************************************

  !> @brief Get the number of rows
  integer function get_nl(S)
    type(csc), intent(in) :: S
    get_nl = S%nl
  end function get_nl

  !> @brief Get the number of columns
  integer function get_nc(S)
    type(csc), intent(in) :: S
    get_nc = S%nc
  end function get_nc

  !> @brief Get the number nonzeros elements
  integer function get_nnz_csc(S)
    type(csc), intent(in) :: S
    get_nnz_csc = S%nnz
  end function get_nnz_csc

  ! ************************************************************************
  ! *                                                                      *
  ! *  Méthodes Calcul de matrices types                                   *
  ! *                                                                      *
  ! ************************************************************************


!!$  !> @brief Calcul de la matrice DF du pb de poisson, u -->
!!$  !> -Laplacien(u) sur le carré (0,1)x(0,1) avec n divisions dans
!!$  !> chaque direction et conditions de Dirichlet.
!!$
!!$  !> @details La matrice est de taille (n-1)^2 et ses coefficients sont
!!$  !>
!!$  !> A = 1/h^2* *(a_{ij}) avec i,j = 1 .. n-1
!!$  !>
!!$  !> et a_{ii} = 4
!!$  !>    a_{ii+1}=a_{ii-1}=a_{ii-(n-1)}=a_{ii+(n-1)}=-1
!!$  !>
!!$  !> avec h = 1/n.
!!$  function matrice_poisson2D(n) result(A)
!!$    type(csc) :: A
!!$    integer, intent(in) :: n
!!$
!!$    real(pr) :: h2
!!$    integer :: i,j,ij,k, nnz
!!$
!!$    ! h^2
!!$    h2 = 1._pr/real(n,pr)**2
!!$
!!$    ! Allocation memoire
!!$    A%nl = (n-1)**2
!!$    A%nc = (n-1)**2
!!$    nnz = 5*(n-1)**2 - 4*(n-1)
!!$    allocate(A%row(A%nl+1), A%col(nnz), A%a(nnz))
!!$
!!$    ! Initialisation
!!$    A%row = 0
!!$    A%col = 0
!!$    A%a = 0._pr
!!$
!!$    ! Construction du 'profil' ou 'graphe' de la matrice
!!$    ! et remplissage de la matrice
!!$    k = 1
!!$    do j = 1,n-1
!!$       do i = 1,n-1
!!$          ij = i+(j-1)*(n-1)
!!$          A%row(ij) = k ! Indice du debut de stockage de la ligne 'ij'
!!$          if (j>1) then
!!$             A%col(k) = ij-(n-1)
!!$             A%a  (k) = -1._pr / h2
!!$             k = k+1
!!$          end if
!!$          if (i>1) then
!!$             A%col(k) = ij-1
!!$             A%a  (k) = -1._pr / h2
!!$             k = k+1
!!$          end if
!!$          A%col(k) = ij
!!$          A%a  (k) = 4._pr / h2
!!$          k = k+1
!!$          if (i<n-1) then
!!$             A%col(k) = ij+1
!!$             A%a  (k) = -1._pr / h2
!!$             k = k+1
!!$          end if
!!$          if (j<n-1) then
!!$             A%col(k) = ij+(n-1)
!!$             A%a  (k) = -1._pr / h2
!!$             k = k+1
!!$          end if
!!$       end do
!!$    end do
!!$    A%row(A%nl+1) = k
!!$
!!$    A%sorted = .true.
!!$
!!$#ifdef DEBUG
!!$    print *,'k :',k
!!$    print *,'nnz :',nnz,' nc :', A%nc,' nl :', A%nl
!!$    print *,'ROW :',A%row
!!$    print *,'COL :',A%col
!!$    print *,' A  :',A%a
!!$    print *,'DIAG:',A%diag
!!$#endif
!!$
!!$  end function matrice_poisson2D

  ! ************************************************************************
  ! *                                                                      *
  ! *  Méthodes "I/O"                                             *
  ! *                                                                      *
  ! ************************************************************************

  ! csc_save
  ! csc_load


  !> @brief affiche la matrice
  !
  !> @param[in] S matrice
  subroutine print_csc(S)
    type(csc), intent(in) :: s

    integer :: j,k,kDeb,kFin,offset

    print "('Matrice <csc> de taille ',2I9)", s%nl,s%nc
    ! Affichage par colonne
    offset = 0
    if (s%c_style) offset = 1
    do j = 1-offset, s%nc-offset
       kDeb = s%col(offset+j)     ! Debut du stockage de la colonne j
       kFin = s%col(offset+j+1)-1 ! Fin du stockage de la colonne j
       do k=offset+kDeb,offset+kFin
          print "('(',I7,',',I7,') --> ',E13.5)", s%row(k),j,s%a(k)
       end do
    end do

  end subroutine print_csc

end module mat_csc
