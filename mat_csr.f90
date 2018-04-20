!> @brief Matrices creuses au format csr
!> @date Time-stamp: <2017-10-26 02:49:47 yves>
!
!> @details C'est le format le plus naturel (efficace ?) pour faire
!>          des produts matrice-vecteur.  Dans ce format, une matrice
!>          est stockée a l'aide de 3 tableaux:
!>   - row : contient l'indice dans les deux autres tableaux des coefficients
!>           de ligne. Il est de taille n+1 pour une matrice n*p
!>   - col : contient les indices de colonne, ligne apres ligne
!>   - a   : contient les coeffcients

!> Les indices de colonne et les valeurs des coefficients de la ligne
!> i sont dans col(row(i):row(i+1)-1) et a(row(i):row(i+1)-1).
!>
!> @todo mettre des macro openmp partout où c'est possible
!> @author Yves Coudière
!> - Code démarré vers 2006
!> @author Charles Pierre
!> - Contributions, debuggages
!> @author Rodolphe Turpault
!> - Contributions diverses
!> @author Adnane Azzouzi
!> - Contribution diverses (pendant son post-doc seulement)
!> @date mai 2013
!> - récupération depuis le svn et mise à plat pour une utilisation
!>   personnelle
module mat_csr

  use prec
  use matrices  ! paramètres communs
  implicit none

  private

  ! Pour la matrcice
  !v     [ 3 0 1 0 0 ]
  !v     [ 0 4 0 0 0 ]
  !v     [ 0 7 5 9 0 ]
  !v     [ 0 0 0 0 2 ]
  !v     [ 0 0 0 6 5 ]
  !v row :: 1 3 4 7 8 10
  !v col :: 1 3 2 2 3 4 5 4 5
  !v  a  :: 3 1 4 7 5 9 2 6 5
  !v        -----------------
  !v        1 2 3 4 5 6 7 8 9
  !> @brief la structure de donnée des matrices csr
  type csr
     private
     integer                             :: nl=0, nc=0
     integer,  dimension(:), allocatable :: row
     integer,  dimension(:), allocatable :: col
     real(pr), dimension(:), allocatable :: a
     ! .true. si les indices de colonne sont par ordre croissant
     logical                             :: sorted=.false.
  end type csr

  !> @brief libère la mémoire de la matrice
  interface deallocate
     module procedure deallocate_csr
  end interface deallocate

  !> @brief allocation mémoire pour une matrice nulle.
  interface zeros
     module procedure zeros_csr, csr_create_from_rows
     module procedure csr_create_from_nnz_per_row, csr_create_from_rows_list
  end interface zeros

  !> @brief allocation mémoire pour une matrice nulle.
  interface allocate
     module procedure zeros_csr, csr_create_from_rows
     module procedure csr_create_from_nnz_per_row, csr_create_from_rows_list
  end interface allocate

  !> @brief construit une matrice csr à partir d'une autre matrice, s
  !> := a. La matrice s est csr mais la matrice a peut être csr ou
  !> m_list.
  interface assignment(=)
     module procedure csr_assign_from_m_list, csr_assign_from_csr
  end interface assignment(=)

  !> @brief Matrix-vector product -- overloads the standard *
  !> operator.
  interface operator(*)
     module procedure csr_times_array
  end interface operator(*)

  !> @brief Matrix-vector product
  interface mat_times_vec
     module procedure csr_times_array
  end interface mat_times_vec

  !> @brief construction de la matrice par ajout ou assignation de
  !>        coefficients
  interface mat_add_value
     module procedure construct_entry_csr, add_mat_csr
  end interface mat_add_value

  !> @brief construction de la matrice par ajout ou assignation de
  !>        coefficients diagonaux
  interface mat_add_diag
     module procedure add_diag_csr
  end interface mat_add_diag

  !> @brief transpose une matrice
  interface transpose
     module procedure transpose_csr
  end interface transpose

  !> @brief affiche la matrice
  interface print
     module procedure print_csr
  end interface print

  !> @brief Renvoie le nombre de coefficients non-nuls
  interface get_nnz
     module procedure get_nnz_csr
  end interface get_nnz

  public :: csr                          ! structure de données publique
  public :: deallocate                   ! destructeur
  public :: zeros, mat_add_value         ! constructeurs
  public :: mat_add_diag                 ! constructeurs
  public :: allocate, assignment(=)      ! constructeurs
  public :: matrice_poisson2D            ! constructeurs
  public :: transpose                    ! opérateurs unitaires
  public :: operator(*)                  ! opérateurs binaires
  public :: get_nl, get_nc, get_nnz      ! interrogations
  public :: print

contains

  ! ************************************************************************
  ! *                                                                      *
  ! *  Destructeur(s)                                                      *
  ! *                                                                      *
  ! ************************************************************************

  subroutine deallocate_csr(s)

    type(csr), intent(inout) :: s

    if (allocated(s%row))  deallocate(s%row)
    if (allocated(s%col))  deallocate(s%col)
    if (allocated(s%a))    deallocate(s%a)
    s%nl = 0
    s%nc = 0
    s%sorted = .false.

  end subroutine deallocate_csr

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
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine zeros_csr(s,n,p, nnz)

    type(csr), intent(inout) :: s
    integer,   intent(in)    :: n,p, nnz

    integer :: i

    s%nl = n
    s%nc = p
    allocate(s%row(n+1), s%col(nnz), s%a(nnz))

    do i = 1,nnz
       s%col(i) = 0
       s%a(i) = 0._pr
    end do

    do i = 1, n+1
       s%row(i) = 1
    end do

  end subroutine zeros_csr

  !> @brief Allocation initiale pour une matrice de taille n*p à
  !> partir du tableau nnz(:) du nombre de coefficients non uls par
  !> ligne
  !
  !> @param[in,out] M matrice
  !> @param[in] nnz(:) taille réservée pour les coefficients. La
  !>                   taille de nnz est le nombre de ligne n.

  !> @details On ne connait pas le nombre de colonnes et on suppose
  !>          donc que la matrice est carrée
  !
  !> @details alloue de la mémoire une matrice n * p avec nnz
  !>          coefficients non nuls par ligne et prépare la structure
  !>          de donnée connaissant exactement le profil de la matrice
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine csr_create_from_nnz_per_row(s,nnz)

    type(csr)            , intent(inout) :: s
    integer, dimension(:), intent(in)    :: nnz

    integer :: nnz_total, i, k

    s%nl = size(nnz)
    s%nc = s%nl
    nnz_total  = sum(nnz(:))
    allocate(s%row(s%nl+1), s%col(nnz_total), s%a(nnz_total))

    do i = 1, nnz_total
       s%col(i) = 0
       s%a(i)   = 0._pr
    end do

    s%sorted = .false.

    ! Construction de s%row pour tenir compte du profil donné en
    ! entrée
    k = 1
    do i = 1, s%nl
       s%row(i) = k
       k = k + nnz(i)
    end do
    s%row(s%nl+1) = k

  end subroutine csr_create_from_nnz_per_row

  !> @brief Allocation initiale pour une matrice de taille n*p à
  !> partir du tableau lignes(:,:) des coefficients non nuls par ligne
  !
  !> @param[in,out] M matrice
  !> @param[in] lignes(:,:), où ligne(:,i) est la liste des
  !>            coefficients non nuls de la ligne i.
  !
  !> @details On ne connait pas le nombre de colonnes et on suppose
  !>          donc que la matrice est carrée
  !
  !> @details alloue de la mémoire une matrice n * p avec nnz
  !>          coefficients non nuls par ligne et prépare la structure
  !>          de donnée connaissant exactement le profil de la matrice
  !
  !> @author Yves Coudière
  !> - 30/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine csr_create_from_rows(s,lignes)

    type(csr),                 intent(inout) :: s
    integer,   dimension(:,:), intent(in)    :: lignes

    integer :: nnz, k, i, j, col, p

#ifdef DEBUG
    if (size(shape(lignes))/=2) then
       print "('csr_create_from_rows: the shape of lines is not (:,:)')"
       return
    end if
#endif

    ! Taille de la matrice : nombre de colonnes de _lignes_
    s%nl = size(lignes,2) ! Une ligne de matrice par colonne du
                          ! tableau lignes
    s%nc = s%nl           ! matrice supposée carrée
    p = size(lignes,1)
    ! Decompte du nombre d'elements non nuls
    nnz = 0
    do i = 1,s%nl
       nnz = nnz + count(lignes(:,i)>0) ! somme des indices positifs
    end do

#ifdef DEBUG
    print "('csr_create_from_rows: nombre d''éléments non nuls : ',I0)", nnz
#endif

    ! Allocation memoire
    allocate(s%row(s%nl+1), s%col(nnz), s%a(nnz))

    do i = 1, nnz
       s%a(i)   = 0._pr ! Les coefficients sont encore inconnus
    end do

    s%sorted = .false.

    ! Construction de la matrice
    k = 1
    do i = 1,s%nl
       s%row(i) = k
       do j = 1,p
          col = lignes(j,i)
          if (col>0) then ! On a un coeff non nul
             s%col(k) = col
             k = k+1
          end if
       end do
    end do
    s%row(s%nl+1) = k

#ifdef DEBUG
    print "('csr_create_from_rows:')"
    print "('nnz : ',I0,' -- nl : ',I0,' -- nc : ',I0)", nnz, s%nl, s%nc
    print "('row : ',10I7)", s%row(:)
    print "('col : ',10I7)", s%col(:)
    print "(' a  : ',10E7.1)", s%a(:)
#endif

  end subroutine csr_create_from_rows


  !> @brief Allocation initiale pour une matrice de taille n*p à
  !> partir du tableau voisins(:) des coefficients non nuls par ligne
  !
  !> @param[in,out] M matrice

  !> @param[in] voisins(:), où voisins(i) est la liste (type liste
  !>            chainée du module list) des coefficients non nuls de
  !>            laligne i
  !
  !> @details On ne connait pas le nombre de colonnes et on suppose
  !>          donc que la matrice est carrée
  !
  !> @details alloue de la mémoire une matrice n * p avec nnz
  !>          coefficients non nuls par ligne et prépare la structure
  !>          de donnée connaissant exactement le profil de la matrice
  !
  !> @author Yves Coudière
  !> - 30/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine csr_create_from_rows_list(s,voisins)

    use list
    type(csr),                  intent(inout) :: s
    type(i_list), dimension(:), intent(in)    :: voisins

    integer :: nnz, i, k, j

    ! Taille de la matrice
    s%nl = size(voisins) ! nombre de lignes
    s%nc = s%nl          ! matrice carrée
    ! Nombre d'elements non nuls
    nnz = 0
    do i = 1,s%nl
       nnz = nnz + get_length(voisins(i))
    end do

#ifdef DEBUG
    print "('csr_create_from_rows_list: nombre d''éléments non nuls : ',I0)", nnz
#endif

    ! Allocation memoire
    allocate(s%row(s%nl+1), s%col(nnz), s%a(nnz))

    do i = 1,nnz
       s%a(i)   = 0._pr
    end do

    s%sorted = .false.

    ! Construction de la matrice
    k = 1
    do i = 1,s%nl
       s%row(i) = k
       do j = 1,get_length(voisins(i)) ! Nombre de voisins de i
          s%col(k) = get_entry(voisins(i),j)
          k = k+1
       end do
    end do
    s%row(s%nl+1) = k

#ifdef DEBUG
    print "('csr_create_from_rows_list:')"
    print "('nnz : ',I0,' -- nl : ',I0,' -- nc : ',I0)", nnz, s%nl, s%nc
    print "('row : ',10I7)", s%row(:)
    print "('col : ',10I7)", s%col(:)
    print "(' a  : ',10E7.1)", s%a(:)
#endif

  end subroutine csr_create_from_rows_list


  !> @brief Allocation et remplissage pour une matrice csr à partir
  !>        d'une matrice list
  !
  !> @param[out] s matrice csr
  !> @param[in] s_list matrice de type m_list
  !
  !> @details alloue de la mémoire une matrice n * p avec nnz
  !>          coefficients non nuls par ligne et prépare la structure
  !>          de donnée connaissant exactement le profil de la matrice
  !>          et aussi les valeurs de coefficients non nuls
  !
  !> @author Yves Coudière
  !> - 30/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine csr_assign_from_m_list(s,s_list)

    use mat_list
    type(csr),    intent(out) :: s
    type(m_list), intent(in)  :: s_list

    integer :: nnz,k,i,j
    integer, dimension(:), allocatable :: col
    real(pr), dimension(:), allocatable :: a

    ! Dimensions
    s%nl = size(s_list,1)
    s%nc = size(s_list,2)
    nnz = get_nnz(s_list)

    ! Allocation memoire
    allocate(s%row(s%nl+1), s%col(nnz), s%a(nnz))

#ifdef DEBUG
    print "('csr_assign_from_m_list: nombre d''éléments non nuls : ',I0)", nnz
#endif

    ! On construit ligne par ligne
    k = 1
    do i = 1,s%nl
       s%row(i) = k
       call get_row(s_list,i, col,a) ! Renvoie la ligne comme deux tableaux
       do j = 1,size(col)
          s%col(k) = col(j)
          s%a(k)   = a(j)
          k = k+1
       end do
       deallocate(col,a)
    end do
    s%row(s%nl+1) = k

    s%sorted = .false.

#ifdef DEBUG
    print "('csr_assign_from_m_list:')"
    print "('nnz : ',I0,' -- nl : ',I0,' -- nc : ',I0)", nnz, s%nl, s%nc
    print "('row : ',10I7)", s%row(:)
    print "('col : ',10I7)", s%col(:)
    print "(' a  : ',7E10.2)", s%a(:)
#endif

  end subroutine csr_assign_from_m_list




!!$  !> @brief groupe 4 matrices en une matrice par blocs 2x2
!!$  !
!!$  !> @param[in] A11, A12, A21, A22 les 4 matrices à assembler
!!$  !> @param[out] B matrice groupée en sortie
!!$  !
!!$  !> @author Yves Coudière
!!$  !> @date 2 mai 2014
!!$  subroutine csr_assign_from_2x2_m_list(B, A11, A12, A21, A22)
!!$
!!$    use mat_list
!!$    type(m_list), intent(in) :: A11,A12,A21,A22
!!$    type(csr), intent(out) :: B
!!$
!!$    integer :: i, nnz,k,j, nl1,nl2, nc1,nc2
!!$    integer, dimension(:), allocatable :: col
!!$    real(pr), dimension(:), allocatable :: a
!!$
!!$    ! Nb de lignes et de colonnes
!!$    nl1 = size(A11,1)
!!$    nl2 = size(A21,1)
!!$    nc1 = size(A11,2)
!!$    nc2 = size(A12,2)
!!$    B%nl = nl1 + nl2
!!$    B%nc = nc1 + nc2
!!$    ! Nb elts on nuls
!!$    nnz = get_nnz(A11) + get_nnz(A12) + get_nnz(A21) + get_nnz(A22)
!!$
!!$    ! Allocation mémoire
!!$    allocate(B%row(B%nl+1), B%col(nnz), B%a(nnz))
!!$
!!$    ! Construction ligne par ligne
!!$    k = 1
!!$    ! Boucle sur les lignes de A11 et A12
!!$    do i = 1,nl1
!!$       B%row(i) = k              ! Début de la ligne i
!!$       ! Partie A_11
!!$       call get_row(A11,i,col,a) ! Renvoie la ligne comme deux tableaux
!!$       do j = 1,size(col)
!!$          B%col(k) = col(j)
!!$          B%a(k) = a(j)
!!$          k = k+1
!!$       end do
!!$       deallocate(col,a)
!!$       ! Partie A_12
!!$       call get_row(A12,i,col,a) ! Renvoie la ligne comme deux tableaux
!!$       do j = 1,size(col)
!!$          B%col(k) = nc1 + col(j)
!!$          B%a(k) = a(j)
!!$          k = k+1
!!$       end do
!!$       deallocate(col,a)
!!$    end do
!!$    ! Boucle sur les lignes de A21 et A22
!!$    do i = 1,nl2
!!$       B%row(nl1+i) = k              ! Début de la ligne i
!!$       ! Partie A_21
!!$       call get_row(A21,i,col,a) ! Renvoie la ligne comme deux tableaux
!!$       do j = 1,size(col)
!!$          B%col(k) = col(j)
!!$          B%a(k) = a(j)
!!$          k = k+1
!!$       end do
!!$       deallocate(col,a)
!!$       ! Partie A_12
!!$       call get_row(A22,i,col,a) ! Renvoie la ligne comme deux tableaux
!!$       do j = 1,size(col)
!!$          B%col(k) = nc1 + col(j)
!!$          B%a(k) = a(j)
!!$          k = k+1
!!$       end do
!!$       deallocate(col,a)
!!$    end do
!!$    B%row(B%nl+1) = k
!!$
!!$    B%sorted = .false.
!!$
!!$#ifdef DEBUG
!!$    print "('csr_assign_from_2x2_m_list:')"
!!$    print "('nnz : ',I0,' -- nl : ',I0,' -- nc : ',I0)", nnz, B%nl, B%nc
!!$    print "('row : ',10I7)", B%row(:)
!!$    print "('col : ',10I7)", B%col(:)
!!$    print "(' a  : ',10E7.1)", B%a(:)
!!$#endif
!!$
!!$  end subroutine csr_assign_from_2x2_m_list


  !> @brief Allocation et remplissage pour une matrice csr à partir
  !>        d'une matrice csr (copie de toute la structure de donnée)
  !
  !> @param[out] s matrice csr
  !> @param[in] a matrice de type m_list
  !
  !> @details alloue de la mémoire une matrice n * p avec nnz
  !>          coefficients non nuls par ligne et prépare la structure
  !>          de donnée connaissant exactement le profil de la matrice
  !>          et aussi les valeurs de coefficients non nuls
  !
  !> @author Yves Coudière
  !> - 4/7/2013 : réécriture depuis la bibliothèque spmlib
  subroutine csr_assign_from_csr(s,a)

    type(csr), intent(in)    :: a
    type(csr), intent(out)   :: s

    integer :: nnz,i

    s%nl = a%nl
    s%nc = a%nc
    s%sorted = a%sorted
    nnz = size(a%col)
    allocate(s%row(s%nl+1), s%col(nnz), s%a(nnz))

    do i = 1,s%nl+1
       s%row(i) = a%row(i)
    end do

    do i = 1,nnz
       s%col(i) = a%col(i)
       s%a(i)   = a%a(i)
    end do

  end subroutine csr_assign_from_csr

  !> @brief ajoute ou affecte un coefficient : M_ij [+]= a
  !
  !> @details si opt vaut mat_set alors c'est une affectation, si opt
  !>          vaut mat_add alors on ajoute : @f$ M(i,j) = M(i,j) +
  !>          a\f$.
  !>
  !>          décale les tableaux s%col et si nécessaire pour créer un
  !>          nouveau coefficient non nul
  !>
  !>          erreur si la place réservée dans col et a n'est pas
  !>          suffisante
  !
  !> @param[in,out] M matrice
  !> @param[in] i,j indices du coefficient dans la matrice
  !> @param[in] a valeur du coefficient
  !> @param opt (optionel) prend la valeur mat_add (défaut) pour
  !>             ajouter le coefficient à la matrice et mat_set pour
  !>             affecter le coefficient à la matrice.
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine construct_entry_csr (s, i, j, val, opt)

    type (csr), intent(inout) :: s     ! Matrice
    integer,    intent(in)    :: i, j  ! coordonnées
    real(pr),   intent(in)    :: val   ! valeur
    integer,    optional      :: opt   ! add[défaut] ou set

    integer :: k, m, actual_opt

#ifdef DEBUG
    if (.not.is_allocated_csr(s)) then
       print "('construct_entry_csr: allocation error')"
       return
    elseif (i<1.or.i>s%nl) then
       print "('construct_entry_csr: index out of bounds, i = ',I0)", i
       return
    elseif (j<1.or.j>s%nc) then
       print "('construct_entry_csr: index out of bounds, j = ',I0)", j
       return
    end if
#endif

    actual_opt = mat_add        ! valeur par défaut
    if (present(opt))  actual_opt = opt ! si opt est passé comme argument


    ! s%row est initialisé à (1...1) pour une matrice nulle

    ! On cherche si l'élément existe déjà
    ! -----------------------------------
    do k = s%row(i),s%row(i+1)-1
       if (s%col(k)==j) then ! C'est gagné, on est sur l'élément i,j
          select case(actual_opt)
          case(mat_add)
             s%a(k) = s%a(k) + val
          case(mat_set)
             s%a(k) = val
          case default
             s%a(k) = s%a(k) + val
          end select
          return ! On sort de la procédure
       end if
    end do

    ! L'élément n'existe pas. Il faut le créer et donc faire un décalage
    ! ------------------------------------------------------------------
    ! On verifie la place disponible
    if (s%row(s%nl+1)>size(s%col)) then
       print "('construct_entry_csr : plus de place dans le tableau&
            & reserve pour la matrice')"
       return
    end if
    ! On peut maintenant décaler
    do m = s%row(s%nl+1),s%row(i+1)+1,-1
       s%col(m) = s%col(m-1)
       s%a(m)   = s%a(m-1)
    end do
    ! Et ajouter le coefficient à sa place (penser que row(i+1)=0 si
    ! toutes les lignes sont vides)
    s%col(s%row(i+1)) = j
    s%a  (s%row(i+1)) = val ! add ou set, ça revient au même
    ! Mise a jour du tableau row
    s%row(i+1:s%nl+1) = s%row(i+1:s%nl+1) + 1

    s%sorted = .false.

  end subroutine construct_entry_csr

  !> @brief ajoute des éléments diagonaux à la matrice
  !
  !> @details affectation suivante: \f$ s(is_i,is_i) = s(is_i,is_i) +
  !> vals(i),\quad i = 1\ldots p \f.
  !
  !> @param[in,out] s la matrice
  !> @param[in] is(:) les indices de lignes sur lesquels ajouter des
  !>            valeurs
  !> @param[in] vals(:) les valeurs à ajouter
  !> @param opt (optionel) prend la valeur mat_add (défaut) pour
  !>             ajouter le coefficient à la matrice et mat_set pour
  !>             affecter le coefficient à la matrice.
  !
  !> @author Yves Coudière
  !> - 29/4/2013 : réécriture depuis la bibliothèque spmllib
  !> - 29/4/2013 : suppression de l'indicateur d'erreur ierr
  subroutine add_diag_csr (s, is, vals, opt)

    type (csr), intent(inout)            :: s    ! Matrice
    integer,    intent(in), dimension(:) :: is   ! Indices
    real(pr),   intent(in), dimension(:) :: vals ! Valeurs
    integer,    optional                 :: opt

    integer :: i,k, actual_opt

#ifdef DEBUG
    if (.not.is_allocated_csr(s)) then
       print "('add_diag_csr: allocation error')"
       return
    elseif (minval(is)<1) then
       print "('add_diag_csr: index out of bounds, q = ',I0)", minval(is)
       return
    elseif (maxval(is)>min(s%nl,s%nc)) then
       print "('add_diag_csr: index out of bounds, q = ',I0)", maxval(is)
       return
    elseif (size(vals)/=size(is)) then
       print "('add_diag_csr: sizes of indices and values differs')"
       return
    end if
#endif

    actual_opt = mat_add               ! valeur par défaut
    if (present(opt)) actual_opt = opt ! si opt est passé comme argument

    do k = 1,size(is)
       i = is(k) ! On travaille sur la ligne i, colonne i
       call construct_entry_csr(s, i,i, vals(k), actual_opt)
    end do

  end subroutine add_diag_csr


  !> @brief ajoute une sous matrice à la matrice
  !
  !> @details affectation suivante: \f$ s(is_k,js_l) = s(is_k,is_l) +
  !> vals(k,l), \quad k = 1\ldots n,\ l = 1 \ldots p \f$
  !
  !> @param[in,out] s la matrice
  !> @param[in] is(:), js(:) les indices de lignes et colonnes sur
  !>            lesquels ajouter des valeurs
  !> @param[in] vals(:,:) les valeurs à ajouter
  !> @param opt (optionel) prend la valeur mat_add (défaut) pour
  !>             ajouter le coefficient à la matrice et mat_set pour
  !>             affecter le coefficient à la matrice.
  !
  !> @author Yves Coudière
  !> - 29/4/2013 : réécriture depuis la bibliothèque spmllib
  !> - 29/4/2013 : suppression de l'indicateur d'erreur ierr
  subroutine add_mat_csr (s, is, js, vals, opt)

    type (csr), intent(inout)              :: s
    integer,    intent(in), dimension(:)   :: is, js
    real(pr),   intent(in), dimension(:,:) :: vals
    integer,    optional                   :: opt

    integer :: k, m, actual_opt

#ifdef DEBUG
    if (.not.is_allocated_csr(s)) then
       print "('add_mat_csr: allocation error')"
       return
    elseif (minval(is)<1.or.maxval(is)>s%nl.or.&
         minval(js)<1.or.maxval(js)>s%nc) then
       print "('add_mat_csr: index out of bounds, q = ',I0)", q
       return
    elseif (size(vals,1)/=size(is).or.size(vals,2)/=size(js)) then
       print "('add_mat_csr: sizes of indices and values differs')"
       return
    end if
#endif

    actual_opt = mat_add               ! valeur par défaut
    if (present(opt)) actual_opt = opt ! si opt est passé comme argument

    ! On procède élément par élément
    do k = 1,size(is)
       do m = 1,size(js)
          call construct_entry_csr(s,is(k),js(m),vals(k,m), actual_opt)
       end do
    end do

  end subroutine add_mat_csr

  ! ************************************************************************
  ! *                                                                      *
  ! *  Méthodes "opérateurs unitaires"                                     *
  ! *                                                                      *
  ! ************************************************************************

  !> @brief transpose une matrice
  !
  !> @param[in,out] M la matrice
  !
  !> @author Yves Coudière
  !> - 29/4/2013 : réécriture depuis la bibliothèque spmllib
  subroutine transpose_csr(M)

    use list
    type(csr), intent(inout) :: M

    ! M est de taille nl x nc et sera de taille nc x nl
    type(i_list), dimension(M%nc) :: lignes
    type(f_list), dimension(M%nc) :: a

    integer :: i,k,j
    integer, dimension(M%nc+1) :: Mrow

    ! On construit les lignes de la transposée en parcourant les
    ! lignes de la matrice M
    do i = 1,M%nl
       do k = M%row(i),M%row(i+1)-1 ! Parcours de la ligne i de M
          call append(lignes(M%col(k)), i) ! lignes(col) := [ lignes(col), i ]
          call append(a(M%col(k)), M%a(k)) ! coeff(col) := [ coeff(col), a ]
       end do
    end do

    ! On remplace M par sa transposée, qui occupe autant de place en
    ! mémoire (donc pas de réallocation).
    k = 1
    do i = 1,M%nc
       Mrow(i) = k
       do j = 1,get_length(lignes(i)) ! Nombre de coeff non nuls de la ligne i
          M%col(k) = get_entry(lignes(i),j)
          M%a(k) = get_entry(a(i),j)
          k = k+1
       end do
    end do

    ! on remplace le vecteur M%row
    Mrow(M%nc+1) = k
    deallocate(M%row)
    allocate(M%row(M%nc+1))
    do i = 1, M%nc+1
       M%row(i)=Mrow(i)
    end do

    M%sorted = .false.

    ! On libere la mémoire intermédiaire
    do i = 1,M%nc
       call deallocate(lignes(i))
       call deallocate(a(i))
    end do

    ! échange nl <-> nc
    i = M%nl
    M%nl = M%nc
    M%nc = i

  end subroutine transpose_csr

  ! csr_permute

  ! ************************************************************************
  ! *                                                                      *
  ! *  Méthodes "opérateurs binaires"                                     *
  ! *                                                                      *
  ! ************************************************************************

  ! csr_add_lambda_id
  ! csr_times_scalar


  !> @brief Computes the matrix-vector product v=s*u. This is an
  !> _orphaned_ OMP subroutine.

  !> @param[in] s matrix
  !> @param[in] u vector that multiply the matrix
  !> @param[out] v output vector
  function csr_times_array(s,u) result(v)

    type(csr),                 intent(in)  :: s
    real(pr), dimension(s%nc), intent(in)  :: u
    real(pr), dimension(s%nl)              :: v ! := s*u

    ! Les variables locales sont _private_ par defaut.
    integer :: i, jDeb, jFin, j

#ifdef DEBUG
    !$OMP SINGLE
    if (size(shape(u))/=1.or.size(shape(v))/=1.or.&
         size(u,1)/=s%nc.or.size(v,1)/=s%nl) then
       print *,'In axCSR: dimension mismatch'
    end if
    !$OMP END SINGLE
    ! Il y a une barriere implicite en fin de clause SINGLE
#endif

    ! Multiplication ligne par ligne
    !$OMP DO SCHEDULE(RUNTIME)
    do i = 1,s%nl
       jDeb = s%row(i)     ! Debut du stockage de la ligne i
       jFin = s%row(i+1)-1 ! Fin du stockage de la ligne i
       v(i) = 0._pr
       do j = jDeb,jFin
          v(i) = v(i) + s%a(j)*u(s%col(j))
       end do
    end do
    !$OMP END DO

  end function csr_times_array


  ! ************************************************************************
  ! *                                                                      *
  ! *  Méthodes interrogation de la matrice                                *
  ! *                                                                      *
  ! ************************************************************************

  !> @brief Get the number of rows
  integer function get_nl(S)
    type(csr), intent(in) :: S
    get_nl = S%nl
  end function get_nl

  !> @brief Get the number of columns
  integer function get_nc(S)
    type(csr), intent(in) :: S
    get_nc = S%nc
  end function get_nc

  !> @brief Get the number nonzeros elements
  integer function get_nnz_csr(S)
    type(csr), intent(in) :: S
    get_nnz_csr = S%row(S%nl+1)-1
  end function get_nnz_csr

  ! ************************************************************************
  ! *                                                                      *
  ! *  Méthodes Calcul de matrices types                                   *
  ! *                                                                      *
  ! ************************************************************************


  !> @brief Calcul de la matrice DF du pb de poisson, u -->
  !> -Laplacien(u) sur le carré (0,1)x(0,1) avec n divisions dans
  !> chaque direction et conditions de Dirichlet.

  !> @details La matrice est de taille (n-1)^2 et ses coefficients sont
  !>
  !> A = 1/h^2* *(a_{ij}) avec i,j = 1 .. n-1
  !>
  !> et a_{ii} = 4
  !>    a_{ii+1}=a_{ii-1}=a_{ii-(n-1)}=a_{ii+(n-1)}=-1
  !>
  !> avec h = 1/n.
  function matrice_poisson2D(n) result(A)
    type(csr) :: A
    integer, intent(in) :: n

    real(pr) :: h2
    integer :: i,j,ij,k, nnz

    ! h^2
    h2 = 1._pr/real(n,pr)**2

    ! Allocation memoire
    A%nl = (n-1)**2
    A%nc = (n-1)**2
    nnz = 5*(n-1)**2 - 4*(n-1)
    allocate(A%row(A%nl+1), A%col(nnz), A%a(nnz))

    ! Initialisation
    A%row = 0
    A%col = 0
    A%a = 0._pr

    ! Construction du 'profil' ou 'graphe' de la matrice
    ! et remplissage de la matrice
    k = 1
    do j = 1,n-1
       do i = 1,n-1
          ij = i+(j-1)*(n-1)
          A%row(ij) = k ! Indice du debut de stockage de la ligne 'ij'
          if (j>1) then
             A%col(k) = ij-(n-1)
             A%a  (k) = -1._pr / h2
             k = k+1
          end if
          if (i>1) then
             A%col(k) = ij-1
             A%a  (k) = -1._pr / h2
             k = k+1
          end if
          A%col(k) = ij
          A%a  (k) = 4._pr / h2
          k = k+1
          if (i<n-1) then
             A%col(k) = ij+1
             A%a  (k) = -1._pr / h2
             k = k+1
          end if
          if (j<n-1) then
             A%col(k) = ij+(n-1)
             A%a  (k) = -1._pr / h2
             k = k+1
          end if
       end do
    end do
    A%row(A%nl+1) = k

    A%sorted = .true.

#ifdef DEBUG
    print *,'k :',k
    print *,'nnz :',nnz,' nc :', A%nc,' nl :', A%nl
    print *,'ROW :',A%row
    print *,'COL :',A%col
    print *,' A  :',A%a
    print *,'DIAG:',A%diag
#endif

  end function matrice_poisson2D

  ! ************************************************************************
  ! *                                                                      *
  ! *  Méthodes "I/O"                                             *
  ! *                                                                      *
  ! ************************************************************************

  ! csr_save
  ! csr_load


  !> @brief affiche la matrice
  !
  !> @param[in] S matrice
  !
  !> @author Yves Coudière
  !> - 14/4/2014 : réécriture depuis la bibliothèque spmlib
  subroutine print_csr(S)
    type(csr), intent(in) :: s

    integer :: i,j,jDeb,jFin

    print "('Matrice <csr> de taille ',2I9)", s%nl,s%nc
    ! Affichage ligne par ligne
    do i = 1,s%nl
       jDeb = s%row(i)     ! Debut du stockage de la ligne i
       jFin = s%row(i+1)-1 ! Fin du stockage de la ligne i
       do j=jDeb,jFin
          print "('(',I7,',',I7,') --> ',E13.5)", i,s%col(j),s%a(j)
       end do
    end do
!!$    write (unit,"('Column numering is sorted ? ',L1)") s%sorted

!!$    case(MAT_BIN)
!!$
!!$       ! Affichage
!!$       write (unit) s%nl,s%nc
!!$       write (unit) s%row
!!$       write (unit) s%col
!!$       write (unit) s%a
!!$       write (unit) associated(s%diag)
!!$       if (associated(s%diag)) write (unit) s%diag
!!$       write (unit) s%sorted
!!$
!!$    case(MAT_MATLAB)
!!$
!!$       ! Affichage ligne par ligne
!!$       do i = 1,s%nl
!!$          jDeb = s%row(i)     ! Debut du stockage de la ligne i
!!$          jFin = s%row(i+1)-1 ! Fin du stockage de la ligne i
!!$          if (jDeb<=jFin) &
!!$               write (unit,"(2I8,E14.6)") (i,s%col(j),s%a(j), j=jDeb,jFin)
!!$       end do
!!$
!!$    end select

  end subroutine print_csr


end module mat_csr
