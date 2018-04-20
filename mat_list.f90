!> @brief Matrices creuses sous forme de tableaux de listes
!> @date Time-stamp: <2017-10-22 13:11:54 yves>
!
!> @details Module de matrices creuses utilisant des listes chainées.
!> C'est le module de base pour remplir des matrices de
!> discrétisation.  Une matrice peut etre rectangulaire (taille n par
!> p)
!
!> @note Pour le calcul il faut ensuite les convertir en un format
!> utile : CSR, CSC...
!
!> @todo
!>  - Mettre à jour les auteurs
!>  - Tester l'implémenter l'opération A := alpha*A + beta*B
!
!> @author Y. Coudière
!> - Code démarré vers 2006
!> @date mai 2013 mise à plat:
!> - changement des noms des procédures (notation avec des _)
!> - documentation doxygen
!> - suppression de tout ce qui semble inutile, quitte à retourner le
!>   chercher plus tard dans la librairie spmlib originale
module mat_list

  use prec
  use matrices
  implicit none

  private

  !> @brief noeud des listes utilisées pour les lignes: couple (j,a)
  type mat_node
     private
     integer :: j = 0  ! indice de colonne
     real(pr) :: a = 0._pr ! coefficient
     type(mat_node), pointer :: next => NULL()
  end type mat_node

  !> @brief structure de donnée pour une ligne de matrice: une liste
  !>        de couples
  type mat_ligne
     private
     integer :: length = 0
     type(mat_node), pointer :: first => NULL()
     type(mat_node), pointer :: last => NULL()
  end type mat_ligne

  !> @brief Une matrice est un tableau des lignes. Chaque ligne est
  !>        une liste de couples (j,a_{ij}) d'éléments non-nuls dans
  !>        la matrice. Les membres du type m_list sont publics.
  type m_list
     private
     integer :: nl = 0 ! nombre de lignes
     integer :: nc = 0 ! nombre de colonnes
     type(mat_ligne), dimension(:), allocatable :: ligne
  end type m_list

  !> @brief libère la mémoire de la matrice
  interface deallocate
     module procedure deallocate_m_list
  end interface deallocate

  !> @brief affiche la matrice
  interface print
     module procedure print_m_list
  end interface print

  !> @brief allocation mémoire pour une matrice nulle de taille n*p
  interface zeros
     module procedure zeros_m_list
  end interface zeros

  !> @brief allocation mémoire pour une matrice nulle de taille n*p
  interface allocate
     module procedure zeros_m_list
  end interface allocate

  !> @brief allocation mémoire et remplissage pour une matrice à
  !>        partir des ses diagonales non nulles
  interface diags
     module procedure diags_m_list, diags_m_list_array
  end interface diags

  !> @brief construction de la matrice par ajout ou assignation de
  !>        coefficients
  interface mat_add_value
     module procedure construct_entry_m_list, construct_array_m_list
  end interface mat_add_value

  !> @brief nombre de coefficients non nuls
  interface get_nnz
     module procedure nnz_m_list
  end interface get_nnz

  !> @brief renvoie un coefficient
  interface get_entry
     module procedure get_entry_m_list
  end interface get_entry

  !> @brief renvoie une ligne entière
  interface get_row
     module procedure get_row_m_list
  end interface get_row

  !> @brief renvoie une colonne entière
  interface get_col
     module procedure get_col_m_list
  end interface get_col

  !> @brief renvoie la norme euclidienne d'une ligne
  interface row_norm
     module procedure row_norm_m_list
  end interface row_norm

  !> @brief renvoie la norme infini d'une matrice
  interface norm
     module procedure norm_m_list
  end interface norm

  !> @brief renvoie le nombre de coefficients non nuls d'une ligne
  interface row_length
     module procedure row_length_m_list
  end interface row_length

  !> @brief renvoie la taille de la matrice
  interface size
     module procedure size_m_list, size_k_m_list
  end interface size

  !> @brief recopie la matrice, A(:,:) = B(:,:)
  interface assignment(=)
     module procedure copy_m_list
  end interface assignment(=)

  public :: m_list                       ! structure de donnée publique
  public :: deallocate                   ! destructeur
  public :: zeros, diags, mat_add_value  ! constructeurs
  public :: allocate, assignment(=)      ! constructeurs
  public :: alpha_a_plus_beta_b          ! constructeurs
  public :: row_length, get_nnz          ! méthodes informatives
  public :: get_entry
  public :: get_row, row_norm, size, norm, get_col
  public :: print

contains

  !> @brief revoie la taille de la matrice
  !>
  !> @param[in] M la matrice
  !> @return size, tableau 2*2 avec (n,p) pour une matrice n*p
  !
  !> @author Yves Coudière
  !> - 30/4/2013 : première version
  function size_m_list(M) result (s)

    type(m_list), intent(in) :: M
    integer, dimension(2) :: s

    s(:) = (/ M%nl, M%nc /)

  end function size_m_list

  !> @brief revoie la taille de la matrice dans la direction k=1 ou 2
  !>
  !> @param[in] M la matrice
  !> @param[in] k = 1 ou 2
  !> @return size taille dans la direction 1 ou 2
  !
  !> @author Yves Coudière
  !> - 30/4/2013 : première version
  function size_k_m_list(M,k) result (s)

    type(m_list), intent(in) :: M
    integer, intent(in) :: k
    integer :: s

    select case(k)
    case(1)
       s = M%nl
    case(2)
       s = M%nc
    case default
       s = 0
       print "('size_k_m_list: la dimension ',I0,' n''existe pas pour cette matrice')", k
  end select

  end function size_k_m_list

  !> @brief renvoie le nombre de coefficients non nuls de la ligne i
  !
  !> @param[in] M matrice
  !> @param[in] i indice de ligne
  !> @return le nombre de coefficients non nuls de la ligne i de la
  !>         matrice M
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  function row_length_m_list(M,i) result(p)

    type(m_list), intent(in) :: M
    integer, intent(in) :: i
    integer :: p

    p = M%ligne(i)%length

  end function row_length_m_list

  !> @brief renvoie le nombre total de coefficients non nuls de la
  !>        matrice
  !
  !> @param[in] M matrice
  !> @return le nombre total de coefficients non nuls de la matrice
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  !> @date 28/4/2013 -- documentation
  function nnz_m_list(M) result(nnz)

    type(m_list), intent(in) :: M
    integer :: nnz

    integer :: i

    nnz = 0
    do i = 1,M%nl
       nnz = nnz + M%ligne(i)%length
    end do

  end function nnz_m_list

  !> @brief libère la mémoire pour une ligne
  !
  !> @param[in,out] Mi ligne de matrice (ie liste chainée des entrées
  !>               de la ligne)
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine deallocate_ligne(Mi)

    type(mat_ligne), intent(inout) :: Mi ! Ligne d'une matrice

    type(mat_node), pointer :: node, nextNode

    if (associated(Mi%first)) then
       node => Mi%first
       nextNode => node%next
       deallocate(Mi%first)
       do while (associated(nextNode))
          node => nextNode
          nextNode => node%next
          deallocate(node)
       end do
    end if

    nullify(Mi%first)
    nullify(Mi%last)
    Mi%length = 0

  end subroutine deallocate_ligne

  !> @brief libère la mémoire pour la matrice
  !
  !> @param[in,out] M matrice
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine deallocate_m_list(M)

    type(m_list), intent(inout) :: M ! Matrice

    integer :: i

    do i = 1,M%nl
       call deallocate_ligne(M%ligne(i)) ! Libère la mémoire de la ligne i
    end do

    !deallocate(M%ligne)
    M%nl = 0
    M%nc = 0

  end subroutine deallocate_m_list

  !> @brief affiche la matrice
  !
  !> @param[in] M matrice
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine print_m_list(M)

    type(m_list), intent(in) :: M

    type(mat_node), pointer :: node
    integer :: i

    print "('Matrice <m_list> de taille ',2I9)", M%nl,M%nc
    do i = 1,M%nl
       ! Lecture de la ligne i
       node => M%ligne(i)%first ! Premier noeud de la ligne
       do while (associated(node))
          print "('(',I7,',',I7,') --> ',E13.5)", i,node%j,node%a
          node => node%next ! Noeud suivant
       end do
    end do

  end subroutine print_m_list

  !> @brief Allocation initiale pour une matrice de taille n*p nulle
  !>        (sans coefficients)
  !
  !> @param[in,out] M matrice
  !> @param[in] n,p taille de la matrice
  !
  !> @details alloue de la mémoire pour un tableau de n listes
  !>          chainées pouvant accueillir les coefficients de la
  !>          matrice. Comme ces listes sont vides, la matrice est
  !>          nulle.
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine zeros_m_list(M,n,p)

    type(m_list), intent(inout) :: M
    integer, intent(in) :: n,p

    M%nl = n
    M%nc = p
    allocate(M%ligne(n))

  end subroutine zeros_m_list

  !> @brief fonction \a diags de matlab pour une construire une
  !>        matrice via ses diagonales non nulles
  !>
  !> @details construction d'une matrice à partir de ses diagonales,
  !>          comme dans la fonction Matlab:
  !>          M = spdiags(diags, vals)
  !>          si diags=[-10,-3,0,5,6,8] et vals = [1.,2.,3.,4.,5.,6.]
  !>          alors on fait s(i,i+diags) = vals(diags) pour 1 <=
  !>          i+diag <= p
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine diags_m_list(M,n,p,diags,vals)

    type(m_list), intent(out) :: M
    integer, intent(in) :: n,p
    integer, dimension(:), intent(in) :: diags
    real(pr), dimension(:), intent(in) :: vals

    integer :: i,j,k
    type(mat_node), pointer :: node

    if (allocated(M%ligne)) then
       print "('diags_m_list: matrice déjà allouée')"
       return
    end if

    M%nl = n
    M%nc = p
    allocate(M%ligne(n))

    ! boucle sur les diagonales
    do k = 1,size(diags)
       ! Boucle sur les lignes
       do i = 1,n
          j = i+diags(k)
          if (0<j.and.j<=p) then
             ! la colonne j rentre bien dans la matrice
             if (M%ligne(i)%length==0) then
                ! C'est le premier coefficient non nul de la ligne
                M%ligne(i)%length = 1
                allocate(M%ligne(i)%first)
                M%ligne(i)%last => M%ligne(i)%first
                M%ligne(i)%first%j = j       ! colonne j = i + diags(k)
                M%ligne(i)%first%a = vals(k) ! coeff de la diagonale k
             else
                ! On ajoute à la fin de la liste
                node => M%ligne(i)%last ! Dernier noeud de la liste
                M%ligne(i)%length = M%ligne(i)%length + 1
                allocate(node%next)
                M%ligne(i)%last => node%next
                node%next%j = j
                node%next%a = vals(k)
             end if
          end if
       end do
    end do

  end subroutine diags_m_list

  !> @brief fonction \a diags de matlab pour une construire une
  !>        matrice via ses diagonales non nulles
  !>
  !> @details construction d'une matrice à partir de ses diagonales,
  !>          comme dans la fonction Matlab:
  !>          M = spdiags(diags, vals)
  !>          si diags=[-10,-3,0,5,6,8] et vals = array(n,n_diags)
  !>          alors on fait: M(i,i+diags) = vals(i,diags) pour 1 <=
  !>          i+diag <= p
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine diags_m_list_array(M,n,p,diags,vals)

    type(m_list), intent(out) :: M
    integer, intent(in) :: n,p
    integer, dimension(:), intent(in) :: diags
    real(pr), dimension(:,:), intent(in) :: vals

    integer :: i,j,k
    type(mat_node), pointer :: node

    if (allocated(M%ligne)) then
       print "('diags_m_list: matrice déjà allouée')"
       return
    end if
    if (size(vals,2)/=n) then
       print "('diags_m_list: le nombre de ligne de vals n''est pas&
            & égal au nombre de lignes de la matrice')"
       return
    end if

    M%nl = n
    M%nc = p
    allocate(M%ligne(n))

    ! boucle sur les diagonales
    do k = 1,size(diags)
       ! Boucle sur les lignes
       do i = 1,n
          j = i+diags(k)
          if (0<j.and.j<=p) then
             ! la colonne j rentre bien dans la matrice
             if (M%ligne(i)%length==0) then
                ! C'est le premier coefficient non nul de la ligne
                M%ligne(i)%length = 1
                allocate(M%ligne(i)%first)
                M%ligne(i)%last => M%ligne(i)%first
                M%ligne(i)%first%j = j         ! colonne j = i + diags(k)
                M%ligne(i)%first%a = vals(i,k) ! coeff de la diagonale k
             else
                ! On ajoute à la fin de la liste
                node => M%ligne(i)%last ! Dernier noeud de la liste
                M%ligne(i)%length = M%ligne(i)%length + 1
                allocate(node%next)
                M%ligne(i)%last => node%next
                node%next%j = j
                node%next%a = vals(i,k)
             end if
          end if
       end do
    end do

  end subroutine diags_m_list_array

  !> @brief ajoute ou affecte un coefficient : M_ij [+]= a
  !
  !> @details si opt vaut mat_set alors c'est une affectation, si opt
  !>          vaut mat_add alors on ajoute : \f$ M(i,j) = M(i,j) +
  !>          a\f$
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
  subroutine construct_entry_m_list(M, i,j,a, opt)

    type(m_list), intent(inout) :: M
    integer, intent(in) :: i,j
    real(pr), intent(in) :: a
    integer, optional :: opt

    integer :: actual_opt
    type(mat_node), pointer :: node
    logical :: flag

    actual_opt = mat_add               ! valeur par défaut
    if (present(opt)) actual_opt = opt ! si opt est passé comme argument

    if (.not.allocated(M%ligne)) then
       print"('Pas de mémoire allouée pour la matrice')"
       return
    end if

    if (.not.associated(M%ligne(i)%first)) then
       ! Rien encore dans la ligne i
       ! ---------------------------
       M%ligne(i)%length = 1
       allocate(M%ligne(i)%first) ! Allocation mémoire pour le noeud
       M%ligne(i)%first%j = j
       M%ligne(i)%first%a = a
       M%ligne(i)%last => M%ligne(i)%first
    else
       ! Il y a déjà des coefficients dans la ligne i
       ! --------------------------------------------
       ! On regarde la colonne j s'y trouve déjà
       node => M%ligne(i)%first
       flag = (node%j==j) ! sommes-nous sur la colonne j ?
       do while (.not.flag)
          if (associated(node%next)) then
             node => node%next
             flag = (node%j==j) ! sommes-nous sur la colonne j ?
          else
             flag = .true. ! Pas de noeud suivant, on s'arrete
          end if
       end do
       if (node%j==j) then
          ! L'entrée (i,j) existe déjà
          ! --------------------------
          select case(actual_opt)
          case(mat_add)
             node%a = node%a + a
          case(mat_set)
             node%a = a
          case default
             node%a = node%a + a
          end select
       else
          ! L'entrée (i,j) n'existe pas encore
          ! ----------------------------------
          ! node est le dernier noeud, node%next => NULL
          M%ligne(i)%length = M%ligne(i)%length + 1
          allocate(node%next) ! Allocation mémoire pour l'élément suivant
          node%next%j = j
          node%next%a = a
          M%ligne(i)%last => node%next ! Dernier élément de la liste
       end if
    end if

  end subroutine construct_entry_m_list

  !> @brief ajoute ou affecte une sous-matrice
  !
  !> @details si opt vaut mat_set alors c'est une affectation, si opt
  !> vaut mat_add alors on ajoute. L'opération par défaut est donc \f$
  !> M(i_k,j_l) = M(i_k,j_l) + a(k,l),\quad k=1\ldots n_i,\ l =
  !> 1\ldots n_j \f$
  !
  !> @param[in,out] M matrice
  !> @param[in] i(:),j(:) indices de ligne et colonne des coefficients
  !> @param[in] a(:,:) valeurs des coefficient. si i(:) et j(:) sont
  !>            de tailles n_i e n_j alors a(:,:) doit être de taille
  !>            n_i * n_j
  !> @param[opt] opt (optionel) prend la valeur mat_add (défaut) pour
  !>             ajouter le coefficient à la matrice et mat_set pour
  !>             affecter le coefficient à la matrice.
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine construct_array_m_list(M,i,j,a,opt)

    type(m_list), intent(inout) :: M
    integer, dimension(:), intent(in) :: i, j
    real(pr), dimension(:,:), intent(in) :: a
    integer, optional :: opt

    integer :: k,l, actual_opt

    actual_opt = mat_add               ! valeur par défaut
    if (present(opt)) actual_opt = opt ! si opt est passé comme argument

    do k = 1,size(i)
       do l = 1,size(j)
          call construct_entry_m_list(M,i(k),j(l),a(k,l),actual_opt)
       end do
    end do

  end subroutine construct_array_m_list

  !> @brief renvoie le coefficient M(i,j)
  !
  !> @param[in] M matrice
  !> @param[in] i,j indices
  !> @return le coefficient \f$M(i,j)\f$
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  function get_entry_m_list(M,i,j) result(x)

    type(m_list), intent(in) :: M
    integer, intent(in) :: i,j
    real(pr) :: x

    type(mat_node), pointer :: node => NULL()

    ! Lecture de la ligne i
    x = 0._pr
    node => M%ligne(i)%first
    do while (associated(node))
       if (node%j==j) then
          x = node%a
          exit
       end if
       node => node%next
    end do

  end function get_entry_m_list

  !> @brief renvoie la ligne M(i,:)
  !
  !> @param[in] M matrice
  !> @param[in] i indice de ligne
  !> @param[out] j(:) indices des colonnes ayant des coefficients dans
  !>             M(i,:)
  !> @param[out] a(:) coefficients de M(i,:) dans l'ordre des j(:)
  !
  !> @details les tableaux j et a sont allocatable et sont alloués au
  !>          cours de la procédure
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  !> @date 28/4/2013 -- arguments j et a passent de pointer à
  !>                    allocatable
  subroutine get_row_m_list(M,i,j,a)

    type(m_list), intent(in) :: M
    integer, intent(in) :: i
    integer, dimension(:), allocatable, intent(out) :: j
    real(pr), dimension(:), allocatable, intent(out) :: a

    integer :: k
    type(mat_node), pointer :: node => NULL()

    ! Allocation
    allocate(j(M%ligne(i)%length), a(M%ligne(i)%length))

    ! Lecture de la ligne i
    k = 0
    node => M%ligne(i)%first
    do while (associated(node))
       k = k+1
       j(k) = node%j
       a(k) = node%a
       node => node%next
    end do

  end subroutine get_row_m_list


  !> @brief renvoie la colonne M(:,j)
  !
  !> @param[in] M matrice
  !> @param[in] j indice de colonne
  !> @param[out] i(:) indices des colonnes ayant des coefficients dans
  !>             M(i,:)
  !> @param[out] a(:) coefficients de M(:,j) dans l'ordre des i(:)
  !
  !> @details les tableaux i et a sont allocatable et sont alloués au
  !>          cours de la procédure
  !
  !> @author Yves Coudière
  subroutine get_col_m_list(M,j,is,as)

    type(m_list), intent(in) :: M
    integer, intent(in) :: j
    integer, dimension(:), allocatable, intent(out) :: is
    real(pr), dimension(:), allocatable, intent(out) :: as

    integer :: k, i
    type(mat_node), pointer :: node => NULL()

    ! Count the number of nonzeros in the column j
    k = 0
    do i = 1,M%nl ! loop through the rows
       node => M%ligne(i)%first ! search for any occurence of the column j
       do while (associated(node))
          if (node%j == j) then
             k = k + 1 ! (i,j) is nonzero, exit the search loop
             exit
          end if
          node => node%next ! go to next node
       end do
    end do

    ! Allocate
    allocate(is(k), as(k))

    ! Actually get column j
    k = 0
    do i = 1,M%nl ! loop through the rows
       node => M%ligne(i)%first ! search for any occurence of the column j
       do while (associated(node))
          if (node%j == j) then
             k = k + 1 ! (i,j) is nonzero, exit the search loop
             is(k) = i
             as(k) = node%a
             exit
          end if
          node => node%next ! go to next node
       end do
    end do

  end subroutine get_col_m_list


  !> @brief norme euclidienne de la ligne M(i,:)
  !
  !> @param[in] M matrice
  !> @param[in] i indice de ligne
  !> @return la norme euclidienne de la ligne, soit \f$ \left(
  !>         \sum_{j=1}^p M(i,j)^2 \right)^{1/2} \f$
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  function row_norm_m_list (M,i) result(x)

    type(m_list), intent(in) :: M
    integer, intent(in) :: i
    real(pr) :: x

    type(mat_node), pointer :: node => NULL()

    x = 0._pr ! valeur si aucune entré dans la ligne

    ! Lecture de la ligne i
    node => M%ligne(i)%first
    do while (associated(node))
       x = x + (node%a)**2
       node => node%next
    end do

    x = sqrt(x)

  end function row_norm_m_list

  !> @brief renvoie la norme infini d'un matrice
  !
  !> @details C'est \f$ \|M\|_\infty = \max_{i=1\ldots n} \sum_{j=1}^n
  !> |a_{ij}| \f$.
  !> @param[in] M matrice
  !> @return la norme \f$ \|M\|_\infty = \max_{i=1\ldots n}
  !> \sum_{j=1}^n |a_{ij}| \f$.
  !> @author Yves Coudière
  !> @date 09/04/2014
  real(pr) function norm_m_list(M) result (x)
    type(m_list), intent(in) :: M

    type(mat_node), pointer :: node => NULL()
    integer :: i
    real(pr) :: x0

    x = 0._pr

    do i = 1,M%nl
       ! Lecture de la ligne i et construction de x0 = sum_j |a_{ij}|.
       x0 = 0._pr
       node => M%ligne(i)%first
       do while (associated(node))
          x0 = x0 + abs(node%a)
          node => node%next
       end do
       nullify(node)
       x = max(x,x0)
    end do

  end function norm_m_list

  !> @brief recopie une matrice A = B <=> "call operator(=) (A,B)"
  !
  !> @param[out] A matrice copie de A
  !> @param[in] B matrice
  !
  !> @note recopie toute la structure de données
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine copy_m_list(A,B) ! operation A:=B (recopie de B dans A)

    type(m_list), intent(out)   :: A
    type(m_list), intent(in)    :: B

    integer :: i
    type(mat_node), pointer :: node => NULL()

    ! allocation mémoire pour la matrice A
    call zeros_m_list(A,B%nl,B%nc)

    ! recopie ligne par ligne
    do i = 1,B%nl
       ! Lecture de la ligne i
       node => B%ligne(i)%first ! Premier noeud de la ligne
       do while (associated(node))
          call construct_entry_m_list(A,i,node%j,node%a)
          node => node%next ! Noeud suivant
       end do
    end do

  end subroutine copy_m_list

  !> @brief Opération \f$ A := \alpha A + \beta B\f$.
  !
  !> @author Yves Coudière
  !> @param[in] alpha coefficient multiplicatif de A en entrée
  !> @param[in,out] A matrice en entrée et modifiée en sortie
  !> @param[in] beta coefficient multiplicatif de B
  !> @param[in] B matrice
  subroutine alpha_a_plus_beta_b(alpha,A,beta,B)
    real(pr),intent(in) :: alpha
    type(m_list), intent(inout) :: A
    real(pr), intent(in) :: beta
    type(m_list), intent(in) :: B

    integer :: i
    type(mat_node), pointer :: node => NULL()

    if (A%nl/=B%nl.or.A%nc/=B%nc) then
       print "('mat_ligne: alpha_a_plus_beta_b -- &
            &les matrices sont de tailles différentes.')"
       stop
    end if

    ! D'abord on calcule alpha*A
    do i = 1,A%nl
       ! Lecture de la ligne i
       node => A%ligne(i)%first ! Premier noeud de la ligne
       do while (associated(node))
          node%a = alpha * node%a
          node => node%next ! Noeud suivant
       end do
    end do

    ! On ajoute beta*B
    do i = 1,B%nl
       ! Lecture de la ligne i de B
       node => B%ligne(i)%first ! Premier noeud de la ligne
       do while (associated(node))
          ! A_{ij} = A_{ij} + beta*B_{ij}
          call mat_add_value(A, i,node%j, beta*node%a)
          node => node%next ! Noeud suivant
       end do
    end do

  end subroutine alpha_a_plus_beta_b

end module mat_list
