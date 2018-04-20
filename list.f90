!> @brief Implémentation simple de quelques listes (entiers et réels)
!> @date Time-stamp: <2016-06-02 15:47:22 yves>
!
!> @details Implémentation de listes chainées pour des entiers et des réels
!
!> @todo Mettre à jour les auteurs
!
!> @author Y. Coudière
!> - Code démarré vers 2006
!> @date mai 2013 mise à plat:
!> - changement des noms des procédures (notation avec des _)
!> - documentation doxygen
!> - suppression de tout ce qui semble inutile, quitte à retourner le
!>   chercher plus tard dans la librairie spmlib originale
module list

  use prec
  implicit none
  
  private

  !> @brief structure de donnée pour une liste d'entiers. On peut
  !>        accéder directement au premier et au dernier élément de la
  !>        liste.
  type i_list
     private
     integer :: length = 0
     type(i_node), pointer :: first => NULL()
     type(i_node), pointer :: last => NULL()
  end type i_list

  !> @brief structure de donnée pour une liste de réels de précision
  !>        \a pr. On peut accéder directement au premier et au dernier
  !>        élément de la liste.
  type f_list
     private
     integer :: length = 0
     type(f_node), pointer :: first => NULL()
     type(f_node), pointer :: last => NULL()
  end type f_list
  
  !> @brief noeud d'une liste d'entier, ie un entier et un pointeur
  !>        vers le noeud suivant
  type i_node
     private
     integer :: value = 0
     type(i_node), pointer :: next => NULL()
  end type i_node

  !> @brief noeud d'une liste de réels, ie un réel (avec \a kind=pr)
  !>        et un pointeur vers le noeud suivant
  type f_node
     private
     real(pr) :: value = 0._pr
     type(f_node), pointer :: next => NULL()
  end type f_node

  !> @brief ajoute une valeur ou un tableau de valeurs à la fin de la
  !>        liste
  interface append
     module procedure append_i,       append_f
     module procedure append_i_array, append_f_array
  end interface

  !> @brief libère la mémoire pour la liste
  interface deallocate
     module procedure deallocate_i,       deallocate_f
     module procedure deallocate_i_array, deallocate_f_array
  end interface

  !> @brief affiche la liste
  interface print
     module procedure print_i, print_f
  end interface

  !> @brief copie une liste dans un tableau
  interface to_array
     module procedure to_array_i, to_array_f
  end interface

  !> @brief test si l'élément appartient à la liste
  interface belongs_to
     module procedure belongs_to_i, belongs_to_f
  end interface

  !> @brief renvoie l'indice de l'élément (partant de 1) dans la liste
  !>        ou 0 s'il ne s'y trouve pas
  interface index
     module procedure index_i, index_f
  end interface

  !> @brief renvoie la valeur de l'entrée numéro \a i dans la liste
  interface get_entry
     module procedure get_entry_i, get_entry_f
  end interface get_entry

  !> @brief modifie la valeur de l'entrée \a i dans la liste
  interface set_entry
     module procedure set_entry_i, set_entry_f
  end interface

  !> @brief renvoie la longeur de la liste
  interface get_length
     module procedure get_length_i, get_length_f
  end interface

  !> @brief permute les éléments d'une liste
  interface permute
     module procedure permute_i, permute_f
  end interface

  public :: i_list, f_list
  public :: append, deallocate, print, belongs_to, get_length
  public :: index, get_entry, set_entry, to_array, permute

contains

  !> @brief renvoie la longueur de la liste
  !
  !> @param[in] list la liste
  !> @return le nombre d'éléments de la liste
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  function get_length_i(list)

    type(i_list), intent(in) :: list
    integer :: get_length_i

    get_length_i = list%length

  end function get_length_i

  !> @brief renvoie la longueur de la liste
  !
  !> @param[in] list la liste
  !> @return le nombre d'éléments de la liste
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  function get_length_f(list)

    type(f_list), intent(in) :: list
    integer :: get_length_f

    get_length_f = list%length

  end function get_length_f

  !> @brief libère la mémoire pour une liste
  !
  !> @param[in,out] list la liste
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine deallocate_i(list)

    type(i_list), intent(inout) :: list

    type(i_node), pointer :: node, next_node

    if (associated(list%first)) then 
       node => list%first 
       next_node => node%next
       deallocate(list%first)
       do while (associated(next_node))
          node => next_node
          next_node => node%next
          deallocate(node)
       end do
    end if
   
    nullify(list%first)
    nullify(list%last) 
    list%length = 0

  end subroutine deallocate_i

  !> @brief libère la mémoire pour un tableau de listes
  !
  !> @param[in,out] list(:) le tableau de listes
  !
  !> @details désalloue chaque item de liste(:)
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine deallocate_i_array(list)

    type(i_list), dimension(:), intent(inout) :: list

    integer :: i

    do i = 1, size(list,1)
       call deallocate_i(list(i))
    end do

  end subroutine deallocate_i_array

  !> @brief libère la mémoire pour une liste
  !
  !> @param[in,out] list la liste
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine deallocate_f(list)

    type(f_list), intent(inout) :: list

    type(f_node), pointer :: node, next_node

    if (associated(list%first)) then 
       node => list%first ! Egalité des pointeurs
       next_node => node%next
       deallocate(list%first)
       do while (associated(next_node))
          node => next_node
          next_node => node%next
          deallocate(node)
       end do
    end if
   
    nullify(list%first)
    nullify(list%last) 
    list%length = 0

  end subroutine deallocate_f

  !> @brief libère la mémoire pour un tableau de listes
  !
  !> @param[in,out] list(:) le tableau de listes
  !
  !> @details désalloue chaque item de liste(:)
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine deallocate_f_array(list)

    type(f_list), dimension(:), intent(inout) :: list

    integer :: i

    do i = 1, size(list,1)
       call deallocate_f(list(i))
    end do

  end subroutine deallocate_f_array

  !> @brief ajoute un élément à la fin de la liste
  !
  !> @param[in,out] list la liste
  !> @param[in] val la valeur à ajouter
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine append_i(list,val)

    type(i_list), intent(inout) :: list
    integer, intent(in) :: val

    if (associated(list%last)) then 
       allocate(list%last%next)
       list%last => list%last%next
       list%last%value = val
       nullify(list%last%next)
       list%length = list%length+1
    else
       ! Première allocation
       list%length = 1
       allocate(list%first)
       list%first%value = val
       nullify(list%first%next)
       list%last => list%first
    end if

  end subroutine append_i

  !> @brief ajoute des éléments à la fin de la liste
  !
  !> @param[in,out] list la liste
  !> @param[in] val(:) les valeurs à ajouter
  !
  !> @details les valeurs du tableau val(:) sont ajoutées dans l'ordre
  !>          dans lequel elles apparaissent dans le tableau
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine append_i_array(list, val)

    type(i_list), intent(inout) :: list
    integer, dimension(:), intent(in) :: val

    integer :: i

    do i = 1,size(val,1)
       call append_i(list, val(i))
    end do

  end subroutine append_i_array
  
  !> @brief ajoute un élément à la fin de la liste
  !
  !> @param[in,out] list la liste
  !> @param[in] val la valeur à ajouter
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine append_f(list,val)

    type(f_list), intent(inout) :: list
    real(pr), intent(in) :: val

    if (associated(list%last)) then 
       allocate(list%last%next)
       list%last => list%last%next
       list%last%value = val
       nullify(list%last%next)
       list%length = list%length+1
    else
       list%length = 1
       allocate(list%first)
       list%first%value = val
       nullify(list%first%next)
       list%last => list%first
    end if

  end subroutine append_f

  
  !> @brief ajoute des éléments à la fin de la liste
  !
  !> @param[in,out] list la liste
  !> @param[in] val(:) les valeurs à ajouter
  !
  !> @details les valeurs du tableau val(:) sont ajoutées dans l'ordre
  !>          dans lequel elles apparaissent dans le tableau
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine append_f_array(list, val)

    type(f_list), intent(inout) :: list
    real(pr), dimension(:), intent(in) :: val

    integer :: i

    do i = 1,size(val,1)
       call append_f(list, val(i))
    end do

  end subroutine append_f_array
  

  !> @brief converti la liste en un tableau
  !
  !> @param[in] list la liste
  !> @return val(list%length) les valeurs à ajouter
  !
  !> @details la variable de sortie doit être allouée par le programme
  !>          appelant
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  function to_array_i(list) result(i_tab)
    
    type(i_list), intent(in) :: list
    integer, dimension(list%length) :: i_tab

    type(i_node), pointer :: node => NULL()
    integer :: i

    node => list%first
    i = 0
    do while (associated(node))
       i = i+1
       i_tab(i) = node%value
       node => node%next
    end do

  end function to_array_i

  !> @brief converti la liste en un tableau
  !
  !> @param[in] list la liste
  !> @return val(list%length) les valeurs à ajouter
  !
  !> @details la variable de sortie doit être allouée par le programme
  !>          appelant
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  function to_array_f(list) result(f_tab)
    
    type(f_list), intent(in) :: list
    real(pr), dimension(list%length) :: f_tab

    type(f_node), pointer :: node => NULL()
    integer :: i

    node => list%first
    i = 0
    do while (associated(node))
       i = i+1
       f_tab(i) = node%value
       node => node%next
    end do

  end function to_array_f

  !> @brief affiche la liste
  !
  !> @param[in] list la liste
  !
  !> @todo ne pas utiliser la conversion en tableau, tronquer
  !>       l'affichage si la liste est lngue
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine print_i(list)
    
    type(i_list), intent(in) :: list

    integer, dimension(list%length) :: aux_list

    aux_list(:) = to_array_i(list)
    print "(10I0)", aux_list(:)

!!$    print '(X,A,I5)', 'List of integers of length ',list%length
!!$    node => list%first
!!$    i = 0
!!$    do while (associated(node))
!!$       i = i+1
!!$       print '(X,I5,A,I5)',i,' --> ', node%value
!!$       node => node%next
!!$    end do

  end subroutine print_i

 
  !> @brief affiche la liste
  !
  !> @param[in] list la liste
  !
  !> @todo tronquer l'affichage si la liste est lngue
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine print_f(list)
    
    type(f_list), intent(in) :: list

    type(f_node), pointer :: node => NULL()
    integer :: i

    print '(1X,A,I5)', 'List of real(pr) of length ',list%length
    node => list%first
    i = 0
    do while (associated(node))
       i = i+1
       print '(1X,I5,A,E15.7)',i,' --> ', node%value
       node => node%next
    end do

  end subroutine print_f

  !> @brief test si un nombre est dans la liste
  !
  !> @param[in] list la liste
  !> @param[in] n le nombre à tester
  !> @return .true. si n est dans la liste, .false. sinon
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  function belongs_to_i(list, n) result(flag)

    logical :: flag
    type(i_list), intent(in) :: list
    integer, intent(in) :: n

    type(i_node), pointer :: node => NULL()

    node => list%first
    flag = .false.
    do while (associated(node).and.(.not.flag))
       flag = (node%value==n)
       node => node%next
    end do
    
  end function belongs_to_i

  !> @brief test si un nombre est dans la liste
  !
  !> @param[in] list la liste
  !> @param[in] x le nombre à tester
  !> @return .true. si n est dans la liste, .false. sinon
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  function belongs_to_f(list, x) result(flag)

    logical :: flag
    type(f_list), intent(in) :: list
    real(pr), intent(in) :: x

    type(f_node), pointer :: node => NULL()

    node => list%first
    flag = .false.
    do while (associated(node).and.(.not.flag))
       flag = (node%value==x)
       node => node%next
    end do
    
  end function belongs_to_f

  !> @brief renvoie l'index d'un nombre dans la liste
  !
  !> @param[in] list la liste
  !> @param[in] n le nombre à tester
  !> @return l'indice d'un nombre dans la liste (>=1) ou 0 s'il n'y
  !>         est pas
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  function index_i(list, n)

    type(i_list), intent(in) :: list
    integer, intent(in) :: n
    integer :: index_i

    type(i_node), pointer :: node => NULL()
    integer :: i

    node => list%first
    i = 0
    index_i = 0
    do while (associated(node).and.index_i==0)
       i = i+1
       if (node%value==n) index_i = i
       node => node%next
    end do

  end function index_i

  !> @brief renvoie l'index d'un nombre dans la liste
  !
  !> @param[in] list la liste
  !> @param[in] x le nombre à tester
  !> @return l'indice d'un nombre dans la liste (>=1) ou 0 s'il n'y
  !>         est pas
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  function index_f(list, x)

    type(f_list), intent(in) :: list
    real(pr), intent(in) :: x
    integer :: index_f

    type(f_node), pointer :: node => NULL()
    integer :: i

    node => list%first
    i = 0
    index_f = 0
    do while (associated(node).and.index_f==0)
       i = i+1
       if (node%value==x) index_f = i
       node => node%next
    end do

  end function index_f

  !> @brief renvoie l'item numéro n de la liste
  !
  !> @param[in] list la liste
  !> @param[in] n item à renvoyer
  !> @return l'item liste(n)
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  function get_entry_i(list, n)

    type(i_list), intent(in) :: list
    integer, intent(in) :: n
    integer :: get_entry_i 

    type(i_node), pointer :: node => NULL()
    integer :: i

    if (1<=n.and.n<=list%length) then 
       node => list%first
       do i = 2,n
          node => node%next
       end do
       get_entry_i = node%value
    else
       ! Liste vide ou entrée en dehors de l'intervalle (1..n)
       print "('get_entry_i: indice hors limites')"
       get_entry_i = 0
    end if

  end function get_entry_i

  !> @brief renvoie l'item numéro n de la liste
  !
  !> @param[in] list la liste
  !> @param[in] n item à renvoyer
  !> @return l'item liste(n)
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  function get_entry_f(list, n)

    type(f_list), intent(in) :: list
    integer, intent(in) :: n
    real(pr) :: get_entry_f

    type(f_node), pointer :: node => NULL()
    integer :: i

    if (1<=n.and.n<=list%length) then 
       node => list%first
       do i = 2,n
          node => node%next
       end do
       get_entry_f = node%value
    else
       get_entry_f = 0._pr
       print "('get_entry_f: indice hors limites')"
    end if

  end function get_entry_f

  !> @brief affecte l'item numéro n de la liste avec la valeur val
  !
  !> @param[in] list la liste
  !> @param[in] n item à renvoyer
  !> @param[in] val la valeur à affecter
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine set_entry_i(list, n, val)
    
    type(i_list), intent(in) :: list
    integer, intent(in) :: n, val

    type(i_node), pointer :: node => NULL()
    integer :: i

    if (1<=n.and.n<=list%length) then 
       ! Si la liste est non vide et si n est dans l'intervalle (1..n)
       node => list%first
       do i = 2,n
          node => node%next
       end do
       node%value = val
    else
       print "('set_entry_i: indice hors limites')"
    end if

  end subroutine set_entry_i


  !> @brief affecte l'item numéro n de la liste avec la valeur val
  !
  !> @param[in] list la liste
  !> @param[in] n item à renvoyer
  !> @param[in] val la valeur à affecter
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine set_entry_f(list, n, val)
    
    type(f_list), intent(in) :: list
    integer, intent(in) :: n
    real(pr), intent(in) :: val

    type(f_node), pointer :: node => NULL()
    integer :: i

    if (1<=n.and.n<=list%length) then 
       ! Si la liste est non vide et si n est dans l'intervalle (1..n)
       node => list%first
       do i = 2,n
          node => node%next
       end do
       node%value = val
    else
       print "('set_entry_f: indice hors limites')"
    end if

  end subroutine set_entry_f

  !> @brief permute les éléments d'une liste
  !
  !> @param[in,out] list la liste
  !> @param[in] perm(:) la permutation
  !
  !> @details effectue la permutation suivante: list(i) := list(perm(i))
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine permute_i(list,perm)

    type(i_list), intent(inout) :: list
    integer, dimension(:), intent(in) :: perm

    integer :: i
    type(i_list) :: new_list

    new_list%length = list%length
    do i = 1,list%length
       call append_i( new_list, get_entry_i(list,perm(i)) )
    end do

    call deallocate_i(list)

    list%length = new_list%length
    list%first => new_list%first
    list%last => new_list%last
    
  end subroutine permute_i


  !> @brief permute les éléments d'une liste
  !
  !> @param[in,out] list la liste
  !> @param[in] perm(:) la permutation
  !
  !> @details effectue la permutation suivante: list(i) := list(perm(i))
  !
  !> @author Yves Coudière
  !> - 26/4/2013 : réécriture depuis la bibliothèque spmlib
  subroutine permute_f(list,perm)

    type(f_list), intent(inout) :: list
    integer, dimension(:), intent(in) :: perm

    integer :: i
    type(f_list) :: new_list

    new_list%length = list%length
    do i = 1,list%length
       call append_f( new_list, get_entry_f(list,perm(i)) )
    end do

    call deallocate_f(list)

    list%length = new_list%length
    list%first => new_list%first
    list%last => new_list%last
    
  end subroutine permute_f


end module list
