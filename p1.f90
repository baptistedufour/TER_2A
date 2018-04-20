!> @brief Assemblage des matrices de u --> u (masse) et u --> - Laplacien(u)
!> avec conditions aux limites de Neuman homogènes, par éléments finis P1
!> Lagrange sur maillage non structuré.
!
!> @date Time-stamp: <2017-10-19 14:30:55 yves>
!> @todo Reprendre toute la doc pour doxygen
!

!> @details L'assemblage est réalisé grâce aux matrices mat_list. L'utilisateur
!> est ensuite libre de convertir ce format en celui qu'il veut utiliser en
!> pratique pour les calculs.
module p1

  use prec
  use fe_mesh_quentin
  use mat_list
  use fonction
  use mat_csc
  use mat_csr
  use matrices
  implicit none

  private

  public :: assemble_mass_matrix, assemble_stiffness_matrix

contains

  !> @brief Assemblage de la matrice de masse
  !
  !> @details Assemblage de la matrice de masse EF P1 avec integration numerique
  !> (condensation). Cette matrice est donc diagonale : c'est la matrice de
  !> a(u,v) = int_O u v dx. Avec une integration numerique, on trouve que
  !>
  !> \f$ a_{ij} = int_{\Omega} u_i v_j dx = \sum_{T} m(T)/(d+1) \sum_{k \in T}
  !>     u_i(x_k) u_j(x_k) = \sum_{T} m(T)/(d+1) \delta_{ik} \delta_{jk} \f$
  !>
  !> ou \f$d\f$ est la dimension et \f$\delta_{ik}\f$ est le symbole de
  !> Kroneker.  On voit donc que \f$ a_{ij} = 0\f$ si \f$i \ne j\f$, et
  !> \f$a_{ii} = \sum_{T} m(T)/(d+1) \sum_{k \in T} \delta_{ik}\f$.
  !>
  !> La condition aux limites de Dirichlet n'entre pas en ligne de compte ici
  !> (on traite un pb de Neumann homogène).
  !>
  !>
  !> @param[in] m le maillage
  !> @param[out] mass matrice de masse assemblée
  subroutine assemble_mass_matrix (m, mass)
    type(mesh), intent(in) :: m
    type(m_list), intent(out) :: mass

    real(pr), dimension(:), allocatable :: mass_elem
    integer :: i, k, dim, n_nodes, n_elts
    real(pr) :: measure
    integer, dimension(:), allocatable :: nodes

    print "('Assemblage de la matrice de masse...')"

    ! Allocation mémoire initiale
    dim = get_dim(m)
    n_nodes = get_n_nodes(m)
    n_elts = get_n_elts(m)
    call allocate(mass,n_nodes,n_nodes)
    allocate(mass_elem(dim+1), nodes(dim+1))

    ! Matrice de masse élémentaire (diagonale)
    mass_elem(:) = 1._pr/real(dim+1, pr)

    ! Assemblage par une boucle sur les elements
    do i = 1,n_elts
       measure = get_measure(m,i)
       nodes(:) = get_nodes(m,i)
       ! mass{JJ} <- mass{JJ} + meas*elemmass{jj}
       do k = 1,dim+1
          call mat_add_value(mass, nodes(k), nodes(k), &
               measure*mass_elem(k))
       end do
    end do

    deallocate(mass_elem, nodes)

  end subroutine assemble_mass_matrix

  !> @brief Assemblage de la matrice de raideur pour un coefficient de diffusion
  !> D(x):=1
  !
  !> @details Assemblage de la matrice de raideur P1, c'est a dire la matrice de
  !> discretisation de \f$ a(u,v) = \int_{\Omega} D(x)\nabla u \cdot \nabla v dx
  !> \f$.
  !
  !> En général, \f$ D(x) \f$ est un tenseur de diffusion et \f$ \Omega \f$ est
  !> un maillage de simplexes 1,2 ou 3D. Les coefficients de la matrice sont
  !> donc \f$ - K_{ij} = \int_\Omega \nabla u_i \cdot D(x) \nabla u_j \f$.
  !
  !> On se place dans le cas de conditions aux limites de Neumann
  !> homogènes, et \f$D(x):=1\f$.
  !
  !>
  !> @param[in] m maillage
  !> @param[out] stiffness matrice de raideur
  subroutine assemble_stiffness_matrix(m, stiffness)
    type(mesh), intent(in)    :: m
    type(m_list), intent(out) :: stiffness

    integer :: dim, n_nodes, n_elts, i, k, j
    integer, dimension(:), allocatable :: nodes
    real(pr), dimension(:,:), allocatable :: gradient, diff, x
    real(pr) :: measure
    real(pr) :: K_11_kj, K_22_kj, K_12_kj, mu, lambda

    print "('Assemblage de la matrice de raideur...')"

    ! Allocation mémoire
    dim = get_dim(m)
    n_nodes = get_n_nodes(m)
    n_elts = get_n_elts(m)
    call allocate(stiffness,2*n_nodes,2*n_nodes)
    allocate( nodes(dim+1), gradient(dim,dim+1), diff(dim,dim), x(dim,dim+1) )

    ! Debut de la boucle d'assemblage
    do i = 1,n_elts
       ! On recupere les sommets
       nodes(:) = get_nodes(m,i)
       ! Et les coordonnées des noeuds
       do k = 1,dim+1
          x(:,k) = get_x(m,nodes(k))
       end do
       ! On calcule les gradients des fonctions de forme dans
       ! l'element de sommets x. Le résultat est dans gradient :
       ! gradient(:,k) est le gradient de la fonction de base du noeud
       ! k
       call compute_shape_functions(x, gradient)
       ! gradient(i,k) := d_i(phi_k)
       ! Taille de l'élément
       measure = get_measure(m,i)
       ! Valeurs moyennes de lambda et mu

       lambda = (1._pr/3._pr)*(get_lambda(x(1,1), x(2,1))+get_lambda(x(1,2), x(2,2))+get_lambda(x(1,3), x(2,3)))
       mu = (1._pr/3._pr)*(get_mu(x(1,1), x(2,1))+get_mu(x(1,2), x(2,2))+get_mu(x(1,3), x(2,3)))
!!$       ! On a 6 matrices élémentaires, 3 en lambda et 3 en mu:
       do k = 1,dim+1
          do j = 1,dim+1

!!$             ! K^{11_lambda}_KJ= lambda*d_1phi_J*d_1phi_I*|T|
            print*, "triangle = ", i, "gradient =", gradient (1,j), gradient(1,k), gradient(2,j), gradient(2,k)
             K_11_kj = gradient(1,j)*gradient(1,k)*measure
             K_12_kj = gradient(1,k)*gradient(2,j)*measure
             K_22_kj = gradient(2,j)*gradient(2,k)*measure


!!$             ! A_KJ += K_kj

!!$          Membre en lambda

             !! a regarder indice dans stifness ou dans mat_add_value
             call mat_add_value(stiffness, nodes(k),nodes(j),lambda*K_11_kj)
             !print*, nodes(k),nodes(j),lambda* K_11_kj, "lamnda*K_11_kj"
             !read*,
             call mat_add_value(stiffness, nodes(k)+n_nodes,nodes(j)+n_nodes,lambda*K_22_kj) !basdroite
             !print*, nodes(k)+n_nodes,nodes(j)+n_nodes,lambda* K_22_kj,"lambda*K_11_kj"
             !read*,

             call mat_add_value(stiffness, nodes(k),nodes(j)+n_nodes,lambda*K_12_kj) !haut droite
             !print*, nodes(k),nodes(j)+n_nodes,lambda* K_12_kj,"lambda*K_12_kj"
             !read*,
             call mat_add_value(stiffness, nodes(j)+n_nodes, nodes(k),lambda*K_12_kj)!bas gauche
             !print*, nodes(j)+n_nodes, nodes(k),lambda* K_12_kj,"lambda*K_21_kj"
             !read*,

!!$          Membre en mu

             call mat_add_value(stiffness, nodes(k),nodes(j),2*mu*K_11_kj)
             !print*, nodes(k),nodes(j),2*mu* K_11_kj, "mu*K_11_kj"
             !read*,
             call mat_add_value(stiffness, nodes(k),nodes(j),mu*K_22_kj)
             !print*, nodes(k),nodes(j),mu*K_22_kj, "mu*K_22_kj"
             !read*,

             call mat_add_value(stiffness, nodes(j),nodes(k)+n_nodes,mu*K_12_kj) !haut droite
             !print*,nodes(j),nodes(k)+n_nodes,mu*K_12_kj
             !read*,
             call mat_add_value(stiffness, nodes(k)+n_nodes,nodes(j),mu*K_12_kj) !bas gauche
             !print*, nodes(k)+n_nodes,nodes(j),mu*K_12_kj
             !read*,
             call mat_add_value(stiffness, nodes(k)+n_nodes,nodes(j)+n_nodes,2*mu*K_22_kj)
             !print*, nodes(k)+n_nodes,nodes(j)+n_nodes,2*mu*K_22_kj
             !read*,
             call mat_add_value(stiffness, nodes(k)+n_nodes,nodes(j)+n_nodes,mu*K_11_kj)
             !print*, nodes(k)+n_nodes,nodes(j)+n_nodes,mu*K_11_kj
             !read*,

          end do

       end do

    end do

    call condition_fixe(stiffness, m)

    !call print(stiffness)

    deallocate(nodes, gradient, diff,  x)

  end subroutine assemble_stiffness_matrix

  subroutine  condition_fixe(stiffness,m)

    type(mesh), intent(in)    :: m
    type(m_list), intent(inout) :: stiffness
    integer :: dim, n_nodes, i, k, j, n_elts
    integer, dimension(:), allocatable :: nodes

    dim = get_dim(m)
    n_nodes = get_n_nodes(m)
    n_elts = get_n_elts(m)
    allocate( nodes(dim+1))
    print*, "nombre d'elements = ", n_elts
    do i = 1,n_elts
        ! On recupere les sommets
       nodes(:) = get_nodes(m,i)
       do k = 1,dim+1
         print*, "dimnesion +1", dim+1
          !print*, get_node_code(m,nodes(k)), nodes(k)
          if (get_node_code(m,nodes(k))==1) then !(m,node(k))
             do j=1,2*n_nodes
                call mat_add_value(stiffness, nodes(k),j,0._pr,mat_set)
                call mat_add_value(stiffness, nodes(k)+n_nodes,j,0._pr,mat_set)
                call mat_add_value(stiffness,j,nodes(k),0._pr,mat_set)
                call mat_add_value(stiffness,j,nodes(k)+n_nodes,0._pr,mat_set)
             end do
             call mat_add_value(stiffness, nodes(k)+n_nodes, nodes(k)+n_nodes,1._pr, mat_set)
             print*,  nodes(k)+n_nodes, nodes(k)+n_nodes,1
             call mat_add_value(stiffness, nodes(k), nodes(k),1._pr, mat_set)
             print*, nodes(k), nodes(k),1.
          end if
       end do
    end do


    call print(stiffness)

    deallocate(nodes)

  end subroutine  condition_fixe



  !> @brief Computes the P1 shape functions in 1 element.
  !
  !> @details The array gradient is of size dim x (dim+1), and gradient(:,j) is
  !> the gradient of the shape function for node j in the element.
  !
  !> @param[in] x coordinates of the vertices of the simplex (the geometrical
  !> element)
  !> @param[out] gradient(:,:) the gradients of the shape functions, tableau de
  !>                taille dx(d+1) ou d=1,2,3 est la dimension d'espace. On a
  !>                gradient(:,k) est le gradient de la fonction de base P1
  !>                associe au noeud k de l'element i.
  subroutine compute_shape_functions(x, gradient)
    real(pr), dimension(:,:), intent(in) :: x
    real(pr), dimension(:,:), intent(inout) :: gradient

    integer :: dim, j
    real(pr), dimension(:,:), allocatable :: a

    dim = size(x,1)

    if (size(gradient,1)/=dim.or.size(gradient,2)/=dim+1 &
         .or.size(x,2)/=dim+1) then
       print "('p1: compute_shape_functions -- problèmes de dimension')"
       stop
    end if

    ! Allocations
    allocate( a(dim,dim) )

    ! Si X_i (i=1,dim+1) sont les noeuds de K, les gradients des
    ! fonctions de forme P1 verifient grad(Phi).dX =
    ! grad(Phi_ref).A^{-1}d(X_ref), ou A est la matrice de passage A =
    ! [X_{1+i}-X_1]. Nous allons donc résoudre A^T grad(Phi) =
    ! grad(Phi_ref)

    ! Gradients dans l'element de reference, en colonne
    ! gradient_1 = (-1,-1,...,-1)^T
    ! gradient_{1+i} = (delta_{ij})^T_{j=1..dim}
    gradient(:,:) = 0._pr
    do j = 1,dim
       gradient(j,1) = -1._pr
       gradient(j,j+1) = 1._pr
    end do

    ! Transpose de la matrice de passage Colonne j de A = x_{1+j}-x_1, j=1,dim
    ! --> Ligne j de A^T = x_{1+j}-x_1, j=1,dim
    do j = 1,dim
       a(j,:) = x(:,1+j) - x(:,1)
    end do

    ! A gradient(:,j) = gradient_ref(:,j) ou A est la matrice ci-dessus
    ! et les gradients sont en colonnes.
    call gauss_method(a,gradient)

    deallocate(a)

  end subroutine compute_shape_functions


  !> @brief Procédure interne pour résoudre des systèmes linéaires
  !> plein de taille 1, 2, 3.
  !
  !> @details Ne sert que pour l'assemblage des matrices de
  !> raideur/rigidité, par la méthode de passage à l'élément de
  !> référence.
  !>
  !> C'est une implémentation simple du pivot de Gauss avec pivot partiel.
  !> @param[in,out] a(n,n) matrice. Elle est modifiée par l'algorithme
  !> de pivot.
  !> @param[in,out] b(n,p) vecteur de p seconds membres en entrée, et
  !> contenant les solutions des p systèmes linéaires en sortie.
  subroutine gauss_method(a,b)
    real(pr), dimension(:,:), intent(inout) :: a
    real(pr), dimension(:,:), intent(inout) :: b

    integer :: n,p,i,j,k, i0
    real(pr) :: pivot
    real(pr), dimension(:), allocatable :: c

    n = size(a,1) ! taille du système linéaire
    p = size(b,2) ! nombre de second membres
    allocate( c(max(n,p)) ) ! vecteur pour échanges de ligne au pivotage

    ! Triangulation de la matrice, et en même temps, opérations sur le
    ! second membre.
    do i = 1,n-1
       ! Recherche du pivot le plus grand
       c(:) = 0._pr
       c(i:n) = abs(a(i:n,i)) ! tous les pivots possibles
       i0 = maxloc(c,1) ! indice de la ligne qui a le plus grand coefficient
       ! Echange des lignes i et i0 dans a
       c(i:n) = a(i,i:n)
       a(i,i:n) = a(i0,i:n)
       a(i0,i:n) = c(i:n)
       ! Echange des lignes i et i0 dans b
       c(1:p) = b(i,1:p)
       b(i,1:p) = b(i0,1:p)
       b(i0,1:p) = c(1:p)
       ! Et on peut continuer
       do j = i+1,n
          pivot = a(j,i)/a(i,i)
          do k = i+1,n   ! Ligne_j := Ligne_j - pivot * Ligne_i
             a(j,k) = a(j,k) - pivot * a(i,k)
          end do
          do k = 1,p     ! b_j := b_j - pivot * b_i
             b(j,k) = b(j,k) - pivot * b(i,k)
          end do
       end do
    end do

    ! Il ne reste qu'à faire la remontée
    do i = n,1,-1
       do k = 1,p
          do j = i+1,n ! b_i = b_i - sum_{i+1}^n a_{ij} b_j
             b(i,k) = b(i,k) - a(i,j)*b(j,k)
          end do
          b(i,k) = b(i,k) / a(i,i)
       end do
    end do

    deallocate(c)

  end subroutine gauss_method

end module p1
