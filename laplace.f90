program laplace_prof

  !> @brief Equation alpha*u - Laplacien(u) = f
  !> @author Yves Coudière
  !> @date Time-stamp: <2017-11-20 17:32:32 yves>
  !
  !> @detail Par elts finis P1-Lagrange, avec le solveur direct LU creux venant
  !> de UMFPACK. Le problème est posé dans un domaine lu dans un fichier de
  !> maillage, on applique des conditions de Neumann homogènes au bord.
  !
  !>
  !> Trois fonctions f (second membre) sont disponibles:
  !> 1. f(x) = 1
  !> 2. f(x) = 1 dans un disque de rayon r_0 et de centre x0, et 0 ailleurs
  !> 3. f(x) = (alpha - dim*(k*pi)**2)*cos(k*pi*x_1)...*cos(k*pi*x_dim)
  !
  !> Pour le cas 3., la solution exacte est connue (si le domaine est le carré
  !> (0,1)^dim), c'est le produit des cos.
  use prec
  use fe_mesh_quentin
  use mat_list
  use mat_csr
  use mat_csc
  use p1
  implicit none

  ! Données du problème, qui sont lues dans le fichier d'entrées
  character(len=256) :: mesh_file
  real(pr) :: alpha, k, r_0,xymin,xymax
  real(pr), dimension(:), allocatable :: x_0
  integer :: choix_rhs,dir_trac

  ! Variables globales
  type(mesh) :: m ! Maillage
  type(m_list) :: K_l, M_l ! Matrices au format de listes
  type(csr) :: M_csr ! Matrice au format CSR
  type(csc) :: K_csc ! Matrice au format CSC
  real(pr), dimension(:), allocatable :: F,Y,X ! Divers vecteurs

  real :: cpu_t1, cpu_t2
  integer :: n_nodes, dim, i


  ! Lecture des données dans un fichier
  call lecture_donnees('data.txt')

  ! Lecture du maillage dans le fichier mesh_file
  call cpu_time(cpu_t1)
  call read_mesh(m, mesh_file)
  call cpu_time(cpu_t2)
  print "(4X,'cpu time for reading the mesh: ',F6.4,' s')", cpu_t2-cpu_t1
  n_nodes = get_n_nodes(m)
  dim = get_dim(m)

  ! Assemblage des matrices sur ce maillage
  call cpu_time(cpu_t1)
  call assemble_mass_matrix(m, M_l)
  call cpu_time(cpu_t2)
  print "(4X,'cpu time for assembling the mass matrix: ',F6.4,' s')", cpu_t2-cpu_t1
  call assemble_stiffness_matrix(m,K_l)
  call cpu_time(cpu_t1)
  print "(4X,'cpu time for assembling the stiffness matrix: ',F6.4,' s')", cpu_t1-cpu_t2
  call alpha_a_plus_beta_b( 1._pr, K_l, alpha, M_l) ! K := alpha*M + K

  ! Changement stockage creux et libération mémoire
  M_csr = M_l ! Matrice de masse
  K_csc = K_l ! Matrice du système linéaire (alpha u - Delta(u))
  call deallocate(M_l)
  call deallocate(K_l)
  call cpu_time(cpu_t2)
  print "(4X,'cpu time for building the final matrices in CSR and CSS storage: ',F6.4,' s')", cpu_t2-cpu_t1

  ! Assemblage du second membre
  allocate(Y(n_nodes),F(n_nodes))
  call compute_rhs_function(F)
  ! Second membre du système linéaire
  Y = M_csr*F

  call deallocate(M_csr) ! La matrice CSR ne sert qu'au produit matrice-vecteur
  call cpu_time(cpu_t1)
  print "(4X,'cpu time for assembling the RHS: ',F6.4,' s')", cpu_t1-cpu_t2
  ! Sauvegarde fichier pour le second membre (pour affichage)
  i = index(mesh_file,'.',back=.true.)
  call write_mesh(m,mesh_file(1:i)//'rhs.vtk',F)
  call cpu_time(cpu_t2)
  print "(4X,'cpu time for writing the RHS function: ',F6.4,' s')", cpu_t2-cpu_t1

  ! Factorisation de la matrice M (csc)
  call sort_cols(K_csc) ! UMFPACK requires to have the column sorted by
  call factorize(K_csc) ! increasing row index.
  call cpu_time(cpu_t1)
  print "(4X,'cpu time for the factorization: ',F6.4,' s')", cpu_t1-cpu_t2

  ! Résolution du système linéaire
  allocate(X(n_nodes))
  call solve(K_csc,X,Y)
  call cpu_time(cpu_t2)
  print "(4X,'cpu time for the resolution: ',F6.4,' s')", cpu_t2-cpu_t1

  ! Write the solution vector X to a vtk file
  call write_mesh(m,mesh_file(1:i)//'solution.vtk',X)

  ! ... to be completed, adapted...

  ! On libère la mémoire
  call deallocate(K_csc)
  call deallocate(m)
  deallocate(F,Y,X)
  if (allocated(x_0)) deallocate(x_0)

contains

  ! Assemblage du second membre (partie supérieure seulement)
  subroutine compute_rhs_function(F)
    real(pr), dimension(:), intent(inout) :: F

    real(pr), dimension(get_dim(m)) :: x
    integer :: i,j

    select case(choix_rhs)
    case(1) ! f(x) = 1.
       do i = 1, get_n_nodes(m)
          F(i) = 1._pr
       end do
    case(2) ! f(x) = 1. si |x-x_0| < r_0, 0 sinon
       do i = 1, get_n_nodes(m)
          x(:) = get_x(m,i)
          if (sqrt(dot_product(x-x_0,x-x_0)) < r_0) then
             F(i) = 1._pr
          else
             F(i) = 0._pr
          end if
       end do
    case(3) ! f(x) = (alpha - dim*(k*pi)**2)*cos(k*pi*x_1)...*cos(k*pi*x_dim)
       do i = 1, get_n_nodes(m)
          x(:) = get_x(m,i)
          F(i) = (alpha - dim*(k*pi)**2)
          do j = 1,dim
             F(i) = F(i) * cos(k*pi*x(j))
          end do
       end do
    case default
       print "('**** ATTENTION : cas non implémenté ****')"
    end select

  end subroutine compute_rhs_function

  ! Données qui sont lues
  ! mesh_file
  ! alpha
  ! choix_rhs
  ! x0, r (cas 2) ou k (cas 3)
  subroutine lecture_donnees(f_name)
    character(len=*), intent(in) :: f_name

    open(unit=10,file=trim(f_name), status='old', action='read')
    read (10,*) mesh_file
    read (10,*) alpha
    read (10,*) choix_rhs
    select case(choix_rhs)
    case(1)
       read(10,*)
    case(2) ! on a besoin de x_0 et r_0
       allocate(x_0(2)) ! seulement en dimension 2 pour l'instant
       read(10,*) x_0(:), r_0
    case(3) ! on a besoin de k
       read(10,*) k
    end select
    read(10,*) dir_trac
    read(10,*) xymin
    read(10,*) xymax

    close(10)

  end subroutine lecture_donnees

end program laplace_prof
