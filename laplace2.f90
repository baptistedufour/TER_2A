program laplace

  !> @author Yves Coudi�re
  !> @date Time-stamp: <2017-11-20 17:32:32 yves>
  !
  !> @detail Par elts finis P1-Lagrange, avec le solveur direct LU creux venant
  !> de UMFPACK. Le probl�me est pos� dans un domaine lu dans un fichier de

  use prec
  use fe_mesh_quentin
  use mat_list
  use mat_csr
  use mat_csc
  use p1
  implicit none

  ! Donn�es du probl�me, qui sont lues dans le fichier d'entr�es
  character(len=256) :: mesh_file
  real(pr) :: alpha, k, r_0, xymin,xymax, valeur_force
  real(pr), dimension(:), allocatable :: x_0
  integer :: choix_rhs, dir_trac

  ! Variables globales
  type(mesh) :: m, m2 ! Maillage
  type(m_list) :: K_l, M_l ! Matrices au format de listes
  type(csr) :: M_csr ! Matrice au format CSR
  type(csc) :: K_csc ! Matrice au format CSC
  real(pr), dimension(:), allocatable :: F,X ! Divers vecteurs

  real :: cpu_t1, cpu_t2
  integer :: n_nodes, dim, i

  ! Lecture des donn�es dans un fichier
  call lecture_donnees('data.txt')

  ! Lecture du maillage dans le fichier mesh_file
  call cpu_time(cpu_t1)
  call read_mesh(m, mesh_file)
  call chg_lbl(m,xymin,xymax,dir_trac)
  call cpu_time(cpu_t2)
  print "(4X,'cpu time for reading the mesh: ',F6.4,' s')", cpu_t2-cpu_t1
  n_nodes = get_n_nodes(m)
  dim = get_dim(m)

  ! Assemblage des matrices sur ce maillage

  call cpu_time(cpu_t2)
  call assemble_stiffness_matrix(m,K_l)
  call cpu_time(cpu_t1)
  print "(4X,'cpu time for assembling the stiffness matrix: ',F6.4,' s')", cpu_t1-cpu_t2

  ! Changement stockage creux et lib�ration m�moire
  M_csr = M_l ! Matrice de masse
  K_csc = K_l ! Matrice du syst�me lin�aire (alpha u - Delta(u))

  call deallocate(M_l)
  call deallocate(K_l)
  call cpu_time(cpu_t2)
  print "(4X,'cpu time for building the final matrices in CSR and CSS storage: ',F6.4,' s')", cpu_t2-cpu_t1

  ! Assemblage du second membre
  allocate(F(2*n_nodes))
  call compute_rhs_function(F, dir_trac, valeur_force)
  ! Second membre du syst�me lin�aire
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

  ! R�solution du syst�me lin�aire
  allocate(X(2*n_nodes))
  X=0_pr
  call solve(K_csc,X,F)
  call cpu_time(cpu_t2)
  print "(4X,'cpu time for the resolution: ',F6.4,' s')", cpu_t2-cpu_t1

  ! Write the solution vector X to a vtk file
  call write_mesh(m,mesh_file(1:i)//'solution.vtk',X)

  !Afficher le maillage déformé
  m2=m
  call strain(m2,X)
  call write_mesh(m2,mesh_file(1:i)//'deplacement.vtk')

  ! ... to be completed, adapted...

  ! On lib�re la m�moire
  call deallocate(K_csc)
  call deallocate(m)

  print*, "X", X
  deallocate(F,X)
  if (allocated(x_0)) deallocate(x_0)

contains

  ! Assemblage du second membre (partie sup�rieure et inferieur )
  subroutine compute_rhs_function(F,dir_trac, valeur_force)
    real(pr), dimension(:), intent(inout) :: F
    real(pr), intent(in) :: valeur_force
    real(pr) :: h
    integer :: i,dim,n_nodes,k
    integer, intent(in) :: dir_trac
    integer:: dir_trac_perp !direction perpendiculaire à la direction de traciotn pour les arrêtes
    integer, dimension(:), allocatable :: nodes
    real(pr), dimension(:,:), allocatable :: x


    dir_trac_perp= dir_trac+(-1)**(dir_trac+1)
    dim = get_dim(m)
    allocate( nodes(dim+1), x(dim,dim+1) )
    n_nodes = get_n_nodes(m)

    F=0_pr

    !Cas où la force est selon x
    if (dir_trac==1) then

      !parcourt tous les elements
      do i = 1, get_n_elts(m)
         !prend les noeux du triangle k
         nodes(:) = get_nodes(m,i)

         !prend les cordonn�e des point du triangle i
         do k = 1,dim+1
            x(:,k) = get_x(m,nodes(k))
         end do

         !test si une arrete du triangle k est sur le bord
         if (((get_node_code(m,nodes(1))==3) .AND. (get_node_code(m,nodes(2))==3)) .OR. &
         ((get_node_code(m,nodes(2))==3) .AND. (get_node_code(m,nodes(3))==3)) &
          .OR. ((get_node_code(m,nodes(1))==3) .AND. (get_node_code(m,nodes(3))==3))) then

          !si le premier point n'est pas sur le bord
          if (get_node_code(m,nodes(1))==0)then
             h=abs(x(dir_trac_perp,2)-x(dir_trac_perp,3))
             print*, x(2,2), x(2,3)
             F(nodes(2))= h*valeur_force
             !F(nodes(2)+n_nodes) = h*valeur_force
             F(nodes(3))= h*valeur_force
             !F(nodes(3)+n_nodes) = h*valeur_force
          end if
          ! le deuxieme point n'est pas sur le bord
          if (get_node_code(m,nodes(2))==0)then
             h=abs(x(dir_trac_perp,1)-x(dir_trac_perp,3))
             F(nodes(3))= h*valeur_force
             !F(nodes(3)+n_nodes) = h*valeur_force
             F(nodes(1))= h*valeur_force
             !F(nodes(1)+n_nodes) = h*valeur_force
          end if
          !le troisieme point n'est pas sur le bord
          if (get_node_code(m,nodes(3))==0)then

             h=abs(x(dir_trac_perp,2)-x(dir_trac_perp,1))
             F(nodes(1))= h*valeur_force
             !F(nodes(1)+n_nodes) = h*valeur_force
             F(nodes(2))= h*valeur_force
             !F(nodes(2)+n_nodes) = h*valeur_force
          end if

       end if

    end do


    !Cas ou la force est selon y
    elseif (dir_trac==2) then

    !parcourt tous les elements
    do i = 1, get_n_elts(m)
       !prend les noeux du triangle k
       nodes(:) = get_nodes(m,i)

       !prend les cordonn�e des point du triangle i
       do k = 1,dim+1
          x(:,k) = get_x(m,nodes(k))
       end do

       !test si une arrete du triangle k est sur le bord
       if (((get_node_code(m,nodes(1))==3) .AND. (get_node_code(m,nodes(2))==3)) .OR. &
       ((get_node_code(m,nodes(2))==3) .AND. (get_node_code(m,nodes(3))==3)) &
        .OR. ((get_node_code(m,nodes(1))==3) .AND. (get_node_code(m,nodes(3))==3))) then

        !si le premier point n'est pas sur le bord
        if (get_node_code(m,nodes(1))==0)then
           h=abs(x(dir_trac_perp,2)-x(dir_trac_perp,3))
           print*, x(2,2), x(2,3)
           !F(nodes(2))= h*valeur_force
           F(nodes(2)+n_nodes) = h*valeur_force
           print*, h,  valeur_force, nodes(2)+n_nodes
           !F(nodes(3))= h*valeur_force
           F(nodes(3)+n_nodes) = h*valeur_force
        end if
        ! le deuxieme point n'est pas sur le bord
        if (get_node_code(m,nodes(2))==0)then
           h=abs(x(dir_trac_perp,1)-x(dir_trac_perp,3))
           !F(nodes(3))= h*valeur_force
           F(nodes(3)+n_nodes) = h*valeur_force
           !F(nodes(1))= h*valeur_force
           F(nodes(1)+n_nodes) = h*valeur_force
        end if
        !le troisieme point n'est pas sur le bord
        if (get_node_code(m,nodes(3))==0)then

           h=abs(x(dir_trac_perp,2)-x(dir_trac_perp,1))
           !F(nodes(1))= h*valeur_force
           F(nodes(1)+n_nodes) = h*valeur_force
           !F(nodes(2))= h*valeur_force
           F(nodes(2)+n_nodes) = h*valeur_force
        end if

     end if

  end do

end if

  end subroutine compute_rhs_function

  ! Donn�es qui sont lues
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
    read(10,*) valeur_force

    close(10)

  end subroutine lecture_donnees

end program laplace
