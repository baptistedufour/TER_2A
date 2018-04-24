!> @brief I/O des maillages pour quelques formats classiques, et type
!> dérivé maillage éléments finis
!
!> @date Time-stamp: <2017-10-22 13:18:07 yves>
!> @todo Reprendre toute la doc pour doxygen
!
!> @details Implements a user-defined mesh type for finite-element
!> computations in 1,2,3D.
!>

!> The internal representation of meshes is based on 2 arrays, of nodes
!> coordinates and the data of the segment, triangles or tetraedra. Additional
!> informations like subdomains and boundaries are coded by integer labels. The
!> module supports reading and writing 2D and 3D meshes from files. File formats
!> supported on reading are currently:
!>   - .mesh format from the Gamma project (Inria) in 2D or 3D,
!>   - .{node,ele} format output from the meshing software triangle, in 2D
!>
!> File format supported on writing are:
!>   - .mesh format,
!>   - .vtk legacy VTK file format.
module fe_mesh_quentin

  use prec
  implicit none

  private

#include "vtkCellType.inc"

  !> Un type de donnee assez simple pour manipuler efficacement les
  !> maillages elements finis : conforme et simpliciaux en 1, 2, 3D.
  !>
  !> Il contient en particulier
  !>   - la dimension
  !>   - la liste des noeuds
  !>   - la liste des elements
  !>   - un code entier par element pour reperer a quel sous domaine
  !>     il appartient
  !>   - la liste des noeuds auxquels on veut imposer ue condition de
  !>     Dirichlet. Les autres sont soumis a une condition de Neumann
  type mesh
     private
     integer :: dim    !> Dimension
     integer :: n_nodes !> Total number of Nodes
     integer :: n_elts  !> Total number of Elements
     integer :: n_faces !> Total number of Faces (usually on ly in 3D meshes)
     !> node(1:dim,1:nnodes) are the coordinates of the nodes.
     real(pr), dimension(:,:), allocatable :: node
     !> node_code(1:dim,1:nnodes) are the codes of the nodes.
     integer, dimension(:), allocatable :: node_code
     !> elt(1:dim+1,1:n_elts) are the indices of the vertices of the elements.
     integer, dimension(:,:), allocatable :: elt
     !> elt_code(1:n_elts) is the code of the elements.
     integer, dimension(:),   allocatable :: elt_code
     !> face(1:dim,n_faces) in 3D only, because faces may be
     !> available in the mesh file
     integer, dimension(:,:), allocatable :: face
     !> face_code are the labels of the faces
     integer, dimension(:), allocatable :: face_code
  end type mesh

  !> @brief Free memory
  interface deallocate
     module procedure deallocate_mesh
  end interface deallocate

  !> @frief Prints some informations on the mesh
  interface my_print
     module procedure print_mesh
  end interface my_print

  public :: mesh
  public :: deallocate, read_mesh, write_mesh, my_print
  public :: get_n_subdomains, get_dim, get_n_elts, get_n_nodes, get_distance
  public :: get_measure, get_nodes, get_elt_code, get_node_code, get_x
  public :: chg_lbl, strain

contains

  !> @brief Get mesh dimension in space
  pure integer function get_dim(m)
    type(mesh), intent(in) :: m

    get_dim = m%dim
  end function get_dim

  !> @brief Get number of nodes
  pure integer function get_n_nodes(m)
    type(mesh), intent(in) :: m

    get_n_nodes = m%n_nodes
  end function get_n_nodes

  !> @brief Get number of elements
  pure integer function get_n_elts(m)
    type(mesh), intent(in) :: m

    get_n_elts = m%n_elts
  end function get_n_elts

  !> @brief Get the number of subdomains, assuming that they are
  !> unmbered from 0 to n-1
  integer function get_n_subdomains(m)
    type(mesh), intent(in) :: m
    get_n_subdomains = maxval(m%elt_code)+1
    if (minval(m%elt_code)<0) then
       print "('fe_mesh: get_n_subdomains -- &
            & the mesh has negative code for some elements. This is not allowed')"
       stop
    end if
  end function get_n_subdomains

  !> @brief Computes the distance between nodes i and j in mesh m.
  pure real(pr) function get_distance(m,i,j)
    type(mesh), intent(in) :: m
    integer, intent(in) :: i,j

    real(pr), dimension(m%dim) :: x,y

    x = m%node(1:m%dim,i)
    y = m%node(1:m%dim,j)

    get_distance = sqrt( dot_product(x-y,x-y) )
  end function get_distance

  !> @brief Get the dim+1 nodes of element number i in mesh m.
  pure function get_nodes(m,i)
    type(mesh), intent(in) :: m
    integer, intent(in) :: i
    integer, dimension(m%dim+1) :: get_nodes

    get_nodes = m%elt(1:m%dim+1,i)
  end function get_nodes

  !> @brief Get the coordinates of vertex number i from mesh m.
  pure function get_x(m,i)
    type(mesh), intent(in) :: m
    integer, intent(in) :: i
    real(pr), dimension(m%dim) :: get_x

    get_x = m%node(1:m%dim,i)
  end function get_x

  !> @brief Returns the measure of an element.
  !> @details Distance in 1D, surface in 2D and volume in 3D.
  pure function get_measure(m,i) result (meas)
    type(mesh), intent(in) :: m
    integer, intent(in)    :: i
    real(pr)               :: meas

    real(pr), dimension(m%dim,m%dim+1) :: x
    real(pr), dimension(m%dim,m%dim) :: a
    integer :: j

    x(1:m%dim,1:m%dim+1) = m%node(1:m%dim, m%elt(1:m%dim+1,i))
    do j = 1,m%dim
       a(:,j) = x(:,1+j) - x(:,1)
    end do

    ! Dans ts les cas, meas = det(A) / fact(dim)
    select case(m%dim)
    case(1)
       meas = abs(a(1,1))
    case(2)
       meas = abs(a(1,1)*a(2,2)-a(2,1)*a(1,2))
       meas = meas/2._pr
    case(3)
       meas = abs(&
            a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1) + a(1,3)*a(2,1)*a(3,2) &
            - a(1,1)*a(2,3)*a(3,2) - a(1,2)*a(2,1)*a(3,3) - a(1,3)*a(2,2)*a(3,1))
       meas = meas/6._pr
    case default
       meas = 0._pr
    end select

  end function get_measure

  !> @brief Get the code of the element i
  pure integer function get_elt_code(m,i)
    type(mesh), intent(in) :: m
    integer, intent(in) :: i
    get_elt_code = m%elt_code(i)
  end function get_elt_code

  !> @brief Get the code of the node i
  pure integer function get_node_code(m,i)
    type(mesh), intent(in) :: m
    integer, intent(in) :: i
    get_node_code = m%node_code(i)
  end function get_node_code

  !> @brief Read a mesh from file file_name
  subroutine read_mesh (m, file_name, scale)
    type(mesh), intent(out)      :: m
    character(len=*), intent(in) :: file_name
    real(pr), optional           :: scale

    integer :: i
    character(len=len(file_name)) :: suffix
    logical :: test

    inquire(file=file_name, exist=test)

    if (.not.test) then
       print "('File ',A,' doesn''t exists')", trim(file_name)
       stop
    end if

    i = index(file_name,'.',back=.true.)

    suffix = file_name(i+1:len(file_name))

    select case(trim(suffix))
    case('mesh')
       call read_mesh_mesh(m,file_name)
    case('node','ele')
       call read_mesh_node_ele(m,file_name)
    case default
       print "('fe_mesh: read_mesh -- unknown file format for file ',A)", &
            trim(file_name)
       stop
    end select

    if (present(scale)) then
       print "('rescaling the mesh by a factor of : ',E13.5)", scale
       do i = 1,m%n_nodes
          m%node(1:m%dim,i) = scale*m%node(1:m%dim,i)
       end do
    end if

  end subroutine read_mesh

  !> @brief Write the mesh to the file file_name, together with the
  !> scalar data u if it is given (vtk only)
  subroutine write_mesh (m, file_name, u)
    type(mesh), intent(in)       :: m
    character(len=*), intent(in) :: file_name
    real(pr), dimension(:), optional :: u

    integer :: i
    character(len=len(file_name)) :: suffix
    logical :: test

    inquire(file=file_name, exist=test)

    if (.not.test) then

       i = index(file_name,'.',back=.true.)

       suffix = file_name(i+1:len(file_name))

       select case(trim(suffix))
       case('mesh')
          call write_mesh_mesh(m,file_name)
       case('vtk')
          if (present(u)) then
             call write_mesh_vtk(m,file_name,u)
          else
             call write_mesh_vtk(m,file_name)
          end if
       case default
          print "('fe_mesh: write_mesh -- unknown file format for file ',A)", &
               trim(file_name)
          stop
       end select

    else

       print "('File ',A,' already exists... PASS')", trim(file_name)

    end if

  end subroutine write_mesh


  !> @brief Lit les fichiers file_name.node et file_name.ele obtenus
  !> avec 'triangle'
  subroutine read_mesh_node_ele(m,file_name)
    type(mesh), intent(out) :: m
    character(len=*), intent(in) :: file_name

    integer :: unit, i_dot, n_attr, n_bound, attr, code, n_per_tr
    integer :: i1,i2,i3, i,j
    real(pr) :: x,y

    i_dot = index(file_name,'.',back=.true.)

    unit = 10
    open(unit, FILE=file_name(1:i_dot)//'node', ACTION='read', STATUS='old')
    print "('lecture du fichier ',A)", file_name(1:i_dot)//'node'

    read(unit,*) m%n_nodes,m%dim, n_attr, n_bound
    allocate(m%node(m%dim, m%n_nodes))
    allocate(m%node_code(m%n_nodes))
    do i = 1,m%n_nodes
       if (n_attr>0.and.n_bound>0) then
          ! Lecture avec attributs et code
          read(unit,*) j,x,y,attr,code
       elseif (n_attr>0.and.n_bound==0) then
          ! Lecture avec attributs et sans code
          read(unit,*) j,x,y,attr
          code = 0 ! pas de code -> 0
       elseif (n_attr==0.and.n_bound>0) then
          ! Lecture sans attribut et avec code
          read(unit,*) j,x,y,code
       else
          ! Lecture sans attribut ni code
          read(unit,*) j,x,y
          code = 0
       end if
       m%node(1,j) = x
       m%node(2,j) = y
       m%node_code(j) = code
    end do

    close(unit)

    open(unit, FILE=file_name(1:i_dot)//'ele', ACTION='read', STATUS='old')
    print "('lecture du fichier ',A)", file_name(1:i_dot)//'ele'

    read (unit,*) m%n_elts,n_per_tr,n_attr
    if (n_per_tr/=m%dim+1) then
       print "('fe_mesh: read_mesh_node_ele -- &
            &probleme de description du maillage, fichier ', A)", &
            file_name(1:i_dot)//'ele'
       stop
    end if
    allocate(m%elt(m%dim+1, m%n_elts), m%elt_code(m%n_elts))
    if (n_attr>0) then
       ! Lecture avec un code par élément.
       do i = 1,m%n_elts
          read (unit,*) j, i1,i2,i3, code
          m%elt(1,j) = i1
          m%elt(2,j) = i2
          m%elt(3,j) = i3
          m%elt_code(j) = code
       end do
    else
       ! Pas de code par élément --> 0
       do i = 1,m%n_elts
          read (unit,*) j, i1,i2,i3
          m%elt(1,j) = i1
          m%elt(2,j) = i2
          m%elt(3,j) = i3
          m%elt_code(j) = 0
       end do
    end if

    close(unit)

  end subroutine read_mesh_node_ele


  !> @brief Reads a mesh from the file file_name
  subroutine read_mesh_mesh (m,file_name)
    character(len=*), intent(in) :: file_name
    type(mesh), intent(out)      :: m

    integer :: unit, i,zero
    character(len=8) :: line

    unit = 10
    open(unit,FILE=trim(file_name), ACTION='read', STATUS='old')
    print "('lecture du fichier ',A)", trim(file_name)

    ! read <dimension>
    line = ' '
    do while (.not.line(1:3)=='Dim')
       read(unit, '(A8)') line
       line = adjustl(line)
    end do
    read (unit, *) m%dim
    if(m%dim<1.or.m%dim>3) then
       print "('fe_mesh: read_mesh_mesh -- &
            &la dimension devrait etre 1,2 ou 3 au lieu de ',I0)", m%dim
       stop
    end if
    m%dim=2 !Dim=3 ds le fichier mais on veut 2

    ! read nodes : <n_nodes>, <node>, <nodecode>
    do while (.not.line(1:4)=='Vert')
       read (unit, '(A8)') line
       line = adjustl(line)
    end do
    read (unit, *) m%n_nodes
    if (m%n_nodes<2) then
       print "('fe_mesh: read_mesh_mesh -- &
            &unable to read build mesh with ',I0,' nodes')", m%n_nodes
       stop
    end if
    m%n_elts = 0
    allocate(m%node(m%dim, m%n_nodes))
    allocate(m%node_code(m%n_nodes))
    do i=1,m%n_nodes
       ! On essaie de lire un code de noeud.
       read (unit,*) m%node(1:m%dim,i), zero, m%node_code(i)


    end do

    ! read simplices : <element>, indices des noeuds et numero de sous
    ! domaine
    select case(m%dim)
    case(1)
       do while (.not.(line(1:5)=='Edges'))
          read (unit, '(A8)') line
          line = adjustl(line)
       end do
    case(2)
       do while (.not.(line(1:5)=='Trian'))
          read (unit, '(A8)') line
          line = adjustl(line)
       end do
    case(3)
       do while (.not.(line(1:5)=='Tetra'))
          read (unit, '(A8)') line
          line = adjustl(line)
       end do
    end select

    read (unit, *) m%n_elts
    if (m%n_elts<m%dim) then
       print "('fe_mesh: read_mesh_mesh -- &
            &Unable to read mesh ',I0,' edges')", m%n_elts
       stop
    end if
    allocate(m%elt(1+m%dim,m%n_elts), m%elt_code(m%n_elts))
    do i=1,m%n_elts
       read (unit, *) m%elt(1:1+m%dim, i), m%elt_code(i)
    end do

    ! read faces, for 3D meshes that have declared triangles


    close(unit)

  end subroutine read_mesh_mesh

  subroutine chg_lbl(m,xymin,xymax,dir)
    type(mesh),intent(inout)::m
    real(pr),intent(in)::xymin,xymax
    integer,intent(in)::dir
    integer::i

    do i=1,m%n_nodes
      ! Si l'élément est sur le bord inférieur code=1, bord supérieur code=3 sinon code=0
      if (abs(m%node(dir,i)-xymin)<=0.000001)  then !bord inférieur: y=-145
        m%node_code(i)=1
      else if (abs(m%node(dir,i)-xymax)<=0.000001) then!bord supérieur y=35
        m%node_code(i)=3
      else
        m%node_code(i)=0
      end if
    end do
  end subroutine


  !> @brief Free the memory.
  subroutine deallocate_mesh(m)
    type(mesh), intent(inout) :: m

    m%dim = 0
    m%n_nodes = 0
    m%n_elts = 0

    deallocate(m%node)
    if (allocated(m%node_code)) deallocate(m%node_code)
    if (allocated(m%elt)) deallocate(m%elt)
    if (allocated(m%elt_code)) deallocate(m%elt_code)
    if (allocated(m%face)) deallocate(m%face)
    if (allocated(m%face_code)) deallocate(m%face_code)

  end subroutine deallocate_mesh


  !> @brief Prints some basic information about the mesh.
  subroutine print_mesh(m)
    type(mesh), intent(in) :: m

    integer :: i

    print "('Maillage en dimension : ',I0)", m%dim
    print "('Nombre de noeuds      : ',I0)", m%n_nodes
    print "('Nombre d''elements     : ',I0)", m%n_elts
    do i = minVal(m%elt_code),maxVal(m%elt_code)
       print "('            code ',I4,' : ',I0)", i,count(m%elt_code==i)
    end do

  end subroutine print_mesh

  !> @brief Writes the mesh m to the specified fortran file using the
  !> .mesh format specification from medit.
  subroutine write_mesh_mesh(m,file_name)
    type(mesh), intent(in)  :: m
    character(len=*), intent(in) :: file_name

    integer :: unit=10, i
    open(unit, FILE=trim(file_name), ACTION='write', STATUS='new')
    print "('writing file ',A)", trim(file_name)

    ! Introduction
    write (unit,*) 'MeshVersionFormatted 1'
    write (unit,*) 'Dimension'
    write (unit,*) m%dim
    write (unit,*) 'Vertices'
    write (unit,*) m%n_nodes
    write (unit,*) (m%node(1:m%dim,i), m%node_code(i), i=1,m%n_nodes)
    select case(m%dim)
    case(1)
       write (unit,*) '# Set of edges'
       write(unit,*) 'Edges'
    case(2)
       write (unit,*) '# Set of triangles'
       write(unit,*) 'Triangles'
    case(3)
       write (unit,*) '# Set of tetrahedra'
       write(unit,*) 'Tetrahedra'
    end select
    write (unit,*) m%n_elts
    do i = 1,m%n_elts
       write(unit,*) m%elt(1:m%dim+1,i),m%elt_code(i)
    end do

  end subroutine write_mesh_mesh

  !> @brief Writes the mesh m to the specified fortran file using the
  !> legacy vtk format
  !> www.vtk.org/pdf/file-formats.pdf
  !> @param[in] m the mesh
  !> @param[in] file_name the file name
  !> @param[optional] a scalar data to append to the mesh file (ie solution)
  subroutine write_mesh_vtk(m,file_name,u)
    type(mesh), intent(in)  :: m
    character(len=*), intent(in) :: file_name
    real(pr), dimension(:), optional :: u

    integer :: j, unit = 10

    open(unit, FILE=trim(file_name), ACTION='write', STATUS='new')
    print "('writing file ',A)", trim(file_name)

    ! Keywords are case insensitive

    ! 1/ File Version identifier
    write (unit,'(A)') '# vtk DataFile Version 3.0.'

    ! 2/ Header
    write (unit,'(A)') 'Written by fe_mesh:write_vtk'

    ! 3/ File Format
    write (unit,'(A)') 'ASCII'
    ! 4/ Dataset structure
    write (unit,'(A)') 'DATASET UNSTRUCTURED_GRID'
    ! 4-1/ POINTS qui ont tjs 3 coordonnees - meme en 1D et 2D
    write (unit,'(A7,I9,1X,A6)') 'POINTS ',m%n_nodes,'float'
    select case(m%dim)
    case(1)
       write (unit,'(3E23.15)') (m%node(1,j),0._pr,0._pr, j=1,m%n_nodes)
    case(2)
       write (unit,'(3E23.15)') (m%node(1:m%dim,j),0._pr, j=1,m%n_nodes)
    case(3)
       write (unit,'(3E23.15)') (m%node(1:m%dim,j), j=1,m%n_nodes)
    end select
    ! 4-2/ CELLS
    write (unit,'(A5,2I9)') 'CELLS', m%n_elts, (1+m%dim+1)*m%n_elts
    select case(m%dim)
    case(1)
       write (unit,'(3I9)') (m%dim+1,m%elt(1:m%dim+1,j)-1, j=1,m%n_elts)
    case(2)
       write (unit,'(4I9)') (m%dim+1,m%elt(1:m%dim+1,j)-1, j=1,m%n_elts)
    case(3)
       do j = 1,m%n_elts
          write (unit,'(5I9)') m%dim+1,m%elt(1:m%dim+1,j)-1
       end do
    end select
    ! 4-1/ CELL TYPES
    write (unit,'(A10,I9)') 'CELL_TYPES', m%n_elts
    select case(m%dim)
    case(1)
       write (unit,'(I2)') (VTK_LINE, j=1,m%n_elts)
    case(2)
       write (unit,'(I2)') (VTK_TRIANGLE, j=1,m%n_elts)
    case(3)
       write (unit,'(I2)') (VTK_TETRA, j=1,m%n_elts)
    end select

    ! 5.1/ Dataset attributes - Les codes des elements
    write (unit,'(A,I9)') 'CELL_DATA',m%n_elts
    write (unit,'(A,1X,A,1X,A,I2)') 'SCALARS','elt_code','int',1
    write (unit,'(A,1X,A)') 'LOOKUP_TABLE','DEFAULT'
    write (unit,'(8I9)') m%elt_code

    ! 5.2/ Dataset -- Codes of the nodes if provided
    if (allocated(m%node_code)) then
       write (unit,'(A,I9)') 'POINT_DATA',m%n_nodes
       write (unit,'(A,1X,A,1X,A,I2)') 'SCALARS','node_code','int', 1
       write (unit,'(A,1X,A)') 'LOOKUP_TABLE','DEFAULT'
       write (unit,'(8(I0,1X))') m%node_code
    end if

    ! 6/ Dataset if provided
    if (present(u)) then
       write (unit,'(A,1X,A,1X,A,I2)') 'SCALARS','ux','float', 1
       write (unit,'(A,1X,A)') 'LOOKUP_TABLE','DEFAULT'
       write (unit,'(1E13.6)') (u(j), j=1,m%n_nodes)

       if(size(u)==2*m%n_nodes) then
         write (unit,'(A,1X,A,1X,A,I2)') 'SCALARS','uy','float', 1
         write (unit,'(A,1X,A)') 'LOOKUP_TABLE','DEFAULT'
         write (unit,'(1E13.6)') (u(j), j=m%n_nodes+1,2*m%n_nodes)
      end if


    end if

    close(unit)

  end subroutine write_mesh_vtk


  ! vector product for vectors in R^3
  function vec_product(u,v)

    real(pr), dimension(3), intent(in) :: u,v
    real(pr), dimension(3) :: vec_product

    vec_product(1) = u(2)*v(3) - u(3)*v(2)
    vec_product(2) = u(3)*v(1) - u(1)*v(3)
    vec_product(3) = u(1)*v(2) - u(2)*v(1)

  end function vec_product

  subroutine strain(m,u)
    type(mesh),intent(inout)::m
    real(pr),dimension(:),intent(in)::u
    integer::i

    do i=1,m%n_nodes
      m%node(1,i)=m%node(1,i)+u(i)
      m%node(2,i)=m%node(2,i)+u(m%n_nodes+i)
    end do
  end subroutine

end module fe_mesh_quentin
