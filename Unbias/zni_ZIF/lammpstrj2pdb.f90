!ITEM: TIMESTEP
!0
!ITEM: NUMBER OF ATOMS
!816
!ITEM: BOX BOUNDS xy xz yz pp pp pp
!-1.9741981549897282e+01 3.0372448059976723e+01 -1.0016673770959233e+01
!-1.3891208734959575e+01 2.8628128422373553e+01 -1.0032159718961333e+01
!2.3768098780503374e-01 2.4811788012194945e+01 -1.4187558312586022e+01
!ITEM: BOX BOUNDS xy xz yz
!xlo_bound xhi_bound xy
!ylo_bound yhi_bound xz
!zlo_bound zhi_bound yz
!ITEM: ATOMS element xs ys zs
!Zn 0.250224 0.125095 0.375126
MODULE vector_module
! Chapman p. 544
! u%x = ARRAY(1)
! u%y = ARRAY(2)
! u%z = ARRAY(3)
 IMPLICIT NONE
 TYPE  :: vector
  SEQUENCE
  REAL :: x
  REAL :: y
  REAL :: z
 END TYPE vector
 CONTAINS
  type (vector) function vector_unitary(v1)
   implicit none
   type(vector), intent(in) :: v1
   real                     :: norma
   norma=absvec(v1)
   vector_unitary%x = v1%x/norma
   vector_unitary%y = v1%y/norma
   vector_unitary%z = v1%z/norma
  end function vector_unitary
  type (vector) function vector_scale(v1, r )
   implicit none
   type(vector), intent(in) :: v1
   real,intent(in)          :: r
   vector_scale%x = v1%x*r
   vector_scale%y = v1%y*r
   vector_scale%z = v1%z*r
  end function vector_scale
!
   TYPE (vector) FUNCTION vector_add (v1,v2)
    implicit none
    TYPE (vector), INTENT(IN) :: v1
    TYPE (vector), INTENT(IN) :: v2
    vector_add%x = v1%x + v2%x
    vector_add%y = v1%y + v2%y
    vector_add%z = v1%z + v2%z
   END FUNCTION vector_add
!
   TYPE (vector) FUNCTION vector_sub (v1,v2)
    IMPLICIT NONE
    TYPE (vector), INTENT(IN) :: v1
    TYPE (vector), INTENT(IN) :: v2
    vector_sub%x = v1%x - v2%x
    vector_sub%y = v1%y - v2%y
    vector_sub%z = v1%z - v2%z
   END FUNCTION vector_sub
!
   TYPE (vector) FUNCTION cross(a,b)
    IMPLICIT NONE
    TYPE (vector), INTENT (in) :: a, b
    cross%x = a%y * b%z - a%z * b%y
    cross%y = a%z * b%x - a%x * b%z
    cross%z = a%x * b%y - a%y * b%x
   END FUNCTION cross
!
   real function dot(a,b)
    implicit none
    type(vector), intent(in) :: a,b
    dot = a%x * b%x + a%y*b%y + a%z*b%z
   end function dot
!
   REAL FUNCTION absvec(a)
    TYPE (vector), INTENT (in) :: a
    absvec = sqrt(a%x**2 + a%y**2 + a%z**2)
   END FUNCTION absvec
!
   real function dist2vectors(a,b)
    type(vector), intent(in) :: a,b
    dist2vectors = absvec( vector_sub(b,a) )
   end function  dist2vectors
!
   real function angle2vectors(a,b)
    ! Vectors are the  r_ij vector
    !                 (r_ij, r_jk)
    type(vector), intent(in) :: a,b
    angle2vectors = acos( dot(a,b)/(absvec(a)*absvec(b)) )
   end function  angle2vectors
!
   real function angle3vectors_jik(r_i,r_j,r_k)
    type(vector), intent(in) :: r_i,r_j,r_k
    type(vector)             :: r_ij,r_ik,r_jk
    r_ij = vector_sub(r_i,r_j)
    r_jk = vector_sub(r_j,r_k)
    r_ik = vector_sub(r_i,r_k)
    angle3vectors_jik = acos(dot(r_ij,r_ik)/(absvec(r_ij)*absvec(r_ik)))
   end function  angle3vectors_jik
!  
   real function dihedral_angle4vectors_ijkl(r_i,r_j,r_k,r_l)
    type(vector), intent(in) :: r_i,r_j,r_k,r_l
    type(vector)             :: r_ij,r_jk,r_kl
    r_ij = vector_sub(r_i,r_j)
    r_jk = vector_sub(r_j,r_k)
    r_kl = vector_sub(r_k,r_l)
    dihedral_angle4vectors_ijkl = acos(dot(cross(r_ij,r_jk),cross(r_jk,r_kl))/(absvec(cross(r_ij,r_jk))*absvec(cross(r_jk,r_kl))))
   end function  dihedral_angle4vectors_ijkl
!
END MODULE vector_module
!
program lammpstrj2pdb
 use iso_fortran_env
 use vector_module
 implicit none
 integer             :: ierr,i,j,k,l
 integer             :: timestep = 0
 integer             :: n_atoms = 0
 integer             :: n_vacs = 0
 real                :: vr(3,3),rv(3,3),cell_0(6)
 real :: xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz,lx,ly,lz,cosa,cosb,cosg
 real :: xlo_bound,ylo_bound,zlo_bound,xhi_bound,yhi_bound,zhi_bound
 real,parameter      :: r_min_criteria_connectivity=0.15
 real                :: r0=r_min_criteria_connectivity
 real                :: ratom(3),rouratom(3),r
 real,allocatable    :: DistanceMatrix(:,:)
 logical,allocatable :: ConnectedAtoms(:,:)
 real,allocatable    :: xcrystal(:,:),xcarte(:,:)
 character(len=2),   allocatable:: label(:)
 character(len=2),   allocatable:: label_for_distance(:,:)
 character(len=2),   allocatable:: label_for_angle(:,:)
 character(len=2),   allocatable:: label_for_dihedral(:,:)
 character(len=80)   :: line
 character(len=40)   :: keyword
 character(len=1)    :: chain = " "
 character(len=3)    :: resname = "MOL"
 integer             :: ires = 0
 integer             :: n_angles_for_study = 0, n_dihedrals_for_study = 0
 integer             :: n_distances_for_study  = 0
 character(len=10)   :: space_group = "P 1"
 real,parameter      :: pi=acos(-1.0)
 real,parameter      :: radtodeg = 180.0/pi
 real,parameter      :: degtorad = pi/180.0
! arguments in line
 character(len=100),dimension(:), allocatable :: args
 integer             :: num_args = 0
 logical             :: remove_shells = .false.
 logical             :: modify_vacants_by_dummy_atoms = .false.
 logical             :: remove_adsorbate = .false.
 logical             :: analyse = .false., analyse_bend_angle =.false., recentre_in_com=.false.
 logical             :: analyse_dihedral_angles=.false., analyse_RDF = .false.
 logical,allocatable :: printable(:)
 type(vector)        :: com
! ------------------------------------------------------------------------------
 num_args = command_argument_count()
 allocate(args(num_args))
 do i = 1, num_args
  call get_command_argument(i,args(i))
 end do
 write(6,'(a,1x,i2)')'Arguments:',command_argument_count()
 write(6,'(a)')( args(i),i=1,num_args)
 if(num_args >= 1) then
  do i=1,num_args
   select case(args(i))
    case("-com")
     recentre_in_com=.true.
    case("-rads")
     remove_adsorbate = .true.
    case ("-rsh","--remove-shells")
     remove_shells = .true.
    case ("-rvac","--remove-vacants")
     modify_vacants_by_dummy_atoms = .true.
    case ("-ad","--analyse_bond_distance")
     analyse=.true.
     analyse_RDF=.true.
     read(args(i+1),'(i1)') n_distances_for_study
     allocate( label_for_distance(2,1:n_distances_for_study))
     do j=1,n_distances_for_study
      write(6,'(3(a2,1x))') args(i+1+j), args(i+2+j)
      read(args(i+1+j),'(a2)') label_for_distance(1,j)
      read(args(i+2+j),'(a2)') label_for_distance(2,j)
      open(122,file="distance.txt")
     end do
    case ("-aa", "--analyse_bend_angle")
     analyse=.true.
     analyse_bend_angle=.true.
     read(args(i+1),'(i1)') n_angles_for_study
     allocate( label_for_angle(1:3,1:n_angles_for_study))
     do j=1,n_angles_for_study 
      write(6,'(3(a2,1x))') args(i+1+j), args(i+2+j), args(i+3+j)
      read(args(i+1+j),'(a2)') label_for_angle(1,j)
      read(args(i+2+j),'(a2)') label_for_angle(2,j)
      read(args(i+3+j),'(a2)') label_for_angle(3,j)
      open(123,file="angles.txt")
     end do
    case("-at")
     analyse=.true.
     analyse_dihedral_angles=.true.
     read(args(i+1),'(i1)') n_dihedrals_for_study
     allocate( label_for_dihedral(1:4,1:n_dihedrals_for_study))
     do j=1,n_dihedrals_for_study
      write(6,'(4(a2,1x))') args(i+1+j), args(i+2+j), args(i+3+j), args(i+4+j)
      read(args(i+1+j),'(a2)') label_for_dihedral(1,j)
      read(args(i+2+j),'(a2)') label_for_dihedral(2,j)
      read(args(i+3+j),'(a2)') label_for_dihedral(3,j)
      read(args(i+4+j),'(a2)') label_for_dihedral(4,j)
      open(124,file="dihedrals.txt")
     end do
    case default
     remove_shells = .false.
     modify_vacants_by_dummy_atoms = .false.
   end select
  end do
 end if
 open(100,file="out.pdb")
 write(100,'(a)')"AUTHOR    Salvador R. G. Balestra"
 write(100,'(a)')"TITLE     Data from LAMMPS simulations (lammpstrj2pdb.f90 code)"
 k=0
 read_lammpstrj:do
  read(5,'(a)',iostat=ierr)line
  if(ierr/=0) exit read_lammpstrj
  if(line(1:5)=="ITEM:")then
   read(line,'(5x,1x,a)') keyword
   select case (keyword)
    case('TIMESTEP')
     read(5,*) timestep
     k=k+1
    case('NUMBER OF ATOMS')
     read(5,*) n_atoms
     allocate(xcrystal(3,n_atoms))
     allocate(xcarte(3,n_atoms))
     allocate(label(n_atoms))
     allocate(printable(n_atoms))
     allocate(ConnectedAtoms(n_atoms,n_atoms))
     allocate(DistanceMatrix(n_atoms,n_atoms))
     DistanceMatrix=0.0
     ConnectedAtoms=.false.
     printable=.true.
    case('BOX BOUNDS xy xz yz pp pp pp')
     read(5,*)xlo_bound,xhi_bound,xy
     read(5,*)ylo_bound,yhi_bound,xz
     read(5,*)zlo_bound,zhi_bound,yz
     xlo=xlo_bound-MIN(0.0,xy,xz,xy+xz)
     xhi=xhi_bound-MAX(0.0,xy,xz,xy+xz)
     ylo=ylo_bound-MIN(0.0,yz)
     yhi=yhi_bound-MAX(0.0,yz)
     zlo=zlo_bound
     zhi=zhi_bound
     lx=xhi-xlo
     ly=yhi-ylo
     lz=zhi-zlo
     cell_0(1)=lx
     cell_0(2)=sqrt(ly*ly+xy*xy)
     cell_0(3)=sqrt(lz*lz+xz*xz+yz*yz)
     cell_0(4)=radtodeg*acos((xy*xz+ly*yz)/(cell_0(2)*cell_0(3)))
     cell_0(5)=radtodeg*acos(xz/cell_0(3))
     cell_0(6)=radtodeg*acos(xy/cell_0(2))
     call cell(rv,vr,cell_0)
    ! ITEM: BOX BOUNDS pp pp pp
    case('BOX BOUNDS pp pp pp')
     read(5,*)xlo_bound,xhi_bound
     read(5,*)ylo_bound,yhi_bound
     read(5,*)zlo_bound,zhi_bound
     xlo=xlo_bound
     xhi=xhi_bound
     ylo=ylo_bound
     yhi=yhi_bound
     zlo=zlo_bound
     zhi=zhi_bound
     xy=0.0 ; xz=0.0 ; yz=0.0
     lx=xhi-xlo
     ly=yhi-ylo
     lz=zhi-zlo
     cell_0(1)=lx
     cell_0(2)=sqrt(ly*ly+xy*xy)
     cell_0(3)=sqrt(lz*lz+xz*xz+yz*yz)
     cell_0(4)=radtodeg*acos((xy*xz+ly*yz)/(cell_0(2)*cell_0(3)))
     cell_0(5)=radtodeg*acos(xz/cell_0(3))
     cell_0(6)=radtodeg*acos(xy/cell_0(2))
     call cell(rv,vr,cell_0)
    case('ATOMS element xs ys zs')
     ! Format from:
     ! ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_A4.pdf
     line(1:80)=" "
     write(line(1:6),'(a6)')"MODEL "
     write(line(7:14),'(i8)')k
     write(100,'(a80)')line(1:80)
     !write(100,'(a6,1x,9(f14.7,1x))')'REMARK',xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz
     write(100,'(a6,1x,a,1x,i6,1x,a,i6)')'REMARK','Atoms:',n_atoms,'Step:',timestep
     write(100,'(a6,1x,a,1x,f14.7)')'REMARK','Volume:',volume(rv)
     line(1:80)=" "
     write(line( 1: 6),'(a)')    "CRYST1"
     write(line( 7:15),'(f9.3)') cell_0(1)
     write(line(16:24),'(f9.3)') cell_0(2)
     write(line(25:33),'(f9.3)') cell_0(3)
     write(line(34:40),'(f7.2)') cell_0(4)
     write(line(41:47),'(f7.2)') cell_0(5)
     write(line(48:54),'(f7.2)') cell_0(6)
     write(line(56:66),'(a10)')  adjustL( space_group )
     write(100,'(a80)')line(1:80)
     call cell(rv,vr,cell_0)
     do i=1,n_atoms
      read(5,*) label(i),(xcrystal(j,i),j=1,3)
      ! crystallyne -> cartesians
      forall ( j=1:3 )
       xcarte(j,i)=rv(j,1)*xcrystal(1,i) + rv(j,2)*xcrystal(2,i) + rv(j,3)*xcrystal(3,i)
      end forall
     end do
     ! -------- select cases for remoe shells and etc.
     if(recentre_in_com)then
      call compute_com(com)
      call recentre(com)
     end if
     if(remove_adsorbate)then
      call to_remove_adsorbate()
     end if
     if(remove_shells)then
       call to_measure_topoly()
       call remove_shells_()
      ! n_printable_atoms
     end if
     if(modify_vacants_by_dummy_atoms)then
       call to_measure_topoly()
       call to_modify_vacancy(n_vacs)
      ! n_printable_atoms
     end if
     if(analyse)then
      call to_measure_topoly()
      call checkatom_type()
      if(analyse_RDF)then
       call write_distances()
      end if
      if(analyse_bend_angle)then
       call write_angles()
      end if
      if ( analyse_dihedral_angles ) then
       call write_dihedrals()
      end if
     end if
     do i=1,n_atoms
      if(printable(i))then
      line(1:80)=" "
      write(line( 1: 6),'(a6)')  adjustL("ATOM  ")
      write(line( 7:11),'(i5)')  n_vacs+i
      write(line(13:16),'(a)')   label(i)
      write(line(17:17),'(a)')   chain
      write(line(18:20),'(a)')   adjustR(label(i))
      write(line(22:22),'(a)')   chain
      write(line(23:26),'(i4)')  ires
      write(line(27:27),'(a)')   chain
      write(line(31:38),'(f8.3)')xcarte(1,i)
      write(line(39:46),'(f8.3)')xcarte(2,i)
      write(line(47:54),'(f8.3)')xcarte(3,i)
      write(line(55:60),'(f6.2)')1.0
      write(line(61:66),'(f6.2)')0.0
      write(line(77:78),'(a2)')  adjustR(label(i))
      write(line(79:80),'(a2)')  "+0"
      write(100,'(a80)')line(1:80)
      ! OLD:
      !write(100,'(a6,i5,1x,a2,3x,a4,1x,i4,4x,3f8.3,2f6.2,10x,a2)') &
      !'ATOM  ',i,label(i),'MOL ',0,(xcarte(j,i),j=1,3),0.0,0.0,label(i)
      end if
     end do
     write(100,'(a6)')'ENDMDL'
     deallocate(xcrystal)
     deallocate(xcarte)
     deallocate(label)
     deallocate(ConnectedAtoms)
     deallocate(DistanceMatrix)
     deallocate(printable)
    case default
     STOP 'Keyword ??'
   end select
  end if
 end do read_lammpstrj
contains
!------------------------------------------------------------------------------
subroutine write_distances()
 implicit none
 integer   :: i,j,k
 real      :: RDF_max = 13.0
 real      :: minim_dist
 do i=1,n_atoms
  if(label(i)==label_for_distance(1,1))then
   minim_dist = RDF_max
   do j=i+1,n_atoms 
    if(label(j)==label_for_distance(2,1).and.DistanceMatrix(i,j)<=minim_dist)then
     minim_dist = DistanceMatrix(i,j)
     k = i
     !write(122,'(f14.7,2(i5,1x))')DistanceMatrix(i,j),i,j
    end if
   end do
   write(122,'(f14.7,2(i5,1x))') minim_dist,i,k
  end if
 end do
end subroutine write_distances
!
subroutine write_angles()
 implicit none
 integer           :: i,j,k
 type(vector)      :: r_i,r_j,r_k
 real, parameter   :: pi = ACOS(-1.0)
 real              :: DEGTORAD = pi/180.0
 do i=1,n_atoms
  if(label(i)==label_for_angle(2,1))then                                               ! i -> Zn
  do j=1,n_atoms
   if(i/=j.and.ConnectedAtoms(i,j).and.label(j)==label_for_angle(1,1))then             ! j -> N
    search_last: do k=j+1,n_atoms
     if(k/=i.and.k/=j.and.ConnectedAtoms(i,k).and.label(k)==label_for_angle(3,1))then  ! k -> N
       r_i%x = xcarte(1,i) 
       r_i%y = xcarte(2,i)
       r_i%z = xcarte(3,i)
       !
       r_j%x = xcarte(1,j) 
       r_j%y = xcarte(2,j)
       r_j%z = xcarte(3,j)
       !
       r_k%x = xcarte(1,k)
       r_k%y = xcarte(2,k)
       r_k%z = xcarte(3,k)
       if( dist2vectors(r_i,r_j) < 3.0 .and. dist2vectors(r_i,r_k) < 3.0 .and.dist2vectors(r_j,r_k) < 4.0) then
       write(6,'(3(a2,1x),3(i5,1x),f14.7)') label(i),label(j),label(k),i,j,k,angle3vectors_jik(r_i,r_j,r_k)/DEGTORAD
       write(123,'(f14.7, 3(i5,1x))') angle3vectors_jik(r_i,r_j,r_k)/DEGTORAD, i,j,k
       end if
       exit search_last
     end if
    end do search_last
   end if
  end do 
  end if
 end do
 return
end subroutine write_angles
!
subroutine write_dihedrals()
 implicit none
 integer           :: i,j,k,l
 type(vector)      :: r_i,r_j,r_k,r_l
 real, parameter   :: pi = ACOS(-1.0)
 real              :: DEGTORAD = pi/180.0
 do i=1,n_atoms
  if(label(i)==label_for_dihedral(1,1))then
  do j=1,n_atoms
   if(i/=j.and.ConnectedAtoms(i,j).and.label(j)==label_for_dihedral(2,1))then
    do k=1,n_atoms
     if(k/=i.and.k/=j.and.ConnectedAtoms(j,k).and.label(k)==label_for_dihedral(3,1))then
      do l=1,n_atoms
       if(l/=i.and.l/=j.and.l/=k.and.ConnectedAtoms(k,l).and.label(l)==label_for_dihedral(4,1))then
       r_i%x = xcarte(1,i)
       r_i%y = xcarte(2,i)
       r_i%z = xcarte(3,i)
       !
       r_j%x = xcarte(1,j)
       r_j%y = xcarte(2,j)
       r_j%z = xcarte(3,j)
       !
       r_k%x = xcarte(1,k)
       r_k%y = xcarte(2,k)
       r_k%z = xcarte(3,k)
       !
       r_l%x = xcarte(1,l)
       r_l%y = xcarte(2,l)
       r_l%z = xcarte(3,l)
       if( dist2vectors(r_i,r_j) < 3.0 .and. dist2vectors(r_j,r_k) < 3.0 .and.dist2vectors(r_k,r_l) < 3.0) then
       write(6,'(4(a2,1x),4(i5,1x),f14.7)') &
        label(i),label(j),label(k),label(l),i,j,k,l,dihedral_angle4vectors_ijkl(r_i,r_j,r_k,r_l)/DEGTORAD
       write(124,'(f14.7, 4(i5,1x))') dihedral_angle4vectors_ijkl(r_i,r_j,r_k,r_l)/DEGTORAD, i,j,k,l
       end if
       end if
      end do
       ! }}
     end if
    end do
   end if
  end do
  end if
 end do
 return
end subroutine write_dihedrals
!
subroutine to_remove_adsorbate()
 implicit NONE
 integer :: i,j,k,l
 do i=1,n_atoms
   if(label(i)=="C ".or.label(i)==" C") printable(i)=.false.
 end do
 RETURN
end subroutine to_remove_adsorbate
!
subroutine to_measure_topoly()
 implicit NONE
 integer :: i,j,k,l,h
 real    :: r, radius = 1.0, rouratom(3), ratom(3)
 real    :: r1,r2
 character(len=2) :: Zlabel
 ! topology:
  do i=1,n_atoms
   k=0
   do j=1,n_atoms
    call checkatom(label(i),r,r1,h,Zlabel)
    call checkatom(label(j),r,r2,h,Zlabel)
    forall (h=1:3)
     rouratom(h)=xcrystal(h,i)
     ratom(h)=xcrystal(h,j)
    end forall
    call make_distances(ratom,rouratom, r)
    DistanceMatrix(i,j) = r
    DistanceMatrix(j,i) = r
    if(i/=j.and.r<=r1+r2+r_min_criteria_connectivity)then
     k=k+1
     ConnectedAtoms(i,j)=.true.
     ConnectedAtoms(j,i)=.true.
    end if
   end do
  end do
  RETURN
end subroutine to_measure_topoly
!
subroutine checkatom_type()
 implicit none
 integer :: i,j
 integer :: cc,nn,hh,zz
 write(6,'(a,1x,i5)')"Checking topology ...", k
 do i=1,n_atoms
  select case(label(i))
   case("C "," C")
    cc = 0 ; nn = 0 ; hh = 0 ; zz = 0
    do j=1,n_atoms
     if(i/=j.and.ConnectedAtoms(i,j))then
      select case(label(j))
       case("C "," C","C1","C2")
        cc=cc+1
       case("N "," N","N1")
        nn=nn+1
       case("H "," H")
        hh=hh+1
       case("Zn")
        zz=zz+1
      end select
     end if
    end do
    if(cc==0.and.nn==2.and.hh==1.and.zz==0)then
     label(i)="C4"
    else if ((cc==1.and.nn==1.and.hh==1.and.zz==0).or.&
             (cc==2.and.nn==1.and.hh==1.and.zz==0))then
     label(i)="C2"
    else if (cc==1.and.nn==2.and.hh==0.and.zz==0)then
     label(i)="C1"
    else if (cc==1.and.nn==0.and.hh==3.and.zz==0)then
     label(i)="C3"
    else
     write(6,*)'C=',cc,'N=',nn,'H=',hh,'Zn=',zz
     STOP "C-atom doesn't recognised"
    end if
   case("N "," N")
    label(i)="N1"
  end select
 end do
 return
end subroutine checkatom_type
!
subroutine compute_com(comc)
 implicit none
 integer                  :: i,j
 character(len=2)         :: ik
 real                     :: m,o,total_mass=0.0
 type(vector),intent(out) :: comc
 type(vector)             :: com
 com%x = 0.0
 com%y = 0.0
 com%z = 0.0
 do i=1,n_atoms
  call checkatom(label(i),m,o,j,ik) 
  total_mass=total_mass+m
  com%x = com%x + m*xcarte(1,i)
  com%y = com%y + m*xcarte(2,i)
  com%z = com%z + m*xcarte(3,i)
 end do
 com%x=com%x/total_mass
 com%y=com%y/total_mass
 com%z=com%z/total_mass
!
 comc%x = vr(1,1)*com%x + vr(1,2)*com%y + vr(1,3)*com%z + 0.5
 comc%y = vr(2,1)*com%x + vr(2,2)*com%y + vr(2,3)*com%z + 0.5
 comc%z = vr(3,1)*com%x + vr(3,2)*com%y + vr(3,3)*com%z + 0.5
!
 com%x = rv(1,1)*comc%x + rv(1,2)*comc%y + rv(1,3)*comc%z 
 com%y = rv(2,1)*comc%x + rv(2,2)*comc%y + rv(2,3)*comc%z 
 com%z = rv(3,1)*comc%x + rv(3,2)*comc%y + rv(3,3)*comc%z 
!
 line(1:80)=" "
 write(line( 1: 6),'(a6)')  adjustL("ATOM  ")
 write(line( 7:11),'(i5)')  1
 write(line(13:16),'(a)')   'Xe  '
 write(line(17:17),'(a)')   ' '
 write(line(18:20),'(a)')   adjustR('Xe')
 write(line(22:22),'(a)')   chain
 write(line(23:26),'(i4)')  ires
 write(line(27:27),'(a)')   chain
 write(line(31:38),'(f8.3)')com%x
 write(line(39:46),'(f8.3)')com%y
 write(line(47:54),'(f8.3)')com%z
 write(line(55:60),'(f6.2)')1.0
 write(line(61:66),'(f6.2)')0.0
 write(line(77:78),'(a2)') 'Xe' 
 write(line(79:80),'(a2)')  "+0"
 write(100,'(a80)')line(1:80)
 return
end subroutine compute_com
!
subroutine recentre(com)
 implicit none
 type(vector),intent(in) :: com
 do i=1,n_atoms
  ! crystallyne -> cartesians
  !xcrystal(1,i)=mod(xcrystal(1,i)-com%x+1.0,1.0)
  !xcrystal(2,i)=mod(xcrystal(2,i)-com%y+1.0,1.0)
  !xcrystal(3,i)=mod(xcrystal(3,i)-com%z+1.0,1.0)
  xcrystal(1,i)=xcrystal(1,i)-com%x
  xcrystal(2,i)=xcrystal(2,i)-com%y
  xcrystal(3,i)=xcrystal(3,i)-com%z
  forall ( j=1:3 )
   xcarte(j,i)=rv(j,1)*xcrystal(1,i) + rv(j,2)*xcrystal(2,i) + rv(j,3)*xcrystal(3,i)
  end forall
 end do
 return
end subroutine recentre
!
subroutine checkatom(Label,m,s,Z,Zlabel)
 implicit none
 character(len=2),intent(in)  :: Label
 real,intent(out)             :: m,s
 integer,intent(out)          :: Z
 character(len=2),intent(out) :: ZLabel
 select case(Label)
  ! 1  2 3 4 5 6 7  8  9  10 11 12
  ! 14 8 1 6 0 2 10 18 36 54 86 11
  case('N   ','N0  ':'N999')
   Z=7 
   m=14.0067-1.00794
   s=0.78
   Zlabel='N '
  case('C   ','C0  ':'C999')
   Z=6
   m=12.0107
   s=0.6
   ZLabel=' C'
  case('H   ','H0  ':'H999')
   Z=1
   m=1.00794
   s=0.620
   Zlabel=' H'
  case('Zn  ','Zn0 ':'Zn99')
   Z=30
   m=65.38-1.00794
   s=1.3
   Zlabel="Zn"
  case('He  ')
   Z=2
   m=1.00794
   s=0.0001
   Zlabel='He'
  case('Ne  ')
   Z=10
   m=1.00794
   s=0.0001
   Zlabel='Ne'
  case default
   write(6,'(a1,a4,a1)')"'",label,"'"
   STOP 'Atom unknowed in checkatom subroutine'
 end select
end subroutine checkatom
!
subroutine to_modify_vacancy(n_vac)
 implicit NONE
 integer             :: i, j, k, l, test
 integer,intent(out) :: n_vac
 integer             :: ii, jj, kk, ll, fou_atoms(1:4)
 real      :: x, y, z, r1(1:3),r2(1:3),r3(1:3)
 real      :: xcrystal_inT4(1:3,0:4), xinbox_inT4(1:3,0:4)
 real      :: vac2vacd = 5.20
 real      :: O2Od     = 4.35
 real      :: A(1:3,1:3),Ainv(1:3,1:3),B(1:3)
 type (vector)      :: x0,x1,x2,x3,c
 logical   :: singular_matrix_flag = .false.
 n_vac = 0
 search_vacancy_i: do i=1,n_atoms
  if((label(i)=="H ".or.label(i)==" H").and.printable(i).eqv..true.)then
   search_vacancy_j: do j=1,n_atoms
    if(i/=j.and.(label(j)=="H ".or.label(j)==" H").and.&
      DistanceMatrix(i,j)<=vac2vacd.and.printable(j).eqv..true.)then
      search_vacancy_k: do k=1,n_atoms
       if(k/=i.and.k/=j.and.(label(k)=="H ".or.label(k)==" H").and.&
         DistanceMatrix(i,k)<=vac2vacd.and.DistanceMatrix(j,k)<=vac2vacd.and.&
         printable(k).eqv..true.)then
         test=0
        search_vacancy_l: do l=1,n_atoms
        check1: if(l/=i.and.l/=j.and.l/=k.and.(label(l)=="H ".or.label(l)==" H").and.&
          DistanceMatrix(i,l)<=vac2vacd.and.DistanceMatrix(j,l)<=vac2vacd.and.&
          DistanceMatrix(k,l)<=vac2vacd.and.printable(l).eqv..true.)then
            ! Search connected Oxygen atoms for check2:
            search_iiO: do ii = 1,n_atoms
             if(ConnectedAtoms(i,ii).and.&
              (label(ii)=="O ".or.label(ii)==" O")) exit search_iiO
            end do search_iiO
            search_jjO: do jj = 1,n_atoms
             if(ConnectedAtoms(j,jj).and.&
              (label(jj)=="O ".or.label(jj)==" O")) exit search_jjO
            end do search_jjO
            search_kkO: do kk = 1,n_atoms
             if(ConnectedAtoms(k,kk).and.&
              (label(kk)=="O ".or.label(kk)==" O")) exit search_kkO
            end do search_kkO
            search_llO: do ll = 1,n_atoms
            if(ConnectedAtoms(l,ll).and.&
             (label(ll)=="O ".or.label(ll)==" O")) exit search_llO
            end do search_llO
            !
            check2: if (DistanceMatrix(ii,jj)<=O2Od.and.&
                        DistanceMatrix(ii,kk)<=O2Od.and.&
                        DistanceMatrix(ii,ll)<=O2Od.and.&
                        DistanceMatrix(jj,kk)<=O2Od.and.&
                        DistanceMatrix(jj,ll)<=O2Od.and.&
                        DistanceMatrix(kk,ll)<=O2Od.and.&
                        ii/=jj.and.ii/=kk.and.ii/=ll.and.jj/=kk.and.jj/=ll.and.kk/=ll)then
            ! All conditions:
            fou_atoms(1)=ii ; fou_atoms(2)=jj ; fou_atoms(3)=kk ; fou_atoms(4)=ll
            n_vac = n_vac + 1
            printable(i)=.false.
            printable(j)=.false.
            printable(k)=.false.
            printable(l)=.false.
            write(6,'(i5,1x,6(f14.7,1x),4(i4,1x),4(a2,1x))')n_vac,DistanceMatrix(i,j),&
              DistanceMatrix(i,k),&
              DistanceMatrix(i,l),&
              DistanceMatrix(j,k),&
              DistanceMatrix(j,l),&
              DistanceMatrix(k,l),i,j,k,l,label(i),label(j),label(k),label(l)
            write(6,'(i5,1x,6(f14.7,1x),4(i4,1x),4(a2,1x))')n_vac,&
              DistanceMatrix(ii,jj),&
              DistanceMatrix(ii,kk),&
              DistanceMatrix(ii,ll),&
              DistanceMatrix(jj,kk),&
              DistanceMatrix(jj,ll),&
              DistanceMatrix(kk,ll),ii,jj,kk,ll,label(ii),label(jj),label(kk),label(ll)
            ! Detect the O atoms:
            label(fou_atoms(1:4))= "Se"
           ! {{ coloco todos los atomos en el mismo sistema de referencia y calculo
           ! las distancias sin PBC }}
            elijo_pivote_4: forall ( kk=1:3 , jj=1:4 )
              xcrystal_inT4(kk,jj) = xcrystal(kk,fou_atoms(jj))
            end forall elijo_pivote_4
            do jj=2,4
             pivoteo_4: forall ( kk=1:3 )
              r1(kk) = xcrystal_inT4(kk,jj)   ! { atomo que cambia }
              r2(kk) = xcrystal_inT4(kk,jj-1) ! { atomo pivote }
             end forall pivoteo_4
             call minimum_image(r1,r2,r3,x)
             forall ( kk=1:3 )
              xcrystal_inT4(kk,jj) = r3(kk)
             end forall
            end do
           ! {{ paso a coordenadas cartesianas
            do ii=1,4
             forall ( jj=1:3 )
              xinbox_inT4(jj,ii)=rv(jj,1)*xcrystal_inT4(1,ii)+&
                                 rv(jj,2)*xcrystal_inT4(2,ii)+&
                                 rv(jj,3)*xcrystal_inT4(3,ii)
             end forall
            end do
            ! The circumcenter of a tetrahedron can be found as intersection
            ! of three bisector planes.
            x0%x=xinbox_inT4(1,1); x0%y=xinbox_inT4(2,1); x0%z=xinbox_inT4(3,1)
            x1%x=xinbox_inT4(1,2); x1%y=xinbox_inT4(2,2); x1%z=xinbox_inT4(3,2)
            x2%x=xinbox_inT4(1,3); x2%y=xinbox_inT4(2,3); x2%z=xinbox_inT4(3,3)
            x3%x=xinbox_inT4(1,4); x3%y=xinbox_inT4(2,4); x3%z=xinbox_inT4(3,4)
            ! B vector:
            B(1)=0.5*(absvec(x1)**2-absvec(x0)**2)
            B(2)=0.5*(absvec(x2)**2-absvec(x0)**2)
            B(3)=0.5*(absvec(x3)**2-absvec(x0)**2)
            ! A-Matrix:
            A(1,1)=x1%x - x0%x; A(2,1)=x1%y - x0%y; A(3,1)=x1%z - x0%z
            A(1,2)=x2%x - x0%x; A(2,2)=x2%y - x0%y; A(3,2)=x2%z - x0%z
            A(1,3)=x3%x - x0%x; A(2,3)=x3%y - x0%y; A(3,3)=x3%z - x0%z
            !
            write(6,'(a)')"B-vector:"
            do ii=1,3
              write(6,'(f14.7)')B(ii)
            end do
            write(6,'(a)')" A-Matrix:"
            do ii=1,3
             write(6,'(3(f14.7,1x))')( A(ii,jj),jj=1,3)
            end do
            ! Inverse: Based on Doolittle LU factorization for Ax=b
            call inverse(A,Ainv,3)
            write(6,'(a)')"Inverse A-Matrix (Doolittle LU factorization):"
            do ii=1,3
             write(6,'(3(f14.7,1x))')( Ainv(ii,jj),jj=1,3)
            end do
            ! check:
            singular_matrix_flag = .false.
            check_singular: do ii=1,3
              do jj=1,3
                if( Ainv(ii,jj) /= Ainv(ii,jj) ) then
                  singular_matrix_flag = .true.
                  exit check_singular
                end if
              end do
            end do check_singular
            if(singular_matrix_flag)then
             call matrixinv(A,Ainv,3)
              write(6,'(a)')'Singular Matrix!'
              write(6,'(a)')"Generalised Pseudo-Inverse A-Matrix (Based on Gauss-Jordan):"
              do ii=1,3
               write(6,'(3(f14.7,1x))')( Ainv(ii,jj),jj=1,3)
              end do
            end if
            ! Use circumcenter
            c%x = Ainv(1,1)*B(1)+Ainv(2,1)*B(2)+Ainv(3,1)*B(3)
            c%y = Ainv(1,2)*B(1)+Ainv(2,2)*B(2)+Ainv(3,2)*B(3)
            c%z = Ainv(1,3)*B(1)+Ainv(2,3)*B(2)+Ainv(3,3)*B(3)
            !
            write(6,'(a)')"Circumcenter:"
            forall ( jj=1:3 )
              xinbox_inT4(jj,0) =rv(jj,1)*c%x+&
                                 rv(jj,2)*c%y+&
                                 rv(jj,3)*c%z
            end forall
            write(6,'(3(f14.7,1x))')c%x,c%y,c%z
            write(6,'(a80)')"============================================================="
            !
           line(1:80)=" "
           write(line( 1: 6),'(a6)')   adjustL("ATOM  ")
           write(line( 7:11),'(i5)')   n_vac
           write(line(13:16),'(a)')    adjustR("Ge")
           write(line(17:17),'(a)')    chain
           write(line(18:20),'(a)')    adjustR("Ge")
           write(line(22:22),'(a)')    chain
           write(line(23:26),'(i4)')   ires
           write(line(27:27),'(a)')    chain
           write(line(31:38),'(f8.3)') c%x
           write(line(39:46),'(f8.3)') c%y
           write(line(47:54),'(f8.3)') c%z
           write(line(55:60),'(f6.2)') 1.0
           write(line(61:66),'(f6.2)') 0.0
           write(line(77:78),'(a2)')   adjustR("Ge")
           write(line(79:80),'(a2)')   "+0"
           write(100,'(a80)')line(1:80)
         else
           test=test+1
           printable(i)=.true.
           printable(j)=.true.
           printable(k)=.true.
           printable(l)=.true.
           label(i)="H "
           label(j)="H "
           label(k)="H "
           label(l)="H "
           label(ii)="O "
           label(jj)="O "
           label(kk)="O "
           label(ll)="O "
           write(6,*)"trying:",test,i,j,k,l,ii,jj,kk,ll
         end if check2
        end if check1
        end do search_vacancy_l
       end if
      end do search_vacancy_k
    end if
  end do search_vacancy_j
  end if
 end do search_vacancy_i
 return
end subroutine to_modify_vacancy
!
subroutine remove_shells_()
  implicit none
  integer :: i,j,k,l
  do i=1,n_atoms
   if(label(i)=="O ".or.label(i)==" O")then
   searching_shell: do j=1,n_atoms
     if( i/=j .and. DistanceMatrix(i,j)<=0.5 .and. &
     (label(j)=="Xe".or.label(j)==" H".or.label(j)=="H ")) then
     printable(j)=.false.
     write(6,'(a,1x,f14.7)')"[Removing shell!]",DistanceMatrix(i,j)
     exit searching_shell
    !else
     !write(6,'(a,1x,f14.7)')"[dev]", DistanceMatrix(i,j)
    end if
   end do searching_shell
  end if
  end do
  return
end subroutine remove_shells_

!
SUBROUTINE make_distances(r2,r1,dist)
IMPLICIT NONE
REAL,    intent(in)  :: r1(3),r2(3)                       ! coordenadas y matriz de cambio
REAL,    intent(out) :: dist
REAL                 :: d_image(1:27),image(3,27)        ! array de distancias
INTEGER              :: k,l,m,n,o,i,j                    ! variables mudas
REAL                 :: atom(3),ouratom(3)               ! coordenadas preparadas
 k=0
 do l=-1,1
  do m=-1,1
     do n=-1,1
        k = k + 1
        ouratom(1) = r1(1)
        ouratom(2) = r1(2)
        ouratom(3) = r1(3)
        atom(1) = r2(1) + l
        atom(2) = r2(2) + m
        atom(3) = r2(3) + n
        d_image(k) = distance(atom,ouratom)
        forall ( i=1:3)
          image(i,k) = atom(i)
        end forall
        !write(6,*)k, d_image(k)
    enddo
  enddo
 enddo
 dist=MINVAL(d_image)
 RETURN
END SUBROUTINE
!
SUBROUTINE minimum_image(r2,r1,r3,dist)
! {{ prepara para el calculo de distancias en una celda triclinica }}
 IMPLICIT NONE
 REAL,    intent(in)  :: r1(3),r2(3)     ! coordenadas y matriz de cambio
 REAL,    intent(out) :: dist,r3(1:3)    ! matriz de distancias N X N
 REAL                 :: d_image(1:27),image(3,27)        ! array de distancias
 REAL                 :: phi
 INTEGER              :: k,l,m,n,o,i,j                    ! variables mudas
 REAL                 :: atom(3),ouratom(3)               ! coordenadas preparadas
  k=0
  do l=-1,1
   do m=-1,1
      do n=-1,1
         k = k + 1
         ouratom(1) = r1(1)
         ouratom(2) = r1(2)
         ouratom(3) = r1(3)
         atom(1) = r2(1) + l
         atom(2) = r2(2) + m
         atom(3) = r2(3) + n
         d_image(k) = distance(atom,ouratom)
         forall ( i=1:3)
           image(i,k) = atom(i)
         end forall
     enddo
   enddo
  enddo
  dist=MINVAL(d_image)
  phi=1000.0
  k=1
  do l=1,27
   if(d_image(l)<=phi)then
      phi=d_image(l) ! seleccionamos el parametro menor
      k=l            ! y el contador correspondiente.
   endif
  enddo
  forall ( l=1:3)
   r3(l)=image(l,k)
  end forall
 RETURN
END SUBROUTINE minimum_image
!
REAL FUNCTION DISTANCE(atom,ouratom)
 IMPLICIT NONE
 INTEGER :: j
 REAL :: atom(3),ouratom(3),per_atom(3),dist(3),o_atom(3),o_ouratom(3)
 FORALL ( j=1:3 )
  o_ouratom(j) = rv(j,1)*ouratom(1)  + rv(j,2)*ouratom(2)  + rv(j,3)*ouratom(3)
  o_atom(j)    = rv(j,1)*atom(1) + rv(j,2)*atom(2) + rv(j,3)*atom(3)
  dist(j) = o_ouratom(j) - o_atom(j)
 END FORALL
 DISTANCE = sqrt(dist(1)*dist(1) + dist(2)*dist(2) + dist(3)*dist(3))
END FUNCTION
!
SUBROUTINE cell(rv,vr,cell_0)
 implicit none
 integer :: i,j
 real, intent(in)  :: cell_0(6)
 real, intent(out) :: rv(3,3),vr(3,3)
 real :: alp,bet
 real :: cosa,cosb,cosg
 real :: gam,sing
 real :: pi,DEGTORAD
 pi = ACOS(-1.0)
 DEGTORAD=pi/180.0
 IF(cell_0(4) == 90.0) THEN
   cosa = 0.0
 ELSE
   ALP=cell_0(4)*degtorad
   COSA=cos(ALP)
 ENDIF
 IF(cell_0(5) == 90.0) THEN
   cosb = 0.0
 ELSE
   bet = cell_0(5)*degtorad
   cosb = cos(bet)
 ENDIF
 IF(cell_0(6) == 90.0) then
   sing = 1.0
   cosg = 0.0
 ELSE
   gam = cell_0(6)*degtorad
   sing = sin(gam)
   cosg = cos(gam)
 ENDIF
 rv(1,1) = cell_0(1)
 rv(1,2) = cell_0(2)*cosg
 rv(1,3) = cell_0(3)*cosb
 rv(2,1) = 0.0
 rv(2,2) = cell_0(2)*sing
 rv(2,3) = cell_0(3)*(cosa - cosb*cosg)/sing
 rv(3,1) = 0.0
 rv(3,2) = 0.0
 rv(3,3) = sqrt( cell_0(3)*cell_0(3) - rv(1,3)*rv(1,3) - rv(2,3)*rv(2,3))
 call inverse(rv,vr,3)
! print*,'Cell:'
! WRITE(*,'(6F14.7)')( cell_0(j), j=1,6 )
! print*,'Box:'
! DO i=1,3
!    WRITE(*,'(F14.7,F14.7,F14.7)')( rv(i,j), j=1,3 )
! ENDDO
! WRITE(*,*)'----------------------------------------'
! WRITE(*,*)'bOX:'
! DO i=1,3
!    WRITE(*,'(F14.7,F14.7,F14.7)')( vr(i,j), j=1,3 )
! ENDDO
 RETURN
END SUBROUTINE cell
!
SUBROUTINE uncell(rv,cell_0)
  implicit none
  real,intent(out) :: cell_0(6)
  real,intent(in)  :: rv(3,3)
  integer  :: i,j
  real     :: temp(6)
  REAL :: radtodeg,PI
  PI=ACOS(-1.0)
  radtodeg=180.0/PI
!
  do i = 1,3
    temp(i) = 0.0
    do j = 1,3
      temp(i) = temp(i) + rv(j,i)**2
    enddo
    temp(i) = sqrt(temp(i))
  enddo
  cell_0(1) = abs(temp(1))
  cell_0(2) = abs(temp(2))
  cell_0(3) = abs(temp(3))
  do i = 1,3
    temp(3+i) = 0.0
  enddo
  do j = 1,3
    temp(4) = temp(4) + rv(j,2)*rv(j,3)
    temp(5) = temp(5) + rv(j,1)*rv(j,3)
    temp(6) = temp(6) + rv(j,1)*rv(j,2)
  enddo
  temp(4) = temp(4)/(temp(2)*temp(3))
  temp(5) = temp(5)/(temp(1)*temp(3))
  temp(6) = temp(6)/(temp(1)*temp(2))
  cell_0(4) = radtodeg*acos(temp(4))
  cell_0(5) = radtodeg*acos(temp(5))
  cell_0(6) = radtodeg*acos(temp(6))
!  Avoid round off errors for 90.0 and 120.0 degrees
  DO i=4,6
     if (abs(cell_0(i) - 90.0 ).lt.0.00001) cell_0(i) = 90.0
     if (abs(cell_0(i) - 120.0).lt.0.00001) cell_0(i) = 120.0
  ENDDO
!
  return
end subroutine uncell
!
real function det_3x3(a)
 !  Result
 !  Input
 real, intent(in) :: a(1:3,1:3)
 det_3x3 =  a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2) &
          -a(1,3)*a(2,2)*a(3,1)-a(1,2)*a(2,1)*a(3,3)-a(1,1)*a(2,3)*a(3,2)
end function
!=======================================================================
subroutine adj(a,b)
  implicit none
!  Adjugate of a 3x3 matrix (the transpose of the cofactor matrix).
!  Result
   real,intent(out) :: b(1:3,1:3)
   real, intent(in) :: a(1:3,1:3)
!
  b(1,1) = a(2,2)*a(3,3) - a(2,3)*a(3,2)
  b(1,2) = a(1,3)*a(3,2) - a(1,2)*a(3,3)
  b(1,3) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
  b(2,1) = a(2,3)*a(3,1) - a(2,1)*a(3,3)
  b(2,2) = a(1,1)*a(3,3) - a(1,3)*a(3,1)
  b(2,3) = a(1,3)*a(2,1) - a(1,1)*a(2,3)
  b(3,1) = a(2,1)*a(3,2) - a(2,2)*a(3,1)
  b(3,2) = a(1,2)*a(3,1) - a(1,1)*a(3,2)
  b(3,3) = a(1,1)*a(2,2) - a(1,2)*a(2,1)
  return
end subroutine adj
!=======================================================================
subroutine inv(a,b)
!  An inverse of a 3x3 matrix.
!  Result
  !real(dp), dimension(3,3) :: inv
!  Input
   real, intent(in) :: a(1:3,1:3)
   real             :: c(1:3,1:3)
   real,intent(out) :: b(1:3,1:3)
   integer :: i,j
!  Locals
   real :: detr
   detr = 1.0/det_3x3(a)
   call adj(a,c)
   do i=1,3
    do j=1,3
     b(i,j) = c(i,j)*detr
   end do
 end do
 return
end subroutine inv
!
subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed
! during the calculation
!===========================================================
 implicit none
 integer n
 real a(n,n), c(n,n)
 real L(n,n), U(n,n), b(n), d(n), x(n)
 real coeff
 integer i, j, k
! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
 L=0.0
 U=0.0
 b=0.0
! step 1: forward elimination
 do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
 end do
! Step 2: prepare L and U matrices
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
 do i=1,n
  L(i,i) = 1.0
 end do
! U matrix is the upper triangular part of A
 do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
 end do
!
! Step 3: compute columns of the inverse matrix C
 do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
 end do
 return
END SUBROUTINE inverse
!
 subroutine matrixinv(a,b,n)
   ! subroutine to calculate the inverse of a matrix
   ! using Gauss-Jordan elimination
   implicit none
   integer                :: i,j,k,l,m,irrow
   integer, intent(in)    :: n
   real,    intent(inout) :: a(1:n,1:n)
   real,    intent(out)   :: b(1:n,1:n)
   real                   :: big,dum
   !
   do i =1,n
    do j=1,n
      b(i,j)=0.0
    end do
    b(i,i)=1.0
   end do
   ! this is the big loop over all the columns of a(n,n)
   ! in case the entry a(i,i) is zero, we need to find a good pivot; this pivot
   ! is chosen as the largest value on the column i from a(j,i) with j = 1,n
   do i=1,n
     big=a(i,i)
     do j=i,n
      if( a(j,i) > big ) then
        big = a(j,i)
        irrow = j
      end if
     end do
     ! interchange lines i with irow for both a() and b() matrices
     if( big > a(i,i)) then
       do k = 1,n
         dum = a(i,k)
         a(i,k) = a(irrow,k)
         a(irrow,k) = dum
         dum = b(i,k)
         b(i,k) = b(irrow,k)
         b(irrow,k)= dum
       end do
     end if
     ! divide all entries in line i from a(i,j) by the value a(i,i);
     ! same operation for the identity matrix
     dum = a(i,i)
     do j = 1,n
       a(i,j) = a(i,j)/dum
       b(i,j) = b(i,j)/dum
     end do
     ! make zero all entries in the column a(j,i); same operation for indent()
     do j = i + 1,n
       dum = a(j,i)
       do k = 1,n
         a(j,k) = a(j,k) - dum*a(i,k)
         b(j,k) = b(j,k) - dum*b(i,k)
       end do
     end do
   end do
   ! Finally... substract appropiate multiple of row j from row j-1
   do i=1,n-1
     do j = i+1,n
       dum = a(i,j)
       do l=1,n
        a(i,l) = a(i,l) - dum*a(j,l)
        b(i,l) = b(i,l) - dum*b(j,l)
       end do
     end do
   end do
   return
 end subroutine matrixinv
!
real function volume(rv)
  implicit none
  real, intent(in)  :: rv(3,3)
  real       :: r1x
  real       :: r1y
  real       :: r1z
  real       :: r2x
  real       :: r2y
  real       :: r2z
  real       :: r3x
  real       :: r3y
  real       :: r3z
  real       :: vol
!
  r1x = rv(1,1)
  r1y = rv(2,1)
  r1z = rv(3,1)
  r2x = rv(1,2)
  r2y = rv(2,2)
  r2z = rv(3,2)
  r3x = rv(1,3)
  r3y = rv(2,3)
  r3z = rv(3,3)
  vol = r1x*(r2y*r3z - r2z*r3y) + r1y*(r3x*r2z - r3z*r2x) + r1z*(r2x*r3y - r2y*r3x)
  volume = abs(vol)
  RETURN
end function
end program lammpstrj2pdb
