module vector_module
implicit none
type  :: vector
 sequence
 real :: x
 real :: y
 real :: z
end type vector
contains
!
 pure type(vector) function array2vector(a)
  implicit none
  real,intent(in)         :: a(1:3)
  array2vector%x  =  a(1) 
  array2vector%y  =  a(2) 
  array2vector%z  =  a(3) 
 end function array2vector
!
 pure function vector2array(v) result(a)
  implicit none
  type(vector),intent(in) :: v
  real :: a(1:3)
  a(1) = v%x
  a(2) = v%y
  a(3) = v%z
 end function vector2array
!
 pure type(vector) function PointCoordinatesAtDistanceFromOriginAndDirection(o,u,d)
  implicit none
  type(vector),intent(in) :: o,u
  real,intent(in)         :: d
  PointCoordinatesAtDistanceFromOriginAndDirection%x = o%x + u%x*d
  PointCoordinatesAtDistanceFromOriginAndDirection%y = o%y + u%y*d
  PointCoordinatesAtDistanceFromOriginAndDirection%z = o%z + u%z*d
 end function PointCoordinatesAtDistanceFromOriginAndDirection
!
 pure type (vector) function vector_unitary(v1)
  implicit none
  type(vector), intent(in) :: v1
  real                     :: norma
  norma=absvec(v1)
  vector_unitary%x = v1%x/norma
  vector_unitary%y = v1%y/norma
  vector_unitary%z = v1%z/norma
 end function vector_unitary
!
 pure type (vector) function vector_scale(v1, r )
  implicit none
  type(vector), intent(in) :: v1
  real,intent(in)          :: r
  vector_scale%x = v1%x*r
  vector_scale%y = v1%y*r
  vector_scale%z = v1%z*r
 end function vector_scale
!
 pure type (vector) function vector_add (v1,v2)
  implicit none
  type (vector), intent(in) :: v1,v2
  vector_add%x = v1%x + v2%x
  vector_add%y = v1%y + v2%y
  vector_add%z = v1%z + v2%z
 end function vector_add
!
 pure type (vector) function vector_sub (v1,v2)
  implicit none
  type (vector), intent(in) :: v1,v2
  vector_sub%x = v1%x - v2%x
  vector_sub%y = v1%y - v2%y
  vector_sub%z = v1%z - v2%z
 end function vector_sub
!
 pure type(vector) function cross(a,b)
  implicit none
  type(vector),intent(in) :: a, b
  cross%x = a%y * b%z - a%z * b%y
  cross%y = a%z * b%x - a%x * b%z
  cross%z = a%x * b%y - a%y * b%x
 end function cross
!
 pure real function dot(a,b)
  implicit none
  type(vector), intent(in) :: a,b
  dot = a%x * b%x + a%y*b%y + a%z*b%z
 end function dot
!
 pure real function absvec(a)
  type (vector), intent (in) :: a
  absvec = sqrt(a%x**2 + a%y**2 + a%z**2)
 end function absvec
!
 pure real function dist2vectors(a,b)
  type(vector), intent(in) :: a,b
  dist2vectors = absvec( vector_sub(b,a) )
 end function  dist2vectors
!
 pure real function angle2vectors(a,b)
  ! vectors are the  r_ij vector
  !                 (r_ij, r_jk)
  type(vector), intent(in) :: a,b
  angle2vectors = acos( dot(a,b)/(absvec(a)*absvec(b)) )
 end function  angle2vectors
!
 pure real function angle3vectors_jik(r_i,r_j,r_k)
  type(vector), intent(in) :: r_i,r_j,r_k
  type(vector)             :: r_ij,r_ik,r_jk
  r_ij = vector_sub(r_i,r_j)
  r_jk = vector_sub(r_j,r_k)
  r_ik = vector_sub(r_i,r_k)
  angle3vectors_jik = acos(dot(r_ij,r_ik)/(absvec(r_ij)*absvec(r_ik)))
 end function  angle3vectors_jik
! 
 pure real function dihedral_angle4vectors_ijkl(r_i,r_j,r_k,r_l)
  type(vector), intent(in) :: r_i,r_j,r_k,r_l
  type(vector)             :: r_ij,r_jk,r_kl
  r_ij = vector_sub(r_i,r_j)
  r_jk = vector_sub(r_j,r_k)
  r_kl = vector_sub(r_k,r_l)
  dihedral_angle4vectors_ijkl = acos(dot(cross(r_ij,r_jk),cross(r_jk,r_kl))/(absvec(cross(r_ij,r_jk))*absvec(cross(r_jk,r_kl))))
 end function  dihedral_angle4vectors_ijkl
end module vector_module

program pdb2cif
 use vector_module
 implicit none
 type(vector), allocatable      :: xcryst(:),xinbox(:),xcryst_all(:)
 type(vector)                   :: ov,uv,v1,v2
 character(len=4), allocatable  :: label(:),label_all(:)
 integer                        :: n_atoms = 0
 integer                        :: n_t = 0
 integer                        :: ierr,i,j,k
 real, parameter                :: pi=acos(-1.0)
 real                           :: atom(3),ouratom(3),cell_0(1:6),r1(3),r2(3),r3(3)
 real                           :: rv(1:3,1:3),vr(1:3,1:3),average(1:3)
 real                           :: s, m
 integer                        :: Z
 character(len=2)               :: Zlabel 
 character(len=80)              :: line,string
 character(len=4)               :: mol
 character(len=6)               :: atomc
 character(len=5)               :: model
 character(len=6)               :: charct(6)
 character(len=4)               :: typ(3)
 character(len=10)              :: spacegroup = "p1"
! ...
 fileopen_pdb: do
  read (5,'(a)',iostat=ierr) line
  if( ierr /= 0 ) exit fileopen_pdb
  if (line(1:5)=='MODEL') then
    n_atoms = 0
  endif
  if (line(1:4)=='ATOM') n_atoms = n_atoms + 1
 end do fileopen_pdb
 allocate(xcryst(n_atoms),stat=ierr)
 allocate(xinbox(n_atoms),stat=ierr)
 allocate(label(n_atoms) ,stat=ierr)
 if( ierr /= 0 ) stop "problema al alicatar variables"
 rewind(5)
 do
  read (5,'(a)') line
  if(line(1:4)=="CRYS") exit
 end do
 read (line,'(6x,3f9.3,3f7.2,1x,a10)') &
 cell_0(1),cell_0(2),cell_0(3),cell_0(4),cell_0(5),cell_0(6),spacegroup
 spacegroup = "p1"
 call cell(rv,vr,cell_0)
 read_coor_pdb: do i=1,n_atoms
    read(5, '(a)', iostat = ierr ) line
    if (ierr/=0) exit read_coor_pdb          ! formatos de lectura para pdb 
       atomc = line(1:6)                     ! leemos
       mol   = line(18:20)
       read(line(31:38),*) atom(1)
       read(line(39:46),*) atom(2)
       read(line(47:54),*) atom(3)
       xinbox(i) = array2vector(atom)
       label(i)= line(13:14)
    call checkatom(label(i),m,s,Z,Zlabel)
    label(i)=Zlabel
    xcryst(i) = array2vector(Box2CrystalCoordinates(vector2array(xinbox(i))))
    if( label(i)(1:2)==" N" ) n_t = n_t + 2
 enddo read_coor_pdb
 read (5,'(a)') line
! add dummy atoms:
 allocate(xcryst_all(n_atoms+n_t),stat=ierr)
 allocate(label_all(n_atoms+n_t),stat=ierr)
 n_t = 0
 do i=1,n_atoms
  xcryst_all(i)=xcryst(i)
  label_all(i)=label(i)
  if(label(i)(1:2)=="Zn")then
   do j=1,n_atoms
    if(label(j)(1:2)==" N")then
    ! put the both atoms in the same image
     call make_distances( .true.,cell_0,vector2array(xcryst(j)),&
                          vector2array(xcryst(i)),rv,atom,s     )
     if (s <= 3.0) then
      write(6,'(i4,1x,i4,1x,a4,1x,a4,1x,a,1x,f14.7)')i,j,label(i),label(j),':',s
      v1 = xinbox(i)
      v2 = array2vector(Crystal2BoxCoordinates(atom))
      ! add He:
      n_t = n_t + 1
      ov = PointCoordinatesAtDistanceFromOriginAndDirection(v1,vector_unitary(vector_sub(v2,v1)),0.9)
      xcryst_all(n_atoms+n_t)=array2vector(Box2CrystalCoordinates(vector2array(ov)))
      label_all(n_atoms+n_t)="He  "
      ! add Ne:
      n_t = n_t + 1
      uv = PointCoordinatesAtDistanceFromOriginAndDirection(v1,vector_unitary(vector_sub(v2,v1)),s-0.5)
      xcryst_all(n_atoms+n_t)=array2vector(Box2CrystalCoordinates(vector2array(uv)))
      label_all(n_atoms+n_t)="Ne  "
      write(6,'(a,1x,3(f14.7,1x),a)')'Zn coordinates: (',v1%x,v1%y,v1%z,')'
      write(6,'(a,1x,3(f14.7,1x),a)')'N  coordinates: (',v2%x,v2%y,v2%z,')'
      write(6,'(a,1x,3(f14.7,1x),a)')'He coordinates: (',ov%x,ov%y,ov%z,')'
      write(6,'(a,1x,3(f14.7,1x),a)')'Ne coordinates: (',uv%x,uv%y,uv%z,')' 
      !check:
      !r2 = Box2CrystalCoordinates(vector2array(uv))
      !r1 = vector2array(xcryst(i)) 
      !call make_distances(.false.,cell_0,r2,r1,rv,r3,s)
      !write(6,'(i4,1x,i4,1x,a4,1x,a4,1x,a,1x,f14.7)')i,0,label(i),'He',':',s
     end if
    end if
   end do
  end if
 end do
!
 call escritura_cif()
 close(5)
 contains
!
 pure function Crystal2BoxCoordinates(r_c) result (r_b)
  implicit none
  real,intent(in)  :: r_c(1:3)
  real             :: r_b(1:3)
  integer          :: j
  forall ( j=1:3 )
   r_b(j) = rv(j,1)*r_c(1)  + rv(j,2)*r_c(2)  + rv(j,3)*r_c(3)
  end forall
 end function Crystal2BoxCoordinates
!
 pure function Box2CrystalCoordinates(r_b) result (r_c)
  implicit none
  real,intent(in)  :: r_b(1:3)
  real             :: r_c(1:3)
  integer          :: j
  forall ( j=1:3 )
   r_c(j) = mod(vr(j,1)*r_b(1)  + vr(j,2)*r_b(2)  + vr(j,3)*r_b(3) + 100.0,1.0)
  end forall
 end function Box2CrystalCoordinates
!
 subroutine escritura_cif()
  implicit none
  integer   :: u,i
  u=1000
  open(u,file="p1_fromPDBFile.cif")
  write(u,'(a)')'data_subtitutions'
  write(u,'(a)')'_audit_creation_method    igor'
  write(u,'(a)')"_audit_author_name 'sponge bob'"
  write(u,'(a,f14.7)')'_cell_length_a',cell_0(1)
  write(u,'(a,f14.7)')'_cell_length_b',cell_0(2)
  write(u,'(a,f14.7)')'_cell_length_c',cell_0(3)
  write(u,'(a,f14.7)')'_cell_angle_alpha',cell_0(4)
  write(u,'(a,f14.7)')'_cell_angle_beta',cell_0(5)
  write(u,'(a,f14.7)')'_cell_angle_gamma',cell_0(6)
  write(u,'(a,f14.7)')'_cell_volume',volume(rv)
  write(u,'(a)')"_symmetry_space_group_name_hall 'p 1'"
  write(u,'(a)')"_symmetry_space_group_name_h-m 'p 1'"
  write(u,'(a)')'_symmetry_int_tables_number 1'
  write(u,'(a)')"_symmetry_equiv_pos_as_xyz 'x,y,z'"
  write(u,'(a)')'loop_'
  write(u,'(a)')'_atom_site_label'
  write(u,'(a)')'_atom_site_fract_x'
  write(u,'(a)')'_atom_site_fract_y'
  write(u,'(a)')'_atom_site_fract_z'
 !write(u,'(a)')'_atom_site_charge' 
  atoms_: do i=1,n_atoms+n_t
   write(u,'(a4,1x,4(f10.8,1x))')label_all(i),xcryst_all(i)%x,xcryst_all(i)%y,xcryst_all(i)%z
  end do atoms_
  close(u)
  return
 end subroutine escritura_cif
!
subroutine cell(rv,vr,cell_0)
 implicit none
 integer :: i,j
 real, intent(in)  :: cell_0(6)
 real, intent(out) :: rv(3,3),vr(3,3)
 real :: alp,bet
 real :: cosa,cosb,cosg
 real :: gam,sing
 real :: pi,degtorad
 pi = acos(-1.0)
 degtorad=pi/180.0
 if(cell_0(4) == 90.0) then
   cosa = 0.0
 else
   alp=cell_0(4)*degtorad
   cosa=cos(alp)
 endif
 if(cell_0(5) == 90.0) then
   cosb = 0.0
 else
   bet = cell_0(5)*degtorad
   cosb = cos(bet)
 endif
 if(cell_0(6) == 90.0) then
   sing = 1.0
   cosg = 0.0
 else
   gam = cell_0(6)*degtorad
   sing = sin(gam)
   cosg = cos(gam)
 endif
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
 return
end subroutine cell
!
subroutine uncell(rv,cell_0)
  implicit none
  real,intent(out) :: cell_0(6)
  real,intent(in)  :: rv(3,3)
  integer  :: i,j
  real     :: temp(6)
  real :: radtodeg,pi
  pi=acos(-1.0)
  radtodeg=180.0/pi
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
  do i=4,6
     if (abs(cell_0(i) - 90.0 ).lt.0.00001) cell_0(i) = 90.0
     if (abs(cell_0(i) - 120.0).lt.0.00001) cell_0(i) = 120.0
  end do
  return
end subroutine uncell
!
subroutine inverse(a,c,n)
 implicit none 
 integer n
 real a(n,n), c(n,n)
 real l(n,n), u(n,n), b(n), d(n), x(n)
 real coeff
 integer i, j, k
 l=0.0
 u=0.0
 b=0.0
 do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      l(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
 end do
 do i=1,n
  l(i,i) = 1.0
 end do
 do j=1,n
  do i=1,j
    u(i,j) = a(i,j)
  end do
 end do
 do k=1,n
  b(k)=1.0
  d(1) = b(1)
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - l(i,j)*d(j)
    end do
  end do
  x(n)=d(n)/u(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-u(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
 end do
 return
end subroutine inverse
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
  return
end function
!
subroutine make_distances(flag,cell_0,r2,r1,rv,r3,dist)
 use vector_module, only: dist2vectors, array2vector 
 implicit none
 real,    intent(in)  :: r1(3),r2(3),rv(3,3),cell_0(6)    ! coordenadas y matriz de cambio
 real,    intent(out) :: dist,r3(1:3)                     ! matriz de distancias n x n
 real                 :: d_image(1:27),image(3,27)        ! array de distancias
 real                 :: distance,rcm(3),phi
 integer              :: k,l,m,n,o,i,j                    ! variables mudas
 real                 :: atom(3),ouratom(3)               ! coordenadas preparadas
 logical              :: flag                             ! out the coordinate of the atom 
! {{ calculamos la matriz de distancias }}
  k=0
  do l=-1,1
   do m=-1,1
      do n=-1,1
         k = k + 1
         ouratom = r1
         atom(1) = r2(1) + real(l)
         atom(2) = r2(2) + real(m)
         atom(3) = r2(3) + real(n)
         d_image(k) = dist2vectors( array2vector(Crystal2BoxCoordinates(ouratom)),&
                                    array2vector(Crystal2BoxCoordinates(atom   )) )
         !d_image(k) = distance(atom,ouratom,rv)
         forall ( i=1:3)       ! r3 sera la imagen de atom a la menor distancia
          image(i,k) = atom(i) ! entre todas las imagenes.
         end forall
     enddo
   enddo
  enddo
  dist=minval(d_image)
  if(flag)then
   phi=9999999.999999 ! initial infinite distance
   k=1                ! selected image
   do l=1,27
    if(d_image(l)<=phi)then
      phi=d_image(l) ! seleccionamos la imagen con la menor distancia
      k=l            ! 
    endif
   enddo
   forall ( l=1:3)
     r3(l)=image(l,k)
   end forall
  else
   r3(1:3)=0.0
  end if
 return
end subroutine make_distances
 !
real function distance(atom,ouratom,rv)
 implicit none
 integer :: j
 real :: atom(3),ouratom(3),per_atom(3),dist(3),o_atom(3),o_ouratom(3)
 real :: rv(3,3)
 forall ( j=1:3 )
  o_ouratom(j) = rv(j,1)*ouratom(1)  + rv(j,2)*ouratom(2)  + rv(j,3)*ouratom(3)
  o_atom(j)    = rv(j,1)*atom(1) + rv(j,2)*atom(2) + rv(j,3)*atom(3)
  dist(j) = o_ouratom(j) - o_atom(j)
 end forall
 distance = sqrt(dist(1)*dist(1) + dist(2)*dist(2) + dist(3)*dist(3))
end function distance
!
 pure subroutine checkatom(Label,m,s,Z,Zlabel)
  implicit none
  character(len=4),intent(in)  :: Label
  real,intent(out)             :: m,s
  integer,intent(out)          :: Z
  character(len=2),intent(out) :: ZLabel
  real,parameter               :: dummy_mass = 1.00794 ! (hydrogen atom mass)
  select case(Label)
   case('O   ','O0  ':'O999')
    Z=8
    m=15.999  ! conventional weight
    s=0.70    ! covalent radious
    Zlabel=' O'
   case('C   ','C0  ':'C999')
    Z=6
    m=12.0107
    s=0.720
    ZLabel=' C'
   case('H   ','H0  ':'H999')
    Z=1
    m=dummy_mass
    s=0.40
    Zlabel=' H'
   case('N   ','N0  ':'N999')
    Z = 7
    m = 14.0067 - dummy_mass
    s = 0.7
    Zlabel=' N'
   case('Zn  ','Zn0 ':'Zn99')
    Z = 30
    m = 65.38 - dummy_mass
    s = 0.6
    Zlabel='Zn'
   case('Cl  ',' Cl ','Cl0 ':'Cl99')
    Z = 17
    m = 35.453
    s = 0.5
    Zlabel = 'Cl'
! dummy atoms:
   case('He  ','He0 ':'He99')
    Z = 2
    m = dummy_mass
    s = (1.470 + 0.1)/2.0
    Zlabel = 'He'
   case('Ne  ','Ne0 ':'Ne99')
    Z = 10
    m = dummy_mass
    s = 0.10
    Zlabel='Ne'
  end select
 end subroutine checkatom
end program pdb2cif
