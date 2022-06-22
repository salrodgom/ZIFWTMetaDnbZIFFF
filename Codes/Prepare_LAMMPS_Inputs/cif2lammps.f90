program zif_cif2gin
 use iso_fortran_env
 implicit none
! locals
 integer             :: i,j,k,l,h,hh,m,n,ierr
 integer             :: ii,jj,kk,ll
 real                :: r = 1.0e12
!parameters
 real,parameter      :: k_B = 8.617332478e-5
 real,parameter      :: r_min_criteria_connectivity=0.15
 INTEGER, PARAMETER  :: muchisimo=100000
! variables
 integer             :: num_args
 integer             :: n_atoms = 0
 real                :: r0=r_min_criteria_connectivity
 real                :: ratom(3),rouratom(3)
 real                :: rv(3,3),vr(3,3),cell_0(1:6) = 0.0
 real                :: xlo_bound,ylo_bound,zlo_bound,xhi_bound,yhi_bound,zhi_bound
 real                :: xy,xz,yz
 real                :: total_Q = 0.0
 character(len=8)    :: spams
 character(len=100)  :: line
 character(len=20)   :: spam,simulation_type="single"
 character(len=80)   :: string_stop_head= "_atom_site_charge" !"_atom_site_charge"
 character(len=100)  :: CIFFilename=" "
 character(len=100)  :: filename=" "
 character(len=4)    :: sc(1:4)
 CHARACTER (LEN=18)  :: string(muchisimo)
 integer,parameter   :: n_atom_types_max=100
 integer             :: n_atom_types = 0
 integer,allocatable :: check_topol_id(:)
 real,allocatable    :: check_topol_distance(:)
 integer,allocatable :: check_topol_degree(:), table_atoms_atom_types(:,:)
 character(len=4)    :: atom_types(n_atom_types_max)
 integer             :: atom_types_histogram(n_atom_types_max)
 character(len=10)   :: atom_types_hybridization(n_atom_types_max)
 ! type
 integer,parameter   :: max_n_componets=2
 type  :: node
  integer           :: element
  integer           :: type_
  character(len=2)  :: label_element
  character(len=4)  :: label
  character(len=4)  :: new_label="Xxxx"
  character(len=6)  :: label_from_CIFFile="Xxxxxx"
  character(len=50) :: topol_label="Xxxx"
  character(len=10) :: hybridization="Unknown"
  integer          :: degree
  real             :: charge
  real             :: radius
  real             :: mass
  integer          :: n_components
  character(5)     :: clabel(max_n_componets)
  real             :: xyzs(1:3,max_n_componets)
  real             :: xyzc(1:3,max_n_componets)
 end type
 ! allocatable
 type(node),allocatable,dimension(:) :: atom
 type(node),allocatable,dimension(:) :: guest_atom
 real,allocatable    :: DistanceMatrix(:,:)
 logical,allocatable :: ConnectedAtoms(:,:)
 integer :: bond_types_max=100, bend_types_max=100, tors_types_max=100, impr_types_max=100
! bonds
 character(len=8),dimension(:), allocatable   :: bond_type_string
 integer,allocatable                          :: bond_type_histogram(:)
 character(len=8)                             :: bond_string
 integer                                      :: n_bond_types = 0
 character(len=80),dimension(:), allocatable  :: bond
 integer                                      :: n_bonds = 0
 integer                                      :: n_bonds_max = 100000
! bends
 character(len=12),dimension(:), allocatable  :: bend_type_string
 integer,allocatable                          :: bend_type_histogram(:)
 character(len=12)                            :: bend_string
 integer                                      :: n_bend_types = 0
 character(len=80),dimension(:), allocatable  :: bend
 integer                                      :: n_bends = 0
 integer                                      :: n_bends_max = 100000
! torsion (dihedrals)
 character(len=16),dimension(:), allocatable  :: tors_type_string
 integer,allocatable                          :: tors_type_histogram(:)
 character(len=16)                            :: tors_string
 integer                                      :: n_tors_types = 0
 character(len=80),dimension(:), allocatable  :: tors
 integer                                      :: n_torss = 0
 integer                                      :: n_torss_max = 100000
! torsion (impropers)
 character(len=16),dimension(:), allocatable  :: impr_type_string
 integer,allocatable                          :: impr_type_histogram(:)
 character(len=16)                            :: impr_string
 integer                                      :: n_impr_types = 0
 character(len=80),dimension(:), allocatable  :: impr
 integer                                      :: n_imprs = 0
 integer                                      :: n_imprs_max = 100000
! arguments in line
 character(len=100),dimension(:), allocatable :: args
 logical                                      :: charges_flag         = .true.
 logical                                      :: modify_topology_flag = .false.
 logical                                      :: fill_with_Ar_flag    = .false.
 logical                                      :: adsorption_fast_atom_saturation_INPUT = .false.
 logical                                      :: input_from_RASPA =.false.
 logical                                      :: input_from_MS = .false.
 logical                                      :: flag_naming
 logical                                      :: flag_lammps = .false.
 ! Element table and Matrix Topology:
 integer,allocatable  :: TopologyMatrix(:,:)
 integer,parameter    :: max_number_of_elements = 12
 integer              :: element_table(1:max_number_of_elements)
 ! C  H  O  N  P  S Zn Cd He Ar Xe Kr ...
 ! 1  2  3  4  5  6  7  8  9 10 11 12 ...
 ! 6  1  8  7 15 16 30 48  2 18 54 36 ...
!C                      H                      O                      N
 element_table(1) = 6 ; element_table(2) = 1 ; element_table(3) = 8 ; element_table(4) = 7
!P                      S                      Zn                     Cd
 element_table(5) =15 ; element_table(6) =16 ; element_table(7) =30 ; element_table(8) =48
!He                     Ne                     Xe
 element_table(9) = 2 ; element_table(10)=10 ; element_table(11)=54 ; element_table(12)=36
 !
 num_args = command_argument_count()
 allocate(args(num_args))
 do i = 1, num_args
  call get_command_argument(i,args(i))
 end do
 write(6,'(a,1x,i2)')'Arguments:',command_argument_count()
 write(6,'(a)')( args(i),i=1,num_args)
 if(num_args==0) then
  call print_help()
  stop
 end if
 do i=1,num_args
  select case(args(i))
   case ('-f','--filled','--filled-with-guest','-Ar')
    adsorption_fast_atom_saturation_INPUT=.true.
   case ('-h','--help')
    call print_help()
    stop
   case ('-c','--cif','-in')
    CIFFilename=args(i+1)
    filename=CIFFilename(1:Clen_trim(CIFFilename)-4)
    write(6,'(a)') filename
   case ('-wq','--no-charges')
    charges_flag=.false.
    if(adsorption_fast_atom_saturation_INPUT)then
     string_stop_head="_atom_site_charge"
    else
     string_stop_head= "_atom_site_fract_z"
    end if
   case ('-S','--search-and-modify-topology')
    modify_topology_flag = .true.
   case ('-F','--fill-with-Ar')
    fill_with_Ar_flag = .true.
   case ("-l","--lammps")
    flag_lammps = .true.
   case ('-R','--from-RASPA')
    input_from_RASPA=.true.
    string_stop_head= "_atom_site_charge"
   case ('-MS')
    input_from_MS=.true.
    string_stop_head="_atom_site_occupancy"
  end select
 end do
 open(100,file=CIFfilename,status='old',iostat=ierr)
 if(ierr/=0) stop 'Error opening CIF file'
 read_cif: do
  read(100,'(a)',iostat=ierr) line
  if(ierr/=0) exit read_cif
  if(line(1:14)=="_cell_length_a")then
   read(line,*)spam,cell_0(1)
   cycle read_cif
  end if
  if(line(1:14)=="_cell_length_b")then
   read(line,*)spam,cell_0(2)
   cycle read_cif
  end if
  if(line(1:14)=="_cell_length_c")then
   read(line,*)spam,cell_0(3)
   cycle read_cif
  end if
  if(line(1:17)=="_cell_angle_alpha")then
   read(line,*)spam,cell_0(4)
   cycle read_cif
  end if
  if(line(1:16)=="_cell_angle_beta")then
   read(line,*)spam,cell_0(5)
   cycle read_cif
  end if
  if(line(1:17)=="_cell_angle_gamma")then
   read(line,*)spam,cell_0(6)
   cycle read_cif
  end if
  if(line(1:)==string_stop_head) exit read_cif
 end do read_cif
 call cell(rv,vr,cell_0)
 n_atoms=0
 read_natoms: do
  read(100,'(a)',iostat=ierr) line
  if(ierr/=0.or.line(1:5)=="loop_".or.len_trim(line)==0) exit read_natoms
  n_atoms=n_atoms+1
 end do read_natoms
!
 allocate(atom(n_atoms))
 atom(1:n_atoms)%new_label="Xxxx"
 allocate(ConnectedAtoms(n_atoms,n_atoms))
 allocate(DistanceMatrix(n_atoms,n_atoms))
 allocate( TopologyMatrix(n_atoms,max_number_of_elements) )
 TopologyMatrix(1:n_atoms,1:max_number_of_elements) = 0
!
 rewind(100)
 write(6,*)'Atoms:',n_atoms
 do
  read(100,'(a)',iostat=ierr) line
  if(ierr/=0) exit
  if(line(1:)==string_stop_head) exit
 end do
 !
 atom_types(1:n_atom_types_max)(1:4)=" "
 atom_types_hybridization(1:n_atom_types_max)(1:10)=" "
 atom_types_histogram=0
 i=0
 read_atoms: do
  read(100,'(a)',iostat=ierr) line
  if(ierr/=0)                     exit read_atoms ! the end of the CIF File
  if(i>=1.and.line(1:5)=="loop_") exit read_atoms ! final loop in CIF files from Material Studio
  if(i>=1.and.len_trim(line)==0)  exit read_atoms ! finals lines in CIF files from RASPA
  i=i+1
  atom(i)%n_components=1
  if(input_from_RASPA)then
!( -R )
!loop_
!_atom_site_label
!_atom_site_type_symbol
!_atom_site_fract_x
!_atom_site_fract_y
!_atom_site_fract_z
!_atom_site_charge
   read(line,*,iostat=ierr)atom(i)%label_from_CIFFile, atom(i)%label, (atom(i)%xyzs(j,1),j=1,3), atom(i)%charge
  elseif(input_from_MS)then
!( -MS )
!loop_
!_atom_site_label
!_atom_site_type_symbol
!_atom_site_fract_x
!_atom_site_fract_y
!_atom_site_fract_z
!_atom_site_U_iso_or_equiv
!_atom_site_adp_type
!_atom_site_occupancy
   read(line,*,iostat=ierr)atom(i)%label_from_CIFFile, atom(i)%label, (atom(i)%xyzs(j,1),j=1,3), xy, spams, xz
  elseif(charges_flag) then
   read(line,*,iostat=ierr)atom(i)%label,(atom(i)%xyzs(j,1),j=1,3)
  else if(adsorption_fast_atom_saturation_INPUT) then
   read(line,*,iostat=ierr)atom(i)%label_from_CIFFile, atom(i)%label, (atom(i)%xyzs(j,1),j=1,3), atom(i)%charge
  else
   read(line,*,iostat=ierr)atom(i)%label,(atom(i)%xyzs(j,1),j=1,3)
   atom(i)%charge=0.0
  end if
  if(ierr/=0) exit read_atoms
  call CheckAtom(atom(i)%label,atom(i)%mass,atom(i)%radius,atom(i)%element,atom(i)%label_element)
  do j=1,3
   atom(i)%xyzc(j,1)=rv(j,1)*atom(i)%xyzs(1,1)+rv(j,2)*atom(i)%xyzs(2,1)+rv(j,3)*atom(i)%xyzs(3,1)
  end do
 end do read_atoms
 close(100)
! topology:
 DistanceMatrix=0.0
 ConnectedAtoms=.false.
 do i=1,n_atoms
  k=0
  do j=1,n_atoms
   forall (h=1:3)
    rouratom(h)=atom(i)%xyzs(h,1)
    ratom(h)=atom(j)%xyzs(h,1)
   end forall
   call make_distances(cell_0,ratom,rouratom,rv,r)
   DistanceMatrix(i,j)=r
   if(r>0.1.and.r<=atom(i)%radius+atom(j)%radius+r_min_criteria_connectivity)then
    k=k+1
    ConnectedAtoms(i,j)=.true.
   end if
  end do
  atom(i)%degree=k
 end do
! recheck:
! No bonds by constrains:
! Zn...OX, Zn...HX, Zn...CX
 do i=1,n_atoms
  do j=1,n_atoms
   if(i/=j.and.ConnectedAtoms(i,j))then
    if(((atom(i)%element==1.and.atom(j)%element==30).or. &
        (atom(i)%element==30.and.atom(j)%element==1)).or.&
!
       ((atom(i)%element==6.and.atom(j)%element==30).or. &
        (atom(i)%element==30.and.atom(j)%element==6)).or.&
!
       ((atom(i)%element==8.and.atom(j)%element==30).or. &
        (atom(i)%element==30.and.atom(j)%element==8)).or.&
!
       ((atom(i)%element==1.and.atom(j)%element==54).or. &
        (atom(i)%element==54.and.atom(j)%element==1)).or.&
! He ... Ne (Zn_s , N_s)
       ((atom(i)%element==10.and.atom(j)%element==2).or. &
        (atom(i)%element==2.and.atom(j)%element==10)).or.&
! 
       ((atom(i)%element==10.and.atom(j)%element==6).or. & ! 
        (atom(i)%element==6.and.atom(j)%element==10)).or.&
!
       ((atom(i)%element==2.and.atom(j)%element==6).or.  & ! H ... He (Zn_s)
        (atom(i)%element==6.and.atom(j)%element==2)).or. &
!
       ((atom(i)%element==1.and.atom(j)%element==2).or.  & ! H ... He (Zn_s)
        (atom(i)%element==2.and.atom(j)%element==1)).or. &
!
       ((atom(i)%element==1.and.atom(j)%element==10).or. & ! H ... Ne (N_s)
        (atom(i)%element==10.and.atom(j)%element==1)).or.&
!
       ((atom(i)%element==7.and.atom(j)%element==30).or. & ! N ... Zn
        (atom(i)%element==30.and.atom(j)%element==7)).or.&
!
       ((atom(i)%element==10.and.atom(j)%element==30).or.& ! Ne ... Zn
       (atom(i)%element==30.and.atom(j)%element==10)).or.&
!
       ((atom(i)%element==7.and.atom(j)%element==1).or.  & ! N ... H
        (atom(i)%element==1.and.atom(j)%element==7)).or. &
!
       ((atom(i)%element==7.and.atom(j)%element==2).or.  &
        (atom(i)%element==2.and.atom(j)%element==7)).or. &
!
       (atom(i)%element==7.and.atom(j)%element==7).or.   &
       (atom(i)%element==1.and.atom(j)%element==1).or.   &
!
       ((atom(i)%element==54.and.atom(j)%element/=36).or. &
        (atom(i)%element==36.and.atom(j)%element/=54))    &
!
       )  then
     ConnectedAtoms(i,j)=.false.
     ConnectedAtoms(j,i)=.false.
     atom(i)%degree=atom(i)%degree-1
     atom(j)%degree=atom(j)%degree-1
    end if
   end if
  end do
 end do
!
 call check_bonds_with_degree(2,4)  ! He dummy atom of Zn
 call check_bonds_with_degree(10,1) ! Ne dummy atom of N
 call check_bonds_with_degree(30,4) ! Zn
 call check_bonds_with_degree(54,1)
 call check_bonds_with_degree(36,1) 
!
 call check_bonds_with_degree(1,1)  ! H
 call check_bonds_with_degree(7,3)  ! N ! resonant
! degree of each node:
 write(6,'(a)') 'Connectivity array for each atom:'
 write(6,'(1000(i2))')(atom(i)%degree,i=1,n_atoms)
! rename atoms with for certain forcefield
 if ( modify_topology_flag ) then
  !
  write(6,'(a4,1x,a3,1x,12(i2,1x))')'Atom','Z',(element_table(i),i=1,max_number_of_elements)
  absolute_scan: do i=1,n_atoms
   do j=1,n_atoms
    if(i/=j.and.ConnectedAtoms(i,j) )then
     do k=1,max_number_of_elements
      if(  atom(j)%element == element_table(k) )  TopologyMatrix(i,k) = TopologyMatrix(i,k) + 1
     end do
    end if
   end do
   write(atom(i)%topol_label,'(a,a,12(i1,a))')atom(i)%label_element,'@',&
    (TopologyMatrix(i,k),'_', k=1,max_number_of_elements )
   atom(i)%topol_label=adjustl(trim(atom(i)%topol_label))
   write(6,'(i4,1x,a,1x,a1,1x,12(i2,1x),a1,1x,a1,1x,a80)')i,atom(i)%label_element,'(',&
    ( TopologyMatrix(i,k),k=1,max_number_of_elements ),&
    ')',':',atom(i)%topol_label
  end do absolute_scan
  !
  scan_for_label_1st_phase: do i=1,n_atoms
  select case(atom(i)%topol_label)
   ! code:
   ! element@C_H_O_N_P_S_Zn_Cd_He_Ne_Xe_Kr
   ! C  H  O  N  P  S Zn Cd He Ar Xe Kr ...
   ! 1  2  3  4  5  6  7  8  9 10 11 12 ...
   ! 6  1  8  7 15 16 30 48  2 18 54 36 ...
   ! Nitrogen:
   case( "N@2_0_0_0_0_0_0_0_0_1_0_0_")
    atom(i)%new_label     = "N1  "
    atom(i)%charge        = 0.0
    atom(i)%hybridization = "N_c_aromatic"
   case("Ne@0_0_0_1_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = "Ne  "
    atom(i)%charge        = -0.42030
    atom(i)%hybridization = "N_s"
   ! Carbons:
   case( "C@1_0_0_2_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = "C1  "
    atom(i)%charge        = 0.43750  !+0.4291  ! 0.43750
    atom(i)%hybridization = "CC_aromatic"
    !
   case( "C@1_1_0_1_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = "C2  "
    atom(i)%charge        = -0.06620 !-0.0839  !-0.06620
    atom(i)%hybridization = "CH_aromatic"
    !
   case( "C@1_3_0_0_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = "C3  "
    atom(i)%charge        = -0.46060 !-0.4526
    atom(i)%hybridization = "C3_sp3"
   ! Oxygens:
   case( "O@0_2_0_0_0_0_0_0_4_0_0_0_") ! water oxygen:
    !  TIP4P/2005 model 
    atom(i)%new_label     = "O1  "
    atom(i)%charge        = -1.1128
    atom(i)%hybridization = "O_w"
   !
   case("Xe@0_0_0_0_0_0_0_0_0_0_0_1_")
    ! fake DMF (negative particle)
    atom(i)%new_label     = "Xe  "
    atom(i)%charge        = -0.35
    atom(i)%hybridization = "-_DMF_fake"
   case("Kr@0_0_0_0_0_0_0_0_0_0_1_0_")
    ! fake DMF (positive particle)
    atom(i)%new_label     = "Kr  "
    atom(i)%charge        = +0.35
    atom(i)%hybridization = "+_DMF_fake"
   ! Metals:
   ! Scenario 2, big central atom
   case("Zn@0_0_0_0_0_0_0_0_4_0_0_0_")
    atom(i)%new_label     = 'Zn  '
    atom(i)%charge        = +0.7072/2.0
    atom(i)%hybridization = "Zn_c"
   case("He@0_0_0_0_0_0_1_0_3_0_0_0_")
    atom(i)%new_label     = 'He  '
    atom(i)%charge        = +0.7072/2.0/4.0
    atom(i)%hybridization = "Zn_s"
   end select
  end do scan_for_label_1st_phase
!
  scan_for_rename_H_atoms: do i = 1,n_atoms
   if ( atom(i)%element == 1) then
    atom(i)%topol_label = "H_   "
    do j=1,n_atoms
     if(j/=i.and.ConnectedAtoms(i,j))then
      if( atom(j)%new_label == "C2  " )      then
       atom(i)%new_label = "H2  "
       atom(i)%charge =  0.11410  !+0.1128
       atom(i)%hybridization = "HC_aromatic"
      else if( atom(j)%new_label == "C4  " ) then
       atom(i)%new_label = "H4  "
       atom(i)%charge =  0.11410  !+0.1128
       atom(i)%hybridization = "HC_aromatic"
      else if( atom(j)%new_label == "C3  " ) then
       atom(i)%new_label = "H3  "
       atom(i)%charge =  0.13810  !+0.131866667
       atom(i)%hybridization = "HC_sp3"
      else if( atom(j)%new_label == "O1  ")then
       atom(i)%new_label = "H1  "
       atom(i)%charge = +0.5564
       atom(i)%hybridization = "H_w"
      else
       write(6,*)'============================================================'
       write(6,*)'unknow H-atom'
       write(6,*)'Bond:', atom(j)%label,'(',atom(j)%new_label,')',atom(j)%degree
       write(6,*)'Distance:',DistanceMatrix(i,j),ConnectedAtoms(i,j)
       write(6,*)'------------------------------------------------------------'
       do jj=1,n_atoms
       if(DistanceMatrix(i,jj)>0.1.and.DistanceMatrix(i,jj)<=2.0)then
        write(6,*)'Around:', atom(jj)%label,'(',atom(jj)%new_label,')',atom(jj)%degree
        write(6,*)'Distance:', DistanceMatrix(i,jj),ConnectedAtoms(i,jj)
        write(6,*)'------------------------------------------------------------'
       end if
       end do
       stop 'unknow H-atom'
      end if
      cycle scan_for_rename_H_atoms
     end if
    end do
   end if
  end do scan_for_rename_H_atoms
  do i=1,n_atoms
   if(atom(i)%new_label/="Xxxx")then
    atom(i)%label=atom(i)%new_label
   end if
  end do
! }}
 end if
! atom types:
 do i = 1,n_atoms
  check_atom_type: if(i==1)then
   n_atom_types=1
   atom_types(n_atom_types)=atom(i)%label
   atom_types_hybridization(n_atom_types)=atom(i)%hybridization 
   atom(i)%type_=n_atom_types
  else
   n_atom_types=n_atom_types+1
   do h=1,n_atom_types-1
    if(atom(i)%label==atom_types(h))then
     n_atom_types=n_atom_types-1
     atom(i)%type_=h
     exit check_atom_type
    end if
   end do
   atom_types(n_atom_types)=atom(i)%label
   atom_types_hybridization(n_atom_types)=atom(i)%hybridization
   atom(i)%type_=n_atom_types
  end if check_atom_type
  ! print all atoms
  write(6,'(i6,1x,a4,1x,3(f10.8,1x),a2,1x,3(f10.6,1x))') &
   i,atom(i)%label,(atom(i)%xyzs(j,1),j=1,3),&
   atom(i)%label_element,atom(i)%charge,atom(i)%radius,atom(i)%mass
 end do
!
 allocate(table_atoms_atom_types( 1:n_atom_types,1:n_atoms))
 do i=1,n_atoms
  do j=1,n_atom_types
   if(atom(i)%label==atom_types(j)) then
    atom_types_histogram(j)=atom_types_histogram(j)+1
    table_atoms_atom_types(j, atom_types_histogram(j) ) = i
   end if
  end do
 end do
!
 open(789,file="index.ndx")
 do i=0,n_atom_types
  select case(i)
   case(0)
    write(789,'(a)')"[ System ]"
    write(789,'(15(i6,1x))')( j, j=1,n_atoms )
   case default
    write(789,'(a1,1x,a2,1x,a1)')"[",atom_types(i)(1:2),"]" 
    write(789,'(15(i6,1x))')( table_atoms_atom_types( i, j),j=1,atom_types_histogram(i) )
  end select
 end do
 close(789)
 deallocate( table_atoms_atom_types )
!
 write(6,'(a)')'=========='
 write(6,'(a,1x,i2)')'atom_types:',n_atom_types
 do i=1,n_atom_types
  write(6,'(i2,1x,a4,1x,a10,1x,i6)')i,atom_types(i), atom_types_hybridization(i),&
   atom_types_histogram(i)
 end do
 open(987,file="atom_types_for_dump.txt")
  write(987,'(100(a4,1x))') ( atom_types(i), i=1,n_atom_types)
 close(987)
 write(6,'(a)')'=========='
 do i=1,n_atoms
  total_Q = total_Q + atom(i)%charge
 end do
 write(6,'(a,1x,f14.7)')'total charge:',total_Q
 if( abs(total_Q) > 1.0 ) write(6,'(a)')'[WARNING] the system has net charge !!'
 write(6,'(a)')'=========='
! lammps:
 if ( flag_lammps ) then
! bonds:
 write(6,'(a)')'Detecting bonds [...]'
 allocate(bond_type_string(bond_types_max))
 allocate(bond_type_histogram(bond_types_max))
 allocate(bond(n_bonds_max))
 bond_type_histogram=0
 bond_string(1:8)="        "
 do i=1,bond_types_max
  do j=1,8
   write(bond_type_string(j:j),'(a1)')" "
  end do
 end do
 forall (i=1:n_bonds_max)
  forall (j=1:80)
   bond(i)(j:j)=" "
  end forall
 end forall
 do_i: do i=1,n_atoms
  do_j: do j=i+1,n_atoms
   if(ConnectedAtoms(i,j))then
    write(bond_string(1:8),'(2a4)') atom(i)%label(1:4),atom(j)%label(1:4)
    forall (l=1:100)
     line(l:l)=" "
    end forall
    if(n_bond_types==0)then
     n_bonds=n_bonds+1
     n_bond_types=n_bond_types+1
     bond_type_string(1)=bond_string
     bond_type_histogram(1)=1
     write(bond(n_bonds)(1:21),'(3i7)')n_bond_types,i,j
     write(bond(n_bonds)(22:24),'(a3)')' # '
     write(bond(n_bonds)(25:),'(a)') bond_string
     cycle do_j
    else
     n_bonds=n_bonds+1
     n_bond_types=n_bond_types+1 ! try
     check_bond: do h=n_bond_types-1,1,-1
      if(bond_string(1:8)==bond_type_string(h)(1:8).or.&
         bond_string(1:8)==bond_type_string(h)(5:8)//bond_type_string(h)(1:4))then
       n_bond_types=n_bond_types-1
       bond_type_histogram(h)=bond_type_histogram(h)+1
       write(bond(n_bonds)(1:21),'(3i7)')h,i,j
       write(bond(n_bonds)(22:24),'(a3)')' # '
       write(bond(n_bonds)(25:),'(a)') bond_string
       cycle do_j
      end if
     end do check_bond
     bond_type_histogram(n_bond_types)=1
     bond_type_string(n_bond_types)=bond_string
     write(bond(n_bonds)(1:21),'(3i7)')n_bond_types,i,j
     write(bond(n_bonds)(22:24),'(a3)')' # '
     write(bond(n_bonds)(25:),'(a)') bond_string
    end if
   end if
  end do do_j
 end do do_i
!
 write(6,*)'Bond types:',n_bond_types
 write(6,*)'bonds:',n_bonds
 do i=1,n_bond_types
  write(6,*) bond_type_string(i),bond_type_histogram(i)
 end do
!
! bends:
 write(6,'(a)')'Detecting bends [...]'
 allocate(bend_type_string(bend_types_max))
 allocate(bend_type_histogram(bend_types_max))
 allocate(bend(n_bends_max))
 bend_type_histogram=0
 bend_string(1:12)="            "
 do i=1,bend_types_max
  do j=1,12
   write(bend_type_string(j:j),'(a1)')" "
  end do
 end do
 forall (i=1:n_bends_max)
  forall (j=1:80)
   bend(i)(j:j)=" "
  end forall
 end forall
 do_i_bend: do i=1,n_atoms ! central atom
  do_j_bend: do j=1,n_atoms
   if(ConnectedAtoms(i,j).and.j/=i)then
    do_k_bend: do k=1,n_atoms
     if(ConnectedAtoms(i,k).and.k/=j.and.k/=i)then
      write(bend_string(1:12),'(3a4)') atom(j)%label(1:4),atom(i)%label(1:4),atom(k)%label(1:4)
      if(n_bend_types==0)then
       n_bends=n_bends+1
       n_bend_types=n_bend_types+1
       bend_type_string(1)=bend_string
       bend_type_histogram(1)=1
       write(bend(n_bends)(1:28),'(4i7)')n_bend_types,j,i,k
       write(bend(n_bends)(29:31),'(a3)')' # '
       write(bend(n_bends)(32:),'(a)') bend_string
       cycle do_k_bend
      else
       n_bends=n_bends+1
       n_bend_types=n_bend_types+1 ! try
       check_bend_type: do h=n_bend_types-1,1,-1
        if(bend_string(1:12)==bend_type_string(h)(1:12).or.&
           bend_string(1:12)==bend_type_string(h)(9:12)//bend_type_string(h)(5:8)//&
           bend_type_string(h)(1:4))then
         n_bend_types=n_bend_types-1
         write(bend(n_bends)(1:28),'(4i7)')h,j,i,k
         write(bend(n_bends)(29:31),'(a3)')' # '
         write(bend(n_bends)(32:),'(a)') bend_string
         check_bend: do l=1,n_bends-1
          read(bend(l)(8:28),'(3i7)')jj,ii,kk
          if((jj==j.and.ii==i.and.kk==k).or.&
             (kk==j.and.ii==i.and.jj==k))then
           n_bends=n_bends-1
           cycle do_k_bend
          end if
         end do check_bend
         bend_type_histogram(h)=bend_type_histogram(h)+1
         cycle do_k_bend
        end if
       end do check_bend_type
       bend_type_histogram(n_bend_types)=1
       bend_type_string(n_bend_types)=bend_string
       write(bend(n_bends)(1:28),'(4i7)')n_bend_types,j,i,k
       write(bend(n_bends)(29:31),'(a3)')' # '
       write(bend(n_bends)(32:),'(a)') bend_string
      end if
     end if
    end do do_k_bend
   end if
  end do do_j_bend
 end do do_i_bend
!
 write(6,*)'Bend types:',n_bend_types
 write(6,*)'bends:',n_bends
 do i=1,n_bend_types
  write(6,*) bend_type_string(i),bend_type_histogram(i)
 end do
!
! dihedrals
 write(6,'(a)')'Detecting dihedral torsions [...]'
 allocate(tors_type_string(tors_types_max))
 allocate(tors_type_histogram(tors_types_max))
 allocate(tors(n_torss_max))
 tors_type_histogram=0
 tors_string(1:16)="                "
 do i=1,tors_types_max
  do j=1,16
   write(tors_type_string(j:j),'(a1)')" "
  end do
 end do
 forall (i=1:n_torss_max)
  forall (j=1:80)
   tors(i)(j:j)=" "
  end forall
 end forall
 do_i_tors: do i=1,n_atoms
  do_j_tors: do j=1,n_atoms
   if(ConnectedAtoms(i,j).and.j/=i)then
    do_k_tors: do k=1,n_atoms
     if(ConnectedAtoms(j,k).and.k/=j.and.k/=i)then
      do_l_tors: do l=1,n_atoms
       if(ConnectedAtoms(k,l).and.l/=k.and.l/=j.and.l/=i)then
        write(tors_string(1:16),'(4a4)')atom(i)%label(1:4),atom(j)%label(1:4),&
                                        atom(k)%label(1:4),atom(l)%label(1:4)
        if(n_tors_types==0)then
         n_torss=n_torss+1
         n_tors_types=n_tors_types+1
         tors_type_string(1)=tors_string
         tors_type_histogram(1)=1
         write(tors(n_torss)(1:35),'(5i7)')n_tors_types,i,j,k,l
         write(tors(n_torss)(36:38),'(a3)')' # '
         write(tors(n_torss)(39:),'(a)') tors_string
         cycle do_l_tors
        else
         n_torss=n_torss+1
         n_tors_types=n_tors_types+1 ! try
         check_tors_type: do h=n_tors_types-1,1,-1
          if(tors_string(1:16)==tors_type_string(h)(1:16).or.&
             tors_string(1:16)==tors_type_string(h)(13:16)//tors_type_string(h)(9:12)//&
                                tors_type_string(h)(5:8)//tors_type_string(h)(1:4))then
           n_tors_types=n_tors_types-1
           write(tors(n_torss)(1:35),'(5i7)')h,i,j,k,l
           write(tors(n_torss)(36:38),'(a3)')' # '
           write(tors(n_torss)(39:),'(a)') tors_string
           check_tors: do m=1,n_torss-1
            read(tors(m)(8:35),'(4i7)')ii,jj,kk,ll
            if((ii==i.and.jj==j.and.kk==k.and.ll==l).or.&
               (ll==i.and.kk==j.and.jj==k.and.ii==l))then
             n_torss=n_torss-1
             cycle do_l_tors
            end if
           end do check_tors
           tors_type_histogram(h)=tors_type_histogram(h)+1
           cycle do_l_tors
          end if
         end do check_tors_type
         tors_type_histogram(n_tors_types)=1
         tors_type_string(n_tors_types)=tors_string
         write(tors(n_torss)(1:35),'(5i7)') n_tors_types,i,j,k,l
         write(tors(n_torss)(36:38),'(a3)')' # '
         write(tors(n_torss)(39:),'(a)') tors_string
        end if
       end if
      end do do_l_tors
     end if
    end do do_k_tors
   end if
  end do do_j_tors
 end do do_i_tors
!
 write(6,*)'Dihedral types:',n_tors_types
 write(6,*)'dihedrals:',n_torss
 do i=1,n_tors_types
  write(6,*) tors_type_string(i),tors_type_histogram(i)
 end do
!
! impropers
 write(6,'(a)')'Detecting improper torsions [...]'
 allocate(impr_type_string(impr_types_max))
 allocate(impr_type_histogram(impr_types_max))
 allocate(impr(n_imprs_max))
 impr_type_histogram=0
 impr_string(1:16)="                "
 do i=1,impr_types_max
  do j=1,16
   write(impr_type_string(j:j),'(a1)')" "
  end do
 end do
 forall (i=1:n_imprs_max)
  forall (j=1:80)
   impr(i)(j:j)=" "
  end forall
 end forall
!   l
!  .|.
!   i - j
!   |
!   k
 do_i_impr: do i=1,n_atoms
  do_j_impr: do j=1,n_atoms
   if(ConnectedAtoms(i,j).and.j/=i)then
    do_k_impr: do k=1,n_atoms
     if(ConnectedAtoms(i,k).and.k/=j.and.k/=i)then
      do_l_impr: do l=1,n_atoms
       if(ConnectedAtoms(i,l).and.l/=k.and.l/=j.and.l/=i)then
        write(impr_string(1:16),'(4a4)')atom(i)%label(1:4),atom(j)%label(1:4),&
                                        atom(k)%label(1:4),atom(l)%label(1:4)
        if(n_impr_types==0)then
         n_imprs=n_imprs+1
         n_impr_types=n_impr_types+1
         impr_type_string(1)=impr_string
         impr_type_histogram(1)=1
         write(impr(n_imprs)(1:35),'(5i7)')n_impr_types,i,j,k,l
         write(impr(n_imprs)(36:38),'(a3)')' # '
         write(impr(n_imprs)(39:),'(a)') impr_string
         cycle do_l_impr
        else
         n_imprs=n_imprs+1
         n_impr_types=n_impr_types+1 ! try
         check_impr_type: do h=n_impr_types-1,1,-1
          !
          sc(1)=impr_type_string(h)(1:4)
          sc(2)=impr_type_string(h)(5:8)
          sc(3)=impr_type_string(h)(9:12)
          sc(4)=impr_type_string(h)(13:16)
          !
          if( impr_string(1:16)== impr_type_string(h)(1:16).or. &
              impr_string(1:16)== sc(1)//sc(2)//sc(3)//sc(4).or.&
              impr_string(1:16)== sc(1)//sc(3)//sc(4)//sc(2).or.&
              impr_string(1:16)== sc(1)//sc(4)//sc(2)//sc(3).or.&
              impr_string(1:16)== sc(1)//sc(4)//sc(3)//sc(2).or.&
              impr_string(1:16)== sc(1)//sc(3)//sc(2)//sc(4).or.&
              impr_string(1:16)== sc(1)//sc(2)//sc(4)//sc(3) ) then
           !
           n_impr_types=n_impr_types-1
           !
           write(impr(n_imprs)(1:35),'(5i7)')h,i,j,k,l
           write(impr(n_imprs)(36:38),'(a3)')' # '
           write(impr(n_imprs)(39:),'(a)') impr_string
           !
           check_impr: do m=1,n_imprs-1
            read(impr(m)(8:35),'(4i7)')ii,jj,kk,ll
            if((ii==i.and.jj==j.and.kk==k.and.ll==l).or.&
               (ii==i.and.kk==j.and.jj==k.and.ll==l).or.&
               (ii==i.and.jj==j.and.ll==k.and.kk==l).or.&
               (ii==i.and.ll==j.and.jj==k.and.kk==l).or.&
               (ii==i.and.kk==j.and.ll==k.and.jj==l).or.&
               (ii==i.and.ll==j.and.kk==k.and.jj==l))then
             n_imprs=n_imprs-1
             cycle do_l_impr
            end if
           end do check_impr
           impr_type_histogram(h)=impr_type_histogram(h)+1
           cycle do_l_impr
          end if
         end do check_impr_type
         impr_type_histogram(n_impr_types)=1
         impr_type_string(n_impr_types)=impr_string
         write(impr(n_imprs)(1:35),'(5i7)') n_impr_types,i,j,k,l
         write(impr(n_imprs)(36:38),'(a3)')' # '
         write(impr(n_imprs)(39:),'(a)') impr_string
        end if
       end if
      end do do_l_impr
     end if
    end do do_k_impr
   end if
  end do do_j_impr
 end do do_i_impr
!
 write(6,*)'Improper types:',n_impr_types
 write(6,*)'impropers:',n_imprs
 do i=1,n_impr_types
  write(6,*) impr_type_string(i),impr_type_histogram(i)
 end do
!
 write(6,'(a)')'Generating outputs:'
!
! output:
 call cellnormal2lammps(cell_0,xlo_bound,ylo_bound,zlo_bound,&
      xhi_bound,yhi_bound,zhi_bound,xy,xz,yz)
 call output_lammps()
 deallocate(bend_type_string)
 deallocate(bend_type_histogram)
 deallocate(bend)
!
 deallocate(impr)
 deallocate(impr_type_string)
 deallocate(impr_type_histogram)
!
 deallocate(tors)
 deallocate(tors_type_string)
 deallocate(tors_type_histogram)
 !
 deallocate(bond_type_string)
 deallocate(bond_type_histogram)
 deallocate(bond)
 end if
 !call output_gulp()
 call output_pdb()
 call output_CIF()
 deallocate(TopologyMatrix)
 deallocate(DistanceMatrix)
 deallocate(ConnectedAtoms)
 deallocate(atom)
 stop 'Done'
 contains
!
 subroutine remove_bond_in_atom( III )
  implicit none
  integer,intent(in) :: III
  integer            :: i,j,k
  logical            :: flag_degree_mode = .true.
  write(6,*)'=============================='
  write(6,*)'Strange ',atom(III)%element,'-atom:',III,atom(III)%degree, atom(III)%label_element
  allocate( check_topol_id(1:atom(III)%degree))
  allocate( check_topol_distance(1:atom(III)%degree) )
  allocate( check_topol_degree(1:atom(III)%degree) )
  j=0
  do i=1,n_atoms
   if(III/=i.and.ConnectedAtoms(III,i))then
    j=j+1
    check_topol_id(j)=i
    check_topol_distance(j)=DistanceMatrix(III,i)
    check_topol_degree(j)=atom(i)%degree
    write(6,*)atom(i)%label,DistanceMatrix(III,i),j,atom(i)%degree
   end if
  end do
  if( flag_degree_mode) then
  fdhgshd: do i=1,atom(III)%degree
   j=check_topol_id(i)
   if( atom(i)%degree == maxval( check_topol_degree ) ) then
    ConnectedAtoms(III,j)=.false.
    ConnectedAtoms(j,III)=.false.
    atom(III)%degree=atom(III)%degree-1
    atom(j)%degree=atom(j)%degree-1
    write(6,*)'Solution (highest atom degree):'
    write(6,*)'remove bond:',III,j,atom(III)%label,atom(j)%label, DistanceMatrix(III,j), maxval( check_topol_distance )
    exit fdhgshd
   end if
  end do fdhgshd
  else
  ghjxxxx: do i=1,atom(III)%degree
   j=check_topol_id(i)
   if( DistanceMatrix(III,j) == maxval( check_topol_distance ) ) then
    ConnectedAtoms(III,j)=.false.
    ConnectedAtoms(j,III)=.false.
    atom(III)%degree=atom(III)%degree-1
    atom(j)%degree=atom(j)%degree-1
    write(6,*)'Solution (remove the farthest atom):'
    write(6,*)'remove bond:',III,j,atom(III)%label,atom(j)%label, DistanceMatrix(III,j), maxval( check_topol_distance )
    exit ghjxxxx
   end if
  end do ghjxxxx
  end if
  deallocate( check_topol_id )
  deallocate( check_topol_distance )
  deallocate( check_topol_degree )
  return
 end subroutine remove_bond_in_atom
!
 subroutine check_bonds_with_degree(ZZZ,KKK)
 implicit none
 integer,intent(in) :: ZZZ,KKK
 integer            :: m = 0
 logical            :: flag_degree_mode = .true.
 do i=1,n_atoms
  do while ( atom(i)%element == ZZZ .and. atom(i)%degree > KKK )
   write(6,*)'=============================='
   write(6,*)'Strange ',ZZZ,'-atom:',i,atom(i)%degree
   h=0
   allocate( check_topol_id(1:atom(i)%degree) )
   allocate( check_topol_distance(1:atom(i)%degree) )
   allocate( check_topol_degree(1:atom(i)%degree) )
   do j=1,n_atoms
    if(i/=j.and.ConnectedAtoms(i,j))then
     h=h+1
     check_topol_id(h)=j
     check_topol_distance(h)=DistanceMatrix(i,j)
     check_topol_degree(h)=atom(j)%degree
     write(6,*)atom(j)%label,DistanceMatrix(i,j),h,r,atom(j)%degree
    end if
   end do
   if( flag_degree_mode) then
    ! highest atom degree
    m = 99999999
    do h=1,atom(i)%degree
     j=check_topol_id(h)
     if( check_topol_degree(h) == maxval( check_topol_degree ) ) then
      m=j
      ConnectedAtoms(i,j)=.false.
      ConnectedAtoms(j,i)=.false.
      atom(i)%degree=atom(i)%degree-1
      atom(j)%degree=atom(j)%degree-1
     end if
    end do
    if ( m == 99999999 ) stop 'infinite loop !'
   else
    ! the farthest atom
    do h=1,atom(i)%degree
     j=check_topol_id(h)
     if( DistanceMatrix(i,j) == maxval( check_topol_distance ) ) then
      m=j
      ConnectedAtoms(i,j)=.false.
      ConnectedAtoms(j,i)=.false.
      atom(i)%degree=atom(i)%degree-1
      atom(j)%degree=atom(j)%degree-1
     end if
    end do
   end if
   if(flag_degree_mode) then
    write(6,*)'Solution (highest atom degree):'
   else
    write(6,*)'Solution (remove the farthest atom):'
   end if
   write(6,*)'remove bond:', i,m,atom(i)%label,atom(m)%label
   deallocate( check_topol_id )
   deallocate( check_topol_distance )
   deallocate( check_topol_degree )
  end do
 end do
 return
 end subroutine check_bonds_with_degree
!
 subroutine search_forcefield(string,interaction,label)
  implicit none
  character(len=80),intent(out):: string
  real                         :: p1,p2,p3,p5
  integer                      :: p0,p4
  character(len=10)            :: type_fff
  character(len=20)            :: mixing_rule
  character(len=4),intent(in)  :: interaction
  character(len=4),intent(in)  :: label(1:4)
  character(len=4)             :: ourlabel(1:4)
  character(len=80)            :: units
  character(len=80)            :: passing
  character(len=80)            :: variables
  integer                      :: u=123,ierr=0
  forall (ii=1:80)
   string(ii:ii)=" "
   units(ii:ii)=" "
   passing(ii:ii)=" "
   variables(ii:ii)=" "
  end forall
  open(u,file='forcefield.lib')
  fff: select case (interaction)
   case('bond')
    initbond: do
     read(u,'(a)') line
     if(line(1:11)=='Bond Coeffs') then
      read(line(12:),'(a)') units
      units=adjustl(trim(units))
      exit initbond
     end if
    end do initbond
    do
     read(u,'(a)') line
     if(line(1:11)=='Angle Coeffs'.or.line(1:3)=="End") then
      string(1:4)="none"
      exit fff
     else if(line(1:1)/="#".or.line(1:1)/=" ")then
      read(line,'(2(a4,1x),a)')ourlabel(1),ourlabel(2),passing
      if((adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(2)))).or.&
         (adjustl(trim(ourlabel(1)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(1)))))then
       read(passing,'(a,a)')type_fff,variables
       if(adjustl(trim(type_fff))/='none')then
        select case(adjustl(trim(type_fff)))
        case("harmonic")
         read(variables,*) p1,p2
         select case(units)
          case('kcal/mol','kcal/mol/A/A')
           p1=kcalmol2ev(p1) ! kcal/mol/A/A
           p2=p2             ! A
          case('K','K/A/A')
           p1=0.5*K2eV(p1)
           p2=p2
         end select
         write(string,*) p1,p2,' # ',type_fff
        case default
         write(6,'(a,1x,a)')'Unknown bond potential:', adjustl(trim(type_fff))
         stop
        end select
       else
        write(string,*) 0.0,0.0,' # ',type_fff
       end if
       exit fff
      end if
     end if
    end do
   case('bend')
    initbend: do
     read(u,'(a)') line
     if(line(1:12)=='Angle Coeffs') then
      read(line(13:),'(a)') units
      units=adjustl(trim(units))
      !write(6,*) units
      exit initbend
     end if
    end do initbend
    do
     read(u,'(a)') line
     if(line(1:15)=='Dihedral Coeffs'.or.line(1:3)=="End") then
      string(1:4)="none"
      exit fff
     end if
     if(line(1:1)/="#".or.line(1:1)/=" ")then
      read(line,'(3(a4,1x),a)')ourlabel(1),ourlabel(2),ourlabel(3),passing
      if((adjustl(trim(ourlabel(2)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(3)))).or.&
         (adjustl(trim(ourlabel(2)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(1)))==adjustl(trim(label(3))).and.&
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(1)))))then
       read(passing,'(a,a)')type_fff,variables
       if(adjustl(trim(type_fff))/='none')then
        select case(adjustl(trim(type_fff)))
        case("harmonic")
         read(variables,*) p1,p2
         select case(units)
          case('kcal/mol','kcal/mol/rad/rad')
           p1=kcalmol2ev(p1) ! kcal/mol/A/A
           p2=p2             ! A
          case('K','K/rad/rad')
           p1=0.5*K2eV(p1)
           p2=p2
         end select
         write(string,*) p1,p2,' # ',type_fff
        case("charmm")
         read(variables,*) p1,p2,p3,p5
         select case(units)
          case('kcal/mol','kcal/mol/rad/rad')
           p1=kcalmol2ev(p1) ! kcal/mol/A/A
           p2=p2
           p3=kcalmol2ev(p3) ! kcal/mol/A/A
           p5=p5
         end select
         write(string,'(4(f14.7,1x),a,a)') p1,p2,p3,p5,' # ',type_fff
        case default
         write(6,'(a,1x,a)')'Unknown bend potential:', adjustl(trim(type_fff))
         stop
        end select
       else
        write(string,*)0.0,0.0,' # ',type_fff
       end if
       exit fff
      end if
     end if
    end do
   case('tors')
    inittors: do
     read(u,'(a)')line
     if(line(1:15)=='Dihedral Coeffs') then
      read(line(16:),'(a)') units
      units=adjustl(trim(units))
      exit inittors
     end if
    end do inittors
    do
     read(u,'(a)') line
     if(line(1:15)=='Improper Coeffs'.or.line(1:3)=="End") then
      string(1:4)="none"
      exit fff
     else if(line(1:1)/="#".or.line(1:1)/=" ")then
      read(line,'(4(a4,1x),a)')(ourlabel(ii),ii=1,4),passing
      if((adjustl(trim(ourlabel(2)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(3))).and.&
          adjustl(trim(ourlabel(4)))==adjustl(trim(label(4)))).or.&
         (adjustl(trim(ourlabel(1)))==adjustl(trim(label(4))).and.&
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(3))).and.&
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(4)))==adjustl(trim(label(1)))))then
       read(passing,'(a,a)')type_fff,variables
       if(adjustl(trim(type_fff))/='none')then
        select case(adjustl(trim(type_fff)))
        case("fourier")
         read(variables,*) p0,p1,p4,p2
         select case(units)
          case('kcal/mol')
           p0=p0
           p1=kcalmol2ev(p1) ! kcal/mol/rad/rad
           p4=p4
           p2=p2             ! deg
          case('K','K/rad/rad')
           p0=p0
           p1=K2eV(p1)
           p2=p2
           p4=p4
         end select
         write(string,*)   p0,p1,p4,p2,' # ',type_fff
        case("charmm")
         read(variables,*) p1,p0,p4,p2
         select case(units)
          case("kcal/mol")
           p1=kcalmol2ev(p1)
           p0=p0
           p4=p4
           p2=p2
         end select
         write(string,'(f14.7,1x,i3,1x,i3,1x,f14.7,a,a)') p1,p0,p4,p2,' # ',type_fff
        case default
         write(6,'(a,1x,a)')'Unknown dihedral potential:', adjustl(trim(type_fff))
         stop
        end select
       else
        write(string,*)1,0.0,1,0.0,' # ',type_fff
       end if
       exit fff
      end if
     end if
    end do
   case('impr')
    initimpr: do
     read(u,'(a)') line
     if(line(1:15)=='Improper Coeffs') then
      read(line(16:),'(a)') units
      units=adjustl(trim(units))
      exit initimpr
     end if
    end do initimpr
    do
     read(u,'(a)',iostat=ierr) line
     if (ierr/=0) exit fff
     if(line(1:11)=='Pair Coeffs'.or.line(1:3)=="End") then
      string(1:4)="none"
      exit fff
     else if(line(1:1)/="#".or.line(1:1)/=" ")then
      read(line,'(4(a4,1x),a)')(ourlabel(ii),ii=1,4),passing
      if((adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&  ! i j k l
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(2))).and.&  ! 1 2 3 4
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(3))).and.&
          adjustl(trim(ourlabel(4)))==adjustl(trim(label(4)))).or.&
         (adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&  ! i j l k
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(2))).and.&  ! 1 2 4 3
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(4))).and.&
          adjustl(trim(ourlabel(4)))==adjustl(trim(label(3)))).or.&
         (adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&  ! i k j l
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(3))).and.&  ! 1 3 2 4
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(4)))==adjustl(trim(label(4)))).or.&
         (adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&  ! i k l j
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(3))).and.&  ! 1 3 4 2
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(4))).and.&
          adjustl(trim(ourlabel(4)))==adjustl(trim(label(2)))).or.&
         (adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&  ! i l j k
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(4))).and.&  ! 1 4 2 3
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(4)))==adjustl(trim(label(3)))).or.&
         (adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&  ! i l k j
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(4))).and.&  ! 1 4 3 2
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(3))).and.&
          adjustl(trim(ourlabel(4)))==adjustl(trim(label(2)))) )then
       read(passing,'(a,a)')type_fff,variables
       if(adjustl(trim(type_fff))/='none')then
        select case(adjustl(trim(type_fff)))
        case("cossq")        ! 0.5*p1*cos(x-p2)^2
         read(variables,*) p1,p2
         select case(units)
          case('kcal/mol')
           p1=kcalmol2ev(p1) ! kcal/mol/rad/rad
           p2=p2             ! [-]
         end select
         write(string,*) p1,p2,' # ',type_fff
        case("cvff")         ! p1*(1+p0*cos(p4*x))
         read(variables,*) p1,p0,p4
         select case(units)
          case('kcal/mol')
           p1=kcalmol2ev(p1) ! kcal/mol/rad/rad
           p0=p0             ! [-]
           p4=p4             ! [-]
          case('K','K/rad/rad')
           p1=K2eV(p1)
           p0=p0
           p4=p4
         end select
         write(string,'(f14.7,1x,i3,1x,i3,a,a)') p1,p0,p4,' # ',type_fff
        case default
         write(6,'(a,1x,a)')'Unknown improper potential:', adjustl(trim(type_fff))
         stop
        end select
       else
        write(string,*) 0.0,0.0,' # ',type_fff
       end if
       exit fff
      end if
     end if
    end do
   case('pair')
    initpair: do
     read(u,'(a)') line
     if(line(1:11)=='Pair Coeffs') then
      read(line(12:),'(a)') units
       units=adjustl(trim(units))
      read(u,'(a)') line
      if(line(1:3)=="mix")then
       read(line(12:),'(a)') mixing_rule
       mixing_rule = adjustl(trim(mixing_rule))
      end if 
      exit initpair
     end if
    end do initpair
    do
     read(u,'(a)',iostat=ierr) line
     if(ierr/=0) exit fff
     if(line(1:3)=="End") then
      string(1:4)="none"
      exit fff
     else if(line(1:1)/="#".or.line(1:1)/=" ")then
      read(line,'(2(a4,1x),a)')ourlabel(1),ourlabel(2),passing
      if((adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(2)))).or.&
         (adjustl(trim(ourlabel(1)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(1)))))then
       read(passing,'(a,a)')type_fff,variables
       if(adjustl(trim(type_fff))/='none')then
        read(variables,*) p1,p2
        select case(units)
         case('kcal/mol')
          p1=kcalmol2ev(p1) ! kcal/mol
          p2=p2             ! A
        end select
        write(string,'(f14.7,1x,f14.7,a,a)') p1,p2,' # ',type_fff
       else
        write(string,'(f14.7,1x,f14.7,a,a)') 0.0,1.0,' # ',type_fff
       end if
       exit fff
      end if
     end if
    end do
  end select fff
  close(u)
  return
 end subroutine search_forcefield
 subroutine output_CIF()
 implicit none
 character(len=120) :: CIFFilenameNew
 integer           :: u=1000
 integer           :: i
 character(len=12)  :: extension="_topol.cif"
 CIFFilenameNew=filename(1:Clen_trim(filename))//extension
 CIFFilenameNew=adjustl(CIFFilenameNew)
 open(u,file=CIFFilenameNew)
 write(u,'(a)')'data_subtitutions'
 write(u,'(a)')'_audit_creation_method    igor'
 write(u,'(a)')"_audit_author_name 'sponge bob'"
 write(u,'(a,5x,f14.7)')'_cell_length_a',cell_0(1)
 write(u,'(a,5x,f14.7)')'_cell_length_b',cell_0(2)
 write(u,'(a,5x,f14.7)')'_cell_length_c',cell_0(3)
 write(u,'(a,5x,f14.7)')'_cell_angle_alpha',cell_0(4)
 write(u,'(a,5x,f14.7)')'_cell_angle_beta',cell_0(5)
 write(u,'(a,5x,f14.7)')'_cell_angle_gamma',cell_0(6)
 write(u,'(a,5x,f14.7)')'_cell_volume',volume(rv)
 write(u,'(a)')"_symmetry_space_group_name_hall 'p 1'"
 write(u,'(a)')"_symmetry_space_group_name_h-m 'p 1'"
 write(u,'(a)')'_symmetry_int_tables_number 1'
 write(u,'(a)')"_symmetry_equiv_pos_as_xyz 'x,y,z'"
 write(u,'(a)')'loop_'
 write(u,'(a)')'_atom_site_label'
 write(u,'(a)')'_atom_site_fract_x'
 write(u,'(a)')'_atom_site_fract_y'
 write(u,'(a)')'_atom_site_fract_z'
 write(u,'(a)')'_atom_site_charge'
  do i=1,n_atoms
   write(u,'(a4,1x,3(f14.7,1x),1x,f14.7)')atom(i)%label,(atom(i)%xyzs(j,1),j=1,3),atom(i)%charge
  end do
 close(u)
 end subroutine output_CIF
 subroutine output_pdb()
  implicit none
  integer           :: u=333,i
  character(len=4)  :: extension=".pdb"
  character(len=100) :: PDBFilename
  PDBFilename=filename(1:Clen_trim(filename))//extension
  PDBFilename=adjustl(PDBfilename)
  open(u,file=PDBFilename)
  write(u,'(a5,i5)')'MODEL',k
  write(u,'(a6,9(f14.7,1x))')'REMARK',&
   xlo_bound,ylo_bound,zlo_bound,xhi_bound,yhi_bound,zhi_bound,xy,xz,yz
  write(u,'(a6,1x,a,1x,i6,1x,a,i6)')'REMARK','# atoms',n_atoms,'# timestep',int(0.0)
  write(u,'(a6,1x,a,1x,f14.7)')'REMARK','Volume:',volume(rv)
  write(u,'(a6,3f9.3,3f7.2)')'CRYST1',(cell_0(j),j=1,6)
  do i=1,n_atoms
   write(u,'(a6,i5,1x,a2,3x,a4,1x,i4,4x,3f8.3,2f6.2,9x,a2)') &
   'ATOM  ',i,atom(i)%label,'MOL ',0,(atom(i)%xyzc(j,1),j=1,3),0.0,0.0,atom(i)%label
  end do
  write(u,'(a6)')'ENDMDL'
  close(u)
 end subroutine output_pdb
 subroutine output_gulp()
  implicit none
  character(len=100) :: GULPFilename
  integer           :: u=444
  integer           :: i
  real              :: mmm,rrr
  integer           :: zzz
  character(len=2)  :: zlz
  character(len=6)  :: atomtype_library
  character(len=4)  :: extension=".gin"
  GULPFilename=filename(1:Clen_trim(filename))//extension
  GULPFilename=adjustl(GULPfilename)
  !adjustl(
  open(u,file=GULPFilename)
  write(u,'(a)')'single conp'
  write(u,'(A)')'cell'
  write(u,'(6(f9.5,1x))') (cell_0(j) , j=1,6)
  write(u,'(A)')'fractional'
  do i=1,n_atoms
   write(u,'(a4,1x,3(f14.7,1x),1x,f14.7)')atom(i)%label,(atom(i)%xyzs(j,1),j=1,3),atom(i)%charge
  end do
   write(u,'(a,1x,i3)')'species',n_atom_types
  do i=1,n_atom_types
   call checkatomtype_gulp(atom_types(i), atomtype_library)
   select case(atom_types(i))
    case("H0  ":"H999")
     write(u,'(a4,1x,a4,1x,a6)') "H_  ", "core", atomtype_library
    case default
     write(u,'(a4,1x,a4,1x,a6)') atom_types(i), "core", atomtype_library
   end select
  end do
  write(u,'(a)') "library GenericZIF"
  write(u,'(a)') "#switch_minimiser rfo gnorm 0.05"
  write(u,'(a)') "stepmx opt 0.1"
  write(u,'(a)') 'dump every 1 optimise.grs'
  close(u)
 end subroutine output_gulp
!
 subroutine output_lammps()
  implicit none
  integer :: day,hour,i4_huge,milli,minute,month,second,year
  integer :: ii,jj
  character(len=10) :: time
  character(len=8)  :: date
  integer           :: u=222
  real              :: mmm,rrr
  integer           :: zzz
  character(len=2)  :: zlz
  character(len=120) :: DataFilename
  character(len=5)  :: extension=".data"
  character(len=4)  :: label(4)
  character(len=80) :: str
  DataFilename=filename(1:Clen_trim(filename))//extension
  DataFilename=adjustl(trim(DataFilename))
  forall (ii=1:4)
   forall (jj=1:4)
    label(ii)(jj:jj)=" "
   end forall
  end forall
  open(u,file=DataFilename)
  call date_and_time (date,time)
  read (date,'(i4,i2,i2)')year,month,day
  read (time,'(i2,i2,i2,1x,i3)')hour,minute,second,milli
  write(u,'(a,1x,i2,1x,i2,1x,i4,1x,i2,a1,i2,a1,i2,a)')&
   'Created on',day,month,year,hour,':',minute,':',second,' using cif2lammps'
  write(u,*)' '
  write(u,*)n_atoms,' atoms'
  write(u,*)n_bonds,' bonds'
  write(u,*)n_bends,' angles'
  write(u,*)n_torss,' dihedrals'
  write(u,*)n_imprs,' impropers'
  write(u,*)' '
  write(u,*)n_atom_types,' atom types'
  write(u,*)n_bond_types,' bond types'
  write(u,*)n_bend_types,' angle types'
  write(u,*)n_tors_types,' dihedral types'
  write(u,*)n_impr_types,' improper types'
  write(u,*)' '
  write(u,*)xlo_bound,xhi_bound,' xlo xhi'
  write(u,*)ylo_bound,yhi_bound,' ylo yhi'
  write(u,*)zlo_bound,zhi_bound,' zlo zhi'
  write(u,*)xy,xz,yz,' xy xz yz'
  write(u,*)' '
  write(u,'(a)')'Masses'
  write(u,*)' '
  do i=1,n_atom_types
   call CheckAtom(atom_types(i),mmm,rrr,zzz,zlz)
   write(u,'(i4,1x,f14.7,1x,a,a)')i,mmm,' # ',atom_types(i)
  end do
  write(u,*)' '
  write(u,'(a)')'Bond Coeffs'
  write(u,*)' '
  do i=1,n_bond_types
   read(bond_type_string(i),'(2a4)')(label(ii),ii=1,2)
   forall (ii=1:80)
    str(ii:ii)=" "
   end forall
   call search_forcefield(str,'bond',label)
   write(u,'(i4,3x,a,a,a)')i,adjustl(str(1:Clen_trim(str))),' # ',bond_type_string(i)
  end do
  write(u,*)' '
  write(u,'(a)')'Angle Coeffs'
  write(u,*)' '
  do i=1,n_bend_types
   read(bend_type_string(i),'(3a4)')(label(ii),ii=1,3)
   forall (ii=1:80)
    str(ii:ii)=" "
   end forall
   call search_forcefield(str,'bend',label)
   write(u,'(i4,3x,a,a,a)')i,adjustl(str(1:Clen_trim(str))),' # ',bend_type_string(i)
  end do
  write(u,*)' '
  write(u,'(a)')'Dihedral Coeffs'
  write(u,*)' '
  do i=1,n_tors_types
   read(tors_type_string(i),'(4a4)')(label(ii),ii=1,4)
   forall (ii=1:80)
    str(ii:ii)=" "
   end forall
   call search_forcefield(str,'tors',label)
   write(u,'(i4,3x,a,a,a)')i,adjustl(str(1:Clen_trim(str))),' # ',tors_type_string(i)
  end do
  write(u,*)' '
  write(u,'(a)')'Improper Coeffs'
  write(u,*)' '
  do i=1,n_impr_types
   read(impr_type_string(i),'(4a4)')(label(ii),ii=1,4)
   forall (ii=1:80)
    str(ii:ii)=" "
   end forall
   call search_forcefield(str,'impr',label)
   write(u,'(i4,3x,a,a,a)')i,adjustl(str(1:Clen_trim(str))),' # ',impr_type_string(i)
  end do
  !write(u,*)' '
  !write(u,'(a)')'Pair Coeffs'
  !write(u,*)' '
  open(654,file="pair.txt")
  do i=1,n_atom_types
   label(1)=atom_types(i)
   label(2)=atom_types(i)
   forall (ii=1:80)
    str(ii:ii)=" "
   end forall
   call search_forcefield(str,'pair',label)
   !write(654,'(i4,i4,1x,a,a3,a4,1x,a4)')i,i,str,' # ', (label(ii),ii=1,2) 
   write(654,'(a4,1x,a,a,1x,i4)')label(1),str,'#',i
  end do
  close(654)
  write(u,*)' '
  write(u,'(a)')'Atoms'
  write(u,*)' '
  do i=1,n_atoms
   write(u,'(i8,1x,i3,1x,i3,1x,f10.7,1x,3(f14.7,1x),a3,a4,1x,a10)') &
    i,1,atom(i)%type_,atom(i)%charge,(atom(i)%xyzc(j,1),j=1,3),&
    ' # ',atom(i)%new_label,atom(i)%hybridization
  end do
  write(u,*)' '
  write(u,'(a)')'Bonds'
  write(u,*)' '
  do i=1,n_bonds
   write(u,'(i8,a)')i,bond(i)
  end do
  write(u,*)' '
  write(u,'(a)')'Angles'
  write(u,*)' '
  do i=1,n_bends
   write(u,'(i8,a)')i,bend(i)
  end do
  write(u,*)' '
  write(u,'(a)')'Dihedrals'
  write(u,*)' '
  do i=1,n_torss
   write(u,'(i8,a)')i,tors(i)
  end do
  write(u,*)' '
  write(u,'(a)')'Impropers'
  write(u,*)' '
  do i=1,n_imprs
   write(u,'(i8,a)')i,impr(i)
  end do
  write(u,*)' '
  close(u)
 end subroutine output_lammps
 subroutine cellnormal2lammps(cell_0,xlo_bound,ylo_bound,zlo_bound,&
                              xhi_bound,yhi_bound,zhi_bound,xy,xz,yz)
  implicit none
  real,intent(in)  :: cell_0(6)
  real,intent(out) :: xlo_bound,ylo_bound,zlo_bound,xhi_bound,yhi_bound,zhi_bound
  real,intent(out) :: xy,xz,yz
  real,parameter :: pi=acos(-1.0)
  real,parameter :: radtodeg = 180.0/pi
  real,parameter :: degtorad = pi/180.0
  real :: xlo,xhi,ylo,yhi,zlo,zhi,lx,ly,lz,cosa,cosb,cosg
  real :: alp,bet,gam,sing
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
  lx=cell_0(1)
  xy=cell_0(2)*cosg
  xz=cell_0(3)*cosb
  ly=sqrt(cell_0(2)*cell_0(2)-xy*xy)
  yz=(cell_0(2)*cell_0(3)*cosa-xy*xz)/ly
  lz=sqrt(cell_0(3)*cell_0(3)-xz*xz-yz*yz)
  xlo_bound=0.0
  ylo_bound=0.0
  zlo_bound=0.0
  xhi_bound=lx
  yhi_bound=ly
  zhi_bound=lz
  return
 end subroutine cellnormal2lammps

 subroutine checkatomtype_gulp(label,newlabel)
  implicit none
  character(len=4),intent(in)  :: label
  character(len=6),intent(out) :: newlabel
  select case(label)
   case("H1  ":"H999")
    newlabel = "H_    "
   case("Zn  ")
    newlabel = "Zn_cor"
   case("He  ")
    newlabel = "Zn_she"
   case("C1  ","C4  ")
    newlabel = "C_R2  "
   case("C2  ","C6  ")
    newlabel = "C_R   "
   case("C5  ")
    newlabel = "C_R3  "
   case("C3  ","C8  ")
    newlabel = "C_3   "
   case("C7  ")
    newlabel = "C_32  "
   case("N1  ")
    newlabel = "N1_cor"
   case("Ne  ")
    newlabel = "N1_she"
   case default
    write(6,'(a1,a4,a1)')"'",label,"'"
    write(6,'(a)')"Atom unknowed"
    newlabel = label
  end select
  return
 end subroutine checkatomtype_gulp
!
 subroutine checkatom(Label,m,s,Z,Zlabel)
  implicit none
  character(len=4),intent(in)  :: Label
  real,intent(out)             :: m,s
  real,parameter               :: dummy_mass = 1.00794 ! (hydrogen atom mass)
  integer,intent(out)          :: Z
  character(len=2),intent(out) :: ZLabel
  select case(Label)
   case('O   ','O0  ':'O999')
    Z=8
    m=15.999  ! conventional weight
    s=0.66    ! covalent radious
    Zlabel=' O'
   case('C   ','C0  ':'C999')
    Z=6
    m=12.0107
    s=0.720
    ZLabel=' C'
   case('H   ','H0  ':'H999')
    Z=1
    m=dummy_mass
    s=0.320
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
    s = 1.0
    Zlabel = 'Cl'
! dummy atoms or fake particles (pseudo-atoms):
   case('He  ','He0 ':'He99') ! shell of the Zn atom
    Z = 2
    m = dummy_mass
    s = (1.470 + 0.1)/2.0
    Zlabel = 'He'
   case('Ne  ','Ne0 ':'Ne99') ! shell of the N atom
    Z = 10
    m = dummy_mass
    s = 0.55
    Zlabel='Ne'
   case('Xe  ','Xe0 ':'Xe99') ! negative particle of the dipole simulating a DMF molecule
    Z = 54
    m = 15.999
    s = 0.505
    Zlabel='Xe'
   case('Kr  ','Kr0 ':'Kr99') ! positive particle ...
    Z = 36
    m = 15.999
    s = 0.505
    Zlabel='Kr'
   case default
    write(6,'(a1,a4,a1)')"'",label,"'"
    STOP 'Atom unknowed'
  end select
 end subroutine checkatom
!
 PURE INTEGER FUNCTION Clen(s)      ! returns same result as LEN unless:
 CHARACTER(*),INTENT(IN) :: s       ! last non-blank char is null
 INTEGER :: i
 Clen = LEN(s)
 i = LEN_TRIM(s)
 IF (s(i:i) == CHAR(0)) Clen = i-1  ! len of C string
 END FUNCTION Clen
!
 PURE INTEGER FUNCTION Clen_trim(s) ! returns same result as LEN_TRIM unless:
 CHARACTER(*),INTENT(IN) :: s       ! last char non-blank is null, if true:
 INTEGER :: i                       ! then len of C string is returned, note:
                                    ! Ctrim is only user of this function
 i = LEN_TRIM(s) ; Clen_trim = i
 IF (s(i:i) == CHAR(0)) Clen_trim = Clen(s)   ! len of C string
 END FUNCTION Clen_trim
!
 SUBROUTINE cell(rv,vr,cell_0)
 implicit none
 integer :: i,j
 real, intent(in)  :: cell_0(6)
 real, intent(out) :: rv(3,3),vr(3,3)
 real, parameter   :: pi = ACOS(-1.0)
 real :: alp,bet
 real :: cosa,cosb,cosg
 real :: gam,sing
 real :: DEGTORAD
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
 WRITE(6,'(a)') 'Cell:'
 WRITE(6,'(6F14.7)')( cell_0(j), j=1,6 )
 WRITE(6,'(a)')'Linear Transformation Operator:'
 DO i=1,3
  WRITE(6,'(F14.7,F14.7,F14.7)')( rv(i,j), j=1,3 )
 ENDDO
 WRITE(6,'(a)')'----------------------------------------'
 WRITE(6,'(a)')'Inverse Linear Transformation Operator:'
 DO i=1,3
  WRITE(6,'(F14.7,F14.7,F14.7)')( vr(i,j), j=1,3 )
 ENDDO
 WRITE(6,'(a)')'----------------------------------------'
 RETURN
 END SUBROUTINE cell
!
 SUBROUTINE uncell(rv,cell_0)
  implicit none
  real,intent(out)   :: cell_0(6)
  real,intent(in)    :: rv(3,3)
  integer            :: i,j
  real               :: temp(6)
  REAL               :: radtodeg
  REAL, PARAMETER    :: pi=ACOS(-1.0)
  radtodeg=180.0/PI
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
  DO i=4,6
     if (abs(cell_0(i) - 90.0 ).lt.0.00001) cell_0(i) = 90.0
     if (abs(cell_0(i) - 120.0).lt.0.00001) cell_0(i) = 120.0
  ENDDO
  RETURN
 END SUBROUTINE uncell
!
 SUBROUTINE inverse(a,c,n)
 implicit none
 integer n
 real a(n,n), c(n,n)
 real L(n,n), U(n,n), b(n), d(n), x(n)
 real coeff
 integer i, j, k
 L=0.0
 U=0.0
 b=0.0
 do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
 end do
 do i=1,n
  L(i,i) = 1.0
 end do
 do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
 end do
 do k=1,n
  b(k)=1.0
  d(1) = b(1)
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
 end do
 RETURN
 END SUBROUTINE inverse
!
subroutine make_dist_matrix(n,cell_0,rv,vr,x,dist_matrix)
 implicit none
 integer,intent(in) :: n
 real,intent(in)    :: cell_0(6),rv(3,3),vr(3,3),x(3,n)
 real,intent(out)   :: dist_matrix(n,n)
 integer            :: i,j,k
 real               :: r1(3),r2(3),s
 DO i=1,n
    dist_matrix(i,i)=0.0
    DO j=i+1,n
       forall ( k=1:3 )
        r1(k)=x(k,i)
        r2(k)=x(k,j)
       end forall
       call make_distances(cell_0,r1,r2,rv,s)
       dist_matrix(i,j)=s
       dist_matrix(j,i)=dist_matrix(i,j)
    END DO
 END DO
 return
end subroutine make_dist_matrix
!
 SUBROUTINE make_distances(cell_0,r2,r1,rv,dist)
 IMPLICIT NONE
 REAL,    intent(in)  :: r1(3),r2(3),rv(3,3),cell_0(6)    ! coordenadas y matriz de cambio
 REAL,    intent(out) :: dist
 REAL                 :: d_image(1:27),image(3,27)        ! array de distancias
 INTEGER              :: k,l,m,n,o,i,j                    ! variables mudas
 REAL                 :: atom(3),ouratom(3)               ! coordenadas preparadas
  k=0
  do concurrent (l=-1:1)
   do concurrent (m=-1:1)
      do concurrent (n=-1:1)
         k = k + 1
         ouratom(1) = r1(1)
         ouratom(2) = r1(2)
         ouratom(3) = r1(3)
         atom(1) = r2(1) + l
         atom(2) = r2(2) + m
         atom(3) = r2(3) + n
         d_image(k) = distance(atom,ouratom,rv)
         forall ( i=1:3)
           image(i,k) = atom(i)
         end forall
     enddo
   enddo
  enddo
  dist=MINVAL(d_image)
  RETURN
 END SUBROUTINE
!
pure real function volume(rv)
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
!
 pure REAL FUNCTION DISTANCE(atom,ouratom,rv)
  IMPLICIT NONE
  INTEGER :: j
  real,intent(in) :: atom(3), ouratom(3), rv(3,3)
  REAL            :: dist(3),o_atom(3),o_ouratom(3)
  FORALL ( j=1:3 )
   o_ouratom(j) = rv(j,1)*ouratom(1)  + rv(j,2)*ouratom(2)  + rv(j,3)*ouratom(3)
   o_atom(j)    = rv(j,1)*atom(1) + rv(j,2)*atom(2) + rv(j,3)*atom(3)
   dist(j) = o_ouratom(j) - o_atom(j)
  END FORALL
  DISTANCE = sqrt(dist(1)*dist(1) + dist(2)*dist(2) + dist(3)*dist(3))
 END FUNCTION
!
 pure real function K2eV(x)
  implicit none
  real,intent(in)    :: x
  real,parameter :: kB_1=11604.52211052
  K2eV=x/kB_1
  return
 end function
!
 pure real function kcalmol2ev(x)
  implicit none
  real,intent(in)    :: x
  kcalmol2ev = 0.04336*x
  return
 end function  kcalmol2ev
!
 real function kjmol2ev(x)
  implicit none
  real,intent(in) :: x
  kjmol2ev=x*(1036.427/100000.0)
  return
 end  function kjmol2ev
!
 real function kjmolnmnm2evAA(x)
  real, intent(in) :: x
  kjmolnmnm2evaa = kjmol2ev(x)/100.0
  return
 end function  kjmolnmnm2evaa
!
 subroutine print_help()
    print '(a)', '  -h, --help   print usage information and exit'
    print '(a)', '  -c, --cif    CIF File input'
    print '(a)', '  -GCMC, --adsorption [1e5] pressure Pa'
    print '(a)', '  -T, --temperatue  -T, [298] temperature in Kelvin'
    print '(a)', '  -fff, --forcefield [WHCJ] WHCJ, HZJ'
 end subroutine print_help
end program zif_cif2gin
