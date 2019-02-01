module PhononUnfoldingModule
implicit none

!----------------------------------
! type define
!----------------------------------
    real(kind=8) :: pi=3.141592653589793
    complex(kind=8) :: I_imag=(0,1)
    integer :: len_characters=100
    integer :: readline_length=300

    type array_car
        integer :: N1, N2, N3
        complex(kind=8), allocatable :: c_1d(:)
        complex(kind=8), allocatable :: c_2d(:,:)
        complex(kind=8), allocatable :: c_3d(:,:,:)
        real(kind=8), allocatable :: r_1d(:)
        real(kind=8), allocatable :: r_2d(:,:)
        real(kind=8), allocatable :: r_3d(:,:,:)
        integer,allocatable :: i_1d(:)
        integer,allocatable :: i_2d(:)
        integer,allocatable :: i_3d(:)
    end type array_car

    type pebu  ! unfolded energy bands
        integer :: n_k ! number of k-points
        integer :: n_b ! number of energy bands for each k-point
        integer :: n_atom_sc ! number of atoms in a supercell
        real(kind=8),allocatable :: r_atom_sc(:,:)
        real(kind=8),allocatable :: r_direct_atom_sc(:,:)
        real(kind=8) :: s_lat(3,3) ! supercell cell
        real(kind=8) :: p_lat(3,3) ! primary cell
        integer :: rpc(3,3)    ! reciprocal primarycell lattice in the direct form of the supercell reciprocal lattice 
        real(kind=8), allocatable :: k(:,:) ! kpoint list, dimension : n_k*3       
        real(kind=8), allocatable :: k_s(:,:) ! corresponding kpoint list in the supercell reciprocal lattice
        real(kind=8), allocatable :: e(:,:,:) ! dimension: n_k*n_b*2, e(i,j,1) is the energy, while e(i,j,2) is the weight
        complex(kind=8), allocatable :: c(:,:,:,:) ! n_k*n_b*n_atom_sc*3 eigen wavefunction coefficient for each wannier function basis

        real(kind=8), allocatable :: klength(:) ! kpoint length for figure x-axis , dimension : n_k
        integer :: n_s             ! number of segments
        real(kind=8), allocatable :: boundary_s(:)   ! segments boundary, dimension: n_s
    end type pebu 



!--------------------------------------------------------
contains
!--subroutines and functions-----------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! lower case 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
character(len=readline_length) function lower_case(str)
character(len=*) :: str
integer :: i, A, Z, a_lower, dif

A=ichar('A')
Z=ichar('Z')
a_lower=ichar('a')
lower_case=str
dif=a_lower-A
do i=1, readline_length
   if( ichar(lower_case(i:i)).ge.A .and. ichar(lower_case(i:i)).le.Z ) lower_case(i:i)=char(ichar(lower_case(i:i))+dif)
enddo

endfunction lower_case


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  function: delete space  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
character(len=readline_length) function delete_space(str)
character(len=*) :: str
integer :: i, flag
delete_space=''
flag=0
do i=1,len_trim(str)
     if(ichar(str(i:i)) .ne. ichar(' ')) then
         flag=flag+1
         delete_space(flag:flag)=str(i:i)
     endif
enddo

endfunction delete_space

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! read integer number from 'file'  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function read_integer(file,name, default_value)
character(len=*):: file
character(len=*):: name
integer :: default_value
integer :: status_file, status_name, i,j, flag
character(len=readline_length) :: temp
flag=0
open(1,file=trim(file),status='old',action='read',iostat=status_file)
if(status_file==0) then
    readloop: do
        read(1,'(a)',iostat=status_name) temp
        if(status_name /=0) exit
        temp=trim(lower_case(temp))
        i=index(temp,'#')
        if(i/=0)temp=temp(1:i-1)
        i=index(temp,'!')
        if(i/=0)temp=temp(1:i-1)
        i=index(temp,trim(name))
        if(i/=0) then
            temp=temp(i+len_trim(name):)
            j=index(temp,'=')
            if(j/=0) temp=temp(j+1:)
            read(temp,*)read_integer
            flag=flag+1
        endif
    end do readloop
else
    open(2,file='check',access='append')
    write(2,*)'err in opening file '//trim(file)
    close(2)
endif


if(flag==0) read_integer=default_value
if(flag>1) then
    open(2,file='check',access='append')
    write(2,*)'More than one lines contain '//trim(name)//' in '//trim(file)
    close(2)
endif
close(1)
endfunction read_integer


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! read real number from 'file'  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind=8) function read_real(file,name, default_value)
character(len=*):: file
character(len=*):: name
real(kind=8) :: default_value
integer :: status_file, status_name, i,j, flag
character(len=readline_length) :: temp

flag=0

open(1,file=trim(file),status='old',action='read',iostat=status_file)
if(status_file==0) then
    readloop: do
        read(1,'(a)',iostat=status_name) temp
        if(status_name /=0) exit
        temp=trim(lower_case(temp))
        i=index(temp,'#')
        if(i/=0)temp=temp(1:i-1)
        i=index(temp,'!')
        if(i/=0)temp=temp(1:i-1)
        i=index(temp,trim(name))
        if(i/=0) then
            temp=temp(i+len_trim(name):)
            j=index(temp,'=')
            if(j/=0) temp=temp(j+1:)
            read(temp,*)read_real
            flag=flag+1
        endif
    end do readloop
else
    open(2,file='check',access='append')
    write(2,*)'err in opening file '//trim(file)
    close(2)
endif

if(flag==0) read_real=default_value
if(flag>1) then
    open(2,file='check',access='append')
    write(2,*)'More than one lines contain '//trim(name)//' in '//trim(file)
    close(2)
endif


close(1)
endfunction read_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! read characters from 'file'  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
character(len=len_characters) function read_characters(file,name, default_value)
character(len=*):: file
character(len=*):: name
character(len=200):: name_lower
character(len=*) :: default_value
integer :: status_file, status_name, i,j, flag
character(len=readline_length) :: temp

flag=0
open(1,file=trim(file),status='old',action='read',iostat=status_file)

if(status_file==0) then
    readloop: do
        read(1,'(a)',iostat=status_name) temp
        if(status_name /=0) exit
        i=index(temp,'#')
        if(i/=0)temp=temp(1:i-1)
        i=index(temp,'!')
        if(i/=0)temp=temp(1:i-1)
        name_lower=trim(lower_case(temp))

        i=index(name_lower,trim(name))
        if(i/=0) then
            temp=temp(i+len_trim(name):)
            j=index(temp,'=')
            if(j/=0) temp=temp(j+1:)
            temp=delete_space(temp)
            read(temp,*)read_characters
            flag=flag+1
        endif
    end do readloop
else
    open(2,file='check',access='append')
      write(2,*)'err in opening file '//trim(file)
    close(2)
endif

if(flag==0) read_characters=default_value
if(flag>1) then
    open(2,file='check',access='append')
    write(2,*)'More than one lines contain '//trim(name)//' in '//trim(file)
    close(2)
endif

close(1)
endfunction read_characters



!!!!!!!!!!!! construct kpoint list in the first Brillouin zone of primary cell and the corresponding
!!!!!!!!!!!! kpoint list in the first Brillouin zone of supercell
subroutine construct_qlist(e)
     implicit none
     type(pebu) :: e
     type(array_car) :: out
     type(array_car) :: cell1_read, cell2_read
     !integer :: n_segment
     integer(kind=8), allocatable :: n_k(:)
     real(kind=8), allocatable :: a(:,:), b(:,:)
     real(kind=8) :: t(3,3),tt(3,3), q1(3), q2(3)
     integer :: i, j, flag
     integer :: ipiv(3)
    !-------------------------------for dgetrf
    integer ::  info,n, lda
    integer, allocatable :: ipiv2(:)
    integer :: lwork
    real(kind=8), allocatable :: work(:)
    !-------------------------------
    
 
     call read_qpoints_setting_primary_cell('input.dat', out)       ! read input.dat 
     e%n_s=out%N1
     
     allocate(n_k(e%n_s), a(e%n_s,3), b(e%n_s,3))
     n_k=out%i_1d
     do i=1, e%n_s
          a(i,:)=out%r_2d(2*i-1,:)
          b(i,:)=out%r_2d(2*i,:)
     enddo
     ! read the primary cell lattice parameters
     call read_primary_cell('input.dat',cell1_read)
     e%p_lat(1,:)=cell1_read%r_2d(1,:)
     e%p_lat(2,:)=cell1_read%r_2d(2,:)
     e%p_lat(3,:)=cell1_read%r_2d(3,:)
     ! read the super cell lattice parameters
     call read_super_cell('input.dat',cell2_read)
     e%s_lat(1,:)=cell2_read%r_2d(1,:)
     e%s_lat(2,:)=cell2_read%r_2d(2,:)
     e%s_lat(3,:)=cell2_read%r_2d(3,:)

    !allocate e%k(:,:)
     e%n_k=sum(n_k)
     allocate(e%k(e%n_k,3), e%k_s(e%n_k,3),e%klength(e%n_k), e%boundary_s(e%n_s+1))
     e%klength=0d0
     e%boundary_s=0d0
!!!!!!!
! calculate inverse primary cell lattice     
!!!!!!
     t=dble(e%p_lat)         ! prepare inverse matrix
!     call getrf(t,ipiv)  ! prepare inverse matrix
     n=3; lda=3
     allocate(ipiv2(2*n))
     call dgetrf(n,n,t,lda,ipiv2,info)

!    call getri(t,ipiv)  ! inverse matrix
     lwork=2*n
     if(.not.allocated(work)) allocate(work(lwork))
     call dgetri(n, t, lda, ipiv2, work, lwork, info)

     


    !!!!!!!!!!
    ! construct kpoint list in primary cell Brillouin zone
    !!!!!!!!!!
     flag=0
     do i=1,e%n_s
         if(flag.lt.0.5) then
            do j=1,n_k(i)
                 flag=flag+1
                   e%k(flag,:)= a(i,:)*(n_k(i)-j)/(n_k(i)-1d0)    +    b(i,:)*(j-1d0)/(n_k(i)-1d0)
            enddo
         elseif( sum(abs(a(i,:)-b(i-1,:))) .lt. 0.0000001 )then
            do j=1,n_k(i)
                 flag=flag+1
                   e%k(flag,:)= a(i,:)*(n_k(i)-j)/(n_k(i)-0d0)    +    b(i,:)*(j-0d0)/(n_k(i)-0d0)
            enddo
         else
            do j=1,n_k(i)
                 flag=flag+1
                   e%k(flag,:)= a(i,:)*(n_k(i)-j)/(n_k(i)-1d0)    +    b(i,:)*(j-1d0)/(n_k(i)-1d0)
            enddo
         endif
     enddo
    ! calculate e%klength for nice figure x-axis
         tt=transpose(t)
         q2=e%k(1,1)*tt(1,:)+e%k(1,2)*tt(2,:)+e%k(1,3)*tt(3,:)
         
     do i=2, e%n_k+1
       q1=q2
       q2=e%k(i,1)*tt(1,:)+e%k(i,2)*tt(2,:)+e%k(i,3)*tt(3,:)
       !write(*,*)q1
       !write(*,*)q2
       !pause
        e%klength(i)=e%klength(i-1)+sqrt((q2(1)-q1(1))**2 + (q2(2)-q1(2))**2 + (q2(3)-q1(3))**2)
     enddo
     do i=1,e%n_s
        e%boundary_s(i+1)=e%klength(sum(n_k(1:i)))
     enddo

     ! construct reciprocal primarycell lattice in the direct form of reciprocal supercell lattice
 ! call DLINRG(3,e%p_lat,3,t,3)
 


     t= matmul(transpose(t), transpose(e%s_lat) )
     e%rpc=anint(t)

     ! make the corresponding kpoint list in the reciprocal supercell cell lattice
     do i=1,e%n_k
         ! e%k_s(i,:)=matmul(e%rpc,e%k(i,:))
          e%k_s(i,:)=matmul(e%k(i,:),e%rpc)
         ! write(*,*)e%k(i,:)
         ! write(*,*)
         ! write(*,*)e%rpc(1,:)
         ! write(*,*)e%rpc(2,:)
         ! write(*,*)e%rpc(3,:)
         ! write(*,*)
         ! write(*,*)e%k_s(i,:)
         ! pause
     enddo

end subroutine







subroutine write_qlist_to_file(e)
     implicit none
     type(pebu) :: e
     integer :: i, j, k
     character(len=20) :: t
     character(len=20) :: flag
     
 flag=trim(read_characters('input.dat','write_q_correspondence','false'))
 
 if(trim(flag).eq.'true'.or.trim(flag).eq.'t' ) then
   open(11,file='Q-points.dat')
   do i=1, e%n_k
          write(11,'(7f15.10)') e%klength(i), e%k(i,:),e%k_s(i,:)
   enddo
   close(11)
 endif

 open(12,file='Q-boundary.dat')
 do i=1,e%n_s+1
    ! write(12,'(I10, f15.10)')i-1, e%boundary_s(i)
     write(12,'(f15.10, f15.3)')e%boundary_s(i), -100.0
     write(12,'(f15.10, f15.3)')e%boundary_s(i), 2000.0
     write(12,*)
 enddo
 close(12)
 
 open(13,file='q-list.dat')
 write(t,*) e%n_k
 write(13,*) adjustl(t) 
 do i=1, e%n_k
        write(13,'(3f15.10)') e%k_s(i,:)
 enddo
 close(13)
 
endsubroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  subroutine: read primary cell kpoints setting in input.uf  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_qpoints_setting_primary_cell(file, out)
character(len=*) :: file
integer, allocatable:: numbers(:)
real(kind=8), allocatable :: Kbegin(:,:), Kend(:,:)
integer :: N, k1, k2, status_file, status_read, status_k1, status_k2
integer :: i,  flag
character(len=readline_length) temp
type(array_car) :: out

open(1,file=trim(file),status='old',action='read',iostat=status_file)
flag=0
status_k1=0
status_k2=0
k1=0
k2=0
if(status_file==0) then
    readloop: do
        flag=flag+1
        read(1,'(a)',iostat=status_read) temp
        if(status_read /=0) exit
        temp=trim(lower_case(temp))
        i=index(temp,'#')
        if(i/=0)temp=temp(1:i-1)
        i=index(temp,'!')
        if(i/=0)temp=temp(1:i-1)
        temp=delete_space(temp)
        i=index(temp,'beginprimarycellqpoint')
        if(i/=0) then
            k1=flag
            status_k1=status_k1+1
        endif
        i=index(temp,'endprimarycellqpoint')
        if(i/=0) then
            k2=flag
            status_k2=status_k2+1
        endif
    end do readloop
else
    open(2,file='check',access='append')
    write(2,*)'err in opening file '//trim(file)
    close(2)
endif
k1=k1+1
k2=k2-1
N=(k2-k1+1)/3
close(1)
out%N1=N; out%N2=3


if(N.ge.1 .and. status_k1.eq.1 .and. status_k2.eq.1) then
    allocate(out%i_1d(out%N1),out%r_2d(out%N1*2,3))
    out%i_1d=0
    out%r_2d=0
    
    open(1,file=trim(file),status='old',action='read',iostat=status_file)
    if(status_file==0) then
        do i=1,k1-1
            read(1,*)
        enddo
        do i=1, out%N1
            !write(*,*) i
            read(1,*)out%i_1d(i)
            !write(*,*)out%i_1d(i)
            read(1,*)out%r_2d(2*i-1,:)
            read(1,*)out%r_2d(2*i  ,:)
            !write(*,*)out%r_2d(2*i-1,:)
            !write(*,*)out%r_2d(2*i  ,:)
        enddo
    else
        open(2,file='check',access='append')
        write(2,*)'err in opening file '//trim(file)
        close(2)
    endif
    close(1)
else
    open(2,file='check',access='append')
    write(2,*)'there is no (or more than 1 ) begin kpoint or end kpoint, or their positions are wrong'
    close(2)
endif

endsubroutine read_qpoints_setting_primary_cell


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  calculate total number of q-points 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function  q_points_number(file)
character(len=*) :: file
type(array_car) :: out
integer i

call read_qpoints_setting_primary_cell(file, out)

q_points_number=sum(out%i_1d)

endfunction







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  calculate total number of atoms in supercell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function atom_number_sc(file)
character(len=*) :: file
type(array_car) :: out

call read_super_cell_atom_positions(file,out)
atom_number_sc=out%N1

endfunction



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  subroutine: read primary cell in input.uf  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_primary_cell(file,out)
character(len=*) :: file
type(array_car) :: out
integer ::N, k1, k2, status_file, status_read, status_k1, status_k2
integer :: i,  flag
character(len=readline_length) temp

open(1,file=trim(file),status='old',action='read',iostat=status_file)
flag=0
status_k1=0
status_k2=0
k1=0
k2=0
if(status_file==0) then
    readloop: do
        flag=flag+1
        read(1,'(a)',iostat=status_read) temp
        if(status_read /=0) exit
        temp=trim(lower_case(temp))
        i=index(temp,'#')
        if(i/=0)temp=temp(1:i-1)
        i=index(temp,'!')
        if(i/=0)temp=temp(1:i-1)
        temp=delete_space(temp)
        i=index(temp,'beginprimarycellvectors')
        if(i/=0) then
            k1=flag
            status_k1=status_k1+1
        endif
        i=index(temp,'endprimarycellvectors')
        if(i/=0) then
            k2=flag
            status_k2=status_k2+1
        endif
    end do readloop
else
    open(2,file='check',access='append')
    write(2,*)'err in opening file '//trim(file)
    close(2)
endif
k1=k1+1
k2=k2-1
N=(k2-k1+1)
close(1)

out%N1=3; out%N2=3
allocate(out%r_2d(3,3))

if(N.eq.3 .and. status_k1.eq.1 .and. status_k2.eq.1) then
    open(1,file=trim(file),status='old',action='read',iostat=status_file)
    if(status_file==0) then
        do i=1,k1-1
            read(1,*)
        enddo
        do i=1, 3
            read(1,*)out%r_2d(i,:)
            !write(*,*)out%r_2d(i,:)
        enddo
    else
        open(2,file='check',access='append')
        write(2,*)'err in opening file '//trim(file)
        close(2)
    endif
    close(1)
else
    open(2,file='check',access='append')
    write(2,*)'there is no (or more than 1 ) begin primary cell vectors or end primary cell vectors, or their positions are wrong'
    close(2)
endif


endsubroutine read_primary_cell


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  subroutine: read super cell from input file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_super_cell(file,out)
character(len=*) :: file
type(array_car) :: out
integer ::N, k1, k2, status_file, status_read, status_k1, status_k2
integer :: i,  flag
character(len=readline_length) temp

open(1,file=trim(file),status='old',action='read',iostat=status_file)
flag=0
status_k1=0
status_k2=0
k1=0
k2=0
if(status_file==0) then
    readloop: do
        flag=flag+1
        read(1,'(a)',iostat=status_read) temp
        if(status_read /=0) exit
        temp=trim(lower_case(temp))
        i=index(temp,'#')
        if(i/=0)temp=temp(1:i-1)
        i=index(temp,'!')
        if(i/=0)temp=temp(1:i-1)
        temp=delete_space(temp)
        i=index(temp,'beginsupercellvectors')
        if(i/=0) then
            k1=flag
            status_k1=status_k1+1
        endif
        i=index(temp,'endsupercellvectors')
        if(i/=0) then
            k2=flag
            status_k2=status_k2+1
        endif
    end do readloop
else
    open(2,file='check',access='append')
    write(2,*)'err in opening file '//trim(file)
    close(2)
endif
k1=k1+1
k2=k2-1
N=(k2-k1+1)
close(1)

out%N1=3; out%N2=3
allocate(out%r_2d(3,3))
if(N.eq.3 .and. status_k1.eq.1 .and. status_k2.eq.1) then
    open(1,file=trim(file),status='old',action='read',iostat=status_file)
    if(status_file==0) then
        do i=1,k1-1
            read(1,*)
        enddo
        do i=1, 3
            read(1,*)out%r_2d(i,:)
            !write(*,*)out%r_2d(i,:)
        enddo
    else
        open(2,file='check',access='append')
        write(2,*)'err in opening file '//trim(file)
        close(2)
    endif
    close(1)
else
    open(2,file='check',access='append')
    write(2,*)'there is no (or more than 1 ) begin super cell vectors or end super cell vectors, or their positions are wrong'
    close(2)
endif


endsubroutine read_super_cell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  subroutine: read super cell from input file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_super_cell_atom_positions(file,out)
character(len=*) :: file
type(array_car) :: out
integer ::N, k1, k2, status_file, status_read, status_k1, status_k2
integer :: i,  flag
character(len=readline_length) temp

open(1,file=trim(file),status='old',action='read',iostat=status_file)
flag=0
status_k1=0
status_k2=0
k1=0
k2=0
if(status_file==0) then
    readloop: do
        flag=flag+1
        read(1,'(a)',iostat=status_read) temp
        if(status_read /=0) exit
        temp=trim(lower_case(temp))
        i=index(temp,'#')
        if(i/=0)temp=temp(1:i-1)
        i=index(temp,'!')
        if(i/=0)temp=temp(1:i-1)
        temp=delete_space(temp)
        i=index(temp,'beginsupercellatompositions')
        if(i/=0) then
            k1=flag
            status_k1=status_k1+1
        endif
        i=index(temp,'endsupercellatompositions')
        if(i/=0) then
            k2=flag
            status_k2=status_k2+1
        endif
    end do readloop
else
    open(2,file='check',access='append')
    write(2,*)'err in opening file '//trim(file)
    close(2)
endif
k1=k1+1
k2=k2-1
N=(k2-k1+1)
close(1)


out%N1=N; out%N2=3
allocate(out%r_2d(out%N1,out%N2))
if(status_k1.eq.1 .and. status_k2.eq.1) then
    open(1,file=trim(file),status='old',action='read',iostat=status_file)
    if(status_file==0) then
        do i=1,k1-1
            read(1,*)
        enddo
        do i=1, out%N1
            read(1,*)out%r_2d(i,:)
            !write(*,*)out%r_2d(i,:)
        enddo
    else
        open(2,file='check',access='append')
        write(2,*)'err in opening file '//trim(file)
        close(2)
    endif
    close(1)
else
    open(2,file='check',access='append')
    write(2,*)'there is no (or more than 1 ) begin super cell atomi positions &
               or end super cell atom positions, or their positions are wrong'
    close(2)
endif


endsubroutine read_super_cell_atom_positions


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  subroutine: read modes from file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_modes(e,file)
implicit none
character(len=*) :: file
character(len=500) :: temp
integer :: status_file, status_read, status_k1, status_k2
type(pebu) :: e
real(kind=8) :: t1, t2, t3, t4, t5, t6  
complex(kind=8) :: a(3)
integer :: i, j, k

e%n_k=q_points_number('input.dat')
e%n_atom_sc=atom_number_sc('input.dat')
e%n_b=e%n_atom_sc*3
allocate(e%e(e%n_k,e%n_b,2),e%c(e%n_k,e%n_b,e%n_atom_sc,3))
if(.not.allocated(e%k_s)) allocate(e%k_s(e%n_k,3))
if(.not.allocated(e%klength))allocate(e%klength(e%n_k))


open(1,file=trim(file),status='old',action='read',iostat=status_file)
! open(3,file='Q-points.dat')
! do i=1,e%n_k
!   read(3,*)e%klength(i)
! enddo
do i=1,e%n_k 
   !write(*,*) 'i=', i
   read(1,*);read(1,*)
   read(1,'(a)')temp
   read(temp(5:),*) e%k_s(i,:)
   !write(*,*)e%k_s(i,:); pause
   read(1,*)
   do j=1, e%n_b    
      !write(*,*)'k=',i, 'band=',j
      read(1,'(a)')temp
      read(temp(43:58),*)e%e(i,j,1)
      do k=1,e%n_atom_sc
          read(1,'(a)')temp
          read(temp(3:73),*)t1, t2, t3, t4, t5, t6
          e%c(i,j,k,1)=cmplx(t1,t2)
          e%c(i,j,k,2)=cmplx(t3,t4)
          e%c(i,j,k,3)=cmplx(t5,t6)
      enddo
   enddo
   read(1,*)
enddo

endsubroutine read_modes


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!    inverse a double real matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine inv_d(a,b)
real(kind=8), allocatable :: a(:,:), b(:,:)
integer :: i, j, n
   !-------------------------------for dgetrf
integer ::  info, lda
integer, allocatable :: ipiv2(:)
integer :: lwork
real(kind=8), allocatable :: work(:)
    !-------------------------------
n=size(a,1); lda=n
allocate(ipiv2(2*n))
if(.not.allocated(b)) allocate(b(n,n))
b=a
call dgetrf(n,n,b,lda,ipiv2,info)
if(info.ne.0) write(*,*)'Err: dgetrf'

lwork=2*n
if(.not.allocated(work)) allocate(work(lwork))
call dgetri(n, b, lda, ipiv2, work, lwork, info)
if(info.ne.0) write(*,*)'Err: dgetri'

end subroutine
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!    Computes the LU factorization of a general m-by-n complex matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine lu_d(a,b)
real(kind=8), allocatable :: a(:,:), b(:,:)
integer :: i, j, n, m
    !-------------------------------for dgetrf
integer ::  info, lda
integer, allocatable :: ipiv2(:)
    !-------------------------------
n=size(a,1); m=size(a,2)
lda=n
allocate(ipiv2(m+n))
if(.not.allocated(b)) allocate(b(n,m))
b=a
call dgetrf(n,m,b,lda,ipiv2,info)

end subroutine lu_d





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  subroutine: Unfolding phonon dispersion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Unfold_Phonon_planewave(e)
implicit none
type(pebu) :: e
type(array_car) cell1_read, cell2_read, r_atoms
integer :: i, j, k, L, ii, jj, kk, iii, jjj, kkk
real(kind=8) :: ikx, iky, ikz, wtclean
integer :: i_k, i_b, i_a
complex(kind=8) ::  t1_c, t_c
real(kind=8) :: t(3,3), q(3), t_sc(3,3), t_pc(3,3),tt_pc(3,3), t_group(3),t_ipc_unit_isc(3,3)
real(kind=8),allocatable :: q_group(:,:), t_new(:,:), t1_new(:,:)
integer :: n1, n2, n3, n_qgroup, n_flag, max_nx=2, max_ny=2, max_nz=2
real(kind=8):: t_r, t1_r
complex(kind=8), allocatable :: phase(:), phase_temp(:)
integer :: ipiv(3)
! read atom position
if(.not.allocated(e%r_atom_sc)) allocate(e%r_atom_sc(e%n_atom_sc,3),e%r_direct_atom_sc(e%n_atom_sc,3))
allocate(phase(e%n_atom_sc),phase_temp(e%n_atom_sc))

     max_nx=read_integer('input.dat','max_qx', 2)
     max_ny=read_integer('input.dat','max_qy', 2)
     max_nz=read_integer('input.dat','max_qz', 2)

     ! read the primary cell lattice parameters
     call read_primary_cell('input.dat',cell1_read)
     e%p_lat(1,:)=cell1_read%r_2d(1,:)
     e%p_lat(2,:)=cell1_read%r_2d(2,:)
     e%p_lat(3,:)=cell1_read%r_2d(3,:)
     ! read the super cell lattice parameters
     call read_super_cell('input.dat',cell2_read)
     e%s_lat(1,:)=cell2_read%r_2d(1,:)
     e%s_lat(2,:)=cell2_read%r_2d(2,:)
     e%s_lat(3,:)=cell2_read%r_2d(3,:)
    
     allocate(t_new(3,3))
     t_new=dble(e%s_lat)     !  prepare inverse matrix
     !call getrf(t,ipiv)  ! prepare inverse matrix
     call lu_d(t_new,t1_new)
     t_r=t1_new(1,1)*t1_new(2,2)*t1_new(3,3) ! det(e%s_lat)
     !call getri(t,ipiv)  ! inverse matrix
     call inv_d(t_new,t_new)
     t_sc=t_new  ! inverse super cell
     t_new=dble(e%p_lat)     !  prepare inverse matrix
     !call getrf(t,ipiv)  ! prepare inverse matrix
     call lu_d(t_new,t1_new) 
     t_r=abs(t_r/(t1_new(1,1)*t1_new(2,2)*t1_new(3,3))) ! det(e%s_lat)/det(e%p_lat)
     !call getri(t,ipiv)  ! inverse matrix
     call inv_d(t_new,t_new)
     t_pc=t_new  ! inverse primary cell
     n_qgroup=anint(t_r) ! number of q points in the primary cell Brillouin zone
   !  write(*,*)n_qgroup
   !  pause
     allocate(q_group(n_qgroup,3)) 
     q_group=0d0

     t_ipc_unit_isc=matmul(e%s_lat,t_pc) ! inverse primary cell in the unit of inverse super cell 

   ! write(*,*)t_ipc_unit_isc(1,:)
   ! write(*,*)t_ipc_unit_isc(2,:)
   ! write(*,*)t_ipc_unit_isc(3,:)
   ! pause
    
     call read_super_cell_atom_positions('input.dat',r_atoms)
     do i=1,e%n_atom_sc
          e%r_atom_sc(i,:)=r_atoms%r_2d(i,:)
          ! write(*,*)i,e%r_atom_sc(i,:)
     enddo
     e%r_direct_atom_sc= matmul(e%r_atom_sc, t_sc)  ! atom positions in a direct form

     
     t=matmul(e%p_lat,t_sc) ! inverse supercell lattice in the unit of inverse primary cell lattice
     t=transpose(t)
    ! write(*,*)'check'
    ! write(*,*)t(1,:)
    ! write(*,*)t(2,:)
    ! write(*,*)t(3,:)
    ! pause
     
     n1=3d0/dsqrt(dot_product(t(1,:),t(1,:)))
     n2=3d0/dsqrt(dot_product(t(2,:),t(2,:)))
     n3=3d0/dsqrt(dot_product(t(3,:),t(3,:)))
         ! write(*,*) n1,n2,n3
         ! pause

     ! find the group of q points that in the first Brillouin zone of primary cell 
     ! The total q point number is n_qgroup, the first q point is (0 0 0).
     n_flag=1
     do i=-n1,n1;do j=-n2,n2;do k=-n3,n3
        t_group=i*t(1,:)+j*t(2,:)+k*t(3,:)
        if(t_group(1).ge.-0.499 .and. t_group(1).lt.0.501 .and. t_group(2).ge.-0.499 .and. t_group(2).lt.0.501 &
                     .and. t_group(3).ge.-0.499 .and. t_group(3).lt.0.501) then
                if(abs(t_group(1)).gt.1e-4 .or. abs(t_group(2)).gt.1e-4 .or. abs(t_group(3)).gt.1e-4) then
              n_flag=n_flag+1
              ! write(*,*) t_group
              q_group(n_flag,:)=t_group
              ! write(*,*)q_group(n_flag,:)
           endif
        endif
     enddo;enddo;enddo
    ! write(*,*)q_group(1,:)
    ! write(*,*)q_group(2,:)
    ! pause
    

     q_group=matmul(matmul(q_group,transpose(t_pc)),transpose(e%s_lat)) ! in the unit of inverse supercell lattice vectors 

        ! q_group=matmul(q_group,transpose(t_pc))
        ! q_group=matmul(q_group,transpose(e%s_lat))
        ! write(*,*)q_group(1,:)
    ! write(*,*)q_group(2,:)
    ! pause

t=transpose(t_sc)
tt_pc=transpose(t_pc)
do i_k=1, e%n_k
! do i_k=1, 1
   write(*,*) 'q=',i_k
!   write(*,*)'q=', e%k_s(i_k,:)
! pause
   q(1)=e%k_s(i_k,1)*t(1,1)+e%k_s(i_k,2)*t(1,2)+e%k_s(i_k,3)*t(1,3)
   q(2)=e%k_s(i_k,1)*t(2,1)+e%k_s(i_k,2)*t(2,2)+e%k_s(i_k,3)*t(2,3)
   q(3)=e%k_s(i_k,1)*t(3,1)+e%k_s(i_k,2)*t(3,2)+e%k_s(i_k,3)*t(3,3)
 ! write(*,*)'q=', q
  ! pause
   ! prepare phase
   phase=0d0
   do i_a=1,e%n_atom_sc
      t1_r=q(1)*e%r_atom_sc(i_a,1)+q(2)*e%r_atom_sc(i_a,2)+q(3)*e%r_atom_sc(i_a,3)
      phase(i_a)=exp(-2*pi*I_imag*t1_r)
      !write(*,*)'atom=',i_a,  phase(i_a)  
      !write(*,'("atom=",I3, 3f8.4, "  (",2f8.4,")")') i_a, e%r_atom_sc(i_a,:),phase(i_a)
   enddo
   ! write(*,*) phase
   ! pause

   do i_b=1,e%n_b
     do i=1, e%n_atom_sc    
       ! write(*,'(6f10.5)')e%c(i_k,i_b,i,1)*phase(i), e%c(i_k,i_b,i,2)*phase(i), e%c(i_k,i_b,i,3)*phase(i)
        e%c(i_k,i_b,i,1)=e%c(i_k,i_b,i,1)*phase(i)
        e%c(i_k,i_b,i,2)=e%c(i_k,i_b,i,2)*phase(i)
        e%c(i_k,i_b,i,3)=e%c(i_k,i_b,i,3)*phase(i)
     enddo 
     ! pause
     t_r=0d0;  t1_r=0d0
     !do i=-max_nx, max_nx; do j=-max_ny, max_ny; do k=-max_nz, max_nz
      do i=-1, 1; do j=-1, 1; do k=0, 0
         
         ! write(*,*) t_pc(1,:)
         ! write(*,*) t_pc(2,:)
         ! write(*,*) t_pc(3,:)
         ! pause
            ikx=i*tt_pc(1,1)+j*tt_pc(2,1)+k*tt_pc(3,1)
        iky=i*tt_pc(1,2)+j*tt_pc(2,2)+k*tt_pc(3,2)
        ikz=i*tt_pc(1,3)+j*tt_pc(2,3)+k*tt_pc(3,3)
        
              !  write(*,*)
        do ii=1,e%n_atom_sc
                     phase_temp(ii)=exp(-I_imag*2*pi*(ikx*e%r_atom_sc(ii,1)+iky*e%r_atom_sc(ii,2)+ikz*e%r_atom_sc(ii,3)))
                    ! write(*,*)phase_temp(ii)
                    ! pause
                  !  write(*,*)e%r_atom_sc(ii,:)
        enddo
           ! write(*,*)
           ! write(*,*) i,j,k
           ! write(*,*)ikx, iky, ikz
           ! write(*,*) phase_temp(:)
           ! pause

       do jj=1,3
          t_c=0d0
          do ii=1,e%n_atom_sc

             t_c=t_c+e%c(i_k,i_b,ii,jj)*phase_temp(ii)
             !write(*,'(4f8.4)')e%c(i_k,i_b,ii,jj), phase_temp(ii)
          enddo
        !  write(*,*)'adds to ', abs(dble(t_c))
        !  pause
         
         t_r=t_r+abs(dble(t_c))
        enddo
        
    do kk=2,n_qgroup
         ! write(*,*)'q_group',kk
                  ikx=i*tt_pc(1,1)+j*tt_pc(2,1)+k*tt_pc(3,1)+q_group(kk,1)*t(1,1)+q_group(kk,2)*t(2,1)+q_group(kk,3)*t(3,1)
          iky=i*tt_pc(1,2)+j*tt_pc(2,2)+k*tt_pc(3,2)+q_group(kk,1)*t(1,2)+q_group(kk,2)*t(2,2)+q_group(kk,3)*t(3,2)
          ikz=i*tt_pc(1,3)+j*tt_pc(2,3)+k*tt_pc(3,3)+q_group(kk,1)*t(1,3)+q_group(kk,2)*t(2,3)+q_group(kk,3)*t(3,3)
          do ii=1,e%n_atom_sc
                       phase_temp(ii)=exp(-I_imag*2*pi*(ikx*e%r_atom_sc(ii,1)+iky*e%r_atom_sc(ii,2)+ikz*e%r_atom_sc(ii,3)))
          enddo

          do jj=1,3
          t1_c=0d0
          do ii=1,e%n_atom_sc
             t1_c=t1_c+e%c(i_k,i_b,ii,jj)*phase_temp(ii)
            ! write(*,'(4f8.4)')e%c(i_k,i_b,ii,jj), phase_temp(ii)
          enddo
         ! write(*,*)'adds to ', abs(dble(t1_c))
          !pause
          t1_r=t1_r+abs(dble(t1_c))
          enddo
        enddo
        !pause
     enddo;enddo;enddo

       ! write(*,'(2I5,f20.5)')i_k, i_b, t_r/(t_r+t1_r)
    ! pause

     e%e(i_k,i_b,2)=t_r/(t_r+t1_r)
   enddo
 
 enddo

wtclean=read_real('input.dat','wtclean',0.0d0)

open(2,file='unfold.dat')
do i_b=1,e%n_b; do i_k=1,e%n_k
   if (e%e(i_k,i_b,2) .gt. wtclean) then
      write(2,'(3f18.10)') e%klength(i_k),e%e(i_k,i_b,1),e%e(i_k,i_b,2)
   endif
enddo
   write(2,*)
enddo
close(2)

endsubroutine Unfold_Phonon_planewave



end module
