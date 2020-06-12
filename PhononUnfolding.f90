!  PhononUnfolding.f90 
!
!  FUNCTIONS:
!   main program for unfold phonon dispersions.
!
!****************************************************************************

program PhononUnfolding
 use PhononUnfoldingModule

 implicit none
 type(pebu) :: eu
 character(len=20) :: calculation
 character(len=20) :: phonon_code
 integer :: i
 integer :: temp

 type(list):: out


 write(*,*) 'Hello, this is PhononUnfolding Code.'
 write(*,*) 'Version: 1.0.20200601'

 calculation=trim(read_characters('input.dat','calculation','uf'))
 phonon_code=trim(read_characters('input.dat','phonon_code','qe'))

 if(trim(calculation).eq.'uf') then
     call construct_qlist(eu)
     if(trim(phonon_code).eq.'qe') then
        call read_modes(eu,'matdyn.modes')
     elseif(trim(phonon_code).eq.'phonopy') then
        call read_modes_phonopy(eu,'qpoints.yaml')
     endif
     call Unfold_Phonon_planewave(eu)   
 elseif(trim(calculation).eq.'qp')then
     write(*,*) 'q-point list for Phonon'
     call construct_qlist(eu)
     call write_qlist_to_file(eu)
 elseif(trim(calculation).eq.'abinitgamma')then
     call read_modes_abinitgamma(eu,'band_Gamma.yaml')
     call Unfold_Phonon_planewave_abinitgamma(eu) 
 elseif(trim(calculation).eq.'abinitgamma_direct')then
     call read_modes_abinitgamma(eu,'band_Gamma.yaml')
     call Unfold_Phonon_gamma(eu)
 elseif(trim(calculation).eq.'test') then
     call construct_qlist(eu)
     call read_modes_abinit(eu,'qpoints.yaml')
     call Unfold_Phonon_abinit(eu)
 endif


write(*,*)'calculation '//trim(calculation)//' is done.'


end program 
