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
 integer :: i
 integer :: temp

 write(*,*) 'Hello, this is PhononUnfolding Code.'
 write(*,*) 'Version: 1.0.20160229'

 calculation=trim(read_characters('input.dat','calculation','uf'))


 if(trim(calculation).eq.'uf') then
     call construct_qlist(eu)
     call read_modes(eu,'matdyn.modes')
     call Unfold_Phonon_planewave(eu)   
 elseif(trim(calculation).eq.'qp')then
     write(*,*) 'q-point list for Phonon'
     call construct_qlist(eu)
     call write_qlist_to_file(eu)
endif


write(*,*)'calculation '//trim(calculation)//' is done.'


end program 
