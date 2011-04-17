!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!=============================================================================== 
! MODULE od_dos  OptaDOS - Density of States
! This is the module that contains all if the DOS routines. It is used through 
! the global calculate_dos subroutine
!-------------------------------------------------------------------------------
! Three other global routines are available, dos_merge, doslin and 
! doslin_sub_cell_corners, these are currently used by the jdos module and 
! routines.
!-------------------------------------------------------------------------------
! Written by: A J Morris Nov - Dec 2010 Modified from LinDOS (CJP+AJM)
!=============================================================================== 
module od_dos_utils
  use od_constants, only : dp
  use od_electronic, only: matrix_weights_array_boundaries
  implicit none

  !-------------------------------------------------------------------------------
  ! D E R I V E D   P R O T O T Y P E S 

  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  ! P U B L I C   V A R I A B L E S
  real(kind=dp), allocatable, public, save :: dos_adaptive(:,:)
  real(kind=dp), allocatable, public, save :: dos_fixed(:,:)
  real(kind=dp), allocatable, public, save :: dos_linear(:,:)

  real(kind=dp), allocatable, public, save :: intdos_adaptive(:,:)
  real(kind=dp), allocatable, public, save :: intdos_fixed(:,:) 
  real(kind=dp), allocatable, public, save :: intdos_linear(:,:)

  real(kind=dp), allocatable, public, save :: E(:)

  real(kind=dp), public, save :: efermi_fixed
  real(kind=dp), public, save :: efermi_adaptive
  real(kind=dp), public, save :: efermi_linear
  !real(kind=dp), public, save :: efermi_quad
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  ! P U B L I C   F U N C T I O N S 
  public :: dos_utils_calculate
  public :: dos_utils_deallocate
  public :: dos_utils_calculate_at_e
  public :: dos_utils_merge          ! Used by od_jdos
  public :: doslin_sub_cell_corners  ! Used by od_jdos
  public :: doslin                   ! Used by od_jdos
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  ! P R I V A T E   V A R I A B L E S 
  private ! unless otherwise indicated

  type(matrix_weights_array_boundaries) :: mw

  real(kind=dp), save                   :: delta_bins ! Width of bins
  logical :: calc_weighted_dos
  !-------------------------------------------------------------------------------

contains

  !=============================================================================== 
  subroutine dos_utils_calculate(matrix_weights, weighted_dos)
    !===============================================================================  
    ! Main routine in dos module, drives the calculation of density of states for
    ! both task : dos and also if it is required elsewhere.
    !------------------------------------------------------------------------------- 
    ! Arguments: matrix_weigths (in) (opt) : LCAO or other weightings for DOS
    !            weighted_dos   (out)(opt) : Output DOS weigthed by matrix_weights
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: mw, E, dos_adaptive, dos_fixed, dos_linear
    ! intdos_adaptive, intdos_fixed, intdos_linear, efermi_fixed, efermi_adaptive
    ! efermi_linear, delta_bins, calc_weighted_dos
    !-------------------------------------------------------------------------------
    ! Modules Used: see below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: One of linear, adaptive or fixed must be .true.
    !-------------------------------------------------------------------------------
    ! Known Worries: (1) If more than one of linear, adaptive or fixed are set it 
    ! uses the most complicated method.
    ! (2) It should be possible to pass optioinal arguments to sub programs as
    ! optional argumnets without checking whether they are there or not. g95 will 
    ! allow this behaviour. gfotran will not.
    !-------------------------------------------------------------------------------
    ! Written by : A J Morris December 2010
    !=============================================================================== 
    use od_io,        only : stdout,io_time,io_error
    use od_comms,     only : on_root,my_node_id
    use od_electronic,only : band_gradient,band_energy, efermi, efermi_castep,nspins, &
         & elec_read_band_gradient,elec_read_band_energy,elec_report_parameters
    use od_parameters,only : linear, adaptive, fixed, quad, compute_band_energy, &
         & compute_efermi,dos,dos_per_volume,fermi_energy,iprint
    use od_cell,         only : cell_volume, nkpoints, cell_calc_lattice, &
         & cell_report_parameters,cell_dist,kpoint_grid_dim
    implicit none

    !-------------------------------------------------------------------------------
    ! I N T E R N A L   V A R I A B L E S
    integer :: ierr, idos, i, ik, is, ib
    real(kind=dp) :: time0, time1
    real(kind=dp),intent(in), allocatable, optional  :: matrix_weights(:,:,:,:)
    real(kind=dp),intent(out),allocatable, optional  :: weighted_dos(:,:,:) ! bins.spins, orbitals
    !-------------------------------------------------------------------------------

    if(.not.(linear.or.adaptive.or.fixed.or.quad)) call io_error (" DOS: No Broadening Set")
    

    calc_weighted_dos=.false.
    if(present(matrix_weights)) calc_weighted_dos=.true.
    
    
    if(calc_weighted_dos.eqv..false.) then ! We are called just to provide dos.
      if(allocated(E)) then
         if(on_root) write(stdout,*) " Already calculated dos, so returning..."
         return  ! The dos has already been calculated previously so just return.       
      endif
    endif

    if(calc_weighted_dos)then
       mw%norbitals=size(matrix_weights,1)
       mw%nbands   =size(matrix_weights,2)
       mw%nkpoints =size(matrix_weights,3)
       mw%nspins   =size(matrix_weights,4)
    end if

    if(on_root) then
       write(stdout,*)
       write(stdout,'(1x,a78)') '+============================================================================+'
       if(calc_weighted_dos) then
          write(stdout,'(1x,a78)') '+============== Weighted Density Of States Calculation ======================+'
       else
          write(stdout,'(1x,a78)')'+=================== Density Of States Calculation ==========================+'
       endif
       write(stdout,'(1x,a78)') '+============================================================================+'
       write(stdout,*)
    endif

    if(calc_weighted_dos) then 
       if(mw%nspins.ne.nspins)     call io_error ("ERROR : DOS :  mw%nspins not equal to nspins.")
       if(mw%nkpoints.ne.nkpoints) call io_error ("ERROR : DOS : mw%nkpoints not equal to nkpoints.")
    endif

    !-------------------------------------------------------------------------------
    ! R E A D   B A N D   G R A D I E N T S 
    ! If we're using one of the more accurate roadening schemes we also need to read in the 
    ! band gradients too
    if(quad.or.linear.or.adaptive) then
      if(.not.allocated(band_gradient)) call elec_read_band_gradient
    endif
    !-------------------------------------------------------------------------------
    ! C A L C U L A T E   D O S 
    ! Now everything is set up, we can perform the dos accumulation in parellel
    time0=io_time()

    call setup_energy_scale

    if(on_root.and.(iprint>1)) write(stdout,*)

    if(fixed)then
       if(calc_weighted_dos)then 
          call calculate_dos("f",dos_fixed,intdos_fixed, matrix_weights=matrix_weights, weighted_dos=weighted_dos) 
          call dos_utils_merge(dos_fixed,weighted_dos=weighted_dos)    
       else
          call calculate_dos("f",dos_fixed,intdos_fixed)
          call dos_utils_merge(dos_fixed) 
       endif
       call dos_utils_merge(intdos_fixed)
    endif
    if(adaptive)then
       if(calc_weighted_dos)then 
          call calculate_dos("a",dos_adaptive,intdos_adaptive, matrix_weights=matrix_weights, weighted_dos=weighted_dos)
          call dos_utils_merge(dos_adaptive,weighted_dos=weighted_dos) 
       else
          call calculate_dos("a",dos_adaptive,intdos_adaptive)
          call dos_utils_merge(dos_adaptive)
       endif
       call dos_utils_merge(intdos_adaptive)
    endif
    if(linear)then
       if(calc_weighted_dos)then 
          call calculate_dos("l",dos_linear, intdos_linear, matrix_weights=matrix_weights, weighted_dos=weighted_dos)
          call dos_utils_merge(dos_linear,weighted_dos=weighted_dos)
       else      
          call calculate_dos("l",dos_linear,intdos_linear)
          call dos_utils_merge(dos_linear)
       endif
       call dos_utils_merge(intdos_linear)
    endif

    if(quad) then
       call io_error("quadratic broadening not implemented") 
       !if(quad)    call merge_dos(dos_quad)
       !if(quad)    call merge_dos(intdos_quad)
    endif

    if(.not.on_root) then
       if(allocated(E)) deallocate(E, stat=ierr)
       if (ierr/=0) call io_error ("cannot deallocate  E")
    endif

    time1=io_time()
    if(on_root)  write(stdout,'(1x,a40,f11.3,a)') 'Time to calculate dos  ',time1-time0,' (sec)'
    !-------------------------------------------------------------------------------


    !-------------------------------------------------------------------------------
    ! F E R M I   E N E R G Y   A N A L Y S I S
    time0=io_time()
    if(on_root)       write(stdout,*)
    if(on_root)   write(stdout,'(1x,a78)')  '+------------------------ Fermi Energy Analysis -----------------------------+'
    if(on_root)   write(stdout,'(1x,a1,a45,f8.4,a3,20x,a1)') "|","Fermi energy from CASTEP :",efermi_castep," eV","|"
    if(on_root)   write(stdout,'(1x,a78)')    '+----------------------------------------------------------------------------+'
       
       if(compute_efermi) then
          if(fixed) then 
             if(on_root)write(stdout,'(1x,a78)') "| From Fixed broadening                                                      | "
             efermi_fixed= calc_efermi_from_intdos(intdos_fixed)
             if(on_root)write(stdout,'(1x,a1,a46,f8.4,a3,19x,a1)')"|", " Fermi energy (Fixed braodening) : ", efermi_fixed,"eV","|"
             if(on_root)write(stdout,'(1x,a78)')    '+----------------------------------------------------------------------------+'
             
          endif
          if(adaptive) then
             if(on_root)write(stdout,'(1x,a78)') "| From Adaptive broadening                                                   | "  
             efermi_adaptive=calc_efermi_from_intdos(intdos_adaptive)
             if(on_root)write(stdout,'(1x,a1,a46,f8.4,a3,19x,a1)')"|", " Fermi energy (Adaptive braodening) : " &
                  , efermi_adaptive,"eV","|"
             if(on_root)write(stdout,'(1x,a78)') &
                  '+----------------------------------------------------------------------------+'
             
          endif
          if(linear) then
             if(on_root)write(stdout,'(1x,a78)') &
                  "| From Linear broadening                                                     | " 
             efermi_linear=calc_efermi_from_intdos(intdos_linear) 
             if(on_root)write(stdout,'(1x,a1,a46,f8.4,a3,19x,a1)')"|", " Fermi energy (Linear braodening) : ",&
                  efermi_linear," eV","|"
             if(on_root)write(stdout,'(1x,a78)')  &
                  '+----------------------------------------------------------------------------+'
             
          endif
       else
          if(on_root) write(stdout,'(1x,a78)')   "| No Fermi energies calculated                                               | "
          ! Use the derived Fermi energy shifts if calculated. It not use the fermi_energy supplied
          ! if not, use the CASTEP Fermi energy
          if(fermi_energy.ne.-990.0_dp)then
             efermi_fixed=fermi_energy
             efermi_adaptive=fermi_energy
             efermi_linear=fermi_energy
          else
             efermi_fixed=efermi_castep
             efermi_adaptive=efermi_castep
             efermi_linear=efermi_castep
          endif
       endif
       
       ! NB If you have asked for more than one type of broadening
       ! If one of your options is linear then all will have the linear efermi
       ! If you have asked for adaptive, but nor linear, you will have the adaptive efermi
       if(fixed) then
          efermi=efermi_fixed  
       endif
       
       if(adaptive) then
          efermi=efermi_adaptive
       endif
       
       if(linear)then
          efermi=efermi_linear
       endif
       
       if(on_root) then
          write(stdout,'(1x,a1,a46,f8.4,a3,19x,a1)')"|", " Fermi energy used : ", efermi,"eV","|"
          E(:)=E(:)-efermi
          band_energy(:,:,:) = band_energy(:,:,:) - efermi
          
          write(stdout,'(1x,a78)')    '+----------------------------------------------------------------------------+'
          time1=io_time()
          write(stdout,'(1x,a40,f11.3,a)') 'Time to calculate Fermi energies ',time1-time0,' (sec)'
          !-------------------------------------------------------------------------------
       end if
       
       !-------------------------------------------------------------------------------
       ! B A N D   E N E R G Y   A N A L Y S I S
       ! Now for a bit of crosschecking  band energies
       ! These should all converge to the same number as the number of bins is increased
       if(compute_band_energy) call compute_band_energies
       !-------------------------------------------------------------------------------
       
       if(on_root) then
          if(dos_per_volume) then
             if(fixed) then
                dos_fixed=dos_fixed/cell_volume   
                intdos_fixed=intdos_fixed/cell_volume
             endif
             if(adaptive) then
                dos_adaptive=dos_adaptive/cell_volume 
                intdos_adaptive=intdos_adaptive/cell_volume
             endif
             if(linear) then
                dos_linear=dos_linear/cell_volume     
                intdos_linear=intdos_linear/cell_volume
             endif
             if(calc_weighted_dos) then
                weighted_dos=weighted_dos/cell_volume
             endif
             ! if(quad) then
             !    dos_quad=dos_quad/cell_volume         
             !    intdos_quad=intdos_quad/cell_volume
             ! endif   
          endif
       endif




    if(on_root) then
       write(stdout,'(1x,a78)')    '+============================================================================+'
       if(calc_weighted_dos) then
          write(stdout,'(1x,a78)')    '+=========== Weighted Density Of States Calculation End =====================+'
       else
          write(stdout,'(1x,a78)')    '+================  Density Of States Calculation End ========================+'
       endif
       write(stdout,'(1x,a78)')    '+============================================================================+'
       write(stdout,*)
    end if


  end subroutine dos_utils_calculate


  !===============================================================================
  function calc_band_energies(dos)
    !===============================================================================  
    ! Function to return the band energy of a DOS by summing DOS*Energy up to 
    ! Fermi energy
    !------------------------------------------------------------------------------- 
    ! Arguments: dos (in) : Density of States array
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: E
    !-------------------------------------------------------------------------------
    ! Modules Used: See below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: none beyond the parent module variables
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: efermi must be set in od_electonic
    !                       E must be set
    !-------------------------------------------------------------------------------
    ! Known Worries: 
    !-------------------------------------------------------------------------------
    ! Written by : A J Morris December 2010
    !=============================================================================== 
    use od_electronic, only : nspins
    use od_parameters, only : dos_nbins
    use od_electronic, only : efermi
    implicit none

    real(kind=dp),intent(in) :: dos(1:dos_nbins,1:nspins)
    real(kind=dp) :: gband, calc_band_energies
    integer :: is, idos

    gband=0.0_dp
    do is=1,nspins
       do idos=1,dos_nbins
          if(E(idos).le.efermi) then
             gband = gband + dos(idos,is)*E(idos)*delta_bins
          else
             exit   
          endif
       end do
    enddo

    calc_band_energies=gband
  endfunction calc_band_energies


  !=============================================================================== 
  function calc_efermi_from_intdos(INTDOS)
    !===============================================================================  
    ! Function to calculate the Fermi energy from an intgrated DOS by serching for 
    ! the bin which contains the correct number of electrons
    !------------------------------------------------------------------------------- 
    ! Arguments: intdos (in) : Density of States array
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: E
    !-------------------------------------------------------------------------------
    ! Modules Used: See below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: E must be set
    !-------------------------------------------------------------------------------
    ! Known Worries: 
    !-------------------------------------------------------------------------------
    ! Written by : A J Morris December 2010 Modified from LinDOS
    !=============================================================================== 
    use od_electronic, only : num_electrons, nspins
    use od_parameters, only : dos_nbins
    use od_io,         only : stdout

    implicit none
    real(kind=dp), intent(in) :: INTDOS(1:dos_nbins,1:nspins)
    real(kind=dp) :: calc_efermi_from_intdos
    real(kind=dp) :: efermi
    real(kind=dp) :: tolerance=0.0001_dp ! Has the effect of forcing the efermi to 
    ! be in the middle of the band gap, and not
    ! the leading edge.
    integer  :: idos, i,j

    do idos=1,dos_nbins
       if(sum(INTDOS(idos,:)).gt.(sum(num_electrons(:))+(tolerance/2.0_dp))) exit
    end do

    efermi = E(idos)

    do i=idos,1,-1
       if(sum(INTDOS(i,:)).lt.(sum(INTDOS(idos-1,:))-(tolerance/2.0_dp))) exit
    end do

    calc_efermi_from_intdos = (efermi + E(i))/2.0_dp

    do j=1,nspins
       write(stdout,'(1x,a1,a20,i1,a20,f12.5,a6,f12.5,5x,a1)') "|","Spin Component : ",j,&
            &" occupation between ", INTDOS(i,j), "and", INTDOS(idos,j),"|"
    enddo

  end function calc_efermi_from_intdos



  !=============================================================================== 
  subroutine compute_band_energies
    !=============================================================================== 
    ! High-level subroutine to compute band energies of the DOS calculated.
    ! Calculates using the band_energies directly and compares with the
    ! function calc_band_energies which does the low level computation on the DOS.
    !------------------------------------------------------------------------------- 
    ! Arguments: None
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: E,dos_fixed,dos_adaptive,dos_linear
    !-------------------------------------------------------------------------------
    ! Modules Used: See below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: E must be set
    !-------------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------------
    ! Written by : A J Morris December 2010 
    !=============================================================================== 
    use od_parameters, only : adaptive, linear, fixed
    use od_electronic, only : electrons_per_state,efermi,nbands,nspins,band_energy
    use od_cell,       only : nkpoints, kpoint_weight,num_kpoints_on_node
    use od_comms,      only : comms_reduce,my_node_id,on_root
    use od_io,         only : stdout,io_time

    implicit none
    real(kind=dp) :: eband
    real(kind=dp) :: time0, time1
    integer :: ik,is,ib

    time0=io_time()
    if(on_root) then
       write(stdout,*)
       write(stdout,'(1x,a78)')    '+--------------------------- Band Energy Analysis ---------------------------+'
       
       if(fixed)then
          write(stdout,'(1x,a1,a46,f12.4,a3,15x,a1)')"|",&
               " Band energy (Fixed broadening)  : ",calc_band_energies(dos_fixed),"eV","|"
       endif
       if(adaptive)then
          write(stdout,'(1x,a1,a46,f12.4,a3,15x,a1)')"|",&
               " Band energy (Adaptive broadening) : ",calc_band_energies(dos_adaptive),"eV","|"
       endif
       if(linear)then 
          write(stdout,'(1x,a1,a46,f12.4,a3,15x,a1)')"|", &
               " Band energy (Linear broadening) : ",calc_band_energies(dos_linear),"eV","|"
       endif
    endif
    
    eband = 0.0_dp
    do ik=1,num_kpoints_on_node(my_node_id)
       do is=1,nspins
          do ib=1,nbands
             if(band_energy(ib,is,ik).le.efermi) eband = eband + band_energy(ib,is,ik)*electrons_per_state*kpoint_weight(ik)
          end do
       end do
    end do
    call comms_reduce(eband,1,'SUM')
    if(on_root) then
       write(stdout,'(1x,a1,a46,f12.4,a3,15x,a1)')"|", " Band energy (From CASTEP) : ", eband,"eV","|"
       write(stdout,'(1x,a78)')    '+----------------------------------------------------------------------------+'
       time1=io_time()
       write(stdout,'(1x,a40,f11.3,a)') 'Time to calculate Band energies ',time1-time0,' (sec)'
    end if

  end subroutine compute_band_energies


 

  !=============================================================================== 
  subroutine setup_energy_scale
    !=============================================================================== 
    ! Sets up all broadening independent DOS concerns. That is basically the energy
    ! scale, E, and the width of its bins, delta_bins
    !------------------------------------------------------------------------------- 
    ! Arguments: None
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: delta_bins
    !-------------------------------------------------------------------------------
    ! Modules Used: See below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: band_energy in od_electronic must be set
    !-------------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------------
    ! Written by : A J Morris December 2010 
    !=============================================================================== 
    use od_parameters, only : dos_nbins,dos_min_energy,dos_max_energy,dos_spacing,iprint
    use od_electronic, only : band_energy
    use od_io,         only : io_error,stdout
    use od_comms!,      only : comms_reduce

    implicit none

    real(kind=dp) :: min_band_energy, max_band_energy 
    integer       :: idos,i,ierr


    ! If we do have dos_min_energy and dos_max_energy set, then we'd better
    ! use them. If not, let's set some sensible values.
    if(dos_min_energy==-huge(dos_min_energy)) then !Do it automatically
       min_band_energy  = minval(band_energy)-5.0_dp
    else
       min_band_energy=dos_min_energy
    endif
    call comms_reduce(min_band_energy,1,'MIN')
    call comms_bcast(min_band_energy,1)
    
    if(dos_max_energy==huge(dos_max_energy)) then !Do it automatically
       max_band_energy  = maxval(band_energy)+5.0_dp
    else
       max_band_energy=dos_max_energy
    endif
    call comms_reduce(max_band_energy,1,'MAX')
    call comms_bcast(max_band_energy,1)

    ! If dos_nbins is set, then we'd better use that
    if(dos_nbins<0) then ! we'll have to work it out
       dos_nbins=abs(ceiling((max_band_energy-min_band_energy)/dos_spacing))
       ! Now modify the max_band_energy
       if(on_root.and.(iprint>2)) then
          write(stdout,*)
          write(stdout,'(1x,a40,f11.3,a13)') 'max_band_energy (before correction) : ', max_band_energy, " <-- DOS Grid"
       endif
       max_band_energy=min_band_energy+dos_nbins*dos_spacing
    endif

    allocate(E(1:dos_nbins),stat=ierr)
    if (ierr/=0) call io_error ("cannot allocate E")

    delta_bins=(max_band_energy-min_band_energy)/real(dos_nbins-1,dp)
    do idos=1,dos_nbins
       E(idos) = min_band_energy+real(idos-1,dp)*delta_bins
    end do
    
    if(on_root.and.(iprint>2))then
       write(stdout,'(1x,a40,e11.5,a13)') 'dos_min_energy : ', dos_min_energy, " <-- DOS Grid" 
       write(stdout,'(1x,a40,e11.5,a13)') 'dos_max_energy : ', dos_max_energy, " <-- DOS Grid"
       write(stdout,'(1x,a40,f11.3,a13)') 'min_band_energy : ', min_band_energy, " <-- DOS Grid"
       write(stdout,'(1x,a40,f11.3,a13)') 'max_band_energy : ', max_band_energy, " <-- DOS Grid"
       write(stdout,'(1x,a40,i11,a13)')' dos_nbins : ', dos_nbins, " <-- DOS Grid"
       write(stdout,'(1x,a40,f11.3,a13)')' dos_spacing : ', dos_spacing, " <-- DOS Grid"
       write(stdout,'(1x,a40,f11.3,a13)')' delta_bins : ', delta_bins, " <-- DOS Grid"
    endif
    
  end subroutine setup_energy_scale


  !=============================================================================== 
  subroutine dos_utils_merge(dos, weighted_dos)
    !===============================================================================
    ! The DOS was calculated accross nodes. Now give them all back to root
    ! and free up the memeory on the slaves
    !------------------------------------------------------------------------------- 
    ! Arguments: dos          (in - slaves) (inout -  root)       : The DOS
    !            weighted_dos (in - slaves) (inout -  root) (opt) : Weighted DOS 
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: mw
    !-------------------------------------------------------------------------------
    ! Modules Used: See below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: None
    !-------------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------------
    ! Written by : A J Morris December 2010 
    !===============================================================================  
    use od_comms!,      only : on_root, comms_reduce
    use od_electronic, only : nspins
    use od_parameters, only : dos_nbins
    use od_io,         only : io_error 

    implicit none
    integer :: idos,ierr
    real(kind=dp),intent(inout), allocatable, optional :: weighted_dos(:,:,:) ! bins.spins, orbitals
    real(kind=dp),allocatable,intent(inout) :: dos(:,:)

    call comms_reduce(dos(1,1),nspins*dos_nbins,"SUM")

    if(present(weighted_dos))  call comms_reduce(weighted_dos(1,1,1),mw%nspins*dos_nbins*mw%norbitals,"SUM")

    if(.not.on_root) then 
       if(allocated(dos)) deallocate(dos,stat=ierr)
       if (ierr/=0) call io_error (" ERROR : dos : merge_dos : cannot deallocate dos")
       if(present(weighted_dos))  then
          if(allocated(weighted_dos)) deallocate(weighted_dos,stat=ierr)
          if (ierr/=0) call io_error (" ERROR : dos : merge_dos : cannot deallocate weighted_dos")
       end if
    endif
  end subroutine dos_utils_merge


  !===============================================================================
  subroutine allocate_dos(dos,intdos,w_dos)
    !===============================================================================
    ! Allocate the dos, intdos and w_dos if necessary. 
    !------------------------------------------------------------------------------- 
    ! Arguments: dos    (inout)       : The Density of States
    !            intdos (inout)       : The Integrated DOS
    !            w_dos  (inout) (opt) : Weighted DOS 
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: mw
    !-------------------------------------------------------------------------------
    ! Modules Used: See below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: None
    !-------------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------------
    ! Written by : A J Morris December 2010 
    !===============================================================================  
    use od_electronic, only : nspins
    use od_parameters, only : dos_nbins
    use od_io,         only : io_error

    implicit none

    real(kind=dp), allocatable, intent(inout)  :: dos(:,:)
    real(kind=dp), allocatable, intent(inout)  :: intdos(:,:)

    real(kind=dp),intent(out),optional,allocatable    :: w_dos(:,:,:) ! bins.spins, orbitals 

    integer :: ierr 

    allocate(dos(dos_nbins,nspins), stat=ierr)
    if(ierr/=0) call io_error("error in allocating dos")
    dos=0.0_dp

    allocate(intdos(dos_nbins,nspins), stat=ierr)
    if(ierr/=0) call io_error("error in allocating intdos")
    intdos=0.0_dp

    if(present(w_dos)) then
       allocate(w_dos(dos_nbins,mw%nspins,mw%norbitals), stat=ierr)
       if(ierr/=0) call io_error("error in allocating weighted_dos")
       w_dos=0.0_dp
    endif
  end subroutine allocate_dos


  subroutine dos_utils_deallocate
    if(allocated(dos_adaptive)) deallocate(dos_adaptive)
    if(allocated(dos_fixed)) deallocate(dos_fixed)
    if(allocated(dos_linear)) deallocate(dos_linear)
    if(allocated(intdos_adaptive)) deallocate(intdos_adaptive)
    if(allocated(intdos_fixed)) deallocate(intdos_fixed)
    if(allocated(intdos_linear)) deallocate(intdos_linear)
    if(allocated(E)) deallocate(E)
   end subroutine dos_utils_deallocate

  !===============================================================================
  subroutine calculate_dos(dos_type, dos, intdos, matrix_weights, weighted_dos)
    !===============================================================================
    ! Once everything is set up this is the main workhorse of the module.
    ! It accumulates the DOS and WDOS be looping over spins, kpoints and bands.
    !------------------------------------------------------------------------------- 
    ! Arguments: dos           (out)       : The Density of States
    !                          (in)        : one of "a", "f", "l", "q" 
    !            intdos        (out)       : The Integrated DOS
    !            matrix_weights(in)  (opt) : The weightings, such as LCAO for the 
    !                                        weighted dos.      
    !            weighted_dos  (out) (opt) : Weighted DOS 
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: delta_bins
    !-------------------------------------------------------------------------------
    ! Modules Used: See below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: None
    !-------------------------------------------------------------------------------
    ! Known Worries: We use local adaptive, fixed and linear variable so that if multiple 
    ! are set, it won't always appear to do the linear scheme.
    !-------------------------------------------------------------------------------
    ! Written by : A J Morris December 2010 Heavliy modified from LinDOS
    !===============================================================================
    use od_constants,  only : H2eV,sqrt_two
    use od_algorithms, only : gaussian, algorithms_erf
    use od_cell,       only : kpoint_grid_dim,nkpoints,kpoint_weight,num_kpoints_on_node
    use od_electronic, only : band_gradient, electrons_per_state, nbands,nspins,band_energy
    use od_parameters, only : adaptive_smearing,fixed_smearing&
         &,finite_bin_correction,iprint,dos_nbins,numerical_intdos 
    use od_io,         only : stdout,io_error
    use od_comms!,         only : my_node_id

    implicit none

    integer :: i,ik,is,ib,idos,ierr,iorb,dunit
    integer :: m,n,o,nn
    real(kind=dp) :: dos_temp, cuml, intdos_accum, width
    real(kind=dp) :: grad(1:3), step(1:3), EV(0:4)

    character(len=1), intent(in)                    :: dos_type
    real(kind=dp),intent(out),allocatable, optional :: weighted_dos(:,:,:)  
    real(kind=dp),intent(in),              optional :: matrix_weights(:,:,:,:)

    real(kind=dp),intent(out),allocatable :: dos(:,:), intdos(:,:)

    logical :: linear,fixed,adaptive

    linear=.false.
    fixed=.false.
    adaptive=.false.

    select case (dos_type)
       case ("l")
          linear=.true.
       case("a")
          adaptive=.true.
       case("f")
          fixed=.true.
       case default
          if (ierr/=0) call io_error (" ERROR : unknown dos_type in calculate_dos ")
    end select
    

    if(linear.or.adaptive) step(:) = 1.0_dp/real(kpoint_grid_dim(:),dp)/2.0_dp
    if(adaptive) adaptive_smearing=adaptive_smearing*sum(step(:))/3
    if(fixed) width=fixed_smearing

    if(calc_weighted_dos)then
       call allocate_dos(dos, intdos, w_dos=weighted_dos)
    else
       call allocate_dos(dos,intdos)
    endif

    do ik=1,num_kpoints_on_node(my_node_id)
       if(iprint>1.and.on_root) then
          if (mod(real(ik,dp),10.0_dp) == 0.0_dp) write(stdout,'(a40,i4,a3,i4,a14,21x,a7)') &
               &"Calculating k-point ", ik, " of", num_kpoints_on_node(my_node_id)," on this node.","<-- DOS"
       endif
       do is=1,nspins
          do ib=1,nbands
             if(linear.or.adaptive) grad(:) = real(band_gradient(ib,ib,:,ik,is),dp)
             if(linear) call doslin_sub_cell_corners(grad,step,band_energy(ib,is,ik),EV)
             if(adaptive) width = sqrt(dot_product(grad,grad))*adaptive_smearing
             ! Hybrid Adaptive -- This way we don't lose weight at very flat parts of the
             ! band. It's a kind of fudge that we wouldn't need if we had infinitely small bins.
             if(finite_bin_correction.and.(width<delta_bins)) width = delta_bins

             intdos_accum=0.0_dp

             do idos=1,dos_nbins
                ! The linear method has a special way to calculate the integrated dos
                ! we have to take account for this here.
                if(linear)then
                   dos_temp=doslin(EV(0),EV(1),EV(2),EV(3),EV(4),E(idos),cuml)*electrons_per_state*kpoint_weight(ik)
                   intdos(idos,is) = intdos(idos,is) + electrons_per_state*kpoint_weight(ik)*cuml
                else
                   dos_temp=gaussian(band_energy(ib,is,ik),width,E(idos))*electrons_per_state*kpoint_weight(ik)
                   if(numerical_intdos) then
                      intdos_accum=intdos_accum+dos_temp
                      intdos(idos,is)=intdos(idos,is)+intdos_accum
                   else ! Do it (semi)-analytically          
                      intdos(idos,is)=intdos(idos,is)+0.5_dp*(1.0_dp+algorithms_erf((E(idos)- &
                           & band_energy(ib,is,ik))/(sqrt_two*width)))*electrons_per_state*kpoint_weight(ik)
                   endif
                endif

                dos(idos,is)=dos(idos,is) + dos_temp

                if(calc_weighted_dos) then
                   if(ik.le.mw%nkpoints) then
                      if(ib.le.mw%nbands) then
                         do iorb=1,mw%norbitals
                            weighted_dos(idos,is,iorb)=weighted_dos(idos,is,iorb)+ &
                                 & dos_temp*matrix_weights(iorb,ib,ik,is)
                         enddo
                      endif
                   endif
                endif
             end do
          end do
       end do
    end do

    if(.not.linear.and.numerical_intdos) intdos=intdos*delta_bins

  end subroutine calculate_dos



  !===============================================================================
  subroutine doslin_sub_cell_corners(grad,step,energy,EigenV)
    !===============================================================================
    ! A helper subroutine for calculated_dos, which is used for the linear 
    ! broadening method. This routine extrapolates the energies at the corner of the 
    ! sub cells by using the gradient at the centre of the cell
    !------------------------------------------------------------------------------- 
    ! Arguments: grad   (in) : The Gradient of the band at the centre of a sub-cell
    !            step   (in) : The directions to the edge of the sub_cell
    !            energy (in) : The Band energy at the centre of the sub cell
    !            EigenV (out): The Energies at the corners of the sub-cell
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: None
    !-------------------------------------------------------------------------------
    ! Modules Used: See below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: None
    !-------------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------------
    ! Written by : A J Morris December 2010 Heavliy modified from LinDOS
    !===============================================================================
    use od_algorithms,only : heap_sort
    use od_constants, only : bohr2ang
    use od_cell,      only : recip_lattice

    implicit none
    integer :: m,n,o,nn,i
    real(kind=dp), intent(in) :: grad(1:3), step(1:3)
    real(kind=dp), intent(out) :: EigenV(0:4) 

    real(kind=dp) :: stepp(1:3),DE(1:8),energy

    nn = 0
    do m=-1,1,2
       do n=-1,1,2
          do o=-1,1,2
             nn=nn+1
             stepp(1) = step(1)*m
             stepp(2) = step(2)*n
             stepp(3) = step(3)*o

             ! Reciprocal lattice in inverse Ang
             stepp = matmul(recip_lattice,stepp)

             ! DE in eV
             DE(nn)=dot_product(grad,stepp)

          end do
       end do
    end do

    ! Yes this is a hack to the sorter to work the right way around.
    DE=-DE
    call heap_sort(8,DE)
    DE=-DE

    ! WHY ARE WE STORING ALL THIS?
    !EV(0,ib,is,ik) = band_energy(ib,is,ik)
    !EV(1,ib,is,ik) = EV(0,ib,is,ik) + DE(5)
    !EV(2,ib,is,ik) = EV(0,ib,is,ik) + DE(6)
    !EV(3,ib,is,ik) = EV(0,ib,is,ik) + DE(7)
    !EV(4,ib,is,ik) = EV(0,ib,is,ik) + DE(8)
    EigenV(0) = Energy
    do i=1,4
       EigenV(i) = EigenV(0) + DE(4+i)
    enddo
  end subroutine doslin_sub_cell_corners

  !===============================================================================
  function doslin(e0,e1,e2,e3,e4,e,int)
    !===============================================================================
    ! Return the DoS contribution for a linear band portion and a cubic cell
    !------------------------------------------------------------------------------- 
    ! Arguments: e0 (in) : Energy at centre of sub cell
    !   e1,e2,e3,e4 (in) : Energies of the four lowest corners of the sub cell
    !            e  (in) : Energy at which DOS is evaluated
    !          int  (out): Integrated DOS contribution for energy, E)
    ! (The function itself returns the DOS couribution for energy, E)
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: None
    !-------------------------------------------------------------------------------
    ! Modules Used: See below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: None
    !-------------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------------
    ! Written by : C J Pickard. From LinDOS. Extra Comments A J Morris Sept 2010
    !===============================================================================
    implicit none

    ! - E X T R A   C O M M E N T S   B Y   A J M
    real (kind=dp), intent(in)  :: e0,e1,e2,e3,e4
    real (kind=dp), intent(in)  :: e
    real (kind=dp), intent(out) :: int

    real (kind=dp) :: doslin
    real (kind=dp) :: alpha,e_use
    logical        :: shift

    ! ** Check the input arguments

    ! The inputs must be in ascending order. These are the four lowest corners of the
    ! cube (sub cell) around the k-point. E_0 is the energy of the k-point which would
    ! in a Gaussian smearing scheme, be the only info we have about this cube of recip
    ! space.
    if(.not.((e4.le.e3).and.(e3.le.e2).and.(e2.le.e1).and.(e1.le.e0))) then
       stop 'doslin: input arguments incorrect'
    end if

    ! ** Return if outside the range

    ! If the energy is below the smallest corner, then return. The CES doesn't cut this
    ! sub cell
    if(e.le.e4) then
       doslin = 0.0_dp ! TRUE
       int    = 0.0_dp ! TRUE
       return
    end if

    ! Since the energy in the sub cell is linearly extrapolated, if the energy is twice as
    ! big as the energy difference between the smallest corner and the middle, then this 
    ! energy is outdside the top of the cell, and whilst this cell doesn't contribute to the
    ! CES, it does to the occupation.
    if(e.ge.(2.0_dp*e0-e4)) then
       doslin = 0.0_dp ! TRUE
       int    = 1.0_dp ! TRUE
       return
    end if

    ! ** Special treatment if all vertices at the same energy
    ! If the CES is perfectly flat in the cell, then we don't want to be dividing by zero.
    ! the below just catches this case and forces alpha to be 1. This is fine as the else 
    ! block below looks for the same problem and the answer comes out correctly to 
    if((e1==e2).and.(e1==e3).and.(e3==e4)) then
       alpha = 1.0_dp
    else
       alpha = 1.0_dp/((e1-e3)*(e1-e4)+(e3-e4)**2/3.0_dp-(e1-e2)**2/3.0_dp+(e1+e2-e3-e4)*(e0-e1)) ! TRUE
    endif

    ! ** Flip if above e0

    ! If e is greater than the energy of the k-point, then we're going to subtract the
    ! contribution from a full cell, rather than add it to an empty one. The extrapolation
    ! is linear, so we can do this fine.
    if(e.gt.e0) then
       e_use = 2.0_dp*e0-e ! TRUE
       shift = .true.
    else
       e_use = e ! TRUE
       shift = .false.
    end if

    ! ** The analytic constributions to the DOS and integrated DOS

    if(e_use.le.e4) then
       ! P O S S I B L E   S P E E D   U P
       ! If we ended up in here, something went wrong as this should already have been trapped.
       doslin = 0.0_dp ! TRUE
       int    = 0.0_dp ! TRUE
    else if(e_use.lt.e3) then
       ! There isn't a problem with divide by zero here. Since if e3=e4 and e_use < e3
       ! we would have been caugh in the above if.
       doslin = (e_use-e4)**2/(e3-e4)/2.0_dp ! TRUE 
       int    = (e_use-e4)**3/(e3-e4)/6.0_dp ! TRUE
    else if(e_use.lt.e2) then
       doslin = (e_use-(e3+e4)/2.0_dp) ! TRUE 
       int    = ((e3-e4)**2/3.0_dp+(e_use-e3)*(e_use-e4))/2.0_dp ! TRUE
    else if(e_use.lt.e1) then
       doslin = (e1+e2-e3-e4)/2.0_dp ! TRUE
       ! Ok, so the IF costs more than the maths. But this way also catches the
       ! divide by zero. Clever!
       if(e1.ne.e2) doslin = doslin - (e1-e_use)**2/(e1-e2)/2.0_dp ! TRUE
       int    = ((e2-e4)*(e_use-e3)+(e1-e3)*(e_use-e2)+(e3-e4)**2/3.0_dp+&
            ((e1-e_use)**3-(e1-e2)**3)/(e1-e2)/3.0_dp)/2.0_dp ! TRUE
    else if(e_use.le.e0) then
       ! Check to see if the band is flat.
       if((e1+e2-e3-e4).gt.0.0_dp) then
          doslin = (e1+e2-e3-e4)/2.0_dp ! TRUE
          int    = ((e1-e3)*(e1-e4)+(e3-e4)**2/3.0_dp-(e1-e2)**2/3.0_dp+&
               (e1+e2-e3-e4)*(e_use-e1))/2.0_dp ! TRUE
       else ! This can only happen if e1=e2=e3=e4,
          ! in this case we stop doing what we were doing and calculate the contirbution from the 
          ! gradient simplistically. grad E = e0-e4. 
          doslin = 1.0_dp/(e0-e4)/2.0_dp ! TRUE
          int    = (e_use-e4)/(e0-e4)/2.0_dp ! SO YES, THIS DOES INTEGRATE FROM THE LINE ABOVE
       end if
    else
       write (*,*) e_use,e
       write (*,*) e0,e1,e2,e3,e4
       stop 'Got here, but not supposed to!'
    end if

    ! ** Normalise

    doslin = doslin*alpha

    if(shift) then
       int    = 1.0_dp - int*alpha ! TRUE
    else
       int    = int*alpha ! TRUE
    end if

    return
  end function doslin
  
   !=============================================================================== 
  subroutine dos_utils_calculate_at_e(energy, matrix_weights, weighted_dos_at_e, dos_at_e)
    !===============================================================================  
    ! Main routine in dos module, drives the calculation of density of states for
    ! both task : dos and also if it is required elsewhere.
    !------------------------------------------------------------------------------- 
    ! Arguments: matrix_weigths (in) (opt) : LCAO or other weightings for DOS
    !            weighted_dos   (out)(opt) : Output DOS weigthed by matrix_weights
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: mw, E, dos_adaptive, dos_fixed, dos_linear
    ! intdos_adaptive, intdos_fixed, intdos_linear, efermi_fixed, efermi_adaptive
    ! efermi_linear, delta_bins, calc_weighted_dos
    !-------------------------------------------------------------------------------
    ! Modules Used: see below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: One of linear, adaptive or fixed must be .true.
    !-------------------------------------------------------------------------------
    ! Known Worries: (1) If more than one of linear, adaptive or fixed are set it 
    ! uses the most complicated method.
    ! (2) It should be possible to pass optioinal arguments to sub programs as
    ! optional argumnets without checking whether they are there or not. g95 will 
    ! allow this behaviour. gfotran will not.
    !-------------------------------------------------------------------------------
    ! Written by : A J Morris December 2010
    !=============================================================================== 
    use od_io,        only : stdout,io_time,io_error
    use od_comms,     only : on_root,my_node_id
    use od_electronic,only : band_gradient,band_energy, efermi, efermi_castep,nspins, &
         & elec_read_band_gradient,elec_read_band_energy,elec_report_parameters
    use od_parameters,only : linear, adaptive, fixed, quad, compute_band_energy, &
         & compute_efermi,dos,dos_per_volume,fermi_energy,iprint
    use od_cell,         only : cell_volume, nkpoints, cell_calc_lattice, &
         & cell_report_parameters,cell_dist,kpoint_grid_dim
    implicit none

    !-------------------------------------------------------------------------------
    ! I N T E R N A L   V A R I A B L E S
    integer :: ierr, idos, i, ik, is, ib
    real(kind=dp) :: time0, time1
    real(kind=dp),intent(in), allocatable, optional  :: matrix_weights(:,:,:,:)
    real(kind=dp),intent(out),allocatable, optional  :: weighted_dos_at_e(:,:) ! spins, orbitals
    real(kind=dp),intent(in) :: energy
    real(kind=dp),intent(out) :: dos_at_e(nspins) ! spins
    !-------------------------------------------------------------------------------

    if(.not.(linear.or.adaptive.or.fixed.or.quad)) call io_error (" DOS: No Broadening Set")
    
    
    calc_weighted_dos=.false.
    if(present(matrix_weights)) calc_weighted_dos=.true.
    
    if(calc_weighted_dos)then
       mw%norbitals=size(matrix_weights,1)
       mw%nbands   =size(matrix_weights,2)
       mw%nkpoints =size(matrix_weights,3)
       mw%nspins   =size(matrix_weights,4)
    end if


    if(calc_weighted_dos) then 
       if(mw%nspins.ne.nspins)     call io_error ("ERROR : DOS :  mw%nspins not equal to nspins.")
       if(mw%nkpoints.ne.nkpoints) call io_error ("ERROR : DOS : mw%nkpoints not equal to nkpoints.")
    endif

    !-------------------------------------------------------------------------------
    ! R E A D   B A N D   G R A D I E N T S 
    ! If we're using one of the more accurate roadening schemes we also need to read in the 
    ! band gradients too
    if(quad.or.linear.or.adaptive) then
      if(.not.allocated(band_gradient)) call elec_read_band_gradient
    endif
    !-------------------------------------------------------------------------------
    ! C A L C U L A T E   D O S 
    ! Now everything is set up, we can perform the dos accumulation in parellel
    
    time0=io_time()

    if(on_root.and.(iprint>1)) write(stdout,*)

 
       if(calc_weighted_dos)then 
          call calculate_dos_at_e(energy,dos_at_e, matrix_weights=matrix_weights, &
&weighted_dos_at_e=weighted_dos_at_e) 
          call dos_utils_merge_at_e(dos_at_e,weighted_dos_at_e=weighted_dos_at_e)    
       else
          call calculate_dos_at_e(energy,dos_at_e)
          call dos_utils_merge_at_e(dos_at_e) 
       endif
    
    time1=io_time()
    if(on_root)  write(stdout,'(1x,a40,f11.3,a)') 'Time to calculate dos  ',time1-time0,' (sec)'
    !-------------------------------------------------------------------------------

  end subroutine dos_utils_calculate_at_e
  
  
  !===============================================================================
  subroutine calculate_dos_at_e(energy, dos_at_e, matrix_weights, weighted_dos_at_e)
    !===============================================================================
    ! Once everything is set up this is the main workhorse of the module.
    ! It accumulates the DOS and WDOS be looping over spins, kpoints and bands.
    !------------------------------------------------------------------------------- 
    ! Arguments: dos           (out)       : The Density of States
    !            intdos        (out)       : The Integrated DOS
    !            matrix_weights(in)  (opt) : The weightings, such as LCAO for the 
    !                                        weighted dos.      
    !            weighted_dos  (out) (opt) : Weighted DOS 
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: delta_bins
    !-------------------------------------------------------------------------------
    ! Modules Used: See below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: None
    !-------------------------------------------------------------------------------
    ! Known Worries: If linear, fixed and adaptive are all set. It produces the 
    ! linear result, no matter what was intended. This would need to be modified to
    ! take an optinal argument, which would force only one of the above to be set
    ! within the subroutine. Since linear, fixed and adaptive are only non-mutually
    ! exculsive when debugging (such things as efermi) this isn't a priority
    !-------------------------------------------------------------------------------
    ! Written by : A J Morris December 2010 Heavliy modified from LinDOS
    !===============================================================================
    use od_constants,  only : H2eV,sqrt_two
    use od_algorithms, only : gaussian
    use od_cell,       only : kpoint_grid_dim,nkpoints,kpoint_weight,num_kpoints_on_node
    use od_electronic, only : band_gradient, electrons_per_state, nbands,nspins,band_energy
    use od_parameters, only : linear,fixed,adaptive,adaptive_smearing,fixed_smearing&
         &,finite_bin_correction,iprint,dos_nbins,numerical_intdos 
    use od_io,         only : stdout,io_error
    use od_comms!,         only : my_node_id

    implicit none

    integer :: i,ik,is,ib,idos,ierr,iorb,dunit
    integer :: m,n,o,nn
    real(kind=dp) :: dos_temp, cuml, intdos_accum, width
    real(kind=dp) :: grad(1:3), step(1:3), EV(0:4)

    real(kind=dp),intent(out),allocatable, optional :: weighted_dos_at_e(:,:)  
    real(kind=dp),intent(in),              optional :: matrix_weights(:,:,:,:)
    real(kind=dp),intent(in) :: energy 
    real(kind=dp),intent(out) :: dos_at_e(nspins)

    if(linear.or.adaptive) step(:) = 1.0_dp/real(kpoint_grid_dim(:),dp)/2.0_dp
    if(adaptive) adaptive_smearing=adaptive_smearing*sum(step(:))/3
    if(fixed) width=fixed_smearing

    weighted_dos_at_e=0.0_dp

    do ik=1,num_kpoints_on_node(my_node_id)
       if(iprint>1.and.on_root) then
          if (mod(real(ik,dp),10.0_dp) == 0.0_dp) write(stdout,'(a40,i4,a3,i4,a14,21x,a7)') &
               &"Calculating k-point ", ik, " of", num_kpoints_on_node(my_node_id)," on this node.","<-- DOS"
       endif
       do is=1,nspins
          do ib=1,nbands
             if(linear.or.adaptive) grad(:) = real(band_gradient(ib,ib,:,ik,is),dp)
             if(linear) call doslin_sub_cell_corners(grad,step,band_energy(ib,is,ik),EV)
             if(adaptive) width = sqrt(dot_product(grad,grad))*adaptive_smearing
             ! Hybrid Adaptive -- This way we don't lose weight at very flat parts of the
             ! band. It's a kind of fudge that we wouldn't need if we had infinitely small bins.
             if(finite_bin_correction.and.(width<delta_bins)) width = delta_bins

             intdos_accum=0.0_dp

                ! The linear method has a special way to calculate the integrated dos
                ! we have to take account for this here.
                if(linear)then
                   dos_temp=doslin(EV(0),EV(1),EV(2),EV(3),EV(4),energy,cuml)*electrons_per_state*kpoint_weight(ik)
                   
                else
                   dos_temp=gaussian(band_energy(ib,is,ik),width,energy)*electrons_per_state*kpoint_weight(ik)
                   
                endif

                dos_at_e(is)=dos_at_e(is)+dos_temp

                if(calc_weighted_dos) then
                   if(ik.le.mw%nkpoints) then
                      if(ib.le.mw%nbands) then
                         do iorb=1,mw%norbitals
                            weighted_dos_at_e(is,iorb)=weighted_dos_at_e(is,iorb)+ &
                                 & dos_temp*matrix_weights(iorb,ib,ik,is)
                         enddo
                      endif
                   endif
                endif
             end do
       end do
    end do

  end subroutine calculate_dos_at_e
  
    !=============================================================================== 
  subroutine dos_utils_merge_at_e(dos, weighted_dos_at_e)
    !===============================================================================
    ! The DOS was calculated accross nodes. Now give them all back to root
    ! and free up the memeory on the slaves
    !------------------------------------------------------------------------------- 
    ! Arguments: dos          (in - slaves) (inout -  root)       : The DOS
    !            weighted_dos (in - slaves) (inout -  root) (opt) : Weighted DOS 
    !-------------------------------------------------------------------------------
    ! Parent Module Varables Used: mw
    !-------------------------------------------------------------------------------
    ! Modules Used: See below
    !-------------------------------------------------------------------------------
    ! Key Internal Variables: None
    !-------------------------------------------------------------------------------
    ! Necessary Conditions: None
    !-------------------------------------------------------------------------------
    ! Known Worries: None
    !-------------------------------------------------------------------------------
    ! Written by : A J Morris December 2010 
    !===============================================================================  
    use od_comms!,      only : on_root, comms_reduce
    use od_electronic, only : nspins
    use od_parameters, only : dos_nbins
    use od_io,         only : io_error 

    implicit none
    integer :: idos,ierr
    real(kind=dp),intent(inout), allocatable, optional :: weighted_dos_at_e(:,:) ! bins.spins, orbitals
    real(kind=dp),intent(inout) :: dos(nspins)

    call comms_reduce(dos(1),nspins,"SUM")

    if(present(weighted_dos_at_e))  call comms_reduce(weighted_dos_at_e(1,1),mw%nspins*mw%norbitals,"SUM")

    if(.not.on_root) then 
       if(present(weighted_dos_at_e))  then
          if(allocated(weighted_dos_at_e)) deallocate(weighted_dos_at_e,stat=ierr)
          if (ierr/=0) call io_error (" ERROR : dos : merge_dos : cannot deallocate weighted_dos")
       end if
    endif
  end subroutine dos_utils_merge_at_e

endmodule od_dos_utils