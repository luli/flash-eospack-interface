!!****if* source/physics/Eos/EosMain/multiTemp/Eospac/eos_eospacData
!!
!! NAME
!!
!!  eos_eospacData
!!
!! 
!! SYNOPSIS
!!
!!  use eos_eospacData
!!
!! DESCRIPTION
!!
!!  This is the data module for the Gamma law Eos implementation with 
!!  multiple fluids/species, for (up to) 3 temperatures.
!!  DEV: Gamma is assumed to be the same for all species !?
!!  It stores all the runtime parameters, and all the unit scope
!!  data. Some of the unit scope data is fecthed by the wrapper layer
!!  from elsewhere in the code and some is local unit data common to
!!  multiple functions in the unit.
!! 
!! PARAMETERS
!!  
!!   These are the runtime parameters used in Multigamma Eos.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory.
!!   You might have over written these values with the flash.par values
!!   for your specific run.  
!!
!!   smalle[Real]  --- the smallest value for energy 
!!
!!***

module eos_eospacData

#include "Flash.h"
#include "Eos_eospac.h"

  use libeospac_interface

  public :: eospacIntVect2Arr, eospacIntArr2Vect

  real, save :: eos_gammaEle
  real, save :: eos_gammam1Ele, eos_gammam1Rad

  real, save :: eos_eMass
  real, save :: eos_eMassInUAmu !electron mass in unified atomic mass units (aka daltons)

  real, save, dimension(NSPECIES) ::  eos_gc, eos_gammam1j, eos_ggprodj, eos_ggprodinvj, eos_gam1invj

  real, save :: eos_maxFactorUp
  real, save :: eos_maxFactorDown

  ! Some of these values can probaly be taked from the Physical constants unit,
  ! but it's simpler to define them statically rather then use
  ! Physicalconstants_get(...) 
  real(EOS_REAL), parameter, private  :: d_ses2cgs = 1._EOS_REAL    ! kg/m3  -> g/cm3 # doens't look like it's needed
  real(EOS_REAL), parameter, private  :: t_ses2cgs = 1._EOS_REAL       ! K      -> K
  real(EOS_REAL), parameter, private  :: p_ses2cgs = 1.0e10_EOS_REAL   ! GPa    -> Barye
  real(EOS_REAL), parameter, private  :: u_ses2cgs = 1.0e10_EOS_REAL   ! MJ/kg  -> erg/g
  real(EOS_REAL), parameter, private  :: s_ses2cgs = 1.0e10_EOS_REAL   ! MJ/kg/K-> erg/g/K

  ! All the *Arr and *Vec objects contain the same data in reshaped arrays.
  ! Both should be synced by the end of eos_initEospac.F90.
  ! They should be accessed in read-only mode in other routines.
  integer(EOS_INTEGER), save :: eos_tableTypeArr(NSPECIES,NEOSTABS)
  integer(EOS_INTEGER), save :: eos_tableTypeVec(NSPECIES*NEOSTABS)

  integer(EOS_INTEGER), save :: eos_matIdArr(NSPECIES,NEOSTABS)
  integer(EOS_INTEGER), save :: eos_matIdVec(NSPECIES*NEOSTABS)

  integer(EOS_INTEGER), save :: eos_tableHandleArr(NSPECIES,NEOSTABS)
  integer(EOS_INTEGER), save :: eos_tableHandleVec(NSPECIES*NEOSTABS)

  integer(EOS_INTEGER), save :: eos_matId(NSPECIES) ! contains SESAME table name

  integer(EOS_INTEGER), save :: eos_tableTypeBase(NEOSTABS) = (/ &
      EOS_Pic_DT, &
      EOS_Uic_DT, &
      EOS_Sic_DT, &
      EOS_Pe_DT, &
      EOS_Ue_DT, &
      EOS_Se_DT,  &
      EOS_Pic_DUic, &
      EOS_T_DUic, &
      EOS_Sic_DUic, &
      EOS_Pe_DUe, &
      EOS_T_DUe, &
      EOS_Se_DUe  &
      /)
  ! This array takes care of the conversions from SESAME units to cgs.
  ! A strong assumption here is that EOS_X_CONVERT, EOS_Y_CONVERT, EOS_F_CONVERT 
  ! defined inside Eospac's eos_Interface.h have
  ! consecutive integer values that start with EOS_X_CONVERT (it is the case
  ! with EOSPAC v6.2.4beta)
  real(EOS_REAL), parameter :: eos_convertSesame2Cgs(NEOSTABS,EOS_X_CONVERT:EOS_X_CONVERT+2) = &
      transpose(reshape( (/ &
         d_ses2cgs, t_ses2cgs, p_ses2cgs, &
         d_ses2cgs, t_ses2cgs, u_ses2cgs, &
         d_ses2cgs, t_ses2cgs, s_ses2cgs, &
         d_ses2cgs, t_ses2cgs, p_ses2cgs, &
         d_ses2cgs, t_ses2cgs, u_ses2cgs, &
         d_ses2cgs, t_ses2cgs, s_ses2cgs, &
         d_ses2cgs, 1./u_ses2cgs, p_ses2cgs, &
         d_ses2cgs, 1./u_ses2cgs, t_ses2cgs, &
         d_ses2cgs, 1./u_ses2cgs, s_ses2cgs, &
         d_ses2cgs, 1./u_ses2cgs, p_ses2cgs, &
         d_ses2cgs, 1./u_ses2cgs, t_ses2cgs, &
         d_ses2cgs, 1./u_ses2cgs, s_ses2cgs /),&
         (/3, NEOSTABS /)))

  real(EOS_REAL), save :: eos_massFraction(NSPECIES) ! temporarly should be read from FLASH

  integer(EOS_INTEGER), save :: eos_errorCode
  integer(EOS_INTEGER), save :: eos_tableHandleErrorCode
  real(EOS_REAL) :: eos_infoVals(NEOSINFO)
  integer(EOS_INTEGER) :: eos_infoItems(NEOSINFO) = (/ &
       EOS_Mean_Atomic_Num, &
       EOS_Mean_Atomic_Mass, &
       EOS_Normal_Density, &
       EOS_X_Convert_Factor, &
       EOS_Y_Convert_Factor,&
       EOS_F_Convert_Factor &
       /)
  logical, save :: eos_optionFlags(EOS_MIN_OPTION_FLAG_VALUE:EOS_NUM_TABLE_OPTIONS+EOS_MIN_OPTION_FLAG_VALUE) = .false.

  character(EOS_MaxErrMsgLen) :: eos_errorMessage


contains
  subroutine eospacIntVect2Arr(vect,arr)
    integer(EOS_INTEGER), intent(in) :: vect(NSPECIES*NEOSTABS)
    integer(EOS_INTEGER), intent(out) :: arr(NSPECIES,NEOSTABS)

    integer :: i,j
    
    do i=1,NEOSTABS
      do j=1,NSPECIES
        arr(j,i) =  vect((i-1)*NSPECIES+j)
      enddo
    enddo

  end subroutine eospacIntVect2Arr

  subroutine eospacIntArr2Vect(arr,vect)
     integer(EOS_INTEGER), intent(out) :: vect(NSPECIES*NEOSTABS)
     integer(EOS_INTEGER), intent(in) :: arr(NSPECIES,NEOSTABS)

     integer :: i,j

     do i=1,NEOSTABS
        do j=1,NSPECIES
           vect((i-1)*NSPECIES+j) = arr(j,i)
        enddo
     enddo

  end subroutine eospacIntArr2Vect

  function eos_indexUnwrap(idx) result(idx_arr)
     integer, intent(in) :: idx
     integer             :: idx_arr(2)

     idx_arr = (/ -1, -1 /)
     ! This should be done with modulo function...
     !idx_arr(1) = modulo(idx, NSPECIES)
     !dx_arr(2) = (idx - idx_arr(1))/NSPECIES + 1
      do i=1,NEOSTABS
          do j=1,NSPECIES
             if (idx == (i-1)*NSPECIES+j) then
                 idx_arr(1) = j
                 idx_arr(2) = i
             endif
          enddo
      enddo

  end function eos_indexUnwrap

  function eos_indexWrap(idx_arr) result(idx)
     integer             :: idx
     integer, intent(in) :: idx_arr(2)

     idx = (idx_arr(2)-1)*NSPECIES+ idx_arr(1)

  end function eos_indexWrap

  function eos_idxVar(var) result(idx)
     integer(EOS_INTEGER), intent(in) :: var
     integer                          :: idx
     idx = 1
     do while (eos_tableTypeBase(idx)/=var.and.idx<=NEOSTABS)
        idx = idx + 1
     end do
  end function eos_idxVar


end module eos_eospacData
