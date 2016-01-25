!!****if* source/physics/Eos/EosMaiGn/multiTemp/Eospac/eos_initEospac
!!
!! NAME
!!
!!  eos_initEospac
!!
!! SYNOPSIS
!!
!!  call eos_initEospac()
!!
!! DESCRIPTION
!!
!!  This routine initializes various scalars and arrays needed
!!  by the EOS unit from the runtime parameters and physical constants.
!!  This version is for use when multiple species
!!  are present. The gamma values for different species are obtained from
!!  the Multispecies unit, initialized in Simulation_initSpecies.F90 .
!!
!! PARAMETERS
!!  
!!   These are the runtime parameters used in Gamma law Eos for multiple
!!   species with different abundances. 
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory. You might over write these values with the 
!!   flash.par values for your specific run.  
!!
!!   eos_singleSpeciesA[Real]  -- Nucleon number for the gas, default 1.00
!!   eos_singleSpeciesZ[Real]  -- Proton number for the gas, default 1.00
!!   eos_tolerance[Real, 1.0e-8] -- Convergence tolerance for the Newton-Rhapson
!!               iterations
!!   eos_maxNewton[Integer, 50] -- Maximum number of Newton-Raphson iterations
!!   eos_forceConstantInput     -- This switch forces the Eos implementation
!!                                 to never modify EINT or PRES in MODE_DENS_EI
!!                                 and MODE_DENS_PRES and similar modes.
!!                                 If this is .false. (the
!!                                 default), calls to Eos may slightly modify
!!                                 these input variables in order to preserve
!!                                 thermodynamic equilibrium.
!!
!!  NOTES
!!
!!  Gamma law Eos defines two mesh-based parameters GAMC_VAR and GAME_VAR in Flash.h
!!
!!***

subroutine eos_initEospac()

  use Simulation_interface, ONLY : Simulation_mapIntToStr
  use Eos_data, ONLY : eos_gasConstant, eos_smalle, eos_eintSwitch, eos_type, &
       eos_tol, eos_maxNewton, eos_forceConstantInput, &
       eos_combinedTempRule, eos_entrEleScaleChoice
  use eos_mgammaData, ONLY : eos_gc, eos_gammam1j, eos_ggprodj, eos_ggprodinvj, eos_gam1invj, &
       eos_eMass, eos_eMassInUAmu, &
       eos_gammaEle, eos_gammam1Ele, eos_gammam1Rad, &
       eos_maxFactorUp, eos_maxFactorDown


  use eos_eospacData, ONLY : eos_tableTypeArr, eos_tableTypeVec,&
       eos_tableTypeBase, eos_matId, eos_massFraction, eos_matIdArr,&
       eos_matIdVec, eos_tableHandleArr, eos_tableHandleVec, &
       eos_errorCode, eos_tableHandleErrorCode, eos_infoVals,&
       eos_infoItems, eos_errorMessage, eos_optionFlags, &
       eospacIntVect2Arr, eospacIntArr2Vect, eos_indexUnwrap, &
       eos_indexWrap, eos_convertSesame2Cgs


  use RuntimeParameters_interface, ONLY : RuntimeParameters_get,&
                                      RuntimeParameters_set
  use Multispecies_interface, ONLY : Multispecies_getProperty,&
                                      Multispecies_setProperty
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get
  use eos_mtInterface, ONLY: eos_initPhysData
  use Logfile_interface, ONLY: Logfile_stampMessage
  use Driver_interface, ONLY : Driver_abortFlash

  use libeospac_interface
 


#include "Flash.h"
#include "Eos.h"
#include "Multispecies.h"
#include "Eos_eospac.h"
#include "constants.h"

  implicit none

  integer(EOS_INTEGER),parameter :: nXYPairs = 4
  real(EOS_REAL) :: X(nXYPairs), Y(nXYPairs), F(nXYPairs), dFx(nXYPairs), dFy(nXYPairs)

  real(EOS_REAL), save :: eos_concInMix(NSPECIES*nXYPairs)
  integer(EOS_INTEGER) ::  i, j
  integer :: spec, specno, matid_int, eosinfono, tabno, infono, optno_rel,&
              optno
  integer :: tabno_arridx(2)
  character(len=20)   :: matid_str
  character(len=300)  :: abortMessage, logMessage
  character(len=4)    :: speciesName
  real                :: propertyVal



  call RuntimeParameters_get("gammaEle", eos_gammaEle)
  call RuntimeParameters_get('eos_tolerance', eos_tol)
  call RuntimeParameters_get('eos_maxNewton', eos_maxNewton)
  call RuntimeParameters_get('eos_maxFactorUp', eos_maxFactorUp)
  call RuntimeParameters_get('eos_maxFactorDown', eos_maxFactorDown)

#ifndef EINT_VAR
  if (eos_eintSwitch > 0.0) then
     call Driver_abortFlash("[Eos_init] eintSwitch is nonzero, but EINT_VAR not defined!")
  end if
#endif

#ifdef OLD_EINTNSWITCH
  call RuntimeParameters_get("eint1Switch",eos_eint1Switch)  
  call RuntimeParameters_get("eint2Switch",eos_eint2Switch)  
  call RuntimeParameters_get("eint3Switch",eos_eint3Switch)  
  if (eos_eint1Switch==-1.0) eos_eint1Switch = eos_eintSwitch
  if (eos_eint2Switch==-1.0) eos_eint2Switch = eos_eintSwitch
  if (eos_eint3Switch==-1.0) eos_eint3Switch = eos_eintSwitch
#endif
  call RuntimeParameters_get("eos_forceConstantInput",eos_forceConstantInput)
  call RuntimeParameters_get("eos_combinedTempRule",eos_combinedTempRule)
  call RuntimeParameters_get("eos_entrEleScaleChoice", eos_entrEleScaleChoice)
  call PhysicalConstants_get("electron mass",eos_eMass) !or value from eos_helmConstData?

  call PhysicalConstants_get("electron mass",eos_eMassInUAmu,unitMass="amu")

  ! EOSPAC options
  call RuntimeParameters_get("eospac_pt_smoothing",eos_optionFlags(EOS_PT_SMOOTHING))
  call RuntimeParameters_get("eospac_smooth",eos_optionFlags(EOS_SMOOTH))
  call RuntimeParameters_get("eospac_create_tzero",eos_optionFlags(EOS_CREATE_TZERO))



  ! Compute radiation constant, maybe other physical constants:
  call eos_initPhysData()

  eos_gammam1Ele = 1.0/(eos_gammaEle-1.0)
  eos_gammam1Rad = 3.0
  

  do specno = 1, NSPECIES
     call Multispecies_getProperty(SPECIES_BEGIN + specno - 1, GAMMA,  eos_gc(specno))
  end do

  ! Note that these are all ARRAYS of size NSPECIES
  eos_gammam1j   = 1. / (eos_gc - 1.)
  eos_ggprodj    = eos_gammam1j * eos_gasConstant
  eos_ggprodinvj = 1. / eos_ggprodj
  eos_gam1invj   = 1. / eos_gammam1j
  eos_type = EOS_EOSPAC

  ! Reading user defined SESAME table id
  do specno=1,NSPECIES
     spec = SPECIES_BEGIN - 1 + specno
     call Multispecies_getProperty (spec , MS_EOSPRESFILE , matid_str)
     read( matid_str, '(i10)' ) matid_int
     eos_matId(specno) = matid_int
  enddo

  do specno=1,NSPECIES
     eos_matIdArr(specno,:) = eos_matId(specno)
  enddo

  do tabno=1,NEOSTABS
     eos_tableTypeArr(:,tabno) = eos_tableTypeBase(tabno)
  enddo

  do specno=1,NSPECIES
     do i=1,nXYPairs
        eos_concInMix(i+(specno-1)*nXYPairs) = eos_massFraction(specno)
     enddo
  enddo
    

  eos_errorCode = EOS_OK
  eos_tableHandleVec(:) = 0

  ! Making sure that all *Vec are contain the same thing as *Arr
  call eospacIntArr2Vect(eos_tableTypeArr, eos_tableTypeVec)
  call eospacIntArr2Vect(eos_matIdArr, eos_matIdVec)
 
  ! Initialize table data objects
  call eos_CreateTables ( NSPECIES*NEOSTABS, eos_tableTypeVec, eos_matIdVec, eos_tableHandleVec, eos_errorCode)
  if (eos_errorCode.NE.EOS_OK) then
     do tabno=1, NSPECIES*NEOSTABS
        eos_tableHandleErrorCode = EOS_OK
        call eos_GetErrorCode (eos_tableHandleVec(tabno), eos_tableHandleErrorCode )
        call eos_GetErrorMessage ( eos_tableHandleErrorCode, eos_errorMessage )
        write(*,998) 'eos_CreateTables ERROR ', eos_tableHandleErrorCode, ': ', &
                     eos_errorMessage(1:(len_trim(eos_errorMessage)-1))
     enddo
  endif
  call eospacIntVect2Arr(eos_tableHandleVec, eos_tableHandleArr)


  ! This part sets table options before loading the tables 
  ! Also changing output units to cgs.
  do optno_rel=1, EOS_NUM_TABLE_OPTIONS 
     optno = EOS_MIN_OPTION_FLAG_VALUE + optno_rel - 1
     do tabno=1, NSPECIES*NEOSTABS
        if (eos_optionFlags(optno))  then
           ! User defined binary options
           call eos_SetOption( eos_tableHandleVec(tabno), optno, &
                                  EOS_NullVal, eos_errorCode )
        else if (optno == EOS_X_CONVERT .or. optno == EOS_X_CONVERT &
               .or. optno == EOS_F_CONVERT) then
           ! Converting from Sesame units to cgs
           tabno_arridx = eos_indexUnwrap(tabno)
           call eos_SetOption( eos_tableHandleVec(tabno), optno,  &
                      eos_convertSesame2Cgs(tabno_arridx(2),optno), &
                      eos_errorCode)
        endif

        if (eos_errorCode.NE.EOS_OK) then
           call eos_GetErrorMessage ( eos_errorCode, eos_errorMessage )
           write(abortMessage,998) 'SetOption ERROR ', eos_errorCode, ': ', &
                        eos_errorMessage(1:(len_trim(eos_errorMessage)-1))
           if (eos_errorCode/=EOS_INVALID_OPTION_FLAG) then
              call Driver_abortFlash(abortMessage)
           else ! Non critical error
              call Logfile_stampMessage(abortMessage)
           endif
        endif
     enddo
  enddo

 
  ! Load data into table data objects
  call eos_LoadTables ( NSPECIES*NEOSTABS, eos_tableHandleVec, eos_errorCode)
  if (eos_errorCode.NE.EOS_OK) then
     call eos_GetErrorMessage ( eos_errorCode, eos_errorMessage )
     write(*,998) 'eos_LoadTables ERROR ', eos_errorCode, ': ', &
                   eos_errorMessage(1:(len_trim(eos_errorMessage)-1))
     do tabno=1, NSPECIES*NEOSTABS
        eos_tableHandleErrorCode = EOS_OK
        call eos_GetErrorCode ( eos_tableHandleVec(tabno), eos_tableHandleErrorCode )
        call eos_GetErrorMessage ( eos_tableHandleErrorCode, eos_errorMessage )
        write(abortMessage,994) 'eos_LoadTables ERROR ', eos_tableHandleErrorCode, ' (TH=', &
                     eos_tableHandleVec(tabno), '): ', &
                     eos_errorMessage(1:(len_trim(eos_errorMessage)-1))
        call Driver_abortFlash(abortMessage)
     enddo
  endif


  ! Set Multispecies properties (A,Z, [dens]) from table metadata and check
  ! that they are consistent with user providied RuntimeParameters
  call Logfile_stampMessage('[Eos eos_initEospac] Warning: material &
          properties will be set from the SESAME database. &
          Corresponding user defined runfile parameters are ignored.')
  do specno=1, NSPECIES

     spec = SPECIES_BEGIN - 1 + specno
     call Simulation_mapIntToStr(spec,speciesName,MAPBLOCK_UNK)

     write(logMessage,'(5X,a,a,a,I5)') 'Material: ', speciesName, &
          ', SESAME table: ', eos_matId(specno) 
     call Logfile_stampMessage(logMessage)
     print *, logMessage
     do infono=1, NEOSINFO
        ! getting info from the 1st table of this species. 
        call eos_GetTableInfo ( eos_tableHandleArr(specno, 1), 1_EOS_INTEGER,&
            eos_infoItems(infono), eos_infoVals(infono), eos_errorCode )
        if (eos_errorCode.NE.EOS_OK) then
           call eos_GetErrorMessage ( eos_errorCode, eos_errorMessage )
           write(abortMessage,998) 'eos_GetTableInfo', &
             eos_errorCode, ': ', &
             eos_errorMessage(1:(len_trim(eos_errorMessage)-1))
           call Driver_abortFlash(abortMessage)
        endif


        if (infono==EOSPAC_INFO_A) then
           call Multispecies_getProperty(spec, A, propertyVal)
           write(logMessage, "(10X,a,1F8.5)"),  'A : ', eos_infoVals(infono)
           call Logfile_stampMessage(logMessage)

           call RuntimeParameters_set('ms_' // speciesName // 'A', &
                                          eos_infoVals(infono))
           call Multispecies_setProperty(spec, A, eos_infoVals(infono))
        else if (infono==EOSPAC_INFO_Z) then
           call Multispecies_getProperty(spec, Z, propertyVal)
           write(logMessage, "(10X,a,1F8.5)"),  'Z : ', eos_infoVals(infono)
           call Logfile_stampMessage(logMessage)

           call RuntimeParameters_set('ms_' // speciesName // 'Z', &
                                          eos_infoVals(infono))
           call Multispecies_setProperty(spec, Z, eos_infoVals(infono))
        !else if (infono>3) then
        !  print *, infono, ':', eos_infoVals(infono)

        !else if (infono==EOSPAC_INFO_DENS) then
        !   call RuntimeParameters_get('sim_rho' // speciesName, propertyVal
        !   if (propertyVal == UNDEFINED_INT) then
        !   write(logMessage, "(10X,a,1F8.5)"),  'Z : ', eos_infoVals(infono)
        !   call Logfile_stampMessage(logMessage)

        !   call RuntimeParameters_set('ms_' // speciesName // 'Z', &
        !                                  eos_infoVals(infono))
        !   call Multispecies_setProperty(spec, Z, eos_infoVals(infono))
        endif
     enddo
  enddo
!
!  !
!  !     interpolate -- errors codes are intentionally produced
!  !
!  X(1) = 300._EOS_REAL
!  X(2) = 60._EOS_REAL
!  X(3) = 800._EOS_REAL
!  X(4) = 80._EOS_REAL
!
!  Y(1) = 20000.0_EOS_REAL
!  Y(2) = 60000.0_EOS_REAL
!  Y(3) = 40000.0_EOS_REAL
!  Y(4) = 2000.0_EOS_REAL
!
!  call eospacIntArr2Vect(tableHandleVec, tableHandleArr)
!
!  do i=1, NEOSTABS
!     write(*,*) ' '
!     call eos_Mix (NSPECIES,  tableHandleArr(:,i), nXYPairs, concInMix, X, Y, F, dFx, dFy, errorCode)
!     if (errorCode.NE.EOS_OK) then
!        call eos_GetErrorMessage ( errorCode, errorMessage )
!        write(*,994) 'eos_Interpolate ERROR ', errorCode, ' (TH=', &
!                     tableHandleVec(i), '): ', &
!                     errorMessage(1:(len_trim(errorMessage)-1))
!     else
!        do j=1, nXYPairs
!            write(*,999) j-1,X(j),Y(j),F(j),dFx(j),dFy(j),errorCode
!        enddo
!     endif
!  enddo
!
!
!  !
!  !     Destroy all data objects
!  !
  !call eos_DestroyAll (eos_errorCode)
  !if (eos_errorCode.NE.EOS_OK) then
  !   do tabno=1, NSPECIES*NEOSTABS
  !      eos_tableHandleErrorCode = EOS_OK
  !      call eos_GetErrorCode ( eos_tableHandleVec(tabno), eos_tableHandleErrorCode )
  !      call eos_GetErrorMessage ( eos_tableHandleErrorCode, eos_errorMessage )
  !      write(*,998) 'eos_DestroyAll ERROR ', eos_tableHandleErrorCode, ': ', &
  !                   eos_errorMessage(1:(len_trim(eos_errorMessage)-1))
  !   enddo
  !endif

994 format (a,i4,a,i1,2a)
995 format (i2,a,f13.6)
997 format (a,:,a,:,i2)
998 format (a,i4,2a)
999 format ('    i=',i2,'    X =',1pe13.6,', Y =',1pe13.6, &
       ', F =',1pe13.6,', dFx =',1pe13.6,', dFy =', &
       1pe13.6,', errorCode: ',i4)
  return

end subroutine eos_initEospac
