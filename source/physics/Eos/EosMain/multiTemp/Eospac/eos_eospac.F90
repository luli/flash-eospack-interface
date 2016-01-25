!!****if* source/physics/Eos/EosMain/multiTemp/Eospac/eos_eospac
!!
!! NAME
!!
!!  eos_eospac
!!
!! SYNOPSIS
!!
!!  call      eos_eospac(integer(IN) :: mode,
!!                integer(IN)                 :: vecLen,
!!                real(INOUT)                 :: eosData(vecLen*EOS_NUM),
!!                integer(IN)                 :: vecBegin,
!!                integer(IN)                 :: vecEnd,
!!             optional,real(IN),dimension(*) :: massFrac(vecLen*NSPECIES),
!!             optional,target, logical(IN)   :: mask(EOS_VARS+1:EOS_NUM),
!!             optional, integer(IN)          :: componentMask(N_EOS_TEMP) )
!!
!! DESCRIPTION
!!
!!  Call Eos implementation routines of (possibly) different types and combine
!!  the results.
!!
!!  This subroutine is called from the Eospac implementation of Eos.F90,
!!  potentially from a Newton-Raphson loop for finding the right temperatures
!!  to satisfy energy requirements.
!!
!!  The model underlying THIS implementation (KW July 2011) is that various
!!  materials that co-occur in the same cells don't notice each other.
!!  Each material is handled as if it filled the cell (fully).
!!  Additive quantities are then added up to form the results returned.
!!
!!  This implementation is typically called for ions and electrons separately,
!!  and not necessarily at the same temperature.  While contributions from
!!  different species in the FLASH sence, i.e., different materials, are
!!  combined here, the caller is responsible for combining contributions
!!  of different components, i.e., ions and electrons (DEV: and radiation?).
!!
!! ARGUMENTS 
!! 
!!  mode :    Selects the mode of operation of the Eos unit.
!!            This routine will normally be called with a MODE_DENS_TEMP* mode.
!!            For current (FLASH4-beta) Multitype Eos.F90 purposes, the modes
!!            used will be MODE_DENS_TEMP_ELE and MODE_DENS_TEMP_ION if the mode
!!            for which Eos is called was MODE_DENS_TEMP_GATHER or
!!            MODE_DENS_EI_GATHER.
!!
!!  vecLen   : number of points (cells) for which the eosData array is sized.
!!             This is the maximum allowed value for vecEnd.
!!
!!  eosData  : This array is the data structure through which variable values are 
!!             passed in and out of the Eos routine. The arrays is sized as 
!!             EOS_NUM*vecLen. EOS_NUM, and individual input and output
!!             Eos variables are defined in Eos.h. The array is organizes such that
!!             the first 1:vecLen entries represent the first Eos variable, vecLen+1:
!!             2*vecLen represent the second Eos variable and so on.
!!             This is the same kind of packaging as for Eos.F90 etc.
!!
!!  vecBegin : Index of first cell in eosData to handle.
!!             Can be used to limit operation to a subrange of cells, untested.
!!             Must be greater than or equal to 1.
!!
!!  vecEnd   : Index of last cell in eosData to handle.
!!             Can be used to limit operation to a subrange of cells, untested.
!!             Must be less than or equal to vecLen.
!!
!!  massFrac : Contains the mass fractions of the species included in
!!             the simulation. The array is sized as NSPECIES*vecLen.
!!
!!  mask     : Mask is a logical array the size of EOS_DERIVS (number
!!              of partial derivatives that can be computed, defined in
!!              Eos.h), where each index represents a specific partial derivative
!!              that can be calculated by the Eos unit. A .true. value in mask 
!!              results in the corresponding derivative being calculated and 
!!              returned. It should preferably be dimensioned as
!!              mask(EOS_VARS+1:EOS_NUM) in the calling routine 
!!              to exactly match the arguments declaration in Eos Unit.
!!             Note that the indexing of mask does not begin at 1, but rather at one past
!!             the number of variables.
!!
!!             An implementation that does not need derivative quantities should
!!             set the mask equal to .false.
!!
!!  componentMask : A 3-element integer vector that contains 1 for componenets to
!!                  consider and 0 for components to ignore in this call, in the
!!                  order (ions, electrons, radiation).
!!                  This kind of duplicates information already present in the
!!                  mode argument (MODE_DENS_TEMP_ION vs. MODE_DENS_TEMP_ELE),
!!                  and may be useless and superfluous as well as redundant.
!!
!! NOTES
!!
!!  The global variable eos_type must be set to EOS_MTYPE in order for this routine to
!!  operate in Multitype mode; otherwise, it will just (pretty uselessly) pass the call
!!  on to the appropriate simpler Eos implementation.
!!
!!  NSPECIES is defined in Flash.h.
!!
!!  EOS_VARS and EOS_NUM  are defined in Eos.h.
!!  
!!  MODE_DENS_TEMP, MODE_DENS_EI, MODE_DENS_PRES, etc. are defined in constants.h.
!!
!!
!! SEE ALSO
!!
!!  Eos.h    defines the variables used.
!!  Eos_wrapped  sets up the data structure.
!!
!!***

subroutine eos_eospac(mode,vecLen,eosData,vecBegin,vecEnd,massFrac,mask,componentMask)

!==============================================================================
  use Multispecies_interface, ONLY : Multispecies_getProperty, &
       Multispecies_getPropertyVector, &
       Multispecies_getSumInv
  use Driver_interface, ONLY : Driver_abortFlash
  use Eos_data, ONLY : eos_meshMe, eos_type, eos_entrEleScaleChoice
  use eos_helmConstData, ONLY: eos_ao3, eos_avo, eos_kerg, eos_kergavo
!!$  use eos_vecData, ONLY:  detRow, dptRow !, dpdRow, dedRow, pelRow, neRow, etaRow, cvRow, cpRow
  use eos_eospacData, ONLY : eos_tableTypeArr, eos_tableTypeVec,&
       eos_tableTypeBase, eos_matId, eos_massFraction, eos_matIdArr,&
       eos_matIdVec, eos_tableHandleArr, eos_tableHandleVec, &
       eos_errorCode, eos_tableHandleErrorCode, eos_infoVals,&
       eos_infoItems, eos_errorMessage, eos_optionFlags, &
       eospacIntVect2Arr, eospacIntArr2Vect, eos_indexUnwrap, &
       eos_indexWrap, eos_idxVar

  use libeospac_interface

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "Multispecies.h"
  integer, INTENT(in) :: mode, vecLen
  real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
  integer, INTENT(in) :: vecBegin,vecEnd
  real, optional, INTENT(in),dimension(:)    :: massFrac
  logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask
  integer, optional, INTENT(in),dimension(N_EOS_TEMP)::componentMask

  integer :: spec, specno
  integer :: ilo,ihi, rowLen
  integer :: passMode
  integer,dimension(N_EOS_TEMP) :: passCMask
  real, dimension(EOS_NUM*vecLen) :: curData, outData
  real, dimension(NSPECIES) :: Avect, Yvect, Zmin_vect
  real, dimension(NSPECIES,vecLen) :: Xmat, Ymat, Fmat !!!, Pmat
  real, dimension(NSPECIES*vecLen) :: Fvec ! reshaped version of Fmat
  real, dimension(vecLen) :: dFx(vecBegin:vecEnd), dFy(vecBegin:vecEnd)
  real :: curA, curZ, curZMin
  real :: curGamma
  real :: abarInv
  integer :: specieStart, specieEnd
  integer :: dens, temp, pres, eint, abar, zbar, ekin
  integer :: gamc
  integer :: tempIon,tempEle,tempRad                                                
  integer :: eintIon,eintEle,eintRad, entrEle, det, dst,cvion,cvele
  integer :: presIon,presEle,presRad, presComp, dpt, dpd
  integer :: i, v, offs
  integer :: subtype
  logical :: arrayBoundHackMode
  logical :: handleRadHere
  integer :: tempToUse
  integer :: tempToPrint
  integer :: speciesEosType, numberExistingSpecies

  logical,save,target, dimension(EOS_VARS+1:EOS_NUM) :: maskInternal
  data maskInternal / EOS_DERIVS * .FALSE. /
  logical,pointer, dimension(:) :: maskPtr

  logical :: presMFrac, presMask
  logical,parameter :: doRad = .TRUE.
  character(len=300)  :: abortMessage
  character(len=15) :: table_strId
  integer(EOS_INTEGER):: xyBounds(vecLen)

  presMask = present(mask)
  presMFrac = present(massFrac)
  
  maskPtr => maskInternal
  if (presMask) then
     maskPtr => mask
  end if
 
  if (lbound(maskPtr,1) == 1) then
     !gfortran does not set lower bound of pointer to EOS_VARS+1.
     !(I think gfortran is doing the correct thing).
     arrayBoundHackMode = .true.

     nullify(maskPtr)
     allocate(maskPtr(EOS_VARS+1:EOS_NUM))
     if (present(mask)) then
        maskPtr(EOS_VARS+1:EOS_NUM) = mask(EOS_VARS+1:EOS_NUM)
     else
        maskPtr(EOS_VARS+1:EOS_NUM) = maskInternal(EOS_VARS+1:EOS_NUM)
     end if
  else
     arrayBoundHackMode = .false.
  end if

  if (present(componentMask)) then
     passCMask = componentMask
  else
     passCMask = (/1,1,1/)
  end if
     

!!$  if (present(vecBegin)) then
     ilo = vecBegin
!!$  else
!!$     ilo = 1
!!$  end if
!!$  if (present(vecEnd)) then
     ihi = vecEnd
!!$  else
!!$     ihi = vecLen
!!$  end if
!!$  rowLen = ihi - ilo + 1
!!$#ifdef DEBUG_EOS
!!$  if (ilo < 1 .OR. ilo > vecLen) then
!!$     print*,'[eos_multiTypeByTemp] ilo is',ilo
!!$     call Driver_abortFlash("[eos_multiTypeByTemp] invalid ilo")
!!$  end if
!!$  if (ihi < 1 .OR. ihi > vecLen) then
!!$     print*,'[eos_multiTypeByTemp] ihi is',ihi
!!$     call Driver_abortFlash("[eos_multiTypeByTemp] invalid ihi")
!!$  end if
!!$  if (rowLen < 0 .OR. rowLen > vecLen) then
!!$     print*,'[eos_multiTypeByTemp] rowLen is',rowLen
!!$     call Driver_abortFlash("[eos_multiTypeByTemp] invalid rowLen")
!!$  end if
!!$#endif
!!$  if (rowLen == 0) then
!!$     print*,'[eos_multiTypeByTemp] rowLen is 0.'
!!$  end if

  pres = (EOS_PRES-1)*vecLen
  dens = (EOS_DENS-1)*vecLen
  temp = (EOS_TEMP-1)*vecLen
  eint = (EOS_EINT-1)*vecLen
  gamc = (EOS_GAMC-1)*vecLen
  abar = (EOS_ABAR-1)*vecLen
  zbar = (EOS_ZBAR-1)*vecLen
  ekin = (EOS_EKIN-1)*vecLen
  tempIon = (EOS_TEMPION-1)*vecLen
  tempEle = (EOS_TEMPELE-1)*vecLen
  tempRad = (EOS_TEMPRAD-1)*vecLen
  eintIon = (EOS_EINTION-1)*vecLen
  eintEle = (EOS_EINTELE-1)*vecLen
  eintRad = (EOS_EINTRAD-1)*vecLen
  presIon = (EOS_PRESION-1)*vecLen
  presEle = (EOS_PRESELE-1)*vecLen
  presRad = (EOS_PRESRAD-1)*vecLen
  entrEle = (EOS_ENTRELE-1)*veclen
  det = (EOS_DET-1)*vecLen
  dpt = (EOS_DPT-1)*vecLen
  dpd = (EOS_DPD-1)*vecLen
  dst = (EOS_DST-1)*vecLen
  cvion = (EOS_CVION-1)*vecLen
  cvele = (EOS_CVELE-1)*vecLen

  rowLen = ihi - ilo + 1
#ifdef DEBUG_EOS
  if (ilo < 1 .OR. ilo > vecLen) then
     print*,'[eos_idealGamma3T] ilo is',ilo
     call Driver_abortFlash("[eos_idealGamma3T] invalid ilo")
  end if
  if (ihi < 1 .OR. ihi > vecLen) then
     print*,'[eos_idealGamma3T] ihi is',ihi
     call Driver_abortFlash("[eos_idealGamma3T] invalid ihi")
  end if
  if (rowLen < 0 .OR. rowLen > vecLen) then
     print*,'[eos_idealGamma3T] rowLen is',rowLen
     call Driver_abortFlash("[eos_idealGamma3T] invalid rowLen")
  end if
#endif
  if (rowLen == 0) then
     print*,'[eos_idealGamma3T] rowLen is 0.'
  end if




  if (presMask) then
     ! nothing special
  end if

  select case (eos_type)
     case(EOS_EOSPAC)
        ! ok...
     case default
        if (eos_meshMe==MASTER_PE) print*,"eos_eospac: only supports EOS type EOS_EOSPAC, called for eosType=", eos_type
        call Driver_abortFlash("eos_eospac: only supports EOS_EOSPAC, called for wrong EOS type!")
  end select

  if (.NOT. presMFrac) then
     call Driver_abortFlash('eos_eospac needs mass fractions!')
  end if
  outData(:) = 0.0
  curData(:) = eosData(:)
  
  do i = ilo,ihi
     specieStart = (i-1)*NSPECIES + 1
     specieEnd = i*NSPECIES
     call Multispecies_getSumInv(A, abarInv ,massFrac(specieStart:specieEnd))
     eosData(abar+i) = 1.0/abarInv
  end do
 
  call Multispecies_getPropertyVector(A,Avect)
  do specno = 1,NSPECIES
     Xmat(specno,ilo:ihi) = massFrac((ilo-1)*NSPECIES+specno:(ihi-1)*NSPECIES+specno:NSPECIES)
     Ymat(specno,ilo:ihi) = Xmat(specno,ilo:ihi) / Avect(specno)
     ! this is apparently the number fraction
     Fmat(specno,ilo:ihi) = Ymat(specno,ilo:ihi) * eosData(abar+ilo:abar+ihi)
     ! number density arranged according to Eospac method 
     Fvec((specno-1)*vecLen+1:specno*vecLen) = Fmat(specno, ilo:ihi)
     !print *, Fvec(specno,ilo:ihi)
 !!!        Pmat(specno,ilo:ihi) = 0.0
     call Multispecies_getProperty(specno,MS_ZMIN,Zmin_vect(specno)) 
  end do

  do v=1,EOS_NUM
    offs=(v-1)*vecLen
    outData(offs+ilo:offs+ihi) = eosData(offs+ilo:offs+ihi)
  end do

  !print *, 'mode', mode, mask(EOS_CV), mask(EOS_DPD), mask(EOS_DPT),&
  !             mask(EOS_DET), mask(EOS_CVION), mask(EOS_CVELE), mask(EOS_GAMC)
  
  994 format (a,i5,4a)
  select case (mode)
  ! density, temperature taken as input
  case (MODE_DENS_TEMP_ION)
!     print *, eosData(dens+ilo:dens+ihi)
     table_strId = 'EOS_Pic_DT'
     call eos_Mix(NSPECIES,  eos_tableHandleArr(:,eos_idxVar(EOS_Pic_DT)), &
                  vecLen, Fvec, &
                  eosData(dens+ilo:dens+ihi), &
                  eosData(tempIon+ilo:tempIon+ihi), &
                  outData(presIon+ilo:presIon+ihi), &
                    dFx, dFy, eos_errorCode)
     if(mask(EOS_DPD)) then
       outData(dpd+ilo:dpd+ihi) = dFx
     end if

     if (eos_errorCode.NE.EOS_OK) goto  70

     table_strId = 'EOS_Uic_DT'
     call eos_Mix(NSPECIES,  eos_tableHandleArr(:,eos_idxVar(EOS_Uic_DT)), &
                  vecLen, Fvec, &
                  eosData(dens+ilo:dens+ihi), &
                  eosData(tempIon+ilo:tempIon+ihi), &
                  outData(eintIon+ilo:eintIon+ihi), &
                    dFx, dFy, eos_errorCode)
      
     if (eos_errorCode.NE.EOS_OK) goto 70

  case (MODE_DENS_TEMP_ELE)
     table_strId = 'EOS_Pe_DT'
     call eos_Mix(NSPECIES,  eos_tableHandleArr(:,eos_idxVar(EOS_Pe_DT)), &
                  vecLen, Fvec, &
                  eosData(dens+ilo:dens+ihi), &
                  eosData(tempEle+ilo:tempEle+ihi), &
                  outData(presEle+ilo:presEle+ihi), &
                    dFx, dFy, eos_errorCode)
     if(mask(EOS_DPD)) then
       outData(dpd+ilo:dpd+ihi) = dFx
     end if

     if (eos_errorCode.NE.EOS_OK) goto 70

     table_strId = 'EOS_Ue_DT'
     call eos_Mix(NSPECIES,  eos_tableHandleArr(:,eos_idxVar(EOS_Ue_DT)), &
                  vecLen, Fvec, &
                  eosData(dens+ilo:dens+ihi), &
                  eosData(tempEle+ilo:tempEle+ihi), &
                  outData(eintEle+ilo:eintEle+ihi), &
                    dFx, dFy, eos_errorCode)
      
     if (eos_errorCode.NE.EOS_OK) goto 70

  case (MODE_DENS_EI_ION)
     table_strId = 'EOS_Pic_DUic'
     !call eos_CheckExtrap(eos_tableHandleArr(1,eos_idxVar(EOS_Pic_DUic)), &
     !             vecLen, & 
     !             eosData(dens+ilo:dens+ihi), &
     !             1.0e-10*eosData(eintIon+ilo:eintIon+ihi), &
     !             xyBounds, &
     !              eos_errorCode)
     !print *, xyBounds
     !if (eos_errorCode.NE.EOS_OK) goto 70
     !call eos_Interpolate(eos_tableHandleArr(1,eos_idxVar(EOS_Pic_DUic)), &
     !             vecLen, & 
     !             eosData(dens+ilo:dens+ihi), &
     !             eosData(eintIon+ilo:eintIon+ihi), &
     !             outData(presIon+ilo:presIon+ihi), &
     !               dFx, dFy, eos_errorCode)
     !if (eos_errorCode.NE.EOS_OK) goto 70

     call eos_Mix(NSPECIES,  eos_tableHandleArr(:,eos_idxVar(EOS_Pic_DUic)), &
                  vecLen, Fvec, &
                  eosData(dens+ilo:dens+ihi), &
                  1.0e-10*eosData(eintIon+ilo:eintIon+ihi), &
                  outData(presIon+ilo:presIon+ihi), &
                    dFx, dFy, eos_errorCode)
     !if(mask(EOS_DPD)) then
     !  outData(dpd+ilo:dpd+ihi) = dFx
     !end if

     if (eos_errorCode.NE.EOS_OK) goto 70

     table_strId = 'EOS_T_DUic'
     call eos_Mix(NSPECIES,  eos_tableHandleArr(:,eos_idxVar(EOS_T_DUic)), &
                  vecLen, Fvec, &
                  eosData(dens+ilo:dens+ihi), &
                  1.0e-10*eosData(eintIon+ilo:eintIon+ihi), &
                  outData(tempIon+ilo:tempIon+ihi), &
                    dFx, dFy, eos_errorCode)
      
     if (eos_errorCode.NE.EOS_OK) goto 70

  case (MODE_DENS_EI_ELE)
     table_strId = 'EOS_Pe_DUe'
     call eos_Mix(NSPECIES,  eos_tableHandleArr(:,eos_idxVar(EOS_Pe_DUe)), &
                  vecLen, Fvec, &
                  eosData(dens+ilo:dens+ihi), &
                  1.0e-10*eosData(eintEle+ilo:eintEle+ihi), &
                  outData(presEle+ilo:presEle+ihi), &
                    dFx, dFy, eos_errorCode)

     if (eos_errorCode.NE.EOS_OK) goto 70

     table_strId = 'EOS_T_DUe'
     call eos_Mix(NSPECIES,  eos_tableHandleArr(:,eos_idxVar(EOS_T_DUe)), &
                  vecLen, Fvec, &
                  eosData(dens+ilo:dens+ihi), &
                  1.0e-10*eosData(eintEle+ilo:eintEle+ihi), &
                  outData(tempEle+ilo:tempEle+ihi), &
                    dFx, dFy, eos_errorCode)
      
     if (eos_errorCode.NE.EOS_OK) goto 70




  case default
        if (eos_meshMe==MASTER_PE) print*,"eos_eospac: eos mode not supported", mode
        call Driver_abortFlash("eos_eospac: eos mode not supported!")

70   continue
     if (eos_errorCode.NE.EOS_OK) then
        call eos_GetErrorMessage ( eos_errorCode, eos_errorMessage )
        write(abortMessage,994) 'eos_eospac/eos_Mix: error', eos_errorCode,&
          ' in table ', table_strId,': ', eos_errorMessage(1:(len_trim(eos_errorMessage)-1))
        call Driver_abortFlash(abortMessage)
     endif
  end select
 
 !!$     print*,'Ymat:',Ymat
 
 !numberExistingSpecies = 0
 ! do spec=SPECIES_BEGIN,SPECIES_END
 !    specno = spec - SPECIES_BEGIN + 1
 !    if (ALL(massFrac((ilo-1)*NSPECIES+specno:(ihi-1)*NSPECIES+specno:NSPECIES).LE. 1.0e-20)) then
 !#ifdef DEBUG_EOS
 !       print*,' skipping species #',specno,' for',&
 !            size(massFrac((ilo-1)*NSPECIES+specno:(ihi-1)*NSPECIES+specno:NSPECIES)),&
 !            ' cells since all massfracs are (nearly) 0.'
 !#endif
 !    else
 !!!$           print*,'NOT Skipping species #',specno,' for',&
 !!!$                size(massFrac((ilo-1)*NSPECIES+specno:(ihi-1)*NSPECIES+specno:NSPECIES)),&
 !!!$                ' cells:',&
 !!!$                massFrac((ilo-1)*NSPECIES+specno:(ihi-1)*NSPECIES+specno:NSPECIES)
 !
 !       numberExistingSpecies = numberExistingSpecies + 1
 !       call Multispecies_getProperty(spec,MS_EOSTYPE,speciesEosType)
 !       do v=1,EOS_NUM
 !          offs=(v-1)*vecLen
 !          curData(offs+ilo:offs+ihi) = eosData(offs+ilo:offs+ihi)
 !       end do
 !       curData(dens+ilo:dens+ihi) = eosData(dens+ilo:dens+ihi) * Xmat(specno,ilo:ihi)
 !
 !       call Multispecies_getProperty(spec,GAMMA,curGamma) !to delete!?
 !       call Multispecies_getProperty(spec,A    ,curA)
 !       call Multispecies_getProperty(spec,Z    ,curZ)
 !       call Multispecies_getProperty(spec,MS_ZMIN,curZMin)
 !       curZ = max(curZ,curZMin)
 !       curData(gamc+ilo:gamc+ihi) = curGamma !to delete!?
 !       curData(abar+ilo:abar+ihi) = curA
 !       curData(zbar+ilo:zbar+ihi) = curZ

 !       call eos_idealGamma3T(mode, vecLen, curData,ilo,ihi, &
 !            eosType=EOS_GAM, material=spec, massFrac=massFrac, mask=mask, componentMask=passCMask)
 !
 !       call eos_rescaleCurData(curData,mode,maskPtr)
 !       call eos_accumData(outData,curData,mode)
 !    end if
 ! end do

  call eos_rescaleOutData(outData,mode)
!
  outData(abar+ilo:abar+ihi) = eosData(abar+ilo:abar+ihi)
  outData(ekin+ilo:ekin+ihi) = eosData(ekin+ilo:ekin+ihi)
  do v=1,EOS_NUM
     offs=(v-1)*vecLen
     if (v>EOS_VARS) then
        if (present(mask)) then
           if(.NOT. mask(v)) cycle
        else
           cycle
        end if
     end if
     eosData(offs+ilo:offs+ihi) = outData(offs+ilo:offs+ihi)
  end do
!!$     detRow(ilo:ihi) = eosData(det+ilo:det+ihi)
!!$     dptRow(ilo:ihi) = eosData(dpt+ilo:dpt+ihi)

#ifdef DEBUG_EOS
9999 format(' E.o.multiTypeByTemp     e:',1P,(7(1x,G20.6)))
     write(*,9999) (eosData(eint+i),eosData(eintIon+i),eosData(eintEle+i), eosData(eintRad+i), eosData(temp+i),&
          eosData(dens+i),eosData(dens+i)/eosData(abar+i)*eos_avo ,i=ilo,ihi)
     select case (mode)
        case(MODE_DENS_TEMP_ION)
           tempToPrint = tempIon
        case(MODE_DENS_TEMP_ELE)
           tempToPrint = tempEle
        case(MODE_DENS_TEMP_RAD)
           tempToPrint = tempRad
        case default
           tempToPrint = temp
        end select
9998 format(' E.o.multiTypeByTemp     p:',1P,5(1x,G20.6))
     write(*,9998) (eosData(pres+i),eosData(presIon+i),eosData(presEle+i), eosData(presRad+i), eosData(tempToPrint+i),&
         i=ilo,ihi)
#if(0)
!!$9997 format(' E.o.multiTypeByTemp     Z:',1P,4(1x,G20.6))
!!$     write(*,9997) (eosData(zbar+i),i=ilo,ihi)
#endif
9996 format(' E.o.multiTypeByTemp     det',1P,4(G20.6,1x))
     write(*,9996) (eosData(det+i),i=ilo,ihi)
#endif

  !print *, eosData(presIon+ilo:presIon+ihi)

  if (arrayBoundHackMode .eqv. .true.) then
     deallocate(maskPtr)
  end if

  return

contains
  subroutine eos_rescaleCurData(curData, mode, mask)
    real,intent(INOUT),  dimension(EOS_NUM*vecLen) :: curData
    integer, intent(in) :: mode
    logical, INTENT(in),dimension(EOS_VARS+1:EOS_NUM) :: mask
    real,               dimension(vecLen) :: presRelevant
    real    :: prefactor(ilo:ihi)
    integer :: ic

!    curData(gamc+ilo:gamc+ihi) = Fmat(specno,ilo:ihi) / (curData(gamc+ilo:gamc+ihi)-1.0)
    if (mode==MODE_DENS_TEMP_ION.OR.mode==MODE_DENS_TEMP_ELE.OR.mode==MODE_DENS_TEMP_RAD.OR.&
         mode==MODE_DENS_TEMP_COMP) then
       presRelevant = 0.0
       do ic=0,N_EOS_TEMP-1
          if (passCMask(ic+1) .ne. 0) then
             presRelevant(ilo:ihi) = presRelevant(ilo:ihi) + &
                  curData(presIon+ic*vecLen+ilo:presIon+ic*vecLen+ihi)
          end if

       end do
    else
       presRelevant(ilo:ihi) = curData(pres+ilo:pres+ihi)
    end if
#ifdef GAMC_AVERAGE_COUPLED
    curData(gamc+ilo:gamc+ihi) = &
         presRelevant(ilo:ihi) / (curData(gamc+ilo:gamc+ihi)-1.0)
#else
    curData(gamc+ilo:gamc+ihi) = &
         presRelevant(ilo:ihi) * curData(gamc+ilo:gamc+ihi)
#endif

!!$    curData(abar+ilo:abar+ihi) = curData(abar+ilo:abar+ihi) * Xmat(specno,ilo:ihi)
    curData(zbar+ilo:zbar+ihi) = curData(zbar+ilo:zbar+ihi) * Ymat(specno,ilo:ihi)

    curData(eint+ilo:eint+ihi) = curData(eint+ilo:eint+ihi) * Xmat(specno,ilo:ihi)
    if (mask(EOS_EINTION)) &
         curData(eintIon+ilo:eintIon+ihi) = curData(eintIon+ilo:eintIon+ihi) * Xmat(specno,ilo:ihi)
    if (mask(EOS_EINTELE)) &
         curData(eintEle+ilo:eintEle+ihi) = curData(eintEle+ilo:eintEle+ihi) * Xmat(specno,ilo:ihi)
    if (mask(EOS_EINTRAD)) &
         curData(eintRad+ilo:eintRad+ihi) = curData(eintRad+ilo:eintRad+ihi) * Xmat(specno,ilo:ihi)
    curData(entrEle+ilo:entrEle+ihi) = curData(entrEle+ilo:entrEle+ihi) * Xmat(specno,ilo:ihi)
!!$    print*,'EntrA:', curData(entrEle+ilo:entrEle+ihi)
    if (mode==MODE_DENS_TEMP_ELE) then
       select case (eos_entrEleScaleChoice)
       case(1)
          prefactor = eos_kerg
       case(2,3)
          prefactor = eosData(zbar+ilo:zbar+ihi) / eosData(abar+ilo:abar+ihi) * eos_kergavo
       case(4,6)
          prefactor = eosData(zbar+ilo:zbar+ihi) / eosData(abar+ilo:abar+ihi) ! Ye
       case(5,7)
          prefactor = 1.0
       case default
          call Driver_abortFlash('eos_eospac - unsupported eos_entrEleScaleChoice value')
       end select
       select case (eos_entrEleScaleChoice)
       case(4,6)
          curData(entrEle+ilo:entrEle+ihi) = curData(entrEle+ilo:entrEle+ihi) + &
            curData(zbar+ilo:zbar+ihi) * log(Xmat(specno,ilo:ihi)) ! includes * Ymat(specno,ilo:ihi)
       case default
          curData(entrEle+ilo:entrEle+ihi) = curData(entrEle+ilo:entrEle+ihi) + &
            prefactor * Fmat(specno,ilo:ihi) * log(Xmat(specno,ilo:ihi))
       end select
!!$       print*,'EntrB:', curData(entrEle+ilo:entrEle+ihi)
    end if

    if (mask(EOS_DET)) &
         curData(det+ilo:det+ihi) = curData(det+ilo:det+ihi) * Xmat(specno,ilo:ihi)
    if (mask(EOS_DST)) &
         curData(dst+ilo:dst+ihi) = curData(dst+ilo:dst+ihi) * Xmat(specno,ilo:ihi)
    if (mask(EOS_CVION)) &
         curData(cvion+ilo:cvion+ihi) = curData(cvion+ilo:cvion+ihi) * Xmat(specno,ilo:ihi)
    if (mask(EOS_CVELE)) &
         curData(cvele+ilo:cvele+ihi) = curData(cvele+ilo:cvele+ihi) * Xmat(specno,ilo:ihi)
  end subroutine eos_rescaleCurData

  subroutine eos_accumData(cumData, curData, mode)
    real,intent(INOUT), dimension(EOS_NUM*vecLen) :: cumData
    real,intent(in),  dimension(EOS_NUM*vecLen) :: curData
    integer, intent(in) :: mode

    integer :: vv,offs

    do vv=1,EOS_NUM
       offs=(vv-1)*vecLen
       if (vv>EOS_VARS) then
          if (present(mask)) then
             if(.NOT. mask(vv)) cycle
          else
             cycle
          end if
       end if
       cumData(offs+ilo:offs+ihi) = cumData(offs+ilo:offs+ihi) + curData(offs+ilo:offs+ihi)
    end do
    cumData(temp+ilo:temp+ihi) = curData(temp+ilo:temp+ihi)
    cumData(tempIon+ilo:tempIon+ihi) = curData(tempIon+ilo:tempIon+ihi)
    cumData(tempEle+ilo:tempEle+ihi) = curData(tempEle+ilo:tempEle+ihi)
    cumData(tempRad+ilo:tempRad+ihi) = curData(tempRad+ilo:tempRad+ihi)
  end subroutine eos_accumData

  subroutine eos_rescaleOutData(outData, mode)
    real,intent(INOUT),  dimension(EOS_NUM*vecLen) :: outData
    integer, intent(in) :: mode
    real,               dimension(vecLen) :: presRelevant
    integer :: ic

    if (mode==MODE_DENS_TEMP_ION.OR.mode==MODE_DENS_TEMP_ELE.OR.mode==MODE_DENS_TEMP_RAD.OR.&
         mode==MODE_DENS_TEMP_COMP) then
       presRelevant = 0.0
       do ic=0,N_EOS_TEMP-1
          if(passCMask(ic+1) .ne. 0) then
             presRelevant(ilo:ihi) = presRelevant(ilo:ihi) + &
                  outData(presIon+ic*vecLen+ilo:presIon+ic*vecLen+ihi)
          end if
       end do
    else
       presRelevant(ilo:ihi) = outData(pres+ilo:pres+ihi)
    end if

!    outData(gamc+ilo:gamc+ihi) = 1.0 +  1.0/outData(gamc+ilo:gamc+ihi)
#ifdef GAMC_AVERAGE_COUPLED
    outData(gamc+ilo:gamc+ihi) = 1.0 +  presRelevant(ilo:ihi)/outData(gamc+ilo:gamc+ihi)
#else
    outData(gamc+ilo:gamc+ihi) = outData(gamc+ilo:gamc+ihi) / presRelevant(ilo:ihi)
#endif

    outData(zbar+ilo:zbar+ihi) = outData(zbar+ilo:zbar+ihi) * eosData(abar+ilo:abar+ihi)

  end subroutine eos_rescaleOutData

end subroutine eos_eospac
