!!****ih* source/physics/Eos/localAPI/eos_eospacInterface
!!
!! NAME
!!     eos_mtInterface
!!
!! SYNOPSIS
!!     use eos_mtInterface
!!
!! DESCRIPTION
!!
!! This is an interface module for internal use of
!! multiTemp Eos implementations.
!!
!!***

module eos_eospacInterface
  implicit none
#include "Eos.h"

  interface
     subroutine eos_eospac(mode,vecLen,eosData,vecBegin,vecEnd,massFrac,mask,componentMask)
       implicit none
       integer, INTENT(in) :: mode, vecLen
       real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
       integer, INTENT(in) :: vecBegin,vecEnd
       real, optional, INTENT(in),dimension(:)    :: massFrac
       logical,optional,target, dimension(EOS_VARS+1:EOS_NUM),INTENT(in)::mask
       integer,optional, dimension(N_EOS_TEMP),INTENT(in)::componentMask
     end subroutine eos_eospac
  end interface

end module eos_eospacInterface
