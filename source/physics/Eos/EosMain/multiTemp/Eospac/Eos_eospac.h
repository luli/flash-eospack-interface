#if 0
This files contains the indices associated with EOSPAC that are used by FLASH.

The values of the constants are distincts from the ones defined in EOSPAC.
The relationship with EOSPAC constants is given bellow:
 FLASH EOSPAC Interface     EOSPAC library
 ----------------------     --------------
    EOSPAC_INFO_Z    : EOS_Mean_Atomic_Num
    EOSPAC_INFO_A    : EOS_Mean_Atomic_Mass
    EOSPAC_INFO_DENS : EOS_Normal_Density
    EOSPAC_X_CFACTOR : EOS_X_Convert_Factor
    EOSPAC_Y_CFACTOR : EOS_Y_Convert_Factor
    EOSPAC_F_CFACTOR : EOS_F_Convert_Factor

    EOSPAC_OPT_SMOOTH: EOS_SMOOTH
    


#endif

#define NEOSTABS 12
#define NEOSINFO 6
#define NEOSOPT 1


#define EOSPAC_INFO_Z 1
#define EOSPAC_INFO_A 2
#define EOSPAC_INFO_DENS 3
#define EOSPAC_INFO_XFACTOR 4
#define EOSPAC_INFO_YFACTOR 5
#define EOSPAC_INFO_FFACTOR 6
