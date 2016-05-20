export g03root=/opt/
export GAUSS_SCRDIR=/home/$USER
if [[ -d "$g03root" ]];
    then gr=$g03root
fi
GAUSS_EXEDIR="$gr/g03/bsd:$gr/g03/private:$gr/g03"
GAUSS_ARCHDIR="$gr/g03/arch"
GMAIN=$GAUSS_EXEDIR
PATH=$PATH:$GAUSS_EXEDIR
LD_LIBRARY_PATH=$GAUSS_EXEDIR
G03BASIS="$gr/g03/basis"
F_ERROPT1="271,271,2,1,2,2,2,2"
#following for sgi debugging
#TRAP_FPE="DEBUG;OVERFL=ABORT;DIVZERO=ABORT;INVALID=ABORT;INT_OVERFL=ABORT"
TRAP_FPE="OVERFL=ABORT;DIVZERO=ABORT;INT_OVERFL=ABORT"
MP_STACK_OVERFLOW="OFF"
# to partially avoid KAI stupidity
KMP_DUPLICATE_LIB_OK="TRUE"
export GAUSS_EXEDIR GAUSS_ARCHDIR PATH GMAIN LD_LIBRARY_PATH F_ERROPT1 TRAP_FPE MP_STACK_OVERFLOW \
  KMP_DUPLICATE_LIB_OK G03BASIS
