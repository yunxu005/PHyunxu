#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module m_MEF90_DefMechAT1REU
#include "finclude/petscdef.h"
   Use m_MEF90
   Use m_MEF90_DefMechAT_class
   implicit none

!!!
!!!  MEF90_DefMechAT1REU_Type: Ambrosio-Tortorelli-type model whose stiffness
!!!  interpolation follows the Reuss (equal-stress) estimate for a two-phase medium
!!!
!!!  a(alpha) = (1-alpha) e_reu(alpha,p0), with
!!!     e_reu(c,p0) = 1/(c p0 + 1 - c)
!!!  where c is the volume fraction of the particle phase (c = alpha here),
!!!  p0 = E1/E0, E0 is the Young modulus of the matrix and E1 that of the particle.
!!!
!!!  The factor (1-alpha) makes a(1) = 0, i.e. the stiffness vanishes completely.
!!!
!!!  w(alpha), cw, aOrder and wOrder are taken from the standard AT1 model.
!!!

   Type, extends(MEF90_DefMechAT_Type)                 :: MEF90_DefMechAT1REU_Type
      PetscReal                                        :: E0
      PetscReal                                        :: E1
      PetscReal                                        :: p0
   Contains
      Procedure, pass(self)                            :: a   => aAT1REU
      Procedure, pass(self)                            :: Da  => DaAT1REU
      Procedure, pass(self)                            :: D2a => D2aAT1REU

      Procedure, pass(self)                            :: w   => wAT1REU
      Procedure, pass(self)                            :: Dw  => DwAT1REU
      Procedure, pass(self)                            :: D2w => D2wAT1REU
   end Type

   interface MEF90_DefMechAT1REU_Type
      module procedure MEF90_DefMechAT1REU_Constructor
   end interface

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_DefMechAT1REU_Constructor"
!!!
!!!  MEF90_DefMechAT1REU_Constructor: the default constructor for a MEF90_DefMechAT1REU_Type
!!!
   Type(MEF90_DefMechAT1REU_Type) Function MEF90_DefMechAT1REU_Constructor(E0,E1)
      PetscReal,Intent(IN)                           :: E0,E1
      MEF90_DefMechAT1REU_Constructor%E0              = E0
      MEF90_DefMechAT1REU_Constructor%E1              = E1
      MEF90_DefMechAT1REU_Constructor%p0              = MEF90_DefMechAT1REU_Constructor%E1 / MEF90_DefMechAT1REU_Constructor%E0
      MEF90_DefMechAT1REU_Constructor%cw              = 2.0_Kr / 3.0_Kr
      MEF90_DefMechAT1REU_Constructor%aorder          = 2
      MEF90_DefMechAT1REU_Constructor%worder          = 1
      MEF90_DefMechAT1REU_Constructor%type            = 'MEF90_DefMechAT1REU'
   End Function MEF90_DefMechAT1REU_Constructor

#undef __FUNCT__
#define __FUNCT__ "aAT1REU"
!!!
!!!  aAT1REU: the stiffness interpolation, a(alpha) = (1-alpha)/(alpha p0 + 1 - alpha)
!!!
   PetscReal function aAT1REU(self,alpha)
      Class(MEF90_DefMechAT1REU_Type),Intent(IN)        :: self
      PetscReal                                        :: alpha

      PetscReal                                        :: den
      PetscLogDouble                                   :: flops
      PetscErrorCode                                   :: ierr

      den = alpha*self%p0 + 1.0_Kr - alpha
      aAT1REU = (1.0_Kr - alpha)/den
      flops = 5.0
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End function aAT1REU

#undef __FUNCT__
#define __FUNCT__ "DaAT1REU"
!!!
!!!  DaAT1REU: derivative of aAT1REU with respect to alpha
!!!           with den = 1 + delta alpha and delta = p0 - 1,
!!!           Da = -p0 / den^2
!!!
   PetscReal function DaAT1REU(self,alpha)
      Class(MEF90_DefMechAT1REU_Type),Intent(IN)        :: self
      PetscReal                                        :: alpha

      PetscReal                                        :: den
      PetscLogDouble                                   :: flops
      PetscErrorCode                                   :: ierr

      den = 1.0_Kr + (self%p0 - 1.0_Kr)*alpha
      DaAT1REU = -self%p0/den**2
      flops = 4.0
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End function DaAT1REU

#undef __FUNCT__
#define __FUNCT__ "D2aAT1REU"
!!!
!!!  D2aAT1REU: second derivative of aAT1REU with respect to alpha
!!!            D2a = 2 p0 delta / den^3
!!!
   PetscReal function D2aAT1REU(self,alpha)
      Class(MEF90_DefMechAT1REU_Type),Intent(IN)        :: self
      PetscReal                                        :: alpha

      PetscReal                                        :: delta,den

      delta = self%p0 - 1.0_Kr
      den   = 1.0_Kr + delta*alpha
      D2aAT1REU = 2.0_Kr*self%p0*delta/den**3
   End function D2aAT1REU

#undef __FUNCT__
#define __FUNCT__ "wAT1REU"
!!!
!!!  wAT1REU: the "w" function, identical to the standard AT1 model, w(alpha) = alpha
!!!
   PetscReal function wAT1REU(self,alpha)
      Class(MEF90_DefMechAT1REU_Type),Intent(IN)        :: self
      PetscReal                                        :: alpha

      wAT1REU = alpha
   End function wAT1REU

#undef __FUNCT__
#define __FUNCT__ "DwAT1REU"
!!!
!!!  DwAT1REU: derivative of wAT1REU, Dw(alpha) = 1
!!!
   PetscReal function DwAT1REU(self,alpha)
      Class(MEF90_DefMechAT1REU_Type),Intent(IN)        :: self
      PetscReal                                        :: alpha

      DwAT1REU = 1.0_Kr
   End function DwAT1REU

#undef __FUNCT__
#define __FUNCT__ "D2wAT1REU"
!!!
!!!  D2wAT1REU: second derivative of wAT1REU, D2w(alpha) = 0
!!!
   PetscReal function D2wAT1REU(self,alpha)
      Class(MEF90_DefMechAT1REU_Type),Intent(IN)        :: self
      PetscReal                                        :: alpha

      D2wAT1REU = 0.0_Kr
   End function D2wAT1REU
End module m_MEF90_DefMechAT1REU
