#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module m_MEF90_DefMechAT1SC
#include "finclude/petscdef.h"
   Use m_MEF90
   Use m_MEF90_DefMechAT_class
   implicit none

!!!
!!!  MEF90_DefMechAT1SC_Type: Ambrosio-Tortorelli-type model whose stiffness
!!!  interpolation follows the self-consistent (SC) estimate for a two-phase medium
!!!
!!!  a(alpha) = (1-alpha) e_sc(alpha,p0), with
!!!     A(c,p0)  = 5(1-c)(1-p0) + 2 p0 - 3
!!!     e_sc(c,p0)= (-A + sqrt(A^2 + 24 p0)) / (6 p0) = 4/(sqrt(A^2+24 p0) + A)
!!!  where c is the volume fraction of the particle phase (c = alpha here),
!!!  p0 = E1/E0, E0 is the Young modulus of the matrix and E1 that of the particle.
!!!
!!!  The factor (1-alpha) makes a(1) = 0, i.e. the stiffness vanishes completely.
!!!
!!!  w(alpha), cw, aOrder and wOrder are taken from the standard AT1 model.
!!!

   Type, extends(MEF90_DefMechAT_Type)                 :: MEF90_DefMechAT1SC_Type
      PetscReal                                        :: E0
      PetscReal                                        :: E1
      PetscReal                                        :: p0
   Contains
      Procedure, pass(self)                            :: a   => aAT1SC
      Procedure, pass(self)                            :: Da  => DaAT1SC
      Procedure, pass(self)                            :: D2a => D2aAT1SC

      Procedure, pass(self)                            :: w   => wAT1SC
      Procedure, pass(self)                            :: Dw  => DwAT1SC
      Procedure, pass(self)                            :: D2w => D2wAT1SC
   end Type

   interface MEF90_DefMechAT1SC_Type
      module procedure MEF90_DefMechAT1SC_Constructor
   end interface

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_DefMechAT1SC_Constructor"
!!!
!!!  MEF90_DefMechAT1SC_Constructor: the default constructor for a MEF90_DefMechAT1SC_Type
!!!
   Type(MEF90_DefMechAT1SC_Type) Function MEF90_DefMechAT1SC_Constructor(E0,E1)
      PetscReal,Intent(IN)                           :: E0,E1
      MEF90_DefMechAT1SC_Constructor%E0              = E0
      MEF90_DefMechAT1SC_Constructor%E1              = E1
      MEF90_DefMechAT1SC_Constructor%p0              = MEF90_DefMechAT1SC_Constructor%E1 / MEF90_DefMechAT1SC_Constructor%E0
      MEF90_DefMechAT1SC_Constructor%cw              = 2.0_Kr / 3.0_Kr
      MEF90_DefMechAT1SC_Constructor%aorder          = 2
      MEF90_DefMechAT1SC_Constructor%worder          = 1
      MEF90_DefMechAT1SC_Constructor%type            = 'MEF90_DefMechAT1SC'
   End Function MEF90_DefMechAT1SC_Constructor

#undef __FUNCT__
#define __FUNCT__ "aAT1SC"
!!!
!!!  aAT1SC: the stiffness interpolation, a(alpha) = (1-alpha) e_sc(alpha,p0)
!!!          with A = -5 delta (1-alpha) + 2 p0 - 3, s = sqrt(A^2 + 24 p0) and
!!!          e_sc = 4/(s + A)
!!!
   PetscReal function aAT1SC(self,alpha)
      Class(MEF90_DefMechAT1SC_Type),Intent(IN)         :: self
      PetscReal                                        :: alpha

      PetscReal                                        :: A,s
      PetscLogDouble                                   :: flops
      PetscErrorCode                                   :: ierr

      A = -5.0_Kr*(self%p0 - 1.0_Kr)*(1.0_Kr - alpha) + 2.0_Kr*self%p0 - 3.0_Kr
      s = sqrt(A*A + 24.0_Kr*self%p0)
      aAT1SC = (1.0_Kr - alpha)*4.0_Kr/(s + A)
      flops = 10.0
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End function aAT1SC

#undef __FUNCT__
#define __FUNCT__ "DaAT1SC"
!!!
!!!  DaAT1SC: derivative of aAT1SC with respect to alpha
!!!          with e_sc = 4/(s+A), e_sc' = -20 delta /(s (s+A)) and
!!!          Da = -e_sc + (1-alpha) e_sc'
!!!
   PetscReal function DaAT1SC(self,alpha)
      Class(MEF90_DefMechAT1SC_Type),Intent(IN)         :: self
      PetscReal                                        :: alpha

      PetscReal                                        :: A,s,sa,delta
      PetscLogDouble                                   :: flops
      PetscErrorCode                                   :: ierr

      delta = self%p0 - 1.0_Kr
      A = -5.0_Kr*delta*(1.0_Kr - alpha) + 2.0_Kr*self%p0 - 3.0_Kr
      s = sqrt(A*A + 24.0_Kr*self%p0)
      sa = s + A
      DaAT1SC = -4.0_Kr/sa - (1.0_Kr - alpha)*20.0_Kr*delta/(s*sa)
      flops = 12.0
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End function DaAT1SC

#undef __FUNCT__
#define __FUNCT__ "D2aAT1SC"
!!!
!!!  D2aAT1SC: second derivative of aAT1SC with respect to alpha
!!!           with e_sc'' = 100 delta^2 / s^3 and
!!!           D2a = -2 e_sc' + (1-alpha) e_sc''
!!!
   PetscReal function D2aAT1SC(self,alpha)
      Class(MEF90_DefMechAT1SC_Type),Intent(IN)         :: self
      PetscReal                                        :: alpha

      PetscReal                                        :: A,s,sa,delta

      delta = self%p0 - 1.0_Kr
      A = -5.0_Kr*delta*(1.0_Kr - alpha) + 2.0_Kr*self%p0 - 3.0_Kr
      s = sqrt(A*A + 24.0_Kr*self%p0)
      sa = s + A
      D2aAT1SC = 40.0_Kr*delta/(s*sa) + (1.0_Kr - alpha)*100.0_Kr*delta**2/s**3
   End function D2aAT1SC

#undef __FUNCT__
#define __FUNCT__ "wAT1SC"
!!!
!!!  wAT1SC: the "w" function, identical to the standard AT1 model, w(alpha) = alpha
!!!
   PetscReal function wAT1SC(self,alpha)
      Class(MEF90_DefMechAT1SC_Type),Intent(IN)         :: self
      PetscReal                                        :: alpha

      wAT1SC = alpha
   End function wAT1SC

#undef __FUNCT__
#define __FUNCT__ "DwAT1SC"
!!!
!!!  DwAT1SC: derivative of wAT1SC, Dw(alpha) = 1
!!!
   PetscReal function DwAT1SC(self,alpha)
      Class(MEF90_DefMechAT1SC_Type),Intent(IN)         :: self
      PetscReal                                        :: alpha

      DwAT1SC = 1.0_Kr
   End function DwAT1SC

#undef __FUNCT__
#define __FUNCT__ "D2wAT1SC"
!!!
!!!  D2wAT1SC: second derivative of wAT1SC, D2w(alpha) = 0
!!!
   PetscReal function D2wAT1SC(self,alpha)
      Class(MEF90_DefMechAT1SC_Type),Intent(IN)         :: self
      PetscReal                                        :: alpha

      D2wAT1SC = 0.0_Kr
   End function D2wAT1SC
End module m_MEF90_DefMechAT1SC
