#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module m_MEF90_DefMechAT1MT
#include "finclude/petscdef.h"
   Use m_MEF90
   Use m_MEF90_DefMechAT_class
   implicit none

!!!
!!!  MEF90_DefMechAT1MT_Type: Ambrosio-Tortorelli-type model whose stiffness
!!!  interpolation follows the Mori-Tanaka estimate for a two-phase medium
!!!
!!!  a(alpha) = (1-alpha) e_mt(alpha,p0) / p0, with
!!!     e_mt(c,p0) = (1 + 5(1-c)(p0-1)/(2 c p0 + 3 + 2(1-c))) / p0
!!!  where c is the volume fraction of the particle phase (c = alpha here),
!!!  p0 = E1/E0, E0 is the Young modulus of the matrix and E1 that of the particle.
!!!
!!!  The factor (1-alpha) makes a(1) = 0, i.e. the stiffness vanishes completely.
!!!
!!!  w(alpha), cw, aOrder and wOrder are taken from the standard AT1 model.
!!!

   Type, extends(MEF90_DefMechAT_Type)                 :: MEF90_DefMechAT1MT_Type
      PetscReal                                        :: E0
      PetscReal                                        :: E1
      PetscReal                                        :: p0
   Contains
      Procedure, pass(self)                            :: a   => aAT1MT
      Procedure, pass(self)                            :: Da  => DaAT1MT
      Procedure, pass(self)                            :: D2a => D2aAT1MT

      Procedure, pass(self)                            :: w   => wAT1MT
      Procedure, pass(self)                            :: Dw  => DwAT1MT
      Procedure, pass(self)                            :: D2w => D2wAT1MT
   end Type

   interface MEF90_DefMechAT1MT_Type
      module procedure MEF90_DefMechAT1MT_Constructor
   end interface

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_DefMechAT1MT_Constructor"
!!!
!!!  MEF90_DefMechAT1MT_Constructor: the default constructor for a MEF90_DefMechAT1MT_Type
!!!
   Type(MEF90_DefMechAT1MT_Type) Function MEF90_DefMechAT1MT_Constructor(E0,E1)
      PetscReal,Intent(IN)                           :: E0,E1
      MEF90_DefMechAT1MT_Constructor%E0              = E0
      MEF90_DefMechAT1MT_Constructor%E1              = E1
      MEF90_DefMechAT1MT_Constructor%p0              = MEF90_DefMechAT1MT_Constructor%E1 / MEF90_DefMechAT1MT_Constructor%E0
      MEF90_DefMechAT1MT_Constructor%cw              = 2.0_Kr / 3.0_Kr
      MEF90_DefMechAT1MT_Constructor%aorder          = 4
      MEF90_DefMechAT1MT_Constructor%worder          = 1
      MEF90_DefMechAT1MT_Constructor%type            = 'MEF90_DefMechAT1MT'
   End Function MEF90_DefMechAT1MT_Constructor

#undef __FUNCT__
#define __FUNCT__ "aAT1MT"
!!!
!!!  aAT1MT: the stiffness interpolation, a(alpha) = (1-alpha) e_mt(alpha,p0)/p0
!!!
   PetscReal function aAT1MT(self,alpha)
      Class(MEF90_DefMechAT1MT_Type),Intent(IN)         :: self
      PetscReal                                        :: alpha

      PetscReal                                        :: p0,den
      PetscLogDouble                                   :: flops
      PetscErrorCode                                   :: ierr

      p0  = self%p0
      den = 2.0_Kr*alpha*p0 + 3.0_Kr + 2.0_Kr*(1.0_Kr - alpha)
      aAT1MT = (1.0_Kr - alpha)*(1.0_Kr + 5.0_Kr*(1.0_Kr - alpha)*(p0 - 1.0_Kr)/den)/p0
      flops = 12.0
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End function aAT1MT

#undef __FUNCT__
#define __FUNCT__ "DaAT1MT"
!!!
!!!  DaAT1MT: derivative of aAT1MT with respect to alpha
!!!          with delta = p0-1, den = 5 + 2 delta alpha,
!!!          f   = (1-alpha) n and n = 5 + 5 delta - 3 delta alpha,
!!!          Da  = (f' den - 2 delta f) / (p0 den^2)
!!!
   PetscReal function DaAT1MT(self,alpha)
      Class(MEF90_DefMechAT1MT_Type),Intent(IN)         :: self
      PetscReal                                        :: alpha

      PetscReal                                        :: delta,den,f,fp
      PetscLogDouble                                   :: flops
      PetscErrorCode                                   :: ierr

      delta = self%p0 - 1.0_Kr
      den   = 5.0_Kr + 2.0_Kr*delta*alpha
      f     = 5.0_Kr + 5.0_Kr*delta - (5.0_Kr + 8.0_Kr*delta)*alpha + 3.0_Kr*delta*alpha**2
      fp    = -5.0_Kr - 8.0_Kr*delta + 6.0_Kr*delta*alpha
      DaAT1MT = (fp*den - 2.0_Kr*delta*f)/(self%p0*den**2)
      flops = 12.0
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End function DaAT1MT

#undef __FUNCT__
#define __FUNCT__ "D2aAT1MT"
!!!
!!!  D2aAT1MT: second derivative of aAT1MT with respect to alpha
!!!           D2a = (f'' den^2 - 4 delta f' den + 8 delta^2 f) / (p0 den^3)
!!!
   PetscReal function D2aAT1MT(self,alpha)
      Class(MEF90_DefMechAT1MT_Type),Intent(IN)         :: self
      PetscReal                                        :: alpha

      PetscReal                                        :: delta,den,f,fp,fpp

      delta = self%p0 - 1.0_Kr
      den   = 5.0_Kr + 2.0_Kr*delta*alpha
      f     = 5.0_Kr + 5.0_Kr*delta - (5.0_Kr + 8.0_Kr*delta)*alpha + 3.0_Kr*delta*alpha**2
      fp    = -5.0_Kr - 8.0_Kr*delta + 6.0_Kr*delta*alpha
      fpp   = 6.0_Kr*delta
      D2aAT1MT = (fpp*den**2 - 4.0_Kr*delta*fp*den + 8.0_Kr*delta**2*f)/(self%p0*den**3)
   End function D2aAT1MT

#undef __FUNCT__
#define __FUNCT__ "wAT1MT"
!!!
!!!  wAT1MT: the "w" function, identical to the standard AT1 model, w(alpha) = alpha
!!!
   PetscReal function wAT1MT(self,alpha)
      Class(MEF90_DefMechAT1MT_Type),Intent(IN)         :: self
      PetscReal                                        :: alpha

      wAT1MT = alpha
   End function wAT1MT

#undef __FUNCT__
#define __FUNCT__ "DwAT1MT"
!!!
!!!  DwAT1MT: derivative of wAT1MT, Dw(alpha) = 1
!!!
   PetscReal function DwAT1MT(self,alpha)
      Class(MEF90_DefMechAT1MT_Type),Intent(IN)         :: self
      PetscReal                                        :: alpha

      DwAT1MT = 1.0_Kr
   End function DwAT1MT

#undef __FUNCT__
#define __FUNCT__ "D2wAT1MT"
!!!
!!!  D2wAT1MT: second derivative of wAT1MT, D2w(alpha) = 0
!!!
   PetscReal function D2wAT1MT(self,alpha)
      Class(MEF90_DefMechAT1MT_Type),Intent(IN)         :: self
      PetscReal                                        :: alpha

      D2wAT1MT = 0.0_Kr
   End function D2wAT1MT
End module m_MEF90_DefMechAT1MT
